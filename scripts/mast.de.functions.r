########## MAST DE functions 

groupify <- function(f, vars) {
	df <- model.frame(f, data=vars)
	as.numeric(as.factor(apply(df, 1, function(a) paste(a, collapse='~'))))
}

subsample_cells <- function(covariates, data, pattern, n_cells) {
	# Short-circuit when no work is needed
	if (n_cells >= nrow(covariates)) {
		return (list(covariates=covariates, data=data))
	}

	# Get groups and their sizes
	groups <- if (length(pattern) > 0) {
		groupify(pattern[[1]], covariates)
	} else {
		rep(1, nrow(covariates))
	}
	Ng <- sort(table(groups))

	# What should be the group size cutoff?
	if (Ng[1] * length(Ng) > n_cells) {
		maxN <- floor(n_cells / length(Ng))
		extraI <- length(Ng) - (n_cells - maxN * length(Ng))
	} else {
		needed <- cumsum(Ng) + (length(Ng) - seq_along(Ng)) * Ng
		cut_i <- rev(which(needed <= n_cells))[1]
		left <- n_cells - needed[cut_i]
		maxNrest <- floor(left / (length(Ng) - cut_i))
		maxN <- Ng[cut_i] + maxNrest
		extraI <- length(Ng) - (left - maxNrest * (length(Ng) - cut_i))
	}

	# Subsample groups
	pattern_rest <- pattern[-1]
	covariates.sub.l <- list()
	data.sub.l <- list()
	for (i in seq_along(Ng)) {
		mask <- groups == names(Ng)[i]
		covi <- covariates[mask,,drop=F]
		datai <- data[,mask,drop=F]
		maxNi <- maxN + ifelse(i > extraI, 1, 0)
		if (Ng[i] <= maxNi) {
			# No sampling needed
		} else if (length(pattern_rest) > 0) {
			# Recurse on the next pattern
			res <- subsample_cells(covi, datai, pattern_rest, maxNi)
			covi <- res$covariates
			datai <- res$data
		} else {
			# Subsample this group
			I <- sample.int(nrow(covi), maxNi)
			covi <- covi[I,,drop=F]
			datai <- datai[,I,drop=F]
		}
		covariates.sub.l <- c(covariates.sub.l, list(covi))
		data.sub.l <- c(data.sub.l, list(datai))
	}
	
	# Re-combine subsampled groups
	covariates.sub <- do.call(rbind, covariates.sub.l)
	data.sub <- do.call(cbind, data.sub.l)
	
	list(covariates=covariates.sub, data=data.sub)
}

mast.de <- function(data, f, covariates, contrasts, ...) {
	# Model fitting and LRT
	fdata <- data.frame(a=rep(0,nrow(data)))[,c(),drop=F]
	sca <- MAST::FromMatrix(data, covariates, fdata)
	gc()
	zlm.res <- zlm(f, sca, force=T, ...)
	gc()
	res <- summary(zlm.res, doLRT=contrasts)$datatable
	
	# Get component information
	res.f <- res[res$component == 'logFC', .(primerid, contrast, coef)]
	res.d <- res[res$component == 'D', .(primerid, contrast, coef, `Pr(>Chisq)`)]
	res.c <- res[res$component == 'C', .(primerid, contrast, coef, `Pr(>Chisq)`)]
	res.h <- res[res$component == 'H', .(primerid, contrast, `Pr(>Chisq)`)]
	
	# Combine results
	res <- merge(res.d, res.c, by=c('primerid', 'contrast'), all=T, suffixes=c('D', 'C'))
	res <- Reduce(function(...) merge(..., by=c('primerid', 'contrast'), all=T), list(res, res.f, res.h))
	res <- data.frame(subset(res, !is.na(`Pr(>Chisq)`)), stringsAsFactors=F)
	
	return (res)
}

mast.de.batched <- function(data, f, covariates, contrasts, batchsize=50000, ...) { # if wish to run batch then set batchsize=500
	res <- list()
	while (nrow(data) %% batchsize == 1) {
		batchsize <- batchsize + 1
	}
	nbatches <- ceiling(nrow(data)/batchsize)
	for (i in seq_len(nbatches)) {
		first <- (i-1) * batchsize + 1
		last <- min(first+batchsize-1,nrow(data))
		cat(sprintf("###########----%s----Running batch %d/%d (genes %d - %d)\n", format(Sys.time(), "%b-%d-%X-%Y"), i, nbatches, first, last))
		flush.console()
		res.sub <- mast.de(data[first:last,,drop=F], f, covariates, contrasts, ...)
		res <- c(res, list(res.sub))
	}
	return (do.call(rbind, res))
}

run.celltype <- function(seur, runlabel, celltype, prefilter_only=F) {
	cat(sprintf("################----------- Working on %s/%s --------#################\n", runlabel, celltype))

	covariates = data.frame(
		n_genes = scale(seur@meta.data$n_genes),
		Layer = factor(seur@meta.data$Layer, levels=c('N', 'E', 'L')),
		Type = factor(seur@meta.data$Type, levels=c('Heal', 'NonI', 'Infl')),
		Donor = factor(seur@meta.data$PubID),
		Channel = factor(seur@meta.data$PubIDSample))
	cv <- covariates[seur@meta.data$anno2==celltype,,drop=F]
	
	if (any(table(cv$Type) < 10)) {
		cat(sprintf("Not enough cells per inflammation condition\n"))
		print(table(cv$Type))
		return()
	}
	
	# Get the data from Seurat
	Idents(seur) <- seur@meta.data$anno2
	data <- as.matrix(GetAssayData(subset(seur, ident=celltype), assay = "RNA", slot = "data"))
	
	# Keep only genes expressing the minimum fraction in the top-level subsample pattern
	cat("Filtering genes by min expression...\n")
	flush.console()
	minexpr_groups <- if (length(subsample.pattern) > 0) {
		groupify(subsample.pattern[[1]], cv)
	} else {
		rep(1, nrow(cv))
	}
	max_grp_expr_frac <- rep(0, nrow(data))
	for (g in unique(minexpr_groups)) {
		max_grp_expr_frac <- pmax(max_grp_expr_frac, rowMeans(data[,minexpr_groups==g,drop=F]>0))
	}
	keep_genes <- max_grp_expr_frac >= min_expression_frac
	cat(sprintf("Filtered out %d genes (%d left, %.1f%%)\n", sum(!keep_genes), sum(keep_genes), 100*mean(keep_genes)))
	flush.console()
	data <- data[keep_genes,,drop=F]
	
	# Subsample if we have too many cells
	if (ncol(data) > max.cells) {
		cat(sprintf("Subsampling %d cells to %d...\n", ncol(data), max.cells))
		flush.console()
		ret <- subsample_cells(cv, data, subsample.pattern, max.cells)
		data <- ret$data
		cv <- ret$covariates
	}

	# Prefilter with a simpler, but anti-conservative model
	cat("Running prefilter...\n")
	flush.console()
	de.res.pre <- mast.de.batched(data, test.formula.prefilter, cv, test.contrasts)
	prefilter_keep <- pmin(pmin(de.res.pre$Pr..Chisq.D, de.res.pre$Pr..Chisq.C), de.res.pre$Pr..Chisq.) < prefilter.P.thresh
	keep_genes <- rownames(data) %in% de.res.pre$primerid[prefilter_keep]
	cat(sprintf("Prefilter eliminated %d genes (%d left, %.1f%%)\n", sum(!keep_genes), sum(keep_genes), 100*mean(keep_genes)))
	flush.console()
	data <- data[keep_genes,,drop=F]
	gc()
	
	if (prefilter_only) {
		res.final <- de.res.pre
	} else {
		# Run the full model
		cat("Running full model fits...\n")
		flush.console()
		de.res <- mast.de.batched(data, test.formula, cv, test.contrasts, method='glmer', ebayes=F, strictConvergence=F, fitArgsD = list(nAGQ = 0))
		
		# Combine the results
		res.final <- merge(de.res, de.res.pre, by=c('primerid','contrast'), suffixes=c("",".pre"), all=TRUE)
	}
	
	# And dump the results
	if (!dir.exists(de.dir)) {
		dir.create(de.dir)
	}
	csvfile <- gsub("[?/]", "!", sprintf("%s.%s.csv", runlabel, celltype))
	write.csv(res.final, file.path(de.dir, csvfile))
}
