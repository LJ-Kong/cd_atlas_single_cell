# Use dirichlet-multinomial regression to find significant changes in cell frequencies during disease
# Fuctions were modified based on Smillie C.S. scripts at: 
# https://github.com/cssmillie/ulcerative_colitis/blob/master/analysis.r
# ---------------------------------------------------------------------------------------------------

library(lmerTest)
library(MASS)
library(egg)
library(DirichletReg)
library(Matrix)
library(data.table)
library(tidyverse)
library(tidyr)
library("RColorBrewer")
library(cowplot)
library(plyr)

dirichlet_regression = function(counts, covariates){

    # Calculate regression
    counts = as.data.frame(as.matrix(counts))
    counts$counts = DR_data(counts)
    data = cbind(counts, covariates)
    fit = DirichReg(counts ~ tissue + condition, data)
    
    # Get p-values
    u = summary(fit)
    pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
    v = names(pvals)
    pvals = matrix(pvals, ncol=length(u$varnames))
    rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
    colnames(pvals) = u$varnames
    fit$pvals = pvals
    
    coefs = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 1]
    v = names(coefs)
    coefs = matrix(coefs, ncol=length(u$varnames))
    rownames(coefs) = gsub('condition', '', v[1:nrow(pvals)])
    colnames(coefs) = u$varnames
    fit$coefs = coefs
    
    fit
}

dirichlet_regression_simple = function(counts, covariates){

    # Dirichlet multinomial regression to detect changes in cell frequencies
    # formula is not quoted, example: counts ~ condition
    # counts is a [samples x cell types] matrix
    # covariates holds additional data to use in the regression
    #
    # Example:
    # counts = do.call(cbind, tapply(seur@data.info$orig.ident, seur@ident, table))
    # covariates = data.frame(condition=gsub('[12].*', '', rownames(counts)))
    # res = dirichlet_regression(counts, covariates, counts ~ condition)
    
    # Calculate regression
    counts = as.data.frame(as.matrix(counts))
    counts$counts = DR_data(counts)
    data = cbind(counts, covariates)
    fit = DirichReg(counts ~ condition, data)
    
    # Get p-values
    u = summary(fit)
    pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
    v = names(pvals)
    pvals = matrix(pvals, ncol=length(u$varnames))
    rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
    colnames(pvals) = u$varnames
    fit$pvals = pvals
    
    coefs = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 1]
    v = names(coefs)
    coefs = matrix(coefs, ncol=length(u$varnames))
    rownames(coefs) = gsub('condition', '', v[1:nrow(pvals)])
    colnames(coefs) = u$varnames
    fit$coefs = coefs
    
    fit
}


cell_proportions_tests <- function(freq, notlayer, simple=F) {
	# Take only samples from a layer appropriate for the data
	freq.layer <- freq[rownames(freq)%in%df$PubIDSample[df$Layer!=notlayer], ]

	# For the dirichlet-multinomial regression, we need to know the disease state for each sample
	# We can get this from the metadata table as follows:
	sample2health <- data.frame(unique(data.frame(sample=df$PubIDSample, health=df$Type)), row.names=1)
	sample2tissue <- data.frame(unique(data.frame(sample=df$PubIDSample, tissue=df$Layer)), row.names=1)
	sample2tissue$tissue <- as.character(sample2tissue$tissue)
	sample2tissue$tissue <- ifelse(sample2tissue$tissue == "N", 1, 0)
	covariates <- data.frame(condition=factor(sample2health[rownames(freq.layer),1], levels=c('Heal', 'NonI', 'Infl')), 
								tissue=factor(sample2tissue[rownames(freq.layer),1]), 
								row.names=rownames(freq.layer))

	# Calculate significant changes using dirichlet multinomial regression
	# This returns a matrix of p-values for each cell type / disease state
	if(simple){
		reg <- dirichlet_regression_simple(counts=freq.layer, covariates=covariates)
	}else{
		reg <- dirichlet_regression(counts=freq.layer, covariates=covariates)
	}
	pvals <- reg$pvals
	colnames(pvals) <- colnames(freq.layer)
	coefs <- reg$coefs
	colnames(coefs) <- colnames(freq.layer)

	pct <- 100*sweep(freq.layer, 1, rowSums(freq.layer), FUN="/")

	qvals <- matrix(p.adjust(pvals, method="fdr"), nrow=nrow(pvals), ncol=ncol(pvals))
	rownames(qvals) <- rownames(pvals)
	colnames(qvals) <- colnames(pvals)
	
	return (list(pct=pct, coefs=coefs, pvals=pvals, qvals=qvals, cov=covariates))
}

# rcolorbrewer "set" palette
set.colors = c(brewer.pal(9, 'Set1'), brewer.pal(7, 'Set2'), brewer.pal(12, 'Set3')[c(3,9,8)], 'violetred4')
set.colors[6] = 'khaki2'
set.colors[8] = 'lightskyblue2'
set.colors = rep(set.colors, 10)

matrix_barplot = function(data, data.cell, group_by=NULL, pvals=NULL, xlab='', ylab='', value='mean', error='se', legend.show=T, legend.title='Groups', colors='Paired', pos='dodge', border=NA,
                          out=NULL, nrow=1.5, ncol=1.5, coord_flip=FALSE, sig_only=F, sig_cut=0.05, do.facet=F, pval_cols=c(Over='red',Under='blue'), lab_num_s=3.5, use_loc=F, no_heal=F){
    
    # Plot barplot of [M x N] matrix
    # x-axis = matrix columns (e.g. cell types)
    # y-axis = matrix values (e.g. frequencies)
    # fill = matrix rows (e.g. samples) or groups (e.g. conditions)
    
    # Arguments:
    # group.by = the group of each row
    # pvals = [G x N] matrix of p-values for each group and column
    # error = sd, se, or none
    
    # Groups (default = rows)
    if(is.null(group_by)){group_by = rownames(data)}
    if(nlevels(group_by) == 0){group_by = as.factor(group_by)}
    
    # Select significant comparisons
    if(sig_only){
		if(use_loc){
			j = pvals <= sig_cut
		}else{
			j = apply(pvals, 2, min) <= sig_cut
		}
		if(sum(j) == 0){return(NULL)}
		data = data[,j,drop=F]
		pvals = pvals[,j,drop=F]
    }
        
    # Construct input data
    names.data = colnames(data)
    data = data.frame(group=group_by, data)
    group_levels = levels(group_by)
    colnames(data)[2:ncol(data)] = names.data
    data = as.data.table(gather_(data, 'x', 'y', setdiff(colnames(data), 'group')))
    
    # Value function
    if(value == 'mean'){vf = mean} else if(value == 'median'){vf = median} else {stop()}
    
    # Error function
    se = function(x, na.rm=T){sd(x, na.rm=na.rm)/sqrt(length(x))}    
    if(error == 'sd'){ef = sd} else if(error == 'se'){ef = se} else {ef = function(x, ...){0}}
    
    # Estimate error bars
    data = data[,.(u=vf(y, na.rm=T), s=ef(y, na.rm=T)),.(group, x)]
    data$merg_group_x <- paste(data$group, data$x, sep=".") # add one extra colname to merge total number of cell
	
	# Summarize the total cells per group
	names.cell <- colnames(data.cell)
	data.cell <- data.frame(group=group_by, data.cell)
	group_levels <- levels(group_by)
	colnames(data.cell)[2:ncol(data.cell)] <- names.cell
	data.cell <- as.data.table(gather_(data.cell, 'x', 'y', setdiff(colnames(data.cell), 'group')))
	data.cell <- as.data.frame(data.cell) %>% group_by(x,group) %>% dplyr::summarize(n=sum(y))
	data.cell$merg_group_x <- paste(data.cell$group, data.cell$x, sep=".")
	
	# Add the cell number to data 
	data <- merge(data, data.cell, all.x=T)
	data$n_over <- data$u < max(data$u, na.rm=T) / 2
	
    # Add p-values 1
    if(!is.null(pvals)){
        pvals = as.data.frame(pvals) %>% rownames_to_column('group') %>% gather(x, pval, -group) %>% as.data.table()
		setkeyv(data, c('x', 'group'))
		setkeyv(pvals, c('x', 'group'))
		data = merge(data, pvals, all=T)
		data$lab1 = ifelse(data$pval <= .001, '**', ifelse(data$pval <= .05, '*', ''))
    }
	

    if(coord_flip){names.cell = rev(names.cell); group_levels=rev(group_levels)}
	# data$x = factor(data$x, levels=names) # original, modified to have a desired order of cell type
    data$x = factor(data$x, levels=as.character(cell_types_anno$anno2))    
    data$group = factor(data$group, levels=group_levels)
	grpidx = paste(data$x, data$group, sep='-')
	baseidx = paste(data$x, levels(data$group)[1], sep='-')
	data$ubase = data$u[match(baseidx, grpidx)]
    
    # Get colors
    if(length(colors) == 1){colors = set.colors[1:length(group_levels)]}
    
    # Plot data
    if(pos == 'stack'){
        p = ggplot(data) + geom_bar(aes(x=x, y=u, fill=group), colour=border, size=.25, stat='identity') + scale_x_discrete(labels = wrap_format(10))
	if(error %in% c('sd', 'se')){p = p + geom_errorbar(aes(x=x, ymin=u-s, ymax=u+s, fill=group), stat='identity', width=.25)}
    } else {
        pos = position_dodge(width=0.7)
        p = ggplot(data) + geom_bar(aes(x=x, y=u, fill=group), colour=border, size=.25, stat='identity', position=position_dodge(width=0.7), width=.7) + scale_x_discrete(labels = wrap_format(10))
	if(error %in% c('sd', 'se')){p = p + geom_errorbar(aes(x=x, ymin=u-s, ymax=u+s, fill=group), stat='identity', position=pos, width=.25)}
    }
	
	if(use_loc & no_heal){
		if(legend.show){
			p = p + theme_bw() +
				scale_fill_manual(values=c(NonI = "#628cd9", Infl = "#d9628c")) + xlab(xlab) + ylab(ylab)
		}else{
			 p = p + theme_bw() +
				scale_fill_manual(values=c(NonI = "#628cd9", Infl = "#d9628c"), name=legend.title) + xlab(xlab) + ylab(ylab) +
				theme(legend.position="none")
		}
	}else if(use_loc){
		if(legend.show){
			p = p + theme_bw() +
				scale_fill_manual(values=c(CO = "#628cd9", TI = "#d9628c")) + xlab(xlab) + ylab(ylab)
		}else{
			 p = p + theme_bw() +
				scale_fill_manual(values=c(CO = "#628cd9", TI = "#d9628c"), name=legend.title) + xlab(xlab) + ylab(ylab) +
				theme(legend.position="none")
		}
	}else{
		if(legend.show){
			p = p + theme_bw() +
				scale_fill_manual(values=c(Heal = "#95d962", NonI = "#628cd9", Infl = "#d9628c")) + xlab(xlab) + ylab(ylab)
		}else{
			 p = p + theme_bw() +
				scale_fill_manual(values=c(Heal = "#95d962", NonI = "#628cd9", Infl = "#d9628c"), name=legend.title) + xlab(xlab) + ylab(ylab) +
				theme(legend.position="none")
		}
	}
    # Facet wrap
    if(do.facet == TRUE){
        p = p + facet_grid(group ~ ., scales='free')
    }

    dy = max(data$u + data$s, na.rm=T)*.01
    if(coord_flip == FALSE){
        p = p + theme(axis.text.x = element_text(angle=55, vjust=1, hjust=1))
		
		if(!is.null(pvals)) {
			p = p + geom_text(aes(x=x, y=u+s+dy, label=lab1, group=group, color=ifelse(u>ubase, 'Over', 'Under')), hjust='center', vjust=0.5, size=5, angle=0, position=pos) +
				scale_color_manual(values=pval_cols) +
				geom_text(aes(x=x, y=ifelse(n_over, u+s+8*dy, u-s-3*dy), label=sprintf("%d", n), group=group, hjust=ifelse(n_over, 0, 1)),
				          vjust=0.4, size=lab_num_s, angle=90, position=pos, color="black")
		}
    } else {
        p = p + coord_flip()
		if(!is.null(pvals)){
			p = p + geom_text(aes(x=x, y=u+s+dy, label=lab1, group=group, color=ifelse(u>ubase, 'Over', 'Under')), hjust='center', vjust=1, size=5, angle=90, position=pos) +
				scale_color_manual(values=pval_cols) + geom_text(aes(x=x, y=0, group=group, label=n), vjust=1, size=lab_num_s, angle=90, position=pos)
		}
    }
    
    # Save plot
    if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}
    p
}
