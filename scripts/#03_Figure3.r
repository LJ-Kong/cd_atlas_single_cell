library(ggplot2)
library(cowplot)
library(dplyr)
library(ggrepel)
library(reshape2)
library(data.table)
library(pheatmap)

source("./data/map_colors_lvls.r")
## to generate de.cmb.rds see #DE_analysis.r
de.cmb <- readRDS("./data/de.cmb.rds")


############################################################
# Figure 3. A-B
############################################################

de.cmb.sub <- de.cmb[de.cmb$FDR.D<0.05,] # subset to significant ones
dim(de.cmb.sub) # 126985     21

## -- Colon
df <- de.cmb.sub[grep("^CO.", de.cmb.sub$X),]
df1 <- df[df$contrast=="TypeInfl", ]
dfgg1c <- dplyr::count(df1, celltype, direct.D,.drop=FALSE)
dfgg1c$loc <- gsub("^(\\w+)\\..*$", "\\1", dfgg1c$celltype)
dfgg1c$direct.D <- factor(dfgg1c$direct.D, levels=c(FALSE, TRUE), labels=c("UP", "DOWN"))
dfgg1c$cy <- gsub("^\\w+\\.", "", dfgg1c$celltype)
dfgg1c$cy <- gsub("!", "/", dfgg1c$cy) # replace "S100A8!9.Monocyte"
dfgg1c$cy <- factor(dfgg1c$cy, levels=rev(as.character(cell_types_anno$anno2)))

df2 <- df[df$contrast=="TypeNonI", ]
dfgg2c <- dplyr::count(df2, celltype, direct.D,.drop=FALSE)
dfgg2c$loc <- gsub("^(\\w+)\\..*$", "\\1", dfgg2c$celltype)
dfgg2c$direct.D <- factor(dfgg2c$direct.D, levels=c(FALSE, TRUE), labels=c("UP", "DOWN"))
dfgg2c$cy <- gsub("^\\w+\\.", "", dfgg2c$celltype)
dfgg2c$cy <- gsub("!", "/", dfgg2c$cy) # replace "S100A8!9.Monocyte"
dfgg2c$cy <- factor(dfgg2c$cy, levels=rev(as.character(cell_types_anno$anno2)))

## -- TI
df <- de.cmb.sub[grep("^TI.", de.cmb.sub$X),]
df1 <- df[df$contrast=="TypeInfl", ]
dfgg1t <- dplyr::count(df1, celltype, direct.D,.drop=FALSE)
dfgg1t$loc <- gsub("^(\\w+)\\..*$", "\\1", dfgg1t$celltype)
dfgg1t$direct.D <- factor(dfgg1t$direct.D, levels=c(FALSE, TRUE), labels=c("UP", "DOWN"))
dfgg1t$cy <- gsub("^\\w+\\.", "", dfgg1t$celltype)
dfgg1t$cy <- gsub("!", "/", dfgg1t$cy) # replace "S100A8!9.Monocyte"
dfgg1t$cy <- factor(dfgg1t$cy, levels=rev(as.character(cell_types_anno$anno2)))


df2 <- df[df$contrast=="TypeNonI", ]
dfgg2t <- dplyr::count(df2, celltype, direct.D,.drop=FALSE)
dfgg2t$loc <- gsub("^(\\w+)\\..*$", "\\1", dfgg2t$celltype)
dfgg2t$direct.D <- factor(dfgg2t$direct.D, levels=c(FALSE, TRUE), labels=c("UP", "DOWN"))
dfgg2t$cy <- gsub("^\\w+\\.", "", dfgg2t$celltype)
dfgg2t$cy <- gsub("!", "/", dfgg2t$cy) # replace "S100A8!9.Monocyte"
dfgg2t$cy <- factor(dfgg2t$cy, levels=rev(as.character(cell_types_anno$anno2)))


legend_title <- "Direction"

library(grid)
# c(Immune.cell = "#ff7e0d", Stromal.cell = "#8c564b", Epithelial.cell = "#c5b1d5")
# strip.background = element_rect(fill=c("#c5b1d5", "#ff7e0d", "#8c564b"))
p1 <- ggplot(dfgg1c)+ labs(title="CO.Infl-CO.Heal", fill=legend_title) + theme_bw() + ylim(0, 2500) + facet_grid(loc ~ ., scales="free_y",space="free") +
		theme(axis.title.x = element_blank(), axis.title.y = element_blank(), strip.text.y = element_text(face="bold")) +
		geom_bar(aes(x=cy, y=n, fill=direct.D), stat="identity", width=0.6, position = position_stack(reverse = T))+ coord_flip()
g <- ggplot_gtable(ggplot_build(p1))
strip_r <- which(grepl('strip-r', g$layout$name))
fills <- c("#ecd4ff", "#ffbe85", "#cf9286")
k <- 1
for (i in strip_r) {
	j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
	g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
	k <- k+1
}
pdf("DEGs_sum_barplot_co.infl-heal_6x7.5_updated.pdf", 6, 7.5)
grid.draw(g)
dev.off()


p2 <- ggplot(dfgg1t)+ labs(title="TI.Infl-TI.Heal", fill=legend_title) + theme_bw() + ylim(0, 800) + facet_grid(loc ~ ., scales="free_y",space="free") +
		theme(axis.title.x = element_blank(), axis.title.y = element_blank(), strip.text.y = element_text(face="bold")) +
		geom_bar(aes(x=cy, y=n, fill=direct.D), stat="identity", width=0.6, position = position_stack(reverse = T))+ coord_flip()
g <- ggplot_gtable(ggplot_build(p2))
strip_r <- which(grepl('strip-r', g$layout$name))
fills <- c("#ecd4ff", "#ffbe85", "#cf9286")
k <- 1
for (i in strip_r) {
	j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
	g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
	k <- k+1
}
pdf("DEGs_sum_barplot_ti.infl-heal_6x7.5_updated.pdf", 6, 7.5)
grid.draw(g)
dev.off()


############################################################
# Figure 3. C
############################################################

# -- plot.ty: bw=between TI and Colon, or within one location
scatter.coef <- function (cmp, conts, thr=7, show.lab=TRUE, plot.ty="bw") {
	df <- de.cmb[grep(cmp, de.cmb$celltype),]
	
	if(plot.ty=="bw"){
		df <- df[df$contrast==conts,]
		df <- df[,c("primerid", "coefD", "FDR.D", "celltype", "location")]
		df <- df[complete.cases(df),]
		
		df <- df[abs(df$coefD)<10, ] # make a cut so there is no crazy point
		var_regex = '^MT|^RP' # remove MT, and RP genes based on HUGO gene names
		df <- df[grep(var_regex, df$primerid, invert=T),]
		
		df$cmb <- paste(df$primerid, df$celltype, sep=".")
		ids <- unique(df$cmb);length(ids) #92517

		dfx <- df[df$location=="TI",]
		dfy <- df[df$location=="CO",]
		dfz <- data.frame(TI=dfx[match(ids, dfx$cmb),]$coefD, CO=dfy[match(ids, dfy$cmb),]$coefD)
		rownames(dfz) <- ids
		dfz <- dfz[complete.cases(dfz),]
		
		cor.res <- cor.test(dfz$TI, dfz$CO, method = "spearman", alternative = "greater")
		cor.info <- sprintf("Spearman cor=%f\n(p.value=%e)", cor.res$estimate, cor.res$p.value)
		cat(cor.info)
		flush.console()
		
		if(conts=="TypeInfl"){
			tmp <- "Inflamed vs. Healthy"
			mtitle <- sprintf("%s\n(%s)", substr(cmp,2,nchar(cmp)-1), tmp)
		}else{
			tmp <- "Non-inflamed vs. Healthy"
			mtitle <- sprintf("%s\n(%s)", substr(cmp,2,nchar(cmp)-1), tmp)
		}
		thrh <- thr # threshold used as distance from diagonal line
		
		if(show.lab){
			p <- ggplot(dfz, aes(x=TI, y=CO)) + xlab("DE coefficient in TI") + ylab("DE coefficient in CO") +
					geom_point(color="black", size=0.55) + theme_bw() + labs(title=mtitle) +
					geom_abline(intercept = 0, slope = 1, color="tan2") + 
					geom_point(data=dfz[dfz$TI^2+dfz$CO^2>thrh^2 & dfz$TI*dfz$CO>0,], color='red', size=1.8) + # correlated
					geom_point(data=dfz[dfz$TI^2+dfz$CO^2>thrh^2 & dfz$TI*dfz$CO<0,], color='blue', size=1.8) + # anti-correlated
					geom_text_repel(box.padding = 0.5, aes(label=ifelse(TI^2+CO^2 > thrh^2, rownames(dfz), "")))
		}else{
			p <- ggplot(dfz, aes(x=TI, y=CO)) + xlab("DE coefficient in TI") + ylab("DE coefficient in CO") +
					geom_point(color="black", size=0.55) + theme_bw(base_size = 12) + labs(title=mtitle) +
					geom_abline(intercept = 0, slope = 1, color="tan2") +
					theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+
					#geom_text(size=6, x=min(dfz$TI)+2, y=max(dfz$CO), label=cor.info)
		}
	}else if(plot.ty=="within"){
		df <- df[df$location==conts,]
		df <- df[,c("primerid", "coefD", "FDR.D", "celltype", "contrast")]
		df <- df[complete.cases(df),]
		
		df <- df[abs(df$coefD)<10, ] # make a cut so there is no crazy point
		var_regex = '^MT|^RP' # remove MT, and RP genes based on HUGO gene names
		df <- df[grep(var_regex, df$primerid, invert=T),]
		
		df$cmb <- paste(df$primerid, df$celltype, sep=".")
		ids <- unique(df$cmb);length(ids) #92517

		dfx <- df[df$contrast=="TypeInfl",]
		dfy <- df[df$contrast=="TypeNonI",]
		dfz <- data.frame(Infl.vs.Heal=dfx[match(ids, dfx$cmb),]$coefD, NonI.vs.Heal=dfy[match(ids, dfy$cmb),]$coefD)
		rownames(dfz) <- ids
		dfz <- dfz[complete.cases(dfz),]
		
		cor.res <- cor.test(dfz$Infl.vs.Heal, dfz$NonI.vs.Heal, method = "spearman", alternative = "greater")
		cor.info <- sprintf("Spearman cor=%f\n(p.value=%e)", cor.res$estimate, cor.res$p.value)
		cat(cor.info)
		flush.console()
		
		if(conts=="TI"){
			mtitle <- sprintf("%s (TI)", substr(cmp,2,nchar(cmp)-1))
		}else{
			mtitle <- sprintf("%s (CO)", substr(cmp,2,nchar(cmp)-1))
		}
		thrh <- thr # threshold used as distance from diagonal line
		
		if(show.lab){
			p <- ggplot(dfz, aes(x=Infl.vs.Heal, y=NonI.vs.Heal)) + xlab("DE coefficient\nInflamed vs. Healthy") + ylab("DE coefficient\nNon-inflamed vs. Healthy") +
					geom_point(color="black", size=0.55) + theme_bw() + labs(title=mtitle) +
					geom_abline(intercept = 0, slope = 1, color="tan2") + 
					geom_point(data=dfz[dfz$Infl.vs.Heal^2+dfz$NonI.vs.Heal^2>thrh^2 & dfz$Infl.vs.Heal*dfz$NonI.vs.Heal>0,], color='red', size=1.8) + # correlated
					geom_point(data=dfz[dfz$Infl.vs.Heal^2+dfz$NonI.vs.Heal^2>thrh^2 & dfz$Infl.vs.Heal*dfz$NonI.vs.Heal<0,], color='blue', size=1.8) + # anti-correlated
					geom_text_repel(box.padding = 0.5, aes(label=ifelse(Infl.vs.Heal^2+NonI.vs.Heal^2 > thrh^2, rownames(dfz), "")))
		}else{
			p <- ggplot(dfz, aes(x=Infl.vs.Heal, y=NonI.vs.Heal)) + xlab("DE coefficient\nInflamed vs. Healthy") + ylab("DE coefficient\nNon-inflamed vs. Healthy") +
					geom_point(color="black", size=0.55) + theme_bw(base_size = 12) + labs(title=mtitle) +
					geom_abline(intercept = 0, slope = 1, color="tan2") +
					theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+ 
					#geom_text(size=6, x=min(dfz$Infl.vs.Heal)+2, y=max(dfz$NonI.vs.Heal), label=cor.info)
		}
	}
	return(p)
}

# coefD in infl-heal in TI vs. CO in three compartments
sc1 <- scatter.coef("^Epithelial.", "TypeInfl", thr=7, show.lab=F) 
# Spearman cor=0.250974 (p.value=0.000000e+00)

sc2 <- scatter.coef("^Stromal.", "TypeInfl", thr=4, show.lab=F) 
# Spearman cor=0.336552 (p.value=0.000000e+00)

sc3 <- scatter.coef("^Immune.", "TypeInfl", thr=4, show.lab=F) 
# Spearman cor=0.211611 (p.value=1.216956e-292)

pvt <- plot_grid(sc1,sc2,sc3,nrow=1)
pdf("fig_scatter.coefD.TI.vs.CO_TypeInfl.new_8.6x3.2.pdf", 8.6, 3.2)
pvt
dev.off()



############################################################
# Figure 3. D
############################################################

comm_cy <- intersect(unique(de.cmb[de.cmb$location=="TI",]$celltype), unique(de.cmb[de.cmb$location=="CO",]$celltype))
length(comm_cy) #48

de.sum.init <- rep(NA, length(comm_cy))
de.sum <- data.frame(ti.de=de.sum.init, co.de=de.sum.init, intersect.pos=de.sum.init, intersect.neg=de.sum.init, intersect.mix=de.sum.init,
			intersect=de.sum.init, union=de.sum.init, intersect_genes=rep("", length(comm_cy)), stringsAsFactors=F)
rownames(de.sum) <- comm_cy
for(i in 1:length(comm_cy)){
	ti.de.pos <- de.cmb[de.cmb$location=="TI" & de.cmb$celltype==comm_cy[i] & de.cmb$contrast=="TypeInfl" & de.cmb$FDR.D<0.05 & de.cmb$coefD>0,]$primerid
	ti.de.neg <- de.cmb[de.cmb$location=="TI" & de.cmb$celltype==comm_cy[i] & de.cmb$contrast=="TypeInfl" & de.cmb$FDR.D<0.05 & de.cmb$coefD<0,]$primerid
	co.de.pos <- de.cmb[de.cmb$location=="CO" & de.cmb$celltype==comm_cy[i] & de.cmb$contrast=="TypeInfl" & de.cmb$FDR.D<0.05 & de.cmb$coefD>0,]$primerid
	co.de.neg <- de.cmb[de.cmb$location=="CO" & de.cmb$celltype==comm_cy[i] & de.cmb$contrast=="TypeInfl" & de.cmb$FDR.D<0.05 & de.cmb$coefD<0,]$primerid
	de.sum$ti.de[i] <-length(ti.de.pos)+length(ti.de.neg)
	de.sum$co.de[i] <-length(co.de.pos)+length(co.de.neg)
	de.sum$intersect.pos[i] <-length(intersect(ti.de.pos, co.de.pos))
	de.sum$intersect.neg[i] <-length(intersect(ti.de.neg, co.de.neg))
	de.sum$intersect.mix[i] <-length(intersect(ti.de.neg, co.de.pos)) + length(intersect(ti.de.pos, co.de.neg))
	de.sum$intersect[i] <- de.sum$intersect.pos[i]+de.sum$intersect.neg[i]+de.sum$intersect.mix[i]
	de.sum$union[i] <-length(union(union(ti.de.pos, co.de.pos), union(ti.de.neg, co.de.neg)))
	de.sum$intersect_genes[i] <- paste(c(intersect(ti.de.pos, co.de.pos),intersect(ti.de.neg, co.de.neg)), collapse=", ")
}

total_genes <- 29756 # from co_ti_cmb.new.loom.rds
de.sum$expected <- de.sum$ti.de * de.sum$co.de / total_genes

de.sum$enrichment.consistent <- (de.sum$intersect.pos + de.sum$intersect.neg) / (de.sum$expected / 2)
de.sum$p_consistent <- ppois(de.sum$intersect.pos + de.sum$intersect.neg, de.sum$co.de * de.sum$ti.de / 2 / total_genes, lower.tail = FALSE)
de.sum$q_consistent <- p.adjust(de.sum$p_consistent, method="fdr")
dim(de.sum) #48 12

de.sum$cy <- gsub("^\\w+\\.", "", rownames(de.sum))
de.sum$compart <- gsub("^(\\w+)\\..*$", "\\1", rownames(de.sum))
#write.csv(de.sum, "Consistency_score.csv")

compart_colors <- c("Immune" = "#ff7e0d", "Stromal" = "#8c564b", "Epithelial" = "#c5b1d5")
p <- ggplot(de.sum, aes(x=reorder(cy,enrichment.consistent), y=enrichment.consistent, fill=compart)) + 
		scale_fill_manual(values=compart_colors) + geom_bar(stat="identity", width=0.8) +
		theme_classic() + coord_flip() + xlab(NULL) + ylab("Consistency score") #"#eeb758"

pdf("Consistency_score_barplot_18x4.pdf", 6.5, 9)
print(p)
dev.off()



############################################################
# Figure 3. E
############################################################

# coefD in infl-heal vs noni-heal in TI in three compartments
sc1 <- scatter.coef("^Epithelial.", "TI", thr=6, show.lab=F, plot.ty="within") 
# Spearman cor=0.669296 (p.value=0.000000e+00)

sc2 <- scatter.coef("^Stromal.", "TI", thr=6, show.lab=F, plot.ty="within") 
# Spearman cor=0.864358 (p.value=0.000000e+00)

sc3 <- scatter.coef("^Immune.", "TI", thr=6, show.lab=F, plot.ty="within") 
# Spearman cor=0.786330 (p.value=0.000000e+00)

pvt <- plot_grid(sc1,sc3,sc2,nrow=1)
pdf("fig_scatter.coefD.TI.within_8.6x3.2.pdf", 8.6, 3.2)
pvt
dev.off()

# coefD in infl-heal vs noni-heal in Colon in three compartments
sc1 <- scatter.coef("^Epithelial.", "CO", thr=8, show.lab=F, plot.ty="within") 
# Spearman cor=0.780434 (p.value=0.000000e+00)

sc2 <- scatter.coef("^Stromal.", "CO", thr=6, show.lab=F, plot.ty="within") 
# Spearman cor=0.721883 (p.value=0.000000e+00)

sc3 <- scatter.coef("^Immune.", "CO", thr=6, show.lab=F, plot.ty="within") 
# Spearman cor=0.823431 (p.value=0.000000e+00)

pvt <- plot_grid(sc1,sc3,sc2,nrow=1)
pdf("fig_scatter.coefD.CO.within_8.6x3.2.pdf", 8.6, 3.2)
pvt
dev.off()



############################################################
# Figure 3. F
############################################################

library(ggplot2)
library(cowplot)
library("RColorBrewer")
library("viridis")
library(dplyr)
library(ggrepel)
library(reshape2)
library(MASS)
library(egg)
library(fgsea)
library(pheatmap)
library(Seurat)


## Download pathways from: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C2
pty <- "./data/c2.cp.kegg.v7.0.symbols.gmt"
pathways.hallmark <- gmtPathways(pty) # files from mysigdb 
pathways.hallmark %>% 
   head() %>% 
   lapply(head)

source("./data/map_colors_lvls.r")

## DE results
de.cmb <- readRDS("./data/de.cmb.rds")

loc <- c("TI", "TI", "CO", "CO")
out_col <- c("TypeInfl", "TypeNonI", "TypeInfl", "TypeNonI")
cmb <- c()
for(j in 1:4){
	
	## different locations?
	de.tests <- de.cmb[de.cmb$location==loc[j], ]

	# go through all celltypes 
	celltys <- unique(de.tests$celltype)
	celltys <- celltys[celltys!="Low_QC"]

	tmp.fgseaRes <- c()
	for(i in 1:length(celltys)){
		cellty <- celltys[i]
		cat(sprintf("Working on: %s\n", cellty))
		flush.console()
		
		#ranklist <- de.tests$coefD[de.tests$celltype==cellty & de.tests$contrast==out_col[j]] # use fold change values
		ranklist <- de.tests$coefD.pre[de.tests$celltype==cellty & de.tests$contrast==out_col[j]] # use pre.coefD!!
		names(ranklist) <- de.tests$primerid[de.tests$celltype==cellty & de.tests$contrast==out_col[j]]
		ranklist <- ranklist[!is.na(ranklist)]
		
		if(length(ranklist)>100){
			fgseaRes <- fgsea(pathways = pathways.hallmark, 
								  stats = sort(ranklist),
								  minSize=3,
								  maxSize=500,
								  nperm=100000)
			
			fgseaRes$cellty <- cellty
			tmp.fgseaRes <- rbind(tmp.fgseaRes, fgseaRes)
		}
	}
	
	tmp.fgseaRes$location <- loc[j]
	tmp.fgseaRes$contrast <- out_col[j]
	if(j==1){
		cmb <- tmp.fgseaRes
	}else{
		cmb <- rbind(cmb, tmp.fgseaRes)
	}
}

dim(cmb) #35784    11
# saveRDS(cmb, "./data/cmb.fgseaRes.coefD.pre.rds")

######################################

cmb$cpart <- substr(cmb$cellty, 1, 3)
cmb$cy <- gsub("^\\w+\\.", "", cmb$cellty)

## pwy.collect returns 2 matrice: NES and padj (at least one significant), 
## and a data.frame contains total number of significant pwy per each direction
pwy.collect <- function(fRes){
	tst <- dcast(data=fRes, formula = pathway~cellty, fun.aggregate=sum, value.var="NES")
	tst <- as.data.frame(tst) # class(tst) -> "data.table" "data.frame"
	annotation_col = data.frame(CellType=substr(colnames(tst)[-1], 1, 3))
	colnames(tst) <- gsub("^\\w+\\.", "", colnames(tst))
	rownames(annotation_col) <- colnames(tst)[-1]
	rownames(tst) <- substr(tst$pathway, 6, nchar(tst$pathway))
	tst <- tst[,-1]
	

	tst2 <- dcast(data=fRes, formula = pathway~cellty, fun.aggregate=sum, value.var="padj")
	tst2 <- as.data.frame(tst2)
	colnames(tst2) <- gsub("^\\w+\\.", "", colnames(tst2))
	rownames(tst2) <- substr(tst2$pathway, 6, nchar(tst2$pathway))
	tst2 <- tst2[,-1]
	tst2[tst2==0] <- 1 # change all 0 padj to 1

	sign_paths <- dplyr::count(fRes[fRes$padj<0.05,],pathway)
	sign_paths$pathway <- substr(sign_paths$pathway, 6, nchar(sign_paths$pathway))
	sign_paths <- sign_paths[order(sign_paths$n, decreasing=T),]
	
	sign_tst <- tst[sign_paths$pathway,]
	sign_tst <- sign(sign_tst)
	
	sign_tst2 <- tst2[sign_paths$pathway,]
	sign_tst2[sign_tst2<0.05] <- -1
	sign_tst2[sign_tst2>=0.05] <- 0
	sign_tst2[sign_tst2==-1] <- 1
	
	sign_cmb <- sign_tst*sign_tst2 # collect NES directions and padj together
	sign_cmb_stat <- data.frame(pos=apply(sign_cmb, 1, function(x) sum(x>0)), neg=apply(sign_cmb, 1, function(x) sum(x<0)))
	rownames(sign_cmb_stat) <- rownames(sign_cmb)
	
	return(list(nes.mat=tst, padj.mat=tst2, collect.stat=sign_cmb_stat))
}

## Epi
fgseaRes <- cmb[cmb$contrast=="TypeInfl"&cmb$location=="TI"&cmb$cpart=="Epi",] #TI, TypeInfl, Epi
res1 <- pwy.collect(fgseaRes)

fgseaRes <- cmb[cmb$contrast=="TypeInfl"&cmb$location=="CO"&cmb$cpart=="Epi",] #CO, TypeInfl, Epi
res2 <- pwy.collect(fgseaRes)

## Str
fgseaRes <- cmb[cmb$contrast=="TypeInfl"&cmb$location=="TI"&cmb$cpart=="Str",] #TI, TypeInfl, Str
res3 <- pwy.collect(fgseaRes)

fgseaRes <- cmb[cmb$contrast=="TypeInfl"&cmb$location=="CO"&cmb$cpart=="Str",] #CO, TypeInfl, Str
res4 <- pwy.collect(fgseaRes)

## Imm
fgseaRes <- cmb[cmb$contrast=="TypeInfl"&cmb$location=="TI"&cmb$cpart=="Imm",] #TI, TypeInfl, Imm
res5 <- pwy.collect(fgseaRes)

fgseaRes <- cmb[cmb$contrast=="TypeInfl"&cmb$location=="CO"&cmb$cpart=="Imm",] #CO, TypeInfl, Imm
res6 <- pwy.collect(fgseaRes)


###############
## Epi, barplot
df1 <- res1$collect.stat
df1$pwy <- rownames(df1)
df1$loc <- "TI" 
df2 <- res2$collect.stat
df2$pwy <- rownames(df2)
df2$loc <- "CO" 
df <- rbind(df1,df2)
df <- melt(df)

df$pct[df$loc=="TI"] <- df$value[df$loc=="TI"]/ncol(res1$nes.mat)
df$pct[df$loc=="CO"] <- df$value[df$loc=="CO"]/ncol(res2$nes.mat)
df$loc <- factor(df$loc, levels=c("TI","CO"))
df$cpart <- "Epi"
df.epi <- df


###############
## Str, barplot
df1 <- res3$collect.stat
df1$pwy <- rownames(df1)
df1$loc <- "TI" 
df2 <- res4$collect.stat
df2$pwy <- rownames(df2)
df2$loc <- "CO" 
df <- rbind(df1,df2)
df <- melt(df)

df$pct[df$loc=="TI"] <- df$value[df$loc=="TI"]/ncol(res3$nes.mat)
df$pct[df$loc=="CO"] <- df$value[df$loc=="CO"]/ncol(res4$nes.mat)
df$loc <- factor(df$loc, levels=c("TI","CO"))
df$cpart <- "Str"
df.str <- df


###############
## Imm, barplot
df1 <- res5$collect.stat
df1$pwy <- rownames(df1)
df1$loc <- "TI" 
df2 <- res6$collect.stat
df2$pwy <- rownames(df2)
df2$loc <- "CO" 
df <- rbind(df1,df2)
df <- melt(df)

df$pct[df$loc=="TI"] <- df$value[df$loc=="TI"]/ncol(res5$nes.mat)
df$pct[df$loc=="CO"] <- df$value[df$loc=="CO"]/ncol(res6$nes.mat)
df$loc <- factor(df$loc, levels=c("TI","CO"))
df$cpart <- "Imm"
df.imm <- df

#############
df.cmb <- rbind(df.epi, df.str, df.imm)

pct.mat <- dcast(data=df.cmb, formula = pwy~cpart+loc, fun.aggregate=sum, value.var="pct")
dim(pct.mat) #133   7
pct.mat <- pct.mat[rowMeans(pct.mat[,-1])>0.1,] # filter out some low percentage pwy
dim(pct.mat) #35  7
pct.mat.f <- dcast(data=df.cmb, formula = pwy~variable+cpart+loc, fun.aggregate=sum, value.var="pct") # add direction back
dim(pct.mat.f) # 133  13
pct.mat.f <- pct.mat.f[pct.mat.f$pwy%in%pct.mat$pwy,] # keep the useful pwy
dim(pct.mat.f) # 35 13

pct.mat.f <- pct.mat.f[!pct.mat.f$pwy%in%c("RIBOSOME"),] # manually remove these uninteresting ones
dim(pct.mat.f) # 34 13

hc <- hclust(dist(pct.mat.f[,-1])) # order the pct
order_pwy <- pct.mat.f$pwy[hc$order]

df.pct.mat.f <- melt(pct.mat.f)
df.pct.mat.f$direction <- factor(substr(df.pct.mat.f$variable, 1, 3), levels=c("neg", "nil", "pos"))
df.pct.mat.f$cpart <- substr(df.pct.mat.f$variable, 5, 7)
df.pct.mat.f$loc <- substr(df.pct.mat.f$variable, 9, 10)
df.pct.mat.f$loc <- factor(df.pct.mat.f$loc, levels=c("TI","CO"))
df.pct.mat.f$pwy <- factor(df.pct.mat.f$pwy, levels=order_pwy)
df.pct.mat.f$cpart[df.pct.mat.f$cpart=="Epi"] <- "Epithelial"
df.pct.mat.f$cpart[df.pct.mat.f$cpart=="Str"] <- "Stromal"
df.pct.mat.f$cpart[df.pct.mat.f$cpart=="Imm"] <- "Immune"


pdf("fig03_pathways_barplot_sum.ti_15x6.pdf", 15, 6)
ggplot(data=df.pct.mat.f) + 
		geom_bar(aes(x=pwy, y=ifelse(direction=="neg",-value, value), fill=direction), stat="identity",
		width=0.6, position = position_stack(reverse = T))+ 
		geom_hline(yintercept=0) +
		facet_grid(~cpart+loc) + coord_flip() + 
		scale_fill_manual(name="Enrichment direction", labels=c(pos="Positive", neg="Negative"), values=c(pos="tomato3", neg="turquoise4")) +
		theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) +
		ylab("Fraction of cell types") + xlab("")
dev.off()
