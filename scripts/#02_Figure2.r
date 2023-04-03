library(Seurat)
library(ggplot2)
library(cowplot)
require("RColorBrewer")
require("viridis")
library(ComplexHeatmap)
library(circlize)
require("png")
library(pheatmap)

library(dplyr)
library(ggrepel)
library(reshape2)
library(lmerTest)
library(MASS)
library(egg)

source("./data/map_colors_lvls.r")

## DE results 
## to generate de.cmb.rds see #DE_analysis.r
de.cmb <- readRDS("./data/de.cmb.rds")

merged.pbmc <- readRDS(file="./data/co_ti_cmb.rds")

############################################################
# Prepare data for Figure 2
############################################################

# tier-1 CD/IBD genes
ibd <- read.table(./data/FineMapping_paper.txt", header=T)
dim(ibd) # 23  2
keys <- as.character(ibd[,2])
length(keys) # 23
setdiff(keys,rownames(merged.pbmc)) #"TNFRSF6B"
keys <- intersect(keys, rownames(merged.pbmc))
length(keys) # 22

# subset to IBD only data
result <- subset(merged.pbmc, features=keys)
dim(result) #22 720633

Idents(result) <- result@meta.data$anno2
gc()
#saveRDS(result, "./data/co_ti_cmb.new.22.IBD.tier1.rds")

# rm full merged.pbmc to save ram
rm(merged.pbmc);gc()

#########################
## IF readin earlier results
result <- readRDS("./data/co_ti_cmb.new.22.IBD.tier1.rds")
dim(result) #22 720633
keys <- rownames(result)
length(keys) #22

meta <- result@meta.data
meta$location <- "TI"
meta$location[meta$Site=="CO"] <- "CO"
dim(meta) #720633     13


# filter out genes that dont have 10% in a cell type
a <- DotPlot(result, features=keys)
tmp <- a$data
setdiff(keys,unique(tmp[tmp$pct.exp>3,]$features.plot))
# [1] "RTEL1"  "IL12B"

keys <- as.character(unique(tmp[tmp$pct.exp>3,]$features.plot))
# subset to these expressing IBD genes
result <- subset(result, features=keys)
dim(result) # 20 720633

cpm_val.sub <- as.matrix(GetAssayData(result, assay = "RNA", slot = "data"))
dim(cpm_val.sub) #20 720633

############################################################
# Figure 2. A
############################################################

# Averaged feature expression by
Idents(result) <- result@meta.data$anno2
full.anno <- unique(result@meta.data$anno2)

# separate for CO and TI
select_loc <- "TI"
hmat.ave.exp <- log1p(AverageExpression(result[,meta$location == select_loc], features = keys, verbose = FALSE)$RNA)
hmat.ave.exp.ti <- hmat.ave.exp

# prepare a tmp df which contains the mismatched cellty between TI and CO
# tmp is set 0
extra_cols <- setdiff(full.anno, colnames(hmat.ave.exp))
tmp <- hmat.ave.exp[,1:length(extra_cols)]
colnames(tmp) <- extra_cols
tmp[,] <- 0
hmat.ave.exp <- cbind(hmat.ave.exp,tmp)


## z-score and order
hmat.ave.exp.z <- scale(t(hmat.ave.exp), center=F)

## Gini coefficients
gini <- function(x) {
	sum(abs(outer(x, x, FUN="-"))) / ((length(x)-1) * 2 * sum(x))
}
gene_gini <- apply(hmat.ave.exp.z, 2, gini)

## return max fraction cell type name
max_cy <- function(x){
	names(x[x==max(x)])
}
max_name <- apply(hmat.ave.exp.z, 2, max_cy)

gini.res <- data.frame(gini=gene_gini, celltype=max_name)
rownames(gini.res) <- names(gene_gini)
gini.res <- gini.res[order(gini.res$gini, decreasing = T),]
#write.csv(gini.res, "./data/Gini.res_20.IBD.tier1_TI.csv")

O <- order(apply(hmat.ave.exp.z, 1, which.max))
hmat.ave.exp.z <- hmat.ave.exp.z[O,]
O <- order(apply(hmat.ave.exp.z, 2, which.max))
hmat.ave.exp.z <- hmat.ave.exp.z[,O]
O <- order(apply(hmat.ave.exp.z, 1, which.max))
hmat.ave.exp.z <- hmat.ave.exp.z[O,]
O <- order(apply(hmat.ave.exp.z, 2, which.max))
hmat.ave.exp.z <- hmat.ave.exp.z[,O]
order.ti <- rownames(hmat.ave.exp.z)
order.ti2 <- colnames(hmat.ave.exp.z)

hmat.ave.exp.z[hmat.ave.exp.z>4] <- 4
ti.hmat.ave.exp.z <- hmat.ave.exp.z

###################
select_loc <- "CO"
hmat.ave.exp <- log1p(AverageExpression(result[,meta$location == select_loc], features = keys, verbose = FALSE)$RNA)
dim(hmat.ave.exp) # 22 55
hmat.ave.exp.co <- hmat.ave.exp

# prepare a tmp df which contains the mismatched cellty between TI and CO
# tmp is set 0
extra_cols <- setdiff(full.anno, colnames(hmat.ave.exp))
tmp <- hmat.ave.exp[,1:length(extra_cols)]
colnames(tmp) <- extra_cols
tmp[,] <- 0
hmat.ave.exp <- cbind(hmat.ave.exp,tmp)

## z-score and order
hmat.ave.exp.z <- scale(t(hmat.ave.exp), center=F)

## Gini coefficients
gini <- function(x) {
	sum(abs(outer(x, x, FUN="-"))) / ((length(x)-1) * 2 * sum(x))
}
gene_gini <- apply(hmat.ave.exp.z, 2, gini)

## return max fraction cell type name
max_cy <- function(x){
	names(x[x==max(x)])
}
max_name <- apply(hmat.ave.exp.z, 2, max_cy)

gini.res <- data.frame(gini=gene_gini, celltype=max_name)
rownames(gini.res) <- names(gene_gini)
gini.res <- gini.res[order(gini.res$gini, decreasing = T),]
#write.csv(gini.res, "./data/Gini.res_20.IBD.tier1_CO.csv")


hmat.ave.exp.z <- hmat.ave.exp.z[order.ti,] # use the same order from TI
hmat.ave.exp.z <- hmat.ave.exp.z[,order.ti2] # use the same order from TI
hmat.ave.exp.z[hmat.ave.exp.z>4] <- 4
co.hmat.ave.exp.z <- hmat.ave.exp.z


threshold <- 1
keep.celltypes <- rowSums(ti.hmat.ave.exp.z > threshold) + rowSums(co.hmat.ave.exp.z > threshold) > 0

ha <- rowAnnotation(cy=cell_types_anno$manu_orig[match(rownames(ti.hmat.ave.exp.z[keep.celltypes,]),cell_types_anno$anno2)], col = list(cy=cell_colors_overall))
ht1 <- Heatmap(ti.hmat.ave.exp.z[keep.celltypes,], column_title = "TI", col=colorRampPalette(c("skyblue","white","tomato"))(100), 
		rect_gp = gpar(col = "#e8e8e8", lwd = 0.5), row_names_side ="left", row_names_max_width = unit(12, "cm"),
		column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 8), heatmap_legend_param=list(title=sprintf("z-score expr", "TI")),
		column_names_rot = 45, cluster_rows=F, cluster_columns=F, show_row_names=T, show_column_names=T, left_annotation=ha)

ht2 <- Heatmap(co.hmat.ave.exp.z[keep.celltypes,], column_title = "CO", col=colorRampPalette(c("skyblue","white","tomato"))(100), 
		rect_gp = gpar(col = "#e8e8e8", lwd = 0.5), row_names_side ="left", row_names_max_width = unit(12, "cm"),
		column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 8), heatmap_legend_param=list(title=sprintf("z-score expr", "CO")),
		column_names_rot = 45, cluster_rows=F, cluster_columns=F, show_row_names=T, show_column_names=T)

ht1+ht2


pdf("AverageExpression_IBD_fraction_of_cells_TI_and_CO_tier1x_9x9.new.pdf", 9, 9)
ht1+ht2
dev.off()


############################################################
# Figure 2. B
############################################################

scatter_ibd <- function(select_loc, cpm_val.other.ti, meta.ti, def1){
	cy.lvl <- cell_types_anno[,c(1,2)]
	colnames(cy.lvl) <- c("anno2", "manu_orig")
	dim(cy.lvl) #66  2

	def1$anno2 <- gsub("^\\w+\\.", "", def1$celltype)
	inf.ids <- paste(def1$primerid[def1$contrast=="TypeInfl"], def1$anno2[def1$contrast=="TypeInfl"], sep="-")# genes.anno2 pairs pass FDR in infl vs. heal
	nni.ids <- paste(def1$primerid[def1$contrast=="TypeNonI"], def1$anno2[def1$contrast=="TypeNonI"], sep="-") # genes.anno2 pairs pass FDR in noni vs. heal

	hmat <- as.matrix(cpm_val.other.ti)

	colnames(hmat) <- meta.ti$anno2
	dim(hmat) #22 430903

	ctype <- unique(colnames(hmat))
	hmat.pctg.heal <- matrix(data = NA, nrow = length(rownames(hmat)), ncol = length(ctype))
	rownames(hmat.pctg.heal) <- rownames(hmat)
	colnames(hmat.pctg.heal) <- ctype
	hmat.pctg.noni <- hmat.pctg.heal
	hmat.pctg.infl <- hmat.pctg.heal

	## seurat has ln(x+1) in normalizated data
	## calculate the non-zero cells percentage per cell type
	for(cy in 1:length(ctype)){
		hmat.pctg.heal[,cy] <- rowMeans(hmat[,(colnames(hmat)==ctype[cy]) & (meta.ti$Type == "Heal")]>0)
		hmat.pctg.noni[,cy] <- rowMeans(hmat[,(colnames(hmat)==ctype[cy]) & (meta.ti$Type == "NonI")]>0)
		hmat.pctg.infl[,cy] <- rowMeans(hmat[,(colnames(hmat)==ctype[cy]) & (meta.ti$Type == "Infl")]>0)
	}

	df <- reshape2::melt(hmat.pctg.heal)
	colnames(df) <- c("gene", "celltype", "heal.pct")
	df$noni.pct <- reshape2::melt(hmat.pctg.noni)[,3]
	df$infl.pct <- reshape2::melt(hmat.pctg.infl)[,3]
	df$label <- sprintf("%s-%s", as.character(df$gene), df$celltype)
	df <- merge(df, cy.lvl, by.x="celltype", by.y="anno2", sort=F)
	dim(df) #1320    7 #1210    7

	thrh <- 0.3 # threshold used as distance from diagonal line
	df$n_h <- abs(df$noni.pct-df$heal.pct)>thrh
	df$i_h <- abs(df$infl.pct-df$heal.pct)>thrh
	df$i_n <- abs(df$infl.pct-df$noni.pct)>thrh
	df$n_h.de <- df$label%in%nni.ids
	df$i_h.de <- df$label%in%inf.ids

	# only keep cell types that have at least 200 cells
	cy_counts <- dplyr::count(meta.ti, anno2,.drop=FALSE)
	cy_counts <- cy_counts[cy_counts$n>200,]
	df <- df[df$celltype%in%unique(cy_counts$anno2),]

	pp1 <- ggplot(df, aes(x=heal.pct, y=noni.pct)) + xlab("Cell fraction expressing in Healthy") + ylab("Cell fraction expressing in Noninflamed") +
			geom_point(aes(fill=as.character(manu_orig), color=ifelse(n_h.de, "DE", as.character(manu_orig)), stroke=ifelse(n_h.de, "DE", "nonDE")), shape=21, size=2, alpha=0.8) + 
			theme_bw() + labs(title=sprintf("Noninflamed vs. Healthy in %s", select_loc, select_loc)) +
			scale_fill_manual(values=c(cell_colors_overall)) +
			scale_color_manual(values=c(c(DE="red"), cell_colors_overall)) +
			scale_discrete_manual("stroke", values=c(DE=0.5, nonDE=0.1)) +
			geom_abline(intercept = thrh, slope = 1, color="grey", linetype="dashed") + 
			geom_abline(intercept = -thrh, slope = 1, color="grey", linetype="dashed") +
			geom_abline(intercept = 0, slope = 1) + 
			geom_text_repel(box.padding = 0.2, segment.size=0.1, size=3, aes(label=ifelse(n_h&(!n_h.de), as.character(gene), ""))) +
			geom_text_repel(segment.size=0.1, point.padding=0.1, size=3, aes(label=ifelse(n_h&n_h.de, as.character(gene), "")), color="red") +
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



	pp2 <- ggplot(df, aes(x=heal.pct, y=infl.pct)) + xlab("Cell fraction expressing in Healthy") + ylab("Cell fraction expressing in Inflamed") +
			geom_point(aes(fill=as.character(manu_orig), color=ifelse(i_h.de, "DE", as.character(manu_orig)), stroke=ifelse(i_h.de, "DE", "nonDE")), shape=21, size=2, alpha=0.8) + 
			theme_bw() + labs(title=sprintf("Inflamed vs. Healthy in %s", select_loc, select_loc)) +
			scale_fill_manual(values=c(cell_colors_overall)) +
			scale_color_manual(values=c(c(DE="red"), cell_colors_overall)) +
			scale_discrete_manual("stroke", values=c(DE=0.5, nonDE=0.1)) +
			geom_abline(intercept = thrh, slope = 1, color="grey", linetype="dashed") + 
			geom_abline(intercept = -thrh, slope = 1, color="grey", linetype="dashed") +
			geom_abline(intercept = 0, slope = 1) + 
			geom_text_repel(box.padding = 0.2, segment.size=0.1, size=3, aes(label=ifelse(i_h&(!i_h.de), as.character(gene), ""))) +
			geom_text_repel(segment.size=0.1, point.padding=0.1, size=3, aes(label=ifelse(i_h&i_h.de, as.character(gene), "")), color="red") + 
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

	pp <- plot_grid(pp2,pp1)
	
	return(pp)
}

####################
select_loc <- "TI"

cpm_val.other.ti <- cpm_val.sub[,meta$location == select_loc]
meta.ti <- meta[meta$location == select_loc,]

def1 <- de.cmb[de.cmb$primerid%in%keys & de.cmb$location==select_loc & de.cmb$FDR.D<0.05,]
def1 <- def1[complete.cases(def1),]
dim(def1) #14 21

pp1 <- scatter_ibd(select_loc, cpm_val.other.ti, meta.ti, def1)

pdf("Scatter_IBD_fraction_of_cells_TI_tier1_updatedx_18x5.4_thr0.3.pdf", 18, 5.4)
print(pp1)
dev.off()

####################
select_loc <- "CO"
cpm_val.other.ti <- cpm_val.sub[,meta$location == select_loc]
meta.ti <- meta[meta$location == select_loc,]

def1 <- de.cmb[de.cmb$primerid%in%keys & de.cmb$location==select_loc & de.cmb$FDR.D<0.05,]
def1 <- def1[complete.cases(def1),]
dim(def1) #215  21

pp2 <- scatter_ibd(select_loc, cpm_val.other.ti, meta.ti, def1)

pdf("Scatter_IBD_fraction_of_cells_CO_tier1_updatedx_18x5.4_thr0.3.pdf", 18, 5.4)
print(pp2)
dev.off()
