############################################################
# Prepare data
############################################################

library(Seurat)
library(ggplot2)
library(cowplot)
require("RColorBrewer")
require("viridis")
library(dplyr)
library(SeuratWrappers)
library(harmony)
library(RColorBrewer)
library(circlize)

result <- readRDS(file="./data/co_ti_cmb.rds")
meta <- result@meta.data
meta$location <- "TI"
meta$location[meta$Site=="CO"] <- "CO"
meta$loc_anno <- paste(meta$location, meta$anno_overall, sep="_")
meta$anno2 <- factor(meta$anno2)
result@meta.data <- meta

## subset to TI stromal cells
Idents(result) <- result@meta.data$loc_anno
result <- subset(result, idents="TI_Stromal cells")
gc()

pbmc <- NormalizeData(result, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- RunFastMNN(object.list = SplitObject(pbmc, split.by = "PubID"))
pbmc <- RunUMAP(pbmc, reduction="mnn", n.neighbors=30L, dims=1:20, seed.use=4)
pbmc <- FindNeighbors(pbmc, reduction="mnn", dims=1:20, k.param=20, force.recalc=T)
pbmc <- FindClusters(pbmc, resolution=2)
#Maximum modularity in 10 random starts: 0.8205
#Number of communities: 31

# saveRDS(pbmc, "./data/TI.STR.seurat.rds")

## further subset to myofiblasts
pbmc <- subset(pbmc, idents=c(6,9,19))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- RunFastMNN(object.list = SplitObject(pbmc, split.by = "PubID"))
pbmc <- RunUMAP(pbmc, reduction="mnn", n.neighbors=30L, dims=1:20, seed.use=4)
pbmc <- FindNeighbors(pbmc, reduction="mnn", dims=1:20, k.param=20, force.recalc=T)
pbmc <- FindClusters(pbmc, resolution=0.8)
#Maximum modularity in 10 random starts: 0.7465
#Number of communities: 9

# saveRDS(pbmc, "./data/TI.MYO.seurat.rds")


############################################################
# Figure 5. A-D
############################################################

library(Seurat)
library(ggplot2)
library(cowplot)
require("RColorBrewer")
require("viridis")
library(dplyr)
library(monocle3)
library(Matrix)
library(SeuratWrappers)
library(pheatmap)
library(magrittr)
library(ggrepel)
set.seed(1234)

myo <- readRDS("./data/TI.MYO.seurat.rds")
meta <- myo@meta.data
myo.cds <- as.cell_data_set(myo)

## Calculate size factors using built-in function in monocle3
myo.cds <- estimate_size_factors(myo.cds)

## Add gene names into CDS
myo.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(myo.cds)

myo.cds <- cluster_cells(cds = myo.cds, reduction_method = "UMAP", k=50, cluster_method="louvain", random_seed=123)
myo.cds <- learn_graph(myo.cds, close_loop=F, use_partition=TRUE)

meta.mat <- colData(myo.cds)
umap_coords <- reducedDims(myo.cds)[["UMAP"]]
root_cells <- rownames(meta.mat)[which.min(umap_coords[,1])]
length(root_cells) # 1

# order cells
myo.cds <- order_cells(myo.cds, reduction_method="UMAP", root_cells=root_cells)
myo.cds@principal_graph_aux$UMAP$pseudotime <- max(myo.cds@principal_graph_aux$UMAP$pseudotime) - myo.cds@principal_graph_aux$UMAP$pseudotime

#############################
# Fig5-C
p1 <- plot_cells(myo.cds, color_cells_by="pseudotime", show_trajectory_graph=F)
pdf("#TI.MYO_monocle3_tree.myo.cds.louvain_4x4_from_c1.pdf", 5, 4)
print(p1)
dev.off()

df <- as.data.frame(reducedDims(myo.cds)[["UMAP"]])
df$monocle3_clusters <- as.character(myo.cds@clusters$UMAP$clusters)
df$celltype[df$monocle3_clusters%in%c("1", "2", "3", "4", "6", "7", "8", "9", "10", "13", "14", "15", "16", "17")] <- "Myofibroblast GREM1+ GREM2+"
df$celltype[df$monocle3_clusters%in%c("12", "18", "11", "5")] <- "Myofibroblast HHIP+ NPNT+"

#############################
# Fig5-B
p1 <- ggplot(data=df) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=celltype), size=1)+ theme_bw() +
		scale_color_manual(values=c("Myofibroblast GREM1+ GREM2+"="#72b3d3", "Myofibroblast HHIP+ NPNT+"="#8c5f0b")) + 
		theme(plot.title=element_text(size=20), panel.grid = element_blank(), axis.text=element_blank(), 
		axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks=element_blank()) +
		guides(color = guide_legend(override.aes = list(size=5)))

pdf("#TI.MYO_monocle3_celltype_UMAP_6.5x4.pdf", 6.5, 4)
print(p1)
dev.off()

# saveRDS(myo.cds, "./data/#TI.MYO.myo.cds.rds")

#############################
# Fig5-A
myo <- readRDS("./data/TI.MYO.seurat.rds")
myo.cds.x <- readRDS("./data/#TI.MYO.myo.cds.rds")

keys <- c('ACTA2', 'ACTG2', 'MYH11', 'HHIP', 'NPNT', 'GREM1', 'GREM2', 'COL18A1', 'COL23A1')
cpm_val.sub <- as.matrix(GetAssayData(myo, assay = "RNA", slot = "data")[keys,]) # original seurat normalized data
cpm_val.sub <- cpm_val.sub[,colnames(myo.cds.x)]
df <- as.data.frame(reducedDims(myo.cds.x)[["UMAP"]])

plotlist <- list()
for(i in 1:length(keys)){
	df$Expr <- cpm_val.sub[i,]
	dfOrd <- df[order(df$Expr),]
	ggp <- ggplot(data=dfOrd) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=Expr), size=0.1) + scale_colour_viridis_c()
	ggp <- ggp + labs(title=keys[i])
	ggp <- ggp + theme_cowplot() + theme(plot.title=element_text(size = 15, face = "bold"), panel.grid = element_blank(), axis.text=element_blank(), 
				axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks=element_blank())
	plotlist <- c(plotlist, list(ggp))
}
p3 <- plot_grid(plotlist=plotlist, nrow=3, ncol=3)
pdf("#TI.MYO_9markers.seu_8x7.pdf", 8, 7)
print(p3)
dev.off()

#############################
# Fig5-D

myo.cds.x@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(myo.cds.x)
res <- graph_test(myo.cds.x, neighbor_graph="principal_graph", cores=5)

# use DE results for removing some low-expressed genes
de.cmb <- readRDS("/mnt/e/Shared/lk/rerun/de-results.2022/de.cmb.rds")

kep_g <- unique(de.cmb$primerid)
kep_g <- grep('^MT|^RP', kep_g, invert=T, value=TRUE) #remove MT or RP genes

q <- res[kep_g,]$q_value
pst <- pseudotime(myo.cds.x)
counts <- exprs(myo.cds.x)[,order(pst)]
lib.sizes <- colSums(counts)
tp10k <- sweep(10000*counts[kep_g,], 2, lib.sizes, FUN="/")
pst <- sort(pst)
pst <- pst + seq_len(length(pst)) * 1e-4

keep <- q < 0.1
tp10k <- tp10k[keep,]
q <- q[keep]

smoothness <- 1.5
t_max <- apply(tp10k, 1, function(x){ pst[which.max(smooth.spline(pst, x,df=3,spar=smoothness)$y)] })

mlq <- -log10(q + 1e-10) - 1
mlq <- mlq + runif(n=length(mlq))*1e-11

tent_slope <- 1000
accept_offset <- 1.0
t_max_d <- outer(t_max, t_max, FUN="-")
mlq_d <- outer(mlq, mlq, FUN="-")

elim <- mlq_d < -abs(t_max_d) * tent_slope - accept_offset
good <- rowSums(elim) == 0

tests <- c("CHMP1A", "TBX3", "RNF168", "GREM1", "GREM2", "HHIP", "NPNT", "COL18A1", "COL23A1")

good2 <- good | names(good)%in%tests
pt.matrix <- tp10k[names(good2[good2])[order(t_max[good2])],]
expr.frac <- rowMeans(pt.matrix>0)
select_n <- 200
n_pick <- select_n
keepg <- union(tests, names(expr.frac)[order(expr.frac, decreasing=T)[1:n_pick]])
while (length(keepg) > select_n) {
	n_pick <- select_n - (length(keepg) - select_n)
	keepg <- union(tests, names(expr.frac)[order(expr.frac, decreasing=T)[1:n_pick]])
}
length(keepg)
pt.matrix <- pt.matrix[rownames(pt.matrix) %in% keepg,]

dim(pt.matrix) #[1]   37 9364

pred_t <- seq(0, max(pst), length.out=100)
pt.matrix <- t(apply(pt.matrix,1,function(x){predict(smooth.spline(pst, x, df=3, spar=smoothness)$fit, pred_t)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){scale(x, center=T)}))
rownames(pt.matrix) <- ifelse(rownames(pt.matrix)%in%tests, rownames(pt.matrix), "")

saveRDS(pt.matrix, "./data/pt.matrix.new.rds")

#################
library(ComplexHeatmap)
library(ggplot2)
library(cowplot)
library("RColorBrewer")

pt.matrix <- readRDS("./data/pt.matrix.new.rds")
tests <- c("CHMP1A", "TBX3", "RNF168", "GREM1", "GREM2", "HHIP", "NPNT", "COL18A1", "COL23A1")

subset.re <- match(tests,rownames(pt.matrix))

pdf("#TI.STR_monocle3_heatmap2.pdf", 7, 8)
Heatmap(pt.matrix, name="Z-scored\nexpression", col=rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(100)), column_names_gp = gpar(fontsize = 10),
cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=F, column_title="Pseudotime")+
rowAnnotation(link = anno_mark(at = subset.re, labels_gp = gpar(fontsize = 18), labels=tests))
dev.off()

############################################################
# Figure 5. E
############################################################

library(ggplot2)
library(cowplot)
require(RColorBrewer)
require(viridis)
library(reshape2)
library(ggrepel)

vp <- read.table("./data/ScreenData_VP.txt", sep="\t", row.names=1, header=T)
vp <- log2(as.matrix(vp))

P <- matrix(NA, nrow(vp)-2, 5)
rownames(P) <- rownames(vp)[-c(1:2)]
lfc <- matrix(NA, nrow(vp)-2, 5)
rownames(lfc) <- rownames(vp)[-c(1:2)]


for (i in 3:nrow(vp)) {
		res <- t.test(vp["Ctrl+TGFb",1:3],vp[i,1:3])
		P[i-2,1] <- res$p.value
		lfc[i-2,1] <- res$estimate[2]-res$estimate[1]
		
		res <- t.test(vp["Ctrl+TGFb",4:6],vp[i,4:6])
		P[i-2,2] <- res$p.value
		lfc[i-2,2] <- res$estimate[2]-res$estimate[1]
		
		res <- t.test(vp["Ctrl+TGFb",7:9],vp[i,7:9])
		P[i-2,3] <- res$p.value
		lfc[i-2,3] <- res$estimate[2]-res$estimate[1]
		
		res <- t.test(vp["Ctrl+TGFb",10:12],vp[i,10:12])
		P[i-2,4] <- res$p.value
		lfc[i-2,4] <- res$estimate[2]-res$estimate[1]
		
		res <- t.test(vp["Ctrl+TGFb",13:15],vp[i,13:15])
		P[i-2,5] <- res$p.value
		lfc[i-2,5] <- res$estimate[2]-res$estimate[1]
}

dff <- data.frame(log2FC=c(lfc[,1],lfc[,2],lfc[,3],lfc[,4],lfc[,5]), p.value=c(P[,1],P[,2],P[,3],P[,4],P[,5]))
dff$FDR <- p.adjust(dff$p.value, method="fdr")
dff$gene <- rep(rownames(P),5)
dff$type <- rep(c("COL1A1", "COL4A1", "COL4A2", "COL5A3", "COL7A1"), each=16)
dff$threshold <- FALSE
dff$threshold[dff$FDR<0.05 & dff$log2FC > 0] <- "Up" 
dff$threshold[dff$FDR<0.05 & dff$log2FC < -0] <- "Down"
dff$toshow <- c(rank(P[,1])<=5, rank(P[,2])<=5, rank(P[,3])<=5, rank(P[,4])<=5, rank(P[,5])<=5) # only show gene names for top5 based on p.values
dff <- dff[dff$type!="COL1A1", ] # remove one screen result

p1 <- ggplot(data=dff) + geom_point(aes(x=log2FC, y=-log10(FDR), colour=threshold)) +
	geom_vline(xintercept=0, linetype = "dotted", color="grey66") +
	geom_hline(yintercept=-log10(0.05), linetype = "dotted", color="red") +
	scale_color_manual(values=c("Up"="#ff4040", "Down"="#4287f5", "FALSE"="black")) +
	theme_bw() + xlab("log2FC") + ylab("-log10(FDR)") + facet_grid(.~type) + xlim(-3.6, 3.6) +
	theme(strip.text.x = element_text(size=14), legend.position="none") + 
	geom_text_repel(data=dff[dff$toshow,], box.padding=0.75, size=4, segment.size = 0.2, 
	aes(x=log2FC, y=-log10(FDR), label=as.character(dff$gene[dff$toshow])))


vp <- read.table("./data/Revision KD qPCR results.txt", sep="\t", row.names=1, header=T)
vp[vp==0] <- NA
vp <- log2(as.matrix(vp))

P <- matrix(NA, nrow(vp)-2, 4)
rownames(P) <- rownames(vp)[-c(1:2)]
lfc <- matrix(NA, nrow(vp)-2, 4)
rownames(lfc) <- rownames(vp)[-c(1:2)]


for (i in 3:nrow(vp)) {
		res <- t.test(vp["Ctrl+TGFb",1:3],vp[i,1:3])
		P[i-2,1] <- res$p.value
		lfc[i-2,1] <- res$estimate[2]-res$estimate[1]
		
		res <- t.test(vp["Ctrl+TGFb",4:6],vp[i,4:6])
		P[i-2,2] <- res$p.value
		lfc[i-2,2] <- res$estimate[2]-res$estimate[1]
		
		if(i!=12 & i!=16){
		res <- t.test(vp["Ctrl+TGFb",7:9],vp[i,7:9])
		P[i-2,3] <- res$p.value
		lfc[i-2,3] <- res$estimate[2]-res$estimate[1]
		}
		
		res <- t.test(vp["Ctrl+TGFb",10:12],vp[i,10:12])
		P[i-2,4] <- res$p.value
		lfc[i-2,4] <- res$estimate[2]-res$estimate[1]
}

dff <- data.frame(log2FC=c(lfc[,1],lfc[,2],lfc[,3],lfc[,4]), p.value=c(P[,1],P[,2],P[,3],P[,4]))
dff$FDR <- p.adjust(dff$p.value, method="fdr")
dff$gene <- rep(rownames(P),4)
dff$type <- rep(c("HHIP", "EDNRB", "MYOCD", "LTBP1"), each=16)
dff$threshold <- FALSE
dff$threshold[dff$FDR<0.05 & dff$log2FC > 0] <- "Up" 
dff$threshold[dff$FDR<0.05 & dff$log2FC < -0] <- "Down"
dff$toshow <- c(rank(P[,1])<=5, rank(P[,2])<=5, rank(P[,3])<=5, rank(P[,4])<=5) # only show gene names for top5 based on p.values

p2 <- ggplot(data=dff) + geom_point(aes(x=log2FC, y=-log10(FDR), colour=threshold)) +
	geom_vline(xintercept=0, linetype = "dotted", color="grey66") +
	geom_hline(yintercept=-log10(0.05), linetype = "dotted", color="red") +
	scale_color_manual(values=c("Up"="#ff4040", "Down"="#4287f5", "FALSE"="black")) +
	theme_bw() + xlab("log2FC") + ylab("-log10(FDR)") + facet_grid(.~type) + #xlim(-3.6, 3.6) +
	theme(strip.text.x = element_text(size=14), legend.position="none") + 
	geom_text_repel(data=dff[dff$toshow,], box.padding=0.35, size=4, segment.size = 0.2, 
	aes(x=log2FC, y=-log10(FDR), label=as.character(dff$gene[dff$toshow])))

pdf("#TI.myo_screens_volcano.pdf", 12.5, 7)
plot_grid(p1, p2, nrow=2)
dev.off()
