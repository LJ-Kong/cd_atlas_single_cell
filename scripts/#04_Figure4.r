library(Seurat)
library(ggplot2)
library(cowplot)
require("RColorBrewer")
require("viridis")
library(dplyr)
library(harmony)
library(RColorBrewer)


source("./data/map_colors_lvls.r")

############################################################
# Prepare EEC data
############################################################

result <- readRDS("./data/co_ti_cmb.rds")
meta <- result@meta.data
meta$location <- "TI"
meta$location[meta$Site=="CO"] <- "CO"
meta$loc_cy <- paste(meta$location, meta$anno2, sep="_")
Idents(result) <- meta$loc_cy

# subset to EEC cells
result.s <- subset(result, idents=c("TI_Enterochromaffin cells", "TI_L cells"))
dim(result.s) #29756  1813
# saveRDS(result.s, "ti_ec.cell.rds")
rm(result); gc()

pbmc <- NormalizeData(result.s) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
pbmc <- RunHarmony(pbmc, group.by.vars = "PubID")
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:20, k.param=20, force.recalc=T)
pbmc <- FindClusters(pbmc, resolution=1)
pbmc <- RunUMAP(pbmc, reduction = "harmony", n.neighbors=40, min.dist = 0.5, dims = 1:20)
pbmc <- subset(pbmc, subset = CHGA>5 & MUC2<5 & FABP1<15 & LYZ<15 & OLFM4<10, slot = 'counts')
dim(pbmc) #29756   802

DimPlot(pbmc, reduction = "umap", label = TRUE, label.size =10, pt.size = 2) + NoLegend()

## donor
df <- pbmc@meta.data[,c("PubID", "seurat_clusters", "Type")]
dfgg <- dplyr::count(df, seurat_clusters, PubID, .drop=FALSE)
ggplot(dfgg)+ labs(title="Donor.cell per cluster") + theme_bw() + scale_fill_manual(values=godsnot_64)+
	theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
	geom_bar(aes(x=seurat_clusters, y=n, fill=PubID), stat="identity", width=0.6)

## disease
dfgg <- dplyr::count(df, seurat_clusters, Type, .drop=FALSE)
celllabels_colors <- c(Heal = "#95d962", NonI = "#628cd9", Infl = "#d9628c")
ggplot(dfgg)+ labs(title="Type.cell per cluster") + theme_bw() +
	scale_fill_manual(values=celllabels_colors) +
	theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
	geom_bar(aes(x=seurat_clusters, y=n, fill=Type), stat="identity", width=0.6)


##################################################################
## cluster2, 4, 9, 10, 12, 14 are contaminant lineages: 
# FABP1-positive enterocytes
# OLFM4-positive stem cells 
# rare MUC2-positive goblet cells 
# LYZ/MMP7-positive Paneth cells
# and several progenitor populations (Figure S3D): APOA1

keep.ids <- c("0", "1", "3", "5", "6", "7", "8", "11", "13", "15", "16", "17")
pbmc <- subset(pbmc, idents=keep.ids)
dim(pbmc) #[1] 29756   706

# re-cluster
pbmc <- NormalizeData(pbmc) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
pbmc <- RunHarmony(pbmc, group.by.vars = "PubID")
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:20, k.param=6, force.recalc=T)
pbmc <- FindClusters(pbmc, resolution=4) # 25 clusters 
pbmc <- RunUMAP(pbmc, reduction="harmony", n.neighbors=80, min.dist=0.8, dims=1:20)

DimPlot(pbmc, reduction = "umap", label = TRUE, label.size =10, pt.size = 2) + NoLegend()

## checking some markers
keys <- c('CHGA', 'CHGB', "TFF1", "NPW", "REG4", "CES1", "TPH1", "KCTD12", "CDKN1C", "SOX4", "NEUROG3", "HES6", "NTS", "PYY", "LYPD8", 
				"GCG", "PRPH", "ISL1", "SST", "VTN", "HHEX", "CCK", "ONECUT3", "DMBT1", "GIP", "MYH10", "STK32A")
cpm_val.sub <- as.matrix(GetAssayData(pbmc, assay = "RNA", slot = "data")[keys,])
df <- as.data.frame(Embeddings(pbmc, reduction = "umap"))
cols.use = c("lightgrey","blue")
plotlist <- list()
for(i in 1:length(keys)){
	df$Expr <- cpm_val.sub[i,]
	dfOrd <- df[order(df$Expr),]
	ggp <- ggplot(data=dfOrd) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=Expr), size=1) + scale_color_viridis(option = "C") #scale_color_gradientn(colours = cols.use)
	ggp <- ggp + labs(title=keys[i])
	ggp <- ggp + theme_cowplot() + theme(plot.title=element_text(size = 20, face = "bold"), panel.grid = element_blank(), axis.text=element_blank(), 
				axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks=element_blank())
	plotlist <- c(plotlist, list(ggp))
}
plot_grid(plotlist=plotlist, nrow=4)


##############################################################################################
## rename clusters
# Idents(pbmc) <- pbmc@meta.data$seurat_clusters
new.cluster.ids <- c("N cells NTS+ PYY+", "EC THP1+ CES1+", "EC THP1+ CES1+", "N cells NTS+ PYY+", 
						"L cells GCG+", "EC REG4+ NPW+", "EC THP1+ CES1+", "EC REG4+ NPW+", "EC THP1+ CES1+", 
						"EC REG4+ NPW+", "EC REG4+ NPW+", "Progenitors NEUROG3+ SOX4+", "EC REG4+ NPW+", 
						"EC REG4+ NPW+", "I cells CCK+", "EC REG4+ NPW+", "EC REG4+ NPW+", "EC REG4+ NPW+", 
						"N cells NTS+ PYY+", "D cells SST+ HHEX+", "EC REG4+ NPW+", "Progenitors NEUROG3+ SOX4+", 
						"EC REG4+ NPW+", "K cells GIP+", "D cells SST+ HHEX+")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
#DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 6, repel = TRUE, pt.size = 2) + NoLegend()
DimPlot(pbmc, reduction="umap", label=F, cols=brewer.pal(8, "Accent"), repel=F, pt.size=2)

pbmc@meta.data$newLabel <- Idents(pbmc) # save new annotation
# saveRDS(pbmc, "./data/EECs706.TI.final.rds")

############################################################
# Figure 4. A-C
############################################################

library(Seurat)
library(ggplot2)
library(cowplot)
require("RColorBrewer")
require("viridis")
library(RColorBrewer)

source("./data/map_colors_lvls.r")
pbmc <- readRDS(file="./data/EECs706.TI.final.rds")

manu_features <- c('CHGA', 'CHGB', "TFF1", "NPW", "REG4", "CES1", "TPH1", "KCTD12", "CDKN1C", "SOX4", "HES6", "NTS", "PYY", "LYPD8", 
				"GCG", "PRPH", "ISL1", "SST", "VTN", "HHEX", "CCK", "ONECUT3", "DMBT1", "GIP", "MYH10", "STK32A")
Idents(pbmc) <- factor(Idents(pbmc), levels=rev(c("K cells GIP+", "I cells CCK+", "D cells SST+ HHEX+", "L cells GCG+", "N cells NTS+ PYY+",
											"Progenitors NEUROG3+ SOX4+", "EC THP1+ CES1+", "EC REG4+ NPW+")))


DimPlot(pbmc, reduction="umap", label=F, cols=brewer.pal(8, "Accent"), repel=F, pt.size=2)
# EECs706.TI_umap_anno_7.2x4.2.pdf

DotPlot(pbmc, features=manu_features, cols = c("lightgrey", "orange")) + RotatedAxis() + guides(color = guide_colorbar(title = 'Scaled Average\nExpression'))
# EECs706.TI_umap_markers_12x6.pdf

pbmc@meta.data$newLabel <- factor(pbmc@meta.data$newLabel, levels=rev(c("K cells GIP+", "I cells CCK+", "D cells SST+ HHEX+", "L cells GCG+", "N cells NTS+ PYY+",
											"Progenitors NEUROG3+ SOX4+", "EC THP1+ CES1+", "EC REG4+ NPW+")))
## donor
df <- pbmc@meta.data[,c("PubID", "newLabel", "Type")]
dfgg <- dplyr::count(df, newLabel, PubID, .drop=FALSE)
p5 <- ggplot(dfgg)+ ylab("Total number of cells") + theme_bw() + scale_fill_manual(values=godsnot_64)+
	theme(axis.title.x = element_blank(), 
	axis.text.x = element_text(angle = -45, hjust = 0)) +
	geom_bar(aes(x=newLabel, y=n, fill=as.character(PubID)), stat="identity", width=0.6)


## disease
dfgg <- dplyr::count(df, newLabel, Type, .drop=FALSE)
celllabels_colors <- c(Heal = "#95d962", NonI = "#628cd9", Infl = "#d9628c")
p6 <- ggplot(dfgg)+ ylab("Total number of cells") + theme_bw() +
	scale_fill_manual(values=celllabels_colors) +
	theme(axis.title.x = element_blank(), 
	axis.text.x = element_text(angle = -45, hjust = 0)) +
	geom_bar(aes(x=newLabel, y=n, fill=Type), stat="identity", width=0.6)

plot_grid(p5, p6, rel_widths=c(1.35,1))
# "EECs760.TI_barplots_11x6.5_noLegend.pdf"
