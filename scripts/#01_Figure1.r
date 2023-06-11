###########################################################
# Raw data was preprocessed using CellRanger and Cumulus by 
# default settings. Fig. 1B-C were generated as part of
# this, the data for which is part of the SCP data.
# To have a quick start please use #00_Prepare_data.r first 
############################################################
# Figure 1. D PCoAs
############################################################

library(ggplot2)
library(cowplot)
library(reshape2)
library(labdsv)

source("./data/map_colors_lvls.r")

# Load metadata
df <- read.table('/data/co_ti_cmb_metadata.txt', sep='\t', header=T, row.names=1, stringsAsFactors=F)
dim(df) #720633     10
length(unique(df$anno2)) #66

propmat1 <- dcast(df, anno2 ~ PubIDSample)

propmat <- matrix(as.numeric(unlist(propmat1[,-1])), nrow=nrow(propmat1)) # convert to numeric
rownames(propmat) <- propmat1$anno2
colnames(propmat) <- colnames(propmat1)[-1]

propmat <- sweep(propmat, 2, colSums(propmat), FUN="/")
D <- dsvdis(t(propmat), index="bray/curtis")
dim(as.matrix(D)) #225 225

pc <- pco(D)

dfc <- df[!duplicated(df$PubIDSample),]
dfc <- dfc[match(colnames(propmat),dfc$PubIDSample),]
dfc$pc1 <- pc$points[,1]
dfc$pc2 <- pc$points[,2]

varexp <- pc$eig / sum(pmax(0, pc$eig))

p1 <- ggplot(data=dfc) +
	geom_point(aes(x=pc1, y=pc2, fill=Type), shape=21, size=3, color="black", stroke=0.2) + theme_bw() +
	scale_fill_manual(values=c(Heal="#95d962", NonI="#628cd9", Infl="#d9628c")) +
	xlab(sprintf("PCo 1 (%.1f%%)", 100*varexp[1])) +
	ylab(sprintf("PCo 2 (%.1f%%)", 100*varexp[2])) + ggtitle("Type") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p2 <- ggplot(data=dfc) +
	geom_point(aes(x=pc1, y=pc2, fill=Layer), shape=21, size=3, color="black", stroke=0.2) + theme_bw() +
	scale_fill_manual(values=c(E="#42cbf5", L="#f5428d", N="#ffde05")) +
	xlab(sprintf("PCo 1 (%.1f%%)", 100*varexp[1])) +
	ylab(sprintf("PCo 2 (%.1f%%)", 100*varexp[2])) + ggtitle("Layer") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p3 <- ggplot(data=dfc) +
	geom_point(aes(x=pc1, y=pc2, fill=Site), shape=21, size=3, color="black", stroke=0.2) + theme_bw() +
	scale_fill_manual(values=c(TI="orange", SB="#fc5d00", CO="purple")) +
	xlab(sprintf("PCo 1 (%.1f%%)", 100*varexp[1])) +
	ylab(sprintf("PCo 2 (%.1f%%)", 100*varexp[2])) + ggtitle("Location") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p4 <- ggplot(data=dfc) +
	geom_point(aes(x=pc1, y=pc2, fill=Chem), shape=21, size=3, color="black", stroke=0.2) + theme_bw() +
	scale_fill_manual(values=c(v1="#edfa5c", v3="#ffa1bb", v2="#5e9900")) +
	xlab(sprintf("PCo 1 (%.1f%%)", 100*varexp[1])) +
	ylab(sprintf("PCo 2 (%.1f%%)", 100*varexp[2])) + ggtitle("Chemistry") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p <- plot_grid(p2, p3, p1, nrow=3)
pdf("fig_Pco_6x15_vertical.pdf", 6, 15)
print(p)
dev.off()

############################################################
# Figure 1. E-F Compositional barplots
############################################################

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
library(scales)

source("./data/map_colors_lvls.r")
source("./scripts/dirichlet_functions.r")

# Load metadata
meta <- read.table('./data/co_ti_cmb_metadata.txt', sep='\t', header=T, row.names=1, stringsAsFactors=F)
dim(meta) #720633     10

###########################
# for TI
df <- meta[meta$Site!="CO",]
dim(df)  # 430903     10
titl_loc <- "TI"

###########################
# create a new column change anno2 from factor to character
df$anno22 <- as.character(df$anno2)

# Count each cell subset in every sample in TI
epi.freq = as.matrix(as.data.frame.matrix(table(df$PubIDSample[df$anno_overall=="Epithelial cells"], df$anno22[df$anno_overall=="Epithelial cells"])))
fib.freq = as.matrix(as.data.frame.matrix(table(df$PubIDSample[df$anno_overall=="Stromal cells"], df$anno22[df$anno_overall=="Stromal cells"])))
imm.freq = as.matrix(as.data.frame.matrix(table(df$PubIDSample[df$anno_overall=="Immune cells"], df$anno22[df$anno_overall=="Immune cells"])))

epi.tests <- cell_proportions_tests(epi.freq, "L")
fib.tests <- cell_proportions_tests(fib.freq, "E")
imm.tests <- cell_proportions_tests(imm.freq, "E")


## Epithelial compartment
tmp <- epi.freq[rownames(epi.tests$pct),colnames(epi.tests$pct)]
p3 <- matrix_barplot(as.matrix(epi.tests$pct), tmp, group_by=epi.tests$cov$condition, pvals=epi.tests$qvals[!grepl("^tissue", rownames(epi.tests$qvals)),], 
		xlab="", ylab="",  colors=set.colors, legend.show=T, sig_only=T, sig_cut=0.05) + ggtitle("Epithelial cells")

## Stromal compartment
tmp <- fib.freq[rownames(fib.tests$pct),colnames(fib.tests$pct)]
p2 <- matrix_barplot(as.matrix(fib.tests$pct), tmp, group_by=fib.tests$cov$condition, pvals=fib.tests$qvals[!grepl("^tissue", rownames(fib.tests$qvals)),], 
		xlab="", ylab="", colors=set.colors, legend.show=F, sig_only=T, sig_cut=0.05) + ggtitle("Stromal cells")

## separate for immune compartment
cond1 <- colnames(imm.tests$pct)%in%cell_types_anno$anno2[cell_types_anno$manu_orig=="T cells"]
tmp <- imm.freq[rownames(imm.tests$pct),colnames(imm.tests$pct)]
tmp <- tmp[,cond1]
p1a <- matrix_barplot(as.matrix(imm.tests$pct[,cond1]), tmp,
		group_by=imm.tests$cov$condition, pvals=imm.tests$qvals[!grepl("^tissue", rownames(imm.tests$qvals)),cond1], 
		xlab="", ylab="", colors=set.colors, legend.show=F, sig_only=T, sig_cut=0.05) + ggtitle("T cells")

cond2 <- colnames(imm.tests$pct)%in%cell_types_anno$anno2[cell_types_anno$manu_orig=="B cells"]
tmp <- imm.freq[rownames(imm.tests$pct),colnames(imm.tests$pct)]
tmp <- tmp[,cond2]
p1b <- matrix_barplot(as.matrix(imm.tests$pct[,cond2]), tmp,
		group_by=imm.tests$cov$condition, pvals=imm.tests$qvals[!grepl("^tissue", rownames(imm.tests$qvals)),cond2], 
		xlab="", ylab=sprintf("Percent of sample (%s)", titl_loc), colors=set.colors, legend.show=F, sig_only=T, sig_cut=0.05) + ggtitle("B cells")

cond3 <- colnames(imm.tests$pct)%in%cell_types_anno$anno2[cell_types_anno$manu_orig=="Myeloid cells"]
tmp <- imm.freq[rownames(imm.tests$pct),colnames(imm.tests$pct)]
tmp <- tmp[,cond3]
p1c <- matrix_barplot(as.matrix(imm.tests$pct[,cond3]), tmp,
		group_by=imm.tests$cov$condition, pvals=imm.tests$qvals[!grepl("^tissue", rownames(imm.tests$qvals)),cond3],
		xlab="", ylab="", colors=set.colors, legend.show=F, sig_only=T, sig_cut=0.05) + ggtitle("Myeloid cells")

p <- ggarrange(p1b, p1a, p1c, p2, p3, nrow=1, widths=c(1.5, 3, 14, 17, 9))
pdf("Cell_composition_barplot_stat.TI_20x4.pdf", 20, 4)
print(p)
dev.off()

###########################
# for CO
df <- meta[meta$Site=="CO",]
dim(df)  # 289730     10
titl_loc <- "CO"

# create a new column change anno2 from factor to character
df$anno22 <- as.character(df$anno2)

# Count each cell subset in every sample in TI
epi.freq = as.matrix(as.data.frame.matrix(table(df$PubIDSample[df$anno_overall=="Epithelial cells"], df$anno22[df$anno_overall=="Epithelial cells"])))
fib.freq = as.matrix(as.data.frame.matrix(table(df$PubIDSample[df$anno_overall=="Stromal cells"], df$anno22[df$anno_overall=="Stromal cells"])))
imm.freq = as.matrix(as.data.frame.matrix(table(df$PubIDSample[df$anno_overall=="Immune cells"], df$anno22[df$anno_overall=="Immune cells"])))

epi.tests <- cell_proportions_tests(epi.freq, "L")
fib.tests <- cell_proportions_tests(fib.freq, "E")
imm.tests <- cell_proportions_tests(imm.freq, "E")


## Epithelial compartment
tmp <- epi.freq[rownames(epi.tests$pct),colnames(epi.tests$pct)]
p3 <- matrix_barplot(as.matrix(epi.tests$pct), tmp, group_by=epi.tests$cov$condition, pvals=epi.tests$qvals[!grepl("^tissue", rownames(epi.tests$qvals)),], 
		xlab="", ylab="",  colors=set.colors, legend.show=T, sig_only=T, sig_cut=0.05) + ggtitle("Epithelial cells")

## Stromal compartment
tmp <- fib.freq[rownames(fib.tests$pct),colnames(fib.tests$pct)]
p2 <- matrix_barplot(as.matrix(fib.tests$pct), tmp, group_by=fib.tests$cov$condition, pvals=fib.tests$qvals[!grepl("^tissue", rownames(fib.tests$qvals)),], 
		xlab="", ylab="", colors=set.colors, legend.show=F, sig_only=T, sig_cut=0.05) + ggtitle("Stromal cells")

## separate for immune compartment
cond1 <- colnames(imm.tests$pct)%in%cell_types_anno$anno2[cell_types_anno$manu_orig=="T cells"]
tmp <- imm.freq[rownames(imm.tests$pct),colnames(imm.tests$pct)]
tmp <- tmp[,cond1]
p1a <- matrix_barplot(as.matrix(imm.tests$pct[,cond1]), tmp,
		group_by=imm.tests$cov$condition, pvals=imm.tests$qvals[!grepl("^tissue", rownames(imm.tests$qvals)),cond1], 
		xlab="", ylab="", colors=set.colors, legend.show=F, sig_only=T, sig_cut=0.05) + ggtitle("T cells")

cond2 <- colnames(imm.tests$pct)%in%cell_types_anno$anno2[cell_types_anno$manu_orig=="B cells"]
tmp <- imm.freq[rownames(imm.tests$pct),colnames(imm.tests$pct)]
tmp <- tmp[,cond2]
p1b <- matrix_barplot(as.matrix(imm.tests$pct[,cond2]), tmp,
		group_by=imm.tests$cov$condition, pvals=imm.tests$qvals[!grepl("^tissue", rownames(imm.tests$qvals)),cond2], 
		xlab="", ylab=sprintf("Percent of sample (%s)", titl_loc), colors=set.colors, legend.show=F, sig_only=T, sig_cut=0.05) + ggtitle("B cells")

cond3 <- colnames(imm.tests$pct)%in%cell_types_anno$anno2[cell_types_anno$manu_orig=="Myeloid cells"]
tmp <- imm.freq[rownames(imm.tests$pct),colnames(imm.tests$pct)]
tmp <- tmp[,cond3]
p1c <- matrix_barplot(as.matrix(imm.tests$pct[,cond3]), tmp,
		group_by=imm.tests$cov$condition, pvals=imm.tests$qvals[!grepl("^tissue", rownames(imm.tests$qvals)),cond3],
		xlab="", ylab="", colors=set.colors, legend.show=F, sig_only=T, sig_cut=0.05) + ggtitle("Myeloid cells")

library(ggpubr)
p <- ggarrange(p1b, p1a, p1c, p2, p3, nrow=1, align="h", widths=c(2,6,3.5,8.5,12.5))
pdf("Cell_composition_barplot_stat.CO_20x4.pdf", 20, 4.5)
print(p)
dev.off()
