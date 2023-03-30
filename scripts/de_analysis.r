############################################################
# Differential expression analysis
############################################################

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggrepel)
library(reshape2)
library(lmerTest)
library(MAST)
library(data.table)

source("mast.de.functions.r")

result <- readRDS(file="./data/co_ti_cmb.rds")
meta <- result@meta.data
meta$location <- "TI"
meta$location[meta$Site=="CO"] <- "CO"
meta$loc_anno <- paste(meta$location, meta$anno_overall, sep="_")
meta$anno2 <- factor(meta$anno2)
result@meta.data <- meta

Idents(result) <- result@meta.data$loc_anno
unique(Idents(result))

########################################################
## DE is done by different compartment per location

###########-- TI 

ti.str <- subset(result, idents="TI_Stromal cells")
ti.epi <- subset(result, idents="TI_Epithelial cells")
ti.imm <- subset(result, idents="TI_Immune cells")

###########-- Colon 

co.str <- subset(result, idents="CO_Stromal cells")
co.imm <- subset(result, idents="CO_Immune cells")
co.epi <- subset(result, idents="CO_Epithelial cells")

########################################################
rm(result);gc()

# Formula for MAST
test.formula <- formula(~n_genes + Layer + Type + (1 | PubID) + (1 | PubIDSample))
test.formula.prefilter <- formula(~n_genes + Layer + Type)
test.contrasts <- c("TypeNonI", "TypeInfl")
subsample.pattern <- list(~Type, ~PubID, ~PubIDSample)
max.cells <- 10000
min_expression_frac <- 0.1

de.dir <- "de-results.out"
prefilter.P.thresh <- 0.05

#################################################################################
## Immune, TI

unique(ti.imm@meta.data$anno2)
for (celltype in unique(ti.imm@meta.data$anno2)) {
	run.celltype(ti.imm, "TI.Immune", celltype)
}

#################################################################################
## Stromal, TI

unique(ti.str@meta.data$anno2)
for (celltype in unique(ti.str@meta.data$anno2)) {
	run.celltype(ti.str, "TI.Stromal", celltype)
}

#################################################################################
## Epi, TI

unique(ti.epi@meta.data$anno2)
for (celltype in unique(ti.epi@meta.data$anno2)) {
	run.celltype(ti.epi, "TI.Epithelial", celltype)
}


################################################################################
#################################################################################
## Immune, CO

unique(co.imm@meta.data$anno2)
for (celltype in unique(co.imm@meta.data$anno2)) {
	run.celltype(co.imm, "CO.Immune", celltype)
}

#################################################################################
## Stromal, CO

unique(co.str@meta.data$anno2)
for (celltype in unique(co.str@meta.data$anno2)) {
	run.celltype(co.str, "CO.Stromal", celltype)
}

#################################################################################
## Epi, CO

unique(co.epi@meta.data$anno2)
# NEED to manual change "Goblet.cells.SPINK4+"
for (celltype in unique(co.epi@meta.data$anno2)) {
	run.celltype(co.epi, "CO.Epithelial", celltype)
}

########################################################
# after running all DE results
# put all DE results together

setwd(de.dir)
de.files <- list.files(pattern="*.csv")

de.cmb <- list()
for (fname in de.files) {
	fn <- file.path(de.dir, fname)
	if (file.exists(fn)) {
		de.res <- read.csv(fn)
		de.res$X <- substr(fname, 1, nchar(fname)-4)
		de.cmb <- c(de.cmb, list(de.res))
	}
}

de.cmb <- do.call(rbind, de.cmb)
de.cmb$celltype <- substr(de.cmb$X, 4, nchar(de.cmb$X)) #cn <- gsub("^\\w+\\.", "", cn)
de.cmb$location <- substr(de.cmb$X, 1, 2)
de.cmb$FDR.D <- p.adjust(ifelse(is.na(de.cmb$Pr..Chisq.D), de.cmb$Pr..Chisq.D.pre, de.cmb$Pr..Chisq.D), method='fdr')
de.cmb$FDR.C <- p.adjust(ifelse(is.na(de.cmb$Pr..Chisq.C), de.cmb$Pr..Chisq.C.pre, de.cmb$Pr..Chisq.C), method='fdr')
de.cmb$FDR <- p.adjust(ifelse(is.na(de.cmb$Pr..Chisq.), de.cmb$Pr..Chisq..pre, de.cmb$Pr..Chisq.), method='fdr')
de.cmb$direct.D <- de.cmb$coefD<0 # FALSE=UP, TRUE=DOWN
saveRDS(./data/de.cmb, "de.cmb.rds")
