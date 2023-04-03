##################################################################################################
# gather & download all the mtx files from SCP and merge to one combined seurat object

library(Seurat)
library(Matrix)

################################
# Colon stromal counts

str.counts <- readMM('CO_STR.scp.raw.mtx')
rownames(str.counts) <- read.table('CO_STR.scp.features.tsv')[,2] # take the genenames
colnames(str.counts) <- readLines('CO_STR.scp.barcodes.tsv')
dim(str.counts) # 28663 39433
str.seur <- CreateSeuratObject(counts=str.counts, min.cells=0, min.features=0, names.field=1, names.delim='\\.')

# readin normalized
str.data <- readMM('CO_STR.scp.matrix.mtx')
rownames(str.data) <- read.table('CO_STR.scp.features.tsv')[,2] # take the genenames
colnames(str.data) <- readLines('CO_STR.scp.barcodes.tsv')
dim(str.data) # 28663 39433
str.data <- CreateSeuratObject(counts=str.data, min.cells=0, min.features=0, names.field=1, names.delim='\\.')

# combine with the raw counts with the noralized data
str.seur <- SetAssayData(object=str.seur,  slot="data", new.data=GetAssayData(str.data, assay = "RNA", slot = "data"), assay= "RNA")
str.seur <- UpdateSeuratObject(str.seur)
dim(str.seur) # 28663 39433

rm(str.counts, str.data); gc()
saveRDS(str.seur, "./data/co.str.seur.rds")


################################
# Colon immune counts

imm.counts <- readMM('CO_IMM.scp.raw.mtx')
rownames(imm.counts) <- read.table('CO_IMM.scp.features.tsv')[,2] # take the genenames
colnames(imm.counts) <- readLines('CO_IMM.scp.barcodes.tsv')
dim(imm.counts) #28663 152509
imm.seur <- CreateSeuratObject(counts=imm.counts, min.cells=0, min.features=0, names.field=1, names.delim='\\.')

# readin normalized
imm.data <- readMM('CO_IMM.scp.matrix.mtx')
rownames(imm.data) <- read.table('CO_IMM.scp.features.tsv')[,2] # take the genenames
colnames(imm.data) <- readLines('CO_IMM.scp.barcodes.tsv')
dim(imm.data) #28663 152509
imm.data <- CreateSeuratObject(counts=imm.data, min.cells=0, min.features=0, names.field=1, names.delim='\\.')

# combine with the raw counts with the noralized data
imm.seur <- SetAssayData(object=imm.seur,  slot="data", new.data=GetAssayData(imm.data, assay = "RNA", slot = "data"), assay= "RNA")
imm.seur <- UpdateSeuratObject(imm.seur)
dim(imm.seur) #28663 152509

rm(imm.counts, imm.data); gc()
saveRDS(imm.seur, "./data/co.imm.seur.rds")

################################
# Colon epithelial counts

epi.counts <- readMM('CO_EPI.scp.raw.mtx')
rownames(epi.counts) <- read.table('CO_EPI.scp.features.tsv')[,2] # take the genenames
colnames(epi.counts) <- readLines('CO_EPI.scp.barcodes.tsv')
dim(epi.counts) #28663 97788
epi.seur <- CreateSeuratObject(counts=epi.counts, min.cells=0, min.features=0, names.field=1, names.delim='\\.')

# readin normalized
epi.data <- readMM('CO_EPI.scp.matrix.mtx')
rownames(epi.data) <- read.table('CO_EPI.scp.features.tsv')[,2] # take the genenames
colnames(epi.data) <- readLines('CO_EPI.scp.barcodes.tsv')
dim(epi.data) #28663 97788
epi.data <- CreateSeuratObject(counts=epi.data, min.cells=0, min.features=0, names.field=1, names.delim='\\.')

# combine with the raw counts with the noralized data
epi.seur <- SetAssayData(object=epi.seur,  slot="data", new.data=GetAssayData(epi.data, assay = "RNA", slot = "data"), assay= "RNA")
epi.seur <- UpdateSeuratObject(epi.seur)
dim(epi.seur) #28663 97788

rm(epi.counts, epi.data); gc()
saveRDS(epi.seur, "./data/co.epi.seur.rds")

#########################################
## combine all obj to make a colon obj

library(Seurat)
library(Matrix)

epi.seur <- readRDS("./data/co.epi.seur.rds")
str.seur <- readRDS("./data/co.str.seur.rds")
imm.seur <- readRDS("./data/co.imm.seur.rds")

co.seur <- merge(epi.seur, str.seur)
co.seur <- merge(co.seur, imm.seur)

dim(co.seur) #  21454 365491
saveRDS(co.seur, "./data/co.seur.rds")


################################
# TI stromal counts

str.counts <- readMM('TI_STR.scp.raw.mtx')
rownames(str.counts) <- read.table('TI_STR.scp.features.tsv')[,2] # take the genenames
colnames(str.counts) <- readLines('TI_STR.scp.barcodes.tsv')
dim(str.counts) # 28923 75695
str.seur <- CreateSeuratObject(counts=str.counts, min.cells=0, min.features=0, names.field=1, names.delim='\\.')

# readin normalized
str.data <- readMM('TI_STR.scp.matrix.mtx')
rownames(str.data) <- read.table('TI_STR.scp.features.tsv')[,2] # take the genenames
colnames(str.data) <- readLines('TI_STR.scp.barcodes.tsv')
dim(str.data) # 28923 75695
str.data <- CreateSeuratObject(counts=str.data, min.cells=0, min.features=0, names.field=1, names.delim='\\.')

# combine with the raw counts with the noralized data
str.seur <- SetAssayData(object=str.seur,  slot="data", new.data=GetAssayData(str.data, assay = "RNA", slot = "data"), assay= "RNA")
str.seur <- UpdateSeuratObject(str.seur)
dim(str.seur) # 28923 75695

rm(str.counts, str.data); gc()
saveRDS(str.seur, "./data/ti.str.seur.rds")


################################
# TI immune counts

imm.counts <- readMM('TI_IMM.scp.raw.mtx')
rownames(imm.counts) <- read.table('TI_IMM.scp.features.tsv')[,2] # take the genenames
colnames(imm.counts) <- readLines('TI_IMM.scp.barcodes.tsv')
dim(imm.counts) #28923 201072
imm.seur <- CreateSeuratObject(counts=imm.counts, min.cells=0, min.features=0, names.field=1, names.delim='\\.')

# readin normalized
imm.data <- readMM('TI_IMM.scp.matrix.mtx')
rownames(imm.data) <- read.table('TI_IMM.scp.features.tsv')[,2] # take the genenames
colnames(imm.data) <- readLines('TI_IMM.scp.barcodes.tsv')
dim(imm.data) #28923 201072
imm.data <- CreateSeuratObject(counts=imm.data, min.cells=0, min.features=0, names.field=1, names.delim='\\.')

# combine with the raw counts with the noralized data
imm.seur <- SetAssayData(object=imm.seur,  slot="data", new.data=GetAssayData(imm.data, assay = "RNA", slot = "data"), assay= "RNA")
imm.seur <- UpdateSeuratObject(imm.seur)
dim(imm.seur) #28923 201072

rm(imm.counts, imm.data); gc()
saveRDS(imm.seur, "./data/ti.imm.seur.rds")

################################
# TI epithelial counts

epi.counts <- readMM('TI_EPI.scp.raw.mtx')
rownames(epi.counts) <- read.table('TI_EPI.scp.features.tsv')[,2] # take the genenames
colnames(epi.counts) <- readLines('TI_EPI.scp.barcodes.tsv')
dim(epi.counts) #28923 154136
epi.seur <- CreateSeuratObject(counts=epi.counts, min.cells=0, min.features=0, names.field=1, names.delim='\\.')

# readin normalized
epi.data <- readMM('TI_EPI.scp.matrix.mtx')
rownames(epi.data) <- read.table('TI_EPI.scp.features.tsv')[,2] # take the genenames
colnames(epi.data) <- readLines('TI_EPI.scp.barcodes.tsv')
dim(epi.data) #28923 154136
epi.data <- CreateSeuratObject(counts=epi.data, min.cells=0, min.features=0, names.field=1, names.delim='\\.')

# combine with the raw counts with the noralized data
epi.seur <- SetAssayData(object=epi.seur,  slot="data", new.data=GetAssayData(epi.data, assay = "RNA", slot = "data"), assay= "RNA")
epi.seur <- UpdateSeuratObject(epi.seur)
dim(epi.seur) #28923 154136

rm(epi.counts, epi.data); gc()
saveRDS(epi.seur, "./data/ti.epi.seur.rds")

#########################################
## combine all obj to make a TI obj

library(Seurat)
library(Matrix)

epi.seur <- readRDS("./data/ti.epi.seur.rds")
str.seur <- readRDS("./data/ti.str.seur.rds")
imm.seur <- readRDS("./data/ti.imm.seur.rds")

ti.seur <- merge(epi.seur, str.seur)
ti.seur <- merge(ti.seur, imm.seur)

dim(ti.seur) #  28923 430903
saveRDS(ti.seur, "./data/ti.seur.rds")


#########################################
## combine Colon and TI objs

library(Seurat)
library(Matrix)

co <- readRDS("./data/co.seur.rds")
ti <- readRDS("./data/ti.seur.rds")
cd <- merge(ti, co)
dim(cd) # 29756 720633
rm(co, ti); gc()

# Load metadata
meta <- read.table('./data/co_ti_cmb_metadata.txt', sep='\t', header=T, row.names=1, stringsAsFactors=F)
dim(meta) #720633     10
length(unique(meta$anno2)) #66

meta.x <- cd@meta.data
# matched and add the extra meta info to seurat obj
meta <- meta[match(rownames(meta.x),rownames(meta)),]
cd@meta.data <- meta
saveRDS(cd, "./data/co_ti_cmb.rds") # final obj
