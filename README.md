## Human CD atlas study between colon and terminal ileum

This repository includes main code to analyse the single cell data from [*'The landscape of immune dysregulation in Crohn’s Disease revealed through single-cell transcriptomic profiling in the ileum and colon'*](https://doi.org/10.1016/j.immuni.2023.01.002), which includes 720,633 quality cells from 71 donors with varying inflammation status and from both terminal ileum (TI) and colon.

### To start:
Download the expression data from the Single Cell Portal (SCP) at: https://singlecell.broadinstitute.org/single_cell/study/SCP1884

Specfically, download these files below:

- CO_STR.scp.raw.mtx, CO_STR.scp.matrix.mtx, CO_STR.scp.features.tsv, CO_STR.scp.barcodes.tsv
- CO_IMM.scp.raw.mtx, CO_IMM.scp.matrix.mtx, CO_IMM.scp.features.tsv, CO_IMM.scp.barcodes.tsv
- CO_EPI.scp.raw.mtx, CO_EPI.scp.matrix.mtx, CO_EPI.scp.features.tsv, CO_EPI.scp.barcodes.tsv
- TI_STR.scp.raw.mtx, TI_STR.scp.matrix.mtx, TI_STR.scp.features.tsv, TI_STR.scp.barcodes.tsv
- TI_IMM.scp.raw.mtx, TI_IMM.scp.matrix.mtx, TI_IMM.scp.features.tsv, TI_IMM.scp.barcodes.tsv
- TI_EPI.scp.raw.mtx, TI_EPI.scp.matrix.mtx, TI_EPI.scp.features.tsv, TI_EPI.scp.barcodes.tsv

Also, extra meta files can be downloaded from: [data folder](https://github.com/LJ-Kong/cd_atlas_single_cell/tree/main/data).
Once all these files were downloaded please use [/scripts/#00_Prepare_data.r](https://github.com/LJ-Kong/cd_atlas_single_cell/blob/main/scripts/%2300_PrepareData.r) to start the analysis. Please also note that all the analyses were done under R base v4.1.2 and Seurat v4.3.0.
