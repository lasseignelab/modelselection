# SCRIPT ORDER
1. _DataDownload.Rmd script with rec3_rse_download.R function script_ - Download data from recount3 with some pre-processing (removing/relabeling incorrect data)
2. _TumPurity.Rmd script_ - Compute what genes are highly associated with tumor purity, renders purity_genes obj utilized in #3
3. _SubsetRSE.Rmd script with rse_to_subset.R function script_ - Create and store individual RSE subset objects to compare model performance between different subset sizes, when tumor purity is in/excluded, and when gene variation rankings are based on GTEx or TCGA
4. _MetaTable.Rmd with barebones_cellline.R function script_ - Create a master matrix of metadata for categorizing data down the line
5. 
