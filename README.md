# Evaluating cancer cell line and patient-derived xenograft recapitulation of tumor and non-diseased tissue gene expression profiles
## Authors: Avery S. Williams, Elizabeth J. Wilk, Jennifer L. Fisher, Brittany N. Lasseigne, Ph.D.

## Overview
> Preclinical models like cancer cell lines and patient-derived xenografts (PDXs) are vital for studying disease mechanisms and evaluating treatment options. It is essential that they accurately recapitulate the disease state of interest to generate results that will translate in the clinic. Prior studies have demonstrated that preclinical models do not recapitulate all biological aspects of human tissues, particularly with respect to the tissue of origin gene expression signatures. Therefore, it is critical to assess how well preclinical model gene expression profiles correlate with human cancer tissues to inform preclinical model selection and data analysis decisions. Here we evaluated how well preclinical models recapitulate human cancer and non-diseased tissue gene expression patterns with respect to the most variable genes, tumor purity, and tissue specificity by using publicly available gene expression profiles across multiple sources. We found that using the full gene set improves correlations between preclinical model and tissue global gene expression profiles, confirmed that GBM PDX global gene expression correlation to GBM tumor global gene expression outperforms GBM cell line to GBM tumor global gene expression correlations, and demonstrated that preclinical models in our study often failed to reproduce tissue-specific expression. While including additional genes for global gene expression comparison between cell lines and tissues decreases the overall correlation, it improves the relative rank between a cell line and its tissue of origin compared to other tissues. Our findings underscore the importance of using the full gene expression set measured when comparing preclinical models and tissues and confirm that tissue-specific patterns are better preserved in GBM PDX models than in GBM cell lines. Future studies can build on these findings to determine the specific pathways and gene sets recapitulated by particular preclinical models to facilitate model selection for a given study design or goal. 

## Data Sets
* RNA-seq gene counts from [recount3](https://rna.recount.bio/)
  * The Cancer Genome Atlas (project_home = "data_sources/tcga") patient tumor tissue projects BRCA, KIRC, LUAD, UCEC, THCA, PRAD, LUSC, HNSC, COAD, LGG, SKCM, LAML, STAD, BLCA, OV, LIHC, KIRP, CESC, SARC, ESCA, PCPG, PAAD, READ, GBM, TGCT, THYM, KICH, MESO, UVM, ACC, UCS, DLBC, and CHOL
  * The Genotype-Tissue Expression (GTEx) project (project_home =	"data_sources/gtex") non-diseased tissue projects BRAIN, SKIN, ESOPHAGUS, BLOOD, BLOOD_VESSEL, ADIPOSE_TISSUE, HEART, MUSCLE, COLON, THYROID, NERVE, LUNG, BREAST, TESTIS, STOMACH, PANCREAS, PITUITARY, ADRENAL_GLAND, PROSTATE, SPLEEN, LIVER, OVARY, SMALL_INTESTINE, SALIVARY_GLAND, VAGINA, UTERUS, KIDNEY, BLADDER, CERVIX_UTERI, FALLOPIAN_TUBE
  * The Human Protein Atlas (HPA) (project_home = "data_sources/sra") non-diseased tissue project ERP003613
  * The Human Protein Atlas (HPA) (project_home = "data_sources/sra") cancer cell line project SRP017465
  * The Cancer Cell Line Encyclopedia (CCLE) (project_home = "data_sources/sra") cancer cell line project SRP186687
  * Scott MC et al. (project_home = "data_sources/sra") osteosarcoma tumor, non-diseased bone, and osteosarcoma cell line project SRP090849
* Tissue-Specific genes from GTEx (download used in this analysis in github [here](https://github.com/lasseignelab/modelselection/blob/main/data/) as gene_set_library_up_crisp.gmt and gene_set_library_dn_crisp.gmt), sourced from [here](https://maayanlab.cloud/Harmonizome/dataset/GTEx+Tissue+Gene+Expression+Profiles)
* Tumor Purity information (download used in this analysis in github [here](https://github.com/lasseignelab/modelselection/blob/main/data/) as TCGA_mastercalls.abs_tables_JSedit.fixed.csv), sourced from [here](https://gdc.cancer.gov/about-data/publications/pancanatlas)
* Problematic cell lines information from The Cellosaurus [here](https://www.cellosaurus.org/search?input=%27problematic%20cell%20line%27)
  
## Docker 
The docker image used to run these analyses (rstudio_modelselection:1.2.1) is available on DockerHub [here](https://hub.docker.com/r/lizzyr/rstudio_modelselection)

## Script Order
1. _DataDownload.Rmd script with rec3_rse_download.R function script_ - Download data from recount3 with some pre-processing (removing/relabeling incorrect data)
2. _TumPurity.Rmd script_ - Compute what genes are highly associated with tumor purity, renders purity_genes obj utilized in #3
3. _SubsetRSE.Rmd script with rse_to_subset.R function script_ - Create and store individual RSE subset objects to compare model performance between different subset sizes, when tumor purity is in/excluded, and when gene variation rankings are based on GTEx or TCGA
4. _MetaTable.Rmd with barebones_cellline.R function script_ - Create a master matrix of metadata for categorizing data down the line
5. _Cor_by_GeneSetSize.Rmd with tissuecormatch.R and tissuecormatch_median.R function scripts_ - Find correlations between cell lines and tissues and how that differs in different gene set sizes
6. _GlobalandVariableCor.Rmd script_ - Comparing model performance for global and most variable gene correlations in tumor and healthy tissues
7. _TissueEnriched_GTEx.Rmd script with tisgenesetcor.R and geom_split_violin.R function scripts_ - Evaluating model performance for GTEx tissue-associated genes in tumor and healthy tissues
8. _Ind_Gene_Pathways.Rmd script with cor_eachgene.R function script_ - Evaluating model performance within individual genes and what pathways these genes belong to

## Lasseigne Lab 
<a href="https://www.lasseigne.org/" target="_blank"><img src="https://www.lasseigne.org/img/main/lablogo.png" width="200" height="200"></a>

## License
[MIT](https://github.com/lasseignelab/modelselection/blob/main/LICENSE)

## Versions 
Available in Docker image (see above)

```
R version 4.2.2 (2022-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

 colorspace_2.1-0            rjson_0.2.21                ellipsis_0.3.2              qvalue_2.30.0               htmlTable_2.4.1            
 XVector_0.38.0              GenomicRanges_1.50.2        base64enc_0.1-3             rstudioapi_0.14             bit64_4.0.5                
 AnnotationDbi_1.60.2        fansi_1.0.4                 xml2_1.3.3                  splines_4.2.2               codetools_0.2-19           
 R.methodsS3_1.8.2           cachem_1.0.7                knitr_1.42                  jsonlite_1.8.4              Formula_1.2-5              
 Rsamtools_2.14.0            cluster_2.1.4               dbplyr_2.3.1                png_0.1-8                   R.oo_1.25.0                
 rentrez_1.2.3               readr_2.1.4                 compiler_4.2.2              httr_1.4.5                  backports_1.4.1            
 Matrix_1.5-3                fastmap_1.1.1               limma_3.54.2                cli_3.6.0                   htmltools_0.5.4            
 prettyunits_1.1.1           tools_4.2.2                 gtable_0.3.1                glue_1.6.2                  GenomeInfoDbData_1.2.9     
 reshape2_1.4.4              dplyr_1.1.0                 rappdirs_0.3.3              doRNG_1.8.6                 Rcpp_1.0.10                
 Biobase_2.58.0              bumphunter_1.40.0           vctrs_0.6.0                 Biostrings_2.66.0           rtracklayer_1.58.0         
 iterators_1.0.14            recount3_1.8.0              xfun_0.37                   stringr_1.5.0               lifecycle_1.0.3            
 restfulr_0.0.15             rngtools_1.5.2              XML_3.99-0.13               zlibbioc_1.44.0             scales_1.2.1               
 BSgenome_1.66.3             VariantAnnotation_1.44.1    hms_1.1.2                   MatrixGenerics_1.10.0       parallel_4.2.2             
 SummarizedExperiment_1.28.0 GEOquery_2.66.0             derfinderHelper_1.32.0      yaml_2.3.7                  curl_5.0.0                 
 memoise_2.0.1               gridExtra_2.3               downloader_0.4              ggplot2_3.4.1               biomaRt_2.54.0             
 recount_1.24.1              rpart_4.1.19                stringi_1.7.12              RSQLite_2.3.0               S4Vectors_0.36.2           
 BiocIO_1.8.0                foreach_1.5.2               checkmate_2.1.0             GenomicFeatures_1.50.4      BiocGenerics_0.44.0        
 filelock_1.0.2              BiocParallel_1.32.5         GenomeInfoDb_1.34.9         rlang_1.1.0                 pkgconfig_2.0.3            
 GenomicFiles_1.34.0         matrixStats_0.63.0          bitops_1.0-7                evaluate_0.20               lattice_0.20-45            
 purrr_1.0.1                 GenomicAlignments_1.34.1    htmlwidgets_1.6.1           bit_4.0.5                   tidyselect_1.2.0           
plyr_1.8.8                  magrittr_2.0.3              R6_2.5.1                    IRanges_2.32.0              generics_0.1.3             
Hmisc_5.0-1                 DelayedArray_0.24.0         DBI_1.1.3                   pillar_1.8.1                foreign_0.8-84             
KEGGREST_1.38.0             RCurl_1.98-1.10             nnet_7.3-18                 tibble_3.2.0                crayon_1.5.2               
derfinder_1.32.0            utf8_1.2.3                  BiocFileCache_2.6.1         tzdb_0.3.0                  rmarkdown_2.20             
progress_1.2.2              locfit_1.5-9.7              grid_4.2.2                  data.table_1.14.8           blob_1.2.3                 
digest_0.6.31               tidyr_1.3.0                 R.utils_2.12.2              stats4_4.2.2                munsell_0.5.0              
sessioninfo_1.2.2  
```
