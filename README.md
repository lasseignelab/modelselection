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
  
## Script Order
1. _DataDownload.Rmd script with rec3_rse_download.R function script_ - Download data from recount3 with some pre-processing (removing/relabeling incorrect data)
2. _TumPurity.Rmd script_ - Compute what genes are highly associated with tumor purity, renders purity_genes obj utilized in #3
3. _SubsetRSE.Rmd script with rse_to_subset.R function script_ - Create and store individual RSE subset objects to compare model performance between different subset sizes, when tumor purity is in/excluded, and when gene variation rankings are based on GTEx or TCGA
4. _MetaTable.Rmd with barebones_cellline.R function script_ - Create a master matrix of metadata for categorizing data down the line
5. _Cor_by_GeneSetSize.Rmd with tissuecormatch.R and tissuecormatch_median.R function scripts_ - Find correlations between cell lines and tissues and how that differs in different gene set sizes
6. _GlobalandVariableCor.Rmd script_ - Comparing model performance for global and most variable gene correlations in tumor and healthy tissues
7. _TissueEnriched_GTEx.Rmd script with tisgenesetcor.R and geom_split_violin.R function scripts_ - Evaluating model performance for GTEx tissue-associated genes in tumor and healthy tissues
8. _Ind_Gene_Pathways.Rmd script with cor_eachgene.R function script_ - Evaluating model performance within individual genes and what pathways these genes belong to

## License
[MIT](https://github.com/lasseignelab/modelselection/blob/main/LICENSE)
