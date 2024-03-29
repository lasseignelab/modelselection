#check for rse, tpm, and raw_counts folders
checkfolders <- function(savefilepath){
  if(!dir.exists(paste0(savefilepath, "rse")) == TRUE){
    dir.create(paste0(savefilepath, "rse"))  
  }
  if(!dir.exists(paste0(savefilepath, "rse/tcga_tissue")) == TRUE){
    dir.create(paste0(savefilepath, "rse/tcga_tissue"))  
  }
  if(!dir.exists(paste0(savefilepath, "rse/gtex_tissue")) == TRUE){
    dir.create(paste0(savefilepath, "rse/gtex_tissue"))
  }
  if(!dir.exists(paste0(savefilepath, "tpm")) == TRUE){
    dir.create(paste0(savefilepath, "tpm"))  
  }
  if(!dir.exists(paste0(savefilepath, "raw_counts")) == TRUE){
    dir.create(paste0(savefilepath, "raw_counts"))  
  }
}



########OTHER SRA FUNCTION
getSRAproj <- function(project, cellosaurusfilepath, savefilepath){
  #check for rse, tpm, and raw_counts folders
  checkfolders(savefilepath)
  
  proj_home <- "data_sources/sra" #project home for recount3
  
  #assign project # by input
  proj <- ifelse(project == "ccle", "SRP186687", 
                 ifelse(project == "hpa_c", "SRP017465", 
                        ifelse(project == "hpa_t", "ERP003613", 
                               ifelse(project == "pdx", "SRP201347", 
                                      ifelse(project %in% c("bone_t", "bone_n", "bone_c"), "SRP090849", stop("project was not one of the available options: ccle, hpa_c, hpa_t, pdx, bone_t, bone_n, bone_c")
                                             
                                      )
                               )
                        )
                 )
  )
  
  #print(paste(proj))
  #create recount3 rse manuals
  rse <- recount3::create_rse_manual(
    project = proj,
    project_home = proj_home,
    organism = "human",
    annotation = "gencode_v26",
    type = "gene"
  )
  #most SRA-derived projects have attributes condensed -- need to be expanded
  if(project != "hpa_c"){
    rse <- expand_sra_attributes(rse)
  }
  
  # quick fix of HPA metadata mislabeling
  if (project == "hpa_c"){
    rse@colData@listData[["sra.sample_name"]][19] <- "PC3_b"
    rse@colData@listData[["sra.sample_name"]][20] <- "RT-4_a"
  }
  
  if (project == "ccle"){
    # missing most of its labeling metadata, including cell line of origin
    problematic_celllines <- read_csv(cellosaurusfilepath, col_names = FALSE)
    rse <- rse[,!rse@colData@listData[["sra_attribute.cell_line"]] %in% problematic_celllines$X2]
    # need to remove some parentheses in some that cause problems down the line
    rse[,grep("SRR8615298", colnames(rse))]@colData@listData[["sra_attribute.cell_line"]] <- "PE/CA-PJ41 clone D2"
    rse[,grep("SRR8615299", colnames(rse))]@colData@listData[["sra_attribute.cell_line"]] <- "PE/CA-PJ34 clone C12"
    rse[,grep("SRR8615493", colnames(rse))]@colData@listData[["sra_attribute.cell_line"]] <- "Ishikawa Heraklio 02 ER-"
    rse[,grep("SRR8615501", colnames(rse))]@colData@listData[["sra_attribute.cell_line"]] <- "SK-N-BE2"
    rse[,grep("SRR8615785", colnames(rse))]@colData@listData[["sra_attribute.cell_line"]] <- "Hs 688A.T"
    # missing most of its labeling metadata, including cell line of origin
    rse <- rse[,grep("SRR8615727", colnames(rse), invert = TRUE)]
  }
  
  if(project == "bone_t"){
    rse <- rse[,grep("osteosarcoma", rse@colData@listData[["sra_attribute.source_name"]])]
  } 
  if(project == "bone_n"){
    rse <- rse[,grep("normal bone", rse@colData@listData[["sra_attribute.source_name"]])]
  } 
  if(project == "bone_c"){
    rse <- rse[,grep("osteosarcoma|normal bone", rse@colData@listData[["sra_attribute.source_name"]], invert = TRUE)]
  }
  
  if (project == "pdx"){
    rse <- rse[,grep("IDH-mutant", rse@colData@listData[["sra_attribute.isolate"]], invert = TRUE)]
  }
  
  #calculate TPM
  assay(rse, "counts") <- recount3::transform_counts(rse)
  assays(rse)$TPM <- recount::getTPM(rse, length_var = "bp_length")
  
  #save RSE 
  saveRDS(rse, paste(savefilepath, "rse/", project, "_rse.rds", sep = ""))
  
  #save TPM and raw counts
  saveRDS(rse@assays@data@listData[["raw_counts"]], paste(savefilepath, "raw_counts/", project, "_rawcounts.rds", sep = ""))
  saveRDS(rse@assays@data@listData[["TPM"]], paste(savefilepath, "tpm/", project, "_tpm.rds", sep = ""))
  
  return(rse)
}




########TCGA FUNCTION
getTCGArse <- function(project, savefilepath){
  #project- a character string with the project name
  #check for rse, tpm, and raw_counts folders
  checkfolders(savefilepath)
  
  #in recount3, TCGA is split up by tissue for projects
  #proj <- c("BRCA","KIRC","LUAD","UCEC","THCA","PRAD","LUSC","HNSC","COAD","LGG","SKCM","LAML", "STAD","BLCA","OV","LIHC",
  #          "KIRP","CESC","SARC","ESCA","PCPG","PAAD","READ","GBM","TGCT","THYM","KICH","MESO","UVM","ACC","UCS","DLBC","CHOL")

  rse <- recount3::create_rse_manual(
    project = project,
    project_home = "data_sources/tcga",
    organism = "human",
    annotation = "gencode_v26",
    type = "gene"
  )
  
    #require(recount3)
    rse <- rse[,grep("Primary", rse@colData@listData[["tcga.cgc_sample_sample_type"]])]
    assay(rse, "counts") <- recount3::transform_counts(rse)
    assays(rse)$TPM <- recount::getTPM(rse, length_var = "bp_length")
  
  saveRDS(rse, paste(savefilepath, "rse/tcga_tissue/", paste(project), "_tcga_rse.rds", sep = ""))
  #return(rse)
  message(paste("Finished RSE pull and save for", project))
  # }
  
}


########GTEX FUNCTION
getGTEXrse <- function(project, savefilepath){
  #project- a character string with the project name
  #check for rse, tpm, and raw_counts folders
  checkfolders(savefilepath)
  
  #in recount3, GTEx is split up by tissue for projects
  #projects <- c("BRAIN", "SKIN", "ESOPHAGUS", "BLOOD", "BLOOD_VESSEL", "ADIPOSE_TISSUE", "HEART", "MUSCLE", "COLON", 
  #          "THYROID", "NERVE", "LUNG", "BREAST", "TESTIS", "STOMACH", "PANCREAS", "PITUITARY", "ADRENAL_GLAND", 
  #          "PROSTATE", "SPLEEN", "LIVER", "OVARY", "SMALL_INTESTINE", "SALIVARY_GLAND", "VAGINA", "UTERUS", 
  #          "KIDNEY", "BLADDER", "CERVIX_UTERI", "FALLOPIAN_TUBE")
  
  rse <- recount3::create_rse_manual(
    project = project,
    project_home = "data_sources/gtex",
    organism = "human",
    annotation = "gencode_v26",
    type = "gene"
  )
      
    #require(recount3)
    rse <- rse[,grep("RNASEQ", rse@colData@listData$gtex.smafrze)]
    assay(rse, "counts") <- recount3::transform_counts(rse)
    assays(rse)$TPM <- recount::getTPM(rse, length_var = "bp_length")
  
  saveRDS(rse, paste(savefilepath, "rse/gtex_tissue/", paste(project), "_gtex_rse.rds", sep = ""))
  #return(rse)
  message(paste("Finished RSE pull and save for", project))
  #}
}


########TPM/RAW_COUNTS FUNCTION
getRawCounts <- function(project, savefilepath){
  rsefiles <- Sys.glob(paste0(savefilepath, "*_", project, "_rse.rds")) #, paste(proj[i]), "_gtex_rse.rds", sep = "")
  #print(rsefiles)
  
  foreach(i = 1:length(rsefiles), .combine=cbind) %dopar% { #run loops in parallel, outputs a combined list
    rse <- readRDS(rsefiles[i])
    raw <-  rse@assays@data@listData[["raw_counts"]]
    
  }
  
}


getTPM <- function(project, savefilepath){
  rsefiles <- Sys.glob(paste0(savefilepath, "*_", project, "_rse.rds")) 
  #print(rsefiles)
  
  foreach(i = 1:length(rsefiles), .combine=cbind) %dopar% { #run loops in parallel, outputs a combined list
    rse <- readRDS(rsefiles[i])
    tpm <-  rse@assays@data@listData[["TPM"]]
    
  }
  
}
