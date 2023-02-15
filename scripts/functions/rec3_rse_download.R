
getgtex <- function(project = "gtex", savefilepath){
  rse <- getgtexrse()
  saverse(rse, savefilepath)
}

saveraw <- function(project = "gtex", savefilepath){
  rsefiles <- Sys.glob(paste0(savefilepath, "rse/*_gtex_rse.rds")) #, paste(proj[i]), "_gtex_rse.rds", sep = "")
  print(rsefiles)
  #rawcounts <- describe(rsefiles)
  foreach(i = 1:length(rsefiles), .combine=cbind, .packages = "bigmemory") %dopar% {
    rse <- readRDS(rsefiles[i])
    raw <-  rse@assays@data@listData[["raw_counts"]]
    #write.csv(raw, file = paste(savefilepath, "raw_counts/gtex_rawcounts.csv", sep = ""))
    #print(class(raw))
   # rawcounts <- attach.big.matrix(raw)
    #tpm <- rse@assays@data@listData[["TPM"]]
#    print(str(raw))
#    print(str(tpm))
#    return(rsefiles)
  
#  saveRDS(tpm, paste(savefilepath, "tpm/", project, "_tpm.rds", sep = ""))
  }
  # now, save RSE as object
#    saveRDS(rse, paste(savefilepath, "rse/", project, "_rse.rds", sep = ""))
#    saveRDS(rse@assays@data@listData[["raw_counts"]], paste(savefilepath, "raw_counts/", project, "_rawcounts.rds", sep = ""))
  
  # also save raw counts and TPM individually
#    saveRDS(rse@assays@data@listData[["TPM"]], paste(savefilepath, "tpm/", project, "_tpm.rds", sep = ""))
#  saveRDS(raw, paste(savefilepath, "raw_counts/", project, "_rawcounts.rds", sep = ""))
  
  #try to make into a bigmemory matrix?
  #as.big.matrix(rawcounts)
  
  #print(str(rawcounts))
#  return(bigmemory::big.matrix(rawcounts))
  #return(rawcounts)
}



savetpm <- function(project = "gtex", savefilepath){
  rsefiles <- Sys.glob(paste0(savefilepath, "rse/*_gtex_rse.rds")) 
  print(rsefiles)

  foreach(i = 1:length(rsefiles), .combine=cbind, .packages = "bigmemory") %dopar% {
    rse <- readRDS(rsefiles[i])
    tpm <-  rse@assays@data@listData[["TPM"]]

  }

}

getgtexrse <- function(project = "gtex", savefilepath){
  proj_home <- "data_sources/gtex"
  proj <- c("BRAIN", "SKIN", "ESOPHAGUS", "BLOOD", "BLOOD_VESSEL", "ADIPOSE_TISSUE", "HEART", "MUSCLE", "COLON", "THYROID", "NERVE", "LUNG", "BREAST", "TESTIS", "STOMACH", "PANCREAS", "PITUITARY", "ADRENAL_GLAND", "PROSTATE", "SPLEEN", "LIVER", "BONE_MARROW", "OVARY", "SMALL_INTESTINE", "SALIVARY_GLAND", "VAGINA", "UTERUS", "KIDNEY", "BLADDER", "CERVIX_UTERI", "FALLOPIAN_TUBE")
  #proj <- c("BLADDER", "CERVIX_UTERI", "FALLOPIAN_TUBE") #subset for testing

  foreach(i = 1:length(proj)) %dopar% {
    require(recount3)
    rse <- recount3::create_rse_manual(
      project = proj[i],
      project_home = proj_home,
      organism = "human",
      annotation = "gencode_v26",
      type = "gene"
      
    )
    require(recount3)
    assay(rse, "counts") <- recount3::transform_counts(rse)
    assays(rse)$TPM <- recount::getTPM(rse, length_var = "bp_length")
    #save RSE for each project
    saveRDS(rse, paste(savefilepath, "rse/", paste(proj[i]), "_gtex_rse.rds", sep = ""))
    #save csv of TPM and rawcounts for each
    #write.csv(rse@assays@data@listData[["TPM"]], paste(savefilepath, "tpm/", paste(proj[i]), "_gtex_tpm.csv", sep = ""))
    #write.csv(rse@assays@data@listData[["raw_counts"]], paste(savefilepath, "raw_counts/", paste(proj[i]), "_gtex_rawcounts.csv", sep = ""))

    return(rse)
  }
#return(rse)
}


rec3_rse_download <- function(project = c("tcga", "gtex", "ccle", "hpa_c", "hpa_t", "pdx", "bone_t", "bone_n", "bone_c"), cellosaurusfilepath, savefilepath){
  # creating a function to specify one of tcga, ccle, pdx, hpa, or gtex and download rse loaded with TPM
  # built to remove/rename specific problematic samples consistently
  # utilizes cellosaurus problematic cell line list in this directory's /data folder
  
  # determining project home
  if (project == "tcga"){
    proj_home <- "data_sources/tcga"
  } else if (project == "gtex"){
    proj_home <- "data_sources/gtex"
  } else {
    proj_home <- "data_sources/sra"
  }
  
  # setting proj object to pull rse
  if (project == "tcga"){
    proj <- c("BRCA","KIRC","LUAD","UCEC","THCA","PRAD","LUSC","HNSC","COAD","LGG","SKCM","LAML", "STAD","BLCA","OV","LIHC","KIRP","CESC","SARC","ESCA","PCPG","PAAD","READ","GBM","TGCT","THYM","KICH","MESO","UVM","ACC","UCS","DLBC","CHOL")
  } else if (project == "gtex"){
    proj <- c("BRAIN", "SKIN", "ESOPHAGUS", "BLOOD", "BLOOD_VESSEL", "ADIPOSE_TISSUE", "HEART", "MUSCLE", "COLON", "THYROID", "NERVE", "LUNG", "BREAST", "TESTIS", "STOMACH", "PANCREAS", "PITUITARY", "ADRENAL_GLAND", "PROSTATE", "SPLEEN", "LIVER", "BONE_MARROW", "OVARY", "SMALL_INTESTINE", "SALIVARY_GLAND", "VAGINA", "UTERUS", "KIDNEY", "BLADDER", "CERVIX_UTERI", "FALLOPIAN_TUBE")
  } else if (project == "ccle"){
    proj <- "SRP186687"
  } else if (project == "hpa_c"){
    proj <- "SRP017465"
  } else if (project == "hpa_t"){
    proj <- "ERP003613"
  } else if (project == "pdx"){
    proj <- "SRP201347"
  } else if (project %in% c("bone_t","bone_n","bone_c")){
    proj <- "SRP090849"
  }else {
    warning("project was not one of the options")
    stop()
  }
  
  # next, can set to pull rse
  if (length(proj) > 1){
    # make list of all rse's, then collapse together
    temp_rse_list <- list()
    for (i in proj){
      temp_rse_list[[paste(i)]] <- recount3::create_rse_manual(
        project = paste(i),
        project_home = proj_home,
        organism = "human",
        annotation = "gencode_v26",
        type = "gene"
      )
      if (project == "tcga"){
        temp_rse_list[[paste(i)]] <- temp_rse_list[[paste(i)]][,grep("Primary", temp_rse_list[[paste(i)]]@colData@listData[["tcga.cgc_sample_sample_type"]])]
      }
    }
    rse <- temp_rse_list[[1]]
    for(i in proj[-1]){
      rse <- cbind(rse, temp_rse_list[[paste(i)]])  
    }
  } else if (length(proj) == 1){
    rse <- recount3::create_rse_manual(
      project = proj,
      project_home = proj_home,
      organism = "human",
      annotation = "gencode_v26",
      type = "gene"
    )
  }
  
  if(proj_home == "data_sources/sra"){
    if(project != "hpa_c"){
      rse <- expand_sra_attributes(rse)
    }
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
  
  if (project == "gtex"){
    rse <- rse[,grep("RNASEQ", rse@colData@listData$gtex.smafrze)]
  }
  
  if (project == "pdx"){
    rse <- rse[,grep("IDH-mutant", rse@colData@listData[["sra_attribute.isolate"]], invert = TRUE)]
  }
  
  if (project %in% c("bone_t","bone_n","bone_c")){
    if(project == "bone_t"){
      rse <- rse[,grep("osteosarcoma", rse@colData@listData[["sra_attribute.source_name"]])]
    } else if(project == "bone_n"){
      rse <- rse[,grep("normal bone", rse@colData@listData[["sra_attribute.source_name"]])]
    } else if(project == "bone_c"){
      rse <- rse[,grep("osteosarcoma|normal bone", rse@colData@listData[["sra_attribute.source_name"]], invert = TRUE)]
    }
  }
  
  # should have assigned an rse object at this point, now need to generate tpm
  assay(rse, "counts") <- transform_counts(rse)
  assays(rse)$TPM <- recount::getTPM(rse, length_var = "bp_length")
  
  # now, save RSE as object
  saveRDS(rse, paste(savefilepath, "rse/", project, "_rse.rds", sep = ""))
  
  # also save raw counts and TPM individually
  saveRDS(rse@assays@data@listData[["raw_counts"]], paste(savefilepath, "raw_counts/", project, "_rawcounts.rds", sep = ""))
  saveRDS(rse@assays@data@listData[["TPM"]], paste(savefilepath, "tpm/", project, "_tpm.rds", sep = ""))
}