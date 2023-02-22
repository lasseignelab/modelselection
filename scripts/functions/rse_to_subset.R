#### function to create needed folders
subsetfolders <- function(savefilepath, subset_num){
  if(!dir.exists(paste0(savefilepath, "rse_subsets/")) == TRUE){
    dir.create(paste0(savefilepath, "rse_subsets/"))  
  }
  for(i in 1:length(subset_num)){
    if(!dir.exists(paste0(savefilepath, "rse_subsets/", subset_num[i])) == TRUE){
      dir.create(paste0(savefilepath, "rse_subsets/", subset_num[i]))  
    }
    purity <- c("purinc", "purexc")
    for(j in 1:length(purity)){
      if(!dir.exists(paste0(savefilepath, "rse_subsets/", subset_num[i], "/", purity[j])) == TRUE){
        dir.create(paste0(savefilepath, "rse_subsets/", subset_num[i], "/", purity[j]))
      }
      vartypes <- c("tcgavar", "gtexvar")
      for(k in 1:length(vartypes)){
        if(!dir.exists(paste0(savefilepath, "rse_subsets/", subset_num[i], "/", purity[j], "/", vartypes[k])) == TRUE){
          dir.create(paste0(savefilepath, "rse_subsets/", subset_num[i], "/", purity[j], "/", vartypes[k]))
        }
        filetype <- c("rse", "tpm")
        for(l in 1:length(filetype)){
          if(!dir.exists(paste0(savefilepath, "rse_subsets/", subset_num[i], "/", purity[j], "/", vartypes[k], "/", filetype[l])) == TRUE){
            dir.create(paste0(savefilepath, "rse_subsets/", subset_num[i], "/", purity[j], "/", vartypes[k],  "/", filetype[l]))
          }
        }
      }
    }
  }
}

rse_to_subset <- function(rse, project, subset_num, var_source, purity_inc, saveRSE, saveTPM, savefilepath){
  # rse is the full input rse to be subsetted
  # project is a character string of how to name the new object
  # subset_num is a numeric of how many of the most variable genes to subset to
  # var_source is either "tcga" or "gtex", whichever we base our variation on
  # purity_exc is TRUE/FALSE, whether to remove purity-related genes from comparison
  # NOTE: should ONLY set purity_exc to FALSE if var_source is TCGA
  # save is TRUE/FALSE, whether to save resultant object to RDS
  
  #check if necessary folders exist, and create them if they don't
  subsetfolders(savefilepath, subset_num)
  
  if(var_source == "tcga"){
    varying_genes <- tcga_varying_genes
  } else if(var_source == "gtex"){
    varying_genes <- gtex_varying_genes
  }
  
  if(purity_inc == TRUE){
    varying_genes <- names(varying_genes[1:subset_num])
    
    output <- rse[rownames(rse) %in% varying_genes,]
    
    
    if(saveRSE == TRUE){
      if(var_source == "tcga"){
        mainDir <- paste(savefilepath, subset_num, "/purinc/tcgavar", sep = "")
          ifelse(!dir.exists(file.path(mainDir, "rse")), dir.create(file.path(mainDir, "rse")), FALSE)
          saveRDS(rse[rownames(rse) %in% varying_genes,], paste(mainDir, "rse/", project, "_rse", "_sub", subset_num, "_purinc_tcgavar.rds", sep = ""))
      } else if(var_source == "gtex"){
        mainDir <- paste(savefilepath, subset_num, "/purinc/gtexvar", sep = "")
          ifelse(!dir.exists(file.path(mainDir, "rse")), dir.create(file.path(mainDir, "rse")), FALSE)
          saveRDS(rse[rownames(rse) %in% varying_genes,], paste(mainDir, "rse/", project, "_rse", "_sub", subset_num, "_purinc_gtexvar.rds", sep = ""))
        
        #saveRDS(rse[rownames(rse) %in% varying_genes,], paste(savefilepath, subset_num, "/purinc/gtexvar/rse/", project, "_rse", "_sub", subset_num, "_purinc_gtexvar.rds", sep = ""))
      }
    }
    if(saveTPM == TRUE){
      if(var_source == "tcga"){
        mainDir <- paste(savefilepath, subset_num, "/purinc/tcgavar", sep = "")
        ifelse(!dir.exists(file.path(mainDir, "tpm")), dir.create(file.path(mainDir, "tpm")), FALSE)
        saveRDS(rse[rownames(rse) %in% varying_genes,], paste(mainDir, "tpm/", project, "_tpm", "_sub", subset_num, "_purinc_tcgavar.rds", sep = ""))
        
#        saveRDS(rse[rownames(rse) %in% varying_genes,]@assays@data@listData[["TPM"]], paste(savefilepath, subset_num, "/purinc/tcgavar/tpm/", project, "_tpm", "_sub", subset_num, "_purinc_tcgavar.rds", sep = ""))
      } else if(var_source == "gtex"){
        mainDir <- paste(savefilepath, subset_num, "/purinc/gtexvar", sep = "")
        ifelse(!dir.exists(file.path(mainDir, "tpm")), dir.create(file.path(mainDir, "tpm")), FALSE)
        saveRDS(rse[rownames(rse) %in% varying_genes,], paste(mainDir, "tpm/", project, "_tpm", "_sub", subset_num, "_purinc_gtexvar.rds", sep = ""))
#        saveRDS(rse[rownames(rse) %in% varying_genes,]@assays@data@listData[["TPM"]], paste(savefilepath, subset_num, "/purinc/gtexvar/tpm/", project, "_tpm", "_sub", subset_num, "_purinc_gtexvar.rds", sep = ""))
      }
    }
    
  } else if(purity_inc == FALSE){
    varying_genes <- names(varying_genes)
    varying_genes <- varying_genes[!varying_genes %in% purity_genes]
    varying_genes <- varying_genes[1:subset_num]
    
    output <-  rse[rownames(rse) %in% varying_genes,]
    
    if(saveRSE == TRUE){
      if(var_source == "tcga"){
        mainDir <- paste(savefilepath, subset_num, "/purexc/tcgavar", sep = "")
        ifelse(!dir.exists(file.path(mainDir, "rse")), dir.create(file.path(mainDir, "rse")), FALSE)
        saveRDS(rse[rownames(rse) %in% varying_genes,], paste(mainDir, "rse/", project, "_rse", "_sub", subset_num, "_purexc_tcgavar.rds", sep = ""))
        
        #saveRDS(rse[rownames(rse) %in% varying_genes,], paste(savefilepath, subset_num, "/purexc/tcgavar/rse/", project, "_rse", "_sub", subset_num, "_purexc_tcgavar.rds", sep = ""))
      } else if(var_source == "gtex"){
        mainDir <- paste(savefilepath, subset_num, "/purexc/gtexvar", sep = "")
        ifelse(!dir.exists(file.path(mainDir, "rse")), dir.create(file.path(mainDir, "rse")), FALSE)
        saveRDS(rse[rownames(rse) %in% varying_genes,], paste(mainDir, "rse/", project, "_rse", "_sub", subset_num, "_purexc_gtexvar.rds", sep = ""))
        
        #saveRDS(rse[rownames(rse) %in% varying_genes,], paste(savefilepath, subset_num, "/purexc/gtexvar/rse/", project, "_rse", "_sub", subset_num, "_purexc_gtexvar.rds", sep = ""))
      }
    }
    if(saveTPM == TRUE){
      if(var_source == "tcga"){
        mainDir <- paste(savefilepath, subset_num, "/purexc/tcgavar", sep = "")
        ifelse(!dir.exists(file.path(mainDir, "tpm")), dir.create(file.path(mainDir, "tpm")), FALSE)
        saveRDS(rse[rownames(rse) %in% varying_genes,], paste(mainDir, "tpm/", project, "_tpm", "_sub", subset_num, "_purexc_tcgavar.rds", sep = ""))
        
        #saveRDS(rse[rownames(rse) %in% varying_genes,]@assays@data@listData[["TPM"]], paste(savefilepath, subset_num, "/purexc/tcgavar/tpm/", project, "_tpm", "_sub", subset_num, "_purexc_tcgavar.rds", sep = ""))
      } else if(var_source == "gtex"){
        mainDir <- paste(savefilepath, subset_num, "/purexc/gtexvar", sep = "")
        ifelse(!dir.exists(file.path(mainDir, "tpm")), dir.create(file.path(mainDir, "tpm")), FALSE)
        saveRDS(rse[rownames(rse) %in% varying_genes,], paste(mainDir, "tpm/", project, "_tpm", "_sub", subset_num, "_purexc_gtexvar.rds", sep = ""))
        
        #saveRDS(rse[rownames(rse) %in% varying_genes,]@assays@data@listData[["TPM"]], paste(savefilepath, subset_num, "/purexc/gtexvar/tpm/", project, "_tpm", "_sub", subset_num, "_purexc_gtexvar.rds", sep = ""))
      }
    }
  }
  
  return(output)
}