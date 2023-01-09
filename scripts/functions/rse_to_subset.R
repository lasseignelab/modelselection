rse_to_subset <- function(rse, project, subset_num, var_source, purity_inc, save, savefilepath){
  # rse is the full input rse to be subsetted
  # project is a character string of how to name the new object
  # subset_num is a numeric of how many of the most variable genes to subset to
  # var_source is either "tcga" or "gtex", whichever we base our variation on
  # purity_exc is TRUE/FALSE, whether to remove purity-related genes from comparison
  # NOTE: should ONLY set purity_exc to FALSE if var_source is TCGA
  # save is TRUE/FALSE, whether to save resultant object to RDS
  
  if(var_source == "tcga"){
    varying_genes <- tcga_varying_genes
  } else if(var_source == "gtex"){
    varying_genes <- gtex_varying_genes
  }
  
  if(purity_inc == TRUE){
    varying_genes <- names(varying_genes[1:subset_num])
    
    output <- rse[rownames(rse) %in% varying_genes,]
    
    if(save == TRUE){
      if(var_source == "tcga"){
        saveRDS(rse[rownames(rse) %in% varying_genes,], paste(savefilepath, subset_num, "/purinc/tcgavar/rse/", project, "_rse", "_sub", subset_num, "_purinc_tcgavar.rds", sep = ""))
      } else if(var_source == "gtex"){
        saveRDS(rse[rownames(rse) %in% varying_genes,], paste(savefilepath, subset_num, "/purinc/gtexvar/rse/", project, "_rse", "_sub", subset_num, "_purinc_gtexvar.rds", sep = ""))
      }
    }
    
  } else if(purity_inc == FALSE){
    varying_genes <- names(varying_genes)
    varying_genes <- varying_genes[!varying_genes %in% purity_genes]
    varying_genes <- varying_genes[1:subset_num]
    
    output <-  rse[rownames(rse) %in% varying_genes,]
    
    if(save == TRUE){
      if(var_source == "tcga"){
        saveRDS(rse[rownames(rse) %in% varying_genes,], paste(savefilepath, subset_num, "/purexc/tcgavar/rse/", project, "_rse", "_sub", subset_num, "_purexc_tcgavar.rds", sep = ""))
      } else if(var_source == "gtex"){
        saveRDS(rse[rownames(rse) %in% varying_genes,], paste(savefilepath, subset_num, "/purexc/gtexvar/rse/", project, "_rse", "_sub", subset_num, "_purexc_gtexvar.rds", sep = ""))
      }
    }  
  }
  
  return(output)
}