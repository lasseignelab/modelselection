# correlation function to correlate group a to b, append b (so this should be tissue) by median, then find which tissue best correlates with each of group a and if that matches origin tissue
# this function is reliant on the tpm being available in the filepath format they would be saved in if rse_to_subset function is used

tissuecormatch <- function(rootfilepath, row_sources = c("bone_c", "ccle", "hpa_c", "pdx"), col_sources = c("tcga", "gtex", "hpa_t", "bone_n", "bone_t"), sub, pur = c(NA, "purinc", "purexc"), var = c(NA, "tcga", "gtex"), app = TRUE, app_by = c(NA, "tissue", "disease"), save = FALSE, savefilepath){
  # rootfilepath: root location where all tpm files are located. requires them to be in the filepath system established by rse_to_subset(), e.g., if filepath to row_sources tpm is ./data/rse_subsets/100/purinc/tcgavar/tpm/gtex_rse_sub100_purinc_tcgavar.rds, this arg should be "./data/rse_subsets/"
  # row_sources: all sources wanted for row
  # col_sources: all sources wanted for col
  # sub: 100, 1000, 5000, or 63856 - which gene subset to use
  # pur: "purinc" or "purexc" - whether to include purity
  # var: "tcga", "gtex" - which to measure variance by
  # app: "tissue" or "disease", whether to average columns across tissues or across diseases (so if control and tumor are both included, disease would separate them while tissue wouldn't)
  # savefilepath: where to generate /corres/ folder. can be the same as rootfilepath and will place /corres in same root as /rse and /tpm
  
  results <- list()
  if(sub == 63856){
    if(length(row_sources) == 1){
      row_tpm <- readRDS(paste(rootfilepath, "63856/tpm/", row_sources, "_tpm.rds", sep = ""))
    }else if(length(row_sources) > 1){
      row_1 <- row_sources[1]
      row_tpm <- readRDS(paste(rootfilepath, "63856/tpm/", row_1, "_tpm.rds", sep = ""))
      for(i in row_sources[-1]){
        row_tpm_temp <- readRDS(paste(rootfilepath, "63856/tpm/", i, "_tpm.rds", sep = ""))
        row_tpm <- cbind(row_tpm, row_tpm_temp)
      }
    }
    # make TPM object including all col sources specified
    if(length(col_sources) == 1){
      col_tpm <- readRDS(paste(rootfilepath, "63856/tpm/", col_sources, "_tpm.rds", sep = ""))
    }else if(length(col_sources) > 1){
      col_1 <- col_sources[1]
      col_tpm <- readRDS(paste(rootfilepath, "63856/tpm/", col_1, "_tpm.rds", sep = ""))
      for(i in col_sources[-1]){
        col_tpm_temp <- readRDS(paste(rootfilepath, "63856/tpm/", i, "_tpm.rds", sep = ""))
        col_tpm <- cbind(col_tpm, col_tpm_temp)
      }
    }
  } else {
    # make TPM object including all row sources specified
    if(length(row_sources) == 1){
      row_tpm <- readRDS(paste(rootfilepath, sub, "/", pur, "/", var, "var/tpm/", row_sources, "_tpm_sub", sub, "_", pur, "_", var, "var.rds", sep = ""))
    }else if(length(row_sources) > 1){
      row_1 <- row_sources[1]
      row_tpm <- readRDS(paste(rootfilepath, sub, "/", pur, "/", var, "var/tpm/", row_1, "_tpm_sub", sub, "_", pur, "_", var, "var.rds", sep = ""))
      for(i in row_sources[-1]){
        row_tpm_temp <- readRDS(paste(rootfilepath, sub, "/", pur, "/", var, "var/tpm/", i, "_tpm_sub", sub, "_", pur, "_", var, "var.rds", sep = ""))
        row_tpm <- cbind(row_tpm, row_tpm_temp)
      }
    }
    # make TPM object including all col sources specified
    if(length(col_sources) == 1){
      col_tpm <- readRDS(paste(rootfilepath, sub, "/", pur, "/", var, "var/tpm/", col_sources, "_tpm_sub", sub, "_", pur, "_", var, "var.rds", sep = ""))
    }else if(length(col_sources) > 1){
      col_1 <- col_sources[1]
      col_tpm <- readRDS(paste(rootfilepath, sub, "/", pur, "/", var, "var/tpm/", col_1, "_tpm_sub", sub, "_", pur, "_", var, "var.rds", sep = ""))
      for(i in col_sources[-1]){
        col_tpm_temp <- readRDS(paste(rootfilepath, sub, "/", pur, "/", var, "var/tpm/", i, "_tpm_sub", sub, "_", pur, "_", var, "var.rds", sep = ""))
        col_tpm <- cbind(col_tpm, col_tpm_temp)
      }
    }
  } 
  # correlate rows and columns
  cor <- cor(row_tpm, col_tpm, method = "spearman")
  
  # append columns if app == TRUE
  if(app == TRUE){
    if(app_by == "tissue"){
      groups <- unique(sample_master_meta[sample_master_meta$rse_sampname %in% colnames(cor),]$tissue_coerced)
      i <- groups[1]
      samples <- sample_master_meta[grep(i, sample_master_meta$tissue_coerced),]$rse_sampname
      cor_appcol <- data.frame(matrixStats::rowMedians(cor[,colnames(cor) %in% samples]))
      
      for(i in groups[-1]){
        samples <- sample_master_meta[grep(i, sample_master_meta$tissue_coerced),]$rse_sampname
        col_2 <- data.frame(matrixStats::rowMedians(cor[,colnames(cor) %in% samples]))
        cor_appcol <- cbind(cor_appcol, col_2)
      }
    }else if(app_by == "disease"){
      groups <- unique(sample_master_meta[sample_master_meta$rse_sampname %in% colnames(cor),]$disease_coerced)
      i <- groups[1]
      samples <- sample_master_meta[grep(i, sample_master_meta$disease_coerced),]$rse_sampname
      cor_appcol <- data.frame(matrixStats::rowMedians(cor[,colnames(cor) %in% samples]))
      
      for(i in groups[-1]){
        samples <- sample_master_meta[grep(i, sample_master_meta$disease_coerced),]$rse_sampname
        col_2 <- data.frame(matrixStats::rowMedians(cor[,colnames(cor) %in% samples]))
        cor_appcol <- cbind(cor_appcol, col_2)
      } 
    }
    colnames(cor_appcol) <- groups
    rownames(cor_appcol) <- rownames(cor)
    results[["fullcorapp"]] <- cor_appcol
    
    # now want to return how many of the rows matched/didn't with their origin tissue
    res_top_first <- colnames(cor_appcol)[max.col(cor_appcol, ties.method = "first")]
    res_top_last <- colnames(cor_appcol)[max.col(cor_appcol, ties.method = "last")]
    if(all.equal(res_top_first, res_top_last) == FALSE){
      warn("ties in max correlation, only first is displayed")
    }
    
    if(app_by == "tissue"){
      res_df <- data.frame(sample = rownames(cor), samp_origin = sample_master_meta[sample_master_meta$rse_sampname %in% rownames(cor),]$tissue_coerced)
    } else if(app_by == "disease"){
      res_df <- data.frame(sample = rownames(cor), samp_origin = sample_master_meta[sample_master_meta$rse_sampname %in% rownames(cor),]$disease_coerced)
    }
    
    res_df <- cbind(res_df, res_top_first)
    
    results[["alltoptissue"]] <- res_df
    
    res_df_sub <- res_df[res_df$samp_origin != res_df$res_top_first,]
    res_df_temp <- res_df[res_df$samp_origin == res_df$res_top_first,]
    
    results[["dif_from_origin"]] <- res_df_sub
    results[["dfo_summary"]] <- data.frame(same_as_origin = length(res_df_temp$sample), different = length(res_df_sub$sample))
    
    print(data.frame(same_as_origin = length(res_df_temp$sample), different = length(res_df_sub$sample)))
    
  } else if(app == FALSE){
    results[["fullcor"]] <- cor
  }
  
  # save
  if(save == TRUE){
    if(app_by == "disease"){
      if(sub != 63856){
        saveRDS(results, paste(savefilepath, sub, "/", pur, "/", var, "var/corres/", "cor_dismatch_", pur, "_", var, "var_res_", sub, ".rds", sep = ""))
      }else if(sub == 63856){
        if("pdx" %in% row_sources){
          saveRDS(results, paste(savefilepath, "63856/corres/", "cor_dismatch_res_all.rds", sep = ""))
        }else{
          saveRDS(results, paste(savefilepath, "63856/corres/", "cor_dismatch_res_allCL.rds", sep = ""))
        }
      }
    }else if(app_by == "tissue"){
      if(sub != 63856){
        saveRDS(results, paste(savefilepath, sub, "/", pur, "/", var, "var/corres/", "cor_tismatch_", pur, "_", var, "var_res_", sub, ".rds", sep = ""))
      }else if(sub == 63856){
        if("pdx" %in% row_sources){
          saveRDS(results, paste(savefilepath, "63856/corres/", "cor_tismatch_res_all.rds", sep = ""))
        }else{
          saveRDS(results, paste(savefilepath, "63856/corres/", "cor_tismatch_res_allCL.rds", sep = ""))
        }
      }
    }
  }
  
  return(results)
}