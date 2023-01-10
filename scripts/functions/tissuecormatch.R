# correlation function to correlate group a to b, append b (so this should be tissue) by median, then find which tissue best correlates with each of group a and if that matches origin tissue

tissuecormatch <- function(row_tpm, col_tpm, sub, pur = c("purinc", "purexc", NA), var = c("tcga", "gtex", NA), app = TRUE, app_by = c("tissue", "disease", NA), save = FALSE, savefilepath){
  # row_sources: all sources wanted for row
  # col_sources: all sources wanted for col
  # sub: 100, 1000, 5000, or 63856 - which gene subset to use
  # pur: "purinc" or "purexc" - whether to include purity
  # var: "tcga", "gtex" - which to measure variance by
  # app: "tissue" or "disease", whether to average columns across tissues or across diseases (so if control and tumor are both included, disease would separate them while tissue wouldn't)
  
  results <- list()
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
      saveRDS(results, paste(savefilepath, sub, "/", pur, "/", var, "var/corres/", "cor_dismatch_", pur, "_", var, "var_res_", sub, ".rds", sep = ""))
    }else if(app_by == "tissue"){
      saveRDS(results, paste(savefilepath, sub, "/", pur, "/", var, "var/corres/", "cor_tismatch_", pur, "_", var, "var_res_", sub, ".rds", sep = ""))
    }
  }
  
  return(results)
}