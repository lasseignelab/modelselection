# after using tissuecormatch() to make correlation matrices
tissuecormatch_median <- function(filepath, subsize = c(100, 1000, 5000, 10000, 63856), pur = c("purinc", "purexc", NA), var = c("tcga", "gtex", NA)){
  # first, pull the corresponding correlation table
  if(subsize == 63856){
    cor_temp <- readRDS(paste(filepath, "63856/corres/cor_dismatch_res_allCL.rds", sep = ""))
  }else{
    cor_temp <- readRDS(paste(filepath, subsize, "/", pur, "/", var, "var/corres/cor_dismatch_", pur, "_", var, "var_res_", subsize, ".rds", sep = ""))
  }  
  fullcor <- cor_temp$fullcorapp # NOTE: requires corres obj to have $fullcorapp
  origintissue <- cor_temp$alltoptissue[,1:2]
  
  # now, find median correlation across each individual tissue
  i <- unique(origintissue$samp_origin)[1]
  samples <- origintissue[grep(i, origintissue$samp_origin),]$sample
  fullcor_sub <- fullcor[rownames(fullcor) %in% samples,colnames(fullcor) %in% i]
  if(length(fullcor_sub != 0)){
    medcor <- median(fullcor_sub)
    tissues <- i
  }
  
  
  for(i in unique(origintissue$samp_origin)[-1]){
    samples <- origintissue[grep(i, origintissue$samp_origin),]$sample
    fullcor_sub <- fullcor[rownames(fullcor) %in% samples,colnames(fullcor) %in% i]
    if(length(fullcor_sub != 0)){
      medcor <- c(medcor, median(fullcor_sub))
      tissues <- c(tissues, i)
    }
  }
  results <- data.frame(median_cor = medcor)
  rownames(results) <- tissues
  return(results)
}