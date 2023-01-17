# correlate based on tissue-specific gene set
tisgenesetcor <- function(x,y,method = "spearman", tissue, onlymatch = TRUE, medtis = TRUE, medmod = FALSE){
  # x = matrix of tpm for preclinical models to cor
  # y = matrix of tpm for tissues to cor
  # method = cor method
  # tissue = tissue to cor by tissue-specific genes, must be an option from gtex_tissues_fin obj, LOWERCASE
  # onlymatch = if TRUE, will only correlate matching preclinical models (x) to specified tissue in TAGs,
  # if FALSE, will correlate all preclinical models (x) to specified tissue in TAGs
  # medtis = if TRUE, will give the median performance of each model (x) across specified tissue in TAGs
  # if FALSE, will output the full correlation matrix of all x to all y (y only in specified tissue)
  # medmod = if TRUE, will take the median of the performance of all models in a given tissue's median correlation to the tissue; thus returns a single value instead of a df matching the length of the model set
  # NOTE, this arg is written to be stacked on top of medtis = TRUE, thus we can only summarize across a tissue and maintain individual models, or summarize across both
  
  # first, pull geneset of interest from gtex_both_tissue_genes
  genes <- gtex_both_tissue_genes[[paste(tissue)]]
  
  # subset tpm matrices to just gene set
  x_sub <- x[substr(rownames(x),1,15) %in% genes,]
  y_sub <- y[substr(rownames(y),1,15) %in% genes,]
  
  # select only x and y samples that match tissue
  if(onlymatch == TRUE){
    if(tissue != "intestine"){
      # first, subset to just samples in each tpm matrix
      x_meta <- sample_master_meta[sample_master_meta$rse_sampname %in% colnames(x_sub),]
      # next, subset by tissue
      x_meta <- x_meta[x_meta$tissue_coerced == stringr::str_to_title(tissue),]
      x_samples <- x_meta$rse_sampname
      
      # first, subset to just samples in each tpm matrix
      y_meta <- sample_master_meta[sample_master_meta$rse_sampname %in% colnames(y_sub),]
      # neyt, subset by tissue
      y_meta <- y_meta[y_meta$tissue_coerced == stringr::str_to_title(tissue),]
      y_samples <- y_meta$rse_sampname
    } else if(tissue == "intestine"){
      # needs to be different because we need both "colon" and "small intestine"
      # first, subset to just samples in each tpm matrix
      x_meta <- sample_master_meta[sample_master_meta$rse_sampname %in% colnames(x_sub),]
      # next, subset by tissue
      x_meta_suba <- x_meta[x_meta$tissue_coerced == "Colon",]
      x_meta_subb <- x_meta[x_meta$tissue_coerced == "Small Intestine",]
      x_samples <- c(x_meta_suba$rse_sampname, x_meta_subb$rse_sampname)
      
      # first, subset to just samples in each tpm matrix
      y_meta <- sample_master_meta[sample_master_meta$rse_sampname %in% colnames(y_sub),]
      # neyt, subset by tissue
      y_meta_suba <- y_meta[y_meta$tissue_coerced == "Colon",]
      y_meta_subb <- y_meta[y_meta$tissue_coerced == "Small Intestine",]
      y_samples <- c(y_meta_suba$rse_sampname, y_meta_subb$rse_sampname)
    }
    
    # now, subset counts matrices to just these samples
    x_sub <- x_sub[,colnames(x_sub) %in% x_samples]
    y_sub <- y_sub[,colnames(y_sub) %in% y_samples]
    
  } else{
    # just want to correlate to gene set in specified tissue still, so just subset y
    if(tissue != "intestine"){
      # first, subset to just samples in each tpm matrix
      y_meta <- sample_master_meta[sample_master_meta$rse_sampname %in% colnames(y_sub),]
      # neyt, subset by tissue
      y_meta <- y_meta[y_meta$tissue_coerced == stringr::str_to_title(tissue),]
      y_samples <- y_meta$rse_sampname
    } else if(tissue == "intestine"){
      # first, subset to just samples in each tpm matrix
      y_meta <- sample_master_meta[sample_master_meta$rse_sampname %in% colnames(y_sub),]
      # next, subset by tissue
      y_meta_suba <- y_meta[y_meta$tissue_coerced == "Colon",]
      y_meta_subb <- y_meta[y_meta$tissue_coerced == "Small Intestine",]
      y_samples <- c(y_meta_suba$rse_sampname, y_meta_subb$rse_sampname)
    }
    
    # now, subset y counts matrix to just these samples
    y_sub <- y_sub[,colnames(y_sub) %in% y_samples]
    
  }
  
  # correlate
  cor <- cor(x_sub, y_sub, method = method)
  
  if(medtis == TRUE){
    cor <- data.frame(cor = MatrixGenerics::rowMedians(cor))
    rownames(cor) <- colnames(x_sub)
  }
  
  if(medmod == TRUE){
    cor <- median(cor)
  }
  
  # return obj
  return(cor)
}