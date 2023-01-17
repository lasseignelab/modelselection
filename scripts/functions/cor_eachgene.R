# NOTE: will return rho = 0 if one or both vectors considered is entirely 0's
cor_eachgene <- function(x, y, method = "spearman"){
  # x and y are each tpm matrices, for my use x is model expression and y is tissue
  
  # first, make an empty dataframe with named rows/cols to fill in after
  ## need to make a loop to make a list of columns comparing every x sample to every y sample
  
  res <- data.frame(matrix(ncol = 2, nrow = length(rownames(x))))
  rownames(res) <- rownames(x)
  colnames(res) <- c("rho", "p.value")
  
  # fill in data frame with correlations
  for(i in rownames(x)){
    # operating under the assumption x and y have the same gene set, which they would with recount3
    if(sum(x[rownames(x) == i,]) == 0 || sum(y[rownames(y) == i,][1:length(na.omit(x[rownames(x) == i,]))]) == 0){
      # above because correlation will be NA if either or both is entirely 0's
      res[rownames(res) == i,1] <- 0
      res[rownames(res) == i,2] <- NA
    }else if(sum(x[rownames(x) == i,], y[rownames(y) == i,]) != 0){
      cortest <- cor.test(na.omit(x[rownames(x) == i,]), na.omit(y[rownames(y) == i,])[1:length(na.omit(x[rownames(x) == i,]))], alternative = "two.sided", method = method, exact = FALSE) # right now it's just picking the first n samples in y matching amount of x samples, assuming y has more samples than x
      res[rownames(res) == i,1] <- cortest[["estimate"]][["rho"]]
      res[rownames(res) == i,2] <- cortest[["p.value"]]
    }
  }
  
  return(res)
}