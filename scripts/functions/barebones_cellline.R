barebones_cellline <- function(x){
  # input a list of cell lines
  # will remove any variation based on dashes, spaces, and case differences, check for underscores etc.
  # pipe operator had issues with gsub()
  x <- gsub(" *-+ *", "", x)
  x <- gsub(" ", "", x)
  x <- gsub("[.]", "", x)
  x <- tolower(x)
  
  return(x)
}