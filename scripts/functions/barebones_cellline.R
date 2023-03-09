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

#Read in all TCGA or GTEx subproject (tissue) RSE's and combine metadata tables
getMetadata <- function(savefilepath){
  #rsefiles <- Sys.glob(paste0(savefilepath, "_rse.rds"))
  rsefiles <- list.files(path = savefilepath, pattern = "_rse.rds", full.names = TRUE)
  print(rsefiles)
  
  foreach(i = 1:length(rsefiles), .combine=rbind) %dopar% { #run loops in parallel, outputs a combined list
    rse <- readRDS(rsefiles[i])
    meta <-  as.data.frame(rse@colData)
  }
}

####FUNCTION altered from jasondubois/sportfish 
ReadRDSFiles <- function(fileDir, envir = .GlobalEnv) {
  
  # pattern for file searching
  p <- ".rds$"
  
  rds <- list.files(path = fileDir, pattern = p)
  
  out <- vapply(rds, FUN = function(.x) {
    
    nm <- sub(pattern = p, replacement = "", x = .x)
    
    # LW alteration to then be able to pull colData
    out <- readRDS(file = file.path(fileDir, .x))
    
    # load in global env
    assign(nm, value = as.data.frame(out@colData), envir = envir)
    
    if (!exists(nm, envir = envir)) return(FALSE)
    
    TRUE
    
  }, FUN.VALUE = logical(1L), USE.NAMES = FALSE)
  
  if (!all(out)) warning("Some `.rds` files not loaded.", call. = FALSE)
  
  spc <- paste0(rep('*', times = nchar(fileDir) + 1), collapse = "")
  
  cat("RDS Files loaded from:", fileDir, spc, rds[out], sep = "\n ", spc)
}
