#!/usr/bin/Rscript

### This function returns the LD matrix computed from the compressed binary file 
### output from REGENIE
###   ld.file: binary file output from Regenie

get.corr.sq.matrix <- function( ld.file = NULL ){

  if(is.null(ld.file) || !file.exists(ld.file))
    stop("Need to pass valid LD file!")
  list.file <- paste0(ld.file, ".snplist")
  if(!file.exists(list.file))
    stop("Cannot open accompagnying snplist file: ", list.file)

  # Variant IDs
  snplist <- read.table(list.file)$V1
  n.snps <- length(snplist)

  rfile <- file(ld.file, "rb")
  # number of samples and snps
  f.info <- readBin(
    con = rfile, 
    what = "integer", 
    n = 2,
  ) 
  if(f.info[2] != n.snps)
    stop("Number of variants from snplist does not match that in correlation file")

  # R^2 as stored using uint16 representation
  rvals <- readBin(
    con = rfile, 
    what = "integer", 
    n = choose(n.snps, 2),
    signed = FALSE,
    size = 2
  ) 
  close(rfile)

  # matrix of R^2 values
  rmat <- matrix(0, n.snps, n.snps)
  rmat[ lower.tri(rmat) ] <- rvals / (2^16 - 1)
  rmat[ upper.tri(rmat) ] <- t(rmat)[ upper.tri(rmat) ]
  diag(rmat) <- 1
  colnames(rmat) <- rownames(rmat) <- snplist

  return(rmat)
}

