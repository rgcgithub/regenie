#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if(!require(data.table)){ install.packages("data.table"); library(data.table) }
  if(!require(tidyverse)){ install.packages("tidyverse"); library(tidyverse) }
  if(!require(optparse)){ install.packages("optparse"); library(optparse) }
})
#########################################
##
## Script used to make sparse GRM file for fastGWA
## for the analyses in the REGENIE 2020 paper
## For more details, visit: https://rgcgithub.github.io/regenie/
##  
#########################################

option_list = list(
  make_option("--step1File", type="character", default="",
    help="bed file prefix for step 1", metavar="string"),
  make_option("--prefix", type="character", default="",
    help="output files prefix", metavar="string"),
  make_option("--partition", type="integer", default=0,
    help="partition number for fastGWA GRM computation"),
  make_option("--npartitions", type="integer", default=250,
    help="number of partitions for fastGWA GRM computation")
  );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

total.partitions <- opt$npartitions
if( opt$partition > total.partitions) stop("Invalid argument")
print(opt)

########### Functions ############
fastGWA.computeGRM <- function(){

  fastGWA.call <- paste0("gcta64 ",
    "--bfile ", bed.file, " ",
    "--thread-num ", parallel::detectCores(), " ",
    "--out ", outprefix, " ",
    "--make-grm-part ", total.partitions, " ", opt$partition)
  print(fastGWA.call)

  tot.time <- system.time({
    system(fastGWA.call)
  })

  write( tot.time, paste0(outprefix, "_time_pt_", opt$partition), ncol = length(tot.time))
}

fastGWA.compute.sparseGRM <- function(){

  ## Combine the partitions
  system( paste0("cat ", outprefix,
      ".part_",total.partitions,"_*.grm.id > ",outprefix, ".grm.id") )
  system( paste0("cat ", outprefix,
      ".part_",total.partitions,"_*.grm.bin > ",outprefix, ".grm.bin") )
  system( paste0("cat ", outprefix,
      ".part_",total.partitions,"_*.grm.N.bin > ",outprefix, ".grm.N.bin") )

  # compute sparse GRM
  fastGWA.call <- paste0("gcta64 ",
    "--grm ", outprefix, " ",
    "--make-bK-sparse 0.05 ",
    "--out ", outprefix, "_sp ")
  system(fastGWA.call)

  if(file.exists( paste0(outprefix, "_sp.grm.id") )){ # in case of error so no need to redo
    system( paste0("rm ", outprefix,".part_",total.partitions,"_*.grm.*") ) 
  }
}

#################################
outprefix <- opt$prefix
bed.file <- opt$step1File

if( opt$partition > 0 ){ # run fastGWA GRM computation partitions

  fastGWA.computeGRM()

} else if(opt$partition == 0){ # finish fastGWA sparse GRM computation

  tot.time <- system.time({
    fastGWA.compute.sparseGRM()
  })
  write( tot.time, paste0(outprefix, "_time_compSparseGRM"), ncol = length(tot.time))

} 

