#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if(!require(data.table)){ install.packages("data.table"); library(data.table) }
  if(!require(tidyverse)){ install.packages("tidyverse"); library(tidyverse) }
  if(!require(optparse)){ install.packages("optparse"); library(optparse) }
})
#########################################
##
## Script used to run REGENIE/BOLT/fastGWA/SAIGE for GWAS
## for the analyses in the REGENIE 2020 paper
## For more details, visit: https://rgcgithub.github.io/regenie/
##  
#########################################
option_list = list(
  make_option("--method", type="character", default="",
    help="method to run in analysis"),
  make_option("--step1File", type="character", default="",
    help="bed file prefix for step 1", metavar="string"),
  make_option("--step2File", type="character", default="",
    help="bgen file (or prefix) for step 2", metavar="string"),
  make_option("--phenoFile", type="character", default="",
    help="phenotype file", metavar="string"),
  make_option("--covarFile", type="character", default="",
    help="covariate file", metavar="string"),
  make_option("--prefix", type="character", default="",
    help="output files prefix", metavar="string"),
  make_option("--pheno", type="integer", default=0,
    help="which phenotype column to run [default is all for regenie]", metavar="number"),
  make_option("--bt", action="store_true", default=FALSE,
    help="run regenie in BT mode"),
  make_option("--lowmem", action="store_true", default=FALSE,
    help="run regenie with lowmem option"),
  make_option("--loocv", action="store_true", default=FALSE,
    help="run regenie with LOOCV"),
  make_option("--grm", type="character", default="",
    help="path to sparse GRM for fastGWA", metavar="string"),
  make_option("--chr", type="integer", default=0,
    help="chromosome to test", metavar="number"),
  make_option("--noapprox", action="store_true", default=FALSE,
    help="use exact Firth"),
  make_option("--spa", action="store_true", default=FALSE,
    help="use SPA"),
  make_option("--skipNull", action="store_true", default=FALSE,
    help="Run step 2")
  );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(!file.exists(opt$phenoFile)) stop("Phenotype file does not exist")
if(!file.exists(opt$covarFile)) stop("Covariate file does not exist")

print(opt)

########### Functions ############
# fit regenie
fit.regenie <- function() {

  mode.rg <- test.type <- phenoCols <- rg.suffix <- rg.suffix2 <- ""

  # get phenotype name
  phenoNames <- fread(pheno.file) %>% select(-FID,-IID) %>% colnames
  if( opt$pheno > 0) {
    if(opt$pheno > length(phenoNames)) stop("Invalid phenotype column provided")
    phenoCols <- paste0("--phenoCol ", phenoNames[opt$pheno])
    rg.suffix <- paste0("_phenoCol", opt$pheno) 
  }

  # BT mode
  if( opt$bt ) mode.rg <- "--bt"
  # cv mode
  if( opt$loocv ){
    rg.suffix <- paste0(rg.suffix, "_loocv" )
    mode.rg <- paste0(mode.rg, " --loocv" )
  }
  # lowmem mode
  if( opt$lowmem ) 
    mode.rg <- paste0(mode.rg, " --lowmem" )
  else
    rg.suffix <- paste0(rg.suffix, "_nowrite" )
  # spa/firth
  if( opt$skipNull & opt$bt ) {
    if(opt$spa){
      mode.rg <- paste0(mode.rg, " --spa" )
      rg.suffix2 <- "_SPA"
    } else if(opt$noapprox){
      mode.rg <- paste0(mode.rg, " --firth" )
      rg.suffix2 <- "_FirthExact"
    } else{
      mode.rg <- paste0(mode.rg, " --firth --approx" )
      rg.suffix2 <- "_Firth"
    }
  }

  if( !opt$skipNull) { # step 1

    # regenie call
    rg.call <- paste0("regenie ",
      "--bed ", bed.file, " ",
      "--phenoFile ", pheno.file, " ", phenoCols, " ",
      "--covarFile ", covar.file, " ",
      "--bsize 1000 ",
      "--step 1 ",
      mode.rg, " ", 
      "--threads ", parallel::detectCores()," ",
      "--out ", outprefix, rg.suffix
    )
    rg.time.suffix <- "_time.step1"

  } else { # step 2 for each chromosome

    # regenie call
    rg.call <- paste0("regenie ",
      "--bgen ", bgen.file, opt$chr, ".bgen ",
      "--phenoFile ", pheno.file, " ", phenoCols, " ",
      "--covarFile ", covar.file, " ",
      "--bsize 400 ",
      "--step 2 ",
      mode.rg, " ", 
      "--threads ", parallel::detectCores()," ",
      "--pred ", outprefix, rg.suffix, "_pred.list ",
      "--out ", outprefix, rg.suffix, rg.suffix2, "_chr", opt$chr
    )
    rg.time.suffix <- paste0("_time.chr", opt$chr)

  }

  cat( cmd <- paste0("/usr/bin/time -v ", rg.call, " > ", track.file, " 2>&1" ) )
  t0 <- system.time(system(cmd))
  write( t0, paste0(outprefix, rg.suffix, rg.suffix2, rg.time.suffix), ncol = length(t0))

  return(NULL)
}

# fit bolt (step 2 for all chromosomes)
fit.bolt <- function() {

  phenoNames <- fread(pheno.file) %>% select(-FID,-IID) %>% colnames
  phenoCol <- phenoNames[opt$pheno]
  covNames <- fread(covar.file) %>% select(-FID,-IID) %>% colnames
  split.bgen.files <- paste0(bgen.file, 1:22, ".bgen")

  bolt.call <- paste0("bolt ",
    "--bfile=", bed.file, "  ",
    "--phenoFile=", pheno.file, " ",
    "--phenoCol=", phenoCol, " ",
    "--covarFile=", covar.file, " ",
    paste0("--qCovarCol=",covNames, collapse=" "), " ",
    "--lmmForceNonInf ",
    "--LDscoresUseChip ",
    "--numThreads=", parallel::detectCores()," ",
    "--statsFile=", outprefix,"_phenoCol", opt$pheno ,".grm ",
    "--predBetasFile=", outprefix,"_phenoCol", opt$pheno ,".loco ",
    paste( paste0("--bgenFile=", split.bgen.files), collapse=" "), " ",
    "--sampleFile=", gsub("bgen", "sample", split.bgen.files[1]), " ",
    "--statsFileBgenSnps=", outprefix, "_phenoCol", opt$pheno ,".test"
  )
  bolt.time.suffix <- paste0("_time")

  cat( cmd <- paste0("/usr/bin/time -v ", bolt.call, " > ", track.file, " 2>&1" ) )
  t0 <- system.time(system(cmd))
  write( t0, paste0(outprefix, "_phenoCol", opt$pheno , bolt.time.suffix), ncol = length(t0))

  return(NULL)
}

# fit saige
fit.saige <- function() {

  if( !opt$skipNull) { ## step 1

    cov.file <- fread(covar.file) %>% select(-FID,-IID)
    ncov <- ncol(cov.file)
    colnames(cov.file) <- paste0("V",1:ncov)
    dp <- fread(pheno.file)
    colnames(dp) <- c("FID","IID", paste0("P",1:(ncol(dp)-2)))

    fwrite(cbind(dp, cov.file), paste0(outprefix, "_phenoCol", opt$pheno, ".yout.SAIGE"), 
      sep=" ", quote=FALSE, na="NA")
    phenoNames <- dp %>% select(-FID,-IID) %>% colnames
    phenoCol <- phenoNames[opt$pheno]

    saige.call <- paste0("step1_fitNULLGLMM.R ",
      "--plinkFile=", bed.file, " ",
      "--phenoFile=", outprefix, "_phenoCol", opt$pheno, ".yout.SAIGE ",
      "--phenoCol=P", opt$pheno, " ",
      "--covarColList=", paste0(paste0("V",1:ncov), collapse=",")," ",
      "--sampleIDColinphenoFile=IID ",
      "--traitType=binary ",
      "--LOCO=TRUE ",
      "--minMAFforGRM=0.0001 ",
      "--nThreads=", parallel::detectCores()," ",
      "--outputPrefix=", outprefix, "_phenoCol", opt$pheno
    )
    saige.time.suffix <- "_time.step1"

  } else { # step 2 for each chromosome

    fread( paste0(bgen.file, opt$chr, ".sample") )  %>%
      select(ID_1) %>%
      slice(-1) %>%
      fwrite(file=paste0(outprefix, "_phenoCol", opt$pheno, ".sample.SAIGE"), quote=FALSE, na="NA", col.names=FALSE)

    saige.call <- paste0("step2_SPAtests.R ",
      "--sampleFile=", outprefix, "_phenoCol", opt$pheno, ".sample.SAIGE ",
      "--LOCO=TRUE ",
      "--chrom=", opt$chr, " ",
      "--bgenFile=", bgen.file, opt$chr,".bgen ",
      "--bgenFileIndex=", bgen.file, opt$chr,".bgen.bgi ",
      "--GMMATmodelFile=", outprefix, "_phenoCol", opt$pheno, ".rda ",
      "--varianceRatioFile=", outprefix, "_phenoCol", opt$pheno, ".varianceRatio.txt ",
      "--SAIGEOutputFile=", outprefix, "_phenoCol", opt$pheno, "_chr", opt$chr, ".test"
    )
    saige.time.suffix <- paste0("_time.chr", opt$chr)
  }

  cat( cmd <- paste0("/usr/bin/time -v ", saige.call, " > ", track.file, " 2>&1" ) )
  t0 <- system.time( system(cmd) )
  write( t0, paste0(outprefix, "_phenoCol", opt$pheno, saige.time.suffix), ncol = length(t0))

  return(NULL)
}

#fastGWA
fit.fastGWA <- function(){

  fread(pheno.file) %>%
    fwrite( paste0(outprefix, ".pheno.fastGWA"), col.names=FALSE, na=NA, sep=" ",quote=FALSE)
  fread(covar.file) %>%
    fwrite( paste0(outprefix, ".covar.fastGWA"), col.names=FALSE, na=NA, sep=" ",quote=FALSE)
  fread(pheno.file) %>% select(FID, IID) %>%
    fwrite(paste0(outprefix, ".sample.fastGWA"), col.names=FALSE, na=NA, sep=" ",quote=FALSE)

  phenoNames <- fread(pheno.file) %>% select(-FID,-IID) %>% colnames
  phenoCol <- phenoNames[opt$pheno]
  data.frame( bgen.names = paste0(bgen.file,1:22, ".bgen") ) %>%
    fwrite(paste0(outprefix, "_phenoCol", opt$pheno, ".bgen.list"),
      quote=F, col.names=F)

  fastGWA.call <- paste0("gcta64 ",
    "--pheno ", outprefix, ".pheno.fastGWA ",
    "--mpheno ", opt$pheno, " ", 
    "--qcovar ", outprefix, ".covar.fastGWA ",
    "--keep ", outprefix, ".sample.fastGWA ",
    "--fastGWA-mlm --h2-limit 2.5 --maf 0 --geno 1 ",
    "--grm-sparse ", opt$grm, " ",
    "--mbgen ", outprefix, "_phenoCol", opt$pheno, ".bgen.list ",
    "--sample ", bgen.file, "1.sample ",
    "--threads ", parallel::detectCores()," ",
    "--out ", outprefix, "_phenoCol", opt$pheno, ".test"
  )

  cat( cmd <- paste0("/usr/bin/time -v ", fastGWA.call, " > ", track.file, " 2>&1" ) )
  t0 <- system.time( system(cmd) )
  write( t0, paste0(outprefix, "_phenoCol", opt$pheno, "_time"), ncol = length(t0))

  return(NULL)
}

#################################
# main paths + files
outprefix <- opt$prefix

pheno.file <- opt$phenoFile
covar.file <- opt$covarFile
bed.file <- opt$step1File
bgen.file <- opt$step2File

## tracking memory usage
track.file <- paste0(outprefix,"_timing")
# all phenos
if( opt$pheno == 0 ){
  if(!opt$skipNull) 
    track.file <- paste0(track.file, "_step1")
  else {
    track.file <- paste0(track.file, "_step2")
    if( opt$chr != 0 ) track.file <- paste0(track.file, "_chr", opt$chr)
  }
} else { # by phenotype
  if(!opt$skipNull) 
    track.file <- paste0(track.file, "_step1_phenoCol", opt$pheno)
  else { 
    track.file <- paste0(track.file, "_step2_phenoCol", opt$pheno)
    if( opt$chr != 0 ) track.file <- paste0(track.file, "_chr", opt$chr)
  }
}

if( opt$method == "regenie") {
  if( !opt$lowmem) track.file <- paste0(track.file, "_nowrite")
  if( opt$loocv) track.file <- paste0(track.file, "_loocv")
  if( opt$skipNull & opt$bt) 
    track.file <- paste0(track.file, ifelse(opt$spa, "_spa", ifelse(opt$noapprox, "_firthexact",  "_firth")))
}

track.file <- paste0(track.file, ".log")
cat( paste0("Tracking memory in file: ", track.file),"\n")

if( opt$method == "regenie") {
  fit.regenie()
} else if( opt$method == "bolt") {
  fit.bolt()
} else if( opt$method == "saige") {
  fit.saige()
} else if( opt$method == "fastgwa") {
  fit.fastGWA()
} else {
  stop("Invalid method")
}

