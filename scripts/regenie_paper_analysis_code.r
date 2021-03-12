
suppressPackageStartupMessages({
  if(!require(data.table)){ install.packages("data.table"); library(data.table) }
  if(!require(optparse)){ install.packages("optparse"); library(optparse) }
})
#########################################
##
## Script used to run REGENIE for GWAS
## for the analyses in the REGENIE 2020 paper
## For more details, visit: https://rgcgithub.github.io/regenie/
##  
#########################################
option_list = list(
  make_option("--step1File", type="character", default="",
    help="bed file prefix for step 1", metavar="string"),
  make_option("--step2File", type="character", default="",
    help="bgen file for step 2", metavar="string"),
  make_option("--phenoFile", type="character", default="",
    help="phenotype file", metavar="string"),
  make_option("--covarFile", type="character", default="",
    help="covariate file", metavar="string"),
  make_option("--prefix", type="character", default="",
    help="output files prefix", metavar="string"),
  make_option("--pheno", type="integer", default=0,
    help="which phenotype column to run [default is all for regenie]", metavar="number"),
  make_option("--bt", action="store_true", default=FALSE,
    help="run in BT mode"),
  make_option("--lowmem", action="store_true", default=FALSE,
    help="run regenie with lowmem option"),
  make_option("--loocv", action="store_true", default=FALSE,
    help="run regenie with LOOCV"),
  make_option("--chr", type="integer", default=0,
    help="chromosome to test", metavar="number"),
  make_option("--noapprox", action="store_true", default=FALSE,
    help="use exact Firth"),
  make_option("--spa", action="store_true", default=FALSE,
    help="use SPA"),
  make_option("--skipNull", action="store_true", default=FALSE,
    help="Run step 2"),
  make_option("--binary", type="character", default="regenie",
    help="path to regenie binary", metavar="string")
  );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if( opt$skipNull & (opt$chr == 0) ) stop("Invalid chromosome number.")
if(!opt$skipNull & !file.exists(paste0(opt$step1File,".bed"))) stop("Step 1 file does not exists")
if(opt$skipNull & !file.exists(opt$step2File)) stop("Step 2 file does not exists")
if(!file.exists(opt$phenoFile)) stop("Phenotype file does not exists")
if(!file.exists(opt$covarFile)) stop("Covariate file does not exists")

cat("##### Running regenie ######\n")
print(opt)

########### Functions ############
# fit regenie
fit.regenie <- function() {

  mode.rg <- test.type <- phenoCols <- rg.suffix <- ""

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
      rg.suffix <- paste0(rg.suffix, "_SPA" )
    } else if(opt$noapprox){
      mode.rg <- paste0(mode.rg, " --firth" )
      rg.suffix <- paste0(rg.suffix, "_FirthExact" )
    } else{
      mode.rg <- paste0(mode.rg, " --firth --approx" )
      rg.suffix <- paste0(rg.suffix, "_Firth" )
    }
  }

  if( !opt$skipNull) { # step 1

    # regenie call
    rg.call <- paste0(opt$binary, " ",
      "--bed ", bed.file, " ",
      "--phenoFile ", pheno.file, " ", phenoCols, " ",
      "--covarFile ", covar.file, " ",
      "--bsize 1000 ",
      "--step 1 ",
      mode.rg, " ", 
      "--out ", outprefix, rg.suffix
    )
    rg.time.suffix <- "_time.step1"

  } else { # step 2

    # regenie call
    rg.call <- paste0(opt$binary, " ",
      "--bgen ", bgen.file, " ",
      "--phenoFile ", pheno.file, " ", phenoCols, " ",
      "--covarFile ", covar.file, " ",
      "--bsize 400 ",
      "--step 2 ",
      mode.rg, " ", 
      "--pred ", outprefix, rg.suffix, "_pred.list ",
      "--out ", outprefix, rg.suffix, "_chr", opt$chr
    )
    rg.time.suffix <- paste0("_time.chr", opt$chr)

  }

  cat( cmd <- paste0("/usr/bin/time -v ", rg.call, " > ", track.file, " 2>&1" ) )
  t0 <- system.time(system(cmd))
  write( t0, paste0(outprefix, rg.suffix, rg.time.suffix), ncol = length(t0))

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
# regenie all phenos
if( opt$pheno == 0 & !opt$skipNull) 
  track.file <- paste0(track.file, "_step1")
if( opt$pheno == 0 & opt$skipNull)
  track.file <- paste0(track.file, "_step2_chr", opt$chr)
# regenie by phenotype
if( opt$pheno != 0 & !opt$skipNull) 
  track.file <- paste0(track.file, "_step1_phenoCol", opt$pheno)
if( opt$pheno != 0 & opt$skipNull) 
  track.file <- paste0(track.file, "_step2_phenoCol", opt$pheno, "_chr", opt$chr)

if( !opt$lowmem) track.file <- paste0(track.file, "_nowrite")
if( opt$loocv) track.file <- paste0(track.file, "_loocv")
if( opt$skipNull & opt$bt) 
  track.file <- paste0(track.file, ifelse(opt$spa, "_spa", ifelse(opt$noapprox, "_firthexact",  "_firth")))

track.file <- paste0(track.file, ".log")
cat( paste0("Tracking memory in file: ", track.file),"\n")

fit.regenie()

