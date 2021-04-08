#########################################
## To prepare step 1 file that will be used in regenie/BOLT/fastGWA/SAIGE in UKB WB samples
##

## array bed
step1file_pre=$1
## imputed bgen
step2file=$2
## list of FID/IID for UKB WB samples with covariate info
sample_keep=$3
# bed region file listing ICLD and low-complexity regions
regions_ignore=$4
outprefix=$5

## subset to samples with array & imputed data
grep -wFf $sample_keep ${step1file_pre}.fam | cut -f1,2 > samples_keep.tmp
grep -wFf samples_keep.tmp <(cat ${step2file}.sample | tr ' ' '\t') | cut -f1,2 > samples_keep.tmp1

## filters for step 1
## samples: %missing < 10%
## variants: MAF>1%, %miss<1%, HWE p>1e-15, no ICLD, no low-complexity regions, light LD pruning
plink2 \
  --write-samples \
  --write-snplist \
  --autosome \
  --bfile $step1file_pre \
  --keep samples_keep.tmp1 \
  --mind 0.1 \
  --geno 0.01 --maf 0.01 --hwe 1e-15 \
  --exclude range $regions_ignore \
  --indep-pairwise 1000 100 0.9 \
  --out $outprefix

rm samples_keep.tmp samples_keep.tmp1

## make bed file applying filters
plink2 \
  --make-bed \
  --bfile $step1file_pre \
  --keep ${outprefix}.id \
  --extract ${outprefix}.prune.in \
  --out $outprefix


