#!/usr/bin/env bash

### REGENIE TEST SCRIPT 
## For version<1.0.5.6, will get error if WITH_GZ is set since option '--gz' did not exist

info_msg="Usage: ./test_bash.sh OPTIONS\n"
info_msg+="  --path  path to Regenie repository\n"
info_msg+="  --gz    Flag to specify compilation was done with Boost Iostream library\n"
if [ "$#" -eq 0 ]; then
  echo -e "$info_msg"; exit 1
fi

# Force only Y1 to be analyzed 
# (single trait runs should produce identical results 
#   removing NA rows of trait from phenotype file)
phenoColArg="--phenoCol Y1"
btArg="--bt"

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --path) REGENIE_PATH="$2"; shift ;;
    --gz) WITH_GZ=1 ;;
    --all) phenoColArg= ;;
    --qt) btArg= ;;
    *) echo -e "Unknown parameter passed: $1.\n$info_msg"; exit 1 ;;
  esac
  shift
done


# quick check src/example folders are present
if [ ! -d "${REGENIE_PATH}/src" ] || [ ! -d "${REGENIE_PATH}/example" ]; then
  echo "ERROR: First input argument must be the directory where Regenie repo was cloned"; exit 1
else
  cd $REGENIE_PATH
fi 

# If compiling was done with Boost Iostream library, use gzipped files as input
if [ "$WITH_GZ" = "1" ]; then
  fsuf=.gz
  arg_gz="--gz"
fi

REGENIE_PATH=$(pwd)/  # use absolute path
regenie_bin=`ls regenie* | head -n 1`
help_msg="Update to most recent REGENIE version (using 'git pull') and re-compile the software (using 'make clean && make')."

if [ ! -f "$regenie_bin" ]; then
  echo "ERROR: Regenie binary cannot be found. Compile the software first using 'make clean && make'"; exit 1
fi



# Run regenie on trait Y1 (missing will be dropped automatically)
echo -e "Regenie on Y1\n=================================="
rgcmd="--step 1 \
  --bed example/example \
  --covarFile example/covariates.txt${fsuf} \
  --phenoFile example/phenotype_bin_wNA.txt \
  $phenoColArg \
  --bsize 100 \
  $btArg \
  $arg_gz \
  --lowmem \
  --lowmem-prefix tmp_rg \
  --out test/fit_bin_out"

# run regenie
./$regenie_bin $rgcmd

# step 2
rgcmd="--step 2 \
  --bed example/example_3chr \
  --covarFile example/covariates.txt \
  --phenoFile example/phenotype_bin_wNA.txt \
  $phenoColArg \
  --bsize 200 \
  $btArg \
  --firth --approx \
  --pThresh 0.01 \
  --pred test/fit_bin_out_pred.list \
  --out test/test_bin_out_firth"

# run regenie
./$regenie_bin $rgcmd


###
# Drop missing from the phenotype file for Y1
###
echo -e "Regenie on Y1 dropping missing from phenotype file\n=================================================="
echo "Making new phenotype/covariate files dropping NAs (test/test_bin_noNA*.txt)..."
grep -v "NA" example/phenotype_bin_wNA.txt > test/test_bin_noNA.txt
grep -wFf <(cut -f1,2 -d ' ' test/test_bin_noNA.txt) example/covariates.txt > test/test_bin_noNA_covs.txt
echo "Comparing original files to that dropping individuals NAs"
if cmp --silent test/test_bin_noNA.txt example/phenotype_bin_wNA.txt && \
 cmp --silent test/test_bin_noNA_covs.txt example/covariates.txt ; then
  echo "Uh oh, files are the same!"; exit 1
else
  echo -e "Files are different; proceeding forward...\n"
fi

rgcmd="--step 1 \
  --bed example/example \
  --covarFile test/test_bin_noNA_covs.txt \
  --phenoFile test/test_bin_noNA.txt \
  $phenoColArg \
  --bsize 100 \
  $btArg \
  $arg_gz \
  --lowmem \
  --lowmem-prefix tmp_rg \
  --out test/fit_bin_out_noNA"

# run regenie
./$regenie_bin $rgcmd

# step 2
rgcmd="--step 2 \
  --bed example/example_3chr \
  --covarFile test/test_bin_noNA_covs.txt \
  --phenoFile test/test_bin_noNA.txt \
  $phenoColArg \
  --bsize 200 \
  $btArg \
  --firth --approx \
  --pThresh 0.01 \
  --pred test/fit_bin_out_noNA_pred.list \
  --out test/test_bin_out_firth_noNA"

# run regenie
./$regenie_bin $rgcmd

# check files are not empty
if [ "`cat test/test_bin_out_firth.regenie | wc -l`" -le 10 ]; then
  echo "Uh oh, result file is empty!"; exit 1
elif cmp --silent test/test_bin_out_firth.regenie test/test_bin_out_firth_noNA.regenie ; then
  echo "Result files are the same"
else
  echo "Uh oh, result files are different"
fi

# cleanup
rm test/test_bin_* test/fit_bin_out*
