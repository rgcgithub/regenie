#!/usr/bin/env bash

### REGENIE TEST SCRIPT for version >= 1.0.6.2
## For version<1.0.6.2, will get error since option '--ref-first' did not exist
## For version<1.0.6.1, will get error since option '--print-pheno' did not exist
## For version<1.0.5.6, will get error if WITH_GZ is set since option '--gz' did not exist

info_msg="Usage: ./test_bash.sh OPTIONS\n"
info_msg+="  --path  path to Regenie repository\n"
info_msg+="  --gz    Flag to specify compilation was done with Boost Iostream library\n"
if [ "$#" -eq 0 ]; then
  echo -e "$info_msg"; exit 1
fi

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --path) REGENIE_PATH="$2"; shift ;;
    --gz) WITH_GZ=1 ;;
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
mntpt=
regenie_bin=`ls regenie* | head -n 1`
help_msg="Update to most recent REGENIE version (using 'git pull') and re-compile the software (using 'make clean && make')."

if [ ! -f "$regenie_bin" ]; then
  echo "ERROR: Regenie binary cannot be found. Compile the software first using 'make clean && make'"; exit 1
fi


# Prepare regenie command to run for Step 1
echo -e "Running step 1 of REGENIE\n=================================="
basecmd="--step 1 \
  --bed ${mntpt}example/example \
  --exclude ${mntpt}example/snplist_rm.txt \
  --covarFile ${mntpt}example/covariates.txt${fsuf} \
  --phenoFile ${mntpt}example/phenotype_bin.txt${fsuf} \
  --remove ${mntpt}example/fid_iid_to_remove.txt \
  --bsize 100 \
  --bt $arg_gz"
rgcmd="$basecmd \
  --lowmem \
  --lowmem-prefix tmp_rg \
  --out ${mntpt}test/fit_bin_out"

# run regenie
./$regenie_bin $rgcmd

## quick check that the correct files have been created
if [ ! -f "${REGENIE_PATH}test/fit_bin_out.log" ] || \
  [ ! -f "${REGENIE_PATH}test/fit_bin_out_pred.list" ] || \
  [ ! -f "${REGENIE_PATH}test/fit_bin_out_1.loco$fsuf" ] || \
  [ ! -f "${REGENIE_PATH}test/fit_bin_out_2.loco$fsuf" ]; then
  echo "Step 1 of REGENIE did not finish successfully. $help_msg"; exit 1
elif [ "`grep \"0.4504\" ${REGENIE_PATH}test/fit_bin_out.log | grep \"min value\"`" = "" ]; then
  echo "Step 1 of REGENIE did not finish successfully. $help_msg"; exit 1
fi

#### Run step 1 splitting across jobs for level 0
njobs=4
echo -e "Re-running step 1 splitting in $njobs jobs"
# pt1 - run regenie before l0
rgcmd="$basecmd \
  --split-l0 ${mntpt}test/fit_bin_parallel,$njobs \
  --out ${mntpt}test/fit_bin_l0"

./$regenie_bin $rgcmd
if [ ! -f "${REGENIE_PATH}test/fit_bin_parallel.master" ]; then
  echo "Step 1 of REGENIE did not finish successfully. $help_msg"; exit 1
fi

# pt2 - run regenie for l0
nj=`seq 1 $njobs`
for job in $nj; do
  rgcmd="$basecmd \
    --run-l0 ${mntpt}test/fit_bin_parallel.master,$job \
    --out ${mntpt}test/fit_bin_l0"

  ./$regenie_bin $rgcmd
  if [ ! -f "${REGENIE_PATH}test/fit_bin_parallel_job${job}_l0_Y1" ]; then
    echo "Step 1 of REGENIE did not finish successfully. $help_msg"; exit 1
  fi
done


# pt3 - run regenie for l1
rgcmd="$basecmd \
  --run-l1 ${mntpt}test/fit_bin_parallel.master \
  --out ${mntpt}test/fit_bin_l1"

./$regenie_bin $rgcmd

if [ ! -f "${REGENIE_PATH}test/fit_bin_l1_1.loco$fsuf" ]; then
  echo "Step 1 of REGENIE did not finish successfully. $help_msg"; exit 1
elif ! cmp --silent \
  "${REGENIE_PATH}test/fit_bin_out_1.loco$fsuf" \
  "${REGENIE_PATH}test/fit_bin_l1_1.loco$fsuf" 
then
  echo "Uh oh, REGENIE did not build successfully. $help_msg"; exit 1
elif ! cmp --silent \
  "${REGENIE_PATH}test/fit_bin_out_2.loco$fsuf" \
  "${REGENIE_PATH}test/fit_bin_l1_2.loco$fsuf" 
then
  echo "Uh oh, REGENIE did not build successfully. $help_msg"; exit 1
fi


# First step 2 command
echo -e "Running step 2 of REGENIE\n=================================="
rgcmd="--step 2 \
  --bgen ${mntpt}example/example.bgen \
  --covarFile ${mntpt}example/covariates.txt${fsuf} \
  --phenoFile ${mntpt}example/phenotype_bin.txt${fsuf} \
  --remove ${mntpt}example/fid_iid_to_remove.txt \
  --bsize 200 \
  --bt \
  --firth --approx \
  --pThresh 0.01 \
  --pred ${mntpt}test/fit_bin_out_pred.list \
  --split $arg_gz \
  --out ${mntpt}test/test_bin_out_firth"

# run regenie
./$regenie_bin $rgcmd

##  do this way so zcat works on OSX
if [ -f ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie.gz ]; then
  ( zcat < ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie.gz ) > ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie
fi

echo "------------------------------------------"
if cmp --silent \
  ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie \
  ${REGENIE_PATH}example/example.test_bin_out_firth_Y1.regenie 
then
  echo -e "Files are identical.\n\nRunning second test...\n"
else
  echo -e "ERROR: Uh oh... Files are different! $help_msg"; exit 1
fi


# Second command
rgcmd="--step 2 \
  --bed ${mntpt}example/example_3chr \
  --ref-first \
  --covarFile ${mntpt}example/covariates.txt${fsuf} \
  --phenoFile ${mntpt}example/phenotype_bin.txt${fsuf} \
  --phenoColList Y2 \
  --bsize 100 \
  --chrList 2,3 \
  --test dominant \
  --ignore-pred \
  --write-samples \
  --print-pheno \
  --out ${mntpt}test/test_out"

# run regenie
./$regenie_bin $rgcmd

# check files
echo "------------------------------------------"
if [ ! -f "${REGENIE_PATH}test/test_out_Y2.regenie.ids" -o -f "${REGENIE_PATH}test/test_out_Y1.regenie.ids" ]
then
  echo "Uh oh, REGENIE did not build successfully. $help_msg"; exit 1
elif (( $(head -n 1 ${REGENIE_PATH}test/test_out_Y2.regenie.ids | cut -f1) != "Y2" )); then
  echo "Uh oh, REGENIE did not build successfully. $help_msg"; exit 1
elif (( $(head -n 1 "${REGENIE_PATH}test/test_out_Y2.regenie.ids" | tr '\t' '\n' | wc -l) != 2 )); then
  echo "Uh oh, REGENIE did not build successfully. $help_msg"; exit 1
elif (( `grep "mog_" "${REGENIE_PATH}test/test_out.regenie" | wc -l` > 0 )); then
  echo "Uh oh, REGENIE did not build successfully. $help_msg"; exit 1
elif (( `grep "ADD" "${REGENIE_PATH}test/test_out.regenie" | wc -l` > 0 )); then
  echo "Uh oh, REGENIE did not build successfully. $help_msg"; exit 1
elif [ "`cut -d ' ' -f1-5 ${REGENIE_PATH}test/test_out.regenie | sed '2q;d'`" != "`grep \"^2\" ${REGENIE_PATH}example/example_3chr.bim | head -n 1 | awk '{print $1,$4,$2,$5,$6}'`" ]; then
  echo "Uh oh, REGENIE did not build successfully. $help_msg"; exit 1
fi


echo -e "Passed.\n\nRunning third test...\n"
# Third command
rgcmd="--step 2 \
  --bed ${mntpt}example/example_3chr \
  --ref-first \
  --extract ${mntpt}test/test_out.snplist \
  --covarFile ${mntpt}example/covariates.txt${fsuf} \
  --phenoFile ${mntpt}example/phenotype_bin.txt${fsuf} \
  --phenoColList Y2 \
  --bsize 100 \
  --test dominant \
  --ignore-pred \
  --out ${mntpt}test/test_out_extract"

grep -v "^1" ${REGENIE_PATH}example/example_3chr.bim | awk '{ print $2 }' > ${REGENIE_PATH}test/test_out.snplist

# run regenie
./$regenie_bin $rgcmd

if cmp --silent \
  ${REGENIE_PATH}test/test_out.regenie \
  ${REGENIE_PATH}test/test_out_extract.regenie 
then
  echo "SUCCESS: REGENIE build passed the tests!"
else
  echo "Uh oh, REGENIE did not build successfully. $help_msg"; exit 1
fi


# file cleanup
rm ${REGENIE_PATH}test/fit_bin* ${REGENIE_PATH}test/test_bin_out_firth* ${REGENIE_PATH}test/test_out*

