#!/usr/bin/env bash

### REGENIE TEST SCRIPT TO USE WITH DOCKER IMAGE ###
#### Working from directory where regenie repo was cloned
if [ "$#" -eq 0 ]; then
  echo "Usage: test.sh <PATH_TO_CLONED_REGENIE_REPO> <DOCKER_IMAGE_TAG>"; exit 1
fi

### Test script for Regenie version >= 1.0.6.2
## For version<1.0.6.1, will get error since option '--ref-first' did not exist
## For version<1.0.6.1, will get error since option '--print-pheno' did not exist
## For version<1.0.5.6, will get error if WITH_GZ is set since option '--gz' did not exist
REGENIE_PATH="$1" 
DOCKER_IMAGE=$2
WITH_GZ=$3

# quick check src/example folders are present
if [ ! -d "${REGENIE_PATH}/src" ] || [ ! -d "${REGENIE_PATH}/example" ]; then
  echo "ERROR: First input argument must be the directory where Regenie repo was cloned"; exit 1
else
  cd $REGENIE_PATH
fi 

# check docker image
if [ -z $DOCKER_IMAGE ]; then
  echo "ERROR: Need to pass docker image tag."; exit 1
elif [[ "$(docker images -q $DOCKER_IMAGE 2> /dev/null)" == "" ]]; then
  echo "ERROR: Image with tag \"${DOCKER_IMAGE}\" does not exist!"; exit 1
fi

# If compiling was done with Boost Iostream library, use gzipped files as input
if [ -z $WITH_GZ ]; then # add suffixes to files
  echo "ERROR: Need to pass indicator for Boost Iostream compilation."; exit 1
elif [ "$WITH_GZ" = "1" ]; then
  fsuf=.gz
  arg_gz="--gz"
fi

# Create test folder to store results and use as mounting point
REGENIE_PATH=$(pwd)/  # use absolute path
# where to mount in container
mntpt=/docker/ 

echo "** Checking docker image \"${DOCKER_IMAGE}\" **"
echo -e "  -> Mounting directory $REGENIE_PATH to /docker/ \n"
echo -e "Running step 1 of REGENIE\n=================================="
# Prepare regenie command to run for Step 1
rgcmd="--step 1 \
  --bed ${mntpt}example/example \
  --exclude ${mntpt}example/snplist_rm.txt \
  --covarFile ${mntpt}example/covariates.txt${fsuf} \
  --phenoFile ${mntpt}example/phenotype_bin.txt${fsuf} \
  --remove ${mntpt}example/fid_iid_to_remove.txt \
  --bsize 100 \
  --bt $arg_gz \
  --lowmem \
  --lowmem-prefix tmp_rg \
  --out ${mntpt}test/fit_bin_out"

docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd

## quick check that the correct files have been created
if [ ! -f ${REGENIE_PATH}test/fit_bin_out.log ] || \
  [ ! -f ${REGENIE_PATH}test/fit_bin_out_pred.list ] || \
  [ ! -f ${REGENIE_PATH}test/fit_bin_out_1.loco$fsuf ] || \
  [ ! -f ${REGENIE_PATH}test/fit_bin_out_2.loco$fsuf ]; then
  echo "Step 1 of REGENIE did not finish successfully. Check the docker image and re-build if needed."; exit 1
elif [ "`grep \"0.456629\" ${REGENIE_PATH}test/fit_bin_out.log | grep \"min value\"`" = "" ]; then
  echo "Step 1 of REGENIE did not finish successfully. Check the docker image and re-build if needed."; exit 1
fi


echo -e "Running step 2 of REGENIE\n=================================="
# First command
rgcmd="--step 2 \
  --bgen ${mntpt}example/example.bgen \
  --with-bgi \
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

docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd

##  do this way so zcat works on OSX
if [ -f ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie.gz ]; then
  ( zcat < ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie.gz ) > ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie
fi

echo "------------------------------------------"
if cmp --silent \
  ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie \
  ${REGENIE_PATH}example/example.test_bin_out_firth_Y1.regenie 
then
  echo -e "Files are identical. Docker image passed 1st test.\n\nSecond test run\n"
else
  echo -e "ERROR: Uh oh... Files are different!\nDocker image did not build successfully."
  exit 1
fi


# Second command
rgcmd="--step 2 \
  --bgen ${mntpt}example/example_3chr.bgen \
  --sample ${mntpt}example/example_3chr.sample \
  --with-bgi \
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

docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd

# check files
echo "------------------------------------------"
if [ ! -f "${REGENIE_PATH}test/test_out_Y2.regenie.ids" -o -f "${REGENIE_PATH}test/test_out_Y1.regenie.ids" ]
then
  echo "Uh oh, docker image did not build successfully"
elif (( $(head -n 1 ${REGENIE_PATH}test/test_out_Y2.regenie.ids | cut -f1) != "Y2" )); then
  echo "Uh oh, docker image did not build successfully"
elif (( $(head -n 1 "${REGENIE_PATH}test/test_out_Y2.regenie.ids" | tr '\t' '\n' | wc -l) != 2 )); then
  echo "Uh oh, docker image did not build successfully"
elif (( `grep "mog_" "${REGENIE_PATH}test/test_out.regenie" | wc -l` > 0 )); then
  echo "Uh oh, docker image did not build successfully"
elif (( `grep "ADD" "${REGENIE_PATH}test/test_out.regenie" | wc -l` > 0 )); then
  echo "Uh oh, docker image did not build successfully"
elif [ "`cut -d ' ' -f1-5 ${REGENIE_PATH}test/test_out.regenie | sed '2q;d'`" != "`grep \"^2\" ${REGENIE_PATH}example/example_3chr.bim | head -n 1 | awk '{print $1,$4,$2,$5,$6}'`" ]; then
  echo "Uh oh, REGENIE did not build successfully"
else

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
    --write-samples \
    --print-pheno \
    --out ${mntpt}test/test_out_extract"

  grep -v "^1" ${REGENIE_PATH}example/example_3chr.bim | awk '{ print $2 }' > ${REGENIE_PATH}test/test_out.snplist

  docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd

  # remove info column from file in run#2
  cut --complement -d ' ' -f7 ${REGENIE_PATH}test/test_out.regenie > ${REGENIE_PATH}test/test_out.regenie.cut

  if cmp --silent \
    ${REGENIE_PATH}test/test_out.regenie.cut \
    ${REGENIE_PATH}test/test_out_extract.regenie 
    then
      echo "SUCCESS: Docker image passed the tests!"
      echo -e "\nYou can run regenie using for example:"
      echo -e "docker run -v <host_path>:<mount_path> $DOCKER_IMAGE regenie <command_options>\n"
    else
      echo "Uh oh, REGENIE did not build successfully. $help_msg"
  fi

fi

# file cleanup
rm ${REGENIE_PATH}test/fit_bin_out* ${REGENIE_PATH}test/test_bin_out_firth* ${REGENIE_PATH}test/test_out*

