#!/usr/bin/env bash

#### Working from directory where regenie repo was cloned
if [ "$#" -eq 0 ]; then
echo "Usage: test.sh <PATH_TO_CLONED_REGENIE_REPO> <DOCKER_IMAGE_TAG>"; exit 1
fi

### Test script for Regenie version >= 1.0.5.6
## In previous versions, will get error if WITH_GZ is set since option '--gz' did not exist
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
elif (( $WITH_GZ == 1 )); then
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
fi


echo -e "Running step 2 of REGENIE\n=================================="
# Prepare regenie command to run for Step 2
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

docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd


## check that compilation was successful
if [ -f ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie.gz ]; then
  # uncompress file (zcat this way should work on OSX)
  ( zcat < ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie.gz ) > ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie
fi

# compare result files to exemplar file
echo "------------------------------------------"
if cmp --silent \
  ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie \
  ${REGENIE_PATH}example/example.test_bin_out_firth_Y1.regenie 
then
  echo "SUCCESS: Docker image passed the tests!"
  echo -e "\nYou can run regenie using:"
  echo -e "docker run -v <host_path>:<mount_path> $DOCKER_IMAGE regenie <command_options>\n"
else
  echo -e "ERROR: Uh oh... Files are different!\nDocker image did not build successfully"
fi

# file cleanup
rm ${REGENIE_PATH}test/fit_bin_out* ${REGENIE_PATH}test/test_bin_out_firth*

