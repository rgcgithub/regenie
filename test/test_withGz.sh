#!/usr/bin/env bash

#### Working from directory where regenie repo was cloned
REGENIE_PATH="$1" 
if [ "$#" -eq 0 ]; then
echo "Usage: test_withGz.sh <PATH_TO_CLONED_REGENIE_REPO>"; exit 1
fi

# quick check src/example folders are present
if [ ! -d "${REGENIE_PATH}/src" ] || [ ! -d "${REGENIE_PATH}/example" ]; then
  echo "ERROR: Input argument must be the directory where Regenie repo was cloned"; exit 1
fi 

cd $REGENIE_PATH


echo -e "Building docker image\n=================================="
dckfile=Dockerfile_gz
docker build -f $dckfile -t regenie_gz:latest .


# Create test folder to store results and use as mounting point
# make sure to use absolute path
REGENIE_PATH=$(pwd)/
# where to mount in container
mntpt=/docker/ 

echo -e "\n\nRunning step 1 of regenie with gzipped-files\n=================================="
# Prepare regenie command to run for Step 1
rgcmd="--step 1 \
  --bed ${mntpt}example/example \
  --exclude ${mntpt}example/snplist_rm.txt \
  --covarFile ${mntpt}example/covariates.txt.gz \
  --phenoFile ${mntpt}example/phenotype_bin.txt.gz \
  --remove ${mntpt}example/fid_iid_to_remove.txt \
  --bsize 100 \
  --bt --lowmem \
  --lowmem-prefix tmp_rg \
  --out ${mntpt}test/fit_bin_out"

docker run -v ${REGENIE_PATH}:${mntpt} --rm regenie_gz:latest $rgcmd

echo -e "Running step 2 of regenie with gzipped-files\n=================================="
# Prepare regenie command to run for Step 2
rgcmd="--step 2 \
  --bgen ${mntpt}example/example.bgen \
  --covarFile ${mntpt}example/covariates.txt.gz \
  --phenoFile ${mntpt}example/phenotype_bin.txt.gz \
  --remove ${mntpt}example/fid_iid_to_remove.txt \
  --bsize 200 \
  --bt \
  --firth --approx \
  --pThresh 0.01 \
  --pred ${mntpt}test/fit_bin_out_pred.list \
  --split \
  --out ${mntpt}test/test_bin_out_firth"

docker run -v ${REGENIE_PATH}:${mntpt} --rm regenie_gz:latest $rgcmd

## quick check that compilation was successful
if cmp --silent \
  ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie \
  ${REGENIE_PATH}example/example.test_bin_out_firth_Y1.regenie 
then
   echo "SUCCESS: Files are identical!"
   echo -e "\nYou can run regenie using:"
   echo "docker run -v <host_path>:<mount_path> --rm regenie_gz:latest <command_options>"
 else
   echo "Uh oh... Files are different!"
fi

# file cleanup
rm ${REGENIE_PATH}test/fit_bin_out* ${REGENIE_PATH}test/test_bin_out_firth*

