#!/usr/bin/env bash

### REGENIE TEST SCRIPT TO USE WITH DOCKER IMAGE ###
# Functions used
help_msg="Check the docker image and re-build if needed."
err_msg="Docker image did not build successfully."
print_err () { 
  echo "$err_msg"; exit 1 
}
print_custom_err () {
  echo "${1} $help_msg"; exit 1 
}

### READ OPTIONS
if [ "$#" -eq 0 ]; then
  echo "Usage: test.sh <PATH_TO_CLONED_REGENIE_REPO> <DOCKER_IMAGE_TAG>"; exit 1
fi

REGENIE_PATH="$1" 
DOCKER_IMAGE=$2
WITH_GZ=$3

# quick check src/example folders are present
if [ ! -d "${REGENIE_PATH}/src" ] || [ ! -d "${REGENIE_PATH}/example" ]; then
  echo "ERROR: First input argument must be the directory where Regenie repo was cloned"; exit 1
fi 
cd $REGENIE_PATH

# check docker image
if [ -z $DOCKER_IMAGE ]; then
  echo "ERROR: Need to pass docker image tag."; exit 1
elif [[ "$(docker images -q $DOCKER_IMAGE 2> /dev/null)" == "" ]]; then
  echo "ERROR: Image with tag \"${DOCKER_IMAGE}\" does not exist!"; exit 1
fi

# If compiling was done with Boost Iostream library, use gzipped files as input
if [ "$WITH_GZ" = "1" ]; then
  fsuf=.gz
  arg_gz="--gz"
fi

# Create test folder to store results and use as mounting point
REGENIE_PATH=$(pwd)/  # use absolute path
# where to mount in container
mntpt=/docker/ 

echo "** Checking docker image \"${DOCKER_IMAGE}\" **"
echo -e "  -> Mounting directory $REGENIE_PATH to /docker/ \n"


echo -e "==>Running step 1 of REGENIE"
# Prepare regenie command to run for Step 1
fail_msg="Step 1 of REGENIE did not finish successfully."
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

docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd

## quick check that the correct files have been created
if [ ! -f "${REGENIE_PATH}test/fit_bin_out.log" ] || \
  [ ! -f "${REGENIE_PATH}test/fit_bin_out_pred.list" ] || \
  [ ! -f "${REGENIE_PATH}test/fit_bin_out_1.loco$fsuf" ] || \
  [ ! -f "${REGENIE_PATH}test/fit_bin_out_2.loco$fsuf" ]; then
  print_custom_err "$fail_msg"
elif [ "`grep \"0.4504\" ${REGENIE_PATH}test/fit_bin_out.log | grep \"min value\"`" = "" ]; then
  print_custom_err "$fail_msg"
fi


#### Run step 1 splitting across jobs for level 0
njobs=4
echo -e "==>Re-running step 1 splitting in $njobs jobs"
# pt1 - run regenie before l0
rgcmd="$basecmd \
  --split-l0 ${mntpt}test/fit_bin_parallel,$njobs \
  --out ${mntpt}test/fit_bin_l0"

docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd
if [ ! -f "${REGENIE_PATH}test/fit_bin_parallel.master" ]; then
  print_custom_err "$fail_msg"
fi

# pt2 - run regenie for l0
nj=`seq 1 $njobs`
for job in $nj; do
  rgcmd="$basecmd \
    --run-l0 ${mntpt}test/fit_bin_parallel.master,$job \
    --out ${mntpt}test/fit_bin_l0"

  docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd
  if [ ! -f "${REGENIE_PATH}test/fit_bin_parallel_job${job}_l0_Y1" ]; then
    print_custom_err "$fail_msg"
  fi
done


# pt3 - run regenie for l1
rgcmd="$basecmd \
  --run-l1 ${mntpt}test/fit_bin_parallel.master \
  --out ${mntpt}test/fit_bin_l1"

docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd

if [ ! -f "${REGENIE_PATH}test/fit_bin_l1_1.loco$fsuf" ]; then
  print_custom_err "$fail_msg"
elif ! cmp --silent \
  "${REGENIE_PATH}test/fit_bin_out_1.loco$fsuf" \
  "${REGENIE_PATH}test/fit_bin_l1_1.loco$fsuf" 
then
  print_custom_err "$fail_msg"
elif ! cmp --silent \
  "${REGENIE_PATH}test/fit_bin_out_2.loco$fsuf" \
  "${REGENIE_PATH}test/fit_bin_l1_2.loco$fsuf" 
then
  print_custom_err "$fail_msg"
fi



##########
##########
#### Step 2
i=1
echo -e "==>Running step 2 of REGENIE; test #$i"
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
  $arg_gz \
  --out ${mntpt}test/test_bin_out_firth"

docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd

##  do this way so zcat works on OSX
if [ -f ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie.gz ]; then
  ( zcat < ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie.gz ) > ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie
fi

if ! cmp --silent \
  ${REGENIE_PATH}test/test_bin_out_firth_Y1.regenie \
  ${REGENIE_PATH}example/example.test_bin_out_firth_Y1.regenie 
then
  print_custom_err "ERROR: Uh oh... Files are different!"
fi


(( i++ ))
echo -e "Files are identical.\n\n==>Running test #$i\n"
# Next test
basecmd="--step 2 \
  --bed ${mntpt}example/example_3chr \
  --ref-first \
  --covarFile ${mntpt}example/covariates.txt${fsuf} \
  --phenoFile ${mntpt}example/phenotype_bin.txt${fsuf} \
  --phenoColList Y2 \
  --bsize 100 \
  --test dominant \
  --ignore-pred"
rgcmd="$basecmd \
  --chrList 2,3 \
  --write-samples \
  --print-pheno \
  --out ${mntpt}test/test_out"

docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd

# check files
if [ ! -f "${REGENIE_PATH}test/test_out_Y2.regenie.ids" -o -f "${REGENIE_PATH}test/test_out_Y1.regenie.ids" ]
then
  print_err
elif (( $(head -n 1 ${REGENIE_PATH}test/test_out_Y2.regenie.ids | cut -f1) != "Y2" )); then
  print_err
elif (( $(head -n 1 "${REGENIE_PATH}test/test_out_Y2.regenie.ids" | tr '\t' '\n' | wc -l) != 2 )); then
  print_err
elif (( `grep "mog_" "${REGENIE_PATH}test/test_out_Y2.regenie" | wc -l` > 0 )); then
  print_err
elif (( `grep "ADD" "${REGENIE_PATH}test/test_out_Y2.regenie" | wc -l` > 0 )); then
  print_err
elif [ "`cut -d ' ' -f1-5 ${REGENIE_PATH}test/test_out_Y2.regenie | sed '2q;d'`" != "`grep \"^2\" ${REGENIE_PATH}example/example_3chr.bim | head -n 1 | awk '{print $1,$4,$2,$5,$6}'`" ]; then
  print_err
fi

(( i++ ))
echo -e "==>Running test #$i"
# Next test
rgcmd="$basecmd \
  --extract ${mntpt}test/test_out.snplist \
  --out ${mntpt}test/test_out_extract"

awk '{if($1!=1) {print $2}}'  ${REGENIE_PATH}example/example_3chr.bim > ${REGENIE_PATH}test/test_out.snplist

docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd

if ! cmp --silent \
  ${REGENIE_PATH}test/test_out_Y2.regenie \
  ${REGENIE_PATH}test/test_out_extract_Y2.regenie 
then
  print_err
fi

(( i++ ))
echo -e "==>Running test #$i"
# First command (V1)
rgcmd="--step 2 \
  --bed ${mntpt}example/example_3chr_masks \
  --covarFile ${mntpt}example/covariates.txt${fsuf} \
  --phenoFile ${mntpt}example/phenotype_bin.txt${fsuf} \
  --remove ${mntpt}example/fid_iid_to_remove.txt \
  --bsize 10 \
  --ignore-pred \
  --htp TEST \
  --out ${mntpt}test/test_out_masks_V1"
# run regenie
docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd

# Second command (V2)
# build masks
awk '{print $4}' ${REGENIE_PATH}example/example_3chr.setlist | tr ',' '\n' > ${REGENIE_PATH}test/tmp1.txt 
rgcmd="--step 2 \
  --ignore-pred \
  --bed ${mntpt}example/example_3chr \
  --extract ${mntpt}test/tmp1.txt \
  --covarFile ${mntpt}example/covariates.txt${fsuf} \
  --phenoFile ${mntpt}example/phenotype_bin.txt${fsuf} \
  --remove ${mntpt}example/fid_iid_to_remove.txt \
  --set-list ${mntpt}example/example_3chr.setlist \
  --anno-file ${mntpt}example/example_3chr.annotations \
  --mask-def ${mntpt}example/example_3chr.masks \
  --write-mask \
  --write-setlist ${mntpt}example/example_3chr.write_sets \
  --bsize 15 \
  --aaf-bins 0.2 \
  --chrList 1,3 \
  --htp TEST \
  --out ${mntpt}test/test_out_masks_V2"

# run regenie
docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd

head ${REGENIE_PATH}test/test_out_masks_V2_Y1.regenie -n 3 | tail -n 2 | cut --complement -f4,5 > ${REGENIE_PATH}test/tmp1.txt
tail -n 1 ${REGENIE_PATH}test/test_out_masks_V2_Y1.regenie | cut --complement -f4,5 >> ${REGENIE_PATH}test/tmp1.txt
cat ${REGENIE_PATH}test/test_out_masks_V1_Y1.regenie | cut --complement -f4,5 > ${REGENIE_PATH}test/tmp2.txt

if ! cmp --silent \
  ${REGENIE_PATH}test/tmp1.txt \
  ${REGENIE_PATH}test/tmp2.txt ; then
  print_err
elif [ ! -f ${REGENIE_PATH}test/test_out_masks_V2_masks.bed ]; then
  print_err
elif [ "$(hexdump -e \"%07_ax\ \"\ 16/1\ \"\ %02x\"\ \"\\n\"  -n 3 ${REGENIE_PATH}test/test_out_masks_V2_masks.bed | head -n 1 | awk '{print $2,$3,$4}' | tr ' ' ',')" != "6c,1b,01" ]; then
  print_err
elif [ "`wc -l ${REGENIE_PATH}test/test_out_masks_V2_masks.{bim,fam} | awk '{print $1}' | head -n 2| paste -sd','`" != "4,494" ]; then
  print_err
elif [ "`cat ${REGENIE_PATH}test/test_out_masks_V2_tmp2.setlist | head -n 1 | tr ',' '\n' | wc -l`" != "2" ]; then
  print_err
fi


(( i++ ))
echo -e "==>Running test #$i"
# build masks
awk '{print $4}' ${REGENIE_PATH}example/example_3chr.setlist | tr ',' '\n' > ${REGENIE_PATH}test/tmp1.txt 
rgcmd="--step 2 \
  --ignore-pred \
  --bed ${mntpt}example/example_3chr \
  --extract ${mntpt}test/tmp1.txt \
  --covarFile ${mntpt}example/covariates.txt${fsuf} \
  --phenoFile ${mntpt}example/phenotype_bin.txt${fsuf} \
  --set-list ${mntpt}example/example_3chr.setlist \
  --anno-file ${mntpt}example/example_3chr.annotations \
  --mask-def ${mntpt}example/example_3chr.masks \
  --mask-lovo SET1,M1,0.2 \
  --htp TEST \
  --out ${mntpt}test/test_out_masks_loo"

# run regenie
docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd

if [ ! -f ${REGENIE_PATH}test/test_out_masks_loo_Y1.regenie ]; then
  print_err
elif [ `cat ${REGENIE_PATH}test/test_out_masks_loo_Y1.regenie | wc -l` != 20 ]; then
  print_err
elif [ `grep "_mog" ${REGENIE_PATH}test/test_out_masks_loo_Y1.regenie | wc -l` != 18 ]; then
  print_err
fi


(( i++ ))
echo -e "==>Running test #$i"
# build masks using set domains
rgcmd="--step 2 \
  --ignore-pred \
  --bed ${mntpt}example/example_3chr \
  --covarFile ${mntpt}example/covariates.txt${fsuf} \
  --phenoFile ${mntpt}example/phenotype_bin.txt${fsuf} \
  --remove ${mntpt}example/fid_iid_to_remove.txt \
  --set-list ${mntpt}example/example_3chr.setlist \
  --anno-file ${mntpt}example/example_3chr.annotationsV2 \
  --mask-def ${mntpt}example/example_3chr.masks \
  --bsize 20 \
  --aaf-bins 0.2 \
  --out ${mntpt}test/test_out_masks_V3"

# run regenie
docker run -v ${REGENIE_PATH}:${mntpt} --rm $DOCKER_IMAGE regenie $rgcmd

if ! [[ "`head -n 1 ${REGENIE_PATH}test/test_out_masks_V3_Y1.regenie`" =~ ^\#\#MASKS.* ]]
then
  print_err
elif [ `grep "SET2.*.M1" ${REGENIE_PATH}test/test_out_masks_V3_Y1.regenie | wc -l` != "4" ]
then
  print_err
fi

echo "SUCCESS: Docker image passed the tests!"
echo -e "\nYou can run regenie using for example:"
echo -e "docker run -v <host_path>:<mount_path> $DOCKER_IMAGE regenie <command_options>\n"
# file cleanup
rm ${REGENIE_PATH}test/fit_bin_* ${REGENIE_PATH}test/test_bin_out_firth* ${REGENIE_PATH}test/test_out* ${REGENIE_PATH}test/tmp[12].txt

