#########################################
## script used to run regenie/BOLT/fastGWA/SAIGE
## for the analyses in the REGENIE 2020 paper
## For more details, visit: https://rgcgithub.github.io/regenie/
##

main_dir=/regenie_paper/
code_dir=${main_dir}scripts/
in_raw_dir=${main_dir}input_raw/
files_dir=${main_dir}input/
out_dir=${main_dir}output/
fig_dir=${main_dir}figures/
mkdir -p $files_dir $out_dir $fig_dir


## Input files needed
## * UKB 500K array data BED file
ukb_array_prefix=${in_raw_dir}/ukb_array
## * UKB 500K HRC imputed data BGEN files split by chromosome (22)
###   - only variants which have minor allele frequency above 0.5% 
###      or have minor allele count above 5 and are annotated as functional
ukb_imputed_prefix=${in_raw_dir}/ukb_imputed_chr
## * List of regions to exclude for step 1 (ICLD and low-complexity regions) 
regions_ignore=${in_raw_dir}/regions.exclude
## * List of white British ancestry samples (FID/IID) 
wb_samples=${in_raw_dir}/wb.samples
## * UKB 500K phenotype files for the QTs and BTs analyzed in the Regenie paper
###   -format is be FID IID followed by the phenotypes analyzed
###   - 4 sets of files for 3 exemplar QTs, 50 multi-trait QTs, 4 exemplar BTs and 50 multi-trait BTs
ukb_pheno_exQT=${files_dir}/ukb_phenos_exQTs.txt
ukb_pheno_mtQT=${files_dir}/ukb_phenos_mtQTs.txt
ukb_pheno_exBT=${files_dir}/ukb_phenos_exBTs.txt
ukb_pheno_mtBT=${files_dir}/ukb_phenos_mtBTs.txt
## * UKB 500K covariate file
###   -format is be FID IID followed by the covariates included in paper
ukb_covars=${files_dir}/ukb_cov.txt


sudo chmod u+x ${code_dir}*
# setup imagequick for tiff to pdf conversion
sudo apt install -y imagemagick-6.q16
cp /etc/ImageMagick-6/policy.xml /etc/ImageMagick-6/policy.xml.cpy
line_change=`grep -n "coder.*PDF"  /etc/ImageMagick-6/policy.xml | head -n 1 | cut -f1 -d':'`
head -n $(( line_change - 1 )) /etc/ImageMagick-6/policy.xml > tmp.policy
sed "${line_change}q;d" /etc/ImageMagick-6/policy.xml | sed 's/none/write/' >> tmp.policy
tail -n +$(( line_change + 1 )) /etc/ImageMagick-6/policy.xml >> tmp.policy
mv tmp.policy /etc/ImageMagick-6/policy.xml
# install parallel for SAIGE step 2
sudo apt install -y parallel


##########
## Prepare Step 1 files
#########

### Apply filters for Step 1 file of array SNPs
array_step1=${files_dir}/ukb_wb_array_step1

${code_dir}prep_files.sh \
  $ukb_array_prefix \
  ${ukb_imputed_prefix}1 \
  $wb_samples \
  $regions_ignore \
  $array_step1

### Create sparse GRM for fastGWA
npartitions=250
for job in $(seq 1 $npartitions); do

  ${code_dir}mk_sparseGRM.r \
    --step1File=$array_step1 \
    --partition=$job \
    --npartitions=$npartitions \
    --prefix=${array_step1}_fastgwa

done

${code_dir}mk_sparseGRM.r \
  --npartitions=$npartitions \
  --prefix=${array_step1}_fastgwa

# Make file with phenotype names
for f in $ukb_pheno_exQT $ukb_pheno_mtQT $ukb_pheno_exBT $ukb_pheno_mtBT ; do 
  head -n 1 $f | tr ' ' '\n' | grep -v "FID\|IID" > ${f}.names
done


##########
## Run GWAS for each data set
#########
nchr=22

##################
### 1. exemplarQTs (analyze one trait at a time)
oprefix=${out_dir}exQTs
npheno=`cat ${ukb_pheno_exQT}.names | wc -l`

for ipheno in $(seq 1 $npheno); do

  ## regenie step 1
  ${code_dir}run_methods.r \
    --method=regenie \
    --phenoFile=$ukb_pheno_exQT \
    --pheno=$ipheno \
    --covarFile=$ukb_covars \
    --step1File=$array_step1 \
    --lowmem \
    --prefix=${oprefix}_regenie

  ## regenie step 2
  for chr in $(seq 1 $nchr); do

    ${code_dir}run_methods.r \
      --method=regenie \
      --skipNull \
      --lowmem \
      --phenoFile=$ukb_pheno_exQT \
      --pheno=$ipheno \
      --covarFile=$ukb_covars \
      --step2File=$ukb_imputed_prefix \
      --chr $chr \
      --prefix=${oprefix}_regenie

    done

  # bolt-lmm
  ${code_dir}run_methods.r \
    --method=bolt \
    --phenoFile=$ukb_pheno_exQT \
    --pheno=$ipheno \
    --covarFile=$ukb_covars \
    --step1File=$array_step1 \
    --step2File=$ukb_imputed_prefix \
    --prefix=${oprefix}_bolt

  # fastgwa
  ${code_dir}run_methods.r \
    --method=fastgwa \
    --phenoFile=$ukb_pheno_exQT \
    --pheno=$ipheno \
    --covarFile=$ukb_covars \
    --grm=${array_step1}_fastgwa_sp \
    --step2File=$ukb_imputed_prefix \
    --prefix=${oprefix}_fastgwa

  done


## make ManP plot
${code_dir}mk_plots_qt.r \
  --loadFuns=${code_dir}std_ffuns.r \
  --manColors=${code_dir}manP.colors \
  --phenoNames=${ukb_pheno_exQT}.names \
  --figfolder=$fig_dir \
  --prefix=$oprefix



##################
### 2. exemplarBTs (analyze one trait at a time)
oprefix=${out_dir}exBTs
npheno=`cat ${ukb_pheno_exBT}.names | wc -l`

for ipheno in $(seq 1 $npheno); do

  ## regenie step 1
  ${code_dir}run_methods.r \
    --method=regenie \
    --bt \
    --phenoFile=$ukb_pheno_exBT \
    --pheno=$ipheno \
    --covarFile=$ukb_covars \
    --step1File=$array_step1 \
    --lowmem \
    --prefix=${oprefix}_regenie


  ## regenie step 2
  for chr in $(seq 1 $nchr); do

    # firth
    ${code_dir}run_methods.r \
      --method=regenie \
      --bt \
      --skipNull \
      --lowmem \
      --phenoFile=$ukb_pheno_exBT \
      --pheno=$ipheno \
      --covarFile=$ukb_covars \
      --step2File=$ukb_imputed_prefix \
      --chr $chr \
      --prefix=${oprefix}_regenie


    # spa
    ${code_dir}run_methods.r \
      --method=regenie \
      --bt --spa \
      --skipNull \
      --lowmem \
      --phenoFile=$ukb_pheno_exBT \
      --pheno=$ipheno \
      --covarFile=$ukb_covars \
      --step2File=$ukb_imputed_prefix \
      --chr $chr \
      --prefix=${oprefix}_regenie

  done

  # bolt-lmm
  ${code_dir}run_methods.r \
    --method=bolt \
    --phenoFile=$ukb_pheno_exBT \
    --pheno=$ipheno \
    --covarFile=$ukb_covars \
    --step1File=$array_step1 \
    --step2File=$ukb_imputed_prefix \
    --prefix=${oprefix}_bolt


  # saige step 1
  ${code_dir}run_methods.r \
    --method=saige \
    --phenoFile=$ukb_pheno_exBT \
    --pheno=$ipheno \
    --covarFile=$ukb_covars \
    --step1File=$array_step1 \
    --prefix=${oprefix}_saige

  # saige step 2 (run in parallel since saige only uses ~1 thread here)
  parallel --jobs $(nproc) --link \
    ${code_dir}run_methods.r \
    --method=saige \
    --skipNull \
    --phenoFile=$ukb_pheno_exBT \
    --pheno=$ipheno \
    --covarFile=$ukb_covars \
    --step2File=$ukb_imputed_prefix \
    --prefix=${oprefix}_saige \
    --chr \
    ::: $(seq 1 $nchr)

done

## make ManP plot
${code_dir}mk_plots_bt.r \
  --loadFuns=${code_dir}std_ffuns.r \
  --manColors=${code_dir}manP.colors \
  --phenoNames=${ukb_pheno_exBT}.names \
  --figfolder=$fig_dir \
  --prefix=$oprefix



##################
### 3. multi-trait QTs
oprefix=${out_dir}mtQTs
npheno=`cat ${ukb_pheno_mtQT}.names | wc -l`

## regenie step 1 (all traits analyzed) K-fold CV
${code_dir}run_methods.r \
  --method=regenie \
  --phenoFile=$ukb_pheno_mtQT \
  --covarFile=$ukb_covars \
  --step1File=$array_step1 \
  --lowmem \
  --prefix=${oprefix}_regenie

## regenie step 1 LOOCV
${code_dir}run_methods.r \
  --method=regenie \
  --phenoFile=$ukb_pheno_mtQT \
  --covarFile=$ukb_covars \
  --step1File=$array_step1 \
  --loocv \
  --lowmem \
  --prefix=${oprefix}_regenie

## regenie step 2
for chr in $(seq 1 $nchr); do

  ${code_dir}run_methods.r \
    --method=regenie \
    --skipNull \
    --lowmem \
    --phenoFile=$ukb_pheno_mtQT \
    --covarFile=$ukb_covars \
    --step2File=$ukb_imputed_prefix \
    --chr $chr \
    --prefix=${oprefix}_regenie

done

for ipheno in $(seq 1 $npheno); do

  # bolt-lmm
  ${code_dir}run_methods.r \
    --method=bolt \
    --phenoFile=$ukb_pheno_mtQT \
    --pheno=$ipheno \
    --covarFile=$ukb_covars \
    --step1File=$array_step1 \
    --step2File=$ukb_imputed_prefix \
    --prefix=${oprefix}_bolt

  # fastgwa
  ${code_dir}run_methods.r \
    --method=fastgwa \
    --phenoFile=$ukb_pheno_mtQT \
    --pheno=$ipheno \
    --covarFile=$ukb_covars \
    --grm=${array_step1}_fastgwa_sp \
    --step2File=$ukb_imputed_prefix \
    --prefix=${oprefix}_fastgwa

done


##################
### 4. multi-trait BTs
oprefix=${out_dir}mtBTs
npheno=`cat ${ukb_pheno_mtBT}.names | wc -l`

## regenie step 1 (all traits analyzed at once) K-fold CV
${code_dir}run_methods.r \
  --method=regenie \
  --bt \
  --phenoFile=$ukb_pheno_mtBT \
  --covarFile=$ukb_covars \
  --step1File=$array_step1 \
  --lowmem \
  --prefix=${oprefix}_regenie

## regenie step 1 LOOCV
${code_dir}run_methods.r \
  --method=regenie \
  --bt \
  --loocv \
  --phenoFile=$ukb_pheno_mtBT \
  --covarFile=$ukb_covars \
  --step1File=$array_step1 \
  --lowmem \
  --prefix=${oprefix}_regenie

## regenie step 2
for chr in $(seq 1 $nchr); do

  # firth
  ${code_dir}run_methods.r \
    --method=regenie \
    --bt \
    --skipNull \
    --lowmem \
    --phenoFile=$ukb_pheno_mtBT \
    --covarFile=$ukb_covars \
    --step2File=$ukb_imputed_prefix \
    --chr $chr \
    --prefix=${oprefix}_regenie


  # spa
  ${code_dir}run_methods.r \
    --method=regenie \
    --bt --spa \
    --skipNull \
    --lowmem \
    --phenoFile=$ukb_pheno_mtBT \
    --covarFile=$ukb_covars \
    --step2File=$ukb_imputed_prefix \
    --chr $chr \
    --prefix=${oprefix}_regenie

done

for ipheno in $(seq 1 $npheno); do

  # bolt-lmm
  ${code_dir}run_methods.r \
    --method=bolt \
    --phenoFile=$ukb_pheno_mtBT \
    --pheno=$ipheno \
    --covarFile=$ukb_covars \
    --step1File=$array_step1 \
    --step2File=$ukb_imputed_prefix \
    --prefix=${oprefix}_bolt


  # saige step 1
  ${code_dir}run_methods.r \
    --method=saige \
    --phenoFile=$ukb_pheno_mtBT \
    --pheno=$ipheno \
    --covarFile=$ukb_covars \
    --step1File=$array_step1 \
    --prefix=${oprefix}_saige

  # saige step 2 (run in parallel since saige only uses ~1 thread here)
  parallel --jobs $(nproc) --link \
    ${code_dir}run_methods.r \
    --method=saige \
    --skipNull \
    --phenoFile=$ukb_pheno_mtBT \
    --pheno=$ipheno \
    --covarFile=$ukb_covars \
    --step2File=$ukb_imputed_prefix \
    --prefix=${oprefix}_saige \
    --chr \
    ::: $(seq 1 $nchr)

done


# reset policies for imagemagick
mv /etc/ImageMagick-6/policy.xml.cpy /etc/ImageMagick-6/policy.xml
