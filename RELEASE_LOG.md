## Changelog for past releases

Version 1.0.7 (Enabled for level 0 models in step 1 to be run in parallel [see [Wiki](https://github.com/rgcgithub/regenie/wiki/Further-parallelization-for-level-0-models-in-Step-1) for details]).

Version 1.0.6.9 (Improved step 2 for BGEN format files not in v1.2 or 8-bit encoding).

Version 1.0.6.8 (New option `--range` to specify a chromosome region of variants to test in step 2).

Version 1.0.6.7 (New option `--print-prs` in step 1 to print the whole genome predictions (i.e. PRS) without using LOCO; 
new flag `--use-prs` in step 2 to use these in the association tests).

Version 1.0.6.6 (Fixed MAC calculation for variants on sex chromosomes when sex information is available in the genotype file).

Version 1.0.6.5 (Enabled options `--extract/--exclude` in step 2).

Version 1.0.6.4 (New option `--minINFO` to filter imputed variants in Step 2; added Regenie binary compiled with Intel MKL (only for x86_64 Linux)).

Version 1.0.6.3 (Improved ridge logistic regression to avoid convergence issues in step 1 with low case-count traits).

Version 1.0.6.2 (New option `--ref-first` to use the first allele for each variant as the reference allele for BGEN or PLINK bed/bim/fam file input [default is to use the last allele as the reference]).

Version 1.0.6.1 (Bug fix: Mach R^2 info score is only printed for PGEN input when dosages are present; added flag `--print-pheno` to write the phenotype name in 1st line of sample IDs file [i.e. when using `--write-samples`]).

Version 1.0.6.0 (Improved logistic regression implementation to address convergence issues with low case counts; add new option `--firth-se` to compute SE using effect size estimate and LRT p-value when using Firth correction).


Version 1.0.5.9 (Fixed bug printing variant ID  when variant with variance = 0 occurs in step 1).

Version 1.0.5.8 (Fixed bug due to input genotype file not sorted by chromosome and one of options `--extract/--exclude/--chr/--chrList` is used).

Version 1.0.5.7 (New option `--with-bgi` to read variant information from a .bgi index file for BGEN input format; added option `--write-samples` to write IDs of samples analyzed for each trait in step 2; added Mach Rsq imputation quality metric in INFO column for step 2 with PGEN input file format).

Version 1.0.5.6 (Enabled output of LOCO predictions files and association result files in gzip compressed format using option `--gz` [requires compiling with Boost Iostream library]; added automatic removal from the analysis of genotyped samples in step 2 not present in the LOCO prediction files from step 1 [done separately for each trait]).

Version 1.0.5.5 (fixed bug when setting the total number of blocks [the bug was introduced in v1.0.5.3 due to `n_blocks` being uninitialized]; addressed bug in step 1 with boost filesystem on some machines due to invalid LC_ALL locale) (Note: can now build docker image using `make docker-build`).

Version 1.0.5.4 (Enable using gzip compressed phenotype/covariate files as input [requires installed Boost Iostream library and setting `HAS_BOOST_IOSTREAM = 1` in Makefile] )

Version 1.0.5.31 (Print out ID of problematic variants with low variance in step 1)

Version 1.0.5.3 (Use cxxopts header-only library to parse command line arguments; changed program options `--p/--c/--b/--o/--1` to `--phenoFile/--covarFile/--bsize/--out/--cc12`, respectively; added options `--lowmem-prefix/--pThresh`)

Version 1.0.5.2 (Changed default behavior to remove individuals who have missing data at all phenotypes in the analysis; absolute paths are written in the predictions list file created in step 1)

Version 1.0.5.1 (Reduced memory usage and computational time when using options to keep/remove genotyped samples from the analysis)

Version 1.0.4.2 (Fixed bug excluding/including variants in step 1 with PGEN input format and improved the implementation of how it's done)

Version 1.0.4.1 (Can specify multiple phenotypes/covariates/chromosomes using comma separated arguments; chromosome names can start with 'chr' in the input genotype file)

Version 1.0.4 (Enabled PLINK 2.0 PGEN format files as input using the PLINK 2.0 PGEN library)

Version 1.0.3 (fixed genotype coding in dominant/recessive test for BGEN input format)

Version 1.0.2 (fixed numerical overflow bug when using option `--chr` in step 2; changed to boost split function to read all input files [either space/tab delimited])

Version 1.0.1 (fixed numerical overflow bug for quantile calculation; added new strategy for fitting null model for approximate Firth test) 

Version 1.0 (22 June 2020): Initial release



