## Documentation
## Getting started

To run **regenie**, use the command `./regenie` on the command line,
followed by options and flags as needed.

To get a full list of options use

```
./regenie --help
```

The directory `examples/` contains some small example files that are
useful when getting started. A test run on a set of binary traits can be achieved by the
following 2 commands.

In **Step 1** the whole genome regression model is fit to the traits, and
a set of genomic predictions are produced as output

```
./regenie \
  --step 1 \
  --bed example/example \
  --exclude example/snplist_rm.txt \
  --covarFile example/covariates.txt \
  --phenoFile example/phenotype_bin.txt \
  --remove example/fid_iid_to_remove.txt \
  --bsize 100 \
  --bt --lowmem \
  --lowmem-prefix tmp_rg \
  --out fit_bin_out
```

In **Step 2** a set of imputed SNPs are tested for association using a
Firth logistic regression model

```
./regenie \
  --step 2 \
  --bgen example/example.bgen \
  --covarFile example/covariates.txt \
  --phenoFile example/phenotype_bin.txt \
  --remove example/fid_iid_to_remove.txt \
  --bsize 200 \
  --bt \
  --firth --approx \
  --pThresh 0.01 \
  --pred fit_bin_out_pred.list \
  --split \
  --out test_bin_out_firth
```

One of the output files from these two commands is included in the `example/` directory and 
you can check they are the same using  the following command (the message should be printed out)

```
cmp test_bin_out_firth_Y1.regenie example/example.test_bin_out_firth_Y1.regenie \
  && echo "Files are identical"
```


## Input 


| Option | Argument | Type | Description|
|---|-------|------|----|
|`--bgen, --bed, --pgen`  | FILE | Required |Input genetic data file. Either BGEN file eg. `file.bgen`, or bed/bim/fam prefix that assumes`file.bed`, `file.bim`, `file.fam` exist, or pgen/pvar/psam prefix that assumes`file.pgen`, `file.pvar`, `file.psam` exist |
|`--sample`  | FILE | Optional |Sample file corresponding to input BGEN file|
|`--with-bgi`  | FLAG | Optional |Specify to use accompanying bgi index file to get variant information (assumes file name is `file.bgen.bgi`|
|`--keep`  | FILE | Optional | Inclusion file that lists individuals to retain in the analysis|
|`--remove`  | FILE | Optional | Exclusion file that lists individuals to remove from the analysis|
|`--extract`  | FILE | Optional | Inclusion file that lists IDs of variants to keep **(only works with option `--step 1`)**|
|`--exclude`  | FILE | Optional | Exclusion file that lists IDs of variants to remove **(only works with option `--step 1`)**|
|`--phenoFile`  | FILE | Required |Phenotypes file|
|`--phenoCol` | STRING | Optional | Use for each phenotype you want to include in the analysis|
|`--phenoColList` | STRING | Optional | Comma separated list of phenotypes to include in the analysis|
|`--covarFile`  | FILE | Optional | Covariates file|
|`--covarCol` | STRING | Optional | Use for each covariate you want to include in the analysis|
|`--covarColList` | STRING | Optional | Comma separated list of covariates to include in the analysis|
|`--pred`  | FILE | Optional  | File containing predictions from Step 1 (see Overview). **This is required for `--step 2`**|


#### Genetic data file format

**regenie** can read BGEN files, bed/bim/fam files or pgen/psam/pvar 
files in Step 1 and Step 2.

The BGEN file format is described
[here](https://www.well.ox.ac.uk/~gav/bgen_format/).

The bed/bim/fam file format is described [here](https://www.cog-genomics.org/plink/1.9/formats).

The pgen/pvar/psam file format is described [here](https://www.cog-genomics.org/plink/2.0/formats#pgen).

Tools useful for genetic data file format conversion are : [PLINK](http://www.cog-genomics.org/plink/), [QCTOOL](https://www.well.ox.ac.uk/~gav/qctool/), [BCFTOOLS](https://samtools.github.io/bcftools/).

Step 2 of **regenie** can be sped up with BGEN files by using v1.2 format with 8 bits encoding 
(genotype file can be generated with [PLINK2](https://www.cog-genomics.org/plink/2.0/) using 
option `--export bgen-1.2 'bits=8'`) as well as passing on an accompanying .bgi index file 
(a useful tool to create such file is bgenix which is part of the BGEN library).

To include X chromosome genotypes in step 1 and/or step 2, males should be coded as diploid 
so that their genotypes are 0/2. This can be done in PLINK by setting the sex of all 
individuals to female before generating the genotype file.
Chromosome values of 23 (for human analyses), X, XY, PAR1 and PAR2 are all acceptable and 
will be collapsed into a single chromosome.


#### Sample inclusion/exclusion file format

```
2 2 
7 7 
.
```

No header. Each line starts with individual FID IID. Space/tab separated.

Samples listed in the file that are not in bgen/bed/pgen file are ignored.

#### Variant inclusion/exclusion file format

```
20
31
.
```

No header. Each line must start with variant ID 
(if there are additional columns, file must be space/tab separated).

Variants listed in this file that are not in bgen/bed/pgen file are ignored.

#### Covariate file format

```
FID IID V1 V2 V3
1 1 1.46837294454993 1.93779743016325 0.152887004505393
2 2 -1.2234390803815 -1.63408619199948 -0.190201446835255
3 3 0.0711531925667286 0.0863906292357564 0.14254739715665
.
```

Line 1 : Header with FID, IID and \(C\) covariate names.

Followed by lines of \(C+2\) values. Space/tab separated.

Each line contains individual FID and IID followed by \(C\) covariate
values.

Samples listed in this file that are not in bgen/bed/pgen file are ignored.
Genotyped samples that are not in this file are removed from the analysis.

No missing values are allowed.

If `--step 2` is specified, then the covariate file should be the same
as that used in Step 1.
#### Phenotype file format

```
FID IID Y1 Y2
1 1 1.64818554321186 2.2765234736685
2 2 -2.67352013711554 -1.53680421614647
3 3 0.217542851471485 0.437289912695016
.
```

Line 1 : Header with FID, IID and \(P\) phenotypes names.

Followed by lines of \(P+2\) values. Space/tab separated. 
Each line contains individual FID and IID followed by P phenotype values
(for binary traits, must be coded as 0=control, 1=case, NA=missing unless using `--1`).

Samples listed in this file that are not in bgen/bed/pgen file are ignored.
Genotyped samples that are not in this file are removed from the analysis.

Missing values must be coded as NA.

With QTs, missing values are mean-imputed in Step 1 and they are dropped when testing each phenotype in Step 2 (unless using `--force-impute`).

With BTs, missing values are mean-imputed in Step 1 when fitting the
level 0 linear ridge regression and 
they are dropped when fitting the level 1 logistic ridge regression for each trait . 
In Step 2, missing values are dropped when testing each trait.

To remove all samples that have missing values at **any** of the \(P\) phenotypes, use option `--strict` in Step 1 and 2. This is also useful when analyzing a single trait to avoid making a new bed/pgen/bgen file just for the complete data set of individuals (so setting the phenotype values of individuals to remove to NA), although `--remove` can also be used in that situation.

#### Predictions file format

Running `--step 1 --out foo` will produce

1. A set of files containing genomic predictions for each phenotype
   from Step 1 (see Output section below).
2. A file called `foo_pred.list` listing the locations of the prediction files.

The file list is needed as an input file when using `--step 2`
via the `--pred` option. 
It has one line per phenotype (in any order) that specifies the name of the phenotype and its
corresponding prediction file name. 
Each phenotype must have exactly one prediction file and phenotype names 
must match with those in the phenotype file.
Phenotypes in this file not included in the analysis are ignored.

Each prediction file contains the genetic predictions for the phenotype (space separated).

Line 1 starts with 'FID_IID' followed by $N$ sample identifiers.
It is followed by 23 lines containing the genetic predictions for each chromosome 
(sex chromosomes are collapsed into chromosome 23).

More specifically, each line has $N+1$ values which are the chromosome number followed by the $N$
leave-one chromosome out (LOCO) predictions for each individual.

Samples in this file not in the bed/pgen/bgen input file are ignored. Genotyped samples not 
present in this file will be ignored in the analysis of the corresponding trait. 

Samples with missing LOCO predictions must have their corresponding phenotype value set to missing.


## Options


| Option | Argument | Type | Description|
|---|-------|------|----|
|`--step`| INT| Required| specify step for the regenie run (see Overview) [argument can be `1` or `2`] |
|`--bt`| FLAG| Optional| specify that traits are binary with 0=control,1=case,NA=missing (default is quantitative)|
|`-1,--cc12`| FLAG| Optional| specify to use 1/2/NA encoding for binary traits (1=control,2=case,NA=missing)|
|`--bsize`| INT| Required| size of the genotype blocks|
|`--cv`| INT| Optional| number of cross validation (CV) folds [default is 5]|
|`--loocv`| FLAG | Optional| flag to use leave-one out cross validation|
|`--lowmem`| FLAG | Optional | flag to reduce memory usage by writing level 0 predictions to disk (details below). This is very useful if the number of traits is large (e.g. greater than 10)|
|`--lowmem-prefix`| FILE PREFIX | Optional | prefix where to temporarily write the level 0 predictions|
|`--nb`| INT| Optional| number of blocks (determined from block size if not provided)|
|`--strict`|FLAG| Optional| flag to removing samples with missing data at any of the phenotypes|
|`--ignore-pred`|FLAG| Optional| skip reading the file specified by `--pred` (corresponds to simple linear/logistic regression)|
|`--split`|FLAG| Optional| flag to split asssociation results into separate files for each trait. 
|`--gz`|FLAG| Optional| flag to output files in compressed gzip format (LOCO prediction files in step 1 and association results files in step 2) **[this only works when compiling with Boost Iostream library (see Install tab)]**. 
|`--force-impute`|FLAG| Optional| flag to keep and impute missing observations for QTs in step 2|
|`--write-samples`|FLAG| Optional| flag to write sample IDs for those kept in the analysis for each trait in step 2|
|`--firth`| FLAG | Optional | specify to use Firth likelihood ratio test (LRT) as fallback for p-values less than threshold|
|`--approx`|FLAG | Optional| flag to use approximate Firth LRT for computational speedup (only works when option `--firth` is used)|
|`--firth-se`| FLAG | Optional | flag to compute SE based on effect size and LRT p-value when using Firth correction (instead of based on Hessian of unpenalized log-likelihood)|
|`--spa`| FLAG | Optional| specify to use Saddlepoint approximation as fallback for p-values less than threshold|
|`--pThresh`| FLOAT | Optional| P-value threshold below which to apply Firth/SPA correction [default is 0.05]
|`--test`| STRING | Optional | specify to carry out dominant or recessive test [default is additive; argument can be `dominant` or `recessive`]|
|`--chr`| INT| Optional| specify which chromosomes to test in step 2 (use for each chromosome to include)|
|`--chrList` | STRING | Optional | Comma separated list of chromosomes to test in step 2|
|`--minMAC`| INT| Optional| flag to specify the minimum minor allele count (MAC) when testing variants [default is 5]. Variants with lower MAC are ignored.|
|`--nauto`| INT| Optional| number of autosomal chromosomes (for non-human studies) [default is 22]|
|`--niter`| INT| Optional| maximum number of iterations for logistic regression [default is 30]|
|`--maxstep-null`| INT| Optional| maximum step size for logistic model with Firth penalty under the null [default is 25]|
|`--maxiter-null`| INT| Optional| maximum number of iterations for logistic model with Firth penalty under the null [default is 1000]|
|`--threads`| INT | Optional| number of computational threads to use [default=all]|
|`--debug`| FLAG | Optional | debug flag (for use by developers)|
|`--verbose`| FLAG | Optional| verbose screen output|
|`--help`| FLAG | Optional| Prints usage and options list to screen|

When step 1 of **regenie** is run in low memory mode (i.e. using `--lowmem`), 
temporary files are created on disk (using `--lowmem-prefix tmp_prefix` determines 
where the files are written [as in `tmp_prefix_l0_Y1`,...,`tmp_prefix_l0_YP` 
for P phenotypes]). If the prefix is not specified, the default is to use the 
prefix specified by `--out` (see below).
These are automatically deleted at the end of the program (unless the run
was not successful in which case the user would need to delete the files)

## Output

| Option | Argument | Type | Description|
|---|-------|------|----|
|`--out`| FILE PREFIX| Required| Output files that depends on `--step`|

A log file `file.log` of the output is generated.

**Using `--step 1 --out file`**

For the \(P\) phenotypes, files `file_1.loco`,...,`file_P.loco` are output with the
per-chromosome LOCO predictions as rows of the files. 
If option `--gz` was used, the files will be compressed in gzip format and have extension `.loco.gz`.

Genotyped individuals specified using option `--remove` are excluded from this file. 
 Hence, this can be used if genotype files in step 1 and 2 have different number of samples 
 (so only keeping samples present in both files).

Individuals with missing phenotype values kept in the analysis 
are included in the file and have their predictions set to missing.

The list of blup files needed for step 2 (association testing) is written to  `file_pred.list`.

**Using`--step 2 --out file`** 

By default, results are written to a single file `file.regenie`, which has one line per
SNP along with a header line.
If option `--gz` was used, the file will be compressed in gzip format and have extension `.regenie.gz`.

The first 7 entries of each row specify chromosome, posistion, ID, reference allele (allele 0), 
alternative allele (allele 1), frequency of the alternative allele, and the test performed 
(additive/dominant/recessive).
With BGEN/PGEN files, the imputation INFO score is also provided (IMPUTE info score for BGEN and Mach Rsq for PGEN).
Allele frequency and INFO score, if applicable, are computed using all 
individuals included in the analysis (so they are the same for all phenotypes).

These are followed by the estimated effect sizes, standard errors, chi-square test statistics 
and \(-\log_{10}\) p-values for each phenotype.

When using option `--split`, the results are written in separate files for
each phenotype
`file_<phenotype1_name>.regenie,...,file_<phenotypeP_name>.regenie` 
with the same format.
If option `--gz` was used, the files will be compressed in gzip format and have extension `.regenie.gz`.
If option `--write-samples` was used, IDs of samples used for each trait will be written in files
`file_<phenotype1_name>.regenie.ids,...,file_<phenotypeP_name>.regenie.ids` (tab separated, no header).
