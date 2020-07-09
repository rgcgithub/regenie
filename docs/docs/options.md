## Documentation
## Getting started

To run **regenie**, use the command `./regenie` on the command line,
followed by options and flags as needed.

To get a full list of options use

```
./regenie --help
```

The directory `/examples` contains some small example files that are
useful when getting started. A test run on a set of binary traits can be achieved by the
following 2 commands.

In **Step 1** the whole genome regression model is fit to the traits, and
a set of genomic predictions are produced as output

```
./regenie \
  --step 1 \
  --bed example/example \
  --exclude example/snplist_rm.txt \
  --c example/covariates.txt \
  --p example/phenotype_bin.txt \
  --remove example/fid_iid_to_remove.txt \
  --b 100 \
  --bt --lowmem \
  --o fit_bin_out
```

In **Step 2** a set of imputed SNPs are tested for association using a
Firth logistic regression model

```
./regenie \
  --step 2 \
  --bgen example/example.bgen \
  --c example/covariates.txt \
  --p example/phenotype_bin.txt \
  --remove example/fid_iid_to_remove.txt \
  --b 200 \
  --bt \
  --firth 0.01 --approx \
  --pred fit_bin_out_pred.list \
  --split \
  --o test_bin_out_firth
```

One of the output files from this command with association results is included in `example/example.test_bin_out_firth_Y1.regenie`.

## Input 


| Option | Argument | Type | Description|
|---|-------|------|----|
|`--bgen, --bed, --pgen`  | FILE | Required |Input genetic data file. Either BGEN file eg. `file.bgen`, or bed/bim/fam prefix that assumes`file.bed`, `file.bim`, `file.fam` exist, or pgen/pvar/psam prefix that assumes`file.pgen`, `file.pvar`, `file.psam` exist |
|`--sample`  | FILE | Optional |Sample file corresponding to input BGEN file|
|`--keep`  | FILE | Optional | Inclusion file that lists individuals to retain in the analysis|
|`--remove`  | FILE | Optional | Exclusion file that lists individuals to remove from the analysis|
|`--extract`  | FILE | Optional | Inclusion file that lists IDs of variants to keep **(only works with option `--step 1`)**|
|`--exclude`  | FILE | Optional | Exclusion file that lists IDs of variants to remove **(only works with option `--step 1`)**|
|`--p`  | FILE | Required |Phenotypes file|
|`--phenoCol` | STRING | Optional | Use for each phenotype you want to include in the analysis|
|`--c`  | FILE | Optional | Covariates file|
|`--covarCol` | STRING | Optional | Use for each covariate you want to include in the analysis|
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
option `--export bgen-1.2 'bits=8'`).

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

Running `--step 1 --o foo` will produce

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

  Samples must be in the bed/pgen/bgen input file and must be included in
  the analysis (otherwise use `--remove`).

  For each phenotype, samples with missing LOCO predictions must have their corresponding 
phenotype value set to missing, and **all samples with non-missing phenotype values must have LOCO
predictions** (otherwise use `--remove`).


## Options


| Option | Argument | Type | Description|
|---|-------|------|----|
|`--step`| INT| Required| specify step for the regenie run (see Overview) [argument can be `1` or `2`] |
|`--bt`| FLAG| Optional| specify that traits are binary with 0=control,1=case,NA=missing (default is quantitative)|
|`--1`| FLAG| Optional| specify to use 1/2/NA encoding for binary traits (1=control,2=case,NA=missing)|
|`--b`| INT| Required| size of the genotype blocks|
|`--cv`| INT| Optional| number of cross validation (CV) folds [default is 5]|
|`--loocv`| FLAG | Optional| flag to use leave-one out cross validation|
|`--lowmem`| FILE PREFIX | Optional | flag to reduce memory usage by writing level 0 predictions to disk (details below). This is very useful if the number of traits is large (e.g. greater than 10)|
|`--nb`| INT| Optional| number of blocks (determined from block size if not provided)|
|`--strict`|FLAG| Optional| flag to removing samples with missing data at any of the phenotypes|
|`--ignore-pred`|FLAG| Optional| skip reading the file specified by `--pred` (corresponds to simple linear/logistic regression)|
|`--split`|FLAG| Optional| flag to split asssociation results into separate files for each trait. 
|`--force-impute`|FLAG| Optional| flag to keep and impute missing observations for QTs in step 2|
|`--firth`| FLOAT | Optional | specify to use Firth likelihood ratio test as fallback for p-values less than the specified threshold [default is 0.05]|
|`--approx`|FLAG | Optional| flag to use approximate Firth LRT for computational speedup (only works when option `--firth` is used)|
|`--spa`| FLOAT | Optional| specify to use Saddlepoint approximation as fallback for p-values less than the specified threshold [default is 0.05]|
|`--test`| STRING | Optional | specify to carry out dominant or recessive test [default is additive; argument can be `dominant` or `recessive`]|
|`--chr`| INT| Optional| specify which chromosomes to test in step 2 (use for each chromosome to include)|
|`--minMAC`| INT| Optional| flag to specify the minimum minor allele count (MAC) when testing variants [default is 5]. Variants with lower MAC are ignored.|
|`--nauto`| INT| Optional| number of autosomal chromosomes (for non-human studies) [default is 22]|
|`--niter`| INT| Optional| maximum number of iterations for logistic regression [default is 30]|
|`--maxstep-null`| INT| Optional| maximum step size for logistic model with Firth penalty under the null [default is 25]|
|`--maxiter-null`| INT| Optional| maximum number of iterations for logistic model with Firth penalty under the null [default is 1000]|
|`--threads`| INT | Optional| number of computational threads to use [default=all]|
|`--debug`| FLAG | Optional | debug flag (for use by developers)|
|`--v`| FLAG | Optional| verbose screen output|
|`--help`| FLAG | Optional| Prints usage and options list to screen|

When step 1 of **regenie** is run in low memory mode (i.e. using `--lowmem prefix`), 
temporary files are created on disk where the prefix argument, if specified, 
determines where the files are written (as in `prefix_l0_Y1`,...,`prefix_l0_YP` 
for P phenotypes). If the prefix argument is omitted, the default is to use the 
prefix specified by `--o` (see below).
These are automatically deleted at the end of the program (unless the run
was not successful in which case the user would need to delete the files)

## Output

| Option | Argument | Type | Description|
|---|-------|------|----|
|`--o`| FILE PREFIX| Required| Output files that depends on `--step`|

A log file `file.log` of the output is generated.

**Using `--step 1 --o file`**

For the \(P\) phenotypes, files `file_1.loco`,...,`file_P.loco` are output with the
per-chromosome LOCO predictions as columns of the files.

Genotyped individuals specified using option `--remove` are excluded from this file. 
 Hence, this can be used if genotype files in step 1 and 2 have different number of samples 
 (so only keeping samples present in both files).

Individuals with missing phenotype values kept in the analysis 
are included in the file and have their predictions set to missing.

The list of blup files needed for step 2 (association testing) is written to  `file_pred.list`.

**Using`--step 2 --o file`** 

By default, results are written to a single file `file.regenie`, which has one line per
SNP along with a header line.

The first 7 entries of each row specify chromosome, posistion, ID, reference allele (allele 0), 
alternative allele (allele 1), frequency of the alternative allele, and the test performed 
(additive/dominant/recessive).
With BGEN files, the imputation INFO score is also provided. 
Allele frequency and INFO score, if applicable, are computed using all 
individuals included in the analysis (so they are the same for all phenotypes).

These are followed by the estimated effect sizes, standard errors, chi-square test statistics 
and \(-\log_{10}\) p-values for each phenotype.

When using option `--split`, the results are written in separate files for
each phenotype
`file_<phenotype1_name>.regenie,...,file_<phenotypeP_name>.regenie` 
with the same format.

