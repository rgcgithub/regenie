[![build](https://github.com/rgcgithub/regenie/actions/workflows/test.yml/badge.svg)](https://github.com/rgcgithub/regenie/actions/workflows/test.yml)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/rgcgithub/regenie?logo=Github)
[![install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg)](https://anaconda.org/bioconda/regenie)
[![Github All Releases](https://img.shields.io/github/downloads/rgcgithub/regenie/total.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**regenie** is a C++ program for whole genome regression modelling of large [genome-wide association studies](https://en.wikipedia.org/wiki/Genome-wide_association_study).

It is developed and supported by a team of scientists at the Regeneron Genetics Center.

The method has the following properties

- It works on quantitative and binary traits, including binary traits with unbalanced case-control ratios
- It can handle population structure and relatedness
- It can process multiple phenotypes at once efficiently
- It is fast and memory efficient 🔥
- For binary traits, it supports Firth logistic regression and an SPA test
- It can perform gene/region-based tests, interaction tests and conditional analyses
- It supports the [BGEN](https://www.well.ox.ac.uk/~gav/bgen_format/), [PLINK](https://www.cog-genomics.org/plink/1.9/formats#bed) bed/bim/fam and [PLINK2](https://www.cog-genomics.org/plink/2.0/formats#pgen) pgen/pvar/psam genetic data formats
- It is ideally suited for implementation in [Apache Spark](https://spark.apache.org/) (see [GLOW](https://projectglow.io/))
- It can be installed with [Conda](https://anaconda.org/bioconda/regenie)

Full documentation for the **regenie** can be found [here](https://rgcgithub.github.io/regenie/).

## Citation 
Mbatchou, J., Barnard, L., Backman, J. et al. Computationally efficient whole-genome regression for quantitative and binary traits. Nat Genet 53, 1097–1103 (2021). https://doi.org/10.1038/s41588-021-00870-7

## License

**regenie** is distributed under an [MIT license](https://github.com/rgcgithub/regenie/blob/master/LICENSE).

## Contact
If you have any questions about regenie please contact

- <jonathan.marchini@regeneron.com>
- <joelle.mbatchou@regeneron.com>

If you want to submit a issue concerning the software please do so
using the **regenie** [Github repository](https://github.com/rgcgithub/regenie/issues).


## Version history
[Version 3.4.1](https://github.com/rgcgithub/regenie/releases/tag/v3.4.1) (Reduction in memory usage for LD computation when writing to text files; Fix bug rejecting valid PVAR files)

[Version 3.4](https://github.com/rgcgithub/regenie/releases/tag/v3.4) (Reduction in memory usage for LD computation with dosages; Minor bug fixes for LD computation; Bug fix for when carriage returns are in optional input files)

[Version 3.3](https://github.com/rgcgithub/regenie/releases/tag/v3.3) (Faster implementation of approximate Firth LRT; New strategy for approximate Firth LRT with ultra-rare variants; Relaxed convergence criterion of Firth LRT from 1E-4 to 2.5E-4)

[Version 3.2.9](https://github.com/rgcgithub/regenie/releases/tag/v3.2.9) (Switch to robust version of ACAT to handle very small p-values; Bug fix for Step1 when sex chromosome was included in the analysis; Allow for 64 domains when using the 4-column annotation file)

[Version 3.2.8](https://github.com/rgcgithub/regenie/releases/tag/v3.2.8) (New option `--bgi` to specify custom index bgi file accompagnying BGEN file; Relax matching criteria between BGEN and index bgi files to use CPRA instead of variant ID)

[Version 3.2.7](https://github.com/rgcgithub/regenie/releases/tag/v3.2.7) (New option `--force-mac-filter` to apply different MAC filter to subset of SNPs; Extend maximum number of domains to 32 for 4-column anno-file; Update PGEN library)

[Version 3.2.6](https://github.com/rgcgithub/regenie/releases/tag/v3.2.6) (Relax tolerance parameter for null unpenalized logistic regression from 1e-8 to 1e-6; Minor bug fixes)

[Version 3.2.5.3](https://github.com/rgcgithub/regenie/releases/tag/v3.2.5.3) (Fix inflation issue when testing main effect of SNP in GxE model; Minor bug fixes)

[Version 3.2.5](https://github.com/rgcgithub/regenie/releases/tag/v3.2.5) (Use pseudo-data representation algorithm as default in step 2 single variant tests; Use ACAT to get SBAT p-value across POS/NEG models; Bug fix for ACATV when set has a single variant with zero weight)

[Version 3.2.4](https://github.com/rgcgithub/regenie/releases/tag/v3.2.4) (Relaxed the requirement on the minimum number of unique values for QTs to 3; Various bug fixes)

[Version 3.2.3](https://github.com/rgcgithub/regenie/releases/tag/v3.2.3) (Address convergence issues in Firth regression; Various bug fixes)

[Version 3.2.2](https://github.com/rgcgithub/regenie/releases/tag/v3.2.2) (New columns in sumstats file (N_CASES/N_CONTROLS) to output the number of cases/controls when using `--af-cc`; Various bug fixes)

[Version 3.2.1](https://github.com/rgcgithub/regenie/releases/tag/v3.2.1) (New option `--lovo-snplist` to only consider a subset of LOVO masks; Improve efficiency of LOVO for large sets to reduce memory usage; Bug fix for SPA with numerical overflow; For SKAT/ACAT tests with Firth correction, don't include SKAT weights when running Firth on single variants)

[Version 3.2](https://github.com/rgcgithub/regenie/releases/tag/v3.2) (Bug fix for SKAT/SKATO when testing on binary traits using Firth/SPA; Switched name of NNLS joint test to SBAT test altering name of corresponding options and applied Bonferroni correction before reporting its p-value [correcting for minP of 2 tests])

[Version 3.1.4](https://github.com/rgcgithub/regenie/releases/tag/v3.1.4) (New option `--par-region` to specify build to determine bounds for chrX PAR regions; new option `--force-qt` to force QT runs for traits with fewer than 10 values [otherwise will throw an error]; phenotype imputation for missing values is now applied after RINTing when using `--apply-rint`; several bug fixes)

[Version 3.1.2](https://github.com/rgcgithub/regenie/releases/tag/v3.1.2) (Reduction in memory usage for SKAT/SKATO tests; Bug fix for LOVO with SKAT/ACAT tests; Improvements for null Firth logistic algorithm to address reported convergence issues)

[Version 3.1.1](https://github.com/rgcgithub/regenie/releases/tag/v3.1.1) (Reduction in memory usage for SKAT/SKATO tests; Improvements for logistic regressions algorithms to address reported convergence issues)

[Version 3.1](https://github.com/rgcgithub/regenie/releases/tag/v3.1) (Fixed bug in SKAT/SKATO tests when applying Firth/SPA correction; Improved SPA implementation by computing both tail probabilities; New option `--set-singletons` to specify variants to consider as singletons for burden masks; New option `--l1-phenoList` to run level 1 models in Step 1 in parallel across phenotypes; Several bug fixes)

[Version 3.0.3](https://github.com/rgcgithub/regenie/releases/tag/v3.0.3) (Skip BTs where null model fit failed; Bug fix for BURDEN-ACAT; Bug fix when nan/inf values are in phenotype/covariate file)

[Version 3.0.1](https://github.com/rgcgithub/regenie/releases/tag/v3.0.1) (Improve ridge logistic regression in Step 1; Add compilation with Cmake)

[Version 3.0](https://github.com/rgcgithub/regenie/releases/tag/v3.0) (New gene-based tests: SKAT, SKATO, ACATV, ACATO and NNLS [Non-Negative Least Square test]; New GxE and GxG interaction testing functionality; New conditional analysis functionality; see [release page](https://github.com/rgcgithub/regenie/releases/tag/v3.0) for minor additions)

For past releases, see [here](RELEASE_LOG.md).
