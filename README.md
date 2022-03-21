[![build](https://github.com/rgcgithub/regenie/actions/workflows/test.yml/badge.svg)](https://github.com/rgcgithub/regenie/actions/workflows/test.yml)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/rgcgithub/regenie?logo=Github)
[![Regenie](https://anaconda.org/bioconda/regenie/badges/installer/conda.svg)](https://anaconda.org/bioconda/regenie)
[![Github All Releases](https://img.shields.io/github/downloads/rgcgithub/regenie/total.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**regenie** is a C++ program for whole genome regression modelling of large [genome-wide association studies](https://en.wikipedia.org/wiki/Genome-wide_association_study).

It is developed and supported by a team of scientists at the Regeneron Genetics Center.

The method has the following properties

- It works on quantitative and binary traits, including binary traits with unbalanced case-control ratios
- It can handle population structure and relatedness
- It can process multiple phenotypes at once efficiently
- It is fast and memory efficient 🔥
- For binary traits it supports Firth logistic regression and an SPA test
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
Version 3.0.1 (Improve ridge logistic regression in Step 1; Add compilation with Cmake)

Version 3.0 (New gene-based tests: SKAT, SKATO, ACATV, ACATO and NNLS [Non-Negative Least Square test]; New GxE and GxG interaction testing functionality; New conditional analysis functionality; see [release page](https://github.com/rgcgithub/regenie/releases/tag/v3.0) for minor additions)

For past releases, see [here](RELEASE_LOG.md).

