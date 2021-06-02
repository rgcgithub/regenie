**regenie** is a C++ program for whole genome regression modelling of large [genome-wide association studies](https://en.wikipedia.org/wiki/Genome-wide_association_study).

It is developed and supported by a team of scientists at the Regeneron Genetics Center.

The method has the following properties

- It works on quantitative and binary traits, including binary traits with unbalanced case-control ratios
- It can process multiple phenotypes at once
- It is fast and memory efficient 🔥
- For binary traits it supports Firth logistic regression and an SPA test
- It can perform gene/region-based burden tests
- It supports the [BGEN](https://www.well.ox.ac.uk/~gav/bgen_format/), [PLINK](https://www.cog-genomics.org/plink/1.9/formats#bed) bed/bim/fam and [PLINK2](https://www.cog-genomics.org/plink/2.0/formats#pgen) pgen/pvar/psam genetic data formats
- It is ideally suited for implementation in [Apache Spark](https://spark.apache.org/) (see [GLOW](https://projectglow.io/))
- It can be installed with [Conda](https://anaconda.org/bioconda/regenie) [![Regenie](https://anaconda.org/bioconda/regenie/badges/installer/conda.svg)](https://anaconda.org/bioconda/regenie)

Full documentation for the **regenie** can be found [here](https://rgcgithub.github.io/regenie/).

## Citation 
Mbatchou, J., Barnard, L., Backman, J. et al. Computationally efficient whole-genome regression for quantitative and binary traits. Nat Genet (2021).[https://doi.org/10.1038/s41588-021-00870-7](https://doi.org/10.1038/s41588-021-00870-7)

## License

**regenie** is distributed under an [MIT license](https://github.com/rgcgithub/regenie/blob/master/LICENSE).

## Contact
If you have any questions about regenie please contact

- <jonathan.marchini@regeneron.com>
- <joelle.mbatchou@regeneron.com>

If you want to submit a issue concerning the software please do so
using the **regenie** [Github repository](https://github.com/rgcgithub/regenie/issues).


## Version history
Version 2.1 (Added interaction testing functionality for GxG and GxE tests [options `--interaction/--interaction-snp, respectively]`; new options `--write-null-firth/--use-null-firth` to store the estimates from approximate Firth null model; new option `--minCaseCount` to filter out BTs with low number of cases from the analysis; new option `--no-split` to enforce output of summary stats to a single file for all traits; added support for tranposed phenotype file format with `--tphenoFile`)

Version 2.0.2 (Bug fix for burden testing with BGEN files not in v1.2 with 8-bit encoding; enabled faster step 2 implementation with Zstd compressed BGEN files in v1.2 with 8-bit encoding)

Version 2.0.1 (New option `--catCovList` to specify categorical covariates; Enabled parameter expansion when specifying select phenotypes/covariates to analyze [e.g. 'PC{1:10}'])

Version 2.0 (Added burden testing functionality for region or gene-based tests [see [website](https://rgcgithub.github.io/regenie/options/#burden-testing) for details]; added sample size column in summary stats output).

For past releases, see [here](RELEASE_LOG.md).

