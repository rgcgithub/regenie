**regenie** is a C++ program for whole genome regression modelling of large [genome-wide association studies](https://en.wikipedia.org/wiki/Genome-wide_association_study).

It is developed and supported by a team of scientists at the Regeneron Genetics Center.

The method has the following properties

- It works on quantitative and binary traits, including binary traits with unbalanced case-control ratios
- It can process multiple phenotypes at once
- It is fast and memory efficient 🔥
- For binary traits it supports Firth logistic regression and an SPA test
- It supports the [BGEN](https://www.well.ox.ac.uk/~gav/bgen_format/), [PLINK](https://www.cog-genomics.org/plink/1.9/formats#bed) bed/bim/fam and [PLINK2](https://www.cog-genomics.org/plink/2.0/formats#pgen) pgen/pvar/psam genetic data formats
- It is ideally suited for implementation in [Apache Spark](https://spark.apache.org/) (see [GLOW](https://projectglow.io/))
- It can be installed with [Conda](https://anaconda.org/bioconda/regenie) [![Regenie](https://anaconda.org/bioconda/regenie/badges/installer/conda.svg)](https://anaconda.org/bioconda/regenie)

Full documentation for the **regenie** can be found [here](https://rgcgithub.github.io/regenie/).

## Citation 
Joelle Mbatchou, Leland Barnard, Joshua Backman, Anthony Marcketta, Jack A. Kosmicki, Andrey Ziyatdinov, Christian Benner, Colm O'Dushlaine, Mathew Barber, Boris Boutkov, Lukas Habegger, Manuel Ferreira, Aris Baras, Jeffrey Reid, Goncalo Abecasis, Evan Maxwell, Jonathan Marchini. (2020) Computationally efficient whole genome regression for quantitative and binary traits [[BioRxiv pre-print]](https://www.biorxiv.org/content/10.1101/2020.06.19.162354v1)

## License

**regenie** is distributed under an [MIT license](https://github.com/rgcgithub/regenie/blob/master/LICENSE).

## Contact
If you have any questions about regenie please contact

- <jonathan.marchini@regeneron.com>
- <joelle.mbatchou@regeneron.com>

If you want to submit a issue concerning the software please do so
using the **regenie** [Github repository](https://github.com/rgcgithub/regenie/issues).


## Version history
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

For past releases, see [here](RELEASE_LOG.md).

