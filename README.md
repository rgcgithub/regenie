**regenie** is a C++ program for whole genome regression modelling of large [genome-wide association studies](https://en.wikipedia.org/wiki/Genome-wide_association_study).

It is developed and supported by a team of scientists at the Regeneron Genetics Center.

The method has the following properties

- It works on quantitative and binary traits, including binary traits with unbalanced case-control ratios
- It can process multiple phenotypes at once
- It is fast and memory efficient 🔥
- For binary traits it supports Firth logistic regression and an SPA test
- It supports the [BGEN](https://www.well.ox.ac.uk/~gav/bgen_format/), [PLINK](https://www.cog-genomics.org/plink/1.9/formats#bed) bed/bim/fam and [PLINK2](https://www.cog-genomics.org/plink/2.0/formats#pgen) pgen/pvar/psam genetic data formats
- It is ideally suited for implementation in [Apache Spark](https://spark.apache.org/) (see [GLOW](https://projectglow.io/))

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



