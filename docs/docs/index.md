## regenie 

**regenie** is a C++ program for whole genome regression modelling of large
[genome-wide association studies](https://en.wikipedia.org/wiki/Genome-wide_association_study).

It is developed and supported by a team of scientists at the Regeneron Genetics Center.

The method has the following properties

- It works on quantitative and binary traits, including binary
traits with unbalanced case-control ratios
- It can handle population structure and relatedness
- It can process multiple phenotypes at once efficiently
- For binary traits, it supports Firth logistic regression and an SPA test
- It can perform gene/region-based tests (Burden, [SBAT](https://doi.org/10.1016/j.ajhg.2024.08.021), SKAT/SKATO, ACATV/ACATO)
- It can perform interaction tests (GxE, GxG) as well as conditional analyses
- Meta-analysis of REGENIE summary statistics can be performed using [REMETA](https://rgcgithub.github.io/remeta/)
- It is fast and memory efficient 🔥
- It supports the [BGEN](https://www.well.ox.ac.uk/~gav/bgen_format/), [PLINK](https://www.cog-genomics.org/plink/1.9/formats#bed) bed/bim/fam and [PLINK2](https://www.cog-genomics.org/plink/2.0/formats#pgen) pgen/pvar/psam genetic data formats
- It is ideally suited for implementation in
  [Apache Spark](https://spark.apache.org/) (see [GLOW](https://projectglow.io/))
- It can be installed with [Conda](https://anaconda.org/bioconda/regenie)

##Citation

Mbatchou, J., Barnard, L., Backman, J. et al. Computationally efficient whole-genome regression for quantitative and binary traits. Nat Genet 53, 1097–1103 (2021). [https://doi.org/10.1038/s41588-021-00870-7](https://doi.org/10.1038/s41588-021-00870-7)


## License 

**regenie** is distributed under an [MIT license](https://github.com/rgcgithub/regenie/blob/master/LICENSE).


## Contact 

If you have any questions about **regenie** please contact

- <jonathan.marchini@regeneron.com>
- <joelle.mbatchou@regeneron.com> 

If you want to submit a issue concerning the software please do so
using the **regenie** [Github repository](https://github.com/rgcgithub/regenie/issues).

<a href="https://github.com/rgcgithub/regenie/issues/new"><button class="btn btn-primary btn-sm" type="submit"><i class="fab fa-github fa-2x"></i> Report Issue</button></a>
<a href="https://github.com/rgcgithub/regenie/issues"><button class="btn btn-primary btn-sm" type="submit"><i class="fab fa-github fa-2x"></i> Active Issues </button></a>

<!--
## Version history

Version 1.0 (22 June 2020): Initial release
-->
