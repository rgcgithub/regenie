##Download

The **regenie** source code is hosted on
[Github](https://github.com/rgcgithub/regenie).

##Installation

1. **regenie** requires the
  [BGEN library](https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk) so
  you will need to download and install that library.
2. In the source code edit the BGEN_PATH variable in the `Makefile`
   to the BGEN library path.
3. On the command line type `make` while in the main source code directory.
    - Compile using [GCC](https://gcc.gnu.org) version >= 5.1
4. This should produce the executable called `regenie`.


##Computing requirements

We have tested **regenie** on Linux 
<!-- Issue when using BGEN library API on MAC (gives different results than on Linux)
and Mac OSX 
-->
computing environment.

### Memory usage
In both Step 1 and Step 2 of a **regenie** run the genetic data file is
read once, in blocks of SNPs, so at no point is the full dataset ever stored in
memory.

**regenie** uses a dimension reduction approach using ridge regression
  to produce a relatively small set of genetic predictors, that are
  then used to fit a whole-genome regression model. These genetic
  predictors are stored in memory by default, and can be relatively
  large if many phenotypes are stored at once.

For example, if there are \(P\) phenotypes, \(M\) SNPs and \(N\) samples, and a
block size of \(B\) SNPs is used with \(R\) ridge parameters,
 then **regenie** needs to store roughly \(N\times M/B\times R\)
doubles per phenotype, which is 8Gb per phenotype when \(M=500,000,
N=400,000, B =1,000,R=5\) and 200Gb in total when \(P=25\).

However, the `--lowmem` option can be used to avoid that memory usage,
at negligible extra computational cost, by writing temporary files to disk.

### Threading

**regenie** can take advantage of multiple cores using threading. The
number of threads can be specified using the `--threads` option.

**regenie** uses the [Eigen library](http://eigen.tuxfamily.org/index.php?title=Main_Page) for 
efficient linear algebra operations and this uses threading where possible.

For PLINK bed/bim/fam files, PLINK2 pgen/pvar/psam files, as well as BGEN v1.2 files with 8-bit encoding (format used for UK Biobank
500K imputed data), step 2 of **regenie** has been optimized by 
using multithreading through [OpenMP](https://www.openmp.org).



