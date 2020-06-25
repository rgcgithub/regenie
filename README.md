**regenie** is a C++ program for whole genome regression modelling of large [genome-wide association studies](https://en.wikipedia.org/wiki/Genome-wide_association_study).

It is developed and supported by a team of scientists at the Regeneron Genetics Center.

The method has the following properties

- It works on quantitative and binary traits, including binary traits with unbalanced case-control ratios
- It can process multiple phenotypes at once
- It is fast and memory efficient 🔥
- For binary traits it supports Firth logistic regression and an SPA test
- It supports the [BGEN](https://www.well.ox.ac.uk/~gav/bgen_format/) and [PLINK](https://www.cog-genomics.org/plink/1.9/formats#bed) bed/bim/fam genetic data formats
- It is ideally suited for implementation in [Apache Spark](https://spark.apache.org/) (see [GLOW](https://projectglow.io/))

Full documentation for the **regenie** can be found [here](https://rgcgithub.github.io/regenie/).

## Frequently asked questions (F.A.Q.)

1. What to do if Step 1 of REGENIE failed for a binary trait when fitting the penalized logsitic regression model? This can occur when the sample size used to fit the model is small.
    - If using K-fold CV, switch to LOOCV (option `--loocv`) to increase the size of the sample used to fit the model
        - LOOCV is now used by default when the sample size is less than 5,000
    - If it is due to quasi-separation (i.e. `Var(Y)=0` occurred in model fitting), either increase the sample size using LOOCV or increase the MAF threshold for variants included in step 1 analysis 

2. What to do if Step 2 of REGENIE fails when fitting the null model for the approximate Firth correction? 
    - We have implemented the same measures as in `logistf` function in R to avoid convergence issues, which include the use of a step size threshold value when performing a Newton step. 
We first try fit the model with a step size threshold that is more liberal (=25) as well as a maximum number of iterations of 1,000 and if convergence fails, we use a more stringent maximum step size value (=5) which slows down convergence as well as increase the number of iterations to 5,000.
    - The user can also specify a maximum step size threshold using `--maxstep-null` (use value <5) as well as increase the maximum number of iterations using `--maxiter-null` (use value >5000). In that case, no retries are perfomed if convergence fails.
        - We recommend to test chromosomes separately (using `--chr`) as these parameters may need to be altered when fitting the null model for each chromosome

3. Can REGENIE be used on data with small sample sizes? 
    - For quantitative traits, we have not obtained issues running REGENIE on small data sets.
    - For binary traits, we have obtained successful runs of REGENIE (step 1 and 2) on data with as little as 300 samples. A few factors to consider:
        - Convergence issues may occur in step 1 (all the more if a trait is highly unbalanced) -> see (1) above
        - Similarly, convergence issues may occur in step 2 when using Firth approximation -> see (2) above 
 

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
Version 1.0.1 (fixed numerical overflow bug for quantile calculation and added new strategy for fitting null model for approximate Firth test [see (2) in FAQ]) 
Version 1.0 (22 June 2020): Initial release



