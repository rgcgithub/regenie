## Frequently asked questions
<br/>

*    <span style="font-size: large; font-style: italic;color:#404040"> What block size to use in step 1? 
</span>

We recommend to use blocks of size 1000 as we have observed that it leads to a reasonable number of ridge predictors 
at level 1 (e.g. 2,500 with 500K SNPs used and the default **regenie** parameters) and have noticed little change in the 
final predictions when varying the block size.

<br/>

*    <span style="font-size: large; font-style: italic;color:#404040"> How many variants to use in step 1? 
</span>

We recommend to use a smaller set of about 500K directly genotyped SNPs in step 1, which should be sufficient to capture genome-wide polygenic effects. Note that using too many SNPs in Step 1 (e.g. >1M) can lead to a high computational burden due to the resulting higher number of predictors in the level 1 models.

<br/>

*    <span style="font-size: large; font-style: italic;color:#404040"> What do I do if I get the error "Uh-oh, SNP XX has low variance (=XX)" in step 1? 
</span>

This is due to variants with very low minor allele count (MAC) being included in step 1. To avoid this, you should use a MAC filter to remove such variants in a pre-processing step before running Regenie.

For example, in PLINK2 you would use the `--mac` option and obtain a list of variants that pass the MAC filter (note that if you are using `--keep/--remove` in Regenie, you should also use it in the PLINK2 command)
```
plink2 \
  --bfile my_bed_file \
  --mac 100 \
  --write-snplist \
  --out snps_pass
```

You would then use the output file in **regenie** as `--extract snps_pass.snplist` (and this would avoid having to make a new genotype file).


 <br/>

*    <span style="font-size: large; font-style: italic;color:#404040"> Can **regenie** be run on small sample sizes? 
</span>

    - For quantitative traits, we have not obtained issues running **regenie** on small data sets.
    - For binary traits, we have obtained successful runs of **regenie** (step 1 and 2) on data sets with as little as 300 samples. A few factors to consider:
        - Convergence issues may occur in step 1 (all the more if a trait is highly unbalanced) \(-\) see below
        - Similarly, convergence issues may occur in step 2 when using Firth approximation \(-\) see below 

 <br/>

*    <span style="font-size: large;font-style: italic; color:#404040"> What to do if Step 1 of **regenie** failed for a binary trait when fitting the penalized logsitic regression model? 
</span>

    - This can occur when the sample size used to fit the model is small and/or if the trait is extremely unbalanced. 
    - If using K-fold CV, switch to LOOCV (option `--loocv`) to increase the size of the sample used to fit the model
        - LOOCV is now used by default when the sample size is less than 5,000
    - If it is due to quasi-separation (i.e. `Var(Y)=0` occurred in model fitting), either increase the sample size using LOOCV or increase the MAF threshold for variants included in step 1 analysis 

<br/>

*    <span style="font-size: large;font-style: italic; color:#404040"> What to do if Step 2 of **regenie** fails when fitting the null model for the approximate Firth correction? 
</span>

    - This can occur when the sample size used to fit the model is small and/or if the trait is extremely unbalanced. 
    - We have implemented the same measures as in the `logistf` function in R to avoid convergence issues, which include the use of a step size threshold when performing a Newton step. 
        - We first try fitting the model with a step size threshold that is more liberal (=25) as well as a maximum number of iterations of 1,000 and if convergence fails, we retry the model fit using a more stringent step size threshold (=5) and a higher threshold for the number of iterations (=5,000), which will slow down convergence.
    - The user can also specify a maximum step size threshold using `--maxstep-null` (use value <5) as well as increase the maximum number of iterations using `--maxiter-null` (use value >5000). In that case, no retries are perfomed if convergence fails.
        - We recommend to test chromosomes separately (using `--chr`) as these parameters may need to be altered when fitting the null model for each chromosome


