## Frequently asked questions
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


