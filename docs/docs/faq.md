## Frequently asked questions (F.A.Q.)

- Can REGENIE be run on small sample sizes? 
    - For quantitative traits, we have not obtained issues running REGENIE on small data sets.
    - For binary traits, we have obtained successful runs of REGENIE (step 1 and 2) on data sets with as little as 300 samples. A few factors to consider:
        - Convergence issues may occur in step 1 (all the more if a trait is highly unbalanced) -> see below
        - Similarly, convergence issues may occur in step 2 when using Firth approximation -> see below 
 

2. What to do if Step 1 of REGENIE failed for a binary trait when fitting the penalized logsitic regression model? This can occur when the sample size used to fit the model is small.
    - If using K-fold CV, switch to LOOCV (option `--loocv`) to increase the size of the sample used to fit the model
        - LOOCV is now used by default when the sample size is less than 5,000
    - If it is due to quasi-separation (i.e. `Var(Y)=0` occurred in model fitting), either increase the sample size using LOOCV or increase the MAF threshold for variants included in step 1 analysis 

3. What to do if Step 2 of REGENIE fails when fitting the null model for the approximate Firth correction? 
    - We have implemented the same measures as in `logistf` function in R to avoid convergence issues, which include the use of a step size threshold value when performing a Newton step. 
We first try fit the model with a step size threshold that is more liberal (=25) as well as a maximum number of iterations of 1,000 and if convergence fails, we use a more stringent maximum step size value (=5) as well as increase the number of iterations to 5,000, which slow down convergence.
    - The user can also specify a maximum step size threshold using `--maxstep-null` (use value <5) as well as increase the maximum number of iterations using `--maxiter-null` (use value >5000). In that case, no retries are perfomed if convergence fails.
        - We recommend to test chromosomes separately (using `--chr`) as these parameters may need to be altered when fitting the null model for each chromosome


