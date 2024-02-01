/* 

   This file is part of the regenie software package.

   Copyright (c) 2020-2024 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

*/

#ifndef TESTMCC_H
#define TESTMCC_H

//-------------------------------------------------------------
// Class MCCResults for MCC test results
//-------------------------------------------------------------
class MCCResults 
{
  public:
    double nl_dbl_dmin = 10.0 * std::numeric_limits<double>::min();

    unsigned int N; // number of samples                
    unsigned int K; // numbr of covariates
    unsigned int q; // number of traits                
    unsigned int M; // number of variables (to test for association)
    Eigen::VectorXd n, n2;

    // sum_x: q x M matrices
    Eigen::MatrixXd sum_x, sum_x_sq;
    Eigen::MatrixXd sum_x2, sum_x2_sq;
    Eigen::MatrixXd sum_x3;
    Eigen::MatrixXd sum_x4;

    Eigen::MatrixXd sum_nx, sum_nx_sq;
    Eigen::MatrixXd sum_nx2, sum_nx2_sq, sum_nx2_cub;
    Eigen::MatrixXd sum_nx3, sum_nx4, sum_nx6;

    // expectations: q x M matrices
    Eigen::MatrixXd A; // observed test statistic
    Eigen::MatrixXd EA, EA2, EA3, EA4;
    // moments
    Eigen::MatrixXd variance, kurtosis, skewness;
    // Beta/Gamma distr.
    Eigen::MatrixXd S1, S2;
    Eigen::MatrixXd S1g, S2g, L;
    MatrixXb SkipBeta, SkipGamma, Skip;
    
    Eigen::MatrixXd R, PvalBeta, PvalGamma, Pval;

    // DKAT
    Eigen::MatrixXd D;
    Eigen::MatrixXd Dm1, Dm2, Dm3;

    Eigen::MatrixXd ShapeD, ScaleD, LocationD;
    Eigen::MatrixXd PvalD;

    //----------
    // Methods
    //----------
    MCCResults();
    MCCResults(unsigned int, unsigned int, unsigned int, unsigned int, const Eigen::VectorXd &n_);
    ~MCCResults();

    void statistics(const MatrixXb &, const Eigen::MatrixXd &, const Eigen::MatrixXd &);
    void expectations(const MatrixXb &, const Eigen::MatrixXd &,
        const Eigen::VectorXd &, const Eigen::VectorXd &,
        const Eigen::VectorXd &, const Eigen::VectorXd &,
        const Eigen::VectorXd &, const Eigen::VectorXd &);
    void moments();
    void distr();

    void dkat(const MatrixXb &, const Eigen::MatrixXd &,
        const Eigen::MatrixXd &, 
        const Eigen::VectorXd &, const Eigen::VectorXd &,
        const Eigen::VectorXd &, const Eigen::VectorXd &,
        const Eigen::VectorXd &, const Eigen::VectorXd &,
        const Eigen::VectorXd &, const Eigen::VectorXd &);
};

//-------------------------------------------------------------
// Class MCC for Moment-Matching Correlation (MCC) test
//-------------------------------------------------------------
class MCC
{
  public:
    /********************
     * Public attributes
     ********************/
    int verbose;

    // dimensions
    unsigned int n_samples; // number of samples
    unsigned int n_traits; // number of traits
    unsigned int n_covariates; // number of covariates

    MatrixXb Mask; // masked samples = samples with missing values on phenotypes; 0 = missing 
    Eigen::VectorXd Neff; // number of non-missing samples (per trait)

    Eigen::MatrixXd Yres; // matrix of residualized traits (missing are set to 0)
    Eigen::MatrixXd Ynorm; // matrix of normalized residualized traits (missing are set to 0)

    // precomputed products
    Eigen::VectorXd sum_y, sum_y_sq;
    Eigen::VectorXd sum_y2, sum_y2_sq;
    Eigen::VectorXd sum_y3;
    Eigen::VectorXd sum_y4;

    // precomputed products for DKAT
    Eigen::VectorXd sum_ny, sum_ny_sq;
    Eigen::VectorXd sum_ny2, sum_ny2_sq, sum_ny2_cub;
    Eigen::VectorXd sum_ny3, sum_ny4, sum_ny6;

    /********************
     * Public methods
     ********************/
    void setup_y(const MatrixXb &, const Eigen::MatrixXd &res, unsigned int _n_covariates = 1);

    MCCResults run(const Eigen::MatrixXd& G);

    MCC();
    ~MCC() { };

  private:
    /********************
     * Private methods
     ********************/
    void check_setup_data();
    void precomp_sumy();
};

#endif
