/*

   This file is part of the regenie software package.

   Copyright (c) 2020-2021 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

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


#ifndef HLM_H
#define HLM_H

#include <LBFGS.h>

// to use heteroskedastic linear model
class HLM
{

  public:

    int first_fit = true; // only run null model once if no blups
    int max_iter = 100; // maximum number of iterations for LBFGS
    int max_iter_retry = 500; // maximum number of iterations for LBFGS
    int linesearch_try = 50; // number of linesearch trials
    int linesearch_retry = 200; // max number of linesearch trials
    int max_step_retry = 1000; // max step size
    int n; // sample size for each trait

    // to fit null
    Eigen::MatrixXd X; // covariates for mean
    Eigen::MatrixXd Vlin, V; // covariates for variance
    Eigen::VectorXd y; // phenotype analyzed
    Eigen::VectorXd alpha;
    Eigen::ArrayXd Vb, Dinv;
    ArrayXb mask;

    // stored est from null
    Eigen::MatrixXd Dinv_sqrt, yres;
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > Px; // for each trait

    void prep_run(struct phenodt const*,struct param const*);
    void get_alpha(Eigen::VectorXd const&);
    void get_beta_approx(Eigen::VectorXd&);
    void store_null_est(int const&);
    void residualize(int const&,Eigen::Ref<Eigen::MatrixXd>,Eigen::Ref<Eigen::MatrixXd>);


    // functor to minimize
    double operator()(Eigen::VectorXd const& beta, Eigen::VectorXd& gradient){

      // get Vb, Dinv and alpha
      get_alpha(beta);
      Eigen::ArrayXd esq = ((y - X * alpha).array()).square();

      // f = -2 ( ll/n + 0.5 log(2pi) )
      double fval = ((Vb + Dinv * esq) * mask.cast<double>()).sum() / n;
      // update gradient
      gradient = V.transpose() * ( ((1 - esq * Dinv) * mask.cast<double>()) / n ).matrix();

      return fval;
    }

    HLM();
    ~HLM();

};


void HLM_fitNull(HLM& nullHLM, struct ests const&,struct phenodt const&,struct in_files const&,struct param const&,mstream&);

#endif
