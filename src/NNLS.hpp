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

#ifndef NNLS_H
#define NNLS_H

#include <vector>
#include <unordered_set>
#include <numeric>
#include <iostream>
#include <fstream>
#include <random>
#include <list>

#include <boost/math/distributions.hpp>
#include <boost/math/special_functions/binomial.hpp> // binomial_coefficient

#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

/****************************************************
 * Declaration of main functions for Joint Burden test 
 * with the prefix "jburden_" 
****************************************************/

// main function for the NNLS test
// - input: residualized y and X (covariates and mean are projected out)
// - output: bhat & p-values for three models
//   -- OLS: y = Xb + e
//   -- NNLS positive: y = Xb + e with b >= 0
//   -- NNLS negative: y = Xb + e with b <= 0
double jburden_test(const Eigen::VectorXd &y, const Eigen::MatrixXd& X,
  int df = 0, double tol = 1e-6, int n_approx = 100, bool 
  strict = false, int verbose = 0);

// compute exact weights for the NNLS test
Eigen::VectorXd jburden_wts(const Eigen::MatrixXd& V, int verbose = 0);
// compute adaptive weights for the NNLS test
int jburden_wts_adapt(const Eigen::MatrixXd& V, Eigen::VectorXd& wts_out,
    int n_approx = 100, bool normalize = true, int verbose = 0);

// compute CDF for MVN
double jburden_pnorm(const Eigen::MatrixXd& A, 
  int maxpts = 25000, double abseps = 1e-3, int verbose = 0);
// the active set algorithm for fitting NNLS
int jburden_fit_nnls(const Eigen::VectorXd &y, const Eigen::MatrixXd& X, 
  Eigen::VectorXd& bhat_out, vector<bool>& selected_out,
  double tol = 1e-6, bool neg = false, int maxit = 1000, int verbose = 0);
// the number of all set of k out of n
int jburden_choose(int n, int k);
// enumerate all sets of k out of n numbers
int jburden_choose(int n, int k);
double jburden_choose_boost(int n, int k);
void jburden_nchoosek(int n, int k, std::list<std::vector<int>> &ll);
void jburden_nchoosek_sample(int n, int k, int s, list<vector<int>> &ll);

struct FitNNLS 
{
  bool executed;
  bool converged;
  int it;
  VectorXd bhat;
  vector<bool> selected;
  double stat;
  double pval;
};

/*
 * TODO: 
 *  - method to check params (e.g. maxit > 0) and write an error message;
 *  - add epsilon parameter for _npd function; use case epsilon = 0 for testing;
 *  - optimize computation of weights by parallelization;
 */
class NNLS 
{
  public:
    int napprox;
    bool normalize;
    double tol;
    int maxit;
    bool strict;
    int verbose;
    string msg_error;

    // 1. OLS
    int p; // number of independent variables in y ~ X model, i.e. p = ncol(X)
    int df;
    MatrixXd XX;
    double sigma2;
    MatrixXd V;
    VectorXd bhat_ols;
    double stat_ols;
    // 2. Positive-definite V
    MatrixXd Vpd;
    // 3. Weights for NNLS test
    int nw;
    VectorXd wts;
    // 4. NNLS model fits
    FitNNLS fit_pos;
    FitNNLS fit_neg;
    // 5. P-value
    bool best_fit; // 1 = pos, 0 = neg
    double pval_min2;

    void set_defaults();
    void run(const Eigen::VectorXd &y, const Eigen::MatrixXd& X, int df = 0);

    void fit_ols(const Eigen::VectorXd &y, const Eigen::MatrixXd& X, int df = 0);
    void compute_weights();
    void fit_nnls(const Eigen::VectorXd &y, const Eigen::MatrixXd& X);
    void fit_nnls_sign(const Eigen::VectorXd &y, const Eigen::MatrixXd& X, bool neg, struct FitNNLS&);
     
    void print_param() 
    { 
      cout << "NNLS parameters: (weights) napprox = " << napprox << ", normalize = " << normalize 
        << "; (model fitting) tol = " << tol << ", maxit = " << maxit 
        << "; (general) verbose = " << verbose
        << endl;
    }
    void print_results() 
    {
      cout << "NNLS results: (model fitting pos/neg) executed = " << fit_pos.executed << "/" << fit_neg.executed 
        << ", converged = " << fit_pos.converged << "/" << fit_neg.converged 
        << "; (inference) pval_min2 = " << pval_min2 
        << "; (general) error message = \"" << msg_error << "\"" 
        << endl;
    }

    // get NNLS results
    FitNNLS* get_best_fit(bool pos) 
    { 
      FitNNLS* ret = pos ? &fit_pos : &fit_neg; 
      return(ret);
    };

    // print info
    string str_bhat_i(unsigned i) 
    {
      ostringstream buffer;
      if(pval_min2 == -1) {
        buffer << "NA";
      } else {
        if(best_fit) buffer << fit_pos.bhat[i];
        else buffer << fit_neg.bhat[i];
      }
      return(buffer.str());
    }

    string str_sel_i(unsigned i) 
    {
      ostringstream buffer;
      if(pval_min2 == -1) {
        buffer << "NA";
      } else {
        if(best_fit) buffer << fit_pos.selected[i];
        else buffer << fit_neg.selected[i];
      }
      return(buffer.str());
    }

    string str_bhat(bool pos) 
    {
      ostringstream buffer;
      for(int i = 0; i < p; i++) {
        if(pos) buffer << fit_pos.bhat[i] << " ";
        else buffer << fit_neg.bhat[i] << " ";
      }
      buffer << endl;
      return(buffer.str());
    }

    string str_selected(bool pos) 
    {
      ostringstream buffer;
      for(int i = 0; i < p; i++) {
        if(pos) buffer << fit_pos.selected[i] << " ";
        else buffer << fit_neg.selected[i] << " ";
      }
      buffer << endl;
      return(buffer.str());
    }

    string str_wts()
    {
      ostringstream buffer;
      for(int i = 0; i < p; i++) {
        buffer << wts[i] << " ";
      }
      buffer << endl;
      return(buffer.str());
    }

    string str_bhat_ols()
    {
      ostringstream buffer;
      for(int i = 0; i < p; i++) {
        buffer << bhat_ols[i] << " ";
      }
      buffer << endl;
      return(buffer.str());
    }

    string str_info()
    {
      ostringstream buffer;
      buffer << "best_fit " << best_fit << endl 
        << "pval_min2 " << pval_min2 << endl 
        << "wts " << str_wts()
        // NNLS pos
        << "nnls_pos" << endl 
        << "stat " << fit_pos.stat << endl << "pval " << fit_pos.pval << endl
        << "selected_pos " << str_selected(true) << "bhat_nnls_pos " << str_bhat(true)
        // NNLS neg
        << "nnls_neg" << endl 
        << "stat " << fit_neg.stat << endl << "pval " << fit_neg.pval << endl
        << "selected_neg " << str_selected(false) << "bhat_nnls_neg " << str_bhat(false)
        // OLS
        << "ols" << endl 
        << "sigma2 " << sigma2 << endl << "stat_ols " << stat_ols << endl
        << "bhat_ols " << str_bhat_ols();

      return(buffer.str());
    };

    NNLS(); 
    NNLS(int napprox_, bool normalize_, double tol_, bool strict_, int verbose_);
    ~NNLS() { };
};


#endif
