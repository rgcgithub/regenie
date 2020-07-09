/* 

   This file is part of the regenie software package.

   Copyright (c) 2020 Joelle Mbatchou & Jonathan Marchini

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

#ifndef TEST_MODELS_H
#define TEST_MODELS_H

struct f_ests {

  Eigen::MatrixXd covs_firth;
  Eigen::MatrixXd beta_null_firth;
  double deviance_logistic;
  double bhat_firth, se_b_firth;

};

struct spa_ests {

  Eigen::MatrixXd SPA_pvals;
  Eigen::ArrayXd Gmod;
  double val_a, val_b, val_c, val_d;
  bool pos_score;

};

// firth
bool fit_firth_logistic(int,int,bool,struct param*,struct phenodt*,struct ests*,struct f_ests*,mstream&);
void fit_firth_logistic_snp(int,int,bool,struct param*,struct phenodt*,struct ests*,struct f_ests*,variant_block*,mstream&);

// spa (multithreading in Eigen)
double solve_K1(const double,const bool,const double,const int,const int,const struct param*,const struct ests*,const struct spa_ests*,const struct geno_block*,mstream&);
double compute_K(const double,const int,const struct ests*,const struct spa_ests*);
double compute_K1(const double,const int,const struct ests*,const struct spa_ests*);
double compute_K2(const double,const int,const struct ests*,const struct spa_ests*);
double compute_K_fast(const double,const double,const int,const int,const struct ests*,const struct spa_ests*,const struct geno_block*);
double compute_K1_fast(const double,const double,const int,const int,const struct ests*,const struct spa_ests*,const struct geno_block*);
double compute_K2_fast(const double,const double,const int,const int,const struct ests*,const struct spa_ests*,const struct geno_block*);
double get_SPA_pvalue(const double,const double,const bool,const double,const int,const int,const struct param*,const struct ests*,const struct spa_ests*,const struct geno_block*,mstream&);

// spa (multithreading in openmp)
double solve_K1_snp(const double,const int,const struct param*,const struct ests*,variant_block*,mstream&);
double compute_K_snp(const double,const struct ests*,variant_block*,const int);
double compute_K1_snp(const double,const struct ests*,variant_block*,const int);
double compute_K2_snp(const double,const struct ests*,variant_block*,const int);
double compute_K_fast_snp(const double,const struct ests*,variant_block*,const int);
double compute_K1_fast_snp(const double,const struct ests*,variant_block*,const int);
double compute_K2_fast_snp(const double,const struct ests*,variant_block*,const int);
double get_SPA_pvalue_snp(const double,const double,const int,struct param*,const struct ests*,variant_block*,mstream&); 

#endif
