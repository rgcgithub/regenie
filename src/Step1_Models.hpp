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

#ifndef MODELS_H
#define MODELS_H

struct ests {

  Eigen::MatrixXd offset_logreg;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> blups;
  Eigen::MatrixXd Gamma_sqrt;
  Eigen::MatrixXd Y_hat_p;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > Xt_Gamma_X_inv;

};

struct ridgel0 {

  Eigen::MatrixXd GGt;
  Eigen::MatrixXd GTY;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > G_folds; // storage for folds at levle 0
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > GtY; // storage for folds at levle 0
  Eigen::MatrixXd GGt_eig_vec, GGt_eig_val;
  Eigen::MatrixXd Wmat;

};

struct ridgel1 {

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > X_folds, XtY; // storage for folds at level 1
  std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > > pred_mat, test_mat;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > test_mat_conc;
  std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > > pred_pheno, test_pheno;
  std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > > pred_pheno_raw, test_pheno_raw;
  std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > > pred_offset, test_offset;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > cumsum_values; // storage of values to compute rsq and values [Sx, Sy, Sx2, Sy2, Sxy]
  std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > > beta_hat_level_1;
  ArrayXb pheno_l1_not_converged;

};


void fit_null_logistic(const int,struct param*,struct phenodt*,struct ests*,mstream&);
double get_logist_dev(const Eigen::ArrayXd& Y, const Eigen::ArrayXd& pi, const ArrayXb& mask);

void ridge_level_0(const int,struct in_files*,struct param*,struct filter*,struct ests*,struct geno_block*,struct phenodt*,std::vector<snp>&,struct ridgel0*,struct ridgel1*,std::vector<MatrixXb>&,mstream&);
void ridge_level_0_loocv(const int,struct in_files*,struct param*,struct filter*,struct ests*,struct geno_block*,struct phenodt*,std::vector<snp>&,struct ridgel0*,struct ridgel1*,mstream&);

void ridge_level_1(struct in_files*,struct param*,struct ridgel1*,mstream&);
bool get_wvec_fold(int ph, Eigen::ArrayXd& etavec, Eigen::ArrayXd& pivec, Eigen::ArrayXd& wvec, const Eigen::ArrayXd& beta, const MatrixXb& masks, const Eigen::MatrixXd& offset, const Eigen::MatrixXd& test_mat);
void ridge_level_1_loocv(struct in_files*,struct param*,struct phenodt*,struct ridgel1*,mstream&);

void ridge_logistic_level_1(struct in_files*,struct param*,struct phenodt*,struct ridgel1*,std::vector<MatrixXb>&,mstream&);
void ridge_logistic_level_1_loocv(struct in_files*,struct param*,struct phenodt*,struct ests*,struct ridgel1*,mstream&);

double compute_log_lik(const double,const double);

#endif

