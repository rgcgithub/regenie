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

#ifndef MODELS_H
#define MODELS_H

#define ETAMINTHR -30.0
#define ETAMAXTHR 30.0

struct ests {

  Eigen::MatrixXd offset_nullreg;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> blups;
  Eigen::MatrixXd ltco_prs;
  Eigen::MatrixXd Gamma_sqrt;
  Eigen::MatrixXd Y_hat_p;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > Xt_Gamma_X_inv;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > X_Gamma;
  Eigen::MatrixXd bhat_start; // for interaction tests

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


void fit_null_logistic(const int&,struct param*,struct phenodt*,struct ests*,struct in_files*,mstream&);
bool fit_logistic(const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const ArrayXb>&,Eigen::ArrayXd&,Eigen::ArrayXd&,Eigen::ArrayXd&,struct param const*,mstream&);
double get_logist_dev(const Eigen::Ref<const Eigen::ArrayXd>& Y, const Eigen::Ref<const Eigen::ArrayXd>& pi, const Eigen::Ref<const ArrayXb>& mask);

void fit_null_poisson(const int&,struct param*,struct phenodt*,struct ests*,struct in_files*,mstream&);
bool fit_poisson(const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const ArrayXb>&,Eigen::ArrayXd&,Eigen::ArrayXd&,Eigen::ArrayXd&,struct param const*,mstream&);
double get_poisson_dev(const Eigen::Ref<const Eigen::ArrayXd>& Y, const Eigen::Ref<const Eigen::ArrayXd>& pi, const Eigen::Ref<const ArrayXb>& mask);

void ridge_level_0(const int&,struct in_files*,struct param*,struct filter*,struct ests*,struct geno_block*,struct phenodt*,std::vector<snp>&,struct ridgel0*,struct ridgel1*,std::vector<MatrixXb>&,mstream&);
void ridge_level_0_loocv(const int,struct in_files*,struct param*,struct filter*,struct ests*,struct geno_block*,struct phenodt*,std::vector<snp>&,struct ridgel0*,struct ridgel1*,mstream&);
void write_l0_file(std::ofstream*,Eigen::MatrixXd&,mstream&);

void set_mem_l1(struct in_files*,struct param*,struct filter*,struct ests*,struct geno_block*,struct phenodt*,struct ridgel1*,std::vector<MatrixXb>&,mstream&);
void ridge_level_1(struct in_files*,struct param*,struct ridgel1*,mstream&);
void ridge_level_1_loocv(struct in_files*,struct param*,struct phenodt*,struct ridgel1*,mstream&);

void ridge_logistic_level_1(struct in_files*,struct param*,struct phenodt*,struct ridgel1*,std::vector<MatrixXb>&,mstream&);
void ridge_logistic_level_1_loocv(struct in_files*,struct param*,struct phenodt*,struct ests*,struct ridgel1*,mstream&);
bool run_log_ridge_loocv(const double&,const int&,const int&,Eigen::ArrayXd&,Eigen::ArrayXd&,Eigen::ArrayXd&,const Eigen::Ref<const Eigen::ArrayXd>&,Eigen::Ref<Eigen::MatrixXd>,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const ArrayXb>&,struct param*,mstream&);
void run_log_ridge_loocv_adam(const int&,const double&,Eigen::ArrayXd&,Eigen::ArrayXd&,Eigen::ArrayXd&,const Eigen::Ref<const Eigen::ArrayXd>&,Eigen::Ref<Eigen::MatrixXd>,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const ArrayXb>&,struct param*,mstream&);

void ridge_poisson_level_1(struct in_files*,struct param*,struct phenodt*,struct ridgel1*,std::vector<MatrixXb>&,mstream&);
void ridge_poisson_level_1_loocv(struct in_files*,struct param*,struct phenodt*,struct ests*,struct ridgel1*,mstream&);
bool run_ct_ridge_loocv(const double&,const int&,const int&,Eigen::ArrayXd&,Eigen::ArrayXd&,const Eigen::Ref<const Eigen::ArrayXd>&,Eigen::Ref<Eigen::MatrixXd>,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const ArrayXb>&,struct param*,mstream&);

bool get_wvec(Eigen::ArrayXd&,Eigen::ArrayXd&,const Eigen::Ref<const ArrayXb>&,const double&);
void get_pvec(Eigen::ArrayXd&,Eigen::ArrayXd&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,double const&);
void get_pvec_poisson(Eigen::ArrayXd&,Eigen::ArrayXd&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,double const&);
double compute_log_lik_bern(const double&,const double&);
double compute_log_lik_poisson(const double&,const double&);

void read_l0(int const&,int const&,struct in_files*,struct param*,struct ridgel1*,mstream&);
void read_l0_chunk(int const&,int const&,int const&,int const&,const std::string&,struct param*,struct ridgel1*,mstream&);

uint64 getSize(std::string const& fname);
#endif

