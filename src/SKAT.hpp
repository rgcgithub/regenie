/* 

   This file is part of the regenie software package.

   Copyright (c) 2020-2022 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

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

#ifndef SKAT_H
#define SKAT_H

// SKAT
void update_vc_gmat(SpMat&,Eigen::ArrayXd&,Eigen::ArrayXd&,ArrayXb&,const int&,const int&,struct param const&,const Eigen::Ref<const ArrayXb>&,Eigen::Ref<Eigen::MatrixXd>,std::vector<variant_block>&,const Eigen::Ref<MatrixXb>);
void compute_vc_masks(SpMat&,Eigen::Ref<Eigen::ArrayXd>,Eigen::Ref<Eigen::ArrayXd>,SpMat&,Eigen::Ref<MatrixXb>,const Eigen::Ref<const Eigen::MatrixXd>&, struct ests const&,struct f_ests const&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const MatrixXb>&,MatrixXb&,std::vector<variant_block>&,const Eigen::Ref<const ArrayXb>&,struct param const&);
void prep_ultra_rare_mask(SpMat&,Eigen::Ref<Eigen::ArrayXd>,Eigen::Ref<Eigen::ArrayXd>,SpMat&,Eigen::Ref<MatrixXb>,MatrixXb&,const Eigen::Ref<const ArrayXb>&,struct param const&);
void compute_vc_masks_qt(SpMat&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const MatrixXb>&,std::vector<variant_block>&,struct param const&);
void compute_vc_masks_qt_fixed_rho(SpMat&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const MatrixXb>&,std::vector<variant_block>&,double const&,double const&,double const&,uint const&,bool const&);
void compute_vc_masks_qt(SpMat&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const MatrixXb>&,std::vector<variant_block>&,const Eigen::Ref<const Eigen::ArrayXd>&,double const&,double const&,uint const&,bool const&);
void compute_vc_mats_qt(Eigen::Ref<Eigen::MatrixXd>,Eigen::Ref<Eigen::MatrixXd>,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<SpMat>,const Eigen::Ref<const Eigen::ArrayXd>&);
void get_single_pvs(Eigen::Ref<Eigen::MatrixXd>,const Eigen::Ref<const Eigen::MatrixXd>&);

void compute_vc_masks_bt(SpMat&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,struct ests const&,struct f_ests const&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const MatrixXb>&,std::vector<variant_block>&,struct param const&);
void compute_vc_masks_bt_fixed_rho(SpMat&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,struct ests const&,struct f_ests const&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const MatrixXb>&,std::vector<variant_block>&,double const&,double const&,double const&,bool const&,uint const&,bool const&,struct param const&);
void compute_vc_masks_bt(SpMat&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,struct ests const&,struct f_ests const&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const MatrixXb>&,std::vector<variant_block>&,const Eigen::Ref<const Eigen::ArrayXd>&,double const&,double const&,bool const&,uint const&,bool const&,struct param const&);
void get_single_pvs_bt(Eigen::Ref<Eigen::ArrayXd>,const Eigen::Ref<const Eigen::ArrayXd>&);
Eigen::MatrixXd get_RsKRs(const Eigen::Ref<const Eigen::MatrixXd>&,const double&,const double&);
Eigen::MatrixXd get_RsKRs(const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,const double&,const double&,const double&);
void get_lambdas(Eigen::VectorXd&,const Eigen::Ref<const Eigen::MatrixXd>&,const double&);
void compute_fixed_skato_p(double&,double&,double const&,double const&,double const&,Eigen::VectorXd&,const double&,bool const&);
void compute_fixed_skato_p(double&,double&,double&,double const&,Eigen::VectorXd&,const double&);
void compute_skat_pv(double&,double&,double const&,Eigen::VectorXd&,const double&);
double get_chisq_mix_pv(double const&,const Eigen::Ref<const Eigen::VectorXd>&);
double get_davies_pv(double const&,Eigen::Ref<Eigen::VectorXd>,bool const&);
double get_kuonen_pv(const double&,const Eigen::Ref<const Eigen::VectorXd>&);
double get_liu_pv(const double&,const Eigen::Ref<const Eigen::VectorXd>&);
double get_tmin_lambda(const double&,const Eigen::Ref<const Eigen::ArrayXd>&);
double get_tmax_lambda(const Eigen::Ref<const Eigen::ArrayXd>&);
void solve_kp(bool&,double&,const double&,const double&,const double&,const Eigen::Ref<const Eigen::ArrayXd>&);
bool valid_bounds(double&,double const&,double const&,const double&,const Eigen::Ref<const Eigen::ArrayXd>&);
double K_lambda(const double&,const Eigen::Ref<const Eigen::ArrayXd>&);
double Kp_lambda(const double&,const Eigen::Ref<const Eigen::ArrayXd>&);
double Kpp_lambda(const double&,const Eigen::Ref<const Eigen::ArrayXd>&);
double get_spa_pv(const double&,const double&,const Eigen::Ref<const Eigen::ArrayXd>&);
void apply_correction_cc(const int&,Eigen::Ref<Eigen::ArrayXd>,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,SpMat const&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::MatrixXd>&,SpMat const&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const ArrayXb>&,struct f_ests const&,struct param const&);
void apply_firth_snp(bool&,double&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const ArrayXb>&,struct param const&);
void get_skato_mom(double&,double&,double&,double&,Eigen::ArrayXd&,const Eigen::Ref<const Eigen::VectorXd>&,double const&,double const&,double const&,const Eigen::Ref<const Eigen::ArrayXd>&,bool const&);
void get_cvals(int const&,Eigen::Ref<Eigen::MatrixXd>,const Eigen::Ref<const Eigen::VectorXd>&);
void get_cvals(Eigen::Ref<Eigen::ArrayXd>,const Eigen::Ref<const Eigen::VectorXd>&);
void get_Qmin(int const&,double&,Eigen::Ref<Eigen::ArrayXd>,const Eigen::Ref<const Eigen::MatrixXd>&);
void get_skato_pv(double &,double&,double const&,int const&,double const&,bool const&);
void print_vc_sumstats(int const&,std::string const&,std::string const&,variant_block*,std::vector<snp> const&,struct in_files const&,struct param const*);

// for numerical integration with skat-o
#ifdef __cplusplus
extern "C"
{
#endif

  extern void dqags_(double f(double*),double*,double*,double*,double*,double*,double*,int*,int*,int*,int*,int*,int*,double*);
  double SKATO_integral_fn(double*);
  double SKATO_integral_fn_liu(double*);

#ifdef __cplusplus
}
#endif

// declare global variables
extern Eigen::ArrayXd flipped_skato_rho;
extern Eigen::ArrayXd skato_Qmin_rho;
extern Eigen::ArrayXd skato_tau;
extern Eigen::VectorXd skato_lambdas;
extern double skato_muQ;
extern double skato_fdavies;
extern double skato_sdQ;
extern double skato_dfQ;
extern double skato_upper;
extern int skato_state; // positive if integration failed
void integrate(double f(double*),double&,int const&,bool const&);

#endif
