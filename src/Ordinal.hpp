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

#ifndef ORDINAL_H
#define ORDINAL_H

//-------------------
// Class FitOrdinal
//-------------------

class FitOrdinal 
{
  public:
    unsigned int cnt_fit;
    unsigned int verbose;
    bool executed, converged;
    bool trace;
    unsigned int it, it2, maxit, maxit2, maxit3, cnt_updates;
    bool strict;
    double stop_value, tol;
    double pseudo_stophalf;

    bool check_step;
    double max_step;

    bool apply_start; // apply starting values?
    bool store_offset; // store offset?
    bool apply_offset; // apply offset?
    Eigen::VectorXd yo; // offset vector (linear predictor Xb)
    Eigen::VectorXd yo_int; // offset intercepts for Multinomial model
    bool exclude_intercepts; // exclude intercepts, e.g., offset is given
    bool exclude_intercepts_offset; // exclude intercepts when applying offset
    Eigen::VectorXd bhat;
    double loglik;

    Eigen::VectorXd Score;
    Eigen::MatrixXd Info;
                   
    std::string response; // response type = [binom, multinom]
    std::string model; // model = [POM: Proportional Odds Model, ACL: Adjacent Category Logit]
    std::string optim; // optimization algorithm = [FisherScoring, WeightHalving, Pseudo]
    bool firth_binom, firth_multinom; // Firth correction
    double firth_mult;

    unsigned int N, Neff; // sample size
    unsigned int Ncov, Nb; // number of covariates in X & number of covariates + all intercepts
    unsigned int  Ncov0, Ncov1; // number of last covariates in X: the effects are fixed to zero
    bool last0;
                    
    unsigned int ncat, ncat1, ncat1sq; // number of categories 
    Eigen::VectorXi Ncat;
                              
    Eigen::VectorXd b0;

    // current results in model fitting
    Eigen::MatrixXd cur_Info;
    Eigen::VectorXd cur_Score, cur_v, cur_b;
    double cur_loglik, cur_dev, prev_dev;

    // intermediate matrices/vectors for model fitting
    // common + multinomial
    Eigen::VectorXd Xb0;
    Eigen::MatrixXd Xb, exp_eta, gamma, PQ;
    Eigen::MatrixXd P;
    Eigen::VectorXd Psum, Pk;
    Eigen::MatrixXd D, Q, V, S, QS, W;
    Eigen::VectorXd WS2, WSS1; // column sums (dim 2) and sqrt of row sums (dim 1)
    Eigen::MatrixXd XW, XWs, XW22;
    Eigen::VectorXd XW1;
    // binomial
    Eigen::VectorXd mub, wb;
    Eigen::MatrixXd XtW;
    // firth
    Eigen::VectorXd ystar;
    Eigen::MatrixXd Ystar;

    void setup_defaults();
    void check_setup_model();
    void check_setup_model_common();
    void check_setup_data();
    void check_setup_data_common();
    void check_setup_data_binom();
    void check_setup_data_multinom();
    void setup_ncov0(unsigned int , bool, bool );

    void setup_restart(const Eigen::VectorXd &);
    void setup_offset_binom(const Eigen::VectorXd &, bool);
    void setup_offset_multinom_pom(const Eigen::VectorXd & , const Eigen::VectorXd & );

    // optimization functions
    bool stop_criterion();
    bool update_par(const VectorXb & , const MatrixXb & , const Eigen::Ref<const Eigen::MatrixXd> & , const Eigen::VectorXd & , bool pseudo = false);
    bool optimize(const VectorXb & , const MatrixXb & , const Eigen::Ref<const Eigen::MatrixXd> & X);
    bool optimize_FisherScoring(const VectorXb & , const MatrixXb & , const Eigen::Ref<const Eigen::MatrixXd> & X);
    bool optimize_FisherScoringPseudo(const VectorXb & , const MatrixXb & , const Eigen::Ref<const Eigen::MatrixXd> & X);
    bool optimize_WeightHalving(const VectorXb & , const MatrixXb & , const Eigen::Ref<const Eigen::MatrixXd> & X);
    bool optimize_WeightHalvingPseudo(const VectorXb & , const MatrixXb & , const Eigen::Ref<const Eigen::MatrixXd> & X);
    void update_fit();

    // multinom (POM only)
    void fit_multinom_pom(const VectorXb & , const MatrixXb & , const Eigen::Ref<const Eigen::MatrixXd> & );
    void setup_start_multinom();
    void setup_par_multinom();
    bool update_par_multinom(const VectorXb & , const MatrixXb & , const Eigen::Ref<const Eigen::MatrixXd> & , const Eigen::VectorXd & , bool pseudo = false);
    double loglik_multinom(const VectorXb &, const MatrixXb &);
    double loglik_multinom_firth(const VectorXb &, const MatrixXb &, const Eigen::LLT<Eigen::MatrixXd> &, bool add = false, double base_loglik = 0.0);
  /* cur_loglik = firth_multinom ? loglik_multinom_firth(Mask, Y, cur_Info) : loglik_multinom(Mask, Y); */
  /* cur_loglik =  loglik_multinom(Mask, Y); */
    double loglik_multinom_firth(const VectorXb &, const MatrixXb &, const Eigen::MatrixXd &);
    /* void update_fit_multinom(FitOrdinal & ); */

    // binom (logistic regression)
    void fit_binom(const VectorXb &, const MatrixXb &, const Eigen::Ref<const Eigen::MatrixXd> &);
    void setup_start_binom();
    void setup_par_binom();
    bool update_par_binom(const VectorXb & , const MatrixXb & , const Eigen::Ref<const Eigen::MatrixXd> & , const Eigen::VectorXd & );
    bool update_par_binom_firth(const VectorXb & , const MatrixXb & , const Eigen::Ref<const Eigen::MatrixXd> & , const Eigen::VectorXd & );
    bool update_par_binom_pseudo(const VectorXb & , const MatrixXb & , const Eigen::Ref<const Eigen::MatrixXd> & , const Eigen::VectorXd & );
    double loglik_binom(const VectorXb &, const MatrixXb &);
    double loglik_binom_firth(const VectorXb &, const MatrixXb &, const Eigen::LLT<Eigen::MatrixXd> &);

    void setup_xy(const Eigen::VectorXi &, const Eigen::MatrixXd& );

    FitOrdinal();
    ~FitOrdinal() { };
};

//-------------------------------------------------------------
// Class Ordinal
//-------------------------------------------------------------
class Ordinal
{
  public:

    std::string response; // response type = [binom, multinom]
    std::string model; // model = [POM: Proportional Odds Model, ACL: Adjacent Category Logit]
    std::string optim; // optimization algorithm = [FisherScoring, WeightHalving]
    bool firth_binom, firth_multinom; // Firth correction

    unsigned int it, it2, maxit, maxit2, maxit3;
    bool strict;
    double tol;
    double pseudo_stophalf;
    bool converged;

    bool check_step;
    double max_step;
                
    bool inc_intercept; // is the interecept term in X included or not
    bool preproc_cov; // center/orth. covariates X?

    unsigned int N, Neff; // sample size
    VectorXb Mask; // masked samples = samples with missing values 
    unsigned int Ncov, Nb; // number of covariates in X & number of covariates + all intercepts
                    
    std::set<int> cat; // ordered (!) set of category levels
    unsigned int ncat, ncat1, ncat1sq; // number of categories 
    Eigen::VectorXi Ncat;
                              
    MatrixXb Y; // N x ncat matrix: Y[i,j] = sample i is in j category
    Eigen::MatrixXd Xcov; // X = U S V'
    Eigen::MatrixXd Xcov1; // [Intercept, X] = [1s, X]
    Eigen::VectorXd b0;

    // intermediate matrices/vectors for model fitting
    // common + multinomial
    Eigen::VectorXd Xb0;
    Eigen::MatrixXd Xb, exp_eta, gamma, PQ;
    Eigen::MatrixXd P;
    Eigen::VectorXd Psum, Pk;
    Eigen::MatrixXd D, Q, V, S, QS, W;
    Eigen::VectorXd WS2, WSS1; // column sums (dim 2) and sqrt of row sums (dim 1)
    Eigen::MatrixXd XW, XWs, XW22;
    Eigen::VectorXd XW1;

    // intermediate matrices/vectors for model fitting
    // binomial
    Eigen::VectorXd yb;
    Eigen::VectorXd mub, wb;
    Eigen::MatrixXd XtW;

    Eigen::MatrixXd cur_Info;
    Eigen::VectorXd cur_Score, cur_v, cur_b;
    double cur_loglik, cur_dev, prev_dev;

    void setup_defaults();
    FitOrdinal setup_fit();
    void setup_xy(const Eigen::VectorXi &, const Eigen::MatrixXd& );

    FitOrdinal fit(const Eigen::VectorXi &y, const Eigen::MatrixXd& X, 
      unsigned int Ncov0 = 0, bool last0 = true);
    void update_par(const Eigen::VectorXd &, const Eigen::MatrixXd &);
    void update_fit_common(FitOrdinal & );
    bool stop_criterion();
    bool optimize(const Eigen::MatrixXd &);
    bool optimize_FisherScoring(const Eigen::MatrixXd &);
    bool optimize_WeightHalving(const Eigen::MatrixXd &);

    // multinom (POM only)
    void fit_multinom_pom(const Eigen::MatrixXd& , FitOrdinal & );
    void setup_start_multinom();
    void setup_par_multinom();
    void update_par_multinom(const Eigen::VectorXd &, const Eigen::MatrixXd &);
    double loglik_multinom();
    void update_fit_multinom(FitOrdinal & );

    // binom (logistic regression)
    void fit_binom(const Eigen::MatrixXd& , FitOrdinal & );
    void setup_start_binom();
    void setup_par_binom();
    void update_par_binom(const Eigen::VectorXd &, const Eigen::MatrixXd &);
    void update_par_binom_firth(const Eigen::VectorXd &, const Eigen::MatrixXd &);
    double loglik_binom();
    void update_fit_binom(FitOrdinal & );

    double test_score(const FitOrdinal &, const Eigen::MatrixXd &, const Eigen::MatrixXd &);
    double test_score(const FitOrdinal &, const Eigen::MatrixXd &);
    double test_score_multinom_pom(const FitOrdinal &, const Eigen::MatrixXd &, const Eigen::MatrixXd &);
    double test_score_binom(const FitOrdinal &, const Eigen::MatrixXd &, const Eigen::MatrixXd &);

    Ordinal();
    ~Ordinal() { };
};

//-------------------
// Class MultiPhen
//-------------------

class MultiPhen 
{
  public:
    unsigned int cnt_fit; // fit index/counter
    unsigned int verbose; // verbose level
    std::string response; // response type = [binom, multinom]
    std::string model; // model = [POM: Proportional Odds Model, ACL: Adjacent Category Logit]
    std::string optim; // optimization algorithm = [FisherScoring, WeightHalving]
    bool firth_binom, firth_multinom; // Firth correction
    double firth_mult; 
    bool reuse_start, reset_start;
    bool approx_offset;
    int mac_approx_offset;
    Eigen::VectorXd b0;
    Eigen::VectorXd yo;
    Eigen::VectorXd yo_int; 
    Eigen::VectorXd w0;
    Eigen::MatrixXd Yres0;

    unsigned int maxit, maxit2, maxit3; // maximum number of interations in Outer/Inner loops of the IRLS
    bool strict;
    double tol;
    double pseudo_stophalf;
    bool check_step;
    double max_step;
    std::string offset_mode;

    bool trace; // trace updates
    unsigned int it, cnt_updates;
                          
    unsigned int N, Neff; // sample size
    VectorXb Mask; // masked samples = samples with missing values 
    unsigned int Nx, Nx1, Ny, Ny1, Ny21; // columns in XY matrix
    /* unsigned int Ncov, Nb; // number of covariates in X & number of covariates + all intercepts */
    /* unsigned int  Ncov0, Ncov1; // number of last covariates in X: the effects are fixed to zero */
                    
    bool pos_intercept_first, pos_phen_first;

    unsigned int ncat, ncat1, ncat1sq; // number of categories 
    Eigen::VectorXi Ncat;
    int Ncat_minor;

    MatrixXb Ym; // N x ncat matrix: Y[i,j] = sample i is in j category. Response for multinomial model.
    Eigen::VectorXd yb;

    bool set_x, set_y;

    bool executed, converged;
    double pval_test;
    double pval_thr;
    std::string test; // test type = [none, nocov_score, nocov_addcov]
    Eigen::VectorXd bhat_y;

    void setup_x(const VectorXb& , const Eigen::MatrixXd& , unsigned int , unsigned int , bool , bool );
    void setup_y(const Eigen::VectorXd & );
    void setup_approx_offset();
    void reset_model();

    FitOrdinal setup_fit(bool , bool , bool use_offset = false);
    FitOrdinal fit(const Eigen::Ref<const Eigen::MatrixXd> & , bool , bool , bool use_res = false);

    void run(const Eigen::VectorXd & , const Eigen::MatrixXd & , unsigned int , unsigned int );
    void run_test_offset(const Eigen::Ref<const Eigen::MatrixXd> &);
    void run_test_qt(const Eigen::Ref<const Eigen::MatrixXd> &);
    void run_test_score(const Eigen::Ref<const Eigen::MatrixXd> & , bool );
    void run_test_lrt(const Eigen::Ref<const Eigen::MatrixXd> & , bool );
    void run_test_addcov(const Eigen::Ref<const Eigen::MatrixXd> & );
    void run_test_add_offset(const Eigen::Ref<const Eigen::MatrixXd> & );

    void run0(const Eigen::VectorXi & , const Eigen::MatrixXd& , const Eigen::MatrixXd& , bool );

    void test0(const Eigen::VectorXi & g, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
      bool firth_binom = false,
      std::string optim = "WeightHalving", double tol = 1e-4,
      unsigned int maxit = 100, bool check_step = true, double max_step = 10.0);
    void test_addcov(const Eigen::VectorXi & g, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
      bool firth_binom = false,
      std::string optim = "WeightHalving", double tol = 1e-4,
      unsigned int maxit = 100, bool check_step = true, double max_step = 10.0);

    // score test
    double test_score(const FitOrdinal & , const VectorXb & , const MatrixXb & , const Eigen::VectorXd & , const Eigen::Ref<const Eigen::MatrixXd> & , bool );
    double test_score_binom(const FitOrdinal & , const VectorXb & , const Eigen::VectorXd & , const Eigen::Ref<const Eigen::MatrixXd> & , const Eigen::Ref<const Eigen::MatrixXd> & );
    double test_score_multinom_pom(const FitOrdinal & , const VectorXb & , const MatrixXb & , const Eigen::Ref<const Eigen::MatrixXd> & , const Eigen::Ref<const Eigen::MatrixXd> & );

    std::string print_sumstats(int const& snp_index, uint32_t const&, 
        std::string const& test_string, std::string const& wgr_string, variant_block* block_info, 
        std::vector<snp> const& snpinfo, struct param const* params);
    std::string print_sum_stats_htp(const variant_block* , struct param const* );

    void setup_defaults();
    MultiPhen();
    MultiPhen(std::string );
    MultiPhen(unsigned int );
    ~MultiPhen() { };
};
 
#endif
