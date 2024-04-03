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

#include "Regenie.hpp"
#include "Files.hpp"
#include "Geno.hpp"
#include "Pheno.hpp"

#include "Ordinal.hpp"

using namespace Eigen;
using namespace std;

//-----------------
// Local functions
//-----------------

Eigen::MatrixXd orth_matrix(const Eigen::MatrixXd & , const MatrixXb &);
void exp_matrix(Eigen::MatrixXd &);
void exp_matrix_ord(Eigen::MatrixXd &);
void exp_vector(Eigen::VectorXd &);
Eigen::VectorXd dlog_vector(const Eigen::VectorXd & );
Eigen::MatrixXd dlog_matrix(const Eigen::MatrixXd & );
bool check_nan(double );

//-------------------
// Class MultiPhen
//-------------------

void MultiPhen::setup_defaults()
{
  // settings
  cnt_fit = 0;
  verbose = 0;
  response = "unknown";
  optim = "WeightHalving";
  firth_binom = false; firth_multinom = false;
  firth_mult = 1.0;
  reuse_start = false; reset_start = false;
  approx_offset = false;
  mac_approx_offset = 0;
  offset_mode = "offset";
  maxit = 150; maxit2 = 10; maxit3 = 10; strict = false;
  tol = 1e-4; pseudo_stophalf = 0.0;
  check_step = true; max_step = 10.0;
  // statuses
  set_x = false; set_y = false;
  // data dimenstions
  N = 0; Neff = 0; // sample size
  /* Ncov = 0, Nb = 0, Ncov0 = 0; Ncov1 = 0; // number of covariates */ 
  ncat = 0, ncat1 = 0, ncat1sq = 0; // number of categories 
  // tests
  pval_thr = 0.1;
  // model fitting results
  executed = false; converged = false;
  trace = false;
  it = 0; cnt_updates = 0;
}

MultiPhen::MultiPhen() 
{
  setup_defaults();
  test = "none";
}

MultiPhen::MultiPhen(std::string _test) 
{
  setup_defaults();
  test = _test;
}

MultiPhen::MultiPhen(unsigned int test_code) 
{
  setup_defaults();

  std::map<unsigned int, std::string> test_map = { {0, "none"}, {1, "cov_score"} };
  test = test_map[test_code];
}

// constructor for FitOrdinal
// - copy model parameters to FitOrdinal object
FitOrdinal MultiPhen::setup_fit(bool inc_cov, bool inc_phen, bool use_offset)
{
  FitOrdinal fit;

  // copy parameters from Ordinal
  fit.verbose = verbose; 
  fit.response = response; // response type = [binom, multinom]
  fit.model = model; // model = [POM: Proportional Odds Model, ACL: Adjacent Category Logit]
  fit.optim = optim; // optimization algorithm = [FisherScoring, WeightHalving]
  fit.firth_binom = firth_binom; fit.firth_multinom = firth_multinom; // Firth correction
  fit.firth_mult = firth_mult; 
      
  fit.maxit = maxit; fit.maxit2 = maxit2; fit.maxit3 = maxit3; fit.strict = strict;
  fit.tol = tol; fit.pseudo_stophalf = pseudo_stophalf;

  fit.check_step = check_step;
  fit.max_step = max_step;

  fit.N = N; fit.Neff = Neff; // samples size
  if(use_offset) {
    // use offset
    if(inc_cov) {
      if(inc_phen) {
        fit.Ncov = Ny; fit.Nb = Ny; // number of covariates
      } else {
        throw std::runtime_error("use offset with covariates only (Ncov = Nb = 0)");
      }
    } else {
      if(inc_phen) {
        fit.Ncov = Ny; fit.Nb = Ny; // number of covariates
      } else {
        throw std::runtime_error("use offset with covariates only (Ncov = Nb = 0)");
      }
    }
  } else {
    // no offset
    if(inc_cov) {
      if(inc_phen) {
        fit.Ncov = Nx + Ny; fit.Nb = ncat1 + Nx + Ny; // number of covariates
      } else {
        fit.Ncov = Nx; fit.Nb = ncat1 + Nx; // number of covariates
      }
    } else {
      if(inc_phen) {
        fit.Ncov = Ny; fit.Nb = ncat1 + Ny; // number of covariates
      } else {
        fit.Ncov = 0; fit.Nb = ncat1; // number of covariates
      }
    }
  }

  fit.ncat = ncat; // number of categories
  fit.ncat1 = ncat1; fit.ncat1sq = ncat1sq;
  fit.Ncat = Ncat;

  fit.cur_dev = 0; fit.prev_dev = 0;
  fit.trace = trace;
  fit.it = 0; fit.it2 = 0; fit.cnt_updates = 0;

  fit.cnt_fit = cnt_fit++;

  return(fit);
}

void MultiPhen::run(const Eigen::VectorXd & g, 
  const Eigen::MatrixXd& XYR, unsigned int n_cov, unsigned int n_phen)
{
  reset_model();

  // check if XYR is set up
  if(!set_x) throw std::runtime_error("run: set_x is false");
  // set y
  setup_y(g); // -> Ym, yb
  if(!set_y) return; // early stop (example #cat = 1 for imputed variant due to rounding)
  setup_approx_offset(); // approx_offset
  // print info
  if(verbose) cout << "MultiPhen: Nx = " << Nx << " Ny = " << Ny << endl;

  // test
  if(test == "none") {
    reset_model();
    // do nothing
  } else if(test == "cov_score_it1") {
    maxit = 1; optim = "FisherScoring";
    run_test_score(XYR, true); // inc_cov = false
  } else if(test == "nocov_score") {
    run_test_score(XYR, false); // inc_cov = false
  } else if(test == "cov_score") {
    run_test_score(XYR, true); // inc_cov = true
  } else if(test == "nocov_lrt") {
    run_test_lrt(XYR, false); // inc_cov = false
  } else if(test == "cov_lrt") {
    run_test_lrt(XYR, true); // inc_cov = true
  } else if(test == "offset") {
    run_test_offset(XYR);
  } else if(test == "nocov_score_addcov") {
    run_test_addcov(XYR);
  } else if(test == "nocov_score_offset") {
    run_test_add_offset(XYR);
  } else {
    throw std::runtime_error("run: unknown test");
  }
}

void MultiPhen::run0(const Eigen::VectorXi & g, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, bool score_lrt)
{
  // set up Ordinal model (no Firth)
  Ordinal ord;
  ord.optim = optim; ord.tol = tol; ord.pseudo_stophalf = pseudo_stophalf; ord.maxit = maxit; ord.maxit2 = maxit2; ord.maxit3 = maxit3; ord.strict = strict;
  ord.check_step = check_step; ord.max_step = max_step;
  ord.firth_binom = false; 

  if(score_lrt) { // Score test
    executed = true; converged = false; pval_test = -1.0;
    FitOrdinal fit;
  
    // fit null model
    fit = ord.fit(g, X);
    if(!fit.converged) { return; }

    /* // run Score test */
    converged = fit.converged;
    pval_test = ord.test_score(fit, Y);
  } else { // LRT
    executed = true; converged = false; pval_test = -1.0;
    FitOrdinal fit0, fit1;
  
    // prepare new matrix of covariates X + Y
    MatrixXd X1(Y.rows(), X.cols() + Y.cols());
    if(X.cols()) X1.leftCols(X.cols()) = X;
    X1.rightCols(Y.cols()) = Y;
  
    // fit null model
    fit0 = ord.fit(g, X);
    if(!fit0.converged) { return; }

    // fit alternative model (Firth)
    fit1 = ord.fit(g, X1);
    if(!fit1.converged) { return; }
    converged = fit1.converged;

    boost::math::chi_squared dist(Y.cols());
    double stat_lrt = 2 * (fit1.loglik - fit0.loglik);
    pval_test = boost::math::cdf(boost::math::complement(dist, stat_lrt));
  }
}

// XYR = [Intercept, X, Y, Inercept, R]
FitOrdinal MultiPhen::fit(const Eigen::Ref<const Eigen::MatrixXd> & XYR, bool inc_cov, bool inc_phen, bool use_res)
{
  if(use_res) throw std::runtime_error("use_res is not implemented yet");

  // initialize defaults settings 
  bool inc_phen_null = false, inc_phen_firth = inc_phen;
  bool use_offset = (inc_phen && approx_offset);
  bool copy_start = (reuse_start && inc_cov && inc_phen && !approx_offset);
  // update settings for Binom: no firth / firth
  if(response == "binom") {
    inc_phen_null = firth_binom && !inc_phen && !approx_offset;
    inc_phen_firth = inc_phen_null ? true : inc_phen;
  } 
  // update settings for Multinom: no firth / firth
  if(response == "multinom") {
    inc_phen_null = firth_multinom && !inc_phen && !approx_offset;
    inc_phen_firth = inc_phen_null ? true : inc_phen;
  } 
  // create a fit object
  FitOrdinal fit = setup_fit(inc_cov, inc_phen_firth, use_offset);
  /* cout << "done MultiPhen setup_fit: response = " << response */ 
  /*   << " Nx = " << Nx << " Ny = " << Ny */ 
  /*   << " inc_cov = " << inc_cov << " inc_phen_firth = " << inc_phen_firth */ 
  /*   << " use_offset = " << use_offset << " fit.Nb = " << fit.Nb << " fit.Ncov = " << fit.Ncov << endl; */

  // reuse starting par. values
  if(copy_start) fit.setup_restart(b0);

  // refine fit for Binom only
  if(response == "binom") {
    // constraint some par. to zero?
    bool reverse_last = firth_binom && !inc_cov && inc_phen_firth;
    bool last0 = !reverse_last;
    if(inc_phen_null) fit.setup_ncov0(Ny, last0, false); // preproc_cov = false
  }
  // refine fit for Multinom only
  if(response == "multinom") {
    // constraint some par. to zero?
    bool reverse_last = firth_multinom && !inc_cov && inc_phen_firth;
    bool last0 = !reverse_last;
    if(inc_phen_null) fit.setup_ncov0(Ny, last0, false); // preproc_cov = false
  }

  // store offset?
  if(!inc_phen && approx_offset) fit.store_offset = true;
  // apply offset?
  if(use_offset) {
    if(response == "binom") fit.setup_offset_binom(yo, false); // decrement_Nb = false
    else if (response == "multinom") fit.setup_offset_multinom_pom(yo, yo_int);
    else throw std::runtime_error("unknown response");
  }

  // do model fitting & control the columns in XYR passed
  if(response == "binom") {
    if(use_offset) { 
      if(inc_phen_firth) {
        /* fit.fit_binom(Mask, Ym, XYR.rightCols(Ny21).leftCols(Ny)); // matrix of phenotypes Y */ 
        fit.fit_binom(Mask, Ym, Yres0); // matrix of phenotypes Y 
      } else throw std::runtime_error("use offset for the null model");
    } else { 
      if(inc_cov) {
        if(inc_phen_firth) fit.fit_binom(Mask, Ym, XYR.leftCols(Nx1 + Ny)); // X + Y + Intercept
        else fit.fit_binom(Mask, Ym, XYR.leftCols(Nx1)); // X + Intercept
      } else {
        if(inc_phen_firth) fit.fit_binom(Mask, Ym, XYR.rightCols(Ny21).leftCols(Ny1)); // matrix of phenotypes Y + Intercept
        else fit.fit_binom(Mask, Ym, XYR.leftCols(1)); // Intercept
      }
    }
  } else if(response == "multinom") {
    if(use_offset) {
      if(inc_phen_firth) fit.fit_multinom_pom(Mask, Ym, XYR.rightCols(Ny21).leftCols(Ny)); // matrix of phenotypes Y 
      else throw std::runtime_error("use offset for the null model");
    } else {
      if(inc_cov) {
        if(inc_phen_firth) fit.fit_multinom_pom(Mask, Ym, XYR.leftCols(Nx1 + Ny).rightCols(Nx + Ny)); // X + Y
        else fit.fit_multinom_pom(Mask, Ym, XYR.leftCols(Nx1).rightCols(Nx)); // X
      } else {
        if(inc_phen_firth) fit.fit_multinom_pom(Mask, Ym, XYR.rightCols(Ny1).leftCols(Ny)); // matrix of phenotypes Y 
        else fit.fit_multinom_pom(Mask, Ym, XYR.leftCols(0)); // 0 columns
      }
    }
  } else {
    throw std::runtime_error("unknown response");
  }

  if(trace) {
    cnt_updates += fit.cnt_updates;
    it += fit.it;
  }

  return(fit);
}

void MultiPhen::run_test_addcov(const Eigen::Ref<const Eigen::MatrixXd> & XYR)
{
  run_test_score(XYR, false); // inc_cov = false
  if(pval_test < pval_thr) {
    run_test_lrt(XYR, true); // inc_cov = true
  }
}

void MultiPhen::run_test_add_offset(const Eigen::Ref<const Eigen::MatrixXd> & XYR)
{
  run_test_score(XYR, false); // inc_cov = false
  if(pval_test < pval_thr) {
    run_test_offset(XYR);
  }
}

void MultiPhen::run_test_offset(const Eigen::Ref<const Eigen::MatrixXd> & XYR)
{
  FitOrdinal null0, null, full;
  VectorXd b0_fit;
  double ll_null, ll_full;
  boost::math::chi_squared dist(Ny);
  double stat_lrt;
  unsigned int i;

  reset_model(); // reset model fit results
                 
  if(response == "binom") {
    executed = true; 

    // fit null model
    null0 = setup_fit(true, false, false); // inc_cov = true, inc_phen = false, use_offset = false
    null0.store_offset = true;
    null0.fit_binom(Mask, Ym, XYR.leftCols(Nx1)); // covariates X + Intercept

    if(trace) { cnt_updates += null0.cnt_updates; it += null0.it; }

    if(!null0.converged) return;

    // store offset/weights from the null model
    yo = null0.yo;
    yo_int = null0.yo;
    yo_int.array() -= null0.bhat(0); // substract intercept bhat
    w0 = null0.wb;

    // residualize phenotypes
    Yres0 = XYR.rightCols(Ny21).leftCols(Ny); // matrix of phenotypes Y 
    ColPivHouseholderQR<MatrixXd> qrXw;
    qrXw.compute(MatrixXd(Nx1, Nx1).setZero().selfadjointView<Lower>().rankUpdate((XYR.leftCols(Nx1).array().colwise() * w0.array().sqrt()).matrix().adjoint()));
    Yres0 -= XYR.leftCols(Nx1).matrix() * qrXw.solve((XYR.leftCols(Nx1).array().colwise() * w0.array()).matrix().transpose() * Yres0);
    for(i = 0; i < Yres0.cols(); i++) {
      Yres0.col(i) = Mask.select(Yres0.col(i), 0.0);
    }

    // extract quantities from null model
    VectorXd mub0 = yo;
    exp_vector(mub0); 
    mub0.array() /= (1.0 + mub0.array()); 

    // fit full model
    if(offset_mode == "offset") {
      // full model: logit(g) = offset + Y beta
      full = setup_fit(false, true, true); // inc_cov = false, inc_phen = true, use_offset = true
      full.Ncov = Ny; full.Nb = Ny; // overwrite Ncov, Nb
      full.setup_offset_binom(yo, false); // decrement_Nb = false
      full.fit_binom(Mask, Ym, Yres0); // Logistic phenotype residuals

      if(!full.converged) return;
      converged = true;

      /* ll_null = 0.0; */ 
      /* ll_null += Ym.col(0).select((1.0 - mub0.array()).log(), 0.0).array().sum(); // controls */
      /* ll_null += Ym.col(1).select(mub0.array().log(), 0.0).array().sum(); // cases */
      ll_null = null.loglik_multinom(Mask, Ym); // depends on Y, P, Pk, Mask
      if(firth_binom) {
        MatrixXd null_Info = Yres0.transpose() * (Yres0.array().colwise() * w0.array()).matrix();
        LLT<MatrixXd> llt_null(null_Info);
        ll_null += llt_null.matrixL().toDenseMatrix().diagonal().array().log().sum();
      }

      ll_full = full.loglik;

      stat_lrt = 2 * (ll_full - ll_null);
      pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt));
    } else if(offset_mode == "offsetcov") {
      if(!firth_binom) throw std::runtime_error("offsetcov for firth_binom only");

      // null model: logit(g) = [offsetcov; Y] [beta0, betaY] wrt betaY = 0
      MatrixXd Yres0_Int(N, Ny1);
      Yres0_Int.leftCols(1) = Mask.select(yo_int, 0.0);
      Yres0_Int.rightCols(Ny) = Yres0;

      null = setup_fit(false, true, true); // inc_cov = false, inc_phen = true, use_offset = true
      null.Ncov = Ny1; null.Nb = Ny1; // overwrite Ncov, Nb
      null.setup_ncov0(Ny, true, false); // last0 = true, preproc_cov = false
      null.fit_binom(Mask, Ym, Yres0_Int); // Logistic phenotype residuals

      if(trace) { cnt_updates += null.cnt_updates; it += null.it; }

      if(!null.converged) return;

      // full model: logit(g) = [offset; Y] beta
      full = setup_fit(false, true, true); // inc_cov = false, inc_phen = true, use_offset = true
      full.Ncov = Ny1; full.Nb = Ny1; // overwrite Ncov, Nb
      full.fit_binom(Mask, Ym, Yres0_Int); // Logistic phenotype residuals

      if(trace) { cnt_updates += full.cnt_updates; it += full.it; }

      if(!full.converged) return;
      converged = true;

      stat_lrt = 2 * (full.loglik - null.loglik);
      pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt));
    } else if(offset_mode == "offsetcov_int") {
      if(!firth_binom) throw std::runtime_error("offsetcov_int for firth_binom only");

      b0_fit.resize(2);
      b0_fit << null0.bhat(0), 1.0;

      // null model: logit(g) = [1, offsetcov; Y] [beta0, betaY] wrt betaY = 0
      MatrixXd Yres0_Int(N, Ny1 + 1);
      Yres0_Int.leftCols(1) = XYR.leftCols(1);
      Yres0_Int.leftCols(2).rightCols(1) = Mask.select(yo_int, 0.0);
      Yres0_Int.rightCols(Ny) = Yres0;

      null = setup_fit(false, true, true); // inc_cov = false, inc_phen = true, use_offset = true
      null.Ncov = Ny1 + 1; null.Nb = Ny1 + 1; // overwrite Ncov, Nb
      null.setup_ncov0(Ny, true, false); // last0 = true, preproc_cov = false
      null.setup_restart(b0_fit);
      null.fit_binom(Mask, Ym, Yres0_Int); // Logistic phenotype residuals

      if(trace) { cnt_updates += null.cnt_updates; it += null.it; }

      if(!null.converged) return;

      // full model: logit(g) = [1, offset; Y] beta
      full = setup_fit(false, true, true); // inc_cov = false, inc_phen = true, use_offset = true
      full.Ncov = Ny1 + 1; full.Nb = Ny1; // overwrite Ncov, Nb
      null.setup_restart(b0_fit);
      full.fit_binom(Mask, Ym, Yres0_Int); // Logistic phenotype residuals

      if(trace) { cnt_updates += full.cnt_updates; it += full.it; }

      if(!full.converged) return;
      converged = true;

      stat_lrt = 2 * (full.loglik - null.loglik);
      pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt));
    } else if(offset_mode == "offset_int") {
      if(!firth_binom) throw std::runtime_error("offset_int for firth_binom only");

      // null model: logit(g) = offset + [1; Y] [beta0, betaY] wrt betaY = 0
      MatrixXd Yres0_Int(N, Ny1);
      Yres0_Int.leftCols(1) = XYR.leftCols(1);
      Yres0_Int.rightCols(Ny) = Yres0;

      null = setup_fit(false, true, true); // inc_cov = false, inc_phen = true, use_offset = true
      null.Ncov = Ny1; null.Nb = Ny1; // overwrite Ncov, Nb
      null.setup_offset_binom(yo_int, false); // decrement_Nb = false
      null.setup_ncov0(Ny, true, false); // last0 = true, preproc_cov = false
      null.fit_binom(Mask, Ym, Yres0_Int); // Logistic phenotype residuals

      if(trace) { cnt_updates += null.cnt_updates; it += null.it; }

      if(!null.converged) return;

      // full model: logit(g) = offset + [1; Y] beta
      full = setup_fit(false, true, true); // inc_cov = false, inc_phen = true, use_offset = true
      full.Ncov = Ny1; full.Nb = Ny1; // overwrite Ncov, Nb
      full.setup_offset_binom(yo_int, false); // decrement_Nb = false
      full.fit_binom(Mask, Ym, Yres0_Int); // Logistic phenotype residuals

      if(trace) { cnt_updates += full.cnt_updates; it += full.it; }

      if(!full.converged) return;
      converged = true;

      stat_lrt = 2 * (full.loglik - null.loglik);
      pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt));
    } else {
      throw std::runtime_error("unknown offset mode");
    }
  } else if(response == "multinom") {
    executed = true; 

    // fit null model
    if(verbose) cout << "fitting initial null model" << endl;
    null = setup_fit(true, false, false); // inc_cov = true, inc_phen = false, use_offset = false
    null.store_offset = true;
    null.fit_multinom_pom(Mask, Ym, XYR.leftCols(Nx1).rightCols(Nx)); // covariates X without Intercept

    if(trace) { cnt_updates += null.cnt_updates; it += null.it; }

    if(!null.converged) return;
    if(verbose) cout << "initial null converged" << endl;

    // store offset/weights from the null model
    yo = null.yo;
    yo_int = null.yo_int;

    // !NB! not residuals
    MatrixXd Yres0 = XYR.rightCols(Ny21).leftCols(Ny); // Phenotypes

    if(offset_mode == "offset") {
      // full model: logit(gamma) = offset + Y betaY
      full = setup_fit(false, true, true); // inc_cov = false, inc_phen = true, use_offset = true
      full.Ncov = Ny; full.Nb = Ny; // overwrite Ncov, Nb
      full.setup_offset_multinom_pom(yo, yo_int); // manually set up offset
      full.exclude_intercepts = true; full.exclude_intercepts_offset = false;
      full.fit_multinom_pom(Mask, Ym, Yres0);

      if(trace) { cnt_updates += full.cnt_updates; it += full.it; }

      if(!full.converged) return;
      converged = true;

      ll_null = null.loglik_multinom(Mask, Ym); // depends on Y, P, Pk, Mask
      if(firth_multinom) {
        MatrixXd null_Info = MatrixXd(Ny, Ny).setZero().selfadjointView<Lower>().
          rankUpdate((Yres0.array().colwise() * null.WSS1.array()).matrix().adjoint());
        LLT<MatrixXd> llt_null(null_Info);
        ll_null += llt_null.matrixL().toDenseMatrix().diagonal().array().log().sum();
      }
                                                                         
      stat_lrt = 2 * (full.loglik - ll_null);
      pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt));
    } else if(offset_mode == "offset_int") {
      if(!firth_multinom) throw std::runtime_error("offset_int for firth_multinom only");

      b0_fit.resize(2);
      b0_fit << yo_int;

      // null model: logit(gamma) = offset + Y betaY wrt betaY = 0
      null = setup_fit(false, true, true); // inc_cov = false, inc_phen = true, use_offset = true
      null.Ncov = Ny; null.Nb = Ny + ncat1; // overwrite Ncov, Nb
      null.setup_offset_multinom_pom(yo, yo_int); // manually set up offset
      null.exclude_intercepts = false; null.exclude_intercepts_offset = true;
      null.setup_ncov0(Ny, true, false); // last0 = true, preproc_cov = false
      null.setup_restart(b0_fit);
      null.fit_multinom_pom(Mask, Ym, Yres0);

      if(trace) { cnt_updates += null.cnt_updates; it += null.it; }

      if(!null.converged) return;
      if(verbose) cout << "null converged" << endl;

      // full model: logit(gamma) = offset + Y betaY
      full = setup_fit(false, true, true); // inc_cov = false, inc_phen = true, use_offset = true
      full.Ncov = Ny; full.Nb = Ny + ncat1; // overwrite Ncov, Nb
      full.setup_offset_multinom_pom(yo, yo_int); // manually set up offset
      full.exclude_intercepts = false; full.exclude_intercepts_offset = true;
      full.setup_restart(b0_fit);
      full.fit_multinom_pom(Mask, Ym, Yres0);

      if(trace) { cnt_updates += full.cnt_updates; it += full.it; }

      if(!full.converged) return;
      converged = true;
      if(verbose) cout << "full converged" << endl;

      stat_lrt = 2 * (full.loglik - null.loglik);
      pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt));
      if(verbose) cout << "pval_test =  " << pval_test << endl;
    } else {
      throw std::runtime_error("unknown offset mode");

      /* // residualize phenotypes */
      /* // !NB! not implemented yet */

      /* // full model */
      /* full = setup_fit(false, true, true); // inc_cov = false, inc_phen = true, use_offset = true */
      /* full.setup_offset_multinom_pom(yo, yo_int); // manually set up offset */
      /* full.exclude_intercepts = true; */
      /* full.Ncov = Ny; full.Nb = Ny; // overwrite Ncov, Nb */
      /* full.fit_multinom_pom(Mask, Ym, XYR.rightCols(Ny21).leftCols(Ny)); // Phenotypes */
      /* /1* if(offset_mode == "offset") { *1/ */
      /* /1*   full = setup_fit(false, true, true); // inc_cov = false, inc_phen = true, use_offset = true *1/ */
      /* /1*   full.setup_offset_multinom_pom(yo, yo_int); // manually set up offset *1/ */
      /* /1*   full.exclude_intercepts = true; *1/ */
      /* /1*   full.Ncov = Ny; full.Nb = Ny; // overwrite Ncov, Nb *1/ */
      /* /1*   full.fit_multinom_pom(Mask, Ym, XYR.rightCols(Ny21).leftCols(Ny)); // Phenotypes *1/ */
      /* /1* } else if(offset_mode == "offset_int") { *1/ */
      /* /1*   full = setup_fit(false, true, true); // inc_cov = false, inc_phen = true, use_offset = true *1/ */
      /* /1*   full.setup_offset_multinom_pom(yo, yo_int); // manually set up offset *1/ */
      /* /1*   full.exclude_intercepts = false; *1/ */
      /* /1*   full.Ncov = Ny; full.Nb = ncat1 + Ny; // overwrite Ncov, Nb *1/ */
      /* /1*   full.fit_multinom_pom(Mask, Ym, XYR.rightCols(Ny21).leftCols(Ny)); // Phenotypes *1/ */
      /* /1* } else { *1/ */
      /* /1*   throw std::runtime_error("unknown offset mode"); *1/ */
      /* /1* } *1/ */

      /* if(trace) { cnt_updates += full.cnt_updates; it += full.it; } */

      /* if(!full.converged) return; */
      /* converged = true; */

      /* stat_lrt = 2 * (full.loglik - null.loglik); */
      /* pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt)); */
    }
  } else {
    throw std::runtime_error("unknown response");
  }

  // store results
  if(converged) {
    bhat_y = full.bhat.tail(Ny);
  }
}

void MultiPhen::run_test_qt(const Eigen::Ref<const Eigen::MatrixXd> & XYR)
{
  reset_model(); // reset model fit results

  if(response == "binom") {
    executed = true; 
    converged = true; 
    VectorXd beta_qt = XYR.leftCols(Nx1).transpose() * yb;
    // residualize
    VectorXd y_qt = yb - XYR.leftCols(Nx1) * beta_qt;
    VectorXd x_qt = XYR.leftCols(Nx1 + 1).rightCols(1);
    // regression
    /* VectorXd bhat_qt = (y_qt.transpose() * x_qt) / x2; */
    /* bhat = (Y.col(i).transpose() * G).array().rowwise() / G2.array().transpose(); */
    /* /1* B.row(i) = bhat; *1/ */
    /* // residuals, s2 */
    /* s2 = (((G.array().rowwise() * bhat.array().transpose()). // predicted yp = X bhat */
    /*   colwise() - Y.col(i).array()). // residuals = y - yp */
    /*   matrix().colwise().squaredNorm()). // residuals^2 */
    /*   array() / (N_data - 1.0); // s2 = residuals^2 / (N - 1) */
    /* Z.row(i) = bhat.array() * (G2.array() / s2.array()).sqrt(); */
    
    /* // regression */
    /* bhat = (Y.col(i).transpose() * G).array().rowwise() / G2.array().transpose(); */
    /* /1* B.row(i) = bhat; *1/ */
    /* // residuals, s2 */
    /* s2 = (((G.array().rowwise() * bhat.array().transpose()). // predicted yp = X bhat */
    /*   colwise() - Y.col(i).array()). // residuals = y - yp */
    /*   matrix().colwise().squaredNorm()). // residuals^2 */
    /*   array() / (N_data - 1.0); // s2 = residuals^2 / (N - 1) */
    /* Z.row(i) = bhat.array() * (G2.array() / s2.array()).sqrt(); */


    /* yb */ 
  /* pval_test = test_score(null, Mask, Ym, yb, XYR, inc_cov); */ 
  } else {
    return;
  }

}

void MultiPhen::run_test_lrt(const Eigen::Ref<const Eigen::MatrixXd> & XYR, bool inc_cov)
{
  reset_model(); // reset MultiPhen model fit results
  executed = true; 
  
  FitOrdinal null, full;

  if(reuse_start & !approx_offset) {
    if(!inc_cov) throw std::runtime_error("reuse_start in not available for inc_cov = false");
    /* if(approx_offset) throw std::runtime_error("reuse_start is not compatible with approx_offset"); */

    // null model: logit(g) = X alpha 
    null = fit(XYR, inc_cov, false); // inc_cov, inc_phen = false
    if(!null.converged) return;

    b0 = null.bhat;

    // full model: logit(g) = X alpha + Y beta
    full = fit(XYR, inc_cov, true); // inc_cov, inc_phen = true
    // give another chance if reuse_start & reset_start
    if(reset_start) {
      reuse_start = false;
      full = fit(XYR, inc_cov, true); // inc_cov, inc_phen = true
    }
    if(!full.converged) return;

    converged = true;
    boost::math::chi_squared dist(Ny);
    double stat_lrt = 2 * (full.loglik - null.loglik);
    pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt));
  } else if(approx_offset && response == "binom") {
    // null model: logit(g) = X alpha 
    null = fit(XYR, inc_cov, false); // inc_cov, inc_phen = false
    if(!null.converged) return;

    // store offset/weights from the null mode 
    yo = null.yo;
    w0 = null.wb;

    Yres0 = XYR.rightCols(Ny21).leftCols(Ny); // matrix of phenotypes Y 
    ColPivHouseholderQR<MatrixXd> qrXw;
    qrXw.compute(MatrixXd(Nx1, Nx1).setZero().selfadjointView<Lower>().rankUpdate((XYR.leftCols(Nx1).array().colwise() * w0.array().sqrt()).matrix().adjoint()));
    Yres0 -= XYR.leftCols(Nx1).matrix() * qrXw.solve((XYR.leftCols(Nx1).array().colwise() * w0.array()).matrix().transpose() * Yres0);

    // full model: logit(g) = X alpha + Y beta
    full = fit(XYR, inc_cov, true); // inc_cov, inc_phen = true
    if(!full.converged) return;
    converged = true;

    // problem: null.mub is not at scale [0, 1]
    /* cout << "null.mub = " << null.mub.head(5).transpose() << endl; */
    VectorXd mub = null.yo;
    exp_vector(mub); // mub <- exp(mub)
    mub.array() /= (1.0 + mub.array()); // mub <- exp(mub) / (1 + exp(mub))
                                        //
    double ll_null = null.loglik_binom(Mask, Ym);
    /* double ll_null = 0.0; */ 
    /* ll_null += Ym.col(0).select((1.0 - mub.array()).log(), 0.0).array().sum(); // controls */
    /* ll_null += Ym.col(1).select(mub.array().log(), 0.0).array().sum(); // cases */
    if(firth_binom) {
      MatrixXd null_Info = Yres0.transpose() * (Yres0.array().colwise() * w0.array()).matrix();
      LLT<MatrixXd> llt_null(null_Info);
      ll_null += llt_null.matrixL().toDenseMatrix().diagonal().array().log().sum();
    }

    double ll_full;
    if(full.firth_binom) {
      LLT<MatrixXd> llt_full(full.Info);
      ll_full = full.loglik_binom_firth(Mask, Ym, llt_full);
    } else {
      ll_full = full.loglik_binom(Mask, Ym);
    }

    boost::math::chi_squared dist(Ny);
    /* double stat_lrt = 2 * (full.loglik - null.loglik); */
    double stat_lrt = 2 * (ll_full - ll_null);
    pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt));
  } else if(approx_offset && response == "multinom") {
    // null model
    null = fit(XYR, inc_cov, false); // inc_cov, inc_phen = false
    if(!null.converged) return;

    // store offset vectors
    yo = null.yo;
    yo_int = null.yo_int;

    // full model: logit(g) = X alpha + Y beta
    full = fit(XYR, inc_cov, true); // inc_cov, inc_phen = true
    if(!full.converged) return;
    converged = true;

    boost::math::chi_squared dist(Ny);
    double stat_lrt = 2 * (full.loglik - null.loglik);
    /* cout << "stat_lrt = " << stat_lrt << " full.loglik = " << full.loglik << " null.loglik = " << null.loglik << endl; */
    pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt));
  } else {
    // null model: logit(g) = X alpha 
    null = fit(XYR, inc_cov, false); // inc_cov, inc_phen = false
    if(!null.converged) return;

    // full model: logit(g) = X alpha + Y beta
    full = fit(XYR, inc_cov, true); // inc_cov, inc_phen = true
    if(!full.converged) return;

    converged = true;
    boost::math::chi_squared dist(Ny);
    double stat_lrt = 2 * (full.loglik - null.loglik);
    /* cout << " lrt = " << stat_lrt << " = " << full.loglik << " - " << null.loglik << endl; */
    pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt));
  }
  // store results
  if(converged) {
    bhat_y = full.bhat.tail(Ny);
  }
}

void MultiPhen::run_test_score(const Eigen::Ref<const Eigen::MatrixXd> & XYR, bool inc_cov)
{
  bool _firth_binom = firth_binom, _firth_multinom = firth_multinom, _approx_offset = approx_offset;
  firth_binom = false; firth_multinom = false; approx_offset = false;
 
  reset_model(); // reset model fit results
  executed = true; 

  FitOrdinal null = fit(XYR, inc_cov, false); // inc_cov, inc_phen = false
  if(!null.converged) { return; }

  converged = true; 
  if(trace) { cnt_updates += null.cnt_updates; it += null.it; }
  pval_test = test_score(null, Mask, Ym, yb, XYR, inc_cov); 

  firth_binom = _firth_binom; firth_multinom = _firth_multinom; approx_offset = _approx_offset;
}

void MultiPhen::setup_x(const VectorXb & _Mask,  const Eigen::MatrixXd& XYR, unsigned int n_cov, unsigned int n_phen, 
    bool _pos_intercept_first, bool _pos_phen_first)
{
  // check
  if(XYR.cols() != 2 + n_cov + 2*n_phen) throw std::runtime_error("setup_x: dimensions XYR");
  if(XYR.rows() != _Mask.size()) throw std::runtime_error("setup_x: dimensions XYR and Mask");
  // extract dimensions from XYR
  N = XYR.rows();
  /* Ncov = n_cov; // Nb = ncat1 + Ncov, where ncat1 depend on g */
  Nx = n_cov; Nx1 = n_cov + 1; Ny = n_phen; Ny1 = n_phen + 1; Ny21 = Ny1 + n_phen;
  pos_intercept_first = _pos_intercept_first;
  pos_phen_first = _pos_phen_first;
  // Mask
  Mask = _Mask; // VectorXb::Constant(N, true);
  Neff = Mask.array().cast<double>().sum();
  // update status
  set_x = true;
}

void MultiPhen::reset_model()
{
  executed = false; converged = false;
  pval_test = -1.0;
  it = 0; cnt_updates = 0;
}

void MultiPhen::setup_approx_offset()
{
  if(!set_y) throw std::runtime_error("setup_approx_offset: set_y is false");

  if(mac_approx_offset == 0) {
    approx_offset = false;
  } else if(mac_approx_offset == 1) {
    approx_offset = true;
  } else if(mac_approx_offset > 1) {
    if(Ncat_minor <= mac_approx_offset) approx_offset = false;
    else approx_offset = true;
  }
}

void MultiPhen::setup_y(const Eigen::VectorXd & _g)
{
  // Eigen::VectorXi g = _g.cast<int>(); // 1.6 -> 1
  Eigen::VectorXi g = _g.array().round().cast<int>(); // 1.6 -> 2

  unsigned int i;
  std::set<int> genotypes; // ordered (!) set of category levels
  set<int>::iterator it_set;

  // checks
  if(N == 0) throw std::runtime_error("setup_y: N == 0");
  if(g.size() != N) throw std::runtime_error("setup_y: g.size() != N");

  // assign category levels 
  for(i = 0; i < g.size(); i++) if(Mask(i)) genotypes.insert(g[i]);

  // check genotypes levels: 0/1 or 0/1/2
  /* for(i = 0, it_set = genotypes.begin(); i < genotypes.size(); i++, it_set++) cout << "genotypes " << i << " = " << *it_set << endl; */
  /* cout << "genotypes.size() = " << genotypes.size() << endl; */
  if(genotypes.size() == 1) {
    /* cerr << "WARNING: number of genotype categories is 1" << endl; */
    return;
  }
  if(!(genotypes.size() == 2 || genotypes.size() == 3)) throw std::runtime_error("setup_y: number of genotype categories must be 2 or 3");

  // assign ncat, ncat1
  ncat = genotypes.size();
  ncat1 = ncat - 1; ncat1sq = ncat1 * ncat1;
  
  // assign response
  if(ncat == 2) response = "binom";
  else if(ncat == 3) response = "multinom";
  else throw std::runtime_error("setup_y: unexpected number of genotype categories");

  // assign Ncov, Nb
  /* Nb = ncat1 + Ncov; */

  // assign Ym
  Ym.resize(N, ncat);
  Ncat = VectorXi::Constant(ncat, 0);
  Ncat_minor = 0;
  int Ncat_max = 0;
  // loop over a a few genotype categories
  for(i = 0, it_set = genotypes.begin(); i < ncat; i++, it_set++) {
    Ym.col(i) = Mask.select(g.array() == *it_set, false);
    /* Ym.col(i) = (g.array() == *it_set); */
    /* Ym.col(i) = Mask.select(Ym.col(i), false); */
    Ncat(i) = Ym.col(i).cast<int>().sum();
    // get the maximum value in Ncat & minor counts in Ncat (all except the maximum)
    if(Ncat(i) > Ncat_max) Ncat_max = Ncat(i);
    Ncat_minor += Ncat(i);
  }
  Ncat_minor -= Ncat_max;
  // assign yb if binomial
  if(response == "binom") {
    yb = Ym.col(1).cast<double>(); // booleans -> 0/1
  }
  // update status
  set_y = true;
}

void MultiPhen::test0(const Eigen::VectorXi & g, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
      bool firth_binom,
      std::string optim, double tol, unsigned int maxit, bool check_step, double max_step)
{
  executed = true;
  converged = false;
  pval_test = -1.0;

  FitOrdinal fit, fit1;
  
  // set up Ordinal model (no Firth)
  Ordinal ord;
  ord.optim = optim; ord.tol = tol; ord.pseudo_stophalf = pseudo_stophalf; ord.maxit = maxit;
  ord.check_step = check_step; ord.max_step = max_step;
  ord.firth_binom = false;
  
  // fit null model
  fit = ord.fit(g, X);
  if(!fit.converged) { return; }

  // run Score test
  converged = fit.converged;
  pval_test = ord.test_score(fit, Y);

  // run LRT test (if needed)
  if(pval_test < pval_thr) {
    pval_test = -1.0;
    converged = false;

    // prepare new matrix of covariates X + Y
    MatrixXd X1(Y.rows(), X.cols() + Y.cols());
    if(X.cols()) {
      X1.leftCols(X.cols()) = X;
    }
    X1.rightCols(Y.cols()) = Y;

    if(firth_binom & (ord.response == "binom")) {
      ord.firth_binom = firth_binom;

      // fit null model (Firth) for LRT
      fit = ord.fit(g, X1, Y.cols());
      if(!fit.converged) { return; }

      // fit alternative model (Firth)
      fit1 = ord.fit(g, X1);
      if(!fit1.converged) { return; }
      converged = fit1.converged;

      boost::math::chi_squared dist(Y.cols());
      double stat_lrt = 2 * (fit1.loglik - fit.loglik);
      pval_test = boost::math::cdf(boost::math::complement(dist, stat_lrt));
    } else {
      // fit alternative model (no Firth)
      fit1 = ord.fit(g, X1);
      if(!fit1.converged) { return; }
      converged = fit1.converged;

      boost::math::chi_squared dist(Y.cols());
      double stat_lrt = 2 * (fit1.loglik - fit.loglik);
      pval_test = boost::math::cdf(boost::math::complement(dist, stat_lrt));
    }
  }
}

void MultiPhen::test_addcov(const Eigen::VectorXi & g, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
      bool firth_binom,
      std::string optim, double tol, unsigned int maxit, bool check_step, double max_step)
{
  executed = true;
  converged = false;
  pval_test = -1.0;

  FitOrdinal fit, fit1;
  MatrixXd X0;
  
  // set up Ordinal model (no Firth)
  Ordinal ord;
  ord.optim = optim; ord.tol = tol; ord.pseudo_stophalf = pseudo_stophalf; ord.maxit = maxit;
  ord.check_step = check_step; ord.max_step = max_step;
  ord.firth_binom = false;
  
  // fit null model
  fit = ord.fit(g, X0);
  if(!fit.converged) { return; }

  // run Score test
  converged = fit.converged;
  pval_test = ord.test_score(fit, Y);

  // run LRT test (if needed)
  if(pval_test < pval_thr) {
    pval_test = -1.0;
    converged = false;

    // null model (with covariates)
    fit = ord.fit(g, X);
    if(!fit.converged) { return; }

    // prepare new matrix of covariates X + Y
    MatrixXd X1(Y.rows(), X.cols() + Y.cols());
    if(X.cols()) {
      X1.leftCols(X.cols()) = X;
    }
    X1.rightCols(Y.cols()) = Y;

    if(firth_binom & (ord.response == "binom")) {
      ord.firth_binom = firth_binom;

      // re-fit null model (Firth)
      fit = ord.fit(g, X);
      if(!fit.converged) { return; }

      // fit alternative model (Firth)
      fit1 = ord.fit(g, X1);
      if(!fit1.converged) { return; }
      converged = fit1.converged;

      boost::math::chi_squared dist(Y.cols());
      double stat_lrt = 2 * (fit1.loglik - fit.loglik);
      pval_test = boost::math::cdf(boost::math::complement(dist, stat_lrt));
    } else {
      // fit alternative model (no Firth)
      fit1 = ord.fit(g, X1);
      if(!fit1.converged) { return; }
      converged = fit1.converged;

      boost::math::chi_squared dist(Y.cols());
      double stat_lrt = 2 * (fit1.loglik - fit.loglik);
      pval_test = boost::math::cdf(boost::math::complement(dist, stat_lrt));
    }
  }
}

//------------------------
// Class FitOrdinal
//------------------------

void FitOrdinal::setup_defaults()
{ 
  cnt_fit = 0;
  verbose = 0;
  // model parameters
  N = 0; Neff = 0; // sample size
  Ncov = 0, Nb = 0, Ncov0 = 0; Ncov1 = 0; // number of covariates 
  ncat = 0, ncat1 = 0, ncat1sq = 0; // number of categories 
  
  firth_binom = false;
  firth_mult = 1.0;

  apply_start = false;
  store_offset = false;
  apply_offset = false;
  exclude_intercepts = false;

  // model fitting results
  executed = false; converged = false;
  trace = false;
  it = 0; maxit = 0; cnt_updates = 0;

}

FitOrdinal::FitOrdinal()
{ 
  setup_defaults();
}

//-----------------------------
//  Class FitOrdinal: Checkers
//-----------------------------

void FitOrdinal::check_setup_model()
{
  if(verbose >= 2) {
    cout << "check_setup_model" << endl;
    cout << " --  N = " << N << " Neff = " << Neff << " Nb = " << Nb << " Ncov = " << Ncov 
      << " Ncov0 = " << Ncov0 << " Ncov1 = " << Ncov1
      << " apply_start = " << apply_start << " apply_offset = " << apply_offset << " store_offset = " << store_offset << " exclude_intercepts = " << exclude_intercepts << "exclude_intercepts_offset = " << exclude_intercepts_offset 
      << " firth_multinom = " << firth_multinom << " firth_binom = " << firth_binom << " firth_mult = " << firth_mult << " check_step = " << check_step << " max_step = " << max_step 
      << " maxit = " << maxit << " maxit2 = " << maxit2 << " maxit3 = " << maxit3
      << endl;
    cout << " Ncat = " << Ncat << endl;
  }

  check_setup_model_common();
}

void FitOrdinal::check_setup_model_common()
{
  if(N == 0) { throw std::runtime_error("check_setup_model: N == 0"); }
  if(Neff == 0) { throw std::runtime_error("check_setup_model: Neff == 0"); }
  if(Nb == 0) { throw std::runtime_error("check_setup_model: Nb == 0"); }
  if(ncat == 0) { throw std::runtime_error("check_setup_model: ncat == 0"); }
}

void FitOrdinal::check_setup_data()
{
  if(verbose >= 2) cout << "check_setup_data\n";

  check_setup_data_common();

  if(response == "multinom") {
    check_setup_data_multinom();
  } else if(response == "binom") {
    check_setup_data_binom();
  } else {
    throw std::runtime_error("unknown response");
  }
}

void FitOrdinal::check_setup_data_common()
{
  /* cur_Score.resize(Nb); */
  /* cur_Info.resize(Nb, Nb); */
  /* cur_v.resize(Nb); cur_b.resize(Nb); */
}

void FitOrdinal::check_setup_data_multinom()
{
}

void FitOrdinal::check_setup_data_binom()
{
  if(mub.size() != N) { throw std::runtime_error("check_setup_model: mub.size() != N"); }
  if(wb.size() != N) { throw std::runtime_error("check_setup_model: wb.size() != N"); }
  /* XtW.resize(Nb, N); */
}

//-----------------------------
//  Class Ordinal: Constructors
//-----------------------------

Ordinal::Ordinal() 
{ 
  setup_defaults();
}

void Ordinal::setup_defaults() 
{
  response = "multinom";
  optim = "WeightHalving";
  firth_binom = false; firth_multinom = false;

  maxit = 100; maxit2 = 7; maxit3 = 25;
  it2 = 0; strict = false;
  tol = 1e-4; pseudo_stophalf = 0.0;

  check_step = false;
  max_step = 10.0;

  preproc_cov = false;

  cur_dev = 0; prev_dev = 0;
}

// constructor for FitOrdinal
// - copy model parameters to FitOrdinal object
FitOrdinal Ordinal::setup_fit()
{
  FitOrdinal fit;

  // copy parameters from Ordinal
  fit.response = response; // response type = [binom, multinom]
  fit.model = model; // model = [POM: Proportional Odds Model, ACL: Adjacent Category Logit]
  fit.optim = optim; // optimization algorithm = [FisherScoring, WeightHalving]
  fit.firth_binom = firth_binom; fit.firth_multinom = firth_multinom; // Firth correction
      
  fit.maxit = maxit; fit.maxit2 = maxit2; fit.maxit3 = maxit3; fit.strict = strict;
  fit.tol = tol; fit.pseudo_stophalf = pseudo_stophalf;

  fit.check_step = check_step;
  fit.max_step = max_step;

  fit.N = N; fit.Neff = Neff; // samples size
  fit.Ncov = Ncov; fit.Nb = Nb; // number of covariates
                    
  fit.ncat = ncat; // number of categories
  fit.ncat1 = ncat1; fit.ncat1sq = ncat1sq;
  fit.Ncat = Ncat;

  fit.cur_dev = 0; fit.prev_dev = 0;
  fit.it = 0; fit.it2 = 0;

  return(fit);
}

//--------------------------
//  MultiPhen: Score Test
//--------------------------

double MultiPhen::test_score(const FitOrdinal & null, 
    const VectorXb & Mask, const MatrixXb & Ym, const Eigen::VectorXd & yb, 
    const Eigen::Ref<const Eigen::MatrixXd> & XYR, bool inc_cov)
{
  double pval;
  if(response == "multinom") {
    if(inc_cov) pval = test_score_multinom_pom(null, Mask, Ym, XYR.leftCols(Nx1).rightCols(Nx), XYR.rightCols(Ny21).leftCols(Ny)); // covariates X (no intercept); phenotypes Y
    else pval = test_score_multinom_pom(null, Mask, Ym, XYR.leftCols(0), XYR.rightCols(Ny21).leftCols(Ny)); // 0 covarites (no intercept); phenotypes Y
  } else if(response == "binom") {
    if(inc_cov) pval = test_score_binom(null, Mask, yb, XYR.leftCols(Nx1).rightCols(Nx), XYR.rightCols(Ny21).leftCols(Ny)); // covariates X (no intercept); phenotypes Y
    else pval = test_score_binom(null, Mask, yb, XYR.leftCols(0), XYR.rightCols(Ny21).leftCols(Ny)); // 0 covarites (no intercept); phenotypes Y
  } else {
    throw std::runtime_error("unknown response");
  }
  return(pval);
}

double MultiPhen::test_score_multinom_pom(const FitOrdinal & null, 
    const VectorXb & Mask, const MatrixXb & Ym, 
    const Eigen::Ref<const Eigen::MatrixXd> & X, const Eigen::Ref<const Eigen::MatrixXd> & G)
{
  unsigned int k;
  unsigned int Ng = G.cols(); // l = ncol(G), p = Nb, m = ncat1;

  // check dimensions
  if(G.cols() == 0) throw std::runtime_error("#cols in G is 0");
  if(G.rows() == 0) throw std::runtime_error("#rows in G is 0");
  if(G.rows() != N) throw std::runtime_error("#rows in G is different from N");
  if(Ym.rows() != N) throw std::runtime_error("#rows in Ym is different from N");
  if(X.cols() != null.Ncov) throw std::runtime_error("#cols in X != null.Ncov (test_score_multinom_pom)");

  // Score vector with l elements
  VectorXd Score1 = ((null.V).transpose() * G).colwise().sum();
  
  // pre-compute
  MatrixXd GW = G.transpose() * null.W;

  // Info matrix
  // 1x1 block of size pxp
  // P = null.Score
  // 1x2 block = W = p x m matrix
  MatrixXd Info1_W(null.Nb, Ng);
  // fill in part 1 of Info1_W: first ncat1 rows
  for(k = 0; k < ncat1; k++) {
    MatrixXd GWs = GW(all, seqN(k, ncat1, ncat1));
    VectorXd GW1 = GWs.rowwise().sum();
    Info1_W.row(k) = GW1.array();
  }
  // fill in part 2 of Info1_W: last Ncov rows
  if(null.Ncov) {
    MatrixXd GW12 = (X.array().colwise() * (null.WSS1).array()).matrix().transpose() * 
      (G.array().colwise() * (null.WSS1).array()).matrix();
    Info1_W(seqN(ncat1, null.Ncov), all) = GW12;
  }

  // Info1_Q
  MatrixXd Info1_Q = MatrixXd(Ng, Ng).setZero().selfadjointView<Lower>().
      rankUpdate((G.array().colwise() * (null.WSS1).array()).matrix().adjoint());

  // Variance matrix V of the scire Score1
  // V = (Q - W' Info0^{-1} W
  LLT<MatrixXd> llt_Info0(null.Info);
  MatrixXd Var_Score1 = (Info1_Q - (Info1_W.transpose() * llt_Info0.solve(Info1_W)));

  // Test statistic = Score1' Var_Score1^{-1} Score1
  LLT<MatrixXd> llt_Var(Var_Score1);
  double stat_score = Score1.transpose() * llt_Var.solve(Score1);
  
  // R: pchisq(stat_score, Ng, lower = FALSE)
  boost::math::chi_squared dist(Ng);
  double pval = boost::math::cdf(boost::math::complement(dist, stat_score));

  return(pval);
}

double MultiPhen::test_score_binom(const FitOrdinal & null, 
    const VectorXb & Mask, const Eigen::VectorXd & yb, 
    const Eigen::Ref<const Eigen::MatrixXd> & X, const Eigen::Ref<const Eigen::MatrixXd> & G)
{
  unsigned int Ng = G.cols(); // l = ncol(G), p = Nb, m = ncat1;

  // check dimensions
  if(G.cols() == 0) throw std::runtime_error("#cols in G is 0");
  if(G.rows() == 0) throw std::runtime_error("#rows in G is 0");
  if(G.rows() != N) throw std::runtime_error("#rows in G is different from N");
  if(yb.size() != N) throw std::runtime_error("size of yb is different from N");
  if(X.cols() != null.Ncov) throw std::runtime_error("#cols in X != null.Ncov (test_score_binom)");

  // Score vector with Ng elements
  VectorXd Score1 = G.transpose() * (yb - null.mub);
  
  // Info matrix 
  // Info1_W
  MatrixXd Info1_W(null.Nb, Ng);
  Info1_W.row(0) = (G.array().colwise() * (null.wb).array()).colwise().sum();
  if(null.Ncov) {
    Info1_W(seqN(1, null.Ncov), all) = X.transpose() * (G.array().colwise() * (null.wb).array()).matrix();
  }

  // Info1_Q
  MatrixXd Info1_Q = MatrixXd(Ng, Ng).setZero().selfadjointView<Lower>().
      rankUpdate((G.array().colwise() * (null.wb).array().sqrt()).matrix().adjoint());

  // Variance matrix V of the scire Score1
  // V = (Q - W' Info0^{-1} W
  LLT<MatrixXd> llt_Info0(null.Info);
  MatrixXd Var_Score1 = (Info1_Q - (Info1_W.transpose() * llt_Info0.solve(Info1_W)));

  // Test statistic = Score1' Var_Score1^{-1} Score1
  LLT<MatrixXd> llt_Var(Var_Score1);
  double stat_score = Score1.transpose() * llt_Var.solve(Score1);
  
  // R: pchisq(stat_score, Ng, lower = FALSE)
  boost::math::chi_squared dist(Ng);
  double pval = boost::math::cdf(boost::math::complement(dist, stat_score));

  return(pval);
}

//------------------
//  Score Test
//------------------

double Ordinal::test_score(const FitOrdinal & null, const Eigen::MatrixXd & G)
{
  double pval;
  if(response == "multinom") {
    pval = test_score_multinom_pom(null, Xcov, G);
  } else if(response == "binom") {
    pval = test_score_binom(null, Xcov, G);
  } else {
    throw std::runtime_error("unknown response");
  }
  return(pval);
}

double Ordinal::test_score(const FitOrdinal & null, const Eigen::MatrixXd & X, const Eigen::MatrixXd & G)
{
  double pval;
  if(response == "multinom") {
    pval = test_score_multinom_pom(null, X, G);
  } else if(response == "binom") {
    pval = test_score_binom(null, X, G);
  } else {
    throw std::runtime_error("unknown response");
  }

  return(pval);
}

double Ordinal::test_score_binom(const FitOrdinal & null, const Eigen::MatrixXd & X, const Eigen::MatrixXd & G)
{
  unsigned int Ng = G.cols(); // l = ncol(G), p = Nb, m = ncat1;

  // check dimensions
  if(G.cols() == 0) throw std::runtime_error("#cols in G is 0");
  if(G.rows() == 0) throw std::runtime_error("#rows in G is 0");
  if(G.rows() != N) throw std::runtime_error("#rows in G is different from N");

  // Score vector with l elements
  VectorXd Score1 = G.transpose() * (yb - null.mub);
  
  // Info matrix 
  // Info1_W
  MatrixXd Info1_W(Nb, Ng);
  Info1_W.row(0) = (G.array().colwise() * (null.wb).array()).colwise().sum();
  if(Ncov) {
    Info1_W(seqN(1, Ncov), all) = X.transpose() * (G.array().colwise() * (null.wb).array()).matrix();
  }
  // Info1_Q
  MatrixXd Info1_Q = MatrixXd(Ng, Ng).setZero().selfadjointView<Lower>().
      rankUpdate((G.array().colwise() * (null.wb).array().sqrt()).matrix().adjoint());

  // Variance matrix V of the scire Score1
  // V = (Q - W' Info0^{-1} W
  LLT<MatrixXd> llt_Info0(null.Info);
  MatrixXd Var_Score1 = (Info1_Q - (Info1_W.transpose() * llt_Info0.solve(Info1_W)));

  // Test statistic = Score1' Var_Score1^{-1} Score1
  LLT<MatrixXd> llt_Var(Var_Score1);
  double stat_score = Score1.transpose() * llt_Var.solve(Score1);
  
  // R: pchisq(stat_score, Ng, lower = FALSE)
  boost::math::chi_squared dist(Ng);
  double pval = boost::math::cdf(boost::math::complement(dist, stat_score));

  return(pval);
}

double Ordinal::test_score_multinom_pom(const FitOrdinal & null, const Eigen::MatrixXd & X, const Eigen::MatrixXd & G)
{
  unsigned int k;
  unsigned int Ng = G.cols(); // l = ncol(G), p = Nb, m = ncat1;

  // check dimensions
  if(G.cols() == 0) throw std::runtime_error("#cols in G is 0");
  if(G.rows() == 0) throw std::runtime_error("#rows in G is 0");
  if(G.rows() != N) throw std::runtime_error("#rows in G is different from N");

  // Score vector with l elements
  VectorXd Score1 = ((null.V).transpose() * G).colwise().sum();
  
  // pre-compute
  MatrixXd GW = G.transpose() * null.W;

  // Info matrix
  // 1x1 block of size pxp
  // P = null.Score
  // 1x2 block = W = p x m matrix
  MatrixXd Info1_W(Nb, Ng);
  // fill in part 1 of Info1_W: first ncat1 rows
  for(k = 0; k < ncat1; k++) {
    MatrixXd GWs = GW(all, seqN(k, ncat1, ncat1));
    VectorXd GW1 = GWs.rowwise().sum();
    Info1_W.row(k) = GW1.array();
  }
  // fill in part 2 of Info1_W: last Ncov rows
  if(Ncov) {
    MatrixXd GW12 = (X.array().colwise() * (null.WSS1).array()).matrix().transpose() * 
      (G.array().colwise() * (null.WSS1).array()).matrix();
    Info1_W(seqN(ncat1, Ncov), all) = GW12;
  }

  // Info1_Q
  MatrixXd Info1_Q = MatrixXd(Ng, Ng).setZero().selfadjointView<Lower>().
      rankUpdate((G.array().colwise() * (null.WSS1).array()).matrix().adjoint());

  // Variance matrix V of the scire Score1
  // V = (Q - W' Info0^{-1} W
  LLT<MatrixXd> llt_Info0(null.Info);
  MatrixXd Var_Score1 = (Info1_Q - (Info1_W.transpose() * llt_Info0.solve(Info1_W)));

  // Test statistic = Score1' Var_Score1^{-1} Score1
  LLT<MatrixXd> llt_Var(Var_Score1);
  double stat_score = Score1.transpose() * llt_Var.solve(Score1);
  
  // R: pchisq(stat_score, Ng, lower = FALSE)
  boost::math::chi_squared dist(Ng);
  double pval = boost::math::cdf(boost::math::complement(dist, stat_score));

  return(pval);
}

//------------------
//  Set up
//------------------

void Ordinal::setup_xy(const Eigen::VectorXi &y, const Eigen::MatrixXd& X)
{
  unsigned int i;
  set<int>::iterator it_set;

  // assign N
  N = y.size();
  // assign category levels 
  for(i = 0; i < y.size(); i++) {
    cat.insert(y[i]);
  }
  // assign ncat, ncat1
  ncat = cat.size();
  ncat1 = ncat - 1;
  ncat1sq = ncat1 * ncat1;
  
  // check if the type response 
  if(ncat == 2) {
    response = "binom";
  }
  
  // assign Ncov, Nb
  Ncov = X.cols();
  Nb = ncat1 + Ncov;
  // assign Mask
  Mask = VectorXb::Constant(N, true);
  // assign Neff
  Neff = Mask.array().cast<double>().sum();
  // process X
  if(Ncov && preproc_cov) {
    Xcov = orth_matrix(X, Mask);
  } else {
    Xcov = X;
  }
  // update Ncov: some colinear columns in X might be removed
  Ncov = Xcov.cols();

  if(Ncov) {
    Xcov1.resize(Xcov.rows(), Xcov.cols() + 1);
    Xcov1.col(0).array() = VectorXd::Ones(Xcov1.rows());
    Xcov1.rightCols(Ncov) = Xcov;
  } else {
    Xcov1 = MatrixXd::Ones(N, 1);
  }

  // assign Y
  Y.resize(N, ncat);
  Ncat = VectorXi::Constant(ncat, 0);
  for(i = 0, it_set = cat.begin(); i < ncat; i++, it_set++) {
    Y.col(i) = Mask.select(y.array() == *it_set, false);
    Ncat(i) = Y.col(i).cast<int>().sum();
  }

  // binom
  if(response == "binom") {
    yb = Y.col(1).cast<double>(); // booleans -> 0/1
  }
}


//------------------
//  Fit
//------------------

void FitOrdinal::setup_restart(const Eigen::VectorXd & _b0)
{
  unsigned Nb0 = _b0.size();
  if(Nb0 == 0) throw std::runtime_error("input b0 has size 0");

  apply_start = true;

  if(Nb0 == Nb) {
    b0 = _b0;
  } else if(Nb0 < Nb) {
    b0.resize(Nb);
    b0.setZero();
    b0.head(Nb0) = _b0;
  } else {
    throw std::runtime_error("Nb0 > Nb");
  }

  /* cout << " _b0 = " << _b0.transpose() << endl; */
  /* cout << " b0 = " << b0.transpose() << endl; */
}

void FitOrdinal::setup_offset_binom(const Eigen::VectorXd & _yo, bool decrement_Nb) 
{
  apply_offset = true;
  exclude_intercepts = true;
  yo = _yo;
  // Intercept is not modeled
  if(decrement_Nb) --Nb; 
}

void FitOrdinal::setup_offset_multinom_pom(const Eigen::VectorXd & _yo, const Eigen::VectorXd & _yo_int)
{
  apply_offset = true;
  exclude_intercepts = true;
  yo = _yo;
  yo_int = _yo_int;
}

void FitOrdinal::setup_ncov0(unsigned int _Ncov0, bool _last0, bool preproc_cov)
{
  Ncov0 = _Ncov0;
  last0 = _last0;

  if(Ncov0) {
    if(preproc_cov) throw std::runtime_error("preproc_cov is on when Ncov0 != 0");

    if(response == "multinom") {
      if(Ncov0 > Ncov) throw std::runtime_error("Ncov0 > Ncov (multinom)");
      Ncov1 = Ncov - Ncov0;
    } else if(response == "binom") {
      if(Ncov0 > Nb) throw std::runtime_error("Ncov0 > Nb (binom)");
      Ncov1 = Nb - Ncov0;
    } else {
      throw std::runtime_error("unknown response");
    }
  }
}

// main function
FitOrdinal Ordinal::fit(const Eigen::VectorXi &y, const Eigen::MatrixXd& X,
    unsigned int Ncov0, bool last0)
{
  // set up X & y
  setup_xy(y, X);

  // fit
  FitOrdinal fit = setup_fit();
  fit.setup_ncov0(Ncov0, last0, preproc_cov);

  if(response == "multinom") {
    fit.fit_multinom_pom(Mask, Y, Xcov);
  } else if(response == "binom") {
    fit.fit_binom(Mask, Y, Xcov1);
  } else {
    throw std::runtime_error("unknown response");
  }

  return(fit);
}

//----------------------------
//  Fit Common (FitOrdinal)
//----------------------------

void FitOrdinal::update_fit()
{
  bhat = cur_b;
  loglik = cur_loglik;
  Score = cur_Score; Info = cur_Info;
}

//----------------------------
//  Fit Multinom (FitOrdinal)
//----------------------------

// main fit multinom function
void FitOrdinal::fit_multinom_pom(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X)
{
  if(verbose >= 2) cout << "fit_multinom_pom\n"; 
  // start values of parameters
  check_setup_model();
  setup_start_multinom();
  // allocate memory for matrices used in the loop
  setup_par_multinom();
  check_setup_data();
  // optim
  converged = optimize(Mask, Y, X);
  // store results into fit
  update_fit();
  if(store_offset) {
    // linear predictor without intercepts
    if(Ncov) Xb0 = X * bhat.tail(Ncov);
    else Xb0.setZero();
    if(apply_offset) Xb0.array() += yo.array();
    yo = Mask.select(Xb0, 0.0); // overwrite offset vector yo if present
    // intercepts
    yo_int = bhat.head(ncat1);
  }
}

// set up starting values
void FitOrdinal::setup_start_multinom()
{
  if(verbose >= 2) cout << "setup_start_multinom\n"; 

  unsigned int i, i_cov, n_nom, n_denom;
  
  if(apply_start) {
    if(b0.size() != Nb) throw std::runtime_error("b0 is not of size Nb");
  } else {
    b0.resize(Nb);
    // initialize intercepts
    if(!exclude_intercepts) {
      // v1
      /* for(i = 0; i < ncat1; i++) { */
      /*   b0[i] = (double)(1 + i); */
      /* } */
      // v2
      for(i = 0, n_nom = 0, n_denom = Neff; i < ncat1; i++) {
        n_nom += Ncat(i);
        n_denom -= Ncat(i);
        b0[i] = log((double)(n_nom)/(double)(n_denom));
      }
      // v3
      /* Eigen::VectorXd Ncat_half(ncat); */
      /* for(i = 0; i < ncat; i++) Ncat_half(i) = (double)(Ncat(i)) + 0.5; */

      /* double n_nom_half, n_denom_half; */
      /* for(i = 0, n_nom_half = 0.0, n_denom_half= (double)(Neff) + ncat*0.5; i < ncat1; i++) { */
      /*   n_nom_half+= Ncat_half(i); */
      /*   n_denom_half -= Ncat_half(i); */
      /*   b0[i] = log(n_nom_half/n_denom_half); */
      /* } */
    }

    // initialize covariate effects
    i_cov = exclude_intercepts ? 0 : ncat1;
    for(i = i_cov; i < Nb; i++) {
      b0[i] = 0.0;
    }
  }

  if(Ncov0) {
    if(last0) b0.tail(Ncov0).setZero();
    else b0.head(Ncov0).setZero();
  }

  if(verbose >= 3) cout << " b0 = " << b0.transpose() << endl;
}

void FitOrdinal::setup_par_multinom()
{
  if(verbose >= 2) cout << "setup_par_multinom\n"; 

  Xb0.resize(N);
  Xb.resize(N, ncat1); exp_eta.resize(N, ncat1); gamma.resize(N, ncat1); PQ.resize(N, ncat1);
  P.resize(N, ncat1); P.setZero();
  Psum.resize(N); Pk.resize(N);

  D.resize(N, ncat1); D.setZero();
  V.resize(N, ncat1); V.setZero();

  Q.resize(N, ncat1sq); Q.setZero();
  S.resize(N, ncat1sq); S.setZero();
  QS.resize(N, ncat1sq); QS.setZero();
  W.resize(N, ncat1sq); W.setZero();

  WS2.resize(ncat1sq); WSS1.resize(N);
  if(Ncov) {
    XW.resize(Ncov, ncat1sq);
    XWs.resize(Ncov, ncat1);
    XW1.resize(Ncov);
    XW22.resize(Ncov, Ncov);
  }
  
  cur_Score.resize(Nb);
  cur_Info.resize(Nb, Nb);
  cur_v.resize(Nb); cur_b.resize(Nb);

  if(Ncov0) {
    if(last0) {
      cur_Score.tail(Ncov0).setZero();
      cur_v.tail(Ncov0).setZero();
      cur_b.tail(Ncov0).setZero();
    } else {
      cur_Score.head(Ncov0).setZero();
      cur_v.head(Ncov0).setZero();
      cur_b.head(Ncov0).setZero();
    }
  }

  if(firth_multinom) {
    Ystar.resize(N, ncat);
  }
}

bool FitOrdinal::update_par_multinom(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X, const Eigen::VectorXd & b, 
    bool pseudo)
{
  if(verbose >= 2) cout << " - update_par_multinom " << (pseudo ? "(pseudo)" : "") << "\n"; 

  unsigned int i, k,  l, m, start;

  VectorXd b_cov;
  if(exclude_intercepts) b_cov = b;
  else b_cov = b.segment(ncat1, Ncov);

  if(Ncov) {
    if(Ncov0) {
      if(last0) Xb0 = X.leftCols(Ncov1) * b_cov.head(Ncov1);
      else throw std::runtime_error("!last0 not implemented yet");
    } else Xb0 = X * b_cov;
  } else Xb0.setZero(); 

  if(apply_offset) Xb0.array() += yo.array(); // offset 

  // update linear predictor Xb with intercepts
  if(exclude_intercepts) for(i = 0; i < ncat1; i++) Xb.col(i).array() = Xb0.array();
  else {
    VectorXd b_int = b.head(ncat1);
    for(i = 0; i < ncat1; i++) Xb.col(i).array() = Xb0.array() + b_int(i);
  }
  if(apply_offset & !exclude_intercepts_offset) for(i = 0; i < ncat1; i++) Xb.col(i).array() += yo_int(i);

  exp_eta = Xb; exp_matrix_ord(exp_eta);
  gamma.array() = exp_eta.array() / (1.0 + exp_eta.array());

  P = gamma;
  for(i = 1; i < ncat1; i++) P.col(i).array() -= gamma.col(i - 1).array();
  Psum = P.rowwise().sum();
  if((Psum.array() >= 1.0).any()) {
    cerr << "WARNING: some elements in Psum >= 1.0" << endl;
    return(false);
  }

  Pk.array() = 1.0 - Psum.array();

  // interim computation of log-lik
  if(!pseudo) {
    cur_loglik = loglik_multinom(Mask, Y); // depends on Y, P, Pk, Mask
    if(verbose > 2) cout << "  -- (iterim) loglik: " << cur_loglik << endl;
    if(check_nan(cur_loglik)) {
      cerr << "WARNING: log-lik is NaN or Inf" << endl;
      return(false);
    }
  }
 
  // D = (Y[, -ncat] / P) - (Y[, ncat] / Pk)
  if(!pseudo) for(i = 0; i < ncat1; i++) D.col(i).array() = Y.col(i).cast<double>().array() / P.col(i).array() - Y.col(ncat1).cast<double>().array() / Pk.array();
  else for(i = 0; i < ncat1; i++) D.col(i).array() = Ystar.col(i).array() / P.col(i).array() - Ystar.col(ncat1).array() / Pk.array();

  // Q = dh / deta
  PQ.array() = gamma.array() * (1.0 - gamma.array());

  for(m = 0; m < ncat1; m++) {
    l = m; start = l * ncat1;
    Q.col(start + m).array() = PQ.col(m).array();
  }
  for(m = 1; m < ncat1; m++) {
    l = m - 1; start = l * ncat1;
    Q.col(start + m).array() = -1.0 * PQ.col(l).array();
  }

  // V
  for(k = 0; k < ncat1; k++) {
    // cols_Q = (k - 1)*ncat1 + seq(ncat1)
    // V[, k] = rowSums(D * Q[, cols_Q])
    V.col(k).array() = (D.array() * Q(all, seqN(k*ncat1, ncat1)).array()).rowwise().sum(); 
  }

  // Sinv or S = diag(1/p) + 1/(1 - sum_p)
  for(k = 0; k < ncat1sq; k++) { // go through all ncat1sq columns
    S.col(k) = Pk.array().inverse();
  }
  for(k = 0; k < ncat1; k++) { // go through ncat1 columns
    l = (ncat1 + 1)*k;
    S.col(l).array() += P.col(k).array().inverse();
  }

  // W = crossprod(Q, Sinv) %*% Q
  // (Q'S)'
  for(l = 0; l < ncat1; l++) {
    for(m = 0; m < ncat1; m++) {
      // col_QS = ncat1 * (l - 1) + m
      // cols_S = ncat1 * (m - 1) + seq(ncat1)
      // cols_Q = ncat1 * (l - 1) + seq(ncat1)
      // QS[, col_QS] = rowSums(Q[, cols_Q] * S[, cols_S])
      k = ncat1 * l + m; // col_QS
      QS.col(k).array() = (Q(all, seqN(l*ncat1, ncat1)).array() * S(all, seqN(ncat1*m, ncat1)).array()).rowwise().sum();
    }
  }
  for(l = 0; l < ncat1; l++) {
    for(m = 0; m < ncat1; m++) {
      // col_W = ncat1 * (l - 1) + m
      // cols_QS = ncat1 * (m - 1) + seq(ncat1)
      // cols_Q = ncat1 * (l - 1) + seq(ncat1)
      // W[, col_W] = rowSums(QS[, cols_QS] * Q[, cols_Q])
      k = ncat1 * l + m; // col_W
      W.col(k).array() = (QS(all, seqN(m*ncat1, ncat1)).array() * Q(all, seqN(ncat1*l, ncat1)).array()).rowwise().sum();
    }
  }

  // Account for miss. via Mask
  for(unsigned int i = 0; i < V.cols(); i++) V.col(i) = Mask.select(V.col(i), 0.0);
  for(unsigned int i = 0; i < W.cols(); i++) W.col(i) = Mask.select(W.col(i), 0.0);

  // Score
  // Score = c(colSums(V), colSums((crossprod(V, X))))
  if(!exclude_intercepts) cur_Score.head(ncat1) = V.colwise().sum();
  if(Ncov) {
    if(Ncov0) {
      if(last0) cur_Score.tail(Ncov).head(Ncov1) = (V.transpose() * X.leftCols(Ncov1)).colwise().sum();
      else throw std::runtime_error("!last0 not implemented yet");
    } else cur_Score.tail(Ncov) = (V.transpose() * X).colwise().sum();
  }

  // Info
  WS2 = W.colwise().sum();
  WSS1 = W.rowwise().sum().array().sqrt();
  // Info 1x1 block: Info[seq(ncat1), seq(ncat1)] = matrix(Wsum, ncat1, ncat1)
  if(!exclude_intercepts) cur_Info(seqN(0, ncat1), seqN(0, ncat1)) = Map<MatrixXd>(WS2.data(), ncat1, ncat1);

  if(Ncov) {
    if(!exclude_intercepts) {
      XW = X.transpose() * W;
      // Info 1x2 & 2x1 blocks
      for(k = 0; k < ncat1; k++) {
        // cols_XW = (k - 1) + seq(1, by = ncat1, length = ncat1)
        // Info[k, seq(ncat, nb)] = rowSums(XW[, cols_XW, drop = FALSE])
        XWs = XW(all, seqN(k, ncat1, ncat1));
        XW1 = XWs.rowwise().sum();
        cur_Info(k, seqN(ncat1, Ncov)).array() = XW1.array();
        cur_Info(seqN(ncat1, Ncov), k).array() = XW1.array();
      }
    }
    // Info 2x2 block
    // Info[seq(ncat, nb), seq(ncat, nb)] = crossprod(sqrt(Wsum1) * X)
    XW22 = MatrixXd(Ncov, Ncov).setZero().selfadjointView<Lower>().
      rankUpdate((X.array().colwise() * WSS1.array()).matrix().adjoint());
    if(exclude_intercepts) {
      cur_Info.array() = XW22.array();
    } else {
      cur_Info(seqN(ncat1, Ncov), seqN(ncat1, Ncov)).array() = XW22.array();
    }
  }

  // solve: v = solve(Info, Score)
  // Firth correction?
  LLT<MatrixXd> llt(cur_Info);

  if(verbose > 2) {
    ColPivHouseholderQR<MatrixXd> qr;
    qr.compute(cur_Info);
    cout << " qr.isInvertible(cur_info) = " << qr.isInvertible() << endl;
  }

  if(!firth_multinom | pseudo) {
    if(Ncov0) {
      if(last0) {
        LLT<MatrixXd> llt1(cur_Info.block(0, 0, ncat1 + Ncov1, ncat1 + Ncov1));
        cur_v.head(ncat1 + Ncov1) = llt1.solve(cur_Score.head(ncat1 + Ncov1));
      } else throw std::runtime_error("!last0 not implemented yet");
    } else cur_v = llt.solve(cur_Score);
  } else {
    /* MatrixXd cur_Info_inv = llt.solve(MatrixXd::Identity(Nb, Nb)); */
    MatrixXd cur_Info_inv;
    if(Ncov0) {
      if(last0) cur_Info_inv = cur_Info.block(0, 0, ncat1 + Ncov1, ncat1 + Ncov1).inverse();
      else throw std::runtime_error("!last0 not implemented yet");
    } else {
      // v1
      /* cur_Info_inv = cur_Info.inverse(); */ 
      // v2
      ColPivHouseholderQR<MatrixXd> qr;
      qr.compute(cur_Info);
      /* if(!qr.isInvertible()) { */
      /*   cerr << "WARNING: Info is not invertible" << endl; */
      /*   return(false); */
      /* } */
      cur_Info_inv = qr.inverse();
    }
    if(verbose > 2) cout << "  -- cur_Info_inv: " << cur_Info_inv.rows() << "x" << cur_Info_inv.cols() << endl;

    MatrixXd diagA(N, ncat1); diagA.setZero();
    if(Ncov0) {
      if(last0) {
        if(Ncov1) diagA = 2 * (X.leftCols(Ncov1) * cur_Info_inv(seqN(ncat1, Ncov1), seqN(0, ncat1))); 
        for(i = 0; i < ncat1; i++) diagA.col(i).array() += cur_Info_inv(i, i);
        if(Ncov1) diagA.array().colwise() += ((X.leftCols(Ncov1) * cur_Info_inv(seqN(ncat1, Ncov1), seqN(ncat1, Ncov1))).array() * X.leftCols(Ncov1).array()).rowwise().sum();
      } else throw std::runtime_error("!last0 not implemented yet");
    } else {
      if(!exclude_intercepts) {
        if(Ncov) diagA = 2 * (X * cur_Info_inv(seqN(ncat1, Ncov), seqN(0, ncat1))); 
        for(i = 0; i < ncat1; i++) diagA.col(i).array() += cur_Info_inv(i, i);
        if(Ncov) diagA.array().colwise() += ((X * cur_Info_inv(seqN(ncat1, Ncov), seqN(ncat1, Ncov))).array() * X.array()).rowwise().sum();
      } else {
        if(Ncov) diagA.array().colwise() += ((X * cur_Info_inv).array() * X.array()).rowwise().sum();
      }
    }
    if(verbose > 2) cout << "  -- diagA: " << diagA.rows() << "x" << diagA.cols() << endl;

    MatrixXd adj_c = 0.5 * diagA.array() * dlog_matrix(Xb).array();

    // adjustment to counts
    MatrixXd adj_a(N, ncat);
    adj_a.leftCols(ncat1) = adj_c; 
    adj_a.col(ncat1) *= 0;
    adj_a.rightCols(ncat1).array() -= adj_c.array();

    // Yadj = Y + adj_a
    Ystar = Y.array().cast<double>(); 
    Ystar.array() += adj_a.array();

    // re-compute D, V (Q doesn't change as it is function of probs.)
    for(i = 0; i < ncat1; i++) D.col(i).array() = Ystar.col(i).array() / P.col(i).array() - Ystar.col(ncat1).array() / Pk.array();
    for(k = 0; k < ncat1; k++) V.col(k).array() = (D.array() * Q(all, seqN(k*ncat1, ncat1)).array()).rowwise().sum(); 
    // account for miss. via Mask
    for(unsigned int i = 0; i < V.cols(); i++) V.col(i) = Mask.select(V.col(i), 0.0);

    // re-compmute Scores
    if(!exclude_intercepts) cur_Score.head(ncat1) = V.colwise().sum();
    if(Ncov) {
      if(Ncov0) {
        if(last0) cur_Score.tail(Ncov).head(Ncov1) = (V.transpose() * X.leftCols(Ncov1)).colwise().sum();
        else throw std::runtime_error("!last0 not implemented yet");
      } else cur_Score.tail(Ncov) = (V.transpose() * X).colwise().sum();
    }
    if(verbose > 2) cout << "  -- cur_Score: " << cur_Score.transpose() << endl;

    if(Ncov0) {
      if(last0) cur_v.head(ncat1 + Ncov1) = cur_Info_inv * cur_Score.head(ncat1 + Ncov1);
      else throw std::runtime_error("!last0 not implemented yet");
    } else cur_v = cur_Info_inv * cur_Score;
    if(verbose > 2) cout << "  -- cur_v: " << cur_v.transpose() << endl;
  }

  if(!pseudo) {
    if(firth_multinom) cur_loglik = loglik_multinom_firth(Mask, Y, llt, true, cur_loglik); // add = true
    cur_dev = -2.0 * cur_loglik;
  }

  // dump
  if(verbose >= 3) {
    bool append = true; // (cnt_updates > 1);
    string name = "ordinal.txt";
    ofstream file;
    if(append) file.open(name.c_str(), ios::out | ios::app);
    else file.open(name.c_str(), ios::out);

    double diff = cur_Score.array().abs().maxCoeff(); 
    file << cnt_fit  << " " << b.size() << " " << cnt_updates << " " << cur_loglik << " " << cur_dev
      << " " << diff; 
    for(unsigned int i = 0; i < b.size(); i++) {
      file << " " << b(i);
    }
    file << endl;
    file.close();

    if(verbose >= 4) {
      const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
      ofstream yfile("y.txt");
      yfile << Y.col(1).format(CSVFormat);
      yfile.close();
      ofstream xfile("x.txt");
      xfile << X.format(CSVFormat);
      xfile.close();
      ofstream mfile("m.txt");
      mfile << Mask.format(CSVFormat);
      mfile.close();

      throw std::runtime_error("verbose level 4: exit after dumping data into y.txt, x.txt and m.txt");
    }
  }

  return(true);
}

double FitOrdinal::loglik_multinom(const VectorXb & Mask, const MatrixXb & Y)
{
  // loglik = colSums(Y*log(cbind(P, Pk))) %>% sum
  double res = 0.0;
  // 1, 2, ..., ncat1 categories
  for(unsigned int i = 0; i < ncat1; i++) {
    res += Mask.select(Y.col(i).select(P.col(i).array().log(), 0.0), 0.0).array().sum();
  }
  // the last ncat category
  res += Mask.select(Y.col(ncat1).select(Pk.array().log(), 0.0), 0.0).array().sum();

  return(res);
}

double FitOrdinal::loglik_multinom_firth(const VectorXb & Mask, const MatrixXb & Y, const LLT<MatrixXd> & llt,
    bool add, double base_loglik)
{
  double res = add ? base_loglik : loglik_multinom(Mask, Y);

  // https://gist.github.com/redpony/fc8a0db6b20f7b1a3f23
  double half_logdet = llt.matrixL().toDenseMatrix().diagonal().array().log().sum();
  /* double half_logdet = log(llt.matrixL().determinant()); */

  res += firth_mult * half_logdet;

  return(res);
}

double FitOrdinal::loglik_multinom_firth(const VectorXb & Mask, const MatrixXb & Y, const MatrixXd & Info)
{
  double res = loglik_multinom(Mask, Y);

  double half_logdet = 0.5*log(Info.determinant());

  res += firth_mult * half_logdet;

  return(res);
}

//----------------------------
//  Fit Binom (FitOrdinal)
//----------------------------

// main fit binom function
void FitOrdinal::fit_binom(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X)
{
  if(verbose >= 2) cout << "fit_binom\n"; 
  /* // start values of parameters */
  check_setup_model();
  setup_start_binom();
  // allocate memory for matrices used in the loop
  setup_par_binom();
  check_setup_data();
  // optim
  converged = optimize(Mask, Y, X);
  // store results into fit
  update_fit();
  // store offset if specified
  if(store_offset) {
    Xb = X * bhat;
    if(apply_offset) Xb.array() += yo.array();
    yo = Mask.select(Xb, 0.0); // overwrite offset vector yo if present
  }
}

// set up starting values
void FitOrdinal::setup_start_binom()
{
  if(verbose >= 2) cout << "setup_start_binom\n";

  unsigned int i;
  
  if(apply_start) {
    if(b0.size() != Nb) throw std::runtime_error("b0 is not of size Nb");
  } else {
    b0.resize(Nb);
    // intercepts
    b0(0) = log((double)(Ncat(1)/(double)(Ncat(0)))); // log(n1/n0)
    // covariates
    for(i = ncat1; i < Nb; i++) {
      b0[i] = 0.0;
    }
  }

  if(Ncov0) {
    if(last0) b0.tail(Ncov0).setZero();
    else b0.head(Ncov0).setZero();
  }

  if(verbose >= 3) cout << " b0 = " << b0.transpose() << endl;

}

void FitOrdinal::setup_par_binom()
{
  if(verbose >= 2) cout << "setup_par_binom\n";

  mub.resize(N); wb.resize(N);
  XtW.resize(Nb, N);

  cur_Score.resize(Nb);
  cur_Info.resize(Nb, Nb);
  cur_v.resize(Nb); cur_b.resize(Nb);

  if(Ncov0) {
    if(last0) {
      cur_Score.tail(Ncov0).setZero();
      cur_v.tail(Ncov0).setZero();
      cur_b.tail(Ncov0).setZero();
    } else {
      cur_Score.head(Ncov0).setZero();
      cur_v.head(Ncov0).setZero();
      cur_b.head(Ncov0).setZero();
    }
  }

  if(firth_binom) {
    ystar.resize(N);
  }
}

bool FitOrdinal::update_par_binom(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X, const Eigen::VectorXd & b)
{
  if(verbose >= 2) cout << " - update_par_binom\n"; 

  mub = X * b;
  if(apply_offset) mub.array() += yo.array();
  exp_vector(mub); // mub <- exp(mub)
  mub.array() /= (1.0 + mub.array()); // mub <- exp(mub) / (1 + exp(mub))

  wb.array() = Mask.select(mub.array() * (1.0 - mub.array()), 1.0);

  // Score: Score = c(sum(y - mu), X' (y - mu)) 
  cur_Score = X.transpose() * Mask.select((Y.col(1).array().cast<double>() - mub.array()), 0.0).matrix();

  // Info
  cur_Info = X.transpose() * (X.array().colwise() * wb.array()).matrix();

  // solve: v = solve(Info, Score)
  LLT<MatrixXd> llt(cur_Info);
  cur_v = llt.solve(cur_Score);

  cur_loglik = loglik_binom(Mask, Y);
  cur_dev = -2.0 * cur_loglik;

  // dump
  if(verbose >= 3) {
    bool append = true; // (cnt_updates > 1);
    string name = "ordinal.txt";
    ofstream file;
    if(append) file.open(name.c_str(), ios::out | ios::app);
    else file.open(name.c_str(), ios::out);

    double diff = cur_Score.array().abs().maxCoeff(); 
    file << cnt_fit << " " << b.size() << " " << cnt_updates << " " << cur_loglik << " " << cur_dev
      << " " << diff; 
    for(unsigned int i = 0; i < b.size(); i++) {
      file << " " << b(i);
    }
    file << endl;
    file.close();
  }

  return(true);
}

bool FitOrdinal::update_par_binom_firth(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X, const Eigen::VectorXd & b)
{
  if(verbose >= 2) cout << " - update_par_binom_firth\n"; 

  /* cout << Ncov0 << " " << last0 << " " << Ncov1 << endl; */
  /* cout << b.size() << " " << X.cols() << " " << X.rows() << endl; */
  if(Ncov0) {
    if(last0) mub = X.leftCols(Ncov1) * b.head(Ncov1);
    else mub = X.rightCols(Ncov1) * b.tail(Ncov1);
  } else {
    mub = X * b;
  }
  if(apply_offset) mub.array() += yo.array();
  exp_vector(mub); // mub <- exp(mub)
  mub.array() /= (1.0 + mub.array()); // mub <- exp(mub) / (1 + exp(mub))

  wb.array() = Mask.select(mub.array() * (1.0 - mub.array()), 1.0);

  // Info = X' W X
  XtW = X.transpose() * wb.array().sqrt().matrix().asDiagonal();
  cur_Info = XtW * XtW.transpose();

  // derive h
  LLT<MatrixXd> llt(cur_Info);
  VectorXd h = (llt.solve(XtW).array() * XtW.array()).colwise().sum();

  // derive pseudo response
  ystar = Y.col(1).cast<double>().array() + firth_mult * h.array() * (0.5 - mub.array());

  // update Score = Ab + Sb
  /* Ab = crossprod(X, h * (0.5 - mu)) */
  /* Sb = crossprod(X, y - mu) */
  // solve: v = solve(Info, Score)

  if(Ncov0) {
    if(last0) {
      LLT<MatrixXd> llt1(cur_Info.block(0, 0, Ncov1, Ncov1));
      cur_Score.head(Ncov1) = (X.leftCols(Ncov1).transpose() * Mask.select(ystar.array() - mub.array(), 0.0).matrix()).array();
      cur_v.head(Ncov1) = llt1.solve(cur_Score.head(Ncov1));
    } else {
      LLT<MatrixXd> llt1(cur_Info.block(Ncov0, Ncov0, Ncov1, Ncov1));
      cur_Score.tail(Ncov1) = (X.rightCols(Ncov1).transpose() * Mask.select(ystar.array() - mub.array(), 0.0).matrix()).array();
      cur_v.tail(Ncov1) = llt1.solve(cur_Score.tail(Ncov1));
    }
  } else {
    cur_Score = (X.transpose() * Mask.select(ystar.array() - mub.array(), 0.0).matrix()).array();
    cur_v = llt.solve(cur_Score);
  }

  cur_loglik = loglik_binom_firth(Mask, Y, llt);
  cur_dev = -2.0 * cur_loglik;

  // dump
  if(verbose >= 3) {
    bool append = true; // (cnt_updates > 1);
    string name = "ordinal.txt";
    ofstream file;
    if(append) file.open(name.c_str(), ios::out | ios::app);
    else file.open(name.c_str(), ios::out);

    double diff = cur_Score.array().abs().maxCoeff(); 
    file << cnt_fit << " " << b.size() << " " << cnt_updates << " " << cur_loglik << " " << cur_dev
      << " " << diff; 
    for(unsigned int i = 0; i < b.size(); i++) {
      file << " " << b(i);
    }
    file << endl;
    file.close();

    if(verbose >= 4) {
      const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
      ofstream yfile("y.txt");
      yfile << Y.col(1).format(CSVFormat);
      yfile.close();
      ofstream xfile("x.txt");
      xfile << X.format(CSVFormat);
      xfile.close();
      ofstream mfile("m.txt");
      mfile << Mask.format(CSVFormat);
      mfile.close();
      
      throw std::runtime_error("verbose level 4: exit after dumping data into y.txt, x.txt and m.txt");
    }
  }

  return(true);
}

bool FitOrdinal::update_par_binom_pseudo(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X, const Eigen::VectorXd & b)
{
  if(verbose >= 2) cout << " - update_par_binom_pseudo\n"; 

  if(Ncov0) {
    if(last0) mub = X.leftCols(Ncov1) * b.head(Ncov1);
    else mub = X.rightCols(Ncov1) * b.tail(Ncov1);
  } else {
    mub = X * b;
  }
  if(apply_offset) mub.array() += yo.array();
  exp_vector(mub); // mub <- exp(mub)
  mub.array() /= (1.0 + mub.array()); // mub <- exp(mub) / (1 + exp(mub))

  wb.array() = Mask.select(mub.array() * (1.0 - mub.array()), 1.0);

  // Info = X' W X
  XtW = X.transpose() * wb.array().sqrt().matrix().asDiagonal();
  cur_Info = XtW * XtW.transpose();

  if(Ncov0) {
    if(last0) {
      LLT<MatrixXd> llt1(cur_Info.block(0, 0, Ncov1, Ncov1));
      cur_Score.head(Ncov1) = X.leftCols(Ncov1).transpose() * Mask.select(ystar.array() - mub.array(), 0.0).matrix();
      cur_v.head(Ncov1) = llt1.solve(cur_Score.head(Ncov1));
    } else {
      LLT<MatrixXd> llt1(cur_Info.block(Ncov0, Ncov0, Ncov1, Ncov1));
      cur_Score.tail(Ncov1) = X.rightCols(Ncov1).transpose() * Mask.select(ystar.array() - mub.array(), 0.0).matrix();
      cur_v.tail(Ncov1) = llt1.solve(cur_Score.tail(Ncov1));
    }
  } else {
    LLT<MatrixXd> llt(cur_Info);
    cur_Score = X.transpose() * Mask.select(ystar.array() - mub.array(), 0.0).matrix();
    cur_v = llt.solve(cur_Score);
  }

  // What is dev/loglik for pseudo response model?
  /* cur_loglik = loglik_binom_firth(Mask, Y, llt); */
  /* cur_dev = -2.0 * cur_loglik; */

  return(true);
}

double FitOrdinal::loglik_binom(const VectorXb & Mask, const MatrixXb & Y)
{
  // loglik = sum(Y*log(p) + (1-Y)*log(1-p)), where p = mu
  double res = 0.0;
  // contols
  res += Mask.select(Y.col(0).select((1.0 - mub.array()).log(), 0.0), 0.0).array().sum();
  // cases
  res += Mask.select(Y.col(1).select(mub.array().log(), 0.0), 0.0).array().sum();

  return(res);
}

double FitOrdinal::loglik_binom_firth(const VectorXb & Mask, const MatrixXb & Y, const LLT<MatrixXd> & llt)
{
  double res = loglik_binom(Mask, Y);

  // https://gist.github.com/redpony/fc8a0db6b20f7b1a3f23
  double half_logdet = llt.matrixL().toDenseMatrix().diagonal().array().log().sum();
  /* double half_logdet = log(llt.matrixL().determinant()); */

  /* cout << "firth_mult = " << firth_mult << endl; */
  res += firth_mult * half_logdet;

  return(res);
}

bool FitOrdinal::stop_criterion()
{
  bool stop;
  // v1: abs. diff. between bhat_cur and bhat_prev
  // stop = (cur_v.norm() < tol);
  // v2: abs. max. value of Scores. 
  // - Example: Regenie Firth model fitting 
  stop_value = cur_Score.array().abs().maxCoeff();
  stop = (stop_value < tol); 
  // v3: relative diff. in deviance. 
  // - Example: glm2::glm.fit2.R
  /* stop = ((abs(cur_dev - prev_dev) / (1.0 + abs(prev_dev))) < tol); // 1.0 to prevent from division by zero */

  if(verbose >= 2) cout << "stop_criterion: " << stop_value << " < " << tol << endl;

  return(stop);
}

bool FitOrdinal::update_par(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X, const Eigen::VectorXd & b, bool pseudo)
{
  bool ret;

  // count
  if(trace) { cnt_updates++; }

  // store prev. values 
  prev_dev = cur_dev;

  // update directions
  if(response == "multinom") {
    if(pseudo) ret = update_par_multinom(Mask, Y, X, b, true); // pseudo = true
    else ret = update_par_multinom(Mask, Y, X, b);
  } else if(response == "binom") {
    if(firth_binom) {
      if(pseudo) ret = update_par_binom_pseudo(Mask, Y, X, b);
      else ret = update_par_binom_firth(Mask, Y, X, b);
    } else {
      ret = update_par_binom(Mask, Y, X, b);
    }
  } else {
    throw std::runtime_error("unknown response");
  }

  return(ret);
}

// optimization loop for binom. response
bool FitOrdinal::optimize(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X)
{
  if(verbose >= 2) cout << "optimize\n"; 

  bool res;

  // special case when firth_multinom is on
  /* if(firth_multinom) { */
  /*   res = optimize_FisherScoring(Mask, Y, X); */
  /* } else */ 
  if(optim == "FisherScoring") {
    res = optimize_FisherScoring(Mask, Y, X);
  } else if(optim == "FisherScoringPseudo") {
    if(firth_binom | firth_multinom) res = optimize_FisherScoringPseudo(Mask, Y, X);
    else res = optimize_FisherScoring(Mask, Y, X);
  } else if(optim == "WeightHalving") {
    res = optimize_WeightHalving(Mask, Y, X);
  } else if(optim == "WeightHalvingPseudo") {
    if(firth_binom | firth_multinom) res = optimize_WeightHalvingPseudo(Mask, Y, X);
    else res = optimize_WeightHalving(Mask, Y, X);
  } else {
    throw std::runtime_error("unknown optimize");
  }

  if(verbose >= 2) cout << "optimize it = " << it << " | cnt_updates = " << cnt_updates << endl;
  if(verbose >= 3) cout << "bhat = " << cur_b.transpose() << endl;
  if(verbose >= 3) cout << "converged = " << converged << endl;

  return(res);
}

bool FitOrdinal::optimize_FisherScoring(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X)
{
  unsigned int i;
  bool res, up = true;
  double ratio_step;

  cur_b = b0;
  for(i = 0; i < maxit; i++) {
    // update directions
    up = update_par(Mask, Y, X, cur_b);
    if(!up) break;

    // check the stopping criteria
    if(stop_criterion()) break;

    // check the  absolute step size to be less than max_step for each entry of step (cur_v2)
    if(check_step) {
      ratio_step = cur_v.array().abs().maxCoeff() / max_step;
      if(verbose >= 3) cout << " cur_v.array().abs().maxCoeff() = " << cur_v.array().abs().maxCoeff() << endl;
      if(verbose >= 3) cout << " ratio_step = " << ratio_step << endl;
      if(ratio_step > 1.0) {
        cur_v.array() /= ratio_step;
      }
      if(verbose >= 3) cout << " cur_v.array().abs().maxCoeff() = " << cur_v.array().abs().maxCoeff() << endl;
    }

    // update parameters
    cur_b += cur_v;

    // check if bhat is nan
    if(cur_b.array().isNaN().any()) { it = i; return false; }

  }
  // assign # iterations performed
  it = i;
  res = (it < maxit) & up;

  // check if any NaN
  if(cur_Score.array().isNaN().any() | cur_b.array().isNaN().any() | isnan(cur_dev)) {
    return false; 
  }

  return(res);
}

bool FitOrdinal::optimize_FisherScoringPseudo(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X)
{
  unsigned int i, i3;
  bool res, up = true;
  double ratio_step;

  if(!(firth_binom | firth_multinom)) throw std::runtime_error("optimize_FisherScoringPseudo is for binomial/multinomial response with Firth correction");

  cur_b = b0;
  for(i = 0; i < maxit; i++) {
    // update directions
    up = update_par(Mask, Y, X, cur_b);
    if(!up) break;

    // check the stopping criteria
    if(stop_criterion()) break;

    // check the  absolute step size to be less than max_step for each entry of step (cur_v2)
    if(check_step) {
      ratio_step = cur_v.array().abs().maxCoeff() / max_step;
      if(verbose >= 3) cout << " cur_v.array().abs().maxCoeff() = " << cur_v.array().abs().maxCoeff() << endl;
      if(verbose >= 3) cout << " ratio_step = " << ratio_step << endl;
      if(ratio_step > 1.0) cur_v.array() /= ratio_step;
      if(verbose >= 3) cout << " cur_v.array().abs().maxCoeff() = " << cur_v.array().abs().maxCoeff() << endl;
    }

    // ystar is derived & stored in update_par_binom_firth

    // Pseudo loop
    for(i3 = 0; i3 < maxit3; i3++) {
      if(verbose >= 3) cout << " - pseudo loop it " << i3 << endl;

      // update directions
      up = update_par(Mask, Y, X, cur_b, true); // pseudo = true
      if(!up) { it = i; return false; }

      // check the stopping criteria
      if(stop_criterion()) break;

      // check step size
      if(check_step) {
        ratio_step = cur_v.array().abs().maxCoeff() / max_step;
        if(ratio_step > 1.0) cur_v.array() /= ratio_step;
      }

      // update parameters
      cur_b += cur_v;

      // check if bhat is nan
      if(cur_b.array().isNaN().any()) { it = i; return false; }
    } // end of pseudo loop 
  } // end of main loop

  // assign # iterations performed
  it = i;
  res = (i < maxit) & up;

  // check if any NaN
  if(cur_Score.array().isNaN().any() | cur_b.array().isNaN().any() | isnan(cur_dev)) {
    return false; 
  }

  return(res);
}

bool FitOrdinal::optimize_WeightHalving(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X)
{
  if(verbose >= 2) cout << "optimize_WeightHalving\n";

  unsigned int i, i2;
  bool res, up = true;

  // declare variables
  VectorXd cur_b2, cur_v2;
  double cur_dev2;
  double denom, ratio_step;

  // initial values for Outer loop
  cur_b = b0;
  up = update_par(Mask, Y, X, cur_b); // get (i) the step size (cur_v); (ii) current value of dev
  if(!up) { it = 0; return false; }

  // Outer loop
  for(i = 1; i < maxit; i++) {
    if(verbose >= 3) cout << " - outer loop it " << i << endl;

    if(verbose >= 3) cout << " - cur_b = " << cur_b.transpose() << endl;
    if(verbose >= 3) cout << " - cur_Score = " << cur_Score.transpose() << " [max = " << cur_Score.array().abs().maxCoeff() << "]" << endl;

    // stopping criteria is checked here in the beginning rather than in the end of loop
    // - reason: update_par was called above
    // - example: i = 0 & no covariates & initial values are proportions --> no optimization inside the Outer/Inner loop is required
    if(stop_criterion()) break;

    // initial values for Inner loop
    cur_b2 = cur_b; cur_v2 = cur_v;
    cur_dev2 = cur_dev;

    denom = 2.0;
    
    // Inner loop (step halving)
    for(i2 = 0; i2 < maxit2; i2++) {
      if(verbose >= 3 && i2 > 0) cout << " - inner loop it " << i2 << endl;

      // update step according to the rule: step(i) = step(initial) / 2^i
      // one exception from the rule: skip halving at the very first iteration, i2 = 0 --> Fisher Scoring at i2 = 0
      if(i2) cur_v2.array() /= denom;

      // check the  absolute step size to be less than max_step for each entry of step (cur_v2)
      if(check_step) {
        ratio_step = cur_v2.array().abs().maxCoeff() / max_step;
        if(ratio_step > 1.0) cur_v2.array() /= ratio_step;
      }

      // update param.
      // - the baseline value (cur_b2) is fixed
      // - the increment step (cur_v2) is reduced at each iteration (see the code line above)
      cur_b = cur_b2 + cur_v2;

      // check if Score is nan
      if(cur_Score.array().isNaN().any()) { it = i; return false; }
      // check if bhat is nan
      if(cur_b.array().isNaN().any()) { it = i; return false; }

      // update Score, Info, loglik, dev
      up = update_par(Mask, Y, X, cur_b);
      if(!up) {it = i; return false; }

      // check if cur_dev is nan
      if(isnan(cur_dev)) { it = i; return false; }

      // stop the inner loop (step halving) if dev. is improved
      if(cur_dev < cur_dev2) break;
    }

    // assign Innter loop iterations
    it2 += i2;

    /* bool strict_WeightHalving = false; */
    if(strict) {
      // check if all Inner loop iterations are used & exit
      if(i2 == maxit2) {
        // let the first iteration (i = 0) go even when convergence failure
        if(i) { it = i; return false; }
      }
    }
  }

  // assign # iterations performed
  it = i;
  res = (i < maxit) & up;

  // check if any NaN
  if(cur_Score.array().isNaN().any() | cur_b.array().isNaN().any() | isnan(cur_dev)) {
    return false; 
  }

  return res;
}

bool FitOrdinal::optimize_WeightHalvingPseudo(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X)
{
  if(verbose >= 2) cout << "optimize_WeightHalvingPseudo\n";

  if(!(firth_binom | firth_multinom)) throw std::runtime_error("optimize_WeightHalvingPseudo is for binomial/multinomial response with Firth correction");

  unsigned int i, i2, i3;
  bool res, stop, up = true;

  // declare variables
  VectorXd cur_b2, cur_v2;
  double cur_dev2;
  double denom, ratio_step;

  // initial values for Outer loop
  cur_b = b0;
  /* update_par(Mask, Y, X, cur_b); // get (i) the step size (cur_v); (ii) current value of dev */

  // Outer loop
  for(i = 1; i < maxit; i++) {
    up = update_par(Mask, Y, X, cur_b); 
    if(!up) break;

    if(verbose >= 3) cout << " - outer loop it " << i << endl;
    if(verbose >= 3) cout << " - cur_b = " << cur_b.transpose() << endl;
    if(verbose >= 3) cout << " - cur_Score = " << cur_Score.transpose() << " [max = " << cur_Score.array().abs().maxCoeff() << "]" << endl;

    // stopping criterion
    if(stop_criterion()) break;

    if(stop_value > pseudo_stophalf) {
      // initial values for Inner loop
      cur_b2 = cur_b; cur_v2 = cur_v;
      cur_dev2 = cur_dev;

      denom = 2.0;
      
      // Inner loop (step halving)
      for(i2 = 0; i2 < maxit2; i2++) {
        if(verbose >= 3 && i2 > 0) cout << " - inner loop it " << i2 << endl;

        // update step according to the rule: step(i) = step(initial) / 2^i
        // one exception from the rule: skip halving at the very first iteration, i2 = 0 --> Fisher Scoring at i2 = 0
        if(i2)  cur_v2.array() /= denom;

        // check the  absolute step size to be less than max_step for each entry of step (cur_v2)
        if(check_step) {
          ratio_step = cur_v2.array().abs().maxCoeff() / max_step;
          if(ratio_step > 1.0) cur_v2.array() /= ratio_step;
        }

        // update param.
        // - the baseline value (cur_b2) is fixed
        // - the increment step (cur_v2) is reduced at each iteration (see the code line above)
        cur_b = cur_b2 + cur_v2;

        // check if Score is nan
        /* if(cur_Score.array().isNaN().any()) { it = i; return false; } */
        // check if bhat is nan
        /* if(cur_b.array().isNaN().any()) { it = i; return false; } */

        // update Score, Info, loglik, dev
        up = update_par(Mask, Y, X, cur_b);
        /* if(!up) { it = i; return false; } */
        if(!up) continue;

        // check if cur_dev is nan
        /* if(isnan(cur_dev)) { it = i; return false; } */

        // stop the inner loop (step halving) if dev.is improved
        if(cur_dev < cur_dev2) {
          break;
        }
      } // end of inner loop

      // assign Inner loop iterations
      it2 += i2;
    } else { // end of condition to start inner loop
      // check step size
      if(check_step) {
        ratio_step = cur_v.array().abs().maxCoeff() / max_step;
        if(ratio_step > 1.0) cur_v.array() /= ratio_step;
      }

      // update parameters
      cur_b += cur_v;
    }

    // ystar is derived & stored in update_par_binom_firth

    // Pseudo loop
    // store initial values before entering Pseudo loop
    cur_b2 = cur_b; cur_v2 = cur_v; cur_dev2 = cur_dev;
    bool loop_pseudo = false;
    for(i3 = 0; i3 < maxit3; i3++) {
      if(verbose >= 3) cout << " - pseudo loop it " << i3 << endl;

      // update directions
      up = update_par(Mask, Y, X, cur_b, true); // pseudo = true
      if(!up) break;

      // check the stopping criteria
      stop = stop_criterion();
      if(check_nan(stop_value)) break;
      if(stop) { loop_pseudo = true; break; }

      // check step size
      if(check_step) {
        ratio_step = cur_v.array().abs().maxCoeff() / max_step;
        if(ratio_step > 1.0) cur_v.array() /= ratio_step;
      }

      // update parameters
      cur_b += cur_v;

      // check if any NaN
      /* if(cur_Score.array().isNaN().any() | cur_b.array().isNaN().any()) return false; */ 
    } // end of pseudo loop (i3)
    // cancel results of Pseudo loop if it failed
    if(!loop_pseudo) {
      cur_b = cur_b2; cur_v = cur_v2; cur_dev = cur_dev2;
    }
  }

  // assign # iterations performed
  it = i;
  res = (i < maxit) & up;

  // check if any NaN
  if(cur_Score.array().isNaN().any() | cur_b.array().isNaN().any() | isnan(cur_dev)) {
    return false; 
  }

  return res;
}

//------------------------------
// Utils
//------------------------------

bool check_nan(double x) 
{
  /* bool res = (boost::math::isnan)(x); */
  bool res = (boost::math::isnan)(x) | !(boost::math::isnormal)(x);
  return(res);
}

void exp_matrix(Eigen::MatrixXd & X)
{
  // See: mu = binomial()$linkinv(eta)
  // - https://github.com/wch/r-source/blob/trunk/src/library/stats/src/family.c
  // - https://stackoverflow.com/a/1566222
  /* double EPSILON = 2.221e-16; */
  double EPSILON = 10 * std::numeric_limits<double>::epsilon();
  double THRESH = 30.0, MTHRESH = -30.0; 
  double INVEPS = 1.0/EPSILON;
  for(unsigned int i = 0; i < X.cols(); i++) {
    X.col(i).array() = (X.col(i).array() < MTHRESH).
      select(EPSILON, (X.col(i).array() > THRESH).
          select(INVEPS, X.col(i).array().exp()));
  }
}

void exp_matrix_ord(Eigen::MatrixXd & X)
{
  double EPSILON = 10 * std::numeric_limits<double>::epsilon();
  double THRESH = 30.0, MTHRESH = -30.0; 
  double INVEPS = 1.0/EPSILON;

  unsigned int ncols = X.cols();
  ArrayXb mask_top = (X.array() > THRESH).rowwise().all(), mask_bottom = (X.array() < MTHRESH).rowwise().all();
  /* cout << "any(mask_top) = " << mask_top.any() << endl; */

  for(unsigned int i = 0; i < ncols; i++) {
    X.col(i).array() = (X.col(i).array() < MTHRESH).
      select(EPSILON, (X.col(i).array() > THRESH).
          select(INVEPS, X.col(i).array().exp()));
  }
  // correction for cases: eta{1,2} > THRESH & eta1 < eta2 
  if(ncols > 1) {
    if(mask_top.any()) { // scaling factor (0.5, 1) for columns 1,2
      for(unsigned int i = 0; i < ncols; i++) {
        double sc = std::pow(0.5, ncols - 1 - i);
        X.col(i).array() = mask_top.select(sc * X.col(i).array(), X.col(i).array());
      }
    }
    /* if(mask_bottom.any()) { */
    /*   for(unsigned int i = 0; i < ncols; i++) { */
    /*     double sc = std::pow(0.5, i); */
    /*     X.col(i).array() = mask_bottom.select(sc * X.col(i).array(), X.col(i).array()); */
    /*   } */
    /* } */
  }
}

void invlogit_matrix(Eigen::MatrixXd & X)
{
  double EPSILON = 10 * std::numeric_limits<double>::epsilon();
  /* double EPSILON = 2.221e-16; */
  double THRESH = 30.0, MTHRESH = -30.0; 
  for(unsigned int i = 0; i < X.cols(); i++) {
    X.col(i).array() = (X.col(i).array() > MTHRESH).
      select(1.0 / (1.0 + EPSILON), (X.col(i).array() < THRESH).
          select(EPSILON / (1.0 + EPSILON), 1.0 - 1.0 / (1.0 + X.col(i).array().exp())));
  }
}

void exp_vector(Eigen::VectorXd & x)
{
  double EPSILON = 2.221e-16;
  double THRESH = 30.0, MTHRESH = -30.0; 
  double INVEPS = 1.0/EPSILON;
  x.array() = (x.array() < MTHRESH).
    select(EPSILON, (x.array() > THRESH).
        select(INVEPS, x.array().exp()));
}

Eigen::VectorXd dlog_vector(const Eigen::VectorXd & x)
{
  // dfun = function(x) exp(x)/(1+exp(x))^2 
  // pfun = function(x) 1/(1+exp(-x)) 
  // ddfun = function(eta) dfun(eta)*(1 - 2*pfun(eta))
  double EPSILON = 2.221e-16;
  double THRESH = 30.0; // , MTHRESH = -30.0; 
  // y = extreme value or exp(x)
  VectorXb mask_extreme = (x.array().abs() > THRESH);
  VectorXd y = mask_extreme.select(EPSILON, x.array().exp());
  for(unsigned int i = 0; i < y.size(); i++) {
    if(!mask_extreme(i)) {
      y(i) = y(i) * (1.0 - y(i)) / pow(y(i) + 1.0, 3);
    } else if(x[i] > THRESH) {
      y(i) *= -1;
    }
  }
  /* y.array() = (x.array() < MTHRESH) */
  /*   select(EPSILON, (x.array() > THRESH). */
  /*       select(-1*EPSILON, */ 
  /*         x.array().exp())); */
  /* y.array() = (x.array() < MTHRESH). */
  /*   select(x, (x.array() > THRESH). */
  /*       select(x, */ 
  /*         (y.array() * (1 - y.array())) / (y.array() + 1).pow(3))); */
  return(y);
}

Eigen::MatrixXd dlog_matrix(const Eigen::MatrixXd & x)
{
  MatrixXd y(x.rows(), x.cols());
  VectorXd ycol(x.rows());
  for(unsigned int i = 0; i < x.cols(); i++) {
    ycol = dlog_vector(x.col(i));
    y.col(i).array() = ycol.array();
  }
  return(y);
}

MatrixXd orth_matrix(const Eigen::MatrixXd & X0, const MatrixXb & Mask)
{

  // check parameters
  if(Mask.cols() != 1) {
    throw std::runtime_error("Mask must have 1 column");
  }
  if(X0.rows() != Mask.rows()) {
    throw std::runtime_error("number of rows different for X0 and Mask");
  }

  // parameters
  double tol_eval = 1e-15; // default in Regenie
                           
  double Neff = Mask.array().cast<double>().sum();

  // step 0. Copy 
  /* MatrixXd X = Mask.select(X0, 0.0); */
  MatrixXd X = X0;
  /* MatrixXd X(nr, nc) */
  /* for(i = 0; i < nc; i++) { */
  /*   X.col(i) = Mask.select(X0.col(i), 0.0); */
  /* } */

  /* // step 1. center (no intercept column in X) */
  ArrayXd means = X.colwise().sum() / Neff;
  X.array().rowwise() -= means.transpose();
  /* ArrayXd means = Mask.select(0.0, X).colwise().sum() / Neff; */
  /* X.array().rowwise() -= means.transpose(); */
  // respect the missingness pattern (Mask)
  /* for(i = 0; i < nc; i++) { */
  /*   X.col(i) = Mask.select(X.col(i), 0.0); */
  /* } */

  /* // step 2. Orthogonalize using SVD */
  MatrixXd XtX = X.transpose() * X;
  SelfAdjointEigenSolver<MatrixXd> es(XtX);
  VectorXd D = es.eigenvalues(); // eigenvalues sorted in increasing order
  MatrixXd V = es.eigenvectors();

  double max_eval = D.tail(1)(0);
  unsigned int nonzero_eval = (D.array() > max_eval * tol_eval).count();
  if(nonzero_eval == 0) {
    throw std::runtime_error("no columns left after thr. EVD");
  }
  ArrayXd sds = D.tail(nonzero_eval).array().sqrt();

  X *= V.rightCols(nonzero_eval);
  X.array().rowwise() /= sds.transpose();

  return X;
}

// print sum_stats
std::string MultiPhen::print_sumstats( int const& isnp, uint32_t const& snp_index,
    string const& test_string, string const& wgr_string, 
    variant_block* block_info, vector<snp> const& snpinfo, 
    struct param const* params)
{
  std::ostringstream buffer_header;
  buffer_header << snpinfo[snp_index].chrom << " " << snpinfo[snp_index].physpos << " "
    << snpinfo[snp_index].ID << " "
    << snpinfo[snp_index].allele1 << " " << snpinfo[snp_index].allele2 << " "
    << block_info->mac1 << " " << block_info->af1 << " " << block_info->ns1;
  string header = buffer_header.str();

  std::ostringstream buffer;
  buffer << header; // write header to buffer

  // output p-values
  if(pval_test == -1.0) {
    buffer << " NA";
  } else {
    buffer << " " << -log10(pval_test);
  }
 
  // model fits results
  buffer << " " << (response == "binom" ? 0 : 1);
  buffer << " " << it << " " << cnt_updates;

  bool firth = (response == "binom" ? firth_binom : false);
  buffer << " " << (firth ? 1 : 0);

  buffer << endl;

  return buffer.str();
}

std::string MultiPhen::print_sum_stats_htp(const variant_block* block_info, struct param const* params)
{
  std::ostringstream buffer;

  bool test_pass = (pval_test != -1.0);

  // bhat & 95 CI
  buffer << "NA\tNA\tNA\t";
  // Pvalue
  if(test_pass) buffer << pval_test << "\t";
  else buffer << "NA\t";

  // print out AF, counts in cases, counts in controls
  unsigned int ph = 0;
  buffer << block_info->af(ph) << "\t";
  // N, N_Ref, N_Het, N_Alt
  buffer << (int) block_info->genocounts.block(0,ph,3,1).sum() << "\t" 
    << (int) block_info->genocounts(0,ph) << "\t" 
    << (int) block_info->genocounts(1,ph) << "\t" 
    << (int) block_info->genocounts(2,ph) << "\t";
  buffer << "NA\tNA\tNA\tNA\t";
  /* double N_total = block_info->genocounts.block(0,ph,3,1).sum(), */
  /*        N_ref = block_info->genocounts(0,ph), */
  /*        N_het = block_info->genocounts(1,ph), */
  /*        N_alt = block_info->genocounts(2,ph); */
  /* double AAF = (N_het + 2*N_alt) / (2*N_total); */

  // info column
  if(test_pass) buffer << "DF=" << Ny;
  else buffer << "DF=0";
  buffer << ";NO_BETA";
  buffer << ";IT=" << it << ";UP=" << cnt_updates;
  buffer << ";BIN=" << (int)(response == "binom");
  if(bhat_y.size() > 0) {
    buffer << ";BHAT=";
    for(unsigned int i = 0; i < bhat_y.size(); i++) {
      if(i) buffer << "," << bhat_y(i);
      else buffer << bhat_y(i);
    }
  }

  buffer << endl;

  return buffer.str();
}

