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
void exp_vector(Eigen::VectorXd &);


//-------------------
// Class MultiPhen
//-------------------

void MultiPhen::setup_defaults()
{
  // settings
  verbose = 0;
  response = "unknown";
  optim = "WeightHalving";
  firth_binom = false;
  firth_mult = 1.0;
  reuse_start = false;
  approx_offset = false;
  maxit = 150; maxit2 = 10; 
  tol = 1e-4;
  check_step = true; max_step = 10.0;
  // statuses
  set_x = false;
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
  fit.firth_binom = firth_binom; // Firth correction
  fit.firth_mult = firth_mult; 
      
  fit.maxit = maxit; fit.maxit2 = maxit2; 
  fit.tol = tol;

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

  return(fit);
}

void MultiPhen::run(const Eigen::VectorXd & g, 
  const Eigen::MatrixXd& XY, unsigned int n_cov, unsigned int n_phen)
{
  // check if XY is set up
  if(!set_x) throw std::runtime_error("run: set_x is false");
  // set y
  setup_y(g); // -> Ym, yb
  // test
  if(test == "none") {
    reset_model();
    // do nothing
  } else if(test == "cov_score_it1") {
    maxit = 1; optim = "FisherScoring";
    run_test_score(XY, true); // inc_cov = false
  } else if(test == "nocov_score") {
    run_test_score(XY, false); // inc_cov = false
  } else if(test == "cov_score") {
    run_test_score(XY, true); // inc_cov = true
  } else if(test == "nocov_lrt") {
    run_test_lrt(XY, false); // inc_cov = false
  } else if(test == "cov_lrt") {
    run_test_lrt(XY, true); // inc_cov = true
  } else if(test == "nocov_score_addcov") {
    run_test_addcov(XY);
  } else {
    throw std::runtime_error("run: unknown test");
  }
}

void MultiPhen::run0(const Eigen::VectorXi & g, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, bool score_lrt)
{
  // set up Ordinal model (no Firth)
  Ordinal ord;
  ord.optim = optim; ord.tol = tol; ord.maxit = maxit; ord.maxit2 = maxit2;
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

// XY = [Intercept, X, Y, Inercept]
FitOrdinal MultiPhen::fit(const Eigen::Ref<const Eigen::MatrixXd> & XY, bool inc_cov, bool inc_phen)
{
  // initialize defaults settings 
  bool inc_phen_null = false, inc_phen_firth = inc_phen;
  bool use_offset = (inc_phen && approx_offset);
  bool copy_start = (reuse_start && inc_cov && inc_phen && !approx_offset);
  // update settings for Binom: no firth / firth
  if(response == "binom") {
    inc_phen_null = firth_binom && !inc_phen && !approx_offset;
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

  // store offset?
  if(!inc_phen && approx_offset) fit.store_offset = true;
  // apply offset?
  if(use_offset) {
    if(response == "binom") fit.setup_offset_binom(yo, false); // decrement_Nb = false
    else if (response == "multinom") fit.setup_offset_multinom_pom(yo, yo_int);
    else throw std::runtime_error("unknown response");
  }

  // do model fitting & control the columns in XY passed
  if(response == "binom") {
    if(use_offset) { 
      if(inc_phen_firth) fit.fit_binom(Mask, Ym, XY.rightCols(Ny1).leftCols(Ny)); // matrix of phenotypes Y 
      else throw std::runtime_error("use offset for the null model");
    } else { 
      if(inc_cov) {
        if(inc_phen_firth) fit.fit_binom(Mask, Ym, XY.leftCols(Nx1 + Ny)); // X + Y + Intercept
        else fit.fit_binom(Mask, Ym, XY.leftCols(Nx1)); // X + Intercept
      } else {
        if(inc_phen_firth) fit.fit_binom(Mask, Ym, XY.rightCols(Ny1)); // matrix of phenotypes Y + Intercept
        else fit.fit_binom(Mask, Ym, XY.leftCols(1)); // Intercept
      }
    }
  } else if(response == "multinom") {
    if(use_offset) {
      if(inc_phen) fit.fit_multinom_pom(Mask, Ym, XY.rightCols(Ny1).leftCols(Ny)); // matrix of phenotypes Y 
      else throw std::runtime_error("use offset for the null model");
    } else {
      if(inc_cov) {
        if(inc_phen) fit.fit_multinom_pom(Mask, Ym, XY.leftCols(Nx1 + Ny).rightCols(Nx + Ny)); // X + Y
        else fit.fit_multinom_pom(Mask, Ym, XY.leftCols(Nx1).rightCols(Nx)); // X
      } else {
        if(inc_phen) fit.fit_multinom_pom(Mask, Ym, XY.rightCols(Ny1).leftCols(Ny)); // matrix of phenotypes Y 
        else fit.fit_multinom_pom(Mask, Ym, XY.leftCols(0)); // 0 columns
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

void MultiPhen::run_test_addcov(const Eigen::Ref<const Eigen::MatrixXd> & XY)
{
  // run score test
  run_test_score(XY, false); // inc_cov = false


  if(pval_test < pval_thr) {
    run_test_lrt(XY, true); // inc_cov = true
  }
}

void MultiPhen::run_test_lrt(const Eigen::Ref<const Eigen::MatrixXd> & XY, bool inc_cov)
{
  reset_model(); // reset MultiPhen model fit results
  executed = true; 

  if(reuse_start) {
    if(!inc_cov) throw std::runtime_error("reuse_start in not available for inc_cov = false");
    if(approx_offset) throw std::runtime_error("reuse_start is not compatible with approx_offset");

    // null model: logit(g) = X alpha 
    FitOrdinal null = fit(XY, inc_cov, false); // inc_cov, inc_phen = false
    if(!null.converged) return;

    b0 = null.bhat;

    // full model: logit(g) = X alpha + Y beta
    FitOrdinal full = fit(XY, inc_cov, true); // inc_cov, inc_phen = true
    if(!full.converged) return;

    converged = true;
    boost::math::chi_squared dist(Ny);
    double stat_lrt = 2 * (full.loglik - null.loglik);
    pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt));
  } else if(approx_offset && response == "binom") {
    // null model: logit(g) = X alpha 
    FitOrdinal null = fit(XY, inc_cov, false); // inc_cov, inc_phen = false
    if(!null.converged) return;

    // store offset vector
    yo = null.yo;

    // full model: logit(g) = X alpha + Y beta
    FitOrdinal full = fit(XY, inc_cov, true); // inc_cov, inc_phen = true
    if(!full.converged) return;
    converged = true;


    // problem: null.mub is not at scale [0, 1]
    /* cout << "null.mub = " << null.mub.head(5).transpose() << endl; */
    VectorXd mub = null.yo;
    exp_vector(mub); // mub <- exp(mub)
    mub.array() /= (1.0 + mub.array()); // mub <- exp(mub) / (1 + exp(mub))
                                        //
    double ll_null = 0.0; //null.loglik_binom(Mask, Ym);
    ll_null += Ym.col(0).select((1.0 - mub.array()).log(), 0.0).array().sum(); // controls
    ll_null += Ym.col(1).select(mub.array().log(), 0.0).array().sum(); // cases
    if(null.firth_binom) {
      /* MatrixXd Y = XY.rightCols(Ny1).leftCols(Ny); // phenotypes Y */
      MatrixXd null_Info = XY.rightCols(Ny1).leftCols(Ny).transpose() * (XY.rightCols(Ny1).leftCols(Ny).array().colwise() * null.wb.array()).matrix();
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
    FitOrdinal null = fit(XY, inc_cov, false); // inc_cov, inc_phen = false
    if(!null.converged) return;

    // store offset vectors
    yo = null.yo;
    yo_int = null.yo_int;

    // full model: logit(g) = X alpha + Y beta
    FitOrdinal full = fit(XY, inc_cov, true); // inc_cov, inc_phen = true
    if(!full.converged) return;
    converged = true;

    boost::math::chi_squared dist(Ny);
    double stat_lrt = 2 * (full.loglik - null.loglik);
    /* cout << "stat_lrt = " << stat_lrt << " full.loglik = " << full.loglik << " null.loglik = " << null.loglik << endl; */
    pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt));
  } else {
    // null model: logit(g) = X alpha 
    FitOrdinal null = fit(XY, inc_cov, false); // inc_cov, inc_phen = false
    if(!null.converged) return;

    // full model: logit(g) = X alpha + Y beta
    FitOrdinal full = fit(XY, inc_cov, true); // inc_cov, inc_phen = true
    if(!full.converged) return;

    converged = true;
    boost::math::chi_squared dist(Ny);
    double stat_lrt = 2 * (full.loglik - null.loglik);
    pval_test = (stat_lrt < 0) ? 1 : boost::math::cdf(boost::math::complement(dist, stat_lrt));
  }
}

void MultiPhen::run_test_score(const Eigen::Ref<const Eigen::MatrixXd> & XY, bool inc_cov)
{
  bool _firth_binom = firth_binom, _approx_offset = approx_offset;
  firth_binom = false; approx_offset = false;
 
  reset_model(); // reset model fit results
  executed = true; 

  FitOrdinal null = fit(XY, inc_cov, false); // inc_cov, inc_phen = false
  if(!null.converged) { return; }

  converged = true; 
  if(trace) { cnt_updates += null.cnt_updates; it += null.it; }
  pval_test = test_score(null, Mask, Ym, yb, XY, inc_cov); 

  firth_binom = _firth_binom;
  approx_offset = _approx_offset;
}

void MultiPhen::setup_x(const VectorXb & _Mask,  const Eigen::MatrixXd& XY, unsigned int n_cov, unsigned int n_phen, 
    bool _pos_intercept_first, bool _pos_phen_first)
{
  // check
  if(XY.cols() != 2 + n_cov + n_phen) throw std::runtime_error("setup_x: dimensions XY");
  if(XY.rows() != _Mask.size()) throw std::runtime_error("setup_x: dimensions XY and Mask");
  // extract dimensions from XY
  N = XY.rows();
  /* Ncov = n_cov; // Nb = ncat1 + Ncov, where ncat1 depend on g */
  Nx = n_cov; Nx1 = n_cov + 1; Ny = n_phen; Ny1 = n_phen + 1;
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
  for(i = 0; i < g.size(); i++) {
    genotypes.insert(g[i]);
  }
  // check genotypes levels: 0/1 or 0/1/2
  if(!(genotypes.size() == 2 || genotypes.size() == 3)) throw std::runtime_error("setup_y: number of genotype categories must be 2 or 3");

  // assign ncat, ncat1
  ncat = genotypes.size();
  ncat1 = ncat - 1;
  ncat1sq = ncat1 * ncat1;
  
  // assign response
  if(ncat == 2) response = "binom";
  else if(ncat == 3) response = "multinom";
  else throw std::runtime_error("setup_y: unexpected number of genotype categories");

  // assign Ncov, Nb
  /* Nb = ncat1 + Ncov; */

  // assign Ym
  Ym.resize(N, ncat);
  Ncat = VectorXi::Constant(ncat, 0);
  // loop over a a few genotype categories
  for(i = 0, it_set = genotypes.begin(); i < ncat; i++, it_set++) {
    Ym.col(i) = Mask.select(g.array() == *it_set, false);
    Ncat(i) = Ym.col(i).cast<int>().sum();
  }
  // assign yb if binomial
  if(response == "binom") {
    yb = Ym.col(1).cast<double>(); // booleans -> 0/1
  }
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
  ord.optim = optim; ord.tol = tol; ord.maxit = maxit;
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
  ord.optim = optim; ord.tol = tol; ord.maxit = maxit;
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
      << " apply_start = " << apply_start << " apply_offset = " << apply_offset << " store_offset = " << store_offset << " exclude_intercepts = " << exclude_intercepts 
      << endl;
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
  firth_binom = false;

  maxit = 100; 
  maxit2 = 7; it2 = 0;
  tol = 1e-4;

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
  fit.firth_binom = firth_binom; // Firth correction
      
  fit.maxit = maxit; fit.maxit2 = maxit2; 
  fit.tol = tol;

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
    const Eigen::Ref<const Eigen::MatrixXd> & XY, bool inc_cov)
{
  double pval;
  if(response == "multinom") {
    if(inc_cov) pval = test_score_multinom_pom(null, Mask, Ym, XY.leftCols(Nx1).rightCols(Nx), XY.rightCols(Ny1).leftCols(Ny)); // covariates X (no intercept); phenotypes Y
    else pval = test_score_multinom_pom(null, Mask, Ym, XY.leftCols(0), XY.rightCols(Ny1).leftCols(Ny)); // 0 covarites (no intercept); phenotypes Y
  } else if(response == "binom") {
    if(inc_cov) pval = test_score_binom(null, Mask, yb, XY.leftCols(Nx1).rightCols(Nx), XY.rightCols(Ny1).leftCols(Ny)); // covariates X (no intercept); phenotypes Y
    else pval = test_score_binom(null, Mask, yb, XY.leftCols(0), XY.rightCols(Ny1).leftCols(Ny)); // 0 covarites (no intercept); phenotypes Y
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
    if(Ncov0 > Ncov) {
      throw std::runtime_error("Ncov0 >= Ncov");
    }
    if(preproc_cov) {
      throw std::runtime_error("preproc_cov is on when Ncov0 != 0");
    }

    if(response == "multinom") {
      Ncov1 = Ncov - Ncov0;
    } else if(response == "binom") {
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
    yo = Xb0; // overwrite offset vector yo if present
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
}

void FitOrdinal::update_par_multinom(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X, const Eigen::VectorXd & b)
{
  if(verbose >= 2) cout << " - update_par_multinom\n"; 
  /* if(verbose >= 3) cout << "   b = " << b.transpose() << endl; */

  unsigned int i, k,  l, m, start;

  VectorXd b_cov;
  if(exclude_intercepts) b_cov = b;
  else b_cov = b.segment(ncat1, Ncov);

  if(Ncov) Xb0 = X * b_cov;
  else Xb0.setZero(); 
  if(apply_offset) Xb0.array() += yo.array(); // offset 

  // update linear predictor Xb with intercepts
  if(exclude_intercepts) {
    for(i = 0; i < ncat1; i++) {
      Xb.col(i).array() = Xb0.array(); 
    }
  } else {
    VectorXd b_int = b.head(ncat1);
    for(i = 0; i < ncat1; i++) {
      Xb.col(i).array() = Xb0.array() + b_int(i);
    }
  }
  if(apply_offset) {
    for(i = 0; i < ncat1; i++) {
      Xb.col(i).array() += yo_int(i);
    }
  }
  exp_eta = Xb; exp_matrix(exp_eta);
  gamma.array() = exp_eta.array() / (1.0 + exp_eta.array());

  P = gamma;
  for(i = 1; i < ncat1; i++) {
    P.col(i).array() -= gamma.col(i - 1).array();
  }
  Psum = P.rowwise().sum();
  
  bool error_psum = (Psum.array() >= 1.0).any();
  if(error_psum) {
    throw std::runtime_error("some elements in Psum >= 1.0");
  }
  Pk.array() = 1.0 - Psum.array();
 
  // D = (Y[, -ncat] / P) - (Y[, ncat] / Pk)
  for(i = 0; i < ncat1; i++) {
    D.col(i).array() = Y.col(i).cast<double>().array() / P.col(i).array() - 
      Y.col(ncat1).cast<double>().array() / Pk.array();
  } 

  // Q = dh / deta
  PQ.array() = gamma.array() * (1.0 - gamma.array());
  for(m = 0; m < ncat1; m++) {
    l = m;
    start = l * ncat1;
    Q.col(start + m).array() = PQ.col(m).array();
  }
  for(m = 1; m < ncat1; m++) {
    l = m - 1;
    start = l * ncat1;
    Q.col(start + m).array() = -1.0 * PQ.col(l).array();
  }

  // V
  for(k = 0; k < ncat1; k++) {
    // R code:
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
  for(unsigned int i = 0; i < V.cols(); i++) {
    V.col(i) = Mask.select(V.col(i), 0.0);
  }
  for(unsigned int i = 0; i < W.cols(); i++) {
    W.col(i) = Mask.select(W.col(i), 0.0);
  }

  // Score
  // Score = c(colSums(V), colSums((crossprod(V, X))))
  if(!exclude_intercepts) cur_Score.head(ncat1) = V.colwise().sum();
  if(Ncov) cur_Score.tail(Ncov) = (V.transpose() * X).colwise().sum();

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
  LLT<MatrixXd> llt(cur_Info);
  cur_v = llt.solve(cur_Score);

  cur_loglik = loglik_multinom(Mask, Y);
  cur_dev = -2.0 * cur_loglik;
}

double FitOrdinal::loglik_multinom(const VectorXb & Mask, const MatrixXb & Y)
{
  // loglik = colSums(Y*log(cbind(P, Pk))) %>% sum
  double res = 0.0;
  // 1, 2, ..., ncat1 categories
  for(unsigned int i = 0; i < ncat1; i++) {
    res += Y.col(i).select(P.col(i).array().log(), 0.0).array().sum();
  }
  // the last ncat category
  res += Y.col(ncat1).select(Pk.array().log(), 0.0).array().sum();

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
    yo = Xb; // overwrite offset vector yo if present
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
}

void FitOrdinal::update_par_binom(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X, const Eigen::VectorXd & b)
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

}

void FitOrdinal::update_par_binom_firth(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X, const Eigen::VectorXd & b)
{
  if(verbose >= 2) cout << " - update_par_binom_firth\n"; 

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

  // update Score = Ab + Sb
  /* Ab = crossprod(X, h * (0.5 - mu)) */
  /* Sb = crossprod(X, y - mu) */
  // solve: v = solve(Info, Score)

  if(Ncov0) {
    if(last0) {
      LLT<MatrixXd> llt1(cur_Info.block(0, 0, Ncov1, Ncov1));
      cur_Score.head(Ncov1) = (X.leftCols(Ncov1).transpose() * Mask.select((Y.col(1).cast<double>().array() - mub.array() + firth_mult * h.array() * (0.5 - mub.array())), 0.0).matrix()).array();
      cur_v.head(Ncov1) = llt1.solve(cur_Score.head(Ncov1));
    } else {
      LLT<MatrixXd> llt1(cur_Info.block(Ncov0, Ncov0, Ncov1, Ncov1));
      cur_Score.tail(Ncov1) = (X.rightCols(Ncov1).transpose() * Mask.select((Y.col(1).cast<double>().array() - mub.array() + firth_mult * h.array() * (0.5 - mub.array())), 0.0).matrix()).array();
      cur_v.tail(Ncov1) = llt1.solve(cur_Score.tail(Ncov1));
    }
  } else {
    cur_Score = (X.transpose() * Mask.select((Y.col(1).cast<double>().array() - mub.array() + firth_mult * h.array() * (0.5 - mub.array())), 0.0).matrix()).array();
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
    file << b.size() << " " << cnt_updates << " " << cur_loglik << " " << cur_dev
      << " " << diff; 
    for(unsigned int i = 0; i < b.size(); i++) {
      file << " " << b(i);
    }
    file << endl;
    file.close();
  }
}

double FitOrdinal::loglik_binom(const VectorXb & Mask, const MatrixXb & Y)
{
  // loglik = sum(Y*log(p) + (1-Y)*log(1-p)), where p = mu
  double res = 0.0;
  // contols
  res += Y.col(0).select((1.0 - mub.array()).log(), 0.0).array().sum();
  // cases
  res += Y.col(1).select(mub.array().log(), 0.0).array().sum();

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
  if(verbose >= 2) cout << "stop_criterion\n";

  bool stop;
  // v1: abs. diff. between bhat_cur and bhat_prev
  // stop = (cur_v.norm() < tol);
  // v2: abs. max. value of Scores. 
  // - Example: Regenie Firth model fitting 
  stop = (cur_Score.array().abs().maxCoeff() < tol); 
  /* stop = (cur_Score.segment(0, Nb - Ncov0).array().abs().maxCoeff() < tol); */ 
  // v3: relative diff. in deviance. 
  // - Example: glm2::glm.fit2.R
  /* stop = ((abs(cur_dev - prev_dev) / (1.0 + abs(prev_dev))) < tol); // 1.0 to prevent from division by zero */

  return(stop);
}

void FitOrdinal::update_par(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X, const Eigen::VectorXd & b)
{
  // count
  if(trace) { cnt_updates++; }

  // store prev. values 
  prev_dev = cur_dev;

  // update directions
  if(response == "multinom") {
    update_par_multinom(Mask, Y, X, b);
  } else if(response == "binom") {
    if(firth_binom) {
      update_par_binom_firth(Mask, Y, X, b);
    } else {
      update_par_binom(Mask, Y, X, b);
    }
  } else {
    throw std::runtime_error("unknown response");
  }
}

// optimization loop for binom. response
bool FitOrdinal::optimize(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X)
{
  if(verbose >= 2) cout << "optimize\n"; 

  bool res;

  if(optim == "FisherScoring") {
    res = optimize_FisherScoring(Mask, Y, X);
  } else if(optim == "WeightHalving") {
    res = optimize_WeightHalving(Mask, Y, X);
  } else {
    throw std::runtime_error("unknown optimize");
  }

  if(verbose >= 2) cout << "optimize it = " << it << " | cnt_updates = " << cnt_updates << endl;
  if(verbose >= 3) cout << "bhat = " << cur_b.transpose() << endl;

  return(res);
}

bool FitOrdinal::optimize_FisherScoring(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X)
{
  unsigned int i;
  bool res, stop;
  double ratio_step;

  cur_b = b0;
  for(i = 0; i < maxit; i++) {
    // update directions
    update_par(Mask, Y, X, cur_b);

    // check the stopping criteria
    stop = stop_criterion();
    if(stop) {
      break;
    }

    // check the  absolute step size to be less than max_step for each entry of step (cur_v2)
    if(check_step) {
      ratio_step = cur_v.array().abs().maxCoeff() / max_step;
      if(ratio_step > 1.0) {
        cur_v.array() /= ratio_step;
      }
    }

    // update parameters
    cur_b += cur_v;

    // check if bhat is nan
    if(cur_b.array().isNaN().any()) {
      it = i;
      return false;
    }

  }
  // assign # iterations performed
  it = i;
  res = (i < maxit);

  return(res);
}

bool FitOrdinal::optimize_WeightHalving(const VectorXb & Mask, const MatrixXb & Y, const Eigen::Ref<const Eigen::MatrixXd> & X)
{
  if(verbose >= 2) cout << "optimize_WeightHalving\n";

  unsigned int i, i2;
  bool stop, res;

  // declare variables
  VectorXd cur_b2, cur_v2;
  double cur_dev2;
  double denom, ratio_step;

  // initial values for Outer loop
  cur_b = b0;
  update_par(Mask, Y, X, cur_b); // get (i) the step size (cur_v); (ii) current value of dev

  // Outer loop
  for(i = 1; i < maxit; i++) {
    if(verbose >= 3) cout << " - outer loop it " << i << endl;

    // stopping criteria is checked here in the beginning rather than in the end of loop
    // - reason: update_par was called above
    // - example: i = 0 & no covariates & initial values are proportions --> no optimization inside the Outer/Inner loop is required
    stop = stop_criterion();
    if(stop) {
      break;
    }

    // initial values for Inner loop
    cur_b2 = cur_b; cur_v2 = cur_v;
    cur_dev2 = cur_dev;

    denom = 2.0;
    
    // Inner loop (step halving)
    for(i2 = 0; i2 < maxit2; i2++) {
      if(verbose >= 3) cout << " - inner loop it " << i2 << endl;

      // update step according to the rule: step(i) = step(initial) / 2^i
      // one exception from the rule: skip halving at the very first iteration, i2 = 0 --> Fisher Scoring at i2 = 0
      if(i2) { 
        cur_v2.array() /= denom;
      }

      // check the  absolute step size to be less than max_step for each entry of step (cur_v2)
      if(check_step) {
        ratio_step = cur_v2.array().abs().maxCoeff() / max_step;
        if(ratio_step > 1.0) {
          cur_v2.array() /= ratio_step;
        }
      }

      // update param.
      // - the baseline value (cur_b2) is fixed
      // - the increment step (cur_v2) is reduced at each iteration (see the code line above)
      cur_b = cur_b2 + cur_v2;

      // check if bhat is nan
      if(cur_b.array().isNaN().any()) {
        it = i;
        return false;
      }

      // update Score, Info, loglik, dev
      update_par(Mask, Y, X, cur_b);

      // stop the inner loop (step halving) if dev. is improved
      if(cur_dev < cur_dev2) {
        break;
      }
    }

    // assign Innter loop iterations
    it2 += i2;

    bool strict_WeightHalving = false;
    if(strict_WeightHalving) {
      // check if all Inner loop iterations are used & exit
      if(i2 == maxit2) {
        // let the first iteration (i = 0) go even when convergence failure
        if(i) {
          // assign & exit
          it = i; 
          return false;
        }
      }
    }
  }

  // assign # iterations performed
  it = i;
  res = (i < maxit);

  return res;
}

//------------------------------
// Utils
//------------------------------

void exp_matrix(Eigen::MatrixXd & X)
{
  // See: mu = binomial()$linkinv(eta)
  // - https://github.com/wch/r-source/blob/trunk/src/library/stats/src/family.c
  // - https://stackoverflow.com/a/1566222
  double EPSILON = 2.221e-16;
  double THRESH = 30.0, MTHRESH = -30.0; 
  double INVEPS = 1.0/EPSILON;
  for(unsigned int i = 0; i < X.cols(); i++) {
    X.col(i).array() = (X.col(i).array() < MTHRESH).
      select(EPSILON, (X.col(i).array() > THRESH).
          select(INVEPS, X.col(i).array().exp()));
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

  buffer << endl;

  return buffer.str();
}
