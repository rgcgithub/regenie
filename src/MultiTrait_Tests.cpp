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
#include "SKAT.hpp" // get_lambdas
#include "Joint_Tests.hpp" // get_acat
#include "MCC.hpp"

#include "MultiTrait_Tests.hpp"

using namespace std;
using namespace Eigen;


double get_fisher_robust(const Eigen::Ref<const ArrayXd>& logp)
{
  double stat = 2.0 * log(10) * logp.sum();
  // get logp from chisq(k)
  double logp_fisher;
  get_logp(logp_fisher, stat, 2 * logp.size());

  return(logp_fisher);
}

//------------------------
// Class MTestsResults
//------------------------

MTestsResults::MTestsResults(unsigned int n_tests_, unsigned int q_, unsigned int M_) 
{
  setup(n_tests_, q_, M_);
}

void MTestsResults::setup(unsigned int n_tests_, unsigned int q_, unsigned int M_) 
{
  // assign
  n_tests = n_tests_;
  q = q_;
  M = M_;
  // allocate
  stats_mt.resize(n_tests);
  logp_mt.resize(n_tests);
  logp_univar.resize(n_tests);
}

MTestsResults::~MTestsResults() { }

//---------------
// Class MTests
//---------------

MTests::MTests() 
{
  setup();
}

MTests::~MTests() { }

void MTests::setup() 
{
  verbose = 0;
  precomp = false;

  mcc_skew_abs = 1.0;
  mcc_z2 = 4; // ~ p-value = 0.05

  n_tests = MULTITRAIT_N_TESTS;
  n_traits = 0;
  Neff0 = 0.0;

  // bayes test parameters
  prior_a0 = 6.0;
  prior_Q0 = 4.0;
  prior_Mbeta0 = 0.0;
  prior_Vbeta0 = 0.02;

  // NNLS
  nnls = NNLS();
}

void MTests::check_setup_no_data() 
{
  if(n_tests != MULTITRAIT_N_TESTS) { throw std::runtime_error("check_setup_no_data: n_tests != MULTITRAIT_N_TESTS"); }
}

void MTests::check_setup_data() 
{
  if(n_traits == 0) { throw std::runtime_error("check_setup_data: n_traits == 0"); }
  if(Mask.cols() == 0) { throw std::runtime_error("check_setup_data: Mask.cols() == 0"); }
  if(Mask.cols() != n_traits) { throw std::runtime_error("check_setup_data: Mask.cols() != n_traits"); }
  if(Neff.size() == 0) { throw std::runtime_error("check_setup_data: Neff.size() == 0"); }
  if(Neff.size() != n_traits) { throw std::runtime_error("check_setup_data: Neff.size() != n_traits"); }
}

void MTests::setup_masks(const MatrixXb &masked_indivs)
{
  // check dimensions
  if(masked_indivs.rows() == 0) { throw std::runtime_error("setup_masks: masked_indivs.rows() == 0"); }
  if(masked_indivs.cols() == 0) { throw std::runtime_error("setup_masks: masked_indivs.cols() == 0"); }

  Mask = masked_indivs;
  Neff = Mask.cast<double>().colwise().sum(); 

  Mask0 = Mask.col(0);
  for(unsigned int i = 1; i < Mask.cols(); i++) {
    Mask0.col(0).array() = Mask0.col(0).array() || Mask.col(i).array();
  }
  Neff0 = Mask0.cast<double>().sum();
}

void MTests::setup_yres(const MatrixXd &res)
{
  // check dimensions
  if(res.rows() == 0) { throw std::runtime_error("setup_yres: res.rows() == 0"); }
  if(res.cols() == 0) { throw std::runtime_error("setup_yres: res.cols() == 0"); }
  // check masks were setup upstream
  if(Neff0 == 0.0) { throw std::runtime_error("setup_yres: Neff0 == 0"); }

  // set up #traits
  n_traits = res.cols();

  // Y: n x q matrix of traits
  Yres = res;
  Y0res = res;
  Y0res.array().colwise() *= Mask0.col(0).array().cast<double>().array();

  // cross-product YtY
  precomp0_YtY = MatrixXd(n_traits, n_traits).setZero().selfadjointView<Lower>().rankUpdate(Y0res.adjoint());
  // Syy = covariance of Y
  precomp0_Syy = precomp0_YtY / (Neff0 - 1.0);
  precomp0_Syy_inv = precomp0_Syy.llt().solve(MatrixXd::Identity(n_traits, n_traits));
  // MANOVA-specific 
  precomp0_ld0= precomp0_YtY.ldlt().vectorD().array().log().sum();
  // Bayes-specific 
  VectorXd Mbeta_0 = VectorXd::Constant(n_traits, prior_Mbeta0);
  MatrixXd Q0 = MatrixXd::Constant(n_traits, n_traits, 0.0);
  Q0.diagonal().array() = prior_Q0;
  double ld = (Q0 + precomp0_YtY).ldlt().vectorD().array().log().sum();
  precomp0_LL_M0 =  0.5* (double)(n_traits) * log(prior_Vbeta0) - 0.5*(Neff0 + prior_a0 + (double)(n_traits) - 1.0) * ld;
  // NNLS
  nnls.ss_weights(precomp0_Syy);
  // Hiearachical Omnibus
  VectorXd lambdas = VectorXd::Zero(n_traits);
  get_lambdas(lambdas, precomp0_Syy, 1e-5);
  precomp0_lambdas_Syy = lambdas;
  // normalize eigen values: l_i = l_i / sum(l_i)
  lambdas /= lambdas.sum();
  precomp0_lambdas_norm_Syy = lambdas;
  // Robust Omnibus
  SelfAdjointEigenSolver<MatrixXd> es(precomp0_Syy);
  PC_Y0res = Y0res * es.eigenvectors();
  PC_Y0res.array().colwise() *= Mask0.col(0).array().cast<double>().array();
  // skewness of PCs
  compute_skew_pc();
  // RINT PCs
  RPC_Y0res = PC_Y0res;
  for(unsigned int i = 0; i < RPC_Y0res.cols(); i++) {
    MatrixXd rpc_col = RPC_Y0res.col(i);
    rint_pheno(rpc_col, Mask0);
    RPC_Y0res.col(i) = rpc_col;
  }
}

//----------------------------------
// Main function to apply tests
//----------------------------------

MTestsResults MTests::run_tests_snp(
    int const& isnp, struct geno_block& gblock,
    const Ref<const MatrixXd>& yres, const Ref<const RowVectorXd>& p_sd_yres, 
    struct param const& params)
{
  check_setup_data();

  MapMatXd Gmat(gblock.Gmat.col(isnp).data(), params.n_samples, 1);

  MTestsResults mt_results(n_tests, n_traits, Gmat.cols());

  assoc_manova(yres, Gmat, mt_results);
  assoc_omnibus0(yres, Gmat, mt_results);
  assoc_bayes(yres, Gmat, mt_results);

  // debug
  if(verbose > 2) {
    /* dump(yres, Gmat); */
    const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
    // dump y
    ofstream file_y("mt.Y.txt");
    file_y << yres.format(CSVFormat);
    file_y.close();
    // dump X
    ofstream file_X("mt.X.txt");
    file_X << Gmat.format(CSVFormat);
    file_X.close();
    // dump masks M
    ofstream file_M("mt.M.txt");
    file_M << Mask.array().cast<double>().matrix().format(CSVFormat);
    file_M.close();
    // dump masks M0
    ofstream file_M0("mt.M0.txt");
    file_M0 << Mask0.array().cast<double>().matrix().format(CSVFormat);
    file_M0.close();
  }

  return(mt_results);
}

MTestsResults MTests::run_tests_snp_precomp(
    int const& isnp, struct geno_block& gblock,
    struct param const& params)
{
  check_setup_data();

  MapMatXd Gmat(gblock.Gmat.col(isnp).data(), params.n_samples, 1);

  MTestsResults mt_results(n_tests, n_traits, Gmat.cols());

  assoc_manova(Y0res, Gmat, mt_results);
  assoc_omnibus0(Y0res, Gmat, mt_results);
  assoc_cpc0(PC_Y0res, Gmat, mt_results);
  assoc_rcpc0(RPC_Y0res, Gmat, mt_results);
  assoc_bayes(Y0res, Gmat, mt_results);
  assoc_nnls0(Y0res, Gmat, mt_results);
  
  // debug
  if(verbose > 2) {
    dump_data(Y0res, Gmat, Mask, Mask0);
  }

  return(mt_results);
}

//----------------------------------
// MANOVA test
//----------------------------------

// Association scan based on MANOVA
void MTests::assoc_manova(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G, MTestsResults& mt_results)
{
  unsigned pos_test = 0;
  unsigned pos_test_npmanova = 5;

  // check dimensions
  if(Y.rows() == 0) { throw std::runtime_error("MANOVA: Y.rows() == 0"); }
  if(Y.cols() == 0) { throw std::runtime_error("MANOVA: Y.cols() == 0"); }
  if(G.rows() == 0) { throw std::runtime_error("MANOVA: G.rows() == 0"); }
  if(G.cols() == 0) { throw std::runtime_error("MANOVA: G.cols() == 0"); }

  // dimensions
  unsigned int N = Y.rows();
  unsigned int q = Y.cols();
  unsigned int M = G.cols();

  double N_data = Neff0;

  // check dimensions
  if(Y.cols() != n_traits) { throw std::runtime_error("MANOVA: Y.cols() != n_traits"); }
  if(G.rows() != N) { throw std::runtime_error("MANOVA: #rows in G != N"); }
  if(M < 1) { throw std::runtime_error("MANOVA: M < 1"); }
  if(q < 1) { throw std::runtime_error("MANOVA: q < 1"); }

  // pre-compute matrix products
  MatrixXd YtY(q, q);
  if(precomp) YtY = precomp0_YtY;
  else YtY = MatrixXd(q, q).setZero().selfadjointView<Lower>().rankUpdate(Y.adjoint());

  VectorXd G2 = G.colwise().squaredNorm();

  // Estimates of beta
  // Bhat, q x M matrix
  MatrixXd Bhat = (Y.transpose() * G).
    array().rowwise() / G2.array().transpose();

  // loop over M variants
  // test statistic for variant i: (q/2 - N + 1) * log(Wi)
  // where Wi = det(E1) / det(E0) is a Wilkstest statistic for MANOVA
  mt_results.stats_mt[pos_test].resize(M);
  mt_results.logp_mt[pos_test].resize(M);
  mt_results.stats_mt[pos_test_npmanova].resize(M);
  mt_results.logp_mt[pos_test_npmanova].resize(M);

  double ld0;
  if(precomp) ld0 = precomp0_ld0;
  else ld0 = YtY.ldlt().vectorD().array().log().sum();
  
  MatrixXd E, zzt;
  VectorXd b;
  double ld1, stat_i, pval_i;
  boost::math::chi_squared dist(q);
  for(unsigned int i = 0; i < M; i++) {
    // 1. MANOVA
    b = Bhat.col(i);
    zzt = b * b.transpose() * G2[i];
    E = YtY - zzt;
    ld1 = E.ldlt().vectorD().array().log().sum();
    // error: ld1 > ld0
    if(ld1 > ld0) { throw std::runtime_error("MANOVA: log(det(E1)) > log(det(E0))"); }
    // MANOVA test statistics & p-value
    stat_i = ((double)(q) / 2 - N_data + 1.0) * (ld1 - ld0);
    pval_i = boost::math::cdf(boost::math::complement(dist, stat_i));
    // store test results
    mt_results.stats_mt[pos_test][i] = stat_i;
    mt_results.logp_mt[pos_test][i] = -log10(pval_i);

    // 2. NPMANOVA
    // NPMANOVA test statistic (pseudo F statistic)
    unsigned int p0 = 0; // no covariates projected out
    double mean_SS_E = E.diagonal().sum() / (N_data - (double)(p0) - 1.0);
    double mean_SS_T = zzt.diagonal().sum(); // df_T = 1, i.e., one SNP is tested (M1)
    stat_i = mean_SS_T / mean_SS_E; // F tilde in notation of PERMANOVA
    // EVD on YtY
    VectorXd lambdas = VectorXd::Zero(q);
    if(precomp) {
      lambdas = precomp0_lambdas_norm_Syy;
    } else {
      MatrixXd Sy(q, q);
      Sy = YtY / (Neff0 - 1.0);
      VectorXd lambdas = VectorXd::Zero(q);
      get_lambdas(lambdas, Sy, 1e-5);
      // normalize eigen values: l_i = l_i / sum(l_i)
      lambdas /= lambdas.sum();
    }
    // NPMANOVA p-value
    // v1
    /* pval_i =  get_chisq_mix_pv(stat_i, lambdas); */
    // v2
    // re-scale so that max lambda is 1 (lambda is sorted)
    double newQ = stat_i / lambdas.tail(1)(0);
    VectorXd newL = lambdas / lambdas.tail(1)(0);
    pval_i = get_kuonen_pv(newQ, newL); // SPA
    if(pval_i <= 0) {// if SPA failed
       pval_i = get_liu_pv(newQ, newL); // only use mod Liu if Davies/SPA failed
    }
    // store test results
    mt_results.stats_mt[pos_test_npmanova][i] = stat_i;
    mt_results.logp_mt[pos_test_npmanova][i] = -log10(pval_i);
  }
  
} 

//----------------------------------------
// NNLS0 test (complete sample overlap)
//----------------------------------------

// Association NNLS0 scan usign summary stats from NNLS0
void MTests::assoc_nnls0(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G, MTestsResults& mt_results)
{
  unsigned pos_test = 3;

  // check bhat
  if(mt_results.zscore_univar.size() == 0) { throw std::runtime_error("assoc_nnls0: bhat"); }

  // dimensions
  unsigned int q = Y.cols();
  unsigned int M = G.cols();

  // pre-compute matrix products
  MatrixXd Sy(q,q);
  if(precomp) {
    Sy = precomp0_Syy; 
  } else {
    Sy = MatrixXd(q, q).setZero().selfadjointView<Lower>().rankUpdate(Y.adjoint());
    Sy /= (Neff0 - 1.0);
  }

  MatrixXd Sy_inv;
  if(precomp) Sy_inv = precomp0_Syy_inv;
  else Sy_inv = Sy.llt().solve(MatrixXd::Identity(q, q));

  // loop over M variants
  // test statistic for variant i: z' S^{-1} z
  mt_results.stats_mt[pos_test].resize(M);
  mt_results.logp_mt[pos_test].resize(M);

  VectorXd z(q);
  NNLS nnls0 = nnls;
  for(unsigned int i = 0; i < M; i++) {
    for(unsigned int j = 0; j < q; j++) {
      z[j] = mt_results.zscore_univar[i][j]; // dim 1 = variants; dim 2 = traits
    }
    
    nnls0.ss_run(z);

    mt_results.stats_mt[pos_test][i] = 0;
    mt_results.logp_mt[pos_test][i] = -log10(nnls0.pval_min2);
  }
} 

//----------------------------------------
// Omnibus0 test (complete sample overlap)
//----------------------------------------

// Association scan based on Omnibus + complete sample overlap
void MTests::assoc_omnibus0(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G, MTestsResults& mt_results)
{
  unsigned pos_test = 1;
  unsigned pos_test_sumz = 4;
  unsigned pos_test_homnibus = 6;

  // check dimensions
  if(Y.rows() == 0) { throw std::runtime_error("Omnibus0: Y.rows() == 0"); }
  if(Y.cols() == 0) { throw std::runtime_error("Omnibus0: Y.cols() == 0"); }
  if(G.rows() == 0) { throw std::runtime_error("Omnibus0: G.rows() == 0"); }
  if(G.cols() == 0) { throw std::runtime_error("Omnibus0: G.cols() == 0"); }

  // dimensions
  unsigned int N = Y.rows();
  unsigned int q = Y.cols();
  unsigned int M = G.cols();

  double N_data = Neff0;

  // check dimensions
  if(Y.cols() != n_traits) { throw std::runtime_error("Omnibus0: Y.cols() != n_traits"); }
  if(G.rows() != N) { throw std::runtime_error("Omnibus0: #rows in G != N"); }
  if(M < 1) { throw std::runtime_error("Omnibus0: M < 1"); }
  if(q < 1) { throw std::runtime_error("Omnibus0: q < 1"); }

  // pre-compute matrix products
  MatrixXd Sy(q,q);
  if(precomp) {
    Sy = precomp0_Syy; 
  } else {
    Sy = MatrixXd(q, q).setZero().selfadjointView<Lower>().rankUpdate(Y.adjoint());
    Sy /= (Neff0 - 1.0);
  }

  MatrixXd Sy_inv;
  if(precomp) Sy_inv = precomp0_Syy_inv;
  else Sy_inv = Sy.llt().solve(MatrixXd::Identity(q, q));
  
  VectorXd lambdas = VectorXd::Zero(q);
  if(precomp) {
      lambdas = precomp0_lambdas_Syy;
  } else {
    VectorXd lambdas = VectorXd::Zero(q);
    get_lambdas(lambdas, Sy, 1e-5);
  }

  VectorXd G2 = G.colwise().squaredNorm();

  // Marginal Z-scores, q x M matrix
  MatrixXd Z(q, M); //, B(q, M);
  VectorXd bhat(M), s2(M);
  for(unsigned i = 0; i < q; i++) {
    bhat = (Y.col(i).transpose() * G).array().rowwise() / G2.array().transpose();
    /* B.row(i) = bhat; */
    // residuals, s2
    s2 = (((G.array().rowwise() * bhat.array().transpose()). // predicted yp = X bhat
      colwise() - Y.col(i).array()). // residuals = y - yp
      matrix().colwise().squaredNorm()). // residuals^2
      array() / (N_data - 1.0); // s2 = residuals^2 / (N - 1)
    Z.row(i) = bhat.array() * (G2.array() / s2.array()).sqrt();
  }

  // loop over M variants
  // test statistic for variant i: z' S^{-1} z
  mt_results.logp_univar[pos_test].resize(M);
  mt_results.zscore_univar.resize(M);
  mt_results.stats_mt[pos_test].resize(M);
  mt_results.logp_mt[pos_test].resize(M);
  mt_results.stats_mt[pos_test_sumz].resize(M);
  mt_results.logp_mt[pos_test_sumz].resize(M);
  mt_results.stats_mt[pos_test_homnibus].resize(M);
  mt_results.logp_mt[pos_test_homnibus].resize(M);

  VectorXd z(q); //, b(q);
  double stat_univar, pval_univar, stat_i, pval_i;
  boost::math::chi_squared dist(q);
  boost::math::chi_squared dist_univar(1);
  for(unsigned int i = 0; i < M; i++) {
    z = Z.col(i);
    /* b = B.col(i); */

    // univar. tests (single traits)
    mt_results.logp_univar[pos_test][i].resize(q); // dim 1 = mt; dim 2 = variants; dim 3 = traits
    mt_results.zscore_univar[i].resize(q); // dim 1 = variants; dim 2 = traits
    for(unsigned int j = 0; j < q; j++) {
      stat_univar = z[j] * z[j];
      pval_univar = boost::math::cdf(boost::math::complement(dist_univar, stat_univar));
      mt_results.logp_univar[pos_test][i][j] = -log10(pval_univar);
      mt_results.zscore_univar[i][j] = z[j];

    }
    // multi-trait test: Omnibus 
    stat_i = z.transpose() * Sy_inv * z;
    /* double stat_i = b.transpose() * Sy_inv * b; */
    pval_i = boost::math::cdf(boost::math::complement(dist, stat_i));

    mt_results.stats_mt[pos_test][i] = stat_i;
    mt_results.logp_mt[pos_test][i] = -log10(pval_i);
    
    // multi-trait test: SumZ with T = sum(Z)^2 / sum(V) 
    stat_i = z.sum();
    stat_i = stat_i*stat_i / Sy.sum();
    pval_i = boost::math::cdf(boost::math::complement(dist_univar, stat_univar));

    mt_results.stats_mt[pos_test_sumz][i] = stat_i;
    mt_results.logp_mt[pos_test_sumz][i] = -log10(pval_i);

    // hOmnibus
    stat_i = z.transpose() * z;
    double newQ = stat_i / lambdas.tail(1)(0);
    VectorXd newL = lambdas / lambdas.tail(1)(0);
    pval_i = get_kuonen_pv(newQ, newL); // SPA
    if(pval_i <= 0) {// if SPA failed
       pval_i = get_liu_pv(newQ, newL); // only use mod Liu if Davies/SPA failed
    }
    mt_results.stats_mt[pos_test_homnibus][i] = pval_i;
    mt_results.logp_mt[pos_test_homnibus][i] = -log10(pval_i);
  }
} 

//----------------------------------------
// CPC0 test (0 = complete sample overlap)
//----------------------------------------

// Association scan based on PCs + complete sample overlap
void MTests::assoc_cpc0(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G, MTestsResults& mt_results)
{
  unsigned pos_test = 7;
  unsigned pos_test_acpc_sumchi2 = 11, pos_test_acpc_fisher = 12, pos_test_acpc_acat = 13;

  // check dimensions
  if(Y.rows() == 0) { throw std::runtime_error("CPC0: Y.rows() == 0"); }
  if(Y.cols() == 0) { throw std::runtime_error("CPC0: Y.cols() == 0"); }
  if(G.rows() == 0) { throw std::runtime_error("CPC0: G.rows() == 0"); }
  if(G.cols() == 0) { throw std::runtime_error("CPC0: G.cols() == 0"); }

  // dimensions
  unsigned int N = Y.rows();
  unsigned int q = Y.cols();
  unsigned int M = G.cols();

  double N_data = Neff0;

  // check dimensions
  if(Y.cols() != n_traits) { throw std::runtime_error("Omnibus0: Y.cols() != n_traits"); }
  if(G.rows() != N) { throw std::runtime_error("Omnibus0: #rows in G != N"); }
  if(M < 1) { throw std::runtime_error("Omnibus0: M < 1"); }
  if(q < 1) { throw std::runtime_error("Omnibus0: q < 1"); }

  // pre-compute matrix products
  VectorXd G2 = G.colwise().squaredNorm();

  // Marginal Z-scores, q x M matrix
  MatrixXd Z(q, M); //, B(q, M);
  VectorXd bhat(M), s2(M);
  for(unsigned i = 0; i < q; i++) {
    bhat = (Y.col(i).transpose() * G).array().rowwise() / G2.array().transpose();
    /* B.row(i) = bhat; */
    // residuals, s2
    s2 = (((G.array().rowwise() * bhat.array().transpose()). // predicted yp = X bhat
      colwise() - Y.col(i).array()). // residuals = y - yp
      matrix().colwise().squaredNorm()). // residuals^2
      array() / (N_data - 1.0); // s2 = residuals^2 / (N - 1)
    Z.row(i) = bhat.array() * (G2.array() / s2.array()).sqrt();
  }

  // loop over M variants
  // test statistic for variant i: stat(CPC) = Sum z^2
  mt_results.stats_mt[pos_test].resize(M);
  mt_results.logp_mt[pos_test].resize(M);
  mt_results.logp_mt[pos_test_acpc_sumchi2].resize(M);
  mt_results.logp_mt[pos_test_acpc_fisher].resize(M);
  mt_results.logp_mt[pos_test_acpc_acat].resize(M);
  mt_results.zscore_cpc.resize(M);
  mt_results.zscore_acpc.resize(M);

  VectorXd z(q); //, b(q);
  double stat_i, pval_i, logp_i;
  boost::math::chi_squared dist(q), dist_univar(1);
  VectorXd logp_univar(q);
  MCC mcc;
  boost::math::chi_squared chisq(1);
  for(unsigned int i = 0; i < M; i++) {
    z = Z.col(i);

    // univar. tests (single traits)
    for(unsigned int j = 0; j < q; j++) {
      logp_univar[j] = -log10(boost::math::cdf(boost::math::complement(dist_univar, z[j]*z[j])));
    }

    // store z-scores
    mt_results.zscore_cpc[i].resize(q); // dim 1 = variants; dim 2 = traits
    Map<VectorXd>(mt_results.zscore_cpc[i].data(), q) = z;

    stat_i = z.transpose() * z;
    pval_i = boost::math::cdf(boost::math::complement(dist, stat_i));

    // store CPC results
    mt_results.stats_mt[pos_test][i] = stat_i;
    logp_i = -log10(pval_i);
    mt_results.logp_mt[pos_test][i] = logp_i;

    // ACPC
    mt_results.zscore_acpc[i].resize(q); // dim 1 = variants; dim 2 = traits
    Map<VectorXd>(mt_results.zscore_acpc[i].data(), q) = z;
    // adjust CPC?
    bool mcc_failed = false;
    double z2j, z2j_adj;
    if(n_skewed_pc) {
      if((z.array().square() > mcc_z2).any()) {
        for(unsigned int j = 0; j < q; j++) {
          z2j = z[j]*z[j];
          if((skew_PC[j] > mcc_skew_abs) && z2j > mcc_z2) { 
            // adjust
            mcc.setup_y(Mask0, PC_Y0res.col(j), 1); // ncov analyzed = 1
            MCCResults mcc_results_i = mcc.run(G.col(i));
            if(mcc_results_i.Skip(0, 0)) {
              mcc_failed = true;
              break;
            }
            logp_univar[j] = -log10(mcc_results_i.Pval(0, 0));
            // adjust Z-score
            z2j_adj = boost::math::quantile(boost::math::complement(chisq, mcc_results_i.Pval(0, 0)));
            mt_results.zscore_acpc[i][j] *= sqrt(z2j_adj / z2j);
          } else { 
            // don't adjust & PC results
            get_logp(logp_univar[j], z2j);
          }
        }
      }
    }
    // check for MCC failure & combine p-values & store
    if(mcc_failed) {
      // !NB! how to encode NA p-values?
      mt_results.logp_mt[pos_test_acpc_sumchi2][i] = -9;
      mt_results.logp_mt[pos_test_acpc_fisher][i] = -9;
      mt_results.logp_mt[pos_test_acpc_acat][i] = -9;
    } else {
      // sum chi2
      stat_i = Map<ArrayXd>(mt_results.zscore_acpc[i].data(), q).square().sum();
      get_logp(logp_i, stat_i, q);
      mt_results.logp_mt[pos_test_acpc_sumchi2][i] = logp_i;
      // Fisher
      logp_i = get_fisher_robust(logp_univar);
      mt_results.logp_mt[pos_test_acpc_fisher][i] = logp_i;
      // ACAT
      logp_i = get_acat(logp_univar);
      mt_results.logp_mt[pos_test_acpc_acat][i] = logp_i;
    }
  }
} 

//------------------------------------------------
// Robust CPC0 test (0 = complete sample overlap)
//------------------------------------------------

// Association scan based on Robust PCs + complete sample overlap
void MTests::assoc_rcpc0(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G, MTestsResults& mt_results)
{
  unsigned pos_test_sumchi2 = 8, pos_test_fisher = 9, pos_test_acat = 10;

  // check dimensions
  if(Y.rows() == 0) { throw std::runtime_error("RCPC0: Y.rows() == 0"); }
  if(Y.cols() == 0) { throw std::runtime_error("RCPC0: Y.cols() == 0"); }
  if(G.rows() == 0) { throw std::runtime_error("RCPC0: G.rows() == 0"); }
  if(G.cols() == 0) { throw std::runtime_error("RCPC0: G.cols() == 0"); }

  // dimensions
  unsigned int N = Y.rows();
  unsigned int q = Y.cols();
  unsigned int M = G.cols();

  double N_data = Neff0;

  // check dimensions
  if(Y.cols() != n_traits) { throw std::runtime_error("Omnibus0: Y.cols() != n_traits"); }
  if(G.rows() != N) { throw std::runtime_error("Omnibus0: #rows in G != N"); }
  if(M < 1) { throw std::runtime_error("Omnibus0: M < 1"); }
  if(q < 1) { throw std::runtime_error("Omnibus0: q < 1"); }

  // pre-compute matrix products
  VectorXd G2 = G.colwise().squaredNorm();

  // Marginal Z-scores, q x M matrix
  MatrixXd Z(q, M); //, B(q, M);
  VectorXd bhat(M), s2(M);
  for(unsigned i = 0; i < q; i++) {
    bhat = (Y.col(i).transpose() * G).array().rowwise() / G2.array().transpose();
    // residuals, s2
    s2 = (((G.array().rowwise() * bhat.array().transpose()). // predicted yp = X bhat
      colwise() - Y.col(i).array()). // residuals = y - yp
      matrix().colwise().squaredNorm()). // residuals^2
      array() / (N_data - 1.0); // s2 = residuals^2 / (N - 1)
    Z.row(i) = bhat.array() * (G2.array() / s2.array()).sqrt();
  }

  // loop over M variants
  // test statistic for variant i: P-value(RCPC) = Fisher(P-values of PC1, PC2, ...)
  /* mt_results.stats_mt[pos_test].resize(M); */
  mt_results.logp_mt[pos_test_sumchi2].resize(M);
  mt_results.logp_mt[pos_test_fisher].resize(M);
  mt_results.logp_mt[pos_test_acat].resize(M);
  mt_results.zscore_rcpc.resize(M);

  VectorXd z(q); 
  double stat_i, logp_i;
  boost::math::chi_squared dist(q), dist_univar(1);
  VectorXd logp_univar(q);
  for(unsigned int i = 0; i < M; i++) {
    z = Z.col(i);
    
    mt_results.zscore_rcpc[i].resize(q); // dim 1 = variants; dim 2 = traits
    Map<VectorXd>(mt_results.zscore_rcpc[i].data(), q) = z;

    // univar. tests (single traits)
    for(unsigned int j = 0; j < q; j++) {
      logp_univar[j] = -log10(boost::math::cdf(boost::math::complement(dist_univar, z[j]*z[j])));
    }

    // Sum Chi2
    stat_i = Map<ArrayXd>(mt_results.zscore_rcpc[i].data(), q).square().sum();
    get_logp(logp_i, stat_i, q);
    mt_results.logp_mt[pos_test_sumchi2][i] = logp_i;
    // Fisher
    logp_i = get_fisher_robust(logp_univar);
    mt_results.logp_mt[pos_test_fisher][i] = logp_i;
    // ACAT
    logp_i = get_acat(logp_univar);
    mt_results.logp_mt[pos_test_acat][i] = logp_i;
  }
} 

//----------------------------------------
// Bayesian test (BF instead of P-value)
//----------------------------------------

// Association scan using BF
void MTests::assoc_bayes(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G, MTestsResults& mt_results)
{
  unsigned pos_test = 2;

  // check dimensions
  if(Y.rows() == 0) { throw std::runtime_error("Bayes: Y.rows() == 0"); }
  if(Y.cols() == 0) { throw std::runtime_error("Bayes: Y.cols() == 0"); }
  if(G.rows() == 0) { throw std::runtime_error("Bayes: G.rows() == 0"); }
  if(G.cols() == 0) { throw std::runtime_error("Bayes: G.cols() == 0"); }

  // dimensions
  unsigned int N = Y.rows();
  unsigned int q = Y.cols();
  unsigned int M = G.cols();

  double N_data = Neff0;

  // check dimensions
  if(Y.cols() != n_traits) { throw std::runtime_error("MANOVA: Y.cols() != n_traits"); }
  if(G.rows() != N) { throw std::runtime_error("MANOVA: #rows in G != N"); }
  if(M < 1) { throw std::runtime_error("MANOVA: M < 1"); }
  if(q < 1) { throw std::runtime_error("MANOVA: q < 1"); }

  // variables
  double ld; // log-determinant of matrices
  double LL_M0; // log-lik. for M0 

  // pre-compute matrix products
  MatrixXd YtY(q, q);
  if(precomp) YtY = precomp0_YtY;
  else YtY = MatrixXd(q, q).setZero().selfadjointView<Lower>().rankUpdate(Y.adjoint());

  VectorXd G2 = G.colwise().squaredNorm();

  // prior p(Sigma) = IW(a0, Q0)
  VectorXd Mbeta_0 = VectorXd::Constant(q, prior_Mbeta0);
  MatrixXd Q0 = MatrixXd::Constant(q, q, 0.0);
  Q0.diagonal().array() = prior_Q0;

  //--------------------
  // null model M0
  //--------------------
  if(precomp) {
    LL_M0 = precomp0_LL_M0;
  } else  {
    ld = (Q0 + YtY).ldlt().vectorD().array().log().sum();
    LL_M0 =  0.5* (double)(q) * log(prior_Vbeta0) - 0.5*(N_data + prior_a0 + (double)(q) - 1.0) * ld;
  }
  
  //--------------------
  // full model M1
  //--------------------

  // posterior parameters
  // p(beta - Mbeta_1 | Sigma) = N(Vbeta_1, Sigma)
  // p(Sigma) = IW(a1, Q1)
  
  // Vector Vbeta_1 of size 1 x M
  VectorXd Vbeta_1 = (G2.array() + 1.0/prior_Vbeta0).inverse();

  // Matrix Mbeta_1 of size q x M
  // 1. Initialize Mbeta_1
  MatrixXd Mbeta_1 = (((Y.transpose() * G).
    // 2. Column-wise update: Mbeta_1 = apply(Mbeta_1, 2, function(x) (x + (Mbeta_0 / Vbeta_0)))
    array().colwise() + (Mbeta_0.array() / prior_Vbeta0)).
    // 3. Row-wise update: Mbeta_1 = apply(Mbeta_1, 1, function(x) (x / (1/Vbeta_0 + scprod_G))) %>% t
    rowwise() * Vbeta_1.transpose().array()).matrix();

  /* double a1 = prior_a0 + (double)(N_data); */

  MatrixXd Q1_common = (Q0.array() + (Mbeta_0.squaredNorm() / prior_Vbeta0) + YtY.array()).matrix();
  MatrixXd Q1 = MatrixXd(q, q);

  // Loop over M variants
  mt_results.stats_mt[pos_test].resize(M);
  mt_results.logp_mt[pos_test].resize(M);

  /* //  0.5*q * log(Vbeta_1) - 0.5*(N + a0 + q - 1) * as.numeric(determinant(Q1, log = TRUE)$modulus) */
  VectorXd LL_M1_base = Vbeta_1.array().log() * (0.5*(double)(q));
  for(unsigned int i = 0; i < M; i++) {
    Q1 = (Q1_common.array() - (Mbeta_1.col(i).squaredNorm() / Vbeta_1[i])).matrix();

    ld = Q1.ldlt().vectorD().array().log().sum();
    double LL_M1 = LL_M1_base[i] + -0.5 * ((double)(N_data) + prior_a0 + (double)(q) - 1.0) * ld;
    double BF = exp(LL_M1 - LL_M0);

    mt_results.logp_mt[pos_test][i] = log10(BF);
  }

} 

//----------------------------------
// Estimate correlations
//----------------------------------

double cor2_mask(const Eigen::ArrayXd &x, const Eigen::ArrayXd &y, 
    const MatrixXb &mask)
{
  size_t n = x.rows();
  size_t n_val = mask.cast<int>().sum();
  if(n_val <= 1) {
    throw std::runtime_error("0 or 1 samples overlapped for a pair of traits");
  }
  if((size_t)(n_val) > n) {
    throw std::runtime_error("unexpected value of n_val (> n)");
  }

  // means for x and y
  double mx = mask.select(x, 0.0).sum() / (double)(n_val);
  double my = mask.select(y, 0.0).sum() / (double)(n_val);

  // Formula for the Pearson correlation, rho(x, y)
  // rho(x, y) = [sum(x_i y_i) - n m_x m_y] / 
  //               [ sqrt{ sum(x_i^2) - n m_x^2 } sqrt{ sum(y_i^2) - n m_y^2 } ] 
  double r = 
      // r_xy = [sum(x_i y_i) - n m_x m_y] / 
      (mask.select(x * y, 0.0).sum() - (double)(n_val)*mx*my) /
        // [ sqrt{ sum(x_i^2) - n m_x^2 } sqrt{ sum(y_i^2) - n m_y^2 } ] 
        sqrt(
          (mask.select(x, 0.0).square().sum() - (double)(n_val)*mx*mx) * 
          (mask.select(y, 0.0).square().sum() - (double)(n_val)*my*my));

  return(r);
}

void MTests::compute_cory(const Ref<const MatrixXd>& Y, const Ref<const MatrixXb>& M)
{
  size_t p = Y.cols();
  if(p == 0) {
    throw std::runtime_error("#columns of Y (p) is zero");
  }

  Ryy = MatrixXd::Constant(p, p, 0.0);
  Ryy.diagonal().array() = 1.0;

  for(size_t i = 0; i < p; i++) {
    for(size_t j = i + 1; j < p; j++) {
      MatrixXb mask = M.col(i) || M.col(j);

      double r = cor2_mask(Y.col(i), Y.col(j), mask);
      Ryy(i, j) = r;
      Ryy(j, i) = r;
    }
  }
}
  
void MTests::compute_ncory(const Ref<const MatrixXb>& M)
{
  size_t p = M.cols();
  if(p == 0) {
    throw std::runtime_error("#columns of M (p) is zero");
  }

  N_Ryy = MatrixXi::Constant(p, p, 0);
  N_Ryy.diagonal().array() = M.cast<int>().colwise().sum(); 

  for(size_t i = 0; i < p; i++) {
    for(size_t j = i + 1; j < p; j++) {
      MatrixXb mask = M.col(i) || M.col(j);

      size_t n_val = mask.cast<int>().sum();
      N_Ryy(i, j) = n_val;
      N_Ryy(j, i) = n_val;
    }
  }
}

void MTests::compute_skew_pc()
{
  if(PC_Y0res.cols() != n_traits) { throw std::runtime_error("MTests::compute_skew_rcpc: PC_Y0res"); }

  // for each RPC, compute skewness
  skew_PC = ArrayXd::Constant(n_traits, 0.0);
  for(unsigned int i = 0; i < n_traits; i++) {
    double skew = skew_pheno(PC_Y0res.col(i), Mask0);
    // debug
    /* cout << "skew = " << skew << endl; */
    skew_PC[i] = abs(skew);
  }
  n_skewed_pc = (skew_PC.array() > mcc_skew_abs).sum();
}

//----------------------------------
// Print results
//----------------------------------

string MTests::print_sumstats(
    const MTestsResults& mt_results,
    int const& isnp, uint32_t const& snp_index,
    string const& test_string, string const& wgr_string, 
    variant_block* block_info, vector<snp> const& snpinfo, 
    struct param const* params)
{
  check_setup_data();

  // check attributes
  if(mt_results.M == 0) { throw std::runtime_error("MTests::print_sumstats: M == 0"); }

  /* string header; */
  std::ostringstream buffer_header;
  buffer_header << snpinfo[snp_index].chrom << " " << snpinfo[snp_index].physpos 
    << " "<< snpinfo[snp_index].ID << " "
    << snpinfo[snp_index].allele1 << " "<< snpinfo[snp_index].allele2 << " " 
    << block_info->mac(0) << " "
    << block_info->af(0) << " " 
    << Neff0;                                
  string header = buffer_header.str();

  std::ostringstream buffer;
  buffer << header; // write header to buffer

  // output p-values
  double max_logp;
  size_t pos_ominibus0 = 1;
  if(mt_results.logp_univar[pos_ominibus0].size() != 1) { throw std::runtime_error("MTests::print_sumstats: logp_univar"); }
  // single-trait p-values
  /* for(size_t s = 0; s < n_traits; s++) { */
  /*   buffer << " " << mt_results.logp_univar[pos_ominibus0][0][s]; // write univar p-value */ 
  /* } */
  // minP / minQ
  max_logp = mt_results.logp_univar[pos_ominibus0][0][0];
  for(size_t s = 1; s < n_traits; s++) {
    double logp_s = mt_results.logp_univar[pos_ominibus0][0][s];
    if(logp_s > max_logp) max_logp = logp_s;
  }
  buffer << " " << max_logp;
  double max_logq = max(0.0, max_logp - log10(n_traits));
  buffer << " " << max_logq;
  // multi-trait test p-values
  for(size_t t = 0; t < n_tests; t++) {
    if(mt_results.logp_mt[t].size() != 1) { throw std::runtime_error("MTests::print_sumstats: logp_mt"); }

    buffer << " " << mt_results.logp_mt[t][0]; // write p-value of test #t to buffer
  }
  // NNLS0 q-value
  size_t pos_nnls0 = 3;
  max_logq = max(0.0, mt_results.logp_mt[pos_nnls0][0] - log10(2.0));
  buffer << " " << max_logq;
  // Single-trait z-scores
  for(size_t j = 0; j < n_traits; j++) {
    buffer << " " << mt_results.zscore_univar[0][j]; // dim 1 = variants; dim 2 = traits
  }
  // PC z-scores
  for(size_t j = 0; j < n_traits; j++) {
    buffer << " " << mt_results.zscore_cpc[0][j]; // dim 1 = variants; dim 2 = traits
  }
  // Robust PC z-scores
  for(size_t j = 0; j < n_traits; j++) {
    buffer << " " << mt_results.zscore_rcpc[0][j]; // dim 1 = variants; dim 2 = traits
  }
  // Adjusted PC z-scores
  for(size_t j = 0; j < n_traits; j++) {
    buffer << " " << mt_results.zscore_acpc[0][j]; // dim 1 = variants; dim 2 = traits
  }

  buffer << endl; // finish writing a line for SNP to buffer

  return(buffer.str());
}


//----------------------------------
// Debug
//----------------------------------

void MTests::dump_data(const MatrixXd &Y, const MatrixXd &X, const MatrixXb &Mask, const MatrixXb &Mask0)
{
  const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
  // dump Y
  ofstream file_y("mt.Y.txt");
  file_y << Y.format(CSVFormat);
  file_y.close();
  // dump X
  ofstream file_X("mt.X.txt");
  file_X << X.format(CSVFormat);
  file_X.close();
  // dump masks M
  ofstream file_M("mt.M.txt");
  file_M << Mask.array().cast<double>().matrix().format(CSVFormat);
  file_M.close();
  // dump masks M0
  ofstream file_M0("mt.M0.txt");
  file_M0 << Mask0.array().cast<double>().matrix().format(CSVFormat);
  file_M0.close();
}

