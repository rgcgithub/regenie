/* 

   This file is part of the regenie software package.

   Copyright (c) 2020-2023 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

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
#include "MultiTrait_Tests.hpp"

using namespace std;
using namespace Eigen;

MTests::MTests() 
{
  setup();
}

MTests::~MTests() { }

void MTests::setup() 
{
  n_tests = MULTITRAIT_N_TESTS;
  logp_tests.resize(n_tests);
  // bayes test parameters
  prior_a0 = 6.0;
  prior_Q0 = 4.0;
  prior_Mbeta0 = 0.0;
  prior_Vbeta0 = 0.02;
}

void MTests::check_setup() 
{
  n_tests = MULTITRAIT_N_TESTS;
}

//----------------------------------
// Main function to apply tests
//----------------------------------

void MTests::apply_tests(
    const int& chrom, const int& block, 
    struct phenodt const* pheno_data, 
    const Eigen::Ref<const Eigen::MatrixXd>& res, 
    struct geno_block const* gblock, std::vector<variant_block>& block_info, 
    vector<string> const& ynames, struct param const* params)
{
  /* cout << "nrow(Gmat) = " << gblock->Gmat.rows() << endl; */
  /* cout << "ncol(Gmat) = " << gblock->Gmat.cols() << endl; */
  assoc_manova(res, gblock->Gmat);
}

void MTests::apply_tests_snp(
    int const& isnp, struct geno_block& gblock,
    const Ref<const MatrixXd>& yres, const Ref<const RowVectorXd>& p_sd_yres, 
    struct param const& params)
{
  check_setup();

  MapMatXd Gmat(gblock.Gmat.col(isnp).data(), params.n_samples, 1);

  assoc_manova(yres, Gmat);
  assoc_omnibus0(yres, Gmat);
  assoc_bayes(yres, Gmat);
  /* cout << "mt.pvals = " << pvals << endl; */
}


//----------------------------------
// MANOVA test
//----------------------------------

// Association scan based on MANOVA
void MTests::assoc_manova(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G)
{
  unsigned pos_test = 0;

  // check dimensions
  if(Y.rows() == 0) { throw std::runtime_error("MANOVA: Y.rows() == 0"); }
  if(Y.cols() == 0) { throw std::runtime_error("MANOVA: Y.cols() == 0"); }
  if(G.rows() == 0) { throw std::runtime_error("MANOVA: G.rows() == 0"); }
  if(G.cols() == 0) { throw std::runtime_error("MANOVA: G.cols() == 0"); }

  // dimensions
  N = Y.rows();
  q = Y.cols();
  M = G.cols();

  double N_data = Neff(0);

  // check dimensions
  if(G.rows() != N) { throw std::runtime_error("MANOVA: #rows in G != N"); }
  if(M < 1) { throw std::runtime_error("MANOVA: M < 1"); }
  if(q < 1) { throw std::runtime_error("MANOVA: q < 1"); }

  // pre-compute matrix products
  MatrixXd YtY(MatrixXd(q, q).setZero().selfadjointView<Lower>().rankUpdate(Y.adjoint()));
  VectorXd G2 = G.colwise().squaredNorm();

  // Estimates of beta
  // Bhat, q x M matrix
  MatrixXd Bhat = (Y.transpose() * G).
    array().rowwise() / G2.array().transpose();

  // loop over M variants
  // test statistic for variant i: (q/2 - N + 1) * log(Wi)
  // where Wi = det(E1) / det(E0) is a Wilkstest statistic for MANOVA
  stats.resize(M);
  pvals.resize(M);
  logp_tests[pos_test].resize(M);

  double ld0 = YtY.ldlt().vectorD().array().log().sum();
  
  MatrixXd E;
  VectorXd b;
  double ld1, stat_i;
  boost::math::chi_squared dist(q);
  for(unsigned int i = 0; i < M; i++) {
    b = Bhat.col(i);
    E = YtY - (b * b.transpose()) * G2[i];

    ld1 = E.ldlt().vectorD().array().log().sum();
    // error: ld1 > ld0
    if(ld1 > ld0) { throw std::runtime_error("MANOVA: log(det(E1)) > log(det(E0))"); }

    stat_i = ((double)(q) / 2 - N_data + 1.0) * (ld1 - ld0);
    stats[i] = stat_i;
    pvals[i] = boost::math::cdf(boost::math::complement(dist, stat_i));

    logp_tests[pos_test][i] = -log10(pvals[i]);
  }
  
} 

//----------------------------------------
// Omnibus0 test (complete sample overlap)
//----------------------------------------

// Association scan based on Omnibus + complete sample overlap
void MTests::assoc_omnibus0(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G)
{
  unsigned pos_test = 1;

  // check dimensions
  if(Y.rows() == 0) { throw std::runtime_error("Omnibus0: Y.rows() == 0"); }
  if(Y.cols() == 0) { throw std::runtime_error("Omnibus0: Y.cols() == 0"); }
  if(G.rows() == 0) { throw std::runtime_error("Omnibus0: G.rows() == 0"); }
  if(G.cols() == 0) { throw std::runtime_error("Omnibus0: G.cols() == 0"); }

  // dimensions
  N = Y.rows();
  q = Y.cols();
  M = G.cols();

  double N_data = Neff(0);

  // check dimensions
  if(G.rows() != N) { throw std::runtime_error("Omnibus0: #rows in G != N"); }
  if(M < 1) { throw std::runtime_error("Omnibus0: M < 1"); }
  if(q < 1) { throw std::runtime_error("Omnibus0: q < 1"); }

  // pre-compute matrix products
  MatrixXd Sy(MatrixXd(q, q).setZero().selfadjointView<Lower>().rankUpdate(Y.adjoint()));
  Sy.array() /= (N_data - 1.0);
  MatrixXd Sy_inv = Sy.llt().solve(MatrixXd::Identity(q, q));

  VectorXd G2 = G.colwise().squaredNorm();

  // Marginal Z-scores, q x M matrix
  MatrixXd Z(q, M);
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
  // test statistic for variant i: z' S^{-1} z
  stats.resize(M);
  pvals.resize(M);
  logp_tests[pos_test].resize(M);

  VectorXd z(q);
  double stat_i;
  boost::math::chi_squared dist(q);
  for(unsigned int i = 0; i < M; i++) {
    z = Z.col(i);
    stat_i = z.transpose() * Sy_inv * z;

    stats[i] = stat_i;
    pvals[i] = boost::math::cdf(boost::math::complement(dist, stat_i));

    logp_tests[pos_test][i] = -log10(pvals[i]);
  }
} 

//----------------------------------------
// Bayesian test (BF instead of P-value)
//----------------------------------------

// Association scan using BF
void MTests::assoc_bayes(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G)
{
  unsigned pos_test = 2;

  // check dimensions
  if(Y.rows() == 0) { throw std::runtime_error("Bayes: Y.rows() == 0"); }
  if(Y.cols() == 0) { throw std::runtime_error("Bayes: Y.cols() == 0"); }
  if(G.rows() == 0) { throw std::runtime_error("Bayes: G.rows() == 0"); }
  if(G.cols() == 0) { throw std::runtime_error("Bayes: G.cols() == 0"); }

  // dimensions
  N = Y.rows();
  q = Y.cols();
  M = G.cols();

  double N_data = Neff(0);

  // check dimensions
  if(G.rows() != N) { throw std::runtime_error("Bayes: #rows in G != N"); }
  if(M < 1) { throw std::runtime_error("Bayes: M < 1"); }
  if(q < 1) { throw std::runtime_error("Bayes: q < 1"); }

  // variables
  double ld; // log-determinant of matrices
  double LL_M0; // log-lik. for M0 

  // pre-compute matrix products
  MatrixXd YtY(MatrixXd(q, q).setZero().selfadjointView<Lower>().rankUpdate(Y.adjoint()));
  VectorXd G2 = G.colwise().squaredNorm();

  // prior p(Sigma) = IW(a0, Q0)
  VectorXd Mbeta_0 = VectorXd::Constant(q, prior_Mbeta0);
  MatrixXd Q0 = MatrixXd::Constant(q, q, 0.0);
  Q0.diagonal().array() = prior_Q0;

  //--------------------
  // null model M0
  //--------------------
  ld = (Q0 + YtY).ldlt().vectorD().array().log().sum();
  LL_M0 =  0.5* (double)(q) * log(prior_Vbeta0) - 0.5*(N_data + prior_a0 + (double)(q) - 1.0) * ld;
  
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
  logp_tests[pos_test].resize(M);

  /* //  0.5*q * log(Vbeta_1) - 0.5*(N + a0 + q - 1) * as.numeric(determinant(Q1, log = TRUE)$modulus) */
  VectorXd LL_M1_base = Vbeta_1.array().log() * (0.5*(double)(q));
  for(unsigned int i = 0; i < M; i++) {
    Q1 = (Q1_common.array() - (Mbeta_1.col(i).squaredNorm() / Vbeta_1[i])).matrix();

    ld = Q1.ldlt().vectorD().array().log().sum();
    double LL_M1 = LL_M1_base[i] + -0.5 * ((double)(N_data) + prior_a0 + (double)(q) - 1.0) * ld;
    double BF = exp(LL_M1 - LL_M0);

    logp_tests[pos_test][i] = log10(BF);
  }

} 

//----------------------------------
// Estimate correlations
//----------------------------------

double cor2_mask(const Eigen::ArrayXd &x, const Eigen::ArrayXd &y, 
    const ArrayXb &mask)
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

  return r;
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
      ArrayXb mask = M.col(i) || M.col(j);

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
      ArrayXb mask = M.col(i) || M.col(j);

      size_t n_val = mask.cast<int>().sum();
      N_Ryy(i, j) = n_val;
      N_Ryy(j, i) = n_val;
    }
  }
}

//----------------------------------
// Print results
//----------------------------------

// print sum_stats
string MTests::print_sumstats(
    int const& isnp, uint32_t const& snp_index,
    string const& test_string, string const& wgr_string, 
    variant_block* block_info, vector<snp> const& snpinfo, 
    struct param const* params)
{
  check_setup();

  // check attributes
  if(M == 0) { throw std::runtime_error("MTests::print_sumstats: M == 0"); }

  /* string header; */
  /* if(!params->htp_out) { header = print_sum_stats_head(snp_index, snpinfo); } */
  std::ostringstream buffer_header;
  buffer_header<< snpinfo[snp_index].chrom << " " << snpinfo[snp_index].physpos 
    << " "<< snpinfo[snp_index].ID << " "
    << snpinfo[snp_index].allele1 << " "<< snpinfo[snp_index].allele2 << " " ;
  string header = buffer_header.str();

  std::ostringstream buffer;
  buffer << header; // write header to buffer

  // output p-values
  /* if(isnp >= pvals.size()) { throw std::runtime_error("MTests::print_sumstats: isnp"); } */
  if(pvals.size() != 1) { throw std::runtime_error("MTests::print_sumstats: isnp"); }

  for(size_t t = 0; t < n_tests; t++) {
    if(logp_tests[t].size() != 1) { throw std::runtime_error("MTests::print_sumstats: logp_tests"); }

    buffer << " " << logp_tests[t][0]; // write p-value of test #t to buffer
  }
  buffer << endl; // finish writing a line for SNP to buffer

  /* // Loop over M SNPs */
  /* std::ostringstream buffer; */
  /* for(size_t i = 0; i < M; i++) { // col 0 = chisq, col 1 = logp */
  /*   bool assoc_ok = (pvals[i] >= 0); */

  /*   if(assoc_ok) { */
  /*     buffer << header << " " << -log10(pvals[i]) << endl; */
  /*     /1* if(params->htp_out) *1/ */ 
  /*     /1* buffer << print_sum_stats_head_htp(snp_index, "trait_set", test_string + wgr_string, snpinfo, params) *1/ */ 
  /*     /1*   << print_sum_stats_htp(-1, -1, stats[i], pvals[i], -1, -1, -1, block_info->genocounts, snp_index, true, 1, params); *1/ */
  /*     /1* else *1/ */ 
  /*     /1*   buffer << header *1/ */ 
  /*     /1*     << print_sum_stats(-1,-1,-1, -1, block_info->ns1, test_string, -1, -1, stats[i], pvals[i], true, 1, params, (i+1)); *1/ */
  /*   } else { // print NA sum stats */
  /*     /1* buffer << print_na_sumstats(i, 1, header, test_string, block_info, *params); *1/ */
  /*   } */
  /* } */
  
  return buffer.str();
}

