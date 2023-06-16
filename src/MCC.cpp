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
/* #include "Files.hpp" */
/* #include "Geno.hpp" */
/* #include "Pheno.hpp" */

#include "MCC.hpp"
/* #include "pgamma/rbase_pgamma.hpp" */
#include <boost/math/distributions/gamma.hpp>

using namespace std;
using namespace Eigen;

// R: pgamma(q, shape = shape, scale = scale, lower = FALSE)
double boost_pgamma(double q, double shape, double scale, bool lower)
{
  double p;
  boost::math::gamma_distribution<> dist(shape, scale);

  if(q < 0) {
    p = 0.99999; // 1.0 - 1e-5;
  } else {
    if(lower) {
      p = boost::math::cdf(dist, q);
    } else {
      p = boost::math::cdf(boost::math::complement(dist, q));
    }
  }

  return(p);
}

//---------------------------------
// MCCResults Constructors
//---------------------------------

MCCResults::MCCResults(unsigned int N_, unsigned int K_, unsigned int q_, unsigned int M_, const VectorXd &n_)
{
  N = N_;
  K = K_;
  q = q_;
  M = M_;

  n = n_.array() - (double)(K);
  n2 = n.array().square();

  // allocate
  sum_x.resize(q, M); sum_x_sq.resize(q, M);
  sum_x2.resize(q, M); sum_x2_sq.resize(q, M);
  sum_x3.resize(q, M); sum_x4.resize(q, M);

  sum_nx.resize(q, M); sum_nx_sq.resize(q, M);
  sum_nx2.resize(q, M); sum_nx2_sq.resize(q, M); sum_nx2_cub.resize(q, M);
  sum_nx3.resize(q, M); sum_nx4.resize(q, M); sum_nx6.resize(q, M);

  A.resize(q, M);
  EA.resize(q, M); EA2.resize(q, M);
  EA3.resize(q, M); EA4.resize(q, M);

  variance.resize(q, M); skewness.resize(q, M); kurtosis.resize(q, M);

  S1.resize(q, M); S2.resize(q, M); 
  S1g.resize(q, M); S2g.resize(q, M); L.resize(q, M); 
  SkipBeta.resize(q, M); SkipGamma.resize(q, M); Skip.resize(q, M);

  R.resize(q, M); 
  PvalBeta.resize(q, M); PvalGamma.resize(q, M); Pval.resize(q, M);
  
  // DKAT
  D.resize(q, M);
  Dm1.resize(q, M); Dm2.resize(q, M); Dm3.resize(q, M);

  PvalD.resize(q, M);
  ShapeD.resize(q, M); ScaleD.resize(q, M); LocationD.resize(q, M);
}

MCCResults::~MCCResults() { }

//---------------------
// MCC Constructors
//---------------------

MCC::MCC() { }

//---------------------
// Set up
//---------------------

void MCC::setup_y(const MatrixXb &masked_indivs, const MatrixXd &res, unsigned int _n_covariates)
{
  // check dimensions
  if(masked_indivs.rows() == 0) { throw std::runtime_error("setup_y: masked_indivs.rows() == 0"); }
  if(masked_indivs.cols() == 0) { throw std::runtime_error("setup_y: masked_indivs.cols() == 0"); }
  if(res.rows() == 0) { throw std::runtime_error("setup_y: res.rows() == 0"); }
  if(res.cols() == 0) { throw std::runtime_error("setup_y: res.cols() == 0"); }
  if(masked_indivs.rows() != res.rows()) { throw std::runtime_error("setup_y: masked_indivs.rows() != res.rows()"); }
  if(masked_indivs.cols() != res.cols()) { throw std::runtime_error("setup_y: masked_indivs.cols() != res.cols()"); }

  n_samples = res.rows();
  n_traits = res.cols();
  n_covariates = _n_covariates;

  Mask = masked_indivs;
  Neff = Mask.cast<double>().colwise().sum(); 
  
  // center & set Y = 0 for missing values
  Yres = res;
  ArrayXd sums = Mask.select(Yres, 0.0).colwise().sum();
  ArrayXd means = sums / Neff.array();
  Yres.array().rowwise() -= means.transpose();
  Yres.array() *= Mask.array().cast<double>().array(); 

  // normalized Y: centered & sum(y^2) = 1
  // don't confuse with scaled Y, for which sum(y^2) = (N - 1)
  Ynorm = Yres; // Ynorm is centered because of centered Yres
  ArrayXd norms = Ynorm.colwise().norm();
  Ynorm.array().rowwise() /= norms.transpose();
  Ynorm.array() *= Mask.array().cast<double>().array(); 

  precomp_sumy();
}

void MCC::check_setup_data() 
{
  if(n_traits == 0) { throw std::runtime_error("check_setup_data: n_traits == 0"); }

  if(Mask.cols() == 0) { throw std::runtime_error("check_setup_data: Mask.cols() == 0"); }
  if(Mask.cols() != n_traits) { throw std::runtime_error("check_setup_data: Mask.cols() != n_traits"); }
  if(Neff.size() == 0) { throw std::runtime_error("check_setup_data: Neff.size() == 0"); }
  if(Neff.size() != n_traits) { throw std::runtime_error("check_setup_data: Neff.size() != n_traits"); }

  if(Yres.cols() == 0) { throw std::runtime_error("check_setup_data: Yres.cols() == 0"); }
  if(Yres.cols() != n_traits) { throw std::runtime_error("check_setup_data: Yres.cols() != n_traits"); }

  if(sum_y.size() == 0) { throw std::runtime_error("check_setup_data: sum_y.size() == 0"); }
}

//---------------------
// Pre-compute
//---------------------

void MCC::precomp_sumy()
{
  sum_y = Yres.colwise().sum();
  sum_y_sq = sum_y.array().square();

  sum_y2 = Yres.array().square().colwise().sum();
  sum_y2_sq = sum_y2.array().square();

  sum_y3 = Yres.array().pow(3).colwise().sum();
  sum_y4 = Yres.array().pow(4).colwise().sum();

  // sums for Ynorm
  sum_ny = Ynorm.colwise().sum();
  sum_ny_sq = sum_ny.array().square();

  sum_ny2 = Ynorm.array().square().colwise().sum();
  sum_ny2_sq = sum_ny2.array().square();
  sum_ny2_cub = sum_ny2.array().pow(3);

  sum_ny3 = Ynorm.array().pow(3).colwise().sum();
  sum_ny4 = Ynorm.array().pow(4).colwise().sum();
  sum_ny6 = Ynorm.array().pow(6).colwise().sum();
}

//---------------------
// Association (X, Y)`
//---------------------

MCCResults MCC::run(const Eigen::MatrixXd& G)
{
  check_setup_data();

  MCCResults mcc_results(n_samples, n_covariates, n_traits, G.cols(), Neff);
  mcc_results.statistics(Mask, G, Yres);
  /* mcc_results.expectations(Mask, G, */
  /*   sum_y, sum_y_sq, sum_y2, sum_y2_sq, sum_y3, sum_y4); */
  /* mcc_results.moments(); */
  /* mcc_results.distr(); */

  mcc_results.dkat(Mask, G, Ynorm, 
    sum_ny, sum_ny_sq, sum_ny2, sum_ny2_sq, sum_ny2_cub, sum_ny3, sum_ny4, sum_ny6);

  return(mcc_results);
}

void MCCResults::statistics(const MatrixXb &Mask, const Eigen::MatrixXd &G, const Eigen::MatrixXd &Y)
{
  unsigned int i, j;
  for(i = 0; i < q; i++) {
    for(j = 0; j < M; j++) {
      A(i, j) = (Mask.col(i).select(G.col(j), 0.0).array() * Y.col(i).array()).sum();
    }
  }
}

void MCCResults::expectations(const MatrixXb &Mask, const Eigen::MatrixXd & G,
    const VectorXd &sum_y, const VectorXd &sum_y_sq, 
    const VectorXd &sum_y2, const VectorXd &sum_y2_sq,
    const VectorXd &sum_y3, const VectorXd &sum_y4)
{
  unsigned int i, j;
  
  VectorXd size1(q), size2(q), size3(q), size4(q), size5(q);
  MatrixXd term1x(q, M), term2x(q, M), term3x(q, M), term4x(q, M), term5x(q, M);
  VectorXd term1y(q), term2y(q), term3y(q), term4y(q), term5y(q);
  
  //------- Part 0: sums of y, y^2, ...
  //-----------------------------------
  // That is precomputed beforehand 
  // and passed by arguments

  //------- Part 1: sums of x, x^2, ...
  //-----------------------------------
  
  // sum_x[i, j] = (G[, j] * Mask[, i]).colSums()
  for(i = 0; i < q; i++) {
    for(j = 0; j < M; j++) {
      sum_x(i, j) = Mask.col(i).select(G.col(j), 0.0).array().sum();
    }
  }
  // sum_x_sq = (sum_x)^2
  sum_x_sq = sum_x.array().square();

  // sum_x2[i, j] = (G[, j] * Mask[, i]).square().colSums()
  for(i = 0; i < q; i++) {
    for(j = 0; j < M; j++) {
      sum_x2(i, j) = Mask.col(i).select(G.col(j), 0.0).array().square().sum();
    }
  }
  // sum_x2_sq = (sum_x2)^2
  sum_x2_sq = sum_x2.array().square();
 
  // sum_x3[i, j] = (G[, j] * Mask[, i]).pow(3).colSums()
  for(i = 0; i < q; i++) {
    for(j = 0; j < M; j++) {
      sum_x3(i, j) = Mask.col(i).select(G.col(j), 0.0).array().pow(3).sum();
    }
  }

  // sum_x4[i, j] = (G[, j] * Mask[, i]).pow(4).colSums()
  for(i = 0; i < q; i++) {
    for(j = 0; j < M; j++) {
      sum_x4(i, j) = Mask.col(i).select(G.col(j), 0.0).array().pow(4).sum();
    }
  }

  //------- Part 2: expectations of sum(x*y), sum(x*y)^2, ...
  //---------------------------------------------------------
  // Expectation 1: E_perm(A) 
  // q x M matrix EA: EA[, j] = (sum_x[, j] * sum_y) / n
  for(j = 0; j < M; j++) {
    EA.col(j) = (sum_x.col(j).array() * sum_y.array()) / n.array();
  }

  // Expectation 2: E_perm(A2) 
  // q x M matrix EA2: EA[, j] = (sum_x2[, i] * sum_y2) / n +
  //    (sum_x_sq[, i] - sum_x2[, i]) * (sum_y_sq - sum_y2) / (n2 - n)
  for(j = 0; j < M; j++) {
    EA2.col(j) = (sum_x2.col(j).array() * sum_y2.array()) / n.array() +
      (sum_x_sq.col(j).array() - sum_x2.col(j).array())  * (sum_y_sq.array() - sum_y2.array()) / (n2.array() - n.array());
  }

  // Expectation 3: E_perm(A3) 
  size1 = n;
  size2 = 3 * n.array() * (n.array() - 1);
  size3 = n2.array()*n.array() - size2.array() - size1.array();
  // fill qxM matrices for terms for X
  term1x = sum_x3;
  term2x = 3*(sum_x2.array() * sum_x.array() - sum_x3.array());
  term3x = sum_x_sq.array() * sum_x.array() - term2x.array() - term1x.array();
  // fill qx1 vectors for terms for Y
  term1y = sum_y3;
  term2y = 3*(sum_y2.array() * sum_y.array() - sum_y3.array());
  term3y = sum_y_sq.array()*sum_y.array() - term2y.array() - term1y.array();
  // NB: expect (all(n > 2))
  for(j = 0; j < M; j++) {
    EA3.col(j) = (term1x.col(j).array() * term1y.array()) / size1.array() + (term2x.col(j).array() * term2y.array()) / size2.array() +
      // if n > 2
      (term3x.col(j).array() * term3y.array()) / size3.array();
  }

  // Expectation 4: E_perm(A^4)
  size1 = n;
  size2 = 4*n.array()*(n.array() - 1);
  size3 = 3*n.array()*(n.array() - 1);
  size4 = 6*n.array()*(n.array() - 1)*(n.array() - 2);
  size5 = n2.array()*n2.array() - 6*n2.array()*n.array() + 11*n2.array() - 6*n.array();
  // fill qxM matrices for terms for X
  term1x = sum_x4;
  term2x = 4*(sum_x3.array() * sum_x.array() - term1x.array());
  term3x = 3*(sum_x2_sq.array() - term1x.array());
  term4x = 6*sum_x2.array()*sum_x_sq.array() - 12*sum_x3.array()*sum_x.array() - 6*sum_x2_sq.array() + 12*sum_x4.array();
  term5x = sum_x_sq.array()*sum_x_sq.array() - term4x.array() - term3x.array() - term2x.array() - term1x.array();
  // fill qx1 vectors for terms for Y
  term1y = sum_y4.array();
  term2y = 4*(sum_y3.array() * sum_y.array() - term1y.array());
  term3y = 3*(sum_y2_sq.array() - term1y.array());
  term4y = 6*sum_y2.array()*sum_y_sq.array() - 12*sum_y3.array()*sum_y.array() - 6*sum_y2_sq.array() + 12*sum_y4.array();
  term5y = sum_y_sq.array()*sum_y_sq.array() - term4y.array() - term3y.array() - term2y.array() - term1y.array();
  // NB: expect all(n > 3)
  for(j = 0; j < M; j++) {
    EA4.col(j) = 
      // if n > 2
      (term1x.col(j).array() * term1y.array()) / size1.array() +
      (term2x.col(j).array() * term2y.array()) / size2.array() +
      (term3x.col(j).array() * term3y.array()) / size3.array() +
      // if n >= 3
      (term4x.col(j).array() * term4y.array()) / size4.array() +
      // if n > 3
      (term5x.col(j).array() * term5y.array()) / size5.array();
  }
  
  //------- Part 3: observed test statistic
  //---------------------------------------------------------
  for(j = 1; j < M; j++) {
  }
}

void MCCResults::moments()
{
  variance = EA2.array() - EA.array().square(); // variance
  // non-centered A: s = (EA3 - 3*EA*V - 3*EA^3) / V^(3/2) 
  skewness = EA3.array() / variance.array().pow(1.5); // skewness
  kurtosis = EA4.array() / variance.array().square() - 3; // kurtosis
}

void MCCResults::distr()
{
  unsigned int i, j;
  
  for(i = 0; i < q; i++) {
    for(j = 0; j < M; j++) {
      double V = variance(i, j);
      double k = kurtosis(i, j), s = skewness(i, j);
      double k2 = k*k, s2 = s*s;
      double s3 = s2*s, s4 = s2*s2;

      double r = A(i, j) / sqrt(variance(i, j) * (n(i) - 1.0));
      R(i, j) = r;
      
      // --- 1. Fit Beta(alpha, beta) distribution
      double alpha_d1 = -k2*s2 + 32*k2 - 84*k*s2 + 96*k + 36*s4 - 180*s2;
      double alpha_d2 = 2*k - 3*s2;
      bool skip = (alpha_d2 >= 0) | (alpha_d2 == 0);
      
      double alpha_t1 = skip ? 0.0 : sqrt(-1.0/alpha_d1);
      double alpha_t2 = (36*s - 18*s3 + 3*k2*s - 3*k*s3 + 24*k*s) * alpha_t1;
      double alpha_t3 = 3*k - 3*s2 + 6;
      double alpha_t4 = -6*s2 + 6*k + 12;

      double alpha1 = (alpha_t3 + alpha_t2 - alpha_t4)/alpha_d2;
      double alpha2 = (alpha_t3 - alpha_t2 - alpha_t4)/alpha_d2;
      double beta1 = -(alpha_t3 + alpha_t2)/alpha_d2;
      double beta2 = -(alpha_t3 - alpha_t2)/alpha_d2;

      bool switch_alpha = (alpha1 <= 0) | (beta1 <= 0);
      double alpha = switch_alpha ? alpha2 : alpha1;
      double beta = switch_alpha ? beta2 : beta1;

      skip = skip | ((alpha < 0) & (beta < 0));
      alpha = skip ? 0.0 : alpha;
      beta = skip ? 0.0 : beta;

      S1(i, j) = alpha;
      S2(i, j) = beta;

      // --- 2. Fit Beta Gamma(alpha, beta, location) distribution
      bool skip_gamma = (abs(s) < nl_dbl_dmin);

      double m2 = V;
      double m3 = abs(s);
      double alpha_gamma = 4.0 / (m3*m3);
      double beta_gamma = sqrt(alpha_gamma / m2);
      double location_gamma = - alpha_gamma * beta_gamma;

      /* double alpha_gamma = skip_gamma ? 0.0 : (4 / s2); */
      /* double beta_gamma = skip_gamma ? 0.0 : sqrt(alpha_gamma / V); */
      /* double location_gamma = skip_gamma ? 0.0 : (-alpha_gamma * beta_gamma); */
  
      S1g(i, j) = alpha_gamma;
      S2g(i, j) = beta_gamma;
      L(i, j) = location_gamma;

      // ---- 3. P-value calculation for Beta
      double rprime;
      double pval_right, pval_left, pval_double;

      bool skip_test = false;
      if(!skip) {
        double alpha_beta = alpha + beta;
        double mean_beta = alpha / alpha_beta;
        double var_beta = (alpha * beta) / (alpha_beta*alpha_beta * (alpha_beta + 1.0));

        double c0 = mean_beta;
        double c1 = sqrt(var_beta * (n(i) - 1.0));
        rprime = c0 + c1 * r;

        skip_test = (rprime < 0) | (rprime > 1);
        if(!skip_test) {
          boost::math::beta_distribution<> dist(alpha, beta);
          double pval_right = boost::math::cdf(boost::math::complement(dist, rprime));
          double pval_left = boost::math::cdf(dist, rprime);
          double pval_double = 2*min(pval_right, pval_left);
          if(pval_double > 1.0) {
            pval_double = 1.0;
          }
          PvalBeta(i, j) = pval_double;
        }
      }
      SkipBeta(i, j) = skip | skip_test;

      // ---- 4. P-value calculation for Gamma
      bool skip_test_gamma = false;
      if(!skip_gamma) {
        // flip the sign of test statistic r?
        bool flip_sign = (s < 0);
        double mult_flip = flip_sign ? -1.0 : 1.0;
        rprime = alpha_gamma * beta_gamma + beta_gamma * beta_gamma * sqrt(V * (n(i) - 1.0)) * (mult_flip * r);

        skip_test_gamma = (rprime < 0);
        if(!skip_test_gamma) {
          /* cout << " s2 = " << s2 << " 1/s2 = " << 1/s2 << " alpha_gamma = " << alpha_gamma << " " << beta_gamma << " " << 1.0 / beta_gamma << endl; */
          /* cout << " rprime = " << rprime << endl; */
          // CDF from Boost library
          /* boost::math::gamma_distribution<> dist_gama(alpha_gamma, 1.0 / beta_gamma); */
          /* pval_right = boost::math::cdf(boost::math::complement(dist_gama, rprime)); */
          /* pval_left = boost::math::cdf(dist_gama, rprime); */
          // CDF from R 
          /* pval_right = rbase_pgamma(rprime, alpha_gamma, beta_gamma , 0, 0); // 0 = lower tail, 0 = log_p */
          /* pval_left = rbase_pgamma(rprime, alpha_gamma, beta_gamma , 1, 0); // 1 = lower tail, 0 = log_p */
          // CDF from Boost
          pval_right = boost_pgamma(rprime, alpha_gamma, 1.0 / beta_gamma , false); // upper tail
          pval_left = boost_pgamma(rprime, alpha_gamma, 1.0 / beta_gamma , true); // lower tail
          /* cout << "rprime|alpha|beta|pval_right " << rprime << "|" << alpha_gamma << "|" << beta_gamma << "|" << pval_right << endl; */
          pval_double = 2*min(pval_right, pval_left);
          if(pval_double > 1.0) {
            pval_double = 1.0;
          }
          PvalGamma(i, j) = pval_double;
        }
      }
      SkipGamma(i, j) = skip_gamma | skip_test_gamma;

      // --- Combine Beta and Gamma results
      // Combine v1: Gamma preferred over Beta
      if(SkipBeta(i, j) & SkipGamma(i, j)) {
        // case 1: both Gamma and Beta fail
        Skip(i, j) = true;
      } else if(!SkipGamma(i, j)) {
        // case 2: Gamma ok
        Skip(i, j) = false;
        Pval(i, j) = PvalGamma(i, j);
      } else if(!SkipBeta(i, j)) {
        // case 3: Gamma fails, Beta ok
        Skip(i, j) = false;
        Pval(i, j) = PvalBeta(i, j);
      }
      
      // Combine v2: Gamma
      /* Skip(i, j) = SkipGamma(i, j); */
      /* Pval(i, j) = PvalGamma(i, j); */
      
      // Combine v3: Beta
      /* Skip(i, j) = SkipBeta(i, j); */
      /* Pval(i, j) = PvalBeta(i, j); */
    } // end of loop for (i, j) = (trait, variant) pair
  }
}

//-----------------
// DKAT
//-----------------

void MCCResults::dkat(const MatrixXb &Mask, const Eigen::MatrixXd & G,
    const Eigen::MatrixXd & Ynorm, 
    const VectorXd &sum_ny, const VectorXd &sum_ny_sq, 
    const VectorXd &sum_ny2, const VectorXd &sum_ny2_sq, const VectorXd &sum_ny2_cub,
    const VectorXd &sum_ny3, const VectorXd &sum_ny4, const VectorXd &sum_ny6)
{
  unsigned int i, j;
  double ni;
  MatrixXd X(N, M);
  ArrayXd means_X(M), norms_X(M);
  double pval_right;
  
  for(i = 0; i < q; i++) {
    ni = n(i); // Neff for trait i
    // create a copy of genotype matrix & use mask for trait i
    X = G;

    // normalize genotype j: center + normalize
    for(j = 0; j < M; j++) {
      means_X(j) = Mask.col(i).select(X.col(j), 0.0).sum() / ni;
    }
    X.array().rowwise() -= means_X.transpose();
    for(j = 0; j < M; j++) {
      norms_X(j) = Mask.col(i).select(X.col(j), 0.0).norm();
    }
    X.array().rowwise() /= norms_X.transpose();
    for(j = 0; j < M; j++) {
      X.col(j) = Mask.col(i).select(X.col(j), 0.0);
    }

    // compute sums
    sum_nx.row(i) = X.colwise().sum();
    sum_nx_sq.row(i) = sum_nx.row(i).array().square();
    sum_nx2.row(i) = X.array().square().colwise().sum();
    sum_nx2_sq.row(i) = sum_nx2.row(i).array().square();
    sum_nx2_cub.row(i) = sum_nx2.row(i).array().pow(3);

    sum_nx3.row(i) = X.array().pow(3).colwise().sum();
    sum_nx4.row(i) = X.array().pow(4).colwise().sum();
    sum_nx6.row(i) = X.array().pow(6).colwise().sum();

    // test statistic D = R^2 (squared Pearson corr.)
    D.row(i) = (X.transpose() * Ynorm.col(i)).array().square();

    // Moment 1 of D: mean
    Dm1.row(i) = sum_nx2.row(i) * sum_ny2(i) / ni;

    // Momemnt 2 of D: variance
    double T = sum_ny2(i), T2 = sum_ny2_sq(i), S2 = sum_ny4(i);
    ArrayXd Ts = sum_nx2.row(i), T2s = sum_nx2_sq.row(i), S2s = sum_nx4.row(i);

    double T_sq = T*T;
    double T_cub = T_sq*T;
    ArrayXd Ts_sq = Ts.square();
    ArrayXd Ts_cub = Ts_sq * Ts;

    double ni_1 = ni - 1.0, ni_2 = ni - 2.0, ni_3 = ni - 3.0; 
    /* double ni1 = ni + 1.0, ni2 = ni + 2.0, ni3 = ni + 3.0, ni4 = ni + 4.0; */
    double ni1 = ni + 1.0, ni4 = ni + 4.0;
    double ni_sq = ni*ni;
    double ni_cub = ni_sq * ni, ni_quad = ni_sq * ni_sq;
    
    ArrayXd temp1 = 2.0 * (ni_1*T2 - T_sq)*(ni_1*T2s - Ts_sq) / (ni_1*ni_1*ni1*ni_2);
    double temp21 = ni*ni1*S2 - ni_1*(T_sq + 2*T2);
    ArrayXd temp22 = ni*ni1*S2s - ni_1*(Ts_sq + 2*T2s);
    double temp23 = ni1*ni*ni_1*ni_2*ni_3;
    ArrayXd temp2 = temp21 * temp22 / temp23;
    Dm2.row(i) = temp1 + temp2;
    
    // Momemnt 3 of D: skewness
    double T3 = sum_ny2_cub(i), S3 = sum_ny6(i);
    ArrayXd T3s = sum_nx2_cub.row(i), S3s = sum_nx6.row(i);
    double U = sum_ny3(i) * sum_ny3(i);
    ArrayXd Us = sum_nx3.row(i).array().square();
    double R = sum_ny2(i) * sum_ny4(i);
    ArrayXd Rs = sum_nx2.row(i).array() * sum_nx4.row(i).array();
    double B = U;
    ArrayXd Bs = Us;

    ArrayXd t1 = ni_sq * ni1 * (ni_sq + 15*ni - 4) * S3 * S3s;
    ArrayXd t2 = 4 * (ni_quad - 8*ni_cub + 19*ni_sq - 4*ni - 16) * U * Us;
    ArrayXd t3 = 24 * (ni_sq - ni - 4) * (U * Bs + B * Us);
    ArrayXd t4 = 6 * (ni_quad - 8*ni_cub + 21*ni_sq - 6*ni - 24) * B * Bs;

    ArrayXd t5 = 12 * (ni_quad - ni_cub - 8*ni_sq + 36*ni - 48) * R * Rs;
    ArrayXd t6 = 12 * (ni_cub - 2*ni_sq + 9*ni - 12) * (T*S2*Rs + R*Ts*S2s); 
    ArrayXd t7 = 3 * (ni_quad - 4*ni_cub - 2*ni_sq + 9*ni - 12) * T*Ts*S2*S2s;

    ArrayXd t81 = (ni_cub - 3*ni_sq - 2*ni + 8) * (R*Us + U*Rs);
    ArrayXd t82 = (ni_cub - 2*ni_sq - 3*ni + 12) * (R*Bs + B*Rs);
    ArrayXd t8 = 24 * (t81 + t82);
    ArrayXd t9 = 12 * (ni_sq - ni + 4) * (T*S2*Us + U*Ts*S2s);
    ArrayXd t10 = 6 * (2*ni_cub - 7*ni_sq - 3*ni + 12) * (T*S2*Bs + B*Ts*S2s);

    ArrayXd t11 = -2*ni*ni_1*(ni_sq - ni + 4) * ((2*U + 3*B)*S3s + (2*Us + 3*Bs)*S3);
    ArrayXd t12 = -3*ni*ni_1*ni_1*ni4 * ((T*S2 + 4*R)*S3s + (Ts*S2s + 4*Rs)*S3);
    ArrayXd t13 = 2*ni*ni_1*ni_2 *((T_cub + 6*T*T2 + 8*T3)*S3s + (Ts_cub + 6*Ts*T2s + 8*T3s)*S3);
    ArrayXd t14 = T_cub * ((ni_cub - 9*ni_sq + 23*ni - 14)*Ts_cub + 6*(ni - 4)*Ts*T2s + 8*T3s);
    ArrayXd t15 = 6*T*T2*((ni - 4)*Ts_cub + (ni_cub - 9*ni_sq + 24*ni - 14)*Ts*T2s + 4*ni_3*T3s);

    ArrayXd t16 = 8*T3*(Ts_cub + 3*ni_3*Ts*T2s + (ni_cub - 9*ni_sq + 26*ni - 22)*T3s);
    ArrayXd t17 = -16*(T_cub*Us + U*Ts_cub) - 6*(T*T2*Us + U*Ts*T2s) * (2*ni_sq - 10*ni + 16);
    ArrayXd t18 = -8*(T3*Us + U*T3s) * (3*ni_sq - 15*ni + 16) - (T_cub*Bs + B*Ts_cub) * (6*ni_sq - 30*ni + 24);
    ArrayXd t19 = -6*(T*T2*Bs + B*Ts*T2s) * (4*ni_sq - 20*ni + 24) - 8*(T3*Bs + B*T3s) * (3*ni_sq - 15*ni + 24);

    ArrayXd t201 = 24*(T_cub*Rs + R*Ts_cub) + 6*(T*T2*Rs + R*Ts*T2s) * (2*ni_sq- 10*ni + 24);
    ArrayXd t202 = 8*(T3*Rs + R*T3s) * (3*ni_sq - 15*ni + 24) + (3*ni_sq - 15*ni + 6) * (T_cub*Ts*S2s + T*S2*Ts_cub);
    ArrayXd t203 = 6*(T*T2*Ts*S2s + Ts*T2s*T*S2) * (ni_sq - 5*ni + 6) + 48*(T3*Ts*S2s + T3s*T*S2);
    ArrayXd t20 = -ni_2 * (t201 + t202 + t203);

    ArrayXd temp31 = t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10 + t11 + t12 + t13 + t14 + t15 + t16 + t17 + t18 + t19 + t20;
    double temp32 = ni * ni_1 * ni_2 * ni_3 * (ni_3 - 1) * (ni_3 - 2);
    ArrayXd mom3 = temp31 / temp32;

    Dm3.row(i) = (mom3.transpose() - 3*Dm1.row(i).array() * Dm2.row(i).array() - Dm1.row(i).array().pow(3)) / Dm2.row(i).array().pow(1.5);
    
    // Parameters of Gamma distribution of D: shape, scale, location
    /* shape=4/m3^2 */
    /* scale=sqrt(m2)*m3/2 */
    /* location=m1-2*sqrt(m2)/m3 */
    ShapeD.row(i).array() = 4.0 / Dm3.row(i).array().square();
    ScaleD.row(i).array() = Dm2.row(i).array().sqrt() * Dm3.row(i).array() / 2.0;
    LocationD.row(i).array() = Dm1.row(i).array() - 2 * Dm2.row(i).array().sqrt() / Dm3.row(i).array();

    for(j = 0; j < M; j++) {
      /* gscale = abs(scale) */ 
      /* ssgn = sign(scale) */
      /* lower.tail = FALSE */
      /* pval = pgamma(ssgn*(Fstar - location), shape = shape, scale = gscale, */
      /*   lower.tail = xor(ssgn < 0, lower.tail), log.p = FALSE) */
      /* cout << "skewness = " << Dm3(i, j) << endl; */
      /* if(abs(Dm3(i, j)) < 1e-6) { */
      /*   /1* cout << "approx" << endl; *1/ */
      /*   boost::math::normal norm(Dm1(i, j), sqrt(Dm2(i, j))); */
      /*   pval_right = boost::math::cdf(boost::math::complement(norm, D(i, j))); */
      /* } else { */
        /* pval_right = rbase_pgamma(D(i, j) - LocationD(i, j), ShapeD(i, j), ScaleD(i, j) , 0, 0); // 0 = lower tail, 0 = log_p */
        pval_right = boost_pgamma(D(i, j) - LocationD(i, j), ShapeD(i, j), ScaleD(i, j) , false);
      /* } */
      PvalD(i, j) = pval_right;

      Skip(i, j) = false;
      Pval(i, j) = PvalD(i, j);
    }
  }
}
