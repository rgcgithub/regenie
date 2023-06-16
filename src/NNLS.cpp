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

#include "NNLS.hpp"
#include "mvtnorm/mvtnorm.h"

using namespace Eigen;
using namespace std;

/****************************************************
 * Compute weights for the NNLS test
 * ic.weights function in the R package ic.infer
****************************************************/

Eigen::MatrixXd _inverse(const Eigen::MatrixXd& V)
{
  int n = V.rows();
  MatrixXd Vinv = V.llt().solve(MatrixXd::Identity(n, n));
  /* MatrixXd Vinv = V.inverse(); */
  return(Vinv);
}

inline void _assign_wts(vector<double> &wts, int index, double value, int verbose) {
  wts[index] = value;
  if(verbose > 1) cout << " w[" << index << "] = " << value << 
    " (equal to zero = " << (value == 0.0) << ")" << endl;
}

void _complement(int n, const vector<int> &s1, vector<int> &s2)
{
  int n_s1 = s1.size();
  s2.reserve(n - n_s1);
  
  std::vector<char> s1_true(n, false);
  for(int i : s1) {
    s1_true[i] = true;
  }

  for(int i = 0; i < n; ++i) {
    if(!s1_true[i]) {
      s2.push_back(i);
    }
  }
}


Eigen::MatrixXd _subset_matrix(const Eigen::MatrixXd& V, 
    vector<int> &rows, vector<int> &cols)
{
  int ns = rows.size(), ps = cols.size(); // nxp dimensions of subset matrix
  Eigen::MatrixXd S(ns, ps);

  for(int i = 0; i < ns; ++i) {
    for(int j = 0; j < ps; ++j) {
      S(i, j) = V(rows[i], cols[j]);
    }
  }

  return(S);
}

inline Eigen::MatrixXd _subset_matrix1(const Eigen::MatrixXd& V, 
    vector<int> &rows, vector<int> &cols)
{
  return(V(rows, cols));
}

Eigen::MatrixXd jburden_subset_matrix(const Eigen::MatrixXd& V, 
    vector<int> &rows, vector<int> &cols, int method)
{
  if(method == 0) {
    return(_subset_matrix(V, rows, cols));
  } else if(method == 1) {
    return(_subset_matrix1(V, rows, cols));
  } else {
    throw std::runtime_error("jburden_subset_matrix: unknown method"); 
  }
}

void _wts_subset(const Eigen::MatrixXd &V, vector<int> &s1, double &w1, double &w2, int index)
{
  int n = V.rows(); // n rows x n colums

  vector<int> s2;
  _complement(n, s1, s2);

  // subset V: V11, V22, V12 & V21
  Eigen::MatrixXd V11 = _subset_matrix(V, s1, s1);
  Eigen::MatrixXd V22 = _subset_matrix(V, s2, s2);
  Eigen::MatrixXd V12 = _subset_matrix(V, s1, s2);
  Eigen::MatrixXd V21 = _subset_matrix(V, s2, s1);

  // matrix V220 = V22 - V21 solve(V11, V12)
  MatrixXd V11inv = _inverse(V11);
  MatrixXd V220 = V22 - V21 * V11inv * V12;

  w1 = jburden_pnorm(V11inv) * jburden_pnorm(V220);

  // needed only if n is odd
  // example 1: skip if n = 4 & index = 1 is the last
  // example 2: pass if n = 5
  if(2*(index + 1) != (n - 2)) {
    // matrix V110 = V11 - V12 solve(V22, V21)
    MatrixXd V22inv = _inverse(V22);
    MatrixXd V110 = V11 - V12 * V22inv * V21;

    w2 = jburden_pnorm(V22inv) * jburden_pnorm(V110);
  } else {
    w2 = 0.0;
  }
}

// main function: compute weights for the NNLS test
Eigen::VectorXd jburden_wts(const Eigen::MatrixXd& V, int verbose)
{
  int n = V.rows(); // n rows x n colums

  // output vertor of weihts: w(n), w(n-1), ..., w(1), w(0)
  int nw = n + 1;  
  vector<double> wts(nw, 0.0);

  // step 1: wts[0] = w(n)
  /* wts[0] = jburden_pnorm(V); */
  // cerr << jburden_pnorm(V) << endl;
  _assign_wts(wts, 0, jburden_pnorm(V), verbose);
  
  // step 2: wts[nw-1] = w(0)
  MatrixXd Vinv = _inverse(V);
  /* wts[nw - 1] = jburden_pnorm(Vinv); */
  _assign_wts(wts, nw - 1, jburden_pnorm(Vinv), verbose);

  // step 3: wts[1] = w(3)
  int nk = floor((n - 2) / 2);
  for(int k = 0; k < nk; ++k) {
    // step 3.1: get all subsets of size (k+1) 
    if(verbose) cout << " k " << k << " / " << nk << endl;

    list<vector<int>> subsets;
    jburden_nchoosek(n, k + 1, subsets);
    int n_subsets = subsets.size();
    if(verbose) cout << "   #subsets = " << n_subsets << endl;

    vector<double> w1(n_subsets);
    vector<double> w2(n_subsets);
    int i = 0;
    list<vector<int>>::iterator it = subsets.begin(); 
    for(; it != subsets.end(); ++it, ++i) {
      double w1_i, w2_i;
      _wts_subset(V, *it, w1_i, w2_i, k);

      w1[i] = w1_i;
      w2[i] = w2_i;
    }

    // assign weights
    /* wts[k + 1] = accumulate(w1.begin(), w1.end(), 0.0); */
    _assign_wts(wts, k + 1, accumulate(w1.begin(), w1.end(), 0.0), verbose);

    /* wts[n - k - 1] = accumulate(w2.begin(), w2.end(), 0.0); */
    _assign_wts(wts, nw - k - 2, accumulate(w2.begin(), w2.end(), 0.0), verbose);
  }
  
  // two sums of odd and even weights
  double sum_odd = 0.0, sum_even = 0.0;
  if(n % 2 == 0) {
    for(int i = 0; i < nw; ++i) {
      if(i % 2 == 0) sum_even += wts[i];
      else sum_odd += wts[i];
    }
  } else {
    for(int i = 0; i < nw; ++i) {
      if(i % 2 == 0) sum_odd += wts[i];
      else sum_even += wts[i];
    }
  }
  double sub_odd = 0.5 - sum_odd, sub_even = 0.5 - sum_even;

  if( n % 2 == 0) {
    // n is even
    if(n % 4 == 0) {
      // example: n = 4
      _assign_wts(wts, n/2, sub_even, verbose);
      _assign_wts(wts, n/2 + 1, sub_odd, verbose);
    } else {
      // example: n = 6
      _assign_wts(wts, n/2, sub_odd, verbose);
      _assign_wts(wts, n/2 + 1, sub_even, verbose);
    }
  } else {
    // n is odd
    if((n + 1) % 4 == 0) {
      // example: n = {3, 7}
      _assign_wts(wts, (n + 1)/2 - 1, sub_even, verbose);
      _assign_wts(wts, (n + 3)/2 - 1, sub_odd, verbose);
    } else {
      // example: n = {5, 9}
      _assign_wts(wts, (n + 1)/2 - 1, sub_odd, verbose);
      _assign_wts(wts, (n + 3)/2 - 1, sub_even, verbose);
    }
  }

  VectorXd ret(nw);
  for(int i = 0; i < nw; ++i) {
    ret[i] = wts[i];
  }

  return(ret);
}

double jburden_pnorm(const Eigen::MatrixXd& A, 
  int maxpts, double abseps, int verbose)
{
  int i, j, k; // iterators
  int n = A.rows(); // n rows x n colums
  int nc = (n*n - n)/2; // number of elements in low-tri part of A
  double* bound = new double[n];
  double* cmat = new double[nc];
  double error, ret; 

  // fill n bounds with 0
  for(i = 0; i < n; ++i) {
    bound[i] = 0.0;
  }
  
  // convert A to correlation matrix 
  // (requirement of pmvnorm)
  /* Eigen::MatrixXd C = wts_cov2cor(A); */
  vector<double> sd(n);
  for(int i = 0; i < n; ++i) {
    sd[i] = sqrt(A(i, i));
  }
  

  // fill low-tri version cmat of correlation matrix A
  /* string name = "C.matirx.tmp"; */
  /* ofstream file(name.c_str()); */
  k = 0; // counter for filled entries in cmat
  for(i = 0; i < n; ++i) {
    for(j = 0; j < i; ++j) {
      cmat[k] = A(i, j) / (sd[i]*sd[j]);
      k += 1;
      /* file << cmat[k] << "\t"; */
    }
    /* file << "\n"; */
  }
  /* file.close(); */
 
  // call C++ function that calls the Fortran code
  ret = pmvnorm_complement(n, maxpts, abseps, bound, cmat, &error);
  //ret = 1.0;

  return(ret);
}

/****************************************************
 * Enumerate all sets of k out of n numbers
 * nchoosek function in the R package ic.infer
****************************************************/

// the number of all set of k out of n
// Implementation in Boost Library: https://www.boost.org/doc/libs/1_56_0/libs/math/doc/html/math_toolkit/factorials/sf_binomial.html
// - Binomial coefficients are calculated using table lookup of factorials where possible.
// - Otherwise, it is implemented in terms of the beta function using the relations.
double jburden_choose_boost(int n, int k)
{
  double ret = boost::math::binomial_coefficient<double>(n, k);
  return ret;
}

// the number of all set of k out of n
int jburden_choose(int n, int k)
{
  if(k == 0) return 1;

  if(k > n /2) return(jburden_choose(n, n - k));
  
  long res = 1;

  for(int i = 1; i <= k; ++i) {
    res *= n - i + 1;
    res /= i;
  }

  return res;
}

// recursive algorithm for nchoosek
void _nchoosek_rec(int size, int left, int index, 
  vector<int> &l, list<vector<int>> &ll)
{
  if(left == 0) { 
    ll.push_back(l);
    return;
  }
  
  for(int i = index; i < size; i++) {
    l.push_back(i);
    _nchoosek_rec(size, left - 1, i + 1, l, ll);
    l.pop_back();
  }
}

void jburden_nchoosek(int n, int k, list<vector<int>> &ll)
{
  vector<int> l;   
  l.reserve(k);
  _nchoosek_rec(n, k, 0, l, ll);
}

/***************************************
 * NLLS model fitting: jburden_fit_nnls
***************************************/

// R: x <- c(TRUE, FALSE); which(x)
vector<int> bool_to_int(vector<bool> &x) 
{ 
  vector<int> y;

  int i = 0;
  for(vector<bool>::iterator it = x.begin(); it != x.end(); ++it, ++i) {
    if(*it) {
      y.push_back(i);
    }
  }
  return(y);
}


// R: x <- c(TRUE, FALSE); max(x)
int sum_bool(vector<bool> &x) 
{ 
  int sum = accumulate(begin(x), end(x), 0);
  return(sum);
}

// R: x <- c(1, 3); subset <- c(TRUE, FALSE); min(x[subset])
double min_vec(VectorXd &x, vector<bool> &subset)
{
  int i = 0;
  vector<bool>::iterator it = subset.begin(); 
  double min;

  // search for first true in subset
  for(; it != subset.end(); ) {
    if(*it) break; 
    ++it;
    ++i;
  }
  // initialize the maximum 
  min = x[i];
  ++it;
  ++i;

  // go over other elements of x
  for(; it != subset.end(); ++it, ++i) {
    if(*it) {
      double val = x[i];
      if(val < min) {
        min = val;
      }
    }
  }

  return(min);
}

// R: x <- c(1, 3); subset <- c(TRUE, FALSE); max(x[subset])
double max_vec(VectorXd &x, vector<bool> &subset)
{
  int i = 0;
  vector<bool>::iterator it = subset.begin(); 
  double max;

  // early stop: all false in subset
  int m = sum_bool(subset);
  if(m == 0) {
    return(0.0);
  }

  // search for first true in subset
  for(; it != subset.end(); ) {
    if(*it) break; 
    ++it;
    ++i;
  }
  // initialize the maximum 
  max = x[i];
  ++it;
  ++i;

  // go over other elements of x
  for(; it != subset.end(); ++it, ++i) {
    if(*it) {
      double val = x[i];
      if(val > max) {
        max = val;
      }
    }
  }

  return(max);
}

// R: x <- c(1, 3); subset <- c(TRUE, FALSE); which.max(x[subset])
int which_max(VectorXd &x, vector<bool> &subset)
{
  int i = 0;
  vector<bool>::iterator it = subset.begin(); 
  double max; // max value
  int m; // index with max value (return)

  // search for first true in subset
  for(; it != subset.end(); ) {
    if(*it) break; 
    ++it;
    ++i;
  }
  // initialize the maximum 
  max = x[i];
  m = i;
  ++it;
  ++i;

  // go over other elements of x
  for(; it != subset.end(); ++it, ++i) {
    if(*it) {
      double val = x[i];
      if(val > max) {
        max = val;
        m = i;
      }
    }
  }

  return(m);
}

// R: x <- c(1, 3); subset <- c(TRUE, FALSE); x[subset]
VectorXd subset_vec(VectorXd &x, vector<bool> subset) 
{
  int m = sum_bool(subset);
  VectorXd y(m);

  int i = 0, j = 0;
  vector<bool>::iterator it = subset.begin(); 
  for(; it != subset.end(); ++it, ++i) {
    if(*it) {
      y[j++] = x[i];
    }
  }

  return(y);
}

// R: X <- matrix(1:9, 3, 3); subset <- 1:2; X[subset, subset]
MatrixXd subset_mat(MatrixXd &X, vector<bool> subset) 
{
  int m = sum_bool(subset);
  MatrixXd Y(m, m);

  vector<int> ind = bool_to_int(subset);

  for(int i = 0; i < m; ++i) {
    int ix = ind[i];
    for(int j = 0; j < m; ++j) {
      int jx = ind[j];
      Y(i, j) = X(ix, jx);
    }
  }

  return(Y);
}

// solve OLS on a subset of variables: y ~ X[, subset] b_subset
VectorXd solve_s(VectorXd &Xty, MatrixXd &XtX, vector<bool> subset)
{
  int p = XtX.rows();
  int m = sum_bool(subset);
  VectorXd b(p); // output vector of betas

  if(m > 0) {
    MatrixXd XtX_subset = subset_mat(XtX, subset);
    VectorXd Xty_subset = subset_vec(Xty, subset);

    LLT<MatrixXd> llt(XtX_subset);
    VectorXd b_subset = llt.solve(Xty_subset);
 
    int i = 0, i_subset = 0;
    vector<bool>::iterator it = subset.begin(); 
    for(; it != subset.end(); ++it, ++i) {
      if(*it) {
        b[i] = b_subset[i_subset];
        ++i_subset;
      } else {
        b[i] = 0.0;
      }
    }
  } else {
    for(int j = 0; j < p; ++j) {
      b[j] = 0.0;
    }
  }

  return(b);
}


/* S <- (s <= tol) & P */
/* alpha <- min(b[S] / (b[S] - s[S])) */
/* b <- b + alpha * (s - b) */
inline void update_b(VectorXd &b, VectorXd &s, vector<bool> &subset, double tol)
{
  // overlap two conditions: s <= tol & subset
  int n = subset.size();
  vector<bool> negs_subset(n);
  for(int i = 0; i < n; i++) {
    negs_subset[i] = subset[i] & (s[i] <= tol);
  }
  
  // convert boolean indices to integer indices
  vector<int> ind = bool_to_int(negs_subset);
  int m = ind.size();

  // initial value of min
  int k = ind[0];
  double min_val = b[k] / (b[k] - s[k]);
  // go over other elements starting from the 2nd index (1)
  for(int i = 1; i < m; ++i) {
    k = ind[i];
    double val = b[k] / (b[k] - s[k]);
    if(val < min_val) {
      min_val = val;
    }
  }

  double alpha = min_val;
  /* VectorXd delta = s - b; */
  b = b + alpha * (s - b);
}

/* l <- which(b <= 0) */
/* P[l] <- FALSE */
/* R[l] <- TRUE */
inline void update_subsets(VectorXd &b, vector<bool> &P, vector<bool> &R, double tol)
{
  int p = b.size();

  for(int i = 0; i < p; ++i) {
    if(b[i] <= tol) {
      P[i] = false;
      R[i] = true;
    }
  }
}

// NNLS model fitting: active set algorithm 
// - summary stat. level data b and V
int jburden_fit_nnls_cprod(const Eigen::VectorXd &Xty_, const Eigen::MatrixXd& XtX_,
  Eigen::VectorXd& bhat_out, vector<bool>& selected_out,
  double tol, bool neg, int maxit, int maxit_inner, int verbose)
{
  VectorXd Xty = Xty_;
  MatrixXd XtX = XtX_;
  const int p = Xty.size();

  if(neg) {
    Xty = -1 * Xty;
  }

  // A. Initialization
  if(verbose > 2) cout << "A. Initialization" << endl;

  vector<bool> P(p, false); // passive set P (empty)
  vector<bool> R(p, true); // active set R (full)

  VectorXd b(p); // vector b of effect sizes
  b.setZero(); // initialize b with zeros

  VectorXd w = Xty; // Lagrange multiplier 
     
  // B. Main loop
  if(verbose > 2) cout << "B. Main loop" << endl;
  
  int cnt_main = 0, cnt_inner = 0;
  int m; // index in R with the maximum w (R subset)
  while((sum_bool(R) > 0) & (max_vec(w, R) > tol) & (cnt_main <= maxit)) {
    if(verbose > 2) cout << " - it (main) " << cnt_main << endl;

    m = which_max(w, R);
    if(verbose > 2) cout << " - max(w) = " << max_vec(w, R) << endl;
    P[m] = true;
    R[m] = false;
    if(verbose > 2) cout << " - sum(R) = " << sum_bool(R) << "; sum(P) = " << sum_bool(P) << endl;

    VectorXd s = solve_s(Xty, XtX, P);
    
    // C. Inner loop
    if(verbose > 2) cout << "C. Inner loop" << endl;

    cnt_inner = 0;
    while(min_vec(s, P) <= tol) {
      if(verbose > 2) cout << " - it (inner) " << cnt_inner << endl;
      
      update_b(b, s, P, tol);
      update_subsets(b, P, R, tol);
      if(verbose > 2) cout << " - sum(R) = " << sum_bool(R) << "; sum(P) = " << sum_bool(P) << endl;

      s = solve_s(Xty, XtX, P);

      cnt_inner++;

      // early break from the inner loop:
      // the set of variables included in the model is empty
      // example: 
      //   - 1st selected variable has b = 1e7, while tol = 1e-6
      //   - so this only variable is removed from set P and P becomes empty.
      if(sum_bool(P) == 0) break;
      // added by JMb (02/16/2022): avoid cases where loops runs indefinitely
      else if(cnt_inner >= maxit_inner) return(-1);
    }
    // early break from the main loop: see above
    if(sum_bool(P) == 0) break;

    if(verbose > 2) cout << "D. Store current results" << endl;
    b = s;
    w = Xty - XtX * b; // X' (y - X b) Lagrange multiplier 

    cnt_main++;
  }
  
  // check convergence
  if(cnt_main >= maxit) {
    return(-1);
  }

  // return
  if(neg) {
    b = -1 * b;
  }
  bhat_out = b;
  selected_out = P; // passive set R (active constraints not applied)

  return(0); // return success
}

// NNLS model fitting: active set algorithm 
// - individual-level data y an X
int jburden_fit_nnls(const Eigen::VectorXd &y, const Eigen::MatrixXd& X,
  Eigen::VectorXd& bhat_out, vector<bool>& selected_out,
  double tol, bool neg, int maxit, int maxit_inner, int verbose)
{
  const int p = X.cols();

  MatrixXd XtX(MatrixXd(p, p).setZero().  
    selfadjointView<Lower>().rankUpdate(X.adjoint()));
  VectorXd Xty = X.adjoint() * y;
  if(neg) {
    Xty = -1 * Xty;
  }

  // A. Initialization
  if(verbose > 2) cout << "A. Initialization" << endl;

  vector<bool> P(p, false); // passive set P (empty)
  vector<bool> R(p, true); // active set R (full)

  VectorXd b(p); // vector b of effect sizes
  b.setZero(); // initialize b with zeros

  VectorXd w = X.adjoint() * y; // Lagrange multiplier 
  if(neg) {
    w = -1 * w;
  }
     
  // B. Main loop
  if(verbose > 2) cout << "B. Main loop" << endl;
  
  int cnt_main = 0, cnt_inner = 0;
  int m; // index in R with the maximum w (R subset)
  while((sum_bool(R) > 0) & (max_vec(w, R) > tol) & (cnt_main <= maxit)) {
    if(verbose > 2) cout << " - it (main) " << cnt_main << endl;

    m = which_max(w, R);
    if(verbose > 2) cout << " - max(w) = " << max_vec(w, R) << endl;
    P[m] = true;
    R[m] = false;
    if(verbose > 2) cout << " - sum(R) = " << sum_bool(R) << "; sum(P) = " << sum_bool(P) << endl;

    VectorXd s = solve_s(Xty, XtX, P);
    
    // C. Inner loop
    if(verbose > 2) cout << "C. Inner loop" << endl;

    cnt_inner = 0;
    while(min_vec(s, P) <= tol) {
      if(verbose > 2) cout << " - it (inner) " << cnt_inner << endl;
      
      update_b(b, s, P, tol);
      update_subsets(b, P, R, tol);
      if(verbose > 2) cout << " - sum(R) = " << sum_bool(R) << "; sum(P) = " << sum_bool(P) << endl;

      s = solve_s(Xty, XtX, P);

      cnt_inner++;

      // early break from the inner loop:
      // the set of variables included in the model is empty
      // example: 
      //   - 1st selected variable has b = 1e7, while tol = 1e-6
      //   - so this only variable is removed from set P and P becomes empty.
      if(sum_bool(P) == 0) break;
      // added by JMb (02/16/2022): avoid cases where loops runs indefinitely
      else if(cnt_inner >= maxit_inner) return(-1);
    }
    // early break from the main loop: see above
    if(sum_bool(P) == 0) break;

    if(verbose > 2) cout << "D. Store current results" << endl;
    b = s;
    w = Xty - XtX * b; // X' (y - X b) Lagrange multiplier 

    cnt_main++;
  }
  
  // check convergence
  if(cnt_main >= maxit) {
    return(-1);
  }

  // return
  if(neg) {
    b = -1 * b;
  }
  bhat_out = b;
  selected_out = P; // passive set R (active constraints not applied)

  return(0); // return success
}

/***************************************
 * NLLS test: jburden_test
***************************************/

// R: pchisq(10, 1, lower = FALSE)
double jburden_pchisq(double q, int df, bool lower)
{
  double p;
  boost::math::chi_squared dist(df);

  if(lower) {
    p = boost::math::cdf(dist, q);
  } else {
    p = boost::math::cdf(boost::math::complement(dist, q));
  }

  return(p);
}

double jburden_pchisq_bar(double x, Eigen::VectorXd& wt)
{
  // order of weigts: w(n), w(n-1), ..., w(1), w(0)
  // - the first weight correspond to df = n
  // - the last weight correspond to df = 0

  int n = wt.size(); // number of weights 
  int df = n - 1; // number of variables tested jointly

  // sum of weights
  double sum_wt = 0;
  for(int i = 0; i < n; ++i) {
    sum_wt += wt[i];
  }

  // sum of (1 - p-value) for non-zero weights
  // pp denotes (1 - p-value)
  double sum_pp = 0.0;
  int i = 0, dfi = df;
  for(; i < (n - 1); ++i, --dfi) {
    double pp_i = jburden_pchisq(x, dfi, false);
    sum_pp += wt[i] * pp_i;

    //debug
    /* cout << " -- jburden_pchisq_bar: i = " << i */ 
    /*   << " pp_i = " << pp_i << " wt[i] = " << wt[i] */ 
    /*   << " sum_pp = " << sum_pp << endl; */
  }

  // output p-value
  // (depreciated) option 1: 
  //  - can produce negative p-values because sum_wt is not exactly 1
  //  - that might happen only when p-values are very small
  /* double p = 1.0 - sum_wt + sum_pp; */ 
  // option 2: 
  // - safe to produce non-negative p-values, as wt >= 0
  double p = sum_pp;
  
  // debug
  /* cout << "jburden_pchisq_bar: sum_wt = " << sum_wt << endl; */
  /* cout << "jburden_pchisq_bar: p = " << p << endl; */

  return(p);
}

/***************************************
 * NNLS adaptive weights
***************************************/

void _sample(int n, int k, vector<int> &l, std::mt19937_64& gen)
{
  unordered_set<int> sample;
  
  for(int i = n - k; i < n; ++i) {
    std::uniform_int_distribution<int> distr(0, i);
    int s = distr(gen);

    if(sample.find(s) == sample.end()) {
      sample.insert(s);
    } else {
      sample.insert(i);
    }
  }

  // copy unordered set (sample) to output vector (l)
  for(unordered_set<int>::iterator it = sample.begin(); it != sample.end(); ++it) {
    l.push_back(*it);
  }
  // sort output vector (l)
  sort(l.begin(), l.end());
}

void jburden_nchoosek_sample(int n, int k, int s, list<vector<int>> &ll, std::mt19937_64& gen)
{
  for(int i = 0; i < s; ++i) {
    vector<int> sample;
    sample.reserve(k);
    _sample(n, k, sample, gen);
    ll.push_back(sample);
  }
}

int jburden_wts_adapt(const Eigen::MatrixXd& V, Eigen::VectorXd& wts_out, std::mt19937_64& gen,
  int n_approx, bool normalize, int verbose)
{
  int n = V.cols(); 
  int nw = n + 1;  

  // output vertor of weights: w(0), w(1), ..., w(n-1), w(n)
  vector<double> wts(nw, 0.0);
  vector<int> ind_exact, ind_approx;
  
  // A. w(n) = wts[nw - 1]
  double wts_n = jburden_pnorm(V);
  if(wts_n < 0) {
    std::cout << "ERROR: computing NNLS weight(n = " << n << ") failed; " 
      << "pnorm(V) returned negative value (" << wts_n << ")\n";
    return(-1);
  }
  _assign_wts(wts, nw - 1, jburden_pnorm(V), verbose);
  ind_exact.push_back(nw - 1);
      
  // B. w(0)
  MatrixXd Vinv = _inverse(V);
  double wts_0 = jburden_pnorm(Vinv);
  if(wts_0 < 0) {
    std::cout << "ERROR: computing NNLS weight(0) failed; " 
      << "pnorm(Vinv) returned negative value (" << wts_0 << ")\n";
    return(-1);
  }
  _assign_wts(wts, 0, wts_0, verbose);
  ind_exact.push_back(0);

  // C. Outer loop over weights: w(1), w(2), ..., w(n-1)
  for(int i = 1; i < (nw - 1); ++i) {
    // NB: choose(40, 20) = 137,846,528,820
    // (depreciated due to overflow) 
    // int n_sets = jburden_choose(n, i);
    double n_sets_numeric = jburden_choose_boost(n, i);
    /* // skip this
    // check if n_set is not overflowed
    int max_int = std::numeric_limits<int>::max();
    bool overflow = (n_sets_numeric > (double)(max_int));
    if(overflow) { throw std::runtime_error("jburden_wts_adapt: integer overflow for #sets"); }
    */

    // approximate?
    bool approx = (n_approx > 0) & ((double)n_approx < n_sets_numeric);
    
    // sets of (n_approx or n_sets) indices from the range [0; n]
    list<vector<int>> sets;
    if(approx) {
      // n_approx randomly sampled sets
      jburden_nchoosek_sample(n, i, n_approx, sets, gen);
      ind_approx.push_back(i);
    } else {
      // all possible sets
      jburden_nchoosek(n, i, sets);
      ind_exact.push_back(i);
    }
    int n_sets_comp = sets.size();

    // print info.
    if(verbose) {
      cout << " w[" << i << "]" << "; #sets = " << n_sets_numeric << 
        "; approx. = " << approx << 
        "; #set to compute = " << n_sets_comp << endl;
    }
    
    // B. Inner loop over the components in the sum:
    //  pnorm(V_alpha'^{-1} * pnorm(V_alpha;alpha') over |alpha| = i
    vector<double> comp(n_sets_comp);

    int j = 0;
    list<vector<int>>::iterator it = sets.begin(); 
    for(; it != sets.end(); ++it, ++j) {
      // indices: alpha
      vector<int> s2 = *it; 
      // indices: complement to alpha
      vector<int> s1; 
      _complement(n, s2, s1);
      
      // subset V: V11, V22, V12 & V21
      Eigen::MatrixXd V11 = _subset_matrix(V, s1, s1);
      Eigen::MatrixXd V22 = _subset_matrix(V, s2, s2);
      Eigen::MatrixXd V12 = _subset_matrix(V, s1, s2);
      Eigen::MatrixXd V21 = _subset_matrix(V, s2, s1);

      // matrix V220 = V22 - V21 solve(V11, V12)
      MatrixXd V11inv = _inverse(V11);
      MatrixXd V220 = V22 - V21 * V11inv * V12;

      double w_comp = jburden_pnorm(V11inv) * jburden_pnorm(V220);
      if(w_comp < 0) {
        std::cout << "ERROR: computing NNLS weight(i =" << i << ")" 
          << " & component(j = " << j << " out of " << n_sets_comp << ") failed; " 
          << "pnorm(V11inv)*pnorm(V220) returned negative value (" << w_comp << ")\n";
        return(-1);
      }
      comp.push_back(w_comp);
    }

    // sum of elements in comp
    double sum_comp = accumulate(comp.begin(), comp.end(), 0.0);

    // assign wts[i] using elements in comp
    if(approx) {
      // 1. Approximated weight
      // No. elements in comp = n_sets_comp
      // No. elements in total = n_sets
      // Approximation by sum of normals: 
      //   wts[i] = n_sets * mu, where mu = sum(comp) / n_sets_comp
      _assign_wts(wts, i, (sum_comp / n_sets_comp) * n_sets_numeric, verbose);
    } else {
      // 2. Exact weight, i.e., n_sets_comp = n_sets
      _assign_wts(wts, i, sum_comp, verbose);
    }
  }

  // debug
  if(verbose > 2) {
    cout << "jburden_wts_adapt: wts (before norm.) = ";
    for(vector<int>::iterator it = ind_approx.begin(); it < ind_approx.end(); ++it) {
      cout << " " << wts[*it] << "; ";
    }
    cout << endl;
  }

  // D. normalize in two steps
  // 1. sum(approximated weights) = 1 - sum(exact weights)
  //    - exact weights are not touched
  // 2. max weight = 1.0 - sum(all weights except max weight)
  //    - force the sum of all weights to be 1.0
  //      (it will be not exactly equal to 1.0 
  //       due to the limited numerical precision)
  if(normalize) {
    // D1. normalize approximate weights
    if(ind_approx.size()) {
      // compute scaling factor = sum_expected / sum_empirical =
      // = (1 - sum_exact) / sum_approx
      double sum_exact = 0.0;
      for(vector<int>::iterator it = ind_exact.begin(); it < ind_exact.end(); ++it) {
        sum_exact += wts[*it];
      }
      double sum_approx = 0.0;
      for(vector<int>::iterator it = ind_approx.begin(); it < ind_approx.end(); ++it) {
        sum_approx += wts[*it];
      }
      double sc_factor = (1 - sum_exact) / sum_approx;

      // update weights
      for(vector<int>::iterator it = ind_approx.begin(); it < ind_approx.end(); ++it) {
        wts[*it] *= sc_factor;
      }
    }

    // D2. normalize all weights
    // find the maximum weights
    vector<double>::iterator it_max = max_element(wts.begin(), wts.end());
    // assign [1 - (sum(wts) - max)] to the max element 
    //  - thus, sum(wts) = 1
    double sum_wts = accumulate(wts.begin(), wts.end(), 0.0);
    *it_max = 1.0 - (sum_wts - *it_max);
  }

  if(verbose > 2) {
    // debug
    cout << "jburden_wts_adapt: wts (after norm.) = ";
    for(vector<int>::iterator it = ind_approx.begin(); it < ind_approx.end(); ++it) {
      cout << " " << wts[*it] << "; ";
    }
    cout << endl;

    // debug
    double sum_wts = accumulate(wts.begin(), wts.end(), 0.0);
    cout << "jburden_wts_adapt: sum_wts = " << sum_wts << endl;
  }

  // return
  //  - returned weights are in the reverse order:
  //    w(n), w(n-1), ..., w(1), w(0)
  wts_out.resize(nw);
  for(int i = 0; i < nw; ++i) {
    wts_out[i] = wts[nw - 1 - i];
  }

  // return success
  return(0);
}

/***************************************
 * NLLS: the main function: 
 * - model fitting
 * - update V matrix & compute weigts
 * - compute p-values using chi2bar 
 * - return min(pval_pos_nnls, pval_neg_nnls)
***************************************/

Eigen::MatrixXd _npd(const Eigen::MatrixXd& X, double eps = 1e-6)
{
  // method 1: diag(X) = diag(X) + eps
  /* V = 0.5 * (V + V.transpose()); // force it symmetric */
  /* double eps = 1e-4; // to be a function argument */
  /* V.diagonal().array() += eps; // force positive-definite */
  
  // method 2: EVD (negative eigenvalues -> 0) + diag(X) = diag(X) + eps
  MatrixXd Y = 0.5 * (X + X.transpose());
  SelfAdjointEigenSolver<MatrixXd> solver(Y);
  VectorXd D = solver.eigenvalues();
  MatrixXd V = solver.eigenvectors();
  VectorXd Dplus = D.cwiseMax(0.0);
  MatrixXd Z = V * Dplus.asDiagonal() * V.transpose();
  Z.diagonal().array() += eps;

  return(Z);
}

// main function for the NNLS test
double jburden_test(const Eigen::VectorXd &y, const Eigen::MatrixXd& X, std::mt19937_64& gen,
  int df, double tol, int n_approx, bool strict, int verbose)
{
  // dimensions
  int n = y.size(), p = X.cols();
  if(df == 0) {
    df = n - p;
  }
  if(verbose) cout << " n = " << n << "; p = " << p << "; df = " << df << endl;

  // OLS solution
  if(verbose) cout << " - OLS" << endl;

  MatrixXd XtX(MatrixXd(p, p).setZero().  
    selfadjointView<Lower>().rankUpdate(X.adjoint()));
  VectorXd Xty = X.adjoint() * y;

  LLT<MatrixXd> llt(XtX);
  VectorXd bhat = llt.solve(Xty);

  VectorXd fitted = X * bhat;
  VectorXd resid = y - fitted;
  double sigma = resid.norm() / sqrt((double)(df));
  double sigma2 = sigma * sigma;
  /* MatrixXd V = sigma2 * llt.matrixL().solve(MatrixXd::Identity(p, p)); */
  MatrixXd V = sigma2 * llt.solve(MatrixXd::Identity(p, p));
  /* MatrixXd V = sigma2 * XtX.inverse(); */

  /* VectorXd Vib = (1.0 / sigma2) * XtX * bhat; */
  /* double stat = bhat.adjoint() * Vib; */
  /* double pval = jburden_pchisq(stat, p, false); */

  // weights for the NNLS test
  // - common for positive NNLS and negative NNLS
  if(verbose) cout << " - NNLS weights" << endl;

  // debug
  if(verbose > 2) {
    VectorXcd evals_complex = V.eigenvalues();
    VectorXd evals = evals_complex.real(); 
    cout << " - jburden_test: " << evals.size() << " eigenvalues = " << evals << endl;

    const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
    string name = "V.matirx.tmp";
    ofstream file(name.c_str());
    file << V.format(CSVFormat);
    file.close();
  }

  // make matrix V positive-definite for computing weights
  //  - not to be used for computing the test statistic
  //  - the significance drops a lot, e.g. LDL-PCSK9, 
  //    (i) -log10(P) = 94 for the corrected V;
  //    (ii) -log10(P) = 288 for the original V
  MatrixXd Vpd = V; //_npd(V);

  // exact weights
  /* VectorXd w = jburden_wts(V); */
  // adaptive weights: the key argument is n_approx
  VectorXd w;
  int ret = jburden_wts_adapt(Vpd, w, gen, n_approx, true, verbose);
  if(ret < 0) {
    if(strict) {
      std::cout << "ERROR: computing NNLS weights failed"
        << " by jburden_wts_adapt\n";
      exit(-1);
    } else {
      return(-1.0); // return p-value = -1
    }
  }
  
  // TEMP
  /* MatrixXd D = MatrixXd::Constant(p, p, 0.0); */
  /* D.diagonal().array() = 1.0; */
  /* VectorXd w = jburden_wts_adapt(D, n_approx, true, verbose); */

  // positive NNLS: b >= 0
  if(verbose) cout << " - NNLS b >= 0" << endl;
  
  VectorXd bhat_pos;
  vector<bool> selected_pos;
  int ret_pos = jburden_fit_nnls(y, X, bhat_pos, selected_pos, tol, false); // maxit, maxit_inner, verbose);
  if(ret_pos < 0) {
    if(strict) {
      std::cout << "ERROR: computing NNLS weights failed"
        << " by jburden_fit_nnls\n";
      exit(-1);
    } else {
      return(-1.0); // return p-value = -1
    }
  }
  
  VectorXd Vib_pos = (1.0 / sigma2) * XtX * bhat_pos;
  double stat_pos = bhat_pos.adjoint() * Vib_pos;
  double pval_pos = jburden_pchisq_bar(stat_pos, w);
  
  // negative NNLS: b <= 0
  if(verbose) cout << " - NNLS b <= 0" << endl;

  VectorXd bhat_neg;
  vector<bool> selected_neg;
  int ret_neg = jburden_fit_nnls(y, X, bhat_neg, selected_neg, tol, true); // maxit, maxit_inner, verbose);
  if(ret_neg < 0) {
    if(strict) {
      std::cout << "ERROR: computing NNLS weights failed"
        << " by jburden_fit_nnls\n";
      exit(-1);
    } else {
      return(-1.0); // return p-value = -1
    }
  }

  VectorXd Vib_neg = (1.0 / sigma2) * XtX * bhat_neg;
  double stat_neg = bhat_neg.adjoint() * Vib_neg;
  double pval_neg = jburden_pchisq_bar(stat_neg, w);

  // return
  double pval_min2 = min(pval_pos, pval_neg);

  return(pval_min2);
}


/***************************************
 * Methods of class NLLS
***************************************/

//-------------------
// NNLS constructors
//-------------------
NNLS::NNLS()
{
  // assign
  napprox = 10;
  normalize = true;
  tol = 1e-6;
  maxit = 1000;
  strict = false;
  verbose = 0;
  // defaults
  set_defaults();
}


NNLS::NNLS(int napprox_, bool normalize_, double tol_, bool strict_, int verbose_)
{
  // assign
  napprox = napprox_;
  normalize = normalize_;
  tol = tol_;
  maxit = 1000;
  maxit_inner = 500;
  strict = strict_;
  verbose = verbose_;
  // defaults
  set_defaults();
}

void NNLS::set_defaults()
{
  p = 0;
  nw = 0;
  msg_error = "";
  fit_pos.executed = false;
  fit_neg.executed = false;
  fit_pos.converged = false;
  fit_neg.converged = false;
  pval_min2 = -1.0;
}

//-------------------
// NNLS weights
//-------------------

void NNLS::compute_weights()
{
  if(verbose) cout << " NNLS: Weights step\n";

  // check previous steps
  if(p == 0) return;
  if(V.cols() != p) return;

  // assign Vpd
  Vpd = V; // _npd(V);

  VectorXd w;

  int ret = jburden_wts_adapt(Vpd, w, *gen, napprox, normalize, verbose);
  if(ret < 0) {
    msg_error = "error in computing NNLS weights";

    if(strict) {
      std::cout << "ERROR: computing NNLS weights failed"
        << " by jburden_wts_adapt\n";
      exit(-1);
    } else {
      return;
    }
  }

  // assign
  wts = w;   
  nw = wts.size();
}

//-------------------
// NNLS fit (X, y)
//-------------------

// OLS step
void NNLS::fit_ols(const Eigen::VectorXd &y, const Eigen::MatrixXd& X, 
  int df)
{
  if(verbose) cout << " NNLS: OLS step\n";

  int n = X.rows();
  p = X.cols();
  if(df == 0) {
    df = n - p;
  }

  MatrixXd XtX(MatrixXd(p, p).setZero().  
    selfadjointView<Lower>().rankUpdate(X.adjoint()));
  VectorXd Xty = X.adjoint() * y;

  LLT<MatrixXd> llt(XtX);
  VectorXd bhat = llt.solve(Xty);

  VectorXd fitted = X * bhat;
  VectorXd resid = y - fitted;
  double ss = resid.array().square().sum();

  // assign sigma2
  sigma2 = ss / df;
  // assign V
  /* V = sigma2 * llt.matrixL().solve(MatrixXd::Identity(p, p)); */
  V = sigma2 * llt.solve(MatrixXd::Identity(p, p));
  // assign XX
  XX = XtX;
  // assign bhat
  bhat_ols = bhat;

  // (optional) OLS test statistic
  VectorXd Vib = (1.0 / sigma2) * XX * bhat_ols;
  stat_ols = bhat_ols.adjoint() * Vib;
  /* double pval = jburden_pchisq(stat, p, false); */
}

// NNLS fit & inference step
void NNLS::pw_calc_pvals()
{
  // assign min p-value if both fits are ok
  if(fit_pos.converged & fit_neg.converged) {
    // re-calculate pval_pos
    VectorXd Vib = (1.0 / sigma2) * XX * fit_pos.bhat;
    fit_pos.stat = fit_pos.bhat.adjoint() * Vib;
    fit_pos.pval = jburden_pchisq_bar(fit_pos.stat, wts);
    // re-calculate pval_neg
    Vib = (1.0 / sigma2) * XX * fit_neg.bhat;
    fit_neg.stat = fit_neg.bhat.adjoint() * Vib;
    fit_neg.pval = jburden_pchisq_bar(fit_neg.stat, wts);
    // pval min2
    pval_min2 = min(fit_pos.pval, fit_neg.pval);
    best_fit = (fit_pos.pval < fit_neg.pval); // 1 = pos, 0 = neg
  }
}

void NNLS::fit_nnls(const Eigen::VectorXd &y, const Eigen::MatrixXd& X)
{
  if(verbose) cout << " NNLS: Fit step\n";

  // NNLS pos: b >= 0
  fit_nnls_sign(y, X, false, fit_pos);
  // NNLS neg: b <= 0
  fit_nnls_sign(y, X, true, fit_neg);

  // assign min p-value if both fits are ok
  if(fit_pos.converged & fit_neg.converged) {
    pval_min2 = min(fit_pos.pval, fit_neg.pval);
    best_fit = (fit_pos.pval < fit_neg.pval); // 1 = pos, 0 = neg
  }
}

void NNLS::fit_nnls_sign(const Eigen::VectorXd &y, const Eigen::MatrixXd& X, bool neg, struct FitNNLS& fit)
{
  // check previous steps
  if(p == 0) return;
  if(XX.cols() != p) return;
  if(nw == 0) return;
  if(wts.size() != p + 1) return;

  fit.executed = true;
  VectorXd bhat;
  vector<bool> selected;
  int ret = jburden_fit_nnls(y, X, bhat, selected, tol, neg, maxit, maxit_inner, verbose);
  if(ret < 0) {
    msg_error = "error in computing NNLS model fit";

    if(strict) {
      std::cout << "ERROR: computing NNLS model fit failed"
        << " by jburden_fit_nnls\n";
      exit(-1);
    } else {
      fit.converged = false;
      return;
    }
  }
 
  // assign results of model fitting
  fit.converged = true;
  fit.bhat = bhat;
  fit.selected = selected;

  // compute test statistic & assign
  VectorXd Vib = (1.0 / sigma2) * XX * bhat;
  fit.stat = bhat.adjoint() * Vib;
  fit.pval = jburden_pchisq_bar(fit.stat, wts);
}

//-------------------
// NNLS fit (b, V)
//-------------------

// NNLS fit & inference step
void NNLS::ss_fit_nnls()
{
  if(verbose) cout << " NNLS: Fit step\n";

  VectorXd Xty = Vinv * bhat_ols;
  MatrixXd XtX = Vinv;

  // NNLS pos: b >= 0
  ss_fit_nnls_sign(Xty, XtX, false, fit_pos);
  // NNLS neg: b <= 0
  ss_fit_nnls_sign(Xty, XtX, true, fit_neg);

  // assign min p-value if both fits are ok
  if(fit_pos.converged & fit_neg.converged) {
    pval_min2 = min(fit_pos.pval, fit_neg.pval);
    best_fit = (fit_pos.pval < fit_neg.pval); // 1 = pos, 0 = neg
  }
}

void NNLS::ss_fit_nnls_sign(const Eigen::VectorXd &Xty, const Eigen::MatrixXd& XtX, bool neg, struct FitNNLS& fit)
{
  // check previous steps
  if(p == 0) return;
  if(nw == 0) return;
  if(wts.size() != p + 1) return;

  fit.executed = true;
  VectorXd bhat;
  vector<bool> selected;
  int ret = jburden_fit_nnls_cprod(Xty, XtX, bhat, selected, tol, neg, maxit, maxit_inner, verbose);
  if(ret < 0) {
    msg_error = "error in computing NNLS model fit";

    if(strict) {
      std::cout << "ERROR: computing NNLS model fit failed"
        << " by jburden_fit_nnls_ss\n";
      exit(-1);
    } else {
      fit.converged = false;
      return;
    }
  }
 
  // assign results of model fitting
  fit.converged = true;
  fit.bhat = bhat;
  fit.selected = selected;

  // compute test statistic & assign
  fit.stat = bhat.adjoint() * Vinv * bhat;
  fit.pval = jburden_pchisq_bar(fit.stat, wts);
}


//-------------------
// NNLS main
//-------------------

void NNLS::pw_weights(int napprox_)
{
  // check dimensions
  if(V.rows() == 0) { throw std::runtime_error("pw_weights: dimensions (nrows = 0)"); }
  if(V.cols() == 0) { throw std::runtime_error("pw_weights: dimensions (ncols = 0)"); }
  if(V.rows() != V.cols()) { throw std::runtime_error("pw_weights: dimensions"); }

  // code copied from compute_weights() with input napprox_ instead of class member napprox
  VectorXd w;

  int ret = jburden_wts_adapt(V, w, *gen, napprox_, normalize, verbose);
  if(ret < 0) { 
    if(strict) { throw std::runtime_error("pw_weights: error in jburden_wts_adapt"); }
    else { return; }
  }

  // assign
  napprox = napprox_;
  wts = w;   
  nw = wts.size();
}

void NNLS::pw_weights(const Eigen::MatrixXd& V_)
{
  // check dimensions
  if(V_.rows() == 0) { throw std::runtime_error("pw_weights: dimensions (nrows = 0)"); }
  if(V_.cols() == 0) { throw std::runtime_error("pw_weights: dimensions (ncols = 0)"); }
  if(V_.rows() != V_.cols()) { throw std::runtime_error("pw_weights: dimensions"); }

  // code copied from compute_weights() with input V_ instead of class member V
  VectorXd w;

  int ret = jburden_wts_adapt(V_, w, *gen, napprox, normalize, verbose);
  if(ret < 0) { 
    if(strict) { throw std::runtime_error("pw_weights: error in jburden_wts_adapt"); }
    else { return; }
  }

  // assign
  wts = w;   
  nw = wts.size();
}

void NNLS::pw_weights(const Eigen::VectorXd& wts_)
{
  // check dimensions
  if(wts_.size() == 0) { throw std::runtime_error("pw_weights: dimensions (input weights size = 0)"); }

  // assign
  wts = wts_;   
  nw = wts.size();
}

void NNLS::ss_weights(const Eigen::MatrixXd& V_)
{
  if(verbose) print_param();

  // assign
  V = V_;
  Vinv = _inverse(V);
  p = V.rows();

  // check dimensions
  if(p == 0) { throw std::runtime_error("ss_set_V: dimensions (p = 0)"); }
  if(V.rows() != V.cols()) { throw std::runtime_error("ss_set_V: dimensions"); }

  compute_weights();
}


void NNLS::ss_run(const Eigen::VectorXd &bhat_)
{
  if(verbose) print_param();

  // assign
  bhat_ols = bhat_;
  p = bhat_ols.size();

  // check dimensions
  if(p == 0) { throw std::runtime_error("ss_run: dimensions (p = 0)"); }
  if(bhat_ols.size() != p) { throw std::runtime_error("ss_run: dimensions (bhat.size() != p)"); }
  if(V.rows() == 0) { throw std::runtime_error("ss_run: dimensions (V.rows() = 0)"); }
  if(Vinv.rows() == 0) { throw std::runtime_error("ss_run: dimensions (Vinv.rows() = 0)"); }
  if(V.rows() != p) { throw std::runtime_error("ss_run: dimensions"); }
  if(Vinv.rows() != p) { throw std::runtime_error("ss_run: dimensions"); }

  // check weights are pre-computed
  if(nw == 0) { throw std::runtime_error("ss_run: weights (nw == 0)"); }
  if(nw != (p + 1)) { throw std::runtime_error("ss_run: weights (nw != p + 1"); }
  if(wts.size() != nw) { throw std::runtime_error("ss_run: weights (wts.size() != nw"); }

  // fit & get p-values
  ss_fit_nnls(); 

  if(verbose) print_results();
}

void NNLS::ss_run(const Eigen::VectorXd &bhat_, const Eigen::MatrixXd& V_)
{
  if(verbose) print_param();

  // assign
  bhat_ols = bhat_;
  V = V_;
  Vinv = _inverse(V);
  p = bhat_ols.size();

  // check dimensions
  if(V.rows() != V.cols()) { throw std::runtime_error("ss_run: dimensions"); }
  if(bhat_ols.size() != V.cols()) { throw std::runtime_error("ss_run: dimensions"); }

  compute_weights();
  ss_fit_nnls(); 

  if(verbose) print_results();
}

void NNLS::pw_run(const Eigen::VectorXd &y, const Eigen::MatrixXd& X, int df)
{
  if(nw == 0) { throw std::runtime_error("pw_run: nw == 0"); }

  if(verbose) print_param();

  fit_ols(y, X, df);
  /* compute_weights(); */
  fit_nnls(y, X);

  if(verbose) print_results();
}

void NNLS::run(const Eigen::VectorXd &y, const Eigen::MatrixXd& X, int df)
{
  if(verbose) print_param();

  fit_ols(y, X, df);
  compute_weights();
  fit_nnls(y, X);

  if(verbose) print_results();
}
