
#include <stdlib.h>
#include <stdio.h>
#include "mvtnorm.h"

// infinity bounds
const static int INFIN_BOUND_NORMAL = 2;        // (..., ...)
const static int INFIN_BOUND_UPPER = 1;         // (..., inf)
const static int INFIN_BOUND_LOWER = 0;         // (-inf, ..)
const static int INFIN_BOUND_LOWER_UPPER = -1;  // (-inf, inf)

// returns the probability between 0 and 1 
// returns a negative value if errors
// -1 = completed with error > abseps
// -2 = n greater 1000 or n < 1
// -3 = matrix not positive semidefinite
double pmvnorm(int* n, int* nu, double* lower, double* upper,
    int* infin,  double* correl,  double* delta,   
    int* maxpts, double* abseps,  double* releps,  
    double* error, double* value, int* inform)
{
  mvtdst_ (n, nu,
      lower, upper, infin, correl, delta,
      maxpts, abseps, releps,
      error, value, inform);
  /* printf ("error = %g, value = %g, inform = %d, abseps = %g \n", *error, *value, *inform, *abseps); */

  switch (*inform) {
    case 0:
      return *value;
    case 1: 
      return -1.0;
    case 2:
      return -2.0;
    case 3:
      return -3.0;
  };

  return *value;
};

// The complement CDF of X ~ MVN(0, Correlation Matrix)
double pmvnorm_complement(int n,
    int maxpts, double abseps,
    double* bound, // all zeros if CDF = P(X <= 0)
    double* cmat, // correlation matrix with only lower-diagonal entires stored: (2,1), (3,1), (3,2)...
    double* error)
{
  int nu_ = 0;
  int maxpts_ = maxpts;
  double abseps_ = abseps;
  double releps_ = 0;

  double* upper = new double[n];
  int* infin = new int[n];
  double* delta = new double[n];

  for (size_t i = 0; i < n; i++) {
    infin[i] = 1; // (-inf, bound]
    upper[i] = 0.0;
    delta[i] = 0.0;
  }

  // values to return
  double value_ = 0.0;
  int inform_ = 0.0;

  double ret = pmvnorm(&n, &nu_, 
      bound, upper, infin, cmat, delta, 
      &maxpts_, &abseps_, &releps_, error, &value_, &inform_);

  delete[] (upper);
  delete[] (infin);
  delete[] (delta);

  return ret;
}
