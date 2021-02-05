#ifndef _MVT_H_
#define _MVT_H_

#ifdef __cplusplus
extern "C" {

  extern void mvtdst_(int* n, int* nu, double* lower, double* upper, int* infin,
      double* correl, double* delta, int* maxpts, double* abseps, double* releps,
      double* error, double* value, int* inform);

#endif

  double pmvnorm(int* n, int* nu, double* lower, double* upper,
      int* infin,  double* correl,  double* delta,   
      int* maxpts, double* abseps,  double* releps,  
      double* error, double* value, int* inform);    

  double pmvnorm_complement(int n, int maxpts, double abseps, double* bound,
      double* cmat, double* error);

#ifdef __cplusplus
}
#endif
#endif /* _MVT_H_ */
