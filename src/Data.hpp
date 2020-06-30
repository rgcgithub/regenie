/* 

   This file is part of the regenie software package.

   Copyright (c) 2020 Joelle Mbatchou & Jonathan Marchini

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

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <random>
#include <map>
#include <fstream>
#include <string>
#include <math.h>       /* exp */
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <thread>
#include <sys/types.h>
#include <sys/stat.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions.hpp>

#include "bgen_to_vcf.hpp"
#include "eigen3.3/Dense"
#include "eigen3.3/StdVector"

using namespace std;
using namespace Eigen;
using namespace boost;

#ifdef __linux__
#include <omp.h>
#endif

#define VERSION_NUMBER "1.0.2"

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

// for the log file
class mstream
{
  public:
    std::ofstream coss;

    template <class S>
      mstream& operator<< (const S& val)
      {
        coss << val;
        std::cout << val;
        return *this;
      };

    mstream& operator<< (ostream& (*pfun)(ostream&))
    {
      pfun(coss);
      pfun(cout);
      return *this;
    };

    mstream(void);
    ~mstream(void);
};

class MeasureTime {

  public:
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
    time_t start_time_info, end_time_info;

    void init() {
      auto start = chrono::system_clock::now(); // wall clock
      start_time_info = chrono::system_clock::to_time_t( start ); 
      begin = std::chrono::steady_clock::now(); // to measure elapsed time
    }

    void stop(){
      auto endtime = chrono::system_clock::now(); 
      end_time_info = chrono::system_clock::to_time_t( endtime ); 
      end = std::chrono::steady_clock::now();
    }

    MeasureTime(void);
    ~MeasureTime(void);
};


struct snp {
  int chrom;
  string ID;
  double genpos; 
  uint32_t physpos;
  string allele1, allele2;
  double MAF;
  uint64 offset;
  bool mask = false;
} ;

struct tally {
  uint64 snp_count = 0;
  uint64 n_failed_tests = 0;
  uint64 n_ignored_snps = 0;
  uint64 n_skipped_snps = 0;
};

struct variant_block {
  double af, info, scale_fac;
  double dif_deviance;
  double val_a, val_b, val_c, val_d; 
  ArrayXd Geno;
  ArrayXd Gmod;
  MatrixXd genocounts;
  vector<bool> test_fail;
  vector<bool> is_corrected; // for firth/spa
  ArrayXd scale_fac_pheno;
  ArrayXd denum;
  ArrayXd stats;
  ArrayXd chisq_val;
  ArrayXd pval_log;
  ArrayXd bhat;
  ArrayXd se_b;
  // for spa
  vector <uint> non_zero_indices;
  bool flipped, pos_score;
  // firth
  MatrixXd beta_null_firth;
  // reset each time
  bool ignored = false;
  bool fastSPA = true;
  uint n_non_zero = 0;
};


struct findID {
  uint64 index;
  bool is_found;
};


class Data {

public:
  // class elements
  string bgen_file = "NULL", sample_file = "NULL", cov_file = "NULL", pheno_file = "NULL", blup_file = "NULL";
  string bedfile="NULL", bimfile="NULL", famfile="NULL", out_file = "NULL", loco_tmp_prefix = "NULL", file_ind_include = "NULL", file_ind_exclude = "NULL", file_snps_include = "NULL", file_snps_exclude = "NULL";
  string err_help = "For list of arguments, run with option --help\n"; // for checks
  mstream sout;
  MeasureTime runtime;
  vector<string> pheno_names;
  bool select_phenos = false; // user specify which phenotype columns to use
  vector<string> pheno_colKeep_names;
  bool select_covs = false; // user specify which covariate columns to use
  vector<string> cov_colKeep_names;
  bool select_chrs = false; // user specify which chromosomes to test
  vector<int> chrKeep_test;
  int n_genofiles = 0; 
  bool bgenSample = false; // user specifies .sample file to go with bgen file
  bool rm_indivs = false; // user specify to remove genotyped samples from analysis
  bool keep_indivs = false; // user specify to keep only select genotyped samples in the analysis
  bool keep_snps = false; // user specify to keep select snps in analysis
  bool rm_snps = false; // user specify to remove snps from analysis
  int run_mode = 0; // running in null model fitting (=1) or association testing (=2)
  int nChrom = 23; // total number of chromosome numbers (sex chromosomes collapsed in chr23)
  int control_code = 0; // (default: 0=controls, 1=cases)
  int block_size = -1; // number of SNPs per block
  int cv_folds = 5; // number of CV folds
  int n_block = -1; // number of blocks to run
  int total_n_block = -1; // number of blocks to run across all chrs
  int n_ridge_l0 = 5; // number of ridge parameters at level 0
  int n_ridge_l1 = 5; // number of ridge parameters at level 1
  bool user_ridge_params_l0 = false; // if user specifies ridge parameters
  bool user_ridge_params_l1 = false; // if user specifies ridge parameters
  bool test_mode = false;
  int test_type = 0; // add=0/dom=1/rec=2 test
  string file_type = "null"; // type of the genotype file format;
  bool streamBGEN = true; // use fast version of step 2 when testing with BGEN v1.2 zlib compressed input
  int min_MAC = 5; // minimum MAC of SNPs in testing mode
  bool use_SPA = false; // saddlepoint approximation to estimate pvalue
  bool firth = false;// firth correction using LRT
  bool firth_approx = false; // approx. to Firth LRT
  bool binary_mode = false; // identifier: QT = false, BT = true
  bool within_sample_l0 = false; // specify to use within-sample predictions as features at level 1 (default is to use out-of-sample predictions)
  bool use_loocv = false; // specify whether to use LOOCV [note: this is not used if method=1]
  bool make_loco = true; // specify whether to compute & ouput LOCO predictions
  bool perchrLOCO = false; // write LOCO predictions for each chromosome
  bool write_l0_pred = false; // specify whether to write level 0 predictions to file to save on RAM
  bool split_by_pheno = false; // specify whether to write testing result in separate phenotype files
  bool htp_out = false; // use HTPv4 format
  string cohort_name = "NULL"; // Name of cohort to add in HTP output
  string model_type, test_string;
  vector < MatrixXd > genocounts;
  bool skip_blups = false;
  bool rm_missing_qt = true; // remove missing individuals when performing test with QTs
  bool strict_mode = false; // remove all missing individuals (at both step 1 and 2)
  bool print_block_betas = false; // print betas from level 0 within each block (for debugging)
  int print_snpcount = 0; // cumulative snp count across blocks
  const string missing_pheno_str = "NA";
  const double missing_value_double = -999;
  double eigen_val_rel_tol = 1e-16;
  int chunk_mb = 1000; // max amount of memory to use with LOOCV
  double numtol = 1e-6;
  double tol = 1e-8; // for logistic regression
  double eps = 10 * std::numeric_limits<double>::epsilon();
  double nl_dbl_dmin = 10.0 * std::numeric_limits<double>::min();
  int niter_max = 30; // max number of iterations for logistic reg.
  vector<double> lambda; // ridge parameters at level 0
  vector<double> tau; // ridge parameters at level 1
  vector<snp> snpinfo;
  vector<int> chr_counts, chr_file_counts;
  map<int, vector<int> > chr_map;
  map <std::string, uint64> FID_IID_to_ind;
  vector < string > snplist_to_keep, snplist_to_rm;
  int Nbed, Nbim;
  BgenParser bgen, bgen_tmp;
  MatrixXd covariates, new_cov;
  MatrixXd phenotypes, phenotypes_raw;
  ArrayXd ind_in_pheno_and_geno;
  ArrayXd ind_in_cov_and_geno;
  ArrayXd ind_in_analysis, ind_ignore;
  MatrixXd masked_indivs;
  vector < MatrixXd > masked_in_folds;
  ArrayXd pheno_l1_not_converged;
  ArrayXd bad_snps;
  MatrixXd blup;
  vector<Matrix<double, Dynamic, Dynamic> > Xt_Gamma_X_inv;
  MatrixXd Gamma_sqrt, Y_hat_p;
  ArrayXd Gmod;
  VectorXd denum_tstat;
  MatrixXd offset_logreg;
  MatrixXd G, G_tmp;
  int n_cov = 0; // number of covariates
  int n_samples = 0; // number of samples
  ArrayXd Neff; // number of non-missing samples (per trait)
  int n_pheno; // number of phenotypes
  uint32_t n_variants = 0; // number of variants in bgen file
  vector<int> cv_sizes;
  vector<Matrix<double, Dynamic, Dynamic> > G_folds, GtY; // storage for folds at levle 0
  vector<Matrix<double, Dynamic, Dynamic> > X_folds, XtY; // storage for folds at level 1
  MatrixXd GGt,GTY;
  MatrixXd GGt_eig_vec, GGt_eig_val, Wmat;
  vector<vector<Matrix<double, Dynamic, Dynamic> > > pred_mat, test_mat;
  vector<Matrix<double, Dynamic, Dynamic> > test_mat_conc;
  vector<vector<Matrix<double, Dynamic, Dynamic> > > pred_pheno, test_pheno;
  vector<vector<Matrix<double, Dynamic, Dynamic> > > pred_pheno_raw, test_pheno_raw;
  vector<vector<Matrix<double, Dynamic, Dynamic> > > pred_offset, test_offset;
  vector<vector<Matrix<double, Dynamic, Dynamic> > > beta_hat_level_1;
  vector<Matrix<double, Dynamic, Dynamic> > predictions;
  bool verbose = false;
  vector<Matrix<double, Dynamic, Dynamic> > cumsum_values; // storage of values to compute rsq and values [Sx, Sy, Sx2, Sy2, Sxy]
  vector<Matrix<double, Dynamic, Dynamic> >  beta_print_out;
  MatrixXd ident_l0, ident_l1;
  ifstream bed_ifstream;
  vector <uchar> inbed;
  vector<string> blup_files;
  vector<int> pheno_index;
  uint threads = 0;
  Matrix<double, Dynamic, Dynamic> blups;
  MatrixXd res, stats, W_hat, snp_afs, snp_info;
  RowVectorXd scale_Y;
  RowVectorXd p_sd_yres;
  VectorXd scale_G; // keep track of sd(Y) (1xP) and sd(G) (M*1)
  double alpha_pvalue = 0.05; // significance threshold above which to use firth correction
  int n_corrected = 0; // to keep track of how many SNPs require correction
  double deviance_logistic = -1, bhat_firth, se_b_firth;
  double numtol_firth = 1e-5; // tolerance level for firth
  int niter_max_firth = 250; // max number of iterations in Firth logistic reg.
  int niter_max_firth_null = 1000; // max number of iterations in Firth logistic reg. null model
  int niter_max_line_search = 5; // max number of iterations for line search in logistic reg.
  bool fix_maxstep_null = false; // if user specifies max step size
  int maxstep = 5; // max step size in penalized logistic regression
  int maxstep_null = 25; // max step size in null penalized logistic regression
  int retry_maxstep_firth=5, retry_niter_firth=5000; // fallback settings for null approx. firth regression
  int pval_converged = false; // keep track of whether SPA/Firth converged
  vector < vector < uint > > non_zero_indices_G;
  MatrixXd covs_firth, beta_null_firth;
  bool fastSPA; // use fast approx. for rare SNPs
  int niter_max_spa = 1000; 
  double tol_spa = pow( std::numeric_limits<double>::epsilon(), 0.25);
  MatrixXd SPA_pvals;
  double val_a, val_b, val_c, val_d;
  bool pos_score;
  vector<bool> snp_flipped;
  uint32_t size1, size2, snp_index_counter, total_chrs_loco;
  uint64 bed_block_size; // prevent overflow
  
  // function definitions
  void run();
  void print_help(bool);
  void start_log(int,char **);
  void print_header(ostream&);
  void read_params_and_check(int argc,char *argv[]);
  void set_ridge_params(int,vector<double>&);
  void print_usage_info();
  void file_read_initialization();
  void prep_bgen();
  void read_bgen_sample(std::vector<string>&);
  void read_bed_bim_fam();
  void read_bim();
  void read_fam();
  void prep_bed(string);
  void set_snps_to_rm();
  void set_snps_to_keep();
  void set_IDs_to_rm();
  void set_IDs_to_keep();
  void pheno_read();
  void covariate_read();
  void getCovBasis();
  double convertDouble(const string&);
  findID getIndivIndex(const string&,const string&);
  void residualize_phenotypes();
  void residualize_genotypes();
  void fit_null_logistic(int);
  void set_blocks();
  void set_folds();
  void setmem();
  void get_G(int,int,int);
  void readChunkFromBGENFileToG(int,int);
  void readChunkFromBedFileToG(int);
  void scale_genotypes(bool);
  void level_0_calculations();
  void ridge_level_0(int);
  void ridge_level_0_loocv(int);
  void ridge_level_1();
  void ridge_level_1_loocv();
  void ridge_logistic_level_1();
  void ridge_logistic_level_1_loocv();
  void calc_cv_matrices(int);
  void read_pheno_and_cov();
  void make_predictions(int,int);
  void make_predictions_loocv(int,int);
  void make_predictions_binary(int,int);
  void make_predictions_binary_loocv(int,int);
  double compute_log_lik(double,double);
  void output();
  void write_predictions(int);

  //void set_test_files();
  void blup_read();
  void blup_read_chr(int);
  void check_bgen();
  void test_snps();
  void set_blocks_for_testing();
  double check_pval(double,int,int,int);
  double run_firth_correction(int,int,int);
  bool fit_firth_logistic(int,int,bool);
  void run_SPA_test(int);
  double get_SPA_pvalue(double,double,bool,bool,int,int);
  double solve_K1(double,bool,bool,int,int);
  double compute_K(double,int,int);
  double compute_K1(double,int,int);
  double compute_K2(double,int,int);
  double compute_K_fast(double,int,int);
  double compute_K1_fast(double,int,int);
  double compute_K2_fast(double,int,int);
  void skip_snps(int);
  // fast streaming functions
  void test_snps_fast();
  void analyze_block(const int&,const int&,tally*,vector<variant_block>&);
  void get_data_blocks(std::istream*,vector<uchar>*);
  void extract_variant_MT(vector<uchar>*,const uint32_t,const uint32_t,snp*,variant_block*);
  void readChunkFromBedFileToG(const int &n_snps, vector<variant_block> &all_snps_info);
  void residualize_geno(variant_block*);
  void run_SPA_test_snp(variant_block*,int,const VectorXd&);
  double get_SPA_pvalue_snp(double,double,variant_block*,int);
  double solve_K1_snp(double,variant_block*,int);
  double compute_K_snp(double,variant_block*,int);
  double compute_K1_snp(double,variant_block*,int);
  double compute_K2_snp(double,variant_block*,int);
  double compute_K_fast_snp(double,variant_block*,int);
  double compute_K1_fast_snp(double,variant_block*,int);
  double compute_K2_fast_snp(double,variant_block*,int);
  void check_pval_snp(variant_block*,int,int);
  void run_firth_correction_snp(int,variant_block*,int);
  void fit_firth_logistic_snp(int,int,bool,variant_block*);


  int chrStrToInt(string);

  Data();
  ~Data();




};



