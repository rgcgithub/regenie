/* 

   This file is part of the regenie software package.

   Copyright (c) 2020-2022 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

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

#ifndef REGENIE_H
#define REGENIE_H


#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <random>
#include <memory>
#include <map>
#include <fstream>
#include <math.h>       /* exp */
#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <sys/types.h>
#include <sys/stat.h>

// if using external LAPACK routines
#ifdef WITH_OPENBLAS
// fix conflict between complex and older boost versions
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include "lapacke.h"
#elif defined(WITH_MKL)
#include "mkl_lapacke.h"
#endif

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wparentheses"
#endif
#include <boost/math/distributions.hpp>

#include "Eigen/Dense"
#include "Eigen/StdVector"
#include <Eigen/SparseCore>

#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#ifdef __linux__
#include <omp.h>
#endif

#define MAXFILELEN 2001

#define BIT_SET(a,b) ((a) |= (1ULL<<(b)))
#define BIT_UNSET(a,b) ((a) &= ~(1ULL << (b)))
#define CHECK_BIT(a,b) ((a) & (1ULL<<(b)))

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;
typedef Eigen::Array<bool,Eigen::Dynamic,1> ArrayXb;
typedef Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> MatrixXb;
typedef Eigen::Map<Eigen::ArrayXd > MapArXd;
typedef Eigen::Map<const Eigen::ArrayXd > MapcArXd;
typedef Eigen::Map<Eigen::MatrixXd > MapMatXd;
typedef Eigen::Map<const Eigen::MatrixXd > MapcMatXd;
typedef Eigen::Map<ArrayXb> MapArXb;
typedef Eigen::Map<const ArrayXb> MapcArXb;
typedef Eigen::Array<uint16_t,Eigen::Dynamic,1> ArrayXt;
typedef Eigen::Array<uint64,Eigen::Dynamic,1> ArrayXui;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseMatrix<bool> SpMatb;
typedef Eigen::SparseVector<double> SpVec;

inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

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
      }

    mstream& operator<< (std::ostream& (*pfun)(std::ostream&))
    {
      pfun(coss);
      pfun(std::cout);
      return *this;
    };

    mstream(void);
    ~mstream(void);
};


class MeasureTime {

  public:
    std::chrono::steady_clock::time_point begin, end;
    std::chrono::high_resolution_clock::time_point ms_begin;
    time_t start_time_info, end_time_info;

    void init() {
      auto start = std::chrono::system_clock::now(); // wall clock
      start_time_info = std::chrono::system_clock::to_time_t( start ); 
      begin = std::chrono::steady_clock::now(); // to measure elapsed time
    }

    void stop(){
      auto endtime = std::chrono::system_clock::now(); 
      end_time_info = std::chrono::system_clock::to_time_t( endtime ); 
      end = std::chrono::steady_clock::now();
    }

    void start_ms() {
      ms_begin = std::chrono::high_resolution_clock::now(); // wall clock
    }

    std::string stop_ms(){
      auto ms_end = std::chrono::high_resolution_clock::now(); // wall clock
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(ms_end - ms_begin);
      std::ostringstream buffer;
      buffer << "done (" << duration.count() << "ms)";
      return buffer.str();
    }

    MeasureTime(void);
    ~MeasureTime(void);
};


struct param {

  std::string err_help = "For list of arguments, run with option --help\n"; // for checks
  std::string webinfo = "For more information, use option '--help' or visit the website: https://rgcgithub.github.io/regenie/"; 

  //////
  // global options
  int run_mode; // running in null model fitting (=1) or association testing (=2)
  bool test_mode = false; // step 1: false; step 2 = true
  int trait_mode = 0; // 0=QT,1=BT,2=CT
  bool strict_mode = false; // remove individuals with any NA
  bool bgenSample = false; // .sample file for bgen file
  bool gzOut = false; // to compress output files (.loco and .regenie files)
  bool transposedPheno = false, tpheno_iid_only = false;
  bool getCorMat = false, cor_out_txt = false, cormat_force_vars = false;
  std::vector<std::string> forced_in_snps;//variant to force in for LD matrix
  std::map<std::string, int> extract_vars_order;//order of variants
  bool condition_snps = false, condition_file = false;
  uint32_t max_condition_vars = 10000;
  int sex_specific = 0; // 0 = all; 1 = male-only; 2=female-only

  // filters 
  bool rm_indivs = false; // user specify to remove genotyped samples from analysis
  bool keep_indivs = false; // user specify to keep only select genotyped samples in the analysis
  bool keep_snps = false, keep_or = false; // user specify to keep select snps in analysis
  bool rm_snps = false, rm_or = false; // user specify to remove snps from analysis
  bool mk_snp_map = false, keep_snp_map = false;
  bool select_phenos = false, select_phenos_rm = false, force_qt_run = false; // user specify which phenotype columns to use
  bool select_covs = false, select_covs_rm = false, cat_cov = false; // user specify which covariate columns to use and if categorical covars present
  int max_cat_levels = 10; // maximum number of categories of categorical covars
  bool select_chrs = false; // user specify which chromosomes to test

  // other global options
  const std::string missing_pheno_str = "NA";
  const double missing_value_double = -999;
  int nChrom = 23; // total number of chromosome numbers (sex chromosomes collapsed in chr23)
  bool CC_ZeroOne = true; // BT: 0/1 encoding?
  int mcc = 10; // minimum case count
  double numtol = 1e-6, qr_tol = 1e-7;
  double numtol_firth = 1e-4; // tolerance level for firth
  double numtol_eps = 10 * std::numeric_limits<double>::epsilon();
  double tol = 1e-8; // for logistic regression
  double eigen_val_rel_tol = 1e-15;
  double nl_dbl_dmin = 10.0 * std::numeric_limits<double>::min();
  int threads = 0, neff_threads = 1;
  bool verbose = false, debug = false;
  bool early_exit = false, l1_full_samples = false, rint = false;
  bool split_l0 = false, run_l0_only = false, run_l1_only = false; // run level 0 in parallel across different jobs
  std::map<std::string, bool> select_pheno_l1;
  int njobs, job_num, parallel_nGeno, parallel_nBlocks, parallel_nSnps;
  int start_block = 1;
  bool use_adam = false, adam_mini = true; // use ADAM for log. reg.
  double adam_alpha = 0.001, adam_beta1 = 0.9, adam_beta2 = 0.999, adam_eps = 1e-7, adam_batch_size = 128;
  std::vector<Eigen::ArrayXi> adam_indices;

  // for input data
  uint32_t n_samples = 0, n_analyzed = 0; // number of samples
  int n_pheno = 0; // number of phenotypes
  int n_cov = 0; // number of covariates
  int ncov, ncov_analyzed, ncov_interaction; // number of linearly independent covariates
  uint32_t n_variants = 0, nvs_stored = 0; // number of variants in genotype file
  std::map <std::string, uint32_t> FID_IID_to_ind;
  std::vector< std::vector<std::string> > FIDvec; // store FID/IID separately (for write-samples option)
  bool with_bgi = false, zlib_compress; // input bgi index file for BGEN format and compression format
  uint BGENbits = 0; // bit-encoding used in BGEN file
  bool ref_first = false; // ordering of REF/ALT alleles in input genotype file
  Eigen::ArrayXi sex; // 0=unknown, 1=male, 2=female
  std::vector<Eigen::ArrayXd> bed_lookup_table; // plink bed lookup table
  ArrayXb pheno_pass;

  // step 1 
  int block_size = -1; // number of SNPs per block
  int cv_folds = 5; // number of CV folds
  int n_block = -1; // number of blocks to run
  int total_n_block = 0; // number of blocks to run across all chrs
  int n_ridge_l0 = 5; // number of ridge parameters at level 0
  int n_ridge_l1 = 5; // number of ridge parameters at level 1
  double alpha_prior = -1; // to set MAF dependent prior on the effect sizes
  int chunk_mb = 1000; // max amount of memory to use with LOOCV
  bool user_ridge_params_l0 = false; // if user specifies ridge parameters
  bool user_ridge_params_l1 = false; // if user specifies ridge parameters
  bool use_loocv = false; // specify whether to use LOOCV [note: this is not used if method=1]
  bool make_loco = true; // specify whether to compute & ouput LOCO predictions
  bool print_prs = false; // specify to print PRS (i.e. no LOCO used)
  bool write_blups = false; // write BLUP predictions for each chromosome
  bool use_rel_path = false; // write relative paths in pred.list file
  bool write_l0_pred = false; // specify whether to write level 0 predictions to file to save on RAM
  bool rm_l0_pred = true; // specify whether to delete written level 0 predictions after level 1
  bool print_block_betas = false; // print betas from level 0 within each block (for debugging)
  bool force_run = false; // if using more than max nvariants in step 1
  int max_step1_variants = 1e6; // prevent users using too many step 1 variants
  int niter_max_ridge = 100, niter_max_ridge_adam = 25; // max number of iterations for ridge logistic reg.
  int niter_max_line_search_ridge = 100; // max number of iterations for line search in ridge logistic reg.
  double l1_ridge_tol = 1e-4; // tolerance level for convergence criteria
  double l1_ridge_eps = 1e-5; // epsilon used to set weights for 0/1 probabilities
  uint32_t print_snpcount = 0; 
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >  beta_print_out;
  Eigen::ArrayXd lambda; // ridge parameters at level 0
  std::vector<Eigen::ArrayXd> tau; // ridge parameters at level 1
  // TO REMOVE
  bool within_sample_l0 = false; // specify to use within-sample predictions as features at level 1 (default is to use out-of-sample predictions)
  Eigen::ArrayXi cv_sizes;


  // step 2
  bool rm_missing_qt = true; // remove missing individuals when performing test with QTs
  std::string file_type; // type of the genotype file format;
  bool streamBGEN = true; //  for BGEN v1.2 with 8-bit encoding
  bool fastMode = true; // use fast version of step 2 
  bool dosage_mode = false; // track if dosages are present for step 2
  bool split_by_pheno = true; // specify whether to write testing result in separate phenotype files
  bool skip_blups = false, blup_cov = false;
  bool with_flip = true; // can flip to minor allele for all variants
  bool use_prs = false; // adjust for whole genome PRS (no LOCO)
  double min_MAC = 5, min_MAC_mask, minHOMs = 0; // minimum MAC of SNPs in testing mode
  bool setMinMAC = false;
  double min_INFO = 0; // minimum INFO score of SNPs (dosages) in testing mode
  bool setMinINFO = false;
  bool write_samples = false; // write sample IDs for each trait
  double alpha_pvalue = 0.05, zcrit, z_thr, chisq_thr; // significance threshold above which to use firth correction
  int test_type = 0; // add=0/dom=1/rec=2 test
  bool w_interaction = false, interaction_cat = false, interaction_snp = false, interaction_prs = false, interaction_file = false, w_ltco = false, print_vcov = false, hlm_vquad = true, int_add_extra_term = false, int_add_esq = false, add_homdev = false; // interaction test
  int interaction_istart = 0, ltco_chr;
  uint64 interaction_snp_offset; // index in genotype file
  bool force_robust = false, force_hc4 = false, no_robust = false; // when using robust SE for rare variants with QTs
  double rareMAC_inter = 1000; // MAC below which to use HLM
  int n_tests_per_variant = 1;
  std::vector<std::string> interaction_lvl_names; // name of levels if using categorical variable for test
  bool gwas_condtl = true;
  std::string condtl_suff;
  // spa
  bool use_SPA = false; // saddlepoint approximation to estimate pvalue
  int niter_max_spa = 1000; 
  double tol_spa = pow( std::numeric_limits<double>::epsilon(), 0.25);
  // firth
  bool firth = false;// firth correction using LRT
  bool firth_approx = false; // approx. to Firth LRT
  bool write_null_firth = false, use_null_firth = false, compute_all_chr = false; // write/use null coefficients from approx. Firth
  int niter_max = 50; // max number of iterations for logistic reg.
  int niter_max_firth = 250, niter_max_firth_adam = 25; // max number of iterations in Firth logistic reg.
  int niter_max_firth_null = 1000; // max number of iterations in Firth logistic reg. null model
  int niter_max_line_search = 25; // max number of iterations for line search in logistic reg.
  int maxstep = 5; // max step size in penalized logistic regression
  int maxstep_null = 25; // max step size in null penalized logistic regression
  bool fix_maxstep_null = false; // if user specifies max step size
  bool back_correct_se = false; // for SE with Firth
  bool print_pheno_name = false; // add phenotype name when writing to file with sample IDs
  bool htp_out = false, af_cc = false; 
  std::string cohort_name; // Name of cohort to add in HTP output
  bool set_range = false;
  int range_chr; 
  double range_min, range_max; // use genomic region to filter variants
  std::string build_code = "hg38"; // to identify chrX PAR region bounds
  uint32_t par1_max_bound, par2_min_bound;

  // snp sets (masks/joint tests)
  bool snp_set = false; 
  bool build_mask = false; 
  std::map<std::string, uint> vc_tests_map = { {"acatv", 0}, {"skat", 1}, {"skato", 2}, {"skato-acat", 3}, {"acato", 4} };
  uint vc_test = 0;
  bool apply_gene_pval_strategy = false;
  std::string genep_mask_sets_file = "";
  std::map<std::string, bool> mask_map;
  double vc_maxAAF = 1; // max AAF for variants in SKAT/ACAT gene-based tests
  bool w_anno_lab = false;
  bool check_mask_files = false, strict_check_burden = false, fail_check = false;
  bool skip_test = false; // skip computing tests
  bool joint_test = false; // for step 2 joint testing
  std::string burden = ""; // type of burden test;
  uint max_set_size = 1000; // maximum number of predictors in joint test
  bool set_select_list = false; // comma separated list of sets given
  bool keep_sets = false, rm_sets = false; // user specify to filter sets in analysis
  bool w_regions = false; // categorize by set regions 
  uint max_cat = 64, nmax_regions = 16; // maximum number of annotations (to fit in uint64)
  std::vector<std::string> mbins; // temporary object to store aaf bins
  bool mask_rule_max = true, mask_rule_comphet = false; // default use max to combine mask
  std::string mask_rule = "max";
  bool set_aaf = false;// for user-given AAFs for building masks
  bool aaf_file_wSingletons = false;//for choosing snps in singleton masks
  bool singleton_carriers = false; // carrier count used to define singletons
  uint64 max_bsize = 0; // number of SNPs per variant set
  bool write_masks = false, write_setlist = false, write_mask_snplist = false; //write masks to bed file
  bool mask_loo = false, mask_lodo = false;
  bool use_max_bsize = false; // set bsize to max set size
  bool p_joint_only = false;
  std::string mask_loo_name, mask_loo_set, mask_loo_region, masks_loo_snpfile; // for LOO with masks
  double mask_loo_aaf;
  bool nnls_out_all = false;
  int nnls_napprox = 10;
  double acat_a1 = 1, acat_a2 = 25, skat_a1 = 1, skat_a2 = 25, skat_tol = 1e-5; // for ACAT & SKAT test
  int skat_collapse_MAC = 10;
  Eigen::ArrayXd skato_rho; // rho parameter from skat-o

  // multi-trait tests 
  bool trait_set = false; 
};

struct geno_file_info {
  std::string file, format;
  bool with_sample = false, with_bgi = false, ref_first = false;
  std::string sample;
};

// for input files
struct in_files {

  std::string bed_prefix;
  std::string pgen_prefix;
  std::string bgen_file, sample_file;
  std::vector<std::string> file_ind_include, file_ind_exclude;
  std::vector<std::string> file_snps_include, file_snps_exclude;
  std::vector<std::string> file_snps_include_or, file_snps_exclude_or;
  std::vector<std::string> file_sets_include, file_sets_exclude;
  std::string sets_include, sets_exclude;
  std::string cov_file, pheno_file;
  std::string loco_tmp_prefix = "";
  std::string split_file;
  std::string out_file;
  std::string blup_list_file;
  std::string null_firth_file;
  std::vector<std::shared_ptr<std::ofstream>> write_preds_files;
  std::map<std::string, std::string> blup_files;
  std::vector<std::string> null_firth_files;
  std::vector<std::string> pheno_names;
  std::vector<int> chr_counts, chr_read;
  uint64 bed_block_size; // prevent overflow
  std::ifstream geno_ifstream;
  std::vector<uchar> inbed;
  std::vector<std::vector<uchar>> bed_data_blocks;
  std::string set_file, new_sets;
  std::string anno_file, anno_labs_file, mask_file, aaf_file;
  std::vector<int> bstart, btot; // for parallel l0
  std::vector<std::string> mprefix; // for parallel l0
  std::string condition_snps_list; // for conditional analyses
  geno_file_info condition_snps_info; 
  geno_file_info interaction_snp_info; 

};

struct filter {

  // to filter phenotype/covariates/genotype
  std::map<std::string, bool> pheno_colKeep_names, pheno_colRm_names, cov_colKeep_names, cov_colRm_names; //cov keep map: true for qVar, false for catVar
  std::map<int, bool> tpheno_colrm;
  uint32_t tpheno_indexCol;
  std::string interaction_cov;
  std::string interaction_cov_null_level;//if categorical for GxE / or coding for GxG
  std::map <int, bool> chrKeep_test;
  std::map <std::string, uint32_t> snpID_to_ind;
  ArrayXb ind_ignore, has_missing, ind_in_analysis;
  uint32_t step1_snp_count = 0;
  std::map <std::string, std::vector<int>> setID_to_ind;//chr,index,is_kept
  std::map <std::string, uint64> condition_snp_names;

};

void start_log(const std::string&,MeasureTime*,mstream&);
template <typename T> 
void print_args(T,std::map<std::string,bool>&,mstream&);

void print_help(bool const&);
void read_params_and_check(int& argc,char *argv[],struct param*,struct in_files*,struct filter*,MeasureTime*,mstream&);
void check_file(std::string const&,std::string const&);
void check_file(std::string const&,std::vector<std::string> const&,std::string const&);
void print_header(std::ostream&);
Eigen::ArrayXd get_unit_params(bool const&,std::string const&,std::vector<std::string> const&,struct param const*,mstream&);
void set_ridge_params(int const&,Eigen::ArrayXd&,mstream&);
void print_usage_info(struct param const*,struct in_files*,mstream&);
int chrStrToInt(const std::string&, const int&);
std::vector<std::string> check_name(std::string const&,mstream&);
void check_build_code(struct param*);
double convertDouble(const std::string&,struct param const*,mstream&);
double convertNumLevel(const std::string&,std::map<std::string,int>&,struct param const*,mstream&);
void check_inter_var(std::string&,std::string&,mstream&);
std::string print_csv(const std::vector<std::string>&);
std::string print_scsv(const std::vector<std::string>&);
template <typename T>
std::string print_sv(const std::vector<T>&,const std::string&);
Eigen::ArrayXi get_true_indices(const Eigen::Ref<const ArrayXb>&);
void get_logp(double&,const double&);
void get_logp(const double&,double&,double&,const double&);
void allocate_mat(Eigen::MatrixXd&,int const&,int const&);
std::string print_mat_dims(Eigen::MatrixXd const&);
int parseLine(char*);
int get_mem();
std::string print_mem();

void set_threads(struct param*);

template <typename KEY, typename VALUE> 
bool in_map(KEY element, std::map<KEY,VALUE> const& emap){
  return emap.find(element) != emap.end();
}

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

template <typename U> 
std::string get_test_list(U const& bit_map, std::map<std::string, U>& srt_map){

  std::vector<std::string> test_list;
  typename std::map <std::string, U>::iterator itr;
  for (itr = srt_map.begin(); itr !=  srt_map.end(); ++itr) {
    // skip nnls_pos and neg
    if(itr->first == "nnls_pos" || itr->first == "nnls_neg") continue;
    if( CHECK_BIT(bit_map, itr->second) ) { // add to test list
      std::string newstr = itr->first;
      std::transform(newstr.begin(), newstr.end(), newstr.begin(), ::toupper);
      test_list.push_back( newstr );
    }
  }
  return print_csv(test_list);
}



#endif
