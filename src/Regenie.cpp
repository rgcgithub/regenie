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

#include "cxxopts.hpp"
#include <regex>
#include <chrono>
#include "Regenie.hpp"
#include "Geno.hpp"
#include "Joint_Tests.hpp"
#include "Step1_Models.hpp"
#include "Step2_Models.hpp"
#include "Files.hpp"
#include "Pheno.hpp"
#include "Masks.hpp"
#include "Data.hpp"


using namespace std;
using namespace Eigen;
using namespace boost;


mstream::mstream(){ }
mstream::~mstream(){ }
MeasureTime::MeasureTime(){ }
MeasureTime::~MeasureTime(){ }


int main( int argc, char** argv ) {

  Data data;
  read_params_and_check(argc, argv, &data.params, &data.files, &data.in_filters, &data.runtime, data.sout);

  data.run();

  data.runtime.stop();

  data.sout << "\nElapsed time : " << std::chrono::duration<double>(data.runtime.end - data.runtime.begin).count() << "s" << endl;
  data.sout << "End time: " << ctime(&data.runtime.end_time_info) << endl; 

}


void print_header(std::ostream& o){

  std::ostringstream oss;
  string vnumber;
  // adjust spacing for version with Boost Iostream library (`.gz` suffix)
#ifndef HAS_BOOST_IOSTREAM
  oss << "  ";
#endif
  oss << "REGENIE v" << VERSION_NUMBER; 
  vnumber = oss.str();

  o << left << std::setw(14) << " " << "|" << std::string(30, '=')<< "|" << endl;
  o << left << std::setw(14) << " " << "|" << left << std::setw(6) << " " <<
    left << std::setw(24) << vnumber << "|" << endl;
  o << left << std::setw(14) << " " << "|" << std::string(30, '=')<< "|\n\n";

  o << "Copyright (c) 2020 Joelle Mbatchou and Jonathan Marchini." << endl;
  o << "Distributed under the MIT License.\n";
#ifdef HAS_BOOST_IOSTREAM
  o << "Compiled with Boost Iostream library.\n";
#endif

  // adding BLAS/LAPACK external routines
#if defined(WITH_MKL)
  o << "Using Intel MKL with Eigen.\n";
#elif defined(WITH_OPENBLAS)
  o << "Using BLAS/LAPACK routines from OpenBLAS with Eigen.\n";
#endif

  o << "\n";
}


void read_params_and_check(int argc, char *argv[], struct param* params, struct in_files* files, struct filter* filters, MeasureTime* mt, mstream& sout) {

  cxxopts::Options AllOptions(argv[0], "");

  AllOptions.add_options()
    ("h,help", "print list of available options")
    ("helpFull", "print list of all available options")
    ;

  // add main options
  AllOptions.add_options("Main")
    ("step", "specify if fitting null model (=1) or association testing (=2)", cxxopts::value<int>(params->run_mode),"INT")
    ("bed", "prefix to PLINK .bed/.bim/.fam files", cxxopts::value<std::string>(files->bed_prefix),"PREFIX")
    ("pgen", "prefix to PLINK2 .pgen/.pvar/.psam files", cxxopts::value<std::string>(files->pgen_prefix),"PREFIX")
    ("bgen", "BGEN file", cxxopts::value<std::string>(files->bgen_file),"FILE")
    ("sample", "sample file corresponding to BGEN file", cxxopts::value<std::string>(files->sample_file),"FILE")
    ("ref-first", "use the first allele as the reference for BGEN or PLINK bed/bim/fam input format [default assumes reference is last]")
    ("keep", "file listing samples to retain in the analysis (no header; starts with FID IID)", cxxopts::value<std::string>(files->file_ind_include),"FILE")
    ("remove", "file listing samples to remove from the analysis (no header; starts with FID IID)", cxxopts::value<std::string>(files->file_ind_exclude),"FILE")
    ("extract", "file with IDs of variants to retain in the analysis", cxxopts::value<std::string>(files->file_snps_include),"FILE")
    ("exclude", "file with IDs of variants to remove from the analysis", cxxopts::value<std::string>(files->file_snps_exclude),"FILE")
    ("p,phenoFile", "phenotype file (header required starting with FID IID)", cxxopts::value<std::string>(files->pheno_file),"FILE")
    ("phenoCol", "phenotype name in header (use for each phenotype to keep)", cxxopts::value< std::vector<std::string> >(filters->pheno_colKeep_names),"STRING")
    ("phenoColList", "comma separated list of phenotype names to keep", cxxopts::value<std::string>(),"STRING,..,STRING")
    ("c,covarFile", "covariate file (header required starting with FID IID)", cxxopts::value<std::string>(files->cov_file),"FILE")
    ("covarCol", "covariate name in header (use for each covariate to keep)", cxxopts::value< std::vector<std::string> >(filters->cov_colKeep_names),"STRING")
    ("covarColList", "comma separated list of covariate names to keep", cxxopts::value<std::string>(),"STRING,..,STRING")
    ("o,out", "prefix for output files", cxxopts::value<std::string>(files->out_file),"PREFIX")
    ("bt", "analyze phenotypes as binary")
    ("1,cc12", "use control=1,case=2,missing=NA encoding for binary traits")
    ("b,bsize", "size of genotype blocks", cxxopts::value<int>(params->block_size),"INT")
    ("cv", "number of cross validation (CV) folds", cxxopts::value<int>(params->cv_folds),"INT(=5)")
    ("loocv", "use leave-one out cross validation (LOOCV)")
    ("l0", "number of ridge parameters to use when fitting models within blocks [evenly spaced in (0,1)]", cxxopts::value<int>(params->n_ridge_l0),"INT(=5)")
    ("l1", "number of ridge parameters to use when fitting model across blocks [evenly spaced in (0,1)]", cxxopts::value<int>(params->n_ridge_l1),"INT(=5)")
    ("lowmem", "reduce memory usage by writing level 0 predictions to temporary files")
    ("lowmem-prefix", "prefix where to write the temporary files in step 1 (default is to use prefix from --out)", cxxopts::value<std::string>(files->loco_tmp_prefix),"PREFIX")
    ("split-l0", "split level 0 across N jobs and set prefix of output files", cxxopts::value<std::string>(),"PREFIX,N")
    ("run-l0", "run level 0 for job K in {1..N} using master file created from '--split-l0'", cxxopts::value<std::string>(),"FILE,K")
    ("run-l1", "run level 1 using master file from '--split-l0'", cxxopts::value<std::string>(files->split_file),"FILE")
    ("keep-l0", "avoid deleting the level 0 predictions written on disk after fitting the level 1 models")
    ("strict", "remove all samples with missingness at any of the traits")
    ("print-prs", "also output polygenic predictions without using LOCO (=whole genome PRS)")
    ("gz", "compress output files (gzip format)")
    ("threads", "number of threads", cxxopts::value<int>(params->threads),"INT")
    ("pred", "file containing the list of predictions files from step 1", cxxopts::value<std::string>(files->blup_file),"FILE")
    ("ignore-pred", "skip reading predictions from step 1 (equivalent to linear/logistic regression with only covariates)")
    ("use-prs", "when using whole genome PRS step 1 output in '--pred'")
    ("write-samples", "write IDs of samples included for each trait (only in step 2)")
    ("minMAC", "minimum minor allele count (MAC) for tested variants", cxxopts::value<double>(params->min_MAC),"FLOAT(=5)")
    ("minINFO", "minimum imputation info score (Impute/Mach R^2) for tested variants", cxxopts::value<double>(params->min_INFO),"DOUBLE(=0)")
    ("split", "split asssociation results into separate files for each trait")
    ("firth", "use Firth correction for p-values less than threshold")
    ("approx", "use approximation to Firth correction for computational speedup")
    ("spa", "use Saddlepoint approximation (SPA) for p-values less than threshold")
    ("pThresh", "P-value threshold below which to apply Firth/SPA correction", cxxopts::value<double>(params->alpha_pvalue),"FLOAT(=0.05)")
    ("chr", "specify chromosome to test in step 2 (use for each chromosome)", cxxopts::value< std::vector<std::string> >(),"STRING")
    ("chrList", "Comma separated list of chromosomes to test in step 2", cxxopts::value<std::string>(),"STRING,..,STRING")
    ("range", "to specify a physical position window for variants to test in step 2", cxxopts::value<std::string>(),"CHR:MINPOS-MAXPOS")
    ("test", "'dominant' or 'recessive' (default is additive test)", cxxopts::value<std::string>(),"STRING")
    ("set-list", "file with sets definition", cxxopts::value<std::string>(files->set_file),"FILE")
    ("extract-sets", "file with IDs of sets to retain in the analysis", cxxopts::value<std::string>(files->file_sets_include),"FILE")
    ("exclude-sets", "file with IDs of sets to remove from the analysis", cxxopts::value<std::string>(files->file_sets_exclude),"FILE")
    ("extract-setlist", "comma separated list of sets to retain in the analysis", cxxopts::value<std::string>(files->file_sets_include),"FILE")
    ("exclude-setlist", "comma separated list of sets to remove from the analysis", cxxopts::value<std::string>(files->file_sets_exclude),"FILE")
    ("anno-file", "file with variant annotations", cxxopts::value<std::string>(files->anno_file),"FILE")
    ("anno-labels", "file with labels to annotations", cxxopts::value<std::string>(files->anno_labs_file),"FILE")
    ("mask-def", "file with mask definitions", cxxopts::value<std::string>(files->mask_file),"FILE")
    ("aaf-file", "file with AAF to use when building masks", cxxopts::value<std::string>(files->aaf_file),"FILE")
    ("aaf-bins", "comma separated list of AAF bins cutoffs for building masks", cxxopts::value<std::string>(),"FLOAT,..,FLOAT")
    ("build-mask", "rule to construct masks, can be 'max', 'sum' or 'comphet' (default is max)", cxxopts::value<std::string>(params->mask_rule),"STRING")
    ("singleton-carrier", "define singletons as variants with a single carrier in the sample")
    ("write-mask", "write masks in PLINK bed/bim/fam format")
    ("mask-lovo", "apply Leave-One-Variant-Out (LOVO) scheme when building masks (<set_name>,<mask_name>,<aaf_cutoff>)", cxxopts::value<std::string>(),"STRING")
    ("skip-test", "skip computing association tests after building masks")
    ;


  // extended options
  AllOptions.add_options("Additional")
    ("v,verbose", "verbose screen output")
    ("setl0", "comma separated list of ridge parameters to use when fitting models within blocks", cxxopts::value<std::string>(), "FLOAT,..,FLOAT")
    ("setl1", "comma separated list of ridge parameters to use when fitting model across blocks", cxxopts::value<std::string>(), "FLOAT,..,FLOAT")
    ("nauto", "number of autosomal chromosomes", cxxopts::value<int>(),"INT")
    ("nb", "number of blocks to use", cxxopts::value<int>(params->n_block),"INT")
    ("force-step1", "run step 1 for more than 1M variants (not recommended)")
    ("niter", "maximum number of iterations for logistic regression", cxxopts::value<int>(params->niter_max),"INT(=30)")
    ("maxstep-null", "maximum step size in null Firth logistic regression", cxxopts::value<int>(params->maxstep_null),"INT(=25)")
    ("maxiter-null", "maximum number of iterations in null Firth logistic regression", cxxopts::value<int>(params->niter_max_firth_null),"INT(=1000)")
    ("force-impute", "keep and impute missing observations when in step 2 (default is to drop missing for each trait)")
    ("firth-se", "Compute SE for Firth based on effect size estimate and LRT p-value")
    ("print-pheno", "Print phenotype name when writing sample IDs to file (only for step 2)")
    ;

  // extra options
  AllOptions.add_options("Extra")
    ("print", "print estimated effect sizes from level 0 and level 1 models")
    //("nostream", "print estimated effect sizes from level 0 and level 1 models")
    ("htp", "output association files in step 2 in HTPv4 format", cxxopts::value<std::string>(params->cohort_name),"STRING")
    ("within", "use within-sample predictions as input when fitting model across blocks in step 1")
    ("early-exit", "Exit program after fitting level 0 models (avoid deleting temporary prediction files from level 0)")
    ("prior-alpha", "alpha value used when speifying the MAF-dependent prior on SNP effect sizes", cxxopts::value<double>(params->alpha_prior),"FLOAT(=-1)")
    ("joint", "comma spearated list of joint tests to perform", cxxopts::value<std::string>(params->burden),"STRING")
    ("joint-only", "only output p-values from joint tests")
    ("write-setlist", "file with list of masks to combine as sets", cxxopts::value<std::string>(files->new_sets),"FILE")
    ("nnls-napprox", "number of random draws to use for approximate NNLS test", cxxopts::value<int>(params->nnls_napprox),"INT(=10)")
    ("nnls-verbose", "To output detailed NNLS test results")
    ("acat-beta", "parameters for Beta(a,b) used for ACAT test statistic", cxxopts::value<std::string>(), "a,b(=1,1)")
    ;


  try
  {
    //AllOptions.parse_positional({"htp"});
    auto vm = AllOptions.parse(argc, argv);
    auto arguments = vm.arguments();

    // help menu
    if (vm.count("help")){
      print_header(std::cout);
      std::cout << AllOptions.help({"", "Main"}) << '\n' << params->webinfo << "\n\n";
      exit(EXIT_SUCCESS);
    } else if (vm.count("helpFull")) {
      print_header(std::cout);
      std::cout << AllOptions.help({"", "Main", "Additional"}) << '\n' << params->webinfo << "\n\n";
      exit(EXIT_SUCCESS);
    } 
    
    if (!vm.count("out")){
      print_header(std::cout);
      std::cout << "ERROR: You must provide an output prefix using '--out'" << '\n' << params->webinfo << "\n\n";
      exit(EXIT_FAILURE);
    }


    // Print output to file and to stdout
    // print command line arguments
    start_log(arguments, files->out_file, mt, sout);
    vector< string > tmp_str_vec;

    if( (vm.count("bgen") + vm.count("bed")  + vm.count("pgen"))  != 1 ){
      sout << "ERROR: You must use either --bed,--bgen or --pgen.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }

    if( vm.count("bgen") ) params->file_type = "bgen";
    if( vm.count("bed") ) params->file_type = "bed";
    if( vm.count("pgen") ) params->file_type = "pgen";
    if( vm.count("sample") ) params->bgenSample = true;
    if( vm.count("ref-first") ) params->ref_first = true;
    if( vm.count("keep") ) params->keep_indivs = true;
    if( vm.count("remove") ) params->rm_indivs = true;
    if( vm.count("extract") ) params->keep_snps = true;
    if( vm.count("exclude") ) params->rm_snps = true;
    if( vm.count("phenoCol") ) params->select_phenos = true;
    if( vm.count("covarCol") ) params->select_covs = true;
    if( vm.count("bt") ) params->binary_mode = true;
    if( vm.count("1") ) params->CC_ZeroOne = false;
    if( vm.count("loocv") ) params->use_loocv = true;
    if( vm.count("strict") ) params->strict_mode = true;
    if( vm.count("print-prs") ) params->print_prs = true;
    if( vm.count("ignore-pred") ) params->skip_blups = true;
    if( vm.count("use-prs") ) params->use_prs = true;
    if( vm.count("force-impute") ) params->rm_missing_qt = false;
    if( vm.count("split") ) params->split_by_pheno = true;
    if( vm.count("approx") ) params->firth_approx = true;
    if( vm.count("nauto") ) params->nChrom = vm["nauto"].as<int>() + 1;
    if( vm.count("maxstep-null") | vm.count("maxiter-null") ) params->fix_maxstep_null = true;
    if( vm.count("firth") ) params->firth = true;
    if( vm.count("spa") ) params->use_SPA = true;
    if( vm.count("minMAC") ) params->setMinMAC = true;
    if( vm.count("minINFO") ) params->setMinINFO = true;
    if( vm.count("htp") ) params->htp_out = params->split_by_pheno = true;
    if( vm.count("v") ) params->verbose = true;
    if( vm.count("range") ) params->set_range = true;
    if( vm.count("print") ) params->print_block_betas = true;
    //if( vm.count("nostream") ) params->streamBGEN = params->fastMode = false;
    if( vm.count("within") ) params->within_sample_l0 = true;
    if( vm.count("write-samples") ) params->write_samples = true;
    if( vm.count("print-pheno") ) params->print_pheno_name = true;
    if( vm.count("early-exit") ) params->early_exit = true;
    if( vm.count("force-step1") ) params->force_run = true;
    if( vm.count("lowmem") ) params->write_l0_pred = true;
    if( vm.count("keep-l0") ) params->rm_l0_pred = false;
    if( vm.count("split-l0") ) params->split_l0 = true;
    if( vm.count("run-l0") ) { params->run_l0_only = params->write_l0_pred = params->keep_snps = true; params->rm_snps = false;}
    if( vm.count("run-l1") ) params->run_l1_only = params->write_l0_pred = true;
    if( vm.count("firth") && vm.count("firth-se") ) params->back_correct_se = true;
    if( vm.count("joint") ) params->joint_test = true;
    if( vm.count("joint-only") ) params->p_joint_only = true;
    if( vm.count("extract-sets") ) params->keep_sets = true;
    if( vm.count("exclude-sets") ) params->rm_sets = true;
    if( vm.count("extract-setlist") ) params->set_select_list = params->keep_sets = true;
    if( vm.count("exclude-setlist") ) params->set_select_list = params->rm_sets = true;
    if( vm.count("aaf-file") ) params->set_aaf = true;
    if( vm.count("singleton-carrier") ) params->singleton_carriers = true;
    if( vm.count("mask-lovo") ) params->mask_loo = true;
    if( vm.count("write-mask") ) params->write_masks = true;
    if( vm.count("write-setlist") ) params->write_setlist = true;
    if( vm.count("skip-test") ) params->skip_test = true;
    if( vm.count("nnls-verbose") ) params->nnls_out_all = true;
    if( vm.count("gz") ) {
# if defined(HAS_BOOST_IOSTREAM)
      // only works when compiled with boost IO library
      params->gzOut = true;
# else
      sout << "WARNING : REGENIE was not compiled with Boost Iostream library so ignoring option '--gz'.\n";
#endif
    }


    if( vm.count("phenoColList") ) {
      params->select_phenos = true;
      tmp_str_vec = string_split(vm["phenoColList"].as<string>(),",");
      filters->pheno_colKeep_names.insert( 
          filters->pheno_colKeep_names.end(),
          std::begin( tmp_str_vec ), 
          std::end( tmp_str_vec)         );
    }
    if( vm.count("covarColList") ) {
      params->select_covs = true;
      tmp_str_vec = string_split(vm["covarColList"].as<string>(),",");
      filters->cov_colKeep_names.insert( 
          filters->cov_colKeep_names.end(),
          std::begin( tmp_str_vec ), 
          std::end( tmp_str_vec)         );
    }
    if( vm.count("chrList") ) {
      params->select_chrs = true;
      tmp_str_vec = string_split(vm["chrList"].as<string>(),",");
      for( size_t ichr = 0; ichr < tmp_str_vec.size(); ichr++)
        filters->chrKeep_test.insert( std::make_pair( chrStrToInt(tmp_str_vec[ichr], params->nChrom), true ) );
    }
    if( vm.count("chr") ) {
      params->select_chrs = true;
      tmp_str_vec = vm["chr"].as<std::vector<string>>();
      for( size_t ichr = 0; ichr < tmp_str_vec.size(); ichr++)
        filters->chrKeep_test.insert( std::make_pair( chrStrToInt(tmp_str_vec[ichr], params->nChrom), true ) );
    }
    if( vm.count("split-l0") ) { // Format: FILE,INT
      tmp_str_vec = string_split(vm["split-l0"].as<string>(),",");
      if(tmp_str_vec.size() != 2 ){
        sout << "ERROR: Wrong format for --split-l0 (must be FILE,INT).\n" << params->err_help;
        exit(EXIT_FAILURE);
      }
      files->split_file = tmp_str_vec[0];
      params->njobs = atoi( tmp_str_vec[1].c_str() );
    }
    if( vm.count("run-l0") ) { // Format: FILE,INT
      tmp_str_vec = string_split(vm["run-l0"].as<string>(),",");
      if(tmp_str_vec.size() != 2 ){
        sout << "ERROR: Wrong format for --run-l0 (must be FILE,INT).\n" << params->err_help;
        exit(EXIT_FAILURE);
      }
      files->split_file = tmp_str_vec[0];
      params->job_num = atoi( tmp_str_vec[1].c_str() );
      if(params->job_num < 1 ){
        sout << "ERROR: Invalid job number for --run-l0 (must be >=1).\n" << params->err_help;
        exit(EXIT_FAILURE);
      }
    }
    if( vm.count("test") ) {
      if( vm["test"].as<string>() == "dominant") params->test_type = 1; 
      else if( vm["test"].as<string>() == "recessive") params->test_type = 2; 
      else {
        sout << "ERROR: Unrecognized argument for option --test, must be either 'dominant' or 'recessive'.\n" << params->err_help;
        exit(EXIT_FAILURE);
      }
    }
    if( vm.count("range") ) { // Format: Chr:min-max
      char tmp_chr[20];
      double p0 = -1, p1 = -1;
      string tmpd = vm["range"].as<string>();

      if(sscanf( tmpd.c_str(), "%[^:]:%lf-%lf", tmp_chr, &p0, &p1 ) != 3
          || (p0 < 0) || (p1 < 0) ){ 
        //cerr << tmp_chr << "\t" << p0 << "\t" << p1 << endl;
        sout << "ERROR: Wrong format for --range (must be CHR:MINPOS-MAXPOS).\n" << params->err_help;
        exit(EXIT_FAILURE);
      }

      tmpd = tmp_chr;
      params->range_chr = chrStrToInt(tmpd, params->nChrom);
      params->range_min = min(p0,p1);
      params->range_max = max(p0,p1);
    }

    if( vm.count("build-mask") ) {
      if( params->mask_rule == "max") params->mask_rule_max = true; 
      else if( params->mask_rule == "sum") params->mask_rule_max = false; 
      else if( params->mask_rule == "comphet") { params->mask_rule_max = false, params->mask_rule_comphet = true; 
      } else {
        sout << "ERROR : Unrecognized argument for option --build-mask (=" << params->mask_rule << ").\n" << params->err_help;
        exit(EXIT_FAILURE);
      }
    }

    if( vm.count("acat-beta") ) {
      tmp_str_vec = string_split(vm["acat-beta"].as<string>(),",");
      params->acat_a1 = convertDouble( tmp_str_vec[0], params, sout);
      params->acat_a2 = convertDouble( tmp_str_vec[1], params, sout);
    }

    if ( params->run_mode == 1 ) params->test_mode = false;
    else if (params->run_mode == 2 ) params->test_mode = true;
    else {
      sout << "ERROR: Specify which mode regenie should be running using option --step.\n" << params->err_help;
      exit(EXIT_FAILURE);
    }

    if(!params->test_mode) {

      // loocv only used with out-of-sample predictions
      if(params->use_loocv && params->within_sample_l0) {
        sout << "WARNING : Option --loocv cannot be used with option --within.\n" ;
        params->use_loocv = false;
      }

      // writing of level 0 predictions only available when using out-of-sample predictions
      if(params->write_l0_pred && params->within_sample_l0){
        sout << "WARNING : Option --lowmem cannot be used with option --within.\n" ;
        params->write_l0_pred = false;
      }

      // user specified ridge parameters to use at l0
      if( vm.count("setl0") ) {
        params->user_ridge_params_l0 = true;
        tmp_str_vec = string_split(vm["setl0"].as<string>(),",");
        for( size_t val = 0; val < tmp_str_vec.size(); val++)
          params->lambda.push_back(convertDouble( tmp_str_vec[val], params, sout));
        std::sort(params->lambda.begin(), params->lambda.end());
        params->lambda.erase( unique( params->lambda.begin(), params->lambda.end() ), params->lambda.end() );
        params->n_ridge_l0 = params->lambda.size();
        // parameters must be less in (0, 1)
        if( std::count_if(params->lambda.begin(), params->lambda.end(), std::bind2nd(std::greater<double>(), 0)) != params->n_ridge_l0 || std::count_if(params->lambda.begin(), params->lambda.end(), std::bind2nd(std::less<double>(), 1)) != params->n_ridge_l0 ){
          sout << "ERROR: You must specify values for --l0 in (0,1).\n" << params->err_help;
          exit(EXIT_FAILURE);
        } 
      } else set_ridge_params(params->n_ridge_l0, params->lambda, params->err_help, sout);

      // user specified ridge parameters to use at l1
      if( vm.count("setl1") ) {
        params->user_ridge_params_l1 = true;
        tmp_str_vec = string_split(vm["setl1"].as<string>(),",");
        for( size_t val = 0; val < tmp_str_vec.size(); val++)
          params->tau.push_back(convertDouble( tmp_str_vec[val], params, sout));
        std::sort(params->tau.begin(), params->tau.end());
        params->tau.erase( unique( params->tau.begin(), params->tau.end() ), params->tau.end() );
        params->n_ridge_l1 = params->tau.size();
        if( std::count_if(params->tau.begin(), params->tau.end(), std::bind2nd(std::greater<double>(), 0)) != params->n_ridge_l1 || std::count_if(params->tau.begin(), params->tau.end(), std::bind2nd(std::less<double>(), 1)) != params->n_ridge_l1 || (params->n_ridge_l1 == 0) ){
          sout << "ERROR: You must specify values for --l1 in (0,1).\n" << params->err_help;
          exit(EXIT_FAILURE);
        }
      } else set_ridge_params(params->n_ridge_l1, params->tau, params->err_help, sout);

      // firth only done in test mode
      if(params->firth) params->firth = false;
      if(params->use_SPA) params->use_SPA = false;
      params->streamBGEN = false;
      params->test_type = 0;
      if( vm.count("range") ) {
        params->set_range = false; 
        sout << "WARNING : Option --range only works for step 2.\n";
      }

    } else if(params->firth && !params->binary_mode) {
      // firth correction is only applied to binary traits
      sout << "WARNING : Option --firth will not be applied (it is only run with binary traits).\n";
      params->firth = false;
    } else if(params->use_SPA && !params->binary_mode) {
      // SPA is only applied to binary traits
      sout << "WARNING : Option --spa will not be applied (it is only run with binary traits).\n";
      params->use_SPA = false;
    }


    if(params->test_mode && params->use_loocv) params->use_loocv = false;

    if( (vm.count("write-samples") || vm.count("write-mask")) && vm.count("bgen") && !vm.count("sample") ){
      sout << "ERROR : must specify sample file (using --sample) if writing sample IDs to file.\n" << params->err_help;
      exit(EXIT_FAILURE);
    }


    if( vm.count("joint") ){

      if( vm.count("test") ){ 
        sout << "ERROR : cannot use --test with --joint.\n" << params->err_help;
        exit(EXIT_FAILURE);
      } else if ( vm.count("nnls-napprox") && params->nnls_napprox < 1 ){
        sout << "ERROR : must pass positive integer for --nnls-napprox\n";
        exit(EXIT_FAILURE);
      }

      params->snp_set = true;
    }

    if( vm.count("anno-file") || vm.count("mask-def") ){

      if(vm.count("anno-labels")) params->w_anno_lab = true;

      if( !(vm.count("anno-file") && vm.count("mask-def")) ){
        sout << "ERROR : must use --anno-file with --mask-def.\n" << params->err_help;
        exit(EXIT_FAILURE);
      } 

      if( (params->test_type > 0) && !params->mask_rule_max && !params->mask_rule_comphet ){
        sout << "ERROR : only additive test allowed when using 'sum' in --build-mask.\n" << params->err_help;
        exit(EXIT_FAILURE);
      } 

      if(params->write_masks && !params->mask_rule_max && !params->mask_rule_comphet ){
        sout << "ERROR : Cannot write masks when using 'sum' in --build-mask.\n";
        exit(EXIT_FAILURE);
      }


      // store aaf bins if given
      if( vm.count("aaf-bins") ) 
        tmp_str_vec = string_split(vm["aaf-bins"].as<string>(),",");
      else tmp_str_vec.resize(0);
      params->mbins = tmp_str_vec;


      if( vm.count("mask-lovo") ) {
        int cstart = 1;
        tmp_str_vec = string_split(vm["mask-lovo"].as<string>(),",");
        if(tmp_str_vec.size() < 3){
          cerr << "ERROR: Wrong format for option --mask-lovo.\n";
          exit(EXIT_FAILURE);
        } else if ( tmp_str_vec.size() == 4 ) {
          params->w_regions = true; cstart++;
        }
        params->mask_loo_set = tmp_str_vec[0];
        if(params->w_regions) params->mask_loo_region = tmp_str_vec[cstart-1];
        params->mask_loo_name = tmp_str_vec[cstart];
        params->mbins.resize(1);
        params->mbins[0] = tmp_str_vec[cstart+1]; // either singleton or AAF cutoff

        if(params->write_masks){
          sout << "WARNING : Cannot use --write-mask with --mask-lovo.\n";
          params->write_masks = false;
        }
      }

      params->snp_set = true;
      params->build_mask = true;

    }

    if(!params->build_mask && params->write_masks) params->write_masks = false;
    if(!params->write_masks && params->skip_test) params->skip_test = false;
    if(!params->write_masks && params->write_setlist) {
      sout << "WARNING : Must use --write-setlist with --write-mask.\n";
      params->write_setlist = false;
    }

    if( params->snp_set && !vm.count("set-list") ){
      sout << "ERROR : must specify set list (using --set-list).\n" << params->err_help;
      exit(EXIT_FAILURE);
    }

    if( params->snp_set && 
        (vm.count("extract-sets")+vm.count("exclude-sets")+vm.count("extract-setlist")+vm.count("exclude-setlist"))>1 
      ){
      sout << "ERROR :Must use only one of --extract-sets/--exclude-sets/--extract-setlist/--exclude-setlist.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    } 

    if(!params->test_mode && params->setMinMAC){
      sout << "WARNING : Option --minMAC only works in step 2 of REGENIE.\n";
      params->setMinMAC = false;
    }
    if(params->test_mode && params->min_MAC < 0.5){
      sout << "ERROR: minimum MAC must be at least 0.5.\n" << params->err_help;
      exit(EXIT_FAILURE);
    }
    if(!params->test_mode && params->setMinINFO){
      sout << "WARNING : Option --minINFO only works in step 2 of REGENIE.\n";
      params->setMinINFO = false;
    }
    if(params->test_mode && (params->min_INFO < 0 || params->min_INFO > 1) ){
      sout << "ERROR: minimum info score must be in [0,1].\n" << params->err_help;
      exit(EXIT_FAILURE);
    }
    if( params->rm_missing_qt && (params->strict_mode || params->binary_mode || !params->test_mode) ) params->rm_missing_qt = false;

    if( !vm.count("bsize") && !params->snp_set ) {
      sout << "ERROR : must specify the block size using '--bsize'.\n" << params->err_help;
      exit(EXIT_FAILURE);
    } else if(vm.count("bsize") && ( params->block_size < 2 )){
      sout << "ERROR : Block size must be at least 2.\n" << params->err_help;
      exit(EXIT_FAILURE);
    }
    if(params->set_aaf && !params->build_mask) params->set_aaf = false;

    // determine number of threads if not specified
    if(params->threads < 1)
      params->threads = std::max(1u, std::thread::hardware_concurrency() - 1); //may return 0 when not able to detect

    // check parallel l0
    if(params->test_mode && 
        (vm.count("split-l0")||vm.count("run-l0")||vm.count("run-l1")) ) {
      sout << "WARNING : Options --split-l0/--run-l0/--run-l1 only work in step 1.\n";
      params->split_l0 = params->run_l0_only = params->run_l1_only = false;
    } else if( vm.count("nb") && 
        (vm.count("split-l0")||vm.count("run-l0")||vm.count("run-l1")) ) {
      sout << "WARNING : Options --split-l0/--run-l0/--run-l1 cannot be used with --nb.\n";
      params->split_l0 = params->run_l0_only = params->run_l1_only = false;
    }
    if( (vm.count("run-l0")||vm.count("run-l1")) && !file_exists(files->split_file) ){
      sout << "ERROR: " << files->split_file  << " doesn't exist.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }

    // set Firth as default if both Firth and SPA are specified
    if(params->use_SPA && params->firth) {
      sout << "WARNING : Only one of --firth/--spa can be used. Only Firth will be used.\n";
      params->use_SPA = false;
    }

    // check firth fallback pvalue threshold
    if(params->firth && ((params->alpha_pvalue < params->nl_dbl_dmin) || (params->alpha_pvalue > 1 - params->numtol)) ){
      sout << "ERROR: Firth fallback p-value threshold must be in (0,1).\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    // check SPA fallback pvalue threshold
    if(params->use_SPA && ((params->alpha_pvalue < params->nl_dbl_dmin) || (params->alpha_pvalue > 1 - params->numtol)) ){
      sout << "ERROR: SPA fallback p-value threshold must be in (0,1).\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(params->firth_approx && !params->firth) params->firth_approx = false;

    // check arguments for logistic regression 
    if(params->binary_mode && (params->niter_max < 1)){
      sout << "ERROR: Invalid argument for --niter (must be positive integer).\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(params->firth && (params->maxstep_null < 1)){
      sout << "ERROR: Invalid argument for --maxstep-null (must be a positive integer).\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(params->firth && (params->niter_max_firth_null < 1)){
      sout << "ERROR: Invalid argument for --maxiter-null (must be a positive integer).\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(params->nChrom < 2){
      sout << "ERROR: Invalid argument for --nauto (must be > 1).\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(params->set_range && (params->range_chr == -1)){
      sout << "ERROR: Unrecognized chromosome in --range.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(params->rm_indivs && params->keep_indivs ){
      sout << "ERROR: Cannot use both --keep and --remove.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(params->rm_snps && params->keep_snps ){
      sout << "ERROR: Cannot use both --extract and --exclude.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }

    if( params->test_mode && params->select_chrs && (filters->chrKeep_test.find(-1) != filters->chrKeep_test.end()) ){
      sout << "ERROR: Invalid chromosome specified by --chr/--chrList.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }

    if(params->test_mode && !params->skip_blups && !vm.count("pred")) {
      sout << "ERROR: You must specify --pred if using --step 2 (otherwise use --ignore-pred).\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }

    if(params->test_mode && (params->file_type == "pgen") && !params->fastMode){
      sout << "ERROR: Cannot use --nostream with PGEN format.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(params->file_type == "bgen") {
      if(!file_exists (files->bgen_file)) {
        sout << "ERROR: " << files->bgen_file  << " doesn't exist.\n" << params->err_help ;
        exit(EXIT_FAILURE);
      }
      string bgifile = files->bgen_file + ".bgi";
      params->with_bgi = file_exists (bgifile) ;
    }
    if(vm.count("covarFile") & !file_exists (files->cov_file)) {
      sout << "ERROR: " << files->cov_file  << " doesn't exist.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(!file_exists (files->pheno_file)) {
      sout << "ERROR: " << files->pheno_file  << " doesn't exist.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if((params->file_type == "bed") & !file_exists (files->bed_prefix + ".bed")) {
      sout << "ERROR: " << files->bed_prefix << ".bed"  << " doesn't exist.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if((params->file_type == "bed") & !file_exists (files->bed_prefix + ".bim")) {
      sout << "ERROR: " << files->bed_prefix << ".bim"  << " doesn't exist.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if((params->file_type == "bed") & !file_exists (files->bed_prefix + ".fam")) {
      sout << "ERROR : " << files->bed_prefix << ".fam"  << " doesn't exist.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(params->snp_set && !file_exists (files->set_file)) {
      sout << "ERROR : " << files->set_file << " doesn't exist.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(params->build_mask && !file_exists (files->anno_file)) {
      sout << "ERROR : " << files->anno_file  << " doesn't exist.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(params->build_mask && vm.count("anno-labels") && !file_exists (files->anno_labs_file)) {
      sout << "ERROR : " << files->anno_labs_file  << " doesn't exist.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(params->build_mask && !file_exists (files->mask_file)) {
      sout << "ERROR : " << files->mask_file  << " doesn't exist.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }
    if(params->set_aaf && !file_exists (files->aaf_file)) {
      sout << "ERROR : " << files->aaf_file  << " doesn't exist.\n" << params->err_help ;
      exit(EXIT_FAILURE);
    }

  } catch (const cxxopts::OptionException& e) {
    print_header(cerr);
    cerr << "ERROR: " << e.what() << endl << params->err_help << endl;
    exit(EXIT_FAILURE);
  }

  return;
}


template <typename T> 
void start_log(T arguments, const string out_file, MeasureTime* mt, mstream& sout){

  string log_name = out_file + ".log";
  sout.coss.open(log_name.c_str(), ios::out | ios::trunc); 
  if (!sout.coss.is_open()) {
    cerr << "ERROR: Cannot write log file '" << log_name << "'\n" ;
    exit(EXIT_FAILURE);
  } 

  mt->init();
  sout << "Start time: " << ctime( &(mt->start_time_info) ) << endl; 
  print_header(sout.coss);
  print_header(cout);
  sout << "Log of output saved in file : " << log_name << endl<< endl;

  // print options
  sout << "Options in effect:" << endl ;
  for(size_t counter=1;counter<arguments.size();counter++){	  
    //if( trimmed_str[0] == '-') sout << "\\" << endl << "  ";
    sout << "  --" << arguments[counter-1].key() << " ";
    if( arguments[counter-1].value() == "true" ) {sout << "\\\n"; continue;}
    sout << arguments[counter-1].value() << " \\" << endl;
  }
  // last option (skip \ at the end)
  sout << "  --" << arguments.back().key() << " ";
  if( arguments.back().value() != "true" ) 
    sout << arguments.back().value();

  sout << "\n\n";

}

void set_ridge_params(int nparams, vector<double>& in_param, const string err_help, mstream& sout){

  if(nparams < 2){
    sout << "ERROR: Number of ridge parameters must be at least 2 (=" << nparams << ").\n" << err_help;
    exit(EXIT_FAILURE);
  } else {
    // endpoints are 0.01 and 0.99 
    double step = 1.0 / ( nparams - 1 );
    double val = step;
    in_param.resize( nparams);

    for( int index_p = 1; index_p < (nparams - 1); index_p++, val += step) in_param[index_p] = val;
    in_param[0] = 0.01;
    in_param[nparams-1] = 0.99;

  }
}

void print_usage_info(struct param* params, struct in_files* files, mstream& sout){

  double total_ram;
  string ram_unit;

  ///// Memory usage
  if(!params->test_mode){
    // Step 1
    // 4P + max( B + PRT, PRT) + #chrs [P:#traits;R=#ridge l0;T=#predictions from l0]
    int t_eff = ( params->write_l0_pred ? 1 : params->total_n_block );
    int p_eff = ( params->write_l0_pred ? 1 : params->n_pheno );
    int b_eff = params->total_n_block;

    total_ram = 4 * params->n_pheno + params->nChrom;
    total_ram += std::max( params->block_size + params->n_pheno * params->n_ridge_l0 * t_eff, p_eff * params->n_ridge_l0 * b_eff );
  } else {
    // Step 2
    // 3P + B
    total_ram = params->n_pheno * 3 + params->block_size; // y, mask, y_resid, g
    if(params->binary_mode) {
      total_ram += 2 * params->n_pheno + params->block_size; // y_raw, gamma_hat, g_resid
      if(params->use_SPA) total_ram += 0.5 * params->block_size; // non_zero_indices of g (4 bytes)
    }
    if((params->file_type == "bed") && params->fastMode) total_ram += params->block_size/4.0/sizeof(double); //for extracting snp_data_block
    if(params->mask_loo) total_ram += params->block_size; // loo masks
  }

  total_ram *= params->n_samples * sizeof(double);
  total_ram += params->nvs_stored * sizeof(struct snp);
  if( params->use_loocv ) total_ram += params->chunk_mb * 1e6; // max amount of memory used for LOO computations involved
  total_ram /= 1024.0 * 1024.0; 
  if( total_ram > 1000 ) {
    total_ram /= 1024.0; 
    ram_unit = "GB";
  } else ram_unit = "MB";

  int ram_int = (int) ceil( total_ram );
  sout << " * approximate memory usage : " << ram_int << ram_unit << endl;

  ///// Disk space usage
  if(!params->test_mode && !params->run_l1_only && params->write_l0_pred){
    if(files->loco_tmp_prefix.empty()) files->loco_tmp_prefix = files->out_file;
    sout << " * writing level 0 predictions to disk" << endl;
    sout << "   -temporary files will have prefix [" << files->loco_tmp_prefix << "_l0_Y]" << endl;
    // N*P*T*R
    int b_eff = params->total_n_block;
    total_ram = params->n_pheno * b_eff * params->n_ridge_l0;
    total_ram *= params->n_samples * sizeof(double);
    total_ram /= 1024.0 * 1024.0; 
    if( total_ram > 1000 ) {
      total_ram /= 1024.0; 
      ram_unit = "GB";
    } else ram_unit = "MB";
    int ram_int = (int) ceil( total_ram );
    sout << "   -approximate disk space needed : " << ram_int << ram_unit << endl;
  }
}

int chrStrToInt(const string chrom, const int nChrom) {

  // if label is chr1, chr2,...
  string s_chr = std::regex_replace(chrom, std::regex(R"(^chr)"), "");

  if (isdigit(s_chr[0])) {
    int chr = atoi(s_chr.c_str());
    if((chr >= 1) && (chr <= nChrom)) return chr;
  } else if ( (s_chr == "X") || (s_chr == "XY") || (s_chr == "PAR1") || (s_chr == "PAR2") ) return nChrom;

  return -1;
}

double convertDouble(const string& val, struct param* params, mstream& sout){

  if(val == params->missing_pheno_str)
    return params->missing_value_double;

  double dval;
  if(sscanf(val.c_str(), "%lf", &dval) != 1){
    sout << "ERROR: Could not convert value to double: " << val << endl;
    exit(EXIT_FAILURE);
  }
  return dval;
}



