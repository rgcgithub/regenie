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

#include "cxxopts/include/cxxopts.hpp"
#include "Regenie.hpp"
#include "Geno.hpp"
#include "Step1_Models.hpp"
#include "Step2_Models.hpp"
#include "Files.hpp"
#include "Pheno.hpp"
#include "Data.hpp"


using namespace std;
using namespace Eigen;
using namespace boost;


mstream::mstream(){ }
mstream::~mstream(){ }
MeasureTime::MeasureTime(){ }
MeasureTime::~MeasureTime(){ }


// This example program reads data from a bgen file specified as the first argument
// and outputs it as a VCF file.
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
  o << "\n";
}


void read_params_and_check(int argc, char *argv[], struct param* params, struct in_files* files, struct filter* filters, MeasureTime* mt, mstream& sout) {

  string webinfo = "For more information, visit the website: https://rgcgithub.github.io/regenie/";

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
    ("with-bgi", "read variant information in BGEN file from a bgi index file")
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
    ("bt", "analyze phenotypes as binary")
    ("1,cc12", "use control=1,case=2,missing=NA encoding for binary traits")
    ("b,bsize", "size of genotype blocks", cxxopts::value<int>(params->block_size),"INT")
    ("cv", "number of cross validation (CV) folds", cxxopts::value<int>(params->cv_folds),"INT(=5)")
    ("loocv", "use leave-one out cross validation (LOOCV)")
    ("l0", "number of ridge parameters to use when fitting models within blocks [evenly spaced in (0,1)]", cxxopts::value<int>(params->n_ridge_l0),"INT(=5)")
    ("l1", "number of ridge parameters to use when fitting model across blocks [evenly spaced in (0,1)]", cxxopts::value<int>(params->n_ridge_l1),"INT(=5)")
    ("lowmem", "reduce memory usage by writing level 0 predictions to temporary files")
    ("lowmem-prefix", "prefix where to write the temporary files in step 1 (default is to use prefix from --out)", cxxopts::value<std::string>(files->loco_tmp_prefix),"PREFIX")
    ("strict", "remove all samples with missingness at any of the traits")
    ("o,out", "prefix for output files", cxxopts::value<std::string>(files->out_file),"PREFIX")
    ("threads", "number of threads", cxxopts::value<int>(params->threads),"INT")
    ("pred", "file containing the list of predictions files from step 1", cxxopts::value<std::string>(files->blup_file),"FILE")
    ("ignore-pred", "skip reading predictions from step 1 (equivalent to linear/logistic regression with only covariates)")
    ("force-impute", "keep and impute missing observations when in step 2 (default is to drop missing for each trait)")
    ("write-samples", "write IDs of samples included for each trait (only in step 2)")
    ("minMAC", "minimum minor allele count (MAC) for tested variants", cxxopts::value<int>(params->min_MAC),"INT(=5)")
    ("split", "split asssociation results into separate files for each trait")
    ("firth", "use Firth correction for p-values less than threshold")
    ("approx", "use approximation to Firth correction for computational speedup")
    ("spa", "use Saddlecxxoptsint approximation (SPA) for p-values less than threshold")
    ("pThresh", "P-value threshold below which to apply Firth/SPA correction", cxxopts::value<double>(params->alpha_pvalue),"FLOAT(=0.05)")
    ("chr", "specify chromosome to test in step 2 (use for each chromosome)", cxxopts::value< std::vector<std::string> >(),"STRING")
    ("chrList", "Comma separated list of chromosomes to test in step 2", cxxopts::value<std::string>(),"STRING,..,STRING")
    ("test", "'dominant' or 'recessive' (default is additive test)", cxxopts::value<std::string>(),"STRING")
    ("gz", "compress output files (gzip format)")
    ;


  // extended options
  AllOptions.add_options("Additional")
    ("v,verbose", "verbose screen output")
    ("setl0", "comma separated list of ridge parameters to use when fitting models within blocks", cxxopts::value<std::string>(), "FLOAT,..,FLOAT")
    ("setl1", "comma separated list of ridge parameters to use when fitting model across blocks", cxxopts::value<std::string>(), "FLOAT,..,FLOAT")
    ("nauto", "number of autosomal chromosomes", cxxopts::value<int>(),"INT")
    ("nb", "number of blocks to use", cxxopts::value<int>(params->n_block),"INT")
    ("niter", "maximum number of iterations for logistic regression", cxxopts::value<int>(params->niter_max),"INT(=30)")
    ("maxstep-null", "maximum step size in null Firth logistic regression", cxxopts::value<int>(params->maxstep_null),"INT(=25)")
    ("maxiter-null", "maximum number of iterations in null Firth logistic regression", cxxopts::value<int>(params->niter_max_firth_null),"INT(=1000)")
    ;

  // extra options
  AllOptions.add_options("Extra")
    ("print", "print estimated effect sizes from level 0 and level 1 models")
    ("nostream", "print estimated effect sizes from level 0 and level 1 models")
    ("htp", "output association files in step 2 in HTPv4 format", cxxopts::value<std::string>(params->cohort_name),"STRING")
    ("within", "use within-sample predictions as input when fitting model across blocks in step 1")
    ("early-exit", "Exit program after fitting level 0 models (avoid deleting temporary prediction files from level 0)")
    ;


  try
  {
    //AllOptions.parse_positional({"htp"});
    auto vm = AllOptions.parse(argc, argv);
    auto arguments = vm.arguments();

    // help menu
    if (vm.count("help")){
      print_header(std::cout);
      std::cout << AllOptions.help({"", "Main"}) << '\n' << webinfo << "\n\n";
      exit(-1);
    } else if (vm.count("helpFull")) {
      print_header(std::cout);
      std::cout << AllOptions.help({"", "Main", "Additional"}) << '\n' << webinfo << "\n\n";
      exit(-1);
    } 
    
    if (!vm.count("out")){
      print_header(std::cout);
      std::cout << "ERROR :You must provide an output prefix using '--out'" << '\n' << webinfo << "\n\n";
      exit(-1);
    }


    // Print output to file and to stdout
    // print command line arguments
    start_log(arguments, files->out_file, mt, sout);
    vector< string > tmp_str_vec;

    if( (vm.count("bgen") + vm.count("bed")  + vm.count("pgen"))  != 1 ){
      sout << "ERROR :You must use either --bed,--bgen or --pgen.\n" << params->err_help ;
      exit(-1);
    }

    if( vm.count("bgen") ) params->file_type = "bgen";
    if( vm.count("bed") ) params->file_type = "bed";
    if( vm.count("pgen") ) params->file_type = "pgen";
    if( vm.count("sample") ) params->bgenSample = true;
    if( vm.count("with-bgi") ) params->with_bgi = true;
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
    if( vm.count("ignore-pred") ) params->skip_blups = true;
    if( vm.count("force-impute") ) params->rm_missing_qt = false;
    if( vm.count("split") ) params->split_by_pheno = true;
    if( vm.count("approx") ) params->firth_approx = true;
    if( vm.count("nauto") ) params->nChrom = vm["nauto"].as<int>() + 1;
    if( vm.count("maxstep-null") | vm.count("maxiter-null") ) params->fix_maxstep_null = true;
    if( vm.count("firth") ) params->firth = true;
    if( vm.count("spa") ) params->use_SPA = true;
    if( vm.count("minMAC") ) params->setMinMAC = true;
    if( vm.count("htp") ) params->htp_out = params->split_by_pheno = true;
    if( vm.count("v") ) params->verbose = true;
    if( vm.count("print") ) params->print_block_betas = true;
    if( vm.count("nostream") ) params->streamBGEN = false;
    if( vm.count("within") ) params->within_sample_l0 = true;
    if( vm.count("write-samples") ) params->write_samples = true;
    if( vm.count("early-exit") ) params->early_exit = true;
    if( vm.count("gz") ) {
# if defined(HAS_BOOST_IOSTREAM)
      // only works when compiled with boost IO library
      params->gzOut = true;
# else
      sout << "WARNING : REGENIE was not compiled with Boost Iostream library so ignoring option '--gz'.\n";
#endif
    }


    if( vm.count("lowmem") ) {
      params->write_l0_pred = true;
    }
    if( vm.count("phenoColList") ) {
      params->select_phenos = true;
      boost::algorithm::split(tmp_str_vec, vm["phenoColList"].as<string>(), is_any_of(","));
      filters->pheno_colKeep_names.insert( 
          filters->pheno_colKeep_names.end(),
          std::begin( tmp_str_vec ), 
          std::end( tmp_str_vec)         );
    }
    if( vm.count("covarColList") ) {
      params->select_covs = true;
      boost::algorithm::split(tmp_str_vec, vm["covarColList"].as<string>(), is_any_of(","));
      filters->cov_colKeep_names.insert( 
          filters->cov_colKeep_names.end(),
          std::begin( tmp_str_vec ), 
          std::end( tmp_str_vec)         );
    }
    if( vm.count("chrList") ) {
      params->select_chrs = true;
      boost::algorithm::split(tmp_str_vec, vm["chrList"].as<string>(), is_any_of(","));
      for( size_t ichr = 0; ichr < tmp_str_vec.size(); ichr++)
        filters->chrKeep_test.push_back( chrStrToInt(tmp_str_vec[ichr], params->nChrom) );
    }
    if( vm.count("chr") ) {
      params->select_chrs = true;
      tmp_str_vec = vm["chr"].as<std::vector<string>>();
      for( size_t ichr = 0; ichr < tmp_str_vec.size(); ichr++)
        filters->chrKeep_test.push_back( chrStrToInt(tmp_str_vec[ichr], params->nChrom) );
    }
    if( vm.count("test") ) {
      if( vm["test"].as<string>() == "dominant") params->test_type = 1; 
      else if( vm["test"].as<string>() == "recessive") params->test_type = 2; 
      else {
        sout << "ERROR : Unrecognized argument for option --test, must be either 'dominant' or 'recessive'.\n" << params->err_help;
        exit(-1);
      }
    }

    if ( params->run_mode == 1 ) params->test_mode = false;
    else if (params->run_mode == 2 ) params->test_mode = true;
    else {
      sout << "ERROR : Specify which mode regenie should be running using option --step.\n" << params->err_help;
      exit(-1);
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
        boost::algorithm::split(tmp_str_vec, vm["setl0"].as<string>(), is_any_of(","));
        for( size_t val = 0; val < tmp_str_vec.size(); val++)
          params->lambda.push_back(convertDouble( tmp_str_vec[val], params, sout));
        std::sort(params->lambda.begin(), params->lambda.end());
        params->lambda.erase( unique( params->lambda.begin(), params->lambda.end() ), params->lambda.end() );
        params->n_ridge_l0 = params->lambda.size();
        // parameters must be less in (0, 1)
        if( std::count_if(params->lambda.begin(), params->lambda.end(), std::bind2nd(std::greater<double>(), 0)) != params->n_ridge_l0 || std::count_if(params->lambda.begin(), params->lambda.end(), std::bind2nd(std::less<double>(), 1)) != params->n_ridge_l0 ){
          sout << "ERROR : You must specify values for --l0 in (0,1).\n" << params->err_help;
          exit(-1);
        } 
      } else set_ridge_params(params->n_ridge_l0, params->lambda, params->err_help, sout);

      // user specified ridge parameters to use at l1
      if( vm.count("setl1") ) {
        params->user_ridge_params_l1 = true;
        boost::algorithm::split(tmp_str_vec, vm["setl1"].as<string>(), is_any_of(","));
        for( size_t val = 0; val < tmp_str_vec.size(); val++)
          params->tau.push_back(convertDouble( tmp_str_vec[val], params, sout));
        std::sort(params->tau.begin(), params->tau.end());
        params->tau.erase( unique( params->tau.begin(), params->tau.end() ), params->tau.end() );
        params->n_ridge_l1 = params->tau.size();
        if( std::count_if(params->tau.begin(), params->tau.end(), std::bind2nd(std::greater<double>(), 0)) != params->n_ridge_l1 || std::count_if(params->tau.begin(), params->tau.end(), std::bind2nd(std::less<double>(), 1)) != params->n_ridge_l1 ){
          sout << "ERROR : You must specify values for --l1 in (0,1).\n" << params->err_help;
          exit(-1);
        }
      } else set_ridge_params(params->n_ridge_l1, params->tau, params->err_help, sout);

      // firth only done in test mode
      if(params->firth) params->firth = false;
      if(params->use_SPA) params->use_SPA = false;
      params->streamBGEN = false;
      params->test_type = 0;

    } else if(params->firth && !params->binary_mode) {
      // firth correction is only applied to binary traits
      params->firth = false;
    } else if(params->use_SPA && !params->binary_mode) {
      // SPA is only applied to binary traits
      params->use_SPA = false;
    }


    if(params->test_mode && params->use_loocv) params->use_loocv = false;
    if(params->test_mode && params->rm_snps) params->rm_snps = false;
    if(params->test_mode && params->keep_snps) params->keep_snps = false;

    if(!params->test_mode && params->setMinMAC){
      sout << "WARNING : Option --minMAC only works in step 2 of REGENIE.\n";
    }
    if(params->test_mode && params->min_MAC < 1){
      sout << "ERROR : minimum MAC must be at least 1.\n" << params->err_help;
      exit(-1);
    }
    if( params->rm_missing_qt && (params->strict_mode || params->binary_mode || !params->test_mode) ) params->rm_missing_qt = false;

    if( !vm.count("bsize") ) {
      sout << "ERROR : must specify the block size using '--bsize'.\n" << params->err_help;
      exit(-1);
    }

    // determine number of threads if not specified
    if(params->threads < 1){
      params->threads = std::thread::hardware_concurrency(); //may return 0 when not able to detect
      if(params->threads < 1) params->threads = 1;
    }

    // set Firth as default if both Firth and SPA are specified
    if(params->use_SPA && params->firth) params->use_SPA = false;

    // check firth fallback pvalue threshold
    if(params->firth && ((params->alpha_pvalue < params->nl_dbl_dmin) || (params->alpha_pvalue > 1 - params->numtol)) ){
      sout << "ERROR :Firth fallback p-value threshold must be in (0,1).\n" << params->err_help ;
      exit(-1);
    }
    // check SPA fallback pvalue threshold
    if(params->use_SPA && ((params->alpha_pvalue < params->nl_dbl_dmin) || (params->alpha_pvalue > 1 - params->numtol)) ){
      sout << "ERROR :SPA fallback p-value threshold must be in (0,1).\n" << params->err_help ;
      exit(-1);
    }
    if(params->firth_approx && !params->firth) params->firth_approx = false;

    // check arguments for logistic regression 
    if(params->binary_mode && (params->niter_max < 1)){
      sout << "ERROR :Invalid argument for --niter (must be positive integer).\n" << params->err_help ;
      exit(-1);
    }
    if(params->firth && (params->maxstep_null < 1)){
      sout << "ERROR :Invalid argument for --maxstep-null (must be a positive integer).\n" << params->err_help ;
      exit(-1);
    }
    if(params->firth && (params->niter_max_firth_null < 1)){
      sout << "ERROR :Invalid argument for --maxiter-null (must be a positive integer).\n" << params->err_help ;
      exit(-1);
    }
    if(params->nChrom < 2){
      sout << "ERROR :Invalid argument for --nauto (must be > 1).\n" << params->err_help ;
      exit(-1);
    }
    if(params->rm_indivs && params->keep_indivs ){
      sout << "ERROR :Cannot use both --keep and --remove.\n" << params->err_help ;
      exit(-1);
    }
    if(params->rm_snps && params->keep_snps ){
      sout << "ERROR :Cannot use both --extract and --exclude.\n" << params->err_help ;
      exit(-1);
    }

    if( params->test_mode && params->select_chrs && std::count( filters->chrKeep_test.begin(), filters->chrKeep_test.end(), -1) ){
      sout << "ERROR :Invalid chromosome specified by --chr/--chrList.\n" << params->err_help ;
      exit(-1);
    }

    if(params->test_mode && !params->skip_blups && !vm.count("pred")) {
      sout << "ERROR :You must specify --pred if using --step 2 (otherwise use --ignore-pred).\n" << params->err_help ;
      exit(-1);
    }

    if(params->test_mode && (params->file_type == "pgen") && !params->streamBGEN){
      sout << "ERROR :Cannot use --nostream with PGEN format.\n" << params->err_help ;
      exit(-1);
    }

    if((params->file_type == "bgen") & (!file_exists (files->bgen_file))) {
      sout << "ERROR : " << files->bgen_file  << " doesn't exist.\n" << params->err_help ;
      exit(-1);
    }
    if(vm.count("covarFile") & !file_exists (files->cov_file)) {
      sout << "ERROR : " << files->cov_file  << " doesn't exist.\n" << params->err_help ;
      exit(-1);
    }
    if(!file_exists (files->pheno_file)) {
      sout << "ERROR : " << files->pheno_file  << " doesn't exist.\n" << params->err_help ;
      exit(-1);
    }
    if((params->file_type == "bed") & !file_exists (files->bed_prefix + ".bed")) {
      sout << "ERROR : " << files->bed_prefix << ".bed"  << " doesn't exist.\n" << params->err_help ;
      exit(-1);
    }
    if((params->file_type == "bed") & !file_exists (files->bed_prefix + ".bim")) {
      sout << "ERROR : " << files->bed_prefix << ".bim"  << " doesn't exist.\n" << params->err_help ;
      exit(-1);
    }
    if((params->file_type == "bed") & !file_exists (files->bed_prefix + ".fam")) {
      sout << "ERROR : " << files->bed_prefix << ".fam"  << " doesn't exist.\n" << params->err_help ;
      exit(-1);
    }

  } catch (const cxxopts::OptionException& e) {
    print_header(cerr);
    cerr << "ERROR: " << e.what() << endl << params->err_help << endl;
    exit(-1);
  }

  return;
}


template <typename T> 
void start_log(T arguments, const string out_file, MeasureTime* mt, mstream& sout){

  string log_name = out_file + ".log";
  sout.coss.open(log_name.c_str(), ios::out | ios::trunc); 
  if (!sout.coss.is_open()) {
    cerr << "ERROR : Cannot write log file '" << log_name << "'\n" ;
    exit(-1);
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
    sout << "ERROR : Number of ridge parameters must be at least 2 (=" << nparams << ").\n" << err_help;
    exit(-1);
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
    total_ram = 4 * params->n_pheno + params->nChrom;
    int t_eff = ( params->write_l0_pred ? 1 : params->total_n_block );
    int p_eff = ( params->write_l0_pred ? 1 : params->n_pheno );
    total_ram += std::max( params->block_size + params->n_pheno * params->n_ridge_l0 * t_eff, p_eff * params->n_ridge_l0 * params->total_n_block );
  } else {
    // Step 2
    // 3P + B
    total_ram = params->n_pheno * 3 + params->block_size; // y, mask, y_resid, g
    if(params->binary_mode) {
      total_ram += 2 * params->n_pheno + params->block_size; // y_raw, gamma_hat, g_resid
      if(params->use_SPA) total_ram += 0.5 * params->block_size; // non_zero_indices of g (4 bytes)
    }
    if((params->file_type == "bed") && params->streamBGEN) total_ram += params->block_size/32.0; //for extracting snp_data_block
  }

  total_ram *= params->n_samples * sizeof(double);
  if( params->use_loocv ) total_ram += params->chunk_mb * 1e6; // max amount of memory used for LOO computations involved
  total_ram /= 1024.0 * 1024.0; 
  if( total_ram > 1000 ) {
    total_ram /= 1024.0; 
    ram_unit = "GB";
  } else ram_unit = "MB";

  int ram_int = (int) ceil( total_ram );
  sout << " * approximate memory usage : " << ram_int << ram_unit << endl;

  ///// Disk space usage
  if(!params->test_mode && params->write_l0_pred){
    if(files->loco_tmp_prefix.empty()) files->loco_tmp_prefix = files->out_file;
    sout << " * writing level 0 predictions to disk" << endl;
    sout << "   -temporary files will have prefix [" << files->loco_tmp_prefix << "_l0_Y]" << endl;
    // N*P*T*R
    total_ram = params->n_pheno * params->total_n_block * params->n_ridge_l0;
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


