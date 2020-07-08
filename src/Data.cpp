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

#include "Data.hpp"

inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}



using namespace Eigen;

using boost::math::normal;
using boost::math::chi_squared;

mstream::mstream(){ }
mstream::~mstream(){ }
MeasureTime::MeasureTime(){ }
MeasureTime::~MeasureTime(){ }

Data::Data() { // @suppress("Class members should be properly initialized")
 }

Data::~Data() {
	// TODO Auto-generated destructor stub
}


void Data::run() {

  if(test_mode) {
    if(streamBGEN) check_bgen();
    if(streamBGEN) test_snps_fast();
    else test_snps();

  } else {
    sout << "Fitting null model" << endl;

    // set number of threads
    setNbThreads(threads);
    // set up file for reading
    file_read_initialization();
    // read phenotype and covariate files
    read_pheno_and_cov();
    // set number of blocks and block size and ridge parameters
    set_blocks();
    // some initializations
    setmem();
    // level 0
    level_0_calculations();
    // level 1 ridge
    if(!binary_mode) // QT
      if(use_loocv) ridge_level_1_loocv();
      else ridge_level_1();
    else // BT
      if(use_loocv) ridge_logistic_level_1_loocv();
      else ridge_logistic_level_1();
    // output results
    output();
  }
}

void Data::print_help( bool help_full ){

  print_header(cout);
  cout << "Options:\n";

  cout << left << std::setw(35) << " --help" << "print list of available options\n";
  cout << left << std::setw(35) << " --helpFull" << "print list of all available options\n";

  cout << left << std::setw(35) << " --step INT"<< "specify if fitting null model (=1) or association testing (=2)\n";
  cout << left << std::setw(35) << " --bed PREFIX" << "prefix to PLINK .bed/.bim/.fam files\n";
  cout << left << std::setw(35) << " --bgen FILE" << "BGEN file\n";
  cout << left << std::setw(35) << " --sample FILE" << "sample file corresponding to BGEN file\n";
  cout << left << std::setw(35) << " --keep FILE" << "file listing samples to retain in the analysis\n" << 
    std::setw(35) << " " << "(no header; starts with FID IID)\n";
  cout << left << std::setw(35) << " --remove FILE" << "file listing samples to remove from the analysis\n" << 
    std::setw(35) << " " << "(no header; starts with FID IID)\n";
  cout << left << std::setw(35) << " --extract FILE" << "file with IDs of variants to retain in the analysis\n";
  cout << left << std::setw(35) << " --exclude FILE" << "file with IDs of variants to remove from the analysis\n";
  cout << left << std::setw(35) << " --p FILE" << "phenotype file (header required starting with FID IID)\n";
  cout << left << std::setw(35) << " --phenoCol STRING"<< "phenotype name in header (use for each phenotype to keep)\n";
  cout << left << std::setw(35) << " --c FILE" << "covariate file (header required starting with FID IID)\n";
  cout << left << std::setw(35) << " --covarCol STRING"<< "covariate name in header (use for each covariate to keep)\n";
  cout << left << std::setw(35) << " --b INT"<< "size of genotype blocks\n";
  cout << left << std::setw(35) << " --cv INT (=5)"<< "number of cross validation (CV) folds\n";
  cout << left << std::setw(35) << " --l0 INT (=5)"<< "number of ridge parameters to use when fitting models\n" << 
    std::setw(35) << " " << "within blocks (evenly spaced in (0,1))\n";
  cout << left << std::setw(35) << " --l1 INT (=5)"<< "number of ridge parameters to use when fitting model\n" <<
    std::setw(35) << " " << "across blocks (evenly spaced in (0,1))\n";
  cout << left << std::setw(35) << " --loocv"<< "use leave-one out cross validation (LOOCV)\n";
  cout << left << std::setw(35) << " --bt"<< "analyze phenotypes as binary (default is quantitative)\n";
  cout << left << std::setw(35) << " --lowmem PREFIX "<< "reduce memory usage by writing level 0 predictions\n" << 
   std::setw(35) << " " << "to temporary files\n";
  cout << left << std::setw(35) << " --threads INT (=ALL)"<< "number of threads\n";
  cout << left << std::setw(35) << " --strict"<< "remove all samples with missingness at any of the traits\n";
  cout << left << std::setw(35) << " --o PREFIX" << "prefix for output files\n";
  cout << left << std::setw(35) << " --pred FILE" << "file containing the list of files with predictions from step 1\n";
  cout << left << std::setw(35) << " --ignore-pred"<< "skip reading predictions from step 1 (equivalent to\n" <<
    std::setw(35) << " " << "linear/logistic regression with only covariates)\n";;
  cout << left << std::setw(35) << " --force-impute" << "keep and impute missing observations when in step 2 (default is\n" <<
    std::setw(35) << " " << "to drop missing for each trait)\n";
  cout << left << std::setw(35) << " --minMAC INT (=5)"<< "minimum minor allele count (MAC) for tested variants\n";
  cout << left << std::setw(35) << " --split" << "split asssociation results into separate files for each trait\n";
  cout << left << std::setw(35) << " --firth FLOAT (=0.05)" << "use Firth correction for p-values less than threshold\n";
  cout << left << std::setw(35) << " --approx" << "use approximation to Firth correction for computational speedup\n";
  cout << left << std::setw(35) << " --spa FLOAT (=0.05)" << "use Saddlepoint approximation (SPA) for p-values less\n" <<
    std::setw(35) << " " << "than threshold\n";
  cout << left << std::setw(35) << " --chr INT" << "specify chromosome to test in step 2 (use for each chromosome)\n";
  cout << left << std::setw(35) << " --test STRING" << "specify to use dominant or recessive test\n";

  if(help_full){
  cout << left << std::setw(35) << " --v" << "verbose screen output\n";
  cout << left << std::setw(35) << " --nb INT"<< "number of blocks to use\n";
  cout << left << std::setw(35) << " --setl0 FLOAT ... FLOAT"<< "specify ridge parameters to use when fitting models within blocks\n";
  cout << left << std::setw(35) << " --setl1 FLOAT ... FLOAT"<< "specify ridge parameters to use when fitting model across blocks\n";
  cout << left << std::setw(35) << " --nauto INT (=22)"<< "number of autosomal chromosomes\n";
  cout << left << std::setw(35) << " --niter INT (=30)"<< "maximum number of iterations for logistic regression\n";
  cout << left << std::setw(35) << " --maxstep-null INT (=25)"<< "maximum step size in null Firth logistic regression\n";
  cout << left << std::setw(35) << " --maxiter-null INT (=25)"<< "maximum number of iterations in null Firth logistic regression\n";
  cout << left << std::setw(35) << " --within" << "use within-sample predictions as input when fitting model\n" <<
    std::setw(35) << " " << "across blocks.\n";
  }

  cout << "\nFor more information, visit the website: https://rgcgithub.github.io/regenie/\n";
  exit(-1);
}

void Data::read_params_and_check(int argc, char *argv[]) {

  int maxargs = argc - 1;

  for(size_t counter=1;counter<argc;counter++){	  
    if(string(argv[counter]) == "--bt") binary_mode = true;
    if(string(argv[counter]) == "--1") control_code = 1;
    if(string(argv[counter]) == "--within") within_sample_l0 = true;
    if(string(argv[counter]) == "--loocv") use_loocv = true;
    if(string(argv[counter]) == "--split") split_by_pheno = true;
    if(string(argv[counter]) == "--strict") strict_mode = true;
    if(string(argv[counter]) == "--force-impute") rm_missing_qt = false;
    if(string(argv[counter]) == "--ignore-pred") skip_blups = true;
    if(string(argv[counter]) == "--print") print_block_betas = true;
    if(string(argv[counter]) == "--approx") firth_approx = true;
    if(string(argv[counter]) == "--v") verbose = true;
    if(string(argv[counter]) == "--nostream") streamBGEN = false;
    if(string(argv[counter]) == "--helpFull") print_help( true );
    if(string(argv[counter]) == "--help") print_help( false );

    if(string(argv[counter]) == "--bgen") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) {
        bgen_file = string(argv[counter+1]);
        file_type = "bgen";
        n_genofiles++;
      }
    }
    if(string(argv[counter]) == "--c") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) cov_file = string(argv[counter+1]);
    }
    if(string(argv[counter]) == "--p") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) pheno_file = string(argv[counter+1]);
    }
    if(string(argv[counter]) == "--pred") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) blup_file = string(argv[counter+1]);
    }
    if(string(argv[counter]) == "--o") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) out_file = string(argv[counter+1]);
    }
    if(string(argv[counter]) == "--lowmem") {
      write_l0_pred = true;
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) loco_tmp_prefix = string(argv[counter+1]);
    }
    if(string(argv[counter]) == "--b") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) block_size = atoi(argv[counter+1]);
    }
    if(string(argv[counter]) == "--nb") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) n_block = atoi(argv[counter+1]);
    }
    if(string(argv[counter]) == "--cv") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) cv_folds = atoi(argv[counter+1]);
    }
    if(string(argv[counter]) == "--l0") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) n_ridge_l0 = atoi(argv[counter+1]);
    }
    if(string(argv[counter]) == "--l1") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) n_ridge_l1 = atoi(argv[counter+1]);
    }
    if(string(argv[counter]) == "--step") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) run_mode = atoi(argv[counter+1]);
    }
    if(string(argv[counter]) == "--nauto") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) nChrom = atoi(argv[counter+1]) + 1;
    }
    if(string(argv[counter]) == "--niter") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) niter_max = atoi(argv[counter+1]);
    }
    if(string(argv[counter]) == "--maxstep-null") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) {
        maxstep_null = atoi(argv[counter+1]);
        fix_maxstep_null = true;
      }
    }
    if(string(argv[counter]) == "--maxiter-null") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) {
        niter_max_firth_null = atoi(argv[counter+1]);
        fix_maxstep_null = true;
      }
    }
    if(string(argv[counter]) == "--minMAC") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) min_MAC = atoi(argv[counter+1]);
    }
    if(string(argv[counter]) == "--bed") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) {
        bedfile = string(argv[counter+1]) + ".bed";
        bimfile = string(argv[counter+1]) + ".bim";
        famfile = string(argv[counter+1]) + ".fam";
        file_type = "bed";
        n_genofiles++;
      }
    }
    if(string(argv[counter]) == "--sample") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) {
        sample_file = string(argv[counter+1]);
        bgenSample = true;
      }
    }
    if(string(argv[counter]) == "--remove") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) {
        rm_indivs = true;
        file_ind_exclude = string(argv[counter+1]);
      }
    }
    if(string(argv[counter]) == "--keep") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) {
        keep_indivs = true;
        file_ind_include = string(argv[counter+1]);
      }
    }
    if(string(argv[counter]) == "--extract") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) {
        keep_snps = true;
        file_snps_include = string(argv[counter+1]);
      }
    }
    if(string(argv[counter]) == "--exclude") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) {
        rm_snps = true;
        file_snps_exclude = string(argv[counter+1]);
      }
    }
    if(string(argv[counter]) == "--phenoCol") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) {
        select_phenos = true;
        pheno_colKeep_names.push_back( string(argv[counter+1]) );
      }
    }
    if(string(argv[counter]) == "--covarCol") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) {
        select_covs = true;
        cov_colKeep_names.push_back( string(argv[counter+1]) );
      }
    }
    if(string(argv[counter]) == "--chr") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) {
        select_chrs = true;
        chrKeep_test.push_back( chrStrToInt(string(argv[counter+1])) );
      }
    }
    if(string(argv[counter]) == "--test") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) {
        if( string(argv[counter+1]) == "dominant") test_type = 1; 
        else if( string(argv[counter+1]) == "recessive") test_type = 2; 
      }
    }
    if(string(argv[counter]) == "--setl0") {
      n_ridge_l0 = 0;
      user_ridge_params_l0 = true;
      for(size_t i = 1; (counter+i) < argc; i++) {
        if( string(argv[counter+i])[0] != '-' ) n_ridge_l0++;
        else break;
      }
      lambda.resize(n_ridge_l0);
      for(size_t i = 0; i < n_ridge_l0; i++) lambda[i] = atof(argv[counter+1+i]);
    }
    if(string(argv[counter]) == "--setl1") {
      n_ridge_l1 = 0;
      user_ridge_params_l1 = true;
      for(size_t i = 1; (counter+i) < argc; i++) {
        if( string(argv[counter+i])[0] != '-' ) n_ridge_l1++;
        else break;
      }
      tau.resize(n_ridge_l1);
      for(size_t i = 0; i < n_ridge_l1; i++) tau[i] = atof(argv[counter+1+i]);
    }
    if(string(argv[counter]) == "--firth") {
      firth = true;
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) alpha_pvalue = atof(argv[counter+1]);;
    }
    if(string(argv[counter]) == "--spa") {
      use_SPA = true;
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) alpha_pvalue = atof(argv[counter+1]);
    }
    if(string(argv[counter]) == "--threads") {
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) threads = atoi(argv[counter+1]);
    }
    if(string(argv[counter]) == "--htp") {
      htp_out = split_by_pheno = true;
      if( (counter < maxargs) && (string(argv[counter+1])[0] != '-') ) cohort_name = string(argv[counter+1]);
    }

  }

  if(out_file == "NULL") {
    print_header(cerr);
    cerr << "ERROR :You must specify an output file with --o.\n" << err_help ;
    exit(-1);
  }
  // Print output to file and to stdout
  // print command line arguments
  start_log(argc, argv);


  if ( run_mode == 1 ) test_mode = false;
  else if (run_mode == 2 ) test_mode = true;
  else {
    sout << "ERROR : Specify which mode regenie should be running using option --step.\n" << err_help;
    exit(-1);
  }

  if(!test_mode) {

    // loocv only used with out-of-sample predictions
    if(use_loocv && within_sample_l0) {
      sout << "WARNING : Option --loocv cannot be used with option --within.\n" ;
      use_loocv = false;
    }

    // writing of level 0 predictions only available when using out-of-sample predictions
    if(write_l0_pred && within_sample_l0){
      sout << "WARNING : Option --lowmem cannot be used with option --within.\n" ;
      write_l0_pred = false;
    }

    // user specified ridge parameters to use
    if( user_ridge_params_l0 ){
      if(n_ridge_l0 < 1){
        sout << "ERROR : Number of ridge parameters must be at least 1 in option --setl0.\n" << err_help;
        exit(-1);
      }
      // parameters must be less in (0, 1)
      if( std::count_if(lambda.begin(), lambda.end(), std::bind2nd(std::greater<double>(), 0)) != n_ridge_l0 || std::count_if(lambda.begin(), lambda.end(), std::bind2nd(std::less<double>(), 1)) != n_ridge_l0 ){
        sout << "ERROR : You must specify values for --l0 in (0,1).\n" << err_help;
        exit(-1);
      } 
    } else set_ridge_params(n_ridge_l0, lambda);

    // user specified ridge parameters to use
    if( user_ridge_params_l1 ){
      if(n_ridge_l1 < 1){
        sout << "ERROR : Number of ridge parameters must be at least 1 in option --setl1.\n" << err_help;
        exit(-1);
      }
      if( std::count_if(tau.begin(), tau.end(), std::bind2nd(std::greater<double>(), 0)) != n_ridge_l1 || std::count_if(tau.begin(), tau.end(), std::bind2nd(std::less<double>(), 1)) != n_ridge_l1 ){
        sout << "ERROR : You must specify values for --l1 in (0,1).\n" << err_help;
        exit(-1);
      }
    } else set_ridge_params(n_ridge_l1, tau);

    // firth only done in test mode
    if(firth) firth = false;
    if(use_SPA) use_SPA = false;
    streamBGEN = false;
    test_type = 0;

  } else if(firth && !binary_mode) {
    // firth correction is only applied to binary traits
    firth = false;
  } else if(use_SPA && !binary_mode) {
    // SPA is only applied to binary traits
    use_SPA = false;
  }

  if(test_mode && rm_snps) rm_snps = false;
  if(test_mode && keep_snps) keep_snps = false;

  if(test_mode && min_MAC < 1){
    sout << "ERROR : minimum MAC must be at least 1.\n" << err_help;
    exit(-1);
  }
  if( rm_missing_qt && (strict_mode || binary_mode || !test_mode) ) rm_missing_qt = false;

  // determine number of threads if not specified
  if(threads < 1){
    threads = std::thread::hardware_concurrency(); //may return 0 when not able to detect
    if(threads < 1) threads = 1;
  }

  // set Firth as default if both Firth and SPA are specified
  if(use_SPA && firth) use_SPA = false;

  // check SPA fallback pvalue threshold
  if(use_SPA && ((alpha_pvalue < nl_dbl_dmin) || (alpha_pvalue > 1 - numtol)) ){
    sout << "ERROR :SPA fallback p-value threshold must be in (0,1).\n" << err_help ;
    exit(-1);
  }
  // check firth fallback pvalue threshold
  if(firth && ((alpha_pvalue < nl_dbl_dmin) || (alpha_pvalue > 1 - numtol)) ){
    sout << "ERROR :Firth fallback p-value threshold must be in (0,1).\n" << err_help ;
    exit(-1);
  }
  if(firth_approx && !firth) firth_approx = false;

  // check arguments for logistic regression 
  if(binary_mode && (niter_max < 1)){
    sout << "ERROR :Invalid argument for --niter (must be positive integer).\n" << err_help ;
    exit(-1);
  }
  if(firth && (maxstep_null < 1)){
    sout << "ERROR :Invalid argument for --maxstep-null (must be a positive integer).\n" << err_help ;
    exit(-1);
  }
  if(firth && (niter_max_firth_null < 1)){
    sout << "ERROR :Invalid argument for --maxiter-null (must be a positive integer).\n" << err_help ;
    exit(-1);
  }
  if(nChrom < 2){
    sout << "ERROR :Invalid argument for --nauto (must be > 1).\n" << err_help ;
    exit(-1);
  }
  if(rm_indivs && keep_indivs ){
    sout << "ERROR :Cannot use both --keep and --remove.\n" << err_help ;
    exit(-1);
  }
  if(rm_snps && keep_snps ){
    sout << "ERROR :Cannot use both --extract and --exclude.\n" << err_help ;
    exit(-1);
  }

  if( test_mode && select_chrs && std::count( chrKeep_test.begin(), chrKeep_test.end(), -1) ){
    sout << "ERROR :Invalid chromosome specified to be tested.\n" << err_help ;
    exit(-1);
  }

  if(test_mode && !skip_blups && blup_file == "NULL") {
    sout << "ERROR :You must specify --pred if using --step 2 (otherwise use --ignore-pred).\n" << err_help ;
    exit(-1);
  }

  if( n_genofiles < 1 ) {
    sout << "ERROR :You must supply an input file using one of --bed or --bgen.\n" << err_help ;
    exit(-1);
  }
  if( n_genofiles > 1 ){
    sout << "ERROR :You must use either --bed or --bgen but not both.\n" << err_help ;
    exit(-1);
  }

  if( bgenSample && (file_type != "bgen") ) {
    sout << "ERROR :You must supply a BGEN file corresponding to the sample file (use `--bgen`).\n" << err_help ;
    exit(-1);
  }
  if(block_size < 1) {
    sout << "ERROR : You must set --b.\n" << err_help ;
    exit(-1);
  }
  if((file_type == "bgen") & (!file_exists (bgen_file))) {
    sout << "ERROR : " << bgen_file  << " doesn't exist.\n" << err_help ;
    exit(-1);
  }
  if((cov_file != "NULL") & (!file_exists (cov_file))) {
    sout << "ERROR : " << cov_file  << " doesn't exist.\n" << err_help ;
    exit(-1);
  }
  if((pheno_file != "NULL") & (!file_exists (pheno_file))) {
    sout << "ERROR : " << pheno_file  << " doesn't exist.\n" << err_help ;
    exit(-1);
  }
  if((file_type == "bed") & (!file_exists (bedfile))) {
    sout << "ERROR : " << bedfile  << " doesn't exist.\n" << err_help ;
    exit(-1);
  }
  if((file_type == "bed") & (!file_exists (famfile))) {
    sout << "ERROR : " << famfile  << " doesn't exist.\n"  << err_help;
    exit(-1);
  }
  if((file_type == "bed") & (!file_exists (bimfile))) {
    sout << "ERROR : " << bimfile  << " doesn't exist.\n" << err_help ;
    exit(-1);
  }

}

void Data::print_header(std::ostream& o){

  o << left << std::setw(14) << " " << "|==================================|" << endl;
  o << left << std::setw(14) << " " << "|           REGENIE v" << left << std::setw(14) << VERSION_NUMBER << "|" << endl;
  o << left << std::setw(14) << " " << "|==================================|" << endl << endl;

  o << "Copyright (c) 2020 Joelle Mbatchou and Jonathan Marchini." << endl;
  o << "Distributed under the MIT License.\n\n";
}

void Data::start_log(int argc, char **argv){
  string log_name = out_file + ".log";
  sout.coss.open(log_name.c_str(), ios::out | ios::trunc); 

  runtime.init();
  sout << "Start time: " << ctime(&runtime.start_time_info) << endl; 
  print_header(sout.coss);
  print_header(cout);
  sout << "Log of output saved in file : " << log_name << endl<< endl;

  // print options
  sout << "Command line arguments:";
  for(size_t counter=1;counter<argc;counter++){	  
    if( string(argv[counter])[0] == '-') sout << endl << "  ";
    sout << string(argv[counter]) << " ";
  }
  sout << endl << endl;

}

void Data::set_ridge_params(int nparams, vector<double> &param){

  if(nparams < 2){
    sout << "ERROR : Number of ridge parameters must be at least 2\n" << err_help;
    exit(-1);
  } else {
    // endpoints are 0.01 and 0.99 
    double step = 1.0 / ( nparams - 1 );
    double val = step;
    param.resize( nparams);

    for( size_t index_p = 1; index_p < (nparams - 1); index_p++, val += step) param[index_p] = val;
    param[0] = 0.01;
    param[nparams-1] = 0.99;

  }
}

void Data::print_usage_info(){

  double total_ram;
  string ram_unit;

  ///// Memory usage
  if(!test_mode){
    // Step 1
    // 4P + max( B + PRT, PRT) + #chrs [P:#traits;R=#ridge l0;T=#predictions from l0]
    total_ram = 4 * n_pheno + nChrom;
    int t_eff = ( write_l0_pred ? 1 : total_n_block );
    int p_eff = ( write_l0_pred ? 1 : n_pheno );
    total_ram += std::max( block_size + n_pheno * n_ridge_l0 * t_eff, p_eff * n_ridge_l0 * total_n_block );
  } else {
    // Step 2
    // 3P + B
    total_ram = n_pheno * 3 + block_size; // y, mask, y_resid, g
    if(binary_mode) {
      total_ram += 2 * n_pheno + block_size; // y_raw, gamma_hat, g_resid
      if(use_SPA) total_ram += 0.5 * block_size; // non_zero_indices of g (4 bytes)
    }
  }

  total_ram *= n_samples * sizeof(double);
  if( use_loocv ) total_ram += chunk_mb * 1e6; // max amount of memory used for LOO computations involved
  total_ram /= 1024.0 * 1024.0; 
  if( total_ram > 1000 ) {
    total_ram /= 1024.0; 
    ram_unit = "GB";
  } else ram_unit = "MB";

  int ram_int = (int) ceil( total_ram );
  sout << " * approximate memory usage : " << ram_int << ram_unit << endl;

  ///// Disk space usage
  if(!test_mode && write_l0_pred){
    if(loco_tmp_prefix == "NULL") loco_tmp_prefix = out_file;
    sout << " * writing level 0 predictions to disk" << endl;
    sout << "   -temporary files will have prefix [" << loco_tmp_prefix << "_l0_Y]" << endl;
    // N*P*T*R
    total_ram = n_pheno * total_n_block * n_ridge_l0;
    total_ram *= n_samples * sizeof(double);
    total_ram /= 1024.0 * 1024.0; 
    if( total_ram > 1000 ) {
      total_ram /= 1024.0; 
      ram_unit = "GB";
    } else ram_unit = "MB";
    int ram_int = (int) ceil( total_ram );
    sout << "   -approximate disk space needed : " << ram_int << ram_unit << endl;
  }
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          read in files
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::file_read_initialization() {

  if( rm_snps ) set_snps_to_rm();
  else if( keep_snps ) set_snps_to_keep();

  // read genotype data
  chr_counts.resize(nChrom, 0.0);
  if(rm_snps || keep_snps) chr_file_counts.resize(nChrom, 0.0);

  if(file_type == "bed") read_bed_bim_fam();
  else prep_bgen();

  if(n_variants == 0){
    sout << "ERROR: No variant left to include in analysis.\n";
    exit(-1);
  }

  //  keep track of samples in geno/pheno/covariates 
  //   and those to ignore
  ind_ignore = ArrayXd::Zero( n_samples );
  ind_in_analysis = ArrayXd::Zero( n_samples );
  ind_in_pheno_and_geno = ArrayXd::Zero( n_samples );
  ind_in_cov_and_geno = ArrayXd::Zero( n_samples );

  if( rm_indivs ) set_IDs_to_rm();
  else if( keep_indivs ) set_IDs_to_keep();
  
}

void Data::prep_bgen(){

  uint32_t nOutofOrder = 0;
  std::vector< int > chr_read ; 
  std::string chromosome, rsid;
  uint32_t position ;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;
  std::vector< string > tmp_ids ;
  snp tmp_snp;

  sout << left << std::setw(20) << " * bgen" << ": [" << bgen_file << "]" << endl;
  // open file and print file info
  bgen_tmp.open( bgen_file ) ;
  bgen_tmp.summarise( sout.coss ) ;
  bgen_tmp.summarise( cerr ) ;

  // retrieve number of samples/snps
  n_samples  = bgen_tmp.number_of_samples();
  if(!(keep_snps||rm_snps)) n_variants = bgen_tmp.number_of_variants();
  else n_variants = 0;

  // get sample IDs (from sample file or directly from bgen file)
  if( bgenSample ) {
    read_bgen_sample(tmp_ids);
  } else {
    bgen_tmp.get_sample_ids(
        [&tmp_ids]( std::string const& id ) { tmp_ids.push_back( id ) ; }
        ) ;
  }

  // check duplicates -- if not, store in map
  for(size_t i = 0; i < n_samples; i++) {
    if (FID_IID_to_ind.find(tmp_ids[i]) != FID_IID_to_ind.end()) {
      sout << "ERROR: Duplicate individual in .bgen file : FID_IID=" << tmp_ids[i] << endl;
      exit(1);
    }
    FID_IID_to_ind.insert( std::make_pair( tmp_ids[i], i ) );
  }

  // get CPRA info for each SNP
  tmp_snp.offset = bgen_tmp.get_position();
  while(bgen_tmp.read_variant( &chromosome, &position, &rsid, &alleles )) {

    bgen_tmp.ignore_probs();
    assert(alleles.size() == 2) ; // only bi-allelic allowed

    tmp_snp.chrom = chrStrToInt(chromosome);
    if (tmp_snp.chrom == -1) {
      sout << "ERROR: Unknown chromosome code in bgen file."<< endl;
      exit(1);
    }
    if( chr_read.empty() || (tmp_snp.chrom != chr_read.back()) ) chr_read.push_back(tmp_snp.chrom);

    tmp_snp.physpos = position;
    tmp_snp.ID = rsid;
    tmp_snp.allele1 = alleles[1];
    tmp_snp.allele2 = alleles[0]; // switch so allele0 is ALT

    // check if variant is in list of retained/ignored SNPs
    if(keep_snps) tmp_snp.mask = (std::find(snplist_to_keep.begin(), snplist_to_keep.end(), tmp_snp.ID) == snplist_to_keep.end());
    else if(rm_snps) tmp_snp.mask = (std::find(snplist_to_rm.begin(), snplist_to_rm.end(), tmp_snp.ID) != snplist_to_rm.end());
    // track how many variants are to be included for each chr
    if(!tmp_snp.mask) {
      chr_counts[tmp_snp.chrom-1]++; 
      if(keep_snps || rm_snps) n_variants++;
    }
    // track how many variants are in file for each chr
    if(keep_snps || rm_snps) chr_file_counts[tmp_snp.chrom-1]++; 

    // check if snps are in order (same chromosome & non-decreasing positions)
    if (!snpinfo.empty() && (tmp_snp.chrom == snpinfo.back().chrom) && ( (tmp_snp.physpos < snpinfo.back().physpos) || (tmp_snp.genpos < snpinfo.back().genpos) )) nOutofOrder++;

    snpinfo.push_back(tmp_snp);

    tmp_snp.offset = bgen_tmp.get_position();
  }

  if (!test_mode && (nOutofOrder > 0)) sout << "WARNING: Total number of snps out-of-order in bgen file : " << nOutofOrder << endl;

  // go through each chromosome in order & save number of snps
  // also save how many are actually read (when using --exclude)
  vector<int> tmp_v;
  tmp_v.resize(3, 0);
  for(int j = 0; j < chr_read.size(); j++){
    int i = chr_read[j];
    tmp_v[0] = chr_counts[i-1];
    if(keep_snps || rm_snps) tmp_v[2] = chr_file_counts[i-1];
    chr_map.insert(pair<int, vector<int> >(i, tmp_v)); 
  }

  if(keep_snps) {
    sout << "   -keeping only variants specified in [" << file_snps_include << "]" << endl;
    sout << "   -number of variants remaining in the analysis = " << n_variants << endl;
    if(snpinfo.size() == n_variants) keep_snps = false;
  } else if(rm_snps) {
    sout << "   -removing variants specified in [" << file_snps_exclude << "]" << endl;
    sout << "   -number of variants remaining in the analysis = " << n_variants << endl;
    if(snpinfo.size() == n_variants) rm_snps = false;
  }

  if( !streamBGEN ) {
    // setup file for reading the genotype probabilities later
    bgen.open( bgen_file ) ;
  }

}

void Data::read_bgen_sample(std::vector<string> &ids){

  ifstream myfile;
  int nline = 0;
  string FID, IID, line, tmp_str;

  sout << "   -sample file: " << sample_file << endl;
  myfile.open (sample_file, ios::in);
  if (!myfile.is_open()) {    
    sout << "ERROR: Cannot open sample file : " << sample_file << endl;
    exit(-1);
  }

  // read fid/iid information
  while (getline (myfile,line)) {
    std::istringstream iss(line);

    if( !(iss >> FID >> IID) ){
      sout << "ERROR: Incorrectly formatted sample file at line" << ids.size() + 1 << endl;
      exit(-1);
    }

    // check first two lines for correct format
    if(nline == 0){
      if( (FID != "ID_1") || (IID != "ID_2") ) {
        sout << "ERROR: Header of the sample file must start with: ID_1 ID_2" << endl;
        exit(1);
      }
    } else if(nline == 1){
      if( (FID != "0") || (IID != "0") ) {
        sout << "ERROR: Second line of sample file must start with: 0 0" << endl;
        exit(1);
      }
    } else {
      tmp_str = FID + "_" + IID;
      ids.push_back(tmp_str);
    }

    nline++;
  }

  if( n_samples != ids.size() ){
    sout << "ERROR: Number of samples in BGEN file does not match that in the sample file." << endl;
    exit(-1);
  }

  myfile.close();
}

void Data::read_bed_bim_fam() {

  read_bim();
  sout << "n_snps = " << snpinfo.size() << endl;

  if(keep_snps) {
    sout << "   -keeping only variants specified in [" << file_snps_include << "]" << endl;
    sout << "   -number of variants remaining in the analysis = " << n_variants << endl;
    if(snpinfo.size() == n_variants) keep_snps = false;
  } else if(rm_snps) {
    sout << "   -removing variants specified in [" << file_snps_exclude << "]" << endl;
    sout << "   -number of variants remaining in the analysis = " << n_variants << endl;
    if(snpinfo.size() == n_variants) rm_snps = false;
  }

  read_fam();
  sout << "n_samples = " << n_samples << endl;
  prep_bed(bedfile);

}


void Data::read_bim() {

  uint32_t nOutofOrder = 0;
  int minChr_read = 0; // enforce that chromosomes in file are sorted
  std::vector< int > chr_read ; 
  std::vector< string > tmp_str_vec ;
  snp tmp_snp; 
  string line;
  ifstream myfile;

  sout << left << std::setw(20) << " * bim" << ": [" << bimfile << "] " << flush;
  myfile.open(bimfile.c_str());
  if (!myfile.is_open()) {    
    sout << "ERROR: Cannot open bim file : " << bimfile << endl;
    exit(1);
  }

  while (getline(myfile, line)) {
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 6 ){
      sout << "ERROR: Incorrectly formatted bim file at line " << snpinfo.size()+1 << endl;
      exit(-1);
    }

    tmp_snp.chrom = chrStrToInt(tmp_str_vec[0]);
    tmp_snp.ID = tmp_str_vec[1];
    tmp_snp.genpos = std::stod( tmp_str_vec[2]);
    tmp_snp.physpos = std::stoul( tmp_str_vec[3],nullptr,0);
    // take ref allele as last
    tmp_snp.allele2 = tmp_str_vec[4];
    tmp_snp.allele1 = tmp_str_vec[5];

    if (tmp_snp.chrom == -1) {
      sout << "ERROR: Unknown chromosome code in bim file at line " << snpinfo.size()+1 << endl;
      exit(1);
    }
    if( chr_read.empty() || (tmp_snp.chrom != chr_read.back() ) ) {
      chr_read.push_back(tmp_snp.chrom);
      if( tmp_snp.chrom <= minChr_read ){
        sout << "ERROR: Chromosomes in .bim file are not in ascending order.\n";
        exit(-1);
      } else minChr_read = tmp_snp.chrom;
    }
    
    // check if variant is in list of ignored SNPs
    if(keep_snps) tmp_snp.mask = (std::find(snplist_to_keep.begin(), snplist_to_keep.end(), tmp_snp.ID) == snplist_to_keep.end());
    else if(rm_snps) tmp_snp.mask = (std::find(snplist_to_rm.begin(), snplist_to_rm.end(), tmp_snp.ID) != snplist_to_rm.end());

    // track how many variants are to be included for chr
    if(!tmp_snp.mask) {
      chr_counts[tmp_snp.chrom-1]++; 
      n_variants++;
    }
    // track how many variants in file for chr
    if(keep_snps || rm_snps) chr_file_counts[tmp_snp.chrom-1]++; 

    // check if snps are in order (same chromosome & non-decreasing positions)
    if (!snpinfo.empty() && (tmp_snp.chrom == snpinfo.back().chrom) && ( (tmp_snp.physpos < snpinfo.back().physpos) || (tmp_snp.genpos < snpinfo.back().genpos) )) nOutofOrder++;

    snpinfo.push_back(tmp_snp);
  }

  if (!test_mode && (nOutofOrder > 0)) sout << "WARNING: Total number of snps out-of-order in bim file : " << nOutofOrder << endl;

  // go through each chromosome in order & save number of snps
  // also save how many are actually read (when using --exclude)
  vector<int> tmp_v;
  tmp_v.resize(3, 0);
  for(int j = 0; j < chr_read.size(); j++){
    int i = chr_read[j];
    tmp_v[0] = chr_counts[i-1];
    if(keep_snps || rm_snps) tmp_v[2] = chr_file_counts[i-1];
    chr_map.insert(pair<int, vector<int> >(i, tmp_v)); 
  }

  myfile.close();

}


void Data::read_fam() {

  int sex, lineread = 0; 
  double yval;
  string line, tmp_id;
  std::vector< string > tmp_str_vec ;
  ifstream myfile;

  sout << left << std::setw(20) << " * fam" << ": [" << famfile << "] ";
  myfile.open(famfile.c_str());
  if (!myfile.is_open()) {    
    sout << "ERROR: Cannot open fam file : " << famfile << endl;
    exit(1);
  }

  while (getline(myfile, line)) {
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 6 ){
      sout << "ERROR: Incorrectly formatted fam file at line " << lineread + 1 << endl;
      exit(-1);
    }

    tmp_id = tmp_str_vec[0] + "_" + tmp_str_vec[1];

    // check duplicates -- if not, store in map
    if (FID_IID_to_ind.find(tmp_id ) != FID_IID_to_ind.end()) {
      sout << "ERROR: Duplicate individual in fam file : FID_IID=" << tmp_id << endl;
      exit(1);
    }
    FID_IID_to_ind.insert( std::make_pair( tmp_id, lineread ) );

    lineread++;
  }
  n_samples = lineread;

  myfile.close();
}

void Data::prep_bed(string bedfile) {

  sout << left << std::setw(20) << " * bed" << ": [" << bedfile << "]" << endl;
  bed_ifstream.open(bedfile.c_str(), std::ios::in | std::ios::binary);
  if (!bed_ifstream.is_open()) {    
    sout << "ERROR: Cannot open bed file : " << bedfile << endl;
    exit(1);
  }

  uchar header[3];
  bed_ifstream.read( reinterpret_cast<char *> (&header[0]), 3);
  if ( (header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01) ) {
    sout << "ERROR: Incorrect magic number in bed file : " << bedfile << endl;
    exit(1);
  }

  // size of genotype block [(n+3)/4 = ceil(n/4.0)]
  bed_block_size = (n_samples+3)>>2;
  inbed.resize( bed_block_size );
}

// snps to retain in step 1 analysis
void Data::set_snps_to_keep() {

  ifstream myfile;
  string line;
  std::vector< string > tmp_str_vec ;

  myfile.open (file_snps_include.c_str(), ios::in);
  if (!myfile.is_open()) {    
    sout << "ERROR: Cannot open file specified by --exclude :" << file_snps_include<< endl;
    exit(-1);
  }

  while( getline (myfile,line) ){
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 1 ){
      sout << "ERROR: Incorrectly formatted file specified by --exclude." << endl;
      exit(-1);
    }
    snplist_to_keep.push_back( tmp_str_vec[0] );
  }

}

// snps to exclude from step 1 analysis
void Data::set_snps_to_rm() {

  ifstream myfile;
  string line;
  std::vector< string > tmp_str_vec ;

  myfile.open (file_snps_exclude.c_str(), ios::in);
  if (!myfile.is_open()) {    
    sout << "ERROR: Cannot open file specified by --exclude :" << file_snps_exclude<< endl;
    exit(-1);
  }

  while( getline (myfile,line) ){
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 1 ){
      sout << "ERROR: Incorrectly formatted file specified by --exclude." << endl;
      exit(-1);
    }
    snplist_to_rm.push_back( tmp_str_vec[0] );
  }

}

void Data::set_IDs_to_keep() {

  ifstream myfile;
  string line;
  std::vector< string > tmp_str_vec ;
  findID person;
  uint64 indiv_index;
  ArrayXd is_included = ArrayXd::Zero(n_samples);

  // track individuals to include -> remaining are ignored
  myfile.open (file_ind_include.c_str(), ios::in);
  if (!myfile.is_open()) {    
    sout << "ERROR: Cannot open file specified by --keep :" << file_ind_include<< endl;
    exit(-1);
  }

  while( getline (myfile,line) ){
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 2 ){
      sout << "ERROR: Incorrectly formatted file specified by --keep." << endl;
      exit(-1);
    }

    person = getIndivIndex(tmp_str_vec[0], tmp_str_vec[1]);
    if(!person.is_found) continue;

    indiv_index = person.index;
    is_included(indiv_index) = 1;
  }

  // check size
  if( is_included.sum() < 1 ) {
    sout << "ERROR: None of the individuals specified by --keep are in the genotype file.\n";
    exit(-1);
  }
  sout << "   -number of genotyped individuals to keep in the analysis = " << is_included.sum() << endl;

  ind_ignore = 1 - is_included; 

}

void Data::set_IDs_to_rm() {

  ifstream myfile;
  string line;
  std::vector< string > tmp_str_vec ;
  findID person;
  uint64 indiv_index;

  // specify individuals to exclude by -1
  myfile.open (file_ind_exclude.c_str(), ios::in);
  if (!myfile.is_open()) {    
    sout << "ERROR: Cannot open file specified by --remove :" << file_ind_exclude<< endl;
    exit(-1);
  }

  while( getline (myfile,line) ){
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 2 ){
      sout << "ERROR: Incorrectly formatted file specified by --remove." << endl;
      exit(-1);
    }

    person = getIndivIndex(tmp_str_vec[0], tmp_str_vec[1]);
    if(!person.is_found) continue;

    indiv_index = person.index;
    ind_ignore(indiv_index) = 1;
  }

  sout << "   -number of genotyped individuals to exclude from the analysis = " << ind_ignore.sum() << endl;

}

void Data::read_pheno_and_cov() {

  if( strict_mode ) sout << " * dropping observations with missing values at any of the phenotypes" << endl;
  else if( !rm_missing_qt  && !binary_mode) sout << " * keeping and mean-imputing missing observations (done for each trait)" << endl;

  // read in phenotype (mean-impute for QT)
  pheno_read();

  // Intercept
  new_cov = MatrixXd::Ones(n_samples, 1);
  if(strict_mode) new_cov.array() *= masked_indivs.col(0).array();

  // read in covariates
  if(cov_file != "NULL") covariate_read();
  else ind_in_cov_and_geno = ArrayXd::Ones( n_samples );

  // mask individuals 
  ind_in_analysis = ind_in_pheno_and_geno * ind_in_cov_and_geno;
  masked_indivs.array().colwise() *= ind_in_analysis;
  if( strict_mode ) ind_in_analysis *= masked_indivs.col(0).array();
  phenotypes.array().colwise() *= ind_in_analysis;
  if(binary_mode) phenotypes_raw.array().colwise() *= ind_in_analysis;
  new_cov.array().colwise() *= ind_in_analysis;
  Neff = masked_indivs.colwise().sum();

  // check sample size
  if( ind_in_analysis.sum() < 1 ) {
    sout << "ERROR: Sample size cannot be < 1\n";
    exit(-1);
  }
  sout << " * number of individuals used in analysis = " << ind_in_analysis.sum() << endl;


  // orthonormal basis
  getCovBasis();

  // compute offset for BT (only in step 1)
  if(binary_mode && !test_mode) fit_null_logistic(0);

  // residualize phenotypes (skipped for BTs when testing)
  if(!test_mode || !binary_mode) residualize_phenotypes();
}

void Data::pheno_read() {

  ifstream myfile;
  string line;
  std::vector< string > tmp_str_vec ;
  findID person;
  uint64 indiv_index;
  vector<bool> pheno_colKeep;
  bool keep_pheno;
  double mean;

  sout << left << std::setw(20) << " * phenotypes" << ": [" << pheno_file << "] ";
  myfile.open (pheno_file.c_str(), ios::in);
  if (!myfile.is_open()) {    
    sout << "ERROR: Cannot open phenotype file : " << pheno_file << endl;
    exit(-1);
  }

  getline (myfile,line); // header
  boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

  // check that FID and IID are first two entries in header
  if( (tmp_str_vec[0] != "FID") || (tmp_str_vec[1] != "IID") ) {
    sout << "ERROR: Header of phenotype file must start with: FID IID" << endl;
    exit(-1);
  }

  // get phenotype names 
  n_pheno = 0;
  for( size_t filecol = 2; filecol < tmp_str_vec.size(); filecol++ ) {
    // default is to use all columns  
    if(!select_phenos){
      pheno_colKeep.push_back( true );
      n_pheno++;
      pheno_names.push_back( tmp_str_vec[filecol] );
    } else {
      // otherwise, check if phenotype is in list of phenotypes to keep
      keep_pheno = std::find(pheno_colKeep_names.begin(), pheno_colKeep_names.end(), tmp_str_vec[filecol]) != pheno_colKeep_names.end();
      pheno_colKeep.push_back( keep_pheno );
      if(keep_pheno){
        n_pheno++;
        pheno_names.push_back( tmp_str_vec[filecol] );
      }
    }
  }

  // check #pheno is > 0
  if(n_pheno < 1){
    sout << "Need at least one phenotype." << endl;
    exit(-1);
  }
  sout << "n_pheno = " << n_pheno << endl;

  // allocate memory
  phenotypes = MatrixXd::Zero(n_samples, n_pheno);
  masked_indivs = MatrixXd::Ones(n_samples, n_pheno);
  // for raw binary traits
  if (binary_mode) {
    phenotypes_raw = MatrixXd::Zero(n_samples, n_pheno);
    if(!test_mode)
      offset_logreg = MatrixXd::Zero(n_samples, n_pheno); 
  }
  VectorXd total, ns;
  total.setZero(n_pheno);
  ns.setZero(n_pheno);

  // read in data
  while( getline (myfile,line) ){
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < pheno_colKeep.size() ){
      sout << "ERROR: Incorrectly formatted phenotype file." << endl;
      exit(-1);
    }

    person = getIndivIndex(tmp_str_vec[0], tmp_str_vec[1]);
    if(!person.is_found) continue;

    indiv_index = person.index;

    // ignore sample if it is in exlusion list
    if( ind_ignore(indiv_index) ) continue;

    // check duplicate
    if( !ind_in_pheno_and_geno(indiv_index) ){
      ind_in_pheno_and_geno( indiv_index ) = 1;
    } else {
      sout << "ERROR: Individual appears more than once in phenotype file: FID=" << tmp_str_vec[0] << " IID=" << tmp_str_vec[1] << endl;
      exit(-1);
    }

    // read phenotypes 
    for(int i_pheno = 0, j = 0; j < pheno_colKeep.size(); j++) {

      if( !pheno_colKeep[j] ) continue;

      phenotypes(indiv_index, i_pheno) = convertDouble(tmp_str_vec[2+j]);

      // for BT, save raw data and force 0/1 values
      if (binary_mode) {
        if((control_code == 1) && (phenotypes(indiv_index, i_pheno) != missing_value_double)) phenotypes(indiv_index, i_pheno) -= control_code; // if using 1/2/NA encoding
        phenotypes_raw(indiv_index, i_pheno) = phenotypes(indiv_index, i_pheno);
        if(fabs(phenotypes_raw(indiv_index, i_pheno)) > numtol && fabs(phenotypes_raw(indiv_index, i_pheno)-1) > numtol ) {
          if(within_sample_l0){
            sout << "ERROR: No missing value allowed in phenotype file with option -within" << endl;
            exit(-1);
          } else if( phenotypes_raw(indiv_index, i_pheno) != missing_value_double ) {
            sout << "ERROR: A phenotype value is not 0/1/NA for individual: FID=" << tmp_str_vec[0] << " IID=" << tmp_str_vec[1] << endl;
            sout << "Use flag '--1' for 1/2/NA encoding [1=control|2=case|NA=missing]." << endl;
            exit(-1);
          }
          phenotypes_raw(indiv_index, i_pheno) = missing_value_double;
          masked_indivs(indiv_index, i_pheno) = 0;
          if( strict_mode ) masked_indivs.row(indiv_index) = MatrixXd::Zero(1, n_pheno);
        }
      }

      if( phenotypes(indiv_index, i_pheno) != missing_value_double ) {
        total(i_pheno) +=  phenotypes(indiv_index, i_pheno);
        ns(i_pheno) +=  1;
      } else {
        if( test_mode && rm_missing_qt ) masked_indivs(indiv_index, i_pheno) = 0;
        if( strict_mode ) masked_indivs.row(indiv_index) = MatrixXd::Zero(1, n_pheno);
      }

      i_pheno++;
    }

  }

  // mask individuals in genotype data but not in phenotype data
  masked_indivs.array().colwise() *= ind_in_pheno_and_geno;

  // check if all individuals have missing/invalid phenotype
  if(masked_indivs.colwise().sum().array().minCoeff() == 0){
    sout << "ERROR: All individuals have missing/invalid phenotype values." << endl;
    exit(-1);
  }

  if(!binary_mode || !test_mode){

    if(!binary_mode){
      // impute missing with mean
      for(size_t i = 0; i < n_samples;i++) 
        for(size_t j = 0; j < n_pheno;j++) {
          if( phenotypes(i,j) != missing_value_double ) {
            phenotypes(i,j) -= total(j) / ns(j);	  
          }  else phenotypes(i,j) = 0.0;
        }
    } else {
      for(size_t j = 0; j < n_pheno; j++) {
        mean = (masked_indivs.col(j).array() == 1).select( phenotypes.col(j).array(), 0).sum() / masked_indivs.col(j).sum();
        phenotypes.col(j).array() = (masked_indivs.col(j).array() == 1).select(phenotypes.col(j).array() - mean, 0);
      }
    }

    // apply masking
    phenotypes.array() *= masked_indivs.array();

    sout <<  "   -reading phenotypes...done" << endl;
  }

  // number of phenotyped individuals 
  sout <<  "   -number of phenotyped individuals = " << ind_in_pheno_and_geno.sum() << endl;

  // check that there cases are present
  if(binary_mode){
    for(size_t j = 0; j < n_pheno;j++) {
      if( ( phenotypes_raw.col(j).array() == 1 ).count() == 0){
        sout << "ERROR: No cases present for phenotype: " << pheno_names[j] << endl; 
        exit(-1);
      }
    }
  }

  Neff = masked_indivs.colwise().sum();
  if(strict_mode) sout << "   -number of individuals remaining with non-missing phenotypes = " << Neff(0) << endl;

  myfile.close();

}

void Data::covariate_read() {
  ifstream myfile;
  string line;
  std::vector< string > tmp_str_vec ;
  findID person;
  uint64 indiv_index;
  vector<bool> cov_colKeep;
  bool keep_cov;

  sout << left << std::setw(20) << " * covariates" << ": [" << cov_file << "] " << flush;
  myfile.open (cov_file.c_str(), ios::in);
  if (!myfile.is_open()) {
    sout << "ERROR: Cannot open covariate file : " << cov_file << endl;
    exit(-1);
  }

  getline (myfile,line);
  boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

  // check that FID and IID are first two entries in header
  if( (tmp_str_vec[0] != "FID") || (tmp_str_vec[1] != "IID") ) {
    sout << "ERROR: Header of covariate file must start with: FID IID" << endl;
    exit(-1);
  }

  // get covariate names 
  n_cov = 0;
  for( size_t filecol = 2; filecol < tmp_str_vec.size(); filecol++ ) {
    // default is to use all columns  
    if(!select_covs){
      cov_colKeep.push_back( true );
      n_cov++;
    } else {
      // otherwise, check if covariate is in list of covariates to keep
      keep_cov = std::find(cov_colKeep_names.begin(), cov_colKeep_names.end(), tmp_str_vec[filecol]) != cov_colKeep_names.end();
      cov_colKeep.push_back( keep_cov );
      if(keep_cov){
        n_cov++;
      }
    }
  }
  // check #covariates is > 0
  if(n_cov < 1){ // only intercept will be included
    sout << "n_cov = " << n_cov << " (+ intercept)" << endl;
    return ;
  }
  sout << "n_cov = " << n_cov << flush;

  // allocate memory 
  covariates.resize(n_samples, n_cov);

  // read in data
  while( getline (myfile,line) ){
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < cov_colKeep.size() ){
      sout << "ERROR: Incorrectly formatted covariate file." << endl;
      exit(-1);
    }

    person = getIndivIndex(tmp_str_vec[0], tmp_str_vec[1]);
    if(!person.is_found) continue;

    indiv_index = person.index;

    // ignore sample if it is in exlusion list
    if( ind_ignore(indiv_index) ) continue;

    // check duplicate
    if( !ind_in_cov_and_geno(indiv_index) ){
      ind_in_cov_and_geno(indiv_index) = 1;
    } else {
      sout << "ERROR: Individual appears more than once in covariate file: FID=" << tmp_str_vec[0] << " IID=" << tmp_str_vec[1] << endl;
      exit(-1);
    }

    // read covariate data and check for missing values
    for(size_t i_cov = 0, j = 0; j < cov_colKeep.size(); j++) {

      if( !cov_colKeep[j] ) continue;

      covariates(indiv_index, i_cov) = convertDouble(tmp_str_vec[2+j]);
      if( covariates(indiv_index, i_cov) == missing_value_double ) {
        sout << "ERROR: Individual has missing value in covariate file: FID=" << tmp_str_vec[0] << " IID=" << tmp_str_vec[1] << endl;
        exit(-1);
      }
      i_cov++;
    }

  }
  myfile.close();

  // mask individuals in genotype data but not in covariate data
  masked_indivs.array().colwise() *= ind_in_cov_and_geno;

  // Append covariates to intercept
  new_cov.resize(n_samples, 1 + n_cov);
  new_cov.col(0) = MatrixXd::Ones(n_samples, 1);
  new_cov.rightCols(n_cov) = covariates;

  // apply masking
  new_cov.array().colwise() *= ind_in_cov_and_geno;
  if(strict_mode) new_cov.array().colwise() *= masked_indivs.col(0).array();

  sout <<  endl;
  sout <<  "   -number of individuals with covariate data = " << ind_in_cov_and_geno.sum() << endl;

}

void Data::getCovBasis(){

    // eigen-decompose cov matrix
    MatrixXd xtx = new_cov.transpose() * new_cov;
    SelfAdjointEigenSolver<MatrixXd> es(xtx);
    VectorXd D = es.eigenvalues();
    MatrixXd V = es.eigenvectors();

    // create basis set
    // eigenvalues sorted in increasing order
    int non_zero_eigen = (D.array() > D.tail(1)(0) * eigen_val_rel_tol).count();
    RowVectorXd vv1 = D.tail(non_zero_eigen).array().sqrt();
    new_cov *= V.rightCols(non_zero_eigen);
    new_cov.array().rowwise() /= vv1.array();

}


double Data::convertDouble(const string &phenoValue){
  if(phenoValue == missing_pheno_str)
    return missing_value_double;

  double pheno_d;
  if(sscanf(phenoValue.c_str(), "%lf", &pheno_d) != 1){
    sout << "ERROR: Could not convert phenotype value to double: " << phenoValue << endl;
    exit(-1);
  }
  return pheno_d;
}

findID Data::getIndivIndex(const string &FID, const string &IID){

  string tmp_str;
  findID indiv;

  // get ID of individual
  tmp_str = FID + "_" + IID;

  // check individual is in genotype data
  if (FID_IID_to_ind.find(tmp_str) == FID_IID_to_ind.end()) {
    if(verbose) sout << "WARNING: Individual not present in genotype data: FID=" << FID << " IID=" << IID << endl;
    indiv.is_found = false;
  } else {
    indiv.is_found = true;
    indiv.index = FID_IID_to_ind[tmp_str];
  }

  return indiv;
}

int Data::chrStrToInt(string chrom) {

  //detect chromosome labels with 'chr' prefix
  regex re("^chr(\\d+)");
  smatch m;

  if (isdigit(chrom[0])) {
    int chr = atoi(chrom.c_str());
    if((chr >= 1) && (chr <= nChrom)) return chr;
  } else if (regex_search(chrom, m, re)){
    string dig = m[1];
    int chr = atoi(dig.c_str());
    if((chr >= 1) && (chr <= nChrom)) return chr;
  } else if ( (chrom == "X") || (chrom == "XY") || (chrom == "PAR1") || (chrom == "PAR2") || (chrom == "chrX") ) return nChrom;

  return -1;
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          adjust for covariates
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::residualize_phenotypes() {
  sout << "   -residualizing and scaling phenotypes...";
  auto t1 = std::chrono::high_resolution_clock::now();

  // residuals (centered) then scale
  MatrixXd beta = phenotypes.transpose() * new_cov;
  phenotypes -= ( (new_cov * beta.transpose()).array() * masked_indivs.array() ).matrix();
  scale_Y = phenotypes.colwise().norm().array() / sqrt(Neff.matrix().transpose().array() - 1);

  // check sd is not 0 
  if(scale_Y.minCoeff() < numtol){
    sout << "ERROR: At least one of the phenotypes has s.d. of 0." << endl;
    exit(-1);
  }
  phenotypes.array().rowwise() /= scale_Y.array();

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;
}

void Data::fit_null_logistic(int chrom) {

  sout << "   -fitting null logistic regression on binary phenotypes..." << flush;

  auto t1 = std::chrono::high_resolution_clock::now();

  int niter_cur;
  double dev_old, dev_new;
  ArrayXd Y1, hvec, score;
  ArrayXd betaold, betanew, etavec, pivec, wvec, zvec;
  MatrixXd  X1, XtW, XtWX;

  for(size_t i = 0; i < n_pheno; ++i ){

    Y1 = phenotypes_raw.col(i).array() * masked_indivs.col(i).array();
    X1 = (new_cov.array().colwise() * masked_indivs.col(i).array()).matrix();

    // starting values
    betaold = ArrayXd::Zero(new_cov.cols());
    betaold(0) = ( 0.5 + Y1.sum()) / ( Neff(i) + 1);
    betaold(0) = log( betaold(0) / (1 - betaold(0) ));
    if(test_mode) betaold(0) -= (blups.col(i).array() * masked_indivs.col(i).array()).mean();
    // compute deviance
    etavec = ( X1 * betaold.matrix()).array();
    if(test_mode) etavec += blups.col(i).array() * masked_indivs.col(i).array();
    pivec = 1 - 1 / (etavec.exp() + 1) ;
    dev_old = -2 * (Y1 * pivec.log() + (1-Y1) * (1-pivec).log() ).sum();

    niter_cur = 0;
    while(niter_cur++ < niter_max){

      // linear predictor = offset + X * beta
      etavec = ( X1 * betaold.matrix()).array();
      if(test_mode) etavec += blups.col(i).array() * masked_indivs.col(i).array();
      // fitted probabilities
      pivec = 1 - 1 / (etavec.exp() + 1) ;
      // diagonal matrix of sqrt( p*(1-p) )
      wvec = ( pivec * (1 - pivec) ).sqrt();

      // check none of the values are 0
      if( ( (masked_indivs.col(i).array() == 1) && (wvec == 0) ).count() > 0 ){
        sout << "ERROR: Zeros occured in Var(Y) during logistic regression! (Check covariates)" << endl;
        exit(-1);
      }
 
      XtW = X1.transpose() * wvec.matrix().asDiagonal();
      XtWX = XtW * XtW.transpose();
      // working vector z = X*beta + (Y-p)/(p*i(1-p)
      zvec = (etavec + (Y1 - pivec) / wvec.square()) * masked_indivs.col(i).array();
      if(test_mode) zvec -= blups.col(i).array() * masked_indivs.col(i).array();
      // parameter estimate
      betanew = ( XtWX ).colPivHouseholderQr().solve( XtW * wvec.matrix().asDiagonal() * zvec.matrix()).array();

      // start step-halving and stop when deviance decreases 
      for( size_t niter_search = 1; niter_search <= niter_max_line_search; niter_search++ ){

        etavec = ( X1 * betanew.matrix()).array();
        if(test_mode) etavec += blups.col(i).array() * masked_indivs.col(i).array();
        pivec = 1 - 1 / (etavec.exp() + 1) ;

        dev_new = -2 * (Y1 * pivec.log() + (1-Y1) * (1-pivec).log() ).sum();
        score = X1.transpose() * (Y1 - pivec).matrix(); 

        if( dev_new < dev_old + numtol ) break;
        // adjusted step size
        betanew = (betaold + betanew) / 2;
      }

      // stopping criterion
      if( score.abs().maxCoeff() < numtol) break;

      betaold = betanew;
      dev_old = dev_new;
    }

    // If didn't converge
    if(niter_cur > niter_max){
      sout << "WARNING: Logistic regression did not converge! (Increase --niter)\n";
      exit(-1);
    }
    // sout << "Converged in "<< niter_cur << " iterations." << endl;

    etavec = (X1 * betanew.matrix()).array();
    if(test_mode){
      etavec += blups.col(i).array() * masked_indivs.col(i).array();
      Y_hat_p.col(i) = (1 - 1 / (etavec.exp() + 1)).matrix() ;
      Gamma_sqrt.col(i) = (Y_hat_p.col(i).array() * (1 - Y_hat_p.col(i).array()) ).sqrt().matrix();
      XtW = ( X1.array().colwise() * Gamma_sqrt.col(i).array()).matrix();
      Xt_Gamma_X_inv[i] = (XtW.transpose() * XtW).colPivHouseholderQr().inverse();
    } else offset_logreg.col(i) = etavec;
  }

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;
}


void Data::residualize_genotypes() {

  sout << "   -residualizing and scaling genotypes..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  // mask missing individuals
  G.array().rowwise() *= ind_in_analysis.matrix().transpose().array();
  if(strict_mode) G.array().rowwise() *= masked_indivs.col(0).matrix().transpose().array();

  // residuals (centered)
  MatrixXd beta = G * new_cov;
  G -= beta * new_cov.transpose();

  // scaling
  if(strict_mode) scale_G = G.rowwise().norm() / sqrt(Neff(0) - 1);
  else scale_G = G.rowwise().norm() / sqrt(ind_in_analysis.sum() - 1);

  // only done in step 1
  if(scale_G.array().minCoeff() < numtol) {
    if(!test_mode) {
      sout << "!! Uh-oh, SNP with low variance.\n" ;
      exit(1);
    } else {
      bad_snps( scale_G.array() < numtol ) =  1;
      for(int i = 0; i < G.rows(); i++){
        // make snps polymorphic
        if(bad_snps(i) == 1) {
          G.row(i).head(20).array() += 1; // +1 to first 20 entries
          scale_G(i) = 1;
          if(verbose) sout << "WARNING: Ignoring SNP with low variance.\n";
        }
      }
    }  
  }

  G.array().colwise() /= scale_G.array();

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;

}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          prepare for level 0
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::set_blocks() {

  total_n_block = 0, total_chrs_loco = 0;
  int blocks_left = n_block;
  map<int, vector<int> >::iterator itr; 
  map<int, vector<int> > m1;

  // compute number of blocks for each chromosome
  for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) { 
    int chrom_nsnps = itr->second[0];
    int nb = ceil(chrom_nsnps * 1.0 / block_size);
    if(n_block > 0) {
      if(blocks_left > 0) {
        int minb = min(nb, blocks_left);
        //sout << << endl;
        itr->second[1] = minb;
        total_n_block += minb;
        blocks_left -= minb;
      }
    } else {
      itr->second[1] = nb;
      total_n_block += nb;
    }

    // track how many chromosome will have blups 
    if(itr->second[1] > 0) total_chrs_loco++;
    m1.insert(pair<int, vector<int> >(itr->first, itr->second));
  }
  chr_map = m1;
  //sout << "#chrs = "<< chr_map.size() << ";#loco chrs = "<< total_chrs_loco << endl;

  if(total_n_block == 0){
    sout << "ERROR: Total number of blocks must be > 0.\n";
    exit(-1);
  }

  // set ridge params 
  for(size_t i = 0; i < n_ridge_l0; i++) lambda[i] =  n_variants / lambda[i];
  for(size_t i = 0; i < n_ridge_l1; i++) { 
    if(!binary_mode)
      tau[i] =  (total_n_block *  n_ridge_l0)  / tau[i];
    else {
      // Assuming input tau[i] is total SNP heritability on the liability scale= m * 3/pi^2 * (1-h2) / h2
      tau[i] =  (total_n_block *  n_ridge_l0) * 3 / (M_PI * M_PI) * (1 - tau[i]) / tau[i];
    }
  }

  // for BTs: check if the sample size is lower than 5K (if so, use loocv)
  if( binary_mode && (ind_in_analysis.sum() < 5000) ) {
    if(!use_loocv){
      sout << "   -WARNING: Sample size is less than 5,000 so using LOOCV instead of " << cv_folds << "-fold CV\n.";
      use_loocv = true;
    }
  }

  /*
  // check block size vs sample size
  if(use_loocv && block_size > ind_in_analysis.sum()){
    sout << "ERROR: Block size must be smaller than the number of samples to perform LOOCV!" << endl;
    exit(-1);
  }
  */
  if(use_loocv) cv_folds = n_samples;

  uint64 neff_folds = use_loocv ? ind_in_analysis.sum() : cv_folds; 

  // summarize block sizes and ridge params
  sout << left << std::setw(20) << " * # threads" << ": [" << threads << "]" << endl;
  sout << left << std::setw(20) << " * block size" << ": [" << block_size << "]" << endl;
  sout << left << std::setw(20) << " * # blocks" << ": [" << total_n_block << "]" << endl;
  sout << left << std::setw(20) << " * # CV folds" << ": [" << neff_folds << "]" << endl;
  sout << left << std::setw(20) << " * ridge data_l0" << ": [" << n_ridge_l0 << " : ";
  for(size_t i = 0; i < n_ridge_l0; i++) sout << n_variants / lambda[i] << " ";
  sout << "]" <<endl;
  sout << left << std::setw(20) << " * ridge data_l1" << ": [" << n_ridge_l1 << " : ";
  for(size_t i = 0; i < n_ridge_l1; i++) {
    if(!binary_mode)
      sout << (total_n_block *  n_ridge_l0)  / tau[i] << " ";
    else 
      sout << (total_n_block *  n_ridge_l0) / ( (total_n_block *  n_ridge_l0) + (M_PI * M_PI) * tau[i] / 3 ) << " ";
  }
  sout << "]" << endl;

  // print approx. amount of memory needed
  print_usage_info();

  // if within sample predictions are used in level 1
  if (within_sample_l0) {
    sout << " * using within-sample predictions from level 0 as features at level 1" << endl;
  }


}


void Data::set_folds() {
  // set up folds
  cv_sizes.resize(cv_folds);

  if(strict_mode && !use_loocv){
    // assign folds accounting for missingness
    uint64 target_size_folds = floor( Neff(0) / cv_folds );
    if( target_size_folds < 1 ){
      sout << "ERROR: Not enough samples are present for " << cv_folds<<"-fold CV.\n";
      exit(-1);
    }

    uint64 n_non_miss = 0, cur_fold = 0, cum_size_folds = 0;
    for(size_t i = 0; i < n_samples; i++){

      if( masked_indivs(i, 0) ) n_non_miss++;
      if( n_non_miss == target_size_folds){
        cv_sizes[cur_fold] = i - cum_size_folds + 1;
        cum_size_folds += cv_sizes[cur_fold];
        n_non_miss = 0, cur_fold++;
      } else if( cur_fold == (cv_folds - 1) ){
        cv_sizes[cur_fold] = n_samples - i;
        break;
      }

      //sout << i << " " << cur_fold << " " << n_non_miss << " " << masked_indivs(i, 0) << " "<< target_size_folds << endl;
    }
  } else {

    // assign folds for individuals in analysis
    if( !use_loocv ){

      uint64 target_size_folds = floor( ind_in_analysis.sum() / cv_folds );
      if( target_size_folds < 1 ){
        sout << "ERROR: Not enough samples are present for " << cv_folds<<"-fold CV.\n";
        exit(-1);
      }

      uint64 n_non_miss = 0, cur_fold = 0, cum_size_folds = 0;
      for(size_t i = 0; i < n_samples; i++){

        if( ind_in_analysis(i) ) n_non_miss++;
        if( n_non_miss == target_size_folds){
          cv_sizes[cur_fold] = i - cum_size_folds + 1;
          cum_size_folds += cv_sizes[cur_fold];
          n_non_miss = 0, cur_fold++;
        } else if( cur_fold == (cv_folds - 1) ){
          cv_sizes[cur_fold] = n_samples - i;
          break;
        }

        //sout << i << " " << cur_fold << " " << n_non_miss << " " << ind_in_analysis(i) << " "<< target_size_folds << endl;
      }
    } else { 
      // loocv
      for(size_t i = 0; i < cv_folds; i++) cv_sizes[i] = 1;
    }

  }


  // check sd(Y) in folds
  if(!use_loocv && binary_mode){
    uint64 cum_size_folds = 0, index_i = 0;
    MatrixXd phenos = ( phenotypes_raw.array() * masked_indivs.array()).matrix();
    for(size_t i = 0; i < (cv_folds - 1); i++) {
      ArrayXd sum = phenos.block(cum_size_folds,0,cv_sizes[i],n_pheno).colwise().sum();
      ArrayXd n_cv = masked_indivs.block(cum_size_folds,0,cv_sizes[i],n_pheno).colwise().sum();
      ArrayXd sd_phenos = (sum/n_cv) * (1 - sum/n_cv);

      if( sd_phenos.minCoeff() < numtol ){
        sout << "ERROR: One of the folds has only cases/controls! Either use smaller #folds (option --cv) or use LOOCV (option --loocv).\n";
        exit(-1);
      }
      cum_size_folds += cv_sizes[i];
    }
  }

  // only used for K-fold CV
  if(!use_loocv){
    G_folds.resize(cv_folds);
    GtY.resize(cv_folds);	  
    if (!within_sample_l0){
      X_folds.resize(cv_folds);
      XtY.resize(cv_folds);	
    }
  }
}


void Data::setmem() {
  sout << " * setting memory..." << flush;

  G.resize(block_size, n_samples);
  set_folds();
  cumsum_values.resize(6);
  predictions.resize(1);
  predictions[0].resize(n_samples, total_chrs_loco);
  if (within_sample_l0) {
    pred_mat.resize(n_pheno);
    pred_pheno.resize(n_pheno);
  } else if(!use_loocv) beta_hat_level_1.resize(n_pheno);
  if(!use_loocv) {
    test_pheno.resize(n_pheno);
    test_mat.resize(n_pheno);
  } else test_mat_conc.resize(n_pheno);
  if(binary_mode){
    if (within_sample_l0) {
      pred_pheno_raw.resize(n_pheno);
      pred_offset.resize(n_pheno);
    }
    test_pheno_raw.resize(n_pheno);
    if(!use_loocv) test_offset.resize(n_pheno);
  }
  masked_in_folds.resize(cv_folds);
  if(print_block_betas) beta_print_out.resize(n_pheno);

  for(size_t i = 0; i < n_pheno; ++i ) {
    if (within_sample_l0) {
      pred_mat[i].resize(cv_folds);
      pred_pheno[i].resize(cv_folds);
    } else if(!use_loocv) beta_hat_level_1[i].resize(cv_folds);
    if(!use_loocv) {
      test_pheno[i].resize(cv_folds);
      test_mat[i].resize(cv_folds);
    } else test_mat_conc[i].resize(n_samples, n_ridge_l0 * ( write_l0_pred ? 1 : total_n_block) );
    if(binary_mode) {
      if (within_sample_l0) {
        pred_pheno_raw[i].resize(cv_folds);
        pred_offset[i].resize(cv_folds);
      }
      test_pheno_raw[i].resize(cv_folds);
      if(!use_loocv) test_offset[i].resize(cv_folds);
    }
    for(size_t j = 0; j < cv_folds; ++j ) {
      if (within_sample_l0) {
        pred_mat[i][j].resize(n_samples - cv_sizes[j], total_n_block * n_ridge_l0);
        pred_pheno[i][j].resize(n_samples - cv_sizes[j], 1);
      } else if(!use_loocv) beta_hat_level_1[i][j].resize(total_n_block * n_ridge_l0, n_ridge_l1);
      if(!use_loocv) {
        test_pheno[i][j].resize(cv_sizes[j], 1);
        test_mat[i][j].resize(cv_sizes[j], n_ridge_l0 * ( write_l0_pred ? 1 : total_n_block));
      }
      if(binary_mode) {
        if (within_sample_l0) {
          pred_pheno_raw[i][j].resize(n_samples - cv_sizes[j], 1);
          pred_offset[i][j].resize(n_samples - cv_sizes[j], 1);
        }
        test_pheno_raw[i][j].resize(cv_sizes[j], 1);
        if(!use_loocv) test_offset[i][j].resize(cv_sizes[j], 1);
      }
      if( i == 0 ) masked_in_folds[j].resize(cv_sizes[j], n_pheno);
    }
  }
  sout << "done" << endl << endl;  
}

void Data::calc_cv_matrices(int bs) {
  
  sout << "   -calc working matrices..." << flush;
  auto t2 = std::chrono::high_resolution_clock::now();

  if(!use_loocv){
    GGt.setZero(bs,bs);
    GTY.setZero(bs,n_pheno);
    uint64 cum_size_folds = 0;

    for( std::size_t i = 0; i < cv_folds; ++i ) {
      GtY[i] = G.block(0, cum_size_folds, bs, cv_sizes[i]) * phenotypes.block(cum_size_folds, 0, cv_sizes[i], n_pheno);
      G_folds[i] = G.block(0, cum_size_folds, bs, cv_sizes[i]) * G.block(0, cum_size_folds, bs, cv_sizes[i]).transpose();
      GGt += G_folds[i];
      GTY += GtY[i];
      cum_size_folds += cv_sizes[i];
    }
  } else { 
    GGt = G * G.transpose();
    GTY = G * phenotypes;
    SelfAdjointEigenSolver<MatrixXd> esG(GGt);
    GGt_eig_vec = esG.eigenvectors();
    GGt_eig_val = esG.eigenvalues();
    Wmat = GGt_eig_vec.transpose() * GTY;
  }

  sout << "done";
  auto t3 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
  sout << " (" << duration.count() << "ms) "<< endl;
}

void Data::level_0_calculations() {

  int block = 0, nread;
  snp_index_counter = 0;
  map<int, vector<int> >::iterator itr; 
  for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) { 
    int chrom = itr->first;
    int chrom_nsnps = itr->second[0];
    int chrom_nb = itr->second[1];
    if(chrom_nb > 0) sout << "Chromosome " << chrom << endl;
    //sout << "Ns="<< chrom_nsnps << ";Nfile="<< itr->second[2] << endl;

    // if all snps in chromosome are excluded
    // stream through file to read next chr 
    if((keep_snps || rm_snps) && (chrom_nb == 0)) {
      skip_snps(itr->second[2]);
      snp_index_counter += itr->second[2];
      continue;
    }

    if(keep_snps || rm_snps) nread = snp_index_counter;
    for(int bb = 0; bb < chrom_nb ; bb++) {

      int bs = block_size;
      if(bb == 0) G.resize(bs, n_samples);
      if((bb +1) * block_size > chrom_nsnps) {
        bs = chrom_nsnps - (bb * block_size) ;
        G.resize(bs,n_samples);
      }

      get_G(block, bs, chrom);

      // residualize and scale genotypes 
      residualize_genotypes();

      // calc working matrices for ridge regressions across folds
      calc_cv_matrices(bs);

      // calc level 0 ridge regressions
      if(use_loocv)
        ridge_level_0_loocv(block);
      else
        ridge_level_0(block);

      block++;
    }

    // if skipping all snps at end of chromosome
    if(keep_snps || rm_snps){
      nread = snp_index_counter - nread;
      //sout << "Nread =" << nread << "; chr=" << chrom << endl;
      if(nread < itr->second[2]) {
        skip_snps( itr->second[2] - nread );
        snp_index_counter += itr->second[2] - nread;
        //sout << "Nskipping=" << itr->second[2] - nread << endl;
      }
    }
  }

  // free up memory not used anymore
  G.resize(0,0); 
  if(write_l0_pred && (n_pheno > 1) ) {
    // free level 0 predictions for (P-1) indices in test_mat
    for(size_t ph = 1; ph < n_pheno; ++ph ) {
      if(!use_loocv){
        for(size_t i = 0; i < cv_folds; ++i ) test_mat[ph][i].resize(0,0);
        test_mat[ph].resize(0);
      } else {
        test_mat_conc[ph].resize(0,0);
      }
    }
  }

}
 
/////////////////////////////////////////////////
/////////////////////////////////////////////////
////      read genotype data in chunks
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::get_G(int block, int bs, int chrom){
  auto t1 = std::chrono::high_resolution_clock::now();

  // prepare vector to store non_zero indices (for SPA)
  if(use_SPA) { 
    non_zero_indices_G.resize(bs);
    for( std::size_t i = 0; i < non_zero_indices_G.size(); ++i ) 
      non_zero_indices_G[i].clear();
  }
  // set counts to 0
  if(htp_out){
    for( std::size_t i = 0; i < bs; ++i ) genocounts[i].setZero();
  }

  sout << " block [" << block+1 << "] : " << flush; 

  if(file_type == "bed") readChunkFromBedFileToG(bs);
  else readChunkFromBGENFileToG(bs, chrom);

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;
}

void Data::readChunkFromBGENFileToG(int bs, int chrom) {

  std::string chromosome, rsid;
  uint32_t position ;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;

  for(int snp = 0; snp < bs; ) {
    bgen.read_variant( &chromosome, &position, &rsid, &alleles );

    // skip probs if snp is in exclude file
    if(keep_snps || rm_snps){
      assert(snpinfo[snp_index_counter].ID == rsid);
      if(snpinfo[snp_index_counter].mask){
        bgen.ignore_probs();
        snp_index_counter++;
        continue;
      }
    }
    //sout << "["<< chrom << "]SNPid stored ("<< snpinfo[snp_index_counter].chrom <<") = " << snpinfo[snp_index_counter].ID<< "/ SNPIDread ("<<chromosome<<")= " << rsid << endl;

    assert(chrStrToInt(chromosome) == chrom);
    bgen.read_probs( &probs ) ;

    int ns = 0, hc_val = 0;
    bool switch_alleles;
    double ds, total = 0, info_num = 0; 
    for( std::size_t i = 0; i < probs.size(); ++i ) {
      ds = 0;
      for( std::size_t j = 1; j < probs[i].size(); ++j ) ds += probs[i][j] * j;

      if(ds != -3) {
        ds = 2 - ds; // switch so that allele0 is ALT
        if( ind_in_analysis(i) ){
          if( !strict_mode || (strict_mode && masked_indivs(i,0)) ){
            total += ds;
            info_num += 4 * probs[i][0] + probs[i][1] - ds * ds;
            ns++;
          }

          // get genotype counts (convert to hardcall)
          if( htp_out ) {
            hc_val = (int) (ds + 0.5); // round to nearest integer (0/1/2)
            if( !binary_mode ) {
              genocounts[snp].row(hc_val) += masked_indivs.row(i);
            } else {
              genocounts[snp].row(hc_val).array() += masked_indivs.row(i).array() * phenotypes_raw.row(i).array();
              genocounts[snp].row(3 + hc_val).array() += masked_indivs.row(i).array() * (1 - phenotypes_raw.row(i).array());
            }
          }
        }
      } 
      G(snp, i) = ds;

    }
    if( test_mode && ((total < min_MAC) || ((2 * ns - total) < min_MAC)) ) bad_snps(snp) = 1;
    //sout << "SNP#" << snp + 1 << "AC=" << total << " BAD="<< bad_snps(snp)<< endl;
    total /= ns;
    if(test_mode) {
      snp_afs(snp, 0) = total / 2;
      if( (snp_afs(snp, 0) == 0) || (snp_afs(snp, 0) == 1) ) snp_info(snp, 0) = 1;
      else snp_info(snp, 0) = 1 - info_num / (2 * ns * snp_afs(snp, 0) * (1 - snp_afs(snp, 0)));
    }

    if(use_SPA) { 
      // switch to minor allele
      switch_alleles = total > 1;
      if( test_type > 0) switch_alleles = false; // skip for DOM/REC test
      if(switch_alleles){
        G.row(snp).array() = ( G.row(snp).array() != -3).select( 2 - G.row(snp).array(), G.row(snp).array() );
        total = 2 - total;
        snp_flipped[snp] = true;
      } else snp_flipped[snp] = false;
    }

    // apply dominant/recessive encoding & recompute mean
    if(test_type > 0){
      for( std::size_t i = 0; i < probs.size(); ++i ) {
        if( (G(snp, i) != -3)  && ind_in_analysis(i) && 
            (!strict_mode || (strict_mode && masked_indivs(i,0))) ){
          if(test_type == 1){ //dominant
            G(snp, i) = probs[i][0] + probs[i][1]; // allele0 is ALT
          } else if(test_type == 2){ //recessive
            G(snp, i) = probs[i][0];
          }
        }
      }
      total = ((G.row(snp).transpose().array()!= -3) && (ind_in_analysis == 1)).select(G.row(snp).transpose().array(), 0).sum() / ns;
      if(total < numtol) bad_snps(snp) = 1;
    }

    // deal with missing data and center SNPs
    for( std::size_t i = 0; i < probs.size(); ++i ) {
      ds = G(snp, i);
      if( use_SPA && ind_in_analysis(i) && ds > 0 ) non_zero_indices_G[snp].push_back(i);

      if(ds != -3 && ind_in_analysis(i) && (!strict_mode || (strict_mode && masked_indivs(i,0)) ) ){
        G(snp, i) -= total;
      } else {
        G(snp, i) = 0;
      }
    }

    snp++;snp_index_counter++;
  }

  if(!verbose) sout << bs << " snps ";
}

void Data::readChunkFromBedFileToG(int bs) {

  int hc, ns, byte_start, bit_start;
  double total; 
  bool switch_alleles;
  // mapping matches the switch of alleles done when reading bim
  const int maptogeno[4] = {2, -3, 1, 0};

  // only for step 1
  for(size_t j = 0; j < bs; ) {
    if(keep_snps || rm_snps){
      if(snpinfo[snp_index_counter].mask){
        bed_ifstream.ignore(bed_block_size);
        snp_index_counter++;
        continue;
      }
    }

    ns = 0, total = 0;
    bed_ifstream.read( reinterpret_cast<char *> (&inbed[0]), bed_block_size);
    for (size_t i = 0; i < n_samples; i++) {
      byte_start = i>>2; // 4 samples per byte
      bit_start = (i&3)<<1; // 2 bits per sample
      hc = maptogeno[ (inbed[byte_start] >> bit_start)&3 ]; 
      G(j, i) = hc;

      if(hc != -3) {
        if( ind_in_analysis(i) ){
          if( !strict_mode || (strict_mode && masked_indivs(i,0)) ){
            total += hc;
            ns++;
          }
        }

        // get genotype counts 
        if( htp_out ) {
          if( !binary_mode ) genocounts[j].row(hc) += masked_indivs.row(i);
          else {
            genocounts[j].row(hc).array() += masked_indivs.row(i).array() * phenotypes_raw.row(i).array();
            genocounts[j].row(3 + hc).array() += masked_indivs.row(i).array() * (1 - phenotypes_raw.row(i).array());
          }
        }
      }
    }

    if( test_mode && ((total < min_MAC) || (( 2 * ns - total) < min_MAC)) )  bad_snps(j) = 1;
    total /= ns;
    if(test_mode) snp_afs(j, 0) = total / 2;

    if(use_SPA) { 
      // switch to minor allele
      switch_alleles = total > 1;
      if( test_type > 0) switch_alleles = false; // skip for DOM/REC test
      if(switch_alleles){
        G.row(j).array() = ( G.row(j).array() != -3).select( 2 - G.row(j).array(), G.row(j).array() );
        total = 2 - total;
        snp_flipped[j] = true;
      } else snp_flipped[j] = false;
    }

    // apply dominant/recessive encoding & recompute mean
    if(test_type > 0){
      if(test_type == 1){ //dominant
        G.row(j).array() = (G.row(j).array() == 2).select(1, G.row(j).array());
      } else if(test_type == 2){ //recessive
        G.row(j).array() = (G.row(j).array() >= 1).select(G.row(j).array() - 1, G.row(j).array());
      }
      total = ((G.row(j).transpose().array() != -3) && (ind_in_analysis == 1)).select(G.row(j).transpose().array(), 0).sum() / ns;
      if(total < numtol) bad_snps(j) = 1;
    }

    //if(j<5) sout << "\nj="<< j+1 << ":" <<  G.row(j).array().head(5);
    // deal with missing data and center SNPs
    for (size_t i = 0; i < n_samples; i++) {
      hc = G(j, i);
      if( use_SPA && (hc > 0) ) non_zero_indices_G[j].push_back(i);

      if(hc != -3 && ind_in_analysis(i) && (!strict_mode || (strict_mode && masked_indivs(i,0)) ) ) {
        G(j, i) -= total;
      } else G(j, i) = 0;

    }

    j++;snp_index_counter++;
  }

  sout << bs << " snps ";
}

void Data::skip_snps(int bs){

  std::string chromosome, rsid;
  uint32_t position;
  std::vector< std::string > alleles ;

  // skip the whole block of snps
  if(file_type == "bed") {
    bed_ifstream.seekg( bs * bed_block_size, ios_base::cur);
  } else {
    for(size_t snp = 0; snp < bs; snp++) {
      bgen.read_variant( &chromosome, &position, &rsid, &alleles );
      bgen.ignore_probs();
    }
  }

}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          level 0 models
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::ridge_level_0(int block) {

  sout << "   -calc level 0 ridge..." << flush;
  auto t2 = std::chrono::high_resolution_clock::now();

  int bs = GGt.rows();
  uint64 low = 0, high = 0;
  int block_eff = write_l0_pred ? 0 : block; // if writing to file
  string op_name, out_pheno;
  ofstream ofile;
  ident_l0 = MatrixXd::Identity(bs,bs);
  MatrixXd ww1, ww2, ww3, beta, pred, vmat, dvec, Xout;
  MatrixXd p_sum = MatrixXd::Zero(n_ridge_l0, n_pheno);
  MatrixXd p_sum2 = MatrixXd::Zero(n_ridge_l0, n_pheno);

  if(!within_sample_l0 && print_block_betas) {
    for(size_t ph = 0; ph < n_pheno; ++ph ){ 
      beta_print_out[ph].resize(n_ridge_l0, bs);
      beta_print_out[ph] = MatrixXd::Zero(n_ridge_l0, bs);
    }
  }

  uint64 i_total = 0, cum_size_folds = 0;
  for(size_t i = 0; i < cv_folds; ++i ) {
    // assign masking within folds
    for(size_t j = 0; j < cv_sizes[i]; ++j) {
      masked_in_folds[i].row(j) = masked_indivs.row(i_total);
      i_total++;
    }

    ww1 = GGt - G_folds[i];
    SelfAdjointEigenSolver<MatrixXd> eig(ww1);
    vmat = eig.eigenvectors();
    dvec = eig.eigenvalues();
    //if(i == 0)sout << ww1 << endl;
    if(i>0) low +=  cv_sizes[i-1];
    high += cv_sizes[i];
    for(size_t j = 0; j < n_ridge_l0; ++j ) {
      ww2 = vmat.transpose() * (GTY - GtY[i]);
      //if(i == 0)sout << ww2 << endl;
      //if(i == 0)sout << "lambda[j] =" << lambda[j] << endl;
      //if(i == 0)sout << ww1.array() + lambda[j]*ident_l0.array() << endl;
      ww3 = (dvec.array() + lambda[j]).inverse().matrix().asDiagonal() * ww2;
      //if(i == 0)sout << beta * (ww1+ lambda[j]*ident_l0)<< endl;
      beta = vmat * ww3;

      // save beta for each phenotype (only when using out-of-sample pred)
      if(!within_sample_l0 && print_block_betas)
        for(size_t ph = 0; ph < n_pheno; ++ph ) {
          beta_print_out[ph].row(j) += beta.col(ph).transpose();
        }

      pred = beta.transpose() * G;
      //if(i == 0)sout << pred.rows() << " " << pred.cols() << endl;
      //if(i == 0)sout << beta << endl;
      if(within_sample_l0) {
        // center and scale predictions
        VectorXd p_mean, p_sd;
        p_mean = pred.array().rowwise().mean();
        //if(i == 0)sout << i << " " << p_mean << endl;
        pred.colwise() -= p_mean;
        p_sd = pred.rowwise().norm() / sqrt(ind_in_analysis.sum() -1);
        //if(i == 0)sout << i << " " << p_sd << endl;
        pred.array().colwise() /= p_sd.array();
      } else {
        p_sum.row(j) += (pred.block(0, cum_size_folds, n_pheno, cv_sizes[i]).array() * masked_in_folds[i].transpose().array()).matrix().rowwise().sum();
        p_sum2.row(j) += (pred.block(0, cum_size_folds, n_pheno, cv_sizes[i]).array() * masked_in_folds[i].transpose().array()).matrix().rowwise().squaredNorm();
      }

      // store predictions
      uint64 kk = 0, jj = 0;
      for(size_t k = 0; k < n_samples; ++k ) {
        if( (k < low) | (k >= high) ) {
          if (within_sample_l0) {
            for(size_t ph = 0; ph < n_pheno; ++ph ) {
              pred_mat[ph][i](kk, block*n_ridge_l0 + j) = pred(ph, k);
              pred_pheno[ph][i](kk, 0) = phenotypes(k, ph);
              if (binary_mode && (block == 0) && (j == 0) ) {
                pred_pheno_raw[ph][i](kk, 0) = phenotypes_raw(k, ph);
                pred_offset[ph][i](kk, 0) = offset_logreg(k, ph);
              }
            }
          }
          kk+=1;
        } else {
          for(size_t ph = 0; ph < n_pheno; ++ph ) {
            test_mat[ph][i](jj, block_eff * n_ridge_l0 + j) = pred(ph, k);	      
            test_pheno[ph][i](jj, 0) = phenotypes(k, ph);
            if (binary_mode && (block == 0) && (j == 0) ) {
              test_pheno_raw[ph][i](jj, 0) = phenotypes_raw(k, ph);
              test_offset[ph][i](jj, 0) = offset_logreg(k, ph);
            }
          }
          jj+=1;
        }
      }
    }
    cum_size_folds += cv_sizes[i];
  }

  // when using only out-of-sample predictions for level 1 input features, center and scale using the whole sample
  if(!within_sample_l0){
    for(size_t ph = 0; ph < n_pheno; ++ph ) {
      RowVectorXd p_mean, p_invsd;
      p_mean = p_sum.col(ph).transpose() / Neff(ph);
      p_invsd = sqrt( (Neff(ph) - 1) / (p_sum2.col(ph).transpose().array() - Neff(ph) * p_mean.array().square()) );

      // scale printed estimates by the sd
      if(print_block_betas){ 
        beta_print_out[ph].array().colwise() *= p_invsd.transpose().array();
      }

      if(write_l0_pred) Xout = MatrixXd(n_samples, n_ridge_l0);

      cum_size_folds = 0;
      for(size_t i = 0; i < cv_folds; ++i ) {
        test_mat[ph][i].block(0, block_eff * n_ridge_l0, cv_sizes[i], n_ridge_l0).rowwise() -= p_mean;
        // mask missing
        test_mat[ph][i].block(0, block_eff * n_ridge_l0, cv_sizes[i], n_ridge_l0).array().colwise() *= masked_in_folds[i].col(ph).array();
        test_mat[ph][i].block(0, block_eff * n_ridge_l0, cv_sizes[i], n_ridge_l0).array().rowwise() *= p_invsd.array();

        if(write_l0_pred) {
          Xout.block(cum_size_folds, 0, cv_sizes[i], n_ridge_l0) = test_mat[ph][i].block(0, block_eff * n_ridge_l0, cv_sizes[i], n_ridge_l0);
          cum_size_folds += cv_sizes[i];
        }
      }

      // write predictions to file if specified
      if(write_l0_pred) {
        out_pheno = loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
        if(block == 0)
          ofile.open(out_pheno.c_str(), ios::out | ios::trunc | ios::binary );
        else
          ofile.open(out_pheno.c_str(), ios::out | ios::app | ios::binary );

        if (!ofile.is_open()) {
          sout << "ERROR : Cannote write temporary file " << out_pheno  << endl ;
          exit(-1);
        } 

        ofile.write( reinterpret_cast<char *> (&Xout(0,0)), Xout.rows() * Xout.cols() * sizeof(double) );
        ofile.close();
        //if(block < 2 && ph == 0 ) sout << endl << "Out " << endl <<  Xout.block(0, 0, 5, Xout.cols()) << endl;
      }

    }
  }

  // if printing betas to file (average over folds) [assume snp IDs are unique]
  //   -> separate file for each block (n_ridge_l0 rows & (2+bs) columns)
  if(!within_sample_l0 && print_block_betas) {
    op_name = out_file + "_block" + to_string(block+1) + ".betas";
    ofile.open(op_name.c_str());

    // Header: [TRAIT PARAM snpID1 ... snpIDk]
    ofile << "TRAIT PARAM " ;
    for(size_t i = 0; i < bs; ++i ) 
      ofile << snpinfo[print_snpcount++].ID << " ";
    ofile << endl;

    // Each line: [pheno# ridge# beta1 ... betak]
    for(size_t ph = 0; ph < n_pheno; ++ph ){ 
      beta_print_out[ph] /= cv_folds;
      for(size_t j = 0; j < n_ridge_l0; ++j ) {
        ofile << ph + 1 << " " <<  j + 1 << " ";
        for(size_t i = 0; i < bs; ++i ) 
          ofile << beta_print_out[ph](j,i) << " ";
        ofile << endl;
      }
    }
    ofile.close();
  }

  sout << "done";
  auto t3 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
  sout << " (" << duration.count() << "ms) "<< endl;

}

void Data::ridge_level_0_loocv(int block) {

  sout << "   -calc level 0 ridge..." << flush;
  auto t2 = std::chrono::high_resolution_clock::now();
  int bs = GGt.rows();
  int block_eff = write_l0_pred ? 0 : block; // if writing to file
  string out_pheno;
  ofstream ofile;
	VectorXd z1, gvec;
	MatrixXd VtG, z2, pred, Xout;
  RowVectorXd p_mean, p_sd;

  /*
  if(bs > n_samples){
    sout << "ERROR: Block size must be smaller than the number of samples to perform LOOCV!";
    exit(-1);
  }
  */


	// make matrix of (eigen-value + lambda)^(-1)
	Map<RowVectorXd> Lmap(lambda.data(), n_ridge_l0);
	MatrixXd dl = GGt_eig_val.asDiagonal() * MatrixXd::Ones(bs, n_ridge_l0);
	dl.rowwise() += Lmap;
	MatrixXd DL_inv = dl.array().inverse().matrix();

  uint64 max_bytes = chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (bs * target_size) matrix ]
  int nchunk = ceil( cv_folds * bs * sizeof(double) * 1.0 / max_bytes );
  if (verbose) sout << nchunk << " chunks..." << flush;
  int chunk, size_chunk, target_size = cv_folds / nchunk;
  int j_start;

  for(chunk = 0; chunk < nchunk; ++chunk ) {
    size_chunk = chunk == nchunk - 1? cv_folds - target_size * chunk : target_size;
    j_start = chunk * target_size;

    VtG = GGt_eig_vec.transpose() * G.block(0, j_start, bs, size_chunk);
    for(size_t i = 0; i < size_chunk; ++i ) {
      z1 = VtG.col(i);
      z2 = DL_inv.array().colwise() * z1.array();
      gvec = z2.transpose() * z1;
      pred = z2.transpose() * Wmat - gvec * phenotypes.row(j_start + i);
      pred.array().colwise() /= 1 - gvec.array();
      for(size_t ph = 0; ph < n_pheno; ++ph )
        test_mat_conc[ph].block(j_start + i, block_eff * n_ridge_l0, 1, n_ridge_l0) = pred.col(ph).transpose();
    }
  }

  // center and scale within the block
  for(size_t ph = 0; ph < n_pheno; ++ph ) { 
    // mask missing first
    test_mat_conc[ph].block(0, block_eff * n_ridge_l0, n_samples, n_ridge_l0).array().colwise() *= masked_indivs.col(ph).array();
    p_mean = test_mat_conc[ph].block(0, block_eff * n_ridge_l0, n_samples, n_ridge_l0).colwise().sum() / Neff(ph);
    //if(i == 0)sout << i << " " << p_mean << endl;
    test_mat_conc[ph].block(0, block_eff * n_ridge_l0, n_samples, n_ridge_l0).rowwise() -= p_mean;
    // mask missing again 
    test_mat_conc[ph].block(0, block_eff * n_ridge_l0, n_samples, n_ridge_l0).array().colwise() *= masked_indivs.col(ph).array();
    p_sd = test_mat_conc[ph].block(0, block_eff * n_ridge_l0, n_samples, n_ridge_l0).colwise().norm() / sqrt(Neff(ph) -1);
    //if(i == 0)sout << i << " " << p_sd << endl;
    test_mat_conc[ph].block(0, block_eff * n_ridge_l0, n_samples, n_ridge_l0).array().rowwise() /= p_sd.array();


    if(write_l0_pred) {
      Xout = test_mat_conc[ph].block(0, 0, n_samples, n_ridge_l0);
      out_pheno = loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
      if(block == 0)
        ofile.open(out_pheno.c_str(), ios::out | ios::trunc | ios::binary );
      else
        ofile.open(out_pheno.c_str(), ios::out | ios::app | ios::binary );

      if (!ofile.is_open()) {
        sout << "ERROR : Cannote write temporary file " << out_pheno  << endl ;
        exit(-1);
      } 

      ofile.write( reinterpret_cast<char *> (&Xout(0,0)), Xout.rows() * Xout.cols() * sizeof(double) );
      ofile.close();
      //if(block < 2 && ph == 0 ) sout << endl << "Out " << endl <<  Xout.block(0, 0, 5, Xout.cols()) << endl;
    }

  }

  sout << "done";
  auto t3 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
  sout << " (" << duration.count() << "ms) "<< endl;

}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          level 1 models
/////////////////////////////////////////////////
/////////////////////////////////////////////////
//
void Data::ridge_level_1() {

  sout << endl << " Level 1 ridge..." << endl << flush;

  int bs_l1 = total_n_block * n_ridge_l0;
  int ph_eff;
  string in_pheno;
  ifstream infile;
  MatrixXd X1, X2, beta_l1, p1, vmat, dvec, dl_inv;
  VectorXd VtX2;
  MatrixXd XtX_sum, XtY_sum;
  ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);
  Map<RowVectorXd> Lmap(tau.data(), n_ridge_l1);

  // to compute Rsq and MSE of predictions
  for (int i = 0; i < 5; i++)
    cumsum_values[i].setZero(n_pheno, n_ridge_l1);

  for(size_t ph = 0; ph < n_pheno; ++ph ) {
    sout << "   -on phenotype " << ph+1 <<" (" << pheno_names[ph] << ")..." << flush;
    auto ts1 = std::chrono::high_resolution_clock::now();
    ph_eff = write_l0_pred ? 0 : ph;

    // read in level 0 predictions from file
    if(write_l0_pred){

      // allocate memory 
      if(ph == 0) {
        for( std::size_t i = 0; i < cv_folds; ++i ) 
          test_mat[ph_eff][i].resize(cv_sizes[i], bs_l1);
      }

      in_pheno = loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
      infile.open(in_pheno.c_str(), ios::in | ios::binary );

      if (!infile.is_open()) {
        sout << "ERROR : Cannote read temporary file " << in_pheno  << endl ;
        exit(-1);
      } 

      // store back values in test_mat
      for( std::size_t m = 0; m < bs_l1; ++m ) 
        for( std::size_t i = 0; i < cv_folds; ++i ) 
          for( std::size_t k = 0; k < cv_sizes[i]; ++k ) 
            infile.read( reinterpret_cast<char *> (&test_mat[ph_eff][i](k,m)), sizeof(double) );

      infile.close();
      //if(ph == 0) sout << endl << "In:\n" << test_mat[ph_eff][0].block(0,0,5,6) << endl;
    }

    // compute XtX and Xty for each fold and cum. sum using test_mat's
    if (!within_sample_l0){
      XtX_sum.setZero(bs_l1, bs_l1);
      XtY_sum.setZero(bs_l1, 1);
      for( std::size_t i = 0; i < cv_folds; ++i ) {
        X_folds[i] = test_mat[ph_eff][i].transpose() * test_mat[ph_eff][i];
        XtY[i]     = test_mat[ph_eff][i].transpose() * test_pheno[ph][i];
        XtX_sum += X_folds[i];
        XtY_sum += XtY[i];
      }
    }
    for(size_t i = 0; i < cv_folds; ++i ) {

      // use either in-sample or out-of-sample predictions
      if (within_sample_l0) {
        X1 = pred_mat[ph][i].transpose() * pred_mat[ph][i];
        X2 = pred_mat[ph][i].transpose() * pred_pheno[ph][i];
      } else{
        X1 = XtX_sum - X_folds[i];
        X2 = XtY_sum - XtY[i];
      }

      SelfAdjointEigenSolver<MatrixXd> eigX1(X1);
      vmat = eigX1.eigenvectors();
      dvec = eigX1.eigenvalues();
      VtX2 = vmat.transpose() * X2;
      // compute solutions for all ridge parameters at once
      // p1 is Nfold x nridge_l1 matrix
      dl_inv = ( (dvec.asDiagonal() *  MatrixXd::Ones(bs_l1, n_ridge_l1)).rowwise() + Lmap).array().inverse().matrix();
      dl_inv.array().colwise() *= VtX2.array();
      beta_l1 = vmat * dl_inv;
      if(!within_sample_l0) beta_hat_level_1[ph][i] = beta_l1;
      p1 = test_mat[ph_eff][i] * beta_l1;

      cumsum_values[0].row(ph) += p1.colwise().sum();
      cumsum_values[1].row(ph).array() += test_pheno[ph][i].array().sum();
      cumsum_values[2].row(ph) += p1.array().square().matrix().colwise().sum();
      cumsum_values[3].row(ph).array() += test_pheno[ph][i].array().square().sum();
      cumsum_values[4].row(ph) += (p1.array().colwise() * test_pheno[ph][i].col(0).array()).matrix().colwise().sum() ;
    }

    sout << "done";
    auto ts2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(ts2 - ts1);
    sout << " (" << duration.count() << "ms) "<< endl;
  }

  sout << endl;
}

void Data::ridge_level_1_loocv() {

  sout << endl << " Level 1 ridge..." << flush;
	    
  int bs_l1 = total_n_block * n_ridge_l0;
  int ph_eff;
  string in_pheno;
  ifstream infile;
  MatrixXd Xmat_chunk, Yvec_chunk, Z1, Z2, dl, dl_inv, xtx;
  VectorXd wvec, zvec;
  RowVectorXd calFactor, pred;

  for (int i = 0; i < 5; i++)
    cumsum_values[i].setZero(n_pheno, n_ridge_l1);

	// make matrix of (eigen-values + tau)^(-1)
	Map<RowVectorXd> Lmap(tau.data(), n_ridge_l1);

  uint64 max_bytes = chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (target_size * bs_l1) matrix ]
  int nchunk = ceil( cv_folds * bs_l1 * sizeof(double) * 1.0 / max_bytes );
  if (verbose) sout << nchunk << " chunks...";
  sout << endl;
  int chunk, size_chunk, target_size = cv_folds / nchunk;
  int j_start;

  for(size_t ph = 0; ph < n_pheno; ++ph ) {
    sout << "   -on phenotype " << ph+1 <<" (" << pheno_names[ph] <<")..." << flush;
    auto ts1 = std::chrono::high_resolution_clock::now();
    ph_eff = write_l0_pred ? 0 : ph;

    // read in level 0 predictions from file
    if(write_l0_pred){

      // allocate memory (re-use same matrix for all traits) 
      if(ph == 0) test_mat_conc[ph_eff].resize(n_samples, bs_l1);

      in_pheno = loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
      infile.open(in_pheno.c_str(), ios::in | ios::binary );

      if (!infile.is_open()) {
        sout << "ERROR : Cannote read temporary file " << in_pheno  << endl ;
        exit(-1);
      } 

      // store back values in test_mat_conc
      infile.read( reinterpret_cast<char *> (&test_mat_conc[ph_eff](0,0)), n_samples * bs_l1 * sizeof(double) );

      infile.close();
      //if(ph == 0) sout << endl << "In:\n" << test_mat_conc[ph_eff].block(0,0,5,6) << endl;
    }

    xtx = test_mat_conc[ph_eff].transpose() * test_mat_conc[ph_eff];
    SelfAdjointEigenSolver<MatrixXd> eigX(xtx);
		dl = eigX.eigenvalues().asDiagonal() * MatrixXd::Ones(bs_l1, n_ridge_l1);
		dl.rowwise() += Lmap;
		dl_inv = dl.array().inverse().matrix();

    zvec = test_mat_conc[ph_eff].transpose() * phenotypes.col(ph);
    wvec = eigX.eigenvectors().transpose() * zvec;

    for(chunk = 0; chunk < nchunk; ++chunk ) {
      size_chunk = chunk == nchunk - 1? cv_folds - target_size * chunk : target_size;
      j_start = chunk * target_size;
      Xmat_chunk = test_mat_conc[ph_eff].block(j_start, 0, size_chunk, bs_l1);
      Yvec_chunk = phenotypes.block(j_start, ph, size_chunk, 1);
      Z1 = (Xmat_chunk * eigX.eigenvectors()).transpose();

      for(size_t i = 0; i < size_chunk; ++i ) {
        Z2 = (dl_inv.array().colwise() * Z1.col(i).array()).matrix();
        calFactor = Z1.col(i).transpose() * Z2;
				pred = wvec.transpose() * Z2;
        pred -=  Yvec_chunk(i, 0) * calFactor;
				pred.array()  /= 1 - calFactor.array();
        //if( ph == 0) sout << pred.head(5) << endl;

        // compute mse and rsq
        cumsum_values[0].row(ph) += pred; // Sx
        // Y is centered so Sy = 0
        cumsum_values[2].row(ph) += pred.array().square().matrix(); // Sx2
        // Y is scaled so Sy2 = n_samples - 1
        cumsum_values[4].row(ph).array() += pred.array() * Yvec_chunk(i,0); // Sxy
      }
    }
    sout << "done";
    auto ts2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(ts2 - ts1);
    sout << " (" << duration.count() << "ms) "<< endl;
  }
	cumsum_values[3].array().colwise() += Neff - 1; // Sy2

  sout << endl;
}

void Data::ridge_logistic_level_1() {

  sout << endl << " Level 1 ridge with logistic regression..." << endl << flush;

  int niter_cur;
  int bs_l1 = total_n_block * n_ridge_l0;
  int ph_eff;
  string in_pheno;
  ifstream infile;
  ArrayXd Y1, W1, p1, score;
  ArrayXd betaold, etavec, pivec, wvec, zvec, betanew, etatest;
  MatrixXd X1, XtW, XtWX, XtWZ;
  pheno_l1_not_converged = ArrayXd::Zero(n_pheno);

  ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);
  for (int i = 0; i < 6; i++)
    cumsum_values[i].setZero(n_pheno, n_ridge_l1);

  for(size_t ph = 0; ph < n_pheno; ++ph ) {
    sout << "   -on phenotype " << ph+1 <<" (" << pheno_names[ph] <<")..." << flush;
    auto ts1 = std::chrono::high_resolution_clock::now();
    ph_eff = write_l0_pred ? 0 : ph;

    // read in level 0 predictions from file
    if(write_l0_pred){

      // allocate memory 
      if(ph == 0) {
        for( std::size_t i = 0; i < cv_folds; ++i ) 
          test_mat[ph_eff][i].resize(cv_sizes[i], bs_l1);
      }

      in_pheno = loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
      infile.open(in_pheno.c_str(), ios::in | ios::binary );

      if (!infile.is_open()) {
        sout << "ERROR : Cannote read temporary file " << in_pheno  << endl ;
        exit(-1);
      } 

      // store back values in test_mat
      for( std::size_t m = 0; m < bs_l1; ++m ) 
        for( std::size_t i = 0; i < cv_folds; ++i ) 
          for( std::size_t k = 0; k < cv_sizes[i]; ++k ) 
            infile.read( reinterpret_cast<char *> (&test_mat[ph_eff][i](k,m)), sizeof(double) );

      infile.close();
      //if(ph == 0) sout << endl << "In:\n" << test_mat[ph_eff][0].block(0,0,5,6) << endl;
    }

    for(size_t i = 0; i < cv_folds; ++i ) {
      if( pheno_l1_not_converged(ph) ) break;

      if( within_sample_l0 ){
        X1 = pred_mat[ph][i];
        Y1 = pred_pheno_raw[ph][i];
        W1 = pred_offset[ph][i]; 
      }

      for(size_t j = 0; j < n_ridge_l1; ++j ) {
        if( pheno_l1_not_converged(ph) ) break;

        // starting values
        betaold = ArrayXd::Zero(bs_l1);

        niter_cur = 0;
        while(niter_cur++ < niter_max){

          if(within_sample_l0) {
            etavec = W1 + (X1 * betaold.matrix()).array();
            pivec = 1 - 1/(etavec.exp() + 1);
            wvec = pivec * (1 - pivec);
            // check none of the values are 0
            if( ( wvec == 0 ).count() > 0 ){
              sout << "ERROR: Zeros occured in Var(Y) during ridge logistic regression! (Try with --loocv)" << endl;
              pheno_l1_not_converged(ph) = 1;
              break;
            }
            zvec = (etavec - W1) + (Y1 - pivec) / wvec;
            XtW = X1.transpose() * wvec.matrix().asDiagonal();
            betanew = (XtW * X1 + tau[j] * ident_l1).colPivHouseholderQr().solve(XtW * zvec.matrix()).array();
            // get the score
            etavec = W1 + (X1 * betanew.matrix()).array();
            pivec = 1 - 1/(etavec.exp() + 1);
            score = (X1.transpose() * (Y1 - pivec).matrix()).array() - tau[j] * betanew; 

          } else {
            
            XtWX = MatrixXd::Zero(bs_l1, bs_l1);
            XtWZ = MatrixXd::Zero(bs_l1, 1);

            for(size_t k = 0; k < cv_folds; ++k ) {
              if( k != i) {
                etavec = (masked_in_folds[k].col(ph).array() == 1).select( (test_offset[ph][k] + test_mat[ph_eff][k] * betaold.matrix()).array() , 0);
                pivec = 1 - 1/(etavec.exp() + 1);
                wvec = (masked_in_folds[k].col(ph).array() == 1).select(pivec * (1 - pivec), 0);
                // check none of the values are 0
                if( ( (masked_in_folds[k].col(ph).array() == 1) &&  (wvec == 0) ).count() > 0 ){
                  sout << "ERROR: Zeros occured in Var(Y) during ridge logistic regression! (Try with --loocv)" << endl;
                  pheno_l1_not_converged(ph) = 1;
                  break;
                }
                zvec = (masked_in_folds[k].col(ph).array() == 1).select((etavec - test_offset[ph][k].array()) + (test_pheno_raw[ph][k].array() - pivec) / wvec, 0);

                XtW = test_mat[ph_eff][k].transpose() * wvec.matrix().asDiagonal();
                XtWX += XtW * test_mat[ph_eff][k];
                XtWZ += XtW * zvec.matrix();
              }
            }
            if( pheno_l1_not_converged(ph) ) break;

            betanew = ((XtWX + tau[j] * ident_l1).llt().solve(XtWZ)).array();

            // get the score
            score = ArrayXd::Zero(betanew.size());
            for(size_t k = 0; k < cv_folds; ++k ) {
              if( k != i) {
                etavec = (masked_in_folds[k].col(ph).array() == 1).select( (test_offset[ph][k] + test_mat[ph_eff][k] * betanew.matrix()).array() , 0);
                pivec = 1 - 1/(etavec.exp() + 1);

                score += (test_mat[ph_eff][k].transpose() * (masked_in_folds[k].col(ph).array() == 1).select(test_pheno_raw[ph][k].array() - pivec, 0).matrix()).array();  
              }
            }
            score -= tau[j] * betanew;

          }

          // stopping criterion
          if( score.abs().maxCoeff() < numtol) break;

          betaold = betanew;
        }

        if(niter_cur > niter_max){
          sout << "WARNING: Penalized logistic regression did not converge! (Increase --niter)\n";
          pheno_l1_not_converged(ph) = 1;
          break;
        }
        //sout << "Converged in "<< niter_cur << " iterations." << endl;
        //sout << score.abs().maxCoeff() << endl;

        etatest = test_offset[ph][i].array() + (test_mat[ph_eff][i] * betanew.matrix()).array();
        p1 = (1 - 1/(etatest.exp() + 1));

        if(!within_sample_l0) beta_hat_level_1[ph][i].col(j) = betanew;
          
        // compute mse
        for(size_t l = 0; l < cv_sizes[i]; l++){
          if(!masked_in_folds[i](l,ph)) continue;
          cumsum_values[0](ph,j) += p1(l); // Sx
          cumsum_values[1](ph,j) += test_pheno_raw[ph][i](l,0); // Sy
          cumsum_values[2](ph,j) += p1(l) * p1(l); // Sx2
          cumsum_values[3](ph,j) += test_pheno_raw[ph][i](l,0) * test_pheno_raw[ph][i](l,0); // Sy2
          cumsum_values[4](ph,j) += p1(l) * test_pheno_raw[ph][i](l,0); // Sxy
          cumsum_values[5](ph,j) += compute_log_lik(test_pheno_raw[ph][i](l,0), p1(l)); // Sxy
        }
      }
    }

    sout << "done";
    auto ts2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(ts2 - ts1);
    sout << " (" << duration.count() << "ms) "<< endl;
  }

  sout << endl;
}

void Data::ridge_logistic_level_1_loocv() {

  sout << endl << " Level 1 ridge with logistic regression..." << flush;

  int niter_cur;
  int bs_l1 = total_n_block * n_ridge_l0;
  int ph_eff;
	double v2, pred, p1;
  string in_pheno;
  ifstream infile;
  ArrayXd betaold, etavec, pivec, wvec, zvec, betanew, score;
  MatrixXd XtWX, XtWZ;
  MatrixXd V1, Xmat_chunk, b_loo;
  VectorXd Yvec_chunk, mask_chunk;
  LLT<MatrixXd> Hinv;
  pheno_l1_not_converged = ArrayXd::Zero(n_pheno);

  ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);
  for (int i = 0; i < 6; i++)
    cumsum_values[i].setZero(n_pheno, n_ridge_l1);

  uint64 max_bytes = chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (bs_l1 * target_size) matrix ]
  int nchunk = ceil( cv_folds * bs_l1 * sizeof(double) * 1.0 / max_bytes );
  if (verbose) 
    sout << nchunk << " chunks..." << endl;
  else 
    sout << endl;
  int chunk, size_chunk, target_size = cv_folds / nchunk;
  int j_start;

  for(size_t ph = 0; ph < n_pheno; ++ph ) {
    sout << "   -on phenotype " << ph+1 <<" (" << pheno_names[ph] <<")..." << flush;
    auto ts1 = std::chrono::high_resolution_clock::now();
    ph_eff = write_l0_pred ? 0 : ph;

    // read in level 0 predictions from file
    if(write_l0_pred){

      // allocate memory (re-use same matrix for all traits) 
      if(ph == 0) test_mat_conc[ph_eff].resize(n_samples, bs_l1);

      in_pheno = loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
      infile.open(in_pheno.c_str(), ios::in | ios::binary );

      if (!infile.is_open()) {
        sout << "ERROR : Cannote read temporary file " << in_pheno  << endl ;
        exit(-1);
      } 

      // store back values in test_mat_conc
      infile.read( reinterpret_cast<char *> (&test_mat_conc[ph_eff](0,0)), n_samples * bs_l1 * sizeof(double) );

      infile.close();
      //if(ph == 0) sout << endl << "In:\n" << test_mat_conc[ph_eff].block(0,0,5,6) << endl;
    }

    for(size_t j = 0; j < n_ridge_l1; ++j ) {
      // starting values
      betaold = ArrayXd::Zero(bs_l1);

      niter_cur = 0;
      while(niter_cur++ < niter_max){

        etavec = (masked_indivs.col(ph).array() == 1).select( (offset_logreg.col(ph) + test_mat_conc[ph_eff] * betaold.matrix()).array(), 0);
        pivec = 1 - 1/(etavec.exp() + 1);
        wvec = (masked_indivs.col(ph).array() == 1).select( pivec * (1 - pivec), 0);
        // check none of the values are 0
        if( ( (masked_indivs.col(ph).array() == 1) &&  (wvec == 0) ).count() > 0 ){
          sout << "ERROR: Zeros occured in Var(Y) during ridge logistic regression! (Try with more common SNPs)" << endl;
          pheno_l1_not_converged(ph) = 1;
          break;
        }
        zvec = (masked_indivs.col(ph).array() == 1).select( (etavec - offset_logreg.col(ph).array()) + (phenotypes_raw.col(ph).array() - pivec) / wvec, 0);
        V1 = test_mat_conc[ph_eff].transpose() * wvec.matrix().asDiagonal();
        XtWX = V1 * test_mat_conc[ph_eff];
        XtWZ = V1 * zvec.matrix();
        Hinv.compute( XtWX + tau[j] * ident_l1 );

        betanew = (Hinv.solve(XtWZ)).array();
        // get the score
        etavec = (masked_indivs.col(ph).array() == 1).select( (offset_logreg.col(ph) + test_mat_conc[ph_eff] * betanew.matrix()).array(), 0);
        pivec = 1 - 1/(etavec.exp() + 1);
        score = ( test_mat_conc[ph_eff].transpose() * (masked_indivs.col(ph).array() == 1).select(phenotypes_raw.col(ph).array() - pivec, 0).matrix()).array() ;
        score -= tau[j] * betanew; 

        if( score.abs().maxCoeff() < numtol) break;

        betaold = betanew;

      }

      if(niter_cur > niter_max){
        sout << "WARNING: Ridge logistic regression did not converge! (Increase --niter)\n";
        pheno_l1_not_converged(ph) = 1;
      }
      if( pheno_l1_not_converged(ph) ) break;

      //sout << "Converged in "<< niter_cur << " iterations." << endl;
      //sout << score.abs().maxCoeff() << endl;

      // compute Hinv 
      etavec = (offset_logreg.col(ph) + test_mat_conc[ph_eff] * betanew.matrix()).array();
      pivec = 1 - 1/(etavec.exp() + 1);
      wvec = (masked_indivs.col(ph).array() == 1).select( pivec * (1 - pivec), 0 );
      zvec = (masked_indivs.col(ph).array() == 1).select( (etavec - offset_logreg.col(ph).array()) + (phenotypes_raw.col(ph).array() - pivec) / wvec, 0);
      V1 = test_mat_conc[ph_eff].transpose() * wvec.matrix().asDiagonal();
      XtWX = V1 * test_mat_conc[ph_eff];
      Hinv.compute( XtWX + tau[j] * ident_l1 );

      // LOOCV estimates
      for(chunk = 0; chunk < nchunk; ++chunk ) {
        size_chunk = chunk == nchunk - 1? cv_folds - target_size * chunk : target_size;
        j_start = chunk * target_size;

        Xmat_chunk = test_mat_conc[ph_eff].block(j_start, 0, size_chunk, bs_l1); // n x k
				Yvec_chunk = phenotypes_raw.block(j_start, ph, size_chunk, 1);
        mask_chunk = masked_indivs.block(j_start, ph, size_chunk,1);

        V1 = Hinv.solve( Xmat_chunk.transpose() ); // k x n
        for(size_t i = 0; i < size_chunk; ++i ) {
          if(!mask_chunk(i)) continue;
          v2 = Xmat_chunk.row(i) * V1.col(i); 
          v2 *= wvec(j_start + i); 
          b_loo = (betanew - V1.col(i).array() * (Yvec_chunk(i) - pivec(j_start + i)) / (1 - v2)).matrix();
					pred = Xmat_chunk.row(i) * b_loo.col(0); 
					pred += offset_logreg(j_start + i, ph);
					p1 = 1 - 1/ ( exp(pred) + 1 );

					// compute mse and rsq
					cumsum_values[0](ph,j) += p1; // Sx
					cumsum_values[1](ph,j) += Yvec_chunk(i); // Sy
					cumsum_values[2](ph,j) += p1 * p1; // Sx2
					cumsum_values[3](ph,j) += Yvec_chunk(i) * Yvec_chunk(i); // Sy2
					cumsum_values[4](ph,j) += p1 * Yvec_chunk(i); // Sxy
					cumsum_values[5](ph,j) += compute_log_lik(Yvec_chunk(i), p1); // Sxy
        }
      }
    }

    sout << "done";
    auto ts2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(ts2 - ts1);
    sout << " (" << duration.count() << "ms) "<< endl;
  }

  sout << endl;
}

double Data::compute_log_lik(double y, double p){
  // negative log likelihood for bernoulli
  double ll;
  ll = - y * log(p) - (1 - y) * log(1-p);
  return(ll);
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          Evaluate level 1 output
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::output() {
  
  int min_index;
  double performance_measure, rsq, sse, ll_avg, min_val;
  string pfile, out_blup_list, pline;
  ofstream outb;
  G.resize(0,0);// free genotype block

  sout << "Output" << endl << "------" << endl;

  if(make_loco || binary_mode){
    out_blup_list = out_file + "_pred.list";
    outb.open(out_blup_list.c_str());
  }

  for(size_t ph = 0; ph < n_pheno; ++ph ) { 
    sout << "phenotype " << ph+1 << " (" << pheno_names[ph] << ") : " ;
    if( make_loco || binary_mode ) {
      // for quantitative traits
      if( !binary_mode ) {
        outb << pheno_names[ph]  << " " << out_file << "_" << ph + 1 << ".loco" << endl;
      } else {
        // for binary traits - check level 1 ridge converged
        if( pheno_l1_not_converged(ph)==0 ) {
          outb << pheno_names[ph]  << " " << out_file << "_" << ph + 1 << ".loco" << endl;
        } else {
          if(write_l0_pred){ // cleanup level 0 predictions
            pfile = loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
            remove(pfile.c_str());
          }
          sout << "Level 1 logistic did not converge. LOCO predictions calculations are skipped.\n\n";
          continue;
        }
      }
    }
    sout << endl;

    min_index = 0;
    min_val = 1e10;

    // determine optimal parameter by cv using: QT: MSE, BT: -loglik
    for(int j = 0; j < n_ridge_l1; ++j ) {
      if(!binary_mode)
        performance_measure = cumsum_values[2](ph, j) + cumsum_values[3](ph,j) - 2 * cumsum_values[4](ph,j);
      else
        performance_measure = cumsum_values[5](ph, j);
      performance_measure /= Neff(ph);

      if( performance_measure < min_val) {
        min_index = j;
        min_val = performance_measure;
      }
    }

    for(int j = 0; j < n_ridge_l1; ++j ) {
      if(!binary_mode)
        sout << "  " << setw(5) << (total_n_block *  n_ridge_l0)  / tau[j] ;
      else 
        sout << "  " << setw(5) << (total_n_block *  n_ridge_l0) / ( (total_n_block *  n_ridge_l0) + (M_PI * M_PI) * tau[j] / 3 );

      // output Rsq and MSE
      rsq = cumsum_values[4](ph,j) - cumsum_values[0](ph,j) * cumsum_values[1](ph,j) / Neff(ph); // num = Sxy - SxSy/n
      rsq = (rsq * rsq) / ((cumsum_values[2](ph,j) - cumsum_values[0](ph,j) * cumsum_values[0](ph,j) / Neff(ph)) * (cumsum_values[3](ph,j) - cumsum_values[1](ph,j) * cumsum_values[1](ph,j) / Neff(ph))); // num^2 / ( (Sx2 - Sx^2/n)* (Sy2 - Sy^2/n) )
      sse = cumsum_values[2](ph, j) + cumsum_values[3](ph,j) - 2 * cumsum_values[4](ph,j); // Sx2 + Sy2 - SxSy
      if(binary_mode) ll_avg = cumsum_values[5](ph, j) / Neff(ph);

      sout << " : ";
      sout << "Rsq = " << rsq << ", MSE = " << sse/Neff(ph);
      if(binary_mode) sout << ", -logLik/N = " << ll_avg;
      if(j == min_index) sout << "<- min value";
      sout << endl;
    }
    
    if(!binary_mode){
      if(use_loocv) make_predictions_loocv(ph, min_index);
      else make_predictions(ph, min_index);
    } else if(use_loocv) make_predictions_binary_loocv(ph, min_index);
    else make_predictions_binary(ph, min_index);


    // delete file used to store l0 predictions
    if(write_l0_pred){
      pfile = loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
      remove(pfile.c_str());
    }

  }

  if(make_loco || binary_mode){
    outb.close();
    sout << "List of blup files written to: [" << out_blup_list << "] " << endl;
  }

}

void Data::make_predictions(int ph, int val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  int bs_l1 = total_n_block * n_ridge_l0;
  int ph_eff = write_l0_pred ? 0 : ph;
  string outname, in_pheno;
  ifstream infile;
  ofstream ofile;
  MatrixXd X1, X2, beta_l1, beta_avg;
  ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);

  // read in level 0 predictions from file
  if(write_l0_pred){

    in_pheno = loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
    infile.open(in_pheno.c_str(), ios::in | ios::binary );

    if (!infile.is_open()) {
      sout << "ERROR : Cannote read temporary file " << in_pheno  << endl ;
      exit(-1);
    } 

    // store back values in test_mat
    for( std::size_t m = 0; m < bs_l1; ++m ) 
      for( std::size_t i = 0; i < cv_folds; ++i ) 
        for( std::size_t k = 0; k < cv_sizes[i]; ++k ) 
          infile.read( reinterpret_cast<char *> (&test_mat[ph_eff][i](k,m)), sizeof(double) );

    infile.close();
  }


  if(within_sample_l0){
    X1 = test_mat[ph_eff][0].transpose() * test_mat[ph_eff][0];
    X2 = test_mat[ph_eff][0].transpose() * test_pheno[ph][0];
    for(size_t i = 1; i < cv_folds; ++i ) {
      X1 += test_mat[ph_eff][i].transpose() * test_mat[ph_eff][i];
      X2 += test_mat[ph_eff][i].transpose() * test_pheno[ph][i];
    }
    beta_l1 = (X1 + tau[val] * ident_l1).llt().solve(X2);
  } else if(print_block_betas) {
    beta_avg = MatrixXd::Zero(bs_l1, 1);
    for(size_t i = 0; i < cv_folds; ++i ) {
      beta_avg += beta_hat_level_1[ph][i].col(val);
    }
    beta_avg /= cv_folds;
  }

  // if specified, write betas to file (open in append mode)
  if(!within_sample_l0 && print_block_betas) {
    outname = out_file + "_level1.betas";
    ofile.open(outname.c_str(), ios::out | ios::app);
    ofile << ph + 1 << " ";
    ofile << beta_avg.transpose() << endl;
    ofile.close();
  }

  // sout << "\nFor tau[" << val <<"] = " << tau[val] << endl <<  beta_l1 << endl ;
  int ctr = 0, chr_ctr = 0; 
  int nn, cum_size_folds;
  map<int, vector<int> >::iterator itr;
  for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) {
    nn = itr->second[1] * n_ridge_l0;
    if(nn > 0) { 
      cum_size_folds = 0;
      for(size_t i = 0; i < cv_folds; ++i ) {
        if(!within_sample_l0) beta_l1 = beta_hat_level_1[ph][i].col(val);
        predictions[0].block(cum_size_folds, chr_ctr, cv_sizes[i], 1) = test_mat[ph_eff][i].block(0, ctr, cv_sizes[i], nn) * beta_l1.block(ctr, 0, nn, 1);
        cum_size_folds += cv_sizes[i];
      }
      chr_ctr++;
      ctr += nn;
    }
  }

  write_predictions(ph);

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl << endl;
}

void Data::make_predictions_loocv(int ph, int val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();
	    
  int bs_l1 = total_n_block * n_ridge_l0;
  int ph_eff = write_l0_pred ? 0 : ph;
  string in_pheno;
  ifstream infile;
  MatrixXd Xmat_chunk, Yvec_chunk,  Z1, Z2, b0, xtx;
  VectorXd w1, w2, Vw2, zvec;
  RowVectorXd calFactor;
	ArrayXd dl_inv;

  uint64 max_bytes = chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (target_size * bs_l1) matrix ]
  int nchunk = ceil( cv_folds * bs_l1 * sizeof(double) * 1.0 / max_bytes );
  if (verbose) sout << nchunk << " chunks..." << flush;
  int chunk, size_chunk, target_size = cv_folds / nchunk;
  int j_start;


  // read in level 0 predictions from file
  if(write_l0_pred){

    in_pheno = loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
    infile.open(in_pheno.c_str(), ios::in | ios::binary );

    if (!infile.is_open()) {
      sout << "ERROR : Cannote read temporary file " << in_pheno  << endl ;
      exit(-1);
    } 

    // store back values in test_mat_conc
    infile.read( reinterpret_cast<char *> (&test_mat_conc[ph_eff](0,0)), n_samples * bs_l1 * sizeof(double) );

    infile.close();
  }

  // fit model on whole data again for optimal ridge param
  xtx = test_mat_conc[ph_eff].transpose() * test_mat_conc[ph_eff];
  SelfAdjointEigenSolver<MatrixXd> eigX(xtx);
	zvec = test_mat_conc[ph_eff].transpose() * phenotypes.col(ph);
	w1 = eigX.eigenvectors().transpose() * zvec;
	dl_inv = (eigX.eigenvalues().array() + tau[val]).inverse();
  w2 = (w1.array() * dl_inv).matrix();
  Vw2 = eigX.eigenvectors() * w2;

  for(chunk = 0; chunk < nchunk; ++chunk ) {
    size_chunk = chunk == nchunk - 1? cv_folds - target_size * chunk : target_size;
    j_start = chunk * target_size;
    Xmat_chunk = test_mat_conc[ph_eff].block(j_start, 0, size_chunk, bs_l1);
    Yvec_chunk = phenotypes.block(j_start, ph, size_chunk, 1);

    Z1 = (Xmat_chunk * eigX.eigenvectors()).transpose();
    Z2 = dl_inv.matrix().asDiagonal() * Z1;
    calFactor = (Z1.array() * Z2.array()).matrix().colwise().sum();
    b0 = eigX.eigenvectors() * Z2;
    b0.array().rowwise() *= (w2.transpose() * Z1 - Yvec_chunk.transpose()).array() / (1 - calFactor.array());
    b0.colwise() += Vw2;

    int ctr = 0, chr_ctr = 0;
    int nn;
    map<int, vector<int> >::iterator itr;
    for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) {
      nn = itr->second[1] * n_ridge_l0;
      if(nn > 0) { 
        predictions[0].block(j_start, chr_ctr, size_chunk, 1) = (test_mat_conc[ph_eff].block(j_start, ctr, size_chunk, nn).array() * b0.block(ctr, 0, nn, size_chunk).transpose().array()).rowwise().sum();
        chr_ctr++;
        ctr += nn;
      }
    }
  }

  write_predictions(ph);
    
  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl << endl;
}

// predictions for binary traits
void Data::make_predictions_binary(int ph, int val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  int bs_l1 = total_n_block * n_ridge_l0;
  int ph_eff = write_l0_pred ? 0 : ph;
  string in_pheno;
  ifstream infile;
  ArrayXd etavec, pivec, wvec, zvec, score;
  MatrixXd betaold, betanew, XtW, XtWX, XtWZ, Xconc;
  ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);

  // read in level 0 predictions from file
  if(write_l0_pred){

    in_pheno = loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
    infile.open(in_pheno.c_str(), ios::in | ios::binary );

    if (!infile.is_open()) {
      sout << "ERROR : Cannote read temporary file " << in_pheno  << endl ;
      exit(-1);
    } 

    // store back values in test_mat
    for( std::size_t m = 0; m < bs_l1; ++m ) 
      for( std::size_t i = 0; i < cv_folds; ++i ) 
        for( std::size_t k = 0; k < cv_sizes[i]; ++k ) 
          infile.read( reinterpret_cast<char *> (&test_mat[ph_eff][i](k,m)), sizeof(double) );

    infile.close();
  }

  // fit model using out-of-sample level 0 predictions from whole data
  if(within_sample_l0){
    betaold = MatrixXd::Zero(bs_l1, 1);

    int niter_cur = 0;
    while(niter_cur++ < niter_max){

      XtWX = MatrixXd::Zero(bs_l1, bs_l1);
      XtWZ = MatrixXd::Zero(bs_l1, 1);

      for(size_t i = 0; i < cv_folds; ++i ) {
        etavec = (test_offset[ph][i] + test_mat[ph_eff][i] * betaold).array();
        pivec = 1 - 1/(etavec.exp() + 1);
        wvec =  pivec * (1 - pivec), 0;
        zvec = (etavec - test_offset[ph][i].array()) + (test_pheno_raw[ph][i].array() - pivec) / wvec;

        XtW = test_mat[ph_eff][i].transpose() * wvec.matrix().asDiagonal();
        XtWX += XtW * test_mat[ph_eff][i];
        XtWZ += XtW * zvec.matrix();
      }
      betanew = (XtWX + tau[val] * ident_l1).llt().solve(XtWZ);
      // compute score
      score = ArrayXd::Zero(betanew.rows());
      for(size_t i = 0; i < cv_folds; ++i ) {
        etavec = (test_offset[ph][i] + test_mat[ph_eff][i] * betanew).array();
        pivec = 1 - 1/(etavec.exp() + 1);
        score += (test_mat[ph_eff][i].transpose() * (test_pheno_raw[ph][i].array() - pivec).matrix()).array();
      }
      score -= tau[val] * betanew.array(); 

      // stopping criterion
      if( score.abs().maxCoeff() < numtol) break;

      betaold = betanew;
    }
  }

  // compute predictor for each chr
  int ctr = 0, chr_ctr = 0; 
  int nn, cum_size_folds;
  map<int, vector<int> >::iterator itr;
  for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) {
    nn = itr->second[1] * n_ridge_l0;
    if(nn > 0) { 
      cum_size_folds = 0;
      for(size_t i = 0; i < cv_folds; ++i ) {
        if(!within_sample_l0) betanew = beta_hat_level_1[ph][i].col(val);
        predictions[0].block(cum_size_folds, chr_ctr, cv_sizes[i], 1) = test_mat[ph_eff][i].block(0, ctr, cv_sizes[i], nn) * betanew.block(ctr, 0, nn, 1);
        cum_size_folds += cv_sizes[i];
      }
      chr_ctr++;
      ctr += nn;
    }
  }

  write_predictions(ph);

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl << endl;
}


void Data::make_predictions_binary_loocv(int ph, int val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  int bs_l1 = total_n_block * n_ridge_l0;
  int ph_eff = write_l0_pred ? 0 : ph;
  double v2;
  string in_pheno;
  ifstream infile;
  MatrixXd XtWX, XtWZ, V1, Xmat_chunk;
  ArrayXd betaold, etavec, pivec, wvec, zvec, betanew, score;
  VectorXd Yvec_chunk;
  LLT<MatrixXd> Hinv;
  MatrixXd beta_final;
  ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);

  uint64 max_bytes = chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (bs_l1 * target_size) matrix ]
  int nchunk = ceil( cv_folds * bs_l1 * sizeof(double) * 1.0 / max_bytes );
  int chunk, size_chunk, target_size = cv_folds / nchunk;
  int j_start;

  // read in level 0 predictions from file
  if(write_l0_pred){

    in_pheno = loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
    infile.open(in_pheno.c_str(), ios::in | ios::binary );

    if (!infile.is_open()) {
      sout << "ERROR : Cannote read temporary file " << in_pheno  << endl ;
      exit(-1);
    } 

    // store back values in test_mat_conc
    infile.read( reinterpret_cast<char *> (&test_mat_conc[ph_eff](0,0)), n_samples * bs_l1 * sizeof(double) );

    infile.close();
  }

  // fit logistic on whole data again for optimal ridge param
  betaold = ArrayXd::Zero(bs_l1);
  int niter_cur = 0;
  while(niter_cur++ < niter_max){
    etavec = (offset_logreg.col(ph) + test_mat_conc[ph_eff] * betaold.matrix()).array();
    pivec = 1 - 1/(etavec.exp() + 1);
    wvec = (masked_indivs.col(ph).array() == 1).select( pivec * (1 - pivec), 0 );
    zvec = (masked_indivs.col(ph).array() == 1).select( (etavec - offset_logreg.col(ph).array()) + (phenotypes_raw.col(ph).array() - pivec) / wvec, 0);
    V1 = test_mat_conc[ph_eff].transpose() * wvec.matrix().asDiagonal();
    XtWX = V1 * test_mat_conc[ph_eff];
    XtWZ = V1 * zvec.matrix();
    Hinv.compute( XtWX + tau[val] * ident_l1 );
    betanew = (Hinv.solve(XtWZ)).array();

    // get the score
    etavec = (offset_logreg.col(ph) + test_mat_conc[ph_eff] * betaold.matrix()).array();
    pivec = 1 - 1/(etavec.exp() + 1);
    score = ( test_mat_conc[ph_eff].transpose() * (masked_indivs.col(ph).array() == 1).select(phenotypes_raw.col(ph).array() - pivec, 0).matrix()).array() ;
    score -= tau[val] * betanew; 

    if( score.abs().maxCoeff() < numtol) break;

    betaold = betanew;
  }
  // compute Hinv 
  etavec = (offset_logreg.col(ph) + test_mat_conc[ph_eff] * betanew.matrix()).array();
  pivec = 1 - 1/(etavec.exp() + 1);
  wvec = pivec * (1 - pivec);
  zvec = (etavec - offset_logreg.col(ph).array()) + (phenotypes_raw.col(ph).array() - pivec) / wvec;
  V1 = test_mat_conc[ph_eff].transpose() * wvec.matrix().asDiagonal();
  XtWX = V1 * test_mat_conc[ph_eff];
  Hinv.compute( XtWX + tau[val] * ident_l1 );

  // loo estimates
  for(chunk = 0; chunk < nchunk; ++chunk ) {
    size_chunk = chunk == nchunk - 1? cv_folds - target_size * chunk : target_size;
    j_start = chunk * target_size;
    if( (chunk == 0) || (chunk == nchunk - 1) ) beta_final.resize(bs_l1, size_chunk);

    Xmat_chunk = test_mat_conc[ph_eff].block(j_start, 0, size_chunk, bs_l1); // n x k
    Yvec_chunk = phenotypes_raw.block(j_start, ph, size_chunk, 1);

    V1 = Hinv.solve( Xmat_chunk.transpose() ); // k x n
    for(size_t i = 0; i < size_chunk; ++i ) {
      v2 = Xmat_chunk.row(i) * V1.col(i); 
      v2 *= wvec(j_start + i); 
      beta_final.col(i) = (betanew - V1.col(i).array() * (Yvec_chunk(i) - pivec(j_start + i)) / (1 - v2)).matrix();
    }

    // compute predictor for each chr
    int ctr = 0, chr_ctr = 0; 
    int nn;
    map<int, vector<int> >::iterator itr;
    for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) {
      nn = itr->second[1] * n_ridge_l0;
      if(nn > 0) { 
        predictions[0].block(j_start, chr_ctr, size_chunk, 1) = ( test_mat_conc[ph_eff].block(j_start, ctr, size_chunk, nn).array() * beta_final.block(ctr, 0, nn, size_chunk).transpose().array() ).matrix().rowwise().sum();
        chr_ctr++;
        ctr += nn;
      }
    }
  }

  write_predictions(ph);

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl << endl;
}

void Data::write_predictions(int ph){
  // output predictions to file
  ofstream ofile;
  map<string, uint64 >::iterator itr_ind;
  string out, id_index;
  uint64 index;

  // for the per chromosome predictions -- only for QT
  if(perchrLOCO) {
    out = out_file + "_" + to_string(ph+1);
    ofile.open(out.c_str());
    sout << "writing file " << out << "..." << flush;

    // enforce all chromosomes are printed
    MatrixXd autosomal_pred = MatrixXd::Zero(predictions[0].rows(), nChrom);

    int chr, nn, chr_ctr = 0;
    map<int, vector<int> >::iterator itr;
    for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) {
      chr  = itr->first;
      nn = itr->second[1];
      if(nn > 0){  
        autosomal_pred.col(chr - 1) = predictions[0].col(chr_ctr);
        ++chr_ctr;
      }
    }

    // header line : FID_IID for all individuals
    ofile << "FID_IID "; 
    for (itr_ind = FID_IID_to_ind.begin(); itr_ind != FID_IID_to_ind.end(); ++itr_ind) {
      id_index = itr_ind->first;
      index = itr_ind->second;

      // check individual is not in exclusion list, otherwise skip
      if( ind_ignore( index ) ) continue;
      ofile << id_index << " ";
    }
    ofile << endl;

    // for each row: print chromosome then blups
    for(chr = 0; chr < nChrom; chr++) {
      ofile << chr + 1 << " ";
      for (itr_ind = FID_IID_to_ind.begin(); itr_ind != FID_IID_to_ind.end(); ++itr_ind) {
        id_index = itr_ind->first;
        index = itr_ind->second;

        // check individual is not in exclusion list, otherwise skip
        if( ind_ignore( index ) ) continue;

        // print blup 
        if( masked_indivs(index, ph) )
          ofile << autosomal_pred(index, chr) << " ";
        else 
          ofile << "NA ";
      }
      ofile << endl;
    }

    ofile.close();
  }

  if(make_loco || binary_mode){
    out = out_file + "_" + to_string(ph+1) + ".loco";
    ofile.open(out.c_str());
    sout << "writing LOCO predictions..." << flush;

    // output LOCO predictions G_loco * beta_loco for each autosomal chr  
    MatrixXd loco_pred (predictions[0].rows(), nChrom);
    loco_pred.colwise() = predictions[0].rowwise().sum();

    int chr, nn, chr_ctr = 0;
    map<int, vector<int> >::iterator itr;
    for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) {
      chr  = itr->first;
      nn = itr->second[1];
      if(nn > 0) { 
        loco_pred.col(chr - 1) -= predictions[0].col(chr_ctr);
        ++chr_ctr;
      }
    }

    // header line : FID_IID for all individuals
    ofile << "FID_IID "; 
    for (itr_ind = FID_IID_to_ind.begin(); itr_ind != FID_IID_to_ind.end(); ++itr_ind) {
      id_index = itr_ind->first;
      index = itr_ind->second;

      // check individual is not in exclusion list, otherwise skip
      if( ind_ignore( index ) ) continue;
      ofile << id_index << " ";
    }
    ofile << endl;

    // print loco predictions for each chromosome
    for(chr = 0; chr < nChrom; chr++) {
      ofile << chr + 1 << " ";
      for (itr_ind = FID_IID_to_ind.begin(); itr_ind != FID_IID_to_ind.end(); ++itr_ind) {
        id_index = itr_ind->first;
        index = itr_ind->second;

        // check individual is not in exclusion list, otherwise skip
        if( ind_ignore( index ) ) continue;

        // print loco prediction 
        if( masked_indivs(index, ph) )
          ofile << loco_pred(index, chr) << " ";
        else 
          ofile << "NA ";
      }
      ofile << endl;
    }

    ofile.close();
  }
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
//// determine if using fast streaming of bgen
/////////////////////////////////////////////////
/////////////////////////////////////////////////

// check if uses Layout 2 (v1.2/1.3) and compressed using zlib's compress() function & check for first SNP if precision for probabilities is 8 bits
void Data::check_bgen(){

  // for non-bgen file input, skip check
  if(file_type != "bgen") return;

  BgenParser bgen_ck;
  bgen_ck.open( bgen_file ) ;
  bool layoutV2 = bgen_ck.get_layout();
  bool zlib_compress = bgen_ck.get_compression();
  if( !layoutV2 || !zlib_compress ){
    streamBGEN = false;
    return;
  }
  uint64 first_snp = bgen_ck.get_position();

  uint minploidy = 0, maxploidy = 0, phasing = 0, bits_prob = 0;
  uint16_t SNPID_size = 0, RSID_size = 0, chromosome_size = 0 , numberOfAlleles = 0 ;
  uint32_t position = 0, allele_size = 0, nindivs = 0;
  string allele, tmp_buffer;

  // check bits only for first snp
  //cout << endl << "Snp1 pos:" << first_snp << endl;
  ifstream bfile;
  bfile.open( bgen_file, ios::in | ios::binary );
  bfile.seekg( first_snp );
  // snpid
  bfile.read( reinterpret_cast<char *> (&SNPID_size), 2 );
  tmp_buffer.resize(SNPID_size);
  bfile.read( reinterpret_cast<char *> (&tmp_buffer[0]), SNPID_size );
  // rsid
  bfile.read( reinterpret_cast<char *> (&RSID_size), 2) ;
  tmp_buffer.resize(RSID_size);
  bfile.read( reinterpret_cast<char *> (&tmp_buffer[0]), RSID_size );
  //cout << "RSID:" << tmp_buffer ;
  // chromosome
  bfile.read( reinterpret_cast<char *> (&chromosome_size), 2 );
  tmp_buffer.resize(chromosome_size);
  bfile.read( reinterpret_cast<char *> (&tmp_buffer[0]), chromosome_size );
  assert( chrStrToInt( tmp_buffer ) > 0 );
  //cout << ",CHR:" << tmp_buffer ;
  // position
  bfile.read( reinterpret_cast<char *> (&position), 4 );
  //cout << ",POS:" << position << endl;
  // number of alleles
  bfile.read( reinterpret_cast<char *> (&numberOfAlleles), 2 );
  assert( numberOfAlleles == 2 ); // only diploid
  //cout << ",Nalleles:" << numberOfAlleles ;
  // alleles
  bfile.read( reinterpret_cast<char *> (&allele_size), 4 );
  tmp_buffer.resize(allele_size);
  bfile.read( reinterpret_cast<char *> (&tmp_buffer[0]), allele_size );
  //cout << ",A0:"<<tmp_buffer ;
  bfile.read( reinterpret_cast<char *> (&allele_size), 4 );
  tmp_buffer.resize(allele_size);
  bfile.read( reinterpret_cast<char *> (&tmp_buffer[0]), allele_size );
  //cout << ",A1:"<<tmp_buffer ;

  // set genotype data block
  vector < uchar > geno_block, geno_block_uncompressed;
  uint32_t size_block = 0, size_block_post_compression = 0;
  bfile.read( reinterpret_cast<char *> (&size_block), 4 );
  bfile.read( reinterpret_cast<char *> (&size_block_post_compression), 4);
  //cout << ",block size:"<<size_block  << ",block size post compress:" << size_block_post_compression << endl;
  geno_block.resize(size_block - 4);
  geno_block_uncompressed.resize(size_block_post_compression);
  bfile.read( reinterpret_cast<char *> (&geno_block[0]), size_block - 4);

  // uncompress the block using zlib
  uLongf dest_size = size_block_post_compression;
  if( (uncompress( &(geno_block_uncompressed[0]), &dest_size, &geno_block[0], size_block - 4) != Z_OK) || (dest_size != size_block_post_compression) ){
    streamBGEN = false;
    return;
  }

// stream to uncompressed block
  uchar *buffer = &geno_block_uncompressed[0];
  // sample size
  std::memcpy(&nindivs, &(buffer[0]), 4);
  //cout << "N:"<< nindivs ;
  assert( nindivs == bgen_ck.number_of_samples() );
  buffer += 4;
  // num alleles
  std::memcpy(&numberOfAlleles, &(buffer[0]), 2);
  //cout << ",allele:"<< numberOfAlleles ;
  assert( numberOfAlleles == 2 );
  buffer += 2;
  // ploidy
  std::memcpy(&minploidy, &(buffer[0]), 1);
  //cout << ",minP:"<< minploidy ;
  assert( minploidy == 2 );
  buffer ++;
  std::memcpy(&maxploidy, &(buffer[0]), 1);
  //cout << ",maxP:"<< maxploidy ;
  assert( maxploidy == 2 );
  buffer ++;

  /* //to identify missing when getting dosages
  vector < uchar > ploidy_n;
  ploidy_n.resize( nindivs );
  std::memcpy(&(ploidy_n[0]), &(buffer[0]), nindivs);
  */
  buffer += nindivs;

  // phasing
  std::memcpy(&phasing, &(buffer[0]), 1);
  //cout << ",phasing:"<< phasing ;
  buffer ++;

  // bits per probability
  std::memcpy(&bits_prob, &(buffer[0]), 1);
  //cout << ",bits:"<< bits_prob ;
  buffer ++;
  if( (phasing != 0) || (bits_prob != 8) ){
    streamBGEN = false;
    return;
  }

  /*
  // get dosages (can compute mean as going along (and identify non-zero entries if SPA is used)
  ArrayXd Geno = ArrayXd::Zero(nindivs);
  bool missing;
  double prob0, prob1;
  // parse genotype probabilities block
  for( uint32_t i = 0; i < nindivs; i++) {
  missing = ((ploidy_n[i]) & 0x80);
  if(missing) {
  Geno(i) = -3;
  buffer+=2;
  continue;
  }
  prob0 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;
  prob1 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;
  Geno(i) = prob1 + 2 * (std::max( 1 - prob0 - prob1, 0.0) );
  }
  cout << "myG*****\n" << Geno.matrix().transpose().array().head(10) << endl;
  */

  return;

}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Testing mode (approx. single-threaded)
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::test_snps() {

  sout << "Association testing mode" << endl;
  chi_squared chisq(1);
  normal nd(0,1);
  double chisq_val, bhat, se_b, pval_log, pval_raw;
  double chisq_thr = quantile(complement(chisq, alpha_pvalue));
  double zcrit = quantile(complement(nd, .025));
  double effect_val, outse_val, outp_val;
  uint64 n_failed_tests = 0;
  uint64 n_ignored_snps = 0;
  uint64 n_skipped_snps = 0;
  bool  has_converged;
  ofstream ofile;
  string out, correction_type;
  vector < string > out_split;
  vector < ofstream > ofile_split;
  MatrixXd WX, GW, sqrt_denum, scaleG_pheno;
  RowVectorXd p_sd;
  std::chrono::high_resolution_clock::time_point t1, t2;

  setNbThreads(threads); // set threads   
  file_read_initialization(); // set up files for reading
  read_pheno_and_cov();   // read phenotype and covariate files  
  blup_read(); // read blups
  set_blocks_for_testing();   // set number of blocks 
  print_usage_info();


  if(test_type == 0) test_string = "ADD";
  else if(test_type == 1) test_string = "DOM";
  else test_string = "REC";
  // write results in one file
  if(!split_by_pheno){
    out = out_file + ".regenie";
    ofile.open(out.c_str());
    // header of output file 
    ofile << "CHROM" << " " << "GENPOS" << " " << "ID" << " " << "ALLELE0" << " " << 
      "ALLELE1" << " " << "A1FREQ" << " " << ((file_type == "bgen")? "INFO ":"") << "TEST" << " ";
    for(size_t i = 0; i < n_pheno; i++) ofile << "BETA.Y" << i + 1 << " " << "SE.Y" << i+1 <<  
      " " << "CHISQ.Y" << i+1 << " " << "LOG10P.Y" << i+1 << " ";  
    ofile << endl;
  } else { 
    // split results in separate files for each phenotype
    out_split.resize( n_pheno );
    ofile_split.resize( n_pheno );

    for(size_t i = 0; i < n_pheno; i++) {
      out_split[i] = out_file + "_" + pheno_names[i] + ".regenie";
      ofile_split[i].open( out_split[i].c_str() );
      // header of output file 
      if(!htp_out){
      ofile_split[i] << "CHROM" << " " << "GENPOS" << " " << "ID" << " " << "ALLELE0" << " " << "ALLELE1" << " " << 
        "A1FREQ" << " " << ((file_type == "bgen")? "INFO ":"") << "TEST" << " " << "BETA" << " " << "SE" << " " << 
        "CHISQ" << " " << "LOG10P" << endl;
      } else {
        ofile_split[i] << "Name" << "\t" << "Chr" << "\t" << "Pos" << "\t" << "Ref" << "	" << "Alt" << "\t" << 
          "Trait" << "\t" << "Cohort" << "\t" << "Model" << "\t" << "Effect" << "\t" << "LCI_Effect" << "\t" << 
          "UCI_Effect" << "\t" << "Pval" << "\t" << "AAF" << "\t" << "Num_Cases"<< "\t" << 
          "Cases_Ref" << "\t" << "Cases_Het" << "\t" << "Cases_Alt" << "\t" << "Num_Controls" << "\t" << 
          "Controls_Ref" << "\t" << "Controls_Het"<< "\t"<< "Controls_Alt" << "\t" << "Info" << endl;
      }
    }
  }
  if(htp_out){
    if(binary_mode & firth) correction_type = "-FIRTH";
    else if(binary_mode & use_SPA) correction_type = "-SPA";
    else if(binary_mode) correction_type = "-LOG";
    else correction_type = "-LR";
    if(skip_blups) model_type = test_string + correction_type;
    else model_type = test_string + "-WGR" + correction_type;
  }


  // set memory for matrices
  Y_hat_p.resize(n_samples, n_pheno);
  Gamma_sqrt.resize(n_samples, n_pheno);
  Xt_Gamma_X_inv.resize(n_pheno);

  sout << " * using minimum MAC of " << min_MAC << " (variants with lower MAC are ignored)" << endl;
  if(firth || use_SPA) {
    sout << " * using ";
    if(firth_approx) sout << "fast ";
    if(firth) sout << "Firth ";
    else sout << "SPA ";
    sout << "correction for logistic regression p-values less than " << alpha_pvalue << endl;
    n_corrected = 0;
  }

  // if testing select chromosomes
  if( select_chrs ) sout << " * user specified to test only on select chromosomes" << endl;
  sout << endl;

  // set covariates for firth
  if(firth){
    covs_firth.resize(n_samples, new_cov.cols() + 1);
    covs_firth.leftCols(new_cov.cols()) << new_cov;
  }


  // start analyzing each chromosome
  int block = 0, chrom, chrom_nsnps, chrom_nb, bs;
  uint64 snp_count = 0;
  map<int, vector<int> >::iterator itr; 
  for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) { 
    chrom = itr->first;
    chrom_nsnps = itr->second[0];
    chrom_nb = itr->second[1];

    if(chrom_nb > 0) {

      // skip chromosome if not in list to keep
      if( select_chrs && !std::count( chrKeep_test.begin(), chrKeep_test.end(), chrom) ) {
        sout << "Chromosome " << chrom << " is skipped" << endl;
        skip_snps( chrom_nsnps );
        block += chrom_nb;
        snp_count += chrom_nsnps;
        n_skipped_snps += chrom_nsnps;
        continue;
      }

      sout << "Chromosome " << chrom << " [" << chrom_nb << " blocks in total]" << endl;
      // read polygenic effect predictions from step 1
      blup_read_chr(chrom);

      // compute phenotype residual (adjusting for BLUP [and covariates for BTs])
      if(!binary_mode){ 

        res = phenotypes - blups;  
        res.array() *= masked_indivs.array(); 

        p_sd = res.colwise().norm(); 
        p_sd.array() /= sqrt(Neff -1); 
        res.array().rowwise() /= p_sd.array();

      } else {

        fit_null_logistic(chrom); // for all phenotypes

        res = phenotypes_raw - Y_hat_p;
        res.array() /= Gamma_sqrt.array();
        res.array() *= masked_indivs.array();

        // if using firth approx., fit null penalized model with only covariates and store the estimates 
        // (used as offset when computing LRT in full model)
        if(firth_approx){
          beta_null_firth.resize(covs_firth.cols(), n_pheno);
          sout << "   -fitting null Firth logistic regression on binary phenotypes..." << flush;
          auto t1 = std::chrono::high_resolution_clock::now();

          for( std::size_t i = 0; i < n_pheno; ++i ) {
            has_converged = fit_firth_logistic(chrom, i, true); 
            if(!has_converged) {
              sout << "ERROR: Firth penalized logistic regression failed to converge for phenotype: " << pheno_names[i] << endl;
              exit(-1);
            }
          }

          sout << "done";
          auto t2 = std::chrono::high_resolution_clock::now();
          auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
          sout << " (" << duration.count() << "ms) "<< endl;
        }

      }
    } else continue;

    // analyze by blocks of SNPs
    for(int bb = 0; bb < chrom_nb ; bb++) {

      bs = block_size;
      if(bb == 0) {
        G.resize(bs,n_samples);
        stats.resize(bs, n_pheno);
        snp_afs.resize(bs, 1);
        if(file_type == "bgen") snp_info.resize(bs, 1);
        if(binary_mode) sqrt_denum.resize(bs, n_pheno);
        else scaleG_pheno.resize(bs, n_pheno);
        if(use_SPA) {
          SPA_pvals.resize(bs, n_pheno);
          snp_flipped.resize(bs);
        }
        if(htp_out) {
          genocounts.resize(bs);
          for( std::size_t j = 0; j < bs; ++j ) genocounts[j].resize(6, n_pheno);
        }
      }
      if((bb +1) * block_size > chrom_nsnps) {
        bs = chrom_nsnps - (bb * block_size) ;
        G.resize(bs,n_samples);
        stats.resize(bs, n_pheno);
        snp_afs.resize(bs, 1);
        if(file_type == "bgen") snp_info.resize(bs, 1);
        if(binary_mode) sqrt_denum.resize(bs, n_pheno);
        else scaleG_pheno.resize(bs, n_pheno);
        if(use_SPA) {
          SPA_pvals.resize(bs, n_pheno);
          snp_flipped.resize(bs);
        }
        if(htp_out) {
          genocounts.resize(bs);
          for( std::size_t j = 0; j < bs; ++j ) genocounts[j].resize(6, n_pheno);
        }
      }
      bad_snps = ArrayXd::Zero(bs);

      // get genotype matrix for block (mean impute)
      get_G(block, bs, chrom);

      if(!binary_mode || firth_approx){
        // residualize and scale genotypes (and identify monomorphic if present)
        // only do it for firth approx. test with BTs (ests. are unchanged for other tests)
        residualize_genotypes();
        n_ignored_snps += bad_snps.sum();
      } else { 
        scale_G = ArrayXd::Ones(bs); // no scaling is applied to the SNPs
        // ensure that snps which failed MAC filter are all polymorphic
        if( (bad_snps==1).count() > 0 ) {
          for(size_t i = 0; i < bs ; i++){
            if(bad_snps(i) == 1) {
              G.row(i).head(10).array() += 1; // +1 to first 10 entries
              if(verbose) sout << "WARNING: Ignoring SNP with low variance.\n";
            }
          }
        }
      }
      

      // perform assoc. testing
      t1 = std::chrono::high_resolution_clock::now();
      sout << "   -computing and writing association test statistics..." << flush;

      if(binary_mode) {

        for( std::size_t i = 0; i < n_pheno; ++i ) {

          // project out covariates from G 
          WX = Gamma_sqrt.col(i).asDiagonal() * new_cov;
          GW = G * Gamma_sqrt.col(i).asDiagonal();
          G_tmp = GW - ((GW * WX) * Xt_Gamma_X_inv[i]) * WX.transpose();
          G_tmp.array().rowwise() *= masked_indivs.col(i).transpose().array();
          denum_tstat = G_tmp.rowwise().squaredNorm(); 
          sqrt_denum.col(i) = denum_tstat.array().sqrt();

          // score test stat for BT
          stats.col(i).array() = (G_tmp * res.col(i)).array() / sqrt_denum.col(i).array();

          if(use_SPA) run_SPA_test(i);
        }

      } else {

        // score test stat for QT
        if( strict_mode )
          stats = (G * res) / sqrt( masked_indivs.col(0).sum() );
        else {
          // compute GtG for each phenotype (different missing patterns)
          for( std::size_t i = 0; i < n_pheno; ++i ) {
            scaleG_pheno.col(i) = ( G.array().rowwise() * masked_indivs.col(i).transpose().array()).square().matrix().rowwise().sum();
          }
          stats = ( (G * res).array() / scaleG_pheno.array().sqrt() ).matrix();
        }

      }

      // write stats to file
      for( std::size_t i = 0; i < bs; ++i ) {

        // don't print anything for monomorphic snps
        if(bad_snps(i) == 1) {
          snp_count++;
          continue;
        }

        if(!split_by_pheno) {
          ofile << (snpinfo[snp_count]).chrom << " " << (snpinfo[snp_count]).physpos << " "<< (snpinfo[snp_count]).ID << " "<< (snpinfo[snp_count]).allele1 << " "<< (snpinfo[snp_count]).allele2 << " " << snp_afs(i, 0) << " " ;
          if(file_type == "bgen") ofile << snp_info(i, 0) << " ";
          ofile << test_string << " ";
        }

        for( std::size_t j = 0; j < n_pheno; ++j ) {
          if(split_by_pheno) {
            if(!htp_out){
              ofile_split[j] << (snpinfo[snp_count]).chrom << " " << (snpinfo[snp_count]).physpos << " "<< (snpinfo[snp_count]).ID << " "<< (snpinfo[snp_count]).allele1 << " "<< (snpinfo[snp_count]).allele2 << " " << snp_afs(i, 0) << " " ;
              if(file_type == "bgen") ofile_split[j] << snp_info(i, 0) << " ";
              ofile_split[j] << test_string << " ";
            } else {
              ofile_split[j] <<  (snpinfo[snp_count]).ID << "\t"<< (snpinfo[snp_count]).chrom << "\t" << (snpinfo[snp_count]).physpos << "\t"<< (snpinfo[snp_count]).allele1 << "\t"<< (snpinfo[snp_count]).allele2 << "\t" << pheno_names[j] << "\t" << cohort_name << "\t" << model_type << "\t";
            }
          }

          chisq_val = stats(i,j) * stats(i,j); 
          // test statistic & pvalue
          if(!use_SPA) {
            pval_log = check_pval(chisq_val, chrom, i, j);
          } else {
            pval_log = SPA_pvals(i, j);
            pval_converged = (pval_log != missing_value_double);
          }

          // summary stats
          if( !binary_mode ){
            // estimate & SE for QT
            if( strict_mode )
              bhat = stats(i,j) * ( scale_Y(j) * p_sd(j)) / ( sqrt(masked_indivs.col(j).sum()) * scale_G(i) ); 
            else
              bhat = stats(i,j) * ( scale_Y(j) * p_sd(j)) / ( sqrt(scaleG_pheno(i,j)) * scale_G(i) ); 
            se_b = bhat / stats(i,j);
          } else {
            // with Firth, get sum. stats from Firth logistic regression
            if( firth && (chisq_val > chisq_thr) && pval_converged ){
              pval_raw = max(nl_dbl_dmin, pow(10, -pval_log)); // to prevent overflow
              chisq_val = quantile(complement(chisq, pval_raw));
              bhat = bhat_firth;
              se_b = se_b_firth;
            } else {
              se_b = 1 / sqrt_denum(i, j);
              // with SPA, calculate test stat based on SPA p-value
              if( use_SPA && (chisq_val > chisq_thr) && pval_converged ){
                pval_raw = max(nl_dbl_dmin, pow(10, -pval_log)); // to prevent overflow
                chisq_val = quantile(complement(chisq, pval_raw));
                bhat = sgn(stats(i,j)) * sqrt(chisq_val);
              } else bhat = stats(i,j);
              bhat *= se_b;
              if( use_SPA && snp_flipped[i] ) bhat *= -1;
            }
            bhat /= scale_G(i);
            se_b /= scale_G(i);
          }

          
          if(!split_by_pheno) 
            ofile << bhat << ' ' << se_b << ' ' << chisq_val << ' ';
          else if(!htp_out) ofile_split[j] << bhat << ' ' << se_b << ' ' << chisq_val << ' ';

          if( pval_converged ) {
            if(!split_by_pheno) ofile << pval_log << ' ';
            else {
              if(!htp_out)  ofile_split[j] << pval_log << ' ';
              else {
                outp_val = max(nl_dbl_dmin, pow(10, - pval_log)); // to prevent overflow
                if(outp_val == 1) outp_val = 1 - 1e-7;
              }
            }
          } else {
            n_failed_tests++;
            //sout << "\nSNP was #" << snp_count+1 << endl;
            if(!split_by_pheno) ofile << "NA" << " ";
            else {
              if(!htp_out)  ofile_split[j] << "NA" << " ";
              else outp_val = -1;
            }
          }

          // for HTPv4 output
          if(htp_out){
            if(!binary_mode){
              ofile_split[j] << bhat << "\t" << (bhat - zcrit * se_b) << "\t" << (bhat + zcrit * se_b) << "\t" << outp_val << "\t";
            } else {
              // for tests that failed, print NA 
              if(outp_val<0){
                ofile_split[j] << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t";
              } else{
                // compute allelic OR
                effect_val = (2*genocounts[i](3,j)+genocounts[i](4,j)+.5)*(2*genocounts[i](2,j)+genocounts[i](1,j)+.5)/(2*genocounts[i](5,j)+genocounts[i](4,j)+.5)/(2*genocounts[i](0,j)+genocounts[i](1,j)+.5);
                // compute SE = log(allelic OR) / zstat
                outse_val = fabs(log(effect_val)) / quantile(complement(nd, outp_val/2 ));
                ofile_split[j] << effect_val << "\t" << effect_val * exp(- zcrit * outse_val) << "\t" << effect_val * exp(zcrit * outse_val) << "\t" << outp_val << "\t";
              }
            }

            // print out AF
            ofile_split[j] << snp_afs(i, 0) << "\t";
            // print counts in cases
            ofile_split[j] << (int) genocounts[i].block(0,j,3,1).sum() << "\t" << (int) genocounts[i](0,j) << "\t" << (int) genocounts[i](1,j) << "\t" << (int) genocounts[i](2,j) << "\t";
            // print counts in controls
            ofile_split[j] << (int) genocounts[i].block(3,j,3,1).sum() << "\t" << (int) genocounts[i](3,j) << "\t" << (int) genocounts[i](4,j) << "\t" << (int) genocounts[i](5,j);

            // info column
            if(outp_val >= 0){
              ofile_split[j] << "\t" << "REGENIE_BETA=" << bhat;
              ofile_split[j] << ";" << "REGENIE_SE=" << se_b;
            } else ofile_split[j] << "\t" << "REGENIE_BETA=NA;REGENIE_SE=NA";
            if(binary_mode) {
              if(outp_val < 0) ofile_split[j]  << ";" << "SE=NA";
              else ofile_split[j]  << ";" << "SE=" << outse_val;
            }
            if(file_type == "bgen") ofile_split[j] << ";" << "INFO=" << snp_info(i,0);
          }

          if(split_by_pheno) ofile_split[j] << endl;
        }
        if(!split_by_pheno) ofile << endl;

        snp_count++;
      }

      sout << "done";
      t2 = std::chrono::high_resolution_clock::now();
      auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
      sout << " (" << duration1.count() << "ms) "<< endl;

      block++;
    }
  }

  sout << endl;
  if(!split_by_pheno){
    sout << "Association results stored in file : " << out << endl;
    ofile.close();
  } else {
    sout << "Association results stored separately for each trait " << ( htp_out ? "(HTPv4 format) " : "" ) << "in files : " << endl;
    for( std::size_t j = 0; j < n_pheno; ++j ) {
      ofile_split[j].close();
      sout << "* [" << out_split[j] << "]" << endl;
    }
    sout << endl;
  }

  if(firth || use_SPA) {
    sout << "Number of tests with ";
    sout << (firth ? "Firth " : "SPA "); 
    sout << "correction : (" << n_corrected << "/" << (snpinfo.size() - n_skipped_snps - n_ignored_snps) * n_pheno << ")" <<  endl;
    sout << "Number of failed tests : (" << n_failed_tests << "/" << n_corrected << ")" << endl;
  }
    sout << "Number of ignored SNPs due to low MAC : " << n_ignored_snps << endl;

} 

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    prep for association test
/////////////////////////////////////////////////
/////////////////////////////////////////////////

// get list of blup files
void Data::blup_read() {
  ifstream blup_list_stream, blupf;
  string line, tmp_pheno;
  std::vector< string > tmp_str_vec ;
  int n_files = 0, tmp_index;
  double in_blup;
  vector<int> read_pheno(n_pheno, 0);

  // allocate memory
  blups = MatrixXd::Zero(n_samples, n_pheno);

  // skip reading if specified by user
  if( skip_blups ) {
    sout << " * no blup predictions given. Simple " << ( binary_mode ? "logistic":"linear" ) << " regression will be performed"<<endl;
    return;
  }

  blup_list_stream.open (blup_file.c_str(), ios::in);
  if (!blup_list_stream.is_open()) {
    sout << "ERROR : " << blup_file  << " not open.\n" ;
    exit(-1);
  } 

  // get list of files containing blups 
  sout << " * loco predictions : [" << blup_file << "] ";
  while (getline(blup_list_stream, line)){
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    // each line contains a phenotype name and the corresponding blup file name
    if( tmp_str_vec.size() != 2 ){
      sout << "ERROR: Incorrectly formatted blup list file : " << blup_file << endl;
      exit(-1);
    }

    // get index of phenotype in phenotype matrix
    vector<string>::iterator it = std::find(pheno_names.begin(), pheno_names.end(), tmp_str_vec[0]);
    if (it == pheno_names.end()) continue; // ignore unrecognized phenotypes

    tmp_index = std::distance(pheno_names.begin(), it);
    pheno_index.push_back(tmp_index);

    // check that phenotype only has one file
    if(read_pheno[tmp_index] != 0){
      sout << "ERROR: Phenotype " << tmp_pheno << " appears more than once in blup list file : " << blup_file << endl;
      exit(1);
    }

    n_files++;
    read_pheno[tmp_index] = 1;
    blup_files.push_back(tmp_str_vec[1]);
  }

  // force all phenotypes in phenotype file to be used 
  if(n_files != n_pheno) {
    sout << "ERROR : Number of files (" << n_files <<")  is not equal to the number of phenotypes.\n" ;
    exit(-1);
  }
  sout << "n_files = " << n_files << endl;
  blup_list_stream.close();

  // read blup file for each phenotype
  for(size_t ph = 0; ph < n_pheno; ph++) {
    int i_pheno = pheno_index[ph];

    sout << "   -file [" <<  blup_files[ph];
    sout << "] for phenotype \'" << pheno_names[i_pheno] << "\'";

    blupf.open (blup_files[ph].c_str(), ios::in);
    if (!blupf.is_open()) {
      sout << "ERROR : " << blup_files[ph]  << " not open.\n" ;
      exit(-1);
    }

    sout << endl;
    blupf.close();
  }

}

void Data::blup_read_chr(int chrom) {
  ifstream blup_list_stream, blupf;
  string line, filename, tmp_pheno;
  std::vector< string > id_strings, tmp_str_vec ;
  int tmp_index;
  double in_blup;
  uint64 indiv_index;

  blups = MatrixXd::Zero(n_samples, n_pheno);

  // skip reading if specified by user
  if( skip_blups ) return;

  sout << "   -reading loco predictions for the chromosome..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();
    
  // read blup file for each phenotype
  for(size_t ph = 0; ph < n_pheno; ph++) {

    int i_pheno = pheno_index[ph];
    ArrayXd read_indiv = ArrayXd::Zero(n_samples);
    blupf.open (blup_files[ph].c_str(), ios::in);

    // check header
    getline (blupf,line);
    boost::algorithm::split(id_strings, line, is_any_of("\t "));
    if( id_strings[0] != "FID_IID") {
      sout << "ERROR: Header of blup file must start with FID_IID." << endl;
      exit(-1);
    }

    // skip to chr
    for (int chr = 1; chr < chrom; chr++) blupf.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    getline (blupf,line);
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    // check number of entries is same as in header
    if(tmp_str_vec.size() != id_strings.size()) {
      sout << "ERROR: blup file for phenotype [" << pheno_names[i_pheno] << "] has different number of entries on line " << chrom + 1 << " compared to the header.\n";
      exit(-1);
    }

    // check starts with chromosome number
    if(chrStrToInt(tmp_str_vec[0]) != chrom) {
      sout << "ERROR: blup file for phenotype [" << pheno_names[i_pheno] << "] start with `" << tmp_str_vec[0]<< "`" <<
        "instead of chromosome number=" << chrom << "." << endl;
      exit(-1);
    }

    // read blup data
    for( size_t filecol = 1; filecol < id_strings.size(); filecol++ ) {

      // ignore sample if it is not in genotype data
      if ( FID_IID_to_ind.find(id_strings[filecol]) == FID_IID_to_ind.end()) continue;
      indiv_index = FID_IID_to_ind[id_strings[filecol]];

      // ignore sample if it is not included in analysis
      if(ind_in_analysis(indiv_index) == 0) continue;

      // check if duplicate
      if( !read_indiv(indiv_index) ){
        read_indiv(indiv_index) = 1;
      } else {
        sout << "ERROR: Individual appears more than once in blup file [" << blup_files[ph] <<"]: FID_IID=" << id_strings[filecol] << endl;
        exit(-1);
      }

      in_blup = convertDouble( tmp_str_vec[filecol] );

      // if blup is NA then individual must be ignored in analysis for the phenotype (ie mask = 0)
      if (in_blup == missing_value_double){
        if( masked_indivs(indiv_index, i_pheno) > 0 ){
          sout << "ERROR: Individual (FID_IID=" << id_strings[filecol] << ") has missing blup prediction at chromosome " << chrom <<" for phenotype " << pheno_names[i_pheno]<< ". ";
          sout << "Either set their phenotype to `NA`, specify to ignore them using option '--remove', or skip reading predictions with option '--ignore-pred'.\n" << err_help ;
          exit(-1);
        };
      } else blups(indiv_index, i_pheno) = in_blup;
    }

    // force all non-masked samples to have loco predictions
    if( (masked_indivs.col(i_pheno).array() == 1).select(read_indiv, 0).sum() < masked_indivs.col(i_pheno).sum() ){
      sout << "ERROR: All samples included in the analysis (for phenotype " << 
        pheno_names[i_pheno]<< ") must have LOCO predictions in file : " << blup_files[ph] << "\n";
      exit(-1);
    }

    blupf.close();
  }

  // mask individuals not in analysis
  //blups.array().colwise() *= ind_in_analysis;
  
  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;
}

void Data::set_blocks_for_testing() {

  total_n_block = 0;
  int blocks_left = n_block;
  map<int, vector<int> >::iterator itr; 
  map<int, vector<int> > m1;
  for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) { 
    //int chrom = itr->first;
    int chrom_nsnps = itr->second[0];
    int nb = ceil((double) chrom_nsnps / block_size);
    if(n_block > 0) {
      if(blocks_left > 0) {
        int minb = min(nb, blocks_left);
        //sout << << endl;
        itr->second[1] = minb;
        total_n_block += minb;
        blocks_left -= minb;
      }
    } else {
      itr->second[1] = nb;
      total_n_block += nb;
    }
    if(itr->second[1] > 0) m1.insert(pair<int, vector<int> >(itr->first, itr->second)); 
  }
  chr_map = m1;

  // summarize block sizes
  sout << left << std::setw(20) << " * # threads" << ": [" << threads << "]" << endl;
  sout << left << std::setw(20) << " * block size" << ": [" << block_size << "]" << endl;
  sout << left << std::setw(20) << " * # blocks" << ": [" << total_n_block << "]" << endl;  
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    functions used for assoc. test
/////////////////////////////////////////////////
/////////////////////////////////////////////////
double Data::check_pval(double tstat, int chrom, int snp,  int ph){

  chi_squared chisq(1);
  double chisq_thr = quantile(chisq, 1 - alpha_pvalue); 
  double pval, logp, LRstat;
  pval_converged = false;

  // if quantitative trait, or firth isn't used, or Tstat < threshold, no correction & return pvalue
  if(!binary_mode || !firth || tstat <= chisq_thr){
    pval = cdf(complement(chisq, tstat));

    // check if pvalue from boost 0 (when tstat is very large)
    // use log10p = log10( 2*f(tstat) ) with f=pdf of chisq(1)
    if(pval == 0) logp = log10(2) - 0.5 * log10( 2 * M_PI * tstat ) - 0.5 * tstat * M_LOG10E ; 
    else logp = log10(pval);

    pval_converged = true;
    logp *= -1;
    return logp;
  }
  n_corrected++;

  // firth
  bhat_firth = 0, se_b_firth = 0;
  LRstat = run_firth_correction(chrom, snp, ph);

  if(pval_converged){
    pval = cdf(complement(chisq, LRstat));

    if(pval == 0) logp = log10(2) - 0.5 * log10( 2 * M_PI * LRstat ) - 0.5 * LRstat * M_LOG10E ; 
    else logp = log10(pval);
    logp *= -1;
    return logp;
  } 

  return missing_value_double;
}

double Data::run_firth_correction(int chrom, int snp,  int ph){

  bool converged_l0, converged_l1; 
  double deviance_l0, deviance_l1, dif_deviance;

  // add tested SNP to set of covariates
  covs_firth.rightCols(1) = G.row(snp).transpose();

  // obtain null deviance (set SNP effect to 0 and compute max. pen. LL)
  if(!firth_approx){
    converged_l0 = fit_firth_logistic(chrom, ph, true);
    if(!converged_l0) return missing_value_double;
  }

  // fit full model and compute deviance
  converged_l1 = fit_firth_logistic(chrom, ph, false);
  if(!converged_l1) return missing_value_double;

  // compute LR stat
  dif_deviance = deviance_logistic;
  //sout << endl << "SNP " << snp+1 << ": Deviance l0 =" <<  deviance_l0 << ", Deviance l1" << deviance_l1 << "(Diff = " << dif_deviance << ")";

  // in case of numerical errors
  if( dif_deviance < 0 ) {
    return missing_value_double;
  }

  pval_converged = true;
  return dif_deviance;
}

bool Data::fit_firth_logistic(int chrom, int ph, bool null_fit) {
  // if firth is used, fit based on penalized log-likelihood

  int niter_cur, niter_search, col_incl;
  double dev_old, denum, mx, deviance_l0 = 0;
  int maxstep_firth = null_fit ? maxstep_null : maxstep;
  int niter_firth = null_fit ? niter_max_firth_null : niter_max_firth;
  ArrayXd Y1, hvec, mod_score;
  ArrayXd betaold, betanew, step_size, etavec, pivec, wvec;
  MatrixXd Xmat, XtW, XtWX;
  ColPivHouseholderQR<MatrixXd> qr, qrX;
  hvec = ArrayXd::Zero(n_samples);
  Y1 = phenotypes_raw.col(ph).array() * masked_indivs.col(ph).array();

  if(firth_approx){
    if(null_fit){
      Xmat = new_cov; // only covariates
    } else {
      Xmat = covs_firth.rightCols(1); // only tested SNP
    }
    col_incl = Xmat.cols();
  } else {
    Xmat = covs_firth; // covariates + tested SNP
    col_incl = Xmat.cols();
    if( null_fit ) col_incl--;
  }

  // mask individuals
  Xmat.array().colwise() *= masked_indivs.col(ph).array();

  // with firth approx. => trial 1: use maxstep_null
  // trial 2 => use fallback options (increase maxstep & niter)
  for( size_t trial = 0; trial < 2; trial++){

    // starting values
    if(null_fit){

      betaold = ArrayXd::Zero(Xmat.cols());
      betaold(0) = ( 0.5 + Y1.sum())  / (Neff(ph) + 1);
      betaold(0) = log( betaold(0) / (1 - betaold(0) ));

      // LOCO prediction is offset
      betaold(0) -= (blups.col(ph).array() * masked_indivs.col(ph).array()).mean();

    } else {

      if(firth_approx) betaold = ArrayXd::Zero(col_incl); // only estimating effect of tested SNP
      else betaold = beta_null_firth.array();

    }
    betanew = ArrayXd::Zero( betaold.size() );

    // get the corresponding deviance
    etavec = (Xmat * betaold.matrix()).array();
    etavec += blups.col(ph).array() * masked_indivs.col(ph).array();
    // covariate effects added as offset in firth approx. (last entry of beta_null_firth = 0)
    if( firth_approx && !null_fit ) etavec += ((new_cov.array().colwise() * masked_indivs.col(ph).array()).matrix() * beta_null_firth.block(0,ph,new_cov.cols(),1)).array(); 
    // fitted probabilities
    pivec = 1 - 1 / (etavec.exp() + 1) ;
    wvec = (masked_indivs.col(ph).array() == 1).select( ( pivec * (1 - pivec) ).sqrt(), 0);
    XtW = Xmat.transpose() * wvec.matrix().asDiagonal();
    XtWX = XtW * XtW.transpose();
    qr.compute(XtWX);
    // use penalized log-lik
    dev_old = (masked_indivs.col(ph).array() == 1).select( (Y1 * pivec.log() + (1-Y1) * (1-pivec).log() ), 0).sum();
    dev_old += 0.5 * qr.logAbsDeterminant();
    dev_old *= -2;

    // at niter=0 (i.e. betaSNP=0) this is null deviance
    if( !null_fit ) deviance_l0 = dev_old;

    // solve S'(beta) = S(beta) + X'(h*(0.5-p)) = 0
    niter_cur = 0;
    while(niter_cur++ < niter_firth){

      ////////  compute step size
      etavec = (Xmat * betaold.matrix()).array();
      etavec += blups.col(ph).array() * masked_indivs.col(ph).array();
      if( firth_approx && !null_fit ) etavec += ((new_cov.array().colwise() * masked_indivs.col(ph).array()).matrix() * beta_null_firth.block(0,ph,new_cov.cols(),1)).array(); 

      pivec = 1 - 1 / (etavec.exp() + 1) ;
      wvec = (masked_indivs.col(ph).array() == 1).select( ( pivec * (1 - pivec) ).sqrt(), 0);
      XtW = Xmat.transpose() * wvec.matrix().asDiagonal();
      XtWX = XtW * XtW.transpose();
      qr.compute(XtWX);
      if(!firth_approx && null_fit )  qrX.compute(XtWX.block(0, 0, col_incl, col_incl));

      // compute diag(H), H = U(U'U)^{-1}U', U = Gamma^(1/2)X
      hvec = (qr.solve(XtW).array() * XtW.array() ).colwise().sum();

      // modified score for beta
      mod_score = (Xmat.leftCols(col_incl).transpose() * (masked_indivs.col(ph).array() == 1).select( Y1 - pivec + hvec * (0.5 - pivec), 0).matrix() ).array();

      // step size
      if(!firth_approx && null_fit )
        step_size = qrX.solve( mod_score.matrix() ).array();
      else
        step_size = qr.solve( mod_score.matrix() ).array();

      // force absolute step size to be less than maxstep for each entry of beta
      mx = step_size.abs().maxCoeff() / maxstep_firth;
      if( mx > 1 ) step_size /= mx;

      // start step-halving and stop when deviance decreases 
      denum = 1;
      for( size_t niter_search = 1; niter_search <= niter_max_line_search; niter_search++ ){

        // adjusted step size
        step_size /= denum;

        ///////// compute corresponding deviance
        betanew.head(col_incl) = betaold.head(col_incl) + step_size;
        etavec = (Xmat * betanew.matrix()).array();
        etavec += blups.col(ph).array() * masked_indivs.col(ph).array();
        if( firth_approx && !null_fit ) etavec += ((new_cov.array().colwise() * masked_indivs.col(ph).array()).matrix() * beta_null_firth.block(0,ph,new_cov.cols(),1)).array(); 

        pivec = 1 - 1 / (etavec.exp() + 1) ;
        wvec = (masked_indivs.col(ph).array() == 1).select( ( pivec * (1 - pivec) ).sqrt(), 0);
        XtW = Xmat.transpose() * wvec.matrix().asDiagonal();
        XtWX = XtW * XtW.transpose();
        qr.compute(XtWX);

        deviance_logistic = (masked_indivs.col(ph).array() == 1).select( (Y1 * pivec.log() + (1-Y1) * (1-pivec).log() ), 0).sum();
        deviance_logistic += 0.5 * qr.logAbsDeterminant();
        deviance_logistic *= -2;

        //sout << "\n["<<niter_cur << " - " << niter_search <<"]  denum =" << denum << ";\n step =" << step_size.matrix().transpose().array() / denum<<"; \nbeta=" << betanew.matrix().transpose().array() << ";\n Lnew= " << deviance_logistic << " vs L0="<< dev_old << ";score="<< mod_score<< endl;
        if( deviance_logistic < dev_old + numtol ) break;
        denum *= 2;
      }

      betaold.head(col_incl) += step_size;
      dev_old = deviance_logistic;

      // stopping criterion using modified score function
      if( mod_score.abs().maxCoeff() < numtol_firth) break;

    }

    if(firth_approx && null_fit){ // only retry for firth approx null model
      if(!fix_maxstep_null) { // don't retry with user-given settings
        if( niter_cur > niter_firth ){ // if failed to converge
          sout << "WARNING: Logistic regression with Firth correction did not converge (maximum step size=" << maxstep_firth <<";maximum number of iterations=" << niter_firth<<").";
          maxstep_firth = retry_maxstep_firth;
          niter_firth = retry_niter_firth;
          sout << "Retrying with fallback parameters: (maximum step size=" << maxstep_firth <<";maximum number of iterations=" << niter_firth<<").\n";
          continue;
        }
      }
    }

    break;
  }

  // If didn't converge
  if(niter_cur > niter_firth){
    if(verbose && !firth_approx) sout << "WARNING: Logistic regression with Firth correction did not converge!\n";
    return false;
  }
  // sout << "\nNiter = " << niter_cur << " : " << mod_score.matrix().transpose() << endl;

  if(null_fit) {
    if(firth_approx) beta_null_firth.block(0,ph,betaold.size(),1) = betaold.matrix();
    else beta_null_firth = betaold.matrix();
  } else {
    // compute beta_hat & SE
    bhat_firth = betaold.tail(1)(0);
    se_b_firth = sqrt( qr.inverse().diagonal().tail(1)(0) );

    // compute LRT test stat. 
    deviance_logistic -= deviance_l0;
    deviance_logistic *= -1;

    if( deviance_logistic < 0 ) return false;
  }

  return true;
}

void Data::run_SPA_test(int ph){

  int bs = G_tmp.rows(), index_j, nnz;
  double pval, logp, zstat, zstat_sq, score_num, tval, limK1_low, limK1_high, root_K1, sumG;
  chi_squared chisq(1);
  double chisq_thr = quantile(chisq, 1 - alpha_pvalue); 
  ArrayXd Gmu;

  for(std::size_t snp = 0; snp < bs; ++snp) {

    // avoid monomorphic snps
    if(bad_snps(snp) == 1){
      SPA_pvals(snp, ph) = missing_value_double;
      continue;
    }

    // check test statistic relative to threshold
    zstat = stats(snp, ph);
    zstat_sq = zstat * zstat;
    if( zstat_sq <= chisq_thr){
      pval = cdf(complement(chisq, zstat_sq));
      if(pval == 0) logp = log10(2) - 0.5 * log10( 2 * M_PI * zstat_sq ) - 0.5 * zstat_sq * M_LOG10E ; 
      else logp = log10(pval);
      logp *= -1;
      SPA_pvals(snp, ph) = logp;
      continue;
    }
    n_corrected++;

    // start SPA - check how many non-zero entries there are
    nnz = 0;
    for( std::size_t j = 0; j < non_zero_indices_G[snp].size(); ++j ) {
      index_j = non_zero_indices_G[snp][j];
      if(!masked_indivs(index_j,ph)) continue;
      nnz++;
    }
    fastSPA = nnz < 0.5 * Neff(ph);

    // compute needed quantities
    val_c = sqrt( denum_tstat(snp) );  // sqrt( G'WG )
    score_num = zstat * val_c;
    Gmod = G_tmp.row(snp).transpose().array() / Gamma_sqrt.col(ph).array() * masked_indivs.col(ph).array();
    Gmu = Gmod * Y_hat_p.col(ph).array();
    val_a = Gmu.sum(); 
    sumG = Gmod.sum();

    if(fastSPA){
      val_b = denum_tstat(snp);
      val_d = 0;
      for( std::size_t j = 0; j < non_zero_indices_G[snp].size(); ++j ) {
        index_j = non_zero_indices_G[snp][j];
        if(!masked_indivs(index_j,ph)) continue;
         val_b -= G_tmp(snp, index_j) * G_tmp(snp, index_j);
         val_d += Gmu(index_j);
      }
    }

    // check if K'(t)= s can be solved 
    limK1_low = (Gmod < 0).select(Gmod, 0 ).sum() - val_a ;
    limK1_high = (Gmod > 0).select(Gmod, 0 ).sum() - val_a ;
    if( score_num < limK1_low || score_num > limK1_high ){
      if(verbose) sout << "WARNING: SPA failed (solution to K'(t)=s is infinite)";
      SPA_pvals(snp, ph) = missing_value_double;
      continue;
    }

    // keep track of whether obs stat is positive
    pos_score = zstat  > 0;
    tval = fabs(zstat);

    // solve K'(t)= tval using a mix of Newton-Raphson and bisection method
    root_K1 = solve_K1(tval, pos_score, fastSPA, snp, ph);
    if( root_K1 == missing_value_double ){
      SPA_pvals(snp, ph) = missing_value_double;
      continue;
    }

    // compute pvalue
    SPA_pvals(snp, ph) = get_SPA_pvalue(root_K1, tval, pos_score, fastSPA, snp, ph);
  }

}

double Data::get_SPA_pvalue(double root, double tval, bool val_pos, bool use_SPAfast, int snp, int ph){

  int lambda = val_pos ? 1 : -1; // if score is negative, adjust K and K''
  double kval, k2val, wval, vval, rval, pvalue;
  normal nd(0,1);

  kval = use_SPAfast ? compute_K_fast(lambda * root, snp, ph) : compute_K(lambda * root, ph, 0);
  k2val = use_SPAfast ? compute_K2_fast(lambda * root, snp, ph) : compute_K2(lambda * root, ph, 0);

  wval = sqrt( 2 * ( root * tval - kval ) );
  vval = root * sqrt( k2val );
  if(vval == 0) return 0; // implies score = 0 so pval = 1

  rval = wval + log( vval / wval ) / wval;
  pvalue = cdf(complement(nd, rval ));
  pvalue *= 2;

  if(pvalue > 1) { // SPA can fail for SNPs with very low counts 
    if(verbose) sout << "WARNING: SPA correction failed (resulted in p-value > 1).\n";
    return missing_value_double;
    //return 0;  
  }

  return -log10( pvalue );
}

double Data::solve_K1(double tval, bool val_pos, bool use_SPAfast, int snp, int ph){

  int niter_cur;
  int lambda = val_pos ? 1 : -1; // if score is negative, adjust K' and K''
  double min_x, max_x, t_old, f_old, t_new = -1, f_new, hess;

  niter_cur = 0;
  min_x = 0, max_x = std::numeric_limits<double>::infinity();
  t_old = 0;
  f_old = use_SPAfast ? compute_K1_fast(lambda * t_old, snp, ph) : compute_K1(lambda * t_old, ph, 0);
  f_old *= lambda;
  f_old -= tval; 

  while( niter_cur++ < niter_max_spa ){

    hess = use_SPAfast ? compute_K2_fast(lambda * t_old, snp, ph) : compute_K2(lambda * t_old, ph, 0);
    t_new = t_old - f_old / hess;
    f_new = use_SPAfast ? compute_K1_fast(lambda * t_new, snp, ph) : compute_K1(lambda * t_new, ph, 0);
    f_new *= lambda;
    f_new -= tval;

    if( fabs( f_new ) < tol_spa ) break;

    // update bounds on root
    if( t_new && (t_new > min_x) && (t_new < max_x) ){
      if( f_new > 0) max_x = t_new;
      else min_x = t_new;
    } else{ // bisection method if t_new went out of bounds and re-compute f_new
      t_new = ( min_x + max_x ) / 2;
      // if( fabs( min_x - t_new ) < tol_spa ) break;
      f_new = use_SPAfast ? compute_K1_fast(lambda * t_new, snp, ph) : compute_K1(lambda * t_new, ph, 0);
      f_new *= lambda;
      f_new -= tval;
    }

    t_old = t_new;
    f_old = f_new;
  }

  // If didn't converge
  if( niter_cur > niter_max_spa ){
    if(verbose) sout << "WARNING: SPA did not converge to root for K'(t)=s.\n";
    return missing_value_double;
  }
    //sout << "#iterations = " << niter_cur << "; f= " << f_new << endl;

  return t_new;
}

double Data::compute_K(double t, int ph, int dummy){

  double val = (1 - Y_hat_p.col(ph).array() + Y_hat_p.col(ph).array() * ( t / val_c * Gmod ).exp() ).log().sum() - t * val_a / val_c;

  return val;
}

double Data::compute_K_fast(double t, int snp, int ph){

  uint64 index_j;
  double val = 0;

  for( std::size_t j = 0; j < non_zero_indices_G[snp].size(); ++j ) {
    index_j = non_zero_indices_G[snp][j];
    val += log( 1 - Y_hat_p(index_j,ph) + Y_hat_p(index_j,ph) * exp( t / val_c * Gmod(index_j)) );
  }
  val += -t * val_d / val_c + t * t / 2 / denum_tstat(snp) * val_b;

  return val;
}

double Data::compute_K1(double t, int ph, int dummy){

  double val = ( ( Gmod * Y_hat_p.col(ph).array() / val_c ) / ( Y_hat_p.col(ph).array() + (1 - Y_hat_p.col(ph).array()) * ( -t / val_c * Gmod).exp() ) ).sum();
 val -= val_a / val_c;

  return val;
}

double Data::compute_K1_fast(double t, int snp, int ph){

  uint64 index_j;
  double val = 0;

  for( std::size_t j = 0; j < non_zero_indices_G[snp].size(); ++j ) {
    index_j = non_zero_indices_G[snp][j];
    val += ( Gmod(index_j) * Y_hat_p(index_j,ph) / val_c ) / ( Y_hat_p(index_j,ph) + (1 - Y_hat_p(index_j,ph)) * exp( -t / val_c * Gmod(index_j)) );
  }
  val += -val_d / val_c + t / denum_tstat(snp) * val_b;

  return val;
}

double Data::compute_K2(double t, int ph, int dummy){

  double val = ( ( Gmod.square() * Gamma_sqrt.col(ph).array().square() / (val_c*val_c) * ( -t / val_c * Gmod).exp()) / ( Y_hat_p.col(ph).array() + (1 - Y_hat_p.col(ph).array()) * ( -t / val_c * Gmod).exp() ).square() ).sum();

  return val;
}
  
double Data::compute_K2_fast(double t, int snp, int ph){

  uint64 index_j;
  double val = 0, denum;

  for( std::size_t j = 0; j < non_zero_indices_G[snp].size(); ++j ) {
    index_j = non_zero_indices_G[snp][j];
    denum = Y_hat_p(index_j,ph) + (1 - Y_hat_p(index_j,ph)) * exp( -t / val_c * Gmod(index_j));
    val += ( Gmod(index_j) * Gmod(index_j) * Gamma_sqrt(index_j,ph) * Gamma_sqrt(index_j,ph) * exp( -t / val_c * Gmod(index_j)) / val_c / val_c ) / (denum * denum);
  }
  val += val_b / denum_tstat(snp);

  return val;
}

 
/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Testing mode (multi-threaded)
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::test_snps_fast() {

  sout << "Association testing mode"; 
  std::chrono::high_resolution_clock::time_point t1, t2;
  normal nd(0,1);
  double zcrit = quantile(complement(nd, .025));
  double effect_val, outse_val, outp_val;

#if defined(_OPENMP)
  omp_set_num_threads(threads); // set threads in OpenMP
  sout << " with multithreading using OpenMP";
#endif
  sout << endl;

  setNbThreads(threads);
  file_read_initialization(); // set up files for reading
  read_pheno_and_cov();   // read phenotype and covariate files  
  blup_read(); // read blups
  set_blocks_for_testing();   // set number of blocks 
  print_usage_info();

  sout << " * using minimum MAC of " << min_MAC << " (variants with lower MAC are ignored)" << endl;
  if(firth || use_SPA) {
    sout << " * using " << (firth_approx ? "fast ": "") << (firth ? "Firth ": "SPA ");
    sout << "correction for logistic regression p-values less than " << alpha_pvalue << endl;
    n_corrected = 0;
  }
  // if testing select chromosomes
  if( select_chrs ) sout << " * user specified to test only on select chromosomes" << endl;
  sout << endl;


  // output files
  ofstream ofile;
  string out, correction_type;
  vector < string > out_split;
  vector < ofstream > ofile_split;
  if(test_type == 0) test_string = "ADD";
  else if(test_type == 1) test_string = "DOM";
  else test_string = "REC";
  if(!split_by_pheno){ // write results in one file
    out = out_file + ".regenie";
    ofile.open(out.c_str());
    // header of output file 
    ofile << "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ " << ((file_type == "bgen")? "INFO ":"") << "TEST ";
    for(size_t i = 0; i < n_pheno; i++) ofile << "BETA.Y" << i + 1 << " SE.Y" << i+1 <<  " CHISQ.Y" << i+1 << " LOG10P.Y" << i+1 << " ";  
    ofile << endl;
  } else {  // split results in separate files for each phenotype
    out_split.resize( n_pheno );
    ofile_split.resize( n_pheno );

    for(size_t i = 0; i < n_pheno; i++) {
      out_split[i] = out_file + "_" + pheno_names[i] + ".regenie";
      ofile_split[i].open( out_split[i].c_str() );
      // header of output file 
      if(!htp_out){
        ofile_split[i] << "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ " << ((file_type == "bgen")? "INFO ":"") << "TEST BETA SE CHISQ LOG10P" << endl;
      } else {
        ofile_split[i] << "Name" << "\t" << "Chr" << "\t" << "Pos" << "\t" << "Ref" << "\t" << "Alt" << "\t" << "Trait" << "\t" << "Cohort" << "\t" << "Model" << "\t" << "Effect" << "\t" << "LCI_Effect" << "\t" << "UCI_Effect" << "\t" << "Pval" << "\t" << "AAF" << "\t" << "Num_Cases"<< "\t" << "Cases_Ref" << "\t" << "Cases_Het" << "\t" << "Cases_Alt" << "\t" << "Num_Controls" << "\t" << "Controls_Ref" << "\t" << "Controls_Het"<< "\t"<< "Controls_Alt" << "\t" << "Info" << endl;
      }
    }
  }

  if(htp_out){
    if(binary_mode & firth) correction_type = "-FIRTH";
    else if(binary_mode & use_SPA) correction_type = "-SPA";
    else if(binary_mode) correction_type = "-LOG";
    else correction_type = "-LR";
    if(skip_blups) model_type = test_string + correction_type;
    else model_type = test_string + "-WGR" + correction_type;
  }



  // start analyzing each chromosome
  int block = 0, chrom, chrom_nsnps, chrom_nb, bs;
  bool has_converged;
  tally snp_tally;
  vector< variant_block > block_info;
  Y_hat_p.resize(n_samples, n_pheno);
  Gamma_sqrt.resize(n_samples, n_pheno);
  Xt_Gamma_X_inv.resize(n_pheno);
  if(firth_approx){ // set covariates for firth
    covs_firth.resize(n_samples, new_cov.cols() + 1);
    covs_firth.leftCols(new_cov.cols()) << new_cov;
  }


  map<int, vector<int> >::iterator itr; 
  for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) { 
    chrom = itr->first;
    chrom_nsnps = itr->second[0];
    chrom_nb = itr->second[1];

    if(chrom_nb > 0) {

      // skip chromosome if not in list to keep
      if( select_chrs && !std::count( chrKeep_test.begin(), chrKeep_test.end(), chrom) ) {
        sout << "Chromosome " << chrom << " is skipped" << endl;
        block += chrom_nb;
        snp_tally.snp_count += chrom_nsnps;
        snp_tally.n_skipped_snps += chrom_nsnps;
        if((file_type == "bed")) bed_ifstream.seekg(chrom_nsnps * bed_block_size, ios_base::cur);
        continue;
      }

      sout << "Chromosome " << chrom << " [" << chrom_nb << " blocks in total]" << endl;
      // read polygenic effect predictions from step 1
      blup_read_chr(chrom);

      // compute phenotype residual (adjusting for BLUP [and covariates for BTs])
      if(!binary_mode){ 

        res = phenotypes - blups;  
        res.array() *= masked_indivs.array(); 

        p_sd_yres = res.colwise().norm(); 
        p_sd_yres.array() /= sqrt(Neff -1); 
        res.array().rowwise() /= p_sd_yres.array();

      } else {

        fit_null_logistic(chrom); // for all phenotypes

        res = phenotypes_raw - Y_hat_p;
        res.array() /= Gamma_sqrt.array();
        res.array() *= masked_indivs.array();

        // if using firth approximation, fit null penalized model with only covariates and store the estimates (to be used as offset when computing LRT in full model)
        if(firth_approx){
          beta_null_firth.resize(covs_firth.cols(), n_pheno);
          sout << "   -fitting null Firth logistic regression on binary phenotypes..." << flush;
          auto t1 = std::chrono::high_resolution_clock::now();

          for( std::size_t i = 0; i < n_pheno; ++i ) {
            has_converged = fit_firth_logistic(chrom, i, true); 
            if(!has_converged) {
              sout << "ERROR: Firth penalized logistic regression failed to converge for phenotype: " << pheno_names[i] << ". ";
              sout << "Try decreasing the maximum step size using `--maxstep-null` (currently=" << (fix_maxstep_null ? maxstep_null : retry_maxstep_firth)<< ") " <<
                "and increasing the maximum number of iterations using `--maxiter-null` (currently=" << (fix_maxstep_null ? niter_max_firth_null : retry_niter_firth) << ").\n";
              exit(-1);
            }
          }
          sout << "done";
          auto t2 = std::chrono::high_resolution_clock::now();
          auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
          sout << " (" << duration.count() << "ms) "<< endl;
        }

      }
    } else continue;

    // analyze by blocks of SNPs
    for(int bb = 0; bb < chrom_nb ; bb++) {

      bs = block_size;
      if(bb == 0) {
        block_info.resize(bs);
      }
      if((bb +1) * block_size > chrom_nsnps) {
        bs = chrom_nsnps - (bb * block_size);
        block_info.resize(bs);
      }

      // read SNP, impute missing & compute association test statistic
      sout << " block [" << block + 1 << "] : " << flush; 
      analyze_block(chrom, bs, &snp_tally, block_info);

      // print the results
      for(size_t isnp = 0; isnp < bs; isnp++) {
        uint64 snpindex = snp_tally.snp_count + isnp;

        if( block_info[isnp].ignored ) continue;

        if(!split_by_pheno) {
          ofile << (snpinfo[snpindex]).chrom << " " << (snpinfo[snpindex]).physpos << " "<< (snpinfo[snpindex]).ID << " "<< (snpinfo[snpindex]).allele1 << " "<< (snpinfo[snpindex]).allele2 << " " << block_info[isnp].af << " " ;
          if(file_type == "bgen") ofile << block_info[isnp].info << " ";
          ofile << test_string << " ";
        }

        for(size_t j = 0; j < n_pheno; ++j) {
          if( (firth || use_SPA) && block_info[isnp].is_corrected[j] ) n_corrected++;

          if(split_by_pheno) {
            if(!htp_out){
              ofile_split[j] << (snpinfo[snpindex]).chrom << " " << (snpinfo[snpindex]).physpos << " "<< (snpinfo[snpindex]).ID << " "<< (snpinfo[snpindex]).allele1 << " "<< (snpinfo[snpindex]).allele2 << " " << block_info[isnp].af << " " ;
              if(file_type == "bgen") ofile_split[j] << block_info[isnp].info << " ";
              ofile_split[j] << test_string << " ";
            } else {
              ofile_split[j] <<  (snpinfo[snpindex]).ID << "\t"<< (snpinfo[snpindex]).chrom << "\t" << (snpinfo[snpindex]).physpos << "\t"<< (snpinfo[snpindex]).allele1 << "\t"<< (snpinfo[snpindex]).allele2 << "\t" << pheno_names[j] << "\t" << cohort_name << "\t" << model_type << "\t";
            }
          }

          if(!split_by_pheno) {
            ofile << block_info[isnp].bhat(j) << ' ' << block_info[isnp].se_b(j)  << ' ' << block_info[isnp].chisq_val(j)  << ' ';
          } else if(!htp_out) ofile_split[j] << block_info[isnp].bhat(j)  << ' ' << block_info[isnp].se_b(j)  << ' ' << block_info[isnp].chisq_val(j)  << ' ';

          if( !block_info[isnp].test_fail[j] ) {
            if(!split_by_pheno) ofile << block_info[isnp].pval_log(j)  << ' ';
            else {
              if(!htp_out)  ofile_split[j] << block_info[isnp].pval_log(j)  << ' ';
              else  {
                outp_val = max(nl_dbl_dmin, pow(10, -block_info[isnp].pval_log(j))); // to prevent overflow
                if(outp_val == 1) outp_val = 1 - 1e-7;
              }
            }
          } else {
            snp_tally.n_failed_tests++;
            if(!split_by_pheno) ofile << "NA" << " ";
            else {
              if(!htp_out) ofile_split[j] << "NA" << " ";
              else outp_val = -1;
            }
          }

          // for HTPv4 output
          if(htp_out){
            if(!binary_mode){
              ofile_split[j] << block_info[isnp].bhat(j) << "\t" << block_info[isnp].bhat(j) - zcrit * block_info[isnp].se_b(j) << "\t" << block_info[isnp].bhat(j) + zcrit * block_info[isnp].se_b(j) << "\t" << outp_val << "\t";
            } else {
              // for tests that failed, print NA 
              if(outp_val<0){
                ofile_split[j] << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t";
              } else{
                // compute allelic OR
                effect_val = (2*block_info[isnp].genocounts(3,j)+block_info[isnp].genocounts(4,j)+.5)*(2*block_info[isnp].genocounts(2,j)+block_info[isnp].genocounts(1,j)+.5)/(2*block_info[isnp].genocounts(5,j)+block_info[isnp].genocounts(4,j)+.5)/(2*block_info[isnp].genocounts(0,j)+block_info[isnp].genocounts(1,j)+.5);
                // compute SE = log(allelic OR) / zstat
                outse_val = fabs(log(effect_val)) / quantile(complement(nd, outp_val/2 ));
                ofile_split[j] << effect_val << "\t" << effect_val * exp(- zcrit * outse_val) << "\t" << effect_val * exp(zcrit * outse_val) << "\t" << outp_val << "\t";
              }
            }

            // print out AF
            ofile_split[j] << block_info[isnp].af << "\t";
            // print counts in cases
            ofile_split[j] << (int) block_info[isnp].genocounts.block(0,j,3,1).sum() << "\t" << (int) block_info[isnp].genocounts(0,j) << "\t" << (int) block_info[isnp].genocounts(1,j) << "\t" << (int) block_info[isnp].genocounts(2,j) << "\t";
            // print counts in controls
            ofile_split[j] << (int) block_info[isnp].genocounts.block(3,j,3,1).sum() << "\t" << (int) block_info[isnp].genocounts(3,j) << "\t" << (int) block_info[isnp].genocounts(4,j) << "\t" << (int) block_info[isnp].genocounts(5,j) ;

            // info column
            if(outp_val >= 0){
              ofile_split[j] << "\t" << "REGENIE_BETA=" << block_info[isnp].bhat(j);
              ofile_split[j] << ";" << "REGENIE_SE=" << block_info[isnp].se_b(j);
            } else ofile_split[j] << "\t" << "REGENIE_BETA=NA;REGENIE_SE=NA";
            if(binary_mode) {
              if(outp_val < 0) ofile_split[j]  << ";" << "SE=NA";
              else ofile_split[j]  << ";" << "SE=" << outse_val;
            }
            if(file_type == "bgen") ofile_split[j] << ";" << "INFO=" << block_info[isnp].info;

          }

          if(split_by_pheno) ofile_split[j] << endl;
        }

        if(!split_by_pheno) ofile << endl;

      }

      snp_tally.snp_count += bs;
      block++;
    }
  }

  sout << endl;

  if(!split_by_pheno){
    sout << "Association results stored in file : " << out << endl;
    ofile.close();
  } else {
    sout << "Association results stored separately for each trait " << ( htp_out ? "(HTPv4 format) " : "" ) << "in files : " << endl;
    for( std::size_t j = 0; j < n_pheno; ++j ) {
      ofile_split[j].close();
      sout << "* [" << out_split[j] << "]" << endl;
    }
    sout << endl;
  }

  if(firth || use_SPA) {
    sout << "Number of tests with " << (firth ? "Firth " : "SPA "); 
    sout << "correction : (" << n_corrected << "/" << (snpinfo.size() - snp_tally.n_skipped_snps - snp_tally.n_ignored_snps) * n_pheno << ")" <<  endl;
    sout << "Number of failed tests : (" << snp_tally.n_failed_tests << "/" << n_corrected << ")" << endl;
  }
  sout << "Number of ignored SNPs due to low MAC : " << snp_tally.n_ignored_snps << endl;

} 

// test SNPs in block
void Data::analyze_block(const int &chrom, const int &n_snps, tally* snp_tally, vector<variant_block> &all_snps_info){
  auto t1 = std::chrono::high_resolution_clock::now();

  const int start = snp_tally->snp_count;
  vector< vector < uchar > > snp_data_blocks;
  vector< uint32_t > insize, outsize;

  if(file_type == "bgen"){
    snp_data_blocks.resize( n_snps );
    insize.resize(n_snps); outsize.resize(n_snps);

    // extract genotype data blocks single-threaded
    uint64 first_snp = snpinfo[start].offset;
    ifstream bfile;
    bfile.open( bgen_file, ios::in | ios::binary );
    bfile.seekg( first_snp );
    for(int isnp = 0; isnp < n_snps; isnp++) {
      get_data_blocks(&bfile, &(snp_data_blocks[isnp]));
      insize[isnp] = size1;
      outsize[isnp] = size2;
    }
    bfile.close();
  } else {
    // read genotypes 
    readChunkFromBedFileToG(n_snps, all_snps_info);
  }

  setNbThreads(1);
  // start openmp for loop
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
  for(int isnp = 0; isnp < n_snps; isnp++) {
    int snp_index = start + isnp;
    chi_squared chisq(1);
    double chisq_thr = quantile(complement(chisq, alpha_pvalue));
    double chisq_val, bhat, se_b, pval_log, pval_raw;

    // to store variant information
    variant_block* block_info = &(all_snps_info[isnp]);

    if(file_type == "bgen"){
      // uncompress and extract the dosages 
      // for QTs (or BTs with firth approx.): project out covariates & scale
      extract_variant_MT(&(snp_data_blocks[isnp]), insize[isnp], outsize[isnp], &snpinfo[snp_index], block_info);
    } else {
    // for QTs (or BTs with firth approx): project out covariates & scale
    residualize_geno(block_info);
    }

    // skip SNP if fails filters
    if( block_info->ignored == true ){
      (snp_tally->n_ignored_snps)++;
      continue;
    }
    block_info->pval_log = ArrayXd::Zero(n_pheno); 
    block_info->bhat = ArrayXd::Zero(n_pheno); 
    block_info->se_b = ArrayXd::Zero(n_pheno); 
    block_info->test_fail.assign(n_pheno, false); 
    block_info->is_corrected.assign(n_pheno, true); 

    
    if(binary_mode) {
      MatrixXd tmpG, WX, GW;
      block_info->stats = ArrayXd::Zero(n_pheno); 
      block_info->denum = ArrayXd::Zero(n_pheno); 

      for( std::size_t i = 0; i < n_pheno; ++i ) {

        // project out covariates from G 
        WX = Gamma_sqrt.col(i).asDiagonal() * new_cov;
        // apply mask
        GW = (block_info->Geno * Gamma_sqrt.col(i).array() * masked_indivs.col(i).array()).matrix();
        tmpG = GW - WX * (Xt_Gamma_X_inv[i] * (WX.transpose() * GW));
        block_info->denum(i) = tmpG.squaredNorm(); 

        // score test stat for BT
        block_info->stats(i) = (tmpG.array() * res.col(i).array()).sum() / sqrt( block_info->denum(i) );

        // check this
        if(use_SPA) run_SPA_test_snp(block_info, i, tmpG);
      }

    } else {

      // score test stat for QT
      if( strict_mode ) {
        block_info->Geno *= ind_in_analysis;
        block_info->stats = (res.array().colwise() * block_info->Geno).matrix().transpose().rowwise().sum() / sqrt( ind_in_analysis.sum() );
      } else {
        // compute GtG for each phenotype (different missing patterns)
        block_info->scale_fac_pheno = (masked_indivs.array().colwise() * block_info->Geno).matrix().colwise().squaredNorm();
        block_info->stats = (block_info->Geno.matrix().transpose() * res).transpose().array() / block_info->scale_fac_pheno.sqrt();
      }

    }

    block_info->chisq_val = block_info->stats.square(); 

    for( std::size_t i = 0; i < n_pheno; ++i ) {

      // test statistic & pvalue
      if(!use_SPA) check_pval_snp(block_info, chrom, i);

      // summary stats
      if( !binary_mode ){
        // estimate & SE for QT
        if( strict_mode )
          block_info->bhat(i) = block_info->stats(i) * ( scale_Y(i) * p_sd_yres(i)) / ( sqrt(masked_indivs.col(i).sum()) * block_info->scale_fac ); 
        else
          block_info->bhat(i) = block_info->stats(i) * ( scale_Y(i) * p_sd_yres(i) ) / ( sqrt(block_info->scale_fac_pheno(i)) * block_info->scale_fac ); 
        block_info->se_b(i) = block_info->bhat(i) / block_info->stats(i);
      } else {
        // with Firth, get sum. stats from Firth logistic regression
        if( firth && block_info->is_corrected[i] && !block_info->test_fail[i] ){
          pval_raw = max(nl_dbl_dmin, pow(10, - block_info->pval_log(i))); // to prevent overflow
          block_info->chisq_val(i) = quantile(complement(chisq, pval_raw));
        } else {
          block_info->se_b(i) = 1 / sqrt(block_info->denum(i));
          // with SPA, calculate test stat based on SPA p-value
          if( use_SPA && block_info->is_corrected[i] && !block_info->test_fail[i] ){
            pval_raw = max(nl_dbl_dmin, pow(10, - block_info->pval_log(i))); // to prevent overflow
            block_info->chisq_val(i) = quantile(complement(chisq, pval_raw));
            block_info->bhat(i) = sgn(block_info->stats(i)) * sqrt(block_info->chisq_val(i));
          } else block_info->bhat(i) = block_info->stats(i);
          block_info->bhat(i) *= block_info->se_b(i);
          if( use_SPA && block_info->flipped ) block_info->bhat(i) *= -1;
        }
        block_info->bhat(i) /= block_info->scale_fac;
        block_info->se_b(i) /= block_info->scale_fac;
      }

    }

    //if( (isnp==0) || (isnp == (n_snps-1)) ) cout << "G"<<isnp+1<<" MAF = " <<  block_info.MAF << endl;
  }

  setNbThreads(threads);

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << "done (" << duration.count() << "ms) "<< endl;
}


void Data::get_data_blocks(std::istream* bfile, vector<uchar>* geno_block){

  uint16_t SNPID_size = 0, RSID_size = 0, chromosome_size = 0 , numberOfAlleles = 0 ;
  uint32_t position = 0, allele_size = 0;
  string tmp_buffer;
  // snpid
  bfile->read( reinterpret_cast<char *> (&SNPID_size), 2 );
  tmp_buffer.resize(SNPID_size);
  bfile->read( reinterpret_cast<char *> (&tmp_buffer[0]), SNPID_size );
  // rsid
  bfile->read( reinterpret_cast<char *> (&RSID_size), 2) ;
  tmp_buffer.resize(RSID_size);
  bfile->read( reinterpret_cast<char *> (&tmp_buffer[0]), RSID_size );
  // chromosome
  bfile->read( reinterpret_cast<char *> (&chromosome_size), 2 );
  tmp_buffer.resize(chromosome_size);
  bfile->read( reinterpret_cast<char *> (&tmp_buffer[0]), chromosome_size );
  // position
  bfile->read( reinterpret_cast<char *> (&position), 4 );
  // number of alleles
  bfile->read( reinterpret_cast<char *> (&numberOfAlleles), 2 );
  // alleles
  bfile->read( reinterpret_cast<char *> (&allele_size), 4 );
  tmp_buffer.resize(allele_size);
  bfile->read( reinterpret_cast<char *> (&tmp_buffer[0]), allele_size );
  bfile->read( reinterpret_cast<char *> (&allele_size), 4 );
  tmp_buffer.resize(allele_size);
  bfile->read( reinterpret_cast<char *> (&tmp_buffer[0]), allele_size );

  // set genotype data block
  bfile->read( reinterpret_cast<char *> (&size1), 4 );
  bfile->read( reinterpret_cast<char *> (&size2), 4);
  geno_block->resize(size1 - 4);
  bfile->read( reinterpret_cast<char *> (&((*geno_block)[0])), size1 - 4);

  return;
}

void Data::extract_variant_MT(vector<uchar>* geno_block, const uint32_t insize, const uint32_t outsize, snp* infosnp, variant_block* snp_data){

  uint minploidy = 0, maxploidy = 0, phasing = 0, bits_prob = 0;
  uint16_t numberOfAlleles = 0 ;
  uint32_t nindivs = 0;
  string tmp_buffer;
  uint64 first_snp = infosnp->offset;

  // reset variant info
  snp_data->ignored = false;
  snp_data->fastSPA = use_SPA;
  snp_data->n_non_zero = 0;

  // set genotype data block
  vector < uchar > geno_block_uncompressed;
  geno_block_uncompressed.resize(outsize);

  // uncompress the block using zlib
  uLongf dest_size = outsize;
  if( (uncompress( &(geno_block_uncompressed[0]), &dest_size, &((*geno_block)[0]), insize) != Z_OK) || (dest_size != outsize) ){
    sout << "ERROR: Failed to decompress genotype data block for variant: " << infosnp->ID << endl;
    exit(-1);
  }

  // stream to uncompressed block
  uchar *buffer = &geno_block_uncompressed[0];
  // sample size
  std::memcpy(&nindivs, &(buffer[0]), 4);
  assert( nindivs == n_samples );
  buffer += 4;
  // num alleles
  std::memcpy(&numberOfAlleles, &(buffer[0]), 2);
  assert( numberOfAlleles == 2 );
  buffer += 2;
  // ploidy
  std::memcpy(&minploidy, &(buffer[0]), 1);
  assert( minploidy == 2 );
  buffer ++;
  std::memcpy(&maxploidy, &(buffer[0]), 1);
  assert( maxploidy == 2 );
  buffer ++;

  //to identify missing when getting dosages
  vector < uchar > ploidy_n;
  ploidy_n.resize( nindivs );
  std::memcpy(&(ploidy_n[0]), &(buffer[0]), nindivs);
  buffer += nindivs;

  // phasing
  std::memcpy(&phasing, &(buffer[0]), 1);
  buffer++;

  // bits per probability
  std::memcpy(&bits_prob, &(buffer[0]), 1);
  buffer++;

  // get dosages (can compute mean as going along (and identify non-zero entries if SPA is used)
  bool missing;
  int ns = 0, hc_val;
  double prob0, prob1, total = 0, ds, info_num = 0;
  snp_data->Geno = ArrayXd::Zero(n_samples);
  snp_data->genocounts = MatrixXd::Zero(6, n_pheno);

  // parse genotype probabilities block
  for(size_t i = 0; i < n_samples; i++) {
    missing = ((ploidy_n[i]) & 0x80);
    if(missing) {
      snp_data->Geno(i) = -3;
      buffer+=2;
      continue;
    }
    prob0 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;
    prob1 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;
    snp_data->Geno(i) = 2 - (prob1 + 2 * (std::max( 1 - prob0 - prob1, 0.0) )); // switch allele0 to ALT

    if( ind_in_analysis(i) ){
      if( !strict_mode || (strict_mode && masked_indivs(i,0)) ){
        total += snp_data->Geno(i);
        info_num += 4 * prob0 + prob1 - snp_data->Geno(i) * snp_data->Geno(i);
        ns++;
      }

      // get genotype counts (convert to hardcall)
      if( htp_out ) {
        hc_val = (int) (snp_data->Geno(i) + 0.5); // round to nearest integer 0/1/2
        if( !binary_mode ) {
          snp_data->genocounts.row(hc_val) += masked_indivs.row(i);
        } else {
          snp_data->genocounts.row(hc_val).array() += masked_indivs.row(i).array() * phenotypes_raw.row(i).array();
          snp_data->genocounts.row(3 + hc_val).array() += masked_indivs.row(i).array() * (1 - phenotypes_raw.row(i).array());
        }
      }

    }
  }

  // check MAC
  if( test_mode && ((total < min_MAC) || ((2 * ns - total) < min_MAC)) ) {
    snp_data->ignored = true;
    return;
  }
  
  total /= ns;
  snp_data->af = total / 2;
  snp_data->info = 1 - info_num / (2 * ns * snp_data->af * (1 - snp_data->af));

  if(use_SPA) { 
    // switch to minor allele
    snp_data->flipped = total > 1;
    if( test_type > 0) snp_data->flipped = false; // skip for DOM/REC test
    if(snp_data->flipped){
      snp_data->Geno = ( snp_data->Geno != -3).select( 2 - snp_data->Geno, snp_data->Geno);
      total = 2 - total;
    }
  }

  // apply dominant/recessive encoding & recompute mean
  if(test_type > 0){
    // go over data block again
    buffer -= 2 * n_samples;
    for(size_t i = 0; i < n_samples; i++) {
      missing = ((ploidy_n[i]) & 0x80);
      if(missing) {
        buffer+=2;
        continue;
      }
      prob0 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;
      prob1 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;

      if(ind_in_analysis(i)){
        if(test_type == 1){ //dominant
          snp_data->Geno(i) = prob0 + prob1; // allele0 is ALT
        } else if(test_type == 2){ //recessive
          snp_data->Geno(i) = prob0; // allele0 is ALT
        }
      }
    }
    total = ((snp_data->Geno != -3) && (ind_in_analysis == 1)).select(snp_data->Geno, 0).sum() / ns;
    if(total < numtol) {
      snp_data->ignored = true;
      return;
    }
  }

  // deal with missing data & prep for spa
  for( size_t i = 0; i < n_samples; ++i ) {
    ds = snp_data->Geno(i);

    // keep track of number of entries filled so avoid using clear
    if( use_SPA && (snp_data->fastSPA) && ind_in_analysis(i) && ds > 0 ) {
      snp_data->n_non_zero++;
      if(snp_data->n_non_zero > (n_samples / 2)) {
        snp_data->fastSPA = false;
      } else {
        if(snp_data->non_zero_indices.size() < snp_data->n_non_zero)
          snp_data->non_zero_indices.push_back(i);
        else
          snp_data->non_zero_indices[snp_data->n_non_zero - 1] = i;
      }
    }

    if(ds != -3 && ind_in_analysis(i) && (!strict_mode || (strict_mode && masked_indivs(i,0)) ) ){
      snp_data->Geno(i) -= total;
    } else {
      snp_data->Geno(i) = 0;
    }
  }

  if(!binary_mode || firth_approx){
    // project out covariates
    MatrixXd beta = new_cov.transpose() * snp_data->Geno.matrix();
    snp_data->Geno -= (new_cov * beta).array();

    // scale
    snp_data->scale_fac = snp_data->Geno.matrix().norm();
    snp_data->scale_fac /= sqrt( (strict_mode ? Neff(0) : ind_in_analysis.sum()) - 1 );

    if( snp_data->scale_fac < numtol ) {
      snp_data->ignored = true;
      return;
    }
    snp_data->Geno /= snp_data->scale_fac;

  } else snp_data->scale_fac = 1;

  return;
}

void Data::readChunkFromBedFileToG(const int &bs, vector<variant_block> &all_snps_info){

  int hc, ns, byte_start, bit_start;
  double total; 
  // mapping matches the switch of alleles done when reading bim
  const int maptogeno[4] = {2, -3, 1, 0};

  for(size_t j = 0; j < bs; j++) {
    variant_block* snp_data = &(all_snps_info[j]);
    snp_data->Geno = ArrayXd::Zero(n_samples);
    snp_data->genocounts = MatrixXd::Zero(6, n_pheno);
    // reset variant info
    snp_data->ignored = false;
    snp_data->fastSPA = use_SPA;
    snp_data->n_non_zero = 0;

    ns = 0, total = 0;
    bed_ifstream.read( reinterpret_cast<char *> (&inbed[0]), bed_block_size);
    for (size_t i = 0; i < n_samples; i++) {
      byte_start = i>>2; // 4 samples per byte
      bit_start = (i&3)<<1; // 2 bits per sample
      hc = maptogeno[ (inbed[byte_start] >> bit_start)&3 ]; 
      snp_data->Geno(i) = hc;

      if(hc != -3) {
        if( ind_in_analysis(i) ){
          if( !strict_mode || (strict_mode && masked_indivs(i,0)) ){
            total += hc;
            ns++;
          }
        }

        // get genotype counts 
        if( htp_out ) {
          if( !binary_mode ) {
            snp_data->genocounts.row(hc) += masked_indivs.row(i);
          } else {
            snp_data->genocounts.row(hc).array() += masked_indivs.row(i).array() * phenotypes_raw.row(i).array();
            snp_data->genocounts.row(3 + hc).array() += masked_indivs.row(i).array() * (1 - phenotypes_raw.row(i).array());
          }
        }

      }
    }

    // check MAC
    if( test_mode && ((total < min_MAC) || (( 2 * ns - total) < min_MAC)) ){
      snp_data->ignored = true;
      continue;
    }

    total /= ns;
    snp_data->af = total / 2;

    if(use_SPA) { 
      // switch to minor allele
      snp_data->flipped = total > 1;
      if( test_type > 0) snp_data->flipped = false; // skip for DOM/REC test
      if(snp_data->flipped){
        snp_data->Geno = ( snp_data->Geno != -3).select( 2 - snp_data->Geno, snp_data->Geno);
        total = 2 - total;
      }
    }

    // apply dominant/recessive encoding & recompute mean
    if(test_type > 0){
      if(test_type == 1){ //dominant
        snp_data->Geno = (snp_data->Geno == 2).select(1, snp_data->Geno);
      } else if(test_type == 2){ //recessive
        snp_data->Geno = (snp_data->Geno >= 1).select(snp_data->Geno - 1, snp_data->Geno);
      }
      total = ((snp_data->Geno != -3) && (ind_in_analysis == 1)).select(snp_data->Geno, 0).sum() / ns;
      if(total < numtol) {
        snp_data->ignored = true;
        continue;
      }
    }


    // deal with missing data & prep for spa
    for( size_t i = 0; i < n_samples; ++i ) {
      hc = snp_data->Geno(i);

      // keep track of number of entries filled so avoid using clear
      if( use_SPA && (snp_data->fastSPA) && ind_in_analysis(i) && hc > 0 ) {
        snp_data->n_non_zero++;
        if(snp_data->n_non_zero > (n_samples / 2)) {
          snp_data->fastSPA = false;
        } else {
          if(snp_data->non_zero_indices.size() < snp_data->n_non_zero)
            snp_data->non_zero_indices.push_back(i);
          else
            snp_data->non_zero_indices[snp_data->n_non_zero - 1] = i;
        }
      }

      if(hc != -3 && ind_in_analysis(i) && (!strict_mode || (strict_mode && masked_indivs(i,0)) ) ){
        snp_data->Geno(i) -= total;
      } else {
        snp_data->Geno(i) = 0;
      }
    }


  }

}

void Data::residualize_geno(variant_block* snp_data){

  if(!binary_mode || firth_approx){
    // project out covariates
    MatrixXd beta = new_cov.transpose() * snp_data->Geno.matrix();
    snp_data->Geno -= (new_cov * beta).array();

    // scale
    snp_data->scale_fac = snp_data->Geno.matrix().norm();
    snp_data->scale_fac /= sqrt( (strict_mode ? Neff(0) : ind_in_analysis.sum()) - 1 );

    if( snp_data->scale_fac < numtol ) {
      snp_data->ignored = true;
      return;
    }
    snp_data->Geno /= snp_data->scale_fac;

  } else snp_data->scale_fac = 1;

}


void Data::run_SPA_test_snp(variant_block* block_info, int ph, const VectorXd& Gtmp){

  int index_j, nnz;
  double pval, logp, zstat, zstat_sq, score_num, tval, limK1_low, limK1_high, root_K1, sumG;
  chi_squared chisq(1);
  double chisq_thr = quantile(chisq, 1 - alpha_pvalue); 
  ArrayXd Gmu;

  // check test statistic relative to threshold
  zstat = block_info->stats(ph);
  zstat_sq = zstat * zstat;
  if( zstat_sq <= chisq_thr){
    pval = cdf(complement(chisq, zstat_sq));
    if(pval == 0) logp = log10(2) - 0.5 * log10( 2 * M_PI * zstat_sq ) - 0.5 * zstat_sq * M_LOG10E ; 
    else logp = log10(pval);
    logp *= -1;
    block_info->pval_log(ph) = logp;
    block_info->is_corrected[ph] = false;
    return;
  }

  // compute needed quantities
  block_info->val_c = sqrt( block_info->denum(ph) );  // sqrt( G'WG )
  score_num = zstat * block_info->val_c;
  block_info->Gmod = Gtmp.array() / Gamma_sqrt.col(ph).array() * masked_indivs.col(ph).array();
  Gmu = block_info->Gmod * Y_hat_p.col(ph).array();
  block_info->val_a = Gmu.sum(); 
  sumG = block_info->Gmod.sum();

  if(block_info->fastSPA){
    block_info->val_b = block_info->denum(ph);
    block_info->val_d = 0;
    for( std::size_t j = 0; j < block_info->n_non_zero; ++j ) {
      index_j = block_info->non_zero_indices[j];
      if(!masked_indivs(index_j,ph)) continue;
      block_info->val_b -= Gtmp(index_j) * Gtmp(index_j);
      block_info->val_d += Gmu(index_j);
    }
  }

  // check if K'(t)= s can be solved 
  limK1_low = (block_info->Gmod < 0).select(block_info->Gmod, 0 ).sum() - block_info->val_a ;
  limK1_high = (block_info->Gmod > 0).select(block_info->Gmod, 0 ).sum() - block_info->val_a ;
  if( score_num < limK1_low || score_num > limK1_high ){
    if(verbose) sout << "WARNING: SPA failed (solution to K'(t)=s is infinite)";
    block_info->test_fail[ph] = true;
    return;
  }

  // keep track of whether obs stat is positive
  block_info->pos_score = zstat  > 0;
  tval = fabs(zstat);

  // solve K'(t)= tval using a mix of Newton-Raphson and bisection method
  root_K1 = solve_K1_snp(tval, block_info, ph);
  if( root_K1 == missing_value_double ){
    block_info->test_fail[ph] = true;
    return;
  }

  // compute pvalue
  block_info->pval_log(ph) = get_SPA_pvalue_snp(root_K1, tval, block_info, ph);

}

double Data::solve_K1_snp(double tval, variant_block* block_info, int ph){

  int niter_cur;
  int lambda = block_info->pos_score ? 1 : -1; // if score is negative, adjust K' and K''
  double min_x, max_x, t_old, f_old, t_new = -1, f_new, hess;

  niter_cur = 0;
  min_x = 0, max_x = std::numeric_limits<double>::infinity();
  t_old = 0;
  f_old = block_info->fastSPA ? compute_K1_fast_snp(lambda * t_old, block_info, ph) : compute_K1_snp(lambda * t_old, block_info, ph);
  f_old *= lambda;
  f_old -= tval; 

  while( niter_cur++ < niter_max_spa ){

    hess = block_info->fastSPA ? compute_K2_fast_snp(lambda * t_old, block_info, ph) : compute_K2_snp(lambda * t_old, block_info, ph);
    t_new = t_old - f_old / hess;
    f_new = block_info->fastSPA ? compute_K1_fast_snp(lambda * t_new, block_info, ph) : compute_K1_snp(lambda * t_new, block_info, ph);
    f_new *= lambda;
    f_new -= tval;

    if( fabs( f_new ) < tol_spa ) break;

    // update bounds on root
    if( t_new && (t_new > min_x) && (t_new < max_x) ){
      if( f_new > 0) max_x = t_new;
      else min_x = t_new;
    } else{ // bisection method if t_new went out of bounds and re-compute f_new
      t_new = ( min_x + max_x ) / 2;
      // if( fabs( min_x - t_new ) < tol_spa ) break;
      f_new = block_info->fastSPA ? compute_K1_fast_snp(lambda * t_new, block_info, ph) : compute_K1_snp(lambda * t_new, block_info, ph);
      f_new *= lambda;
      f_new -= tval;
    }

    t_old = t_new;
    f_old = f_new;
  }

  // If didn't converge
  if( niter_cur > niter_max_spa ){
    if(verbose) sout << "WARNING: SPA did not converge to root for K'(t)=s.\n";
    return missing_value_double;
  }
    //sout << "#iterations = " << niter_cur << "; f= " << f_new << endl;

  return t_new;
}

double Data::compute_K_snp(double t, variant_block* block_info, int ph){

  double val = (1 - Y_hat_p.col(ph).array() + Y_hat_p.col(ph).array() * ( t / block_info->val_c * block_info->Gmod ).exp() ).log().sum() - t * block_info->val_a / block_info->val_c;

  return val;
}

double Data::compute_K_fast_snp(double t, variant_block* block_info, int ph){

  uint64 index_j;
  double val = 0;

  for( std::size_t j = 0; j < block_info->n_non_zero; ++j ) {
    index_j = block_info->non_zero_indices[j];
    val += log( 1 - Y_hat_p(index_j,ph) + Y_hat_p(index_j,ph) * exp( t / block_info->val_c * block_info->Gmod(index_j)) );
  }
  val += -t * block_info->val_d / block_info->val_c + t * t / 2 / block_info->denum(ph) * block_info->val_b;

  return val;
}

double Data::compute_K1_snp(double t, variant_block* block_info, int ph){

  double val = ( ( block_info->Gmod * Y_hat_p.col(ph).array() / block_info->val_c ) / ( Y_hat_p.col(ph).array() + (1 - Y_hat_p.col(ph).array()) * ( -t / block_info->val_c * block_info->Gmod).exp() ) ).sum();
 val -= block_info->val_a / block_info->val_c;

  return val;
}

double Data::compute_K1_fast_snp(double t, variant_block* block_info,  int ph){

  uint64 index_j;
  double val = 0;

  for( std::size_t j = 0; j < block_info->n_non_zero; ++j ) {
    index_j = block_info->non_zero_indices[j];
    val += ( block_info->Gmod(index_j) * Y_hat_p(index_j,ph) / block_info->val_c ) / ( Y_hat_p(index_j,ph) + (1 - Y_hat_p(index_j,ph)) * exp( -t / block_info->val_c * block_info->Gmod(index_j)) );
  }
  val += -block_info->val_d / block_info->val_c + t / block_info->denum(ph) * block_info->val_b;

  return val;
}

double Data::compute_K2_snp(double t, variant_block* block_info, int ph){

  double val = ( ( block_info->Gmod.square() * Gamma_sqrt.col(ph).array().square() / (block_info->val_c*block_info->val_c) * ( -t / block_info->val_c * block_info->Gmod).exp()) / ( Y_hat_p.col(ph).array() + (1 - Y_hat_p.col(ph).array()) * ( -t / block_info->val_c * block_info->Gmod).exp() ).square() ).sum();

  return val;
}
  
double Data::compute_K2_fast_snp(double t, variant_block* block_info, int ph){

  uint64 index_j;
  double val = 0, denum;

  for( std::size_t j = 0; j < block_info->n_non_zero; ++j ) {
    index_j = block_info->non_zero_indices[j];
    denum = Y_hat_p(index_j,ph) + (1 - Y_hat_p(index_j,ph)) * exp( -t / block_info->val_c * block_info->Gmod(index_j));
    val += ( block_info->Gmod(index_j) * block_info->Gmod(index_j) * Gamma_sqrt(index_j,ph) * Gamma_sqrt(index_j,ph) * exp( -t / block_info->val_c * block_info->Gmod(index_j)) / block_info->val_c / block_info->val_c ) / (denum * denum);
  }
  val += block_info->val_b / block_info->denum(ph);

  return val;
}

double Data::get_SPA_pvalue_snp(double root, double tval, variant_block* block_info, int ph){

  int lambda = block_info->pos_score ? 1 : -1; // if score is negative, adjust K and K''
  double kval, k2val, wval, vval, rval, pvalue;
  normal nd(0,1);

  kval = block_info->fastSPA ? compute_K_fast_snp(lambda * root, block_info, ph) : compute_K_snp(lambda * root, block_info, ph);
  k2val = block_info->fastSPA ? compute_K2_fast_snp(lambda * root, block_info, ph) : compute_K2_snp(lambda * root, block_info, ph);

  wval = sqrt( 2 * ( root * tval - kval ) );
  vval = root * sqrt( k2val );
  if(vval == 0) return 0; // implies score = 0 so pval = 1

  rval = wval + log( vval / wval ) / wval;
  pvalue = cdf(complement(nd, rval ));
  pvalue *= 2;

  if(pvalue > 1) { // SPA can fail for SNPs with very low counts 
    block_info->test_fail[ph] = true;
    if(verbose) sout << "WARNING: SPA correction failed (resulted in p-value > 1).\n";
    return missing_value_double;
    //return 0;  
  }

  return -log10( pvalue );
}
 

void Data::check_pval_snp(variant_block* block_info, int chrom, int ph){

  chi_squared chisq(1);
  double chisq_thr = quantile(chisq, 1 - alpha_pvalue); 
  double pval, logp, LRstat;
  double tstat = block_info->chisq_val(ph);

  // if quantitative trait, or firth isn't used, or Tstat < threshold, no correction done
  if(!binary_mode || !firth || tstat <= chisq_thr){
    pval = cdf(complement(chisq, tstat));

    // check if pvalue from boost 0 (when tstat is very large)
    if(pval == 0) logp = log10(2) - 0.5 * log10( 2 * M_PI * tstat ) - 0.5 * tstat * M_LOG10E ; 
    else logp = log10(pval);

    logp *= -1;
    block_info->pval_log(ph) = logp;
    block_info->is_corrected[ph] = false;
    return;
  }

  // firth
  run_firth_correction_snp(chrom, block_info, ph);
  if(block_info->test_fail[ph]) return;

  LRstat = block_info->dif_deviance;
  pval = cdf(complement(chisq, LRstat));
  if(pval == 0) logp = log10(2) - 0.5 * log10( 2 * M_PI * LRstat ) - 0.5 * LRstat * M_LOG10E ; 
  else logp = log10(pval);
  logp *= -1;

  block_info->pval_log(ph) = logp;
  return;
}


void Data::run_firth_correction_snp(int chrom, variant_block* block_info, int ph){

  // obtain null deviance (set SNP effect to 0 and compute max. pen. LL)
  if(!firth_approx){
    fit_firth_logistic_snp(chrom, ph, true, block_info);
    if(block_info->test_fail[ph]) return ;
  }

  // fit full model and compute deviance
  fit_firth_logistic_snp(chrom, ph, false, block_info);

  return ;
}


void Data::fit_firth_logistic_snp(int chrom, int ph, bool null_fit, variant_block* block_info) {
  // if firth is used, fit based on penalized log-likelihood

  int niter_cur, niter_search, col_incl;
  double dl, dev_old, denum, mx, deviance_l0 = 0;
  int maxstep_firth = null_fit ? maxstep_null : maxstep;
  int niter_firth = null_fit ? niter_max_firth_null : niter_max_firth;
  ArrayXd mod_score, betaold, betanew, step_size, etavec, pivec, wvec, hvec;
  MatrixXd Xmat, XtW, XtWX;
  ColPivHouseholderQR<MatrixXd> qr, qrX;

  ArrayXd Y1 = phenotypes_raw.col(ph).array() * masked_indivs.col(ph).array();

  if(firth_approx){
    if(null_fit){
      Xmat = new_cov; // only covariates
    } else {
      Xmat = block_info->Geno.matrix(); // only tested SNP
    }
    col_incl = Xmat.cols();
  } else {
    Xmat.resize(n_samples, new_cov.cols() + 1); // covariates + tested SNP
    Xmat << new_cov, block_info->Geno.matrix();
    col_incl = Xmat.cols();
    if( null_fit ) col_incl--;
  }
  // mask individuals
  Xmat.array().colwise() *= masked_indivs.col(ph).array();

  // starting values
  if(null_fit){

    // start at logit^-1(mean(Y))-mean(offset)
    betaold = ArrayXd::Zero(Xmat.cols()); // last entry in exact Firth is kept at 0
    betaold(0) = ( 0.5 + Y1.sum())  / (masked_indivs.col(ph).sum() + 1);
    betaold(0) = log( betaold(0) / (1 - betaold(0) ));

    // LOCO prediction is offset
    betaold(0) -= (blups.col(ph).array() * masked_indivs.col(ph).array()).mean();

  } else {
    // start at 0
    if(firth_approx) betaold = ArrayXd::Zero(Xmat.cols()); 
    // start at estimate from null fit
    else betaold = block_info->beta_null_firth.col(0);
  }
  betanew = ArrayXd::Zero( betaold.size() );

  // get the corresponding deviance
  etavec = (Xmat * betaold.matrix()).array();
  etavec += blups.col(ph).array() * masked_indivs.col(ph).array();
  // covariate effects added as offset in firth approx.
  if( firth_approx && !null_fit ) etavec += ((new_cov.array().colwise() * masked_indivs.col(ph).array()).matrix() * beta_null_firth.block(0,ph,new_cov.cols(),1)).array(); 
  // fitted probabilities
  pivec = 1 - 1 / (etavec.exp() + 1) ;
  wvec = (masked_indivs.col(ph).array() == 1).select( ( pivec * (1 - pivec) ).sqrt(), 0);
  XtW = Xmat.transpose() * wvec.matrix().asDiagonal();
  XtWX = XtW * XtW.transpose();
  qr.compute(XtWX.block(0, 0, col_incl, col_incl));
  // use penalized log-lik
  dev_old = (masked_indivs.col(ph).array() == 1).select( (Y1 * pivec.log() + (1-Y1) * (1-pivec).log() ), 0).sum();
  dev_old += 0.5 * qr.logAbsDeterminant();
  dev_old *= -2;

  // at niter=0 (i.e. betaSNP=0), this is null deviance
  if( !null_fit ) deviance_l0 = dev_old;

  // solve S'(beta) = S(beta) + X'(h*(0.5-p)) = 0
  niter_cur = 0;
  while(niter_cur++ < niter_firth){

    ////////  compute step size
    etavec = (Xmat * betaold.matrix()).array();
    etavec += blups.col(ph).array() * masked_indivs.col(ph).array();
    if( firth_approx && !null_fit ) etavec += ((new_cov.array().colwise() * masked_indivs.col(ph).array()).matrix() * beta_null_firth.block(0,ph,new_cov.cols(),1)).array(); 

    pivec = 1 - 1 / (etavec.exp() + 1) ;
    wvec = (masked_indivs.col(ph).array() == 1).select( ( pivec * (1 - pivec) ).sqrt(), 0);
    XtW = Xmat.transpose() * wvec.matrix().asDiagonal();
    XtWX = XtW * XtW.transpose();
    qr.compute(XtWX);
    if( !firth_approx && null_fit ) qrX.compute(XtWX.block(0, 0, col_incl, col_incl));
    // compute diag(H), H = U(U'U)^{-1}U', U = Gamma^(1/2)X
    hvec = (qr.solve(XtW).array() * XtW.array() ).colwise().sum();

    // modified score for beta
    mod_score = (Xmat.leftCols(col_incl).transpose() * (masked_indivs.col(ph).array() == 1).select( Y1 - pivec + hvec * (0.5 - pivec), 0).matrix() ).array();

    // step size
    if( !firth_approx && null_fit )
      step_size = qrX.solve( mod_score.matrix() ).array();
    else
      step_size = qr.solve( mod_score.matrix() ).array();

    // force absolute step size to be less than maxstep for each entry of beta
    mx = step_size.abs().maxCoeff() / maxstep_firth;
    if( mx > 1 ) step_size /= mx;

    // start step-halving and stop when deviance decreases 
    denum = 1;
    for( size_t niter_search = 1; niter_search <= niter_max_line_search; niter_search++ ){

      // adjusted step size
      step_size /= denum;

      ///////// compute corresponding deviance
      betanew.head(col_incl) = betaold.head(col_incl) + step_size;
      etavec = (Xmat * betanew.matrix()).array();
      etavec += blups.col(ph).array() * masked_indivs.col(ph).array();
      if( firth_approx && !null_fit ) etavec += ((new_cov.array().colwise() * masked_indivs.col(ph).array()).matrix() * beta_null_firth.block(0,ph,new_cov.cols(),1)).array(); 

      pivec = 1 - 1 / (etavec.exp() + 1) ;
      wvec = (masked_indivs.col(ph).array() == 1).select( ( pivec * (1 - pivec) ).sqrt(), 0);
      XtW = Xmat.transpose() * wvec.matrix().asDiagonal();
      XtWX = XtW * XtW.transpose();
      qr.compute(XtWX);

      dl = (masked_indivs.col(ph).array() == 1).select( (Y1 * pivec.log() + (1-Y1) * (1-pivec).log() ), 0).sum();
      dl += 0.5 * qr.logAbsDeterminant();
      dl *= -2;

      //sout << "\n["<<niter_cur << " - " << niter_search <<"]  denum =" << denum << ";\n step =" << step_size.matrix().transpose().array() / denum<<"; \nbeta=" << betanew.matrix().transpose().array() << ";\n Lnew= " << dl << " vs L0="<< dev_old << ";score="<< mod_score<< endl;
      if( dl < dev_old + numtol ) break;
      denum *= 2;
    }

    betaold.head(col_incl) += step_size;
    dev_old = dl;

    // stopping criterion using modified score function
    if( mod_score.abs().maxCoeff() < numtol_firth) break;

  }

  // If didn't converge
  if(niter_cur > niter_firth){
    if(verbose) sout << "WARNING: Logistic regression with Firth correction did not converge!\n";
    block_info->test_fail[ph] = true;
    return ;
  }
  // sout << "\nNiter = " << niter_cur << " : " << mod_score.matrix().transpose() << endl;


  if(null_fit) {
    if(!firth_approx) block_info->beta_null_firth = betaold.matrix();
  } else {
    // compute beta_hat & SE
    block_info->bhat(ph) = betaold.tail(1)(0);
    block_info->se_b(ph) = sqrt( qr.inverse().diagonal().tail(1)(0) );

    dl -= deviance_l0;
    dl *= -1;
    if( dl < 0 ) {
      block_info->test_fail[ph] = true;
      return ;
    }
    block_info->dif_deviance = dl;
  }

  return ;
}
