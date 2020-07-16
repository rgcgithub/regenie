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

#include "Regenie.hpp"
#include "Geno.hpp"
#include "Step1_Models.hpp"
#include "Step2_Models.hpp"
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


void print_help( bool help_full ){

  print_header(cout);
  cout << "Options:\n";

  cout << left << std::setw(35) << " --help" << "print list of available options\n";
  cout << left << std::setw(35) << " --helpFull" << "print list of all available options\n";

  cout << left << std::setw(35) << " --step INT"<< "specify if fitting null model (=1) or association testing (=2)\n";
  cout << left << std::setw(35) << " --bed PREFIX" << "prefix to PLINK .bed/.bim/.fam files\n";
  cout << left << std::setw(35) << " --pgen PREFIX" << "prefix to PLINK2 .pgen/.pvar/.psam files\n";
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
  cout << left << std::setw(35) << " --phenoColList STRING"<< "Comma separated list of phenotype names to keep\n";
  cout << left << std::setw(35) << " --c FILE" << "covariate file (header required starting with FID IID)\n";
  cout << left << std::setw(35) << " --covarCol STRING"<< "covariate name in header (use for each covariate to keep)\n";
  cout << left << std::setw(35) << " --covarColList STRING"<< "Comma separated list of covariate names to keep\n";
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
  cout << left << std::setw(35) << " --chrList STRING"<< "Comma separated list of chromosomes to test in step 2\n";
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

void print_header(std::ostream& o){

  o << left << std::setw(14) << " " << "|==================================|" << endl;
  o << left << std::setw(14) << " " << "|           REGENIE v" << left << std::setw(14) << VERSION_NUMBER << "|" << endl;
  o << left << std::setw(14) << " " << "|==================================|" << endl << endl;

  o << "Copyright (c) 2020 Joelle Mbatchou and Jonathan Marchini." << endl;
  o << "Distributed under the MIT License.\n\n";
}


void read_params_and_check(int argc, char *argv[], struct param* params, struct in_files* files, struct filter* filters, MeasureTime* mt, mstream& sout) {

  int maxargs = argc - 1;
  vector< string > tmp_str_vec;
  std::string trimmed_str, trimmed_str_arg;

  for(size_t counter=1; counter<argc; counter++){	  
    trimmed_str = string(argv[counter]);  // trim this
    trimmed_str.erase(std::remove_if(trimmed_str.begin(), trimmed_str.end(), ::isspace), trimmed_str.end());

    if(trimmed_str == "--bt") params->binary_mode = true;
    if(trimmed_str == "--1") params->CC_ZeroOne = false;
    if(trimmed_str == "--within") params->within_sample_l0 = true;
    if(trimmed_str == "--loocv") params->use_loocv = true;
    if(trimmed_str == "--split") params->split_by_pheno = true;
    if(trimmed_str == "--strict") params->strict_mode = true;
    if(trimmed_str == "--force-impute") params->rm_missing_qt = false;
    if(trimmed_str == "--ignore-pred") params->skip_blups = true;
    if(trimmed_str == "--print") params->print_block_betas = true;
    if(trimmed_str == "--approx") params->firth_approx = true;
    if(trimmed_str == "--v") params->verbose = true;
    if(trimmed_str == "--nostream") params->streamBGEN = false;
    if(trimmed_str == "--helpFull") print_help( true );
    if(trimmed_str == "--help") print_help( false );

    if((trimmed_str == "--bgen") && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        files->bgen_file = trimmed_str_arg;
        params->file_type = "bgen";
        params->n_genofiles++;
      }
    }

    if((trimmed_str == "--pgen")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        files->pgen_prefix = trimmed_str_arg;
        params->file_type = "pgen";
        params->n_genofiles++;
      }
    }

    if((trimmed_str == "--bed")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        files->bedfile = trimmed_str_arg + ".bed";
        files->bimfile = trimmed_str_arg + ".bim";
        files->famfile = trimmed_str_arg + ".fam";
        params->file_type = "bed";
        params->n_genofiles++;
      }
    }

    if((trimmed_str == "--sample")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        files->sample_file = trimmed_str_arg;
        params->bgenSample = true;
      }
    }

    if((trimmed_str == "--remove")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        files->file_ind_exclude = trimmed_str_arg;
        params->rm_indivs = true;
      }
    }

    if((trimmed_str == "--keep")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        files->file_ind_include = trimmed_str_arg;
        params->keep_indivs = true;
      }
    }

    if((trimmed_str == "--extract")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        files->file_snps_include = trimmed_str_arg;
        params->keep_snps = true;
      }
    }

    if((trimmed_str == "--exclude")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        files->file_snps_exclude = trimmed_str_arg;
        params->rm_snps = true;
      }
    }

    if((trimmed_str == "--phenoCol")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        params->select_phenos = true;
        filters->pheno_colKeep_names.push_back( trimmed_str_arg );
      }
    }

    if((trimmed_str == "--phenoColList")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        params->select_phenos = true;
        // split csv
        boost::algorithm::split(tmp_str_vec, trimmed_str_arg, is_any_of(","));
        filters->pheno_colKeep_names.insert( 
            filters->pheno_colKeep_names.end(),
            std::begin( tmp_str_vec ), 
            std::end( tmp_str_vec)         );
      }
    }

    if((trimmed_str == "--covarCol")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        params->select_covs = true;
        filters->cov_colKeep_names.push_back( trimmed_str_arg );
      }
    }

    if((trimmed_str == "--covarColList")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        params->select_covs = true;
        // split csv
        boost::algorithm::split(tmp_str_vec, trimmed_str_arg, is_any_of(","));
        filters->cov_colKeep_names.insert( 
            filters->cov_colKeep_names.end(),
            std::begin( tmp_str_vec ), 
            std::end( tmp_str_vec )      );
      }
    }

    if((trimmed_str == "--c")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') files->cov_file = trimmed_str_arg;
    }

    if((trimmed_str == "--p")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') files->pheno_file = trimmed_str_arg;
    }

    if((trimmed_str == "--pred")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') files->blup_file = trimmed_str_arg;
    }

    if((trimmed_str == "--o")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') files->out_file = trimmed_str_arg;
    }

    if(trimmed_str == "--lowmem"){
      params->write_l0_pred = true;
      if(counter < maxargs){ 
        trimmed_str_arg = string(argv[counter+1]);  // trim this
        trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
        if(trimmed_str_arg[0] != '-') files->loco_tmp_prefix = trimmed_str_arg;
      }
    }

    if((trimmed_str == "--b")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') params->block_size = atoi(argv[counter+1]);
    }

    if((trimmed_str == "--nb")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') params->n_block = atoi(argv[counter+1]);
    }

    if((trimmed_str == "--cv")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') params->cv_folds = atoi(argv[counter+1]);
    }

    if((trimmed_str == "--l0")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') params->n_ridge_l0 = atoi(argv[counter+1]);
    }

    if((trimmed_str == "--l1")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') params->n_ridge_l1 = atoi(argv[counter+1]);
    }

    if((trimmed_str == "--step")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') params->run_mode = atoi(argv[counter+1]);
    }

    if((trimmed_str == "--nauto")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') params->nChrom = atoi(argv[counter+1]) + 1;
    }

    if((trimmed_str == "--niter")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') params->niter_max = atoi(argv[counter+1]);
    }

    if((trimmed_str == "--maxstep-null")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        params->maxstep_null = atoi(argv[counter+1]);
        params->fix_maxstep_null = true;
      }
    }

    if((trimmed_str == "--maxiter-null")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        params->niter_max_firth_null = atoi(argv[counter+1]);
        params->fix_maxstep_null = true;
      }
    }

    if((trimmed_str == "--minMAC")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        params->setMinMAC = true;
        params->min_MAC = atoi(argv[counter+1]);
      }
    }

    if((trimmed_str == "--chr")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        params->select_chrs = true;
        filters->chrKeep_test.push_back( chrStrToInt(trimmed_str_arg, params->nChrom) );
      }
    }

    if((trimmed_str == "--chrList")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        params->select_chrs = true;
        // split csv
        boost::algorithm::split(tmp_str_vec, trimmed_str_arg, is_any_of(","));
        for( size_t ichr = 0; ichr < tmp_str_vec.size(); ichr++)
          filters->chrKeep_test.push_back( chrStrToInt(tmp_str_vec[ichr], params->nChrom) );
      }
    }

    if((trimmed_str == "--test")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') {
        if( trimmed_str_arg == "dominant") params->test_type = 1; 
        else if( trimmed_str_arg == "recessive") params->test_type = 2; 
        else {
          print_header(cerr);
          cerr << "ERROR : Unrecognized argument, must be either 'dominant' or 'recessive' for option `--test`.\n" << params->err_help ;
          exit(-1);
        }
      }
    }

    if((trimmed_str == "--setl0")  && (counter < maxargs)){
      params->n_ridge_l0 = 0;
      // get number of ridge parameters
      for(size_t i = 1; (counter+i) < argc; i++) {
        trimmed_str_arg = string(argv[counter+i]);  // trim this
        trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
        if(trimmed_str_arg[0] != '-') params->n_ridge_l0++;
        else break;
      }
      // read them in
      params->lambda.resize(params->n_ridge_l0);
      for(size_t i = 1; i <= params->n_ridge_l0; i++) params->lambda[i] = atof(argv[counter+i]);
      params->user_ridge_params_l0 = true;
    }

    if((trimmed_str == "--setl1")  && (counter < maxargs)){
      params->n_ridge_l1 = 0;
      // get number of ridge parameters
      for(size_t i = 1; (counter+i) < argc; i++) {
        trimmed_str_arg = string(argv[counter+1]);  // trim this
        trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
        if(trimmed_str_arg[0] != '-') params->n_ridge_l1++;
        else break;
      }
      // read them in
      params->tau.resize(params->n_ridge_l1);
      for(size_t i = 1; i <= params->n_ridge_l1; i++) params->tau[i] = atof(argv[counter+i]);
      params->user_ridge_params_l1 = true;
    }

    if(trimmed_str == "--firth"){
      params->firth = true;
      if(counter < maxargs){ 
        trimmed_str_arg = string(argv[counter+1]);  // trim this
        trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
        if(trimmed_str_arg[0] != '-') params->alpha_pvalue = atof(argv[counter+1]);;
      }
    }

    if(trimmed_str == "--spa"){
      params->use_SPA = true;
      if(counter < maxargs){ 
        trimmed_str_arg = string(argv[counter+1]);  // trim this
        trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
        if(trimmed_str_arg[0] != '-') params->alpha_pvalue = atof(argv[counter+1]);
      }
    }

    if((trimmed_str == "--threads")  && (counter < maxargs)){
      trimmed_str_arg = string(argv[counter+1]);  // trim this
      trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
      if(trimmed_str_arg[0] != '-') params->threads = atoi(argv[counter+1]);
    }

    if(trimmed_str == "--htp"){
      params->htp_out = params->split_by_pheno = true;
      if(counter < maxargs){ 
        trimmed_str_arg = string(argv[counter+1]);  // trim this
        trimmed_str_arg.erase(std::remove_if(trimmed_str_arg.begin(), trimmed_str_arg.end(), ::isspace), trimmed_str_arg.end());
        if(trimmed_str_arg[0] != '-') params->cohort_name = trimmed_str_arg;
      }
    }

  }

  if(files->out_file == "NULL") {
    print_header(cerr);
    cerr << "ERROR :You must specify an output file with --o.\n" << params->err_help ;
    exit(-1);
  }

  // Print output to file and to stdout
  // print command line arguments
  start_log(argc, argv, files->out_file, mt, sout);


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

    // user specified ridge parameters to use
    if( params->user_ridge_params_l0 ){
      if(params->n_ridge_l0 < 1){
        sout << "ERROR : Number of ridge parameters must be at least 1 in option --setl0.\n" << params->err_help;
        exit(-1);
      }
      // parameters must be less in (0, 1)
      if( std::count_if(params->lambda.begin(), params->lambda.end(), std::bind2nd(std::greater<double>(), 0)) != params->n_ridge_l0 || std::count_if(params->lambda.begin(), params->lambda.end(), std::bind2nd(std::less<double>(), 1)) != params->n_ridge_l0 ){
        sout << "ERROR : You must specify values for --l0 in (0,1).\n" << params->err_help;
        exit(-1);
      } 
    } else set_ridge_params(params->n_ridge_l0, params->lambda, params->err_help, sout);

    // user specified ridge parameters to use
    if( params->user_ridge_params_l1 ){
      if(params->n_ridge_l1 < 1){
        sout << "ERROR : Number of ridge parameters must be at least 1 in option --setl1.\n" << params->err_help;
        exit(-1);
      }
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

  // determine number of threads if not specified
  if(params->threads < 1){
    params->threads = std::thread::hardware_concurrency(); //may return 0 when not able to detect
    if(params->threads < 1) params->threads = 1;
  }

  // set Firth as default if both Firth and SPA are specified
  if(params->use_SPA && params->firth) params->use_SPA = false;

  // check SPA fallback pvalue threshold
  if(params->use_SPA && ((params->alpha_pvalue < params->nl_dbl_dmin) || (params->alpha_pvalue > 1 - params->numtol)) ){
    sout << "ERROR :SPA fallback p-value threshold must be in (0,1).\n" << params->err_help ;
    exit(-1);
  }
  // check firth fallback pvalue threshold
  if(params->firth && ((params->alpha_pvalue < params->nl_dbl_dmin) || (params->alpha_pvalue > 1 - params->numtol)) ){
    sout << "ERROR :Firth fallback p-value threshold must be in (0,1).\n" << params->err_help ;
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
    sout << "ERROR :Invalid chromosome specified to be tested.\n" << params->err_help ;
    exit(-1);
  }

  if(params->test_mode && !params->skip_blups && files->blup_file == "NULL") {
    sout << "ERROR :You must specify --pred if using --step 2 (otherwise use --ignore-pred).\n" << params->err_help ;
    exit(-1);
  }

  if( params->n_genofiles < 1 ) {
    sout << "ERROR :You must supply an input file using one of --bed, --pgen, or --bgen.\n" << params->err_help ;
    exit(-1);
  }
  if( params->n_genofiles > 1 ){
    sout << "ERROR :You must use either --bed,--bgen or --pgen.\n" << params->err_help ;
    exit(-1);
  }

  if(params->test_mode && (params->file_type == "pgen") && !params->streamBGEN){
    sout << "ERROR :Cannot use --nostream with PGEN format.\n" << params->err_help ;
    exit(-1);
  }
  if( params->bgenSample && (params->file_type != "bgen") ) {
    sout << "ERROR :You must supply a BGEN file corresponding to the sample file (use `--bgen`).\n" << params->err_help ;
    exit(-1);
  }
  if(params->block_size < 1) {
    sout << "ERROR : You must set --b.\n" << params->err_help ;
    exit(-1);
  }
  if((params->file_type == "bgen") & (!file_exists (files->bgen_file))) {
    sout << "ERROR : " << files->bgen_file  << " doesn't exist.\n" << params->err_help ;
    exit(-1);
  }
  if((files->cov_file != "NULL") & (!file_exists (files->cov_file))) {
    sout << "ERROR : " << files->cov_file  << " doesn't exist.\n" << params->err_help ;
    exit(-1);
  }
  if((files->pheno_file != "NULL") & (!file_exists (files->pheno_file))) {
    sout << "ERROR : " << files->pheno_file  << " doesn't exist.\n" << params->err_help ;
    exit(-1);
  }
  if((params->file_type == "bed") & (!file_exists (files->bedfile))) {
    sout << "ERROR : " << files->bedfile  << " doesn't exist.\n" << params->err_help ;
    exit(-1);
  }
  if((params->file_type == "bed") & (!file_exists (files->famfile))) {
    sout << "ERROR : " << files->famfile  << " doesn't exist.\n"  << params->err_help;
    exit(-1);
  }
  if((params->file_type == "bed") & (!file_exists (files->bimfile))) {
    sout << "ERROR : " << files->bimfile  << " doesn't exist.\n" << params->err_help ;
    exit(-1);
  }

}

void start_log(int argc, char **argv, const string out_file, MeasureTime* mt, mstream& sout){

  string trimmed_str;
  string log_name = out_file + ".log";
  sout.coss.open(log_name.c_str(), ios::out | ios::trunc); 

  mt->init();
  sout << "Start time: " << ctime( &(mt->start_time_info) ) << endl; 
  print_header(sout.coss);
  print_header(cout);
  sout << "Log of output saved in file : " << log_name << endl<< endl;

  // print options
  sout << "Command line arguments:";
  for(size_t counter=1;counter<argc;counter++){	  
    trimmed_str = string(argv[counter]);  // trim this
    trimmed_str.erase(std::remove_if(trimmed_str.begin(), trimmed_str.end(), ::isspace), trimmed_str.end());

    if( trimmed_str[0] == '-') sout << endl << "  ";
    sout << trimmed_str << " ";
  }
  sout << endl << endl;

}

void set_ridge_params(int nparams, vector<double>& in_param, const string err_help, mstream& sout){

  if(nparams < 2){
    sout << "ERROR : Number of ridge parameters must be at least 2\n" << err_help;
    exit(-1);
  } else {
    // endpoints are 0.01 and 0.99 
    double step = 1.0 / ( nparams - 1 );
    double val = step;
    in_param.resize( nparams);

    for( size_t index_p = 1; index_p < (nparams - 1); index_p++, val += step) in_param[index_p] = val;
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
    if(files->loco_tmp_prefix == "NULL") files->loco_tmp_prefix = files->out_file;
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


