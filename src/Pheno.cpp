/* 

   This file is part of the regenie software package.

   Copyright (c) 2020-2021 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

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
#include "Files.hpp"
#include "Pheno.hpp"

using namespace std;
using namespace Eigen;
using namespace boost;


void read_pheno_and_cov(struct in_files* files, struct param* params, struct filter* filters, struct phenodt* pheno_data, struct ests* m_ests, mstream& sout) {

  ArrayXb ind_in_pheno_and_geno = ArrayXb::Constant( params->n_samples, false );
  ArrayXb ind_in_cov_and_geno = ArrayXb::Constant( params->n_samples, files->cov_file.empty());

  if(params->getCorMat){ // intitiate values for phenotype info

    params->n_pheno = 1;
    ind_in_pheno_and_geno = ArrayXb::Constant( params->n_samples, true );
    pheno_data->phenotypes = MatrixXd::Constant( params->n_samples, 1, 1);
    pheno_data->masked_indivs = MatrixXb::Constant( params->n_samples, 1, true);
    pheno_data->Neff = pheno_data->masked_indivs.cast<double>().colwise().sum();

  } else { // read in phenotype (mean-impute for QT)

    pheno_read(params, files, filters, pheno_data, ind_in_pheno_and_geno, sout);
    if(params->binary_mode && !params->test_mode)
      m_ests->offset_logreg = MatrixXd::Zero(params->n_samples, params->n_pheno); 

  }

  // Intercept
  pheno_data->new_cov = MatrixXd::Ones(params->n_samples, 1);
  if(params->strict_mode) pheno_data->new_cov.array() *= pheno_data->masked_indivs.col(0).array().cast<double>();;

  // read in covariates
  if(!files->cov_file.empty()) covariate_read(params, files, filters, pheno_data, ind_in_cov_and_geno, sout);

  // mask individuals 
  filters->ind_in_analysis = ind_in_pheno_and_geno && ind_in_cov_and_geno;
  pheno_data->masked_indivs.array().colwise() *= filters->ind_in_analysis;
  if( params->strict_mode ) 
    filters->ind_in_analysis = filters->ind_in_analysis && pheno_data->masked_indivs.col(0).array();
  pheno_data->phenotypes.array().colwise() *= filters->ind_in_analysis.cast<double>();
  if(params->binary_mode) pheno_data->phenotypes_raw.array().colwise() *= filters->ind_in_analysis.cast<double>();
  pheno_data->new_cov.array().colwise() *= filters->ind_in_analysis.cast<double>();

  // identify individuals with at least one phenotpe with NA
  filters->has_missing = !(pheno_data->masked_indivs.array().rowwise().all());
  //for(int i = 0; i <5; i++) cerr << std::boolalpha << filters->has_missing(i) << endl;

  // check sample size
  if( filters->ind_in_analysis.cast<int>().sum() < 1 ) {
    sout << "ERROR: Sample size cannot be < 1\n";
    exit(EXIT_FAILURE);
  }
  sout << " * number of individuals used in analysis = " << filters->ind_in_analysis.cast<int>().sum() << endl;

}

void pheno_read(struct param* params, struct in_files* files, struct filter* filters, struct phenodt* pheno_data, ArrayXb& ind_in_pheno_and_geno, mstream& sout) {

  uint32_t indiv_index;
  bool all_miss;
  double mean;
  string line;
  std::vector< string > tmp_str_vec;
  ArrayXb keep_cols;
  findID person;
  Files fClass;

  sout << left << std::setw(20) << " * phenotypes" << ": [" << files->pheno_file << "] ";
  fClass.openForRead(files->pheno_file, sout);
  fClass.readLine(line);
  check_str(line); // remove carriage returns at the end of line if any

  // check that FID and IID are first two entries in header
  tmp_str_vec = string_split(line,"\t ");
  if( (tmp_str_vec[0] != "FID") || (tmp_str_vec[1] != "IID") ) {
    sout << "ERROR: Header of phenotype file must start with: FID IID" << endl;
    exit(EXIT_FAILURE);
  }

  // get phenotype names 
  keep_cols = ArrayXb::Constant(tmp_str_vec.size() - 2, true);
  for(int i = 0; i < keep_cols.size(); i++ ) {
    if(params->select_phenos) // check if keeping pheno
      keep_cols(i) = in_map(tmp_str_vec[i+2], filters->pheno_colKeep_names);

    if(keep_cols(i)) files->pheno_names.push_back( tmp_str_vec[i+2] );
  }
  params->n_pheno = keep_cols.count();

  // check #pheno
  if(params->n_pheno < 1){
    sout << "ERROR: Need at least one phenotype." << endl;
    exit(EXIT_FAILURE);
  } 
  sout << "n_pheno = " << params->n_pheno << endl;
  params->strict_mode = params->n_pheno == 1; // drop all missing observations

  // how missingness is handles
  if( params->strict_mode ) sout << "   -dropping observations with missing values at any of the phenotypes" << endl;
  else if( !params->rm_missing_qt  && !params->binary_mode) sout << "   -keeping and mean-imputing missing observations (done for each trait)" << endl;


  // allocate memory
  pheno_data->phenotypes = MatrixXd::Zero(params->n_samples, params->n_pheno);
  pheno_data->masked_indivs = MatrixXb::Constant(params->n_samples, params->n_pheno, true);
  if(params->binary_mode)  
    pheno_data->phenotypes_raw = MatrixXd::Zero(params->n_samples, params->n_pheno);

  VectorXd total, ns;
  total.setZero(params->n_pheno);
  ns.setZero(params->n_pheno);

  // read in data
  while( fClass.readLine(line) ){
    tmp_str_vec = string_split(line,"\t ");

    if( (int)tmp_str_vec.size() != (2+keep_cols.size()) ){
      sout << "ERROR: Incorrectly formatted phenotype file." << endl;
      exit(EXIT_FAILURE);
    }

    person = getIndivIndex(tmp_str_vec[0], tmp_str_vec[1], params, sout);
    if(!person.is_found) continue;

    indiv_index = person.index;

    // check duplicate
    if( !ind_in_pheno_and_geno(indiv_index) ){
      ind_in_pheno_and_geno( indiv_index ) = true;
    } else {
      sout << "ERROR: Individual appears more than once in phenotype file: FID=" << tmp_str_vec[0] << " IID=" << tmp_str_vec[1] << endl;
      exit(EXIT_FAILURE);
    }

    // read phenotypes 
    all_miss = true;
    for(int j = 0, i_pheno = 0; j < keep_cols.size(); j++) {

      if( !keep_cols(j) ) continue;

      pheno_data->phenotypes(indiv_index, i_pheno) = convertDouble(tmp_str_vec[2+j], params, sout);

      // for BT, save raw data and force 0/1 values
      if (params->binary_mode) {

        if(!params->CC_ZeroOne && (pheno_data->phenotypes(indiv_index, i_pheno) != params->missing_value_double)) 
          pheno_data->phenotypes(indiv_index, i_pheno) -= 1; // if using 1/2/NA encoding

        pheno_data->phenotypes_raw(indiv_index, i_pheno) = pheno_data->phenotypes(indiv_index, i_pheno);

        if(fabs(pheno_data->phenotypes_raw(indiv_index, i_pheno)) > params->numtol && fabs(pheno_data->phenotypes_raw(indiv_index, i_pheno)-1) > params->numtol ) {

          if(params->within_sample_l0){
            sout << "ERROR: No missing value allowed in phenotype file with option -within" << endl;
            exit(EXIT_FAILURE);
          } else if( pheno_data->phenotypes_raw(indiv_index, i_pheno) != params->missing_value_double ) {
            sout << "ERROR: A phenotype value is not "<<
              (params->CC_ZeroOne ? "0/1/NA" : "1/2/NA") <<
              " for individual: FID=" << tmp_str_vec[0] << 
              " IID=" << tmp_str_vec[1] << 
              " Y=" << tmp_str_vec[2+j] << endl;
            //sout << "Use flag '--1' for 1/2/NA encoding [1=control|2=case|NA=missing]." << endl;
            exit(EXIT_FAILURE);
          }

          pheno_data->phenotypes_raw(indiv_index, i_pheno) = params->missing_value_double;
          pheno_data->masked_indivs(indiv_index, i_pheno) = false;
          if( params->strict_mode ) pheno_data->masked_indivs.row(indiv_index) = MatrixXb::Constant(1, params->n_pheno, false);

        }
      }

      if( pheno_data->phenotypes(indiv_index, i_pheno) != params->missing_value_double ) {
        total(i_pheno) +=  pheno_data->phenotypes(indiv_index, i_pheno);
        ns(i_pheno) +=  1;
        all_miss = false;
      } else {
        if( params->test_mode && params->rm_missing_qt ) pheno_data->masked_indivs(indiv_index, i_pheno) = false;
        if( params->strict_mode ) pheno_data->masked_indivs.row(indiv_index) = MatrixXb::Constant(1, params->n_pheno, false);
      }

      i_pheno++;
    }

    if( all_miss ) ind_in_pheno_and_geno( indiv_index ) = false; // if individual has no phenotype data at all

  }

  // mask individuals in genotype data but not in phenotype data
  pheno_data->masked_indivs.array().colwise() *= ind_in_pheno_and_geno;

  // check if all individuals have missing/invalid phenotype
  if(pheno_data->masked_indivs.cast<int>().colwise().sum().array().minCoeff() == 0){
    sout << "ERROR: All individuals have missing/invalid phenotype values." << endl;
    exit(EXIT_FAILURE);
  }

  if(!params->binary_mode || !params->test_mode){

    if(!params->binary_mode){
      // impute missing with mean
      for(size_t i = 0; i < params->n_samples;i++) 
        for(int j = 0; j < params->n_pheno;j++) {
          if( pheno_data->phenotypes(i,j) != params->missing_value_double ) {
            pheno_data->phenotypes(i,j) -= total(j) / ns(j);	  
          }  else pheno_data->phenotypes(i,j) = 0.0;
        }
    } else {
      for(int j = 0; j < params->n_pheno; j++) {
        mean = (pheno_data->masked_indivs.col(j).array()).select( pheno_data->phenotypes.col(j).array(), 0).sum() / pheno_data->masked_indivs.col(j).cast<double>().sum();
        pheno_data->phenotypes.col(j).array() = (pheno_data->masked_indivs.col(j).array()).select(pheno_data->phenotypes.col(j).array() - mean, 0);
      }
    }

    // apply masking
    pheno_data->phenotypes.array() *= pheno_data->masked_indivs.array().cast<double>();

  }

  // number of phenotyped individuals 
  sout <<  "   -number of phenotyped individuals = " << ind_in_pheno_and_geno.cast<int>().sum() << endl;

  // check that there cases are present
  if(params->binary_mode){
    for(int j = 0; j < params->n_pheno;j++) {
      if( ( pheno_data->phenotypes_raw.col(j).array() == 1 ).count() == 0){
        sout << "ERROR: No cases present for phenotype: " << files->pheno_names[j] << endl; 
        exit(EXIT_FAILURE);
      }
    }
  }

  pheno_data->Neff = pheno_data->masked_indivs.cast<double>().colwise().sum();
  if(params->strict_mode) sout << "   -number of individuals remaining with non-missing phenotypes = " << pheno_data->Neff(0) << endl;

  fClass.closeFile();

}

void covariate_read(struct param* params, struct in_files* files, struct filter* filters, struct phenodt* pheno_data, ArrayXb& ind_in_cov_and_geno, mstream& sout) {

  int nc_cat = 0;
  uint32_t indiv_index;
  ArrayXb keep_cols;
  string line;
  std::vector< string > tmp_str_vec, covar_names;
  std::vector< std::map<std::string,int> > categories;
  findID person;
  Files fClass;

  sout << left << std::setw(20) << " * covariates" << ": [" << files->cov_file << "] " << flush;
  fClass.openForRead(files->cov_file, sout);
  fClass.readLine(line);
  check_str(line); // remove carriage returns at the end of line if any

  // check header
  tmp_str_vec = string_split(line,"\t ");
  if( (tmp_str_vec[0] != "FID") || (tmp_str_vec[1] != "IID") ) {
    sout << "ERROR: Header of covariate file must start with: FID IID" << endl;
    exit(EXIT_FAILURE);
  }

  // get covariate names 
  keep_cols = ArrayXb::Constant(tmp_str_vec.size() - 2, true);
  for(int i = 0; i < keep_cols.size(); i++) {

    if(!params->select_covs && !in_map(tmp_str_vec[i+2], filters->cov_colKeep_names)) // in case specified as categorical
      filters->cov_colKeep_names[tmp_str_vec[i+2]] = true;
    else keep_cols(i) = in_map(tmp_str_vec[i+2], filters->cov_colKeep_names);

    if(keep_cols(i)){
      covar_names.push_back( tmp_str_vec[i+2] );
      nc_cat += !filters->cov_colKeep_names[ tmp_str_vec[i+2] ];
    }
  }
  categories.resize(nc_cat);

  // check all covariates specified are in the file
  params->n_cov = keep_cols.count(); 
  if( (int)filters->cov_colKeep_names.size() != params->n_cov ) {
    sout << "ERROR: Not all covariates specified are found in the covariate file.\n";
    exit(EXIT_FAILURE);
  }

  // check #covariates is > 0
  if(params->n_cov < 1){ // only intercept will be included
    sout << "n_cov = " << params->n_cov << " (+ intercept)" << endl;
    ind_in_cov_and_geno = ArrayXb::Constant( params->n_samples, true );
    return ;
  }

  // allocate memory 
  pheno_data->new_cov = MatrixXd::Zero(params->n_samples, 1 + params->n_cov);
  pheno_data->new_cov.col(0) = MatrixXd::Ones(params->n_samples, 1);

  // read in data
  while( fClass.readLine(line) ){
    tmp_str_vec = string_split(line,"\t ");

    if( (int)tmp_str_vec.size() != (keep_cols.size()+2) ){
      sout << "ERROR: Incorrectly formatted covariate file." << endl;
      exit(EXIT_FAILURE);
    }

    person = getIndivIndex(tmp_str_vec[0], tmp_str_vec[1], params, sout);
    if(!person.is_found) continue;

    indiv_index = person.index;

    // check duplicate
    if( !ind_in_cov_and_geno(indiv_index) ){
      ind_in_cov_and_geno(indiv_index) = true;
    } else {
      sout << "ERROR: Individual appears more than once in covariate file: FID=" << tmp_str_vec[0] << " IID=" << tmp_str_vec[1] << endl;
      exit(EXIT_FAILURE);
    }

    // read covariate data and check for missing values
    for(int i_cov = 0, i_cat = 0, j = 0; j < keep_cols.size(); j++) {

      if( !keep_cols(j) ) continue;

      // if quantitative
      if( filters->cov_colKeep_names[ covar_names[i_cov] ] )
        pheno_data->new_cov(indiv_index, 1 + i_cov) = convertDouble(tmp_str_vec[2+j], params, sout);
      else // if categorical, convert to numeric category
        pheno_data->new_cov(indiv_index, 1 + i_cov) = convertNumLevel(tmp_str_vec[2+j], categories[i_cat++], params, sout);

      if( pheno_data->new_cov(indiv_index, 1 + i_cov) == params->missing_value_double ) { // ignore individual
        ind_in_cov_and_geno(indiv_index) = false;
        break;
      }

      i_cov++;
    }

  }
  //cerr << endl<<pheno_data->new_cov.block(0,0,5,pheno_data->new_cov.cols())<< endl;

  // add dummy variables if needed
  if(nc_cat > 0){
    int n_dummies;
    int n_add = check_categories(covar_names, categories, params, filters, sout) - nc_cat; // new columns to add (or remove if single category)
    MatrixXd full_covarMat (pheno_data->new_cov.rows(), pheno_data->new_cov.cols() + n_add);

    // copy intercept column
    full_covarMat.col(0) = pheno_data->new_cov.col(0);

    for(int i = 1, j = 1; i < pheno_data->new_cov.cols(); i++){
      //cerr << i << "/" << pheno_data->new_cov.cols()-1 << " - " << covar_names[i-1] << " is q: " << std::boolalpha << filters->cov_colKeep_names[ covar_names[i-1] ] << endl;

      if( filters->cov_colKeep_names[ covar_names[i-1] ] ) // copy col
        full_covarMat.col(j++) = pheno_data->new_cov.col(i);
      else { 
        n_dummies = pheno_data->new_cov.col(i).maxCoeff();
        if( n_dummies > 0 ){
          full_covarMat.block(0, j, full_covarMat.rows(), n_dummies) = get_dummies(pheno_data->new_cov.col(i));
          j+= n_dummies;
        }
      }
    }

    pheno_data->new_cov = full_covarMat;
    params->n_cov = full_covarMat.cols() - 1; // ignore intercept
  }
  //cerr << endl<<pheno_data->new_cov.block(0,0,5,pheno_data->new_cov.cols())<< endl;

  // mask individuals in genotype data but not in covariate data
  pheno_data->masked_indivs.array().colwise() *= ind_in_cov_and_geno;

  // apply masking
  pheno_data->new_cov.array().colwise() *= ind_in_cov_and_geno.cast<double>();
  if(params->strict_mode) pheno_data->new_cov.array().colwise() *= pheno_data->masked_indivs.col(0).array().cast<double>();

  sout << "n_cov = " << params->n_cov << endl;
  sout <<  "   -number of individuals with covariate data = " << ind_in_cov_and_geno.cast<int>().sum() << endl;

  fClass.closeFile();

}

int check_categories(vector<std::string>& covar, vector<std::map<std::string,int>>& categories, struct param* params, struct filter* filters, mstream& sout){

  int ntotal = 0;

  for(size_t i = 0, j = 0; i < covar.size(); i++){

    // skip qCovar
    if( filters->cov_colKeep_names[ covar[i] ] ) continue;

    if((int)categories[j].size() > params->max_cat_levels){
      sout << "ERROR: Too many categories for covariate: " << covar[i] << ". Either use '--maxCatLevels' or combine categories.\n";
      exit(EXIT_FAILURE);
    }

    ntotal += categories[j++].size() - 1; // add K-1 dummy vars
  }

  return ntotal;
}

MatrixXd get_dummies(const Eigen::Ref<const Eigen::MatrixXd>& numCov) {

  int index, nvars;
  nvars = numCov.maxCoeff();

  MatrixXd dout = MatrixXd::Zero(numCov.rows(), nvars);

  for(int i = 0; i < numCov.rows(); i++){
    if(numCov(i,0) == 0) continue; // will go to intercept

    index = numCov(i,0) - 1;
    dout(i, index) = 1;
  }

  return dout;

}

// Adjust for covariates (incl. intercept)
// in step 2, also read blups and check
void prep_run (struct in_files* files, struct param* params, struct phenodt* pheno_data, struct ests* m_ests, mstream& sout){

  // for step 2, check blup files
  if (params->test_mode && !params->getCorMat){
    // individuals not in blup file will have their phenotypes masked
    blup_read(files, params, pheno_data, m_ests, sout);
    if(params->write_samples) write_ids(files, params, pheno_data, sout);
  }

  // compute N for each trait
  pheno_data->Neff = pheno_data->masked_indivs.cast<double>().colwise().sum();

  // orthonormal basis
  getCovBasis(pheno_data->new_cov, params);

  // compute offset for BT (only in step 1)
  if(params->binary_mode && !params->test_mode) fit_null_logistic(0, params, pheno_data, m_ests, sout);

  // residualize phenotypes (skipped for BTs when testing)
  if( !params->getCorMat && (!params->test_mode || !params->binary_mode) ) residualize_phenotypes(params, pheno_data, files->pheno_names, sout);

}

// get list of blup files
void blup_read(struct in_files* files, struct param* params, struct phenodt* pheno_data, struct ests* m_ests, mstream& sout) {

  int n_files = 0, tmp_index, n_masked_prior, n_masked_post;
  uint32_t indiv_index;
  double blup_val;
  string line, tmp_pheno;
  std::vector< string > tmp_str_vec, tmp_prs_vec;
  vector<bool> read_pheno(params->n_pheno, false);
  MatrixXb blupf_mask;
  Files fClass;

  // allocate memory
  m_ests->blups = MatrixXd::Zero(params->n_samples, params->n_pheno);

  // skip reading if specified by user
  if( params->skip_blups ) {
    sout << " * no step 1 predictions given. Simple " << ( params->binary_mode ? "logistic":"linear" ) << " regression will be performed" <<endl;
    return;
  }

  sout << " * " << (params->use_prs ? "PRS" : "LOCO") << " predictions : [" << files->blup_file << "] ";
  fClass.openForRead(files->blup_file, sout);

  // get list of files containing blups
  while (fClass.readLine(line)){
    tmp_str_vec = string_split(line,"\t ");

    // each line contains a phenotype name and the corresponding blup file name
    if( tmp_str_vec.size() != 2 ){
      sout << "ERROR: Incorrectly formatted blup list file : " << files->blup_file << endl;
      exit(EXIT_FAILURE);
    }

    // get index of phenotype in phenotype matrix
    vector<string>::iterator it = std::find(files->pheno_names.begin(), files->pheno_names.end(), tmp_str_vec[0]);
    if (it == files->pheno_names.end()) continue; // ignore unrecognized phenotypes

    tmp_index = std::distance(files->pheno_names.begin(), it);
    files->pheno_index.push_back(tmp_index);

    // check that phenotype only has one file
    if(read_pheno[tmp_index]){
      sout << "ERROR: Phenotype \'" << tmp_pheno << "\' appears more than once in blup list file." << endl;
      exit(EXIT_FAILURE);
    }

    n_files++;
    read_pheno[tmp_index] = true;
    files->blup_files.push_back(tmp_str_vec[1]);
  }

  // force all phenotypes in phenotype file to be used
  if(n_files != params->n_pheno) {
    sout << "ERROR : Number of files (" << n_files <<")  is not equal to the number of phenotypes.\n" ;
    exit(EXIT_FAILURE);
  }
  sout << "n_files = " << n_files << endl;
  fClass.closeFile();

  // read blup file for each phenotype
  for(int ph = 0; ph < params->n_pheno; ph++) {
    int i_pheno = files->pheno_index[ph];

    sout << "   -file [" <<  files->blup_files[ph];
    sout << "] for phenotype \'" << files->pheno_names[i_pheno] << "\'\n";
    fClass.openForRead(files->blup_files[ph], sout);

    // mask all individuals not present in .loco file
    blupf_mask = MatrixXb::Constant(params->n_samples, 1, false);
    n_masked_prior = pheno_data->masked_indivs.col(ph).cast<int>().sum();
    // only read first line which has FID_IID
    fClass.readLine(line);
    tmp_str_vec = string_split(line,"\t ");

    if( tmp_str_vec[0] != "FID_IID") {
      sout << "ERROR: Header of blup file must start with FID_IID (=" << tmp_str_vec[0] << ").\n";
      exit(EXIT_FAILURE);
    }

    // read second line to check for missing predictions
    fClass.readLine(line);
    tmp_prs_vec = string_split(line,"\t ");

    if( params->use_prs  && (tmp_prs_vec[0] != "0") ){ // read in second line
      sout << "ERROR: Second line must start with 0 (=" << tmp_prs_vec[0] << ").\n";
      exit(EXIT_FAILURE);
    }

    for (size_t i = 1; i < tmp_str_vec.size(); i++){
      // ignore sample if it is not in genotype data
      if (!in_map(tmp_str_vec[i], params->FID_IID_to_ind)) continue;
      indiv_index = params->FID_IID_to_ind[tmp_str_vec[i]];
      blup_val = convertDouble(tmp_prs_vec[i], params, sout);

      // ignore samples where prediction is NA
      blupf_mask( indiv_index , 0 ) = (blup_val != params->missing_value_double);
      //cerr << tmp_str_vec[i] << "\t" << std::boolalpha << blupf_mask( indiv_index , 0 ) << endl; 
      if (!blupf_mask( indiv_index , 0 )) continue;

      if( params->use_prs ) 
        m_ests->blups(indiv_index, ph) = blup_val;
    }

    // mask samples not in file
    pheno_data->masked_indivs.col(ph).array() = pheno_data->masked_indivs.col(ph).array() && blupf_mask.col(0).array();
    n_masked_post = pheno_data->masked_indivs.col(ph).cast<int>().sum();

    if( n_masked_post < n_masked_prior ){
      sout << "    + " << n_masked_prior - n_masked_post <<
        " individuals with missing LOCO predictions will be ignored for the trait\n";
    }

    // check not everyone is masked
    if( n_masked_post < 1 ){
      sout << "ERROR: none of the individuals remaining in the analysis are in the LOCO predictions file from step 1.\n Either re-run step 1 including individuals in current genotype file or use option '--ignore-pred'.\n";
    exit(EXIT_FAILURE);
    }

    fClass.closeFile();
  }

}


// write ids of samples included in step 2 (done for each trait)
void write_ids(struct in_files* files, struct param* params, struct phenodt* pheno_data, mstream& sout){

  uint32_t index;
  map<string, uint32_t >::iterator itr_ind;
  string idfile;
  Files fout;

    sout << " * user specified to write sample IDs for each trait"<<endl;

  for( int ph = 0; ph < params->n_pheno; ph++){
    idfile = files->out_file + "_" + files->pheno_names[ph] + ".regenie.ids";
    fout.openForWrite(idfile, sout);

    // print phenotype name on 1st line (ensure 2 column format)
    if( params->print_pheno_name ) fout << files->pheno_names[ph] << "\tNA\n"; 

    // go through map and check if individual is not masked
    for (itr_ind = params->FID_IID_to_ind.begin(); itr_ind != params->FID_IID_to_ind.end(); ++itr_ind) {

      index = itr_ind->second;

      if( !pheno_data->masked_indivs(index, ph) ) continue;
      fout << params->FIDvec[index][0] << "\t" << params->FIDvec[index][1] << endl;

    }

    fout.closeFile();
  } 

  if(!params->write_masks) params->FIDvec.clear();
}


void getCovBasis(MatrixXd& new_cov,struct param* params){

  // eigen-decompose cov matrix
  MatrixXd xtx = new_cov.transpose() * new_cov;
  SelfAdjointEigenSolver<MatrixXd> es(xtx);
  VectorXd D = es.eigenvalues();
  MatrixXd V = es.eigenvectors();

  // create basis set
  // eigenvalues sorted in increasing order
  int non_zero_eigen = (D.array() > D.tail(1)(0) * params->eigen_val_rel_tol).count();
  RowVectorXd vv1 = D.tail(non_zero_eigen).array().sqrt();
  new_cov *= V.rightCols(non_zero_eigen);
  new_cov.array().rowwise() /= vv1.array();

  params->ncov = non_zero_eigen; // save number of lin. indep. covars.

}


void residualize_phenotypes(struct param* params, struct phenodt* pheno_data, const std::vector<std::string>& pheno_names, mstream& sout) {
  sout << "   -residualizing and scaling phenotypes...";
  auto t1 = std::chrono::high_resolution_clock::now();

  // residuals (centered) then scale
  MatrixXd beta = pheno_data->phenotypes.transpose() * pheno_data->new_cov;
  pheno_data->phenotypes -= ( (pheno_data->new_cov * beta.transpose()).array() * pheno_data->masked_indivs.array().cast<double>() ).matrix();
  pheno_data->scale_Y = pheno_data->phenotypes.colwise().norm().array() / sqrt(pheno_data->Neff.matrix().transpose().array() - params->ncov);

  // check sd is not 0 
  MatrixXd::Index minIndex;
  if(pheno_data->scale_Y.minCoeff(&minIndex) < params->numtol){
    sout << "ERROR: Phenotype \'" << pheno_names[minIndex] << "\' has sd=0." << endl;
    exit(EXIT_FAILURE);
  }
  pheno_data->phenotypes.array().rowwise() /= pheno_data->scale_Y.array();

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;
}

void check_str(string& mystring ){

  // check there are no '\r' at the end
  if( !mystring.empty() && mystring[mystring.size() - 1] == '\r' )
    mystring.erase(mystring.size() - 1);

}

