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
#include "Pheno.hpp"

using namespace std;
using namespace Eigen;
using namespace boost;


void read_pheno_and_cov(struct in_files* files, struct param* params, struct filter* filters, struct phenodt* pheno_data, struct ests* m_ests, mstream& sout) {


  ArrayXb ind_in_pheno_and_geno = ArrayXb::Constant( params->n_samples, false );
  ArrayXb ind_in_cov_and_geno = ArrayXb::Constant( params->n_samples, files->cov_file == "NULL" );

  // read in phenotype (mean-impute for QT)
  pheno_read(params, files, filters, pheno_data, ind_in_pheno_and_geno, sout);
  if(params->binary_mode && !params->test_mode)
    m_ests->offset_logreg = MatrixXd::Zero(params->n_samples, params->n_pheno); 

  // Intercept
  pheno_data->new_cov = MatrixXd::Ones(params->n_samples, 1);
  if(params->strict_mode) pheno_data->new_cov.array() *= pheno_data->masked_indivs.col(0).array().cast<double>();;

  // read in covariates
  if(files->cov_file != "NULL") covariate_read(params, files, filters, pheno_data, ind_in_cov_and_geno, sout);

  // mask individuals 
  filters->ind_in_analysis = ind_in_pheno_and_geno * ind_in_cov_and_geno;
  pheno_data->masked_indivs.array().colwise() *= filters->ind_in_analysis;
  if( params->strict_mode ) filters->ind_in_analysis *= pheno_data->masked_indivs.col(0).array();
  pheno_data->phenotypes.array().colwise() *= filters->ind_in_analysis.cast<double>();
  if(params->binary_mode) pheno_data->phenotypes_raw.array().colwise() *= filters->ind_in_analysis.cast<double>();
  pheno_data->new_cov.array().colwise() *= filters->ind_in_analysis.cast<double>();
  pheno_data->Neff = pheno_data->masked_indivs.cast<double>().colwise().sum();

  // check sample size
  if( filters->ind_in_analysis.cast<int>().sum() < 1 ) {
    sout << "ERROR: Sample size cannot be < 1\n";
    exit(-1);
  }
  sout << " * number of individuals used in analysis = " << filters->ind_in_analysis.cast<int>().sum() << endl;


  // orthonormal basis
  getCovBasis(pheno_data->new_cov, params);

  // compute offset for BT (only in step 1)
  if(params->binary_mode && !params->test_mode) fit_null_logistic(0, params, pheno_data, m_ests, sout);

  // residualize phenotypes (skipped for BTs when testing)
  if(!params->test_mode || !params->binary_mode) residualize_phenotypes(params, pheno_data, sout);

}

void pheno_read(struct param* params, struct in_files* files, struct filter* filters, struct phenodt* pheno_data, ArrayXb& ind_in_pheno_and_geno, mstream& sout) {

  uint32_t indiv_index;
  bool keep_pheno;
  double mean;
  string line;
  std::vector< string > tmp_str_vec;
  vector<bool> pheno_colKeep;
  findID person;
  ifstream myfile;

  sout << left << std::setw(20) << " * phenotypes" << ": [" << files->pheno_file << "] ";
  myfile.open (files->pheno_file.c_str(), ios::in);
  if (!myfile.is_open()) {    
    sout << "ERROR: Cannot open phenotype file : " << files->pheno_file << endl;
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
  params->n_pheno = 0;
  for( size_t filecol = 2; filecol < tmp_str_vec.size(); filecol++ ) {
    // default is to use all columns  
    if(!params->select_phenos){
      pheno_colKeep.push_back( true );
      params->n_pheno++;
      files->pheno_names.push_back( tmp_str_vec[filecol] );
    } else {
      // otherwise, check if phenotype is in list of phenotypes to keep
      keep_pheno = std::find(filters->pheno_colKeep_names.begin(), filters->pheno_colKeep_names.end(), tmp_str_vec[filecol]) != filters->pheno_colKeep_names.end();
      pheno_colKeep.push_back( keep_pheno );
      if(keep_pheno){
        params->n_pheno++;
        files->pheno_names.push_back( tmp_str_vec[filecol] );
      }
    }
  }

  // check #pheno is > 0
  if(params->n_pheno < 1){
    sout << "ERROR: Need at least one phenotype." << endl;
    exit(-1);
  }
  sout << "n_pheno = " << params->n_pheno << endl;

  // if n_pheno = 1, drop all missing observations
  if(params->n_pheno == 1) params->strict_mode = true; 

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
  while( getline (myfile,line) ){
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < pheno_colKeep.size() ){
      sout << "ERROR: Incorrectly formatted phenotype file." << endl;
      exit(-1);
    }

    person = getIndivIndex(tmp_str_vec[0], tmp_str_vec[1], params, sout);
    if(!person.is_found) continue;

    indiv_index = person.index;

    // check duplicate
    if( !ind_in_pheno_and_geno(indiv_index) ){
      ind_in_pheno_and_geno( indiv_index ) = true;
    } else {
      sout << "ERROR: Individual appears more than once in phenotype file: FID=" << tmp_str_vec[0] << " IID=" << tmp_str_vec[1] << endl;
      exit(-1);
    }

    // read phenotypes 
    for(int i_pheno = 0, j = 0; j < pheno_colKeep.size(); j++) {

      if( !pheno_colKeep[j] ) continue;

      pheno_data->phenotypes(indiv_index, i_pheno) = convertDouble(tmp_str_vec[2+j], params, sout);

      // for BT, save raw data and force 0/1 values
      if (params->binary_mode) {

        if(!params->CC_ZeroOne && (pheno_data->phenotypes(indiv_index, i_pheno) != params->missing_value_double)) 
          pheno_data->phenotypes(indiv_index, i_pheno) -= 1; // if using 1/2/NA encoding

        pheno_data->phenotypes_raw(indiv_index, i_pheno) = pheno_data->phenotypes(indiv_index, i_pheno);

        if(fabs(pheno_data->phenotypes_raw(indiv_index, i_pheno)) > params->numtol && fabs(pheno_data->phenotypes_raw(indiv_index, i_pheno)-1) > params->numtol ) {

          if(params->within_sample_l0){
            sout << "ERROR: No missing value allowed in phenotype file with option -within" << endl;
            exit(-1);
          } else if( pheno_data->phenotypes_raw(indiv_index, i_pheno) != params->missing_value_double ) {
            sout << "ERROR: A phenotype value is not 0/1/NA for individual: FID=" << tmp_str_vec[0] << " IID=" << tmp_str_vec[1] << endl;
            sout << "Use flag '--1' for 1/2/NA encoding [1=control|2=case|NA=missing]." << endl;
            exit(-1);
          }

          pheno_data->phenotypes_raw(indiv_index, i_pheno) = params->missing_value_double;
          pheno_data->masked_indivs(indiv_index, i_pheno) = false;
          if( params->strict_mode ) pheno_data->masked_indivs.row(indiv_index) = MatrixXb::Constant(1, params->n_pheno, false);

        }
      }

      if( pheno_data->phenotypes(indiv_index, i_pheno) != params->missing_value_double ) {
        total(i_pheno) +=  pheno_data->phenotypes(indiv_index, i_pheno);
        ns(i_pheno) +=  1;
      } else {
        if( params->test_mode && params->rm_missing_qt ) pheno_data->masked_indivs(indiv_index, i_pheno) = false;
        if( params->strict_mode ) pheno_data->masked_indivs.row(indiv_index) = MatrixXb::Constant(1, params->n_pheno, false);
      }

      i_pheno++;
    }

  }

  // mask individuals in genotype data but not in phenotype data
  pheno_data->masked_indivs.array().colwise() *= ind_in_pheno_and_geno;

  // check if all individuals have missing/invalid phenotype
  if(pheno_data->masked_indivs.cast<int>().colwise().sum().array().minCoeff() == 0){
    sout << "ERROR: All individuals have missing/invalid phenotype values." << endl;
    exit(-1);
  }

  if(!params->binary_mode || !params->test_mode){

    if(!params->binary_mode){
      // impute missing with mean
      for(size_t i = 0; i < params->n_samples;i++) 
        for(size_t j = 0; j < params->n_pheno;j++) {
          if( pheno_data->phenotypes(i,j) != params->missing_value_double ) {
            pheno_data->phenotypes(i,j) -= total(j) / ns(j);	  
          }  else pheno_data->phenotypes(i,j) = 0.0;
        }
    } else {
      for(size_t j = 0; j < params->n_pheno; j++) {
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
    for(size_t j = 0; j < params->n_pheno;j++) {
      if( ( pheno_data->phenotypes_raw.col(j).array() == 1 ).count() == 0){
        sout << "ERROR: No cases present for phenotype: " << files->pheno_names[j] << endl; 
        exit(-1);
      }
    }
  }

  pheno_data->Neff = pheno_data->masked_indivs.cast<double>().colwise().sum();
  if(params->strict_mode) sout << "   -number of individuals remaining with non-missing phenotypes = " << pheno_data->Neff(0) << endl;

  myfile.close();

}

void covariate_read(struct param* params, struct in_files* files, struct filter* filters, struct phenodt* pheno_data, ArrayXb& ind_in_cov_and_geno, mstream& sout) {

  uint32_t indiv_index;
  bool keep_cov;
  vector<bool> cov_colKeep;
  string line;
  std::vector< string > tmp_str_vec ;
  findID person;
  ifstream myfile;

  sout << left << std::setw(20) << " * covariates" << ": [" << files->cov_file << "] " << flush;
  myfile.open (files->cov_file.c_str(), ios::in);
  if (!myfile.is_open()) {
    sout << "ERROR: Cannot open covariate file : " << files->cov_file << endl;
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
  params->n_cov = 0;
  for( size_t filecol = 2; filecol < tmp_str_vec.size(); filecol++ ) {
    // default is to use all columns  
    if(!params->select_covs){
      cov_colKeep.push_back( true );
      params->n_cov++;
    } else {
      // otherwise, check if covariate is in list of covariates to keep
      keep_cov = std::find(filters->cov_colKeep_names.begin(), filters->cov_colKeep_names.end(), tmp_str_vec[filecol]) != filters->cov_colKeep_names.end();
      cov_colKeep.push_back( keep_cov );
      if(keep_cov) params->n_cov++;
    }
  }

  // check #covariates is > 0
  if(params->n_cov < 1){ // only intercept will be included
    sout << "n_cov = " << params->n_cov << " (+ intercept)" << endl;
    ind_in_cov_and_geno = ArrayXb::Constant( params->n_samples, true );
    return ;
  }
  sout << "n_cov = " << params->n_cov << flush;

  // allocate memory 
  pheno_data->new_cov = MatrixXd::Zero(params->n_samples, 1 + params->n_cov);
  pheno_data->new_cov.col(0) = MatrixXd::Ones(params->n_samples, 1);

  // read in data
  while( getline (myfile,line) ){
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < cov_colKeep.size() ){
      sout << "ERROR: Incorrectly formatted covariate file." << endl;
      exit(-1);
    }

    person = getIndivIndex(tmp_str_vec[0], tmp_str_vec[1], params, sout);
    if(!person.is_found) continue;

    indiv_index = person.index;

    // check duplicate
    if( !ind_in_cov_and_geno(indiv_index) ){
      ind_in_cov_and_geno(indiv_index) = true;
    } else {
      sout << "ERROR: Individual appears more than once in covariate file: FID=" << tmp_str_vec[0] << " IID=" << tmp_str_vec[1] << endl;
      exit(-1);
    }

    // read covariate data and check for missing values
    for(size_t i_cov = 0, j = 0; j < cov_colKeep.size(); j++) {

      if( !cov_colKeep[j] ) continue;

      pheno_data->new_cov(indiv_index, 1 + i_cov) = convertDouble(tmp_str_vec[2+j], params, sout);
      if( pheno_data->new_cov(indiv_index, 1 + i_cov) == params->missing_value_double ) {
        sout << "ERROR: Individual has missing value in covariate file: FID=" << tmp_str_vec[0] << " IID=" << tmp_str_vec[1] << endl;
        exit(-1);
      }
      i_cov++;
    }

  }
  myfile.close();

  // mask individuals in genotype data but not in covariate data
  pheno_data->masked_indivs.array().colwise() *= ind_in_cov_and_geno;

  // apply masking
  pheno_data->new_cov.array().colwise() *= ind_in_cov_and_geno.cast<double>();
  if(params->strict_mode) pheno_data->new_cov.array().colwise() *= pheno_data->masked_indivs.col(0).array().cast<double>();

  sout <<  endl;
  sout <<  "   -number of individuals with covariate data = " << ind_in_cov_and_geno.cast<int>().sum() << endl;

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

}


void residualize_phenotypes(struct param* params, struct phenodt* pheno_data, mstream& sout) {
  sout << "   -residualizing and scaling phenotypes...";
  auto t1 = std::chrono::high_resolution_clock::now();

  // residuals (centered) then scale
  MatrixXd beta = pheno_data->phenotypes.transpose() * pheno_data->new_cov;
  pheno_data->phenotypes -= ( (pheno_data->new_cov * beta.transpose()).array() * pheno_data->masked_indivs.array().cast<double>() ).matrix();
  pheno_data->scale_Y = pheno_data->phenotypes.colwise().norm().array() / sqrt(pheno_data->Neff.matrix().transpose().array() - 1);

  // check sd is not 0 
  if(pheno_data->scale_Y.minCoeff() < params->numtol){
    sout << "ERROR: At least one of the phenotypes has sd=0." << endl;
    exit(-1);
  }
  pheno_data->phenotypes.array().rowwise() /= pheno_data->scale_Y.array();

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;
}


double convertDouble(const string& phenoValue, struct param* params, mstream& sout){

  if(phenoValue == params->missing_pheno_str)
    return params->missing_value_double;

  double pheno_d;
  if(sscanf(phenoValue.c_str(), "%lf", &pheno_d) != 1){
    sout << "ERROR: Could not convert value to double: " << phenoValue << endl;
    exit(-1);
  }
  return pheno_d;
}



