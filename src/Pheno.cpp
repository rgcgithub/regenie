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

#include <unordered_set>
#include "Regenie.hpp"
#include "Files.hpp"
#include "Geno.hpp"
#include "Step1_Models.hpp"
#include "Pheno.hpp"

using namespace std;
using namespace Eigen;
using namespace boost;
using boost::math::normal;


void read_pheno_and_cov(struct in_files* files, struct param* params, struct filter* filters, struct phenodt* pheno_data, struct ests* m_ests, struct geno_block* gblock, mstream& sout) {

  ArrayXb ind_in_pheno_and_geno = ArrayXb::Constant( params->n_samples, false );
  ArrayXb ind_in_cov_and_geno = ArrayXb::Constant( params->n_samples, files->cov_file.empty());

  if(params->getCorMat){ // intitiate values for phenotype info

    params->n_pheno = 1;
    ind_in_pheno_and_geno = ArrayXb::Constant( params->n_samples, true );
    pheno_data->phenotypes = MatrixXd::Constant( params->n_samples, 1, 1);
    pheno_data->masked_indivs = MatrixXb::Constant( params->n_samples, 1, true);
    pheno_data->Neff = pheno_data->masked_indivs.cast<double>().colwise().sum();

  } else { // read in phenotype (mean-impute for QT)

    if(params->transposedPheno)
      tpheno_read(params, files, filters, pheno_data, ind_in_pheno_and_geno, sout);
    else
      pheno_read(params, files, filters, pheno_data, ind_in_pheno_and_geno, sout);

    if(params->trait_mode && !params->test_mode)
      m_ests->offset_nullreg = MatrixXd::Zero(params->n_samples, params->n_pheno); 

  }

  // used for step 2 if using firth and it failed
  set_pheno_pass(files, params);

  // Intercept
  pheno_data->new_cov = MatrixXd::Ones(params->n_samples, 1);

  // read in covariates
  if(!files->cov_file.empty()) covariate_read(params, files, filters, pheno_data, ind_in_cov_and_geno, sout);
  if(params->condition_snps)
      extract_condition_snps(params, files, filters, pheno_data, gblock, ind_in_cov_and_geno, sout);
  if(params->w_interaction){
    if(params->interaction_snp) // if doing GxG interaction
      extract_interaction_snp(params, files, filters, pheno_data, gblock, ind_in_cov_and_geno, sout);
    else if(params->interaction_prs) // if doing GxPRS interaction
      extract_interaction_prs(params, files, filters, pheno_data, ind_in_cov_and_geno, sout);
    if(params->gwas_condtl){ // append to new_cov
      pheno_data->new_cov.conservativeResize(pheno_data->new_cov.rows(), pheno_data->new_cov.cols() + pheno_data->interaction_cov.cols());
      pheno_data->new_cov.rightCols(pheno_data->interaction_cov.cols()) = pheno_data->interaction_cov;
      params->n_cov = pheno_data->new_cov.cols() - 1;// ignore intercept
    }
  }
  //cerr << endl<<pheno_data->new_cov.topRows(5) << endl;
  //if(params->w_interaction) cerr << endl<<pheno_data->interaction_cov.topRows(5)<< endl;

  // mask individuals 
  filters->ind_in_analysis = ind_in_pheno_and_geno && ind_in_cov_and_geno;
  setMasks(params, filters, pheno_data, sout);
  sout << " * number of individuals used in analysis = " << params->n_analyzed << endl;

  // apply rint
  if(params->rint) {
    sout << "   -applying RINT to all phenotypes\n";
    apply_rint(pheno_data, params);
  }

  // impute missing
  pheno_impute_miss(pheno_data, filters->ind_in_analysis, params);

  // print case-control counts per trait
  if(params->trait_mode==1)
    print_cc_info(params, files, pheno_data, sout);
  else
    print_info(params, files, pheno_data, sout);

}

void pheno_read(struct param* params, struct in_files* files, struct filter* filters, struct phenodt* pheno_data, Ref<ArrayXb> ind_in_pheno_and_geno, mstream& sout) {

  uint32_t indiv_index;
  bool all_miss;
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
  if( tmp_str_vec.size() < 2 ) 
    throw "header of phenotype file has too few columns.";
  else if( (tmp_str_vec[0] != "FID") || (tmp_str_vec[1] != "IID") ) 
    throw "header of phenotype file must start with: FID IID.";

  // check pheno with preds 
  if(params->test_mode && !params->getCorMat) check_blup(files, params, sout);

  // get phenotype names 
  keep_cols = ArrayXb::Constant(tmp_str_vec.size() - 2, true);
  for(int i = 0; i < keep_cols.size(); i++ ) {
    if(params->select_phenos_rm) // check if should skip phenotypes
      keep_cols(i) = !in_map(tmp_str_vec[i+2], filters->pheno_colRm_names);
    if(!keep_cols(i)) continue;
    if(params->select_phenos) // check if keeping pheno
      keep_cols(i) = in_map(tmp_str_vec[i+2], filters->pheno_colKeep_names);
    if(params->test_mode && !params->skip_blups && keep_cols(i)) // check phenotype had prs from step 1
      keep_cols(i) = has_blup(tmp_str_vec[i+2], files->blup_files, params, sout);

    if(keep_cols(i)) files->pheno_names.push_back( tmp_str_vec[i+2] );
  }
  params->n_pheno = keep_cols.count();

  // check #pheno
  if(params->n_pheno < 1)
    throw "need at least one phenotype.";

  sout << "n_pheno = " << params->n_pheno << endl;
  params->strict_mode |= (params->n_pheno == 1); // drop all missing observations

  // how missingness is handles
  if( params->strict_mode ) sout << "   -dropping observations with missing values at any of the phenotypes" << endl;
  else if( !params->rm_missing_qt  && (params->trait_mode==0)) sout << "   -keeping and mean-imputing missing observations (done for each trait)" << endl;


  // allocate memory
  pheno_data->phenotypes = MatrixXd::Zero(params->n_samples, params->n_pheno);
  pheno_data->masked_indivs = MatrixXb::Constant(params->n_samples, params->n_pheno, true);
  if(params->trait_mode)  
    pheno_data->phenotypes_raw = MatrixXd::Zero(params->n_samples, params->n_pheno);

  // read in data
  while( fClass.readLine(line) ){
    tmp_str_vec = string_split(line,"\t ");

    if( (int)tmp_str_vec.size() != (2+keep_cols.size()) )
      throw "incorrectly formatted phenotype file.";

    person = getIndivIndex(tmp_str_vec[0], tmp_str_vec[1], params, sout);
    if(!person.is_found) continue;

    indiv_index = person.index;

    // check duplicate
    if( !ind_in_pheno_and_geno(indiv_index) ){
      ind_in_pheno_and_geno( indiv_index ) = true;
    } else 
      throw "individual appears more than once in phenotype file: FID=" + tmp_str_vec[0] + " IID=" + tmp_str_vec[1] ;

    // read phenotypes 
    all_miss = true;
    for(int j = 0, i_pheno = 0; j < keep_cols.size(); j++) {

      if( !keep_cols(j) ) continue;

      pheno_data->phenotypes(indiv_index, i_pheno) = convertDouble(tmp_str_vec[2+j], params, sout);

      // for non-QT, save raw data
      if (params->trait_mode) {

        if((params->trait_mode==1) && !params->CC_ZeroOne && (pheno_data->phenotypes(indiv_index, i_pheno) != params->missing_value_double)) 
          pheno_data->phenotypes(indiv_index, i_pheno) -= 1; // if using 1/2/NA encoding for BTs

        pheno_data->phenotypes_raw(indiv_index, i_pheno) = pheno_data->phenotypes(indiv_index, i_pheno);

        // for BTs check 0/1/NA values
        if( (params->trait_mode==1) && (pheno_data->phenotypes_raw(indiv_index, i_pheno)!= 0) && 
            (pheno_data->phenotypes_raw(indiv_index, i_pheno)!= 1) ) {

          if(params->within_sample_l0)
            throw "no missing value allowed in phenotype file with option -within";
          else if( pheno_data->phenotypes_raw(indiv_index, i_pheno) != params->missing_value_double ) {
            std::string msg = (params->CC_ZeroOne ? "0/1/NA" : "1/2/NA");
            throw "a phenotype value is not " + msg + " for individual: FID=" + tmp_str_vec[0] + " IID=" + tmp_str_vec[1] + " Y=" + tmp_str_vec[2+j];
          }

          pheno_data->masked_indivs(indiv_index, i_pheno) = false;

        } else if( (params->trait_mode==2) && (pheno_data->phenotypes_raw(indiv_index, i_pheno)<0) ) { // CT check non-neg

          if(params->within_sample_l0)
            throw "no missing value allowed in phenotype file with option -within";
          else if( pheno_data->phenotypes_raw(indiv_index, i_pheno) != params->missing_value_double ) {
            throw "a phenotype value is <0 for individual: FID=" + tmp_str_vec[0] + " IID=" + tmp_str_vec[1] + " Y=" + tmp_str_vec[2+j];
          }
          pheno_data->masked_indivs(indiv_index, i_pheno) = false;

        }

      }

      if( pheno_data->phenotypes(indiv_index, i_pheno) != params->missing_value_double )
        all_miss = false;
      else {
        if( params->test_mode && params->rm_missing_qt ) pheno_data->masked_indivs(indiv_index, i_pheno) = false;
        if( params->strict_mode ) {
          pheno_data->masked_indivs.row(indiv_index) = MatrixXb::Constant(1, params->n_pheno, false);
          all_miss = true;
          break; // skip rest of the row
        }
      }

      i_pheno++;
    }

    if( all_miss ) ind_in_pheno_and_geno( indiv_index ) = false; // if individual has no phenotype data at all
  }

  // mask individuals in genotype data but not in phenotype data
  pheno_data->masked_indivs.array().colwise() *= ind_in_pheno_and_geno;

  // check if all individuals have missing/invalid phenotype
  int mInd;
  ArrayXi nobs_per_trait = pheno_data->masked_indivs.array().colwise().count().cast<int>();
  if((nobs_per_trait == 0).all())
    throw "all individuals have missing/invalid values for all traits." ;
  if(nobs_per_trait.minCoeff(&mInd) == 0)
    throw "all individuals have missing/invalid values for phenotype '" + files->pheno_names[mInd] + "'." ;

  // ignore traits with fewer than the specified minimum case count
  if(params->trait_mode==1)
    rm_phenoCols(ind_in_pheno_and_geno, files, params, pheno_data, sout); 

  // number of phenotyped individuals 
  sout <<  "   -number of phenotyped individuals " <<
   (params->strict_mode ? "with no missing data" : "" ) << 
   " = " << ind_in_pheno_and_geno.count() << endl;

  fClass.closeFile();

}

// in transposed format
void tpheno_read(struct param* params, struct in_files* files, struct filter* filters, struct phenodt* pheno_data, Ref<ArrayXb> ind_in_pheno_and_geno, mstream& sout) {

  uint32_t nid;
  string line, yname;
  std::vector< string > header, tmp_str_vec;
  map<int,uint32_t> indiv_index;
  map<int,uint32_t>::iterator itr;
  Files fClass;

  sout << left << std::setw(20) << " * phenotypes" << ": [" << files->pheno_file << "] ";
  fClass.openForRead(files->pheno_file, sout);
  fClass.readLine(line);
  check_str(line); // remove carriage returns at the end of line if any
  header = string_split(line,"\t ");

  // parse first line
  size_t ncols_file = header.size();
  for(size_t i=0; i < ncols_file; i++ ){
    if(in_map((int)(i+1), filters->tpheno_colrm)) continue;
    else if((i+1) == filters->tpheno_indexCol) continue;
    else { // get index of individuals in genotype file

      if(params->tpheno_iid_only) // assume FID=IID
        line = header[i] + "_" + header[i];
      else
        line = header[i];

      if (!in_map(line, params->FID_IID_to_ind)) continue;
      nid = params->FID_IID_to_ind[line];
      // check duplicate
      if( !ind_in_pheno_and_geno(nid) ){
        ind_in_pheno_and_geno( nid ) = true;
      } else 
        throw "individual appears more than once in phenotype file: ID=" + header[i];

      indiv_index[i] = nid;
    }

  }

  // check sample size
  if(indiv_index.size() == 0)
    throw "no individuals in phenotype file have genetic data.";

  // check pheno with preds 
  if(params->test_mode && !params->getCorMat) check_blup(files, params, sout);

  params->n_pheno = 0;
  int icol = 0;
  // for each trait
  while( fClass.readLine(line) ){

    tmp_str_vec = string_split(line,"\t ");
    if( tmp_str_vec.size() != ncols_file )
      throw "incorrectly formatted phenotype file.";

    // check trait name
    yname = tmp_str_vec[ filters->tpheno_indexCol - 1 ]; 
    if(params->select_phenos_rm && in_map(yname, filters->pheno_colRm_names)) continue;
    if(params->select_phenos && !in_map(yname, filters->pheno_colKeep_names)) continue;
    if(params->test_mode && !params->getCorMat && !has_blup(yname, files->blup_files, params, sout)) continue;
    files->pheno_names.push_back( yname );
    //cerr << params->n_pheno << " " << yname << endl;

    // resize matrices
    params->n_pheno++;
    pheno_data->phenotypes.conservativeResize(params->n_samples, params->n_pheno);
    pheno_data->phenotypes.rightCols(1).array() = 0;
    pheno_data->masked_indivs.conservativeResize(params->n_samples, params->n_pheno);
    pheno_data->masked_indivs.rightCols(1).array() = true;
    if(params->trait_mode) {
      pheno_data->phenotypes_raw.conservativeResize(params->n_samples, params->n_pheno);
      pheno_data->phenotypes_raw.rightCols(1).array() = 0;
    }

    // read in phenotype data
    for (itr = indiv_index.begin(); itr != indiv_index.end(); ++itr) {

      nid = itr->second;
      pheno_data->phenotypes(nid, icol) = convertDouble(tmp_str_vec[itr->first], params, sout);
      
      if (params->trait_mode) { // for nonQT, save raw data

        if((params->trait_mode==1) && !params->CC_ZeroOne && (pheno_data->phenotypes(nid, icol) != params->missing_value_double)) 
          pheno_data->phenotypes(nid, icol) -= 1; // if using 1/2/NA encoding for BTs

        pheno_data->phenotypes_raw(nid, icol) = pheno_data->phenotypes(nid, icol);

        if((params->trait_mode==1) &&  (pheno_data->phenotypes_raw(nid, icol)!= 0) && 
            (pheno_data->phenotypes_raw(nid, icol)!= 1) ) { // force 0/1/NA for BTs

          if(params->within_sample_l0)
            throw "no missing value allowed in phenotype file with option -within";
          else if( pheno_data->phenotypes_raw(nid, icol) != params->missing_value_double ){
            std::string msg = (params->CC_ZeroOne ? "0/1/NA" : "1/2/NA");
            throw "a phenotype value is not " + msg + " for individual: ID=" + header[itr->first] + " Y=" + tmp_str_vec[itr->first];
          }

          pheno_data->masked_indivs(nid, icol) = false;
        } else if((params->trait_mode==2) && (pheno_data->phenotypes_raw(nid, icol)<0) ) { // force non-neg for CTs

          if(params->within_sample_l0)
            throw "no missing value allowed in phenotype file with option -within";
          else if( pheno_data->phenotypes_raw(nid, icol) != params->missing_value_double ){
            throw "a phenotype value is <0 for individual: ID=" + header[itr->first] + " Y=" + tmp_str_vec[itr->first];
          }

          pheno_data->masked_indivs(nid, icol) = false;
        }
      }

      if(params->test_mode && params->rm_missing_qt && (pheno_data->phenotypes(nid, icol) == params->missing_value_double) ) 
        pheno_data->masked_indivs(nid, icol) = false;

    }
    icol++;
  }

  // check #pheno
  if(params->n_pheno < 1)
    throw "need at least one phenotype.";

  sout << "n_pheno = " << params->n_pheno << endl;

  params->strict_mode |= (params->n_pheno == 1); // drop all missing observations
  // how missingness is handled
  if( params->strict_mode ) {
    sout << "   -dropping observations with missing values at any of the phenotypes" << endl;
    pheno_data->masked_indivs.array().colwise() *= pheno_data->masked_indivs.array().rowwise().all();
  } else if( !params->rm_missing_qt  && (params->trait_mode==0)) 
    sout << "   -keeping and mean-imputing missing observations (done for each trait)" << endl;
  
  ind_in_pheno_and_geno = pheno_data->masked_indivs.array().rowwise().any(); // if individual has no phenotype data at all
  // mask individuals in genotype data but not in phenotype data
  pheno_data->masked_indivs.array().colwise() *= ind_in_pheno_and_geno;


  // check if all individuals have missing/invalid phenotype
  int mInd;
  ArrayXi nobs_per_trait = pheno_data->masked_indivs.array().colwise().count().cast<int>();
  if((nobs_per_trait == 0).all())
    throw "all individuals have missing/invalid values for all traits." ;
  if(nobs_per_trait.minCoeff(&mInd) == 0)
    throw "all individuals have missing/invalid values for phenotype '" + files->pheno_names[mInd] + "'." ;

  // ignore traits with fewer than the specified minimum case count
  if(params->trait_mode==1)
    rm_phenoCols(ind_in_pheno_and_geno, files, params, pheno_data, sout); 

  // number of phenotyped individuals 
  sout <<  "   -number of phenotyped individuals " <<
   (params->strict_mode ? "with no missing data" : "" ) << 
   " = " << ind_in_pheno_and_geno.count() << endl;

  fClass.closeFile();

}

// remove phenotypes with low case counts
void rm_phenoCols(Ref<ArrayXb> sample_keep, struct in_files* files, struct param* params, struct phenodt* pheno_data, mstream& sout) { 

  ArrayXb colrm = (pheno_data->phenotypes_raw.array() == 1).colwise().count() < params->mcc;
  int npass = (!colrm).count(); 
  //sout << npass << endl;

  if(npass == 0) 
    throw "all phenotypes have less than " + to_string( params->mcc ) +  " cases.";
  else if( colrm.count() == 0 ) 
    return;

  std::vector< string > tmp_str_vec;
  MatrixXd ynew (params->n_samples, npass);
  MatrixXb mnew (params->n_samples, npass);

  sout << "   -removing phenotypes with fewer than " << params->mcc << " cases\n";

  for(int i = 0, j = 0; i < params->n_pheno; i++ ) {

    if(colrm(i)) {
      sout << "    +WARNING: Phenotype '" << files->pheno_names[i] << "' has too few cases so it will be ignored.\n";
      continue;
    } 

    //cerr << j+1 << ":" << files->pheno_names[i] << endl;
    tmp_str_vec.push_back( files->pheno_names[i] );
    ynew.col(j) = pheno_data->phenotypes_raw.col(i);
    mnew.col(j) = pheno_data->masked_indivs.col(i);
    j++;

  }
  //sout << pheno_data->masked_indivs.colwise().count().array() << "\n\n" << mnew.colwise().count().array() << endl;

  files->pheno_names = tmp_str_vec;
  pheno_data->phenotypes = ynew;
  if(params->trait_mode) pheno_data->phenotypes_raw = ynew;
  pheno_data->masked_indivs = mnew;
  params->n_pheno = pheno_data->masked_indivs.cols();

  // remove samples with no phenotype values
  sample_keep = sample_keep && (pheno_data->masked_indivs.rowwise().count().array() > 0);
  sout << "    + n_pheno = " << npass << endl;

}

void covariate_read(struct param* params, struct in_files* files, struct filter* filters, struct phenodt* pheno_data, Ref<ArrayXb> ind_in_cov_and_geno, mstream& sout) {

  int nc_cat = 0, np_inter = 0;
  uint32_t indiv_index;
  ArrayXb keep_cols;
  ArrayXd inter_cov_column;
  MatrixXd inter_cov_matrix;
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
  if( (tmp_str_vec[0] != "FID") || (tmp_str_vec[1] != "IID") ) 
    throw "header of covariate file must start with: FID IID.";

  // get covariate names 
  keep_cols = ArrayXb::Constant(tmp_str_vec.size() - 2, true);
  for(int i = 0; i < keep_cols.size(); i++) {

    if(params->select_covs_rm) // check if should skip covariates
      keep_cols(i) = !in_map(tmp_str_vec[i+2], filters->cov_colRm_names);
    if(!keep_cols(i)) continue;

    if(!params->select_covs && !in_map(tmp_str_vec[i+2], filters->cov_colKeep_names)) // in case specified as categorical
      filters->cov_colKeep_names[tmp_str_vec[i+2]] = true;
    else keep_cols(i) = in_map(tmp_str_vec[i+2], filters->cov_colKeep_names);

    if(keep_cols(i)){
      covar_names.push_back( tmp_str_vec[i+2] );
      nc_cat += !filters->cov_colKeep_names[ tmp_str_vec[i+2] ];
      // with interaction test
      if(params->w_interaction && !params->interaction_snp && !params->interaction_prs && (filters->interaction_cov == tmp_str_vec[i+2]) ) {
        np_inter = 1;
        params->interaction_cat = !filters->cov_colKeep_names[ tmp_str_vec[i+2] ];
      }
    }
  }
  categories.resize(nc_cat);

  // check all covariates specified are in the file
  params->n_cov = keep_cols.count(); 
  if( (int)filters->cov_colKeep_names.size() != params->n_cov ) 
    throw "not all covariates specified are found in the covariate file.";

  if(params->w_interaction && !params->interaction_snp && !params->interaction_prs && (np_inter != 1))
    throw "cannot find the interaction covariate specified in the covariate file.";

  // check #covariates is > 0
  if(params->n_cov < 1){ // only intercept will be included
    sout << "n_cov = " << params->n_cov << " (+ intercept)" << endl;
    ind_in_cov_and_geno = true;
    return ;
  }
  sout << "n_cov = " << params->n_cov << endl;

  // allocate memory 
  pheno_data->new_cov = MatrixXd::Zero(params->n_samples, 1 + params->n_cov - np_inter);
  pheno_data->new_cov.col(0) = MatrixXd::Ones(params->n_samples, 1);
  if(params->w_interaction && !params->interaction_snp && !params->interaction_prs) inter_cov_column.resize(params->n_samples);

  // read in data
  while( fClass.readLine(line) ){
    tmp_str_vec = string_split(line,"\t ");

    if( (int)tmp_str_vec.size() != (keep_cols.size()+2) )
      throw "incorrectly formatted covariate file.";

    person = getIndivIndex(tmp_str_vec[0], tmp_str_vec[1], params, sout);
    if(!person.is_found) continue;

    indiv_index = person.index;

    // check duplicate
    if( !ind_in_cov_and_geno(indiv_index) )
      ind_in_cov_and_geno(indiv_index) = true;
    else 
      throw "individual appears more than once in covariate file: FID=" + tmp_str_vec[0] + " IID=" + tmp_str_vec[1];

    // read covariate data and check for missing values
    for(int i_cov = 0, i_col = 0, i_cat = 0, j = 0; j < keep_cols.size(); j++) {

      if( !keep_cols(j) ) continue;

      // interaction covariate
      if(params->w_interaction && !params->interaction_snp && !params->interaction_prs && (covar_names[i_cov] == filters->interaction_cov)){

        if( filters->cov_colKeep_names[ covar_names[i_cov] ] ) // if quantitative
          inter_cov_column(indiv_index) = convertDouble(tmp_str_vec[2+j], params, sout);
        else{ // set null category if interaction test and base level is specified
          if( (categories[i_cat].size() == 0) && filters->interaction_cov_null_level.size() > 0 )
            categories[i_cat][filters->interaction_cov_null_level] = 0;
          inter_cov_column(indiv_index) = convertNumLevel(tmp_str_vec[2+j], categories[i_cat++], params, sout);
        }

        if( inter_cov_column(indiv_index) == params->missing_value_double ) { // ignore individual
          ind_in_cov_and_geno(indiv_index) = false;
          break;
        }

      } else { // regular covariate

        if( filters->cov_colKeep_names[ covar_names[i_cov] ] ) // quantitative
          pheno_data->new_cov(indiv_index, 1 + i_col) = convertDouble(tmp_str_vec[2+j], params, sout);
        else // categorical so convert to numerical
          pheno_data->new_cov(indiv_index, 1 + i_col) = convertNumLevel(tmp_str_vec[2+j], categories[i_cat++], params, sout);

        if( pheno_data->new_cov(indiv_index, 1 + i_col) == params->missing_value_double ) { // ignore individual
          ind_in_cov_and_geno(indiv_index) = false;
          break;
        }

        i_col++;
      }

      i_cov++;
    }

  }
  //cerr << endl<<pheno_data->new_cov.block(0,0,3,pheno_data->new_cov.cols())<< endl;
  //if(params->w_interaction && !params->interaction_snp && !params->interaction_prs) cerr << inter_cov_column.head(3);

  // mask individuals in genotype data but not in covariate data
  pheno_data->new_cov.array().colwise() *= ind_in_cov_and_geno.cast<double>();
  if(inter_cov_column.size() > 0) inter_cov_column *= ind_in_cov_and_geno.cast<double>();

  // add dummy variables if needed
  if(nc_cat > 0){
    int n_dummies;
    int n_add = check_categories(covar_names, categories, params, filters, sout) - nc_cat + params->interaction_cat; // new columns to add (or remove if single category) & ignore interaction cov if categorical

    MatrixXd full_covarMat (pheno_data->new_cov.rows(), pheno_data->new_cov.cols() + n_add);
    // copy intercept column
    full_covarMat.col(0) = pheno_data->new_cov.col(0);

    for(int i = 0, raw_col = 1, full_col = 1, icat = -1; i < params->n_cov; i++){
      n_dummies = 1;

      if( filters->cov_colKeep_names[ covar_names[i] ] ) { // qCovar so copy column

        if(params->w_interaction && !params->interaction_snp && !params->interaction_prs && (covar_names[i] == filters->interaction_cov)) { 
          inter_cov_matrix = inter_cov_column.matrix(); continue; 
        } else
          full_covarMat.col(full_col) = pheno_data->new_cov.col(raw_col);

      } else { // cCovar

        icat++;

        if(params->w_interaction && !params->interaction_snp && !params->interaction_prs && (covar_names[i] == filters->interaction_cov)) { 

          np_inter = inter_cov_column.maxCoeff(); // get number of dummies to use

          if( np_inter == 0 ) // too few categories
            throw "interacting covariate '" + covar_names[i] + "' only has a single category.";
          else if((inter_cov_column>0).all())
            throw "no individual has baseline level specified for the interacting covariate (=" + filters->interaction_cov_null_level + ")";

          inter_cov_matrix = get_dummies(inter_cov_column);
          extract_names(params->interaction_lvl_names, categories[icat]); // save levels

          continue;

        } else {

          n_dummies = pheno_data->new_cov.col(raw_col).maxCoeff();
          if( n_dummies > 0 )
            full_covarMat.block(0, full_col, full_covarMat.rows(), n_dummies) = get_dummies(pheno_data->new_cov.col(raw_col).array());

        }

      }
      //cerr << i << " " << raw_col << " " << icat << " " << covar_names[i]  << endl;

      raw_col++;
      full_col += n_dummies;
    }

    pheno_data->new_cov = full_covarMat;
    params->n_cov = pheno_data->new_cov.cols() - 1; // ignore intercept
  }

  if(params->w_interaction && !params->interaction_snp && !params->interaction_prs) // save inter cov
    pheno_data->interaction_cov = filters->cov_colKeep_names[filters->interaction_cov] ? inter_cov_column.matrix() : inter_cov_matrix;

  sout <<  "   -number of individuals with covariate data = " << ind_in_cov_and_geno.count() << endl;
  if(params->w_interaction) {
    sout <<  "   -testing for interaction with "
      << (params->interaction_snp? "variant " : "")
      << "'" << filters->interaction_cov << "'\n";

    if((params->trait_mode==0) && !params->no_robust) {
      sout <<  "    +using " << 
        (params->force_robust && params->force_hc4? "HC4 robust SE" : "" ) << 
        (params->force_robust && !params->force_hc4? "HC3 robust SE" : "" ) << 
        (!params->force_robust ? "HLM model" : "" ) << 
        " when testing variants with MAC below " << params->rareMAC_inter << endl;
      if(!params->force_robust && !params->rint)
        sout <<  "    +WARNING: HLM should be used with RINTed traits (otherwise use option --apply-rint)\n"; 
    }

  }

  fClass.closeFile();

}

void setMasks(struct param* params, struct filter* filters, struct phenodt* pheno_data, mstream& sout){

  // mask samples
  if( params->strict_mode ) // keep if non-missing for all traits
    filters->ind_in_analysis = filters->ind_in_analysis && pheno_data->masked_indivs.array().rowwise().all();
  else // keep if non-missing for any trait
    filters->ind_in_analysis = filters->ind_in_analysis && pheno_data->masked_indivs.array().rowwise().any();

  // individuals kept in the analysis
  pheno_data->masked_indivs.array().colwise() *= filters->ind_in_analysis;

  // mask Y and X matrices
  pheno_data->phenotypes.array().colwise() *= filters->ind_in_analysis.cast<double>();
  if(params->trait_mode) 
    pheno_data->phenotypes_raw.array().colwise() *= filters->ind_in_analysis.cast<double>();
  pheno_data->new_cov.array().colwise() *= filters->ind_in_analysis.cast<double>();
  if( params->w_interaction ) 
    pheno_data->interaction_cov.array().colwise() *= filters->ind_in_analysis.cast<double>();

  // identify individuals masked for at least 1 trait
  filters->has_missing = !(pheno_data->masked_indivs.array().rowwise().all());
  //for(int i = 0; i <5; i++) cerr << std::boolalpha << filters->has_missing(i) << endl;

  // check sample size
  params->n_analyzed = filters->ind_in_analysis.count();
  if( params->n_analyzed < 1 ) 
    throw "sample size cannot be < 1.";
  pheno_data->Neff = pheno_data->masked_indivs.colwise().count().cast<double>();
  //sout << pheno_data->Neff << endl;

}


void print_cc_info(struct param* params, struct in_files* files, struct phenodt* pheno_data, mstream& sout){

  ArrayXd yvec;
  params->pheno_counts = MatrixXi::Constant(files->pheno_names.size(), 2, 0);

  // go through each trait and print number of cases and controls
  sout << " * case-control counts for each trait:\n";

  for (size_t i = 0; i < files->pheno_names.size(); i++){
    if( !params->pheno_pass(i) ) continue;

    params->pheno_counts(i, 0) = pheno_data->masked_indivs.col(i).select( pheno_data->phenotypes_raw.col(i).array(), 0).sum();
    params->pheno_counts(i, 1) = pheno_data->masked_indivs.col(i).select( 1 - pheno_data->phenotypes_raw.col(i).array(), 0).sum();
    sout << "   - '" << files->pheno_names[i] << "': " <<
      params->pheno_counts(i, 0) << " cases and " << params->pheno_counts(i, 1) << " controls\n";

  }
}

void print_info(struct param* params, struct in_files* files, struct phenodt* pheno_data, mstream& sout){

  params->pheno_counts = MatrixXi::Constant(files->pheno_names.size(), 2, 0);

  // go through each trait and print number of samples used
  sout << " * number of observations for each trait:\n";
  for (size_t i = 0; i < files->pheno_names.size(); i++)
    if( params->pheno_pass(i) ) {
      if((params->trait_mode == 0) && !params->force_qt_run) 
        check_nvals(i, files->pheno_names[i], params, pheno_data); // check there is not a binary trait
      params->pheno_counts(i, 0) = pheno_data->masked_indivs.col(i).count();
      sout << "   - '" << files->pheno_names[i] << "': " << params->pheno_counts(i, 0) << " observations\n";
    }
}

void check_nvals(int const& i_pheno, string const& pheno, struct param const* params, struct phenodt const* pheno_data){

  map<double, bool> uniq_vals;
  size_t n_min_vals = 10; // o.w. analyze as bt or ct

  for(size_t i = 0; i < params->n_samples; i++){
    if(!pheno_data->masked_indivs(i, i_pheno)) continue;
    if(!in_map(pheno_data->phenotypes(i, i_pheno), uniq_vals)) {
      uniq_vals[pheno_data->phenotypes(i, i_pheno)] = true;
      if(uniq_vals.size() > n_min_vals) return; // more than 2 values
    }
  }

  // if here then there are not more than 2 unique values
  throw "phenotype '" + pheno + "' has very few unique values (=" + to_string(uniq_vals.size()) + "). If you really want to analyze it as a QT, use `--force-qt`.";

}


void extract_interaction_snp(struct param* params, struct in_files* files, struct filter* filters, struct phenodt* pheno_data, struct geno_block* gblock, Ref<ArrayXb> ind_in_cov_and_geno, mstream& sout) {

  bool mean_impute = false;
  pheno_data->interaction_cov.resize(params->n_samples, 1);
  MapArXd Gcov (pheno_data->interaction_cov.col(0).data(), params->n_samples, 1);

  // read snp
  if(params->interaction_file) {// from external file
    Gcov = extract_from_genofile("interaction", mean_impute, ind_in_cov_and_geno, filters, files, params, sout).col(0).array();
  } else { // from input file
    read_snp(mean_impute, params->interaction_snp_offset, Gcov, ind_in_cov_and_geno, filters->ind_ignore, files, gblock->pgr, params, true);
    /*
       cerr << params->interaction_snp_offset << " " << ind_in_cov_and_geno.count() <<  "\n\n"
       << endl << pheno_data->interaction_cov.topRows(5)
       << endl; exit(-1);
       */
  }

  // apply coding
  code_snp(pheno_data->interaction_cov, ind_in_cov_and_geno, params->interaction_snp_offset, filters, files, params, sout);
  if(params->debug) cerr << pheno_data->interaction_cov.topRows(5) << endl;

}


void extract_condition_snps(struct param* params, struct in_files* files, struct filter* filters, struct phenodt* pheno_data, struct geno_block* gblock, Ref<ArrayXb> ind_in_cov_and_geno, mstream& sout) {

  bool mean_impute = true;
  int count = 0;
  std::map <std::string, uint64>::iterator itr;
  MatrixXd Gcov;

  if(params->condition_file) {

    sout << "    +conditioning on variants in [" << files->condition_snps_list << "]\n";
    Gcov = extract_from_genofile("conditional", mean_impute, ind_in_cov_and_geno, filters, files, params, sout);

  } else { // just read the snps

    sout << "    +conditioning on variants in [" << files->condition_snps_list << "] n_used = " << filters->condition_snp_names.size() << endl;

    Gcov.resize(params->n_samples, filters->condition_snp_names.size());
    for (itr = filters->condition_snp_names.begin(); itr != filters->condition_snp_names.end(); ++itr, count++) 
      read_snp(mean_impute, itr->second, Gcov.col(count).array(), ind_in_cov_and_geno, filters->ind_ignore, files, gblock->pgr, params, true);

  }
  
  // Add to covariates
  pheno_data->new_cov.conservativeResize( pheno_data->new_cov.rows(), pheno_data->new_cov.cols() + Gcov.cols());
  pheno_data->new_cov.rightCols(Gcov.cols()) = Gcov;

  //cerr << Gcov.topRows(5) << "\n\n" << pheno_data->new_cov.topRows(5) << "\n\n";

}


int check_categories(vector<std::string>& covar, vector<std::map<std::string,int>>& categories, struct param* params, struct filter* filters, mstream& sout){

  int ntotal = 0, n_levels = 0;

  for(size_t i = 0, j = 0; i < covar.size(); i++){

    // skip qCovar
    if( filters->cov_colKeep_names[ covar[i] ] ) continue;

    n_levels = categories[j++].size();

    if( n_levels > params->max_cat_levels) // too many categories
      throw "too many categories for covariate: " + covar[i] + " (=" + to_string( n_levels ) + "). Either use '--maxCatLevels' or combine categories.";
    else if( n_levels == 1) // too few categories
      sout << "WARNING: covariate ' " << covar[i] << "' only has a single category so it will be ignored\n";

    // for interaction test
    if(params->w_interaction && !params->interaction_snp && (covar[i] == filters->interaction_cov)) // skip it
      continue;
     
    ntotal += n_levels - 1; // add K-1 dummy vars
  }

  return ntotal;
}

MatrixXd get_dummies(const Eigen::Ref<const Eigen::ArrayXd>& numCov) {

  int index, nvars = numCov.maxCoeff();
  MatrixXd dout = MatrixXd::Zero(numCov.size(), nvars);
  //cerr << dout.rows() << "***" << dout.cols() << "--min=" << numCov.minCoeff() << endl;

  for(int i = 0; i < numCov.size(); i++){
    if(numCov(i) == 0) continue; // will go to intercept

    index = numCov(i) - 1;
    dout(i, index) = 1;
  }

  return dout;

}

// check if need to add E^2 (i.e. if single column and not 0/1)
bool add_square_term(const Eigen::Ref<const Eigen::MatrixXd>& X) {

  if(X.cols() > 1) // categorical
    return false;

  std::unordered_set<double> vals(X.col(0).data(), X.col(0).data() + X.rows());
  if(vals.size() > 2) // not dichotomous
    return true;

  if(vals.find(0) != vals.end()) // one of the values is 0
    return false;

  return true;
}

void extract_names(vector<string>& names, map<string,int>& map_names){

  map<string, int >::iterator itr;
  names.resize( map_names.size() - 1); // ignore 0 category

  for (itr = map_names.begin(); itr != map_names.end(); ++itr) {
    if(itr->second == 0) continue; // ignore 0 category
    names[ itr->second - 1 ] = itr->first; // map_names has values in 0-{K-1}
  }
  //for(auto i:names) cerr << i << "\n";

}

// Adjust for covariates (incl. intercept)
// in step 2, also read blups and check
void prep_run (struct in_files* files, struct filter* filters, struct param* params, struct phenodt* pheno_data, struct ests* m_ests, mstream& sout){

  // for step 2, check blup files
  if (params->test_mode && !params->getCorMat){
    // individuals not in blup file will have their phenotypes masked
    blup_read(files, params, pheno_data, m_ests, sout);
    if(params->write_samples) write_ids(files, params, pheno_data, sout);
  }

  // compute N for each trait
  setMasks(params, filters, pheno_data, sout);

  // for interaction test with BTs, add E^2 to covs
  if( (params->trait_mode==1) && params->w_interaction && params->gwas_condtl ) {
    pheno_data->new_cov.conservativeResize( pheno_data->new_cov.rows(), pheno_data->new_cov.cols() + pheno_data->interaction_cov.cols());
    pheno_data->new_cov.rightCols(pheno_data->interaction_cov.cols()) = pheno_data->interaction_cov.array().square().matrix();
  }

  // orthonormal basis (save number of lin. indep. covars.)
  params->ncov = getBasis(pheno_data->new_cov, params);
  if(params->ncov > (int)params->n_samples)
    throw "number of covariates is larger than sample size!";

  // compute offset for nonQT (only in step 1)
  if((params->trait_mode==1) && !params->test_mode) fit_null_logistic(false, 0, params, pheno_data, m_ests, files, sout);
  else if((params->trait_mode==2) && !params->test_mode) fit_null_poisson(0, params, pheno_data, m_ests, files, sout);

  // with interaction test, remove colinear columns
  if( params->w_interaction ) {
    // apply QR decomp
    QRcheck(pheno_data->interaction_cov, params);
    //cerr << pheno_data->interaction_cov.topRows(3) << "\n\n";
    params->ncov_interaction = pheno_data->interaction_cov.cols();
    params->n_tests_per_variant += 3; // marginal + inter + joint
    params->int_add_extra_term = add_square_term(pheno_data->interaction_cov);
    params->add_homdev = params->add_homdev && params->int_add_extra_term && (pheno_data->interaction_cov.array()>=1.5).any();

    if(!params->gwas_condtl) { // include main effects of Xinter

      params->int_add_esq = (params->trait_mode==1) && !params->add_homdev && params->int_add_extra_term;
      if(params->int_add_esq){ // use G_E and G_E^2
        params->interaction_istart = 2 * params->ncov_interaction;
        pheno_data->interaction_cov_res.resize(pheno_data->interaction_cov.rows(), pheno_data->interaction_cov.cols() * 2);
        pheno_data->interaction_cov_res << pheno_data->interaction_cov, pheno_data->interaction_cov.array().square().matrix();
      } else if(params->add_homdev){ // only with additive coding (add hom. correction term)
        pheno_data->interaction_homdev = (pheno_data->interaction_cov.array()>=1.5).cast<double>().matrix();
        params->interaction_istart = 2 * params->ncov_interaction;
        pheno_data->interaction_cov_res.resize(pheno_data->interaction_cov.rows(), pheno_data->interaction_cov.cols() * 2);
        pheno_data->interaction_cov_res << pheno_data->interaction_cov, pheno_data->interaction_homdev;
      } else { // use only G_E
        params->interaction_istart = params->ncov_interaction;
        // keep original and residualized version
        pheno_data->interaction_cov_res = pheno_data->interaction_cov;
      }
      //cerr << pheno_data->interaction_cov_res.topRows(5) << "\n\n";

      // remove covariate effects
      if(!residualize_matrix(pheno_data->interaction_cov_res, pheno_data->scl_inter_X, pheno_data->new_cov.leftCols( params->ncov + (params->blup_cov && (params->trait_mode == 1) ? -1 : 0)), params->n_analyzed, params->numtol))
        throw "Var=0 for the interaction risk factor.";

    }
    //cerr << pheno_data->interaction_cov.topRows(3) << "\n\n";

    // for interaction tests with QTs using HLM - keep raw Y
    if(params->trait_mode==0)
      pheno_data->phenotypes_raw = pheno_data->phenotypes;

    // allocate per thread if using OpenMP
    pheno_data->Hmat.resize(params->neff_threads);
    pheno_data->scf_i.resize(params->neff_threads);
    
  }

  // residualize phenotypes (skipped for nonQTs when testing)
  if( !params->getCorMat && (!params->test_mode || (params->trait_mode==0)) ) 
    residualize_phenotypes(params, pheno_data, files->pheno_names, sout);

  // if using step 1 preds as covariate
  check_cov_blup(pheno_data, params);

  // store indices for ADAM
  if(params->use_adam && params->adam_mini){
    params->adam_indices.resize(params->n_pheno);
    for(int ph = 0; ph < params->n_pheno; ph++){
      if( !params->pheno_pass(ph) ) continue;
      params->adam_indices[ph].resize(pheno_data->masked_indivs.col(ph).count(),1);
      for(size_t i = 0, j = 0; i < params->n_samples; i++)
        if(pheno_data->masked_indivs(i,ph)) params->adam_indices[ph](j++) = i;
    }
  }
}

// get list of phenotypes in pred file
void check_blup(struct in_files* files, struct param* params, mstream& sout) {

  string line, tmp_pheno;
  std::vector< string > tmp_str_vec;
  Files fClass;

  // skip reading if specified by user
  if( params->skip_blups && !params->interaction_prs ) return;

  fClass.openForRead(files->blup_list_file, sout);

  while (fClass.readLine(line)){
    tmp_str_vec = string_split(line,"\t ");

    // each line contains a phenotype name and the corresponding blup file name
    if( tmp_str_vec.size() != 2 )
      throw "step 1 list file is not in the right format : " + files->blup_list_file;

    if(in_map(tmp_str_vec[0], files->blup_files))
      throw "phenotype \'" + tmp_str_vec[0] + "\' appears more than once in step 1 list file.";

    files->blup_files[ tmp_str_vec[0] ] = tmp_str_vec[1];
  }

  fClass.closeFile();
}

bool has_blup(string const& yname, map<string,string> const& y_read, struct param const* params, mstream& sout) {

  if( ( params->skip_blups && !params->interaction_prs ) || in_map(yname, y_read)) return true;

  sout << "WARNING: No step 1 file provided for phenotype '" << yname << "' so it will be ignored.\n";
  return false;

}

// get list of blup files
void blup_read(struct in_files* files, struct param* params, struct phenodt* pheno_data, struct ests* m_ests, mstream& sout) {

  int n_masked_prior, n_masked_post;
  uint32_t indiv_index;
  double blup_val;
  string yfile, line;
  std::vector< string > tmp_str_vec, tmp_prs_vec;
  ArrayXd full_prs;
  ArrayXb blupf_mask, all_miss_pheno;
  Files fClass;

  // allocate memory
  m_ests->blups = MatrixXd::Zero(params->n_samples, params->n_pheno);

  // skip reading if specified by user
  if( params->skip_blups && !params->interaction_prs) {
    string mode;
    if(params->trait_mode==0) mode = "linear";
    else if(params->trait_mode==1) mode = "logistic";
    else if(params->trait_mode==2) mode = "poisson";
      sout << " * no step 1 predictions given. Simple " << mode << " regression will be performed" <<endl;
    return;
  } else if(params->interaction_prs) return;

  sout << " * " << (params->use_prs ? "PRS" : "LOCO") << " predictions : [" << files->blup_list_file << "]\n";
  all_miss_pheno = ArrayXb::Constant(params->n_pheno, false);

  // allocate memory for LTCO 
  if(params->w_ltco) m_ests->ltco_prs = MatrixXd::Zero(params->n_samples, params->n_pheno);

  // read blup file for each phenotype
  for(int ph = 0; ph < params->n_pheno; ph++) {

    if( !params->pheno_pass(ph) ) continue;

    yfile = files->blup_files[ files->pheno_names[ph] ];
    sout << "   -file [" << yfile  << "] for phenotype '" << files->pheno_names[ph] << "'\n";

    fClass.openForRead(yfile, sout);

    // to mask all individuals not present in .loco file
    blupf_mask = ArrayXb::Constant(params->n_samples, false);
    n_masked_prior = pheno_data->masked_indivs.col(ph).count();

    // read first line which has FID_IID
    fClass.readLine(line);
    tmp_str_vec = string_split(line,"\t ");

    if( tmp_str_vec[0] != "FID_IID") 
      throw "header of blup file must start with FID_IID (=" + tmp_str_vec[0] + ")";

    // read second line to check for missing predictions
    fClass.readLine(line);
    tmp_prs_vec = string_split(line,"\t ");

    if( params->use_prs && (tmp_prs_vec[0] != "0") )
      throw "second line must start with 0 (=" + tmp_prs_vec[0] + ").";

    for (size_t i = 1; i < tmp_str_vec.size(); i++){
      // ignore sample if it is not in genotype data
      if (!in_map(tmp_str_vec[i], params->FID_IID_to_ind)) continue;
      indiv_index = params->FID_IID_to_ind[tmp_str_vec[i]];
      blup_val = convertDouble(tmp_prs_vec[i], params, sout);

      // ignore samples where prediction is NA
      blupf_mask( indiv_index ) = (blup_val != params->missing_value_double);
      //cerr << tmp_str_vec[i] << "\t" << std::boolalpha << blupf_mask( indiv_index ) << endl; 
      if (!blupf_mask( indiv_index )) continue;
      
      if( params->use_prs ) 
        m_ests->blups(indiv_index, ph) = blup_val;
    }

    // mask samples not in file
    pheno_data->masked_indivs.col(ph).array() = pheno_data->masked_indivs.col(ph).array() && blupf_mask;
    n_masked_post = pheno_data->masked_indivs.col(ph).count();

    // check not everyone is masked
    all_miss_pheno(ph) = n_masked_post < 1;
    if( all_miss_pheno(ph) ) {
      fClass.closeFile();
      continue;
    }

    if( n_masked_post < n_masked_prior ){
      sout << "    + " << n_masked_prior - n_masked_post <<
        " individuals with missing LOCO predictions will be ignored for the trait\n";
    }

    if(params->w_ltco){ // go through each line and sum up the loco prs

      bool chr_ltco;
      int nchr_file = 0;
      double ds;
      full_prs = ArrayXd::Zero( params->n_samples );

      // Re-open file (since skipped 2nd row)
      fClass.closeFile();
      fClass.openForRead(yfile, sout);
      fClass.ignoreLines(1); // skip first row

      while( fClass.readLine(line) ){

        tmp_prs_vec = string_split(line,"\t ");
        if( tmp_prs_vec.size() != tmp_str_vec.size() )
          throw "number of entries for chromosome " + tmp_prs_vec[0] + 
            " does not match with that in header (" + 
            to_string(tmp_prs_vec.size()) + " vs " + to_string(tmp_str_vec.size()) + ")";

        // check if it is the LTCO chromosome
        chr_ltco = (chrStrToInt(tmp_prs_vec[0], params->nChrom) == params->ltco_chr);

        for (size_t i = 1; i < tmp_str_vec.size(); i++){
          // ignore sample if it is not in genotype data
          if (!in_map(tmp_str_vec[i], params->FID_IID_to_ind)) continue;
          indiv_index = params->FID_IID_to_ind[tmp_str_vec[i]];
          if(!pheno_data->masked_indivs(indiv_index,ph)) continue;

          ds = convertDouble(tmp_prs_vec[i], params, sout);
          if(chr_ltco) m_ests->ltco_prs(indiv_index, ph) = - ds;
          full_prs(indiv_index) += ds;
        }

        nchr_file++;
      }

      if( nchr_file != params->nChrom )
        throw "incorrectly formatted file.";

      m_ests->ltco_prs.col(ph).array() += full_prs / (nchr_file - 1);
      //cerr << m_ests->ltco_prs.col(ph).head(5)<<endl;

    }

    fClass.closeFile();
  }

    check_phenos(all_miss_pheno, files->pheno_names, files->out_file + "_" + "pheno_all_miss.txt", sout);

}

void extract_interaction_prs(struct param* params, struct in_files* files, struct filter* filters, struct phenodt* pheno_data, Ref<ArrayXb> ind_in_cov_and_geno, mstream& sout) {

  pheno_data->interaction_cov.resize(params->n_samples, 1);
  MapArXd PRS (pheno_data->interaction_cov.col(0).data(), params->n_samples, 1);

  // read prs
  read_prs(PRS, files, params, ind_in_cov_and_geno, sout);
  if(params->debug) cerr << "full PRS head:" << PRS.matrix().transpose().array().head(5) << endl;

  // disable adjusting for PRS
  params->use_prs = false;
  params->skip_blups = true;

}

void read_prs(Ref<ArrayXd> full_prs, struct in_files* files, struct param* params, Ref<ArrayXb> ind_in_cov_and_geno, mstream& sout) {

  int n_masked_prior, n_masked_post;
  uint32_t indiv_index;
  double blup_val;
  string yfile, line;
  std::vector< string > tmp_str_vec, tmp_prs_vec;
  ArrayXb blupf_mask;
  Files fClass;

  sout << "    + extracting full PRS using " << (params->use_prs ? "PRS" : "LOCO") << " predictions : [" << files->blup_list_file << "]\n";
  if(params->n_pheno > 1) throw "option '--interaction-prs' only works with a single phenotype"; 

  // read blup file for each phenotype
  for(int ph = 0; ph < params->n_pheno; ph++) {
    if( !params->pheno_pass(ph) ) continue;

    yfile = files->blup_files[ files->pheno_names[ph] ];
    sout << "    -file [" << yfile  << "] for phenotype '" << files->pheno_names[ph] << "'\n";

    fClass.openForRead(yfile, sout);

    // to mask all individuals not present in .loco file
    blupf_mask = ArrayXb::Constant(params->n_samples, false);
    n_masked_prior = ind_in_cov_and_geno.count();
    full_prs = 0;

    // read first line which has FID_IID
    fClass.readLine(line);
    tmp_str_vec = string_split(line,"\t ");

    if( tmp_str_vec[0] != "FID_IID") 
      throw "header of blup file must start with FID_IID (=" + tmp_str_vec[0] + ")";

    // read second line to check for missing predictions
    fClass.readLine(line);
    tmp_prs_vec = string_split(line,"\t ");

    if( params->use_prs && (tmp_prs_vec[0] != "0") )
      throw "second line must start with 0 (=" + tmp_prs_vec[0] + ").";

    for (size_t i = 1; i < tmp_str_vec.size(); i++){
      // ignore sample if it is not in genotype data
      if (!in_map(tmp_str_vec[i], params->FID_IID_to_ind)) continue;
      indiv_index = params->FID_IID_to_ind[tmp_str_vec[i]];
      blup_val = convertDouble(tmp_prs_vec[i], params, sout);

      // ignore samples where prediction is NA
      blupf_mask( indiv_index ) = (blup_val != params->missing_value_double);
      //cerr << tmp_str_vec[i] << "\t" << std::boolalpha << blupf_mask( indiv_index ) << endl; 
      if (!blupf_mask( indiv_index )) continue;

      if( params->use_prs ) full_prs(indiv_index) = blup_val;
    }

    // mask samples not in file
    ind_in_cov_and_geno = ind_in_cov_and_geno && blupf_mask;
    n_masked_post = ind_in_cov_and_geno.count();

    // check not everyone is masked
    if( n_masked_post < 1 )
      throw "none of the samples have step 1 predictions in file."; 

    if( n_masked_post < n_masked_prior ){
      sout << "     * " << n_masked_prior - n_masked_post <<
        " individuals with missing LOCO predictions will be ignored for the trait\n";
    }

    if(!params->use_prs){ // go through each line and sum up the loco prs to get the full PRS

      int nchr_file = 0;
      double ds;

      // Re-open file (since skipped 2nd row)
      fClass.closeFile();
      fClass.openForRead(yfile, sout);
      fClass.ignoreLines(1); // skip first row

      while( fClass.readLine(line) ){

        tmp_prs_vec = string_split(line,"\t ");
        if( tmp_prs_vec.size() != tmp_str_vec.size() )
          throw "number of entries for chromosome " + tmp_prs_vec[0] + 
            " does not match with that in header (" + 
            to_string(tmp_prs_vec.size()) + " vs " + to_string(tmp_str_vec.size()) + ")";

        for (size_t i = 1; i < tmp_str_vec.size(); i++){
          // ignore sample if it is not in genotype data
          if (!in_map(tmp_str_vec[i], params->FID_IID_to_ind)) continue;
          indiv_index = params->FID_IID_to_ind[tmp_str_vec[i]];
          if(!ind_in_cov_and_geno(indiv_index)) continue;

          ds = convertDouble(tmp_prs_vec[i], params, sout);
          full_prs(indiv_index) += ds;
        }

        nchr_file++;
      }

      if( nchr_file != params->nChrom )
        throw "incorrectly formatted file (not enough chromosomes).";

      full_prs /= (nchr_file - 1);

    }

    fClass.closeFile();
  }

}

void check_phenos(Ref<ArrayXb> pheno_ind, vector<string> const& pheno_names, string const& fname, mstream& sout){

  // print phenotype names where ind is true and exit
  if(!pheno_ind.any()) return;

  Files fout;
  fout.openForWrite(fname, sout);

  for(int i = 0; i < pheno_ind.size(); i++)
    if(pheno_ind(i))
      fout << pheno_names[i] << "\n";
  fout.closeFile();

  throw "Problematic phenotypes found (all individuals have missing step 1 predictions); names written to [" + fname + "].\n" +
    "You can use `--phenoExcludeList` to ignore these traits.";

}


// write ids of samples included in step 2 (done for each trait)
void write_ids(struct in_files const* files, struct param* params, struct phenodt const* pheno_data, mstream& sout){

  uint32_t index;
  map<string, uint32_t >::iterator itr_ind;
  string idfile;
  Files fout;

    sout << " * user specified to write sample IDs for each trait"<<endl;

  for( int ph = 0; ph < params->n_pheno; ph++){

    if( !params->pheno_pass(ph) ) continue;

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


int getBasis(MatrixXd& X,struct param const* params){

  // eigen-decompose NxK matrix
  MatrixXd xtx = X.transpose() * X;
  SelfAdjointEigenSolver<MatrixXd> es(xtx);
  VectorXd D = es.eigenvalues();
  MatrixXd V = es.eigenvectors();

  // create basis set
  // eigenvalues sorted in increasing order
  int non_zero_eigen = (D.array() > D.tail(1)(0) * params->eigen_val_rel_tol).count();
  RowVectorXd vv1 = D.tail(non_zero_eigen).array().sqrt();
  X *= V.rightCols(non_zero_eigen);
  X.array().rowwise() /= vv1.array();

  return non_zero_eigen;
}

void QRcheck(MatrixXd& mat, struct param* params){

  vector<string> new_names;

  // find set of linearly independent cols
  ColPivHouseholderQR<MatrixXd> qrA(mat);
  qrA.setThreshold(params->qr_tol); 
  int indCols = qrA.rank();

  if(indCols == 0)
    throw "rank of matrix is 0.";
  else if ( indCols < mat.cols() ){
    ArrayXi colKeep = qrA.colsPermutation().indices();
    std::vector<int> new_indices;

    // keep only linearly independent columns
    MatrixXd tmpM (mat.rows(), indCols);

    for(int i = 0; i < indCols; i++){
      tmpM.col(i) = mat.col( colKeep(i) );
      if(params->interaction_cat) new_names.push_back( params->interaction_lvl_names[ colKeep(i) ]);
    }

    mat = tmpM;
    if(params->interaction_cat) params->interaction_lvl_names = new_names;
  } 

  // check no columns has sd = 0
  check_sd(mat, params->n_analyzed - params->ncov, params->numtol);

}

void check_sd(const Eigen::Ref<const Eigen::MatrixXd>& mat, int const& n, double const& numtol){

  RowVectorXd mu = mat.colwise().mean();
  VectorXd sd = (mat.rowwise() - mu).colwise().norm().array() / sqrt(n);

  if(sd.minCoeff() < numtol)
    throw "one of the columns has sd=0";

}

// only in step 2
void check_cov_blup(struct phenodt* pheno_data, struct param* params) {

  params->ncov_analyzed = params->ncov + (int) params->blup_cov;

  // for BTs, add extra column for LOCO PRS
  if((params->trait_mode==1) && params->blup_cov){
    pheno_data->new_cov.conservativeResize( pheno_data->new_cov.rows(), pheno_data->new_cov.cols() + 1);
    params->ncov = pheno_data->new_cov.cols();
    params->ncov_analyzed = params->ncov;
  }

}

void residualize_phenotypes(struct param const* params, struct phenodt* pheno_data, const std::vector<std::string>& pheno_names, mstream& sout) {
  sout << "   -residualizing and scaling phenotypes...";
  auto t1 = std::chrono::high_resolution_clock::now();

  // residuals (centered) then scale
  MatrixXd beta = pheno_data->phenotypes.transpose() * pheno_data->new_cov;
  pheno_data->phenotypes -= ( (pheno_data->new_cov * beta.transpose()).array() * pheno_data->masked_indivs.array().cast<double>() ).matrix();
  pheno_data->scale_Y = pheno_data->phenotypes.colwise().norm().array() / sqrt(pheno_data->Neff.matrix().transpose().array() - params->ncov);

  // set sd for phenotypes which are ignored to 1
  pheno_data->scale_Y = params->pheno_pass.select(pheno_data->scale_Y.transpose().array(), 1).matrix().transpose();

  // check sd is not 0 
  MatrixXd::Index minIndex;
  if(pheno_data->scale_Y.minCoeff(&minIndex) < params->numtol)
    throw "phenotype \'" + pheno_names[minIndex] + "\' has sd=0.";

  pheno_data->phenotypes.array().rowwise() /= pheno_data->scale_Y.array();

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;
}

bool residualize_matrix(MatrixXd& mat, ArrayXd& scf, const Eigen::Ref<const Eigen::MatrixXd>& X, size_t const& N, double const& numtol) {

  // residuals (centered) 
  MatrixXd beta = mat.transpose() * X;
  mat -= X * beta.transpose();

  scf = mat.colwise().norm().array() / sqrt(N - X.cols());

  // check sd is not 0 
  if(scf.minCoeff() < numtol)
    return false;

  // scale
  mat.array().rowwise() /= scf.matrix().transpose().array();

  return true;
}

void get_lm_resid(MatrixXd& res, const Eigen::Ref<const Eigen::MatrixXd>& step1_preds, const Eigen::Ref<const Eigen::MatrixXd>& Y) {

  ArrayXd beta = (1 / step1_preds.colwise().squaredNorm().array()) * (step1_preds.array() * Y.array()).colwise().sum();
  res = Y - step1_preds * beta.matrix().asDiagonal();

}

void apply_QR(MatrixXd& mat, struct param const* params, bool const& scale){

  // find set of linearly independent cols
  ColPivHouseholderQR<MatrixXd> qrA(mat);
  qrA.setThreshold(params->qr_tol); 
  int indCols = qrA.rank();

  if(indCols == 0)
    throw "rank of matrix is 0.";
  else if ( indCols < mat.cols() ){
    ArrayXi colKeep = qrA.colsPermutation().indices();
    std::vector<int> new_indices;

    // keep only linearly independent columns
    MatrixXd tmpM (mat.rows(), indCols);

    for(int i = 0; i < indCols; i++)
      tmpM.col(i) = mat.col( colKeep(i) );

    mat = tmpM;
  } 

  if(scale)
    rescale_mat(mat, params);

}

void rescale_mat(Ref<MatrixXd> mat, struct param const* params){

    RowVectorXd mu = mat.colwise().sum() / params->n_samples;
    mat.rowwise() -= mu;

    // check sd is not 0 
    ArrayXd scf = mat.colwise().norm().array() / sqrt(params->n_samples - 1);
    if(scf.minCoeff() < params->numtol)
      throw "sd = 0 occurred";

    // scale
    mat.array().rowwise() /= scf.matrix().transpose().array();

}

void check_str(string& mystring ){

  // check there are no '\r' at the end
  if( !mystring.empty() && mystring[mystring.size() - 1] == '\r' )
    mystring.erase(mystring.size() - 1);

}

void pheno_impute_miss(struct phenodt* pheno_data, const Eigen::Ref<const ArrayXb>& ind_in_analysis, struct param const* params){

  if((params->trait_mode==0) || !params->test_mode){
    double total, ns;

    // for each trait, impute missing with mean
    for(int j = 0; j < params->n_pheno; j++)
      if( params->pheno_pass(j) ){

        MapArXd Y (pheno_data->phenotypes.col(j).data(), params->n_samples, 1);
        MapArXb mask (pheno_data->masked_indivs.col(j).data(), params->n_samples, 1);

        if(params->trait_mode==0){ // impute missing with mean
          total = ( Y != params->missing_value_double ).select(Y, 0).sum();
          ns = ( ind_in_analysis && (Y != params->missing_value_double) ).count();
          Y = ( Y != params->missing_value_double ).select(Y, total / ns);
        } else { // mask tracks missingness
          total = mask.select(Y, 0).sum() / mask.count();
          Y = mask.select(Y, total);
        }
      }

    // apply masking
    pheno_data->phenotypes.array() *= pheno_data->masked_indivs.array().cast<double>();

  }

}

void apply_rint(struct phenodt* pheno_data, struct param const* params){

  // for each trait, apply rank-inverse normal transformation
  for(int ph = 0; ph < params->n_pheno; ph++)
    if( params->pheno_pass(ph) )
      rint_pheno(pheno_data->phenotypes.col(ph), (pheno_data->phenotypes.col(ph).array() != params->missing_value_double) && pheno_data->masked_indivs.col(ph).array());

}

void set_pheno_pass(struct in_files const* files, struct param* params){

  bool select_phenos = params->select_pheno_l1.size() > 0;
  params->pheno_pass = ArrayXb::Constant(params->n_pheno, false);
  params->pheno_fail_nullreg = ArrayXb::Constant(params->n_pheno, false);
  for(int ph = 0; ph < params->n_pheno; ph++)
    if( select_phenos )
      params->pheno_pass(ph) = in_map( files->pheno_names[ph], params->select_pheno_l1 );
    else
      params->pheno_pass(ph) = true;

  // sanity check
  if((!params->pheno_pass).all())
    throw "none of the specified phenotypes for level 1 were found.\n";

}

void rint_pheno(Ref<MatrixXd> Y, const Eigen::Ref<const ArrayXb>& mask){

  int nvals = mask.count();
  vector<rank_pair> yvals;
  yvals.resize(nvals);

  // get the index for each value
  for(int i = 0, j = 0; i < Y.rows(); i++){
    if(!mask(i)) continue;
    yvals[j].val = Y(i,0);
    yvals[j].index = i;
    j++;
  }

  // sort by values keeping track of index
  std::sort(yvals.begin(), yvals.end(), cmp_rank_pair);

  // take care of ties
  int n_eq;
  for(int i = 0; i < nvals; i+=n_eq){
    n_eq = 1;
    while(((i+n_eq) < nvals) && (yvals[i+n_eq].val == yvals[i].val)) n_eq++;
    for(int j = 0; j < n_eq; j++)
      yvals[i+j].val = (i+1) + (n_eq-1)/2.0;
  }

  // apply INT with the ranks
  double kc = 3/8.0, rint_val;
  normal nd(0,1);
  //cerr << Y.block(0,0,6,1)<<endl;
  for(auto const& ypair : yvals){
    rint_val = (ypair.val - kc) / (nvals - 2 * kc + 1);
    Y( ypair.index, 0 ) = quantile(nd, rint_val);
  }
  //cerr << endl << endl << Y.block(0,0,6,1)<<endl;

}

bool cmp_rank_pair(struct rank_pair& a, struct rank_pair& b) {
  return a.val < b.val;
}
