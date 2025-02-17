/*

   This file is part of the regenie software package.

   Copyright (c) 2020-2024 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

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
#include "Files.hpp"
#include "Geno.hpp"
#include "Joint_Tests.hpp"
#include "survival_data.hpp"
#include "cox_ridge.hpp"
#include "cox_score.hpp"
#include "cox_firth.hpp"
#include "Step1_Models.hpp"
#include "Step2_Models.hpp"
#include "HLM.hpp"
#include "Pheno.hpp"
#include "MultiTrait_Tests.hpp"
#include "Ordinal.hpp"
#include "SKAT.hpp"
#include "Masks.hpp"
#include "Data.hpp"

using namespace std;
using namespace Eigen;
using namespace boost;
using boost::math::beta_distribution;
using boost::math::chi_squared;


// null models
void fit_null_logistic(bool const& silent, const int& chrom, struct param* params, struct phenodt* pheno_data, struct ests* m_ests, struct in_files* files, mstream& sout, bool const& save_betas) {

  if(!silent) sout << "   -fitting null logistic regression on binary phenotypes..." << flush;
  if(params->test_mode) params->pheno_pass = true;

  auto t1 = std::chrono::high_resolution_clock::now();
  ArrayXd betaold, etavec, pivec, loco_offset, wvec;
  MatrixXd XtW;
  if(params->w_interaction || params->firth || (params->use_SPA && params->vc_test) || params->write_null_firth) m_ests->bhat_start.resize(pheno_data->new_cov.cols(), params->n_pheno);
  if(params->w_interaction) m_ests->offset_nullreg.resize(pheno_data->new_cov.rows(), params->n_pheno);
  betaold = ArrayXd::Zero(pheno_data->new_cov.cols());

  for(int i = 0; i < params->n_pheno; ++i ){

    if( !params->pheno_pass(i) ) continue;

    MapArXd Y (pheno_data->phenotypes_raw.col(i).data(), pheno_data->phenotypes_raw.rows());
    MapArXb mask (pheno_data->masked_indivs.col(i).data(), pheno_data->masked_indivs.rows());

    if(params->blup_cov) {
      pheno_data->new_cov.rightCols(1) = (m_ests->blups.col(i).array() * mask.cast<double>()).matrix();
      loco_offset = ArrayXd::Zero(Y.size(), 1);
    } else if(params->test_mode) 
      loco_offset = m_ests->blups.col(i).array() * mask.cast<double>();
    else loco_offset = ArrayXd::Zero(Y.size(), 1);

    // starting values
    betaold = 0;
    if(params->print_cov_betas) {betaold(0) = (0.5 + Y.sum()) / (pheno_data->Neff(i) + 1); betaold(0) = log( betaold(0) / (1 - betaold(0))) - loco_offset.mean();}
    get_pvec(etavec, pivec, betaold, loco_offset, pheno_data->new_cov, params->numtol_eps);

    // check if model converged
    if(!(fit_logistic(Y, pheno_data->new_cov, loco_offset, mask, pivec, etavec, betaold, params, sout, true, params->numtol) || fit_logistic(Y, pheno_data->new_cov, loco_offset, mask, pivec, etavec, betaold, params, sout, false, params->numtol))) {
      bool skip_pheno = true;

      // if not, get starting values by omitting loco offset (instead of at 0)
      if(params->test_mode && !params->skip_blups && !params->blup_cov){
        if(!silent) sout << "\n     WARNING: logistic regression did not converge for phenotype '" << files->pheno_names[i] <<"'. Retrying using starting beta from model without LOCO offset.";
        ArrayXd loco_dummy = ArrayXd::Zero(Y.size(), 1);
        pivec = ( 0.5 + Y ) / 2;
        etavec = mask.select( log(pivec/ (1-pivec)), 0);
        betaold = 0;
        if(params->print_cov_betas) {
          betaold(0) = (0.5 + Y.sum()) / (pheno_data->Neff(i) + 1); betaold(0) = log( betaold(0) / (1 - betaold(0)));
          get_pvec(etavec, pivec, betaold, loco_dummy, pheno_data->new_cov, params->numtol_eps);
        }
        if( fit_logistic(Y, pheno_data->new_cov, loco_dummy, mask, pivec, etavec, betaold, params, sout, true, params->numtol) || fit_logistic(Y, pheno_data->new_cov, loco_dummy, mask, pivec, etavec, betaold, params, sout, false, params->numtol) ){ 
          get_pvec(etavec, pivec, betaold, loco_dummy, pheno_data->new_cov, params->numtol_eps);
          skip_pheno = !(fit_logistic(Y, pheno_data->new_cov, loco_offset, mask, pivec, etavec, betaold, params, sout, true, params->numtol) || fit_logistic(Y, pheno_data->new_cov, loco_offset, mask, pivec, etavec, betaold, params, sout, false, params->numtol));
        }
      }

      if(skip_pheno){
        params->pheno_pass(i) = false; // phenotype will be ignored
        params->pheno_fail_nullreg(i) = true;
        if(!silent) sout << "\n     WARNING: logistic regression did not converge for phenotype '" << files->pheno_names[i] <<"'.";
        continue;
        // throw "logistic regression did not converge for phenotype " + files->pheno_names[i] + ". Perhaps increase --niter or check the covariates.";
      }
    } 
    
    if( !silent && (mask && (pivec < params->numtol_eps || pivec > 1 - params->numtol_eps)).any() )
      sout << "\n     WARNING: Fitted probabilities numerically 0/1 occurred (phenotype '" << files->pheno_names[i] <<"').";

    if(params->test_mode){
      if(save_betas && params->print_cov_betas) {
        params->cov_betas.col(i) = betaold;
        // get se
        get_wvec(pivec, wvec, mask, params->l1_ridge_eps);
        MatrixXd XWsqrt = ( pheno_data->new_cov.array().colwise() * (wvec.sqrt() * mask.cast<double>()) ).matrix();
        MatrixXd xtx_inv = ( XWsqrt.transpose() * XWsqrt ).colPivHouseholderQr().inverse();
        params->xtx_inv_diag.col(i).array() = xtx_inv.diagonal().array().sqrt();
        continue;
      }
      m_ests->Y_hat_p.col(i) = pivec.matrix() ;
      get_wvec(pivec, wvec, mask, params->l1_ridge_eps);
      m_ests->Gamma_sqrt.col(i) = wvec.sqrt().matrix();
      m_ests->Gamma_sqrt_mask.col(i) = (m_ests->Gamma_sqrt.col(i).array() * mask.cast<double>()).matrix();
      m_ests->X_Gamma[i] = m_ests->Gamma_sqrt_mask.col(i).asDiagonal() * pheno_data->new_cov;
      getBasis(m_ests->X_Gamma[i], params);
      if(params->w_interaction || params->firth || (params->use_SPA && params->vc_test)) m_ests->bhat_start.col(i) = betaold.matrix();
      if(params->w_interaction) m_ests->offset_nullreg.col(i) = etavec;
    } else {
      m_ests->offset_nullreg.col(i) = etavec;
      if(params->write_null_firth) m_ests->bhat_start.col(i) = betaold.matrix();
    }

    /*
     Files fstar;
     fstar.openForWrite("offsets.txt", sout);
     fstar << etavec;
     fstar.closeFile();
     */

  }

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  if(!silent) sout << "done (" << duration.count() << "ms) "<< endl;

}

bool fit_logistic(const Ref<const ArrayXd>& Y1, const Ref<const MatrixXd>& X1, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, ArrayXd& pivec, ArrayXd& etavec, ArrayXd& betavec, struct param const* params, mstream& sout, bool const& check_hs_dev, double const& numtol) {

  int niter_cur = 0, niter_search;
  double dev_old, dev_new=0, diff_dev;
  ArrayXd score, betanew, wvec, zvec;
  MatrixXd XtW, XtWX;

  dev_old = get_logist_dev(Y1, pivec, mask);
  //cerr << dev_old << endl << pivec.head(5)<<"\n--\n";

  while(niter_cur++ < params->niter_max){

    // p*(1-p) and check for zeroes
    if( get_wvec(pivec, wvec, mask, params->numtol_eps) ){
      if(params->verbose) sout << "ERROR: Zeros occurred in Var(Y) during logistic regression.\n";
      return false;
    }

    XtW = X1.transpose() * mask.select(wvec,0).matrix().asDiagonal();
    XtWX = XtW * X1;

    // working vector z = X*beta + (Y-p)/(p*(1-p))
    zvec = mask.select(etavec - offset + (Y1 - pivec) / wvec, 0);

    // parameter estimate
    betanew = ( XtWX ).colPivHouseholderQr().solve( XtW * zvec.matrix() ).array();

    // start step-halving
    for( niter_search = 1; niter_search <= params->niter_max_line_search; niter_search++ ){

      get_pvec(etavec, pivec, betanew, offset, X1, params->numtol_eps);
      dev_new = get_logist_dev(Y1, pivec, mask);

      if(params->debug) cerr << "HS#" << niter_search << setprecision(16) << ": p in (" << pivec.minCoeff() << "," << pivec.maxCoeff() << "); dev " << dev_old << "->" << dev_new << " \n";
      if( mask.select((pivec > 0) && (pivec < 1), true).all() && (!check_hs_dev || (dev_new < dev_old)) ) break;

      // adjust step size
      betanew = (betavec + betanew) / 2;

    }
    if( niter_search > params->niter_max_line_search ) return false; // step-halving failed

    // stopping criterion
    score = X1.transpose() * mask.select(Y1 - pivec, 0).matrix();
    if( score.abs().maxCoeff() < numtol ) break; // prefer for score to be below tol

    diff_dev = abs(dev_new - dev_old)/(0.1 + abs(dev_new));
    if(params->debug) cerr << "#" << niter_cur << ": score_max=" << score.abs().maxCoeff() << ";dev_diff=" << setprecision(16) << diff_dev << "\n";

    betavec = betanew;
    dev_old = dev_new;
  }
  if(params->debug) cerr << "Log. reg iter#" << niter_cur << ": beta=" << betanew.matrix().transpose() << "; score_max=" << score.abs().maxCoeff() << ";dev_diff=" << 
   setprecision(16) << diff_dev << "\n";

  // If didn't converge (check frac. change in deviances)
  if( ((diff_dev == 0) || (diff_dev >= numtol)) && (niter_cur > params->niter_max) )
    return false;

  betavec = betanew;

  return true;
}

// poisson models
void fit_null_poisson(const int& chrom, struct param* params, struct phenodt* pheno_data, struct ests* m_ests, struct in_files* files, mstream& sout, bool const& save_betas) {

  sout << "   -fitting null poisson regression..." << flush;

  auto t1 = std::chrono::high_resolution_clock::now();
  ArrayXd betaold, etavec, pivec, loco_offset, wvec;
  if(params->w_interaction) m_ests->bhat_start.resize(pheno_data->new_cov.cols(), params->n_pheno);

  for(int i = 0; i < params->n_pheno; ++i ){

    if( !params->pheno_pass(i) ) continue;

    MapArXd Y (pheno_data->phenotypes_raw.col(i).data(), pheno_data->phenotypes_raw.rows());
    MapArXb mask (pheno_data->masked_indivs.col(i).data(), pheno_data->masked_indivs.rows());

    if(params->test_mode) loco_offset = m_ests->blups.col(i).array() * mask.cast<double>();
    else loco_offset = ArrayXd::Zero(Y.size(), 1);

    // starting values
    pivec = Y + 1e-1; // avoid 0
    etavec = mask.select( log(pivec), 0);
    betaold = ArrayXd::Zero(pheno_data->new_cov.cols());
    betaold(0) = etavec.mean() - loco_offset.mean();

    if(!fit_poisson(Y, pheno_data->new_cov, loco_offset, mask, pivec, etavec, betaold, params, sout)){
      params->pheno_pass(i) = false; // phenotype will be ignored
      sout << "\n     WARNING: poisson regression did not converge for phenotype '" << files->pheno_names[i] <<"'.";
      continue;
      // throw "poisson regression did not converge for phenotype " + files->pheno_names[i] + ". Perhaps increase --niter?";
    }
    else if( (mask && pivec < params->numtol_eps).any() )
      sout << "\n     WARNING: Fitted rates numerically 0 occurred (phenotype #" << files->pheno_names[i] <<").";

    if(params->test_mode){
      if(save_betas && params->print_cov_betas) {
        params->cov_betas.col(i) = betaold;
        // get se
        MatrixXd XWsqrt = ( pheno_data->new_cov.array().colwise() * mask.select(pivec,0).sqrt() ).matrix();
        MatrixXd xtx_inv = ( XWsqrt.transpose() * XWsqrt ).colPivHouseholderQr().inverse();
        params->xtx_inv_diag.col(i).array() = xtx_inv.diagonal().array().sqrt();
        continue;
      }
      m_ests->Y_hat_p.col(i) = pivec.matrix() ;
      m_ests->Gamma_sqrt.col(i) = pivec.sqrt().matrix();
      m_ests->X_Gamma[i] = ( pheno_data->new_cov.array().colwise() * (m_ests->Gamma_sqrt.col(i).array() * mask.cast<double>()) ).matrix();
      getBasis(m_ests->X_Gamma[i], params);
      if(params->w_interaction) m_ests->bhat_start.col(i) = betaold.matrix();
    } else m_ests->offset_nullreg.col(i) = etavec;

    /*
     Files fstar;
     fstar.openForWrite("offsets.txt", sout);
     fstar << etavec;
     fstar.closeFile();
     */

  }

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;

}

bool fit_poisson(const Ref<const ArrayXd>& Y1, const Ref<const MatrixXd>& X1, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, ArrayXd& pivec, ArrayXd& etavec, ArrayXd& betavec, struct param const* params, mstream& sout) {

  bool dev_conv = false;
  int niter_cur = 0;
  double dev_old, dev_new=0;
  ArrayXd score, betanew, zvec;
  MatrixXd XtW, XtWX;

  dev_old = get_poisson_dev(Y1, pivec, mask);
  if(params->debug) cerr << "\nstarting deviance: " << dev_old << endl;

  while(niter_cur++ < params->niter_max){

    // check for zeroes
    if( (mask && (pivec == 0)).any() ){
      if(params->verbose) sout << "ERROR: Zeros occurred in Var(Y) during poisson regression.\n";
      return false;
    }

    XtW = X1.transpose() * mask.select(pivec,0).matrix().asDiagonal();
    XtWX = XtW * X1;

    // working vector z = X*beta + (Y-p)/w
    zvec = mask.select(etavec - offset + (Y1 - pivec) / pivec, 0);

    // parameter estimate
    betanew = ( XtWX ).colPivHouseholderQr().solve( XtW * zvec.matrix() ).array();

    // start step-halving
    for( int niter_search = 1; niter_search <= params->niter_max_line_search; niter_search++ ){

      get_pvec_poisson(etavec, pivec, betanew, offset, X1, params->numtol_eps);
      dev_new = get_poisson_dev(Y1, pivec, mask);
      if( !(mask && (pivec == 0)).any() ) break;

      // adjust step size
      betanew = (betavec + betanew) / 2;

    }
    if(params->debug && (niter_cur%5==0)) cerr << "#" << niter_cur << " -> " << 
      betanew.matrix().transpose() << " : " << dev_new << "\n";

    score = X1.transpose() * mask.select(Y1 - pivec, 0).matrix();

    // stopping criterion
    dev_conv = (abs(dev_new - dev_old)/(0.1 + abs(dev_new))) < params->tol;
    if( score.abs().maxCoeff() < params->tol ) break;

    betavec = betanew;
    dev_old = dev_new;
  }

  // If didn't converge
  if( !dev_conv && (niter_cur > params->niter_max) )
    return false;

  betavec = betanew;

  if(params->debug) cerr << "Final (" << niter_cur << ") : " << betavec.matrix().transpose() << " : " << score.abs().maxCoeff() << "\n";
  return true;
}

// Cox regression null model
void fit_null_cox(bool const& silent, const int& chrom, struct param* params, struct phenodt* pheno_data, struct ests* m_ests, struct in_files* files, mstream& sout, bool const& save_betas){
  if(!silent) sout << "   -fitting null cox regression on time-to-event phenotypes..." << flush;
  
  Eigen::VectorXd loco_offset;
  auto t1 = std::chrono::high_resolution_clock::now();

  for (const auto& entry: files->t2e_map) {
    const std::string& time_name = entry.first;
    const std::string& event_name = entry.second;
    // find time column index
    std::vector<std::string>::iterator it_time = std::find(files->pheno_names.begin(), files->pheno_names.end(), time_name);
    int time_index = std::distance(files->pheno_names.begin(), it_time);
    Eigen::VectorXd ph_time = pheno_data->phenotypes_raw.col(time_index);
    MapArXb mask (pheno_data->masked_indivs.col(time_index).data(), pheno_data->masked_indivs.rows());

    // find event column index
    std::vector<std::string>::iterator it_event = std::find(files->pheno_names.begin(), files->pheno_names.end(), event_name);
    int event_index = std::distance(files->pheno_names.begin(), it_event);
    Eigen::VectorXd ph_event = pheno_data->phenotypes_raw.col(event_index);
    
    if(params->blup_cov) {
      pheno_data->new_cov.rightCols(1) = (m_ests->blups.col(time_index).array() * mask.cast<double>()).matrix();
      loco_offset = Eigen::VectorXd::Zero(ph_time.size());
    } else if(params->test_mode) 
      loco_offset = (m_ests->blups.col(time_index).array() * mask.cast<double>()).matrix();
    else loco_offset = Eigen::VectorXd::Zero(ph_time.size());

    survival_data survivalNullData;
    survivalNullData.setup(ph_time, ph_event, mask, !params->test_mode);

    cox_ridge coxRidge_null_lamb0(survivalNullData, pheno_data->new_cov, loco_offset, mask, 0, params->niter_max, params->niter_max_line_search, params->numtol_cox);
    coxRidge_null_lamb0.fit(survivalNullData, pheno_data->new_cov, loco_offset, mask);

    if (params->test_mode) {
      cox_mle coxMLE;
      coxMLE.setup(survivalNullData, pheno_data->new_cov, loco_offset, mask, params->niter_max, params->niter_max_line_search, params->numtol_cox, false, coxRidge_null_lamb0.beta, coxRidge_null_lamb0.eta);
      coxMLE.fit(survivalNullData, pheno_data->new_cov, loco_offset, mask);

      if (coxMLE.converge == false) {
        cox_firth cox_null_model;
        cox_null_model.setup(survivalNullData, pheno_data->new_cov, loco_offset, pheno_data->new_cov.cols(), params->niter_max, params->niter_max_line_search, params->numtol_cox, 0, params->numtol_beta_cox, params->maxstep_null, false, false);
        cox_null_model.fit(survivalNullData, pheno_data->new_cov, loco_offset);

        coxMLE.setup(survivalNullData, pheno_data->new_cov, loco_offset, mask, params->niter_max, params->niter_max_line_search, params->numtol_cox, false, cox_null_model.beta, cox_null_model.eta);
        coxMLE.fit(survivalNullData, pheno_data->new_cov, loco_offset, mask);
      }

      if (coxMLE.converge == false) {
        params->pheno_pass(time_index) = false; // phenotype will be ignored
        params->pheno_fail_nullreg(time_index) = true;
        if(!silent) sout << "\n     WARNING: step2 cox null regression did not converge for phenotype '" << time_name <<"'.";
        continue;
      } else {
        coxMLE.cox_test_prep(survivalNullData, pheno_data->new_cov, loco_offset, mask);
        m_ests->survival_data_pheno[time_index] = survivalNullData;
        m_ests->cox_MLE_NULL[time_index] = coxMLE;
      }

      if(save_betas && params->print_cov_betas) {
        params->cov_betas.col(time_index) = coxMLE.beta;
        params->xtx_inv_diag.col(time_index).array() = coxMLE.XtWX.diagonal().array().sqrt();
      }
    } else {
      if (coxRidge_null_lamb0.converge == false) {
        // try cox firth, without firth
        cox_firth cox_null_model;
        cox_null_model.setup(survivalNullData, pheno_data->new_cov, loco_offset, pheno_data->new_cov.cols(), params->niter_max, params->niter_max_line_search, params->numtol_cox, params->numtol_cox_stephalf, params->numtol_beta_cox, params->maxstep_null, false, false);
        cox_null_model.fit(survivalNullData, pheno_data->new_cov, loco_offset);

        if (cox_null_model.converge == false) {
          cox_null_model.setup(survivalNullData, pheno_data->new_cov, loco_offset, pheno_data->new_cov.cols(), params->niter_max, params->niter_max_line_search, params->numtol_cox, 0, params->numtol_beta_cox, params->maxstep_null, false, false);
          cox_null_model.fit(survivalNullData, pheno_data->new_cov, loco_offset);
        }

        if (cox_null_model.converge == false) {
          params->pheno_pass(time_index) = false; // phenotype will be ignored
          params->pheno_fail_nullreg(time_index) = true;
          if(!silent) sout << "\n     WARNING: step1 cox null regression did not converge for phenotype '" << time_name <<"'.";
          continue;
        } else {
          m_ests->offset_nullreg.col(time_index) = cox_null_model.eta;
        }
      } else {
        m_ests->offset_nullreg.col(time_index) = coxRidge_null_lamb0.eta;
      }
      
    }
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  if(!silent) sout << "done (" << duration.count() << "ms) "<< endl;
}

double getCoxLambdaMax(const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& gradient) {
    Eigen::VectorXd g = (Xmat.transpose() * gradient).array().abs();
    double lambda_max = g.maxCoeff() / 1e-3;
    return lambda_max;
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          level 0 models
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void ridge_level_0(const int& block, struct in_files* files, struct param* params, struct filter* filters, struct ests* m_ests, struct geno_block* Gblock, struct phenodt* pheno_data, vector<snp>& snpinfo, struct ridgel0* l0, struct ridgel1* l1, vector<MatrixXb>& masked_in_folds, mstream& sout) {

  sout << "   -calc level 0 ridge..." << flush;
  auto t2 = std::chrono::high_resolution_clock::now();

  int bs = l0->GGt.rows();
  int block_eff = params->write_l0_pred ? 0 : block;
  string op_name, out_pheno;
  ofstream ofile;

  MatrixXd ww1, ww2, beta, pred, vmat, dvec, Xout;
  MatrixXd ident_l0 = MatrixXd::Identity(bs, bs);
  MatrixXd p_sum = MatrixXd::Zero(params->n_ridge_l0, params->n_pheno);
  MatrixXd p_sum2 = MatrixXd::Zero(params->n_ridge_l0, params->n_pheno);

  if(!params->within_sample_l0 && params->print_block_betas) {
    for(int ph = 0; ph < params->n_pheno; ++ph )
      if( params->pheno_pass(ph) )
        params->beta_print_out[ph] = MatrixXd::Zero(params->n_ridge_l0, bs);
  }

  uint32_t cum_size_folds = 0;
  for(int i = 0; i < params->cv_folds; ++i ) {
    // assign masking within folds
    masked_in_folds[i] = pheno_data->masked_indivs.block(cum_size_folds, 0, params->cv_sizes(i), pheno_data->masked_indivs.cols());

    ww1 = l0->GGt - l0->G_folds[i];
    SelfAdjointEigenSolver<MatrixXd> eig(ww1);
    vmat = eig.eigenvectors();
    dvec = eig.eigenvalues();
    //if(i == 0)sout << ww1 << endl;
    ww2 = vmat.transpose() * (l0->GTY - l0->GtY[i]);

    for(int j = 0; j < params->n_ridge_l0; ++j ) {

      // b = U (D+sI)^(-1) U^t GtY
      beta = vmat * (dvec.array() + params->lambda(j)).inverse().matrix().asDiagonal() * ww2;

      // save beta for each phenotype (only when using out-of-sample pred)
      if(!params->within_sample_l0 && params->print_block_betas)
        for(int ph = 0; ph < params->n_pheno; ++ph ) 
          if( params->pheno_pass(ph) )
            params->beta_print_out[ph].row(j) += beta.col(ph).transpose();

      // out-of-sample predictions (mask missing)
      pred = ( (beta.transpose() * Gblock->Gmat.block(0, cum_size_folds, bs, params->cv_sizes(i))).array()  * masked_in_folds[i].transpose().array().cast<double>() ).matrix();
      p_sum.row(j) += pred.rowwise().sum();
      p_sum2.row(j) += pred.rowwise().squaredNorm();

      // store predictions
      for(int ph = 0; ph < params->n_pheno; ++ph ) {
        if( !params->pheno_pass(ph) ) continue;
        if (params->trait_mode != 3) {
          l1->test_mat[ph][i].col(block_eff * params->n_ridge_l0 + j) = pred.row(ph).transpose();
        } else {
          l1->test_mat_conc[ph].block(cum_size_folds, block_eff * params->n_ridge_l0 + j, params->cv_sizes(i), 1) = pred.row(ph).transpose();
        }
        if((block == 0) && (j == 0)) { // same for all blocks & ridge params
          if (params->trait_mode != 3) {
            l1->test_pheno[ph][i].col(0) = pheno_data->phenotypes.block(cum_size_folds, ph, params->cv_sizes(i), 1);
            if (params->trait_mode != 0) {
              l1->test_pheno_raw[ph][i].col(0) = pheno_data->phenotypes_raw.block(cum_size_folds, ph, params->cv_sizes(i), 1);
              l1->test_offset[ph][i].col(0) = m_ests->offset_nullreg.block(cum_size_folds, ph, params->cv_sizes(i), 1);
            }
          } else {
            l1->fold_id[ph](seqN(cum_size_folds, params->cv_sizes(i))) = Eigen::VectorXi::Constant(params->cv_sizes(i), i);
          }
        }
      }
    }
    cum_size_folds += params->cv_sizes(i);
  }
  if(params->debug && (block < 5) && params->trait_mode != 3) {
    if(params->test_l0)
      cerr << "Ymat (Y1):\n" << (l1->test_pheno[0][0].topRows(5) - l1->top_snp_pgs[0].topRows(5)) << endl;
    else
      cerr << "Ymat (Y1):\n" << l1->test_pheno[0][0].topRows(5) << endl;
    cerr << "Wmat (Y1):\n" << l1->test_mat[0][0].topRows(5).middleCols(block_eff * params->n_ridge_l0, params->n_ridge_l0) << endl;
  }

  // center and scale using the whole sample
  for(int ph = 0; ph < params->n_pheno; ++ph ) {

    if( !params->pheno_pass(ph) ) continue;

    RowVectorXd p_mean, p_invsd;
    p_mean = p_sum.col(ph).transpose() / pheno_data->Neff(ph);
    p_invsd = sqrt( (pheno_data->Neff(ph) - 1) / (p_sum2.col(ph).transpose().array() - pheno_data->Neff(ph) * p_mean.array().square()) );

    // scale printed estimates by the sd
    if(params->print_block_betas)
      params->beta_print_out[ph].array().colwise() *= p_invsd.transpose().array();

    if(params->write_l0_pred) Xout = MatrixXd::Zero(params->n_samples, params->n_ridge_l0);

    cum_size_folds = 0;
    for(int i = 0; i < params->cv_folds; ++i ) {
      if( params->trait_mode != 3){
        l1->test_mat[ph][i].block(0, block_eff * params->n_ridge_l0, params->cv_sizes(i), params->n_ridge_l0).rowwise() -= p_mean;
        l1->test_mat[ph][i].block(0, block_eff * params->n_ridge_l0, params->cv_sizes(i), params->n_ridge_l0).array().rowwise() *= p_invsd.array();
      } else {
        l1->test_mat_conc[ph].block(cum_size_folds, block_eff * params->n_ridge_l0, params->cv_sizes(i), params->n_ridge_l0).rowwise() -= p_mean;
        l1->test_mat_conc[ph].block(cum_size_folds, block_eff * params->n_ridge_l0, params->cv_sizes(i), params->n_ridge_l0).array().rowwise() *= p_invsd.array();
      }

      if(params->write_l0_pred) {
        if (params->trait_mode != 3) {
          Xout.block(cum_size_folds, 0, params->cv_sizes(i), params->n_ridge_l0) = l1->test_mat[ph][i].block(0, block_eff * params->n_ridge_l0, params->cv_sizes(i), params->n_ridge_l0);
        } else {
          Xout.block(cum_size_folds, 0, params->cv_sizes(i), params->n_ridge_l0) = l1->test_mat_conc[ph].block(cum_size_folds, block_eff * params->n_ridge_l0, params->cv_sizes(i), params->n_ridge_l0);
        }
      }
      cum_size_folds += params->cv_sizes(i);
    }

    // write predictions to file if specified
    if(params->write_l0_pred) {
      write_l0_file(files->write_preds_files[ph].get(), Xout, sout);
      //if(block ==0 && ph == 0 ) sout << endl << "Out " << endl <<  Xout.block(0, 0, 3, 3) << endl;
    }

  }

  // if printing betas to file (average over folds) [assume snp IDs are unique]
  //   -> separate file for each block (params->n_ridge_l0 rows & (2+bs) columns)
  if(!params->within_sample_l0 && params->print_block_betas) {
    op_name = files->out_file + "_block" + to_string(block+1) + ".betas";
    openStream(&ofile, op_name, std::ios::out, sout);

    // Header: [TRAIT PARAM snpID1 ... snpIDk]
    ofile << "TRAIT PARAM " ;
    for(int i = 0; i < bs; ++i )
      ofile << snpinfo[params->print_snpcount++].ID << " ";
    ofile << endl;

    // Each line: [pheno# ridge# beta1 ... betak]
    for(int ph = 0; ph < params->n_pheno; ++ph ){

      if( !params->pheno_pass(ph) ) continue;

      params->beta_print_out[ph] /= params->cv_folds;
      for(int j = 0; j < params->n_ridge_l0; ++j ) {
        ofile << ph + 1 << " " <<  j + 1 << " ";
        for(int i = 0; i < bs; ++i )
          ofile << params->beta_print_out[ph](j,i) << " ";
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

void ridge_level_0_loocv(const int block, struct in_files* files, struct param* params, struct filter* filters, struct ests* m_ests, struct geno_block* Gblock, struct phenodt* pheno_data, vector<snp>& snpinfo, struct ridgel0* l0, struct ridgel1* l1, mstream& sout) {

  sout << "   -calc level 0 ridge..." << flush;
  auto t2 = std::chrono::high_resolution_clock::now();
  int bs = l0->GGt_eig_val.size();
  int block_eff = params->write_l0_pred ? 0 : block; // if writing to file
  string out_pheno;
  ofstream ofile;
  VectorXd z1, gvec;
  MatrixXd VtG, z2, pred, Xout;
  RowVectorXd p_mean, p_sd;

  /*
     if(bs > params->n_samples)
     throw "block size must be smaller than the number of samples to perform LOOCV!";
     */

  uint64 max_bytes = params->chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (bs * target_size) matrix ]
  int nchunk = ceil( params->cv_folds * bs * sizeof(double) * 1.0 / max_bytes );
  if (params->verbose) sout << nchunk << " chunks..." << flush;
  int chunk, size_chunk, target_size = params->cv_folds / nchunk;
  int j_start;

  // make matrix of (eigen-value + lambda)^(-1)
  MatrixXd DL_inv = ( l0->GGt_eig_val.rowwise().replicate(params->n_ridge_l0).array().rowwise() + params->lambda.matrix().transpose().array() ).inverse().matrix(); // kxR

  if(!params->test_l0 || (l0->nspns_picked_block.maxCoeff() == 0) || (params->n_pheno == 1)){

    if(params->print_block_betas) // assumes P=1
      l1->beta_snp_step1.middleRows(filters->step1_snp_count, bs) = l0->GGt_eig_vec * (DL_inv.array().colwise() * l0->Wmat.col(0).array()).matrix();// K x R

    for(chunk = 0; chunk < nchunk; ++chunk ) {
      size_chunk = chunk == nchunk - 1? params->cv_folds - target_size * chunk : target_size;
      j_start = chunk * target_size;

      if((params->n_pheno == 1) && l0->subset_l0_snps_gmat)
        VtG = l0->GGt_eig_vec.transpose() * Gblock->Gmat(l0->indices_gmat_keep, seqN(j_start, size_chunk));
      else
        VtG = l0->GGt_eig_vec.transpose() * Gblock->Gmat(all, seqN(j_start, size_chunk));
      for(int i = 0; i < size_chunk; ++i ) {
        z1 = VtG.col(i); // Kx1
        z2 = DL_inv.array().colwise() * z1.array(); // K x R
        gvec = z2.transpose() * z1; // R x 1
        if(params->test_l0)
          pred = z2.transpose() * l0->Wmat - gvec * l0->ymat_res.row(j_start + i);
        else
          pred = z2.transpose() * l0->Wmat - gvec * pheno_data->phenotypes.row(j_start + i);
        pred.array().colwise() /= 1 - gvec.array(); // R x P

        for(int ph = 0; ph < params->n_pheno; ++ph )
          if( params->pheno_pass(ph) )
            l1->test_mat_conc[ph].block(j_start + i, block_eff * params->n_ridge_l0, 1, params->n_ridge_l0) = pred.col(ph).transpose();
      }
    } 

  } else for(chunk = 0; chunk < nchunk; ++chunk ) {
    size_chunk = chunk == nchunk - 1? params->cv_folds - target_size * chunk : target_size;
    j_start = chunk * target_size;
    VtG = l0->GGt_eig_vec.transpose() * Gblock->Gmat(l0->indices_gmat_keep, seqN(j_start, size_chunk)); // k x N
    MatrixXd gamma_rho = (DL_inv.transpose() * VtG.array().square().matrix()).transpose(); // N x R
    for(int ph = 0; ph < params->n_pheno; ++ph ) {
      if(!params->pheno_pass(ph) ) continue;
      Ref<MatrixXd> X_l1 = l1->test_mat_conc[ph].block(j_start, block_eff * params->n_ridge_l0, size_chunk, params->n_ridge_l0); // NxR
      X_l1 = VtG.transpose() * (DL_inv.array().colwise() * l0->Wmat.col(ph).array()).matrix() - (gamma_rho.array().colwise() * l0->ymat_res.col(ph).segment(j_start, size_chunk).array()).matrix(); // N x R
      X_l1.array() /= (1 - gamma_rho.array());
    }
  }

  if(params->debug && (block < 5)) {
    if(params->test_l0)
      cerr << "Ymat (Y1-Y5):\n" << l0->ymat_res.topLeftCorner(5,min(5, params->n_pheno)) << endl;
    else
      cerr << "Ymat (Y1-Y5):\n" << pheno_data->phenotypes.topLeftCorner(5,min(5, params->n_pheno)) << endl;
    cerr << "Wmat (Y1):\n" << l1->test_mat_conc[0].topRows(5).middleCols(block_eff * params->n_ridge_l0, params->n_ridge_l0) << endl;
  }

  // center and scale within the block
  for(int ph = 0; ph < params->n_pheno; ++ph ) {
    if( !params->pheno_pass(ph) ) continue;
    // mask missing first
    l1->test_mat_conc[ph].middleCols(block_eff * params->n_ridge_l0, params->n_ridge_l0).array().colwise() *= pheno_data->masked_indivs.col(ph).array().cast<double>();
    p_mean = l1->test_mat_conc[ph].middleCols(block_eff * params->n_ridge_l0, params->n_ridge_l0).colwise().sum() / pheno_data->Neff(ph);
    //if(i == 0)sout << i << " " << p_mean << endl;
    l1->test_mat_conc[ph].middleCols(block_eff * params->n_ridge_l0, params->n_ridge_l0).rowwise() -= p_mean;
    // mask missing again
    l1->test_mat_conc[ph].middleCols(block_eff * params->n_ridge_l0, params->n_ridge_l0).array().colwise() *= pheno_data->masked_indivs.col(ph).array().cast<double>();
    p_sd = l1->test_mat_conc[ph].middleCols(block_eff * params->n_ridge_l0, params->n_ridge_l0).colwise().norm() / sqrt(pheno_data->Neff(ph) -1);
    //if(i == 0)sout << i << " " << p_sd << endl;
    l1->test_mat_conc[ph].middleCols(block_eff * params->n_ridge_l0, params->n_ridge_l0).array().rowwise() /= p_sd.array();


    if(params->write_l0_pred) {
      Xout = l1->test_mat_conc[ph].leftCols(params->n_ridge_l0);
      write_l0_file(files->write_preds_files[ph].get(), Xout, sout);
      //if(block < 2 && ph == 0 ) sout << endl << "Out " << endl <<  Xout.topLeftCorner(5, Xout.cols()) << endl;
    }

    if(params->print_block_betas) {
      l1->beta_snp_step1.middleRows(filters->step1_snp_count, bs) *= (1/p_sd.array()).matrix().asDiagonal();
      //cerr << "Gb:\n"<<(Gblock->Gmat.transpose() * l1->beta_snp_step1.middleRows(filters->step1_snp_count, bs)).topRows(5) << "\n\n" << "W(centered):\n"<<
       // l1->test_mat_conc[0].block(0, block_eff * params->n_ridge_l0, 5, params->n_ridge_l0) << "\n\n";
    }

  }

  sout << "done";
  auto t3 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
  sout << " (" << duration.count() << "ms) "<< endl;

}

void write_l0_file(ofstream* ofs, MatrixXd& Xout, mstream& sout){

  ofs->write( reinterpret_cast<char *> (&Xout(0,0)), Xout.size() * sizeof(double) );
  if( ofs->fail() )
    throw "cannot successfully write temporary level 0 predictions to disk";

}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          level 1 models
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void set_mem_l1(struct in_files* files, struct param* params, struct filter* filters, struct ests* m_ests, struct geno_block* Gblock, struct phenodt* pheno_data, struct ridgel1* l1, vector<MatrixXb>& masked_in_folds, mstream& sout){ // when l0 was run in parallel

  if(!(params->use_loocv || params->trait_mode == 3)) {
    uint32_t cum_size_folds = 0;
    for(int i = 0; i < params->cv_folds; ++i ) {
      // assign masking within folds
      masked_in_folds[i] = pheno_data->masked_indivs.middleRows(cum_size_folds, params->cv_sizes(i));
      // store predictions
      for(int ph = 0; ph < params->n_pheno; ++ph ) {
        if( !params->pheno_pass(ph) ) continue;
        l1->test_pheno[ph][i] = pheno_data->phenotypes.block(cum_size_folds, ph, params->cv_sizes(i), 1);
        if (params->trait_mode) {
          l1->test_pheno_raw[ph][i] = pheno_data->phenotypes_raw.block(cum_size_folds, ph, params->cv_sizes(i), 1);
          l1->test_offset[ph][i] = m_ests->offset_nullreg.block(cum_size_folds, ph, params->cv_sizes(i), 1);
        }
      }
      cum_size_folds += params->cv_sizes(i);
    }
  } else if (params->trait_mode == 3) {
    uint32_t cum_size_folds = 0;
    for(int i = 0; i < params->cv_folds; ++i ) {
      for(int ph = 0; ph < params->n_pheno; ++ph ) {
        if( !params->pheno_pass(ph) ) continue;
        l1->fold_id[ph](seqN(cum_size_folds, params->cv_sizes(i))) = Eigen::VectorXi::Constant(params->cv_sizes(i), i);
      }
      cum_size_folds += params->cv_sizes(i);
    }
  }
}

void ridge_level_1(struct in_files* files, struct param* params, struct phenodt* pheno_data, struct ridgel1* l1, mstream& sout) {

  sout << endl << " Level 1 ridge..." << endl << flush;

  string in_pheno;
  ifstream infile;
  MatrixXd X1, X2, beta_l1, p1, vmat, dvec, dl_inv, XtX_tau;
  VectorXd VtX2;
  MatrixXd XtX_sum, XtY_sum;

  // to compute Rsq and MSE of predictions
  for (int i = 0; i < 5; i++){
    l1->cumsum_values[i].setZero(params->n_pheno, params->n_ridge_l1);
    if(params->test_l0) l1->cumsum_values_full[i].setZero(params->n_pheno, params->n_ridge_l1);
  }

  for(int ph = 0; ph < params->n_pheno; ++ph ) {
    if( !params->pheno_pass(ph) ) continue; // should not happen for qts

    sout << "   -on phenotype " << ph+1 <<" (" << files->pheno_names[ph] << ")..." << flush;
    auto ts1 = std::chrono::high_resolution_clock::now();
    int ph_eff = params->write_l0_pred ? 0 : ph;
    int bs_l1 = params->total_n_block * params->n_ridge_l0;

    // read in level 0 predictions from file
    if(params->write_l0_pred)
      read_l0(ph, ph_eff, files, params, l1, sout);
    check_l0(ph, ph_eff, params, l1, pheno_data, sout);
    bs_l1 = l1->test_mat[ph_eff][0].cols();
    bool use_simple_ridge = (l1->ridge_param_mult == 1).all();

    // compute XtX and Xty for each fold and cum. sum using test_mat's
    if (!params->within_sample_l0){
      XtX_sum.setZero(bs_l1, bs_l1);
      XtY_sum.setZero(bs_l1, 1);
      for( int i = 0; i < params->cv_folds; ++i ) {
        l1->X_folds[i] = l1->test_mat[ph_eff][i].transpose() * l1->test_mat[ph_eff][i];
        l1->XtY[i]     = l1->test_mat[ph_eff][i].transpose() * l1->test_pheno[ph][i];
        XtX_sum += l1->X_folds[i];
        XtY_sum += l1->XtY[i];
      }
    }

    uint32_t cum_size_folds = 0;
    for(int i = 0; i < params->cv_folds; ++i ) {

      // use either in-sample or out-of-sample predictions
      if (params->within_sample_l0) { // DEPRECATED
        X1 = l1->pred_mat[ph][i].transpose() * l1->pred_mat[ph][i];
        X2 = l1->pred_mat[ph][i].transpose() * l1->pred_pheno[ph][i];
      } else{
        X1 = XtX_sum - l1->X_folds[i];
        X2 = XtY_sum - l1->XtY[i];
      }

      if(use_simple_ridge){
        SelfAdjointEigenSolver<MatrixXd> eigX1(X1);
        vmat = eigX1.eigenvectors();
        dvec = eigX1.eigenvalues();
        VtX2 = vmat.transpose() * X2;
      // compute solutions for all ridge parameters at once
        dl_inv = ( dvec.rowwise().replicate(params->n_ridge_l1) + l1->ridge_param_mult.matrix() * params->tau[ph].matrix().transpose() ).array().inverse().matrix();
        dl_inv.array().colwise() *= VtX2.array();
        beta_l1 = vmat * dl_inv;
      } else { // need to compute seperately for each parameter
        beta_l1.resize(bs_l1, params->n_ridge_l1);
        for(int j = 0; j < params->n_ridge_l1; ++j) {
          XtX_tau = X1;
          XtX_tau.diagonal().array() += params->tau[ph](j) * l1->ridge_param_mult;
          SelfAdjointEigenSolver<MatrixXd> eigMat(XtX_tau);
          beta_l1.col(j) = eigMat.eigenvectors() * (1/eigMat.eigenvalues().array()).matrix().asDiagonal() * eigMat.eigenvectors().transpose() * X2;
        }
      }
      if(!params->within_sample_l0) l1->beta_hat_level_1[ph][i] = beta_l1;
      // p1 is Nfold x nridge_l1 matrix
      p1 = l1->test_mat[ph_eff][i] * beta_l1;
      l1->cumsum_values[0].row(ph) += p1.colwise().sum();
      l1->cumsum_values[1].row(ph).array() += l1->test_pheno[ph][i].array().sum();
      l1->cumsum_values[2].row(ph) += p1.array().square().matrix().colwise().sum();
      l1->cumsum_values[3].row(ph).array() += l1->test_pheno[ph][i].array().square().sum();
      l1->cumsum_values[4].row(ph) += (p1.array().colwise() * l1->test_pheno[ph][i].col(0).array()).matrix().colwise().sum() ;
      if(params->test_l0){ // pred = p1 + top_snp_pgs; Y is res pheno
        p1.colwise() += l1->top_snp_pgs[0].col(ph).segment(cum_size_folds, params->cv_sizes(i));
        l1->cumsum_values_full[0].row(ph) += p1.colwise().sum();
        l1->cumsum_values_full[1].row(ph).array() += pheno_data->phenotypes.block(cum_size_folds, ph, params->cv_sizes(i), 1).array().sum();
        l1->cumsum_values_full[2].row(ph) += p1.array().square().matrix().colwise().sum();
        l1->cumsum_values_full[3].row(ph).array() += pheno_data->phenotypes.block(cum_size_folds, ph, params->cv_sizes(i), 1).array().square().sum();
        l1->cumsum_values_full[4].row(ph) += (p1.array().colwise() * pheno_data->phenotypes.block(cum_size_folds, ph, params->cv_sizes(i), 1).col(0).array()).matrix().colwise().sum() ;
      }

      cum_size_folds += params->cv_sizes(i);
    }

    sout << "done";
    auto ts2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(ts2 - ts1);
    sout << " (" << duration.count() << "ms) "<< endl;
  }

  sout << endl;
}


void ridge_level_1_loocv(struct in_files* files, struct param* params, struct phenodt* pheno_data, struct ridgel1* l1, mstream& sout) {

  sout << endl << " Level 1 ridge..." << flush;

  int bs_l1 = params->total_n_block * params->n_ridge_l0;
  int ph_eff;
  string in_pheno;
  ifstream infile;
  MatrixXd XH_chunk, Z1, Z2, dl, dl_inv, xtx, tmpMat;
  VectorXd wvec, zvec, Yvec, tmpVec, calFactor, pred;

  for (int i = 0; i < 5; i++){
    l1->cumsum_values[i].setZero(params->n_pheno, params->n_ridge_l1);
    if(params->test_l0) l1->cumsum_values_full[i].setZero(params->n_pheno, params->n_ridge_l1);
  }
  if(params->test_l0) l1->cumsum_values_full[3].array().colwise() += pheno_data->Neff - params->ncov; 
  else l1->cumsum_values[3].array().colwise() += pheno_data->Neff - params->ncov; // Sy2

  uint64 max_bytes = params->chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (target_size * bs_l1) matrix ]
  int nchunk = ceil( params->cv_folds * bs_l1 * sizeof(double) * 1.0 / max_bytes );
  if (params->verbose) sout << nchunk << " chunks...";
  sout << endl;
  int chunk, size_chunk, target_size = params->cv_folds / nchunk;
  int j_start;

  for(int ph = 0; ph < params->n_pheno; ++ph ) {
    if( !params->pheno_pass(ph) ) continue; // should not happen for qts
    sout << "   -on phenotype " << ph+1 <<" (" << files->pheno_names[ph] <<")..." << flush;
    auto ts1 = std::chrono::high_resolution_clock::now();
    ph_eff = params->write_l0_pred ? 0 : ph;

    // read in level 0 predictions from file
    if(params->write_l0_pred)
      read_l0(ph, ph_eff, files, params, l1, sout);
    check_l0(ph, ph_eff, params, l1, pheno_data, sout);
    bs_l1 = l1->test_mat_conc[ph_eff].cols();
    bool use_simple_ridge = (l1->ridge_param_mult == 1).all();
    if(params->test_l0)
      Yvec = pheno_data->phenotypes.col(ph) - l1->top_snp_pgs[0].col(ph);
    else
      Yvec = pheno_data->phenotypes.col(ph);

    xtx = l1->test_mat_conc[ph_eff].transpose() * l1->test_mat_conc[ph_eff]; // kxk
    SelfAdjointEigenSolver<MatrixXd> eigX(xtx);
    zvec = eigX.eigenvectors().transpose() * (l1->test_mat_conc[ph_eff].transpose() * Yvec); // kx1

    for(chunk = 0; chunk < nchunk; ++chunk) {
      size_chunk = chunk == nchunk - 1? params->cv_folds - target_size * chunk : target_size;
      j_start = chunk * target_size;
      Ref<VectorXd> Y_chunk = Yvec.segment(j_start, size_chunk);
      Ref<MatrixXd> X_chunk = l1->test_mat_conc[ph_eff].middleRows(j_start, size_chunk);

      tmpMat = X_chunk * eigX.eigenvectors(); // N_c x k

      for(int j = 0; j < params->n_ridge_l1; ++j) { // compute seperately for each parameter
        if(use_simple_ridge)
          tmpVec = (1/(eigX.eigenvalues().array() + params->tau[ph](j))).matrix(); // kx1
        else
          tmpVec = (1/(eigX.eigenvalues().array() + params->tau[ph](j) * l1->ridge_param_mult)).matrix();
        calFactor = tmpMat.array().square().matrix() * tmpVec; //N_cx1
        pred = tmpMat * (tmpVec.array() * zvec.array()).matrix() - (calFactor.array() * Y_chunk.array()).matrix();
        pred.array() /= (1 - calFactor.array());

        // compute mse and rsq
        l1->cumsum_values[0](ph, j) += pred.sum(); // Sx
                                                   // Y is centered so Sy = 0
        l1->cumsum_values[2](ph, j) += pred.squaredNorm(); // Sx2
                                                                    // Y is scaled so Sy2 = params->n_samples - ncov
        l1->cumsum_values[4](ph, j) += pred.dot(Y_chunk); // Sxy
        if(params->test_l0){ // pred = p1 + top_snp_pgs; Y is res pheno
          if(j == 0) l1->cumsum_values[3].row(ph).array() += Y_chunk.squaredNorm(); // (Y-PGS) is not standardized
          pred += l1->top_snp_pgs[0].col(ph).segment(j_start, size_chunk);
          l1->cumsum_values_full[0](ph, j) += pred.sum(); // Sx
          l1->cumsum_values_full[2](ph, j) += pred.squaredNorm(); // Sx2
          l1->cumsum_values_full[4](ph, j) += pred.dot(pheno_data->phenotypes.col(ph).segment(j_start, size_chunk)); // Sxy
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


// Logistic models
void ridge_logistic_level_1(struct in_files* files, struct param* params, struct phenodt* pheno_data, struct ridgel1* l1, vector<MatrixXb>& masked_in_folds, mstream& sout) {

  sout << endl << " Level 1 ridge with logistic regression..." << endl << flush;

  int niter_cur;
  int ph_eff;
  string in_pheno;
  ifstream infile;

  ArrayXd Y1, W1, p1, score;
  ArrayXd betaold, etavec, pivec, wvec, zvec, betanew, etatest;
  MatrixXd X1, XtW, XtWX, XtWZ;
  l1->pheno_l1_not_converged = ArrayXb::Constant(params->n_pheno, false);

  for (int i = 0; i < 6; i++)
    l1->cumsum_values[i].setZero(params->n_pheno, params->n_ridge_l1);

  for(int ph = 0; ph < params->n_pheno; ++ph ) {
    if( !params->pheno_pass(ph) ) continue;
    sout << "   -on phenotype " << ph+1 <<" (" << files->pheno_names[ph] <<")..." << flush;
    auto ts1 = std::chrono::high_resolution_clock::now();
    ph_eff = params->write_l0_pred ? 0 : ph;
    int bs_l1 = params->total_n_block * params->n_ridge_l0;

    // read in level 0 predictions from file
    if(params->write_l0_pred)
      read_l0(ph, ph_eff, files, params, l1, sout);
    check_l0(ph, ph_eff, params, l1, pheno_data, sout);
    bs_l1 = l1->test_mat[ph_eff][0].cols();
    MatrixXd ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);

    for(int i = 0; i < params->cv_folds; ++i ) {
      if( l1->pheno_l1_not_converged(ph) ) break;

      if( params->within_sample_l0 ){
        X1 = l1->pred_mat[ph][i];
        Y1 = l1->pred_pheno_raw[ph][i];
        W1 = l1->pred_offset[ph][i];
      }

      // starting values for each trait
      betaold = betanew = ArrayXd::Zero(bs_l1);

      for(int j = 0; j < params->n_ridge_l1; ++j ) {
        if( l1->pheno_l1_not_converged(ph) ) break;

        niter_cur = 0;
        // use warm starts (i.e. set final beta of previous ridge param 
        // as initial beta for current ridge param)
        betaold = betanew;

        while(niter_cur++ < params->niter_max_ridge){

          if(params->within_sample_l0) {
            etavec = W1 + (X1 * betaold.matrix()).array();
            pivec = 1 - 1/(etavec.exp() + 1);
            wvec = pivec * (1 - pivec);
            // check none of the values are 0
            if( ( wvec == 0 ).count() > 0 ){
              sout << "ERROR: Zeros occurred in Var(Y) during ridge logistic regression! (Try with --loocv)" << endl;
              l1->pheno_l1_not_converged(ph) = true;
              break;
            }
            zvec = (etavec - W1) + (Y1 - pivec) / wvec;
            XtW = X1.transpose() * wvec.matrix().asDiagonal();
            XtWX = params->tau[ph](j) * l1->ridge_param_mult.matrix().asDiagonal();
            XtWX += XtW * X1;
            betanew = XtWX.colPivHouseholderQr().solve(XtW * zvec.matrix()).array();
            // get the score
            etavec = W1 + (X1 * betanew.matrix()).array();
            pivec = 1 - 1/(etavec.exp() + 1);
            score = (X1.transpose() * (Y1 - pivec).matrix()).array() - params->tau[ph](j) * l1->ridge_param_mult * betanew;

          } else {

            XtWX = params->tau[ph](j) * l1->ridge_param_mult.matrix().asDiagonal();
            XtWZ = MatrixXd::Zero(bs_l1, 1);

            for(int k = 0; k < params->cv_folds; ++k ) {
              if( k != i) {

                // get w=p*(1-p) and check none of the values are 0
                get_pvec(etavec, pivec, betaold, l1->test_offset[ph][k].array(), l1->test_mat[ph_eff][k], params->numtol_eps);
                if( get_wvec(pivec, wvec, masked_in_folds[k].col(ph).array(), params->l1_ridge_eps) ){
                  sout << "ERROR: Zeros occurred in Var(Y) during ridge logistic regression! (Try with --loocv)" << endl;
                  l1->pheno_l1_not_converged(ph) = true;
                  break;
                }

                zvec = masked_in_folds[k].col(ph).array().select((etavec - l1->test_offset[ph][k].array()) + (l1->test_pheno_raw[ph][k].array() - pivec) / wvec, 0);

                XtW = l1->test_mat[ph_eff][k].transpose() * masked_in_folds[k].col(ph).array().select(wvec,0).matrix().asDiagonal();
                XtWX += XtW * l1->test_mat[ph_eff][k];
                XtWZ += XtW * zvec.matrix();
              }
            }
            if( l1->pheno_l1_not_converged(ph) ) break;

            betanew = XtWX.llt().solve(XtWZ).array();

            // start step-halving
            for( int niter_search = 1; niter_search <= params->niter_max_line_search_ridge; niter_search++ ){

              bool invalid_wvec = false;

              for(int k = 0; k < params->cv_folds; ++k ) {
                if( k != i) {
                  // get w=p*(1-p) and check none of the values are 0
                  get_pvec(etavec, pivec, betanew, l1->test_offset[ph][k].array(), l1->test_mat[ph_eff][k], params->numtol_eps);
                  invalid_wvec = get_wvec(pivec, wvec, masked_in_folds[k].col(ph).array(), params->l1_ridge_eps);
                  if( invalid_wvec ) break; // do another halving
                }
              }

              if( !invalid_wvec ) break;

              // halve step size
              betanew = (betaold + betanew) / 2;

            }

            // compute score
            score = ArrayXd::Zero(bs_l1);
            for(int k = 0; k < params->cv_folds; ++k ) {
              if( k != i) {
                // get w=p*(1-p) and check none of the values are 0
                get_pvec(etavec, pivec, betanew, l1->test_offset[ph][k].array(), l1->test_mat[ph_eff][k], params->numtol_eps);
                if( get_wvec(pivec, wvec, masked_in_folds[k].col(ph).array(), params->l1_ridge_eps) ){
                  sout << "ERROR: Zeros occurred in Var(Y) during ridge logistic regression! (Try with --loocv)" << endl;
                  l1->pheno_l1_not_converged(ph) = true;
                  break;
                }
                score += (l1->test_mat[ph_eff][k].transpose() * masked_in_folds[k].col(ph).array().select(l1->test_pheno_raw[ph][k].array() - pivec, 0).matrix()).array();
              }
            }
            score -= params->tau[ph](j) * l1->ridge_param_mult * betanew;


          }

          // stopping criterion
          if( (score.abs().maxCoeff() < params->l1_ridge_tol) || l1->pheno_l1_not_converged(ph)) break;

          betaold = betanew;
        }

        //cerr << "\nFold=" << i << " tau = " << params->tau[ph](j) << " beta=" << betanew.matrix().transpose().array() << endl;
        //if(i==1) exit(EXIT_FAILURE);

        if(niter_cur > params->niter_max_ridge){
          sout << "WARNING: Penalized logistic regression did not converge! (Increase --niter)\n";
          l1->pheno_l1_not_converged(ph) = true;
          break;
        } else if(l1->pheno_l1_not_converged(ph)) break;
        //sout << "Converged in "<< niter_cur << " iterations. Score max = " << score.abs().maxCoeff() << endl;


        etatest = l1->test_offset[ph][i].array() + (l1->test_mat[ph_eff][i] * betanew.matrix()).array();
        p1 = (1 - 1/(etatest.exp() + 1));

        if(!params->within_sample_l0) l1->beta_hat_level_1[ph][i].col(j) = betanew;


        // compute mse
        for(int l = 0; l < params->cv_sizes(i); l++){
          if(!masked_in_folds[i](l,ph)) continue;

          // if p is within eps of 0/1, set to eps/1-eps
          if( p1(l) < params->l1_ridge_eps ) p1(l) = params->l1_ridge_eps;
          else if( p1(l) > (1-params->l1_ridge_eps) ) p1(l) = 1 - params->l1_ridge_eps;

          l1->cumsum_values[0](ph,j) += p1(l); // Sx
          l1->cumsum_values[1](ph,j) += l1->test_pheno_raw[ph][i](l,0); // Sy
          l1->cumsum_values[2](ph,j) += p1(l) * p1(l); // Sx2
          l1->cumsum_values[3](ph,j) += l1->test_pheno_raw[ph][i](l,0) * l1->test_pheno_raw[ph][i](l,0); // Sy2
          l1->cumsum_values[4](ph,j) += p1(l) * l1->test_pheno_raw[ph][i](l,0); // Sxy
          l1->cumsum_values[5](ph,j) += compute_log_lik_bern(l1->test_pheno_raw[ph][i](l,0), p1(l)); // -LL
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


void ridge_logistic_level_1_loocv(struct in_files* files, struct param* params, struct phenodt* pheno_data, struct ests* m_ests, struct ridgel1* l1, mstream& sout) {

  sout << endl << " Level 1 ridge with logistic regression..." << flush;

  int ph_eff, bs_l1 = params->total_n_block * params->n_ridge_l0;
  double p1;
  string in_pheno;
  ifstream infile;
  ofstream ofile;

  ArrayXd beta, pivec, wvec, v2, pred;
  MatrixXd XtWX, V1, b_loo;
  LLT<MatrixXd> Hinv;
  l1->pheno_l1_not_converged = ArrayXb::Constant(params->n_pheno, false);
  for (int i = 0; i < 6; i++)
    l1->cumsum_values[i].setZero(params->n_pheno, params->n_ridge_l1);

  uint64 max_bytes = params->chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (bs_l1 * target_size) matrix ]
  int nchunk = ceil( params->cv_folds * bs_l1 * sizeof(double) * 1.0 / max_bytes );
  int j_start, chunk, size_chunk, target_size = params->cv_folds / nchunk;
  sout << (params->verbose ? to_string(nchunk) + " chunks..." : "" ) << endl;
  if(params->debug && (params->n_pheno == 1) )
    openStream(&ofile, files->out_file + "_beta.l1.txt", ios::out, sout);

  for(int ph = 0; ph < params->n_pheno; ++ph ) {
    if( !params->pheno_pass(ph) ) continue;

    sout << "   -on phenotype " << ph+1 << " (" << files->pheno_names[ph] <<")..." << flush;
    auto ts1 = std::chrono::high_resolution_clock::now();
    ph_eff = params->write_l0_pred ? 0 : ph;

    // read in level 0 predictions from file
    if(params->write_l0_pred)
      read_l0(ph, ph_eff, files, params, l1, sout);
    check_l0(ph, ph_eff, params, l1, pheno_data, sout);
    bs_l1 = l1->test_mat_conc[ph_eff].cols();

    MapArXd Y (pheno_data->phenotypes_raw.col(ph).data(), pheno_data->phenotypes_raw.rows());
    MapMatXd X (l1->test_mat_conc[ph_eff].data(), pheno_data->phenotypes_raw.rows(), bs_l1);
    MapArXd offset (m_ests->offset_nullreg.col(ph).data(), pheno_data->phenotypes_raw.rows());
    MapArXb mask (pheno_data->masked_indivs.col(ph).data(), pheno_data->masked_indivs.rows());

    // starting values for each trait
    beta = ArrayXd::Zero(bs_l1);
    for(int j = 0; j < params->n_ridge_l1; ++j ) {
      if(params->debug) cerr << "Ridge param #" << j+1 << " (=" << params->tau[ph](j) << ")\n";

      // using warm starts (i.e. set final beta of previous ridge param 
      // as initial beta for current ridge param)
      if( params->use_adam ) // run ADAM to get close to max
        run_log_ridge_loocv_adam(ph, params->tau[ph](j), l1->ridge_param_mult, beta, pivec, wvec, Y, X, offset, mask, params, sout);

      if(!run_log_ridge_loocv(params->tau[ph](j), l1->ridge_param_mult, target_size, nchunk, beta, pivec, wvec, Y, X, offset, mask, params, sout)){
        sout << "WARNING: Ridge logistic regression did not converge! (Increase --niter)\n";
        l1->pheno_l1_not_converged(ph) = true;
        break;
      }
      if(params->debug && (params->n_pheno == 1) ){
        IOFormat Fmt(FullPrecision, DontAlignCols, " ", "\n", "", "","","");
        ofile << beta.matrix().transpose().format(Fmt) << "\n";
      }

      // compute Hinv
      // zvec = (pheno_data->masked_indivs.col(ph).array()).select( (etavec - m_ests->offset_nullreg.col(ph).array()) + (pheno_data->phenotypes_raw.col(ph).array() - pivec) / wvec, 0);
      XtWX = params->tau[ph](j) * l1->ridge_param_mult.matrix().asDiagonal(); // compute XtWX in chunks
      for(chunk = 0; chunk < nchunk; ++chunk){
        size_chunk = ( chunk == nchunk - 1 ? params->cv_folds - target_size * chunk : target_size );
        j_start = chunk * target_size;

        Ref<MatrixXd> Xmat_chunk = X.block(j_start, 0, size_chunk, bs_l1); // n x k
        Ref<ArrayXd> w_chunk = wvec.segment(j_start, size_chunk);
        Ref<ArrayXb> mask_chunk = mask.segment(j_start, size_chunk);

        XtWX.noalias() += Xmat_chunk.transpose() * mask_chunk.select(w_chunk,0).matrix().asDiagonal() * Xmat_chunk;
      }
      Hinv.compute( XtWX );

      // LOOCV estimates
      for(chunk = 0; chunk < nchunk; ++chunk ) {
        size_chunk = ( chunk == nchunk - 1 ? params->cv_folds - target_size * chunk : target_size );
        j_start = chunk * target_size;

        Ref<MatrixXd> Xmat_chunk = X.block(j_start, 0, size_chunk, bs_l1); // n x k
        Ref<ArrayXd> Yvec_chunk = Y.segment(j_start, size_chunk);
        Ref<ArrayXb> mask_chunk = mask.segment(j_start, size_chunk);
        Ref<ArrayXd> p_chunk = pivec.segment(j_start, size_chunk);
        Ref<ArrayXd> w_chunk = wvec.segment(j_start, size_chunk);
        Ref<ArrayXd> off_chunk = offset.segment(j_start, size_chunk);

        V1 = Hinv.solve( Xmat_chunk.transpose() ); // k x n
        v2 = (Xmat_chunk.array() * V1.transpose().array()).rowwise().sum() * w_chunk;
        b_loo.resize(beta.size(), size_chunk);
        b_loo.array().colwise() = beta;
        b_loo -= V1 * ((Yvec_chunk - p_chunk)/(1-v2)).matrix().asDiagonal();
        pred = (Xmat_chunk.array() * b_loo.transpose().array()).rowwise().sum() + off_chunk;

        for(int i = 0; i < size_chunk; ++i ) {
          if(!mask_chunk(i)) continue;
          p1 = 1 - 1/ ( exp(pred(i)) + 1 );

          // if p is within eps of 0/1, set to eps/1-eps
          if( p1 < params->l1_ridge_eps ) p1 = params->l1_ridge_eps;
          else if( p1 > (1-params->l1_ridge_eps) ) p1 = 1 - params->l1_ridge_eps;

          // compute mse and rsq
          l1->cumsum_values[0](ph,j) += p1; // Sx
          l1->cumsum_values[1](ph,j) += Yvec_chunk(i,0); // Sy
          l1->cumsum_values[2](ph,j) += p1 * p1; // Sx2
          l1->cumsum_values[3](ph,j) += Yvec_chunk(i,0) * Yvec_chunk(i,0); // Sy2
          l1->cumsum_values[4](ph,j) += p1 * Yvec_chunk(i,0); // Sxy
          l1->cumsum_values[5](ph,j) += compute_log_lik_bern(Yvec_chunk(i,0), p1); // -LL
        }
      }

    }

    sout << "done";
    auto ts2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(ts2 - ts1);
    sout << " (" << duration.count() << "ms) "<< endl;
  }

  if(params->debug && (params->n_pheno == 1) )
    ofile.close();

  sout << endl;
}

bool run_log_ridge_loocv(const double& lambda, const Ref<const ArrayXd>& ridge_param_mult, const int& target_size, const int& nchunk, ArrayXd& betaold, ArrayXd& pivec, ArrayXd& wvec, const Ref<const ArrayXd>& Y, Ref<MatrixXd> X, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, struct param* params, mstream& sout) {

  bool dev_conv = false;
  int bs_l1 = X.cols();
  int niter_cur = 0, j_start, chunk, size_chunk;
  double fn_start = 0, fn_end = 0;
  ArrayXd etavec, betanew, score, step_size, vweights;
  MatrixXd XtWX;
  LLT<MatrixXd> Hinv;

  //// get starting values
  // get w=p*(1-p) and check none of the values are 0
  get_pvec(etavec, pivec, betaold, offset, X, params->numtol_eps);
  // get -2f(b)
  fn_start = get_logist_dev(Y, pivec, mask) + lambda * (ridge_param_mult * betaold).square().sum();
  if( get_wvec(pivec, wvec, mask, params->l1_ridge_eps) ){
    sout << "ERROR: Zeros occurred in Var(Y) during ridge logistic regression.\n";
    return false;
  }
  // get the score
  score = ( X.transpose() * mask.select(Y - pivec, 0).matrix()).array() ;
  score -= lambda * ridge_param_mult * betaold;
  // for convergence check
  //vweights = (X.array().square().matrix().transpose() * mask.select(wvec,0).matrix()).array();
  //vweights /= mask.select(wvec,0).sum();

  while(niter_cur++ < params->niter_max_ridge) {

    // get step size
    XtWX = lambda * ridge_param_mult.matrix().asDiagonal(); // compute XtWX in chunks
    for(chunk = 0; chunk < nchunk; ++chunk ) {
      size_chunk = ( chunk == nchunk - 1 ? params->cv_folds - target_size * chunk : target_size );
      j_start = chunk * target_size;

      Ref<MatrixXd> Xmat_chunk = X.block(j_start, 0, size_chunk, bs_l1); // n x k
      Ref<ArrayXd> w_chunk = wvec.segment(j_start, size_chunk); // nx1
      Ref<const ArrayXb> mask_chunk = mask.segment(j_start, size_chunk); //nx1

      XtWX.noalias() += Xmat_chunk.transpose() * mask_chunk.select(w_chunk,0).matrix().asDiagonal() * Xmat_chunk;
    }
    Hinv.compute( XtWX );
    step_size = Hinv.solve(score.matrix()).array();
    //cerr << (step_size * score).sum()/(.01 + abs(fn_end)) << "\n\n";

    // check f(b)
    for( int niter_search = 1; niter_search <= params->niter_max_line_search; niter_search++ ){

      betanew = betaold + step_size;

      // get w=p*(1-p) and check none of the values are 0
      get_pvec(etavec, pivec, betanew, offset, X, params->numtol_eps);
      // -2f(b)
      fn_end = get_logist_dev(Y, pivec, mask) + lambda * (ridge_param_mult * betanew).square().sum();
      if( get_wvec(pivec, wvec, mask, params->l1_ridge_eps) ){
        sout << "ERROR: Zeros occurred in Var(Y) during ridge logistic regression.\n";
        return false;
      }

      if(params->debug) cerr << "#" << niter_cur << "(#" << niter_search << "): " << setprecision(16) << fn_start << "->" << fn_end << "\n";

      if( fn_end < (fn_start + params->numtol) ) break;
      // adjusted step size
      step_size /= 2;
    }

    // get the score
    score = ( X.transpose() * mask.select(Y - pivec, 0).matrix()).array() ;
    score -= lambda * ridge_param_mult * betanew;
    if(params->debug) cerr << "#"<< niter_cur << ": score max = " << score.abs().maxCoeff() << ";dev_diff=" << setprecision(16) << abs(fn_end - fn_start)/(.01 + abs(fn_end)) << "\n";

    dev_conv = (abs(fn_end - fn_start)/(.01 + abs(fn_end)) < params->tol); // fractional change - same as glm
    if( score.abs().maxCoeff() < params->l1_ridge_tol ) 
      //( (vweights * (betanew - betaold).square()).maxCoeff() < params->tol)  // convergence criteria from glmnet
      break; // prefer for score to be below tol

    betaold = betanew;
    fn_start = fn_end;
    //vweights = (X.array().square().matrix().transpose() * mask.select(wvec,0).matrix()).array();
    //vweights /= mask.select(wvec,0).sum();
  }

  if( !dev_conv && (niter_cur > params->niter_max_ridge) )
    return false;

  betaold = betanew;
  return true;
}

// Ridge logistic with ADAM using mini batch
void run_log_ridge_loocv_adam(const int& ph, const double& lambda, const Ref<const ArrayXd>& ridge_param_mult, ArrayXd& betavec, ArrayXd& pivec, ArrayXd& wvec, const Ref<const ArrayXd>& Y, Ref<MatrixXd> X, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, struct param* params, mstream& sout) {

  int niter_cur = 0, index;
  double p_alpha = params->adam_alpha, p_beta1 = params->adam_beta1, p_beta2 = params->adam_beta2, p_eps = params->adam_eps, p_alpha_t;
  double eta, phat;
  //cerr << p_alpha << " " << p_beta1 << " " << p_beta2 << " " << p_eps << endl;
  std::uniform_int_distribution<> d(0, mask.count() - 1);
  std::mt19937 gen;
  ArrayXd etavec, gradient_f, mt, vt, step_size;

  // starting values for ADAM params
  mt = vt = betavec * 0;
  gradient_f.resize( betavec.size() );

  while(niter_cur++ < params->niter_max_ridge_adam) {

    gradient_f = lambda * ridge_param_mult * betavec;

    if(params->adam_mini){ // ADAM using mini-batch (only non-masked samples)

      for (int i = 0; i < params->adam_batch_size; i++){
        index = params->adam_indices[ph](d(gen));
        eta = offset(index) + X.row(index) * betavec.matrix();
        phat = 1 - 1/(exp(eta) + 1);
        gradient_f -= X.row(index).transpose().array() * (Y(index)-phat); 
      }
      gradient_f /= params->adam_batch_size;

    } else {

      get_pvec(etavec, pivec, betavec, offset, X, params->numtol_eps);
      gradient_f -= ( X.transpose() * mask.select(Y - pivec, 0).matrix()).array() ;

    }
    //if(niter_cur%100 == 1) sout << "At iteration #"<< niter_cur << "; score max = " << gradient_f.abs().maxCoeff() << endl;

    mt = p_beta1 * mt + (1 - p_beta1) * gradient_f;
    vt = p_beta2 * vt + (1 - p_beta2) * gradient_f.square();
    p_alpha_t = p_alpha * sqrt(1 - pow(p_beta2, niter_cur)) / (1 - pow(p_beta1, niter_cur));
    step_size = p_alpha_t * mt / (vt.sqrt() + p_eps);

    if( step_size.abs().maxCoeff() < params->numtol ) break;

    betavec -= step_size;

  }

  if(params->verbose) sout << "ADAM took "<< niter_cur << " iterations (score max = " << gradient_f.abs().maxCoeff() << ")...";

}

// Poisson regression
void ridge_poisson_level_1(struct in_files* files, struct param* params, struct phenodt* pheno_data, struct ridgel1* l1, vector<MatrixXb>& masked_in_folds, mstream& sout) {

  sout << endl << " Level 1 ridge with poisson regression..." << endl << flush;

  int niter_cur;
  int ph_eff;
  string in_pheno;
  ifstream infile;

  ArrayXd Y1, W1, p1, score;
  ArrayXd betaold, etavec, pivec, wvec, zvec, betanew, etatest;
  MatrixXd X1, XtW, XtWX, XtWZ;
  l1->pheno_l1_not_converged = ArrayXb::Constant(params->n_pheno, false);

  for (int i = 0; i < 6; i++)
    l1->cumsum_values[i].setZero(params->n_pheno, params->n_ridge_l1);

  for(int ph = 0; ph < params->n_pheno; ++ph ) {
    if( !params->pheno_pass(ph) ) continue;
    sout << "   -on phenotype " << ph+1 <<" (" << files->pheno_names[ph] <<")..." << flush;
    auto ts1 = std::chrono::high_resolution_clock::now();
    ph_eff = params->write_l0_pred ? 0 : ph;
    int bs_l1 = params->total_n_block * params->n_ridge_l0;

    // read in level 0 predictions from file
    if(params->write_l0_pred)
      read_l0(ph, ph_eff, files, params, l1, sout);
    check_l0(ph, ph_eff, params, l1, pheno_data, sout);
    bs_l1 = l1->test_mat[ph_eff][0].cols();
    MatrixXd ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);

    for(int i = 0; i < params->cv_folds; ++i ) {
      if( l1->pheno_l1_not_converged(ph) ) break;

      // starting values for each trait
      betaold = betanew = ArrayXd::Zero(bs_l1);

      for(int j = 0; j < params->n_ridge_l1; ++j ) {
        if( l1->pheno_l1_not_converged(ph) ) break;

        niter_cur = 0;
        // use warm starts (i.e. set final beta of previous ridge param 
        // as initial beta for current ridge param)
        betaold = betanew;

        while(niter_cur++ < params->niter_max_ridge){

          XtWX = params->tau[ph](j) * l1->ridge_param_mult.matrix().asDiagonal();
          XtWZ = MatrixXd::Zero(bs_l1, 1);

          for(int k = 0; k < params->cv_folds; ++k ) {
            if( k != i) {

              // get w=p*(1-p) and check none of the values are 0
              get_pvec_poisson(etavec, pivec, betaold, l1->test_offset[ph][k].array(), l1->test_mat[ph_eff][k], params->numtol_eps);
              if( (masked_in_folds[k].col(ph).array() && (pivec == 0) ).any() ){
                sout << "ERROR: Zeros occurred in Var(Y) during ridge poisson regression! (Try with --loocv)" << endl;
                l1->pheno_l1_not_converged(ph) = true;
                break;
              }

              zvec = masked_in_folds[k].col(ph).array().select((etavec - l1->test_offset[ph][k].array()) + (l1->test_pheno_raw[ph][k].array() - pivec) / pivec, 0);

              XtW = l1->test_mat[ph_eff][k].transpose() * masked_in_folds[k].col(ph).array().select(pivec,0).matrix().asDiagonal();
              XtWX += XtW * l1->test_mat[ph_eff][k];
              XtWZ += XtW * zvec.matrix();
            }
          }
          if( l1->pheno_l1_not_converged(ph) ) break;

          betanew = XtWX.llt().solve(XtWZ).array();

          // start step-halving
          for( int niter_search = 1; niter_search <= params->niter_max_line_search_ridge; niter_search++ ){

            bool invalid_pvec = false;

            for(int k = 0; k < params->cv_folds; ++k ) {
              if( k != i) {
                // get w=p*(1-p) and check none of the values are 0
                get_pvec_poisson(etavec, pivec, betanew, l1->test_offset[ph][k].array(), l1->test_mat[ph_eff][k], params->numtol_eps);
                invalid_pvec = (masked_in_folds[k].col(ph).array() && (pivec == 0) ).any();
                if( invalid_pvec ) break; // do another halving
              }
            }

            if( !invalid_pvec ) break;

            // halve step size
            betanew = (betaold + betanew) / 2;

          }

          // compute score
          score = ArrayXd::Zero(bs_l1);
          for(int k = 0; k < params->cv_folds; ++k ) {
            if( k != i) {
              // get w=p*(1-p) and check none of the values are 0
              get_pvec_poisson(etavec, pivec, betanew, l1->test_offset[ph][k].array(), l1->test_mat[ph_eff][k], params->numtol_eps);
              if( (masked_in_folds[k].col(ph).array() && (pivec == 0) ).any() ){
                sout << "ERROR: Zeros occurred in Var(Y) during ridge logistic regression! (Try with --loocv)" << endl;
                l1->pheno_l1_not_converged(ph) = true;
                break;
              }
              score += (l1->test_mat[ph_eff][k].transpose() * masked_in_folds[k].col(ph).array().select(l1->test_pheno_raw[ph][k].array() - pivec, 0).matrix()).array();
            }
          }
          score -= params->tau[ph](j) * l1->ridge_param_mult * betanew;

          // stopping criterion
          if( (score.abs().maxCoeff() < params->l1_ridge_tol) || l1->pheno_l1_not_converged(ph)) break;

          betaold = betanew;
        }

        if(params->debug && (ph==0) && (i==0)) cerr << "\nFold=" << i << " tau = " << params->tau[ph](j) << 
          " beta(" << niter_cur << ")=" << betanew.matrix().transpose().array() << " : " << score.abs().maxCoeff() << endl;

        if(niter_cur > params->niter_max_ridge){
          sout << "WARNING: Penalized poisson regression did not converge! (Increase --niter)\n";
          l1->pheno_l1_not_converged(ph) = true;
          break;
        } else if(l1->pheno_l1_not_converged(ph)) break;
        //sout << "Converged in "<< niter_cur << " iterations. ;

        get_pvec_poisson(etatest, p1, betanew, l1->test_offset[ph][i].array(), l1->test_mat[ph_eff][i], params->numtol_eps);
        l1->beta_hat_level_1[ph][i].col(j) = betanew;

        // compute mse
        for(int l = 0; l < params->cv_sizes(i); l++){
          if(!masked_in_folds[i](l,ph)) continue;

          // if p is within eps of 0, set to eps
          if( p1(l) < params->l1_ridge_eps ) p1(l) = params->l1_ridge_eps;

          l1->cumsum_values[0](ph,j) += p1(l); // Sx
          l1->cumsum_values[1](ph,j) += l1->test_pheno_raw[ph][i](l,0); // Sy
          l1->cumsum_values[2](ph,j) += p1(l) * p1(l); // Sx2
          l1->cumsum_values[3](ph,j) += l1->test_pheno_raw[ph][i](l,0) * l1->test_pheno_raw[ph][i](l,0); // Sy2
          l1->cumsum_values[4](ph,j) += p1(l) * l1->test_pheno_raw[ph][i](l,0); // Sxy
          l1->cumsum_values[5](ph,j) += compute_log_lik_poisson(l1->test_pheno_raw[ph][i](l,0), p1(l)); // -LL
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

void ridge_poisson_level_1_loocv(struct in_files* files, struct param* params, struct phenodt* pheno_data, struct ests* m_ests, struct ridgel1* l1, mstream& sout) {

  sout << endl << " Level 1 ridge with poisson regression..." << flush;

  int ph_eff, bs_l1 = params->total_n_block * params->n_ridge_l0;
  double v2, pred, p1;
  string in_pheno;
  ifstream infile;

  ArrayXd beta, pivec;
  MatrixXd XtWX, V1, b_loo;
  LLT<MatrixXd> Hinv;
  l1->pheno_l1_not_converged = ArrayXb::Constant(params->n_pheno, false);
  for (int i = 0; i < 6; i++)
    l1->cumsum_values[i].setZero(params->n_pheno, params->n_ridge_l1);

  uint64 max_bytes = params->chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (bs_l1 * target_size) matrix ]
  int nchunk = ceil( params->cv_folds * bs_l1 * sizeof(double) * 1.0 / max_bytes );
  int j_start, chunk, size_chunk, target_size = params->cv_folds / nchunk;
  sout << (params->verbose ? to_string(nchunk) + " chunks..." : "" ) << endl;

  for(int ph = 0; ph < params->n_pheno; ++ph ) {
    if( !params->pheno_pass(ph) ) continue;

    sout << "   -on phenotype " << ph+1 << " (" << files->pheno_names[ph] <<")..." << flush;
    auto ts1 = std::chrono::high_resolution_clock::now();
    ph_eff = params->write_l0_pred ? 0 : ph;

    // read in level 0 predictions from file
    if(params->write_l0_pred)
      read_l0(ph, ph_eff, files, params, l1, sout);
    check_l0(ph, ph_eff, params, l1, pheno_data, sout);
    bs_l1 = l1->test_mat[ph_eff][0].cols();
    MatrixXd ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);

    MapArXd Y (pheno_data->phenotypes_raw.col(ph).data(), pheno_data->phenotypes_raw.rows());
    MapMatXd X (l1->test_mat_conc[ph_eff].data(), pheno_data->phenotypes_raw.rows(), bs_l1);
    MapArXd offset (m_ests->offset_nullreg.col(ph).data(), pheno_data->phenotypes_raw.rows());
    MapArXb mask (pheno_data->masked_indivs.col(ph).data(), pheno_data->masked_indivs.rows());

    // starting values for each trait
    beta = ArrayXd::Zero(bs_l1);
    for(int j = 0; j < params->n_ridge_l1; ++j ) {

      // using warm starts (i.e. set final beta of previous ridge param 
      // as initial beta for current ridge param)
      if( params->use_adam ) // run ADAM to get close to max
        throw "not yet implemented"; //run_ct_ridge_loocv_adam(ph, params->tau[ph](j), beta, pivec, Y, X, offset, mask, params, sout);

      if(!run_ct_ridge_loocv(params->tau[ph](j), l1->ridge_param_mult, target_size, nchunk, beta, pivec, Y, X, offset, mask, params, sout)){
        sout << "WARNING: Ridge poisson regression did not converge! (Increase --niter)\n";
        l1->pheno_l1_not_converged(ph) = true;
        break;
      }

      // compute Hinv
      XtWX = MatrixXd::Zero(bs_l1, bs_l1);
      for(chunk = 0; chunk < nchunk; ++chunk){
        size_chunk = ( chunk == nchunk - 1 ? params->cv_folds - target_size * chunk : target_size );
        j_start = chunk * target_size;

        Ref<MatrixXd> Xmat_chunk = X.block(j_start, 0, size_chunk, bs_l1); // n x k
        Ref<MatrixXd> w_chunk = pivec.matrix().block(j_start, 0, size_chunk,1);
        Ref<MatrixXb> mask_chunk = mask.matrix().block(j_start, 0, size_chunk,1);

        XtWX += Xmat_chunk.transpose() * mask_chunk.array().select(w_chunk.array(),0).matrix().asDiagonal() * Xmat_chunk;
      }
      Hinv.compute( XtWX + params->tau[ph](j) * ident_l1 );

      // LOOCV estimates
      for(chunk = 0; chunk < nchunk; ++chunk ) {
        size_chunk = ( chunk == nchunk - 1 ? params->cv_folds - target_size * chunk : target_size );
        j_start = chunk * target_size;

        Ref<MatrixXd> Xmat_chunk = X.block(j_start, 0, size_chunk, bs_l1); // n x k
        Ref<MatrixXd> Yvec_chunk = Y.matrix().block(j_start, 0, size_chunk, 1);
        Ref<MatrixXb> mask_chunk = mask.matrix().block(j_start, 0, size_chunk,1);

        V1 = Hinv.solve( Xmat_chunk.transpose() ); // k x n
        for(int i = 0; i < size_chunk; ++i ) {
          if(!mask_chunk(i,0)) continue;
          v2 = Xmat_chunk.row(i) * V1.col(i);
          v2 *= pivec(j_start + i);
          b_loo = (beta - V1.col(i).array() * (Yvec_chunk(i,0) - pivec(j_start + i)) / (1 - v2)).matrix();
          pred = Xmat_chunk.row(i) * b_loo.col(0);
          pred += offset(j_start + i);
          p1 = exp(pred);

          // if p is within eps of 0, set to eps
          if( p1 < params->l1_ridge_eps ) p1 = params->l1_ridge_eps;

          // compute mse and rsq
          l1->cumsum_values[0](ph,j) += p1; // Sx
          l1->cumsum_values[1](ph,j) += Yvec_chunk(i,0); // Sy
          l1->cumsum_values[2](ph,j) += p1 * p1; // Sx2
          l1->cumsum_values[3](ph,j) += Yvec_chunk(i,0) * Yvec_chunk(i,0); // Sy2
          l1->cumsum_values[4](ph,j) += p1 * Yvec_chunk(i,0); // Sxy
          l1->cumsum_values[5](ph,j) += compute_log_lik_poisson(Yvec_chunk(i,0), p1); // -LL
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

bool run_ct_ridge_loocv(const double& lambda, const Ref<const ArrayXd>& ridge_param_mult, const int& target_size, const int& nchunk, ArrayXd& betaold, ArrayXd& pivec, const Ref<const ArrayXd>& Y, Ref<MatrixXd> X, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, struct param* params, mstream& sout) {

  int bs_l1 = X.cols();
  int niter_cur = 0, j_start, chunk, size_chunk;
  ArrayXd etavec, zvec, betanew, score;
  MatrixXd XtWX, XtWZ, V1;
  LLT<MatrixXd> Hinv;
  MatrixXd ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);


  while(niter_cur++ < params->niter_max_ridge) {

    // get p and check none of the values are 0
    get_pvec_poisson(etavec, pivec, betaold, offset, X, params->numtol_eps);
    if( (mask && (pivec == 0)).any() ){
      sout << "ERROR: Zeros occurred in Var(Y) during ridge logistic regression.\n";
      return false;
    }
    zvec = mask.select( (etavec - offset) + (Y - pivec) / pivec, 0);

    // compute XtWX and XtWZ in chunks
    XtWX = lambda * ridge_param_mult.matrix().asDiagonal(); 
    XtWZ = MatrixXd::Zero(bs_l1, 1);
    for(chunk = 0; chunk < nchunk; ++chunk ) {
      size_chunk = ( chunk == nchunk - 1 ? params->cv_folds - target_size * chunk : target_size );
      j_start = chunk * target_size;

      Ref<MatrixXd> Xmat_chunk = X.block(j_start, 0, size_chunk, bs_l1); // n x k
      Ref<MatrixXd> w_chunk = pivec.matrix().block(j_start, 0, size_chunk, 1);
      Ref<MatrixXd> z_chunk = zvec.matrix().block(j_start, 0, size_chunk,1);
      Ref<const MatrixXb> mask_chunk = mask.matrix().block(j_start, 0, size_chunk,1);

      V1 = Xmat_chunk.transpose() * mask_chunk.array().select(w_chunk.array(),0).matrix().asDiagonal();
      XtWX += V1 * Xmat_chunk;
      XtWZ += V1 * z_chunk;
    }
    Hinv.compute( XtWX );
    betanew = Hinv.solve(XtWZ).array();

    get_pvec_poisson(etavec, pivec, betanew, offset, X, params->numtol_eps);
    if( (mask && (pivec == 0)).any() ){
      sout << "ERROR: Zeros occurred in Var(Y) during ridge logistic regression.\n";
      return false;
    }

    // get the score
    score = ( X.transpose() * mask.select(Y - pivec, 0).matrix()).array() ;
    score -= lambda * ridge_param_mult * betanew;

    if( score.abs().maxCoeff() < params->l1_ridge_tol ) break;

    betaold = betanew;
  }

  if(niter_cur > params->niter_max_ridge) 
    return false;
  //sout << "Converged in "<< niter_cur << " iterations. Score max = " << score.abs().maxCoeff() << endl;

  betaold = betanew;
  return true;
}

bool get_wvec(ArrayXd& pivec, ArrayXd& wvec, const Ref<const ArrayXb>& mask, const double& tol){

  wvec = mask.select(pivec*(1-pivec), 1); // get_pvec below causes wvec to be at least 9.999800003056407e-06
  /*
  wvec = ArrayXd::Ones( mask.size() );// set all entries to 1
  // avoid 0 weights by setting w to eps when p is within eps of 0/1
  // (strategy used in glmnet)
  for (int i = 0; i < mask.size(); i++){
    if( !mask(i) ) continue;

    if( pivec(i) < tol) {
      pivec(i) = 0;
      wvec(i) = tol;
    } else if ( pivec(i) > (1-tol) ){
      pivec(i) = 1;
      wvec(i) = tol;
    } else wvec(i) = pivec(i) * (1-pivec(i));

  }
  //wvec = masks.col(ph).array().select(pivec * (1 - pivec), 1);
  */

  return (wvec == 0).any();
}

void get_pvec(ArrayXd& etavec, ArrayXd& pivec, const Ref<const ArrayXd>& beta, const Ref<const ArrayXd>& offset, const Ref<const MatrixXd>& Xmat, double const& eps){

  etavec = offset + (Xmat * beta.matrix()).array();
  get_pvec(pivec, etavec, eps);

}

void get_pvec(ArrayXd& etavec, ArrayXd& pivec, const double& beta, const Ref<const ArrayXd>& offset, const Ref<const VectorXd>& Xmat, double const& eps){

  etavec = offset + Xmat.array() * beta;
  get_pvec(pivec, etavec, eps);

}

void get_pvec(ArrayXd& pivec, const Ref<const ArrayXd>& etavec, double const& eps){

  // strategy used in glm
  pivec = (etavec > ETAMAXTHR).select( 1 /(1+eps),
      (etavec < ETAMINTHR).select( eps/(1+eps), 1 - 1/(etavec.exp() + 1) ) );
  //cerr << setprecision(16) << etavec.head(5) << "\n" << pivec.head(5) << "\n";

}

// for firth (ok if wvec=0)
void get_wvec(ArrayXd& pivec, ArrayXd& wvec, const Ref<const ArrayXb>& mask){
  wvec = mask.select(pivec*(1-pivec), 1);
}

void get_pvec_poisson(ArrayXd& etavec, ArrayXd& pivec, const Ref<const ArrayXd>& beta, const Ref<const ArrayXd>& offset, const Ref<const MatrixXd>& Xmat, double const& eps){

  etavec = offset + (Xmat * beta.matrix()).array();
  pivec = etavec.exp(); // lambda = E(Y)

}

double get_logist_dev(const Ref<const ArrayXd>& Y, const Ref<const ArrayXd>& pi, const Ref<const ArrayXb>& mask){

  double dev = 0;

  for( int i = 0; i < Y.size(); i++)
    if(mask(i)) dev += compute_log_lik_bern(Y(i), pi(i));

  return 2 * dev; // -2 log.lik
}

double get_poisson_dev(const Ref<const ArrayXd>& Y, const Ref<const ArrayXd>& pi, const Ref<const ArrayXb>& mask){

  double dev = 0;

  for(int i = 0; i < Y.size(); i++)
    if(mask(i)) dev += compute_log_lik_poisson(Y(i), pi(i));

  return 2 * dev; // -2 log.lik
}


double compute_log_lik_bern(const double& y, const double& p){
  // negative log likelihood for bernoulli
  return -((y==0) ? log(1-p) : log(p));
}

double compute_log_lik_poisson(const double& y, const double& p){
  // negative log likelihood for poisson
  return -(y * log(p) - p); // ignore log(y!) constant
}

double y_log_ypi(const double& y, const double& p){
  return (y != 0) ? y * log(y/p) : 0;
}

double get_deviance_logistic(const Ref<const ArrayXd>& Y, const Ref<const ArrayXd>& pi, const Ref<const ArrayXd>& wt, const Ref<const ArrayXb>& mask){

  double dev = 0;
  for(int i = 0; i < Y.size(); i++)
    if(mask(i)) {
      dev += 2 * wt(i) * ( y_log_ypi(Y(i), pi(i)) + y_log_ypi(1 - Y(i), 1 - pi(i)) );
      //cerr << "(" << wt(i) << ", " << Y(i) << ", " << pi(i) << ") - " << i << " -> " << dev << "\n";
      //if(i>10) exit(-1);
    }

  return dev;
}

void test_assoc_block(int const& chrom, int const& block, struct ridgel0& l0,struct ridgel1& l1, struct geno_block* Gblock, struct phenodt* pheno_data, snp const* snpinfo, struct param const& params, mstream& sout){

  sout << "   -extracting highly associated SNPs..." << flush;
  auto t2 = std::chrono::high_resolution_clock::now();

  bool run_algo = params.l0_snp_pval_thr > 0;
  ArrayXi ind_g_l1;
  MatrixXd tmpM;
  l0.picked_top_snp = MatrixXb::Constant(l0.GGt.cols(), params.n_pheno, false);

  if(run_algo){
    for(int ph = 0; ph < params.n_pheno; ph++)
      apply_iter_cond(chrom, block, ph, l0, l1, Gblock, snpinfo, params);
    sout << "number selected across phenotypes = [ " << l0.nspns_picked_block.matrix().transpose() << " ]...";
  }

  // discard variants that are picked across all traits
  ArrayXb rm_var = l0.picked_top_snp.rowwise().all();
  l0.subset_l0_snps_gmat = rm_var.any();
  l0.indices_gmat_keep = get_true_indices(!rm_var); // resize matrices to remove the picked SNPs
  if(l0.subset_l0_snps_gmat) {
    if(params.use_loocv){
      SelfAdjointEigenSolver<MatrixXd> esG(l0.GGt(l0.indices_gmat_keep,l0.indices_gmat_keep));
      l0.GGt_eig_vec = esG.eigenvectors();
      l0.GGt_eig_val = esG.eigenvalues();
      l0.Wmat = l0.GGt_eig_vec.transpose() * l0.GTY(l0.indices_gmat_keep, all);
    } else {
      l0.GGt = MatrixXd::Zero(l0.indices_gmat_keep.size(),l0.indices_gmat_keep.size());
      l0.GTY = MatrixXd::Zero(l0.indices_gmat_keep.size(),params.n_pheno);
      uint32_t cum_size_folds = 0;
      for(int i = 0; i < params.cv_folds; ++i ) {
        tmpM = l0.G_folds[i](l0.indices_gmat_keep,l0.indices_gmat_keep); l0.G_folds[i] = tmpM;
        l0.GGt += l0.G_folds[i];
        tmpM = l0.GtY[i](l0.indices_gmat_keep,all); l0.GtY[i] = tmpM;
        l0.GTY += l0.GtY[i];
        cum_size_folds += params.cv_sizes(i);
      }
      tmpM = Gblock->Gmat(l0.indices_gmat_keep, all); Gblock->Gmat = tmpM;
    }
  } else if(params.use_loocv){ // loocv
    SelfAdjointEigenSolver<MatrixXd> esG(l0.GGt);
    l0.GGt_eig_vec = esG.eigenvectors();
    l0.GGt_eig_val = esG.eigenvalues();
    l0.Wmat = l0.GGt_eig_vec.transpose() * l0.GTY;
  }

  sout << "done";
  auto t3 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
  sout << " (" << duration.count() << "ms) "<< endl;

}

void read_l0(int const& ph, int const& ph_eff, struct in_files* files, struct param* params, struct ridgel1* l1, mstream& sout){

  int start, np, bs_l1_tot = params->total_n_block * params->n_ridge_l0;
  string fin;

  // resize matrix
  if(params->use_loocv || params->trait_mode == 3)
    l1->test_mat_conc[ph_eff].resize(params->n_samples, bs_l1_tot);
  else for( int i = 0; i < params->cv_folds; ++i )
    l1->test_mat[ph_eff][i].resize(params->cv_sizes(i), bs_l1_tot);

  // all blocks in same file
  if(!params->run_l1_only){

    start = 0;
    np = bs_l1_tot;
    fin = files->loco_tmp_prefix;

    read_l0_chunk(ph, ph_eff, start, np, fin, params, l1, sout);

  } else { // blocks in separate file

    for(size_t i = 0; i < files->bstart.size(); i++){

      start = files->bstart[i] * params->n_ridge_l0;
      np = files->btot[i] * params->n_ridge_l0;
      fin = files->mprefix[i];

      read_l0_chunk(ph, ph_eff, start, np, fin, params, l1, sout);
    }

  }
}

// read in l0 predictors in columns [start,start+np)
void read_l0_chunk(int const& ph, int const& ph_eff, int const& start, int const& np, const string& prefix, struct param* params, struct ridgel1* l1, mstream& sout){

  string in_pheno = prefix + "_l0_Y" + to_string(ph+1);
  ifstream infile;
  openStream(&infile, in_pheno, ios::in | ios::binary, sout);

  if( getSize(in_pheno) != (sizeof(double) * params->n_samples * np ))
    throw "file " + in_pheno + " is not the right size." ;
  //cerr << in_pheno << "  " << getSize(in_pheno) << endl;

  // store back values in test_mat
  if(params->use_loocv || params->trait_mode == 3) {

    infile.read( reinterpret_cast<char *> (&l1->test_mat_conc[ph_eff](0, start)), params->n_samples * np * sizeof(double) );

    //if(ph == 0) sout << endl << "In:\n" << l1->test_mat_conc[ph_eff].block(0,0,5,6) << endl;
  } else {

    int nt = 0;

    for( int m = start; nt < np; nt++, m++ )
      for( int i = 0; i < params->cv_folds; ++i )
        for( int k = 0; k < params->cv_sizes(i); ++k )
          infile.read( reinterpret_cast<char *> (&l1->test_mat[ph_eff][i](k,m)), sizeof(double) );

  //if(start==0) cerr << endl <<l1->test_mat[ph_eff][0].block(0,0,3,3) << endl;

  }

  infile.close();
  
}

void check_l0(int const& ph, int const& ph_eff, struct param* params, struct ridgel1* l1, struct phenodt const* pheno_data, mstream& sout, bool const& silent_mode){

  l1->l0_colkeep.col(ph).array() = true;
  // adjustment to chr map for level 1 (number of level 0 predictors per chromosome)
  l1->chrom_map_ndiff.setZero();
  // to apply ridge parameter
  l1->ridge_param_mult = ArrayXd::Ones(l1->l0_colkeep.rows());

  // select blocks to retain min shrinkage
  if(params->select_l0){

    int N = l1->l0_pv_block.rows();
    MapArXd pv_arr (l1->l0_pv_block.col(ph).data(), N);
    vector<double> quantile_vec; quantile_vec.resize(N);
    VectorXd::Map(&quantile_vec[0], N) = pv_arr;
    if(!silent_mode) sout << "checking level 0 predictors for strength of associations with phenotype..." << flush;

    // sort p-values
    auto const Q1 = N / 4;
    auto const Q2 = N / 2;
    auto const Q3 = Q1 + Q2;
    double thr, beta_q, conf_alpha = 0.05/N; // to be conservative
    std::sort(quantile_vec.begin(), quantile_vec.end());
    if(!silent_mode && params->debug) cout << "[ Q1="<< quantile_vec[Q1] << ", Q2=" << quantile_vec[Q2] << ", Q3=" << quantile_vec[Q3] << " ]..." << flush;

    if (params->rm_l0_pct != 0){ // get the threshold to use & select blocks

      if(!silent_mode) sout << "removing the least " << params->rm_l0_pct << "% significant block level 0 predictors..." << flush;
      int N_rm = N * (params->rm_l0_pct/100);
      if(N_rm >= N) throw "invalid proportion specified for --rm-l0-pct";
      thr = quantile_vec[N_rm-1];
      for( int i = 0; i < N; ++i ) 
        if(pv_arr(i) <= thr) { // throw out all level 0 predictors from block
          l1->l0_colkeep.block(i * params->n_ridge_l0, ph, params->n_ridge_l0, 1).array() = false;
          // subtract from chr_map
          l1->chrom_map_ndiff( l1->chrom_block(i) - 1 ) += params->n_ridge_l0;
        }

    } else {

      // with U(0,1) independent p-values under H0, kth order statistic is Beta(k, N-k) 
      thr = quantile_vec[N-1] + 1e-6;
      for( int i = 0; i < (N-1); ++i ) {
        beta_distribution<>  bd(i+1, N - i - 1);
        beta_q = -log10( quantile(bd, conf_alpha/2.0) );
        if((quantile_vec[N-i-1] < beta_q) || (i==(N-2))){
          if(i>0) thr = quantile_vec[N-i];
          break;
        }
      }

      int Ntop = ceil(0.05*N);
      // threshold = median + 3*IQR
      // or 3 SD above mean?
      if((pv_arr >= thr).count() >= Ntop){
        double iqr = quantile_vec[Q3] - quantile_vec[Q1];
        thr = max(thr, quantile_vec[Q3] + 3 * iqr); 
      }

      // if too many are selected, keep top 1%
      if((pv_arr >= thr).count() >= Ntop){
        thr = quantile_vec[N-Ntop]; 
      }

      if(!silent_mode) sout << (pv_arr >= thr).count() << "/" << N << " blocks selected (Upper bound = " << thr << ")..." << flush;

      if( (pv_arr >= thr).any() )
        // go through each block
        for( int i = 0; i < N; ++i ) 
          if(pv_arr(i) >= thr) { // only use last ridge parameter (min shrinkage)
            l1->l0_colkeep.block(i * params->n_ridge_l0, ph, params->n_ridge_l0 - 1, 1).array() = false;
            l1->ridge_param_mult( (i+1) * params->n_ridge_l0 - 1 ) = 0;
            // subtract from chr_map
            l1->chrom_map_ndiff( l1->chrom_block(i) - 1 ) += params->n_ridge_l0 - 1;
          }
    }
  }

  // subset columns
  if(!l1->l0_colkeep.col(ph).all()){
    ArrayXi ind_keep = get_true_indices(l1->l0_colkeep.col(ph));

    if(params->use_loocv || params->trait_mode == 3){

      //cout << "\n\n" <<  ind_keep.matrix().transpose() << "\n" << l1->test_mat_conc[ph_eff].rows() << "," << l1->test_mat_conc[ph_eff].cols() << "\n";
      // update over row chunks to keep memory usage low
      uint64 max_bytes = params->chunk_mb * 1e6;
      int nchunk = ceil( l1->test_mat_conc[ph_eff].rows() * l1->test_mat_conc[ph_eff].cols() * sizeof(double) * 1.0 / max_bytes );
      int chunk, size_chunk, j_start, target_size = l1->test_mat_conc[ph_eff].rows() / nchunk, ncols_keep = l1->l0_colkeep.col(ph).count();

      for(chunk = 0; chunk < nchunk; ++chunk ) {
        size_chunk = chunk == nchunk - 1? l1->test_mat_conc[ph_eff].rows() - target_size * chunk : target_size;
        j_start = chunk * target_size;
        MatrixXd tmpMat = l1->test_mat_conc[ph_eff](seqN(j_start, size_chunk), ind_keep);
        l1->test_mat_conc[ph_eff].block(j_start, 0, size_chunk, ncols_keep) = tmpMat;
      }

      l1->test_mat_conc[ph_eff].conservativeResize(l1->test_mat_conc[ph_eff].rows(), ncols_keep);
      //cout << "\n\n" <<  ph_eff << ":" << l1->test_mat_conc[ph_eff].rows() << "," << l1->test_mat_conc[ph_eff].cols() << "\n";

    } else
      for( int i = 0; i < params->cv_folds; ++i ) {
        MatrixXd tmpMat = l1->test_mat[ph_eff][i](all, ind_keep);
        l1->test_mat[ph_eff][i] = tmpMat;
      }

    ArrayXd tmpvec = l1->ridge_param_mult( ind_keep );
    l1->ridge_param_mult = tmpvec;
  }

  // set ridge params
  if(!silent_mode) {
    if(params->trait_mode == 2){
      double rate = pheno_data->phenotypes_raw.col(ph).sum() / pheno_data->Neff(ph); // masked entries are 0
      params->tau[ph] = l1->l0_colkeep.col(ph).count() / (1 + params->tau[ph] / (rate * (1 - params->tau[ph]))).log();
      //cerr << endl << params->tau[i].matrix().transpose() << endl;
    } else if(params->trait_mode == 3){
      if (params->t2e_l1_pi6) {
        params->tau[ph] = l1->l0_colkeep.col(ph).count() * (1 - params->tau[ph]) / params->tau[ph];
        // Assuming input tau is total SNP heritability on the liability scale= m * 6/pi^2 * (1-h2) / h2
        params->tau[ph] *= 6 / (M_PI * M_PI);
      } else {
        Eigen::VectorXd index = Eigen::VectorXd::LinSpaced(params->n_ridge_l1, 0, params->n_ridge_l1 - 1);
        params->tau[ph] = ((index.array() / (params->n_ridge_l1 - 1)) * log(1e-6) + log(pheno_data->cox_max_tau[ph])).exp();
      }
    }else {
      params->tau[ph] = l1->l0_colkeep.col(ph).count() * (1 - params->tau[ph]) / params->tau[ph];
      // Assuming input tau is total SNP heritability on the liability scale= m * 3/pi^2 * (1-h2) / h2
      if(params->trait_mode == 1) params->tau[ph] *= 3 / (M_PI * M_PI);
    }
  }

}

uint64 getSize(string const& fname){
  
  struct stat stat_buf;
  int rc = stat(fname.c_str(), &stat_buf);

  return ( rc == 0 ? stat_buf.st_size : 0);

}

//////////// dev functions
void apply_iter_cond(int const& chrom, int const& block, int const& ph, struct ridgel0& l0, struct ridgel1& l1, struct geno_block* Gblock, snp const* snpinfo, struct param const& params){

  chi_squared chisq(1);
  int bs = l0.GGt.rows(), maxIndex = 0;
  if(bs == 1) return;
  double chisq_thr = quantile(complement(chisq, params.l0_snp_pval_thr)), r2_thr = 0.9, ss_y, ss_x1;
  double beta_top_snp, ggt_diag = l0.GGt(0,0); // G is residualized & scaled so GtG=N-K
  vector<int> top_indices;
  MapArXd GtY (l0.GTY.col(ph).data(), bs);
  ArrayXd block_top_pgs = ArrayXd::Zero(Gblock->Gmat.cols());
  ArrayXd snp_pgs, chisq_v, bvec, v_beta, bstart, v_y;
  MatrixXd tmpM, X2tX1_X1tX1_inv, X1tX1_inv_X1ty, LDmat = (l0.GGt.array() / (params.n_samples - params.ncov_analyzed) ).square().matrix();
  // initial values
  l0.nspns_picked_block(ph) = 0;
  bstart = GtY/ggt_diag; bvec = bstart;
  ss_y = l0.ymat_res.col(ph).squaredNorm();
  v_y = (ss_y - bvec * GtY) / (ggt_diag - l0.nspns_picked(ph) - 1);
  v_beta = v_y/ggt_diag;
  ArrayXi tmpVi, indices_start(bs); std::iota(indices_start.begin(), indices_start.end(), 0);
  
  for(int itr = 1; itr < bs; itr++){

    if( (!l0.picked_top_snp.array().col(ph)).count() == 1 ) break; // must have at least one non-picked SNP

    // run marginal tests
    chisq_v = bvec.square() / v_beta / v_y;
    if(chisq_v.maxCoeff(&maxIndex) < chisq_thr) break;

    // keep track of top SNP
    top_indices.push_back(indices_start(maxIndex));
    l0.picked_top_snp.array().col(ph)(indices_start(maxIndex)) = true;
    // get beta and update pgs
    beta_top_snp = bvec(maxIndex);
    if(params.debug) cout << "round " << itr << " - top SNP '" << snpinfo[indices_start(maxIndex)].ID << "'" << ": chisq=" << chisq_v(maxIndex) << "/beta=" << beta_top_snp << "/vy=" << v_y(maxIndex) << "\n";
    l0.nspns_picked_block(ph)++;
    snp_pgs = Gblock->Gmat.row(indices_start(maxIndex)).transpose().array() * beta_top_snp;
    block_top_pgs += snp_pgs; // update top_snp_pgs
    // ignore snps in high LD with top SNP
    l0.picked_top_snp.col(ph).array() = ( LDmat.col(indices_start(maxIndex)).array() > r2_thr ).select(true, l0.picked_top_snp.col(ph).array());
    // track indices of snps not picked
    indices_start = get_true_indices(!l0.picked_top_snp.array().col(ph));

    // quantities needed for beta &v(beta)
    /* // with eigen decomp
    SelfAdjointEigenSolver<MatrixXd> eig_x1tx1(l0.GGt(top_indices, top_indices));
    tmpM = eig_x1tx1.eigenvectors() * (1/eig_x1tx1.eigenvalues().array().sqrt()).matrix().asDiagonal();
    X2tX1_X1tX1_inv = l0.GGt(indices_start, top_indices) * tmpM;
    X1tX1_inv_X1ty = tmpM.transpose() * GtY(top_indices).matrix();
    // get bvec conditional on picked snps
    bvec = bstart(indices_start) - (X2tX1_X1tX1_inv * X1tX1_inv_X1ty).array() / ggt_diag;
    //if(params.debug) cerr << bvec.head(5).matrix().transpose() << "\n";
    v_y = ((l0.ymat_res.col(ph) - block_top_pgs.matrix()).squaredNorm() - bvec * GtY(indices_start)) / (ggt_diag - itr - 1);
    v_beta = (ggt_diag - X2tX1_X1tX1_inv.array().square().rowwise().sum()) / ggt_diag / ggt_diag;  
*/
    // switch to cholesky (X1tX1 should always be pd)
    if(itr == 1){
      X2tX1_X1tX1_inv = l0.GGt(indices_start,top_indices) / ggt_diag;
      ss_x1 = GtY(top_indices).square()(0) / ggt_diag;
    } else if(itr <5) {
      MatrixXd inv_X1tX1 = l0.GGt(top_indices, top_indices).inverse();
      X2tX1_X1tX1_inv = l0.GGt(indices_start,top_indices) * inv_X1tX1;
      //if(params.debug) cerr << (inv_X1tX1 * GtY(top_indices).matrix()).transpose() << "\n";
      ss_x1 = GtY(top_indices).matrix().transpose() * inv_X1tX1 * GtY(top_indices).matrix();
    } else {
      const LLT<MatrixXd> llt_X1tX1 = l0.GGt(top_indices, top_indices).llt();
      X2tX1_X1tX1_inv = llt_X1tX1.solve(l0.GGt(top_indices, indices_start)).transpose();
      //if(params.debug) cerr << llt_X1tX1.solve(GtY(top_indices).matrix()).transpose() << "\n";
      ss_x1 = (GtY(top_indices).matrix().transpose() * llt_X1tX1.solve(GtY(top_indices).matrix()))(0,0);
    }
    // get bvec conditional on picked snps
    bvec = bstart(indices_start) - (X2tX1_X1tX1_inv * GtY(top_indices).matrix()).array() / ggt_diag;
    //if(params.debug) cerr << bvec.head(5).matrix().transpose() << "\n";
    v_y = (ss_y - ss_x1 - bvec * GtY(indices_start)) / (ggt_diag - l0.nspns_picked(ph) - l0.nspns_picked_block(ph) - 1);
    v_beta = (ggt_diag - (X2tX1_X1tX1_inv.array() * l0.GGt(indices_start,top_indices).array()).rowwise().sum()) / ggt_diag / ggt_diag;  
  }

  if(params.debug) cout << "max Tchisq in l0 block (removing top SNPs) = " << chisq_v(maxIndex) << "\n";
  if(bs>1){ // if at least one snp was picked
    l1.top_snp_pgs[chrom].col(ph).array() += block_top_pgs; // update top_snp_pgs (sum per chromosome)
    l1.top_snp_pgs[0].col(ph).array() += block_top_pgs; // update top_snp_pgs (sum across chr)
    l0.ymat_res.col(ph).array() -= block_top_pgs;
    l0.nspns_picked(ph) += l0.nspns_picked_block(ph);
    if(params.use_loocv){ // loocv
      GtY = (Gblock->Gmat * l0.ymat_res.col(ph)).array();  
    } else {
      uint32_t cum_size_folds = 0; GtY = 0;
      for(int i = 0; i < params.cv_folds; ++i ) {
        l0.GtY[i].col(ph) = Gblock->Gmat.middleCols(cum_size_folds, params.cv_sizes(i)) * l0.ymat_res.col(ph).segment(cum_size_folds, params.cv_sizes(i));
        GtY += l0.GtY[i].col(ph).array();
        cum_size_folds += params.cv_sizes(i);
      }
    }
  }
}

void ridge_cox_level_1(struct in_files* files, struct param* params, struct phenodt* pheno_data, struct ridgel1* l1, struct ests* m_ests, mstream& sout) {
  sout << endl << " Level 1 ridge with cox regression..." << endl << flush;
  
  int ph_eff, l0_idx;
  int time_index, event_index;
  Eigen::VectorXd ph_time, ph_event;
  string in_pheno;
  ifstream infile;
  l1->pheno_l1_not_converged = ArrayXb::Constant(params->n_pheno, false);

  for (int i = 0; i < 6; i++)
    l1->cumsum_values[i].setZero(params->n_pheno, params->n_ridge_l1);

  for (const auto& entry: files->t2e_map) {
    const std::string& time_name = entry.first;
    const std::string& event_name = entry.second;
    // find time column index
    std::vector<std::string>::iterator it_time = std::find(files->pheno_names.begin(), files->pheno_names.end(), time_name);
    time_index = std::distance(files->pheno_names.begin(), it_time);
    ph_time = pheno_data->phenotypes_raw.col(time_index);
    // find event column index
    std::vector<std::string>::iterator it_event = std::find(files->pheno_names.begin(), files->pheno_names.end(), event_name);
    event_index = std::distance(files->pheno_names.begin(), it_event);
    ph_event = pheno_data->phenotypes_raw.col(event_index);

    params->pheno_pass(event_index) = false;
    
    sout << "   -on phenotype " << time_name << "..." << flush;

    auto ts1 = std::chrono::high_resolution_clock::now();
    ph_eff = params->write_l0_pred ? 0 : time_index;
    l0_idx = params->t2e_event_l0 ? event_index : time_index;
    if(params->write_l0_pred)
      read_l0(l0_idx, ph_eff, files, params, l1, sout);

    MapArXb mask (pheno_data->masked_indivs.col(time_index).data(), pheno_data->masked_indivs.rows());

    // find max lambda for each t2e trait
    survival_data survivalNullData;
    survivalNullData.setup(ph_time, ph_event, mask, true);
    // initialize at lambda 0, to find lambda_max
    cox_ridge coxRidge_null_lamb0(survivalNullData, l1->test_mat_conc[ph_eff], m_ests->offset_nullreg(all, time_index), mask, 0, params->niter_max, params->niter_max_line_search, params->numtol_cox);
    coxRidge_null_lamb0.coxGrad(survivalNullData);
    Eigen::VectorXd gradient = coxRidge_null_lamb0.get_gradient();
    double lambda_max = getCoxLambdaMax(l1->test_mat_conc[ph_eff], gradient);
    pheno_data->cox_max_tau[time_index] = lambda_max;

    check_l0(time_index, ph_eff, params, l1, pheno_data, sout);

    for(int i = 0; i < params->cv_folds; ++i ) {
      ArrayXb fold_train_mask = (l1->fold_id[time_index].array() != i) && mask;
      ArrayXb fold_test_mask = (l1->fold_id[time_index].array() == i) && mask;
      
      survival_data survivalData_fold;
      survivalData_fold.setup(ph_time, ph_event, fold_train_mask, true);
      cox_ridge_path coxRidgePath_fold(survivalData_fold, l1->test_mat_conc[ph_eff], m_ests->offset_nullreg.col(time_index), fold_train_mask, params->n_ridge_l1, 1e-4, params->tau[time_index], params->niter_max_ridge, params->niter_max_line_search_ridge, params->l1_ridge_tol, true);
      coxRidgePath_fold.fit(survivalData_fold, l1->test_mat_conc[ph_eff], m_ests->offset_nullreg.col(time_index), fold_train_mask);

      if (!coxRidgePath_fold.converge.all()) {
        l1->pheno_l1_not_converged(time_index) = true;
      }
      l1->beta_hat_level_1[time_index][i] = coxRidgePath_fold.beta_mx;

      // prediction (eta), and compute deviance on test set
      survival_data survivalData_test;
      survivalData_test.setup(ph_time, ph_event, fold_test_mask, true);
      for (int l = 0; l < params->tau[time_index].size(); ++l) {
        cox_ridge coxRidge_test(survivalData_test, l1->test_mat_conc[ph_eff], m_ests->offset_nullreg.col(time_index), fold_test_mask, params->tau[time_index](l), params->niter_max_ridge, params->niter_max_line_search_ridge, params->l1_ridge_tol, false, coxRidgePath_fold.beta_mx.col(l));
        l1->cumsum_values[5](time_index, l) += coxRidge_test.get_null_deviance();
      }
    }
    sout << "done";
    auto ts2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(ts2 - ts1);
    sout << " (" << duration.count() << "ms) "<< endl;
  }
  sout << endl;
}
