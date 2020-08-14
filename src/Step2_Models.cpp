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
#include "Files.hpp"
#include "Pheno.hpp"
#include "Data.hpp"

using namespace std;
using namespace Eigen;
using namespace boost;
using boost::math::normal;


// Firth
bool fit_firth_logistic(int chrom, int ph, bool null_fit, struct param* params, struct phenodt* pheno_data, struct ests* m_ests, struct f_ests* fest, mstream& sout) {
  // if firth is used, fit based on penalized log-likelihood

  int niter_cur, col_incl;
  double dev_old, dev_new=0, denum, mx, deviance_l0 = 0;
  int maxstep_firth = null_fit ? params->maxstep_null : params->maxstep;
  int niter_firth = null_fit ? params->niter_max_firth_null : params->niter_max_firth;

  ArrayXd Y1, hvec, mod_score;
  ArrayXd betaold, betanew, step_size, etavec, pivec, wvec, loco_offset, covar_offset;
  MatrixXd Xmat, XtW, XtWX;
  ColPivHouseholderQR<MatrixXd> qr, qrX;
  hvec = ArrayXd::Zero(params->n_samples);
  Y1 = pheno_data->phenotypes_raw.col(ph).array() * pheno_data->masked_indivs.col(ph).array().cast<double>();

  if(params->firth_approx){
    if(null_fit){
      Xmat = pheno_data->new_cov; // only covariates
    } else {
      Xmat = fest->covs_firth.rightCols(1); // only tested SNP
    }
    col_incl = Xmat.cols();
  } else {
    Xmat = fest->covs_firth; // covariates + tested SNP
    col_incl = Xmat.cols();
    if( null_fit ) col_incl--;
  }

  // mask individuals
  Xmat.array().colwise() *= pheno_data->masked_indivs.col(ph).array().cast<double>();
  loco_offset = m_ests->blups.col(ph).array() * pheno_data->masked_indivs.col(ph).array().cast<double>(); 
  // covariate effects added as offset in firth approx. (last entry of beta_null_firth = 0)
  if( params->firth_approx && !null_fit ) covar_offset = ((pheno_data->new_cov.array().colwise() * pheno_data->masked_indivs.col(ph).array().cast<double>()).matrix() * fest->beta_null_firth.block(0,ph, pheno_data->new_cov.cols(),1)).array(); 

  // with firth approx. => trial 1: use maxstep_null
  // trial 2 => use fallback options (increase maxstep & niter)
  for( size_t trial = 0; trial < 2; trial++){

    // starting values
    if(null_fit){

      betaold = ArrayXd::Zero(Xmat.cols());
      betaold(0) = ( 0.5 + Y1.sum())  / (pheno_data->Neff(ph) + 1);
      betaold(0) = log( betaold(0) / (1 - betaold(0) ));

      // LOCO prediction is offset
      betaold(0) -= loco_offset.mean();

    } else {

      if(params->firth_approx) betaold = ArrayXd::Zero(col_incl); // only estimating effect of tested SNP
      else betaold = fest->beta_null_firth.array();

    }
    betanew = ArrayXd::Zero( betaold.size() );

    // get the corresponding deviance
    etavec = (Xmat * betaold.matrix()).array();
    etavec += loco_offset;
    if( params->firth_approx && !null_fit ) etavec += covar_offset; 
    // fitted probabilities
    pivec = 1 - 1 / (etavec.exp() + 1) ;
    wvec = (pheno_data->masked_indivs.col(ph).array()).select( ( pivec * (1 - pivec) ).sqrt(), 0);
    XtW = Xmat.transpose() * wvec.matrix().asDiagonal();
    XtWX = XtW * XtW.transpose();
    qr.compute(XtWX);
    // use penalized log-lik
    dev_old = (pheno_data->masked_indivs.col(ph).array()).select( (Y1 * pivec.log() + (1-Y1) * (1-pivec).log() ), 0).sum();
    dev_old += 0.5 * qr.logAbsDeterminant();
    dev_old *= -2;

    // at niter=0 (i.e. betaSNP=0) this is null deviance
    if( !null_fit ) deviance_l0 = dev_old;

    // solve S'(beta) = S(beta) + X'(h*(0.5-p)) = 0
    niter_cur = 0;
    while(niter_cur++ < niter_firth){

      ////////  compute step size
      etavec = (Xmat * betaold.matrix()).array();
      etavec += loco_offset;
      if( params->firth_approx && !null_fit ) etavec += covar_offset; 

      pivec = 1 - 1 / (etavec.exp() + 1) ;
      wvec = (pheno_data->masked_indivs.col(ph).array()).select( ( pivec * (1 - pivec) ).sqrt(), 0);
      XtW = Xmat.transpose() * wvec.matrix().asDiagonal();
      XtWX = XtW * XtW.transpose();
      qr.compute(XtWX);
      if(!params->firth_approx && null_fit )  qrX.compute(XtWX.block(0, 0, col_incl, col_incl));

      // compute diag(H), H = U(U'U)^{-1}U', U = Gamma^(1/2)X
      hvec = (qr.solve(XtW).array() * XtW.array() ).colwise().sum();

      // modified score for beta
      mod_score = (Xmat.leftCols(col_incl).transpose() * (pheno_data->masked_indivs.col(ph).array()).select( Y1 - pivec + hvec * (0.5 - pivec), 0).matrix() ).array();

      // step size
      if(!params->firth_approx && null_fit )
        step_size = qrX.solve( mod_score.matrix() ).array();
      else
        step_size = qr.solve( mod_score.matrix() ).array();

      // force absolute step size to be less than maxstep for each entry of beta
      mx = step_size.abs().maxCoeff() / maxstep_firth;
      if( mx > 1 ) step_size /= mx;

      // start step-halving and stop when deviance decreases 
      denum = 1;
      for( size_t niter_search = 1; niter_search <= params->niter_max_line_search; niter_search++ ){

        // adjusted step size
        step_size /= denum;

        ///////// compute corresponding deviance
        betanew.head(col_incl) = betaold.head(col_incl) + step_size;
        etavec = (Xmat * betanew.matrix()).array();
        etavec += loco_offset;
        if( params->firth_approx && !null_fit ) etavec += covar_offset; 

        pivec = 1 - 1 / (etavec.exp() + 1) ;
        wvec = (pheno_data->masked_indivs.col(ph).array()).select( ( pivec * (1 - pivec) ).sqrt(), 0);
        XtW = Xmat.transpose() * wvec.matrix().asDiagonal();
        XtWX = XtW * XtW.transpose();
        qr.compute(XtWX);

        dev_new = (pheno_data->masked_indivs.col(ph).array()).select( (Y1 * pivec.log() + (1-Y1) * (1-pivec).log() ), 0).sum();
        dev_new += 0.5 * qr.logAbsDeterminant();
        dev_new *= -2;

        //sout << "\n["<<niter_cur << " - " << niter_search <<"]  denum =" << denum << ";\n step =" << step_size.matrix().transpose().array() / denum<<"; \nbeta=" << betanew.matrix().transpose().array() << ";\n Lnew= " << dev_new << " vs L0="<< dev_old << ";score="<< mod_score<< endl;
        if( dev_new < dev_old + params->numtol ) break;
        denum *= 2;
      }

      betaold.head(col_incl) += step_size;
      dev_old = dev_new;

      // stopping criterion using modified score function
      if( mod_score.abs().maxCoeff() < params->numtol_firth) break;

    }

    if(params->firth_approx && null_fit){ // only retry for firth approx null model
      if(!params->fix_maxstep_null) { // don't retry with user-given settings
        if( niter_cur > niter_firth ){ // if failed to converge
          sout << "WARNING: Logistic regression with Firth correction did not converge (maximum step size=" << maxstep_firth <<";maximum number of iterations=" << niter_firth<<").";
          maxstep_firth = params->retry_maxstep_firth;
          niter_firth = params->retry_niter_firth;
          if(trial == 0) sout << "Retrying with fallback parameters: (maximum step size=" << maxstep_firth <<";maximum number of iterations=" << niter_firth<<").\n";
          continue;
        }
      }
    }

    break;
  }

  // If didn't converge
  if(niter_cur > niter_firth){
    if(params->verbose && !params->firth_approx) sout << "WARNING: Logistic regression with Firth correction did not converge!\n";
    return false;
  }
  // sout << "\nNiter = " << niter_cur << " : " << mod_score.matrix().transpose() << endl;

  if(null_fit) {
    if(params->firth_approx) fest->beta_null_firth.block(0,ph,betaold.size(),1) = betaold.matrix();
    else fest->beta_null_firth = betaold.matrix();
  } else {
    // compute beta_hat & SE
    fest->bhat_firth = betaold.tail(1)(0);
    fest->se_b_firth = sqrt( qr.inverse().diagonal().tail(1)(0) );

    // compute LRT test stat. 
    fest->deviance_logistic = dev_new - deviance_l0;
    fest->deviance_logistic *= -1;

    if( fest->deviance_logistic < 0 ) return false;
  }

  return true;
}



void fit_firth_logistic_snp(int chrom, int ph, bool null_fit, struct param* params, struct phenodt* pheno_data, struct ests* m_ests, struct f_ests* fest, variant_block* block_info, mstream& sout) {
  // if firth is used, fit based on penalized log-likelihood

  int niter_cur, col_incl;
  double dl=0, dev_old, denum, mx, deviance_l0 = 0;
  int maxstep_firth = null_fit ? params->maxstep_null : params->maxstep;
  int niter_firth = null_fit ? params->niter_max_firth_null : params->niter_max_firth;
  ArrayXd mod_score, betaold, betanew, step_size, etavec, pivec, wvec, hvec, loco_offset, covar_offset;
  MatrixXd Xmat, XtW, XtWX;
  ColPivHouseholderQR<MatrixXd> qr, qrX;

  ArrayXd Y1 = pheno_data->phenotypes_raw.col(ph).array() * pheno_data->masked_indivs.col(ph).array().cast<double>();

  if(params->firth_approx){
    if(null_fit){
      Xmat = pheno_data->new_cov; // only covariates
    } else {
      Xmat = block_info->Geno.matrix(); // only tested SNP
    }
    col_incl = Xmat.cols();
  } else {
    Xmat = MatrixXd::Zero(params->n_samples, pheno_data->new_cov.cols() + 1); // covariates + tested SNP
    Xmat << pheno_data->new_cov, block_info->Geno.matrix();
    col_incl = Xmat.cols();
    if( null_fit ) col_incl--;
  }
  // mask individuals
  Xmat.array().colwise() *= pheno_data->masked_indivs.col(ph).array().cast<double>();
  loco_offset = m_ests->blups.col(ph).array() * pheno_data->masked_indivs.col(ph).array().cast<double>(); 
  // covariate effects added as offset in firth approx.
  if( params->firth_approx && !null_fit ) covar_offset = ((pheno_data->new_cov.array().colwise() * pheno_data->masked_indivs.col(ph).array().cast<double>()).matrix() * fest->beta_null_firth.block(0,ph,pheno_data->new_cov.cols(),1)).array(); 

  // starting values
  if(null_fit){

    // start at logit^-1(mean(Y))-mean(offset)
    betaold = ArrayXd::Zero(Xmat.cols()); // last entry in exact Firth is kept at 0
    betaold(0) = ( 0.5 + Y1.sum())  / (pheno_data->Neff(ph) + 1);
    betaold(0) = log( betaold(0) / (1 - betaold(0) ));

    // LOCO prediction is offset
    betaold(0) -= loco_offset.mean();

  } else {
    // start at 0
    if(params->firth_approx) betaold = ArrayXd::Zero(Xmat.cols()); 
    // start at estimate from null fit
    else betaold = block_info->beta_null_firth.col(0);
  }
  betanew = ArrayXd::Zero( betaold.size() );

  // get the corresponding deviance
  etavec = (Xmat * betaold.matrix()).array();
  etavec += loco_offset;
  if( params->firth_approx && !null_fit ) etavec += covar_offset; 
  // fitted probabilities
  pivec = 1 - 1 / (etavec.exp() + 1) ;
  wvec = (pheno_data->masked_indivs.col(ph).array()).select( ( pivec * (1 - pivec) ).sqrt(), 0);
  XtW = Xmat.transpose() * wvec.matrix().asDiagonal();
  XtWX = XtW * XtW.transpose();
  qr.compute(XtWX.block(0, 0, col_incl, col_incl));
  // use penalized log-lik
  dev_old = (pheno_data->masked_indivs.col(ph).array()).select( (Y1 * pivec.log() + (1-Y1) * (1-pivec).log() ), 0).sum();
  dev_old += 0.5 * qr.logAbsDeterminant();
  dev_old *= -2;

  // at niter=0 (i.e. betaSNP=0), this is null deviance
  if( !null_fit ) deviance_l0 = dev_old;

  // solve S'(beta) = S(beta) + X'(h*(0.5-p)) = 0
  niter_cur = 0;
  while(niter_cur++ < niter_firth){

    ////////  compute step size
    etavec = (Xmat * betaold.matrix()).array();
    etavec += loco_offset;
    if( params->firth_approx && !null_fit ) etavec += covar_offset; 

    pivec = 1 - 1 / (etavec.exp() + 1) ;
    wvec = (pheno_data->masked_indivs.col(ph).array()).select( ( pivec * (1 - pivec) ).sqrt(), 0);
    XtW = Xmat.transpose() * wvec.matrix().asDiagonal();
    XtWX = XtW * XtW.transpose();
    qr.compute(XtWX);
    if( !params->firth_approx && null_fit ) qrX.compute(XtWX.block(0, 0, col_incl, col_incl));
    // compute diag(H), H = U(U'U)^{-1}U', U = Gamma^(1/2)X
    hvec = (qr.solve(XtW).array() * XtW.array() ).colwise().sum();

    // modified score for beta
    mod_score = (Xmat.leftCols(col_incl).transpose() * (pheno_data->masked_indivs.col(ph).array()).select( Y1 - pivec + hvec * (0.5 - pivec), 0).matrix() ).array();

    // step size
    if( !params->firth_approx && null_fit )
      step_size = qrX.solve( mod_score.matrix() ).array();
    else
      step_size = qr.solve( mod_score.matrix() ).array();

    // force absolute step size to be less than maxstep for each entry of beta
    mx = step_size.abs().maxCoeff() / maxstep_firth;
    if( mx > 1 ) step_size /= mx;

    // start step-halving and stop when deviance decreases 
    denum = 1;
    for( size_t niter_search = 1; niter_search <= params->niter_max_line_search; niter_search++ ){

      // adjusted step size
      step_size /= denum;

      ///////// compute corresponding deviance
      betanew.head(col_incl) = betaold.head(col_incl) + step_size;
      etavec = (Xmat * betanew.matrix()).array();
      etavec += loco_offset;
      if( params->firth_approx && !null_fit ) etavec += covar_offset; 

      pivec = 1 - 1 / (etavec.exp() + 1) ;
      wvec = (pheno_data->masked_indivs.col(ph).array()).select( ( pivec * (1 - pivec) ).sqrt(), 0);
      XtW = Xmat.transpose() * wvec.matrix().asDiagonal();
      XtWX = XtW * XtW.transpose();
      qr.compute(XtWX);

      dl = (pheno_data->masked_indivs.col(ph).array()).select( (Y1 * pivec.log() + (1-Y1) * (1-pivec).log() ), 0).sum();
      dl += 0.5 * qr.logAbsDeterminant();
      dl *= -2;

      //sout << "\n["<<niter_cur << " - " << niter_search <<"]  denum =" << denum << ";\n step =" << step_size.matrix().transpose().array() / denum<<"; \nbeta=" << betanew.matrix().transpose().array() << ";\n Lnew= " << dl << " vs L0="<< dev_old << ";score="<< mod_score<< endl;
      if( dl < dev_old + params->numtol ) break;
      denum *= 2;
    }

    betaold.head(col_incl) += step_size;
    dev_old = dl;

    // stopping criterion using modified score function
    if( mod_score.abs().maxCoeff() < params->numtol_firth) break;

  }

  // If didn't converge
  if(niter_cur > niter_firth){
    if(params->verbose) sout << "WARNING: Logistic regression with Firth correction did not converge!\n";
    block_info->test_fail[ph] = true;
    return ;
  }
  // sout << "\nNiter = " << niter_cur << " : " << mod_score.matrix().transpose() << endl;


  if(null_fit) {
    if(!params->firth_approx) block_info->beta_null_firth = betaold.matrix();
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




// SPA (MT in Eigen)
double solve_K1(const double tval, const bool use_SPAfast, const double denum_t, const int snp, const int ph, const struct param* params, const struct ests* m_ests, const struct spa_ests* s_est, const struct geno_block* gblock, mstream& sout){

  int niter_cur;
  const int lambda = s_est->pos_score ? 1 : -1; // if score is negative, adjust K' and K''
  double min_x, max_x, t_old, f_old, t_new = -1, f_new, hess;

  niter_cur = 0;
  min_x = 0, max_x = std::numeric_limits<double>::infinity();
  t_old = 0;
  f_old = use_SPAfast ? compute_K1_fast(lambda * t_old, denum_t, snp, ph, m_ests, s_est, gblock) : compute_K1(lambda * t_old, ph, m_ests, s_est);
  f_old *= lambda;
  f_old -= tval; 

  while( niter_cur++ < params->niter_max_spa ){

    hess = use_SPAfast ? compute_K2_fast(lambda * t_old, denum_t, snp, ph, m_ests, s_est, gblock) : compute_K2(lambda * t_old, ph, m_ests, s_est);
    t_new = t_old - f_old / hess;
    f_new = use_SPAfast ? compute_K1_fast(lambda * t_new, denum_t, snp, ph, m_ests, s_est, gblock) : compute_K1(lambda * t_new, ph, m_ests, s_est);
    f_new *= lambda;
    f_new -= tval;

    if( fabs( f_new ) < params->tol_spa ) break;

    // update bounds on root
    if( t_new && (t_new > min_x) && (t_new < max_x) ){
      if( f_new > 0) max_x = t_new;
      else min_x = t_new;
    } else{ // bisection method if t_new went out of bounds and re-compute f_new
      t_new = ( min_x + max_x ) / 2;
      // if( fabs( min_x - t_new ) < params->tol_spa ) break;
      f_new = use_SPAfast ? compute_K1_fast(lambda * t_new, denum_t, snp, ph, m_ests, s_est, gblock) : compute_K1(lambda * t_new, ph, m_ests, s_est);
      f_new *= lambda;
      f_new -= tval;
    }

    t_old = t_new;
    f_old = f_new;
  }

  // If didn't converge
  if( niter_cur > params->niter_max_spa ){
    if(params->verbose) sout << "WARNING: SPA did not converge to root for K'(t)=s.\n";
    return params->missing_value_double;
  }
  //sout << "#iterations = " << niter_cur << "; f= " << f_new << endl;

  return t_new;
}

double compute_K(const double t, const int ph, const struct ests* m_ests, const struct spa_ests* s_est ){

  double val = (1 - m_ests->Y_hat_p.col(ph).array() + m_ests->Y_hat_p.col(ph).array() * ( t / s_est->val_c * s_est->Gmod ).exp() ).log().sum() - t * s_est->val_a / s_est->val_c;

  return val;
}

double compute_K_fast(const double t, const double denum_t, const int snp, const int ph, const struct ests* m_ests, const struct spa_ests* s_est, const struct geno_block* gblock){

  uint32_t index_j;
  double val = 0;

  for( std::size_t j = 0; j < gblock->non_zero_indices_G[snp].size(); ++j ) {
    index_j = gblock->non_zero_indices_G[snp][j];
    val += log( 1 - m_ests->Y_hat_p(index_j,ph) + m_ests->Y_hat_p(index_j,ph) * exp( t / s_est->val_c * s_est->Gmod(index_j)) );
  }
  val += -t * s_est->val_d / s_est->val_c + t * t / 2 / denum_t * s_est->val_b;

  return val;
}

double compute_K1(const double t, const int ph, const struct ests* m_ests, const struct spa_ests* s_est ){

  double val = ( ( s_est->Gmod * m_ests->Y_hat_p.col(ph).array() / s_est->val_c ) / ( m_ests->Y_hat_p.col(ph).array() + (1 - m_ests->Y_hat_p.col(ph).array()) * ( -t / s_est->val_c * s_est->Gmod).exp() ) ).sum();
  val -= s_est->val_a / s_est->val_c;

  return val;
}

double compute_K1_fast(const double t, const double denum_t, const int snp, const int ph, const struct ests* m_ests, const struct spa_ests* s_est, const struct geno_block* gblock ){

  uint32_t index_j;
  double val = 0;

  for( std::size_t j = 0; j < gblock->non_zero_indices_G[snp].size(); ++j ) {
    index_j = gblock->non_zero_indices_G[snp][j];
    val += ( s_est->Gmod(index_j) * m_ests->Y_hat_p(index_j,ph) / s_est->val_c ) / ( m_ests->Y_hat_p(index_j,ph) + (1 - m_ests->Y_hat_p(index_j,ph)) * exp( -t / s_est->val_c * s_est->Gmod(index_j)) );
  }
  val += -s_est->val_d / s_est->val_c + t / denum_t * s_est->val_b;

  return val;
}

double compute_K2(const double t, const int ph, const struct ests* m_ests, const struct spa_ests* s_est ){

  double val = ( ( s_est->Gmod.square() * m_ests->Gamma_sqrt.col(ph).array().square() / (s_est->val_c*s_est->val_c) * ( -t / s_est->val_c * s_est->Gmod).exp()) / ( m_ests->Y_hat_p.col(ph).array() + (1 - m_ests->Y_hat_p.col(ph).array()) * ( -t / s_est->val_c * s_est->Gmod).exp() ).square() ).sum();

  return val;
}

double compute_K2_fast(const double t, const double denum_t, const int snp, const int ph, const struct ests* m_ests, const struct spa_ests* s_est, const struct geno_block* gblock ){

  uint32_t index_j;
  double val = 0, denum;

  for( std::size_t j = 0; j < gblock->non_zero_indices_G[snp].size(); ++j ) {
    index_j = gblock->non_zero_indices_G[snp][j];
    denum = m_ests->Y_hat_p(index_j,ph) + (1 - m_ests->Y_hat_p(index_j,ph)) * exp( -t / s_est->val_c * s_est->Gmod(index_j));
    val += ( s_est->Gmod(index_j) * s_est->Gmod(index_j) * m_ests->Gamma_sqrt(index_j,ph) * m_ests->Gamma_sqrt(index_j,ph) * exp( -t / s_est->val_c * s_est->Gmod(index_j)) / s_est->val_c / s_est->val_c ) / (denum * denum);
  }
  val += s_est->val_b / denum_t;

  return val;
}

double get_SPA_pvalue(const double root, const double tval, const bool use_SPAfast, const double denum_t, const int snp, const int ph, const struct param* params, const struct ests* m_ests, const struct spa_ests* s_est, const struct geno_block* gblock, mstream& sout){

  const int lambda = s_est->pos_score ? 1 : -1; // if score is negative, adjust K' and K''
  double kval, k2val, wval, vval, rval, pvalue;
  normal nd(0,1);

  kval = use_SPAfast ? compute_K_fast(lambda * root, denum_t, snp, ph, m_ests, s_est, gblock) : compute_K(lambda * root, ph, m_ests, s_est);
  k2val = use_SPAfast ? compute_K2_fast(lambda * root, denum_t, snp, ph, m_ests, s_est, gblock) : compute_K2(lambda * root, ph, m_ests, s_est);

  wval = sqrt( 2 * ( root * tval - kval ) );
  vval = root * sqrt( k2val );
  if(vval == 0) return 0; // implies score = 0 so pval = 1

  rval = wval + log( vval / wval ) / wval;
  pvalue = cdf(complement(nd, rval ));
  pvalue *= 2;

  if(pvalue > 1) { // SPA can fail for SNPs with very low counts 
    if(params->verbose) sout << "WARNING: SPA correction failed (resulted in p-value > 1).\n";
    return params->missing_value_double;
    //return 0;  
  }

  return -log10( pvalue );
}



// SPA (MT in OpenMP)
double solve_K1_snp(const double tval, const int ph, const struct param* params, const struct ests* m_ests, variant_block* block_info, mstream& sout){

  int niter_cur;
  int lambda = block_info->pos_score ? 1 : -1; // if score is negative, adjust K' and K''
  double min_x, max_x, t_old, f_old, t_new = -1, f_new, hess;

  niter_cur = 0;
  min_x = 0, max_x = std::numeric_limits<double>::infinity();
  t_old = 0;
  f_old = block_info->fastSPA ? compute_K1_fast_snp(lambda * t_old, m_ests, block_info, ph) : compute_K1_snp(lambda * t_old, m_ests, block_info, ph);
  f_old *= lambda;
  f_old -= tval; 

  while( niter_cur++ < params->niter_max_spa ){

    hess = block_info->fastSPA ? compute_K2_fast_snp(lambda * t_old, m_ests, block_info, ph) : compute_K2_snp(lambda * t_old, m_ests, block_info, ph);
    t_new = t_old - f_old / hess;
    f_new = block_info->fastSPA ? compute_K1_fast_snp(lambda * t_new, m_ests, block_info, ph) : compute_K1_snp(lambda * t_new, m_ests, block_info, ph);
    f_new *= lambda;
    f_new -= tval;

    if( fabs( f_new ) < params->tol_spa ) break;

    // update bounds on root
    if( t_new && (t_new > min_x) && (t_new < max_x) ){
      if( f_new > 0) max_x = t_new;
      else min_x = t_new;
    } else{ // bisection method if t_new went out of bounds and re-compute f_new
      t_new = ( min_x + max_x ) / 2;
      // if( fabs( min_x - t_new ) < params->tol_spa ) break;
      f_new = block_info->fastSPA ? compute_K1_fast_snp(lambda * t_new, m_ests, block_info, ph) : compute_K1_snp(lambda * t_new, m_ests, block_info, ph);
      f_new *= lambda;
      f_new -= tval;
    }

    t_old = t_new;
    f_old = f_new;
  }

  // If didn't converge
  if( niter_cur > params->niter_max_spa ){
    if(params->verbose) sout << "WARNING: SPA did not converge to root for K'(t)=s.\n";
    return params->missing_value_double;
  }
  //sout << "#iterations = " << niter_cur << "; f= " << f_new << endl;

  return t_new;
}

double compute_K_snp(const double t, const struct ests* m_ests, variant_block* block_info, const int ph){

  double val = (1 - m_ests->Y_hat_p.col(ph).array() + m_ests->Y_hat_p.col(ph).array() * ( t / block_info->val_c * block_info->Gmod ).exp() ).log().sum() - t * block_info->val_a / block_info->val_c;

  return val;
}

double compute_K_fast_snp(const double t, const struct ests* m_ests, variant_block* block_info, const int ph){

  uint32_t index_j;
  double val = 0;

  for( std::size_t j = 0; j < block_info->n_non_zero; ++j ) {
    index_j = block_info->non_zero_indices[j];
    val += log( 1 - m_ests->Y_hat_p(index_j,ph) + m_ests->Y_hat_p(index_j,ph) * exp( t / block_info->val_c * block_info->Gmod(index_j)) );
  }
  val += -t * block_info->val_d / block_info->val_c + t * t / 2 / block_info->denum(ph) * block_info->val_b;

  return val;
}

double compute_K1_snp(const double t, const struct ests* m_ests, variant_block* block_info, const int ph){

  double val = ( ( block_info->Gmod * m_ests->Y_hat_p.col(ph).array() / block_info->val_c ) / ( m_ests->Y_hat_p.col(ph).array() + (1 - m_ests->Y_hat_p.col(ph).array()) * ( -t / block_info->val_c * block_info->Gmod).exp() ) ).sum();
  val -= block_info->val_a / block_info->val_c;

  return val;
}

double compute_K1_fast_snp(const double t, const struct ests* m_ests, variant_block* block_info, const int ph){

  uint32_t index_j;
  double val = 0;

  for( std::size_t j = 0; j < block_info->n_non_zero; ++j ) {
    index_j = block_info->non_zero_indices[j];
    val += ( block_info->Gmod(index_j) * m_ests->Y_hat_p(index_j,ph) / block_info->val_c ) / ( m_ests->Y_hat_p(index_j,ph) + (1 - m_ests->Y_hat_p(index_j,ph)) * exp( -t / block_info->val_c * block_info->Gmod(index_j)) );
  }
  val += -block_info->val_d / block_info->val_c + t / block_info->denum(ph) * block_info->val_b;

  return val;
}

double compute_K2_snp(const double t, const struct ests* m_ests, variant_block* block_info, const int ph){

  double val = ( ( block_info->Gmod.square() * m_ests->Gamma_sqrt.col(ph).array().square() / (block_info->val_c*block_info->val_c) * ( -t / block_info->val_c * block_info->Gmod).exp()) / ( m_ests->Y_hat_p.col(ph).array() + (1 - m_ests->Y_hat_p.col(ph).array()) * ( -t / block_info->val_c * block_info->Gmod).exp() ).square() ).sum();

  return val;
}

double compute_K2_fast_snp(const double t, const struct ests* m_ests, variant_block* block_info, const int ph){

  uint32_t index_j;
  double val = 0, denum;

  for( std::size_t j = 0; j < block_info->n_non_zero; ++j ) {
    index_j = block_info->non_zero_indices[j];
    denum = m_ests->Y_hat_p(index_j,ph) + (1 - m_ests->Y_hat_p(index_j,ph)) * exp( -t / block_info->val_c * block_info->Gmod(index_j));
    val += ( block_info->Gmod(index_j) * block_info->Gmod(index_j) * m_ests->Gamma_sqrt(index_j,ph) * m_ests->Gamma_sqrt(index_j,ph) * exp( -t / block_info->val_c * block_info->Gmod(index_j)) / block_info->val_c / block_info->val_c ) / (denum * denum);
  }
  val += block_info->val_b / block_info->denum(ph);

  return val;
}

double get_SPA_pvalue_snp(const double root, const double tval, const int ph, struct param* params, const struct ests* m_ests, variant_block* block_info, mstream& sout){

  int lambda = block_info->pos_score ? 1 : -1; // if score is negative, adjust K and K''
  double kval, k2val, wval, vval, rval, pvalue;
  normal nd(0,1);

  kval = block_info->fastSPA ? compute_K_fast_snp(lambda * root, m_ests, block_info, ph) : compute_K_snp(lambda * root, m_ests, block_info, ph);
  k2val = block_info->fastSPA ? compute_K2_fast_snp(lambda * root, m_ests, block_info, ph) : compute_K2_snp(lambda * root, m_ests, block_info, ph);

  wval = sqrt( 2 * ( root * tval - kval ) );
  vval = root * sqrt( k2val );
  if(vval == 0) return 0; // implies score = 0 so pval = 1

  rval = wval + log( vval / wval ) / wval;
  pvalue = cdf(complement(nd, rval ));
  pvalue *= 2;

  if(pvalue > 1) { // SPA can fail for SNPs with very low counts 
    block_info->test_fail[ph] = true;
    if(params->verbose) sout << "WARNING: SPA correction failed (resulted in p-value > 1).\n";
    return params->missing_value_double;
    //return 0;  
  }

  return -log10( pvalue );
}


