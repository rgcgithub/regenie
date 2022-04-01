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

#include "Regenie.hpp"
#include "Files.hpp"
#include "Geno.hpp"
#include "Pheno.hpp"
#include "NNLS.hpp"
#include "Joint_Tests.hpp"
#include "Step1_Models.hpp"
#include "Step2_Models.hpp"
#include "SKAT.hpp"

#include "qf/qfc.h"

using namespace Eigen;
using namespace std;
using namespace boost;
using boost::math::normal;
using boost::math::chi_squared;
using boost::math::non_central_chi_squared;
using boost::math::beta_distribution;

// numerical integration using quadpack
// global variable for SKAT-O if used
ArrayXd flipped_skato_rho = ArrayXd::Zero(1);
ArrayXd skato_Qmin_rho = ArrayXd::Zero(1);
ArrayXd skato_tau = ArrayXd::Zero(1);
VectorXd skato_lambdas = VectorXd::Zero(1);
double skato_muQ = 0;
double skato_sd = 0;
int skato_state = 0;

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Functions for SKAT/SKAT-O
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void update_vc_gmat(SpMat& mat, ArrayXd& weights, ArrayXd& weights_acat, ArrayXb& ur_ind, int const& start, int const& bs, struct param const& params, const Ref<const ArrayXb>& in_analysis, Ref<MatrixXd> Gmat, vector<variant_block> &all_snps_info, Ref<MatrixXb> Jmat){

  beta_distribution<>  dist(params.skat_a1, params.skat_a2);

  if(params.mask_loo){ // update dimensions
    mat.resize(mat.rows(), bs + Jmat.cols());
    mat.setZero();
    weights = ArrayXd::Zero(bs + Jmat.cols(), 1);
    weights_acat = weights;
    if(params.debug) cerr << "Updating VC gmat...";
  }

#if defined(_OPENMP)
      setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < bs; ++i) {

    MapArXd Gvec (Gmat.col(i).data(), Gmat.rows(), 1);
    double maf;

    if(Jmat.row(start + i).any()){ // if variant is in at least one mask
      // flip if af is above 0.5
      if( all_snps_info[i].af1 > 0.5 ) Gvec = (Gvec == -3).select(-3, 2 - Gvec);
      maf = min(all_snps_info[i].af1, 1 - all_snps_info[i].af1);

      // check if ultra-rare (if so set to 0)
      if(ur_ind(start + i)){
        Jmat.row(start+i).array() = false; // ignore variant for all sets
        Gvec = 0; // don't store the variant
        continue;
      }

      // impute missing with mean
      Gvec = (Gvec == -3).select(2 * maf, Gvec);
      // mask individuals
      Gvec *= in_analysis.cast<double>();
      // store SKAT weight
      weights(start + i) = pdf(dist, maf);
      weights_acat(start + i) = weights(start + i) * weights(start + i) * maf * (1-maf); // for acatv
    } else Gvec = 0; // otherwise set the column to 0

  }
#if defined(_OPENMP)
      setNbThreads(params.threads);
#endif

  mat.middleCols(start, bs) = Gmat.sparseView();
}

void compute_vc_masks(SpMat& mat, Ref<ArrayXd> weights, Ref<ArrayXd> weights_acat, SpMat& vc_rare_mask, Ref<MatrixXb> vc_rare_non_miss, const Ref<const MatrixXd>& X, struct ests const& m_ests, struct f_ests const& fest, const Ref<const MatrixXd>& yres,  const Ref<const MatrixXd>& yraw, const Ref<const MatrixXb>& masked_indivs, MatrixXb& Jmat, vector<variant_block> &all_snps_info, const Ref<const ArrayXb>& in_analysis, struct param const& params){

  prep_ultra_rare_mask(mat, weights, weights_acat, vc_rare_mask, vc_rare_non_miss, Jmat, in_analysis, params);

  if(params.trait_mode==0)
    compute_vc_masks_qt(mat, weights, weights_acat, X, yres, Jmat, all_snps_info, params);
  else if(params.trait_mode==1)
    compute_vc_masks_bt(mat, weights, weights_acat, X, m_ests, fest, yres, yraw, masked_indivs, Jmat, all_snps_info, params);
  else throw "not yet implemented";

}

void prep_ultra_rare_mask(SpMat& mat, Ref<ArrayXd> weights, Ref<ArrayXd> weights_acat, SpMat& rare_mask_mat, Ref<MatrixXb> rare_mask_non_miss, MatrixXb& Jmat, const Ref<const ArrayXb>& in_analysis, struct param const& params){

  int nsets = Jmat.cols();
  int bs = mat.cols() - nsets; // contains single variants + ultra-rare masks
  // check if ultra-rare mask is used in any of the sets
  if(rare_mask_mat.nonZeros() == 0) {
    if(params.debug) cerr << "No ultra-rare variants (MAC <= " << params.skat_collapse_MAC << ") present in any of the sets.\n";
    return;
  }

  ArrayXd gv;
  boost::math::beta_distribution<>  dist(params.skat_a1, params.skat_a2);
  double mean, maf;

  for(int iset = 0; iset < nsets; iset++){

    // set entries of individuals not included in the analysis to 0
    gv = rare_mask_mat.col(iset).cwiseProduct(in_analysis.matrix().cast<double>());

    // check if any UR variants were included
    Jmat(bs+iset, iset) = Jmat.col(iset).any() && (gv>0).count();
    if(!Jmat(bs+iset, iset)) continue;
    rare_mask_non_miss.col(iset).array() = rare_mask_non_miss.col(iset).array() && in_analysis; 

    // compute mean
    mean = gv.sum() / rare_mask_non_miss.col(iset).count();
    maf = min(mean/2, 1 - mean/2);

    // store SKAT weight
    weights(bs+iset) = pdf(dist, maf);
    weights_acat(bs+iset) = weights(bs+iset) * weights(bs+iset) * maf * (1-maf);

    if(params.debug) cerr << "set #" << iset+1 << "; rare_mask [mu,nZ,w,w_a] = [" << mean << "," << (gv>0).count() << "," << weights(bs+iset) << "," << weights_acat(bs+iset) << "]";

    // impute missing entries which were set to 0
    gv = (!in_analysis || rare_mask_non_miss.col(iset).array()).select(gv, mean);
    mat.col(bs + iset) = gv.matrix().sparseView();

    if(params.debug) cerr << "\n";

  }

}

/////////////////////
/////////////////////
///// QTs
/////////////////////
/////////////////////
void compute_vc_masks_qt(SpMat& mat, const Ref<const ArrayXd>& weights, const Ref<const ArrayXd>& weights_acat, const Ref<const MatrixXd>& X, const Ref<const MatrixXd>& yres, const Ref<const MatrixXb>& Jmat, vector<variant_block> &all_snps_info, struct param const& params){

  if(params.skato_rho.size() == 1)
    compute_vc_masks_qt_fixed_rho(mat, weights, weights_acat, X, yres, Jmat, all_snps_info, params.skato_rho(0), params.skat_tol, params.nl_dbl_dmin, params.vc_test, params.debug);
  else
    compute_vc_masks_qt(mat, weights, weights_acat, X, yres, Jmat, all_snps_info, params.skato_rho, params.skat_tol, params.nl_dbl_dmin, params.vc_test, params.debug);

}

// for a given rho value
void compute_vc_masks_qt_fixed_rho(SpMat& mat, const Ref<const ArrayXd>& weights, const Ref<const ArrayXd>& weights_acat, const Ref<const MatrixXd>& X, const Ref<const MatrixXd>& yres, const Ref<const MatrixXb>& Jmat, vector<variant_block> &all_snps_info, double const& rho, double const& skat_lambda_tol, double const& nl_dbl_dmin, uint const& vc_test, bool const& debug){

  bool with_acatv = CHECK_BIT(vc_test,0);
  bool with_skat = (vc_test>>1)&15;
  int jcol, n_pheno = yres.cols(), bs = mat.cols(), nnz;
  double q;
  ArrayXd D, w_mask;
  VectorXd lambdas;
  MatrixXd Qs, Qb, Kv, Svals, Kmat, Rsqrt, sum_stats, pvals, pv_mask;
  SelfAdjointEigenSolver<MatrixXd> es;

  Svals.resize(n_pheno, bs); // PxM
  Kmat.resize(bs, bs); // MxM

  // get score stats & kernel matrices
  compute_vc_mats_qt(Svals, Kmat, X, yres, mat, weights);
  mat.resize(0,0); // not needed anymore

  // SKAT for all masks & traits
  Qs = Svals.array().square().matrix(); // PxM
  // if using acat-v, get single variant p-values
  if(with_acatv) {
    pvals.resize(bs, n_pheno); // MxP
    get_single_pvs(pvals, (Qs * (weights!=0).select(1/Kmat.diagonal().array(),1).matrix().asDiagonal()).transpose()); 
  }
  Qs *= Jmat.cast<double>(); // P x Km

  // burden
  Qb = (Svals * Jmat.cast<double>()).array().square().matrix();
  if(debug) cerr << "Q_SKAT for all masks:\n" << Qs << "\nQ_BURDEN for all masks:\n" << Qb << "\n";

  // for now don't parallelize this as it causes issues with qfc lib
  // but should be ok since dimensions don't depend on N
  for(size_t imask = 0; imask < all_snps_info.size(); imask++){

    variant_block* block_info = &(all_snps_info[imask]);
    if(block_info->sum_stats_vc.size()>0) block_info->sum_stats_vc.clear();
    if(block_info->skip_for_vc) continue;
    sum_stats = MatrixXd::Constant(n_pheno, 2, -1); // chisq & logp

    // get index of mask in Jmat
    jcol = block_info->col_jmat_skat;
    if(jcol < 0) continue; // this should not happen though
    MapcArXb Jvec (Jmat.col(jcol).data(), Jmat.rows(), 1);
    nnz = Jvec.count();
    if(nnz == 0) continue;

    // get kernel matrix
    SpMat Jstar (nnz, Kmat.rows()); // KmxM
    Jstar.reserve(nnz);
    for(int i = 0, j = 0; i < Jvec.size(); i++)
      if(Jvec(i)) Jstar.insert(j++, i) = 1;

    // ACAT-V 
    if(with_acatv){
      pv_mask = Jstar * pvals; // KmxP
      w_mask = (Jstar * weights_acat.matrix()).array(); // Kmx1
      if(debug) cerr << "SV log10p:\n" << pv_mask.col(0).transpose() << "\nWsq:\n" << w_mask.matrix().transpose() << "\n\n";
      for(int ph = 0; ph < n_pheno; ph++)
        get_logp(get_acat(pv_mask.col(ph).array(), w_mask), sum_stats(ph, 1), sum_stats(ph, 0), nl_dbl_dmin);
      block_info->sum_stats_vc["ACATV"] = sum_stats;
      sum_stats.array() = -1; // reset
    }
    if(!with_skat) continue;

    // R_rho sqrt
    es.compute( MatrixXd::Ones(nnz, nnz) ); // psd matrix
    D = 1 - rho + rho * es.eigenvalues().array();
    Rsqrt = es.eigenvectors() * D.max(0).sqrt().matrix().asDiagonal() * es.eigenvectors().transpose(); // 0 ev is ok but not neg due to numerical error

    Kv = Rsqrt * Jstar * (Kmat * Jstar.transpose()) * Rsqrt;

    // get eigen values
    get_lambdas(lambdas, Kv, skat_lambda_tol);
    if(lambdas.size() == 0) continue;
    //cerr << lambdas.matrix().transpose() << "\n\n";

    // compute test statistic & p-value
    for(int ph = 0; ph < n_pheno; ph++) {
      q = (1 - rho) * Qs(ph, jcol) + rho * Qb(ph, jcol);
      if( (rho == 1) || (lambdas.size() == 1) ){ // burden or single variant
        sum_stats(ph, 0) = q / lambdas.tail(1)(0);
        get_logp(sum_stats(ph, 1), sum_stats(ph, 0)); 
      } else compute_skat_pv(sum_stats(ph, 1), sum_stats(ph, 0), q, lambdas, nl_dbl_dmin);
    }
    if( (sum_stats.col(1).array() >= 0).any() ){
      string test_name = (rho > 0 ? "SKAT-RHO" : "SKAT");
      block_info->sum_stats_vc[test_name] = sum_stats;
    }

  }

}

void compute_vc_masks_qt(SpMat& mat, const Ref<const ArrayXd>& weights, const Ref<const ArrayXd>& weights_acat, const Ref<const MatrixXd>& X, const Ref<const MatrixXd>& yres, const Ref<const MatrixXb>& Jmat, vector<variant_block> &all_snps_info, const Ref<const ArrayXd>& rho_vec, double const& skat_lambda_tol, double const& nl_dbl_dmin, uint const& vc_test, bool const& debug){

  bool with_acatv = CHECK_BIT(vc_test,0);
  bool with_omnibus = (vc_test>>2)&7; // any of the omnibus tests
  bool with_skato_int = CHECK_BIT(vc_test,2);
  bool with_skato_acat = CHECK_BIT(vc_test,3);
  bool with_acato = CHECK_BIT(vc_test,4);
  int jcol, n_pheno = yres.cols(), bs = mat.cols(), nnz, nrho = rho_vec.size();
  double ztz, q, minp;
  ArrayXd D, w_mask, p_acato;
  VectorXd lambdas, z_mean, zmZ;
  MatrixXd XtG, ZtZ, ZtUZ, ZtIUZ, Qs, Qb, Qopt, Kv, Svals, Kmat, Rsqrt, cvals, sum_stats, pvals, pv_mask;
  MatrixXd pvs_skato, chisq_skato, pvs_skato_acat, chisq_skato_acat, pvs_acato, chisq_acato;
  SpMat Gw;
  SelfAdjointEigenSolver<MatrixXd> es;

  Svals.resize(n_pheno, bs); // PxM
  Kmat.resize(bs, bs); // MxM
  cvals.resize(nrho, 5);
  skato_Qmin_rho.resize(nrho, 1);
  if(with_acato) p_acato.resize(nrho+1);
  flipped_skato_rho = 1 - rho_vec;

  // project covariates
  Gw = mat  * weights.matrix().asDiagonal(); // NxM
  XtG = X.transpose() * Gw; // CxM

  // get score stats
  // need to use Gresid (must have yres centered)
  Svals = yres.transpose() * Gw - (yres.transpose() * X) * XtG; // PxM
  mat.resize(0,0); // not needed anymore

  // get kernel & resid matrices
  Kmat = Gw.transpose() * Gw - XtG.transpose() * XtG; // ZtZ

  // SKAT for all masks & traits
  Qs = Svals.array().square().matrix(); // PxM
  // if using acat-v, get single variant p-values
  if(with_acatv) {
    pvals.resize(bs, n_pheno); // MxP
    get_single_pvs(pvals, (Qs * (weights!=0).select(1/Kmat.diagonal().array(),1).matrix().asDiagonal()).transpose()); 
  }
  Qs *= Jmat.cast<double>(); // P x Km

  // burden
  Qb = (Svals * Jmat.cast<double>()).array().square().matrix();
  if(debug) cerr << "Q_SKAT for all masks:\n" << Qs << "\nQ_BURDEN for all masks:\n" << Qb << "\n";

  // for now don't parallelize this as it causes issues with qfc lib
  // but should be ok since dimensions don't depend on N
  for(size_t imask = 0; imask < all_snps_info.size(); imask++){

    variant_block* block_info = &(all_snps_info[imask]);
    if(block_info->sum_stats_vc.size()>0) block_info->sum_stats_vc.clear();
    if(block_info->skip_for_vc) continue;
    pvs_skato = MatrixXd::Constant(n_pheno, nrho, -1);
    chisq_skato = MatrixXd::Constant(n_pheno, nrho, -1);
    if(with_skato_acat){
      pvs_skato_acat = MatrixXd::Constant(n_pheno, 1, -1);
      chisq_skato_acat = MatrixXd::Constant(n_pheno, 1, -1);
    }
    if(with_acato){
      pvs_acato = MatrixXd::Constant(n_pheno, 1, -1);
      chisq_acato = MatrixXd::Constant(n_pheno, 1, -1);
    }
    sum_stats = MatrixXd::Constant(n_pheno, 2, -1); // chisq & logp

    // get index of mask in Jmat
    jcol = block_info->col_jmat_skat;
    if(jcol < 0) continue; // this should not happen though
    MapcArXb Jvec (Jmat.col(jcol).data(), Jmat.rows(), 1);
    nnz = Jvec.count();
    if(nnz == 0) continue;

    // to get variants in mask
    SpMat Jstar (nnz, Kmat.rows()); // KmxM
    Jstar.reserve(nnz);
    for(int i = 0, j = 0; i < Jvec.size(); i++)
      if(Jvec(i)) Jstar.insert(j++, i) = 1;
    ZtZ = Jstar * (Kmat * Jstar.transpose()); // KmxKm

    // ACAT-V 
    if(with_acatv){
      pv_mask = Jstar * pvals; // KmxP
      w_mask = (Jstar * weights_acat.matrix()).array(); // Kmx1
      if(debug) cerr << "#in mask=" << pv_mask.rows() << "\nSV log10p:\n" << pv_mask.col(0).transpose() << "\nWsq:\n" << w_mask.matrix().transpose() << "\n\n";
      for(int ph = 0; ph < n_pheno; ph++)
        get_logp(get_acat(pv_mask.col(ph).array(), w_mask), sum_stats(ph, 1), sum_stats(ph, 0), nl_dbl_dmin);
      block_info->sum_stats_vc["ACATV"] = sum_stats;
      sum_stats.array() = -1; // reset
    }
    if(!with_omnibus) continue;

    z_mean = (Gw * Jvec.matrix().cast<double>()).rowwise().sum(); // Nx1
    z_mean -= X * (XtG * Jvec.matrix().cast<double>()).rowwise().sum(); 
    z_mean /= nnz;
    ztz = z_mean.squaredNorm();
    zmZ = ZtZ.rowwise().mean(); // Kmx1
    ZtUZ = zmZ * zmZ.transpose() / ztz; // KmxKm
    ZtIUZ = ZtZ - ZtUZ;
    get_lambdas(skato_lambdas, ZtIUZ, skat_lambda_tol);
    if(skato_lambdas.size() == 0) continue;
    if(debug) cerr << "L:\n" << skato_lambdas.transpose() << "\n";

    if(nnz > 1)
      get_skato_mom(skato_muQ, skato_sd, skato_tau, nnz, skato_lambdas, ztz, zmZ.squaredNorm(), ZtIUZ, ZtUZ, rho_vec, debug);

    Qopt = Qs.col(jcol) * flipped_skato_rho.matrix().transpose() + Qb.col(jcol) * rho_vec.matrix().transpose(); // P x Nrho
    if(debug) cerr << "Q:\n" << Qopt.row(0) << "\n";

    es.compute( MatrixXd::Ones(nnz, nnz) ); // psd matrix
    for(int j = 0; j < nrho; j++){

      // get kernel matrix
      D = flipped_skato_rho(j) + rho_vec(j) * es.eigenvalues().array();
      Rsqrt = es.eigenvectors() * D.max(0).sqrt().matrix().asDiagonal() * es.eigenvectors().transpose(); // 0 ev is ok but check not neg due to numerical error
      Kv = Rsqrt * ZtZ * Rsqrt;

      // get eigen values
      get_lambdas(lambdas, Kv, skat_lambda_tol);
      if(lambdas.size() == 0) break; // SKAT & SKAT-O failed

      // needed for skato (M>1)
      if(nnz > 1)  get_cvals(j, cvals, lambdas);

      for(int ph = 0; ph < n_pheno; ph++) {
        q = Qopt(ph, j);
        if( (rho_vec(j) == 1) || (lambdas.size() == 1) ){ // burden
          chisq_skato(ph, j) = q / lambdas.tail(1)(0);
          get_logp(pvs_skato(ph, j), chisq_skato(ph, j)); 
        } else compute_skat_pv(pvs_skato(ph, j), chisq_skato(ph, j), q, lambdas, nl_dbl_dmin);
      }

      // store SKAT results
      if(rho_vec(j)==0) {
        sum_stats.col(0) = chisq_skato.col(j);
        sum_stats.col(1) = pvs_skato.col(j);
      }

      if(nnz == 1) break; // sum stats same for all rhos
    }

    if((sum_stats.col(1).array() >= 0).any())
      block_info->sum_stats_vc["SKAT"] = sum_stats;

    if(nnz == 1) { // same p for all tests
      if((pvs_skato.col(0).array() < 0).all()) continue; // go to next mask set
      sum_stats.col(0) = chisq_skato.col(0);
      sum_stats.col(1) = pvs_skato.col(0);
      if(with_acato)
      block_info->sum_stats_vc["ACATO"] = sum_stats;
      if(with_skato_acat)
      block_info->sum_stats_vc["SKATO-ACAT"] = sum_stats;
      if(with_skato_int)
      block_info->sum_stats_vc["SKATO"] = sum_stats;
      continue;
    }

    // for each phenotype, check all SKATO pvs are defined
    if((pvs_skato.array() < 0).rowwise().any().all()) continue; // go to next mask set if no phenotype with all pvals defined

    // Get minimum of p-values and corresponding chisq quantile for each rho
    for(int ph = 0; ph < n_pheno; ph++) {
      if( (pvs_skato.row(ph).array() < 0).any() ) continue;

      if(with_skato_acat){
        if(debug) cerr << "skato-acat logp=" << pvs_skato.row(ph) <<"\n";
        get_logp(get_acat(pvs_skato.row(ph).array()), pvs_skato_acat(ph,0), chisq_skato_acat(ph, 0), nl_dbl_dmin);
      } 
      if(with_acato){ // include acatv pvalue
        p_acato(0) = block_info->sum_stats_vc["ACATV"](ph, 1);
        p_acato.tail(nrho) = pvs_skato.row(ph).transpose().array();
        if(debug) cerr << "acato logp=" << p_acato.matrix().transpose() <<"\n";
        get_logp(get_acat(p_acato), pvs_acato(ph,0), chisq_acato(ph, 0), nl_dbl_dmin);
      }
      if(with_skato_int){
        minp = pow(10, -pvs_skato.row(ph).maxCoeff());
        get_Qmin(nrho, minp, skato_Qmin_rho, cvals);
        if(debug) cerr << "Qmin=" << skato_Qmin_rho.matrix().transpose() << "\nlogp=" << pvs_skato.row(ph) <<"\n";
        // numerical integration
        get_skato_pv(pvs_skato(ph,0), chisq_skato(ph, 0), minp, nrho, nl_dbl_dmin, debug);
      }

    }

    // SKATO-ACAT
    if(with_skato_acat && !(pvs_skato_acat.col(0).array() < 0).all()){
      sum_stats.col(0) = chisq_skato_acat.col(0);
      sum_stats.col(1) = pvs_skato_acat.col(0);
      block_info->sum_stats_vc["SKATO-ACAT"] = sum_stats;
    }
    // ACATO
    if(with_acato && !(pvs_acato.col(0).array() < 0).all()){
      sum_stats.col(0) = chisq_acato.col(0);
      sum_stats.col(1) = pvs_acato.col(0);
      block_info->sum_stats_vc["ACATO"] = sum_stats;
    }
    // SKATO
    if(with_skato_int && !(pvs_skato.col(0).array() < 0).all()){
      sum_stats.col(0) = chisq_skato.col(0);
      sum_stats.col(1) = pvs_skato.col(0);
      block_info->sum_stats_vc["SKATO"] = sum_stats;
    }

  }

}

void compute_vc_mats_qt(Ref<MatrixXd> Svals, Ref<MatrixXd> Kmat, const Ref<const MatrixXd>& X, const Ref<const MatrixXd>& yres, Ref<SpMat> Gmat, const Ref<const ArrayXd>& weights){

  // project covariates
  MatrixXd GtX = Gmat.transpose() * X; // MxK

  // get test statistics (PxM matrix)
  // need to use Gresid (must have yres centered)
  Svals = ( yres.transpose() * Gmat - (yres.transpose() * X) * GtX.transpose() ) * weights.matrix().asDiagonal();

  // get kernel matrix (MxM matrix)
  Kmat = weights.matrix().asDiagonal() * (Gmat.transpose() * Gmat - GtX * GtX.transpose()) * weights.matrix().asDiagonal();

}

void get_single_pvs(Ref<MatrixXd> pvals, const Ref<const MatrixXd>& chisq_vals){
  for(int ph = 0; ph < pvals.cols(); ph++)
    for(int isnp = 0; isnp < pvals.rows(); isnp++)
      get_logp(pvals(isnp, ph), chisq_vals(isnp, ph));
  //cerr << pvals.leftCols(1) << "\n\n" << chisq_vals.leftCols(1) << endl;
}

/////////////////////
/////////////////////
///// BTs
/////////////////////
/////////////////////

void compute_vc_masks_bt(SpMat& mat, const Ref<const ArrayXd>& weights, const Ref<const ArrayXd>& weights_acat, const Ref<const MatrixXd>& X, struct ests const& m_ests, struct f_ests const& fest, const Ref<const MatrixXd>& yres, const Ref<const MatrixXd>& yraw, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXb>& Jmat, vector<variant_block> &all_snps_info, struct param const& params){

  if(params.skato_rho.size() == 1)
    compute_vc_masks_bt_fixed_rho(mat, weights, weights_acat, X, m_ests, fest, yres, yraw, masked_indivs, Jmat, all_snps_info, params.skato_rho(0), params.skat_tol, params.nl_dbl_dmin, params.firth || params.use_SPA, params.vc_test, params.debug, params);
  else 
    compute_vc_masks_bt(mat, weights, weights_acat, X, m_ests, fest, yres, yraw, masked_indivs, Jmat, all_snps_info, params.skato_rho, params.skat_tol, params.nl_dbl_dmin, params.firth || params.use_SPA, params.vc_test, params.debug, params);

}

void compute_vc_masks_bt_fixed_rho(SpMat& mat, const Ref<const ArrayXd>& weights, const Ref<const ArrayXd>& weights_acat, const Ref<const MatrixXd>& X, struct ests const& m_ests, struct f_ests const& fest, const Ref<const MatrixXd>& yres, const Ref<const MatrixXd>& yraw, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXb>& Jmat, vector<variant_block> &all_snps_info, double const& rho, double const& skat_lambda_tol, double const& nl_dbl_dmin, bool const& apply_correction, uint const& vc_test, bool const& debug, struct param const& params){

  bool with_acatv = CHECK_BIT(vc_test,0);
  bool with_skat = (vc_test>>1)&15;
  int jcol, n_pheno = yres.cols();
  double q;
  VectorXd lambdas, Qs, Qb;
  ArrayXd Svals, Rvec, Rw, D, w_mask, pvals;
  MatrixXd Kv, Kraw, Kmat, GtWX, Rsqrt, sum_stats, pv_mask;
  SpMat GWs;
  SelfAdjointEigenSolver<MatrixXd> es;

  MatrixXd pvs_m = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);
  MatrixXd chisq_m = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);
  MatrixXd pvs_m_a = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);
  MatrixXd chisq_m_a = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);
  sum_stats = MatrixXd::Constant(n_pheno, 2, -1); // chisq & logp

  for(int ph = 0; ph < n_pheno; ph++) { 
    if( !params.pheno_pass(ph) ) continue;

    MapcArXd Y (yraw.col(ph).data(), yraw.rows());
    MapcArXb mask (masked_indivs.col(ph).data(), yraw.rows());
    MapcArXd Wsqrt (m_ests.Gamma_sqrt.col(ph).data(), yraw.rows());
    MapcMatXd XWsqrt (m_ests.X_Gamma[ph].data(), yraw.rows(), m_ests.X_Gamma[ph].cols());
    MapcArXd phat (m_ests.Y_hat_p.col(ph).data(), yraw.rows());

    // multiply by sqrt(p(1-p)) and mask entries (NxM)
    GWs = (Wsqrt * mask.cast<double>()).matrix().asDiagonal() * mat; // NxM

    // get score stats for all variants (Mx1)
    Svals = (GWs.transpose() * yres.col(ph)).array();

    // kernel matrix for all variants (MxM)
    GtWX = XWsqrt.transpose() * GWs ; // KxM
    Kraw = GWs.transpose() * GWs - GtWX.transpose() * GtWX;

    // apply firth/spa corrections (set R=0 if failed)
    Rvec = (weights > 0).cast<double>();
    //cerr << weights.matrix().transpose() << "\n--\n" << Rvec.matrix().transpose() << "\n\n";
    if(apply_correction)
      apply_correction_cc(ph, Rvec, Svals, Kraw.diagonal().array(), mat, GtWX, XWsqrt, GWs, Wsqrt, phat, Y, mask, fest, params);
    if(with_acatv) {
      pvals.resize(Svals.size()); // Mx1
      get_single_pvs_bt(pvals, (Rvec!=0).select((Svals/Rvec).square() / Kraw.diagonal().array(),1)); 
    }

    // SKAT for all masks (Kmx1)
    Qs = Jmat.transpose().cast<double>() * (Rvec != 0).select(Svals * weights, 0).square().matrix();
    // burden
    Qb = (Jmat.transpose().cast<double>() * (Rvec != 0).select(Svals * weights, 0).matrix()).array().square().matrix();
    if(debug) cerr << "Q_SKAT for all masks:\n" << Qs.transpose() << "\nQ_BURDEN for all masks:\n" << Qb.transpose() << "\n";

    // apply correction factor
    Rw = Rvec * weights;
    Kmat = Rw.matrix().asDiagonal() * Kraw * Rw.matrix().asDiagonal();
    //cerr << "\n--\n" << Svals.matrix().transpose() << "\n--\n" << Rvec.matrix().transpose() << "\n\n" << Kmat.diagonal() << "\n\n";

    for(size_t imask = 0; imask < all_snps_info.size(); imask++){

      variant_block* block_info = &(all_snps_info[imask]);
      if(block_info->skip_for_vc) continue;

      // get index of mask in Jmat
      jcol = block_info->col_jmat_skat;
      if(jcol < 0) continue; // this should not happen though
      MapcArXb Jvec (Jmat.col(jcol).data(), Jmat.rows(), 1);
      int npass = (Jvec && (Rvec != 0)).count();
      if(npass == 0) continue;

      // extract rows/columns
      SpMat Jstar (npass, Kmat.rows()); // KmxM
      Jstar.reserve(npass);
      for(int i = 0, j = 0; i < Jvec.size(); i++)
        if(Jvec(i) && (Rvec(i) != 0)) Jstar.insert(j++, i) = 1;

      // ACAT-V 
      if(with_acatv){
        pv_mask = Jstar * pvals.matrix(); // Kmx1
        w_mask = (Jstar * weights_acat.matrix()).array(); // Kmx1
        if(debug) cerr << "SV log10p:\n" << pv_mask.col(0).transpose() << "\nWsq:\n" << w_mask.matrix().transpose() << "\n\n";
        get_logp(get_acat(pv_mask.col(0).array(), w_mask), pvs_m_a(ph, imask), chisq_m_a(ph, imask), nl_dbl_dmin);
      }
      if(!with_skat) continue;

      // R_rho sqrt
      es.compute( MatrixXd::Ones(npass, npass) ); // psd matrix
      D = 1 - rho + rho * es.eigenvalues().array();
      Rsqrt = es.eigenvectors() * D.max(0).sqrt().matrix().asDiagonal() * es.eigenvectors().transpose(); // 0 ev is ok but not neg due to numerical error

      // get kernel matrix
      Kv = Rsqrt * Jstar * (Kmat * Jstar.transpose()) * Rsqrt;

      // get eigen values
      get_lambdas(lambdas, Kv, skat_lambda_tol);
      if(lambdas.size() == 0) continue;
      //cerr << lambdas << "\n\n";

      // compute SKAT
      q = ((1 - rho) * Qs(jcol) + rho * Qb(jcol));
      if( (rho == 1) || (lambdas.size() == 1) ){ // burden or single variant
        chisq_m(ph, imask) = q / lambdas.tail(1)(0);
        get_logp(pvs_m(ph, imask), chisq_m(ph, imask)); 
      } else compute_skat_pv(pvs_m(ph, imask), chisq_m(ph, imask), q, lambdas, nl_dbl_dmin);
      //cerr << "\n\n" << p_skat << "\n\n";

    }
  }

  // store results
  string test_name = (rho > 0 ? "SKAT-RHO" : "SKAT");
  for(size_t imask = 0; imask < all_snps_info.size(); imask++){
    variant_block* block_info = &(all_snps_info[imask]);
    if(block_info->sum_stats_vc.size()>0) block_info->sum_stats_vc.clear();
    if((pvs_m_a.col(imask).array() >= 0).any()){// acatv
      sum_stats.col(0) = chisq_m_a.col(imask);
      sum_stats.col(1) = pvs_m_a.col(imask);
      block_info->sum_stats_vc["ACATV"] = sum_stats;
    }
    if((pvs_m.col(imask).array() >= 0).any()){
      sum_stats.col(0) = chisq_m.col(imask);
      sum_stats.col(1) = pvs_m.col(imask);
      block_info->sum_stats_vc[test_name] = sum_stats;
    }
  }

  mat.resize(0,0); // not needed anymore

}

void compute_vc_masks_bt(SpMat& mat, const Ref<const ArrayXd>& weights, const Ref<const ArrayXd>& weights_acat, const Ref<const MatrixXd>& X, struct ests const& m_ests, struct f_ests const& fest, const Ref<const MatrixXd>& yres, const Ref<const MatrixXd>& yraw, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXb>& Jmat, vector<variant_block> &all_snps_info, const Ref<const ArrayXd>& rho_vec, double const& skat_lambda_tol, double const& nl_dbl_dmin, bool const& apply_correction, uint const& vc_test, bool const& debug, struct param const& params){

  bool with_acatv = CHECK_BIT(vc_test,0);
  bool with_omnibus = (vc_test>>2)&7;
  bool with_skato_int = CHECK_BIT(vc_test,2);
  bool with_skato_acat = CHECK_BIT(vc_test,3);
  bool with_acato = CHECK_BIT(vc_test,4);
  int jcol, n_pheno = yres.cols(), nrho = rho_vec.size();
  double q, ztz, minp;
  VectorXd lambdas, Qs, Qb, Qopt, z_mean, zmZ;
  ArrayXd Svals, Rvec, Rw, D, pvs_skato, chisq_skato, w_mask, pvals, p_acato;
  MatrixXd ZtZ, ZtUZ, ZtIUZ, Kraw, Kmat, Kv, GtWX, Rsqrt, cvals, sum_stats, pv_mask;
  SpMat GWs;
  SelfAdjointEigenSolver<MatrixXd> es;

  MatrixXd pvs_m = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);//skat
  MatrixXd chisq_m = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);
  MatrixXd pvs_m_o = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);//skato
  MatrixXd chisq_m_o = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);
  MatrixXd pvs_m_o_acat = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);//skato-acat
  MatrixXd chisq_m_o_acat = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);
  MatrixXd pvs_m_acato = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);//acato
  MatrixXd chisq_m_acato = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);
  MatrixXd pvs_m_a = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);//acatv
  MatrixXd chisq_m_a = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);
  sum_stats = MatrixXd::Constant(n_pheno, 2, -1); // chisq & logp
  pvs_skato.resize(nrho);
  chisq_skato.resize(nrho);
  if(with_acato) p_acato.resize(nrho+1);
  cvals.resize(nrho, 5);
  skato_Qmin_rho.resize(nrho, 1);
  flipped_skato_rho = 1 - rho_vec;

  for(int ph = 0; ph < n_pheno; ph++) { 
    if( !params.pheno_pass(ph) ) continue;

    MapcArXd Y (yraw.col(ph).data(), yraw.rows());
    MapcArXb mask (masked_indivs.col(ph).data(), yraw.rows());
    MapcArXd Wsqrt (m_ests.Gamma_sqrt.col(ph).data(), yraw.rows());
    MapcMatXd XWsqrt (m_ests.X_Gamma[ph].data(), yraw.rows(), m_ests.X_Gamma[ph].cols());
    MapcArXd phat (m_ests.Y_hat_p.col(ph).data(), yraw.rows());

    // multiply by sqrt(p(1-p)) and mask entries (NxM)
    GWs = (Wsqrt * mask.cast<double>()).matrix().asDiagonal() * mat; // NxM

    // get score stats for all variants (Mx1)
    Svals = (GWs.transpose() * yres.col(ph)).array();

    // kernel matrix for all variants (MxM)
    GtWX = XWsqrt.transpose() * GWs ; // CxM
    Kraw = GWs.transpose() * GWs - GtWX.transpose() * GtWX; // ZtZ

    // apply firth/spa corrections (set R=0 if failed)
    Rvec = (weights > 0).cast<double>();
    //cerr << weights.matrix().transpose() << "\n--\n" << Rvec.matrix().transpose() << "\n\n";
    if(apply_correction)
      apply_correction_cc(ph, Rvec, Svals, Kraw.diagonal().array(), mat, GtWX, XWsqrt, GWs, Wsqrt, phat, Y, mask, fest, params);
    if(with_acatv) {
      pvals.resize(Svals.size()); // Mx1
      get_single_pvs_bt(pvals, (Rvec!=0).select((Svals/Rvec).square() / Kraw.diagonal().array(),1)); 
    }

    // SKAT for all masks (Kmx1)
    Qs = Jmat.transpose().cast<double>() * (Rvec != 0).select(Svals * weights, 0).square().matrix();
    // burden
    Qb = (Jmat.transpose().cast<double>() * (Rvec != 0).select(Svals * weights, 0).matrix()).array().square().matrix();
    if(debug) cerr << "Q_SKAT for all masks:\n" << Qs.transpose() << "\nQ_BURDEN for all masks:\n" << Qb.transpose() << "\n";

    // apply correction factor
    Rw = Rvec * weights;
    Kmat = Rw.matrix().asDiagonal() * Kraw * Rw.matrix().asDiagonal();
    //cerr << "\n--\n" << Svals.matrix().transpose() << "\n--\n" << Rvec.matrix().transpose() << "\n\n" << Kmat.diagonal() << "\n\n";

    for(size_t imask = 0; imask < all_snps_info.size(); imask++){

      variant_block* block_info = &(all_snps_info[imask]);
      if(block_info->skip_for_vc) continue;

      // get index of mask in Jmat
      jcol = block_info->col_jmat_skat;
      if(jcol < 0) continue; // this should not happen though
      MapcArXb Jvec (Jmat.col(jcol).data(), Jmat.rows(), 1);
      int npass = (Jvec && (Rvec != 0)).count();
      if(npass == 0) continue;

      // extract rows/columns
      SpMat Jstar (npass, Kmat.rows()); // KmxM
      Jstar.reserve(npass);
      for(int i = 0, j = 0; i < Jvec.size(); i++)
        if(Jvec(i) && (Rvec(i) != 0)) Jstar.insert(j++, i) = 1;

      // ACAT-V 
      if(with_acatv){
        pv_mask = Jstar * pvals.matrix(); // Kmx1
        w_mask = (Jstar * weights_acat.matrix()).array(); // Kmx1
        if(debug) cerr << "SV log10p:\n" << pv_mask.col(0).transpose() << "\nWsq:\n" << w_mask.matrix().transpose() << "\n\n";
        get_logp(get_acat(pv_mask.col(0).array(), w_mask), pvs_m_a(ph, imask), chisq_m_a(ph, imask), nl_dbl_dmin);
      }
      if(!with_omnibus) continue;

      ZtZ = Jstar * (Kmat * Jstar.transpose()); // KmxKm

      z_mean = (GWs * (Jvec && (Rvec != 0)).select(weights, 0).matrix()).rowwise().sum(); // Nx1
      z_mean -= XWsqrt * (GtWX * (Jvec && (Rvec != 0)).select(weights, 0).matrix()).rowwise().sum(); 
      z_mean /= npass;
      ztz = z_mean.squaredNorm();
      zmZ = ZtZ.rowwise().mean(); // Kmx1
      ZtUZ = zmZ * zmZ.transpose() / ztz;
      ZtIUZ = ZtZ - ZtUZ;
      get_lambdas(skato_lambdas, ZtIUZ, skat_lambda_tol);
      if(skato_lambdas.size() == 0) continue;
      if(debug) cerr << "L:\n" << skato_lambdas.transpose() << "\n";

      if(npass > 1)
        get_skato_mom(skato_muQ, skato_sd, skato_tau, npass, skato_lambdas, ztz, zmZ.squaredNorm(), ZtIUZ, ZtUZ, rho_vec, debug);

      Qopt = (Qs(jcol) * flipped_skato_rho + Qb(jcol) * rho_vec).matrix(); // Nrho x 1
      if(debug) cerr << "Q:\n" << std::setprecision(10) << Qopt.transpose() << "\n";

      es.compute( MatrixXd::Ones(npass, npass) ); // psd matrix
      for(int j = 0; j < nrho; j++){

        // get kernel matrix
        D = flipped_skato_rho(j) + rho_vec(j) * es.eigenvalues().array();
        Rsqrt = es.eigenvectors() * D.max(0).sqrt().matrix().asDiagonal() * es.eigenvectors().transpose(); // 0 ev is ok but check not neg due to numerical error
        Kv = Rsqrt * ZtZ * Rsqrt;

        // get eigen values
        get_lambdas(lambdas, Kv, skat_lambda_tol);
        if(lambdas.size() == 0) continue;
        //if(rho_vec(j) >0.9) cerr << "rho=" << rho_vec(j) << "\nL:"<<lambdas.matrix().transpose() << "\n\nQ=" << Qopt.col(j);

        // needed for skato (M>1)
        if(npass > 1)  get_cvals(j, cvals, lambdas);

        q = Qopt(j);
        if((rho_vec(j) == 1) || (lambdas.size() == 1) ){ // burden
          chisq_skato(j) = q / lambdas.tail(1)(0);
          get_logp(pvs_skato(j), chisq_skato(j)); 
        } else compute_skat_pv(pvs_skato(j), chisq_skato(j), q, lambdas, nl_dbl_dmin);

        // store SKAT results
        if(rho_vec(j)==0) {
          pvs_m(ph, imask) = pvs_skato(j);
          chisq_m(ph, imask) = chisq_skato(j);
        }

        if(npass == 1) break; // sum stats same for all rhos
      }

      if(npass == 1) { // same p for all tests
        if(pvs_skato(0) < 0) continue; // go to next mask set
        if(with_skato_int){
          pvs_m_o(ph, imask) = pvs_skato(0);
          chisq_m_o(ph, imask) = chisq_skato(0);
        }
        if(with_skato_acat){
          pvs_m_o_acat(ph, imask) = pvs_skato(0);
          chisq_m_o_acat(ph, imask) = chisq_skato(0);
        }
        if(with_acato){
          pvs_m_acato(ph, imask) = pvs_skato(0);
          chisq_m_acato(ph, imask) = chisq_skato(0);
        }
        continue;
      }

      // check pvs
      if((pvs_skato < 0).any()) continue;

      if(with_skato_acat){
        if(debug) cerr << "skato-acat logp=" << pvs_skato.matrix().transpose() <<"\n";
        get_logp(get_acat(pvs_skato), pvs_m_o_acat(ph, imask), chisq_m_o_acat(ph, imask), nl_dbl_dmin);
      } 
      if(with_acato){ // include acatv pvalue
        p_acato(0) = pvs_m_a(ph, imask);
        p_acato.tail(nrho) = pvs_skato;
        if(debug) cerr << "acato logp=" << p_acato.matrix().transpose() <<"\n";
        get_logp(get_acat(p_acato), pvs_m_acato(ph, imask), chisq_m_acato(ph, imask), nl_dbl_dmin);
      }
      if(with_skato_int){
        minp = pow(10, -pvs_skato.maxCoeff());
        get_Qmin(nrho, minp, skato_Qmin_rho, cvals);
        if(debug) cerr << "Qmin=" << skato_Qmin_rho.matrix().transpose() << "\nlogp=" << pvs_skato.matrix().transpose() <<"\n";

        // numerical integration
        get_skato_pv(pvs_m_o(ph, imask), chisq_m_o(ph, imask), minp, nrho, nl_dbl_dmin, debug);
      }

    }
  }

  // store skat results
  for(size_t imask = 0; imask < all_snps_info.size(); imask++){
    variant_block* block_info = &(all_snps_info[imask]);
    if(block_info->sum_stats_vc.size()>0) block_info->sum_stats_vc.clear();
    if((pvs_m_a.col(imask).array() >= 0).any()){// acatv
      sum_stats.col(0) = chisq_m_a.col(imask);
      sum_stats.col(1) = pvs_m_a.col(imask);
      block_info->sum_stats_vc["ACATV"] = sum_stats;
    }
    if((pvs_m.col(imask).array() >= 0).any()){//skat
      sum_stats.col(0) = chisq_m.col(imask);
      sum_stats.col(1) = pvs_m.col(imask);
      block_info->sum_stats_vc["SKAT"] = sum_stats;
    }
    if((pvs_m_o_acat.col(imask).array() >= 0).any()){// skato-acat
      sum_stats.col(0) = chisq_m_o_acat.col(imask);
      sum_stats.col(1) = pvs_m_o_acat.col(imask);
      block_info->sum_stats_vc["SKATO-ACAT"] = sum_stats;
    }
    if((pvs_m_acato.col(imask).array() >= 0).any()){//acato
      sum_stats.col(0) = chisq_m_acato.col(imask);
      sum_stats.col(1) = pvs_m_acato.col(imask);
      block_info->sum_stats_vc["ACATO"] = sum_stats;
    }
    if((pvs_m_o.col(imask).array() >= 0).any()){//skato
      sum_stats.col(0) = chisq_m_o.col(imask);
      sum_stats.col(1) = pvs_m_o.col(imask);
      block_info->sum_stats_vc["SKATO"] = sum_stats;
    }
  }

  mat.resize(0,0); // not needed anymore

}

// correcting for high cc imbalance
void apply_correction_cc(int const& ph, Ref<ArrayXd> Rvec, const Ref<const ArrayXd>& Tstats_sqrt, const Ref<const ArrayXd>& varT, SpMat const& Gsparse, const Ref<const MatrixXd>& GtWX, const Ref<const MatrixXd>& XWsqrt, SpMat const& GWs, const Ref<const ArrayXd>& Wsqrt, const Ref<const ArrayXd>& phat, const Ref<const ArrayXd>& Y, const Ref<const ArrayXb>& mask, struct f_ests const& fest, struct param const& params){

  int npass = Rvec.sum();

  // loop over the markers
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < Rvec.size(); i++) {

    bool test_fail = false;
    double chisq, pv;

    if(Rvec(i) == 0) continue;

    // if Tstat < threshold, no correction done
    double tstat_cur = Tstats_sqrt(i) / sqrt(varT(i));
    if(fabs(tstat_cur) <= params.z_thr) continue;

    MatrixXd Gres = - XWsqrt * GtWX.col(i); // get genotypic residuals
    Gres += GWs.col(i);

    // use SPA by default
    run_SPA_test_snp(chisq, pv, tstat_cur, varT(i), true, Gsparse.col(i), Gres.array(), phat, Wsqrt, mask, test_fail, params.tol_spa, params.niter_max_spa, params.missing_value_double);
    //cerr << "wSPA\n" ;

    // use firth as fallback if spa failed
    if(params.firth && test_fail){
      apply_firth_snp(test_fail, chisq, Gres.cwiseQuotient(Wsqrt.matrix()), Y, fest.cov_blup_offset.col(ph).array(), mask, params);
      //cerr << "wFirth\n" ;
    }

    if(params.debug) 
      cerr << "skat in: " << tstat_cur << " [=" << Tstats_sqrt(i) << " /sqrt(" << varT(i) << ")] -> " << chisq << endl;

    if( test_fail || (chisq == 0) ) { // set R to 0 for variant
      Rvec(i) = 0;
      continue;
    }

    Rvec(i) = fabs(tstat_cur) / sqrt(chisq);

  }
#if defined(_OPENMP)
  setNbThreads(params.threads);
#endif

  int npass_post = (Rvec > 0).count();
  if(params.verbose && (npass_post < npass)) cerr << "WARNING: correction failed for " << npass - npass_post << "/" << npass << " variants.";

}

// firth wrapper 
void apply_firth_snp(bool& fail, double& lrt, const Ref<const MatrixXd>& Gvec, const Ref<const ArrayXd>& Y, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, struct param const& params) {

  double dev0 = 0, dev;
  ArrayXd betaold, se, etavec, pivec;
  betaold = ArrayXd::Zero(1); // start at 0

  fail = !fit_firth_nr(dev0, Y, Gvec, offset, mask, pivec, etavec, betaold, se, 1, dev, true, lrt, params.maxstep, params.niter_max_firth, params.numtol_firth, &params);

}

void get_single_pvs_bt(Ref<ArrayXd> pvals, const Ref<const ArrayXd>& chisq_vals){
  for(int isnp = 0; isnp < pvals.rows(); isnp++)
    get_logp(pvals(isnp), chisq_vals(isnp));
  //cerr << pvals.matrix().transpose() << "\n\n" << chisq_vals.matrix().transpose() << endl;
}


/////////////////////
/////////////////////
///// General
/////////////////////
/////////////////////

void get_lambdas(VectorXd& lambdas, const Ref<const MatrixXd>& K, const double& tol){

  if(K.rows() == 1) { // K is scalar (e.g. single variant in set)
    lambdas = VectorXd::Constant(1, K(0,0));
    return;
  }

  // eigenvalues sorted in increasing order
  SelfAdjointEigenSolver<MatrixXd> esK(K, EigenvaluesOnly);
  // ignore zero eigen-values (matrix is psd)
  //int nonzero = (esK.eigenvalues().array() > esK.eigenvalues().tail(1)(0) * tol).count();
  // use filter strategy in R SKAT
  int nng = (esK.eigenvalues().array() >= 0).count();
  int nonzero = (esK.eigenvalues().array() > ( (esK.eigenvalues().array() >= 0).select(esK.eigenvalues().array(),0).sum() / nng * tol) ).count();
  lambdas = esK.eigenvalues().tail(nonzero);
}

void compute_skat_pv(double& pval, double& chival, double const& Q, VectorXd& lambdas, const double& tol){

  double pv;

  pv = get_chisq_mix_pv(Q, lambdas);
  //cerr << "mixture pv = " << pv << "\n";

  if(pv <= 0) { // spa also failed
    pval = -1; chival = -1;
  } else // take log10 and get chisq quantile
    get_logp(pv, pval, chival, tol);

}

// returns p-value or -1
double get_chisq_mix_pv(double const& q, const Ref<const VectorXd>& lambdas){

  double pv;

  // re-scale so that max lambda is 1 (lambda is sorted)
  double newQ = q / lambdas.tail(1)(0);
  VectorXd newL = lambdas / lambdas.tail(1)(0);

  // exact
  pv = get_davies_pv(newQ, newL);
  if(pv <= 0) pv = get_kuonen_pv(newQ, newL);
  if(pv <= 0) pv = get_liu_pv(newQ, newL);

  return pv;

}

// return 1-F(x) for chisq mixture
double get_davies_pv(double const& q, Ref<VectorXd> lambdas){

  int k = lambdas.size(), ifault = 0; // p & error
  double cdf, pv;
  ArrayXd nc = ArrayXd::Constant(k, 0); // ncp
  ArrayXi df = ArrayXi::Constant(k, 1); // df
  ArrayXd tr = ArrayXd::Constant(7, 0); // params for qf

  // use default lim/acc values from CompQuadForm R package and SKAT resp.
  try {
    cdf = qf(lambdas.data(), nc.data(), df.data(), k, 0, q, 1e4, 1e-6, tr.data(), &ifault); 
    pv = 1 - cdf;
  } catch (...){
    return -1;
  }
  //cerr << "Davies p=" << pv << "\n";

  if((ifault != 0) || (pv <= 0) || (pv > 1))
    return -1;

  return pv;
}

// return 1-F(x) for chisq mixture
double get_kuonen_pv(const double& q, const Ref<const VectorXd>& L){

  bool success = false;
  double pv, t_root = -1;
  MapcArXd lambdas (L.data(), L.size(), 1);
  //cerr << "q=" << q << "\ntop lambdas=" << L.tail(6) << "\n\n";

  // lambdas are sorted in increasing order (from eigen)
  double tmin = get_tmin_lambda(q, lambdas);
  double tmax = get_tmax_lambda(lambdas);
  //cerr << "(" << tmin << "," << tmax << ")\n\n";
  if(tmax < tmin) return -1;

  solve_kp(success, t_root, q, tmin, tmax, lambdas);
  if(!success) return -1;

  pv = get_spa_pv(t_root, q, lambdas);
  //cerr << "SPA p=" << pv << "\n";

  if((pv <= 0) || (pv > 1)) return -1;

  return pv;
}

double get_tmin_lambda(const double& q, const Ref<const ArrayXd>& lambdas){
  if(lambdas(0) < 0) // not applicable here since matrix is psd
    return 1 / (2 * lambdas(0));
  else if(q > lambdas.sum())
    return 0;
  else
    return -0.5 * lambdas.size() / q;
}

double get_tmax_lambda(const Ref<const ArrayXd>& lambdas){
  //return 1 / (2 * lambdas.tail(1)(0));
  return 0.5 - 1e-8; // lambdas are re-scaled so max=1
}

void solve_kp(bool& success, double& t_new,const double& q,const double& tmin,const double& tmax, const Ref<const ArrayXd>& lambdas){

  int niter_cur = 0, niter_max = 1e3;
  double min_x, max_x, t_old, f_old, f_new, hess, tol = 1e-8;

  min_x = tmin, max_x = tmax;
  t_old = min_x;
  // check sign switches
  if(!valid_bounds(f_old, min_x, tmax, q, lambdas)) {
    success = false;
    return;
  }

  while( niter_cur++ < niter_max ){

    hess = Kpp_lambda(t_old,lambdas);
    t_new = t_old - f_old / hess;
    f_new = Kp_lambda(t_new,lambdas) - q;

    //cerr << "#" << niter_cur << ": t=" << t_old << "->" << t_new << " f(t)=" << f_new << "; bounds = (" << min_x << "," << max_x << ")\n";
    if( fabs( f_new ) < tol ) break;

    // update bounds on root
    if( (t_new > min_x) && (t_new < max_x) ){
      if( f_new > 0) max_x = t_new;
      else min_x = t_new;
    } else { // bisection method if t_new went out of bounds and re-compute f_new
      t_new = ( min_x + max_x ) * 0.5;
      f_new = Kp_lambda(t_new,lambdas) - q;
      if(f_new <= 0) min_x = t_new; // reduce interval
      else max_x = t_new;
    }

    t_old = t_new;
    f_old = f_new;
  }

  // If didn't converge
  success = niter_cur <= niter_max;
  //cerr << "#iterations = " << niter_cur << "; f= " << f_new << endl;

}

bool valid_bounds (double& fmin, double const& tmin, double const& tmax, const double& q, const Ref<const ArrayXd>& lambdas){ 

  fmin = Kp_lambda(tmin,lambdas) - q;
  double fmax = Kp_lambda(tmax,lambdas) - q;

  return ((fmin<=0) && (fmax>=0));
}

double K_lambda (const double& t, const Ref<const ArrayXd>& lambdas){ 
  return -0.5 * (1 - 2 * t * lambdas).log().sum();
}

double Kp_lambda (const double& t, const Ref<const ArrayXd>& lambdas){ 
  return (lambdas / (1 - 2 * t * lambdas)).sum();
}

double Kpp_lambda (const double& t, const Ref<const ArrayXd>& lambdas){ 
  return (( 2 * lambdas.square()) / (1 - 2 * t * lambdas).square()).sum();
}

double get_spa_pv(const double& root,const double& q, const Ref<const ArrayXd>& lambdas){

  double u,w,r;
  normal nd(0,1);

  w = sgn(root) * sqrt( 2 * (q * root - K_lambda(root, lambdas)));
  u = root * sqrt(Kpp_lambda(root, lambdas));
  if( fabs(u) < 1e-4 ) return -1;

  r = w + log(u/w) / w;
  return cdf(complement(nd, r));

}

double get_liu_pv(double const& q, const Ref<const VectorXd>& lambdas){

  ArrayXd cvals(6);
  get_cvals(cvals, lambdas);
  //cerr << "cvals liu=" << cvals.matrix().transpose() << endl;
  
  double pv;
  double tstar = (q - cvals(0)) * cvals(1);
  double val = tstar * cvals(3) + cvals(2);

  if(val < 0) return -1;

  // 0 ncp gives strange behavior with non_central_chi_squared (returns -cdf instead of 1-cdf)
  if(cvals(5) == 0) pv = cdf(complement(chi_squared(cvals(4)), val));
  else  pv = cdf(complement(non_central_chi_squared(cvals(4), cvals(5)), val));
  //cerr << "pv liu=" << val << " -> " << pv << endl;

  if((pv <= 0) || (pv > 1)) return -1;
  return pv;
}

void get_skato_mom(double& mu, double& sd, ArrayXd& tau, int const& nnz, const Ref<const VectorXd>& lambdas, double const& ztz, double const& zmZ_nsq, const Ref<const MatrixXd>& ZtIUZ, const Ref<const MatrixXd>& ZtUZ, const Ref<const ArrayXd>& rho, bool const& debug){

  double s2_xi, s2Q;

  mu = lambdas.sum();
  s2_xi = 4 * (ZtIUZ.array() * ZtUZ.array()).sum();
  s2Q = 2 * lambdas.squaredNorm() + s2_xi;
  sd = sqrt( (s2Q - s2_xi) / s2Q );
  if(debug) cerr << "[muQ, scFac]=" << mu << " " << sd << "\n";
  tau = nnz * nnz * rho * ztz + (1-rho) * zmZ_nsq / ztz;
  if(debug) cerr << "tau=" << tau.matrix().transpose() << "\n";

}

void get_cvals(int const& irho, Ref<MatrixXd> cvals, const Ref<const VectorXd>& lambdas){

  double s1, s1_sq, s2, a, dlt;

  cvals(irho, 0) = lambdas.sum();
  cvals(irho, 1) = lambdas.squaredNorm();
  cvals(irho, 2) = lambdas.array().pow(3).sum();
  cvals(irho, 3) = lambdas.array().pow(4).sum();
  s1 = cvals(irho, 2) / cvals(irho, 1) / sqrt(cvals(irho, 1));
  s1_sq = s1 * s1;
  s2 = cvals(irho, 3) / (cvals(irho, 1) * cvals(irho, 1));
  if(s1_sq <= s2)
    cvals(irho, 4) = 1 / s2;
  else {
    a = 1 / (s1 - sqrt(s1_sq - s2));
    dlt = (s1 * a - 1) * a * a;
    cvals(irho, 4) = a * a - 2 * dlt;
  }

}

void get_cvals(Ref<ArrayXd> cvals, const Ref<const VectorXd>& lambdas){

  // cvals = [muQ, invsQ, muX, sX, df, ncp]
  double c1, c2, c3, c4, s1, s1_sq, s2, df, ncp, a;

  c1 = lambdas.sum();
  c2 = lambdas.squaredNorm();
  c3 = lambdas.array().pow(3).sum();
  c4 = lambdas.array().pow(4).sum();
  s1 = c3 / c2 / sqrt(c2);
  s1_sq = s1 * s1;
  s2 = c4 / (c2 * c2);
  if(s1_sq <= s2) {
    df = 1 / s2;
    a = sqrt(df);
    ncp = 0;
  } else {
    a = 1 / (s1 - sqrt(s1_sq - s2));
    ncp = (s1 * a - 1) * a * a;
    df = a * a - 2 * ncp;
  }

  cvals(0) = c1; //muQ
  cvals(1) = 1 / sqrt(2 * c2); //invsQ
  cvals(2) = df + ncp; // muX
  cvals(3) = sqrt(2) * a; // sX
  cvals(4) = df;
  cvals(5) = ncp;

}

void get_Qmin(int const& nrho, double& pmin, Ref<ArrayXd> Qmin, const Ref<const MatrixXd>& cvals){
  for(int j = 0; j < nrho; j++){
    chi_squared chisq( cvals(j, 4) );
    Qmin(j) = cvals(j, 0) + (quantile(complement(chisq, pmin)) - cvals(j, 4)) * sqrt(cvals(j, 1)/cvals(j, 4)) ;
  }
}

double SKATO_integral_fn(double* x){ // variables used beside x are global

  double val = ((skato_Qmin_rho - skato_tau * (*x)) / flipped_skato_rho).minCoeff();
  double S, dlt;
  chi_squared chisq1( 1 );

  // get first term in integral
  if( val > (skato_muQ * 1e4) ) S = 0; // value check from SKAT R package
  else {
    dlt = (val - skato_muQ) * skato_sd + skato_muQ;
    if(dlt <= 0) return 0; // cdf=0
    else{
      // davies/kuonen returns -1 if fails
      S = get_chisq_mix_pv(dlt, skato_lambdas);

      if(S <= 0) { // failed
        skato_state = 1; 
        return 0;
      }
    }
  }

  if(S >= 1) return 0; // cdf=0

  return (1-S) * pdf(chisq1, *x);

}

// for skato num int
void integrate(double f(double*), double& pv, bool const& debug){

  int neval, ierror, ilimit = 1000, last; 
  int lenw = 4 * ilimit;
  double lower = 0, upper = 40, epsabs = 1e-25, epsrel = 0.0001220703, result, abserr; // rel. eps = eps^.25 in R
  VectorXi iwork = VectorXi::Zero(ilimit);
  VectorXd work = VectorXd::Zero(lenw);
  skato_state = 0;

  dqags_(f, &lower, &upper, &epsabs, &epsrel, &result, &abserr, &neval, &ierror, &ilimit, &lenw, &last, iwork.data(), work.data());
  if(debug) {
    cerr << "Niter=" << neval << ";integral=" << result << "Abs.error=" << abserr << ";rel.error=" << epsrel << 
      ";fail="<< skato_state << "/" << ierror << "\n";
    for(int i = 1; i < 5; i++) {lower=0.2*i;cerr << "g(" << lower << ")=" << SKATO_integral_fn(&lower) << " ";}
  }

  if ((skato_state != 0) || (ierror != 0))  pv = -1; 
  else pv = 1 - result;

} 

void get_skato_pv(double &logp, double& chisq, double const& minp, int const& nrhos, double const& nl_dbl_dmin, bool const& debug){

  double a, p_bc = minp * nrhos;

  integrate(SKATO_integral_fn, a, debug);
  if(debug) cerr << "\nSKATO p=" << a << " (minP="<< minp <<"; Bonf=" << p_bc << ")\n";

  if( p_bc < a ) a = p_bc; // bonferroni corrected p
  else if( (a <= 0) && (p_bc <= 1) ) a = p_bc; // if integrate function failed

  if(a <= 0) {logp = -1; return;} // if pmin=0
  
  get_logp(a, logp, chisq, nl_dbl_dmin); 

}

// print sum_stats
void print_vc_sumstats(int const& snp_index, string const& test_string, string const& wgr_string, variant_block* block_info, vector<snp> const& snpinfo, struct in_files const& files, struct param const* params){

  int print_index;
  string header;
  std::map <std::string, MatrixXd>::iterator itr;
  if(!params->htp_out) header = print_sum_stats_head(snp_index, snpinfo);

  for (itr = block_info->sum_stats_vc.begin(); itr != block_info->sum_stats_vc.end(); ++itr) {
    // for each pheno
    for(int i = 0; i < params->n_pheno; i++) { // col 0 = chisq, col 1 = logp

      // make sure results for test are all on same line
      print_index = params->split_by_pheno ? i : 0;

      if(itr->second(i, 1) >= 0) {
        std::ostringstream buffer;

        if(params->htp_out) 
          buffer << print_sum_stats_head_htp(snp_index, files.pheno_names[i], test_string + wgr_string + "-" + itr->first, snpinfo, params) << print_sum_stats_htp(-1, -1, itr->second(i, 0), itr->second(i, 1), -1, -1, -1, block_info->genocounts, i, true, 1, params);
        else 
          buffer << (!params->split_by_pheno && (i>0) ? "" : header) << print_sum_stats(-1,-1,-1, -1, (params->split_by_pheno ? block_info->ns(i) : block_info->ns1), test_string + "-" + itr->first, -1, -1, itr->second(i, 0), itr->second(i, 1), true, 1, params, (i+1));

        block_info->sum_stats[print_index].append( buffer.str() );
      } else if(!params->split_by_pheno) // print NA sum stats
        block_info->sum_stats[print_index].append( print_na_sumstats(i, 1, header, test_string + "-" + itr->first, block_info, *params) );

    }
  }

}

