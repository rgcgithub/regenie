/*  

    This file is part of the regenie software package.

    Copyright (c) 2020-2023 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

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
double skato_fdavies = 0;
double skato_sdQ = 0;
double skato_dfQ = 0;
double skato_upper = 0;
int skato_state = 0;
// for LOVO with BTs
MatrixXd vc_Rvec_start;

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Functions for SKAT/SKAT-O
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void update_vc_gmat(SpMat& mat, ArrayXd& weights, ArrayXd& weights_acat, ArrayXb& ur_ind, int const& start, int const& bs, struct param const& params, const Ref<const ArrayXb>& in_analysis, Ref<MatrixXd> Gmat, vector<variant_block> &all_snps_info, Ref<MatrixXb> Jmat){

  beta_distribution<>  dist(params.skat_a1, params.skat_a2);

  /*if(params.mask_loo){ // update dimensions
    mat.resize(mat.rows(), bs + Jmat.cols());
    mat.setZero();
    weights = ArrayXd::Zero(bs + Jmat.cols(), 1);
    weights_acat = weights;
    if(params.debug) cerr << "Updating VC gmat...";
  }*/

#if defined(_OPENMP)
      setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < bs; ++i) {

    MapArXd Gvec (Gmat.col(i).data(), Gmat.rows(), 1);
    double maf;

    if(Jmat.row(start + i).any()){ // if variant is in at least one mask
      // check if ultra-rare (if so set to 0)
      if(ur_ind(start + i)){
        Jmat.row(start+i).array() = false; // ignore variant for all sets
        Gvec = 0; // don't store the variant
        continue;
      }

      // flip if af is above 0.5
      if( all_snps_info[i].af1 > 0.5 ) Gvec = (Gvec == -3).select(-3, 2 - Gvec);
      maf = min(all_snps_info[i].af1, 1 - all_snps_info[i].af1);

      // impute missing with mean
      Gvec = (Gvec == -3).select(2 * maf, Gvec);
      // mask individuals
      Gvec *= in_analysis.cast<double>();
      // store SKAT weight
      if(!params.vc_with_weights){
        weights(start + i) = pdf(dist, maf);
        weights_acat(start + i) = weights(start + i) * weights(start + i) * maf * (1-maf); // for acatv
      } else if(params.vc_multiply_weights){
        weights(start + i) *= pdf(dist, maf);
        weights_acat(start + i) = weights(start + i) * weights(start + i) * maf * (1-maf); // for acatv
      }
    } else Gvec = 0; // otherwise set the column to 0

  }
#if defined(_OPENMP)
      setNbThreads(params.threads);
#endif

  mat.middleCols(start, bs) = Gmat.sparseView();
}

// with lovo
void update_vc_gmat(SpMat& mat, ArrayXd& weights, ArrayXd& weights_acat, SpMat const& Gmat_sp, const Ref<const ArrayXb>& ur_ind, const Ref<const ArrayXb>& to_flip, const Ref<const ArrayXd>& mafs, const Ref<const ArrayXb>& in_analysis, struct param const& params){

  beta_distribution<>  dist(params.skat_a1, params.skat_a2);
  int bs = Gmat_sp.cols(), bsize = floor(1e9/8.0/params.n_samples), start = 0;
  int nchunks = ceil( bs * 1.0 / bsize );

  for (int j = 0; j < nchunks; ++j) {
    if( j == (nchunks-1) ) bsize = bs - j * bsize;
    MatrixXd Gmat = MatrixXd::Zero(Gmat_sp.rows(), bsize); // no mt with spmat

#if defined(_OPENMP)
    setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < bsize; ++i) {
      if(ur_ind(start + i)) continue; // check if ultra-rare

      MapArXd Gvec (Gmat.col(i).data(), Gmat.rows(), 1);
      Gvec = Gmat_sp.col(start + i);
      if(!(Gvec>0).any()) continue; // not used

      // flip if af is above 0.5
      if( to_flip(start + i) ) Gvec = (Gvec == -3).select(Gvec, 2 - Gvec);

      // impute missing with mean
      Gvec = (Gvec == -3).select(2 * mafs(start + i), Gvec);
      // mask individuals
      Gvec *= in_analysis.cast<double>();
      // store SKAT weight
      if(!params.vc_with_weights){
        weights(start + i) = pdf(dist, mafs(start + i));
        weights_acat(start + i) = weights(start + i) * weights(start + i) * mafs(start + i) * (1-mafs(start + i)); // for acatv
      } else if(params.vc_multiply_weights){
        double v_pdf = pdf(dist, mafs(start + i));
        weights(start + i) *= v_pdf;
        weights_acat(start + i) = weights(start + i) * weights(start + i) * mafs(start + i) * (1-mafs(start + i)); // for acatv
      }

    }
#if defined(_OPENMP)
    setNbThreads(params.threads);
#endif
    mat.middleCols(start, bsize) = Gmat.sparseView();
    start += bsize;
  }

}

bool get_custom_weights(string const& setname, Ref<ArrayXd> weights, Ref<ArrayXd> weights_acat, vector<snp>& snpinfo, vector<uint64> const& indices){

  // load custom user weights
  for(size_t i = 0; i < indices.size(); i++){
    if(!in_map(setname, snpinfo[ indices[i] ].set_weight)) // this shouldn't happen
      throw "no custom weight found for variant " + snpinfo[ indices[i] ].ID;
    weights(i) = snpinfo[ indices[i] ].set_weight[setname];
  }
  double sum_w = weights.sum(); // make weights sum to 1
  if(sum_w == 0) return false;

  weights /= sum_w;
  weights_acat = weights;
;

  return true;
}

bool get_custom_weights(string const& setname, Ref<ArrayXd> weights, vector<snp>& snpinfo, vector<uint64> const& offsets){

  // load custom user weights
  for(size_t i = 0; i < offsets.size(); i++){
    if(!in_map(setname, snpinfo[ offsets[i] ].set_weight)) // this shouldn't happen
      throw "no custom weight found for variant " + snpinfo[ offsets[i] ].ID;
    weights(i) = snpinfo[ offsets[i] ].set_weight[setname];
  }
  double sum_w = weights.sum(); // make weights sum to 1
  if(sum_w == 0) return false;

  weights /= sum_w;
  return true;
}

// with lovo
bool get_custom_weights(string const& setname, Ref<ArrayXd> weights, vector<snp>& snpinfo, const Ref<const ArrayXi>& indices, vector<uint64> const& offsets){

  // load custom user weights
  for(int i = 0; i < indices.size(); i++){
    if(!in_map(setname, snpinfo[ offsets[indices(i)] ].set_weight)) // this shouldn't happen
      throw "no custom weight found for variant " + snpinfo[ offsets[indices(i)] ].ID;
    weights(i) = snpinfo[ offsets[indices(i)] ].set_weight[setname];
  }
  double sum_w = weights.sum(); // make weights sum to 1
  if(sum_w == 0) return false;

  weights /= sum_w;
  return true;
}

void compute_vc_masks(SpMat& mat, Ref<ArrayXd> weights, Ref<ArrayXd> weights_acat, SpMat& vc_rare_mask, Ref<MatrixXb> vc_rare_non_miss, const Ref<const MatrixXd>& X, struct ests const& m_ests, struct f_ests const& fest, const Ref<const MatrixXd>& yres,  const Ref<const MatrixXd>& yraw, const Ref<const MatrixXb>& masked_indivs, MatrixXb& Jmat, vector<variant_block> &all_snps_info, const Ref<const ArrayXb>& in_analysis, struct param const& params, struct remeta_sumstat_writer& remeta_sumstats){

  prep_ultra_rare_mask(mat, weights, weights_acat, vc_rare_mask, vc_rare_non_miss, Jmat, in_analysis, params);

  //if(params.debug) check_sizes(mat, vc_rare_mask, Jmat);

  if(params.trait_mode==0)
    compute_vc_masks_qt(mat, weights, weights_acat, X, yres, Jmat, all_snps_info, params, remeta_sumstats);
  else if(params.trait_mode==1)
    compute_vc_masks_bt(mat, weights, weights_acat, X, m_ests, fest, yres, yraw, masked_indivs, Jmat, all_snps_info, params, remeta_sumstats);
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
    Jmat(bs+iset, iset) = (gv>0).any();
    //cerr << "#" << iset << " " << Jmat.col(iset).any() << "/"<< (gv>0).count() << "\n";
    if(!Jmat(bs+iset, iset)) continue;
    rare_mask_non_miss.col(iset).array() = rare_mask_non_miss.col(iset).array() && in_analysis; 

    // compute mean
    mean = gv.sum() / rare_mask_non_miss.col(iset).count();
    maf = min(mean/2, 1 - mean/2);

    if(params.vc_with_weights && !params.vc_multiply_weights){ // weights already incorporated when taking max of weighted geno
      weights(bs+iset) = weights_acat(bs+iset) = 1;
    } else { // use default SKAT/ACAT weight
      weights(bs+iset) = pdf(dist, maf);
      weights_acat(bs+iset) = weights(bs+iset) * weights(bs+iset) * maf * (1-maf);
    }

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
void compute_vc_masks_qt(SpMat& mat, const Ref<const ArrayXd>& weights, const Ref<const ArrayXd>& weights_acat, const Ref<const MatrixXd>& X, const Ref<const MatrixXd>& yres, const Ref<const MatrixXb>& Jmat, vector<variant_block> &all_snps_info, struct param const& params, struct remeta_sumstat_writer& remeta_sumstats){

  if(params.skato_rho.size() == 1)
    compute_vc_masks_qt_fixed_rho(mat, weights, weights_acat, X, yres, Jmat, all_snps_info, params.skato_rho(0), params.skat_tol, params.nl_dbl_dmin, params.vc_test, params.debug, params, remeta_sumstats);
  else
    compute_vc_masks_qt(mat, weights, weights_acat, X, yres, Jmat, all_snps_info, params.skato_rho, params.skat_tol, params.nl_dbl_dmin, params.vc_test, params.debug, params, remeta_sumstats);

}

// for a given rho value
void compute_vc_masks_qt_fixed_rho(SpMat& mat, const Ref<const ArrayXd>& weights, const Ref<const ArrayXd>& weights_acat, const Ref<const MatrixXd>& X, const Ref<const MatrixXd>& yres, const Ref<const MatrixXb>& Jmat, vector<variant_block> &all_snps_info, double const& rho, double const& skat_lambda_tol, double const& nl_dbl_dmin, uint const& vc_test, bool const& debug, struct param const& params, struct remeta_sumstat_writer& remeta_sumstats){

  bool with_acatv = CHECK_BIT(vc_test,0);
  bool with_skat = (vc_test>>1)&15;
  int jcol, n_pheno = yres.cols(), nnz;
  double c1 = sqrt(1 - rho);
  ArrayXd D;
  VectorXd lambdas;
  MatrixXd Qs, Qb, Svals, Kmat, sum_stats, pvals;

  ArrayXi snp_indices = get_true_indices(Jmat.rowwise().any());
  int bs = snp_indices.size(); // subset to snps included in at least 1 skat mask
  if( !(weights(snp_indices) > 0).any() ) return;

  // slice sparse matrix (cannot use indexing)
  SpMat Jstar (Jmat.rows(), bs); // Mall x M
  Jstar.reserve(bs);
  MatrixXd weights_ordered;
  if(params.remeta_save_ld) {
    weights_ordered.resize(bs, 1);
  }
  for(int i = 0; i < bs; i++) {
    Jstar.insert(snp_indices(i), i) = weights(snp_indices(i));
    if(params.remeta_save_ld) {
        weights_ordered(i) = weights(snp_indices(i));
    }
  }
  SpMat mat2 = mat * Jstar; // mat should be pretty sparse since major-ref
  mat.setZero(); mat.resize(0,0); mat.data().squeeze(); // not needed anymore
  Jstar.setZero(); Jstar.resize(0,0); Jstar.data().squeeze(); // not needed anymore

  // get score stats & kernel matrices
  Svals.resize(n_pheno, bs); // PxM
  Kmat.resize(bs, bs); // MxM
  compute_vc_mats_qt(Svals, Kmat, X, yres, mat2);
  mat2.setZero(); mat2.resize(0,0); mat2.data().squeeze(); // not needed anymore

#ifdef WITH_HTSLIB
  if(params.remeta_save_ld && remeta_sumstats.skat_snplist->size() > 0) {
    MatrixXd weight_inv = weights_ordered.array()
                                  .inverse()
                                  .matrix()
                                  .asDiagonal();
    MatrixXd unweighted_Kmat = weight_inv * Kmat * weight_inv;
    for(int i = 0; i < n_pheno; ++i) {
      if(remeta_sumstats.sparsity_threshold > 0) {
        remeta_sumstats.skat_matrix_writers[i].write_matrix_sparse(
          unweighted_Kmat,
          *remeta_sumstats.gene_name,
          *remeta_sumstats.skat_snplist,
          remeta_sumstats.sparsity_threshold
        );
      } else {
        remeta_sumstats.skat_matrix_writers[i].write_matrix_dense(
          unweighted_Kmat,
          *remeta_sumstats.gene_name,
          *remeta_sumstats.skat_snplist
        );      
      }
    }
  }
#endif

  // SKAT for all masks & traits
  compute_skat_q(Qs, Qb, Svals, Kmat, pvals, weights(snp_indices) != 0, Jmat(snp_indices, all), with_acatv, debug);

  // for now don't parallelize this as it causes issues with qfc lib
  // but should be ok since dimensions don't depend on N
  for(size_t imask = 0; imask < all_snps_info.size(); imask++){

    variant_block* block_info = &(all_snps_info[imask]);
    if(debug) cerr << "Mask : " << block_info->mask_name << "\n";
    if(block_info->sum_stats_vc.size()>0) block_info->sum_stats_vc.clear();
    if(block_info->skip_for_vc) continue;
    sum_stats = MatrixXd::Constant(n_pheno, 2, -1); // chisq & logp

    // get index of mask in Jmat
    jcol = block_info->col_jmat_skat;
    if(jcol < 0) continue; // this should not happen though
    MapcArXb Jvec (Jmat.col(jcol).data(), Jmat.rows(), 1);
    nnz = Jvec.count();
    if(debug) cerr << "#sites in mask=" << nnz << "\n";
    if(nnz == 0) continue;

    // subset to variants kept in mask
    ArrayXi m_indices = get_true_indices(Jvec(snp_indices)); // across markers kept in skat tests
    ArrayXi mall_indices = snp_indices(m_indices); // across all markers in set
    if(debug) cerr <<"W(skat):\n" << weights(mall_indices).head(min(20,nnz)).matrix().transpose() << "\n";

    // ACAT-V 
    if(with_acatv && (weights_acat(mall_indices) > 0).any()){
      for(int ph = 0; ph < n_pheno; ph++)
        get_acatv_pv( ph, pvals(m_indices, ph), weights_acat(mall_indices), sum_stats(ph, 1), sum_stats(ph, 0), nl_dbl_dmin, debug); 
      block_info->sum_stats_vc["ACATV"] = sum_stats;
      sum_stats.array() = -1; // reset
    }
    if(!with_skat) continue;

    // get eigen values of Rsqrt*V*Rsqrt
    //if(debug) cerr << "Kmat:\n" << Kmat(m_indices, m_indices) << "\nrho_Kmat:\n" << get_RsKRs(Kmat(m_indices, m_indices), rho, c1) << "\n";
    get_lambdas(lambdas, get_RsKRs(Kmat(m_indices, m_indices), rho, c1), skat_lambda_tol);
    if(lambdas.size() == 0) continue;
    if(debug) cerr << "L:" << lambdas.head(min(150, (int) lambdas.size())).transpose() << "\n";

    // compute test statistic & p-value
    for(int ph = 0; ph < n_pheno; ph++)
      compute_fixed_skato_p(sum_stats(ph, 1), sum_stats(ph, 0), Qs(ph, jcol), Qb(ph, jcol), rho, lambdas, nl_dbl_dmin, debug);

    if( (sum_stats.col(1).array() >= 0).any() ){
      string test_name = (rho > 0 ? "SKAT-RHO" : "SKAT");
      block_info->sum_stats_vc[test_name] = sum_stats;
    }

  }

}

void compute_vc_masks_qt(SpMat& mat, const Ref<const ArrayXd>& weights, const Ref<const ArrayXd>& weights_acat, const Ref<const MatrixXd>& X, const Ref<const MatrixXd>& yres, const Ref<const MatrixXb>& Jmat, vector<variant_block> &all_snps_info, const Ref<const ArrayXd>& rho_vec, double const& skat_lambda_tol, double const& nl_dbl_dmin, uint const& vc_test, bool const& debug, struct param const& params, struct remeta_sumstat_writer& remeta_sumstats){

  bool with_acatv = CHECK_BIT(vc_test,0);
  bool with_omnibus = (vc_test>>2)&7; // any of the omnibus tests
  bool with_skato_int = CHECK_BIT(vc_test,2);
  bool with_skato_acat = CHECK_BIT(vc_test,3);
  bool with_acato = CHECK_BIT(vc_test,4);
  int jcol, n_pheno = yres.cols(), nnz, nrho = rho_vec.size();
  double minp, gamma1, gamma2, gamma3, tmpv, log10_nl_dbl_dmin = -log10(nl_dbl_dmin);
  ArrayXd D, p_acato, flip_rho_sqrt;
  VectorXd lambdas;
  MatrixXd Qs, Qb, Qopt, Svals, Kmat, cvals, sum_stats, pvals, r_outer_sum;
  MatrixXd pvs_skato, chisq_skato, pvs_skato_acat, chisq_skato_acat, pvs_acato, chisq_acato;

  cvals.resize(nrho, 5);
  skato_Qmin_rho.resize(nrho, 1);
  if(with_acato) p_acato.resize(nrho+1);
  flipped_skato_rho = 1 - rho_vec;
  flip_rho_sqrt = flipped_skato_rho.sqrt();

  ArrayXi snp_indices = get_true_indices(Jmat.rowwise().any());
  int bs = snp_indices.size(); // subset to snps included in at least 1 skat mask
  if( !(weights(snp_indices) > 0).any() ) return;

  // slice sparse matrix (cannot use indexing)
  SpMat Jstar (Jmat.rows(), bs); // Mall x M
  Jstar.reserve(bs);
  MatrixXd weights_ordered;
  if(params.remeta_save_ld) {
    weights_ordered.resize(bs, 1);
  }
  for(int i = 0; i < bs; i++) {
    Jstar.insert(snp_indices(i), i) = weights(snp_indices(i));
    if(params.remeta_save_ld) {
        weights_ordered(i) = weights(snp_indices(i));
    }
  }
    
  SpMat mat2 = mat * Jstar; // mat should be pretty sparse since major-ref
  mat.setZero(); mat.resize(0,0); mat.data().squeeze(); // not needed anymore
  Jstar.setZero(); Jstar.resize(0,0); Jstar.data().squeeze(); // not needed anymore

  // get score stats & kernel matrices
  Svals.resize(n_pheno, bs); // PxM
  Kmat.resize(bs, bs); // MxM
  compute_vc_mats_qt(Svals, Kmat, X, yres, mat2);
  mat2.setZero(); mat2.resize(0,0); mat2.data().squeeze(); // not needed anymore

#ifdef WITH_HTSLIB
  if(params.remeta_save_ld && remeta_sumstats.skat_snplist->size() > 0) {
    MatrixXd weight_inv = weights_ordered.array()
                                  .inverse()
                                  .matrix()
                                  .asDiagonal();
    MatrixXd unweighted_Kmat = weight_inv * Kmat * weight_inv;
    for(int i = 0; i < n_pheno; ++i) {
      if(remeta_sumstats.sparsity_threshold > 0) {
        remeta_sumstats.skat_matrix_writers[i].write_matrix_sparse(
          unweighted_Kmat,
          *remeta_sumstats.gene_name,
          *remeta_sumstats.skat_snplist,
          remeta_sumstats.sparsity_threshold
        );
      } else {
        remeta_sumstats.skat_matrix_writers[i].write_matrix_dense(
          unweighted_Kmat,
          *remeta_sumstats.gene_name,
          *remeta_sumstats.skat_snplist
        );      
      }
    }
  }
#endif

  // SKAT for all masks & traits
  compute_skat_q(Qs, Qb, Svals, Kmat, pvals, weights(snp_indices) != 0, Jmat(snp_indices, all), with_acatv, debug);

  // for now don't parallelize this as it causes issues with qfc lib
  // but should be ok since dimensions don't depend on N
  for(size_t imask = 0; imask < all_snps_info.size(); imask++){

    variant_block* block_info = &(all_snps_info[imask]);
    if(debug) cerr << "Mask : " << block_info->mask_name << "\n";
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
    if(debug) cerr << "#sites in mask=" << nnz << "\n";
    if(nnz == 0) continue;

    // subset to variants kept in mask
    ArrayXi m_indices = get_true_indices(Jvec(snp_indices)); // across markers kept in skat tests
    ArrayXi mall_indices = snp_indices(m_indices); // across all markers in set

    // ACAT-V 
    if(with_acatv && (weights_acat(mall_indices) > 0).any()){
      for(int ph = 0; ph < n_pheno; ph++)
        get_acatv_pv( ph, pvals(m_indices, ph), weights_acat(mall_indices), sum_stats(ph, 1), sum_stats(ph, 0), nl_dbl_dmin, debug); 
      block_info->sum_stats_vc["ACATV"] = sum_stats;
      sum_stats.array() = -1; // reset
    }
    if(!with_omnibus) continue;

    // get eigen values of Zt(I-U)Z
    if(!get_ztz_evals(Kmat(m_indices, m_indices), r_outer_sum, gamma1, gamma2, gamma3, skat_lambda_tol, debug)) continue;

    if(nnz > 1){
      get_skato_mom(skato_muQ, skato_fdavies, skato_sdQ, skato_dfQ, skato_tau, skato_lambdas, gamma1, gamma2, gamma3, rho_vec, debug);
      if(skato_sdQ < 0) continue; // failed
    }
    Qopt = Qs.col(jcol) * flipped_skato_rho.matrix().transpose() + Qb.col(jcol) * rho_vec.matrix().transpose(); // P x Nrho
    if(debug) cerr << "Q:" << Qopt.row(0) << "\n";

    for(int j = 0; j < nrho; j++){

      // get eigen values of Rsqrt*ZtZ*Rsqrt
      get_lambdas(lambdas, get_RsKRs(Kmat(m_indices, m_indices), r_outer_sum, gamma1, rho_vec(j), flip_rho_sqrt(j)), skat_lambda_tol);
      //if(debug) cerr << "rho:" << rho_vec(j) << "-> L:" << lambdas.transpose() << "\n";
      if(lambdas.size() == 0) {
        if(debug) cerr << "all eigen values are 0 for rho = " << rho_vec(j) << "\n";
        break; // SKAT & SKAT-O failed
      }
      // needed for skato (M>1)
      if(nnz > 1)  get_cvals(j, cvals, lambdas);

      for(int ph = 0; ph < n_pheno; ph++) 
        compute_fixed_skato_p(pvs_skato(ph, j), chisq_skato(ph, j), Qopt(ph, j), rho_vec(j), lambdas, nl_dbl_dmin);

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
    if((pvs_skato.array() < 0).rowwise().any().all()) {
      if(debug) cerr << "Some SKATO pvalues failed to compute for all the phenotypes (p = " << pvs_skato << " )\n";
      continue; // go to next mask set if no phenotype with all pvals defined
    }

    // Get minimum of p-values and corresponding chisq quantile for each rho
    for(int ph = 0; ph < n_pheno; ph++) {
      if( (pvs_skato.row(ph).array() < 0).any() ) {
        if(debug) cerr << "One of the SKATO pvalues failed for phenotype (p = " << pvs_skato.row(ph) << " )\n";
        continue;
      }

      if(with_skato_acat){
        if(debug) cerr << "skato-acat logp=" << pvs_skato.row(ph) <<"\n";
        pvs_skato_acat(ph,0) = get_acat(pvs_skato.row(ph).array());
        get_chisq_stat_pv(tmpv, chisq_skato_acat(ph, 0), pvs_skato_acat(ph,0), nl_dbl_dmin, log10_nl_dbl_dmin);
      } 
      if(with_acato){ // include acatv pvalue
        p_acato(0) = block_info->sum_stats_vc["ACATV"](ph, 1);
        p_acato.tail(nrho) = pvs_skato.row(ph).transpose().array();
        if(debug) cerr << "acato logp=" << p_acato.matrix().transpose() <<"\n";
        pvs_acato(ph,0) = get_acat(p_acato);
        get_chisq_stat_pv(tmpv, chisq_acato(ph, 0), pvs_acato(ph,0), nl_dbl_dmin, log10_nl_dbl_dmin);
      }
      if(with_skato_int){
        minp = max(nl_dbl_dmin, pow(10, -pvs_skato.row(ph).maxCoeff())); // prevent underflow
        get_Qmin(nrho, minp, skato_Qmin_rho, cvals);
        if(debug) cerr << "Qmin=" << skato_Qmin_rho.matrix().transpose() << "\nminP=" << minp <<"; logp=" << pvs_skato.row(ph) <<"\n";
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

void compute_vc_mats_qt(Ref<MatrixXd> Svals, Ref<MatrixXd> Kmat, const Ref<const MatrixXd>& X, const Ref<const MatrixXd>& yres, Ref<SpMat> GW){

  // project covariates
  MatrixXd WGtX = GW.transpose() * X; // MxK

  // get test statistics (PxM matrix)
  // need to use Gresid (must have yres centered)
  Svals = yres.transpose() * GW - (yres.transpose() * X) * WGtX.transpose(); // 2nd term: PxK * KxM

  // get kernel matrix (MxM matrix)
  Kmat = -WGtX * WGtX.transpose();
  Kmat += GW.transpose() * GW;

}

void compute_skat_q(MatrixXd& Qs, MatrixXd& Qb, Ref<MatrixXd> Svals, const Ref<const MatrixXd>& Kmat, MatrixXd& pvals, const Ref<const ArrayXb>& mask_w, const Ref<const MatrixXb>& Jmat, bool const& w_acatv, bool const& debug){

  Qs = Svals.array().square().matrix(); // PxM
  // if using acat-v, get single variant p-values
  if(w_acatv) {
    pvals.resize(Svals.cols(), Svals.rows()); // MxP
    get_single_pvs(pvals, (Qs * mask_w.select(1/Kmat.diagonal().array(),1).matrix().asDiagonal()).transpose()); 
  }
  Qs *= Jmat.cast<double>(); // P x Km

  // burden
  Qb = (Svals * Jmat.cast<double>()).array().square().matrix();
  if(debug) cerr << "Q_SKAT for all masks:\n" << Qs << "\nQ_BURDEN for all masks:\n" << Qb << "\n";

}

void get_acatv_pv(int const& ph, const Ref<const MatrixXd>& pvals, const Ref<const ArrayXd>& weights, double& logp, double& chisq, double const& nl_dbl_dmin, bool const& debug){

  if(debug && (ph==0)) {
    int bs = pvals.rows(), nmax = min(150, bs); // avoid over-printing 
    cerr << "SV log10p:\n" << pvals.col(0).transpose().array().head(nmax) << 
    "\nWsq:\n" << weights.head(nmax).matrix().transpose() << "\n\n";
  }

  double tmpv, log10_nl_dbl_dmin = -log10(nl_dbl_dmin);
  logp = get_acat(pvals.array(), weights);
  get_chisq_stat_pv(tmpv, chisq, logp, nl_dbl_dmin, log10_nl_dbl_dmin);

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

void compute_vc_masks_bt(SpMat& mat, const Ref<const ArrayXd>& weights, const Ref<const ArrayXd>& weights_acat, const Ref<const MatrixXd>& X, struct ests const& m_ests, struct f_ests const& fest, const Ref<const MatrixXd>& yres, const Ref<const MatrixXd>& yraw, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXb>& Jmat, vector<variant_block> &all_snps_info, struct param const& params, struct remeta_sumstat_writer& remeta_sumstats){

  if(params.skato_rho.size() == 1)
    compute_vc_masks_bt_fixed_rho(mat, weights, weights_acat, X, m_ests, fest, yres, yraw, masked_indivs, Jmat, all_snps_info, params.skato_rho(0), params.skat_tol, params.nl_dbl_dmin, params.firth || params.use_SPA, params.vc_test, params.debug, params, remeta_sumstats);
  else 
    compute_vc_masks_bt(mat, weights, weights_acat, X, m_ests, fest, yres, yraw, masked_indivs, Jmat, all_snps_info, params.skato_rho, params.skat_tol, params.nl_dbl_dmin, params.firth || params.use_SPA, params.vc_test, params.debug, params, remeta_sumstats);

}

void compute_vc_masks_bt_fixed_rho(SpMat& mat, const Ref<const ArrayXd>& weights, const Ref<const ArrayXd>& weights_acat, const Ref<const MatrixXd>& X, struct ests const& m_ests, struct f_ests const& fest, const Ref<const MatrixXd>& yres, const Ref<const MatrixXd>& yraw, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXb>& Jmat, vector<variant_block> &all_snps_info, double const& rho, double const& skat_lambda_tol, double const& nl_dbl_dmin, bool const& apply_correction, uint const& vc_test, bool const& debug, struct param const& params, struct remeta_sumstat_writer& remeta_sumstats){

  bool with_acatv = CHECK_BIT(vc_test,0);
  bool with_skat = (vc_test>>1)&15;
  int jcol, n_pheno = yres.cols();
  double c1 = sqrt(1 - rho);
  VectorXd lambdas, Qs, Qb;
  ArrayXb masked_sites;
  ArrayXd Svals, Rvec_sqrt, pvals;
  MatrixXd Kmat, GtWX, sum_stats;
  SpMat GWs;

  MatrixXd pvs_m = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);
  MatrixXd chisq_m = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);
  MatrixXd pvs_m_a = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);
  MatrixXd chisq_m_a = MatrixXd::Constant( n_pheno, all_snps_info.size(), -1);
  sum_stats = MatrixXd::Constant(n_pheno, 2, -1); // chisq & logp

  ArrayXi snp_indices = get_true_indices(Jmat.rowwise().any());
  int bs = snp_indices.size(); // subset to snps included in at least 1 skat mask

  // slice sparse matrix (cannot use indexing)
  SpMat Jstar (Jmat.rows(), bs); // M x Mall
  Jstar.reserve(bs);
  MatrixXd weights_ordered;
  if(params.remeta_save_ld) {
    weights_ordered.resize(bs, 1);
  }
  for(int i = 0; i < bs; i++) {
    Jstar.insert(snp_indices(i), i) = weights(snp_indices(i));
    if(params.remeta_save_ld) {
        weights_ordered(i) = weights(snp_indices(i));
    }
  }
  SpMat mat2 = mat * Jstar; // mat should be pretty sparse since major-ref
  mat.setZero(); mat.resize(0,0); mat.data().squeeze(); // not needed anymore
  Jstar.setZero(); Jstar.resize(0,0); Jstar.data().squeeze(); // not needed anymore
  Svals.resize(bs, 1); // Mx1
  Kmat.resize(bs, bs); // MxM

  for(int ph = 0; ph < n_pheno; ph++) { 
    if( !params.pheno_pass(ph) ) continue;

    MapcArXd Y (yraw.col(ph).data(), yraw.rows());
    MapcArXb mask (masked_indivs.col(ph).data(), yraw.rows());
    MapcArXd Wsqrt (m_ests.Gamma_sqrt.col(ph).data(), yraw.rows());
    MapcMatXd XWsqrt (m_ests.X_Gamma[ph].data(), yraw.rows(), m_ests.X_Gamma[ph].cols());
    MapcArXd phat (m_ests.Y_hat_p.col(ph).data(), yraw.rows());

    // get score stats & kernel matrices
    compute_vc_mats_bt(Svals, Kmat, XWsqrt, Wsqrt * mask.cast<double>(), yres.col(ph), mat2, GWs, GtWX);

    // apply firth/spa corrections (set R=0 if failed)
    masked_sites = (weights(snp_indices) > 0);
    Rvec_sqrt = masked_sites.cast<double>();
    if(apply_correction)
      correct_vcov(ph, snp_indices, weights(snp_indices), masked_sites, Rvec_sqrt, Svals, Kmat, mat2, GtWX, XWsqrt, GWs, Wsqrt, phat, Y, mask, fest, params);

  #ifdef WITH_HTSLIB
    if(params.remeta_save_ld && remeta_sumstats.skat_snplist->size() > 0) {
      MatrixXd weight_inv = weights_ordered.array()
                                    .inverse()
                                    .matrix()
                                    .asDiagonal();
      MatrixXd unweighted_Kmat = weight_inv * Kmat * weight_inv;
      if(remeta_sumstats.sparsity_threshold > 0) {
        remeta_sumstats.skat_matrix_writers[ph].write_matrix_sparse(
          unweighted_Kmat,
          *remeta_sumstats.gene_name,
          *remeta_sumstats.skat_snplist,
          remeta_sumstats.sparsity_threshold
        );
      } else {
        remeta_sumstats.skat_matrix_writers[ph].write_matrix_dense(
          unweighted_Kmat,
          *remeta_sumstats.gene_name,
          *remeta_sumstats.skat_snplist
        );
      }
    }
  #endif

    if(with_acatv) {
      pvals.resize(Svals.size()); // Mx1
      get_single_pvs_bt(pvals, masked_sites.select(Svals.square() / Kmat.diagonal().array(),1)); 
    }

    // SKAT for all masks (Kmx1)
    compute_skat_q(Qs, Qb, masked_sites.select(Svals, 0), Kmat, Jmat(snp_indices, all), debug);

    for(size_t imask = 0; imask < all_snps_info.size(); imask++){

      variant_block* block_info = &(all_snps_info[imask]);
      if(debug) cerr << "Mask : " << block_info->mask_name << "\n";
      if(block_info->skip_for_vc) continue;

      // get index of mask in Jmat
      jcol = block_info->col_jmat_skat;
      if(jcol < 0) continue; // this should not happen though
      MapcArXb Jvec (Jmat.col(jcol).data(), Jmat.rows(), 1);
      int npass = (Jvec(snp_indices) && masked_sites).count();
      if(debug) cerr << "#sites in mask=" << npass << "\n";
      if(npass == 0) continue;

      // subset to variants kept in mask
      ArrayXi m_indices = get_true_indices(Jvec(snp_indices) && masked_sites); // across markers kept in skat tests
      ArrayXi mall_indices = snp_indices(m_indices); // across all markers in set

      // ACAT-V 
      if(with_acatv)
        get_acatv_pv( ph, pvals(m_indices).matrix(), weights_acat(mall_indices), pvs_m_a(ph, imask), chisq_m_a(ph, imask), nl_dbl_dmin, debug); 

      if(!with_skat) continue;

      // correct using burden test
      double rfrac = 1;
      if(apply_correction && !params.skip_cf_burden && (npass > 1)) {// no need if M=1

        // to slice sparse matrix (cannot use indexing)
        SpMat Jtmp (bs, npass); // Mall x Mpass
        Jtmp.reserve(npass);
        for(int i = 0; i < npass; i++)
          Jtmp.insert(m_indices(i), i) = 1;

        if(!correct_vcov_burden(ph, rfrac, Qb(jcol), Kmat(m_indices, m_indices).sum(), GtWX(all, m_indices), XWsqrt, GWs * Jtmp, Wsqrt, phat, Y, mask, fest.cov_blup_offset, params))
          continue; // failed to correct with burden mask
      }

    // get eigen values of Rsqrt*V*Rsqrt
      get_lambdas(lambdas, get_RsKRs(rfrac * Kmat(m_indices, m_indices), rho, c1), skat_lambda_tol);
      if(lambdas.size() == 0) continue;
      //cerr << lambdas << "\n\n";

      // compute SKAT
      compute_fixed_skato_p(pvs_m(ph, imask), chisq_m(ph, imask), Qs(jcol), Qb(jcol), rho, lambdas, nl_dbl_dmin, debug);

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

}

void compute_vc_masks_bt(SpMat& mat, const Ref<const ArrayXd>& weights, const Ref<const ArrayXd>& weights_acat, const Ref<const MatrixXd>& X, struct ests const& m_ests, struct f_ests const& fest, const Ref<const MatrixXd>& yres, const Ref<const MatrixXd>& yraw, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXb>& Jmat, vector<variant_block> &all_snps_info, const Ref<const ArrayXd>& rho_vec, double const& skat_lambda_tol, double const& nl_dbl_dmin, bool const& apply_correction, uint const& vc_test, bool const& debug, struct param const& params, struct remeta_sumstat_writer& remeta_sumstats){

  bool with_acatv = CHECK_BIT(vc_test,0);
  bool with_omnibus = (vc_test>>2)&7;
  bool with_skato_int = CHECK_BIT(vc_test,2);
  bool with_skato_acat = CHECK_BIT(vc_test,3);
  bool with_acato = CHECK_BIT(vc_test,4);
  int jcol, n_pheno = yres.cols(), nrho = rho_vec.size();
  double minp, gamma1, gamma2, gamma3, tmpv, log10_nl_dbl_dmin = -log10(nl_dbl_dmin);
  VectorXd lambdas, Qs, Qb, Qopt;
  ArrayXb masked_sites;
  ArrayXd Svals, Rvec_sqrt, pvs_skato, chisq_skato, pvals, p_acato, flip_rho_sqrt;
  MatrixXd Kmat, GtWX, cvals, sum_stats, r_outer_sum;
  SpMat GWs;

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
  flip_rho_sqrt = flipped_skato_rho.sqrt();

  ArrayXi snp_indices = get_true_indices(Jmat.rowwise().any());
  int bs = snp_indices.size(); // subset to snps included in at least 1 skat mask

  // slice sparse matrix (cannot use indexing)
  SpMat Jstar (Jmat.rows(), bs); // M x Mall
  Jstar.reserve(bs);
  MatrixXd weights_ordered;
  if(params.remeta_save_ld) {
    weights_ordered.resize(bs, 1);
  }
  for(int i = 0; i < bs; i++) {
    Jstar.insert(snp_indices(i), i) = weights(snp_indices(i));
    if(params.remeta_save_ld) {
        weights_ordered(i) = weights(snp_indices(i));
    }
  }
  SpMat mat2 = mat * Jstar; // mat should be pretty sparse since major-ref
  mat.setZero(); mat.resize(0,0); mat.data().squeeze(); // not needed anymore
  Jstar.setZero(); Jstar.resize(0,0); Jstar.data().squeeze(); // not needed anymore
  Svals.resize(bs, 1); // Mx1
  Kmat.resize(bs, bs); // MxM

  for(int ph = 0; ph < n_pheno; ph++) { 
    if( !params.pheno_pass(ph) ) continue;

    MapcArXd Y (yraw.col(ph).data(), yraw.rows());
    MapcArXb mask (masked_indivs.col(ph).data(), yraw.rows());
    MapcArXd Wsqrt (m_ests.Gamma_sqrt.col(ph).data(), yraw.rows());
    MapcMatXd XWsqrt (m_ests.X_Gamma[ph].data(), yraw.rows(), m_ests.X_Gamma[ph].cols());
    MapcArXd phat (m_ests.Y_hat_p.col(ph).data(), yraw.rows());

    // get score stats & kernel matrices
    compute_vc_mats_bt(Svals, Kmat, XWsqrt, Wsqrt * mask.cast<double>(), yres.col(ph), mat2, GWs, GtWX);

    // apply firth/spa corrections (set R=0 if failed)
    masked_sites = (weights(snp_indices) > 0);
    Rvec_sqrt = masked_sites.cast<double>();
    if(apply_correction) {
      correct_vcov(ph, snp_indices, weights(snp_indices), masked_sites, Rvec_sqrt, Svals, Kmat, mat2, GtWX, XWsqrt, GWs, Wsqrt, phat, Y, mask, fest, params);
    }

    #ifdef WITH_HTSLIB
      if(params.remeta_save_ld && remeta_sumstats.skat_snplist->size() > 0) {
        MatrixXd weight_inv = weights_ordered.array()
                                      .inverse()
                                      .matrix()
                                      .asDiagonal();
        MatrixXd unweighted_Kmat = weight_inv * Kmat * weight_inv;

        if(remeta_sumstats.sparsity_threshold > 0) {
          remeta_sumstats.skat_matrix_writers[ph].write_matrix_sparse(
            unweighted_Kmat,
            *remeta_sumstats.gene_name,
            *remeta_sumstats.skat_snplist,
            remeta_sumstats.sparsity_threshold
          );
        } else {
          remeta_sumstats.skat_matrix_writers[ph].write_matrix_dense(
            unweighted_Kmat,
            *remeta_sumstats.gene_name,
            *remeta_sumstats.skat_snplist
          );
        }
      }
    #endif

    if(with_acatv) {
      pvals.resize(Svals.size()); // Mx1
      get_single_pvs_bt(pvals, masked_sites.select(Svals.square() / Kmat.diagonal().array(),1)); 
    }

    // SKAT for all masks (Kmx1)
    compute_skat_q(Qs, Qb, masked_sites.select(Svals, 0), Kmat, Jmat(snp_indices, all), debug);

    for(size_t imask = 0; imask < all_snps_info.size(); imask++){

      variant_block* block_info = &(all_snps_info[imask]);
      if(debug) cerr << "Mask : " << block_info->mask_name << "\n";
      if(block_info->skip_for_vc) continue;

      // get index of mask in Jmat
      jcol = block_info->col_jmat_skat;
      if(jcol < 0) continue; // this should not happen though
      MapcArXb Jvec (Jmat.col(jcol).data(), Jmat.rows(), 1);
      int npass = (Jvec(snp_indices) && masked_sites).count();
      if(debug) cerr << "#sites in mask=" << npass << "\n";
      if(npass == 0) continue;

      // subset to variants kept in mask
      ArrayXi m_indices = get_true_indices(Jvec(snp_indices) && masked_sites); // across markers kept in skat tests
      ArrayXi mall_indices = snp_indices(m_indices); // across all markers in set

      // ACAT-V 
      if(with_acatv)
        get_acatv_pv( ph, pvals(m_indices).matrix(), weights_acat(mall_indices), pvs_m_a(ph, imask), chisq_m_a(ph, imask), nl_dbl_dmin, debug); 
      if(!with_omnibus) continue;

      // correct using burden test
      double rfrac = 1;
      if(apply_correction && !params.skip_cf_burden && (npass > 1)) {// no need if M=1

        // to slice sparse matrix (cannot use indexing)
        SpMat Jtmp (bs, npass); // Mall x Mpass
        Jtmp.reserve(npass);
        for(int i = 0; i < npass; i++)
          Jtmp.insert(m_indices(i), i) = 1;

        if(!correct_vcov_burden(ph, rfrac, Qb(jcol), Kmat(m_indices, m_indices).sum(), GtWX(all, m_indices), XWsqrt, GWs * Jtmp, Wsqrt, phat, Y, mask, fest.cov_blup_offset, params))
          continue; // failed to correct with burden mask
      }
      if(apply_correction) block_info->cf_burden(ph) = rfrac;

      // get eigen values of Zt(I-U)Z
      if(!get_ztz_evals(rfrac * Kmat(m_indices, m_indices), r_outer_sum, gamma1, gamma2, gamma3, skat_lambda_tol, debug)) continue;

      if(npass > 1){
        get_skato_mom(skato_muQ, skato_fdavies, skato_sdQ, skato_dfQ, skato_tau, skato_lambdas, gamma1, gamma2, gamma3, rho_vec, debug);
        if(skato_sdQ < 0) continue; // failed
      }

      Qopt = (Qs(jcol) * flipped_skato_rho + Qb(jcol) * rho_vec).matrix(); // Nrho x 1
      if(debug) cerr << "Q:\n" << std::setprecision(10) << Qopt.transpose() << "\n";

      for(int j = 0; j < nrho; j++){

        // get eigen values of Rsqrt*ZtZ*Rsqrt
        get_lambdas(lambdas, get_RsKRs(rfrac * Kmat(m_indices, m_indices), r_outer_sum, gamma1, rho_vec(j), flip_rho_sqrt(j)), skat_lambda_tol);
        if(lambdas.size() == 0) continue;
        //if(rho_vec(j) >0.9) cerr << "rho=" << rho_vec(j) << "\nL:"<<lambdas.matrix().transpose() << "\n";

        // needed for skato (M>1)
        if(npass > 1)  get_cvals(j, cvals, lambdas);

        compute_fixed_skato_p(pvs_skato(j), chisq_skato(j), Qopt(j), rho_vec(j), lambdas, nl_dbl_dmin);

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
        pvs_m_o_acat(ph, imask) = get_acat(pvs_skato);
        get_chisq_stat_pv(tmpv, chisq_m_o_acat(ph, imask), pvs_m_o_acat(ph, imask), nl_dbl_dmin, log10_nl_dbl_dmin);
      } 
      if(with_acato){ // include acatv pvalue
        p_acato(0) = pvs_m_a(ph, imask);
        p_acato.tail(nrho) = pvs_skato;
        if(debug) cerr << "acato logp=" << p_acato.matrix().transpose() <<"\n";
        pvs_m_acato(ph, imask) = get_acat(p_acato);
        get_chisq_stat_pv(tmpv, chisq_m_acato(ph, imask), pvs_m_acato(ph, imask), nl_dbl_dmin, log10_nl_dbl_dmin);
      }
      if(with_skato_int){
        minp = max(nl_dbl_dmin, pow(10, -pvs_skato.maxCoeff())); // prevent underflow
        get_Qmin(nrho, minp, skato_Qmin_rho, cvals);
        if(debug) cerr << "Qmin=" << skato_Qmin_rho.matrix().transpose() << "\nminP=" << minp <<"; logp=" << pvs_skato.matrix().transpose() <<"\n";
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

}

void compute_vc_mats_bt(Ref<ArrayXd> Svals, Ref<MatrixXd> Kmat, const Ref<const MatrixXd>& XWsqrt, const Ref<const ArrayXd>& Wsqrt, const Ref<const MatrixXd>& yres, Ref<SpMat> Gmat, SpMat& GWs, MatrixXd& GtWX){

  // multiply by sqrt(p(1-p)) and mask entries (NxM)
  GWs = Wsqrt.matrix().asDiagonal() * Gmat; // NxM
  GtWX = XWsqrt.transpose() * GWs ; // CxM

  // get score stats for all variants (Mx1)
  // yres is Wsqrt^{-1}(Y-pi)
  Svals = (GWs.transpose() * yres).array();

  // kernel matrix for all variants (MxM)
  Kmat = - GtWX.transpose() * GtWX;
  Kmat += GWs.transpose() * GWs; // ZtZ

}

void compute_skat_q(VectorXd& Qs, VectorXd& Qb, const Ref<const ArrayXd>& Svals, Ref<MatrixXd> Kmat, const Ref<const MatrixXb>& Jmat, bool const& debug){

    Qs = Jmat.transpose().cast<double>() * Svals.square().matrix();
    // burden
    Qb = (Jmat.transpose().cast<double>() * Svals.matrix()).array().square().matrix();

    if(debug) cerr << "Q_SKAT for all masks:\n" << Qs.transpose() << "\nQ_BURDEN for all masks:\n" << Qb.transpose() << "\n";

}

void correct_vcov(int const& ph, const Ref<const ArrayXi>& indices, const Ref<const ArrayXd>& weights, Ref<ArrayXb> masked_sites, Ref<ArrayXd> Rvec_sqrt, const Ref<const ArrayXd>& score_stats, Ref<MatrixXd> Kmat, SpMat const& Gsparse, const Ref<const MatrixXd>& GtWX, const Ref<const MatrixXd>& XWsqrt, SpMat const& GWs, const Ref<const ArrayXd>& Wsqrt, const Ref<const ArrayXd>& phat, const Ref<const ArrayXd>& Y, const Ref<const ArrayXb>& mask, struct f_ests const& fest, struct param const& params){

  apply_correction_cc(ph, indices, weights, Rvec_sqrt, score_stats, Kmat.diagonal().array(), Gsparse, GtWX, XWsqrt, GWs, Wsqrt, phat, Y, mask, fest, params, true);
  if(params.debug) {
    int bs = Rvec_sqrt.size();
    cerr << "Rsqrt:" << Rvec_sqrt.head(min(150, bs)).matrix().transpose() << "\n";
  }

  // apply correction factor
  Kmat = Rvec_sqrt.matrix().asDiagonal() * Kmat * Rvec_sqrt.matrix().asDiagonal();
  masked_sites = (Rvec_sqrt > 0);
}

// when using lovo with bts
void check_cc_correction(SpMat& Gsparse, const Ref<const ArrayXd>& weights, const Ref<const MatrixXd>& X, struct ests const& m_ests, struct f_ests const& fest, const Ref<const MatrixXd>& yres, const Ref<const MatrixXd>& yraw, const Ref<const MatrixXb>& masked_indivs, struct param const& params){

  // if no correction, return 
  if(!(params.firth || params.use_SPA)) return;

  vc_Rvec_start.resize(weights.size(), params.n_pheno);
  vc_Rvec_start.array().colwise() = (weights > 0).cast<double>();

  ArrayXi indices;
  ArrayXd Svals, varS;
  MatrixXd GtWX, Kmat;
  SpMat GWs;
  Svals.resize(weights.size(), 1); // Mx1
  Kmat.resize(weights.size(), weights.size()); // MxM

  SpMat mat2 = Gsparse * weights.matrix().asDiagonal(); // include weights to G

  // loop over each trait
  for(int ph = 0; ph < params.n_pheno; ph++) { 
    if( !params.pheno_pass(ph) ) continue;

    MapcArXd Y (yraw.col(ph).data(), yraw.rows());
    MapcArXb mask (masked_indivs.col(ph).data(), yraw.rows());
    MapcArXd Wsqrt (m_ests.Gamma_sqrt.col(ph).data(), yraw.rows());
    MapcMatXd XWsqrt (m_ests.X_Gamma[ph].data(), yraw.rows(), m_ests.X_Gamma[ph].cols());
    MapcArXd phat (m_ests.Y_hat_p.col(ph).data(), yraw.rows());

    // get test stats with no correction
    compute_vc_mats_bt(Svals, Kmat, XWsqrt, Wsqrt * mask.cast<double>(), yres.col(ph), mat2, GWs, GtWX);
    varS = Kmat.diagonal().array();

    // apply correction and store Rvecs
    apply_correction_cc(ph, indices, weights, vc_Rvec_start.col(ph).array(), Svals, varS, mat2, GtWX, XWsqrt, GWs, Wsqrt, phat, Y, mask, fest, params, false);
  }

  if(params.debug) {
    int bs = vc_Rvec_start.rows();
    cerr << "Rsqrt_start:" << vc_Rvec_start.block(0, 0, min(150, bs), 1).transpose() << "\n";
  }

}

// correcting for high cc imbalance
void apply_correction_cc(int const& ph, const Ref<const ArrayXi>& indices, const Ref<const ArrayXd>& weights, Ref<ArrayXd> Rvec, const Ref<const ArrayXd>& score_stats, const Ref<const ArrayXd>& var_score, SpMat const& Gsparse, const Ref<const MatrixXd>& GtWX, const Ref<const MatrixXd>& XWsqrt, SpMat const& GWs, const Ref<const ArrayXd>& Wsqrt, const Ref<const ArrayXd>& phat, const Ref<const ArrayXd>& Y, const Ref<const ArrayXb>& mask, struct f_ests const& fest, struct param const& params, bool const& check_rvec_start){

  bool use_rvec_start = check_rvec_start && (vc_Rvec_start.size() > 0);
  int npass = Rvec.sum();

  // loop over the markers
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < Rvec.size(); i++) {
    if(Rvec(i) == 0) continue;
    
    // if already pre-computed
    if( use_rvec_start ){
      int ix = indices(i); // match indices (in mask vs in set)
      if( ix < vc_Rvec_start.rows() ) { // ignore ur masks
        Rvec(i) = vc_Rvec_start(ix, ph);
        continue;
      }
    }

    bool test_fail = true;
    double chisq, pv, corrected_var;

    // if Tstat < threshold, no correction done (R=1)
    double tstat_cur = score_stats(i) / sqrt(var_score(i));
    if(fabs(tstat_cur) <= params.z_thr) continue;

    MatrixXd Gres = - XWsqrt * GtWX.col(i); // get genotypic residuals
    Gres += GWs.col(i);

    if(params.use_SPA){ // SPA
      run_SPA_test_snp(chisq, pv, tstat_cur, var_score(i), true, Gsparse.col(i), Gres.array(), phat, Wsqrt, mask, test_fail, params.tol_spa, params.niter_max_spa, params.missing_value_double, params.nl_dbl_dmin);
    } else if(params.firth) { // Firth
      // remove skat weights as it can lead to different model fit for ur masks
      apply_firth_snp(test_fail, chisq, Gres.cwiseQuotient(Wsqrt.matrix()) / weights(i), Y, fest.cov_blup_offset.col(ph).array(), mask, params);
    }

    if( test_fail || (chisq == 0) ) { // set R to 0 for variant
      Rvec(i) = 0;
      continue;
    }

    if(params.debug) cerr << "uncorrected: " << tstat_cur * tstat_cur << " [=(" << score_stats(i) << ")^2/" << var_score(i) << "] -> " << chisq << endl;

    corrected_var = score_stats(i) * score_stats(i) / chisq;
    Rvec(i) = sqrt(corrected_var / var_score(i));

  }
#if defined(_OPENMP)
  setNbThreads(params.threads);
#endif

  int npass_post = (Rvec > 0).count();
  if(npass_post < npass) cerr << "WARNING: Firth/SPA correction failed for " << npass - npass_post << "/" << npass << " variants.";

}

// firth wrapper 
/*
void apply_firth_snp(bool& fail, double& lrt, const Ref<const MatrixXd>& Gvec, const Ref<const ArrayXd>& Y, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, struct param const& params) {

  double dev0 = 0, dev;
  ArrayXd betaold, se, etavec, pivec;
  betaold = ArrayXd::Zero(1); // start at 0

  fail = !fit_firth_pseudo(dev0, Y, Gvec, offset, mask, pivec, etavec, betaold, se, 1, dev, true, lrt, params.maxstep, params.niter_max_firth/2, params.numtol_firth, &params);

  if(!fail) return;

  betaold = 0; // start at 0
  fail = !fit_firth_nr(dev0, Y, Gvec, offset, mask, pivec, etavec, betaold, se, 1, dev, true, lrt, params.maxstep, params.niter_max_firth/2, params.numtol_firth, &params);

}
*/
void apply_firth_snp(bool& fail, double& lrt, const Ref<const MatrixXd>& Gvec, const Ref<const ArrayXd>& Y, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, struct param const& params) {

  double dev0 = 0, bstart = 0, betaold = bstart, se;
  ArrayXi index_carriers;

  // get dev0
  ArrayXd pivec, wvec;
  get_pvec(pivec, offset, params.numtol_eps);
  dev0 = get_logist_dev(Y, pivec, mask);
  get_wvec(pivec, wvec, mask);
  dev0 -= log( mask.select(Gvec.array().square() * wvec, 0).sum() );

  fail = fit_firth_pseudo(dev0, Y, Gvec.col(0), offset, mask, index_carriers, betaold, se, lrt, params.maxstep, params.niter_max_firth/2, params.numtol_firth, &params); // try pseudo

  if(!fail) return;

  betaold = bstart; // start at 0
  fail = !fit_firth(dev0, Y, Gvec, offset, mask, index_carriers, betaold, se, lrt, params.maxstep, params.niter_max_firth/2, params.numtol_firth, &params); // try NR (slower)

}

bool correct_vcov_burden(int const& ph, double& rfrac, double const& qb, double const& var_qb, const Ref<const MatrixXd>& GtWX, const Ref<const MatrixXd>& XWsqrt, SpMat const& GWs, const Ref<const ArrayXd>& Wsqrt, const Ref<const ArrayXd>& phat, const Ref<const ArrayXd>& Y, const Ref<const ArrayXb>& mask, const Ref<const MatrixXd>& offset, struct param const& params){

  if(qb == 0) return true; // no need to apply it

  // check T_burden
  double tstat_cur = sqrt(qb / var_qb);
  if(fabs(tstat_cur) <= params.z_thr) return true; // no need to apply it

  SpVec g_burden; // not needed since not using fastSPA
  bool test_fail = true;
  double chisq, pv;

  // get residuals for burden mask
  VectorXd g_res = GWs * VectorXd::Ones(GWs.cols()) - XWsqrt * GtWX.rowwise().sum(); // get mask residuals

  if( params.use_SPA ){ // use SPA
    run_SPA_test_snp(chisq, pv, tstat_cur, var_qb, false, g_burden, g_res.array(), phat, Wsqrt, mask, test_fail, params.tol_spa, params.niter_max_spa, params.missing_value_double, params.nl_dbl_dmin);
    /*if(params.debug && !test_fail)
      cerr << "SPA // uncorrected: " << tstat_cur * tstat_cur << " -> " << chisq <<
      ";logp="<< pv << ";rfrac=" << tstat_cur * tstat_cur / chisq << "\n";*/

  } else if( params.firth ){ // use firth
    apply_firth_snp(test_fail, chisq, g_res.cwiseQuotient(Wsqrt.matrix()), Y, offset.col(ph).array(), mask, params);
    /*if(params.debug && !test_fail)
      cerr << "Firth // uncorrected: " << tstat_cur * tstat_cur << " -> " << chisq <<
      ";logp="<< pv << ";rfrac=" << tstat_cur * tstat_cur / chisq << "\n";*/
  }

  if( test_fail || (chisq == 0) ) {
    if(params.debug) cerr << "WARNING: failed to correct T_burden.";
    return false;
  }

  // need to make variance bigger (so bigger p-values) so take max
  rfrac = max(1.0, tstat_cur * tstat_cur / chisq );
  if(params.debug) cerr << "T_burden=" << tstat_cur << ";R_factor_burden:" << rfrac << "\n";

  return true;
}

void get_single_pvs_bt(Ref<ArrayXd> pvals, const Ref<const ArrayXd>& chisq_vals){
  //cerr << pvals.matrix().transpose() << "\n\n" << chisq_vals.matrix().transpose() << endl;
  for(int isnp = 0; isnp < pvals.rows(); isnp++)
    get_logp(pvals(isnp), chisq_vals(isnp));
}


/////////////////////
/////////////////////
///// General
/////////////////////
/////////////////////

Eigen::MatrixXd get_RsKRs(const Ref<const MatrixXd>& K, const double& rho, const double& c1){

  int m = K.rows(); // M
  double c2 = sqrt( 1 - rho + m * rho), gamma1;

  VectorXd b = K.rowwise().sum(); // Mx1
  gamma1 = b.sum();

  return (
      (1-rho) * K.array() + 
      c1 * (c2-c1)/m * (b.rowwise().replicate(m) + b.transpose().colwise().replicate(m)).array() + // last term is outer sum of b
      (c2-c1)/m * (c2-c1)/m * gamma1
      ).matrix();
}

Eigen::MatrixXd get_RsKRs(const Ref<const MatrixXd>& K, const Ref<const MatrixXd>& b_outer_sum, const double& gamma1, const double& rho, const double& c1){

  int m = K.rows(); // M
  double c2 = sqrt( 1 - rho + m * rho);

  return (
      (1-rho) * K.array() + 
      c1 * (c2-c1)/m * b_outer_sum.array() +
      (c2-c1)/m * (c2-c1)/m * gamma1
      ).matrix();
}

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

void compute_fixed_skato_p(double& pval, double& chival, double const& Qs, double const& Qb, double const& rho, VectorXd& lambdas, const double& tol, bool const& debug){

  double q = (1 - rho) * Qs + rho * Qb;
  if(debug) cerr << "Q:" << q << "\n";

  if( (rho == 1) || (lambdas.size() == 1) ){ // burden or single variant
    chival = q / lambdas.tail(1)(0);
    get_logp(pval, chival); 
  } else compute_skat_pv(pval, chival, q, lambdas, tol);

}

void compute_fixed_skato_p(double& pval, double& chival, double& q, double const& rho, VectorXd& lambdas, const double& tol){

  if( (rho == 1) || (lambdas.size() == 1) ){ // burden or single variant
    chival = q / lambdas.tail(1)(0);
    get_logp(pval, chival); 
  } else compute_skat_pv(pval, chival, q, lambdas, tol);

}

void compute_skat_pv(double& logp, double& chival, double const& Q, VectorXd& lambdas, const double& tol){
  // use log10P directly to handle small pvalues
  logp = get_chisq_mix_logp(Q, lambdas, chival);
}

// returns p-value or -1
double get_chisq_mix_pv(double const& q, const Ref<const VectorXd>& lambdas){

  double pv, pv_davies_thr = 1e-5; // davies can be unreliable if pv is too small

  // re-scale so that max lambda is 1 (lambda is sorted)
  double newQ = q / lambdas.tail(1)(0);
  VectorXd newL = lambdas / lambdas.tail(1)(0);
  //cerr << "Qval= " << newQ << "\n";
  // exact
  pv = get_davies_pv(newQ, newL, false);
  //cerr << "davies: " << pv << "\n";

  // if failed or is very low, use SPA
  if(pv <= pv_davies_thr){ 
    pv = get_kuonen_pv(newQ, newL); // SPA
    //cerr << "kuonen: " << pv << "\n";
    if(pv <= 0) {// if SPA failed
      pv = get_davies_pv(newQ, newL, true); // use Davies with stringent parameters
      //cerr << "davies strict: " << pv << "\n";
      if(pv <= 0) {
        pv = get_liu_pv(newQ, newL); // only use mod Liu if Davies/SPA failed
        //cerr << "liu: " << pv << "\n";
      }
    }
  }

  return pv;

}

// get log10p or -1
double get_chisq_mix_logp(double const& q, const Ref<const VectorXd>& lambdas, double& chival){

  double logp, pv, pv_davies_thr = 1e-5; // davies can be unreliable if pv is too small
  double nl_dbl_dmin = 10.0 * std::numeric_limits<double>::min();
  double log10_nl_dbl_dmin = -log10(nl_dbl_dmin);

  // re-scale so that max lambda is 1 (lambda is sorted)
  double newQ = q / lambdas.tail(1)(0);
  VectorXd newL = lambdas / lambdas.tail(1)(0);
  //cerr << "Qval= " << newQ << "\n";
  // exact
  pv = get_davies_pv(newQ, newL, false);
  //cerr << "davies: " << pv << "\n";

  // if failed or is very low, use SPA
  if(pv <= pv_davies_thr){ 
    pv = get_kuonen_pv(newQ, newL); // SPA
    //cerr << "kuonen: " << pv << "\n";

    if(pv <= 0) {// if SPA failed
      pv = get_davies_pv(newQ, newL, true); // use Davies with stringent parameters
      //cerr << "davies strict: " << pv << "\n";

      if(pv <= 0) {
        logp = get_liu_pv(newQ, newL, chival); // only use mod Liu if Davies/SPA failed
        // get corresponding test stat for chisq(1)
        get_chisq_stat_pv(pv, chival, logp, nl_dbl_dmin, log10_nl_dbl_dmin);
        //cerr << "liu: " << logp << "\n";
      } else get_logp(pv, logp, chival, nl_dbl_dmin);

    } else get_logp(pv, logp, chival, nl_dbl_dmin);

  } else get_logp(pv, logp, chival, nl_dbl_dmin); 

  if(logp < 0) chival = -1;

  return logp;

}

// return 1-F(x) for chisq mixture
double get_davies_pv(double const& q, Ref<VectorXd> lambdas, bool const& force_stringent){

  // use default lim/acc values from CompQuadForm R package and SKAT resp.
  int k = lambdas.size(), ifault = 0, lim = 1e4; // p & error
  double cdf, pv, acc1 = 1e-6;
  if(force_stringent){ lim=1e6; acc1 = 1e-9;}
  ArrayXd nc = ArrayXd::Constant(k, 0); // ncp
  ArrayXi df = ArrayXi::Constant(k, 1); // df
  ArrayXd tr = ArrayXd::Constant(7, 0); // params for qf

  try {
    cdf = qf(lambdas.data(), nc.data(), df.data(), k, 0, q, lim, acc1, tr.data(), &ifault); 
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

double get_liu_pv(double const& q, const Ref<const VectorXd>& lambdas, const bool& lax){

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

  if(!lax && ((pv <= 0) || (pv > 1))) return -1;
  else if(lax && ((pv < 0) || (pv > 1))) return -1;
  return pv;
}


double get_liu_pv(double const& q, const Ref<const VectorXd>& lambdas, double& chival){

  ArrayXd cvals(6);
  get_cvals(cvals, lambdas);
  //cerr << "cvals liu=" << cvals.matrix().transpose() << endl;
  
  double pv, logpv;
  double tstar = (q - cvals(0)) * cvals(1);
  double val = tstar * cvals(3) + cvals(2);

  //cerr << "liu val = " << val << " ";
  if(val < 0) {
    chival = -1;
    return -1;
  }

  // 0 ncp gives strange behavior with non_central_chi_squared (returns -cdf instead of 1-cdf)
  if(cvals(5) == 0) get_logp(logpv, val, cvals(4));
  else  {
    pv = cdf(complement(non_central_chi_squared(cvals(4), cvals(5)), val));
    logpv = ( ((pv <= 0) || (pv > 1)) ? -1.0 : -log10(pv) );
  }

  //cerr << "; params = ( "  << std::setprecision(10) << cvals(4) << ", " << cvals(5) << ") -> logpv liu = " << logpv << "\n";

  if(logpv < 0) chival = -1;
  chival = val;
  return logpv;
}

bool get_ztz_evals(const Ref<const MatrixXd>& Kmat, MatrixXd& outer_sum, double& gamma1, double& gamma2, double& gamma3, double const& skat_lambda_tol, bool const& debug){

  VectorXd ZtZ_rsum = Kmat.rowwise().sum();//Kmx1
  outer_sum = ZtZ_rsum.rowwise().replicate(Kmat.cols()) + ZtZ_rsum.transpose().colwise().replicate(Kmat.rows());
  gamma1 = ZtZ_rsum.sum();
  gamma2 = ZtZ_rsum.squaredNorm();
  gamma3 = ZtZ_rsum.dot( Kmat * ZtZ_rsum);
  get_lambdas(skato_lambdas, Kmat - ZtZ_rsum * (ZtZ_rsum/gamma1).transpose(), skat_lambda_tol);
  if(skato_lambdas.size() == 0) return false;
  if(debug) {
    int bs = skato_lambdas.size();
    cerr << "L:\n" << skato_lambdas.head(min(150,bs)).transpose() << "\n";
  }

  return true;
}

void get_skato_mom(double& mu, double& sc_fac, double& sd, double& df, ArrayXd& tau, const Ref<const VectorXd>& lambdas, double const& gamma1, double const& gamma2, double const& gamma3, const Ref<const ArrayXd>& rho, bool const& debug){

  double v0, ve, vq;

  mu = lambdas.sum();
  v0 = 2 * lambdas.squaredNorm();
  ve = 4 * (gamma3/gamma1 - gamma2*gamma2/gamma1/gamma1);
  vq = v0 + ve;
  if(vq < 0){sd = -1; return;}
  sd = sqrt(vq);
  sc_fac = sqrt( v0 / vq );
  df = 0.5 * 0.5 * v0 * v0 / lambdas.array().pow(4).sum();
  if(debug) cerr << "[muQ, scFac, sd, df, v0, vq]= [" << mu << " " << sc_fac << " " << sd << " " << df << " " << v0 << " " << vq << " ]\n";
  tau = gamma1 * rho + gamma2/gamma1 * (1-rho);
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
  skato_upper = ((Qmin + flipped_skato_rho * skato_muQ * (1 - skato_fdavies) / skato_fdavies)/skato_tau).minCoeff();
}

double SKATO_integral_fn(double* x){ // variables used beside x are global

  double val = ((skato_Qmin_rho - skato_tau * (*x)) / flipped_skato_rho).minCoeff();
  double S, dlt;
  chi_squared chisq1( 1 );

  if(skato_state == 1) return 0; // skip if failed for other x values
  if(*x == 0) {skato_state = 1; return 0;} // failed

  // get first term in integral (1-cdf)
  if( val > (skato_muQ * 1e4) ) S = 0; // value check from SKAT R package
  else {
    dlt = (val - skato_muQ) * skato_fdavies + skato_muQ;
    if(dlt <= 0) S = 1;
    else{
      S = get_chisq_mix_pv(dlt, skato_lambdas);
      //cerr << *x << " " << S << " " << val << " " << skato_muQ << " " << skato_sdQ << " " << skato_dfQ << endl;

      if(S <= 0) { // failed
        skato_state = 1; 
        return 0;
      } else if(S >= 1) S = 1;
    }
  }

  return S * pdf(chisq1, *x);

}

double SKATO_integral_fn_liu(double* x){ // variables used beside x are global

  double val = ((skato_Qmin_rho - skato_tau * (*x)) / flipped_skato_rho).minCoeff();
  double S, dlt;
  chi_squared chisq1( 1 );

  if(skato_state == 1) return 0; // skip if failed for other x values
  if(*x == 0) {skato_state = 1; return 0;} // failed

  chi_squared chisqL( skato_dfQ );

  // get first term in integral
  dlt = (val - skato_muQ) / skato_sdQ * sqrt(2*skato_dfQ) + skato_dfQ;
  if(dlt<0) return 0; // cdf=0
  S = cdf(complement(chisqL, dlt));
  //cerr << *x << " " << S << " " << val << " " << skato_muQ << " " << skato_sdQ << " " << skato_dfQ << endl;

  return S * pdf(chisq1, *x);

}


// for skato num int
void integrate(double f(double*), double& pv, int const& subd, bool const& debug){

  int neval, ierror, ilimit = subd, last; 
  int lenw = 4 * ilimit;
  double lower = 0, upper = skato_upper, epsabs = 1e-25, epsrel = pow(std::numeric_limits<double>::epsilon(), .25), result, abserr;
  VectorXi iwork = VectorXi::Zero(ilimit);
  VectorXd work = VectorXd::Zero(lenw);
  skato_state = 0;

  dqags_(f, &lower, &upper, &epsabs, &epsrel, &result, &abserr, &neval, &ierror, &ilimit, &lenw, &last, iwork.data(), work.data());
  if(ierror != 0) skato_state = 1;
  if(debug) {
    cerr << "Niter=" << neval << ";integral=" << result << "Abs.error=" << abserr << ";rel.error=" << epsrel <<  ";fail="<< skato_state << "/" << ierror << "\n";
    if(skato_state == 0) for(int i = 1; i < 6; i++) {lower=skato_upper*0.2*i;cerr << "g(" << lower << ")=" << f(&lower) << " ";}
  }

  if (skato_state != 0)  pv = -1; 
  else pv = result;

} 

void get_skato_pv(double &logp, double& chisq, double const& minp, int const& nrhos, double const& nl_dbl_dmin, bool const& debug){

  double a, p_bc = minp * nrhos;
  chi_squared chisq1( 1 );
  double tstar = cdf(complement(chisq1, skato_upper)); 

  if(minp >= (1 - std::numeric_limits<float>::epsilon())) {logp = 0; chisq=0; return;}

  integrate(SKATO_integral_fn, a, 1000, debug);
  if(debug) cerr << "SKATO p=" << (skato_state == 0 ? (a+tstar) : -1) << "=" << a << "+" << tstar  << " (minP="<< minp <<"; Bonf=" << p_bc << ")\n";
  if(skato_state == 0) a += tstar; // add s(q*) to integral

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
          buffer << print_sum_stats_head_htp(snp_index, files.pheno_names[i], test_string + wgr_string + "-" + itr->first, snpinfo, params) << print_sum_stats_htp(-1, -1, itr->second(i, 0), itr->second(i, 1), -1, -1, -1, block_info->genocounts, i, true, 1, params, params->missing_value_double, -1, ( (params->firth || params->use_SPA) && ((itr->first == "SKATO-ACAT") || (itr->first == "SKATO")) ) ? block_info->cf_burden(i): -1.0, params->missing_value_double);
        else 
          buffer << (!params->split_by_pheno && (i>0) ? "" : header) << print_sum_stats(-1,-1,-1, -1, params->pheno_counts.row(i).sum(), params->pheno_counts(i, 0), params->pheno_counts(i, 1), test_string + "-" + itr->first, -1, -1, itr->second(i, 0), itr->second(i, 1), true, 1, params, (i+1));

        block_info->sum_stats[print_index].append( buffer.str() );
      } else if(!params->split_by_pheno) // print NA sum stats
        block_info->sum_stats[print_index].append( print_na_sumstats(i, 1, header, test_string + "-" + itr->first, block_info, *params) );

    }
  }

}

void check_sizes(SpMat const& Gmat_sp, SpMat const& Gmat_sp_urm, 
    //const Ref<const MatrixXd>& Gmat, 
    const Ref<const MatrixXb>& Jmat){

  cerr << "Printing sizes of SKAT objects\n" <<
    "-Gsparse = " << sizeof(double) * Gmat_sp.nonZeros() / 1024.0 / 1024.0 << "MB\n" <<
    "-Jmat = " << sizeof(Jmat(0,0)) * Jmat.size() / 1024.0 / 1024.0 << "MB\n" <<
    "-Gmat_sp_ur = " << sizeof(double) * Gmat_sp_urm.nonZeros() / 1024.0 / 1024.0 << "MB\n";

}
