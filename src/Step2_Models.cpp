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
#include "Files.hpp"
#include "Geno.hpp"
#include "Joint_Tests.hpp"
#include "Step1_Models.hpp"
#include "Step2_Models.hpp"
#include "HLM.hpp"
#include "Pheno.hpp"
#include "Masks.hpp"
#include "Data.hpp"

using namespace std;
using namespace Eigen;
using namespace boost;
using boost::math::normal;
using boost::math::chi_squared;


void blup_read_chr(bool const& silent, int const& chrom, struct ests& m_ests, struct in_files& files, struct filter const& filters, struct phenodt const& pheno_data, struct param& params, mstream& sout) {

  string line, filename, tmp_pheno;
  std::vector< string > id_strings, tmp_str_vec ;
  double in_blup;
  uint32_t indiv_index;
  Files blupf;

  // skip reading if specified by user or if PRS is given (same for all chromosomes)
  if( params.use_prs || params.skip_blups ) return;

  m_ests.blups = MatrixXd::Zero(params.n_samples, params.n_pheno);

  if(!silent) sout << "   -reading loco predictions for the chromosome..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  // read blup file for each phenotype
  for(int ph = 0; ph < params.n_pheno; ph++) {

    int i_pheno = files.pheno_index[ph];
    ArrayXb read_indiv = ArrayXb::Constant(params.n_samples, false);
    blupf.openForRead(files.blup_files[ph], sout);

    // check header
    blupf.readLine(line);
    id_strings = string_split(line,"\t ");
    if( id_strings[0] != "FID_IID") 
      throw "header of blup file must start with FID_IID.";

    // skip to chr
    blupf.ignoreLines(chrom-1);

    blupf.readLine(line);
    tmp_str_vec = string_split(line,"\t ");

    // check number of entries is same as in header
    if(tmp_str_vec.size() != id_strings.size()) 
      throw "blup file for phenotype [" + files.pheno_names[i_pheno] + "] has different number of entries on line " + to_string( chrom + 1 ) + " compared to the header (=" + to_string( tmp_str_vec.size() ) + " vs " + to_string( id_strings.size() ) + ").";

    // check starts with chromosome number
    if(chrStrToInt(tmp_str_vec[0], params.nChrom) != chrom) 
      throw "blup file for phenotype [" + files.pheno_names[i_pheno] + "] starts with `" +  tmp_str_vec[0]  + "`"
        + "instead of chromosome number=" + to_string( chrom ) + ".";

    // read blup data
    for( size_t filecol = 1; filecol < id_strings.size(); filecol++ ) {

      // ignore sample if it is not in genotype data
      if (!in_map(id_strings[filecol], params.FID_IID_to_ind)) continue;
      indiv_index = params.FID_IID_to_ind[id_strings[filecol]];

      // ignore sample if it is not included in analysis
      if(!filters.ind_in_analysis(indiv_index)) continue;

      // ignore sample if it is masked for the trait (prs will be 0)
      if(!pheno_data.masked_indivs(indiv_index,i_pheno)) continue;

      // check if duplicate
      if( !read_indiv(indiv_index) )
        read_indiv(indiv_index) = true;
      else 
        throw "individual appears more than once in blup file [" + files.blup_files[ph] + "]: FID_IID=" + id_strings[filecol];

      in_blup = convertDouble( tmp_str_vec[filecol], &params, sout);

      // if blup is NA then individual must be ignored in analysis for the phenotype (ie mask = 0)
      if (in_blup == params.missing_value_double)
        throw "individual has missing predictions (FID_IID=" + id_strings[filecol] + ";chr=" + to_string( chrom ) + ";phenotype=" + files.pheno_names[i_pheno] + 
          "). Either ignore these individuals using option '--remove', or skip reading predictions with option '--ignore-pred'.\n" + params.err_help ;
      else if(params.w_ltco && (chrom != params.ltco_chr)) // use ltco
        m_ests.blups(indiv_index, i_pheno) = in_blup - m_ests.ltco_prs(indiv_index, i_pheno);
      else // loco
        m_ests.blups(indiv_index, i_pheno) = in_blup;
    }

    // force all non-masked samples to have loco predictions
    //   -> this should not occur as masking of absent samples is done in blup_read() function
    if( (pheno_data.masked_indivs.col(i_pheno).array() && read_indiv).count() < pheno_data.masked_indivs.col(i_pheno).count() )
      throw "all samples included in the analysis (for phenotype " +
        files.pheno_names[i_pheno] + ") must have LOCO predictions in file : " + files.blup_files[ph] ;

    //cerr << m_ests.blups.col(i_pheno).head(5)<<endl;

    blupf.closeFile();
  }

  if(silent) return;

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;

}


/*
 * // tried this way for step 2 but it is slower than per SNP analysis
// marginal score test done for all variants/traits
void compute_score(vector<uint64> const& indices, int const& chrom, string const& test_string, string const& model_type, const Ref<const MatrixXd>& yres, const Ref<const RowVectorXd>& p_sd_yres, struct param const& params, struct phenodt& pheno_data, struct geno_block& gblock, vector<variant_block>& all_snps_info, vector<snp> const& snpinfo, struct ests const& m_ests, struct f_ests& fest, struct in_files const& files, mstream& sout){

  if(params.binary_mode)
    throw "not for BT";//compute_score_bt();
  else
    compute_score_qt(indices, chrom, test_string, model_type, yres, p_sd_yres, params, pheno_data, gblock, all_snps_info, snpinfo, files);

}


void compute_score_qt(vector<uint64> const& indices, int const& chrom, string const& test_string, string const& model_type, const Ref<const MatrixXd>& yres, const Ref<const RowVectorXd>& p_sd_yres, struct param const& params, struct phenodt& pheno_data, struct geno_block& gblock, vector<variant_block>& all_snps_info, vector<snp> const& snpinfo, struct in_files const& files){

  normal nd(0,1);
  double zcrit = quantile(complement(nd, .025));
  MatrixXd stats, bhat, scale_fac_pheno;

  if( params.strict_mode ) {

    gblock.Gmat.array().colwise() *= pheno_data.masked_indivs.col(0).array().cast<double>();
    stats = yres.transpose() * gblock.Gmat; // PxM
    stats.array() /= sqrt( params.n_analyzed - params.ncov );
    // estimate
    bhat = (stats.array().colwise() * (pheno_data.scale_Y.array() * p_sd_yres.array()).matrix().transpose().array()).matrix() / sqrt(params.n_analyzed - params.ncov); // need to divide by scale_fac for G block_info->scale_fac

  } else {

    // compute GtG for each phenotype (different missing patterns)
    scale_fac_pheno = pheno_data.masked_indivs.cast<double>().transpose() * gblock.Gmat.array().square().matrix(); // PxM, each element is GtG for each markerxtrait respecting missingness pattern
    stats = ((yres.transpose() * gblock.Gmat).array() / scale_fac_pheno.array().sqrt()).matrix(); // PxM
    // estimate
    bhat = ((stats.array().colwise() * (pheno_data.scale_Y.array() * p_sd_yres.array()).matrix().transpose().array()) / scale_fac_pheno.array().sqrt()).matrix(); // need to divide by scale_fac for G block_info->scale_fac

  }
  //cerr << stats.block(0,0,5,1).array().square() << endl;

  for(size_t isnp = 0; isnp < indices.size(); isnp++){

    variant_block* block_info = &(all_snps_info[isnp]);

    string tmpstr; // for sum stats
    if(!params.htp_out) tmpstr = print_sum_stats_head(indices[isnp], snpinfo);

    // beta
    block_info->bhat = bhat.col(isnp).array() / block_info->scale_fac;
    // SE
    block_info->se_b = block_info->bhat / stats.col(isnp).array();

    // get test statistic
    block_info->chisq_val = stats.col(isnp).array().square();

    for( int i = 0; i < params.n_pheno; ++i ) {

      if( block_info->ignored_trait(i) ) 
        continue;

      // get  pvalue
      get_logp(block_info->pval_log(i), block_info->chisq_val(i));

      if(!params.p_joint_only)
        block_info->sum_stats[i] = print_sum_stats_line(indices[isnp], i, tmpstr, test_string, model_type, block_info, snpinfo, files, params);

    }
  }

}
*/

// marginal score test for each snp
void compute_score(int const& isnp, int const& snp_index, int const& chrom, int const& thread_num, string const& test_string, string const& model_type, const Ref<const MatrixXd>& yres, const Ref<const RowVectorXd>& p_sd_yres, struct param const& params, struct phenodt& pheno_data, struct geno_block& gblock, variant_block* block_info, vector<snp> const& snpinfo, struct ests const& m_ests, struct f_ests& fest, struct in_files const& files, mstream& sout){

  if(params.binary_mode)
    compute_score_bt(isnp, snp_index, chrom, thread_num, test_string, model_type, yres, params, pheno_data, gblock, block_info, snpinfo, m_ests, fest, files, sout);
  else
    compute_score_qt(isnp, snp_index, thread_num, test_string, model_type, yres, p_sd_yres, params, pheno_data, gblock, block_info, snpinfo, files, sout);

}

// score test stat for QT
void compute_score_qt(int const& isnp, int const& snp_index, int const& thread_num, string const& test_string, string const& model_type, const Ref<const MatrixXd>& yres, const Ref<const RowVectorXd>& p_sd_yres, struct param const& params, struct phenodt& pheno_data, struct geno_block& gblock, variant_block* block_info, vector<snp> const& snpinfo, struct in_files const& files, mstream& sout){

  double gsc = block_info->flipped ? (4 * params.n_samples + block_info->scale_fac) : block_info->scale_fac;
  string tmpstr; // for sum stats
  MapArXd Geno (gblock.Gmat.col(isnp).data(), params.n_samples, 1);
  data_thread* dt_thr = &(gblock.thread_data[thread_num]);

  if( params.strict_mode ) {

    if(params.skip_blups && dt_thr->is_sparse) // Gsparse is on raw scale (must have yres centered)
      dt_thr->stats = (yres.transpose() * dt_thr->Gsparse.cwiseProduct(pheno_data.masked_indivs.col(0).cast<double>()) / gsc) / (sqrt( params.n_analyzed - params.ncov ));
    else
      dt_thr->stats = (yres.transpose() * (Geno * pheno_data.masked_indivs.col(0).cast<double>().array()).matrix()) / sqrt( params.n_analyzed - params.ncov );

    // estimate
    dt_thr->bhat = dt_thr->stats * ( pheno_data.scale_Y.array() * p_sd_yres.array()).matrix().transpose().array() / ( sqrt(params.n_analyzed - params.ncov) * gsc );

  } else {

    // compute GtG for each phenotype (different missing patterns)
    dt_thr->scale_fac_pheno = pheno_data.masked_indivs.transpose().cast<double>() * Geno.square().matrix();
    dt_thr->stats = (yres.transpose() * Geno.matrix()).array() / dt_thr->scale_fac_pheno.sqrt();

    // estimate
    dt_thr->bhat = dt_thr->stats * ( pheno_data.scale_Y.array() * p_sd_yres.array() ).matrix().transpose().array() / ( sqrt(dt_thr->scale_fac_pheno) * gsc );

  }

  // SE
  dt_thr->se_b = dt_thr->bhat / dt_thr->stats;

  // get test statistic
  dt_thr->chisq_val = dt_thr->stats.square();

  if(!params.htp_out) tmpstr = print_sum_stats_head(snp_index, snpinfo);

  for( int i = 0; i < params.n_pheno; ++i ) {

    if( block_info->ignored_trait(i) ) 
      continue;

    if(block_info->flipped) dt_thr->bhat(i) *= -1;

    // get pvalue
    get_logp(dt_thr->pval_log(i), dt_thr->chisq_val(i));

    if(!params.p_joint_only)
      block_info->sum_stats[i] = print_sum_stats_line(snp_index, i, tmpstr, test_string, model_type, block_info, dt_thr, snpinfo, files, params);

  }

}

void compute_score_bt(int const& isnp, int const& snp_index, int const& chrom, int const& thread_num, string const& test_string, string const& model_type, const Ref<const MatrixXd>& yres, struct param const& params, struct phenodt& pheno_data, struct geno_block& gblock, variant_block* block_info, vector<snp> const& snpinfo, struct ests const& m_ests, struct f_ests& fest, struct in_files const& files, mstream& sout){

  string tmpstr; 
  MatrixXd GW;
  SpVec GWs;
  data_thread* dt_thr = &(gblock.thread_data[thread_num]);

  // header snp info for sum stats
  if(!params.htp_out) tmpstr = print_sum_stats_head(snp_index, snpinfo);

  // genetype for marker
  MapArXd Geno (gblock.Gmat.col(isnp).data(), params.n_samples, 1);

  for( int i = 0; i < params.n_pheno; ++i ) {

    if( block_info->ignored_trait(i) ) 
      continue;
    MapArXb mask (pheno_data.masked_indivs.col(i).data(), params.n_samples, 1);

    // project out covariates from G
    if(dt_thr->is_sparse) {
      GWs = dt_thr->Gsparse.cwiseProduct((m_ests.Gamma_sqrt.col(i).array() * mask.cast<double>()).matrix());
      dt_thr->Gres = -m_ests.X_Gamma[i] * (m_ests.Xt_Gamma_X_inv[i] * (m_ests.X_Gamma[i].transpose() * GWs));
      dt_thr->Gres += GWs;
    } else {
      GW = (Geno * m_ests.Gamma_sqrt.col(i).array() * mask.cast<double>()).matrix();
      dt_thr->Gres = GW - m_ests.X_Gamma[i] * (m_ests.Xt_Gamma_X_inv[i] * (m_ests.X_Gamma[i].transpose() * GW));
    }

    // denominator
    dt_thr->denum(i) = dt_thr->Gres.squaredNorm();
    if( dt_thr->denum(i) < params.numtol ){
      block_info->ignored_trait(i) = true;
      continue;
    }

    // score test stat for BT
    if(dt_thr->is_sparse) 
      dt_thr->stats(i) = GWs.cwiseProduct(yres.col(i)).sum() / sqrt( dt_thr->denum(i) );
    else
      dt_thr->stats(i) = dt_thr->Gres.cwiseProduct(yres.col(i)).sum() / sqrt( dt_thr->denum(i) );

    // use firth/spa
    check_pval_snp(block_info, dt_thr, chrom, i, isnp, pheno_data, gblock, m_ests, fest, params, sout);

    dt_thr->bhat(i) /= block_info->scale_fac;
    dt_thr->se_b(i) /= block_info->scale_fac;
    if(block_info->flipped) dt_thr->bhat(i) *= -1;

    // print sum stats
    if(!params.p_joint_only)
      block_info->sum_stats[i] = print_sum_stats_line(snp_index, i, tmpstr, test_string, model_type, block_info, dt_thr, snpinfo, files, params);

  }


}

// Firth
bool fit_firth_logistic(int const& chrom, int const& ph, bool const& null_fit, struct param const* params, struct phenodt* pheno_data, struct ests const* m_ests, struct f_ests* fest, mstream& sout) {
  // if firth is used, fit based on penalized log-likelihood

  bool success, set_start = true;
  int col_incl;
  int maxstep = null_fit ? params->maxstep_null : params->maxstep;
  int niter = null_fit ? params->niter_max_firth_null : params->niter_max_firth;
  double dev, lrt;

  ArrayXd betaold, se, etavec, pivec, offset;
  MatrixXd Xmat;

  MapArXd Y (pheno_data->phenotypes_raw.col(ph).data(), pheno_data->phenotypes_raw.rows());
  MapArXb mask (pheno_data->masked_indivs.col(ph).data(), pheno_data->masked_indivs.rows());

  if(params->firth_approx){
    if(null_fit)
      Xmat = pheno_data->new_cov; // only covariates
    else 
      Xmat = fest->covs_firth.rightCols(1); // only tested SNP
    col_incl = Xmat.cols();
  } else {
    Xmat = fest->covs_firth; // covariates + tested SNP
    col_incl = Xmat.cols();
    if( null_fit ) col_incl--; // remove SNP column
  }

  offset = m_ests->blups.col(ph).array(); 
  // covariate effects added as offset in firth approx. (last entry of beta_null_firth = 0)
  if( params->firth_approx && !null_fit ) offset += (pheno_data->new_cov * fest->beta_null_firth.block(0,ph, pheno_data->new_cov.cols(),1)).array(); 

  // with firth approx. => trial 1: use maxstep_null
  // trial 2 => use fallback options (increase maxstep & niter)
  for(int trial = 0; trial < 2; trial++){

    // starting values
    if(set_start){
      if(null_fit){

        if(params->firth_approx && params->use_null_firth){
          betaold = fest->beta_null_firth.col(ph).head(Xmat.cols()).array();
        } else {
          betaold = ArrayXd::Zero(Xmat.cols());
          betaold(0) = ( 0.5 + mask.select(Y,0).sum())  / (pheno_data->Neff(ph) + 1);
          betaold(0) = log( betaold(0) / (1 - betaold(0) ));
          // LOCO prediction is offset
          betaold(0) -= mask.select(offset,0).mean();
        }

      } else {

        if(params->firth_approx) betaold = ArrayXd::Zero(col_incl); // only estimating effect of tested SNP
        else betaold = fest->beta_null_firth.array();

      }
    }

    success = fit_firth(ph, Y, Xmat, offset, mask, pivec, etavec, betaold, se, col_incl, dev, !null_fit, lrt, maxstep, niter, params, sout);

    if(params->firth_approx && null_fit){ // only retry for firth approx null model
      if(!params->fix_maxstep_null) { // don't retry with user-given settings
        if( !success ){ // if failed to converge
          sout << "WARNING: Logistic regression with Firth correction did not converge (maximum step size=" << maxstep <<";maximum number of iterations=" << niter <<").\n";
          maxstep = params->retry_maxstep_firth;
          niter = params->retry_niter_firth;
          if(trial == 0) sout << "Retrying with fallback parameters: (maximum step size=" << maxstep <<";maximum number of iterations=" << niter<<").\n";
          if(params->use_adam) set_start = false;
          continue;
        }
      }
    }

    break;
  }

  // If didn't converge
  if(!success){
    if(params->verbose && !params->firth_approx) sout << "WARNING: Logistic regression with Firth correction did not converge!\n";
    return false;
  }
  // sout << "\nNiter = " << niter_cur << " : " << mod_score.matrix().transpose() << endl;
  //cerr << betaold << endl;

  if(null_fit) {
    if(params->firth_approx) fest->beta_null_firth.block(0,ph,betaold.size(),1) = betaold.matrix();
    else fest->beta_null_firth = betaold.matrix();
  } else {
    // compute beta_hat
    fest->bhat_firth = betaold.tail(1)(0);
    // compute SE based on Hessian for unpenalized LL
    if(!params->back_correct_se)
      fest->se_b_firth = se.tail(1)(0);

    // compute LRT test stat. 
    fest->deviance_logistic = lrt;
    if( fest->deviance_logistic < 0 ) return false;
  }

  return true;
}


void fit_null_firth(bool const& silent, int const& chrom, struct f_ests* firth_est, struct phenodt* pheno_data, struct ests const* m_ests, struct in_files* files, struct param const* params, mstream& sout){

  auto t1 = std::chrono::high_resolution_clock::now();
  IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","");
  if(!silent) sout << "   -fitting null Firth logistic regression on binary phenotypes..." << flush;

  if(params->use_null_firth) // get starting values from file
    get_beta_start_firth(chrom, firth_est, files, params, sout);
  //cerr << firth_est->beta_null_firth.topRows(3) << "\n\n\n";

  const auto n_pheno = params->n_pheno;
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic) shared(m_ests, firth_est, sout)
#endif
  for( int i = 0; i < n_pheno; ++i ) {
    bool has_converged = fit_firth_logistic(chrom, i, true, params, pheno_data, m_ests, firth_est, sout);
    if(!has_converged) {
      string msg1 = to_string( (params->fix_maxstep_null ? params->maxstep_null : params->retry_maxstep_firth) );
      string msg2 = to_string( (params->fix_maxstep_null ? params->niter_max_firth_null : params->retry_niter_firth) );
      throw "Firth penalized logistic regression failed to converge for phenotype: " + files->pheno_names[i] + "." +
        " Try decreasing the maximum step size using `--maxstep-null` (currently=" + msg1 +  ") " +
        "and increasing the maximum number of iterations using `--maxiter-null` (currently=" + msg2 + ").";
    }

    if(params->write_null_firth)
      (*firth_est->firth_est_files[i]) << chrom << " " << firth_est->beta_null_firth.block(0,i,params->ncov,1).transpose().format(Fmt) << endl;
    //cerr << firth_est->beta_null_firth.topRows(3) << "\n\n\n";
  }
#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

  if(silent) return;

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;

}



void fit_firth_logistic_snp(int const& chrom, int const& ph, int const& isnp, bool const& null_fit, struct param const* params, struct phenodt* pheno_data, struct ests const* m_ests, struct f_ests const* fest, const Ref<const MatrixXd>& Gvec, variant_block* block_info, data_thread* dt_thr, mstream& sout) {
  // if firth is used, fit based on penalized log-likelihood

  bool success;
  int col_incl;
  int maxstep = null_fit ? params->maxstep_null : params->maxstep;
  int niter = null_fit ? params->niter_max_firth_null : params->niter_max_firth;
  double dev, lrt;

  ArrayXd betaold, se, etavec, pivec, offset;
  MatrixXd Xmat;

  MapArXd Y (pheno_data->phenotypes_raw.col(ph).data(), pheno_data->phenotypes_raw.rows());
  MapArXb mask (pheno_data->masked_indivs.col(ph).data(), pheno_data->masked_indivs.rows());

  if(params->firth_approx){
    if(null_fit){
      Xmat = pheno_data->new_cov; // only covariates
    } else {
      Xmat = Gvec; // only tested SNP
    }
    col_incl = Xmat.cols();
  } else {
    Xmat = MatrixXd::Zero(params->n_samples, pheno_data->new_cov.cols() + 1); // covariates + tested SNP
    Xmat << pheno_data->new_cov, Gvec;
    col_incl = Xmat.cols();
    if( null_fit ) col_incl--;
  }

  offset = m_ests->blups.col(ph).array(); 
  // covariate effects added as offset in firth approx.
  if( params->firth_approx && !null_fit ) offset += (pheno_data->new_cov * fest->beta_null_firth.block(0,ph,pheno_data->new_cov.cols(),1)).array(); 

  // starting values
  if(null_fit){

    // start at logit^-1(mean(Y))-mean(offset)
    betaold = ArrayXd::Zero(Xmat.cols()); // last entry in exact Firth is kept at 0
    betaold(0) = ( 0.5 + mask.select(Y,0).sum())  / (pheno_data->Neff(ph) + 1);
    betaold(0) = log( betaold(0) / (1 - betaold(0) ));

    // LOCO prediction is offset
    betaold(0) -= mask.select(offset,0).mean();

  } else {
    // start at 0
    if(params->firth_approx) betaold = ArrayXd::Zero(col_incl); 
    // start at estimate from null fit
    else betaold = dt_thr->beta_null_firth.col(0);
  }

  success = fit_firth(ph, Y, Xmat, offset, mask, pivec, etavec, betaold, se, col_incl, dev, !null_fit, lrt, maxstep, niter, params, sout);

  // If didn't converge
  if(!success){
    if(params->verbose) sout << "WARNING: Logistic regression with Firth correction did not converge!\n";
    block_info->test_fail(ph) = true;
    return ;
  }
  // sout << "\nNiter = " << niter_cur << " : " << mod_score.matrix().transpose() << endl;


  if(null_fit) {
    if(!params->firth_approx) dt_thr->beta_null_firth = betaold.matrix();
  } else {
    // compute beta_hat
    dt_thr->bhat(ph) = betaold.tail(1)(0);
    // compute SE based on Hessian for unpenalized LL
    if(!params->back_correct_se)
      dt_thr->se_b(ph) = se.tail(1)(0);

    if( lrt < 0 ) {
      block_info->test_fail(ph) = true;
      return ;
    }
    dt_thr->dif_deviance = lrt;
  }

  return ;
}

// use NR or ADAM for Firth
bool fit_firth(int const& ph, const Ref<const ArrayXd>& Y1, const Ref<const MatrixXd>& X1, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, ArrayXd& pivec, ArrayXd& etavec, ArrayXd& betavec, ArrayXd& sevec, int const& cols_incl, double& dev, bool const& comp_lrt, double& lrt, int const& maxstep_firth, int const& niter_firth, struct param const* params, mstream& sout) {

  double dev0 = 0;

  // get starting beta from ADAM ( compute and save null deviance )
  if(!comp_lrt && params->use_adam) 
    fit_firth_adam(ph, dev0, Y1, X1, offset, mask, pivec, etavec, betavec, sevec, cols_incl, dev, comp_lrt, lrt, params, sout);

  return fit_firth_nr(dev0, Y1, X1, offset, mask, pivec, etavec, betavec, sevec, cols_incl, dev, comp_lrt, lrt, maxstep_firth, niter_firth, params, sout);

}

// fit based on penalized log-likelihood using NR
bool fit_firth_nr(double& dev0, const Ref<const ArrayXd>& Y1, const Ref<const MatrixXd>& X1, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, ArrayXd& pivec, ArrayXd& etavec, ArrayXd& betavec, ArrayXd& sevec, int const& cols_incl, double& dev, bool const& comp_lrt, double& lrt, int const& maxstep_firth, int const& niter_firth, struct param const* params, mstream& sout) {
  // fit with first cols_incl columns of X1 (non-used entries of betavec should be 0)
  // else assuming using all columns 

  bool use_offset = Y1.size() == offset.size();
  int niter_cur = 0, nc = X1.cols();
  double dev_old=0, dev_new=0, denum, mx;

  ArrayXd hvec, mod_score;
  ArrayXd betanew, step_size, wvec;
  MatrixXd XtW, XtWX;
  ColPivHouseholderQR<MatrixXd> qr, qrX;

  // solve S'(beta) = S(beta) + X'(h*(0.5-p)) = 0
  betanew = betavec * 0;
  while(niter_cur++ < niter_firth){

    // update quantities
    etavec = (X1 * betavec.matrix()).array();
    if(use_offset) etavec += offset;
    pivec = 1 - 1 / (etavec.exp() + 1) ;
    wvec = mask.select( ( pivec * (1 - pivec) ).sqrt(), 0);
    XtW = X1.transpose() * wvec.matrix().asDiagonal();
    XtWX = XtW * XtW.transpose();
    qr.compute(XtWX);
    // compute deviance
    dev_old = get_logist_dev(Y1, pivec, mask) - qr.logAbsDeterminant();
    if(comp_lrt && (niter_cur == 1)) // at first iter (i.e. betaSNP=0) this is null deviance (if not using ADAM)
      dev0 = dev_old;
    // compute diag(H), H = U(U'U)^{-1}U', U = Gamma^(1/2)X
    hvec = (qr.solve(XtW).array() * XtW.array() ).colwise().sum();
    // compute modified score & step size
    if(cols_incl < nc) { 
      qrX.compute(XtWX.block(0, 0, cols_incl, cols_incl));
      mod_score = (X1.leftCols(cols_incl).transpose() * mask.select( Y1 - pivec + hvec * (0.5 - pivec), 0).matrix() ).array();
      step_size = qrX.solve( mod_score.matrix() ).array();
    } else {
      mod_score = (X1.transpose() * mask.select( Y1 - pivec + hvec * (0.5 - pivec), 0).matrix() ).array();
      step_size = qr.solve( mod_score.matrix() ).array();
    }

    // stopping criterion using modified score function
    if( mod_score.abs().maxCoeff() < params->numtol_firth ) break;

    // force absolute step size to be less than maxstep for each entry of beta
    mx = step_size.abs().maxCoeff() / maxstep_firth;
    if( mx > 1 ) step_size /= mx;

    // start step-halving and stop when deviance decreases 
    denum = 1;
    for( int niter_search = 1; niter_search <= params->niter_max_line_search; niter_search++ ){

      // adjusted step size
      step_size /= denum;

      ///////// compute corresponding deviance
      if(cols_incl < nc) 
        betanew.head(cols_incl) = betavec.head(cols_incl) + step_size;
      else 
        betanew = betavec + step_size;
      etavec = (X1 * betanew.matrix()).array();
      if(use_offset) etavec += offset;

      pivec = 1 - 1 / (etavec.exp() + 1) ;
      wvec = mask.select( ( pivec * (1 - pivec) ).sqrt(), 0);
      XtW = X1.transpose() * wvec.matrix().asDiagonal();
      XtWX = XtW * XtW.transpose();
      qr.compute(XtWX);
      dev_new = get_logist_dev(Y1, pivec, mask) - qr.logAbsDeterminant();

      //sout << "\n["<<niter_cur << " - " << niter_search <<"]  denum =" << denum << ";\n step =" << step_size.matrix().transpose().array() / denum<<"; \nbeta=" << betanew.matrix().transpose().array() << ";\n Lnew= " << dev_new << " vs L0="<< dev_old << ";score="<< mod_score<< endl;
      if( dev_new < dev_old + params->numtol ) break;
      denum *= 2;
    }

    if(cols_incl < nc)  
      betavec.head(cols_incl) += step_size;
    else
      betavec += step_size;
    dev_old = dev_new;

    //if(niter_cur>1) sout << "\nNiter = " << niter_cur << " (beta = " << betanew.matrix().transpose() << ") : " << mod_score.matrix().transpose() << endl;

  }

  // If didn't converge
  if(niter_cur > niter_firth) return false;
  // sout << "\nNiter = " << niter_cur << " : " << mod_score.matrix().transpose() << endl;

  dev = dev_new;
  if( comp_lrt ) {
    lrt = dev0 - dev_new;
    if(lrt < 0) return false;

    sevec = qr.inverse().diagonal().array().sqrt();
  }

  return true;
}


// fit based on penalized log-likelihood using ADAM
bool fit_firth_adam(int const& ph, double& dev0, const Ref<const ArrayXd>& Y1, const Ref<const MatrixXd>& X1, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, ArrayXd& pivec, ArrayXd& etavec, ArrayXd& betavec, ArrayXd& sevec, int const& cols_incl, double& dev, bool const& comp_lrt, double& lrt, struct param const* params, mstream& sout) {
  // fit with first cols_incl columns of X1 (non-used entries of betavec should be 0)
  // else assuming using all columns 

  bool use_offset = Y1.size() == offset.size();
  bool force_batch_adam = comp_lrt && params->adam_mini; // force batch adam for 1st iteration to get dev0
  int niter_cur = 0, index;
  double p_alpha = params->adam_alpha, p_beta1 = params->adam_beta1, p_beta2 = params->adam_beta2, p_eps = params->adam_eps, p_alpha_t;

  std::uniform_int_distribution<> d(0, mask.count() - 1);
  std::mt19937 gen;
  ArrayXd hvec, gradient_f, wvec;
  ArrayXd mt, vt, step_size, Ytmp, offset_tmp;
  MatrixXd XtW, XtWX, Xtmp;
  ColPivHouseholderQR<MatrixXd> qr;

  // starting values for ADAM params
  mt = vt = betavec.head(cols_incl) * 0;

  // for mini-batch
  if(params->adam_mini){
    Xtmp.resize(params->adam_batch_size, X1.cols());
    Ytmp.resize(params->adam_batch_size);
    if(use_offset) offset_tmp.resize(params->adam_batch_size);
  }

  // minimize f=-2*pen.LL using ADAM
  while(niter_cur++ < params->niter_max_firth_adam){

    if(params->adam_mini && !force_batch_adam){ // ADAM using mini-batch

      for (int i = 0; i < params->adam_batch_size; i++) {
        index = params->adam_indices[ph](d(gen));
        Xtmp.row(i) = X1.row(index);
        Ytmp(i) = Y1(index);
        if(use_offset) offset_tmp(i) = offset(index);
      }
      // update quantities
      etavec = (Xtmp * betavec.matrix()).array();
      if(use_offset) etavec += offset_tmp;
      // fitted probabilities
      pivec = 1 - 1 / (etavec.exp() + 1) ;
      wvec = pivec * (1 - pivec);
      XtW = Xtmp.transpose() * wvec.matrix().asDiagonal();
      XtWX = XtW * XtW.transpose();
      qr.compute(XtWX);
      // compute diag(H), H = U(U'U)^{-1}U', U = Gamma^(1/2)X
      hvec = (qr.solve(XtW).array() * XtW.array() ).colwise().sum();
      // compute gradient of f
      gradient_f = - (Xtmp.leftCols(cols_incl).transpose() * (Ytmp - pivec + hvec * (0.5 - pivec)).matrix() ).array();

    } else {

      // update quantities
      etavec = (X1 * betavec.matrix()).array();
      if(use_offset) etavec += offset;
      // fitted probabilities
      pivec = 1 - 1 / (etavec.exp() + 1) ;
      wvec = mask.select( ( pivec * (1 - pivec) ).sqrt(), 0);
      XtW = X1.transpose() * wvec.matrix().asDiagonal();
      XtWX = XtW * XtW.transpose();
      qr.compute(XtWX);
      // at first iter (i.e. betaSNP=0) this is null deviance
      if(comp_lrt && (niter_cur == 1))
        dev0 = get_logist_dev(Y1, pivec, mask) - qr.logAbsDeterminant();
      // compute diag(H), H = U(U'U)^{-1}U', U = Gamma^(1/2)X
      hvec = (qr.solve(XtW).array() * XtW.array() ).colwise().sum();
      // compute gradient of f
      gradient_f = - (X1.leftCols(cols_incl).transpose() * mask.select( Y1 - pivec + hvec * (0.5 - pivec), 0).matrix() ).array();

    }

    //if(niter_cur>1 && (niter_cur%5000)==0) sout << "\nNiter = " << niter_cur << " (beta = " << betavec.matrix().transpose() << ") : " << gradient_f.matrix().transpose() << endl;

    mt = p_beta1 * mt + (1 - p_beta1) * gradient_f;
    vt = p_beta2 * vt + (1 - p_beta2) * gradient_f.square();
    p_alpha_t = p_alpha * sqrt(1 - pow(p_beta2, niter_cur)) / (1 - pow(p_beta1, niter_cur));
    step_size = p_alpha_t * mt / (vt.sqrt() + p_eps);

    // stopping criterion
    if( step_size.abs().maxCoeff() < params->numtol) break;

    betavec.head(cols_incl) -= step_size;

  }
  if(params->verbose) sout << "ADAM took "<< niter_cur << " iterations (score max = " << gradient_f.abs().maxCoeff() << ")...";

  return (niter_cur <= params->niter_max_firth_adam);
}


string get_firth_est_allChr(struct in_files& files, struct filter const& filters, struct ests& m_ests, struct f_ests& fest, struct phenodt& pheno_data, struct param& params, mstream& sout){

  sout << "   -computing and storing Firth estimates for all chromosomes..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  // go through each chromosome
  for(int chr = 1; chr <= params.nChrom; chr++){

    if(params.verbose) sout << "chr" << chr <<"..." << flush;

    // read the prs
    blup_read_chr(true, chr, m_ests, files, filters, pheno_data, params, sout);

    // run null firth for each trait and write estimates to file
    fit_null_firth(true, chr, &fest, &pheno_data, &m_ests, &files, &params, sout);

  }

  string fname = print_null_firth_info(files, fest, params);

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;

  sout << "    +written to file [" << fname << "]\n";

  return fname;

}

string print_null_firth_info(struct in_files const& files, struct f_ests& fest, struct param const& params){

    string path_firth, firth_filename, out_firth_list = files.out_file + "_firth.list";
    
    ofstream outf;
    outf.open(out_firth_list);
    if(outf.fail())
      throw "cannot write file : " + out_firth_list;

    for( int j = 0; j < params.n_pheno; ++j ) {
      firth_filename = files.out_file + "_" + to_string(j + 1) + ".firth" + (params.gzOut ? ".gz" : "");
      path_firth = get_fullpath(firth_filename);
      outf << files.pheno_names[j]  << " " <<  path_firth << endl;
      fest.firth_est_files[j]->closeFile();
    }

    outf.close();

    return out_firth_list;
}

void check_beta_start_firth(struct in_files& files, struct param const& params, mstream& sout){


  int tmp_index;
  string line;
  std::vector< string > tmp_str_vec;
  ArrayXb read = ArrayXb::Constant(params.n_pheno, false);
  Files fClass;

  // get list of files containing blups
  files.null_firth_files.resize(params.n_pheno);
  fClass.openForRead(files.null_firth_file, sout);

  while (fClass.readLine(line)){
    tmp_str_vec = string_split(line,"\t ");

    // each line contains a phenotype name and the corresponding blup file name
    if( tmp_str_vec.size() != 2 )
      throw "incorrectly formatted blup list file : " + files.null_firth_file;


    // get index of phenotype in phenotype matrix
    vector<string>::iterator it = std::find(files.pheno_names.begin(), files.pheno_names.end(), tmp_str_vec[0]);
    if (it == files.pheno_names.end()) continue; // ignore unrecognized phenotypes

    tmp_index = std::distance(files.pheno_names.begin(), it);
    files.null_firth_files[tmp_index] = tmp_str_vec[1];

    // check that phenotype only has one file
    if(read(tmp_index))
      throw "phenotype \'" + tmp_str_vec[0] + "\' appears more than once in file.";
    else if(!file_exists(tmp_str_vec[1]))
      throw "file " + tmp_str_vec[1] + " cannot be opened.";

    read(tmp_index) = true;
  }

  // force all phenotypes in phenotype file to be used
  if(read.count() != params.n_pheno) 
    throw "number of files (" + to_string( read.count() ) + ")  is not equal to the number of phenotypes." ;

}

void get_beta_start_firth(int const& chrom, struct f_ests* firth_est, struct in_files* files, struct param const* params, mstream& sout){

  int npar, nmax = firth_est->beta_null_firth.rows();
  double in_beta;
  string line;
  std::vector< string > tmp_str_vec ;
  Files fClass;

  // for each phenotype get b0
  for( int i = 0; i < params->n_pheno; ++i ) {

    bool chr_found = false;
    fClass.openForRead(files->null_firth_files[i], sout);

    while(fClass.readLine(line)){
      tmp_str_vec = string_split(line,"\t ");
      if(tmp_str_vec.size() == 0)
        throw "error reading null firth estimates file";
      else if(chrStrToInt(tmp_str_vec[0], params->nChrom) == chrom) {
        chr_found = true; break;
      }
    }

    //cerr << std::boolalpha << chr_found << endl;
    if(!chr_found) continue; // use 0 as start

    npar = tmp_str_vec.size();
    if((npar-1) > nmax) 
      throw "file has more predictors than included in analysis (=" + to_string(npar) + " vs " + to_string(nmax) + ")";

    for(int j = 1; j < npar; j++ ) {
      in_beta = convertDouble( tmp_str_vec[j], params, sout);
      if (in_beta == params->missing_value_double)
        throw "no missing values allowed in file";
      firth_est->beta_null_firth(j-1,i) = in_beta;
    }

    fClass.closeFile();
  }

}


void check_pval_snp(variant_block* block_info, data_thread* dt_thr, int const& chrom, int const& ph, int const& isnp, struct phenodt& pheno_data, struct geno_block& gblock, struct ests const& m_ests, struct f_ests& fest, struct param const& params, mstream& sout){

  // if firth isn't used, or Tstat < threshold, no correction done
  if(!block_info->is_corrected(ph) || (fabs(dt_thr->stats(ph)) <= params.z_thr)){
    get_sumstats(false, ph, dt_thr);
    block_info->is_corrected(ph) = false;
    return;
  }

  if(params.firth){ // firth
    
    run_firth_correction_snp(chrom, ph, isnp, gblock, block_info, dt_thr, pheno_data, m_ests, fest, params, sout);
    if(block_info->test_fail(ph)) {
      get_sumstats(true, ph, dt_thr);
      return;
    }

    dt_thr->chisq_val(ph) = dt_thr->dif_deviance;
    get_logp(dt_thr->pval_log(ph), dt_thr->chisq_val(ph));

    // compute SE from beta & pvalue
    if( params.back_correct_se && (dt_thr->chisq_val(ph) > 0) )
      dt_thr->se_b(ph) = fabs(dt_thr->bhat(ph)) / sqrt(dt_thr->chisq_val(ph));

  } else if(params.use_SPA) { // spa

    run_SPA_test_snp(block_info, dt_thr, ph, pheno_data.masked_indivs.col(ph).array(), m_ests, params, sout);
    if(block_info->test_fail(ph)) {
      get_sumstats(true, ph, dt_thr);
      return;
    }

    dt_thr->se_b(ph) = 1 / sqrt(dt_thr->denum(ph));
    dt_thr->bhat(ph) = sgn(dt_thr->stats(ph)) * sqrt(dt_thr->chisq_val(ph)) * dt_thr->se_b(ph);

  }

}

void get_sumstats(bool const& no_pv, int const& ph, data_thread* dt_thr) {

  // beta & se
  dt_thr->se_b(ph) = 1 / sqrt(dt_thr->denum(ph));
  dt_thr->bhat(ph) = dt_thr->stats(ph) * dt_thr->se_b(ph);
  if(no_pv) return;

  // chisq & lpv
  dt_thr->chisq_val(ph) = pow(dt_thr->stats(ph), 2);
  get_logp(dt_thr->pval_log(ph), dt_thr->chisq_val(ph));
}

void run_firth_correction_snp(int const& chrom, int const& ph, int const& isnp, struct geno_block& gblock, variant_block* block_info, data_thread* dt_thr, struct phenodt& pheno_data, struct ests const& m_ests, struct f_ests& fest, struct param const& params, mstream& sout){

  if(!params.firth_approx){ // exact firth
    // obtain null deviance (set SNP effect to 0 and compute max. pen. LL)
    fit_firth_logistic_snp(chrom, ph, isnp, true, &params, &pheno_data, &m_ests, &fest, gblock.Gmat.col(isnp), block_info, dt_thr, sout);
    if(block_info->test_fail(ph)) return ;
    // fit full model and compute deviance
    fit_firth_logistic_snp(chrom, ph, isnp, false, &params, &pheno_data, &m_ests, &fest, gblock.Gmat.col(isnp), block_info, dt_thr, sout);
  } else // approx firth - only fit full model
  fit_firth_logistic_snp(chrom, ph, isnp, false, &params, &pheno_data, &m_ests, &fest, dt_thr->Gres.cwiseQuotient(m_ests.Gamma_sqrt.col(ph)), block_info, dt_thr, sout);

}



void run_SPA_test_snp(variant_block* block_info, data_thread* dt_thr, int const& ph, const Ref<const ArrayXb>& mask, struct ests const& m_ests, struct param const& params, mstream& sout){

  int index_j;
  double score_num, tval, limK1_low, limK1_high, root_K1;
  ArrayXd Gmu;

  // compute needed quantities
  dt_thr->val_c = sqrt( dt_thr->denum(ph) );  // sqrt( G'WG )
  score_num = dt_thr->stats(ph) * dt_thr->val_c;
  dt_thr->Gmod = dt_thr->Gres.array() / m_ests.Gamma_sqrt.col(ph).array() * mask.cast<double>();
  Gmu = dt_thr->Gmod * m_ests.Y_hat_p.col(ph).array();
  dt_thr->val_a = Gmu.sum();

  if(dt_thr->fastSPA){
    dt_thr->val_b = dt_thr->denum(ph);
    dt_thr->val_d = 0;
    for (SparseVector<double>::InnerIterator it(dt_thr->Gsparse); it; ++it) {
      index_j = it.index();
      if(!mask(index_j)) continue;
      dt_thr->val_b -= dt_thr->Gres(index_j) * dt_thr->Gres(index_j);
      dt_thr->val_d += Gmu(index_j);
    }
  }

  // check if K'(t)= s can be solved
  limK1_low = (dt_thr->Gmod < 0).select(dt_thr->Gmod, 0 ).sum() - dt_thr->val_a ;
  limK1_high = (dt_thr->Gmod > 0).select(dt_thr->Gmod, 0 ).sum() - dt_thr->val_a ;
  if( score_num < limK1_low || score_num > limK1_high ){
    if(params.verbose) sout << "WARNING: SPA failed (solution to K'(t)=s is infinite)";
    block_info->test_fail(ph) = true;
    return;
  }

  // keep track of whether obs stat is positive
  dt_thr->pos_score = dt_thr->stats(ph) > 0;
  tval = fabs(dt_thr->stats(ph));

  // solve K'(t)= tval using a mix of Newton-Raphson and bisection method
  root_K1 = solve_K1_snp(tval, ph, &params, &m_ests, dt_thr, mask, sout);
  if( root_K1 == params.missing_value_double ){
    block_info->test_fail(ph) = true;
    return;
  }

  // compute pvalue
  get_SPA_pvalue_snp(root_K1, tval, ph, &params, &m_ests, block_info, dt_thr, mask, sout);

}



// SPA (MT in OpenMP)
double solve_K1_snp(const double& tval, const int& ph, const struct param* params, const struct ests* m_ests, data_thread* dt_thr, const Ref<const ArrayXb>& mask, mstream& sout){

  int niter_cur;
  int lambda = dt_thr->pos_score ? 1 : -1; // if score is negative, adjust K' and K''
  double min_x, max_x, t_old, f_old, t_new = -1, f_new, hess;

  niter_cur = 0;
  min_x = 0, max_x = std::numeric_limits<double>::infinity();
  t_old = 0;
  f_old = dt_thr->fastSPA ? compute_K1_fast_snp(lambda * t_old, m_ests, dt_thr, mask, ph) : compute_K1_snp(lambda * t_old, m_ests, dt_thr, mask, ph);
  f_old *= lambda;
  f_old -= tval; 

  while( niter_cur++ < params->niter_max_spa ){

    hess = dt_thr->fastSPA ? compute_K2_fast_snp(lambda * t_old, m_ests, dt_thr, mask, ph) : compute_K2_snp(lambda * t_old, m_ests, dt_thr, mask, ph);
    t_new = t_old - f_old / hess;
    f_new = dt_thr->fastSPA ? compute_K1_fast_snp(lambda * t_new, m_ests, dt_thr, mask, ph) : compute_K1_snp(lambda * t_new, m_ests, dt_thr, mask, ph);
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
      f_new = dt_thr->fastSPA ? compute_K1_fast_snp(lambda * t_new, m_ests, dt_thr, mask, ph) : compute_K1_snp(lambda * t_new, m_ests, dt_thr, mask, ph);
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

double compute_K_snp(const double& t, const struct ests* m_ests, data_thread* dt_thr, const Ref<const ArrayXb>& mask, const int& ph){

  double val = mask.select( ( 1 - m_ests->Y_hat_p.col(ph).array() + m_ests->Y_hat_p.col(ph).array() * ( t / dt_thr->val_c * dt_thr->Gmod ).exp() ).log(), 0).sum() - t * dt_thr->val_a / dt_thr->val_c;

  return val;
}

double compute_K_fast_snp(const double& t, const struct ests* m_ests, data_thread* dt_thr, const Ref<const ArrayXb>& mask, const int& ph){

  uint32_t index_j;
  double val = 0;

  for (SparseVector<double>::InnerIterator it(dt_thr->Gsparse); it; ++it) {
    index_j = it.index();
    if(!mask(index_j)) continue;

    val += log( 1 - m_ests->Y_hat_p(index_j,ph) + m_ests->Y_hat_p(index_j,ph) * exp( t / dt_thr->val_c * dt_thr->Gmod(index_j)) );
  }
  val += -t * dt_thr->val_d / dt_thr->val_c + t * t / 2 / dt_thr->denum(ph) * dt_thr->val_b;

  return val;
}

double compute_K1_snp(const double& t, const struct ests* m_ests, data_thread* dt_thr, const Ref<const ArrayXb>& mask, const int& ph){

  double val = mask.select( ( dt_thr->Gmod * m_ests->Y_hat_p.col(ph).array() / dt_thr->val_c ) / ( m_ests->Y_hat_p.col(ph).array() + (1 - m_ests->Y_hat_p.col(ph).array()) * ( -t / dt_thr->val_c * dt_thr->Gmod).exp() ), 0).sum();
  val -= dt_thr->val_a / dt_thr->val_c;

  return val;
}

double compute_K1_fast_snp(const double& t, const struct ests* m_ests, data_thread* dt_thr, const Ref<const ArrayXb>& mask, const int& ph){

  uint32_t index_j;
  double val = 0;

  for (SparseVector<double>::InnerIterator it(dt_thr->Gsparse); it; ++it) {
    index_j = it.index();
    if(!mask(index_j)) continue;

    val += ( dt_thr->Gmod(index_j) * m_ests->Y_hat_p(index_j,ph) / dt_thr->val_c ) / ( m_ests->Y_hat_p(index_j,ph) + (1 - m_ests->Y_hat_p(index_j,ph)) * exp( -t / dt_thr->val_c * dt_thr->Gmod(index_j)) );
  }
  val += -dt_thr->val_d / dt_thr->val_c + t / dt_thr->denum(ph) * dt_thr->val_b;

  return val;
}

double compute_K2_snp(const double& t, const struct ests* m_ests, data_thread* dt_thr, const Ref<const ArrayXb>& mask, const int& ph){

  double val = mask.select( ( dt_thr->Gmod.square() * m_ests->Gamma_sqrt.col(ph).array().square() / (dt_thr->val_c*dt_thr->val_c) * ( -t / dt_thr->val_c * dt_thr->Gmod).exp()) / ( m_ests->Y_hat_p.col(ph).array() + (1 - m_ests->Y_hat_p.col(ph).array()) * ( -t / dt_thr->val_c * dt_thr->Gmod).exp() ).square(), 0).sum();

  return val;
}

double compute_K2_fast_snp(const double& t, const struct ests* m_ests, data_thread* dt_thr, const Ref<const ArrayXb>& mask, const int& ph){

  uint32_t index_j;
  double val = 0, denum;

  for (SparseVector<double>::InnerIterator it(dt_thr->Gsparse); it; ++it) {
    index_j = it.index();
    if(!mask(index_j)) continue;

    denum = m_ests->Y_hat_p(index_j,ph) + (1 - m_ests->Y_hat_p(index_j,ph)) * exp( -t / dt_thr->val_c * dt_thr->Gmod(index_j));
    val += ( dt_thr->Gmod(index_j) * dt_thr->Gmod(index_j) * m_ests->Gamma_sqrt(index_j,ph) * m_ests->Gamma_sqrt(index_j,ph) * exp( -t / dt_thr->val_c * dt_thr->Gmod(index_j)) / dt_thr->val_c / dt_thr->val_c ) / (denum * denum);
  }
  val += dt_thr->val_b / dt_thr->denum(ph);

  return val;
}

void get_SPA_pvalue_snp(const double& root, const double& tval, const int& ph, struct param const* params, const struct ests* m_ests, variant_block* block_info, data_thread* dt_thr, const Ref<const ArrayXb>& mask, mstream& sout){

  int lambda = dt_thr->pos_score ? 1 : -1; // if score is negative, adjust K and K''
  double kval, k2val, wval, vval, rval;
  normal nd(0,1);

  kval = dt_thr->fastSPA ? compute_K_fast_snp(lambda * root, m_ests, dt_thr, mask, ph) : compute_K_snp(lambda * root, m_ests, dt_thr, mask, ph);
  k2val = dt_thr->fastSPA ? compute_K2_fast_snp(lambda * root, m_ests, dt_thr, mask, ph) : compute_K2_snp(lambda * root, m_ests, dt_thr, mask, ph);

  wval = sqrt( 2 * ( root * tval - kval ) );
  vval = root * sqrt( k2val );
  if(vval == 0) {
    dt_thr->chisq_val(ph) = 0;
    dt_thr->pval_log(ph) = 0;
    return;
  }

  rval = wval + log( vval / wval ) / wval;

  if(rval < 0) { // SPA can fail for SNPs with very low counts and give p-value>1
    block_info->test_fail(ph) = true;
    if(params->verbose) sout << "WARNING: SPA correction failed (resulted in p-value > 1).\n";
    return;
  }

  dt_thr->chisq_val(ph) = pow(rval, 2);
  get_logp(dt_thr->pval_log(ph), dt_thr->chisq_val(ph));
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Functions for sum stats output
/////////////////////////////////////////////////
/////////////////////////////////////////////////

//// header line
std::string print_header_output(struct param const* params){

  if(params->split_by_pheno)
    return print_header_output_single(params);
  else
    return print_header_output_all(params);

}

std::string print_header_output_all(struct param const* params){

  int i;
  std::ostringstream buffer;

  buffer << "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ " << 
    ( params->af_cc ? "A1FREQ_CASES A1FREQ_CONTROLS ":"") <<
    ( !params->build_mask && params->dosage_mode ? "INFO ":"") 
    << "N TEST";

  for(i = 1; i < params->n_pheno; i++) 
    buffer << " BETA.Y" << i << " SE.Y" << i << " CHISQ.Y" << i << " LOG10P.Y" << i;
  // last phenotype
  buffer << " BETA.Y" << i << " SE.Y" << i <<  " CHISQ.Y" << i << " LOG10P.Y" << i <<  " EXTRA\n";

  return buffer.str();
}

std::string print_header_output_single(struct param const* params){

  std::ostringstream buffer;

  buffer << "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ " << 
    ( params->af_cc ? "A1FREQ_CASES A1FREQ_CONTROLS ":"") <<
    ( !params->build_mask && params->dosage_mode ? "INFO ":"") << 
    "N TEST BETA SE CHISQ LOG10P EXTRA\n";

  return buffer.str();
}

std::string print_header_output_htp(){

  std::ostringstream buffer;

  buffer << "Name" << "\t" << "Chr" << "\t" << "Pos" << "\t" << "Ref" << "\t" << "Alt" << "\t" << "Trait" << "\t" << "Cohort" << "\t" << "Model" << "\t" << "Effect" << "\t" << "LCI_Effect" << "\t" << "UCI_Effect" << "\t" << "Pval" << "\t" << "AAF" << "\t" << "Num_Cases"<< "\t" << "Cases_Ref" << "\t" << "Cases_Het" << "\t" << "Cases_Alt" << "\t" << "Num_Controls" << "\t" << "Controls_Ref" << "\t" << "Controls_Het"<< "\t"<< "Controls_Alt" << "\t" << "Info\n";

  return buffer.str();
}

//// header info for each snp
std::string print_sum_stats_head(const int& snp_count, vector<snp> const& snpinfo){

  std::ostringstream buffer;

  buffer << snpinfo[snp_count].chrom << " " << snpinfo[snp_count].physpos << " "<< snpinfo[snp_count].ID << " "<< snpinfo[snp_count].allele1 << " "<< snpinfo[snp_count].allele2 << " " ;

  return buffer.str();
}

std::string print_sum_stats_head_htp(const int& snp_count, const string& pheno_name, const string& model, vector<snp> const& snpinfo, struct param const* params){

  std::ostringstream buffer;

  buffer << snpinfo[snp_count].ID << "\t"<< snpinfo[snp_count].chrom << "\t" << snpinfo[snp_count].physpos << "\t"<< snpinfo[snp_count].allele1 << "\t"<< snpinfo[snp_count].allele2 << "\t" << pheno_name << "\t" << params->cohort_name << "\t" << model << "\t";

  return buffer.str();
}


// print sum stats row per snp
std::string  print_sum_stats_line(int const& snp_index, int const& i, string const& tmpstr, string const& test_string, string const& model_type, variant_block* block_info, data_thread* dt_thr, vector<snp> const& snpinfo, struct in_files const& files, struct param const& params){

  std::ostringstream buffer;

  if(params.htp_out) 
    buffer <<  print_sum_stats_head_htp(snp_index, files.pheno_names[i], model_type, snpinfo, &params) << print_sum_stats_htp(dt_thr->bhat(i), dt_thr->se_b(i), dt_thr->chisq_val(i), dt_thr->pval_log(i), block_info->af(i), block_info->info(i), block_info->mac(i), block_info->genocounts, i, !block_info->test_fail(i), 1, &params);
  else  
    buffer << (!params.split_by_pheno && (i>0) ? "" : tmpstr) << print_sum_stats((params.split_by_pheno ? block_info->af(i) : block_info->af1),block_info->af_case(i),block_info->af_control(i), (params.split_by_pheno ? block_info->info(i) : block_info->info1), (params.split_by_pheno ? block_info->ns(i) : block_info->ns1), test_string, dt_thr->bhat(i), dt_thr->se_b(i), dt_thr->chisq_val(i), dt_thr->pval_log(i), !block_info->test_fail(i), 1, &params, (i+1));

  return buffer.str();
}


//// test info for each snp
std::string print_sum_stats(const double& af, const double& af_case, const double& af_control, const double& info, const int& n, const string& model, const double& beta, const double& se, const double& chisq, const double& pv, const bool& test_pass, const int& df, struct param const* params, int const& ipheno){

  if(params->split_by_pheno)
    return print_sum_stats_single(af, af_case, af_control, info, n, model, beta, se, chisq, pv, test_pass, df, params);
  else
    return print_sum_stats_all(af, af_case, af_control, info, n, model, beta, se, chisq, pv, test_pass, df, params, ipheno);
}

// native format - all phenos
std::string print_sum_stats_all(const double& af, const double& af_case, const double& af_control, const double& info, const int& n, const string& model, const double& beta, const double& se, const double& chisq, const double& pv, const bool& test_pass, const int& df, struct param const* params, int const& ipheno){

  std::ostringstream buffer;

  // AF N INFO TEST
  if(ipheno == 1) {
    buffer << af ;
    if( params->af_cc ) buffer << " " << af_case << " " << af_control;
    if(!params->build_mask && params->dosage_mode) buffer << " " << info;
    buffer << " " << n << " " << model ;
  }

  // BETA SE
  if( se < 0 ) buffer << " NA NA";
  else buffer << ' ' << beta << ' ' << se;

  // CHISQ PV
  if((chisq < 0) || !test_pass) 
    buffer << " NA NA";
  else 
    buffer << ' ' << chisq << ' ' << pv;

  // extra column
  if(ipheno == params->n_pheno) {
    if(params->joint_test) buffer << " DF=" << df << endl;
    else buffer << " NA\n";
  }

  return buffer.str();
}

// native format - single pheno
std::string print_sum_stats_single(const double& af, const double& af_case, const double& af_control, const double& info, const int& n, const string& model, const double& beta, const double& se, const double& chisq, const double& pv, const bool& test_pass, const int& df, struct param const* params){

  std::ostringstream buffer;

  // AF N INFO TEST
  buffer << af << " " ;
  if( params->af_cc ) buffer << " " << af_case << " " << af_control;
  if(!params->build_mask && params->dosage_mode) buffer << info << " ";
  buffer << n << " " << model << " " ;

  // BETA SE
  if( se < 0 ) buffer << "NA NA";
  else buffer << beta << ' ' << se;

  // CHISQ PV
  if((chisq < 0) || !test_pass) 
    buffer << " NA NA";
  else 
    buffer << ' ' << chisq << ' ' << pv;

  // extra column
  vector<string> extraCol;
  if(!test_pass) extraCol.push_back("TEST_FAIL");
  if(params->joint_test) extraCol.push_back("DF=" + to_string(df));
  buffer << " " << (extraCol.size() > 0 ? print_scsv(extraCol) : "NA") << endl;

  return buffer.str();
}


std::string print_sum_stats_htp(const double& beta, const double& se, const double& chisq, const double& lpv, const double& af, const double& info, const double& mac, const Ref<const MatrixXd>& genocounts, const int& ph, const bool& test_pass, const int& df, struct param const* params) {

  std::ostringstream buffer;
  bool print_beta = test_pass && (se>=0);
  bool print_pv = test_pass && (chisq>=0);
  double outp_val = -1, effect_val, outse_val;


  if(print_pv) {
    outp_val = max(params->nl_dbl_dmin, pow(10, -lpv)); // to prevent overflow
    if(outp_val == 1) outp_val = 1 - 1e-7;
  } 

  // Effect / CI bounds / Pvalue columns
  if(print_pv && !print_beta)
    buffer << "NA\tNA\tNA\t" << outp_val << "\t";
  else if(!params->binary_mode || (params->firth && test_pass) ){ // qt or firth

    if(!params->binary_mode) // QT
      buffer << beta << "\t" << (beta - params->zcrit * se) << "\t" << (beta + params->zcrit * se) << "\t";
    else // BT (on OR scale)
      buffer << exp(beta) << "\t" << exp(beta - params->zcrit * se) << "\t" << exp(beta + params->zcrit * se) << "\t"; 

    if(print_pv) buffer << outp_val << "\t";
    else buffer << "NA\t";

  } else { // spa/logistic

    if(print_pv) { // for spa or uncorrected logistic score test
      // compute allelic OR
      effect_val = (2*genocounts(3,ph)+genocounts(4,ph)+.5)*(2*genocounts(2,ph)+genocounts(1,ph)+.5)/(2*genocounts(5,ph)+genocounts(4,ph)+.5)/(2*genocounts(0,ph)+genocounts(1,ph)+.5);
      // compute SE = log(allelic OR) / zstat
      outse_val = fabs(log(effect_val)) / sqrt(chisq);
      buffer << effect_val << "\t" << effect_val * exp(- params->zcrit * outse_val) << "\t" << effect_val * exp(params->zcrit * outse_val) << "\t" << outp_val << "\t";
    } else if(!print_beta) 
      buffer << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t";
    else // used in interaction tests
      buffer << exp(beta) << "\t" << exp(beta - params->zcrit * se) << "\t" << exp(beta + params->zcrit * se) << "\tNA\t"; 

  }

  // print out AF
  buffer << af << "\t";

  // print counts in cases
  buffer << (int) genocounts.block(0,ph,3,1).sum() << "\t" << (int) genocounts(0,ph) << "\t" << (int) genocounts(1,ph) << "\t" << (int) genocounts(2,ph) << "\t";
  // print counts in controls
  if(params->binary_mode){
    buffer << (int) genocounts.block(3,ph,3,1).sum() << "\t" << (int) genocounts(3,ph) << "\t" << (int) genocounts(4,ph) << "\t" << (int) genocounts(5,ph);
  } else buffer << "NA\tNA\tNA\tNA";

  // info column
  vector<string> infoCol;
  if(print_beta){
    if(params->binary_mode && test_pass){
      infoCol.push_back( "REGENIE_BETA=" + to_string(beta) );
      infoCol.push_back( "REGENIE_SE=" + to_string(se) );
      // SPA/uncorrected logistic => also print SE from allelic OR
      if(print_pv && !params->firth) infoCol.push_back( "SE=" + to_string(outse_val) );
    } else if(params->binary_mode){
      infoCol.push_back( "REGENIE_BETA=NA" );
      infoCol.push_back( "REGENIE_SE=NA");
      // SPA/uncorrected logistic => also print SE from allelic OR
      if(print_pv && !params->firth) infoCol.push_back( "SE=" + to_string(outse_val) );
    } else infoCol.push_back( "REGENIE_SE=" + to_string(se) );// fot QTs
  }
  // info score
  if(!params->build_mask && params->dosage_mode) infoCol.push_back( "INFO=" + to_string(info) );
  // mac
  infoCol.push_back( "MAC=" + to_string(mac) );
  // df
  if(params->joint_test) infoCol.push_back("DF=" + to_string(df));
  // print info column
  buffer << "\t" << print_scsv(infoCol) << endl;

  return buffer.str();
}



//// print summary of step 2 run
std::string print_summary(Files* ofile, string const& out, std::vector<std::shared_ptr<Files>>& ofile_split, std::vector< string >const& out_split, int const& n_corrected, struct tally const& snp_tally, struct in_files const& files, struct f_ests& fest, struct param const& params){

  std::ostringstream buffer;

  if(!params.skip_test) {

    if(!params.split_by_pheno){
      buffer << "\nAssociation results stored in file : " << out << endl;
      buffer << " + dictionary with trait names in file : " << files.out_file << ".regenie.Ydict\n";
      ofile->closeFile();
    } else {
      buffer << "\nAssociation results stored separately for each trait " << ( params.htp_out ? "(HTPv4 format) " : "" ) << "in files : \n";
      for( int j = 0; j < params.n_pheno; ++j ) {
        buffer << "* [" << out_split[j] << "]\n";
        ofile_split[j]->closeFile();
      }
    }
    buffer << endl;

    if(params.firth || params.use_SPA) {
      buffer << "Number of tests with " << (params.firth ? "Firth " : "SPA ");
      buffer << "correction : " << n_corrected <<  endl;
      buffer << "Number of failed tests : (" << snp_tally.n_failed_tests << "/" << n_corrected << ")\n";
    }

  }

  buffer << "Number of ignored tests due to low MAC ";
  if( params.setMinINFO ) buffer << "or info score ";
  buffer << ": " << snp_tally.n_ignored_snps * params.n_tests_per_variant * params.n_pheno + snp_tally.n_ignored_tests * params.n_tests_per_variant << endl;

  if(params.write_masks)
    buffer << "\nMasks written to : [" << files.out_file << "_masks.{bed,bim,fam}]\n";

  if(params.write_null_firth){ // store file names with null ests
    buffer << "List of files with null Firth estimates written to: [" 
      << print_null_firth_info(files, fest, params) << "]\n";
  }

  return buffer.str();

}
