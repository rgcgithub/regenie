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
#include "cox_score.hpp"
#include "cox_firth.hpp"
#include "Step1_Models.hpp"
#include "Step2_Models.hpp"
#include "HLM.hpp"
#include "Pheno.hpp"
#include "MultiTrait_Tests.hpp"
#include "Ordinal.hpp"
#include "Masks.hpp"
#include "Data.hpp"
#include "MCC.hpp"

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
    if( !params.pheno_pass(ph) ) continue;

    filename = files.blup_files[ files.pheno_names[ph] ];
    ArrayXb read_indiv = ArrayXb::Constant(params.n_samples, false);
    blupf.openForRead(filename, sout);

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
      throw "blup file for phenotype '" + files.pheno_names[ph] + "' has different number of entries on line " + to_string( chrom + 1 ) + " compared to the header (=" + to_string( tmp_str_vec.size() ) + " vs " + to_string( id_strings.size() ) + ").";

    // check starts with chromosome number
    if(chrStrToInt(tmp_str_vec[0], params.nChrom) != chrom) 
      throw "blup file for phenotype '" + files.pheno_names[ph] + "' starts with `" +  tmp_str_vec[0]  + "`"
        + "instead of chromosome number=" + to_string( chrom ) + ".";

    // read blup data
    for( size_t filecol = 1; filecol < id_strings.size(); filecol++ ) {

      // ignore sample if it is not in genotype data
      if (!in_map(id_strings[filecol], params.FID_IID_to_ind)) continue;
      indiv_index = params.FID_IID_to_ind[id_strings[filecol]];

      // ignore sample if it is not included in analysis
      if(!filters.ind_in_analysis(indiv_index)) continue;

      // ignore sample if it is masked for the trait (prs will be 0)
      if(!pheno_data.masked_indivs(indiv_index,ph)) continue;

      // check if duplicate
      if( !read_indiv(indiv_index) )
        read_indiv(indiv_index) = true;
      else 
        throw "individual appears more than once in blup file [" + filename + "]: FID_IID=" + id_strings[filecol];

      in_blup = convertDouble( tmp_str_vec[filecol], &params, sout);

      // if blup is NA then individual must be ignored in analysis for the phenotype (ie mask = 0)
      if (in_blup == params.missing_value_double)
        throw "individual has missing predictions (FID_IID=" + id_strings[filecol] + ";chr=" + to_string( chrom ) + ";phenotype=" + files.pheno_names[ph] + 
          "). Either ignore these individuals using option '--remove', or skip reading predictions with option '--ignore-pred'.\n" + params.err_help ;
      else if(params.w_ltco && (chrom != params.ltco_chr)) // use ltco
        m_ests.blups(indiv_index, ph) = in_blup - m_ests.ltco_prs(indiv_index, ph);
      else // loco
        m_ests.blups(indiv_index, ph) = in_blup;
    }

    // force all non-masked samples to have loco predictions
    //   -> this should not occur as masking of absent samples is done in blup_read() function
    if( (pheno_data.masked_indivs.col(ph).array() && read_indiv).count() < pheno_data.masked_indivs.col(ph).count() )
      throw "all samples included in the analysis (for phenotype " +
        files.pheno_names[ph] + ") must have LOCO predictions in file : " + filename;

    //cerr << m_ests.blups.col(ph).head(5)<<endl;

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

  if(params.trait_mode)
    throw "not for nonQTs";//compute_score_bt();
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

  if(params.trait_mode==1)
    compute_score_bt(isnp, snp_index, chrom, thread_num, test_string, model_type, yres, params, pheno_data, gblock, block_info, snpinfo, m_ests, fest, files, sout);
  else if(params.trait_mode==2)
    compute_score_ct(isnp, snp_index, chrom, thread_num, test_string, model_type, yres, params, pheno_data, gblock, block_info, snpinfo, m_ests, fest, files, sout);
  else if(params.trait_mode==3)
    compute_score_cox(isnp, snp_index, chrom, thread_num, test_string, model_type, params, pheno_data, gblock, block_info, snpinfo, m_ests, fest, files, sout);
  else if(params.trait_mode==0) {
    if(params.mcc_test) {
      compute_score_qt_mcc(isnp, snp_index, thread_num, test_string, model_type, yres, p_sd_yres, params, pheno_data, gblock, block_info, snpinfo, files, sout);
    } else {
      compute_score_qt(isnp, snp_index, thread_num, test_string, model_type, yres, p_sd_yres, params, pheno_data, gblock, block_info, snpinfo, files, sout);
    }
  }
}

// MCC test stat for QT 
void compute_score_qt_mcc(int const& isnp, int const& snp_index, int const& thread_num, string const& test_string, string const& model_type, const Ref<const MatrixXd>& yres, const Ref<const RowVectorXd>& p_sd_yres, struct param const& params, struct phenodt& pheno_data, struct geno_block& gblock, variant_block* block_info, vector<snp> const& snpinfo, struct in_files const& files, mstream& sout){

  double gsc = block_info->flipped ? (4 * params.n_samples + block_info->scale_fac) : block_info->scale_fac;
  string tmpstr; // for sum stats
  MapArXd Geno (gblock.Gmat.col(isnp).data(), params.n_samples, 1);
  data_thread* dt_thr = &(gblock.thread_data[thread_num]);

  if( params.strict_mode ) {
    double n_sq = sqrt( params.n_analyzed - params.ncov_analyzed );
    if(params.skip_blups && dt_thr->is_sparse) // Gsparse is on raw scale (must have yres centered)
      dt_thr->stats = (yres.transpose() * dt_thr->Gsparse.cwiseProduct(pheno_data.masked_indivs.col(0).cast<double>()) / gsc) / n_sq;
    else
      dt_thr->stats = (yres.transpose() * (Geno * pheno_data.masked_indivs.col(0).cast<double>().array()).matrix()) / n_sq;

    if(params.htp_out)
      dt_thr->scores = dt_thr->stats * n_sq * gsc;

    // estimate
    dt_thr->bhat = dt_thr->stats * ( pheno_data.scale_Y.array() * p_sd_yres.array()).matrix().transpose().array() / ( n_sq * gsc );
  } else {
    // compute GtG for each phenotype (different missing patterns)
    dt_thr->scale_fac_pheno = pheno_data.masked_indivs.transpose().cast<double>() * Geno.square().matrix();
    dt_thr->stats = (yres.transpose() * Geno.matrix()).array() / dt_thr->scale_fac_pheno.sqrt();

    if(params.htp_out)
      dt_thr->scores = dt_thr->stats * dt_thr->scale_fac_pheno.sqrt() * gsc;

    // estimate
    dt_thr->bhat = dt_thr->stats * ( pheno_data.scale_Y.array() * p_sd_yres.array() ).matrix().transpose().array() / ( sqrt(dt_thr->scale_fac_pheno) * gsc );
  }

  // SE
  dt_thr->se_b = dt_thr->bhat / dt_thr->stats;

  // get test statistic
  dt_thr->chisq_val = dt_thr->stats.square();

  // (1) MCC if mcc_apply_thr == false; (2) Score -> MCC if Pval(Score) < mcc_thr
  MCC mcc;
  boost::math::chi_squared chisq(1);
  double chisq_val_adj;

  if(!params.mcc_apply_thr) {
    // (1) only MCC
    mcc.setup_y(pheno_data.masked_indivs, yres, params.ncov_analyzed);
    MCCResults mcc_results = mcc.run(Geno);
    // store MCC results into dt_thr
    for( int i = 0; i < params.n_pheno; ++i ) {
      if(mcc_results.Skip(i, 0)) {
        dt_thr->pval_log(i) = -1;
        block_info->test_fail(i) = true;
      } else {
        dt_thr->pval_log(i) = -log10(mcc_results.Pval(i, 0));
        // adjust SE
        chisq_val_adj = boost::math::quantile(boost::math::complement(chisq, mcc_results.Pval(i, 0)));
        dt_thr->se_b(i) *= sqrt(dt_thr->chisq_val(i) / chisq_val_adj);
      }
    }
  } else {
    // (2) Score -> MCC
    for( int i = 0; i < params.n_pheno; ++i ) {
      get_logp(dt_thr->pval_log(i), dt_thr->chisq_val(i));
      // check for skewness of phenotype i
      if(dt_thr->pval_log(i) > params.mcc_thr_nlog10 && pheno_data.mcc_Y[i]) {
        mcc.setup_y(pheno_data.masked_indivs.col(i), yres.col(i), params.ncov_analyzed);
        MCCResults mcc_results_i = mcc.run(Geno);
        if(mcc_results_i.Skip(0, 0)) {
          dt_thr->pval_log(i) = -1;
          block_info->test_fail(i) = true;
        } else {
          dt_thr->pval_log(i) = -log10(mcc_results_i.Pval(0, 0));
          // adjust SE
          chisq_val_adj = boost::math::quantile(boost::math::complement(chisq, mcc_results_i.Pval(0, 0)));
          dt_thr->se_b(i) *= sqrt(dt_thr->chisq_val(i) / chisq_val_adj);
        }
      }
    }
  }

  if(!params.htp_out) tmpstr = print_sum_stats_head(snp_index, snpinfo);

  for( int i = 0; i < params.n_pheno; ++i ) {

    if( !params.pheno_pass(i) || block_info->ignored_trait(i) ) {
      if(!params.p_joint_only && !params.split_by_pheno)
        block_info->sum_stats[i].append( print_na_sumstats(i, 1, tmpstr, test_string, block_info, params) );
      continue;
    }
    if(block_info->flipped) dt_thr->bhat(i) *= -1;

    // get MCC pvalue
    /* get_logp(dt_thr->pval_log(i), dt_thr->chisq_val(i)); */
    /* if(mcc_results.Skip(i, 0)) { */
    /*   dt_thr->pval_log(i) = -1; */
    /*   block_info->test_fail(i) = true; */
    /* } else { */
    /*   dt_thr->pval_log(i) = -log10(mcc_results.Pval(i, 0)); */
    /* } */

    if(!params.p_joint_only)
      block_info->sum_stats[i].append( print_sum_stats_line(snp_index, i, tmpstr, test_string, model_type, block_info, dt_thr, snpinfo, files, params) );

  }

}
// score test stat for QT
void compute_score_qt(int const& isnp, int const& snp_index, int const& thread_num, string const& test_string, string const& model_type, const Ref<const MatrixXd>& yres, const Ref<const RowVectorXd>& p_sd_yres, struct param const& params, struct phenodt& pheno_data, struct geno_block& gblock, variant_block* block_info, vector<snp> const& snpinfo, struct in_files const& files, mstream& sout){

  double gsc = block_info->flipped ? (4 * params.n_samples + block_info->scale_fac) : block_info->scale_fac;
  string tmpstr; // for sum stats
  MapArXd Geno (gblock.Gmat.col(isnp).data(), params.n_samples, 1);
  data_thread* dt_thr = &(gblock.thread_data[thread_num]);
  double n_sq = sqrt( params.n_analyzed - params.ncov_analyzed );

  if( params.strict_mode ) {

    if(params.skip_blups && dt_thr->is_sparse) // Gsparse is on raw scale (must have yres centered)
      dt_thr->stats = (yres.transpose() * dt_thr->Gsparse.cwiseProduct(pheno_data.masked_indivs.col(0).cast<double>()) / gsc) / n_sq;
    else
      dt_thr->stats = (yres.transpose() * (Geno * pheno_data.masked_indivs.col(0).cast<double>().array()).matrix()) / n_sq;

    if(params.htp_out) {
      dt_thr->scores = dt_thr->stats * n_sq * gsc;
      dt_thr->skat_var = (gsc * n_sq)*(gsc * n_sq);
    }

    // estimate
    dt_thr->bhat = dt_thr->stats * ( pheno_data.scale_Y.array() * p_sd_yres.array()).matrix().transpose().array() / ( n_sq * gsc );

  } else {

    // compute GtG for each phenotype (different missing patterns)
    dt_thr->scale_fac_pheno = pheno_data.masked_indivs.transpose().cast<double>() * Geno.square().matrix();
    dt_thr->stats = (yres.transpose() * Geno.matrix()).array() / dt_thr->scale_fac_pheno.sqrt();

    if(params.htp_out) {
      dt_thr->scores = dt_thr->stats * dt_thr->scale_fac_pheno.sqrt() * gsc;
      dt_thr->skat_var = (gsc * n_sq)*(gsc * n_sq);
    }

    // estimate
    dt_thr->bhat = dt_thr->stats * ( pheno_data.scale_Y.array() * p_sd_yres.array() ).matrix().transpose().array() / ( sqrt(dt_thr->scale_fac_pheno) * gsc );

  }

  // SE
  dt_thr->se_b = dt_thr->bhat / dt_thr->stats;

  // get test statistic
  dt_thr->chisq_val = dt_thr->stats.square();

  if(!params.htp_out) tmpstr = print_sum_stats_head(snp_index, snpinfo);

  for( int i = 0; i < params.n_pheno; ++i ) {

    if( !params.pheno_pass(i) || block_info->ignored_trait(i) ) {
      if(!params.p_joint_only && !params.split_by_pheno)
        block_info->sum_stats[i].append( print_na_sumstats(i, 1, tmpstr, test_string, block_info, params) );
      continue;
    }
    if(block_info->flipped) dt_thr->bhat(i) *= -1;

    // get pvalue
    get_logp(dt_thr->pval_log(i), dt_thr->chisq_val(i));

    if(!params.p_joint_only)
      block_info->sum_stats[i].append( print_sum_stats_line(snp_index, i, tmpstr, test_string, model_type, block_info, dt_thr, snpinfo, files, params) );

  }

}

void compute_score_bt(int const& isnp, int const& snp_index, int const& chrom, int const& thread_num, string const& test_string, string const& model_type, const Ref<const MatrixXd>& yres, struct param const& params, struct phenodt& pheno_data, struct geno_block& gblock, variant_block* block_info, vector<snp> const& snpinfo, struct ests const& m_ests, struct f_ests& fest, struct in_files const& files, mstream& sout){

  string tmpstr; 
  VectorXd GW, XtWG;
  SpVec GWs;
  data_thread* dt_thr = &(gblock.thread_data[thread_num]);

  // header snp info for sum stats
  if(!params.htp_out) tmpstr = print_sum_stats_head(snp_index, snpinfo);

  // genotype for marker
  MapArXd Geno (gblock.Gmat.col(isnp).data(), params.n_samples, 1);

  for( int i = 0; i < params.n_pheno; ++i ) {

    if( !params.pheno_pass(i) || block_info->ignored_trait(i) ){
      if(!params.p_joint_only && !params.split_by_pheno)
        block_info->sum_stats[i].append( print_na_sumstats(i, 1, tmpstr, test_string, block_info, params) );
      continue;
    }

    MapArXb mask (pheno_data.masked_indivs.col(i).data(), params.n_samples, 1);
    MapcArXd Wsqrt (m_ests.Gamma_sqrt.col(i).data(), params.n_samples, 1);
    MapcMatXd XWsqrt (m_ests.X_Gamma[i].data(), params.n_samples, m_ests.X_Gamma[i].cols());

    // project out covariates from G
    if(dt_thr->is_sparse) {
      GWs = dt_thr->Gsparse.cwiseProduct( (Wsqrt * mask.cast<double>()).matrix() );
      XtWG = XWsqrt.transpose() * GWs;
      dt_thr->Gres = -XWsqrt * XtWG;
      dt_thr->Gres += GWs;
    } else {
      GW = (Geno * Wsqrt * mask.cast<double>()).matrix();
      dt_thr->Gres = GW - XWsqrt * (XWsqrt.transpose() * GW);
    }

    // denominator
    if(dt_thr->is_sparse) 
      dt_thr->denum(i) = GWs.squaredNorm() - XtWG.squaredNorm();
    else
      dt_thr->denum(i) = dt_thr->Gres.squaredNorm();
    if( dt_thr->denum(i) < params.numtol ){
      block_info->ignored_trait(i) = true;
      if(!params.p_joint_only && !params.split_by_pheno)
        block_info->sum_stats[i].append( print_na_sumstats(i, 1, tmpstr, test_string, block_info, params) );
      continue;
    }
    // score test stat for BT
    if(dt_thr->is_sparse) 
      dt_thr->stats(i) = GWs.dot(yres.col(i)) / sqrt( dt_thr->denum(i) );
    else
      dt_thr->stats(i) = dt_thr->Gres.col(0).dot(yres.col(i)) / sqrt( dt_thr->denum(i) );

    if(params.htp_out) {
      dt_thr->scores(i) = dt_thr->stats(i) * sqrt( dt_thr->denum(i) );
      dt_thr->skat_var(i) = dt_thr->denum(i);
    }

    /*
    if(params.debug) {
      cerr << "\ny:\n" << yres.col(i).topRows(2) << endl;
      cerr << "\nGresid:\n" << dt_thr->Gres.topRows(2) << endl;
      if(dt_thr->is_sparse) cerr << "\nsum(GW)=" << GWs.sum() << endl;
      cerr << "\nscore=" << dt_thr->Gres.col(0).dot(yres.col(i)) << " var(score)=" << dt_thr->Gres.squaredNorm() << endl;
    }
    */

    // use firth/spa
    check_pval_snp(block_info, dt_thr, chrom, i, isnp, pheno_data, gblock, m_ests, fest, params, sout);

    dt_thr->bhat(i) /= block_info->scale_fac;
    dt_thr->se_b(i) /= block_info->scale_fac;
    if(block_info->flipped) dt_thr->bhat(i) *= -1;

    // print sum stats
    if(!params.p_joint_only)
      block_info->sum_stats[i].append( print_sum_stats_line(snp_index, i, tmpstr, test_string, model_type, block_info, dt_thr, snpinfo, files, params) );

  }


}


// poisson
void compute_score_ct(int const& isnp, int const& snp_index, int const& chrom, int const& thread_num, string const& test_string, string const& model_type, const Ref<const MatrixXd>& yres, struct param const& params, struct phenodt& pheno_data, struct geno_block& gblock, variant_block* block_info, vector<snp> const& snpinfo, struct ests const& m_ests, struct f_ests& fest, struct in_files const& files, mstream& sout){

  string tmpstr; 
  MatrixXd GW;
  SpVec GWs;
  data_thread* dt_thr = &(gblock.thread_data[thread_num]);

  // header snp info for sum stats
  if(!params.htp_out) tmpstr = print_sum_stats_head(snp_index, snpinfo);

  // genetype for marker
  MapArXd Geno (gblock.Gmat.col(isnp).data(), params.n_samples, 1);

  for( int i = 0; i < params.n_pheno; ++i ) {

    if( !params.pheno_pass(i) || block_info->ignored_trait(i) ) {
      if(!params.p_joint_only && !params.split_by_pheno)
        block_info->sum_stats[i].append( print_na_sumstats(i, 1, tmpstr, test_string, block_info, params) );
      continue;
    }
    MapArXb mask (pheno_data.masked_indivs.col(i).data(), params.n_samples, 1);
    MapcArXd Wsqrt (m_ests.Gamma_sqrt.col(i).data(), params.n_samples, 1);
    MapcMatXd XWsqrt (m_ests.X_Gamma[i].data(), params.n_samples, m_ests.X_Gamma[i].cols());

    // project out covariates from G
    if(dt_thr->is_sparse) {
      GWs = dt_thr->Gsparse.cwiseProduct( (Wsqrt * mask.cast<double>()).matrix() );
      dt_thr->Gres = -XWsqrt * (XWsqrt.transpose() * GWs);
      dt_thr->Gres += GWs;
    } else {
      GW = (Geno * Wsqrt * mask.cast<double>()).matrix();
      dt_thr->Gres = GW - XWsqrt * (XWsqrt.transpose() * GW);
    }

    // denominator
    dt_thr->denum(i) = dt_thr->Gres.squaredNorm();
    if( dt_thr->denum(i) < params.numtol ){
      block_info->ignored_trait(i) = true;
      if(!params.p_joint_only && !params.split_by_pheno)
        block_info->sum_stats[i].append( print_na_sumstats(i, 1, tmpstr, test_string, block_info, params) );
      continue;
    }
    // score test stat for CT
    if(dt_thr->is_sparse) 
      dt_thr->stats(i) = GWs.dot(yres.col(i)) / sqrt( dt_thr->denum(i) );
    else
      dt_thr->stats(i) = dt_thr->Gres.col(0).dot(yres.col(i)) / sqrt( dt_thr->denum(i) );

    if(params.debug) {
      cerr << "\ny:\n" << yres.col(i).topRows(2) << endl;
      cerr << "\nGresid:\n" << dt_thr->Gres.topRows(2) << endl;
      if(dt_thr->is_sparse) cerr << "\nsum(GW)=" << GWs.sum() << endl;
      cerr << "\nscore=" << dt_thr->Gres.col(0).dot(yres.col(i)) << " var(score)=" << dt_thr->Gres.squaredNorm() << endl;
    }

    // apply correction
    //check_pval_snp(block_info, dt_thr, chrom, i, isnp, pheno_data, gblock, m_ests, fest, params, sout);
    get_sumstats(false, i, dt_thr);

    dt_thr->bhat(i) /= block_info->scale_fac;
    dt_thr->se_b(i) /= block_info->scale_fac;
    if(block_info->flipped) dt_thr->bhat(i) *= -1;

    // print sum stats
    if(!params.p_joint_only)
      block_info->sum_stats[i].append( print_sum_stats_line(snp_index, i, tmpstr, test_string, model_type, block_info, dt_thr, snpinfo, files, params) );

  }


}

void compute_score_cox(int const& isnp, int const& snp_index, int const& chrom, int const& thread_num, string const& test_string, string const& model_type, struct param const& params, struct phenodt& pheno_data, struct geno_block& gblock, variant_block* block_info, vector<snp> const& snpinfo, struct ests const& m_ests, struct f_ests& fest, struct in_files const& files, mstream& sout){

  string tmpstr; 
  data_thread* dt_thr = &(gblock.thread_data[thread_num]);

  Eigen::VectorXd sqrtWG;
  SpVec sqrtWGs;
  Eigen::VectorXd RGammaG;
  Eigen::VectorXd UhalfG;
  Eigen::VectorXd XtWG;
  Eigen::VectorXd XtUG;
  Eigen::VectorXd XtVG;
  double T;

  // header snp info for sum stats
  if(!params.htp_out) tmpstr = print_sum_stats_head(snp_index, snpinfo);

  // genotype for marker
  MapArXd Geno (gblock.Gmat.col(isnp).data(), params.n_samples, 1);
  
  for( int i = 0; i < params.n_pheno; ++i ) {
    if( !params.pheno_pass(i) || block_info->ignored_trait(i) ){
      if(!params.p_joint_only && !params.split_by_pheno)
        block_info->sum_stats[i].append( print_na_sumstats(i, 1, tmpstr, test_string, block_info, params) );
      continue;
    }
    MapArXb mask(pheno_data.masked_indivs.col(i).data(), params.n_samples, 1);

    // score stat
    if (dt_thr->is_sparse) {
      dt_thr->Gres = dt_thr->Gsparse - m_ests.cox_MLE_NULL[i].X1_X1WX1inv * (dt_thr->Gsparse.transpose() * m_ests.cox_MLE_NULL[i].WX1).transpose();
    } else {
      dt_thr->Gres = Geno.matrix() - m_ests.cox_MLE_NULL[i].X1_X1WX1inv * (Geno.matrix().transpose() * m_ests.cox_MLE_NULL[i].WX1).transpose();
    }
    T = (dt_thr->Gres.array() * m_ests.cox_MLE_NULL[i].residual.array() * mask.cast<double>()).sum();
    
    dt_thr->denum(i) = m_ests.cox_MLE_NULL[i].res_var * (dt_thr->Gres.array()).pow(2).sum();
    
    if (params.coxscore_exact) {
      sqrtWG = dt_thr->Gres.array() * (m_ests.cox_MLE_NULL[i].mu.array().sqrt()) * mask.cast<double>();
      RGammaG = cumulativeSum_reverse2( m_ests.survival_data_pheno[i].R.transpose() * (m_ests.cox_MLE_NULL[i].w_exp_eta.array() * (m_ests.survival_data_pheno[i].permute_mtx * dt_thr->Gres).array()).matrix());
      UhalfG = m_ests.cox_MLE_NULL[i].Dhalf.array() * RGammaG.array();

      XtWG = m_ests.cox_MLE_NULL[i].sqrtWX.transpose() * sqrtWG;
      XtUG = m_ests.cox_MLE_NULL[i].UhalfX.transpose() * UhalfG;
      XtVG = XtWG - XtUG;
      dt_thr->denum(i) = sqrtWG.squaredNorm() - UhalfG.squaredNorm() - (XtVG.array() * (m_ests.cox_MLE_NULL[i].cov_inv * XtVG).array()).sum();
    }
    dt_thr->stats(i) = T/sqrt(dt_thr->denum(i));

    if(params.htp_out) {
      dt_thr->scores(i) = dt_thr->stats(i) * sqrt( dt_thr->denum(i) );
      dt_thr->skat_var(i) = dt_thr->denum(i);
    }

    // use firth/spa
    check_pval_snp(block_info, dt_thr, chrom, i, isnp, pheno_data, gblock, m_ests, fest, params, sout);

    dt_thr->bhat(i) /= block_info->scale_fac;
    dt_thr->se_b(i) /= block_info->scale_fac;
    if(block_info->flipped) dt_thr->bhat(i) *= -1;

    // print sum stats
    if(!params.p_joint_only) {
      block_info->sum_stats[i].append( print_sum_stats_line(snp_index, i, tmpstr, test_string, model_type, block_info, dt_thr, snpinfo, files, params) );
    }
  }
}

// Cox Null firth model
void fit_null_firth_cox(bool const& silent, int const& chrom, struct f_ests* firth_est, struct phenodt* pheno_data, struct ests const* m_ests, struct in_files* files, struct param* params, mstream& sout){

  auto t1 = std::chrono::high_resolution_clock::now();
  ArrayXb has_converged = params->pheno_pass; // if null log reg converged
  IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","");

  if(!silent) sout << "   -fitting null Firth cox regression on time-to-event phenotypes..." << flush;

  // fit null firth (in parallel for MT mode)
#if defined(_OPENMP)
  if((params->n_pheno>2) && !params->blup_cov) setNbThreads(1); // for < 3, mt in eigen should be similar
#pragma omp parallel for schedule(dynamic) if((params->n_pheno>2) && !params->blup_cov)
#endif
  for( int i = 0; i < params->n_pheno; ++i ) {
    if( !params->pheno_pass(i) ) continue;
    if(params->blup_cov) // add step 1 predictions as a covariate (skip multithreading)
      pheno_data->new_cov.rightCols(1) = m_ests->blups.col(i);
    
    Eigen::VectorXd offset;
    if(params->blup_cov) offset = Eigen::VectorXd::Zero(m_ests->blups.rows()); // if step 1 is covariate
    else offset = m_ests->blups.col(i).array();

    cox_firth cox_firth_null;
    cox_firth_null.setup(m_ests->survival_data_pheno[i], pheno_data->new_cov, offset, pheno_data->new_cov.cols(), params->niter_max_firth_null, params->niter_max_line_search, params->numtol_cox, params->numtol_beta_cox, params->maxstep_null, !params->cox_nofirth, false, m_ests->cox_MLE_NULL[i].beta);
    cox_firth_null.fit(m_ests->survival_data_pheno[i], pheno_data->new_cov, offset);
    if( !cox_firth_null.converge ){ // if failed to converge
      cerr << "WARNING: Cox regression with Firth correction did not converge (maximum step size=" << params->maxstep_null <<";maximum number of iterations=" << params->niter_max_firth_null <<").";

      cerr << "Retrying with fallback parameters: (maximum step size=" << params->maxstep_null/5 <<";maximum number of iterations=" << params->niter_max_firth_null*5 <<").\n";

      cox_firth_null.setup(m_ests->survival_data_pheno[i], pheno_data->new_cov, offset, pheno_data->new_cov.cols(), params->niter_max_firth_null*5, params->niter_max_line_search, params->numtol_cox, params->numtol_beta_cox, params->maxstep_null/5, !params->cox_nofirth, false);
      cox_firth_null.fit(m_ests->survival_data_pheno[i], pheno_data->new_cov, offset);
    }
    has_converged(i) = cox_firth_null.converge;
    if(!has_converged(i)) continue;
    firth_est->cov_blup_offset.col(i) = cox_firth_null.eta;
    firth_est->beta_null_firth.col(i) = cox_firth_null.beta;
    
    if(params->write_null_firth)
      (*firth_est->firth_est_files[i]) << chrom << " " << cox_firth_null.beta.transpose().format(Fmt) << endl;

  }
#if defined(_OPENMP)
  if((params->n_pheno>2) && !params->blup_cov) setNbThreads(params->threads);
#endif

  // check if some did not converge
  if(!has_converged.any()) { //  none passed

    string msg1 = to_string( params->maxstep_null / 5 );
    string msg2 = to_string( params->niter_max_firth_null * 5 );
    throw "Firth penalized cox regression failed to converge for all phenotypes."
      " Try decreasing the maximum step size using `--maxstep-null` (currently=" + msg1 +  ") "
      "and increasing the maximum number of iterations using `--maxiter-null` (currently=" + msg2 + ").";

  } else if( ((!has_converged) && (params->pheno_pass || params->pheno_fail_nullreg)).any() ) { // some phenotypes failed (at null reg or null firth) - write their names to file

    ArrayXb pheno_flagged = (!has_converged) && (params->pheno_pass || params->pheno_fail_nullreg);
    Files outf;
    string failed_file = files->out_file + "_failedNullFirth_chr" + to_string(chrom) + ".list";
    outf.openForWrite( failed_file, sout);
    for( int i = 0; i < params->n_pheno; ++i )
      if(pheno_flagged(i))
        outf << files->pheno_names[i] << endl;
    outf.closeFile();
    sout << "WARNING: null Firth failed for " << pheno_flagged.count() << " phenotypes (list of traits written to '" << failed_file << "' and these will be skipped)\n";
    params->pheno_pass = has_converged;

  }
  if(params->blup_cov)
    pheno_data->new_cov.rightCols(1).array() = 0;

  if (silent) return;

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;
}

// firth cov + test snp
void fit_firth_cox_snp(int const& chrom, int const& ph, int const& isnp, struct param const* params, struct phenodt* pheno_data, struct ests const* m_ests, struct f_ests const* fest, const Ref<const MatrixXd>& Gvec, variant_block* block_info, data_thread* dt_thr, mstream& sout) {
  // if firth is used, fit based on penalized log-likelihood
  int col_incl;
  double lrt;
  VectorXd beta0;
  MatrixXd Xmat;

  Xmat = MatrixXd::Zero(params->n_samples, pheno_data->new_cov.cols() + 1); // covariates + tested SNP
  Xmat << pheno_data->new_cov, Gvec;
  col_incl = Xmat.cols();

  // null model
  cox_firth cox_firth_null;
  cox_firth_null.setup(m_ests->survival_data_pheno[ph], Xmat, m_ests->blups.col(ph), col_incl-1, params->niter_max_firth, params->niter_max_line_search, params->numtol_cox, params->numtol_beta_cox, params->maxstep_null, !params->cox_nofirth, false);
  cox_firth_null.fit(m_ests->survival_data_pheno[ph], Xmat, m_ests->blups.col(ph));

  // if (!cox_firth_null.converge) {
  //   cox_firth_null.setup(m_ests->survival_data_pheno[ph], Xmat, m_ests->blups.col(ph), col_incl-1, params->niter_max_firth*5, params->niter_max_line_search_cox, params->numtol_cox, params->maxstep, !params->cox_nofirth, false);
  //   cox_firth_null.fit(m_ests->survival_data_pheno[ph], Xmat, m_ests->blups.col(ph));
  // }

  if (!cox_firth_null.converge) {
    if(params->verbose) cerr << "WARNING: Cox regression with Firth correction null model did not converge!\n";
      block_info->test_fail(ph) = true;
      return ;
  }

  // test model
  cox_firth cox_firth_test;
  cox_firth_test.setup(m_ests->survival_data_pheno[ph], Xmat, m_ests->blups.col(ph), col_incl, params->niter_max_firth, params->niter_max_line_search, params->numtol_cox, params->numtol_beta_cox, params->maxstep, !params->cox_nofirth, false, cox_firth_null.beta);
  cox_firth_test.fit(m_ests->survival_data_pheno[ph], Xmat, m_ests->blups.col(ph));

  // if (!cox_firth_test.converge) {
  //   cox_firth_test.setup(m_ests->survival_data_pheno[ph], Xmat, m_ests->blups.col(ph), col_incl, params->niter_max_firth*5, params->niter_max_line_search_cox, params->numtol_cox, params->maxstep, true, false);
  //   cox_firth_test.fit(m_ests->survival_data_pheno[ph], Xmat, m_ests->blups.col(ph));
  // }

  if (!cox_firth_test.converge) {
    if(params->verbose) cerr << "WARNING: Cox regression with Firth correction did not converge!\n";
      block_info->test_fail(ph) = true;
      return ;
  }

  dt_thr->bhat(ph) = cox_firth_test.beta.tail(1)(0);
  if(!params->back_correct_se)
    dt_thr->se_b(ph) = cox_firth_test.qrsd.inverse().diagonal().array().sqrt().tail(1)(0);

  lrt = 2*(cox_firth_test.loglike.tail(1)(0) - cox_firth_null.loglike.tail(1)(0));
  if( lrt < 0 ) {
    block_info->test_fail(ph) = true;
    return ;
  }
  dt_thr->dif_deviance = lrt;
  return ;
}

// for approx firth testing step
void fit_firth_cox_snp_fast(int const& chrom, int const& ph, int const& isnp, struct param const* params, struct phenodt* pheno_data, struct ests const* m_ests, struct f_ests const* fest, const Ref<const VectorXd>& Gvec, variant_block* block_info, data_thread* dt_thr, mstream& sout) {
  // if firth is used, fit based on penalized log-likelihood
  double lrt;

  // // For rare variants, set entries in Gvec for non-carriers to 0
  // int mac_thr_sparse = (params->skip_fast_firth ? 0 : 50), i = 0, index_j;
  // ArrayXi index_carriers;
  // if(dt_thr->is_sparse && (block_info->mac(ph) < mac_thr_sparse)) {
  //   index_carriers.resize(dt_thr->Gsparse.nonZeros());
  //   for (SpVec::InnerIterator it(dt_thr->Gsparse); it; ++it) {
  //     index_j = it.index();
  //     // check for small entries in G (eg with imputed data)
  //     if(mask(index_j) && (it.value() > 1e-4)) index_carriers(i++) = index_j;
  //   }
  //   index_carriers.conservativeResize(i);
  // }

  cox_firth cox_firth_test;
  cox_firth_test.setup(m_ests->survival_data_pheno[ph], Gvec, fest->cov_blup_offset.col(ph), 1, params->niter_max_firth, params->niter_max_line_search, params->numtol_cox, params->numtol_beta_cox, params->maxstep, !params->cox_nofirth, false);
  cox_firth_test.fit_1(m_ests->survival_data_pheno[ph], Gvec, fest->cov_blup_offset.col(ph));

  // if(!cox_firth_test.converge){
  //   cox_firth_test.setup(m_ests->survival_data_pheno[ph], Gvec, fest->cov_blup_offset.col(ph), 1, params->niter_max_firth*5, params->niter_max_line_search_cox, params->numtol_cox*10, params->maxstep, true, false);
  //   cox_firth_test.fit_1(m_ests->survival_data_pheno[ph], Gvec, fest->cov_blup_offset.col(ph));
  // }
  if(!cox_firth_test.converge){
    if(params->verbose) cerr << "WARNING: Cox regression with Firth correction did not converge!\n";
    block_info->test_fail(ph) = true;
    return ;
  }

  // compute beta_hat
  dt_thr->bhat(ph) = cox_firth_test.beta(0);
  // compute SE based on Hessian for unpenalized LL
  if(!params->back_correct_se)
    dt_thr->se_b(ph) = sqrt(1/cox_firth_test.second_der_1);

  lrt = 2*(cox_firth_test.loglike.tail(1)(0) - cox_firth_test.loglike(0));
  if( lrt < 0 ) {
    block_info->test_fail(ph) = true;
    return ;
  }
  dt_thr->dif_deviance = lrt;
  return ;
}


// Firth (currently only used for null approximate firth)
bool fit_approx_firth_null(int const& chrom, int const& ph, struct phenodt const* pheno_data, struct ests const* m_ests, Ref<ArrayXd> betavec, struct param* params, bool const& save_se) {

  bool success, set_start = true;
  int col_incl;
  int maxstep = params->maxstep_null;
  int niter = params->niter_max_firth_null;
  double tol = 50*params->numtol;
  double dev, lrt;

  ArrayXd betaold, se, etavec, pivec, offset;

  MapcArXd Y (pheno_data->phenotypes_raw.col(ph).data(), pheno_data->phenotypes_raw.rows());
  MapcMatXd Xmat (pheno_data->new_cov.data(), pheno_data->new_cov.rows(), pheno_data->new_cov.cols());
  MapcArXb mask (pheno_data->masked_indivs.col(ph).data(), pheno_data->masked_indivs.rows());
  col_incl = Xmat.cols();

  if(params->blup_cov) offset = ArrayXd::Zero(m_ests->blups.rows()); // if step 1 is covariate
  else offset = m_ests->blups.col(ph).array();

  // with firth approx. => trial 1: use maxstep_null
  // trial=1+ => start at 0 (update maxstep & niter)
  // trial=2+ => use fallback options (update maxstep & niter)
  for(int trial = 0; trial < 3; trial++){

    // starting values
    if(set_start){
        if(params->use_null_firth || (trial == 0) ){ // use saved est or those from unpenalized log. reg
          betaold = betavec.head(Xmat.cols());
        } else {// set to 0 if null firth failed
          betaold = 0;
          //betaold(0) = ( 0.5 + mask.select(Y,0).sum())  / (pheno_data->Neff(ph) + 1);
          //betaold(0) = log( betaold(0) / (1 - betaold(0) ));
          // LOCO prediction is offset
          betaold(0) -= mask.select(offset,0).mean();
        }
    }

    success = fit_firth(ph, Y, Xmat, offset, mask, pivec, etavec, betaold, se, col_incl, dev, false, lrt, maxstep, niter, tol, params);

    if(!params->fix_maxstep_null) { // don't retry with user-given settings
      if( !success ){ // if failed to converge
        cerr << "WARNING: Logistic regression with Firth correction did not converge (maximum step size=" << maxstep <<";maximum number of iterations=" << niter <<").";

        // try fitting pseudo-data representation with IRLS
        double dev0 = 0;
        if(
            fit_firth_pseudo(dev0, Y, Xmat, offset, mask, pivec, etavec, betaold, se, col_incl, dev, false, lrt, maxstep, niter, tol, params)
          ){
          success = true;
          break;
        }

        if( trial == 1 ){
          maxstep /= 5;
          niter *= 5;
          if(params->debug) cerr << "Retrying with fallback parameters: (maximum step size=" << maxstep <<";maximum number of iterations=" << niter<<").\n";
        }
        if(params->use_adam) set_start = false;
        continue;
      }
    }

    break;
  }

  // If didn't converge
  if(!success)
    return false;

  betavec.head(betaold.size()) = betaold;
  if(save_se && params->print_cov_betas) { // get se
    ArrayXd wvec;
    get_wvec(pivec, wvec, mask);
    MatrixXd XWsqrt = ( Xmat.array().colwise() * (wvec.sqrt() * mask.cast<double>()) ).matrix();
    MatrixXd xtx_inv = ( XWsqrt.transpose() * XWsqrt ).colPivHouseholderQr().inverse();
    params->xtx_inv_diag.col(ph).array() = xtx_inv.diagonal().array().sqrt();
  }
  return true;

}

// Approximate null firth model
void fit_null_firth(bool const& silent, int const& chrom, struct f_ests* firth_est, struct phenodt* pheno_data, struct ests const* m_ests, struct in_files* files, struct param* params, mstream& sout){

  auto t1 = std::chrono::high_resolution_clock::now();
  ArrayXb has_converged = params->pheno_pass; // if null log reg converged
  IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","");

  if(!silent && params->firth) sout << "   -fitting null Firth logistic regression on binary phenotypes..." << flush;

  // get starting values
  if(params->use_null_firth) // saved in file
    get_beta_start_firth(chrom, firth_est, files, params, sout);
  else // from null log. reg.
    get_beta_start_firth(firth_est, m_ests);

  // fit null firth (in parallel for MT mode)
#if defined(_OPENMP)
  if((params->n_pheno>2) && !params->blup_cov) setNbThreads(1); // for < 3, mt in eigen should be similar
#pragma omp parallel for schedule(dynamic) if((params->n_pheno>2) && !params->blup_cov)
#endif
  for( int i = 0; i < params->n_pheno; ++i ) {
    if( !params->pheno_pass(i) ) continue;

    if(params->blup_cov) // add step 1 predictions as a covariate (skip multithreading)
      pheno_data->new_cov.rightCols(1) = m_ests->blups.col(i);

    MapArXd bvec (firth_est->beta_null_firth.col(i).data(), firth_est->beta_null_firth.rows());
    has_converged(i) = fit_approx_firth_null(chrom, i, pheno_data, m_ests, bvec, params);
    if(!has_converged(i)) continue; // cannot use break

    if(params->test_mode){
      firth_est->cov_blup_offset.col(i) = pheno_data->new_cov * bvec.head(pheno_data->new_cov.cols()).matrix(); // store offset used for approx firth
     if(!params->blup_cov) firth_est->cov_blup_offset.col(i) += m_ests->blups.col(i); // if offset  
    }

    if(params->write_null_firth)
      (*firth_est->firth_est_files[i]) << chrom << " " << bvec.head(params->ncov).matrix().transpose().format(Fmt) << endl;

  }
#if defined(_OPENMP)
  if((params->n_pheno>2) && !params->blup_cov) setNbThreads(params->threads);
#endif

  // check if some did not converge
  if(!has_converged.any()) { //  none passed

    string msg1 = to_string( params->maxstep_null / (params->fix_maxstep_null ? 1 : 5) );
    string msg2 = to_string( params->niter_max_firth_null * (params->fix_maxstep_null ? 1 : 5) );
    throw "Firth penalized logistic regression failed to converge for all phenotypes."
      " Try decreasing the maximum step size using `--maxstep-null` (currently=" + msg1 +  ") "
      "and increasing the maximum number of iterations using `--maxiter-null` (currently=" + msg2 + ").";

  } else if( ((!has_converged) && (params->pheno_pass || params->pheno_fail_nullreg)).any() ) { // some phenotypes failed (at null reg or null firth) - write their names to file

    ArrayXb pheno_flagged = (!has_converged) && (params->pheno_pass || params->pheno_fail_nullreg);
    Files outf;
    string failed_file = files->out_file + "_failedNullFirth_chr" + to_string(chrom) + ".list";
    outf.openForWrite( failed_file, sout);
    for( int i = 0; i < params->n_pheno; ++i )
      if(pheno_flagged(i))
        outf << files->pheno_names[i] << endl;
    outf.closeFile();
    sout << "WARNING: null Firth failed for " << pheno_flagged.count() << " phenotypes (list of traits written to '" << failed_file << "' and these will be skipped)\n";
    params->pheno_pass = has_converged;

  }
  if(params->blup_cov)
    pheno_data->new_cov.rightCols(1).array() = 0;

  if(silent || !params->firth) return;

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
  double tol = null_fit ? (10*params->numtol) : params->numtol_firth;
  double dev, lrt, dev0 = 0;

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

  // covariate effects added as offset in firth approx.
  if( params->firth_approx && !null_fit ) offset = fest->cov_blup_offset.col(ph).array(); 
  else offset = m_ests->blups.col(ph).array(); 

  // starting values
  if(null_fit){

    betaold = ArrayXd::Zero(Xmat.cols()); // last entry in exact Firth is kept at 0
    if(params->firth_approx){
      // start at logit^-1(mean(Y))-mean(offset)
      betaold(0) = ( 0.5 + mask.select(Y,0).sum())  / (pheno_data->Neff(ph) + 1);
      betaold(0) = log( betaold(0) / (1 - betaold(0) ));

      // LOCO prediction is offset
      betaold(0) -= mask.select(offset,0).mean();
    } else betaold.head(col_incl) = params->cov_betas.col(ph); // start at estimates from null firth with no snp column
      
  } else {
    // start at 0
    if(params->firth_approx) betaold = ArrayXd::Zero(col_incl); 
    // start at estimate from null fit
    else betaold = dt_thr->beta_null_firth.col(0);
  }

  success = fit_firth_pseudo(dev0, Y, Xmat, offset, mask, pivec, etavec, betaold, se, col_incl, dev, !null_fit, lrt, maxstep, niter/2, tol, params); // try pseudo

  // If didn't converge
  if(!success){
    if(!null_fit) { // reset beta
      if(params->firth_approx) betaold = ArrayXd::Zero(col_incl); 
      else betaold = dt_thr->beta_null_firth.col(0);
    } else if (fabs(betaold(0))>1e12) {
      if(params->firth_approx) betaold = 0;
      else betaold.head(col_incl) = params->cov_betas.col(ph);
  }
    success = fit_firth(ph, Y, Xmat, offset, mask, pivec, etavec, betaold, se, col_incl, dev, !null_fit, lrt, maxstep, niter/2, tol, params); // try NR (slower)

    if(!success){
      if(params->verbose) cerr << "WARNING: Logistic regression with Firth correction did not converge!\n";
      block_info->test_fail(ph) = true;
      return ;
    }
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

// for approx firth testing step
void fit_firth_logistic_snp_fast(int const& chrom, int const& ph, int const& isnp, bool const& null_fit, struct param const* params, struct phenodt* pheno_data, struct ests const* m_ests, struct f_ests const* fest, const Ref<const VectorXd>& Gvec, variant_block* block_info, data_thread* dt_thr, mstream& sout) {
  // if firth is used, fit based on penalized log-likelihood

  uint fit_state;
  int maxstep = params->maxstep;
  int niter = params->niter_max_firth, niter_pseudo = min(niter/2, 50), niter_nr = niter/2;
  double tol = params->numtol_firth;
  double lrt, dev0 = 0;

  double bstart = 0, betaold, se;
  ArrayXd offset;

  MapArXd Y (pheno_data->phenotypes_raw.col(ph).data(), pheno_data->phenotypes_raw.rows());
  MapArXb mask (pheno_data->masked_indivs.col(ph).data(), pheno_data->masked_indivs.rows());
  
  // For rare variants, set entries in Gvec for non-carriers to 0
  int mac_thr_sparse = (params->skip_fast_firth ? 0 : 50), i = 0, index_j;
  ArrayXi index_carriers;
  if(dt_thr->is_sparse && (block_info->mac(ph) < mac_thr_sparse)) {
    index_carriers.resize(dt_thr->Gsparse.nonZeros());
    for (SpVec::InnerIterator it(dt_thr->Gsparse); it; ++it) {
      index_j = it.index();
      // check for small entries in G (eg with imputed data)
      if(mask(index_j) && (it.value() > 1e-4)) index_carriers(i++) = index_j;
      }
    index_carriers.conservativeResize(i);
    niter_pseudo = niter/2;
  }

  // warm starts using estimates ignoring covariates
  if( params->htp_out && (block_info->genocounts(2,ph) == 0) && (block_info->genocounts(5,ph) == 0) )
    bstart = log( (block_info->genocounts(1,ph) + 0.5) * (block_info->genocounts(3,ph) + 0.5) / (block_info->genocounts(0,ph) + 0.5) / (block_info->genocounts(4,ph) + 0.5) );
  betaold = bstart;

  // covariate effects added as offset in firth approx.
  offset = fest->cov_blup_offset.col(ph).array(); 

  // get dev0
  ArrayXd pivec, wvec, Gvec_mask;
  get_pvec(pivec, offset, params->numtol_eps);
  dev0 = get_logist_dev(Y, pivec, mask);
  if((index_carriers.size() > 0)) { // bug fix to use the right deviance fn if using approximate penalty based on carrier status
    get_pvec(pivec, offset(index_carriers), params->numtol_eps);
    get_wvec(pivec, wvec, mask(index_carriers));
    Gvec_mask = Gvec(index_carriers);
  } else {
    get_wvec(pivec, wvec, mask);
    Gvec_mask = mask.select(Gvec.array(),0);
  }
  dev0 -= log( (Gvec_mask.square() * wvec).sum() );

  // fit state =
  //  0 - fit was successful
  //  1 - too slow convergence
  //  2 - diff_beta increased
  //  3 - fitted p = 0
  //  4 - lrt < 0
  fit_state = fit_firth_pseudo(dev0, Y, Gvec, offset, mask, index_carriers, betaold, se, lrt, maxstep, niter_pseudo, tol, params); // try pseudo

  // If didn't converge, try again with NR at 0
  if(fit_state && (bstart != 0) && index_carriers.size()) {
    if(params->debug) cerr << "WARNING: Pseudo-firth did not converge (" << fit_state << "; LRT = " << lrt << "; dev0 = " << dev0 << ") !\n";
    betaold = 0;
    fit_state = !fit_firth(dev0, Y, Gvec, offset, mask, index_carriers, betaold, se, lrt, maxstep, 100, tol, params); // try NR (slower)
  }

  // If didn't converge, try with NR
  if(fit_state){
    if(params->debug) cerr << "WARNING: NR-firth did not converge (" << fit_state << "; LRT = " << lrt << ") !\n";
    betaold = bstart; 
    fit_state = !fit_firth(dev0, Y, Gvec, offset, mask, index_carriers, betaold, se, lrt, maxstep, niter_nr, tol, params); // try NR (slower)
  }

  if(fit_state){
    if(params->verbose) cerr << "WARNING: Logistic regression with Firth correction did not converge (" << fit_state << "; LRT = " << lrt << ") !\n";
    block_info->test_fail(ph) = true;
    return ;
  }
  // sout << "\nNiter = " << niter_cur << " : " << mod_score.matrix().transpose() << endl;

  // compute beta_hat
  dt_thr->bhat(ph) = betaold;
  // compute SE based on Hessian for unpenalized LL
  if(!params->back_correct_se)
    dt_thr->se_b(ph) = se;

  if( lrt < 0 ) {
    block_info->test_fail(ph) = true;
    return ;
  }
  dt_thr->dif_deviance = lrt;

  return ;
}

// use NR or ADAM for Firth
bool fit_firth(int const& ph, const Ref<const ArrayXd>& Y1, const Ref<const MatrixXd>& X1, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, ArrayXd& pivec, ArrayXd& etavec, ArrayXd& betavec, ArrayXd& sevec, int const& cols_incl, double& dev, bool const& comp_lrt, double& lrt, int const& maxstep_firth, int const& niter_firth, double const& tol, struct param const* params) {

  double dev0 = 0;

  // get starting beta from ADAM ( compute and save null deviance )
  if(!comp_lrt && params->use_adam) 
    fit_firth_adam(ph, dev0, Y1, X1, offset, mask, pivec, etavec, betavec, sevec, cols_incl, dev, comp_lrt, lrt, params);

  return fit_firth_nr(dev0, Y1, X1, offset, mask, pivec, etavec, betavec, sevec, cols_incl, dev, comp_lrt, lrt, maxstep_firth, niter_firth, tol, params);

}

// fit based on penalized log-likelihood using NR
bool fit_firth_nr(double& dev0, const Ref<const ArrayXd>& Y1, const Ref<const MatrixXd>& X1, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, ArrayXd& pivec, ArrayXd& etavec, ArrayXd& betavec, ArrayXd& sevec, int const& cols_incl, double& dev, bool const& comp_lrt, double& lrt, int const& maxstep_firth, int const& niter_firth, double const& tol, struct param const* params) {
  // fit with first cols_incl columns of X1 (non-used entries of betavec should be 0)
  // else assuming using all columns 

  int niter_cur = 0, niter_search, nc = X1.cols(), n_score_inc = 0;
  double dev_old=0, dev_new=0, denum, mx, bdiff = 1, score_max_new, score_max_old = 1e16;

  ArrayXd hvec, mod_score;
  ArrayXd betanew, step_size, wvec;
  MatrixXd XtW, XtWX;
  ColPivHouseholderQR<MatrixXd> qr, qrX;

  if(params->debug) cerr << "\nFirth starting beta = " << betavec.matrix().transpose() << "\n";

  // solve S'(beta) = S(beta) + X'(h*(0.5-p)) = 0
  betanew = betavec * 0;
  while(niter_cur++ < niter_firth){

    // update quantities
    get_pvec(etavec, pivec, betavec, offset, X1, params->numtol_eps);
    dev_old = get_logist_dev(Y1, pivec, mask);
    get_wvec(pivec, wvec, mask);
    XtW = X1.transpose() * wvec.sqrt().matrix().asDiagonal();
    XtWX = XtW * XtW.transpose();
    qr.compute(XtWX);
    // compute deviance
    dev_old -= qr.logAbsDeterminant();
    if(comp_lrt && (niter_cur == 1)) // at first iter (i.e. betaSNP=0)
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
    // edit 5.31.12 for edge cases with approx Firth
    score_max_new = mod_score.abs().maxCoeff();
    if( ( score_max_new < tol) && (niter_cur >= 2) ) break;

    // try to catch convergence failures early
    if(!comp_lrt){
      if( score_max_new > score_max_old ) n_score_inc++; // track consecutive increases
      else n_score_inc = 0;
      if( n_score_inc > 10 ) return false;
    }

    // force absolute step size to be less than maxstep for each entry of beta
    mx = step_size.abs().maxCoeff() / maxstep_firth;
    if( mx > 1 ) step_size /= mx;

    // start step-halving and stop when deviance decreases 
    denum = 2;
    for( niter_search = 1; niter_search <= params->niter_max_line_search; niter_search++ ){

      // adjusted step size
      if(niter_search > 1) step_size /= denum;

      ///////// compute corresponding deviance
      if(cols_incl < nc) 
        betanew.head(cols_incl) = betavec.head(cols_incl) + step_size;
      else 
        betanew = betavec + step_size;

      get_pvec(etavec, pivec, betanew, offset, X1, params->numtol_eps);
      dev_new = get_logist_dev(Y1, pivec, mask);
      get_wvec(pivec, wvec, mask);
      XtW = X1.transpose() * wvec.sqrt().matrix().asDiagonal();
      XtWX = XtW * XtW.transpose();
      qr.compute(XtWX);
      dev_new -= qr.logAbsDeterminant();

      if(params->debug){
        if(niter_search == 1) bdiff = step_size.abs().maxCoeff();
        cerr << "["<<niter_cur << ":" << niter_search <<"] L1=" << setprecision(16)<< dev_new << "/L0="<< dev_old<< "\n";
      }
      if( dev_new < dev_old ) break;
    }

    if( niter_search > params->niter_max_line_search ) {
      if( comp_lrt ) step_size(0) += 1e-6;
      else return false; // step-halving failed
    }

    if(params->debug) cerr << "[" << niter_cur <<setprecision(16)<< "] beta.head=(" << betanew.head(min(5,cols_incl)).matrix().transpose() << "...); beta_diff.max=" << bdiff << "; score.max=" << score_max_new << "\n";


    if(cols_incl < nc)  
      betavec.head(cols_incl) += step_size;
    else
      betavec += step_size;
    dev_old = dev_new;
    score_max_old = score_max_new;

  }
  if(params->debug) cerr << "Ni=" << niter_cur<<setprecision(16) << "; beta.head=(" << betavec.head(min(15,cols_incl)).matrix().transpose() << "); score.max=" << mod_score.abs().maxCoeff() << "\n";

  // If didn't converge
  if( niter_cur > niter_firth ) return false;

  dev = dev_new;
  if( comp_lrt ) {
    lrt = dev0 - dev_new;
    if(lrt < 0) return false;

    sevec = qr.inverse().diagonal().array().sqrt();
  }

  return true;
}

// using pseudo-data representation with unpenalized logistic (strategy from brglm)
bool fit_firth_pseudo(double& dev0, const Ref<const ArrayXd>& Y1, const Ref<const MatrixXd>& X1, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, ArrayXd& pivec, ArrayXd& etavec, ArrayXd& betavec, ArrayXd& sevec, int const& cols_incl, double& dev, bool const& comp_lrt, double& lrt, int const& maxstep_firth, int const& niter_firth, double const& tol, struct param const* params) {
  // fit with first cols_incl columns of X1 (non-used entries of betavec should be 0)
  // else assuming using all columns 

  int niter_cur = 0, niter_log = 0, niter_search, niter_max = 25, nc = X1.cols(), niter_score_max_unchanged = 0;
  double dev_new=0, mx, maxstep = (comp_lrt && cols_incl == 1) ? 5 : maxstep_firth;
  double bdiff=1e16, bdiff_new=1e16;
  double score_max_old = 1e16, score_max_new;
  //double dev_log0, dev_log1=0;

  ArrayXd hvec, mod_score, ystar, score;
  ArrayXd betanew, step_size, wvec, zvec;
  MatrixXd XtW, XtWX;
  ColPivHouseholderQR<MatrixXd> qr, qrX;

  if(params->debug) cerr << "\nPseudo-firth starting beta = " << betavec.matrix().transpose() << "\n";

  betanew = betavec * 0;
  while(niter_cur++ < niter_firth){

    // update quantities
    get_pvec(etavec, pivec, betavec, offset, X1, params->numtol_eps);
    dev_new = get_logist_dev(Y1, pivec, mask);
    get_wvec(pivec, wvec, mask);
    XtW = X1.transpose() * wvec.sqrt().matrix().asDiagonal();
    XtWX = XtW * XtW.transpose();
    qr.compute(XtWX);
    // compute deviance
    dev_new -= qr.logAbsDeterminant();
    if(comp_lrt && (niter_cur == 1)) // at first iter (i.e. betaSNP=0)
      dev0 = dev_new;

    // compute diag(H), H = U(U'U)^{-1}U', U = Gamma^(1/2)X
    hvec = (qr.solve(XtW).array() * XtW.array() ).colwise().sum();
    // compute pseudo-response
    ystar = Y1 + hvec * (0.5 - pivec); 
    // compute modified score & step size
    if(cols_incl < nc) { 
      qrX.compute(XtWX.block(0, 0, cols_incl, cols_incl));
      mod_score = (X1.leftCols(cols_incl).transpose() * mask.select(ystar - pivec, 0).matrix() ).array();
    } else {
      mod_score = (X1.transpose() * mask.select(ystar - pivec, 0).matrix() ).array();
    }

    // stopping criterion using modified score function
    // edit 5.31.12 for edge cases with approx Firth
    score_max_new = mod_score.abs().maxCoeff();
    if( (score_max_new < tol) && (niter_cur >= 2) ) {
      if(params->debug) cerr << "stopping criterion met (" << score_max_new << " < " << tol << ")\n";
      break;
    }
    if(params->debug) cerr << "[" << niter_cur <<setprecision(16)<< "] beta.head=(" << betavec.head(min(5,cols_incl)).matrix().transpose() << "...); score.max=" << score_max_new << "\n";
    // to catch convergence failure sooner
    if( (niter_cur > 2) && (fabs(betavec(0)) > 1e13) ) return false;
    if(niter_score_max_unchanged > 3) return false;
    if( (niter_cur > 50) && ((score_max_new > 1000) || (betavec.abs().maxCoeff() > 1e12)) ) return false;

    // fit unpenalized logistic on transformed Y
    niter_log = 0;
    bdiff = 1e16;
    //dev_log0 = std::numeric_limits<double>::max();
    while(niter_log++ < niter_max){
      // p*(1-p) and check for zeroes
      if( get_wvec(pivec, wvec, mask, params->numtol_eps) ){
        if(params->debug) cerr << "WARNING: pseudo-firth gave fitted p=0 in logistic reg step\n";
        return false;
      }
      XtW = X1.leftCols(cols_incl).transpose() * mask.select(wvec,0).matrix().asDiagonal();
      XtWX = XtW * X1.leftCols(cols_incl);
      // working vector z = X*beta + (Y-p)/(p*(1-p))
      zvec = mask.select(etavec - offset + (ystar - pivec) / wvec, 0);
      // parameter estimate
      betanew.head(cols_incl) = ( XtWX ).colPivHouseholderQr().solve( XtW * zvec.matrix() ).array();

    // force absolute step size to be less than maxstep for each entry of beta
      if(comp_lrt && (cols_incl == 1)){ // only do this when testing each SNP
        step_size = betanew.head(cols_incl) - betavec.head(cols_incl);
        bdiff_new = fabs(step_size(0));
        if(bdiff_new > bdiff) { // step size should get smaller closer to soln
          if(params->debug) cerr << "WARNING: bdiff in pseudo-firth increased (" << bdiff << " -> " << bdiff_new << ")\n";
          return false; 
        }
        mx = bdiff_new / maxstep;
        if( mx > 1 ) {
          betanew.head(cols_incl) = betavec.head(cols_incl) + step_size / mx;
          if(params->debug) cerr << "step = " << step_size(0) << " -- mx = " << mx << " -- beta = " << betanew(0) << "\n";
        }
      }

      // skip step-halving
      for( niter_search = 1; niter_search <= params->niter_max_line_search; niter_search++ ){
        get_pvec(etavec, pivec, betanew, offset, X1, params->numtol_eps);
        //dev_log1 = get_deviance_logistic((ystar + 0.5 * hvec)/(1+hvec), pivec, 1 + hvec, mask);
        //if(params->debug) cerr << "[[" << niter_log << " - " << niter_search <<setprecision(16) << "]] D0=" << dev_log0 << " -> D1=" << dev_log1 << "\n";
        //if( dev_log1 < dev_log0 ) break;
        break;
        // adjust step size
        //betanew = (betavec + betanew) / 2;
      }
      /*if( niter_search > params->niter_max_line_search ){
        if(params->debug) cerr << "step halving failed in pseudo-firth log. reg step\n";
        return false; // step-halving failed
      }*/

      // stopping criterion
      score = X1.leftCols(cols_incl).transpose() * mask.select(ystar - pivec, 0).matrix();
      if( score.abs().maxCoeff() < tol ) break; // prefer for score to be below tol

      if(params->debug) cerr << "[[" << niter_log <<setprecision(16) << "]] beta.head=(" << betanew.head(min(5,cols_incl)).matrix().transpose() << "...); bdiff=" << bdiff_new << "; score.max=" << score.abs().maxCoeff() << "\n";

      betavec = betanew;
      if(comp_lrt && (cols_incl == 1)) bdiff = bdiff_new;
      //dev_log0 = dev_log1;
    }
    if( niter_log > params->niter_max ) return false;

    betavec = betanew;
    if(score_max_new < score_max_old) {
      score_max_old = score_max_new;
      niter_score_max_unchanged = 0;
    } else niter_score_max_unchanged++;
  }

  if(params->debug) cerr << "Ni=" << niter_cur<<setprecision(16) << "; beta.head=(" << betavec.head(min(15,cols_incl)).matrix().transpose() << "); score.max=" << mod_score.abs().maxCoeff() << "\n";

  // If didn't converge
  if( niter_cur > niter_firth ) return false;

  dev = dev_new;
  if( comp_lrt ) {
    lrt = dev0 - dev_new;
    if(lrt < 0) return false;

    sevec = qr.inverse().diagonal().array().sqrt();
  }

  return true;
}

// for approx firth testing step
uint fit_firth_pseudo(double const& dev0, const Ref<const ArrayXd>& Y1, const Ref<const VectorXd>& Gvec, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, const Ref<const ArrayXi>& index_carriers, double& betavec, double& sevec, double& lrt, int const& maxstep_firth, int const& niter_firth, double const& tol, struct param const* params) {

  bool fastFirth = index_carriers.size() > 0;
  int niter_cur = 0, niter_log = 0, niter_max = 25;
  double dev_new=0, dev_non_carriers = 0, mx, maxstep = 5;
  double bdiff=1e16, bdiff_new=1e16;
  //double dev_log0, dev_log1=0;

  double score, betanew = 0, step_size, XtWX = 0, beta_itr_14 = 0;
  ArrayXd hvec, ystar, etavec, pivec, wvec, XtWX_diag, Gvec_mask, Gvec_sq;

  if(fastFirth) {
    get_pvec(etavec, pivec, betavec, offset, Gvec, params->numtol_eps);
    dev_new = get_logist_dev(Y1, pivec, mask);
    dev_non_carriers = dev_new - get_logist_dev(Y1(index_carriers), pivec(index_carriers), mask(index_carriers));
    Gvec_mask = Gvec(index_carriers);
  } else Gvec_mask = mask.select(Gvec.array(),0);
  Gvec_sq = Gvec_mask.square();

  if(params->debug) cerr << "\nPseudo-firth (fast) starting beta = " << betavec << "\n";

  while(niter_cur++ < niter_firth){

    // update quantities
    if(fastFirth) {
      get_pvec(etavec, pivec, betavec, offset(index_carriers), Gvec(index_carriers), params->numtol_eps);
      dev_new = dev_non_carriers + get_logist_dev(Y1(index_carriers), pivec, mask(index_carriers));
      get_wvec(pivec, wvec, mask(index_carriers));
    } else {
      get_pvec(etavec, pivec, betavec, offset, Gvec, params->numtol_eps);
      dev_new = get_logist_dev(Y1, pivec, mask);
      get_wvec(pivec, wvec, mask);
    }
    XtWX_diag = Gvec_sq * wvec;
    XtWX = XtWX_diag.sum();
    // compute deviance
    dev_new -= log(XtWX);

    // compute diag(H), H = U(U'U)^{-1}U', U = Gamma^(1/2)X
    hvec = XtWX_diag / XtWX;
    // compute pseudo-response 
    ystar = (fastFirth ? Y1(index_carriers) : Y1) + hvec * (0.5 - pivec); 
    // compute modified score & step size
    score = (Gvec_mask * (ystar - pivec)).sum();

    // stopping criterion using modified score function
    // edit 5.31.12 for edge cases with approx Firth
    if( (fabs(score) < tol) && (niter_cur >= 2) ) {
      if(params->debug) cerr << "stopping criterion met (|" << score << "| < " << tol << ")\n";
      break;
    }
    if(params->debug) cerr << "[" << niter_cur <<setprecision(16)<< "] beta.head=(" << betavec << "...); score=" << score << "\n";

    // check for change in beta at iteration 15 (if too large, try with NR)
    if(niter_cur == 14) beta_itr_14 = betavec;
    if((niter_cur == 15) && (fabs(betavec - beta_itr_14) > .1)) return 1;

    // fit unpenalized logistic on transformed Y
    niter_log = 0;
    bdiff = 1e16;
    //dev_log0 = std::numeric_limits<double>::max();
    while(niter_log++ < niter_max){

      // force absolute step size to be less than maxstep for each entry of beta
      step_size = score / XtWX;
      bdiff_new = fabs(step_size);
      if(bdiff_new > bdiff) { // step size should get smaller closer to soln
        if(params->debug) cerr << "WARNING: bdiff in pseudo-firth increased (" << bdiff << " -> " << bdiff_new << ")\n";
        return 2; 
      }
      mx = bdiff_new / maxstep;

      // parameter estimate
      if( mx > 1 ) {
        betanew = betavec + step_size / mx;
        if(params->debug) cerr << "step = " << step_size << " -- mx = " << mx << " -- beta = " << betanew << "\n";
      } else betanew = betavec + step_size;

      // compute score at new beta
      if(fastFirth) get_pvec(etavec, pivec, betanew, offset(index_carriers), Gvec(index_carriers), params->numtol_eps); 
      else get_pvec(etavec, pivec, betanew, offset, Gvec, params->numtol_eps);
      score = (Gvec_mask * (ystar - pivec)).sum();
      if( fabs(score) < tol ) break; // prefer for score to be below tol

      if(params->debug) cerr << "[[" << niter_log <<setprecision(16) << "]] beta=(" << betanew << "...); bdiff=" << bdiff_new << "; score=" << score << "\n";

      // p*(1-p) and check for zeroes
      if( get_wvec(pivec, wvec, (fastFirth ? mask(index_carriers) : mask), params->numtol_eps) ) {
        if(params->debug) cerr << "WARNING: pseudo-firth gave fitted p=0 in logistic reg step\n";
        return 3;
      }
      XtWX_diag = Gvec_sq * wvec;
      XtWX = XtWX_diag.sum();

      betavec = betanew;
      bdiff = bdiff_new;
      //dev_log0 = dev_log1;
    }
    if( niter_log > params->niter_max ) return 1;

    betavec = betanew;
  }

  if(params->debug) cerr << "Ni=" << niter_cur<<setprecision(16) << "; beta=(" << betavec << "); score=" << score << "\n";

  // If didn't converge
  if( niter_cur > niter_firth ) return 1;

  lrt = dev0 - dev_new;
  if(lrt < 0) return 4;

  sevec = sqrt(1/XtWX);

  return 0;
}

// for approx firth testing step (using NR)
bool fit_firth(double const& dev0, const Ref<const ArrayXd>& Y1, const Ref<const VectorXd>& X1, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, const Ref<const ArrayXi>& index_carriers, double& betavec, double& sevec, double& lrt, int const& maxstep_firth, int const& niter_firth, double const& tol, struct param const* params) {

  bool fastFirth = index_carriers.size() > 0;
  int niter_cur = 0, niter_search;
  double dev_old=0, dev_new=0, dev_non_carriers = 0, denum, mx;
  double bdiff=1e16;

  double score, betanew = 0, step_size, XtWX = 0;
  ArrayXd hvec, etavec, pivec, wvec, XtWX_diag, Gvec_mask, Gvec_sq;
 
  get_pvec(etavec, pivec, betavec, offset, X1, params->numtol_eps);
  dev_old = get_logist_dev(Y1, pivec, mask);
  if(fastFirth) {
    get_pvec(etavec, pivec, betavec, offset(index_carriers), X1(index_carriers), params->numtol_eps);
    dev_non_carriers = dev_old - get_logist_dev(Y1(index_carriers), pivec, mask(index_carriers));
    get_wvec(pivec, wvec, mask(index_carriers));
    Gvec_mask = X1(index_carriers);
  } else {
    get_wvec(pivec, wvec, mask);
    Gvec_mask = mask.select(X1.array(),0);
  }
  Gvec_sq = Gvec_mask.square();

  // solve S'(beta) = S(beta) + X'(h*(0.5-p)) = 0
  // starting values
  if(params->debug) cerr << "\nFirth starting beta = " << betavec << "\n";
  // compute deviance
  XtWX_diag = Gvec_sq * wvec;
  XtWX = XtWX_diag.sum();
  dev_old -= log(XtWX);

  while(niter_cur++ < niter_firth){

    // compute diag(H), H = U(U'U)^{-1}U', U = Gamma^(1/2)X
    hvec = XtWX_diag / XtWX;
    // compute modified score
    score = (Gvec_mask * ((fastFirth ? Y1(index_carriers) : Y1) - pivec + hvec * (0.5 - pivec))).sum();
    // stopping criterion using modified score function
    // edit 5.31.12 for edge cases with approx Firth
    if( (fabs(score) < tol) && (niter_cur >= 2) ) break;

    // force absolute step size to be less than maxstep for each entry of beta
    step_size = score / XtWX;
    bdiff = fabs(step_size);
    mx = bdiff / maxstep_firth;
    if( mx > 1 ) step_size /= mx;

    // start step-halving and stop when deviance decreases 
    denum = 2;
    for( niter_search = 1; niter_search <= params->niter_max_line_search; niter_search++ ){

      // adjusted step size
      if(niter_search > 1) step_size /= denum;

      betanew = betavec + step_size;

      ///////// compute corresponding deviance
      if(fastFirth) {
        get_pvec(etavec, pivec, betanew, offset(index_carriers), X1(index_carriers), params->numtol_eps); 
        dev_new = dev_non_carriers + get_logist_dev(Y1(index_carriers), pivec, mask(index_carriers));
      } else {
        get_pvec(etavec, pivec, betanew, offset, X1, params->numtol_eps);
        dev_new = get_logist_dev(Y1, pivec, mask);
      }
      get_wvec(pivec, wvec, (fastFirth ? mask(index_carriers) : mask));
      XtWX_diag = Gvec_sq * wvec;
      XtWX = XtWX_diag.sum();
      dev_new -= log(XtWX);

      if(params->debug) cerr << "["<<niter_cur << ":" << niter_search <<"] L1=" << setprecision(16)<< dev_new << "/L0="<< dev_old<< "\n";
      if( dev_new < dev_old ) break;
    }

    if( niter_search > params->niter_max_line_search ) step_size += 1e-6;

    if(params->debug) cerr << "[" << niter_cur <<setprecision(16)<< "] beta=(" << betanew << "...); beta_diff.max=" << bdiff << "; score=" << score << "\n";

    betavec += step_size;
    dev_old = dev_new;

  }
  if(params->debug) cerr << "Ni=" << niter_cur<<setprecision(16) << "; beta=(" << betavec << "); score=" << score << "\n";

  // If didn't converge
  if( niter_cur > niter_firth ) return false;

  lrt = dev0 - dev_new;
  if(lrt < 0) return false;

  sevec = sqrt(1/XtWX);

  return true;
}

// fit based on penalized log-likelihood using ADAM
bool fit_firth_adam(int const& ph, double& dev0, const Ref<const ArrayXd>& Y1, const Ref<const MatrixXd>& X1, const Ref<const ArrayXd>& offset, const Ref<const ArrayXb>& mask, ArrayXd& pivec, ArrayXd& etavec, ArrayXd& betavec, ArrayXd& sevec, int const& cols_incl, double& dev, bool const& comp_lrt, double& lrt, struct param const* params) {
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

    if(params->debug && (niter_cur>1) && (niter_cur%100==0) ) cerr << "\nNiter = " << niter_cur << " (beta = " << betavec.matrix().transpose() << ") : " << gradient_f.matrix().transpose() << endl;

    mt = p_beta1 * mt + (1 - p_beta1) * gradient_f;
    vt = p_beta2 * vt + (1 - p_beta2) * gradient_f.square();
    p_alpha_t = p_alpha * sqrt(1 - pow(p_beta2, niter_cur)) / (1 - pow(p_beta1, niter_cur));
    step_size = p_alpha_t * mt / (vt.sqrt() + p_eps);

    // stopping criterion
    if( step_size.abs().maxCoeff() < params->numtol) break;

    betavec.head(cols_incl) -= step_size;

  }
  if(params->debug) cerr << "ADAM took "<< niter_cur << " iterations (score max = " << gradient_f.abs().maxCoeff() << ")...";

  return (niter_cur <= params->niter_max_firth_adam);
}


string get_firth_est_allChr(struct in_files& files, struct filter const& filters, struct ests& m_ests, struct f_ests& fest, struct phenodt& pheno_data, struct param& params, mstream& sout){

  sout << "   -computing and storing null Firth estimates for all chromosomes..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  // go through each chromosome
  for(int chr = 1; chr <= params.nChrom; chr++){

    if(params.verbose) sout << "chr" << chr <<"..." << flush;

    // read the prs
    blup_read_chr(true, chr, m_ests, files, filters, pheno_data, params, sout);

    if (params.trait_mode == 1) {
      // run null logistic regression to get the starting values for firth
      fit_null_logistic(true, chr, &params, &pheno_data, &m_ests, &files, sout); // for all phenotypes

      // run null firth for each trait and write estimates to file
      fit_null_firth(true, chr, &fest, &pheno_data, &m_ests, &files, &params, sout);
    } else if (params.trait_mode == 3) {
      // run null cox regression to get the starting values for firth
      fit_null_cox(true, chr, &params, &pheno_data, &m_ests, &files, sout); // for all phenotypes

      // run null firth for each trait and write estimates to file
      fit_null_firth_cox(true, chr, &fest, &pheno_data, &m_ests, &files, &params, sout);
    }
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
      if( !params.pheno_pass(j) ) continue;
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
  files.null_firth_files.assign(params.n_pheno, "");
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

  // // force all phenotypes in phenotype file to be used
  // if(read.count() != params.n_pheno) 
  //   throw "number of valid step 1 files (" + to_string( read.count() ) + ")  is not equal to the number of phenotypes." ;

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

    // if file has not been given, use 0 as start
    if(files->null_firth_files[i] == "") continue;
    if( !params->pheno_pass(i) ) continue;

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

void get_beta_start_firth(struct f_ests* firth_est, struct ests const* m_ests){
  // get b0 from null logistic regression
  firth_est->beta_null_firth.topRows(m_ests->bhat_start.rows()) = m_ests->bhat_start;
}

void check_pval_snp(variant_block* block_info, data_thread* dt_thr, int const& chrom, int const& ph, int const& isnp, struct phenodt& pheno_data, struct geno_block& gblock, struct ests const& m_ests, struct f_ests& fest, struct param const& params, mstream& sout){

  // if firth isn't used, or Tstat < threshold, no correction done
  if(!block_info->is_corrected(ph) || (fabs(dt_thr->stats(ph)) <= params.z_thr)){
    get_sumstats(false, ph, dt_thr);
    dt_thr->cal_factor(ph) = 1;
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

    run_SPA_test(block_info->test_fail(ph), ph, dt_thr, pheno_data.masked_indivs.col(ph).array(), m_ests, params);
    if(block_info->test_fail(ph)) {
      get_sumstats(true, ph, dt_thr);
      return;
    }

    dt_thr->se_b(ph) = 1 / sqrt(dt_thr->denum(ph));
    dt_thr->bhat(ph) = sgn(dt_thr->stats(ph)) * sqrt(dt_thr->chisq_val(ph)) * dt_thr->se_b(ph);

  }

  //if(params.debug) cerr << "uncorrected: " << dt_thr->stats(ph) * dt_thr->stats(ph) <<  "] -> " << dt_thr->chisq_val(ph) << endl;
  dt_thr->cal_factor(ph) =  dt_thr->chisq_val(ph) == 0 ? 0 : dt_thr->stats(ph) * dt_thr->stats(ph) / dt_thr->chisq_val(ph);

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
    if (params.trait_mode == 1) {
      // obtain null deviance (set SNP effect to 0 and compute max. pen. LL)
      fit_firth_logistic_snp(chrom, ph, isnp, true, &params, &pheno_data, &m_ests, &fest, gblock.Gmat.col(isnp), block_info, dt_thr, sout);
      if(block_info->test_fail(ph)) return ;
      // fit full model and compute deviance
      fit_firth_logistic_snp(chrom, ph, isnp, false, &params, &pheno_data, &m_ests, &fest, gblock.Gmat.col(isnp), block_info, dt_thr, sout);
    } else if (params.trait_mode == 3) {
      if (pheno_data.new_cov.cols() == 0) {
        fit_firth_cox_snp_fast(chrom, ph, isnp, &params, &pheno_data, &m_ests, &fest, dt_thr->Gres, block_info, dt_thr, sout);
      } else {
        fit_firth_cox_snp(chrom, ph, isnp, &params, &pheno_data, &m_ests, &fest, gblock.Gmat.col(isnp), block_info, dt_thr, sout);
      }
    }
  } else { // approx firth - only fit full model
    if (params.trait_mode == 1) {
      fit_firth_logistic_snp_fast(chrom, ph, isnp, false, &params, &pheno_data, &m_ests, &fest, dt_thr->Gres.cwiseQuotient(m_ests.Gamma_sqrt.col(ph)), block_info, dt_thr, sout);
    } else if (params.trait_mode == 3) {
      fit_firth_cox_snp_fast(chrom, ph, isnp, &params, &pheno_data, &m_ests, &fest, dt_thr->Gres, block_info, dt_thr, sout);
    }
  }
}

void run_SPA_test(bool& test_fail, int const& ph, data_thread* dt_thr, const Ref<const ArrayXb>& mask, struct ests const& m_ests, struct param const& params){
  run_SPA_test_snp(dt_thr->chisq_val(ph), dt_thr->pval_log(ph), dt_thr->stats(ph), dt_thr->denum(ph), dt_thr->fastSPA, dt_thr->Gsparse, dt_thr->Gres.array(), m_ests.Y_hat_p.col(ph).array(), m_ests.Gamma_sqrt.col(ph).array(), mask, test_fail, params.tol_spa, params.niter_max_spa, params.missing_value_double, params.nl_dbl_dmin);
}

void run_SPA_test_snp(double& chisq, double& pv, const double& stats, const double& denum, bool const& fastSPA, SpVec const& Gsparse, const Ref<const ArrayXd>& Gres, const Ref<const ArrayXd>& phat, const Ref<const ArrayXd>& Gamma_sqrt, const Ref<const ArrayXb>& mask, bool& test_fail, const double& tol, const double& niter_max, const double& missing_value_double, const double& nl_dbl_dmin){

  int index_j;
  double score_num, tval, limK1_low, limK1_high, root_K1, pval1, pval2;
  spa_data spa_df;
  ArrayXd Gmu;

  // compute needed quantities
  spa_df.val_c = sqrt( denum );  // sqrt( G'WG )
  score_num = stats * spa_df.val_c;
  spa_df.Gmod = Gres / Gamma_sqrt * mask.cast<double>();
  Gmu = spa_df.Gmod * phat;
  spa_df.val_a = Gmu.sum();
  spa_df.fastSPA = fastSPA;

  if(spa_df.fastSPA){
    spa_df.val_b = denum;
    spa_df.val_d = 0;
    for (SpVec::InnerIterator it(Gsparse); it; ++it) {
      index_j = it.index();
      if(!mask(index_j)) continue;
      spa_df.val_b -= Gres(index_j) * Gres(index_j);
      spa_df.val_d += Gmu(index_j);
    }
  }

  // check if K'(t)= s can be solved
  limK1_low = (spa_df.Gmod < 0).select(spa_df.Gmod, 0 ).sum() - spa_df.val_a ;
  limK1_high = (spa_df.Gmod > 0).select(spa_df.Gmod, 0 ).sum() - spa_df.val_a ;
  if( score_num < limK1_low || score_num > limK1_high ){
    //if(params.verbose) sout << "WARNING: SPA failed (solution to K'(t)=s is infinite)";
    test_fail = true;
    return;
  }

  tval = stats >= 0 ? -stats : stats;

  // 1.for T
  spa_df.pos_score = true;
  // solve K'(t)= tval using a mix of Newton-Raphson and bisection method
  root_K1 = solve_K1_snp(tval, denum, Gsparse, phat, Gamma_sqrt, spa_df, mask, tol, niter_max, missing_value_double);
  if( root_K1 == missing_value_double ){
    test_fail = true;
    return;
  }
  // compute pvalue (one tail)
  get_SPA_pvalue_snp(root_K1, tval, pval1, test_fail, denum, Gsparse, phat, Gamma_sqrt, spa_df, mask);
  if(test_fail) {return;}

  // 2.for -T
  spa_df.pos_score = false;
  // solve K'(t)= tval using a mix of Newton-Raphson and bisection method
  root_K1 = solve_K1_snp(tval, denum, Gsparse, phat, Gamma_sqrt, spa_df, mask, tol, niter_max, missing_value_double);
  if( root_K1 == missing_value_double ){
    test_fail = true;
    return;
  }
  // compute pvalue (other tail)
  get_SPA_pvalue_snp(root_K1, tval, pval2, test_fail, denum, Gsparse, phat, Gamma_sqrt, spa_df, mask);
  if(test_fail) {return;}

  // get quantile
  //cerr << "\nstats: " << stats << ":" << pval1 << " " << pval2 << "\n";
  if( (pval1 + pval2) > 1 ){
    test_fail = true;
    return;
  }
  get_logp(pval1+pval2, pv, chisq, nl_dbl_dmin);

}



// SPA (MT in OpenMP)
double solve_K1_snp(const double& tval, const double& denum, SpVec const& Gsparse, const Ref<const ArrayXd>& phat, const Ref<const ArrayXd>& Gamma_sqrt, struct spa_data& spa_df, const Ref<const ArrayXb>& mask, double const& tol, int const& niter_max, double const& missing_value_double){

  int niter_cur;
  int lambda = spa_df.pos_score ? 1 : -1; // if score is negative, adjust K' and K''
  double min_x, max_x, t_old, f_old, t_new = -1, f_new, hess;

  niter_cur = 0;
  if(tval >=0){min_x = 0, max_x = std::numeric_limits<double>::max();}
  else{min_x = std::numeric_limits<double>::lowest(), max_x = 0;}
  t_old = 0;
  f_old = spa_df.fastSPA ? compute_K1_fast_snp(lambda * t_old, spa_df.val_b, spa_df.val_c, spa_df.val_d, denum, Gsparse, spa_df.Gmod, phat, mask) : compute_K1_snp(lambda * t_old, spa_df.val_a, spa_df.val_c, spa_df.Gmod, phat, mask);
  f_old *= lambda;
  f_old -= tval; 

  while( niter_cur++ < niter_max ){

    hess = spa_df.fastSPA ? compute_K2_fast_snp(lambda * t_old, spa_df.val_b, spa_df.val_c, spa_df.val_d, denum, Gsparse, spa_df.Gmod, phat, Gamma_sqrt, mask) : compute_K2_snp(lambda * t_old, spa_df.val_a, spa_df.val_c, spa_df.Gmod, phat, Gamma_sqrt, mask);
    if(hess == 0) return missing_value_double;
    t_new = t_old - f_old / hess;
    f_new = spa_df.fastSPA ? compute_K1_fast_snp(lambda * t_new, spa_df.val_b, spa_df.val_c, spa_df.val_d, denum, Gsparse, spa_df.Gmod, phat, mask) : compute_K1_snp(lambda * t_new, spa_df.val_a, spa_df.val_c, spa_df.Gmod, phat, mask);
    f_new *= lambda;
    f_new -= tval;

    if( fabs( f_new ) < tol ) break;

    // update bounds on root
    if( t_new && (t_new > min_x) && (t_new < max_x) ){
      if( f_new > 0) max_x = t_new;
      else min_x = t_new;
    } else{ // bisection method if t_new went out of bounds and re-compute f_new
      t_new = ( min_x + max_x ) / 2;
      // if( fabs( min_x - t_new ) < params->tol_spa ) break;
      f_new = spa_df.fastSPA ? compute_K1_fast_snp(lambda * t_new, spa_df.val_b, spa_df.val_c, spa_df.val_d, denum, Gsparse, spa_df.Gmod, phat, mask) : compute_K1_snp(lambda * t_new, spa_df.val_a, spa_df.val_c, spa_df.Gmod, phat, mask);
      f_new *= lambda;
      f_new -= tval;
      // reduce bounds based on new value
      if(f_new <= 0) min_x = t_new;
      else max_x = t_new;
    }

    t_old = t_new;
    f_old = f_new;
  }

  // If didn't converge
  if( niter_cur > niter_max ){
    //if(params->verbose) sout << "WARNING: SPA did not converge to root for K'(t)=s.\n";
    return missing_value_double;
  }
  //sout << "#iterations = " << niter_cur << "; f= " << f_new << endl;

  return t_new;
}

double compute_K_snp(const double& t, const double& a, const double& c, const Ref<const ArrayXd>& Gmod, const Ref<const ArrayXd>& phat, const Ref<const ArrayXb>& mask){
  double val = mask.select( ( 1 - phat + phat * ( t / c * Gmod ).exp() ).log(), 0).sum() - t * a / c;

  return val;
}

double compute_K_fast_snp(const double& t, const double& b, const double& c, const double& d, const double& denum, SpVec const& Gsparse, const Ref<const ArrayXd>& Gmod, const Ref<const ArrayXd>& phat, const Ref<const ArrayXb>& mask){

  uint32_t index_j;
  double val = 0;

  for (SpVec::InnerIterator it(Gsparse); it; ++it) {
    index_j = it.index();
    if(!mask(index_j)) continue;

    val += log( 1 - phat(index_j) + phat(index_j) * exp( t / c * Gmod(index_j)) );
  }
  val += -t * d / c + t * t / 2 / denum * b;

  return val;
}

double compute_K1_snp(const double& t, const double& a, const double& c, const Ref<const ArrayXd>& Gmod, const Ref<const ArrayXd>& phat, const Ref<const ArrayXb>& mask){

  double val = mask.select( ( Gmod * phat / c ) / ( phat + (1 - phat) * ( -t / c * Gmod ).exp() ), 0).sum();
  val -= a / c;

  return val;
}

double compute_K1_fast_snp(const double& t, const double& b, const double& c, const double& d, const double& denum, SpVec const& Gsparse, const Ref<const ArrayXd>& Gmod, const Ref<const ArrayXd>& phat, const Ref<const ArrayXb>& mask){

  uint32_t index_j;
  double val = 0;

  for (SpVec::InnerIterator it(Gsparse); it; ++it) {
    index_j = it.index();
    if(!mask(index_j)) continue;

    val += ( Gmod(index_j) * phat(index_j) / c ) / ( phat(index_j) + (1 - phat(index_j)) * exp( -t / c * Gmod(index_j)) );
  }
  val += -d / c + t / denum * b;

  return val;
}

double compute_K2_snp(const double& t, const double& a, const double& c, const Ref<const ArrayXd>& Gmod, const Ref<const ArrayXd>& phat, const Ref<const ArrayXd>& Gamma_sqrt, const Ref<const ArrayXb>& mask){

  ArrayXd Vexp = -t / c * Gmod;
  if((mask && (Vexp > MAX_EXP_LIM)).any()) { return 0; }
  double val = mask.select( ( Gmod.square() * Gamma_sqrt.square() / (c*c) * Vexp.exp()) / ( phat + (1 - phat) * Vexp.exp() ).square(), 0).sum();

  return val;
}

double compute_K2_fast_snp(const double& t, const double& b, const double& c, const double& d, const double& denum, SpVec const& Gsparse, const Ref<const ArrayXd>& Gmod, const Ref<const ArrayXd>& phat, const Ref<const ArrayXd>& Gamma_sqrt, const Ref<const ArrayXb>& mask){

  uint32_t index_j;
  double val = 0, denum_v, vexp;

  for (SpVec::InnerIterator it(Gsparse); it; ++it) {
    index_j = it.index();
    if(!mask(index_j)) continue;
    vexp = -t / c * Gmod(index_j);
    if(vexp > MAX_EXP_LIM) { return 0; }
    denum_v = phat(index_j) + (1 - phat(index_j)) * exp( vexp );
    val += ( Gmod(index_j) * Gmod(index_j) * Gamma_sqrt(index_j) * Gamma_sqrt(index_j) * exp( vexp ) / (c*c) ) / (denum_v * denum_v);
    //cerr << "phat:" << phat(index_j) << "; t:"<< t << "; c:" << c << ";G:"<< Gmod(index_j)<< " ;denum:" << denum_v <<"\n";
  }
  val += b / denum;

  return val;
}

void get_SPA_pvalue_snp(const double& root, const double& tval, double& pv, bool& test_fail, const double& denum, SpVec const& Gsparse, const Ref<const ArrayXd>& phat, const Ref<const ArrayXd>& Gamma_sqrt, struct spa_data& spa_df, const Ref<const ArrayXb>& mask){

  int lambda = spa_df.pos_score ? 1 : -1; // if score is negative, adjust K and K''
  double kval, k2val, wval, vval, rval;
  normal nd(0,1);

  kval = spa_df.fastSPA ? compute_K_fast_snp(lambda * root, spa_df.val_b, spa_df.val_c, spa_df.val_d, denum, Gsparse, spa_df.Gmod, phat, mask) : compute_K_snp(lambda * root, spa_df.val_a, spa_df.val_c, spa_df.Gmod, phat, mask);
  k2val = spa_df.fastSPA ? compute_K2_fast_snp(lambda * root, spa_df.val_b, spa_df.val_c, spa_df.val_d, denum, Gsparse, spa_df.Gmod, phat, Gamma_sqrt, mask) : compute_K2_snp(lambda * root, spa_df.val_a, spa_df.val_c, spa_df.Gmod, phat, Gamma_sqrt, mask);
  if(k2val == 0) {
    test_fail = true;
    return;
  }

  wval = sgn(root) * sqrt( 2 * ( root * tval - kval ) );
  vval = root * sqrt( k2val );
  //cerr << " root:" << root << " kval:" << kval << " k2val:" << k2val << " wval:" << wval << " vval:" << vval << " ";
  if(vval == 0) { // root is 0 so s=0 (K'(0)=0)
    pv = 0.5;
  } else {
    rval = wval + log( vval / wval ) / wval;
    pv = cdf(nd, rval); // one-sided
  }
  test_fail = false;
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
    if(params->trait_set) 
      return print_header_output_all_multitrait(params);
    else if(params->multiphen)
      return print_header_output_all_multiphen(params);
    else
      return print_header_output_all(params);
}

std::string print_header_output_all_multiphen(struct param const* params){

  std::ostringstream buffer;

  buffer << "CHROM GENPOS ID ALLELE0 ALLELE1 MAC A1FREQ N LOG10P MULTINOM IT UP FIRTH";
  buffer << endl;

  return buffer.str();
}

std::string print_header_output_all_multitrait(struct param const* params){

  std::ostringstream buffer;

  buffer << "CHROM GENPOS ID ALLELE0 ALLELE1 MAC A1FREQ N";
  // p-values for single-trait tests
  /* for(int i = 0; i < params->n_pheno; i++) { */
  /*   buffer << "LOG10P.Y0" << i+1 << " "; */
  /* } */
  buffer << " LOG10P.MINP0 LOG10Q.MINP0";
  // p-values for multi-trait tests
  buffer << " LOG10P.MANOVA LOG10P.OMNIBUS0 LOG10BF.BAYES LOG10P.NNLS0 LOG10P.SUMZ0 LOG10P.NPMANOVA LOG10P.HOMNIBUS0 LOG10P.CPC0"
      << " LOG10P.RCPC0SUMCHI2 LOG10P.RCPC0FISHER LOG10P.RCPC0ACAT" 
      << " LOG10P.ACPC0SUMCHI2 LOG10P.ACPC0FISHER LOG10P.ACPC0ACAT" 
      << " LOG10Q.NNLS0";
  // z-scores for single-trait models
  for(int i = 0; i < params->n_pheno; i++) {
    buffer << " " << "Z.Y0" << i+1;
  }
  // z-scores for PCs
  for(int i = 0; i < params->n_pheno; i++) {
    buffer << " " << "Z.PC0" << i+1;
  }
  // z-scores for Robust PCs
  for(int i = 0; i < params->n_pheno; i++) {
    buffer << " " << "Z.RPC0" << i+1;
  }
  // z-scores for Adjusted PCs
  for(int i = 0; i < params->n_pheno; i++) {
    buffer << " " << "Z.APC0" << i+1;
  }
  buffer << endl;

  return buffer.str();
}

std::string print_header_output_all(struct param const* params){

  int i;
  std::ostringstream buffer;

  buffer << "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ " << 
    ( params->af_cc ? "A1FREQ_CASES A1FREQ_CONTROLS ":"") <<
    ( !params->build_mask && params->dosage_mode ? "INFO ":"") 
    << "N " <<
    ( params->af_cc ? "N_CASES N_CONTROLS ":"") <<
    "N_RR N_RA N_AA " << // across all analyzed samples (with dosages then hardcounts)
    "TEST";

  for(i = 0; i < params->n_pheno; i++) 
    buffer << " BETA.Y" << i+1 << " SE.Y" << i+1 << " CHISQ.Y" << i+1 << " LOG10P.Y" << i+1;
  // end of line
  buffer << " EXTRA\n";

  return buffer.str();
}

std::string print_header_output_single(struct param const* params){

  std::ostringstream buffer;

  buffer << "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ " << 
    ( params->af_cc ? "A1FREQ_CASES A1FREQ_CONTROLS ":"") <<
    ( !params->build_mask && params->dosage_mode ? "INFO ":"") << 
    "N " <<
    ( params->af_cc ? "N_CASES N_CONTROLS ":"") <<
    "TEST BETA SE CHISQ LOG10P EXTRA\n";

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

  if(params.htp_out) {
    buffer <<  print_sum_stats_head_htp(snp_index, files.pheno_names[i], model_type, snpinfo, &params) << print_sum_stats_htp(dt_thr->bhat(i), dt_thr->se_b(i), dt_thr->chisq_val(i), dt_thr->pval_log(i), block_info->af(i), block_info->info(i), block_info->mac(i), block_info->genocounts, i, !block_info->test_fail(i), 1, &params, dt_thr->scores(i), dt_thr->cal_factor(i), -1, dt_thr->skat_var(i));
  } else {
    buffer << (!params.split_by_pheno && (i>0) ? "" : tmpstr) << print_sum_stats((params.split_by_pheno ? block_info->af(i) : block_info->af1), block_info->af_case(i), block_info->af_control(i), block_info->n_rr, block_info->n_aa, (params.split_by_pheno ? block_info->info(i) : block_info->info1), (params.split_by_pheno ? block_info->ns(i) : block_info->ns1), block_info->ns_case(i), block_info->ns_control(i), test_string, dt_thr->bhat(i), dt_thr->se_b(i), dt_thr->chisq_val(i), dt_thr->pval_log(i), !block_info->test_fail(i), 1, &params, (i+1));
  }

  return buffer.str();
}


//// test info for each snp
std::string print_sum_stats(const double& af, const double& af_case, const double& af_control, const int& n_rr, const int& n_aa, const double& info, const int& n, const int& ns_case, const int& ns_control, const string& model, const double& beta, const double& se, const double& chisq, const double& pv, const bool& test_pass, const int& df, struct param const* params, int const& ipheno){

  if(params->split_by_pheno)
    return print_sum_stats_single(af, af_case, af_control, info, n, ns_case, ns_control, model, beta, se, chisq, pv, test_pass, df, params);
  else
    return print_sum_stats_all(af, af_case, af_control, n_rr, n_aa, info, n, ns_case, ns_control, model, beta, se, chisq, pv, test_pass, df, params, ipheno);
}

// native format - all phenos
std::string print_sum_stats_all(const double& af, const double& af_case, const double& af_control, const int& n_rr, const int& n_aa, const double& info, const int& n, const int& ns_case, const int& ns_control, const string& model, const double& beta, const double& se, const double& chisq, const double& pv, const bool& test_pass, const int& df, struct param const* params, int const& ipheno){

  std::ostringstream buffer;
  bool print_afs = (af >= 0), print_info = (info >= 0), print_se = (se >= 0), print_genoc = (n_rr >= 0);
  bool print_pv = (chisq>=0) && test_pass;

  // AF N INFO TEST
  if(ipheno == 1) {
    if(print_afs) buffer << af ;
    else buffer << "NA" ;
    if( params->af_cc ){
      if(print_afs) buffer << " " << af_case << " " << af_control;
      else buffer << " NA NA";
    }
    if(!params->build_mask && params->dosage_mode) {
      if(print_info) buffer << " " << info;
      else buffer << " NA";
    }
    buffer << " " << n ;
    if( params->af_cc )  buffer << " NA NA";
    if(print_genoc) buffer << " " << n_rr << " " << n - n_rr - n_aa << " " << n_aa ;
    else buffer << " NA NA NA";
    buffer << " " << model ;
  }

  // BETA SE
  if(print_se) buffer << ' ' << beta << ' ' << se;
  else buffer << " NA NA";

  // CHISQ PV
  if(print_pv) buffer << ' ' << chisq << ' ' << pv;
  else buffer << " NA NA";

  // extra column
  if(ipheno == params->n_pheno) {
    if(params->joint_test && (df<0)) buffer << " DF=NA\n";
    else if(params->joint_test) buffer << " DF=" << df << endl;
    else buffer << " NA\n";
  }

  return buffer.str();
}

std::string print_na_sumstats(int const& ph, int const& df, string const& header, string const& model, variant_block const* block_info, struct param const& params){
  return ( ( ph==0 ? header : "" ) + print_sum_stats_all(block_info->af1, block_info->af_case(ph), block_info->af_control(ph), block_info->n_rr,  block_info->n_aa, block_info->info1, block_info->ns1, block_info->ns_case(ph), block_info->ns_control(ph), model, -1, -1, -1, -1, false, df, &params, ph + 1) ); // pheno index is 1-based
}

// native format - single pheno
std::string print_sum_stats_single(const double& af, const double& af_case, const double& af_control, const double& info, const int& n, const int& ns_case, const int& ns_control, const string& model, const double& beta, const double& se, const double& chisq, const double& pv, const bool& test_pass, const int& df, struct param const* params){

  std::ostringstream buffer;
  bool print_afs = (af >= 0), print_info = (info >= 0), print_se = (se >= 0);
  bool print_pv = (chisq>=0) && test_pass;

  // AF N INFO TEST
  if(print_afs) buffer << af << " " ;
  else buffer << "NA " ;
  if( params->af_cc ){
    if(print_afs) buffer << af_case << " " << af_control << " ";
    else buffer << "NA NA ";
  }
  if(!params->build_mask && params->dosage_mode) {
    if(print_info) buffer << info << " ";
    else buffer << "NA ";
  }
  buffer << n ;
  if( params->af_cc )  buffer << " " << ns_case << " " << ns_control;
  buffer << " " << model << " ";

  // BETA SE
  if(print_se) buffer << beta << ' ' << se;
  else buffer << "NA NA";

  // CHISQ PV
  if(print_pv) buffer << ' ' << chisq << ' ' << pv;
  else buffer << " NA NA";

  // extra column
  vector<string> extraCol;
  if(!test_pass) extraCol.push_back("TEST_FAIL");
  if(params->joint_test && (df<0)) extraCol.push_back("DF=NA");
  else if(params->joint_test) extraCol.push_back("DF=" + to_string(df));
  buffer << " " << (extraCol.size() > 0 ? print_scsv(extraCol) : "NA") << endl;

  return buffer.str();
}


std::string print_sum_stats_htp(const double& beta, const double& se, const double& chisq, const double& lpv, const double& af, const double& info, const double& mac, const Ref<const MatrixXi>& genocounts, const int& ph, const bool& test_pass, const int& df, struct param const* params, const double& score, const double& cal_factor, const double& cal_factor_burden, const double& skat_var) {

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
  else if((params->trait_mode!=1) || ((params->trait_mode==1) && params->firth && test_pass) ){ // non-bt or firth

    if(params->trait_mode==0) // QT
      buffer << beta << "\t" << (beta - params->zcrit * se) << "\t" << (beta + params->zcrit * se) << "\t";
    else // BT (on OR scale) or CT
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
  if(af>=0)
    buffer << af << "\t";
  else
    buffer << "NA\t";

  if(mac>0) {

    // print counts in cases
    buffer << genocounts.block(0,ph,3,1).sum() << "\t" << genocounts(0,ph) << "\t" << genocounts(1,ph) << "\t" << genocounts(2,ph) << "\t";
    // print counts in controls
    if(params->trait_mode==1 || params->trait_mode==3)
      buffer << genocounts.block(3,ph,3,1).sum() << "\t" << genocounts(3,ph) << "\t" << genocounts(4,ph) << "\t" << genocounts(5,ph);
    else buffer << "NA\tNA\tNA\tNA";

  } else { // for skat/acat-type tests
    buffer << params->pheno_counts(ph, 0) << "\tNA\tNA\tNA\t";
    if(params->trait_mode==1)
      buffer << params->pheno_counts(ph, 1) << "\tNA\tNA\tNA"; 
    else buffer << "NA\tNA\tNA\tNA";
  }

  // info column
  vector<string> infoCol;
  if(print_beta){
    if(params->trait_mode && test_pass){
      infoCol.push_back( "REGENIE_BETA=" + to_string(beta) );
      infoCol.push_back( "REGENIE_SE=" + to_string(se) );
      // SPA/uncorrected logistic => also print SE from allelic OR
      if((params->trait_mode==1) && print_pv && !params->firth) infoCol.push_back( "SE=" + to_string(outse_val) );
    } else if(params->trait_mode){
      infoCol.push_back( "REGENIE_BETA=NA" );
      infoCol.push_back( "REGENIE_SE=NA");
      // SPA/uncorrected logistic => also print SE from allelic OR
      if((params->trait_mode==1) && print_pv && !params->firth) infoCol.push_back( "SE=" + to_string(outse_val) );
    } else infoCol.push_back( "REGENIE_SE=" + to_string(se) );// fot QTs
  }
  // info score
  if(!params->build_mask && params->dosage_mode && (info >= 0) ) infoCol.push_back( "INFO=" + to_string(info) );
  // mac
  if(mac>=0) infoCol.push_back( "MAC=" + to_string(mac) );
  // score test statistic 
  if(score != params->missing_value_double) infoCol.push_back( "SCORE=" + to_string(score) );
  if(skat_var != params->missing_value_double) infoCol.push_back("SKATV=" + to_string(skat_var*abs(cal_factor)));
  //if(cal_factor != -1) infoCol.push_back( "CF=" + to_string(cal_factor) );
  if(cal_factor_burden != -1) infoCol.push_back( "CF_BURDEN=" + to_string(cal_factor_burden) );
  // df
  if(params->joint_test) infoCol.push_back("DF=" + to_string(df));
  // log10P
  infoCol.push_back( "LOG10P=" + (print_pv ? to_string(lpv) : "NA") );
  // indicator for no beta printed (joint or vc tests)
  if(se<0) infoCol.push_back( "NO_BETA" );
  // print info column
  buffer << "\t" << (infoCol.size() > 0 ? print_scsv(infoCol) : "NA") << "\n";

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

    int n_corrected_cox = 0;
    if (params.trait_mode == 3) n_corrected_cox = n_corrected / 2;
    if(params.firth || params.use_SPA) {
      buffer << "Number of tests with " << (params.firth ? "Firth " : "SPA ");
      if (params.trait_mode == 3) {
        buffer << "correction : " << n_corrected_cox <<  endl;
        buffer << "Number of failed tests : (" << snp_tally.n_failed_tests << "/" << n_corrected_cox << ")\n";
      } else {
        buffer << "correction : " << n_corrected <<  endl;
        buffer << "Number of failed tests : (" << snp_tally.n_failed_tests << "/" << n_corrected << ")\n";
      }
      
    }

  }

  buffer << "Number of ignored tests due to low MAC ";
  if( params.setMinINFO ) buffer << "or info score ";
  if(params.trait_mode == 3) {
    buffer << ": " << snp_tally.n_ignored_snps/2 * params.n_tests_per_variant * params.n_pheno/2 + snp_tally.n_ignored_tests/2 * params.n_tests_per_variant << endl;
  } else {
    buffer << ": " << snp_tally.n_ignored_snps * params.n_tests_per_variant * params.n_pheno + snp_tally.n_ignored_tests * params.n_tests_per_variant << endl;
  }

  if(params.write_masks)
    buffer << "\nMasks written to : [" << files.out_file << "_masks.{bed,bim,fam}]\n";

  if(params.write_null_firth){ // store file names with null ests
    buffer << "List of files with null Firth estimates written to: [" 
      << print_null_firth_info(files, fest, params) << "]\n";
  }

  return buffer.str();

}
