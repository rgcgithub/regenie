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
#include "Step1_Models.hpp"
#include "Step2_Models.hpp"
#include "Pheno.hpp"
#include "HLM.hpp"
#include "Interaction.hpp"

using namespace std;
using namespace Eigen;
using namespace boost;
using boost::math::normal;
using boost::math::chi_squared;

void get_interaction_terms(const int& isnp, const int& thread, struct phenodt* pheno_data, struct geno_block* gblock, variant_block* snp_data, HLM& nullHLM, struct param const* params, mstream& sout){

  if(snp_data->skip_int) return; 

  data_thread* dt_thr = &(gblock->thread_data[thread]);
  MatrixXd iMat;
  SpMat Gdiag;

  if(dt_thr->is_sparse)
    Gdiag = gblock->Gmat.col(isnp).asDiagonal();

  // if rare, use HLM
  if((params->trait_mode==0) && !params->no_robust && !params->force_robust && (snp_data->mac < params->rareMAC_inter).any()){

    // get (G, G*E)
    allocate_mat(pheno_data->Hmat[thread], nullHLM.Vlin.rows(), params->interaction_istart + nullHLM.Vlin.cols());
    if(dt_thr->is_sparse)
      pheno_data->Hmat[thread].rightCols(nullHLM.Vlin.cols()) = Gdiag * nullHLM.Vlin;
    else
      pheno_data->Hmat[thread].rightCols(nullHLM.Vlin.cols()) = (nullHLM.Vlin.array().colwise() * gblock->Gmat.col(isnp).array()).matrix();

    // add main effects for G_E if specified (not re-scaled)
    if(!params->gwas_condtl) {
      pheno_data->Hmat[thread].leftCols(pheno_data->interaction_cov.cols()) = pheno_data->interaction_cov;
      if(params->add_homdev)
        pheno_data->Hmat[thread].col(pheno_data->interaction_cov.cols()) = pheno_data->interaction_homdev;
    }

    snp_data->fitHLM = true;
    return;

  }

  if(dt_thr->is_sparse) 
    iMat = Gdiag * pheno_data->interaction_cov;
  else
    iMat = (pheno_data->interaction_cov.array().colwise() * gblock->Gmat.col(isnp).array()).matrix();
  //if(isnp==0) cerr << iMat.topRows(10) << endl; 

  // remove covariate effects
  snp_data->skip_int = !residualize_matrix(iMat, pheno_data->scf_i[thread], pheno_data->new_cov.leftCols( params->ncov + (params->blup_cov && (params->trait_mode == 1) ? -1 : 0)), params->n_analyzed, params->numtol);
  if(snp_data->skip_int) return;

  // start filling matrix with C*G terms (if not condtl, also add residual of G_E)
  pheno_data->Hmat[thread].resize(pheno_data->interaction_cov.rows(), params->ncov_interaction + params->interaction_istart + 1);
  pheno_data->Hmat[thread].rightCols(iMat.cols()) = iMat;
  if(!params->gwas_condtl) pheno_data->Hmat[thread].leftCols(pheno_data->interaction_cov_res.cols()) = pheno_data->interaction_cov_res;

}


/// Interaction tests
void apply_interaction_tests(const int& index, const int& isnp, const int& thread, const Ref<const MatrixXd>& res, const Ref<const RowVectorXd>& sd_yres, string const& model_type, string const& test_string, struct phenodt* pheno_data, HLM& nullHLM, struct filter const* filters, struct in_files* files, struct geno_block* gblock, variant_block* snp_data, vector<snp> const& snpinfo, struct ests* m_ests, struct f_ests* fest, struct param const* params, mstream& sout){

  if(snp_data->skip_int) return;

  if(params->trait_mode==1)
    apply_interaction_tests_bt(index, isnp, thread, model_type, test_string, pheno_data, filters, files, gblock, snp_data, snpinfo, m_ests, fest, params, sout);
  else if((params->trait_mode==0) && snp_data->fitHLM)
    apply_interaction_tests_HLM(index, isnp, thread, res, sd_yres, model_type, test_string, pheno_data, nullHLM, filters, files, gblock, snp_data, snpinfo, params, sout);
  else if(params->trait_mode==0)
    apply_interaction_tests_qt(index, isnp, thread, res, sd_yres, model_type, test_string, pheno_data, filters, files, gblock, snp_data, snpinfo, params, sout);

}

void apply_interaction_tests_qt(const int& index, const int& isnp, const int& thread, const Ref<const MatrixXd>& res, const Ref<const RowVectorXd>& sd_yres, string const& model_type, string const& test_string, struct phenodt* pheno_data, struct filter const* filters, struct in_files* files, struct geno_block* gblock, variant_block* snp_data, vector<snp> const& snpinfo, struct param const* params, mstream& sout){

  int beg = params->interaction_istart;
  string df_str = to_string(1+params->ncov_interaction);

  // fill rest of matrix [ G, C*G ]
  pheno_data->Hmat[thread].col(beg) = gblock->Gmat.col(isnp);
  //cerr << pheno_data->Hmat[thread].topRows(3) << endl; exit(-1);

  // pre-compute Z = (M^tM)^(-1) for all phenos 
  SelfAdjointEigenSolver<MatrixXd> esM(pheno_data->Hmat[thread].transpose() * pheno_data->Hmat[thread]);
  if( esM.eigenvalues().minCoeff() < params->numtol ) return;
  MatrixXd Z = esM.eigenvectors() * esM.eigenvalues().cwiseInverse().asDiagonal() * esM.eigenvectors().transpose();

  // get leverage h = diag( M * Z * M^t )
  VectorXd hvec = ((pheno_data->Hmat[thread] * Z).array() * pheno_data->Hmat[thread].array()).matrix().rowwise().sum();

  // estimates for all traits
  ArrayXd hc3, hc4;
  MatrixXd tau = Z * pheno_data->Hmat[thread].transpose() * res;
  MatrixXd e_sq = ((res - pheno_data->Hmat[thread] * tau).array().square() * pheno_data->masked_indivs.array().cast<double>()).matrix();
  if(!params->no_robust){
    hc3 = (1 - hvec.array()).square(); 
    if(params->force_hc4)
      hc4 = (1 - hvec.array()).pow( (pheno_data->Hmat[thread].rows() * hvec.array() / pheno_data->Hmat[thread].cols()).min(4) );
    //if(isnp==0) cerr << "tau=\n" << tau << "\n\n ei=\n" << e_sq.topRows(3) << endl;
  }

  chi_squared chisqI(params->ncov_interaction);
  chi_squared chisqK(params->ncov_interaction+1);
  double logp, gscale, tstat, sehat;
  string head = "", stmp;
  MatrixXd Vmat, Vinv;
  ArrayXd iscale, cscale;

  // for output
  if(!params->htp_out) head = print_sum_stats_head(index, snpinfo);

  for(int i = 0; i < params->n_pheno; ++i ) {
    if( !params->pheno_pass(i) ) continue;

    if( snp_data->ignored_trait(i) ) continue;

    std::ostringstream buffer;

    MapArXd bhat (tau.col(i).data(), tau.rows(), 1);
    gscale = pheno_data->scale_Y(i) * sd_yres(i) / snp_data->scale_fac;
    iscale = pheno_data->scale_Y(i) * sd_yres(i) / pheno_data->scf_i[thread];
    if(!params->gwas_condtl) cscale = pheno_data->scale_Y(i) * sd_yres(i) / pheno_data->scl_inter_X;

    // using sandwich estimator
    if(params->no_robust) // model-based
      Vmat = e_sq.col(i).sum() / (pheno_data->Neff(i) - params->ncov_analyzed - Z.cols()) * Z; // s^2*(XtX)^-1
    else if(params->force_hc4 && (snp_data->mac(i) <= params->rareMAC_inter)) // HC4
      Vmat = Z * pheno_data->Hmat[thread].transpose() * (e_sq.col(i).array() / hc4).matrix().asDiagonal() * pheno_data->Hmat[thread] * Z;
    else // HC3
      Vmat = Z * pheno_data->Hmat[thread].transpose() * (e_sq.col(i).array() / hc3).matrix().asDiagonal() * pheno_data->Hmat[thread] * Z;
    //if(index==500) {cerr << "\nZ:\n" << Z << "\nV=\n" << Vmat ; exit(-1);}

    // print cov(beta) (rescale)
    if(params->print_vcov && !params->gwas_condtl){
      Files fout;
      fout.openForWrite(files->out_file + "_" + files->pheno_names[i] + "_" + filters->interaction_cov + "_" + snpinfo[index].ID + ".vcov", sout);
      MatrixXd scvec (pheno_data->Hmat[thread].cols(), 1);
      scvec.col(0).array().head( cscale.size() ) = cscale;
      scvec(beg, 0) = gscale;
      scvec.col(0).array().tail( iscale.size() ) = iscale;
      IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","");
      fout << (scvec.col(0).asDiagonal() * Vmat * scvec.col(0).asDiagonal()).format(Fmt); 
      fout.closeFile();
    }

    ///////////////////////
    // print main effect of G_E
    if(beg > 0){
      for(int j = 0; j < beg; j++){ 
        tstat = bhat(j) * bhat(j) / Vmat(j,j);
        sehat = sqrt(Vmat(j,j)) * cscale(j);
        get_logp(logp, tstat);
        if(params->interaction_cat)
          stmp="-INT_" + filters->interaction_cov + "=" + params->interaction_lvl_names[j];
        else if(params->add_homdev && (j != 0))
          stmp="-INT_" + filters->interaction_cov + "-HOM"; // G_E>=1.5
        else
          stmp="-INT_" + filters->interaction_cov; // single cov

        // print sum_stats
        if(params->htp_out) 
          buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bhat(j) * cscale(j), sehat, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1, params);
        else 
          buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, bhat(j) * cscale(j), sehat, tstat, logp, true, 1, params, (i+1));
      }
    }

    ///////////////////////
    //////  marginal test
    // T, beta, se & pv
    tstat = bhat(beg) * bhat(beg) / Vmat(beg,beg);
    sehat = sqrt(Vmat(beg,beg)) * gscale;
    get_logp(logp, tstat);
    stmp="-INT_SNP";

    // print sum_stats
    if(params->htp_out) 
      buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bhat(beg) * gscale, sehat, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1, params);
    else 
      buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, bhat(beg) * gscale, sehat, tstat, logp, true, 1, params, (i+1));


    ///////////////////////
    //////  interaction tests
    if(params->ncov_interaction > 1){ 

      // print effects for each interaction term
      for(int j = 0; j < params->ncov_interaction; j++){ 
        tstat = bhat(beg+1+j) * bhat(beg+1+j) / Vmat(beg+1+j,beg+1+j);
        sehat = sqrt(Vmat(beg+1+j,beg+1+j)) * iscale(j);
        get_logp(logp, tstat);
        stmp="-INT_SNPx" + filters->interaction_cov + "=" + params->interaction_lvl_names[j];
        // print sum_stats
        if(params->htp_out) 
          buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bhat(beg+1+j) * iscale(j), sehat, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1, params);
        else 
          buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, bhat(beg+1+j) * iscale(j), sehat, tstat, logp, true, 1, params, (i+1));
      }

      // joint test for interaction terms
      // T, beta, se & pv
      Vinv = Vmat.block(beg+1,beg+1,params->ncov_interaction,params->ncov_interaction).inverse();
      tstat = fabs( (bhat.tail(params->ncov_interaction).matrix().transpose() * Vinv * bhat.tail(params->ncov_interaction).matrix()).sum() );
      logp = max(params->nl_dbl_dmin, cdf(complement(chisqI, tstat)));
      logp = -log10( logp );
      stmp="-INT_SNPx" + filters->interaction_cov;
      // print sum_stats
      if(params->htp_out) 
        buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(-1, -1, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, params->ncov_interaction, params);
      else 
        buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, -1, -1, tstat, logp, true, params->ncov_interaction, params, (i+1));

    } else {
      // T, beta, se & pv
      tstat = bhat(beg+1) * bhat(beg+1) / Vmat(beg+1,beg+1);
      sehat = sqrt(Vmat(beg+1,beg+1)) * iscale(0);
      get_logp(logp, tstat);
      if(params->interaction_cat)
        stmp="-INT_SNPx" + filters->interaction_cov + "=" + params->interaction_lvl_names[0];
      else
        stmp="-INT_SNPx" + filters->interaction_cov; // single cov

      // print sum_stats
      if(params->htp_out) 
        buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bhat(beg+1) * iscale(0), sehat, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1, params);
      else 
        buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, bhat(beg+1) * iscale(0), sehat, tstat, logp, true, 1, params, (i+1));
    }

    ///////////////////////
    //////  joint test for G and C*G
    Vinv = Vmat.block(beg,beg,params->ncov_interaction+1,params->ncov_interaction+1).inverse();

    // T & pv
    tstat = fabs( (bhat.tail(params->ncov_interaction+1).matrix().transpose() * Vinv * bhat.tail(params->ncov_interaction+1).matrix()).sum() );
    logp = max(params->nl_dbl_dmin, cdf(complement(chisqK, tstat)));
    logp = -log10( logp );
    stmp="-INT_" + df_str + "DF";

    // print sum_stats
    if(params->htp_out) 
      buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(-1, -1, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1+params->ncov_interaction, params);
    else 
      buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, -1, -1, tstat, logp, true, 1+params->ncov_interaction, params, (i+1));

    //if(isnp==0 & i==0) cerr << endl << buffer.str() << endl;
    snp_data->sum_stats[i].append( buffer.str() );

  }

}


void apply_interaction_tests_HLM(const int& index, const int& isnp, const int& thread, const Ref<const MatrixXd>& res, const Ref<const RowVectorXd>& sd_yres, string const& model_type, string const& test_string, struct phenodt* pheno_data, HLM& nullHLM, struct filter const* filters, struct in_files* files, struct geno_block* gblock, variant_block* snp_data, vector<snp> const& snpinfo, struct param const* params, mstream& sout){

  int beg = params->interaction_istart;
  string df_str = to_string(1+params->ncov_interaction);

  //cerr << pheno_data->Hmat[thread].topRows(20) << endl; exit(-1);

  chi_squared chisqI(params->ncov_interaction);
  chi_squared chisqK(params->ncov_interaction+1);
  double logp, tstat, sehat;
  string head = "", stmp;
  ArrayXd bhat;
  MatrixXd Xres, Vmat, Vinv;
  allocate_mat(Xres, pheno_data->Hmat[thread].rows(), pheno_data->Hmat[thread].cols());

  // for output
  if(!params->htp_out) head = print_sum_stats_head(index, snpinfo);

  for(int i = 0; i < params->n_pheno; ++i ) {
    if( !params->pheno_pass(i) ) continue;

    if( snp_data->ignored_trait(i) ) continue;

    std::ostringstream buffer;

    // get the residuals using null HLM model
    nullHLM.residualize(i, pheno_data->Hmat[thread], Xres);

    // OLS (V is different for each trait) - sigma^2=1
    SelfAdjointEigenSolver<MatrixXd> esM(Xres.transpose() * Xres);
    if( esM.eigenvalues().minCoeff() < params->numtol ) return;
    Vmat = esM.eigenvectors() * esM.eigenvalues().cwiseInverse().asDiagonal() * esM.eigenvectors().transpose();
    bhat = Vmat * (Xres.transpose() * nullHLM.yres.col(i));
    //cerr << "\n" << bhat << "\n\n" << Vmat.array().sqrt().matrix() << "\n\n"; exit(-1);

    // print cov(beta) (rescale)
    if(params->print_vcov && !params->gwas_condtl){
      Files fout;
      fout.openForWrite( files->out_file + "_" + files->pheno_names[i] + "_" + filters->interaction_cov + "_" + snpinfo[index].ID + ".vcov", sout);
      IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","");
      fout << Vmat.format(Fmt); 
      fout.closeFile();
    }

    ///////////////////////
    // print main effect of G_E
    if(beg > 0){
      for(int j = 0; j < beg; j++){ 
        tstat = bhat(j) * bhat(j) / Vmat(j,j);
        sehat = sqrt(Vmat(j,j));
        get_logp(logp, tstat);
        if(params->interaction_cat)
          stmp="-INT_" + filters->interaction_cov + "=" + params->interaction_lvl_names[j];
        else if(params->add_homdev && (j != 0))
          stmp="-INT_" + filters->interaction_cov + "-HOM"; // G_E>=1.5
        else
          stmp="-INT_" + filters->interaction_cov; // single cov

        // print sum_stats
        if(params->htp_out) 
          buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bhat(j), sehat, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1, params);
        else 
          buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, bhat(j), sehat, tstat, logp, true, 1, params, (i+1));
      }
    }

    ///////////////////////
    //////  marginal test
    // T, beta, se & pv
    tstat = bhat(beg) * bhat(beg) / Vmat(beg,beg);
    sehat = sqrt(Vmat(beg,beg));
    get_logp(logp, tstat);
    stmp="-INT_SNP";

    // print sum_stats
    if(params->htp_out) 
      buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bhat(beg), sehat, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1, params);
    else 
      buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, bhat(beg), sehat, tstat, logp, true, 1, params, (i+1));


    ///////////////////////
    //////  interaction tests
    if(params->ncov_interaction > 1){ 

      // print effects for each interaction term
      for(int j = 0; j < params->ncov_interaction; j++){ 
        tstat = bhat(beg+1+j) * bhat(beg+1+j) / Vmat(beg+1+j,beg+1+j);
        sehat = sqrt(Vmat(beg+1+j,beg+1+j));
        get_logp(logp, tstat);
        stmp="-INT_SNPx" + filters->interaction_cov + "=" + params->interaction_lvl_names[j];
        // print sum_stats
        if(params->htp_out) 
          buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bhat(beg+1+j), sehat, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1, params);
        else 
          buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, bhat(beg+1+j), sehat, tstat, logp, true, 1, params, (i+1));
      }

      // joint test for interaction terms
      // T, beta, se & pv
      Vinv = Vmat.block(beg+1,beg+1,params->ncov_interaction,params->ncov_interaction).inverse();
      tstat = fabs( (bhat.tail(params->ncov_interaction).matrix().transpose() * Vinv * bhat.tail(params->ncov_interaction).matrix()).sum() );
      logp = max(params->nl_dbl_dmin, cdf(complement(chisqI, tstat)));
      logp = -log10( logp );
      stmp="-INT_SNPx" + filters->interaction_cov;
      // print sum_stats
      if(params->htp_out) 
        buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(-1, -1, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, params->ncov_interaction, params);
      else 
        buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, -1, -1, tstat, logp, true, params->ncov_interaction, params, (i+1));

    } else {
      // T, beta, se & pv
      tstat = fabs( bhat(beg+1) * bhat(beg+1) / Vmat(beg+1,beg+1) );
      sehat = sqrt(Vmat(beg+1,beg+1));
      get_logp(logp, tstat);
      if(params->interaction_cat)
        stmp="-INT_SNPx" + filters->interaction_cov + "=" + params->interaction_lvl_names[0];
      else
        stmp="-INT_SNPx" + filters->interaction_cov; // single cov

      // print sum_stats
      if(params->htp_out) 
        buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bhat(beg+1), sehat, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1, params);
      else 
        buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, bhat(beg+1), sehat, tstat, logp, true, 1, params, (i+1));
    }

    ///////////////////////
    //////  joint test for G and C*G
    Vinv = Vmat.block(beg,beg,params->ncov_interaction+1,params->ncov_interaction+1).inverse();

    // T & pv
    tstat = fabs( (bhat.tail(params->ncov_interaction+1).matrix().transpose() * Vinv * bhat.tail(params->ncov_interaction+1).matrix()).sum() );
    logp = max(params->nl_dbl_dmin, cdf(complement(chisqK, tstat)));
    logp = -log10( logp );
    stmp="-INT_" + df_str + "DF";

    // print sum_stats
    if(params->htp_out) 
      buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(-1, -1, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1+params->ncov_interaction, params);
    else 
      buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, -1, -1, tstat, logp, true, 1+params->ncov_interaction, params, (i+1));

    //cerr << endl << buffer.str() << endl; exit(-1);
    snp_data->sum_stats[i].append( buffer.str() );
  }

}



void apply_interaction_tests_bt(const int& index, const int& isnp, const int& thread, string const& model_type, string const& test_string, struct phenodt* pheno_data, struct filter const* filters, struct in_files* files,struct geno_block* gblock, variant_block* snp_data, vector<snp> const& snpinfo, struct ests* m_ests, struct f_ests* firth_est, struct param const* params, mstream& sout){

  int beg = params->interaction_istart;
  int np = pheno_data->Hmat[thread].cols() - beg;
  string df_str = to_string(np);

  // fill rest of matrix [ G, C*G ]
  // remove covariate effects (not done for marginal test)
  // use projection from linear regression so only needs to be done once for all traits
  // with step 1 preds as cov, ignore last column of X
  residualize_geno(isnp, thread, snp_data, true, pheno_data->new_cov.leftCols( params->ncov + (params->blup_cov ? -1 : 0)), gblock, params);
  pheno_data->Hmat[thread].col(beg) = gblock->Gmat.col(isnp);
  //cerr << pheno_data->Hmat[thread].topRows(3) << endl;

  chi_squared chisqI(np-1);
  chi_squared chisqK(np);
  double logp, tstat, sehat;
  double lpfirth = -log10( params->alpha_pvalue );
  string head = "", stmp;

  // for output
  if(!params->htp_out) head = print_sum_stats_head(index, snpinfo);

  // for logistic regression
  ArrayXd bhat, etavec, pivec, hvec;
  MatrixXd WX, Vmat, XWX_inv, V_robust;

  for(int i = 0; i < params->n_pheno; ++i ){
    if( !params->pheno_pass(i) ) continue;

    if( snp_data->ignored_trait(i) ) continue;
      
    std::ostringstream buffer, buffer_int;

    MapArXd Y (pheno_data->phenotypes_raw.col(i).data(), pheno_data->phenotypes_raw.rows());
    MapArXb mask (pheno_data->masked_indivs.col(i).data(), pheno_data->masked_indivs.rows());
    MapcArXd offset (m_ests->offset_nullreg.col(i).data(), m_ests->offset_nullreg.rows());

    // starting values
    bhat = ArrayXd::Zero(np+beg);
    etavec = mask.select(offset, 0);
    get_pvec(pivec, etavec, params->numtol_eps);

    if(!fit_logistic(Y, pheno_data->Hmat[thread], offset, mask, pivec, etavec, bhat, params, sout))
      continue; // no results for trait
    /*else if( (mask && (pivec < params->numtol_eps || pivec > 1 - params->numtol_eps)).count() > 0 )
      sout << "\n     WARNING: Fitted probabilities numerically 0/1 occurred (phenotype #" << files->pheno_names[i] <<").";*/
    //cerr << bhat << endl;

    // get cov(beta)
    WX = ( pheno_data->Hmat[thread].array().colwise() * mask.select(pivec * (1 - pivec), 0).sqrt() ).matrix();
    SelfAdjointEigenSolver<MatrixXd> esM(WX.transpose() * WX);
    if( esM.eigenvalues().minCoeff() < params->numtol ) continue;
    if( !params->force_robust && 
        (params->no_robust || (snp_data->mac(i) <= params->rareMAC_inter)) ) // using model-based estimator [robust is inflated for ultra-rare]
      Vmat = esM.eigenvectors() * esM.eigenvalues().cwiseInverse().asDiagonal() * esM.eigenvectors().transpose();
    else { // robust se
      XWX_inv = esM.eigenvectors() * esM.eigenvalues().cwiseInverse().asDiagonal() * esM.eigenvectors().transpose();
      hvec = ((WX * XWX_inv).array() * WX.array()).rowwise().sum();
      V_robust = pheno_data->Hmat[thread].transpose() * mask.select((Y - pivec)/(1-hvec), 0).square().matrix().asDiagonal() * pheno_data->Hmat[thread];
      Vmat = XWX_inv * V_robust * XWX_inv;
      if(params->debug) cerr << "h:\n" << hvec.minCoeff() << " - " << hvec.maxCoeff() << "\n\nb:\n" << bhat << "\n\nV:\n" << Vmat << "\n\nXWXinv:\n" << XWX_inv << "\n\nVh:\n"<< V_robust << "\n\n";
    }
    if( Vmat.diagonal().minCoeff() < 0 ) continue; // if robust SE computation fails
    if(snp_data->flipped) bhat *= -1;

    ///////////////////////
    //////  interaction tests
    bool use_firth = false;
    if(params->ncov_interaction > 1){ 

      // print effects for each interaction term
      for(int j = 0; j < params->ncov_interaction; j++){
        // T, beta, se & pv
        tstat = fabs( bhat(beg+1+j) * bhat(beg+1+j) / Vmat(beg+1+j,beg+1+j) );
        sehat = sqrt(Vmat(beg+1+j,beg+1+j));
        get_logp(logp, tstat);

        // if firth, check pvalue <= thresh
        use_firth = params->firth && (logp >= lpfirth);
        if(use_firth) break;

        stmp="-INT_SNPx" + filters->interaction_cov + "=" + params->interaction_lvl_names[j];

        // print sum_stats
        if(params->htp_out) 
          buffer_int << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bhat(beg+1+j)/pheno_data->scf_i[thread](j), sehat/pheno_data->scf_i[thread](j), tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1, params);
        else 
          buffer_int << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, bhat(beg+1+j)/pheno_data->scf_i[thread](j), sehat/pheno_data->scf_i[thread](j), tstat, logp, true, 1, params, (i+1));
      }

      // switch to firth
      if(use_firth){
        //cerr << i+1 << "   " << snpinfo[index].ID << endl;
        snp_data->sum_stats[i].append( 
            apply_interaction_tests_firth(index, isnp, thread, i, model_type, test_string, pheno_data, filters, files, gblock, snp_data, snpinfo, m_ests, firth_est, params, sout)
            );
        continue; // go to next trait
      }

      /// joint test for interaction
      // T & pv
      tstat = fabs( (bhat.tail(params->ncov_interaction).matrix().transpose() * Vmat.block(beg+1,beg+1,params->ncov_interaction,params->ncov_interaction).inverse() * bhat.tail(params->ncov_interaction).matrix()).sum() );
      logp = max(params->nl_dbl_dmin, cdf(complement(chisqI, tstat)));
      logp = -log10( logp );
      stmp="-INT_SNPx" + filters->interaction_cov;

      // print sum_stats
      if(params->htp_out) 
        buffer_int << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(-1, -1, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, np-1, params);
      else 
        buffer_int << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, -1, -1, tstat, logp, true, np-1, params, (i+1));

    } else {

      // T, beta, se & pv
      tstat = fabs( bhat(beg+1) * bhat(beg+1) / Vmat(beg+1,beg+1) );
      sehat = sqrt(Vmat(beg+1,beg+1));
      get_logp(logp, tstat);

      // if firth, check pvalue <= thresh
      use_firth = params->firth && (logp >= lpfirth);
      // switch to firth
      if(use_firth){
        snp_data->sum_stats[i].append( 
            apply_interaction_tests_firth(index, isnp, thread, i, model_type, test_string, pheno_data, filters, files, gblock, snp_data, snpinfo, m_ests, firth_est, params, sout)
            );
        continue; // go to next trait
      }

      if(params->interaction_cat)
        stmp="-INT_SNPx" + filters->interaction_cov + "=" + params->interaction_lvl_names[0];
      else
        stmp="-INT_SNPx" + filters->interaction_cov; // single cov

      // print sum_stats
      if(params->htp_out) 
        buffer_int << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bhat(beg+1)/pheno_data->scf_i[thread](0), sehat/pheno_data->scf_i[thread](0), tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1, params);
      else 
        buffer_int << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, bhat(beg+1)/pheno_data->scf_i[thread](0), sehat/pheno_data->scf_i[thread](0), tstat, logp, true, 1, params, (i+1));

    }

    // print cov(beta) (rescale)
    if(params->print_vcov && !params->gwas_condtl){
      Files fout;
      fout.openForWrite( files->out_file + "_" + files->pheno_names[i] + "_" + filters->interaction_cov + "_" + snpinfo[index].ID + ".vcov", sout);
      MatrixXd scvec (pheno_data->Hmat[thread].cols(), 1);
      scvec.col(0).array().head( pheno_data->scl_inter_X.size() ) = 1/pheno_data->scl_inter_X;
      scvec(beg, 0) = 1/snp_data->scale_fac;
      scvec.col(0).array().tail( pheno_data->scf_i[thread].size() ) = 1/pheno_data->scf_i[thread];
      IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","");
      fout << (scvec.col(0).asDiagonal() * Vmat * scvec.col(0).asDiagonal()).format(Fmt); 
      fout.closeFile();
    }

    ///////////////////////
    // print main effect of G_E
    if(beg > 0){
      for(int j = 0; j < beg; j++){ 
        tstat = bhat(j) * bhat(j) / Vmat(j,j);
        sehat = sqrt(Vmat(j,j));
        get_logp(logp, tstat);
        if(params->interaction_cat)
          stmp="-INT_" + filters->interaction_cov + "=" + params->interaction_lvl_names[j];
        else if(params->int_add_esq && (j != 0))
          stmp="-INT_" + filters->interaction_cov + "^2"; // G_E^2
        else if(params->add_homdev && (j != 0))
          stmp="-INT_" + filters->interaction_cov + "-HOM"; // G_E>=1.5
        else
          stmp="-INT_" + filters->interaction_cov; // single cov

        // print sum_stats
        if(params->htp_out) 
          buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bhat(j)/pheno_data->scl_inter_X(j), sehat/pheno_data->scl_inter_X(j), tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1, params);
        else 
          buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, bhat(j)/pheno_data->scl_inter_X(j), sehat/pheno_data->scl_inter_X(j), tstat, logp, true, 1, params, (i+1));
      }
    }

    ///////////////////////
    //////  marginal test
    // T, beta, se & pv
    tstat = bhat(beg) * bhat(beg) / Vmat(beg,beg);
    sehat = sqrt(Vmat(beg,beg));
    get_logp(logp, tstat);
    stmp="-INT_SNP";

    // print sum_stats
    if(params->htp_out) 
      buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bhat(beg)/ snp_data->scale_fac, sehat/ snp_data->scale_fac, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, 1, params);
     else 
      buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, bhat(beg)/ snp_data->scale_fac, sehat/ snp_data->scale_fac, tstat, logp, true, 1, params, (i+1));

    // add interaction test results
    buffer << buffer_int.str();

    ///////////////////////
    //////  joint test for G and C*G
    // T & pv
    if(beg!=0)
      tstat = fabs( (bhat.tail(params->ncov_interaction+1).matrix().transpose() * Vmat.block(beg,beg,params->ncov_interaction+1,params->ncov_interaction+1).inverse() * bhat.tail(params->ncov_interaction+1).matrix()).sum() );
    else
      tstat = fabs( (bhat.matrix().transpose() * Vmat.inverse() * bhat.matrix()).sum() );
    logp = max(params->nl_dbl_dmin, cdf(complement(chisqK, tstat)));
    logp = -log10( logp );
    stmp="-INT_" + df_str + "DF";

    // print sum_stats
    if(params->htp_out) 
      buffer << print_sum_stats_head_htp(index, files->pheno_names[i], model_type + stmp, snpinfo, params) << print_sum_stats_htp(-1, -1, tstat, logp, snp_data->af(i), snp_data->info(i), snp_data->mac(i), snp_data->genocounts, i, true, np, params);
    else 
      buffer << (!params->split_by_pheno && (i>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(i) : snp_data->af1),snp_data->af_case(i),snp_data->af_control(i), (params->split_by_pheno ? snp_data->info(i) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(i) : snp_data->ns1),snp_data->ns_case(i),snp_data->ns_control(i), test_string + stmp, -1, -1, tstat, logp, true, np, params, (i+1));

    //if(isnp==0 & i==0) cerr << endl << buffer.str() << endl;
    snp_data->sum_stats[i].append( buffer.str() );

  }

}

std::string apply_interaction_tests_firth(const int& index, const int& isnp, const int& thread, const int& ipheno, string const& model_type, string const& test_string, struct phenodt* pheno_data, struct filter const* filters, struct in_files const* files, struct geno_block* gblock, variant_block* snp_data, vector<snp> const& snpinfo, struct ests const* m_ests, struct f_ests const* firth_est, struct param const* params, mstream& sout){

  int beg = params->interaction_istart;
  int np = pheno_data->Hmat[thread].cols() - beg;
  string df_str = to_string(np);
  snp_data->is_corrected_inter(ipheno) = true;

  chi_squared chisqI(np-1);
  chi_squared chisqK(np);
  double logp, tstat;
  double bsign = snp_data->flipped ? -1 : 1;
  string head = "", stmp;

  // for output
  if(!params->htp_out) head = print_sum_stats_head(index, snpinfo);

  // for firth regression
  double dev, dev_s, se_val;
  ArrayXd bhat, bhat_s, se, se_s, etavec, pivec;
  std::ostringstream buffer, buffer_joint;

  MapArXd Y (pheno_data->phenotypes_raw.col(ipheno).data(), pheno_data->phenotypes_raw.rows());
  MapArXb mask (pheno_data->masked_indivs.col(ipheno).data(), pheno_data->masked_indivs.rows());
  MapcArXd offset (firth_est->cov_blup_offset.col(ipheno).data(), firth_est->cov_blup_offset.rows());

  // fit full model with G & C*G
  // starting values
  bhat = ArrayXd::Zero(np+beg);
  if(beg!=0) {// fit null model with only G_E
    if(!fit_firth(ipheno, Y, pheno_data->Hmat[thread], offset, mask, pivec, etavec, bhat, se, beg, dev_s, false, tstat, params->maxstep_null, params->niter_max_firth_null, params->numtol_firth, params))
      return "";
  }
  if(!fit_firth(ipheno, Y, pheno_data->Hmat[thread], offset, mask, pivec, etavec, bhat, se, np+beg, dev, true, tstat, params->maxstep, params->niter_max_firth, params->numtol_firth, params))
    return "";
  //if(isnp==0 & ipheno==0)cerr << bhat << "\n\n" << se << endl;exit(-1);

  // print cov(beta) (rescale)
  if(params->print_vcov && !params->gwas_condtl){
    Files fout;
    fout.openForWrite( files->out_file + "_" + files->pheno_names[ipheno] + "_" + filters->interaction_cov + "_" + snpinfo[index].ID + ".vcov", sout);
    ArrayXd wvec = mask.select( ( pivec * (1 - pivec) ).sqrt(), 0);
    MatrixXd XtW = pheno_data->Hmat[thread].transpose() * wvec.matrix().asDiagonal();
    ColPivHouseholderQR<MatrixXd> qr(XtW * XtW.transpose());
    MatrixXd scvec (pheno_data->Hmat[thread].cols(), 1);
    scvec.col(0).array().head( pheno_data->scl_inter_X.size() ) = 1/pheno_data->scl_inter_X;
    scvec(beg, 0) = 1/snp_data->scale_fac;
    scvec.col(0).array().tail( pheno_data->scf_i[thread].size() ) = 1/pheno_data->scf_i[thread];
    IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","");
    fout << (scvec.col(0).asDiagonal() * qr.inverse() * scvec.col(0).asDiagonal()).format(Fmt); 
    fout.closeFile();
  }

  ///////////////////////
  //////  GxG: G_E main effect
  if(!params->gwas_condtl){
    for(int j = 0; j < beg; j++){
      if(params->interaction_cat)
        stmp="-INT_" + filters->interaction_cov + "=" + params->interaction_lvl_names[j];
      else if(params->int_add_esq && (j != 0))
        stmp="-INT_" + filters->interaction_cov + "^2"; // G_E^2
      else if(params->add_homdev && (j != 0))
        stmp="-INT_" + filters->interaction_cov + "-HOM"; // G_E>=1.5
      else
        stmp="-INT_" + filters->interaction_cov; // single cov

      // print sum_stats
      if(params->htp_out) 
        buffer << print_sum_stats_head_htp(index, files->pheno_names[ipheno], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bhat(j)/pheno_data->scl_inter_X(j), se(j)/pheno_data->scl_inter_X(j), -1, -1, snp_data->af(ipheno), snp_data->info(ipheno), snp_data->mac(ipheno), snp_data->genocounts, ipheno, true, 1, params);
      else 
        buffer << (!params->split_by_pheno && (ipheno>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(ipheno) : snp_data->af1),snp_data->af_case(ipheno),snp_data->af_control(ipheno), (params->split_by_pheno ? snp_data->info(ipheno) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(ipheno) : snp_data->ns1),snp_data->ns_case(ipheno),snp_data->ns_control(ipheno), test_string + stmp, bhat(j)/pheno_data->scl_inter_X(j), se(j)/pheno_data->scl_inter_X(j), -1, -1, true, 1, params, (ipheno+1));
    }
  }

  /////////////// joint test
  // pv
  if(beg!=0) tstat = dev_s - dev;
  if(tstat < 0) return "";
  logp = max(params->nl_dbl_dmin, cdf(complement(chisqK, tstat)));
  logp = -log10( logp );
  stmp="-INT_" + df_str + "DF";

  // print sum_stats
  if(params->htp_out) 
    buffer_joint << print_sum_stats_head_htp(index, files->pheno_names[ipheno], model_type + stmp, snpinfo, params) << print_sum_stats_htp(-1, -1, tstat, logp, snp_data->af(ipheno), snp_data->info(ipheno), snp_data->mac(ipheno), snp_data->genocounts, ipheno, true, np, params);
  else 
    buffer_joint << (!params->split_by_pheno && (ipheno>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(ipheno) : snp_data->af1),snp_data->af_case(ipheno),snp_data->af_control(ipheno), (params->split_by_pheno ? snp_data->info(ipheno) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(ipheno) : snp_data->ns1),snp_data->ns_case(ipheno),snp_data->ns_control(ipheno), test_string + stmp, -1, -1, tstat, logp, true, np, params, (ipheno+1));

  // get lrt values for each test
  ///////////////////////
  //////  marginal test
  pheno_data->Hmat[thread].col(beg).swap( pheno_data->Hmat[thread].rightCols(1) ); // put G in last column
  bhat_s = bhat;
  bhat_s(beg) = bhat.tail(1)(0);
  bhat_s.tail(1)(0) = 0;

  if(!fit_firth(ipheno, Y, pheno_data->Hmat[thread], offset, mask, pivec, etavec, bhat_s, se_s, beg+np-1, dev_s, false, tstat, params->maxstep, params->niter_max_firth, params->numtol_firth, params))
    return "";
  // T, beta, se & pv
  tstat = dev_s - dev;
  if(tstat < 0) return "";
  get_logp(logp, tstat);
  stmp="-INT_SNP";
  if( params->back_correct_se && (tstat > 0) )
    se_val = fabs(bhat(beg)) / sqrt(tstat);
  else
    se_val = se(beg);

  // print sum_stats
  if(params->htp_out) 
    buffer << print_sum_stats_head_htp(index, files->pheno_names[ipheno], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bsign * bhat(beg)/snp_data->scale_fac, se_val/snp_data->scale_fac, tstat, logp, snp_data->af(ipheno), snp_data->info(ipheno), snp_data->mac(ipheno), snp_data->genocounts, ipheno, true, 1, params);
  else 
    buffer << (!params->split_by_pheno && (ipheno>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(ipheno) : snp_data->af1),snp_data->af_case(ipheno),snp_data->af_control(ipheno), (params->split_by_pheno ? snp_data->info(ipheno) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(ipheno) : snp_data->ns1),snp_data->ns_case(ipheno),snp_data->ns_control(ipheno), test_string + stmp, bsign * bhat(beg)/ snp_data->scale_fac, se_val/ snp_data->scale_fac, tstat, logp, true, 1, params, (ipheno+1));

  ///////////////////////
  //////  interaction tests
  pheno_data->Hmat[thread].col(beg).swap( pheno_data->Hmat[thread].rightCols(1) ); // put back G in correct column
  if(params->ncov_interaction > 1){

    // print b/se for each interaction term (from full model)
    for(int j = 0; j < params->ncov_interaction; j++){
      stmp="-INT_SNPx" + filters->interaction_cov + "=" + params->interaction_lvl_names[j];
      // print sum_stats
      if(params->htp_out) 
        buffer << print_sum_stats_head_htp(index, files->pheno_names[ipheno], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bsign * bhat(beg+1+j)/pheno_data->scf_i[thread](j), se(beg+1+j)/pheno_data->scf_i[thread](j), -1, -1, snp_data->af(ipheno), snp_data->info(ipheno), snp_data->mac(ipheno), snp_data->genocounts, ipheno, true, 1, params);
      else 
        buffer << (!params->split_by_pheno && (ipheno>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(ipheno) : snp_data->af1),snp_data->af_case(ipheno),snp_data->af_control(ipheno), (params->split_by_pheno ? snp_data->info(ipheno) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(ipheno) : snp_data->ns1),snp_data->ns_case(ipheno),snp_data->ns_control(ipheno), test_string + stmp, bsign * bhat(beg+1+j)/pheno_data->scf_i[thread](j), se(beg+1+j)/pheno_data->scf_i[thread](j), -1, -1, true, 1, params, (ipheno+1));
    }

    /// joint test for interaction
    bhat_s = bhat;
    bhat_s.tail(params->ncov_interaction) = 0;

    // run firth
    if(!fit_firth(ipheno, Y, pheno_data->Hmat[thread], offset, mask, pivec, etavec, bhat_s, se_s, beg+1, dev_s, false, tstat, params->maxstep, params->niter_max_firth, params->numtol_firth, params))
      return "";
    // pv
    tstat = dev_s - dev;
    if(tstat < 0) return "";
    logp = max(params->nl_dbl_dmin, cdf(complement(chisqI, tstat)));
    logp = -log10( logp );
    stmp="-INT_SNPx" + filters->interaction_cov; // single cov

    // print sum_stats
    if(params->htp_out) 
      buffer << print_sum_stats_head_htp(index, files->pheno_names[ipheno], model_type + stmp, snpinfo, params) << print_sum_stats_htp(-1, -1, tstat, logp, snp_data->af(ipheno), snp_data->info(ipheno), snp_data->mac(ipheno), snp_data->genocounts, ipheno, true, np-1, params);
    else 
      buffer << (!params->split_by_pheno && (ipheno>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(ipheno) : snp_data->af1),snp_data->af_case(ipheno),snp_data->af_control(ipheno), (params->split_by_pheno ? snp_data->info(ipheno) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(ipheno) : snp_data->ns1),snp_data->ns_case(ipheno),snp_data->ns_control(ipheno), test_string + stmp, -1, -1, tstat, logp, true, np-1, params, (ipheno+1));

  } else { // single interaction term

    bhat_s = bhat;
    bhat_s.tail(1) = 0;

    // run firth
    if(!fit_firth(ipheno, Y, pheno_data->Hmat[thread], offset, mask, pivec, etavec, bhat_s, se_s, beg+np-1, dev_s, false, tstat, params->maxstep, params->niter_max_firth, params->numtol_firth, params))
      return "";
    // T, beta, se & pv
    tstat = dev_s - dev;
    if(tstat < 0) return "";
    get_logp(logp, tstat);

    if(params->interaction_cat)
      stmp="-INT_SNPx" + filters->interaction_cov + "=" + params->interaction_lvl_names[0];
    else
      stmp="-INT_SNPx" + filters->interaction_cov; // single cov

    if( params->back_correct_se && (tstat > 0) )
      se_val = fabs(bhat(beg+1)) / sqrt(tstat);
    else
      se_val = se(beg+1);

    // print sum_stats
    if(params->htp_out) 
      buffer << print_sum_stats_head_htp(index, files->pheno_names[ipheno], model_type + stmp, snpinfo, params) << print_sum_stats_htp(bsign * bhat(beg+1)/pheno_data->scf_i[thread](0), se(beg+1)/pheno_data->scf_i[thread](0), tstat, logp, snp_data->af(ipheno), snp_data->info(ipheno), snp_data->mac(ipheno), snp_data->genocounts, ipheno, true, 1, params);
    else 
      buffer << (!params->split_by_pheno && (ipheno>0) ? "" : head) << print_sum_stats((params->split_by_pheno ? snp_data->af(ipheno) : snp_data->af1),snp_data->af_case(ipheno),snp_data->af_control(ipheno), (params->split_by_pheno ? snp_data->info(ipheno) : snp_data->info1), (params->split_by_pheno ? snp_data->ns(ipheno) : snp_data->ns1),snp_data->ns_case(ipheno),snp_data->ns_control(ipheno), test_string + stmp, bsign * bhat(beg+1)/pheno_data->scf_i[thread](0), se(beg+1)/pheno_data->scf_i[thread](0), tstat, logp, true, 1, params, (ipheno+1));

  }

  //if(isnp==0 & ipheno==0) cerr << endl << buffer.str() << endl; exit(-1);

  snp_data->test_fail_inter[ipheno] = false;
  return buffer.str() + buffer_joint.str();
}

