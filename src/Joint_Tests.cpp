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
#include "Pheno.hpp"
#include "NNLS.hpp"
#include "Joint_Tests.hpp"

using namespace std;
using namespace Eigen;
using boost::math::fisher_f;
using boost::math::chi_squared;
using boost::math::cauchy;

JTests::JTests() { // @suppress("Class members should be properly initialized")
}

JTests::~JTests() {
  // TODO Auto-generated destructor stub
}

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v){
  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}


bool JTests::get_test_info(const struct param* params, string const& test_string, mstream& sout){

  bool with_flip = true; // allow to flip to minor alleles
  std::vector< string > tmp_str_vec ;
  test_list = 0ULL;
  // for tests requiring QR decomp
  BIT_SET(qr_tests, joint_tests_map["ftest"]);
  BIT_SET(qr_tests, joint_tests_map["gates"]);
  BIT_SET(qr_tests, joint_tests_map["sbat"]);

  test_pfx = test_string + "-";
  burden_str = test_pfx + "BURDEN-";
  if(params->htp_out){
    if(params->skip_blups) burden_model = burden_str;
    else { test_pfx = test_string + "-WGR-"; burden_model = test_pfx + "BURDEN-"; }
  }

  // activate tests chosen by user
  tmp_str_vec = string_split(params->burden,",");
  for( auto const& input_test: tmp_str_vec ){

    if( input_test == "minp" || input_test == "gates" )
      BIT_SET(test_list, joint_tests_map[input_test]);
    else if( input_test == "ftest" ){
      if(params->trait_mode) sout << "WARNING: Joint F-test only for QTs.\n";
      else BIT_SET(test_list, joint_tests_map[input_test]);
    } else if( input_test == "sbat" ){
      if(params->trait_mode) sout << "WARNING: Joint SBAT test only for QTs.\n";
      else { 
        BIT_SET(test_list, joint_tests_map[input_test]); 
        nnls_napprox = params->nnls_napprox; 
        nnls_verbose_out = params->nnls_out_all;
      }
      with_flip = false;
    } else if( input_test == "acat" ){
      BIT_SET(test_list, joint_tests_map[input_test]);
      valid_snp_mode = !params->build_mask || (params->build_mask && params->mask_rule_max) ;
      acat_a1 = params->acat_a1;
      acat_a2 = params->acat_a2;
    } else throw "unrecognized joint test (='" + input_test + "').";

  }

  ncovars = params->ncov_analyzed;
  rng_rd = params->rng_rd;
  debug_mode = params->debug;
  nnls_adaptive = params->nnls_adaptive && CHECK_BIT(test_list,joint_tests_map["sbat"]);
  nnls_mt_weights = params->nnls_mt_weights && CHECK_BIT(test_list,joint_tests_map["sbat"]);
  apply_single_p = params->apply_gene_pval_strategy;

  if(apply_single_p) {
    check_class_genep(params->genep_mask_sets_file, params->mask_map);
    int nclass = genep_all_masks + gene_p_tests.size();
    if(nclass == 0) throw "No valid mask groups were specified for GENE_P strategy";
    else sout << " * number of mask groups run through gene-p strategy = " << nclass << "\n";
  }

  if(nnls_verbose_out) { // header
    string fname = out_file_prefix + "_sbat.info";
    ofstream file_info(fname, std::ios_base::out);
    file_info << "MASK_GROUP ID SEL BETA_SBAT BETA_OLS\n";
    file_info.close();
  }

  return with_flip;
}

vector<string> JTests::apply_joint_test(const int& chrom, const int& block, struct phenodt const* pheno_data, const Eigen::Ref<const Eigen::MatrixXd>& res, struct geno_block const* gblock, std::vector<variant_block>& block_info, vector<string> const& ynames, struct param const* params){

  int print_index, bs = setinfo[chrom - 1][block].snp_indices.size();
  vector<string> out_str(params->n_pheno);
  vector<vector<string>> sum_stats_str;
  sum_stats_str.resize(joint_tests_map.size());
  for (size_t i = 0; i < joint_tests_map.size(); i++)
    sum_stats_str[i].resize(params->n_pheno); // store sum stats for each test/phenotype
  std::map <std::string, std::map <std::string, bool>>::iterator itr;
  if( nnls_mt_weights ) prep_nnls_weights(gblock->Gmat.cols());

  for(int ph = 0; ph < params->n_pheno; ++ph) {

    std::map<std::string, double> overall_p;
    reset_vals();

    MapcMatXd yres (res.col(ph).data(), res.rows(), 1);
    string pheno_name = ynames[ph];

    // keep track of this when not splitting sum stats file
    bool run_tests = params->pheno_pass(ph) && (nvars > 0) && !set_vars(bs, ph, block_info);
    bool print_stats = run_tests || !params->split_by_pheno;

    if( CHECK_BIT(test_list,joint_tests_map["minp"]) ) { // minP
      if(run_tests) compute_minp();
      if(print_stats) sum_stats_str[joint_tests_map["minp"]][ph] = print_output(joint_tests_map["minp"], ph+1, chrom, block, pheno_name, params);
    } 
    if( CHECK_BIT(test_list,joint_tests_map["acat"]) ) { // ACAT
      if(run_tests && ((apply_single_p && genep_all_masks) || !apply_single_p)) {
        compute_acat(bs, ph, block_info);
        if(apply_single_p && genep_all_masks && (plog >= 0)) overall_p["BURDEN-ACAT"] = plog;
      }
      if(print_stats && ((apply_single_p && genep_all_masks) || !apply_single_p))
        sum_stats_str[joint_tests_map["acat"]][ph] = print_output(joint_tests_map["acat"], genep_all_sfx, ph+1, chrom, block, pheno_name, params);
    } 

    // check other test
    if( test_list & qr_tests ) {

      if(run_tests) compute_qr_G(pheno_data->masked_indivs.col(ph), gblock);

      if( CHECK_BIT(test_list,joint_tests_map["ftest"]) ) { // F-test
        if(run_tests) compute_ftest(pheno_data->masked_indivs.col(ph), yres); 
        if(print_stats) sum_stats_str[joint_tests_map["ftest"]][ph] = print_output(joint_tests_map["ftest"], ph+1, chrom, block, pheno_name, params);
      } 
      if( CHECK_BIT(test_list,joint_tests_map["gates"]) ) { // GATES
        if(run_tests) compute_gates(ph, block_info);
        if(print_stats) sum_stats_str[joint_tests_map["gates"]][ph] = print_output(joint_tests_map["gates"], ph+1, chrom, block, pheno_name, params);
      } 
      if( CHECK_BIT(test_list,joint_tests_map["sbat"]) ) { // SBAT (NNLS)
        if(run_tests && ((apply_single_p && genep_all_masks) || !apply_single_p))
          compute_nnls(pheno_data->masked_indivs.col(ph), yres, (genep_all_sfx == "" ? "ALL" : genep_all_sfx));
        else reset_vals();

        if( apply_single_p && genep_all_masks && valid_pval(pval_nnls_pos) && valid_pval(pval_nnls_neg)) overall_p["SBAT"] = plog;
        if( ((apply_single_p && genep_all_masks) || !apply_single_p) && print_stats) // default output
          sum_stats_str[joint_tests_map["sbat"]][ph] = print_output(joint_tests_map["sbat"], genep_all_sfx, ph+1, chrom, block, pheno_name, params);
        if((apply_single_p && genep_all_masks) || (!apply_single_p && nnls_verbose_out)) {
          // verbose output with NNLS pos & neg split into two
          // 1. NNLS pos
          if(run_tests && valid_pval(pval_nnls_pos)) get_pv(pval_nnls_pos);
          else reset_vals();
          if(print_stats) sum_stats_str[joint_tests_map["sbat_pos"]][ph] = print_output(joint_tests_map["sbat_pos"], genep_all_sfx, ph+1, chrom, block, pheno_name, params);
          // 2. NNLS neg
          if(run_tests && valid_pval(pval_nnls_neg)) get_pv(pval_nnls_neg);
          else reset_vals();
          if(print_stats) sum_stats_str[joint_tests_map["sbat_neg"]][ph] = print_output(joint_tests_map["sbat_neg"], genep_all_sfx, ph+1, chrom, block, pheno_name, params);
        }
      }
    }

    // should at least have burden-acat p-value
    if(apply_single_p) {
      if(run_tests) run_single_p_acat(bs, chrom, block, ph, pheno_name, block_info, overall_p, gblock, yres, pheno_data->masked_indivs.col(ph), sum_stats_str, params);
      else if(!params->split_by_pheno) { // when printing to single file and test failed
        if( genep_all_masks )
          sum_stats_str[joint_tests_map["gene_p"]][ph] = print_gene_output("GENE_P" + (genep_all_sfx == "" ? "" : "_" + genep_all_sfx), "", ph+1, chrom, block, pheno_name, params);
        for (itr = gene_p_tests.begin(); itr !=  gene_p_tests.end(); ++itr) 
          sum_stats_str[joint_tests_map["gene_p" + itr->first]][ph] = print_gene_output("GENE_P_" + itr->first, "", ph+1, chrom, block, pheno_name, params);
      }
    }

  }

  // store sum stats (if single file, store at index 0)
  for(size_t i = 0; i < joint_tests_map.size(); ++i)
    for(int ph = 0; ph < params->n_pheno; ++ph){
      print_index = params->split_by_pheno ? ph : 0;
      out_str[print_index].append( sum_stats_str[i][ph] );
    }

  return out_str;
}


// determine if marginal test failed 
bool JTests::set_vars(const int& bs, const int& ph, std::vector<variant_block> const& block_info){

  good_vars = ArrayXb::Constant(bs, false);
  log10pv = ArrayXd::Zero(bs);

  //if(debug_mode) cerr << "checking burden masks in set...";
  for(int isnp = 0; isnp < bs; isnp++){
    good_vars(isnp) = !block_info[isnp].ignored && !block_info[isnp].ignored_trait(ph) && !block_info[isnp].test_fail(ph);
    if(!good_vars(isnp)) continue;
    log10pv(isnp) = block_info[isnp].pval_log(ph);
  }
  nvars = good_vars.count();

  return (nvars == 0);
}


void JTests::compute_minp(){

  df_test = good_vars.count();
  if( df_test == 0 ) {reset_vals(); return;}

  // get minimum p-value (on log scale)
  get_pv( pow(10, -(log10pv.maxCoeff())) );

}


void JTests::compute_acat(const int& bs, const int& ph, const vector<variant_block>& block_info){

  double v_maf, tmpd;
  boost::math::beta_distribution<>  dist(acat_a1, acat_a2);

  df_test = good_vars.count();
  //if(debug_mode) cerr << "# burden masks for joint acat test = " << df_test << "\n";
  if( df_test == 0 ) {reset_vals(); return;}

  // make array of weights
  ArrayXd wts = ArrayXd::Zero(bs);
  for(int isnp = 0; isnp < bs; isnp++) {
    if( !good_vars(isnp) ) continue;
    // compute weights
    if( valid_snp_mode && !apply_single_p) {// sqrt(w)=dbeta(maf,a1,a2)*sqrt(maf*(1-maf))
      v_maf = min( block_info[isnp].af(ph), 1 - block_info[isnp].af(ph) );
      //cerr << v_maf << endl;
      tmpd = pdf( dist, v_maf );
      wts(isnp) = v_maf * (1-v_maf) * tmpd * tmpd;
    } else wts(isnp) = 1; // assume weight=1
  }
  //if(debug_mode) cerr << "done building acat weights\n";

  // get ACAT test stat
  get_chisq(get_acat(log10pv, wts));

}

double get_acat_robust(const Eigen::Ref<const ArrayXd>& logpvals, const Eigen::Ref<const ArrayXd>& weights){ // robust to low pvalues

  // if single pval, return log10p
  int n_pv = ((weights!=0) && (logpvals >= 0)).count();
  if(n_pv == 0) return -1;
  else if(n_pv == 1) return logpvals.maxCoeff();

  cauchy dc(0,1);
  double lpv_thr = 15, lpval_out;

  // split pvals by thr
  int n_A = ((weights!=0) && (logpvals >= lpv_thr)).count(); // very small pvals
  int n_B = ((weights!=0) && (logpvals >= 0) && (logpvals < lpv_thr)).count();
  double wsum = (logpvals >= 0).select(weights, 0).sum();
  double l_TA = 0, TB = 0;

  // T_A
  if(n_A > 0){ // compute on log scale to handle the very small pvalues
    ArrayXi vind = get_true_indices((weights!=0) && (logpvals >= lpv_thr));
    ArrayXd lp = logpvals( vind ), ws = weights( vind ) / wsum;
    ArrayXd zvec = lp * log(10) + ws.log() - log(M_PI);
    double zmax = zvec.maxCoeff();
    l_TA = zmax + log( (zvec - zmax).exp().sum() );
  }
  // T_B (can be negative)
  if(n_B > 0){
    ArrayXi vind = get_true_indices((weights!=0) && (logpvals >= 0) && (logpvals < lpv_thr));
    ArrayXd pv = pow(10, -logpvals(vind)).min(0.999); // avoid pvalues of 1
    ArrayXd ws = weights( vind ) / wsum; 
    TB = ( ws * tan( M_PI * (0.5 - pv)) ).sum(); 
  }

  // T_ACAT = TA + TB
  if(n_A == 0){ // avoid computing log(TB) as TB can be negative
    lpval_out = ( TB >= 8886111 ? -log(TB) - log(M_PI) : log(cdf(complement(dc, TB))) );
  } else if ((n_B == 0) || (TB == 0)){
    lpval_out = ( l_TA >= 16 ? -l_TA - log(M_PI) : log(cdf(complement(dc, exp(l_TA)))) );
  } else {
    double lsum; // get sum on log scale
    if(TB < 0){
      double l_abs_TB = log(fabs(TB));
      if(l_abs_TB < l_TA)
        lsum = l_TA + log1p(-exp(l_abs_TB - l_TA));
      else { // compute log(-Tacat)
        lsum = l_abs_TB + log1p(-exp(l_TA - l_abs_TB)); 
        lpval_out = ( lsum >= 16 ? log1p(-exp(-lsum-log(M_PI))) : log(cdf(complement(dc, -exp(lsum)))) );
        return -lpval_out/log(10);
      }
    } else {
      double l_TB = log(TB);
      lsum = fmax(l_TA, l_TB) + log1p(exp(-fabs(l_TB - l_TA)));
    } 
    lpval_out = ( lsum >= 16 ? -lsum - log(M_PI) : log(cdf(complement(dc, exp(lsum) ))) );
  }

  // return log10P
  return -lpval_out/log(10);
}

double get_acat_robust(const Eigen::Ref<const ArrayXd>& logpvals){
  ArrayXd wts = ArrayXd::Constant(logpvals.size(), 1); // uniform weights
  return get_acat_robust(logpvals, wts);
}

double get_acat(const Eigen::Ref<const ArrayXd>& logpvals, const Eigen::Ref<const ArrayXd>& weights){
  double logp = get_acat_robust(logpvals, weights);
  return logp;
}

double get_acat(const Eigen::Ref<const ArrayXd>& logpvals){ // uniform weights
  double logp = get_acat_robust(logpvals);
  return logp;
}

/*
double get_acat(const Eigen::Ref<const ArrayXd>& logpvals, const Eigen::Ref<const ArrayXd>& weights){

  cauchy dc(0,1);
  double tol = 10.0 * std::numeric_limits<double>::min(), pv_thr = 1e-15;

  // if single pval, return pval
  if(logpvals.size() == 1) {
    if((logpvals(0) >= 0) && (weights(0) != 0)) return pow(10, -logpvals(0));
    else return -1;
  }

  // use approx for small p-values (from ACAT R package)
  ArrayXd pvals = ((weights!=0) && (logpvals >= 0)).select( pow(10, -logpvals) , 0.5).max(tol).min(0.999); // to prevent underflow/overflow
  //cerr << "log10pv=" << logpvals.matrix().transpose() << "\npv=" << pvals.matrix().transpose() << "\nw=" << weights.matrix().transpose() << "\n";
  double acat = (pvals > pv_thr).select( weights * tan( M_PI * (0.5 - pvals)), (weights / pvals) / M_PI).sum();
  double wsum = (logpvals >= 0).select(weights, 0).sum();
  //cerr << std::setprecision(10) << "acat num=" << acat << " denum=" << wsum << endl;

  return cdf(complement(dc, acat/wsum ));
}

double get_acat(const Eigen::Ref<const ArrayXd>& logpvals){ // uniform weights

  cauchy dc(0,1);
  double tol = 10.0 * std::numeric_limits<double>::min(), pv_thr = 1e-15;

  // if single pval, return pval
  if(logpvals.size() == 1) return pow(10, -logpvals(0));

  // use approx for small p-values (from ACAT R package)
  ArrayXd pvals = pow(10, -logpvals).max(tol).min(0.999); // to prevent underflow/overflow
  double acat = (pvals > pv_thr).select( tan( M_PI * (0.5 - pvals)), (1.0 / pvals) / M_PI).sum();
  double wsum = logpvals.size();
  //cerr << std::setprecision(10) << "acat num=" << acat << " denum=" << wsum << endl;

  return cdf(complement(dc, acat/wsum ));
}
*/

void JTests::compute_qr_G(const Eigen::Ref<const MatrixXb>& mask, struct geno_block const* gblock){

  ArrayXi colkeep;
  MatrixXd Gnew;
  indices_vars.resize(0);

  // filter out bad variants
  Gnew = MatrixXd::Zero( gblock->Gmat.rows(), good_vars.count() );
  for(int i = 0, j = 0; i < gblock->Gmat.cols(); i++){
    if(!good_vars(i)) continue;
    Gnew.col(j++) = gblock->Gmat.col(i);
    indices_vars.push_back(i);
  }
  Gnew.array().colwise() *= mask.col(0).array().cast<double>();

  // find set of linearly independent cols
  ColPivHouseholderQR<MatrixXd> qrA(Gnew);
  qrA.setThreshold(qr_tol); 
  df_test = qrA.rank();

  if(df_test == 0) return;
  else if ( df_test < good_vars.count() ){
    colKeep = qrA.colsPermutation().indices();
    //ArrayXi tmp1(df_test);tmp1 << 0,1,3,4,5,6,7,9,10;colKeep = tmp1;
    //cerr << qr_tol << " -> " << colKeep.matrix().transpose().array() << endl;
    std::vector<int> new_indices;

    // keep only linearly independent columns
    Gtmp.resize(gblock->Gmat.rows(), df_test);

    for(int i = 0; i < df_test; i++){
      Gtmp.col(i) = Gnew.col( colKeep(i,0) );
      new_indices.push_back( indices_vars[ colKeep(i,0) ] );
    }

    indices_vars = new_indices;

  } else Gtmp = Gnew;

  /*
  // check min eigenvalue
  MatrixXd gtg = Gtmp.transpose() * Gtmp;
  SelfAdjointEigenSolver<MatrixXd> es(gtg, false);
  cerr << es.eigenvalues().head(2) << endl << endl;
  */

  //cerr << Gtmp.block(0,0,5,5) << endl;
}


void JTests::compute_ftest(const Eigen::Ref<const MatrixXb>& mask, const Eigen::Ref<const Eigen::MatrixXd>& ymat){

  if( df_test == 0 ) {
    reset_vals();
    return;
  }

  int ns = mask.col(0).array().count() - ncovars;
  int df_ur = ns - df_test;

  if( df_ur <= 0 ) {
    reset_vals();
    return;
  }

  double ss_r, ss_m, tmpd;
  ArrayXd y_tmp;
  MatrixXd bhat, GtG;
  fisher_f dist(df_test, df_ur);

  y_tmp = ymat.col(0).array() *  mask.col(0).array().cast<double>();
  GtG = Gtmp.transpose() * Gtmp;
  LLT<MatrixXd> lltOfA(GtG);
  bhat = lltOfA.solve(Gtmp.transpose() * y_tmp.matrix()) ; // vector
  //cerr << bhat << endl;

  ArrayXd yhat = (Gtmp * bhat).array();

  // SSM
  ss_m = yhat.square().sum();
  // SSR
  ss_r = ns - ss_m; // y is standardized

  // Ftest
  zval = (ss_m / df_test) / (ss_r / df_ur);
  //cerr << "DF1=" << df_test  << ";DF2=" << df_ur << ";SS1=" << ss_r << ";SS2=" << ss_m << ";F=" << zval << endl;

  // pvalue
  if(zval >= 0) {
    tmpd = cdf(complement(dist, zval)); 
    get_pv( tmpd );
  } else reset_vals();

}


void JTests::compute_nnls(const Eigen::Ref<const MatrixXb>& mask, const Eigen::Ref<const Eigen::MatrixXd>& ymat, string const& input_class){

  pval_nnls_pos = -1; pval_nnls_neg = -1;
  if( df_test == 0 ) {
    reset_vals();
    return;
  }

  int ns = mask.col(0).array().cast<int>().sum() - ncovars;
  int df_ur = ns - df_test, adapt_napprox = 2;

  if( df_ur <= 0 ) {
    reset_vals();
    return;
  }
  
  double pval_min2, adapt_thr = 1e-3; 
  VectorXd y_tmp = ymat.col(0).array() * mask.col(0).array().cast<double>();
  
  // (depreciated) compute NNLS p-value by function
  /* double pval_min2_old = jburden_test(y_tmp, Gtmp, df_ur, nnls_tol, nnls_napprox, true, 3); */

  // initialize an object of NNLS class & pass parameters
  NNLS nnls(nnls_napprox, nnls_normalize, nnls_tol, nnls_strict, nnls_verbose);
  nnls.gen = rng_rd;

  if(nnls_adaptive){ // run the NNLS test using adaptive strategy

    nnls.napprox = adapt_napprox;

    int p = Gtmp.cols();
    MatrixXd XtX(MatrixXd(p, p).setZero().selfadjointView<Lower>().rankUpdate(Gtmp.adjoint()));
    MatrixXd XtX_inv = XtX.llt().solve(MatrixXd::Identity(p, p));

    nnls.pw_weights(XtX_inv);
    if(nnls.nw > 0) nnls.pw_run(y_tmp, Gtmp, df_ur);

    if((nnls.pval_min2 >= 0) && (nnls.pval_min2 < adapt_thr)) {
      nnls.pw_weights(nnls_napprox);
      // note: no need to refit the NNLS model; done above by nnls.pw_run(Y, X)
      nnls.pw_calc_pvals();
    }

  } else if (nnls_mt_weights && nnls_weights[input_class][Gtmp.cols()-1].size() ) { // use pre-computed weights
    nnls.wts = nnls_weights[input_class][Gtmp.cols()-1];
    nnls.nw = nnls_weights[input_class][Gtmp.cols()-1].size();
    nnls.pw_run(y_tmp, Gtmp, df_ur);
  } else  // run the NNLS test: model fitting and inference
    nnls.run(y_tmp, Gtmp, df_ur);

  // get the final p-value = min(NNLS with b>=0, NNLS with b<=0)
  // -1 value means that NNLS failed; check the error message
  pval_min2 = nnls.pval_min2; 
  // get additional p-values & assign
  // to be used downstream if NNLS pos/neg are reported separately
  pval_nnls_pos = nnls.fit_pos.pval;
  pval_nnls_neg = nnls.fit_neg.pval;

  if(valid_pval(pval_min2)) {

    // store wts
    if (nnls_mt_weights && !nnls_weights[input_class][Gtmp.cols()-1].size() ) // store weights used
      nnls_weights[input_class][Gtmp.cols()-1] = nnls.wts;

    // pvalue
    //if(apply_single_p) { // compute pval_min2 using ACAT
    ArrayXd nnls_lpvs(2); 
    nnls_lpvs << -log10(max(nl_dbl_dmin, pval_nnls_pos)), -log10(max(nl_dbl_dmin, pval_nnls_neg)); // avoid 0 p-values (TBD: switch to log10p)
    get_chisq(get_acat(nnls_lpvs)); 
    pval_min2 = pval;
    //} else pval_min2 = min(1.0, 2 * pval_min2); // apply bonferroni correction

    // print extra NNLS information if requested
    if(nnls_verbose_out) {
      string fname = out_file_prefix + "_sbat.info";
      ofstream file_info(fname, std::ios_base::out | std::ios_base::app);

      // v1: tabular format
      // variant results per line
      for(int i = 0; i < df_test; i++) {
        unsigned int k = indices_vars[i]; 
        file_info << input_class << " "; // to accomodate for gene_p classes
        file_info << variant_names[k] << " "; // variant ID
        file_info << nnls.str_sel_i(i) << " "; // NNLS selection status
        file_info << nnls.str_bhat_i(i) << " "; // NNLS beta
        file_info << nnls.bhat_ols[i] << "\n"; // OLS beta
      }


      // v2: non-tabular format
      /* // write line with variant names */
      /* for(unsigned int i = 0; i < df_test; i++) { */
      /*   file_info << variant_names[ indices_vars[i] ] << " "; */
      /* } */
      /* file_info << endl; */

      /* // write nnls info */ 
      /* file_info << nnls.str_info(); */

      file_info.close();
    } 

  } else reset_vals();

}

void JTests::compute_gates(const int& ph, const std::vector<variant_block>& block_info){

  int gcol;
  double p_gates, m_e, p_i, m_ei;

  if( df_test == 0 ) {
    reset_vals();
    return;
  } else if( df_test == 1 ){
    p_gates = pow(10, -log10pv(indices_vars[0]));
    if( p_gates >= 0) get_pv( p_gates );
    else reset_vals();
    return;
  }

  // Sort p-values
  vector<double> pvals, sorted_pv;
  MatrixXd tmpG = Gtmp;

  pvals.resize(df_test);
  for(int i = 0; i < df_test; i++) {
    //cerr << log10pv(indices_vars[i]) ;
    pvals[i] = pow(10, -log10pv(indices_vars[i]));
    //cerr << "\t" << pvals[i] << endl;
  }

  gcol = 0;
  for (auto i: sort_indexes(pvals)) {
    sorted_pv.push_back(pvals[i]);
    tmpG.col(gcol++) = Gtmp.col(i);
    //cout << i << " " << pvals[i] << endl;
  }

  // Compute corr(G)
  MatrixXd GtG, corP;
  GtG = tmpG.transpose() * tmpG / scale_denum;
  //cerr << endl << scale_denum << endl << GtG.block(0,0,4,4) << endl << endl ;

  corP.resize( GtG.rows(), GtG.cols() ); 
  //approximation used from Postgwas R package in gene2p.R (Milan Hiersche et al, 2013)
  corP.array() = 0.7723 * GtG.array().pow(6) -
    1.5659 * GtG.array().pow(5) +
    1.201 * GtG.array().pow(4) -
    0.2355 * GtG.array().pow(3) +
    0.2184 * GtG.array().pow(2) +
    0.6086 * GtG.array();
  //cerr << corP.block(0,0,4,4) << endl;
  //cerr<<endl<<sorted_pv[0] << " " << sorted_pv[1] << endl << endl;

  // Me
  m_e = get_me(corP);
  //cerr << m_e << "\n\n";

  // get Gates p-value
  p_gates = 1;
  for(int i = 0; i < df_test; i++) {
    m_ei = get_me( corP.block(0,0,i+1,i+1) );
    //if(i<2) cerr << m_ei << endl;
    p_i = m_e * sorted_pv[i] / m_ei;
    if(p_i < p_gates) p_gates = p_i;
  }

  // pvalue
  if( p_gates >= 0) get_pv( p_gates );
  else reset_vals();
  
}

double JTests::get_me(const Ref<const MatrixXd>& ldmat){

  int ncols = ldmat.cols();
  double m_e; 

  if(ncols == 1) return 1;

  // Get eigen values
  SelfAdjointEigenSolver<MatrixXd> es(ldmat, Eigen::EigenvaluesOnly);
  ArrayXd D = es.eigenvalues().array();
  //cerr << D.head(5) << endl;
  m_e = ncols - ( D > 1 ).select(D - 1, 0).sum();

  return m_e;
}

void JTests::run_single_p_acat(int const& bs, const int& chrom, const int& block, int const& ph, const string& pheno_name, std::vector<variant_block>& block_info, std::map<std::string, double>& overall_p, struct geno_block const* gblock, const Eigen::Ref<const Eigen::MatrixXd>& yres, const Eigen::Ref<const MatrixXb>& mask, vector<vector<string>>& sum_stats_str, struct param const* params){

  double max_logp = -1, pv; 
  string mname, max_logp_mask = "";
  vector<string> keep_tests = { "ACATV", "SKATO-ACAT" };
  ArrayXd pvals_gene;

  // gene_p combining all masks
  if( genep_all_masks ){
    vector<double> acatv_acat, skato_acat;
    // get SKATO/ACATV p-values as well as top mask
    for(int imask = 0; imask < bs; imask++){
      mname = block_info[imask].mask_name;
      if( (log10pv(imask) > max_logp) && (log10pv(imask) > 0) ){ // check strongest signal also in burden-only test
        max_logp_mask = mname;
        max_logp = log10pv(imask);
      }
      if(block_info[imask].skip_for_vc) continue; 
      for(auto const& extract_test : keep_tests)
        if(in_map(extract_test, block_info[imask].sum_stats_vc)){
          pv = block_info[imask].sum_stats_vc[extract_test](ph,1); 
          if(pv>=0){
            if(pv>max_logp){
              max_logp_mask = mname;
              max_logp = pv;
            }
            if(extract_test == "ACATV") acatv_acat.push_back( pv );
            else if(extract_test == "SKATO-ACAT") skato_acat.push_back( pv );
          }
        }
    }
    // compute acat for acatv & skato
    if(acatv_acat.size() > 0){
      df_test = acatv_acat.size();
      ArrayXd pvals_arr = MapArXd( acatv_acat.data(), df_test); 
      get_chisq(get_acat(pvals_arr));
      if(plog>=0) overall_p["ACATV-ACAT"] = plog;
      sum_stats_str[joint_tests_map["acatv_acat"]][ph] = print_gene_output(test_pfx + "ACATV-ACAT" + (genep_all_sfx == "" ? "" : "_" + genep_all_sfx), "", ph+1, chrom, block, pheno_name, params);
    } else if(!params->split_by_pheno){
      reset_vals();
      sum_stats_str[joint_tests_map["acatv_acat"]][ph] = print_gene_output(test_pfx + "ACATV-ACAT" + (genep_all_sfx == "" ? "" : "_" + genep_all_sfx), "", ph+1, chrom, block, pheno_name, params);
    }
    if(skato_acat.size() > 0){
      df_test = skato_acat.size();
      ArrayXd pvals_arr = MapArXd( skato_acat.data(), df_test); 
      get_chisq(get_acat(pvals_arr));
      if(plog>=0) overall_p["SKATO-ACAT"] = plog;
      sum_stats_str[joint_tests_map["skato_acat"]][ph] = print_gene_output(test_pfx + "SKATO-ACAT" + (genep_all_sfx == "" ? "" : "_" + genep_all_sfx), "", ph+1, chrom, block, pheno_name, params);
    } else if(!params->split_by_pheno){
      reset_vals();
      sum_stats_str[joint_tests_map["skato_acat"]][ph] = print_gene_output(test_pfx + "SKATO-ACAT" + (genep_all_sfx == "" ? "" : "_" + genep_all_sfx), "", ph+1, chrom, block, pheno_name, params);
    }
    // combine all p-values and pass through acat
    if(overall_p.size()>0){
      map_to_vec(df_test, overall_p, pvals_gene);
      get_chisq(get_acat(pvals_gene));
      sum_stats_str[joint_tests_map["gene_p"]][ph] = print_gene_output("GENE_P" + (genep_all_sfx == "" ? "" : "_" + genep_all_sfx), max_logp_mask, ph+1, chrom, block, pheno_name, params);
    } else if(!params->split_by_pheno){
      reset_vals();
      sum_stats_str[joint_tests_map["gene_p"]][ph] = print_gene_output("GENE_P" + (genep_all_sfx == "" ? "" : "_" + genep_all_sfx), "", ph+1, chrom, block, pheno_name, params);
    }
  }

  // go through each set of masks
  std::map <std::string, std::map <std::string, bool>>::iterator itr;
  for (itr = gene_p_tests.begin(); itr !=  gene_p_tests.end(); ++itr) {

    std::map <std::string, double> overall_p_set;
    good_vars = false;
    max_logp = -1;
    max_logp_mask = "";
    bool get_top_mask = itr->second.size() > 1;
    vector<double> acatv_acat, skato_acat;
    if(params->debug) cerr << itr->first << ":\n";

    // identify all the masks in the set
    for(int imask = 0; imask < bs; imask++){
      mname = block_info[imask].mask_name;
      good_vars(imask) = in_map(mname, itr->second) && !block_info[imask].ignored && !block_info[imask].ignored_trait(ph) && !block_info[imask].test_fail(ph);
      //if(params->debug) cerr << mname << " " << std::boolalpha << in_map(mname, itr->second) << " && " << (!block_info[imask].ignored && !block_info[imask].ignored_trait(ph) && !block_info[imask].test_fail(ph)) << " -> " << good_vars(imask) << "\n";
      if(!good_vars(imask)) continue;
      //cerr << mname << ":" << log10pv(imask) << "\n";
      if( get_top_mask && (log10pv(imask) > max_logp) && (log10pv(imask) > 0) ){ // check strongest signal also in burden-only test
        max_logp_mask = mname;
        max_logp = log10pv(imask);
      }
      if(block_info[imask].skip_for_vc) continue; 
      for(auto const& extract_test : keep_tests)
        if( in_map(extract_test, block_info[imask].sum_stats_vc) ){
          pv = block_info[imask].sum_stats_vc[extract_test](ph,1); 
          if(pv>=0){
            if(get_top_mask && (pv > max_logp)){
              max_logp_mask = mname;
              max_logp = pv;
            }
            if(extract_test == "ACATV") acatv_acat.push_back( pv );
            else if(extract_test == "SKATO-ACAT") skato_acat.push_back( pv );
          }
        }
    }

    if(good_vars.any()){
      if(params->debug) cerr << "M=" << good_vars.count() << "\n";

      // run acat
      compute_acat(bs, ph, block_info);
      if(plog >= 0) overall_p_set["BURDEN-ACAT"] = plog;
      sum_stats_str[joint_tests_map["acat" + itr->first]][ph].append(print_output(joint_tests_map["acat"], itr->first, ph+1, chrom, block, pheno_name, params));

      // run nnls
      if( CHECK_BIT(test_list, joint_tests_map["sbat"]) ) {
        compute_qr_G(mask, gblock);
        compute_nnls(mask, yres, itr->first);
        if(valid_pval(pval_nnls_pos) && valid_pval(pval_nnls_neg)) {
          overall_p_set["SBAT"] = plog;
          sum_stats_str[joint_tests_map["sbat" + itr->first]][ph] = print_output(joint_tests_map["sbat"], itr->first, ph+1, chrom, block, pheno_name, params);
          get_pv(pval_nnls_pos);
          sum_stats_str[joint_tests_map["sbat_pos" + itr->first]][ph] = print_output(joint_tests_map["sbat_pos"], itr->first, ph+1, chrom, block, pheno_name, params);
          get_pv(pval_nnls_neg);
          sum_stats_str[joint_tests_map["sbat_neg" + itr->first]][ph] = print_output(joint_tests_map["sbat_neg"], itr->first, ph+1, chrom, block, pheno_name, params);
        } else if(!params->split_by_pheno) {
          reset_vals();
          sum_stats_str[joint_tests_map["sbat" + itr->first]][ph] = print_output(joint_tests_map["sbat"], itr->first, ph+1, chrom, block, pheno_name, params);
          sum_stats_str[joint_tests_map["sbat_pos" + itr->first]][ph] = print_output(joint_tests_map["sbat_pos"], itr->first, ph+1, chrom, block, pheno_name, params);
          sum_stats_str[joint_tests_map["sbat_neg" + itr->first]][ph] = print_output(joint_tests_map["sbat_neg"], itr->first, ph+1, chrom, block, pheno_name, params);
        }
      }

      // compute acat for acatv & skato
      if(acatv_acat.size() > 0){
        df_test = acatv_acat.size();
        ArrayXd pvals_arr = MapArXd( acatv_acat.data(), df_test); 
        get_chisq(get_acat(pvals_arr));
        if(plog>=0) overall_p_set["ACATV-ACAT"] = plog;
        sum_stats_str[joint_tests_map["acatv_acat" + itr->first]][ph] = print_gene_output(test_pfx + "ACATV-ACAT_" + itr->first, "", ph+1, chrom, block, pheno_name, params);
      } else if(!params->split_by_pheno){
        reset_vals();
        sum_stats_str[joint_tests_map["acatv_acat" + itr->first]][ph] = print_gene_output(test_pfx + "ACATV-ACAT_" + itr->first, "", ph+1, chrom, block, pheno_name, params);
      }
      if(skato_acat.size() > 0){
        df_test = skato_acat.size();
        ArrayXd pvals_arr = MapArXd( skato_acat.data(), df_test); 
        get_chisq(get_acat(pvals_arr));
        if(plog>=0) overall_p_set["SKATO-ACAT"] = plog;
        sum_stats_str[joint_tests_map["skato_acat" + itr->first]][ph] = print_gene_output(test_pfx + "SKATO-ACAT_" + itr->first, "", ph+1, chrom, block, pheno_name, params);
      } else if(!params->split_by_pheno){
        reset_vals();
        sum_stats_str[joint_tests_map["skato_acat" + itr->first]][ph] = print_gene_output(test_pfx + "SKATO-ACAT_" + itr->first, "", ph+1, chrom, block, pheno_name, params);
      }

      // apply acat to all p
      if(overall_p_set.size()>0){
        map_to_vec(df_test, overall_p_set, pvals_gene);
        get_chisq(get_acat(pvals_gene));
        sum_stats_str[joint_tests_map["gene_p" + itr->first]][ph] = print_gene_output("GENE_P_" + itr->first, max_logp_mask, ph+1, chrom, block, pheno_name, params);
      } else if(!params->split_by_pheno){
        reset_vals();
        sum_stats_str[joint_tests_map["gene_p" + itr->first]][ph] = print_gene_output("GENE_P_" + itr->first, "", ph+1, chrom, block, pheno_name, params);
      }

    } else if(!params->split_by_pheno){
      reset_vals();
      // print NA for all tests for that phenotype
      sum_stats_str[joint_tests_map["acat" + itr->first]][ph].append(print_output(joint_tests_map["acat"], itr->first, ph+1, chrom, block, pheno_name, params));
      if( CHECK_BIT(test_list, joint_tests_map["sbat"]) ) {
          sum_stats_str[joint_tests_map["sbat" + itr->first]][ph] = print_output(joint_tests_map["sbat"], itr->first, ph+1, chrom, block, pheno_name, params);
          sum_stats_str[joint_tests_map["sbat_pos" + itr->first]][ph] = print_output(joint_tests_map["sbat_pos"], itr->first, ph+1, chrom, block, pheno_name, params);
          sum_stats_str[joint_tests_map["sbat_neg" + itr->first]][ph] = print_output(joint_tests_map["sbat_neg"], itr->first, ph+1, chrom, block, pheno_name, params);
      }
      sum_stats_str[joint_tests_map["acatv_acat" + itr->first]][ph] = print_gene_output(test_pfx + "ACATV-ACAT_" + itr->first, "", ph+1, chrom, block, pheno_name, params);
      sum_stats_str[joint_tests_map["skato_acat" + itr->first]][ph] = print_gene_output(test_pfx + "SKATO-ACAT_" + itr->first, "", ph+1, chrom, block, pheno_name, params);
      sum_stats_str[joint_tests_map["gene_p" + itr->first]][ph] = print_gene_output("GENE_P_" + itr->first, "", ph+1, chrom, block, pheno_name, params);
    }
  }

}

string JTests::print_output(const int& ttype, const int& ipheno, const int& chrom, const int& block, const string& pheno_name, struct param const* params){

  if(!params->htp_out) return print_sum_stats(test_names[ttype], ipheno, chrom, block, params);
  else return print_sum_stats_htp(test_names[ttype], chrom, block, pheno_name, ipheno, params);

}

string JTests::print_output(const int& ttype, string const& tsuf, const int& ipheno, const int& chrom, const int& block, const string& pheno_name, struct param const* params){

  if(!params->htp_out) return print_sum_stats(test_names[ttype] + (tsuf == "" ? "" : "_" + tsuf), ipheno, chrom, block, params);
  else return print_sum_stats_htp(test_names[ttype] + (tsuf == "" ? "" : "_" + tsuf), chrom, block, pheno_name, ipheno, params);

}

// normal regenie format
std::string JTests::print_sum_stats(const string& tname, const int& ipheno, const int& chrom, const int& block, struct param const* params){

  std::ostringstream buffer;

  // chr pos id a0 a1 af
  if(params->split_by_pheno || ipheno == 1) {
    buffer << setinfo[chrom - 1][block].chrom << " " << setinfo[chrom - 1][block].physpos << " " << setinfo[chrom - 1][block].ID << " NA NA NA " ;
    if( params->af_cc ) buffer << "NA NA ";
    // info
    if(!params->build_mask && params->dosage_mode) buffer << "NA ";
    // n
    if(params->split_by_pheno) buffer << params->pheno_counts.row(ipheno-1).sum() << " ";
    else buffer << "NA ";
    if( params->af_cc ) {
      if(params->split_by_pheno) buffer << params->pheno_counts(ipheno-1, 0) << " " << params->pheno_counts(ipheno-1, 1) << " ";
      else buffer << "NA NA ";
    }
    // test
    buffer << burden_str << tname << " ";
  }

  //beta se
  buffer << "NA NA ";

  // chisq
  if( zval != -9 ) buffer << zval << " ";
  else buffer << "NA ";

  // pval
  if( plog != -9 )  buffer << plog << " ";
  else buffer << "NA ";

  //df (print it out only if split by pheno)
  if(params->split_by_pheno || (ipheno == params->n_pheno)) {
    if(params->split_by_pheno && (plog != -9))  buffer << "DF=" << df_test << endl;
    else buffer << "DF=NA\n";
  }

  reset_vals();
  return buffer.str();
}


// htpv4 format
std::string JTests::print_sum_stats_htp(const string& tname, const int& chrom, const int& block, const string& yname, const int& ipheno, struct param const* params){

  std::ostringstream buffer;
  bool test_pass = (pval != -9);
  const string cohort = params->cohort_name;

  if( pval == 1) pval = 1 - 1e-7;

  // SNP info
  buffer << setinfo[chrom - 1][block].ID << "\t"<< setinfo[chrom - 1][block].chrom << "\t" << setinfo[chrom - 1][block].physpos << "\tref\tset\t";
  // trait, cohort, test name
  buffer << yname << "\t" << params->cohort_name << "\t" << burden_model << tname;

  // bhat & 95 CI
  buffer << "\tNA\tNA\tNA\t" ;
  // Pvalue
  if(test_pass) buffer << pval << "\t";
  else buffer << "NA\t";

  // print out AF, counts in cases, counts in controls
  buffer << "NA\t" << params->pheno_counts(ipheno-1, 0) << "\tNA\tNA\tNA\t";
  if(params->trait_mode == 1) buffer << params->pheno_counts(ipheno-1, 1);
  else buffer << "NA";
  buffer << "\tNA\tNA\tNA\t";

  // info column
  if(test_pass) buffer << "DF=" << df_test;
  else buffer << "DF=0";
  // log10P
  if(test_pass) buffer << ";LOG10P=" << plog;
  else buffer << ";LOG10P=NA";

  buffer << ";NO_BETA\n";

  reset_vals();
  return buffer.str();

}

// single gene p
string JTests::print_gene_output(const string& mname, const string& max_name, const int& ipheno, const int& chrom, const int& block, const string& pheno_name, struct param const* params){

  if(!params->htp_out) return print_sum_stats_gene(mname, max_name, ipheno, chrom, block, params);
  else return print_sum_stats_htp_gene(mname, max_name, chrom, block, pheno_name, ipheno, params);

}

// normal regenie format
std::string JTests::print_sum_stats_gene(const string& mname, const string& max_name, const int& ipheno, const int& chrom, const int& block, struct param const* params){

  std::ostringstream buffer;

  if(params->split_by_pheno || ipheno == 1) {
    // chr pos id a0 a1 af
    buffer << setinfo[chrom - 1][block].chrom << " " << setinfo[chrom - 1][block].physpos << " " << setinfo[chrom - 1][block].ID << " NA NA NA " ;
    if( params->af_cc ) buffer << "NA NA ";
    // info
    if(!params->build_mask && params->dosage_mode) buffer << "NA ";
    // n
    if(params->split_by_pheno) buffer << params->pheno_counts.row(ipheno-1).sum() << " ";
    else buffer << "NA ";
    if( params->af_cc ) {
      if(params->split_by_pheno) buffer << params->pheno_counts(ipheno-1, 0) << " " << params->pheno_counts(ipheno-1, 1) << " ";
      else buffer << "NA NA ";
    }
    // test
    buffer << mname << " ";
  }

  //beta se
  buffer << "NA NA ";

  // chisq
  if( zval != -9 ) buffer << zval << " ";
  else buffer << "NA ";

  // pval
  if( plog != -9 )  buffer << plog << " ";
  else buffer << "NA ";

  //df (print it out only if split by pheno)
  if(params->split_by_pheno || (ipheno == params->n_pheno)) {
    if(params->split_by_pheno && (plog != -9))  buffer << "DF=" << df_test;
    else buffer << "DF=NA";

    // top signal (only if split files)
    if( params->split_by_pheno && (max_name != "") )
      buffer << ";STRONGEST_MASK=" << max_name;

    buffer << endl;
  }

  reset_vals();
  return buffer.str();
}


// htpv4 format
std::string JTests::print_sum_stats_htp_gene(const string& mname, const string& max_name, const int& chrom, const int& block, const string& yname, const int& ipheno, struct param const* params){

  std::ostringstream buffer;
  bool test_pass = (pval != -9);
  const string cohort = params->cohort_name;

  if( pval == 1) pval = 1 - 1e-7;

  // SNP info
  buffer << setinfo[chrom - 1][block].ID << "\t"<< setinfo[chrom - 1][block].chrom << "\t" << setinfo[chrom - 1][block].physpos << "\tref\tset\t";
  // trait, cohort, test name
  buffer << yname << "\t" << params->cohort_name << "\t" << mname;

  // bhat & 95 CI
  buffer << "\tNA\tNA\tNA\t" ;
  // Pvalue
  if(test_pass) buffer << pval << "\t";
  else buffer << "NA\t";

  // print out AF, counts in cases, counts in controls
  buffer << "NA\t" << params->pheno_counts(ipheno-1, 0) << "\tNA\tNA\tNA\t";
  if(params->trait_mode == 1) buffer << params->pheno_counts(ipheno-1, 1);
  else buffer << "NA";
  buffer << "\tNA\tNA\tNA\t";

  // info column
  if(test_pass) buffer << "DF=" << df_test;
  else buffer << "DF=0";
  // top signal
  if(max_name != "") buffer << ";STRONGEST_MASK=" << max_name;
  // log10P
  if(test_pass) buffer << ";LOG10P=" << plog;
  else buffer << ";LOG10P=NA";
  buffer << ";NO_BETA\n";

  reset_vals();
  return buffer.str();

}

void JTests::get_variant_names(int const& chrom, int const& block, vector<snp> const& snpinfo){

  // only for NNLS (for now)
  if( !(CHECK_BIT(test_list, joint_tests_map["sbat"]) && nnls_verbose_out) ) return;

  vector<uint64> *indices =  &(setinfo[chrom - 1][block].snp_indices);
  int bs = indices->size();
  variant_names.resize(bs);

  for(int i = 0; i < bs; i++)
    variant_names[i] = snpinfo[ indices->at(i) ].ID;

}

void JTests::check_class_genep(string const& mask_set_file, std::map<std::string, bool> const& mask_map){

  if(mask_set_file == ""){
    vector<string> tmp_vec = {"M1"};//plof only
    add_class("M1", tmp_vec, mask_map);
    tmp_vec[0] = "pLoF";//plof only
    add_class("pLoF", tmp_vec, mask_map);
    tmp_vec[0] = "LoF";//plof only
    add_class("LoF", tmp_vec, mask_map);

  } else {

    genep_all_masks = false;
    string line;
    std::vector< string > tmp_str_vec;
    ifstream mfile;
    mfile.open( mask_set_file, ios::in );

    while( getline(mfile, line) ){
      tmp_str_vec = string_split(line," \t");
      if(tmp_str_vec.size() < 2)
        throw "invalid line = '" + line + "'";
      add_class(tmp_str_vec[0], string_split(tmp_str_vec[1],","), mask_map);
    }
    mfile.close();
    
  }

  // allocate entry in map for each test on set of masks
  // make sure results are printed in the order of each test then list of mask groups
  std::map <std::string, std::map <std::string, bool>>::iterator itr;
  vector<string> tests = {"acat", "sbat", "acatv_acat", "skato_acat", "gene_p"};
  for (auto const& test_name : tests)
    for (itr = gene_p_tests.begin(); itr !=  gene_p_tests.end(); ++itr)
      if( (test_name == "sbat") && CHECK_BIT(test_list, joint_tests_map["sbat"]) ) {
          joint_tests_map[test_name + itr->first] = joint_tests_map.size();
          joint_tests_map[test_name + "_pos" + itr->first] = joint_tests_map.size();
          joint_tests_map[test_name + "_neg" + itr->first] = joint_tests_map.size();
      } else if (test_name != "sbat") joint_tests_map[test_name + itr->first] = joint_tests_map.size();

}

void JTests::add_class(string const& sfx_test, vector<string> const& mask_vec, std::map<std::string, bool> const& mask_map){
  
  std::map <std::string, bool> tmp_map;// keep track of masks in genep class

  // check it is a valid mask
  for(auto const& mask : mask_vec) 
    if(in_map(mask, mask_map))
      tmp_map[mask] = true;

  if(in_map(sfx_test, gene_p_tests))
    throw "GENE_P_'" + sfx_test + "' has already been defined (check for duplicates in the `--rgc-gene-def` file).";

  // check it has at least one mask
  if( tmp_map.size() > 0 ){
    if(tmp_map.size() == mask_map.size()) {
      genep_all_masks = true;
      genep_all_sfx = sfx_test;
    } else gene_p_tests[sfx_test] = tmp_map;
  }

}

void JTests::prep_nnls_weights(int const& max_cols){
  
  std::vector<VectorXd> tmp_v(max_cols);
  std::map <std::string, std::map <std::string, bool>>::iterator itr;

  // for each gene-p class, fill with empty weight vectors
  for (itr = gene_p_tests.begin(); itr !=  gene_p_tests.end(); ++itr)  
    nnls_weights[itr->first] = tmp_v;
  nnls_weights[ (genep_all_sfx == "" ? "ALL" : genep_all_sfx) ] = tmp_v;

}

void map_to_vec(int& nvals, std::map <std::string, double>& map_pvals, ArrayXd& vec_pvals){
  std::map <std::string, double>::iterator itr;
  nvals = map_pvals.size();
  vec_pvals.resize( nvals );
  int i = 0;
  for (itr = map_pvals.begin(); itr !=  map_pvals.end(); ++itr) 
    vec_pvals(i++) = itr->second;
}

void JTests::reset_vals(){
  pval = -9, plog = -9, zval = -9;
}

void JTests::get_pv(const double& pv){

  if(pv < 0) {reset_vals(); return;}
  chi_squared chisq(1);

  pval = max(nl_dbl_dmin, pv); // to prevent underflow
  zval = quantile(complement(chisq, pval)); // chisq stat
  plog = -log10(pval); // -log10p

}

void JTests::get_chisq(const double& lpv){

  if(lpv < 0) {reset_vals(); return;}

  get_chisq_stat_pv(pval, zval, lpv, nl_dbl_dmin, log10_nl_dbl_dmin);
  plog = lpv;

}


bool valid_pval(double const& pv){
  return (pv >= 0) && (pv <= 1);
}
