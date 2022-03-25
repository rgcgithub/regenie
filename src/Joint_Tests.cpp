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
  test_list = 0u;

  burden_str = test_string + "-BURDEN-";
  if(params->htp_out){
    if(params->skip_blups) burden_model = burden_str;
    else burden_model = test_string + "-WGR-BURDEN-";
  }

  // activate tests chosen by user
  tmp_str_vec = string_split(params->burden,",");
  for( auto const& input_test: tmp_str_vec ){

    if( input_test == "minp" || input_test == "gates" )
      BIT_SET(test_list, joint_tests_map[input_test]);
    else if( input_test == "ftest" ){
      if(params->trait_mode) sout << "WARNING: Joint F-test only for QTs.\n";
      else BIT_SET(test_list, joint_tests_map[input_test]);
    } else if( input_test == "nnls" ){
      if(params->trait_mode) sout << "WARNING: Joint NNLS test only for QTs.\n";
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

  ncovars = params->ncov;
  apply_single_p = params->apply_gene_pval_strategy;

  return with_flip;
}

string JTests::apply_joint_test(const int& chrom, const int& block, const int& ph, struct phenodt const* pheno_data, const Eigen::Ref<const Eigen::MatrixXd>& yres, struct geno_block const* gblock, std::vector<variant_block>& block_info, const string& pheno_name, struct param const* params){

  int bs = setinfo[chrom - 1][block].snp_indices.size();
  std::map<std::string, double> overall_p;
  std::ostringstream buffer;
  reset_vals();

  if(nvars == 0 || set_vars(bs, ph, block_info)) 
    return ""; // if no variants in set passed filters
  
  if( CHECK_BIT(test_list,joint_tests_map["minp"]) ) { // minP
    compute_minp();
    buffer << print_output(joint_tests_map["minp"], ph+1, chrom, block, pheno_name, params);
  } 
  if( CHECK_BIT(test_list,joint_tests_map["acat"]) ) { // ACAT
    compute_acat(bs, ph, block_info);
    if(apply_single_p && (plog >= 0)) overall_p["BURDEN-ACAT"] = plog;
    buffer << print_output(joint_tests_map["acat"], ph+1, chrom, block, pheno_name, params);
  } 

  // check other test
  if( test_list & qr_tests ) {

    compute_qr_G(pheno_data->masked_indivs.col(ph), gblock);

    if( CHECK_BIT(test_list,joint_tests_map["ftest"]) ) { // F-test
      compute_ftest(pheno_data->masked_indivs.col(ph), yres); 
      buffer << print_output(joint_tests_map["ftest"], ph+1, chrom, block, pheno_name, params);
    } 
    if( CHECK_BIT(test_list,joint_tests_map["gates"]) ) { // GATES
      compute_gates(ph, block_info);
      buffer << print_output(joint_tests_map["gates"], ph+1, chrom, block, pheno_name, params);
    } 
    if( CHECK_BIT(test_list,joint_tests_map["nnls"]) ) { // NNLS
      compute_nnls(pheno_data->masked_indivs.col(ph), yres); 

      if(!nnls_verbose_out) { 
        // default output
        buffer << print_output(joint_tests_map["nnls"], ph+1, chrom, block, pheno_name, params);
      } else {
        // verbose output with NNLS pos & neg split into two 
        // 1. NNLS pos (test code 5)
        if((pval_nnls_pos >= 0) && (pval_nnls_pos <= 1)) get_pv(pval_nnls_pos);
        else reset_vals();
        buffer << print_output(joint_tests_map["nnls_pos"], ph+1, chrom, block, pheno_name, params);
        // 2. NNLS neg (test code 6)
        if((pval_nnls_neg >= 0) && (pval_nnls_neg <= 1)) get_pv(pval_nnls_neg);
        else reset_vals();
        buffer << print_output(joint_tests_map["nnls_neg"], ph+1, chrom, block, pheno_name, params);
      }

      if( apply_single_p &&
          (pval_nnls_pos >= 0) && (pval_nnls_pos <= 1) &&
          (pval_nnls_neg >= 0) && (pval_nnls_neg <= 1) 
        ) {
        get_pv(pval_nnls_pos);overall_p["NNLS_POS"] = plog;
        get_pv(pval_nnls_neg);overall_p["NNLS_NEG"] = plog;
      }

    }
  }

  // should at least have burden-acat p-value
  if(apply_single_p) buffer << run_single_p_acat(bs, chrom, block, ph, pheno_name, block_info, overall_p, gblock, yres, pheno_data->masked_indivs.col(ph), params);

  return buffer.str();
}


// determine if marginal test failed 
bool JTests::set_vars(const int& bs, const int& ph, std::vector<variant_block> const& block_info){

  log10pv.resize(bs);
  good_vars.resize(bs);

  for(int isnp = 0; isnp < bs; isnp++){
    log10pv(isnp) = block_info[isnp].pval_log(ph);
    good_vars(isnp) = !block_info[isnp].ignored && !block_info[isnp].ignored_trait(ph) && !block_info[isnp].test_fail(ph);
  }
  nvars = good_vars.count();

  return (nvars == 0);
}


void JTests::compute_minp(){

  df_test = good_vars.count();
  if( df_test == 0 ) {reset_vals(); return;}

  // get minimum p-value (on log scale)
  double logpmax = good_vars.select(log10pv, 0).maxCoeff();
  get_pv( pow(10, -logpmax) );

}


void JTests::compute_acat(const int& bs, const int& ph, const vector<variant_block>& block_info){

  double v_maf, tmpd;
  boost::math::beta_distribution<>  dist(acat_a1, acat_a2);

  df_test = good_vars.count();
  if( df_test == 0 ) {reset_vals(); return;}

  // make array of weights
  ArrayXd wts = ArrayXd::Zero(bs);
  for(int isnp = 0, j = 0; isnp < bs; isnp++) {
    if( !good_vars(isnp) ) continue;
    // compute weights
    if( valid_snp_mode ) {// sqrt(w)=dbeta(maf,a1,a2)*sqrt(maf*(1-maf))
      v_maf = min( block_info[isnp].af(ph), 1 - block_info[isnp].af(ph) );
      //cerr << v_maf << endl;
      tmpd = pdf( dist, v_maf );
      wts(j++) = v_maf * (1-v_maf) * tmpd * tmpd;
    } else wts(j++) = 1; // assume weight=1
  }

  // get ACAT test stat
  get_pv( get_acat(log10pv, wts) );

}

double get_acat(const Eigen::Ref<const ArrayXd>& logpvals, const Eigen::Ref<const ArrayXd>& weights){

  cauchy dc(0,1);
  double tol = 10.0 * std::numeric_limits<double>::min(), pv_thr = 1e-15;

  // if single pval, return pval
  if(logpvals.size() == 1) return pow(10, -logpvals(0));

  // use approx for small p-values (from ACAT R package)
  ArrayXd pvals = (weights!=0).select( pow(10, -logpvals) , 0.5).max(tol); // to prevent underflow
  //cerr << "log10pv=" << logpvals.matrix().transpose() << "\npv=" << pvals.matrix().transpose() << "\nw=" << weights.matrix().transpose() << "\n";
  double acat = (pvals > pv_thr).select( weights * tan( M_PI * (0.5 - pvals)), (weights / pvals) / M_PI).sum();
  double wsum = weights.sum();
  //cerr << std::setprecision(10) << "acat num=" << acat << " denum=" << wsum << endl;

  return cdf(complement(dc, acat/wsum ));
}

double get_acat(const Eigen::Ref<const ArrayXd>& logpvals){ // uniform weights

  cauchy dc(0,1);
  double tol = 10.0 * std::numeric_limits<double>::min(), pv_thr = 1e-15;

  // if single pval, return pval
  if(logpvals.size() == 1) return pow(10, -logpvals(0));

  // use approx for small p-values (from ACAT R package)
  ArrayXd pvals = pow(10, -logpvals).max(tol); // to prevent underflow
  double acat = (pvals > pv_thr).select( tan( M_PI * (0.5 - pvals)), (1.0 / pvals) / M_PI).sum();
  double wsum = logpvals.size();
  //cerr << std::setprecision(10) << "acat num=" << acat << " denum=" << wsum << endl;

  return cdf(complement(dc, acat/wsum ));
}

void JTests::compute_qr_G(const Eigen::Ref<const MatrixXb>& mask, struct geno_block const* gblock){

  ArrayXi colkeep;
  MatrixXd Gnew;
  indices_vars.resize(0);

  // filter out bad variants
  Gnew = MatrixXd::Zero( gblock->Gmat.rows(), nvars );
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
  else if ( df_test < nvars ){
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


void JTests::compute_nnls(const Eigen::Ref<const MatrixXb>& mask, const Eigen::Ref<const Eigen::MatrixXd>& ymat){

  if( df_test == 0 ) {
    reset_vals();
    return;
  }

  int ns = mask.col(0).array().cast<int>().sum() - ncovars;
  int df_ur = ns - df_test;
  double pval_min2; 
  
  VectorXd y_tmp = ymat.col(0).array() * mask.col(0).array().cast<double>();
  
  // (depreciated) compute NNLS p-value by function
  /* double pval_min2_old = jburden_test(y_tmp, Gtmp, df_ur, nnls_tol, nnls_napprox, true, 3); */

  // initialize an object of NNLS class & pass parameters
  NNLS nnls(nnls_napprox, nnls_normalize, nnls_tol, nnls_strict, nnls_verbose);
  // run the NNLS test: model fitting and inference
  nnls.run(y_tmp, Gtmp, df_ur);
  // get the final p-value = min(NNLS with b>=0, NNLS with b<=0)
  // -1 value means that NNLS failed; check the error message
  pval_min2 = nnls.pval_min2; 
  // get additional p-values & assign
  // to be used downstream if NNLS pos/neg are reported separately
  pval_nnls_pos = nnls.fit_pos.pval;
  pval_nnls_neg = nnls.fit_neg.pval;

  // pvalue
  if((pval_min2 >= 0) & (pval_min2 <= 1)) get_pv( pval_min2 );
  else reset_vals();

  // print extra NNLS information if requested
  if(nnls_verbose_out) {
    string fname = out_file_prefix + "_nnls.info";
    ofstream file_info(fname);
    
    // v1: tabular format
    // header
    file_info << "ID SEL BETA_NNLS BETA_OLS " << endl;
    // variant results per line
    for(int i = 0; i < df_test; i++) {
      unsigned int k = indices_vars[i]; 
      file_info << variant_names[k] << " "; // variant ID
      file_info << nnls.str_sel_i(i) << " "; // NNLS selection status
      file_info << nnls.str_bhat_i(i) << " "; // NNLS beta
      file_info << nnls.bhat_ols[i] << " "; // OLS beta
      file_info << endl; // end of line
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

std::string JTests::run_single_p_acat(int const& bs, const int& chrom, const int& block, int const& ph, const string& pheno_name, std::vector<variant_block>& block_info, std::map<std::string, double>& overall_p, struct geno_block const* gblock, const Eigen::Ref<const Eigen::MatrixXd>& yres, const Eigen::Ref<const MatrixXb>& mask, struct param const* params){

  int ip;
  double max_logp = -1, pv; 
  string mname, max_logp_mask = "NA";
  vector<string> keep_tests = { "ACATV", "SKATO-ACAT" };
  ArrayXd pvals_gene;
  std::map <std::string, double> overall_p_m1;
  std::map <std::string, double>::iterator itr;
  std::ostringstream buffer;

  // overall
  good_vars = false;
  for(int isnp = 0; isnp < bs; isnp++){

    if(block_info[isnp].skip_for_vc) continue; 
    mname = block_info[isnp].mask_name;

    // identify M1 masks
    good_vars(isnp) = (mname == "M1") || (mname == "pLoF");

    for(auto const& extract_test : keep_tests)
      if(in_map(extract_test, block_info[isnp].sum_stats_vc)){
        pv = block_info[isnp].sum_stats_vc[extract_test](ph,1); 
        if(pv>=0){
          overall_p[extract_test + "." + mname] = pv;
          if(pv>max_logp){
            max_logp_mask = mname;
            max_logp = pv;
          }
          if(good_vars(isnp)) overall_p_m1[extract_test + "." + mname] = pv; // have m1 signal
        }
      }

  }

  // combine all p-values and pass through acat
  if(overall_p.size()>0){
    df_test = overall_p.size();
    pvals_gene.resize( df_test );
    ip = 0;
    for (itr = overall_p.begin(); itr !=  overall_p.end(); ++itr) 
      pvals_gene(ip++) = itr->second;
    get_pv( get_acat(pvals_gene) );
    buffer << print_gene_output("GENE_P", max_logp_mask, ph+1, chrom, block, pheno_name, params);
  }

  // for M1 masks
  if(good_vars.any()){
    // run acat only on M1 signals
    compute_acat(bs, ph, block_info);
    if(plog >= 0) overall_p_m1["BURDEN-ACAT"] = plog;
    // run nnls only on M1 signals
    if( CHECK_BIT(test_list, joint_tests_map["nnls"]) ) {
      compute_qr_G(mask, gblock);
      pval_nnls_pos = -1; pval_nnls_neg = -1;
      compute_nnls(mask, yres); 
      if((pval_nnls_pos >= 0) && (pval_nnls_pos <= 1) &&
          (pval_nnls_neg >= 0) && (pval_nnls_neg <= 1) ) {
        get_pv(pval_nnls_pos);overall_p_m1["NNLS_POS"] = plog;
        get_pv(pval_nnls_neg);overall_p_m1["NNLS_NEG"] = plog;
      }; 
    }
    // apply acat to all p for M1
    if(overall_p_m1.size()>0){
      df_test = overall_p_m1.size();
      pvals_gene.resize( df_test );
      ip = 0;
      for (itr = overall_p_m1.begin(); itr !=  overall_p_m1.end(); ++itr) 
        pvals_gene(ip++) = itr->second;
      get_pv( get_acat(pvals_gene) );
      buffer << print_gene_output("GENE_M1_P", "", ph+1, chrom, block, pheno_name, params);
    }
  }

  return buffer.str();
}

string JTests::print_output(const int& ttype, const int& ipheno, const int& chrom, const int& block, const string& pheno_name, struct param const* params){

  if(!params->htp_out) return print_sum_stats(ttype, ipheno, chrom, block, params);
  else return print_sum_stats_htp(ttype, chrom, block, pheno_name, params);

}

// normal regenie format
std::string JTests::print_sum_stats(const int& ttype, const int& ipheno, const int& chrom, const int& block, struct param const* params){

  std::ostringstream buffer;

  // chr pos id a0 a1 af
  if(params->split_by_pheno || ipheno == 1) {
    buffer << setinfo[chrom - 1][block].chrom << " " << setinfo[chrom - 1][block].physpos << " " << setinfo[chrom - 1][block].ID << " NA NA NA " ;
    if( params->af_cc ) buffer << " NA NA ";
    // info
    if(!params->build_mask && params->dosage_mode) buffer << "NA ";
    // n test
    buffer << "NA " << burden_str << test_names[ttype];
  }

  //beta se
  buffer << " NA NA ";

  // chisq
  if( zval != -9 ) buffer << zval << " ";
  else buffer << "NA ";

  // pval
  if( plog != -9 )  buffer << plog << " ";
  else buffer << "NA ";

  //df (print it out only if split by pheno)
  if(params->split_by_pheno || (ipheno == params->n_pheno)) {
    if(params->split_by_pheno && (plog != -9))  buffer << " DF=" << df_test << endl;
    else buffer << "DF=NA\n";
  }

  reset_vals();
  return buffer.str();
}


// htpv4 format
std::string JTests::print_sum_stats_htp(const int& ttype, const int& chrom, const int& block, const string& yname, struct param const* params){

  std::ostringstream buffer;
  bool test_pass = (pval != -9);
  const string cohort = params->cohort_name;

  if( pval == 1) pval = 1 - 1e-7;

  // SNP info
  buffer << setinfo[chrom - 1][block].ID << "\t"<< setinfo[chrom - 1][block].chrom << "\t" << setinfo[chrom - 1][block].physpos << "\tNA\tNA\t";
  // trait, cohort, test name
  buffer << yname << "\t" << params->cohort_name << "\t" << burden_model << test_names[ttype];

  // bhat & 95 CI
  buffer << "\tNA\tNA\tNA\t" ;
  // Pvalue
  if(test_pass) buffer << pval << "\t";
  else buffer << "NA\t";

  // print out AF, counts in cases, counts in controls
  buffer << "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t";

  // info column
  if(test_pass) buffer << "DF=" << df_test;
  else buffer << "DF=0";

  buffer << ";NO_BETA\n";

  reset_vals();
  return buffer.str();

}

// single gene p
string JTests::print_gene_output(const string& mname, const string& max_name, const int& ipheno, const int& chrom, const int& block, const string& pheno_name, struct param const* params){

  if(!params->htp_out) return print_sum_stats_gene(mname, max_name, ipheno, chrom, block, params);
  else return print_sum_stats_htp_gene(mname, max_name, chrom, block, pheno_name, params);

}

// normal regenie format
std::string JTests::print_sum_stats_gene(const string& mname, const string& max_name, const int& ipheno, const int& chrom, const int& block, struct param const* params){

  std::ostringstream buffer;

  // chr pos id a0 a1 af
  if(params->split_by_pheno || ipheno == 1) {
    buffer << setinfo[chrom - 1][block].chrom << " " << setinfo[chrom - 1][block].physpos << " " << setinfo[chrom - 1][block].ID << " NA NA NA " ;
    if( params->af_cc ) buffer << " NA NA ";
    // info
    if(params->dosage_mode) buffer << "NA ";
    // n test
    buffer << "NA " << mname;
  }

  //beta se
  buffer << " NA NA ";

  // chisq
  if( zval != -9 ) buffer << zval << " ";
  else buffer << "NA ";

  // pval
  if( plog != -9 )  buffer << plog << " ";
  else buffer << "NA ";

  //df (print it out only if split by pheno)
  if(params->split_by_pheno || (ipheno == params->n_pheno)) {
    if(params->split_by_pheno && (plog != -9))  buffer << " DF=" << df_test;
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
std::string JTests::print_sum_stats_htp_gene(const string& mname, const string& max_name, const int& chrom, const int& block, const string& yname, struct param const* params){

  std::ostringstream buffer;
  bool test_pass = (pval != -9);
  const string cohort = params->cohort_name;

  if( pval == 1) pval = 1 - 1e-7;

  // SNP info
  buffer << setinfo[chrom - 1][block].ID << "\t"<< setinfo[chrom - 1][block].chrom << "\t" << setinfo[chrom - 1][block].physpos << "\tNA\tNA\t";
  // trait, cohort, test name
  buffer << yname << "\t" << params->cohort_name << "\t" << mname;

  // bhat & 95 CI
  buffer << "\tNA\tNA\tNA\t" ;
  // Pvalue
  if(test_pass) buffer << pval << "\t";
  else buffer << "NA\t";

  // print out AF, counts in cases, counts in controls
  buffer << "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t";

  // info column
  if(test_pass) buffer << "DF=" << df_test;
  else buffer << "DF=0";

  // top signal
  if(max_name != "") buffer << ";STRONGEST_MASK=" << max_name;

  buffer << ";NO_BETA\n";
  reset_vals();
  return buffer.str();

}

void JTests::get_variant_names(int const& chrom, int const& block, vector<snp> const& snpinfo){

  // only for NNLS (for now)
  if( !CHECK_BIT(test_list,3) || !nnls_verbose_out ) return;

  int bs = setinfo[chrom - 1][block].snp_indices.size();
  uint32_t index;
  variant_names.resize(bs);

  for(int i = 0; i < bs; i++){
    index = setinfo[chrom - 1][block].snp_indices[ i ];
    variant_names[i] = snpinfo[index].ID;
  }

}

void JTests::reset_vals(){
  pval = -9, plog = -9, zval = -9;
}

void JTests::get_pv(const double& pv){

  chi_squared chisq(1);

  pval = max(nl_dbl_dmin, pv); // to prevent underflow
  zval = quantile(complement(chisq, pval)); // chisq stat
  plog = -log10(pval); // -log10p

}

