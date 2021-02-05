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
#include "Pheno.hpp"
#include "NNLS.hpp"
#include "Files.hpp"
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


void JTests::get_test_info(const struct param* params, string test_string, mstream& sout){

  std::vector< string > tmp_str_vec ;
  test_list = 0u;

  burden_str = test_string + "-BURDEN-";
  if(params->htp_out){
    if(params->skip_blups) burden_model = burden_str;
    else burden_model = test_string + "-WGR-BURDEN-";
  }

  // activate tests chosen by user
  tmp_str_vec = string_split(params->burden,",");
  for( size_t i = 0; i < tmp_str_vec.size(); i++){
    switch ( test_ind(tmp_str_vec[i]) ){

      case 0:
        BIT_SET(test_list, 0);
        break;

      case 1:
        if(params->binary_mode) sout << "WARNING: Joint F-test only for QTs.\n";
        else BIT_SET(test_list, 1);
        break;

      case 2:
        BIT_SET(test_list, 2);
        break;

      case 3:
        if(params->binary_mode) sout << "WARNING: Joint NNLS test only for QTs.\n";
        else { 
          BIT_SET(test_list, 3); 
          nnls_napprox = params->nnls_napprox; 
          nnls_verbose_out = params->nnls_out_all;
        }
        break;

      case 4:
        BIT_SET(test_list, 4);
        valid_snp_mode = !params->build_mask || (params->build_mask && params->mask_rule_max) ;
        acat_a1 = params->acat_a1;
        acat_a2 = params->acat_a2;
        break;

      default:
        sout << "ERROR: Unrecognized joint test (=" << tmp_str_vec[i] << ")" <<endl;
        exit(EXIT_FAILURE);

    }
  }

  ncovars = params->ncov;

}

string JTests::apply_joint_test(const int chrom, const int block, const int ph, struct phenodt* pheno_data, const Eigen::Ref<const Eigen::MatrixXd>& yres, struct geno_block* gblock, std::vector<variant_block>& block_info, const string pheno_name, struct param* params){

  int bs = setinfo[chrom - 1][block].snp_indices.size();
  std::ostringstream buffer;
  reset_vals();

  if(nvars == 0 || set_vars(bs, ph, block_info)) 
    return ""; // if no variants in set passed filters
  
  if( CHECK_BIT(test_list,0) ) { // minP
    compute_minp(bs, ph, block_info);
    buffer << print_output(0, chrom, block, pheno_name, params);
  } 
  if( CHECK_BIT(test_list,4) ) { // ACAT
    compute_acat(bs, ph, block_info);
    buffer << print_output(4, chrom, block, pheno_name, params);
  } 

  // check other test
  if( ((int) (test_list&qr_tests)) > 0 ) {

    compute_qr_G(pheno_data->masked_indivs.col(ph), gblock);

    if( CHECK_BIT(test_list,1) ) { // F-test
      compute_ftest(pheno_data->masked_indivs.col(ph), yres); 
      buffer << print_output(1, chrom, block, pheno_name, params);
    } 
    if( CHECK_BIT(test_list,2) ) { // GATES
      compute_gates(ph, block_info);
      buffer << print_output(2, chrom, block, pheno_name, params);
    } 
    if( CHECK_BIT(test_list,3) ) { // NNLS
      compute_nnls(pheno_data->masked_indivs.col(ph), yres); 
      if(!nnls_verbose_out) { 
        // default output
        buffer << print_output(3, chrom, block, pheno_name, params);
      } else {
        // verbose output with NNLS pos & neg split into two 
        // 1. NNLS pos (test code 5)
        if((pval_nnls_pos >= 0) & (pval_nnls_pos <= 1)) get_pv(pval_nnls_pos);
        else reset_vals();
        buffer << print_output(5, chrom, block, pheno_name, params);
        // 2. NNLS neg (test code 6)
        if((pval_nnls_neg >= 0) & (pval_nnls_neg <= 1)) get_pv(pval_nnls_neg);
        else reset_vals();
        buffer << print_output(6, chrom, block, pheno_name, params);
      }
    }
  }

  return buffer.str();
}


// determine if marginal test failed 
bool JTests::set_vars(const int bs, const int ph, std::vector<variant_block>& block_info){

  int ngood = 0;
  good_vars.resize(bs);

  for(int isnp = 0; isnp < bs; isnp++) {
    good_vars[isnp] = !block_info[isnp].ignored && !block_info[isnp].ignored_trait(ph) && !block_info[isnp].test_fail[ph];
    if( good_vars[isnp] ) ngood++;
  }
  nvars = ngood;

  return (ngood == 0);
}


void JTests::compute_minp(const int bs, const int ph, const std::vector<variant_block>& block_info){

  df_test = 0; 

  // get minimum p-value
  for(int isnp = 0; isnp < bs; isnp++) {
    if( !good_vars[isnp] ) continue;

    if(block_info[isnp].pval_log(ph) > plog) 
      plog = block_info[isnp].pval_log(ph);

    df_test++;
  }

  if( df_test == 0 ) reset_vals();
  else get_pv( pow(10, -plog) );

}


void JTests::compute_acat(const int bs, const int ph, const vector<variant_block>& block_info){

  double acat = 0, tval, v_p, wt, wsum = 0, v_maf, tmpd;
  df_test = 0;
  boost::math::beta_distribution<>  dist(acat_a1, acat_a2);
  cauchy dc(0,1);
    

  // get ACAT test stat
  for(int isnp = 0; isnp < bs; isnp++) {
    if( !good_vars[isnp] ) continue;

    // compute weights
    if( valid_snp_mode ) {// sqrt(w)=dbeta(maf,a1,a2)*sqrt(maf*(1-maf))
      v_maf = min( block_info[isnp].af(ph), 1 - block_info[isnp].af(ph) );
      //cerr << v_maf << endl;
      tmpd = pdf( dist, v_maf );
      wt = v_maf * (1-v_maf) * tmpd * tmpd;
    } else wt = 1; // assume weight=1

    // use approx for small p-values (from ACAT R package)
    v_p = max(nl_dbl_dmin, pow(10, -block_info[isnp].pval_log(ph))); // to prevent underflow
    if(v_p > acat_pv_thr)
      tval = wt * tan( (0.5 - v_p) * M_PI );
    else 
      tval = (wt / v_p) / M_PI ;
    //cerr << "m:" << v_maf << ";p:" << v_p << ";w:" << wt << ";t:" << tval << endl;

    acat += tval;
    wsum += wt;

    df_test++;
  }

  if( df_test == 0 ) reset_vals();
  else get_pv( cdf(complement(dc, acat/wsum )) );

}


void JTests::compute_qr_G(const Eigen::Ref<const MatrixXb>& mask, struct geno_block* gblock){

  ArrayXi colkeep;
  MatrixXd Gnew;
  indices_vars.resize(0);

  // filter out bad variants
  Gnew = MatrixXd::Zero( gblock->Gmat.rows(), nvars );
  for(int i = 0, j = 0; i < gblock->Gmat.cols(); i++){
    if(!good_vars[i]) continue;
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

void JTests::compute_gates(const int ph, const std::vector<variant_block>& block_info){

  int gcol;
  double p_gates, m_e, p_i, m_ei;

  if( df_test == 0 ) {
    reset_vals();
    return;
  } else if( df_test == 1 ){
    p_gates = pow(10, -block_info[ indices_vars[0] ].pval_log(ph));
    if( p_gates >= 0) get_pv( p_gates );
    else reset_vals();
    return;
  }

  // Sort p-values
  vector<double> pvals, sorted_pv;
  MatrixXd tmpG = Gtmp;

  pvals.resize(df_test);
  for(int i = 0; i < df_test; i++) {
    //cerr << block_info[ indices_vars[i] ].pval_log(ph) ;
    pvals[i] = pow(10, -block_info[ indices_vars[i] ].pval_log(ph));
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

string JTests::print_output(const int ttype, const int chrom, const int block, const string pheno_name, struct param* params){

  if(!params->htp_out) return print_sum_stats(ttype, chrom, block, params);
  else return print_sum_stats_htp(ttype, chrom, block, pheno_name, params);

}


// normal regenie format
std::string JTests::print_sum_stats(const int ttype, const int chrom, const int block, struct param* params){

  std::ostringstream buffer;

  buffer << setinfo[chrom - 1][block].chrom << " " << setinfo[chrom - 1][block].physpos << " " << setinfo[chrom - 1][block].ID << " NA NA NA " ;
  if(params->dosage_mode) buffer << "NA ";
  buffer << "NA " << burden_str << test_names[ttype] << " NA NA ";

  if( zval != -9 ) buffer << zval << " ";
  else buffer << "NA ";

  if( plog != -9 )  buffer << plog << " " << df_test;
  else buffer << "NA 0";

  buffer << endl;

  reset_vals();
  return buffer.str();
}


// htpv4 format
std::string JTests::print_sum_stats_htp(const int ttype, const int chrom, const int block, const string yname, struct param* params){

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

  buffer << endl;

  reset_vals();
  return buffer.str();

}

void JTests::get_variant_names(int chrom, int block, vector<snp>& snpinfo){

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

void JTests::get_pv(const double pv){

  chi_squared chisq(1);

  pval = max(nl_dbl_dmin, pv); // to prevent underflow
  zval = quantile(complement(chisq, pval)); // chisq stat
  plog = -log10(pval); // -log10p

}

int JTests::test_ind(const string test_str){

  if(test_str == "minp") return 0;
  else if(test_str == "ftest") return 1;
  else if(test_str == "gates") return 2;
  else if(test_str == "nnls") return 3;
  else if(test_str == "acat") return 4;

  return -1;
}
