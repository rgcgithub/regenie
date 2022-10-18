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
#include "Joint_Tests.hpp"
#include "Step1_Models.hpp"
#include "Step2_Models.hpp"
#include "Masks.hpp"


using namespace std;
using namespace Eigen;

GenoMask::GenoMask() { // @suppress("Class members should be properly initialized")
}

GenoMask::~GenoMask() {
  // TODO Auto-generated destructor stub
}

void GenoMask::prep_run(struct param& params, struct in_files const& files){

  params.min_MAC_mask = params.min_MAC; // for association tests
  params.min_MAC = 0.5; // set this so can retain singletons (0.5 for dosages)
  take_max = params.mask_rule_max;
  take_comphet = params.mask_rule_comphet;
  w_loo = params.mask_loo;
  w_lodo = params.mask_lodo;
  w_vc_tests = params.vc_test;
  vc_aaf = params.vc_maxAAF;
  vc_collapse_MAC = params.skat_collapse_MAC;
  w_vc_cust_weights = params.vc_with_weights;
  write_masks = params.write_masks;
  write_snplist = params.write_mask_snplist;
  force_singleton = params.aaf_file_wSingletons;
  verbose = params.verbose || params.debug;

  if(!take_max && !take_comphet) params.htp_out = false; // due to genocounts with sum rule
  if(write_masks) gfile_prefix = files.out_file + "_masks";

}

void GenoMask::setBins(struct param* params, mstream& sout){

  vector<double> tmpv;

  // check if not singleton when using LOO
  if(!params->mask_loo || ( params->mask_loo && (params->mbins[0] != "singleton")) ){

    if(params->mbins.size() >= 1){
      // convert them to double
      for( size_t i = 0; i < params->mbins.size(); i++){
        if(params->mbins[i] == "all")
          tmpv.push_back( 1 );
        else
          tmpv.push_back( convertDouble( params->mbins[i], params, sout) );
      }
    } else tmpv.push_back( default_aaf );

    if(w_vc_tests) tmpv.push_back( vc_aaf );

    // sort and retain unique values
    std::sort(tmpv.begin(), tmpv.end());
    tmpv.erase( unique( tmpv.begin(), tmpv.end() ), tmpv.end() );
    // store in eigen object
    aafs = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(tmpv.data(), tmpv.size());

    // check validity
    if( ( (aafs.array() < minAAF) || (aafs.array() > 1) ).any() && !w_vc_tests )
      throw "must specify values for --aaf-bins in [" + to_string( minAAF ) + ", 1]";
    else if( (aafs.array()>= 0.5).any() && !w_vc_tests )
      sout << "WARNING: For computational efficiency, it is recommended that AAF cutoffs < 0.5\n";

    max_aaf = aafs.tail(1)(0);
  }

  n_aaf_bins = aafs.size();

  if(n_aaf_bins > max_aaf_bins)
    throw "Number of AAF bins (=" + to_string( n_aaf_bins ) + ") above maximum (=" + to_string( max_aaf_bins ) + ")\n";

  IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "", "", "","","");
  if(!params->mask_loo)
    sout << left << std::setw(20) << " * aaf cutoffs" << ": [ " << n_aaf_bins << " : "
      << aafs.transpose().format(Fmt) << " ] + singletons\n";
  else {
    sout << left << std::setw(20) << " * aaf cutoff" << ": " ;
    if(n_aaf_bins == 0) sout << "singleton";
    else sout << aafs(0);
    sout << endl;
  }

  n_aaf_bins++; // add singleton bin

  nmasks_total = n_aaf_bins * masks.size(); // total number of masks

  w_regions = params->w_regions;
  if(w_regions) base_masks = masks;
}

void GenoMask::prepMasks(int const& ntotal, const string& setID) {

  maskinfo tmp_region_mask;
  std::map <std::string, uint16_t>::iterator itr;

  // make new set of masks if using set regions
  if(w_regions){ 
    masks.resize(0);
    // go through each original mask and create region specific mask
    for(size_t i = 0; i < base_masks.size(); i++ ){
      tmp_region_mask = base_masks[i];
      for (itr = regions[setID].begin(); itr != regions[setID].end(); ++itr) { // make region mak
        if(w_lodo){ // LODO scheme
          tmp_region_mask.region_name = "LODO_" + itr->first + ".";
          tmp_region_mask.region = (65535 & ~itr->second); // unset bits for region [2-byte]
          masks.push_back(tmp_region_mask);
        } else {
          tmp_region_mask.region_name = itr->first + ".";
          tmp_region_mask.region = itr->second;
          masks.push_back(tmp_region_mask);
        }
      }
      if(!w_loo){// add mask across all regions
        tmp_region_mask = base_masks[i];
        tmp_region_mask.region |= 65535; //set all 16 bits
        masks.push_back(tmp_region_mask);
      }
    }
    nmasks_total = n_aaf_bins * masks.size();
    if(write_masks) reset_gvec();
  } 

  Gtmp = MatrixXd::Constant(ntotal, nmasks_total, -3);
  colset = ArrayXb::Constant( nmasks_total, false );
  if(!take_max) {
    non_missing = MatrixXb::Constant(ntotal, nmasks_total, false);
    if(!take_comphet) nsites = ArrayXi::Constant(nmasks_total, 0);
  }

  if(write_setlist) {
    list_masks.resize(nmasks_total);
    for(size_t i = 0;i < list_masks.size(); i++)
      list_masks[i].resize(0);
  }
  if(write_snplist) {
    list_snps.resize(nmasks_total);
    for(size_t i = 0;i < list_snps.size(); i++)
      list_snps[i].resize(0);
  }

}

void GenoMask::updateMasks(int const& start, int const& bs, struct param* params, struct filter* filters, const Ref<const MatrixXb>& masked_indivs, struct geno_block* gblock, const Ref<const ArrayXd>& vc_weights, vector<variant_block> &all_snps_info, vset& setinfo, vector<snp>& snpinfo, mstream& sout){

  // identify which snps are in each mask
  set_snp_masks(start, bs, all_snps_info, setinfo, snpinfo, sout);
  // identify which snps are in each aaf bin
  set_snp_aafs(start, bs, params->set_aaf, all_snps_info, setinfo, snpinfo, sout);

  MatrixXb Jmat, ur_miss;
  MatrixXd rare_mask_tmp;
  SpMat ur_sp_mat;
  ArrayXi ur_indices;
  if(w_vc_tests) {
    Jmat = MatrixXb::Constant(bs, nmasks_total, false);
    if(setinfo.ultra_rare_ind.segment(start, bs).any()) {
      int n_ur = setinfo.ultra_rare_ind.segment(start, bs).count();
      rare_mask_tmp = setinfo.vc_rare_mask; // not safe to update SpMat in parallel (not many columns)
      ur_indices = ArrayXi::Constant(bs, -1);
      ur_sp_mat.resize(params->n_samples, n_ur);
      ur_miss.resize(params->n_samples, n_ur);

      // store the ur variants in spmat & keep track of index/missingness
      for(int i = 0, j = 0; i < bs; i++){
        if(!setinfo.ultra_rare_ind(start+i)) continue;
        MapArXd garr (gblock->Gmat.col(i).data(), params->n_samples, 1);
        // flip if necessary
        if(all_snps_info[start+i].af1 > 0.5) ur_sp_mat.col(j) = (garr == -3).select(0, 2 - garr).matrix().sparseView();
        else ur_sp_mat.col(j) = (garr < 0).select(0, garr).matrix().sparseView();
        // if using custom user weights, rescale before collapsing ur variants 
        ur_sp_mat.col(j) *= vc_weights(start+i);
        ur_miss.col(j) = (garr >= 0);
        // store the index
        ur_indices(i) = j++;
      }
    }

  }

  // update each mask 
#if defined(_OPENMP)
  setNbThreads(1);
  // use MT in both loops
#pragma omp parallel for schedule(dynamic) collapse(2)
#endif
  for(size_t i = 0; i < masks.size(); i++){
    for(int j = 0; j < n_aaf_bins; j++){

      int index_start = i * n_aaf_bins + j;
      ArrayXb colkeep = keepmask.col(i).array() && keepaaf.col(j).array();
      if(!take_max && !take_comphet) nsites(index_start) += colkeep.count();
      if(write_snplist) append_snplist(index_start, colkeep, start, setinfo, snpinfo);

      if(w_vc_tests && (j > 0) && ( aafs(j-1) == vc_aaf )) // track variants in mask
        Jmat.col(index_start) = colkeep.matrix();

      // ignore variants in previous AAF categories (accumulation is in next loop)
      if(j>0) colkeep = colkeep && !keepaaf.col(j-1).array(); 

      // if there are no variants included, continue
      if( colkeep.count() == 0 ) continue;
      //if(i==2 && j==1) cerr << i << " " << j << " " << colkeep.count() << " ";

      // update mask
      MapArXd maskvec (Gtmp.col(index_start).data(), params->n_samples, 1);

      if(take_max) {

        SpVec gv, mv;
        mv = maskvec.matrix().sparseView();
        for(int k = 0; k < colkeep.size(); k++){
          if(!colkeep(k)) continue;
          gv = (vc_weights(start+k) * gblock->Gmat.col(k)).sparseView();
          mv = gv.cwiseMax(mv);
        }
        maskvec = MatrixXd(mv).array();

      } else {

        int l;
        double ds;
        SpVec gv;

        for(int k = 0; k < colkeep.size(); k++){
          if(!colkeep(k)) continue;
          gv = (vc_weights(start+k) * gblock->Gmat.col(k)).sparseView();

          // sum rule (ignore -3)
          for (SparseVector<double>::InnerIterator it(gv); it; ++it) {
            l = it.index();
            ds = it.value();

            if( !filters->ind_in_analysis(l) || (ds < 0)) continue;

            if( maskvec(l) < 0 ) maskvec(l) = ds;
            else maskvec(l) += ds;
          }

          // for genotype counts, identify when (-3) is 0
          non_missing.col(index_start).array() = non_missing.col(index_start).array() || ( gblock->Gmat.col(k).array() >= 0 );
        }
       //if(i==2 && j==1) cout << (maskvec == -3).select(0,maskvec).sum() << endl;

      }

      // get ultra-rare mask if using VC test (take max)
      if(w_vc_tests && setinfo.ultra_rare_ind.segment(start, bs).any() && ((j == 0) || ( aafs(j-1) <= vc_aaf )) ) {
        SpVec mv = setinfo.vc_rare_mask.col(index_start);
        for(int k = 0; k < colkeep.size(); k++){
          if(!colkeep(k) || !setinfo.ultra_rare_ind(start+k)) continue;
          mv = mv.cwiseMax( ur_sp_mat.col(ur_indices(k)) );
          setinfo.vc_rare_mask_non_missing.col(index_start).array() = setinfo.vc_rare_mask_non_missing.col(index_start).array() || ur_miss.col( ur_indices(k) ).array();
        }
        rare_mask_tmp.col(index_start) = mv;
      } 

    }
  }

#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

  if(w_vc_tests) {
    setinfo.Jmat.middleRows(start, bs) = Jmat;
    if(setinfo.ultra_rare_ind.segment(start, bs).any()) setinfo.vc_rare_mask = rare_mask_tmp.sparseView();
  }

}

void GenoMask::apply_rule(SpVec& out_mask, SpVec const& Gvec, const Ref<const ArrayXb>& in_analysis, bool const& force_max) { 

  int l;
  double ds;

  if(take_max || force_max) { // max rule to combine variants across sites
    out_mask = Gvec.cwiseMax(out_mask);
  } else { // sum rule (ignore missing) 
    MatrixXd out_mask_vec = out_mask;
    for (SpVec::InnerIterator it(Gvec); it; ++it) {
      l = it.index();
      ds = it.value();

      if( !in_analysis(l) || (ds == -3)) continue;

      if( out_mask_vec(l,0) == -3 ) out_mask_vec(l,0) = ds;
      else if(ds > 0) out_mask_vec(l,0) += ds;
    }
    out_mask = out_mask_vec.col(0).sparseView();
  }

}

void GenoMask::apply_rule(Ref<ArrayXd> out_mask, SpVec const& Gvec, const Ref<const ArrayXb>& in_analysis, bool const& force_max) { 

  int l;
  double ds;

  if(take_max || force_max) { // max rule to combine variants across sites
    SpVec tmpv = out_mask.matrix().sparseView();
    out_mask = Gvec.cwiseMax(tmpv);
    out_mask = in_analysis.select(out_mask, -3);;
  } else { // sum rule (ignore missing) 
    for (SpVec::InnerIterator it(Gvec); it; ++it) {
      l = it.index();
      ds = it.value();

      if( !in_analysis(l) || (ds<0)) continue;

      if( out_mask(l) < 0 ) out_mask(l) = ds;
      else if(ds > 0) out_mask(l) += ds;
    }
  }

}

void GenoMask::apply_rule(Ref<ArrayXd> maskvec, const Ref<const MatrixXd>& Gmat, const Ref<const ArrayXb>& in_analysis, bool const& force_max) { 

  if(take_max || force_max) { // max rule to combine variants across sites
    maskvec = in_analysis.select(maskvec.max(Gmat.rowwise().maxCoeff().array()), maskvec);
  } else { // sum rule (ignore missing) 
    ArrayXb non_miss_G = in_analysis && (Gmat.array() >= 0).rowwise().any();
    maskvec = non_miss_G.select( maskvec.max(0) + (Gmat.array() >= 0).select(Gmat.array(), 0).rowwise().sum(), maskvec);
  }

}

// should only be called once
void GenoMask::collapse_mask_chunk(const Ref<const ArrayXi>& indices, SpMat const& Gmat_sp, const Ref<const ArrayXb>& is_ultra_rare, const Ref<const ArrayXb>& to_flip, const Ref<const ArrayXd>& vc_weights, Ref<ArrayXd> out_mask, Ref<ArrayXd> out_ur_mask, const Ref<const ArrayXb>& in_analysis){ 

  int nkept = indices.size(), icol;
  double weight;
  if(nkept == 0) return;

  // collapse variants
  for(int i = 0; i < nkept; i++){
    icol = indices(i);
    weight = vc_weights(icol); // apply custom user weight to ur variant

    if( w_vc_tests && is_ultra_rare(icol) ){ // collapse into a rare mask
      if( to_flip(icol) ){ // need to flip 
        ArrayXd Gvec = Gmat_sp.col(icol);
        Gvec = (in_analysis && (Gvec >=0 )).select(2 - Gvec, Gvec);
        SpVec G_flip = Gvec.matrix().sparseView();
        apply_rule(out_ur_mask, weight * G_flip, in_analysis, true);
      } else
        apply_rule(out_ur_mask, weight * Gmat_sp.col(icol), in_analysis, true);
    }  
    // for lovo mask
    apply_rule(out_mask, weight * Gmat_sp.col(icol), in_analysis, false);

  }

}

void GenoMask::updateMasks_loo(const Ref<const ArrayXi>& indices_chunk, bool const& comp_full_mask, SpMat const& Gmat_sp, const Ref<const ArrayXb>& is_ultra_rare, const Ref<const ArrayXb>& to_flip, const Ref<const ArrayXd>& vc_weights, const Ref<const ArrayXd>& excl_vars_mask, const Ref<const ArrayXd>& excl_vars_ur_mask, const Ref<const ArrayXb>& in_analysis, vset& setinfo, int const& nthreads){

  bool with_ur = is_ultra_rare(indices_chunk).any() || (excl_vars_ur_mask > 0).any();
  int bs = indices_chunk.size();
  Gtmp.resize(Gmat_sp.rows(), nmasks_total); // incl. full mask

  vector<SpVec> rare_mask_tmp;
  ArrayXb g_ur_nmiss;
  SpVec g_ur_start;

  // store matrix as dense (should be fairly small)
  SpMat Jstar (Gmat_sp.cols(), bs); // Mall x M
  Jstar.reserve(bs);
  for(int i = 0; i < bs; i++)
    Jstar.insert(indices_chunk(i), i) = 1;
  MatrixXd Gmat_d = Gmat_sp * Jstar;

  if(w_vc_tests && with_ur) { // for skat tests and ur vars present
    setinfo.vc_rare_mask.resize(Gtmp.rows(), Gtmp.cols());
    setinfo.vc_rare_mask_non_missing.resize(Gtmp.rows(), Gtmp.cols());
    g_ur_nmiss = excl_vars_ur_mask >= 0;
    g_ur_start = excl_vars_ur_mask.max(0).matrix().sparseView();
    setinfo.vc_rare_mask_non_missing.array().colwise() = g_ur_nmiss;
    rare_mask_tmp.assign(Gtmp.cols(), g_ur_start);
  }
  Gtmp.array().colwise() = excl_vars_mask;

  // start openmp for loop
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(int i = 0; i < bs; i++){ // generate lovo mask

    MapArXd maskvec (Gtmp.col(i).data(), Gtmp.rows(), 1);
    ArrayXb in_loo_mask = ArrayXb::Constant(bs, true); in_loo_mask(i) = false;
    ArrayXi in_loo_indices = get_true_indices(in_loo_mask);

    if(in_loo_mask.any()){
      // LOVO mask
      if(take_max){ // max rule to combine variants across sites
        maskvec = in_analysis.select(maskvec.max((Gmat_d(all, in_loo_indices) * vc_weights(indices_chunk(in_loo_indices)).matrix().asDiagonal()).rowwise().maxCoeff().array()), maskvec);
      } else {
        ArrayXb non_miss_G = in_analysis && (Gmat_d(all, in_loo_indices).array() >= 0).rowwise().any();
        maskvec = non_miss_G.select( maskvec.max(0) + (Gmat_d(all, in_loo_indices).array() >= 0).select((Gmat_d(all, in_loo_indices) * vc_weights(indices_chunk(in_loo_indices)).matrix().asDiagonal()).array(), 0).rowwise().sum(), maskvec);
      }
    }

    if( comp_full_mask && (i == (bs-1)) ) { // compute full mask
      Gtmp.rightCols(1) = Gtmp.col(i); 
      // add to burden mask for full set
      apply_rule(Gtmp.rightCols(1).array(), (vc_weights(indices_chunk(i)) * Gmat_d.col(i)).sparseView(), in_analysis, false);
    }

  }
#if defined(_OPENMP)
  setNbThreads(nthreads);
#endif

  if(w_vc_tests && with_ur) { // for skat tests and ur vars present

    // get spmat for ur variants only
    int n_ur = is_ultra_rare(indices_chunk).count();
    ArrayXi ur_indices = ArrayXi::Constant(bs, -1);
    MatrixXb ur_miss; if(n_ur>0) ur_miss.resize(Gmat_d.rows(), n_ur);
    SpMat ur_sp_mat; if(n_ur>0) ur_sp_mat.resize(Gmat_d.rows(), n_ur);
    SpVec ur_mask_all = g_ur_start;

    // store the ur variants in spmat & keep track of index/missingness
    int m = 0;
    for (int i = 0; i < bs; i++ ) {
      if(!is_ultra_rare(indices_chunk(i))) continue;
      MapArXd garr (Gmat_d.col(i).data(), Gmat_d.rows(), 1);
      // flip if necessary (missing set to 0)
      if(to_flip(indices_chunk(i))) ur_sp_mat.col(m) = (garr<0).select(0, 2 - garr).matrix().sparseView();
      else ur_sp_mat.col(m) = garr.max(0).matrix().sparseView();
      ur_sp_mat.col(m) *= vc_weights(indices_chunk(i)); // apply custom user weight to ur variant
      // track missingness
      ur_miss.col(m) = (garr >= 0);
      // track max across sites
      ur_mask_all = ur_mask_all.cwiseMax(ur_sp_mat.col(m));
      // store the index
      ur_indices(i) = m++;
    }

    // start openmp for loop
#if defined(_OPENMP)
    setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
    for(int i = 0; i < bs; i++){

      ArrayXb in_loo_mask = ArrayXb::Constant(bs, true); in_loo_mask(i) = false;

      if(n_ur>0){
        // if not ur site, set to max across ur sites
        if(!is_ultra_rare(indices_chunk(i))) {
          rare_mask_tmp[i] = ur_mask_all;
          setinfo.vc_rare_mask_non_missing.col(i).array() = setinfo.vc_rare_mask_non_missing.col(i).array() || ur_miss.array().rowwise().any();
        } else // otherwise, take max using lovo
          for (auto const& j : get_true_indices( in_loo_mask && is_ultra_rare(indices_chunk) )) {
            rare_mask_tmp[i] = rare_mask_tmp[i].cwiseMax(ur_sp_mat.col(ur_indices(j)));
            setinfo.vc_rare_mask_non_missing.col(i).array() = setinfo.vc_rare_mask_non_missing.col(i).array() || ur_miss.col(ur_indices(j)).array() ;
          }
      }

      if( comp_full_mask && (i == (bs-1)) ) { // compute full mask
        rare_mask_tmp.back() = rare_mask_tmp[i];
        setinfo.vc_rare_mask_non_missing.rightCols(1) = setinfo.vc_rare_mask_non_missing.col(i);
        if(is_ultra_rare(indices_chunk(i))) { // if variant is UR
          rare_mask_tmp.back() = rare_mask_tmp.back().cwiseMax(ur_sp_mat.col(ur_indices(i)));
          setinfo.vc_rare_mask_non_missing.rightCols(1).array() = setinfo.vc_rare_mask_non_missing.rightCols(1).array() || ur_miss.col(ur_indices(i)).array();
        }
      }

    }
#if defined(_OPENMP)
    setNbThreads(nthreads);
#endif

    for(size_t j = 0; j < rare_mask_tmp.size(); j++)
      setinfo.vc_rare_mask.col(j) = rare_mask_tmp[j];

  }

}


void GenoMask::tally_masks(struct param const* params, struct filter const* filters, const Ref<const MatrixXb>& masked_indivs, SpMat& vc_rare_mask, MatrixXb& vc_rare_mask_non_missing){

  MatrixXd rare_mask_tmp;
  if(w_vc_tests) rare_mask_tmp = vc_rare_mask; // not safe to update SpMat in parallel (not many columns)

#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  // get mask by aggregating across increasing AAF categories
  for(size_t i = 0; i < masks.size(); i++){
    // first mask AAF
    int index_mask = i * n_aaf_bins;
    double ds;
    bool column_set = (Gtmp.col(index_mask).array() >= 0).any();
    colset(index_mask) = column_set;

    // don't parallelize inner loop (cumulative updates)
    for(int j = 1; j < n_aaf_bins; j++){

      int index_start = index_mask + j;
      // check if there are variants in mask
      if(colset(index_start-1) || (Gtmp.col(index_start).array() >= 0).any()) 
        colset(index_start) = true;

      if( !colset(index_start) ) continue;

      // add mask from previous AAF category
      MapArXd maskvec (Gtmp.col(index_start).data(), params->n_samples, 1);

      if(take_max) {

        SpMat gv, mv;
        // sparse of current
        mv = maskvec.matrix().sparseView();
        // sparse of previous
        gv = Gtmp.col(index_start-1).sparseView();
        // aggregate to current
        maskvec = MatrixXd(gv.cwiseMax(mv)).array();

      } else {

        SpVec gv;

        // add previous
        gv = Gtmp.col(index_start-1).sparseView();

        // sum rule (ignore -3)
        for (SparseVector<double>::InnerIterator it(gv); it; ++it) {
          int l = it.index();
          ds = it.value();

          if( !filters->ind_in_analysis(l) || (ds < 0)) continue;

          if( maskvec(l) < 0 ) maskvec(l) = ds;
          else maskvec(l) += ds;

        }

        // for genotype counts, identify when (-3) is 0
        non_missing.col(index_start).array() = non_missing.col(index_start).array() || non_missing.col(index_start-1).array();

      }

      // update ultra-rare mask
      if( w_vc_tests && ( aafs(j-1) <= vc_aaf ) && rare_mask_tmp.col(index_start-1).nonZeros() ){
        rare_mask_tmp.col(index_start) = rare_mask_tmp.col(index_start).cwiseMax(rare_mask_tmp.col(index_start-1));
        vc_rare_mask_non_missing.col(index_start).array() = vc_rare_mask_non_missing.col(index_start).array() || vc_rare_mask_non_missing.col(index_start-1).array();
        if( aafs(j-1) == vc_aaf ){
         for( int k = 1; k <= j; k++){ // remove data not needed anymore
          rare_mask_tmp.col(index_start-k).array() = 0;
          vc_rare_mask_non_missing.col(index_start-k).array() = false;
         }
        }
      }
      
    }
  }

#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

  if(w_vc_tests) vc_rare_mask = rare_mask_tmp.sparseView();

}

void GenoMask::computeMasks(struct param* params, struct filter* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& ymat, struct geno_block* gblock, vector<variant_block> &all_snps_info, vset& setinfo, vector<snp>& snpinfo, mstream& sout){

  // check size of vblock vector
  if(((int) all_snps_info.size()) < nmasks_total) all_snps_info.resize( nmasks_total );

  // tally masks
  tally_masks(params, filters, masked_indivs, setinfo.vc_rare_mask, setinfo.vc_rare_mask_non_missing);

  ArrayXb in_bed = colset;

  // finish building each mask 
#if defined(_OPENMP)
  setNbThreads(1);
  // use MT in both loops
#pragma omp parallel for schedule(dynamic) collapse(2)
#endif
  for(size_t i = 0; i < masks.size(); i++){
    for(int j = 0; j < n_aaf_bins; j++){

      int index_start = i * n_aaf_bins + j;

      // check variants were included in mask
      if(!colset(index_start)) continue;

      // compute mask
      buildMask(index_start, setinfo.chrom, params, filters, masked_indivs, ymat, &all_snps_info[index_start]);

      colset(index_start) = !all_snps_info[index_start].ignored;

    }
  }
#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

  n_mask_pass = colset.count();
  if(verbose && (!colset).any()) sout << "WARNING: " << (nmasks_total - n_mask_pass) << "/" << nmasks_total << " masks fail MAC filter and will be skipped...";
  //cerr << endl << Gtmp.block(0,0,10,5) << endl;

  // reset indices
  setinfo.snp_indices.resize(n_mask_pass); 
  snpinfo.resize(params->n_variants + n_mask_pass); 
  // update Gmat
  gblock->Gmat.resize(gblock->Gmat.rows(), n_mask_pass);
  vector<variant_block> tmp_snp_info; 
  snp tmpsnp;

  if(n_mask_pass == 0){
    all_snps_info = tmp_snp_info;
    return;
  }

  // store masks for testing (ignore those that failed filters)
  int k = 0;
  for(size_t i = 0; i < masks.size(); i++){
    for(int j = 0; j < n_aaf_bins; j++){

      int index_start = i * n_aaf_bins + j;

      std::ostringstream buffer;

      // mask + aaf
      buffer <<  masks[i].name << ".";
      if(j==0) buffer << "singleton";
      else if(aafs(j-1)==1) buffer << "all";
      else buffer << aafs(j-1);

      // save in snpinfo
      tmpsnp.chrom = setinfo.chrom;
      tmpsnp.ID = setinfo.ID + "." + masks[i].region_name + buffer.str();
      tmpsnp.physpos = setinfo.physpos;
      tmpsnp.allele1 = "ref";
      tmpsnp.allele2 = buffer.str();

      if(write_masks && in_bed(index_start)) {
        write_genovec(index_start);
        write_genobim(tmpsnp);
        if(write_setlist) append_setlist(index_start, tmpsnp.ID);
      }

      if(!colset(index_start)) continue;
      if(write_snplist) make_snplist(index_start, tmpsnp.ID);

      // update snpinfo
      tmpsnp.offset = params->n_variants + k; // new index in snpinfo vec.
      snpinfo[ params->n_variants + k ] = tmpsnp;

      // save index
      setinfo.snp_indices[k] = tmpsnp.offset;

      // store mask in G
      gblock->Gmat.col(k) = Gtmp.col(index_start);
      if(w_vc_tests) {
        all_snps_info[index_start].col_jmat_skat = index_start;
        all_snps_info[index_start].skip_for_vc = (j == 0) || ( aafs(j-1) != vc_aaf );
        all_snps_info[index_start].mask_name = masks[i].name;
      }
      tmp_snp_info.push_back(all_snps_info[index_start]);
      k++;

    }
  }

  // update written set list files with new set
  if(write_setlist) make_setlist(setinfo.ID, setinfo.chrom, setinfo.physpos);

  all_snps_info = tmp_snp_info;

}

void GenoMask::computeMasks_loo(const Ref<const ArrayXi>& indices_chunk, bool const& comp_full_mask, struct param* params, struct filter* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& ymat, struct geno_block* gblock, vector<variant_block> &all_snps_info, vset& setinfo, vector<snp>& snpinfo, mstream& sout){

  // check size of vblock vector
  all_snps_info.resize( nmasks_total );

  // finish building each mask 
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(int i = 0; i < nmasks_total; i++) // compute mask
    buildMask(i, setinfo.chrom, params, filters, masked_indivs, ymat, &all_snps_info[i]);
#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

  n_mask_pass = 0;
  for(int i = 0; i < nmasks_total; i++) // check if failed
    n_mask_pass += !all_snps_info[i].ignored;
  if(verbose && (n_mask_pass < nmasks_total)) sout << "WARNING: " << n_mask_pass << "/" << nmasks_total << " masks fail MAC filter...";

  // reset indices
  setinfo.snp_indices.resize(n_mask_pass); 
  snpinfo.resize(params->n_variants + n_mask_pass); 
  // update Gmat
  gblock->Gmat.resize(Gtmp.rows(), n_mask_pass);
  vector<variant_block> tmp_snp_info; 
  snp tmpsnp;

  // store masks for testing (ignore those that failed filters)
  int k = 0;
  std::ostringstream buffer;
  buffer <<  masks[0].name << "." ;
  if(n_aaf_bins == 1) buffer << "singleton";
  else if(aafs(0)==1) buffer << "all";
  else buffer << aafs(0);

  for(int i = 0; i < nmasks_total; i++){

    if(all_snps_info[i].ignored) continue;

    // update snpinfo
    tmpsnp.chrom = setinfo.chrom;
    if(comp_full_mask && (i == (nmasks_total-1) )) { // full mask
      tmpsnp.ID = setinfo.ID + "." + masks[0].region_name + buffer.str();
      tmpsnp.physpos = setinfo.physpos;
    } else { // loo mask
      tmpsnp.ID = setinfo.ID + "." + masks[0].region_name + buffer.str() + "_" + snpinfo[ setinfo.snp_indices[ indices_chunk(i) ] ].ID;
      tmpsnp.physpos = snpinfo[ setinfo.snp_indices[ indices_chunk(i) ] ].physpos;
    }
    tmpsnp.allele1 = "ref";
    tmpsnp.allele2 = buffer.str();
    tmpsnp.offset = params->n_variants + k; // new index in snpinfo vec.
    snpinfo[ params->n_variants + k ] = tmpsnp;

    // save index
    setinfo.snp_indices[k] = tmpsnp.offset;

    // store mask in G
    gblock->Gmat.col(k) = Gtmp.col(i);
    if(w_vc_tests) {
      all_snps_info[i].col_jmat_skat = i;
      all_snps_info[i].skip_for_vc = false;
    }
    tmp_snp_info.push_back(all_snps_info[i]);
    k++;

  }

  all_snps_info = tmp_snp_info;

}

void GenoMask::set_snp_masks(int const& start, int const& bs, vector<variant_block> const &all_snps_info, vset const& setinfo, vector<snp>& snpinfo, mstream& sout){

  uint64 res;
  int res2 = 1;
  keepmask = MatrixXb::Constant(bs, masks.size(), true);

  // go through each mask
  for(size_t i = 0; i < masks.size(); i++){

    // get snps who match with mask
    for(int j = 0; j < bs; j++){
      if(all_snps_info[j].ignored){
        keepmask(j, i) = false;
        continue;
      }

      // check if bit is set for at least one of the categories in mask
      // bitwise AND should return positive value
      res = (snpinfo[ setinfo.snp_indices[start + j] ].anno[setinfo.ID].id & masks[i].id);
      if(w_regions) res2 = (int)(snpinfo[ setinfo.snp_indices[start + j] ].anno[setinfo.ID].regionid & masks[i].region);
      keepmask(j, i) = (res > 0) && (res2 > 0);

      //cerr << snpinfo[ setinfo.snp_indices[j] ].ID << " " <<  snpinfo[ setinfo.snp_indices[j] ].anno[setinfo.ID].id << "\t" << masks[i].id << endl;
    }

  // cerr << i+1<<"npass="<< keepmask.col(i).count()<<"\t";
  }

}

void GenoMask::set_snp_aafs(int const& start, int const& bs, const bool& aaf_given, vector<variant_block> const &all_snps_info, vset& setinfo, vector<snp> const& snpinfo, mstream& sout){

  double upper;
  ArrayXb colkeep = ArrayXb::Constant( bs, true );// these will be nested
  keepaaf = MatrixXb::Constant(bs, n_aaf_bins, true);

  // go through each aaf cutoff (also includes singletons)
  for(int i = (n_aaf_bins-1); i >= 0; i--){
    if(i>0) upper = aafs(i-1);

    // get snps who match with mask
    for(int j = 0; j < bs; j++){
      if(all_snps_info[j].ignored || !colkeep(j) ){
        colkeep(j) = false;
        continue;
      }
       
      if( i == 0 ) colkeep(j) = force_singleton ? snpinfo[ setinfo.snp_indices[start + j] ].force_singleton : all_snps_info[j].singleton;
      else if(aaf_given) colkeep(j) = (snpinfo[ setinfo.snp_indices[start + j] ].aaf <= upper);
      else colkeep(j) = (all_snps_info[j].af1 <= upper);

      if(w_vc_tests) setinfo.ultra_rare_ind(start + j) = all_snps_info[j].mac1 <= vc_collapse_MAC;
      //cerr << snpinfo[ setinfo.snp_indices[start + j] ].aaf  << " " << all_snps_info[j].af1 << endl;
      //if(i==0 && all_snps_info[j].singleton) cerr << snpinfo[ setinfo.snp_indices[start + j] ].ID << endl;
    }

    keepaaf.col(i) = colkeep.matrix();
    //cerr << i+1 << "/" << n_aaf_bins << ": " <<  upper << "--"<< keepaaf.col(i).count()<< endl;

  }

}

bool GenoMask::check_in_lovo_mask(const Ref<const ArrayXd>& Geno, struct filter const& filters, string const& setinfo_ID, snp& snp_info, bool& is_ur, bool& to_flip, double& maf, int const& chrom, struct param const* params){

  bool is_singleton = false, in_lovo_set = false;
  bool non_par = in_non_par(chrom, snp_info.physpos, params);
  int nmales = 0, lval;
  double total, aaf, mac = 0, n_nonmiss = 0, mval;

  if(((Geno < -3) || (Geno > 2)).any())
    throw "out of bounds genotype value for variant '" + snp_info.ID + "'";

  ArrayXb keep_index = filters.ind_in_analysis && (Geno != -3.0);
  total = keep_index.select(Geno,0).sum();
  n_nonmiss = keep_index.count();

  // for MAC, check sex for non-PAR chrX
  if(non_par) {
    for (int i = 0, index = 0; i < filters.ind_ignore.size(); i++) {
      // skip samples that were ignored from the analysis
      if( filters.ind_ignore(i) ) continue;
      if( keep_index(index) ){
        // compute MAC using 0.5*g for males for variants on sex chr non-PAR (males coded as diploid) - sex is 1 for males and 0 o.w.
        lval = (params->sex(i) == 1);
        mval = Geno(index) * 0.5 * (2 - lval);
        // check if not 0/2
        if( !params->dosage_mode && (lval == 1) && (Geno(index) == 1) )
          cerr << "WARNING: genotype is 1 for a male on chrX at " << snp_info.ID << " (males should coded as diploid).";
        mac += mval;
        nmales += lval;
      }
      index++;
    }
  }

  // compute AAF and AAC
  aaf = total / n_nonmiss / 2.0;
  maf = min(aaf, 1 - aaf);
  if(!non_par) mac = total;

  // check if singleton
  if(!params->singleton_carriers) is_singleton = ( ((int)(mac+0.5)) == 1 );
  else is_singleton = (filters.ind_in_analysis && (Geno > 0)).count() == 1;

  // compute MAC (nmales=0 in auto/par)
  mac = min(mac, 2 * n_nonmiss - nmales - mac);

  if(mac < params->min_MAC)
    return false;

  // check if bit is set for at least one of the categories in mask
  // bitwise AND should return positive value
  uint64 res = snp_info.anno[setinfo_ID].id & masks[0].id;
  if( res == 0 ) return false;

  // check aaf/singleton
  if( n_aaf_bins == 1 ) in_lovo_set = force_singleton ? snp_info.force_singleton : is_singleton;
  else if(params->set_aaf) in_lovo_set = (snp_info.aaf <= aafs(0));
  else in_lovo_set = aaf <= aafs(0);

  if(!in_lovo_set) return false;

  if(w_vc_tests) {
    is_ur = mac <= vc_collapse_MAC;
    if(aaf > 0.5) to_flip = true;
  }
  //cerr << snp_info.aaf  << " " << aaf << endl;
  //if(i==0 && is_singleton) cerr << "singleton:" << snp_info.ID << "\n";

  return true;

}


void GenoMask::buildMask(int const& isnp, int const& chrom, struct param const* params, struct filter const* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& ymat, variant_block* snp_data){

  int hc_val, lval, nmales = 0;
  double ds, total = 0, mac = 0, mval, sum_pos;

  MapArXd maskvec (Gtmp.col(isnp).data(), params->n_samples, 1);
  // reset variant info
  prep_snp_stats(snp_data, params);

  // if comphet rule, threshold to 2
  if(take_comphet) maskvec = maskvec.min(2);

  // if dosages were given and writing to PLINK bed, convert dosages to hardcalls
  if(params->dosage_mode && write_masks) maskvec = maskvec.round();

  // get counts
  for (int i = 0, index = 0; i < filters->ind_ignore.size(); i++) {

    // skip samples that were ignored from the analysis
    if( filters->ind_ignore(i) ) continue;
    ds = 0;

    if( filters->ind_in_analysis(index) ){

      ds = maskvec(index);
      // distinguish missing from 0 for sum rule
      if(!w_loo && !take_max && (ds < 0) && non_missing(index,isnp)) 
        ds = 0;

      if( ds >= 0 ){
        lval = 0, mval = ds;
        if(params->test_mode && (chrom == params->nChrom)) {
          lval = (params->sex(i) == 1);
          mval = ds * 0.5 * (2 - lval);
        }
        total += ds;
        mac += mval;
        nmales += lval;
        snp_data->ns1++;

        // counts by trait
        if(filters->has_missing(index)) update_trait_counts(index, ds, mval, lval, 0, snp_data, masked_indivs);

        // get genotype counts (convert to hardcall)
        if( params->htp_out && (take_max || take_comphet) ) {
          if(params->test_mode && (chrom == params->nChrom) && (lval>0)) 
            hc_val = (ds < 1 ? 0 : 2);
          else
            hc_val = (int) (ds + 0.5); // round to nearest integer (0/1/2)
          update_genocounts(params->trait_mode==1, index, hc_val, snp_data->genocounts, masked_indivs, ymat);
        } else if( params->af_cc )
            update_af_cc(index, ds, snp_data, masked_indivs, ymat);

      }
    }

    // force masked entries to be 0
    maskvec(index++) = ds;
  }
  //cerr << maskvec.matrix().transpose().array().head(5) << endl << endl << maskvec.mean()<<endl;
  if(write_masks) make_genovec(isnp, maskvec, filters);

  // check MAC
  if(chrom != params->nChrom) mac = total; // use MAC assuming diploid coding
  // get counts by trait 
  snp_data->mac += mac; // aac
  snp_data->ns += snp_data->ns1; // ns

  // only do this when masks is in [0,2]
  if(take_max || take_comphet){
    // get counts by trait 
    snp_data->nmales += nmales; // nmales

    if(chrom != params->nChrom) {
      mac = min( mac, 2 * snp_data->ns1 - mac );
      snp_data->mac = snp_data->mac.min( 2 * snp_data->ns.cast<double>() - snp_data->mac );
    } else {
      mac = min(mac, 2 * snp_data->ns1 - nmales - mac); // males are 0/1
      snp_data->mac = snp_data->mac.min( 2 * snp_data->ns.cast<double>() - snp_data->nmales.cast<double>() - snp_data->mac );
    }

    if(mac < params->min_MAC_mask) { // don't do this with sum mask
      snp_data->ignored = true; return;
    }
  }
  snp_data->ignored_trait = snp_data->mac < params->min_MAC_mask;

  // get counts by trait 
  snp_data->af += total;

  if(params->af_cc){
    snp_data->af_control = snp_data->af - snp_data->af_case;
    snp_data->af_case /= 2 * snp_data->ns_case.cast<double>();
    snp_data->af_control /= 2 * (snp_data->ns - snp_data->ns_case).cast<double>();
  }

  total /= snp_data->ns1;
  snp_data->af1 = total / 2; // all traits
  snp_data->af /= 2 * snp_data->ns.cast<double>(); // single trait

  if(!take_max && !take_comphet) {
    snp_data->af1 /= nsites(isnp); // take average AAF across sites for sum rule
    snp_data->af /= nsites(isnp); 
    if(params->af_cc){
      snp_data->af_case /= nsites(isnp);;
      snp_data->af_control /= nsites(isnp);;
    }
  }

  if(params->use_SPA) {
    // switch to minor allele
    snp_data->flipped = ((!take_max && !take_comphet) || (params->test_type > 0)) ? false : (total > 1); // skip for DOM/REC test

    if(snp_data->flipped){
      maskvec = ( maskvec >= 0 ).select( 2 -  maskvec, maskvec);
      total = 2 - total;
    }
  }

  // apply dominant/recessive encoding & recompute mean
  if(params->test_type > 0){
    // convert to hard call if it is in dosage form
    if(params->dosage_mode) maskvec = maskvec.round();

    if(params->test_type == 1){ //dominant
      maskvec = (maskvec == 2).select(1, maskvec);
    } else if(params->test_type == 2){ //recessive
      maskvec = (maskvec >= 1).select(maskvec - 1, maskvec);
    }

    sum_pos = ((maskvec >= 0) && filters->ind_in_analysis).select(maskvec, 0).sum();
    if((params->test_type == 2) && (sum_pos < params->minHOMs)) { // filter on homALT carriers
      snp_data->ignored = true;
      return;
    }

    total = sum_pos / snp_data->ns1;
    if(total < params->numtol) {
      snp_data->ignored = true;
      return;
    }
  }

  // impute missing
  mean_impute_g(total, maskvec, filters->ind_in_analysis);

}



// compute MAF from AAF
void GenoMask::get_mafs(int const& bs, ArrayXd& mafvec, vector<variant_block> const &all_snps_info){

  for(int j = 0; j < bs; j++){
    mafvec(j) = min( all_snps_info[j].af1, 1 - all_snps_info[j].af1 );
  }

}


void GenoMask::write_info(struct param* params, struct filter const* filters, mstream& sout){

  // write fam file
  write_famfile(params, filters, sout);

  // prepare ofstream for bim file
  string fname = gfile_prefix + ".bim";
  openStream(&outfile_bim, fname, std::ios::out, sout);

  // write magic number to bed file
  uchar header[3] = {0x6c, 0x1b, 0x01};
  fname = gfile_prefix + ".bed";
  openStream(&outfile_bed, fname, std::ios::out | std::ios::binary, sout);
  outfile_bed.write( reinterpret_cast<char*> (&header[0]), sizeof(uchar) * 3);

  // number of bytes [=ceil(N/4.0)]
  gblock_size = (filters->ind_in_analysis.count() + 3) >> 2;
  reset_gvec();

  // track number of bits empty 0/2/4/6
  int nbits_left = 2 * (( gblock_size * 4 ) - filters->ind_in_analysis.count());
  // set last bits to 0 (use this uchar and apply '&' to last byte)
  last_byte_correction_factor = (1 << (8 - nbits_left)) - 1; 

}

// write to fam
void GenoMask::write_famfile(struct param* params, struct filter const* filters, mstream& sout){

  const string fname = gfile_prefix + ".fam";
  Files out;
  out.openForWrite(fname, sout);

  // columns: FID IID FA MO SEX
  for (int i = 0, index = 0; i < filters->ind_ignore.size(); i++) {
    if( filters->ind_ignore(i) ) continue;

    if( filters->ind_in_analysis(index) ){
      out << 
        params->FIDvec[index][0] << "\t" <<
        params->FIDvec[index][1] << "\t" <<
        "0\t0\t" << params->sex(i) << "\t-9\n";
    }

    index++;
  }

  out.closeFile();

  params->FIDvec.clear();
}

void GenoMask::reset_gvec(){
  gvec.resize(nmasks_total);
  for(int i = 0; i < nmasks_total; i++) 
    gvec[i].resize(gblock_size);
}

// convert to bits
void GenoMask::make_genovec(int const& isnp, Ref<const ArrayXd> mask, struct filter const* filters){

  int byte, bit_start, hc;
  setAllBitsOne(isnp);

  for(int i = 0, index = 0; i < mask.size(); i++){

    if( !filters->ind_in_analysis(i) ) continue;

    // round to nearest int
    hc = (int) (mask(i) + 0.5); 

    // using 'ref-last':
    //  00 -> hom. alt
    //  10 -> missing
    //  01 -> het
    //  11 -> hom. ref
    //  
    //  so ignore mask=0 since gvec is initialized to 11 for everyone
    if(hc == 0) {
      index++;
      continue;
    }
    byte = index >> 2;
    bit_start = (index & 3) <<1; 
    set_gvalue(isnp, byte, bit_start, hc);
    index++;
  }

  // set trailing bits to 0
  gvec[isnp][gblock_size-1] &= last_byte_correction_factor;

}

void GenoMask::setAllBitsZero(int const& isnp){
  std::fill(gvec[isnp].begin(), gvec[isnp].end(), 0ULL);
}
void GenoMask::setAllBitsOne(int const& isnp){
  std::fill(gvec[isnp].begin(), gvec[isnp].end(), ~0ULL);
}
void GenoMask::set_gvalue(int const& isnp, int const& byte, int const& bit_start, int const& val){
  // initial value is : 11
  if(val < 0) BIT_UNSET(gvec[isnp][byte], bit_start + 1);  // set to 10
  else if(val == 1) BIT_UNSET(gvec[isnp][byte], bit_start); // set to 01 
  else if(val == 2) gvec[isnp][byte] &= ~(3<<bit_start); // set to 00
}
void GenoMask::write_genovec(int const& isnp){

  outfile_bed.write( reinterpret_cast<char*> (&gvec[isnp][0]), gblock_size);

}

// get list of indices for each mask (across all AAF bins)
void GenoMask::build_map(map<string,vector<int>>& mask_map){

  for(size_t i = 0; i < masks.size(); i++){
    vector<int> myints;
    for(int j = 0; j < n_aaf_bins; j++){
      // collect indices
      int index_start = i * n_aaf_bins + j;
      myints.push_back(index_start);
    }
    // insert in map
    mask_map[ masks[i].name ] = myints;
  }

}

std::string GenoMask::build_header(){

  std::ostringstream buffer;
  size_t const nmask = mask_out.size();

  // header = ##MASKS=<Mask1="X,X";Mask2="X,X";...;MaskK="X,X">
  buffer << "##MASKS=<";
  for(size_t i = 0; i < nmask; i++)
    buffer << mask_out[i][0] << "=\"" << mask_out[i][1] << "\"" << ((i+1) < nmask ? ";" : "");

  buffer << ">\n";

  return buffer.str();
}

// prep to write list for variants in each set
void GenoMask::prep_snplist(const std::string& prefix, mstream& sout){

  string outfile = prefix + "_masks.snplist";
  sout << " * writing list of variants for each mask in file [" << outfile << "]\n";
  snplist_out.openForWrite(outfile, sout);
  list_snps.resize(nmasks_total);
}

// prep to write set list files
void GenoMask::prep_setlists(const std::string& fin, const std::string& prefix, mstream& sout){

  int nfiles = 0, lineread = 0;
  string line;
  std::vector< string > tmp_str_vec, suffix;
  Files myfile;

  myfile.openForRead(fin, sout);
  sout << "   +writing new set list files using [" << fin << "] ";

  setfiles_index.resize(nmasks_total);
  map<string,vector<int>> mask_map;
  build_map(mask_map);

  while( myfile.readLine(line) ){

    lineread++;
    tmp_str_vec = string_split(line,"\t ,");
    // file suffix + list of masks to include
    if( tmp_str_vec.size() < 2 )
      throw "line " + to_string( lineread ) + " has too few entries." ;

    // get index of masks
    vector<int> mindices;
    for(size_t i = 1; i < tmp_str_vec.size(); i++){
      if (!in_map(tmp_str_vec[i], mask_map)) continue;
      mindices.insert(mindices.end(), mask_map[ tmp_str_vec[i] ].begin(), mask_map[ tmp_str_vec[i] ].end());
    }
    // sort and remove duplicates
    std::sort(mindices.begin(), mindices.end());
    mindices.erase( unique( mindices.begin(), mindices.end() ), mindices.end() );

    // check at least one mask
    if(mindices.size() == 0) continue;

    suffix.push_back(tmp_str_vec[0]);
    for(size_t i = 0; i < mindices.size(); i++)
      setfiles_index[ mindices[i] ].push_back(nfiles);
    nfiles++;

    //cerr << suffix.back() << " -> " << mindices.size() << endl;
  }

  if(nfiles < 1) 
    throw "all set list files have unknown masks.";

  sout << " n_files = " << nfiles << endl;
  write_setlist = true;

  // open file for writing
  setfiles.resize(nfiles);
  for(int i = 0; i < nfiles; i++) {
    line = prefix + "_" + suffix[i] + ".setlist";
    setfiles[i] = std::make_shared<Files>();
    setfiles[i]->openForWrite( line, sout );
  }
  list_masks.resize(nfiles);

}


void GenoMask::write_genobim(struct snp const& tsnp){

  // write mask info to bim file using ref-last
  // CHR ID 0 BP ALT REF 
  outfile_bim << tsnp.chrom << "\t" << tsnp.ID << "\t0\t" << tsnp.physpos << "\t" << tsnp.allele2 << "\t" << tsnp.allele1 << endl;

}

void GenoMask::append_snplist(int const& imask, ArrayXb const& colkeep, int const& start, vset const& setinfo, vector<snp> const& snpinfo){

  // add snps
  if( colkeep.count() == 0 ) return;

  for(int k = 0; k < colkeep.size(); k++){
    if(!colkeep(k)) continue;
    list_snps[imask].push_back( snpinfo[ setinfo.snp_indices[start + k] ].ID );
  }
}

void GenoMask::make_snplist(int const& imask, string const& mask_name){
  // add snplist
  if( list_snps[imask].size() > 0 )
    snplist_out << mask_name << "\t" << print_csv( list_snps[imask] ) << endl;
}


void GenoMask::append_setlist(int const& imask, string const& mname){
  // add mask name
  for(size_t i = 0; i < setfiles_index[imask].size(); i++)
    list_masks[ setfiles_index[imask][i] ].push_back(mname);
}

void GenoMask::make_setlist(string const& sname, int const& chr, uint32_t const& pos){
  for(size_t i = 0; i < setfiles.size(); i++)
    if( list_masks[i].size() > 0 )// add set name + masks
      (*setfiles[i]) << sname << " " << chr << " " << pos << " " << print_csv( list_masks[i] ) << endl;
}

void GenoMask::closeFiles(){
  outfile_bim.close();
  outfile_bed.close();
  if(write_setlist){
    for(size_t i = 0; i < setfiles.size(); i++) 
      setfiles[i]->closeFile();
  }
  if(write_snplist) snplist_out.closeFile();
}

// for LOVO
ArrayXi check_lovo_snplist(const Ref<const ArrayXi>& indices, vector<uint64> const& offsets, vector<snp> const& snpinfo, string const& masks_loo_snpfile){

  int bs = indices.size();
  ArrayXi lovo_masks;

  if(masks_loo_snpfile == ""){
    lovo_masks = ArrayXi::LinSpaced(bs, 0, bs-1);;
    return lovo_masks;
  }

  string line;
  map<string, bool> comp_lovo_mask;
  vector<int> lovo_masks_vec;
  ifstream infile;

  // get list of all variants for which to compute lovo mask
  infile.open(masks_loo_snpfile);
  while(getline(infile, line))
    comp_lovo_mask[ line ] = true;
  infile.close();

  if(comp_lovo_mask.size() == 0)
    throw "no variants were specified in the '--lovo-snplist' file.";

  // check variants which are in map
  for(int i = 0; i < indices.size(); i++)
    if(in_map(snpinfo[ offsets[indices(i)] ].ID, comp_lovo_mask))
      lovo_masks_vec.push_back( i );

  if(lovo_masks_vec.size() == 0)
    throw "none of the genotyped variants are present in the '--lovo-snplist' file.";

  lovo_masks = ArrayXi::Map(lovo_masks_vec.data(), lovo_masks_vec.size());
  return lovo_masks;
}


// make seq(0,n-1) removing i-th entry
ArrayXi get_index_vec_loo(int const& i, int const& n){

  ArrayXi iseq ( n -1 );
  for(int j = 0, k = 0; j < n; j++){
    if(j != i) iseq(k++) = j;
  }

  return iseq;
}
