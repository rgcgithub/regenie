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


void GenoMask::updateMasks(int const& start, int const& bs, struct param* params, struct filter* filters, const Ref<const MatrixXb>& masked_indivs, struct geno_block* gblock, vector<variant_block> &all_snps_info, vset& setinfo, vector<snp>& snpinfo, mstream& sout){

  // identify which snps are in each mask
  set_snp_masks(start, bs, all_snps_info, setinfo, snpinfo, sout);
  // identify which snps are in each aaf bin
  set_snp_aafs(start, bs, params->set_aaf, all_snps_info, setinfo, snpinfo, sout);

  MatrixXb Jmat;
  MatrixXd rare_mask_tmp;
  if(w_vc_tests) {
    Jmat = MatrixXb::Constant(bs, nmasks_total, false);
    if(setinfo.ultra_rare_ind.any()) rare_mask_tmp = setinfo.vc_rare_mask; // not safe to update SpMat in parallel (not many columns)
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

        SpMat gv, mv;
        mv = maskvec.matrix().sparseView();
        for(int k = 0; k < colkeep.size(); k++){
          if(!colkeep(k)) continue;
          gv = gblock->Gmat.col(k).sparseView();
          mv = gv.cwiseMax(mv);
        }
        maskvec = MatrixXd(mv).array();

      } else {

        int l;
        double ds;
        SpVec gv;

        for(int k = 0; k < colkeep.size(); k++){
          if(!colkeep(k)) continue;
          gv = gblock->Gmat.col(k).sparseView();

          // sum rule (ignore -3)
          for (SparseVector<double>::InnerIterator it(gv); it; ++it) {
            l = it.index();
            ds = it.value();

            if( !filters->ind_in_analysis(l) || (ds == -3)) continue;

            if( maskvec(l) == -3 ) maskvec(l) = ds;
            else maskvec(l) += ds;
          }

          // for genotype counts, identify when (-3) is 0
          non_missing.col(index_start).array() = non_missing.col(index_start).array() || ( gblock->Gmat.col(k).array() >= 0 );
        }
       //if(i==2 && j==1) cout << (maskvec == -3).select(0,maskvec).sum() << endl;

      }

      // get ultra-rare mask if using VC test (take max)
      if(w_vc_tests && setinfo.ultra_rare_ind.segment(start, bs).any() && ((j == 0) || ( aafs(j-1) <= vc_aaf )) ) {
        SpMat gv, mv;
        mv = setinfo.vc_rare_mask.col(index_start);
        for(int k = 0; k < colkeep.size(); k++){
          if(!colkeep(k) || !setinfo.ultra_rare_ind(start+k)) continue;
          MapArXd garr (gblock->Gmat.col(k).data(), params->n_samples, 1);
          // flip if necessary
          if(all_snps_info[start+k].af1 > 0.5) gv = (garr == -3).select(0, 2 - garr).matrix().sparseView();
          else gv = garr.matrix().sparseView();
          mv = gv.cwiseMax(mv);
          setinfo.vc_rare_mask_non_missing.col(index_start).array() = setinfo.vc_rare_mask_non_missing.col(index_start).array() || (garr != -3);
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


// should only be called once
void GenoMask::updateMasks_loo(int const& start, int const& bs, struct param const* params, struct filter const* filters, const Ref<const MatrixXb>& masked_indivs, struct geno_block* gblock, vector<variant_block> &all_snps_info, vset& setinfo, vector<snp>& snpinfo, mstream& sout){

  // identify which snps are in each mask
  set_snp_masks(start, bs, all_snps_info, setinfo, snpinfo, sout);
  // identify which snps are in each aaf bin
  set_snp_aafs(start, bs, params->set_aaf, all_snps_info, setinfo, snpinfo, sout);

  colset = keepmask.col(0).array() && keepaaf.col(n_aaf_bins - 1).array(); // take last aaf

  int nkept = colset.count();
  if(nkept == 0) return;

  nmasks_total = nkept; // number of LOO masks
  Gtmp = MatrixXd::Constant(params->n_samples, nkept + 1, -3); // add full mask
  if(!take_max && !take_comphet) nsites = nkept - 1; // each loo mask has (n-1) sites included for AAF calculation

  MatrixXb Jmat;
  vector<SpMat> rare_mask_tmp;
  if(w_vc_tests) {
    Jmat = MatrixXb::Constant(bs, nkept+1, false);
    setinfo.Jmat = MatrixXb::Constant(bs+nkept+1, nkept+1, false); // last K rows are for ultra rare masks
    if(setinfo.ultra_rare_ind.any()) {
      rare_mask_tmp.resize(Jmat.cols()); // not safe to update SpMat in parallel
      setinfo.vc_rare_mask.resize(params->n_samples, Jmat.cols());
      setinfo.vc_rare_mask_non_missing = MatrixXb::Constant(params->n_samples, Jmat.cols(), false);
    }
  }

  // update mask using LOO
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(int i = 0; i < colset.size(); i++){
    if(!colset(i)) continue;

    ArrayXb colkeep_loo = colset;
    colkeep_loo(i) = false; // mask snp

    bool has_non_missing;
    double ds;
    int ix = colset.head(i).count(); // new index among unmasked snps
    MapArXd maskvec (Gtmp.col(ix).data(), Gtmp.rows(), 1);

    for(int k = 0; k < Gtmp.rows(); k++){
      if( !filters->ind_in_analysis(k) ) continue;

      if(take_max){ // max rule to combine variants across sites
        ds = max( maskvec(k), (colkeep_loo).select(gblock->Gmat.row(k).transpose().array(),-3).maxCoeff());
      } else {

        // sum rule (ignore missing)
        ds = 0;
        for(int l = 0; l < colkeep_loo.size(); l++)
          if(colkeep_loo(l) && (gblock->Gmat(k,l) != -3)) {
            has_non_missing = true;
            ds += gblock->Gmat(k,l);
          }

        if(maskvec(k) != -3) ds += maskvec(k);
        else if( (ds == 0) && !has_non_missing ) ds = -3;

      }

      maskvec(k) = ds;

    }

    if(w_vc_tests) // track variants in mask
      Jmat.col(ix) = colkeep_loo.matrix();

  }

#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif


  // for vc tests
  if(w_vc_tests && setinfo.ultra_rare_ind.any()) {
    // flip to minor
    for(int i = 0; i < colset.size(); i++){
      if(!colset(i) || !setinfo.ultra_rare_ind(i) ||
          (all_snps_info[start+i].af1 <= 0.5) ) continue;
      MapArXd garr (gblock->Gmat.col(i).data(), params->n_samples, 1);
      garr = (garr == -3).select(-3, 2 - garr);
    }

    // get ultra-rare masks
#if defined(_OPENMP)
    setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
    for(int i = 0; i < colset.size(); i++){
      if(!colset(i)) continue;
      ArrayXb colkeep_loo = colset;
      colkeep_loo(i) = false; // mask snp
      ArrayXd garr = ArrayXd::Zero(params->n_samples,1);

      double ds;
      int ix = colset.head(i).count(); // new index among unmasked snps

      for(int k = 0; k < Gtmp.rows(); k++){
        if( !filters->ind_in_analysis(k) ) continue;
        ds = (colkeep_loo && setinfo.ultra_rare_ind).select(gblock->Gmat.row(k).transpose().array(),-3).maxCoeff();
        setinfo.vc_rare_mask_non_missing(k, ix) = (ds>=0);
        if(ds > 0) garr(k) = ds;
      }
      if((garr>0).any())
        rare_mask_tmp[ix] = garr.matrix().sparseView();
    }
#if defined(_OPENMP)
    setNbThreads(params->threads);
#endif

    int jx; // new index among unmasked snps
    for(int i = 0; i < colset.size(); i++){
      jx = colset.head(i).count(); // new index among unmasked snps
      if(colset(i) && setinfo.vc_rare_mask_non_missing.col(jx).any()){
        setinfo.vc_rare_mask.col(jx) = rare_mask_tmp[jx];
      }
    }

  }

  MatrixXd tmp_v;
  // compute full mask
  for(int i = 0; i < colset.size(); i++){
    if(!colset(i)) continue;

    bool has_non_missing;
    double ds;
    int ix = colset.head(i).count(); // new index among unmasked snps
    MapArXd maskvec (Gtmp.rightCols(1).data(), Gtmp.rows(), 1); // in last column
    maskvec = Gtmp.col(ix); // start from LOO mask of 1st unmasked site
    if(w_vc_tests && setinfo.ultra_rare_ind.any()){ // for ur mask in vc tests
      tmp_v = setinfo.vc_rare_mask.col(ix);
      setinfo.vc_rare_mask_non_missing.rightCols(1) = setinfo.vc_rare_mask_non_missing.col(ix);
    }

    for(int k = 0; k < Gtmp.rows(); k++){
      if( !filters->ind_in_analysis(k) ) continue;

      if(take_max) ds = max( maskvec(k), gblock->Gmat(k, i));
      else {
        // sum rule (ignore missing)
        ds = 0;
        if(gblock->Gmat(k,i) != -3) {
          has_non_missing = true;
          ds += gblock->Gmat(k,i);
        }

        if(maskvec(k) != -3) ds += maskvec(k);
        else if( (ds == 0) && !has_non_missing ) ds = -3;
      }
      maskvec(k) = ds;

      // get ultra-rare mask if using VC test (take max)
      if(w_vc_tests && setinfo.ultra_rare_ind(i)){
        ds = gblock->Gmat(k, i);
        setinfo.vc_rare_mask_non_missing.rightCols(1)(k) = setinfo.vc_rare_mask_non_missing(k, ix) || (ds>=0);
        if(all_snps_info[i].af1 > 0.5) ds = ds == -3 ? 0 : 2 - ds;
        tmp_v(k, 0) = max(ds, tmp_v(k, 0));
      }
    }
    // cerr << endl << Gtmp.col(ix).head(3) << "\n\n" << gblock->Gmat.col(i).head(3) << "\n\n" << maskvec.head(3) << endl;

    if(w_vc_tests){ // track variants in mask
      Jmat.rightCols(1) = Jmat.col(ix);
      Jmat.rightCols(1)(i) = true; 
      if(setinfo.ultra_rare_ind.any()) 
        setinfo.vc_rare_mask.rightCols(1) = tmp_v.sparseView();
    }

    break;
  }

  if(!take_max && !take_comphet) {
    nsites.conservativeResize( nkept + 1);
    nsites(nkept) = nkept; // to compute AAF with sum for full mask
  }
  if(w_vc_tests) 
    setinfo.Jmat.topRows(Jmat.rows()) = Jmat;

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
    bool column_set = (Gtmp.col(index_mask).array() != -3).count() > 0;
    colset(index_mask) = column_set;

    // don't parallelize inner loop (cumulative updates)
    for(int j = 1; j < n_aaf_bins; j++){

      int index_start = index_mask + j;
      // check if there are variants in mask
      if(colset(index_start-1) || ((Gtmp.col(index_start).array() != -3).count() > 0)) 
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

          if( !filters->ind_in_analysis(l) || (ds == -3)) continue;

          if( maskvec(l) == -3 ) maskvec(l) = ds;
          else maskvec(l) += ds;

        }

        // for genotype counts, identify when (-3) is 0
        non_missing.col(index_start).array() = non_missing.col(index_start).array() || non_missing.col(index_start-1).array();

      }

      // update ultra-rare mask
      if( w_vc_tests && ( aafs(j-1) <= vc_aaf ) && rare_mask_tmp.col(index_start-1).nonZeros() ){
        rare_mask_tmp.col(index_start) = rare_mask_tmp.col(index_start).cwiseMax(rare_mask_tmp.col(index_start-1));
        vc_rare_mask_non_missing.col(index_start).array() = vc_rare_mask_non_missing.col(index_start).array() || vc_rare_mask_non_missing.col(index_start-1).array();
        if( aafs(j-1) < vc_aaf ){ // remove data as not needed anymore
          rare_mask_tmp.col(index_start-1) *= 0;
          vc_rare_mask_non_missing.col(index_start-1).array() = false;
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

void GenoMask::computeMasks_loo(struct param* params, struct filter* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& ymat, struct geno_block* gblock, vector<variant_block> &all_snps_info, vset& setinfo, vector<snp>& snpinfo, mstream& sout){

  // check size of vblock vector
  all_snps_info.resize( nmasks_total + 1 ); // add full mask

  ArrayXb in_bed = colset;
  vector<uint64> old_indices = setinfo.snp_indices;

  // finish building each mask 
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(int i = 0; i <= colset.size(); i++){

    int index_start;

    if( (i < colset.size()) && !colset(i) ) continue;
    index_start = in_bed.head(i).count();

    // compute mask
    buildMask(index_start, setinfo.chrom, params, filters, masked_indivs, ymat, &all_snps_info[index_start]);

    if(i < colset.size()) colset(i) = !all_snps_info[index_start].ignored;

  }
#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

  n_mask_pass = colset.count();
  if(verbose && ((!colset).count() > 0)) sout << "WARNING: " << (!colset).count() << "/" << nmasks_total << " masks fail MAC filter...";
  //cerr << "Npass="<< n_mask_pass << "/" << nmasks_total <<endl;

  // reset indices (add full mask)
  setinfo.snp_indices.resize(n_mask_pass+1); 
  snpinfo.resize(params->n_variants + n_mask_pass+1); 
  // update Gmat
  gblock->Gmat.resize(gblock->Gmat.rows(), n_mask_pass+1);
  vector<variant_block> tmp_snp_info; 
  snp tmpsnp;

  // store masks for testing (ignore those that failed filters)
  int k = 0;
  for(int i = 0; i <= colset.size(); i++){

    int index_start = in_bed.head(i).count();

    if(i < colset.size() && !colset(i)) 
      continue;
    else if(i == colset.size() && all_snps_info[index_start].ignored)
      continue;

    std::ostringstream buffer;
    buffer <<  masks[0].name << "." ;
    if(n_aaf_bins == 1) buffer << "singleton";
    else if(aafs(0)==1) buffer << "all";
    else buffer << aafs(0);

    // update snpinfo
    tmpsnp.chrom = setinfo.chrom;
    if(i < colset.size()) { // loo mask
      tmpsnp.ID = setinfo.ID + "." + masks[0].region_name + buffer.str() + "_" + snpinfo[ old_indices[i] ].ID;
      tmpsnp.physpos = snpinfo[ old_indices[i] ].physpos;
    } else { // full mask
      tmpsnp.ID = setinfo.ID + "." + masks[0].region_name + buffer.str();
      tmpsnp.physpos = setinfo.physpos;
    }
    tmpsnp.allele1 = "ref";
    tmpsnp.allele2 = buffer.str();
    tmpsnp.offset = params->n_variants + k; // new index in snpinfo vec.
    snpinfo[ params->n_variants + k ] = tmpsnp;

    // save index
    setinfo.snp_indices[k] = tmpsnp.offset;

    // store mask in G
    gblock->Gmat.col(k) = Gtmp.col(index_start);
    if(w_vc_tests) {
      all_snps_info[index_start].col_jmat_skat = index_start;
      all_snps_info[index_start].skip_for_vc = false;
    }
    tmp_snp_info.push_back(all_snps_info[index_start]);
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
      if(!w_loo && !take_max && (ds == -3) && non_missing(index,isnp)) 
        ds = 0;

      if( ds != -3 ){
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
      maskvec = ( maskvec != -3.0 ).select( 2 -  maskvec, maskvec);
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

    sum_pos = ((maskvec != -3) && filters->ind_in_analysis).select(maskvec, 0).sum();
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

