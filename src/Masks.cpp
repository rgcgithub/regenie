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
#include "Files.hpp"
#include "Masks.hpp"

#include <Eigen/SparseCore>

using namespace std;
using namespace Eigen;

GenoMask::GenoMask() { // @suppress("Class members should be properly initialized")
}

GenoMask::~GenoMask() {
  // TODO Auto-generated destructor stub
}

void GenoMask::setBins(struct param* params, mstream& sout){

  int n_vals;

  // check if not singleton when using LOO
  if(!params->mask_loo || ( params->mask_loo && (params->mbins[0] != "singleton")) ){
    if(params->mbins.size() >= 1){
      // convert them to double
      for( size_t i = 0; i < params->mbins.size(); i++)
        aafs.push_back( convertDouble( params->mbins[i], params, sout) );
      // sort and retain unique values
      std::sort(aafs.begin(), aafs.end());
      aafs.erase( unique( aafs.begin(), aafs.end() ), aafs.end() );

      n_vals = aafs.size();

      // check valid (in (0,1))
      if( std::count_if(aafs.begin(), aafs.end(), std::bind2nd(std::greater<double>(), minAAF)) != n_vals || std::count_if(aafs.begin(), aafs.end(), std::bind2nd(std::less<double>(), 1-minAAF)) != n_vals ){
        sout << "ERROR : You must specify values for --aaf-bins in (" << minAAF << "," << 1-minAAF << ")\n" << params->err_help;
        exit(-1);
      } else if(std::count_if(aafs.begin(), aafs.end(), std::bind2nd(std::greater<double>(), 0.5)) > 0)
        sout << "WARNING: For computational efficiency, it is recommended that AAF cutoffs <0.5\n";

      max_aaf = *max_element(aafs.begin(), aafs.end());
    } else aafs.push_back(default_aaf); // default
  }

  n_aaf_bins = aafs.size();

  if(n_aaf_bins > max_aaf_bins){
    sout << "ERROR: Number of AAF bins (=" << n_aaf_bins << ") above maximum (="<< max_aaf_bins << ")\n";
    exit(-1);
  }

  if(!params->mask_loo){
    sout << left << std::setw(20) << " * aaf cutoffs" << ": [" << n_aaf_bins << " : "; 
    for(int i = 0; i < n_aaf_bins; i++) sout << aafs[i] << " ";
    sout << "] + singletons\n";
  } else {
    sout << left << std::setw(20) << " * aaf cutoff" << ": " ;
    if(n_aaf_bins == 0) sout << "singleton";
    else sout << aafs[0];
    sout << endl;
  }

  n_aaf_bins++; // add singleton bin

  nmasks_total = n_aaf_bins * masks.size(); // total number of masks

  w_regions = params->w_regions;
  w_loo = params->mask_loo;

  if(w_regions) base_masks = masks;
}

void GenoMask::prepMasks(const int ntotal, const string& setID) {

  maskinfo tmp_region_mask;
  std::map <std::string, uchar>::iterator itr;

  // make new set of masks if using set regions
  if(w_regions){ 
    masks.resize(0);
    for(int i = 0; i < base_masks.size(); i++ ){
      tmp_region_mask = base_masks[i];
      for (itr = regions[setID].begin(); itr != regions[setID].end(); ++itr) {
        // add region info
        tmp_region_mask.region_name = itr->first + ".";
        tmp_region_mask.region = itr->second;
        masks.push_back(tmp_region_mask);
      }
      if(!w_loo){// add mask across all regions
        tmp_region_mask = base_masks[i];
        tmp_region_mask.region |= 255; //set all 8 bits
        masks.push_back(tmp_region_mask);
      }
    }
    nmasks_total = n_aaf_bins * masks.size();
  } 

  Gtmp = MatrixXd::Constant(ntotal, nmasks_total, -3);
  colset = ArrayXb::Constant( nmasks_total, false );
  if(!take_max) {
    non_missing = MatrixXb::Constant(ntotal, nmasks_total, false);
    if(!take_comphet) nsites = ArrayXi::Constant(nmasks_total, 0);
  }

  if(write_setlist) {
    for(size_t i = 0;i < list_masks.size(); i++)
      list_masks[i].resize(0);
  }

}


void GenoMask::updateMasks(const int start, const int bs, struct param* params, struct filter* filters, const Ref<const MatrixXb>& masked_indivs, struct geno_block* gblock, vector<variant_block> &all_snps_info, vset& setinfo, vector<snp>& snpinfo, mstream& sout){

  // identify which snps are in each mask
  set_snp_masks(start, bs, all_snps_info, setinfo, snpinfo, sout);
  // identify which snps are in each aaf bin
  set_snp_aafs(start, bs, params->set_aaf, all_snps_info, setinfo, snpinfo, sout);

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

      // ignore variants in previous AAF categories (accumulation is in next loop)
      if(j>0) colkeep = colkeep && !keepaaf.col(j-1).array(); 

      // if there are no variants included, continue
      if( colkeep.count() == 0 ) continue;
      //if(i==2 && j==1) cerr << i << " " << j << " " << colkeep.count() << " ";

      // update mask
      MapArXd maskvec (Gtmp.col(index_start).data(), params->n_samples, 1);

      if(take_max) {

        SparseMatrix<double> gv, mv;
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
        SparseVector<double> gv;

        for(int k = 0; k < colkeep.size(); k++){
          if(!colkeep(k)) continue;
          gv = gblock->Gmat.col(k).sparseView();

          // sum rule (ignore -3)
          for (SparseVector<double>::InnerIterator it(gv); it; ++it) {
            l = it.index();
            ds = it.value();

            if( !filters->ind_in_analysis(l) || (ds == -3)) continue;
            if( !params->strict_mode || (params->strict_mode && masked_indivs(l,0)) ){
              if( maskvec(l) == -3 ) maskvec(l) = ds;
              else maskvec(l) += ds;
            }
          }

          // for genotype counts, identify when (-3) is 0
          non_missing.col(index_start).array() = non_missing.col(index_start).array() || ( gblock->Gmat.col(k).array() >= 0 );
        }
       //if(i==2 && j==1) cout << (maskvec == -3).select(0,maskvec).sum() << endl;

      }

    }
  }

#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

}


// should only be called once
void GenoMask::updateMasks_loo(const int start, const int bs, struct param* params, struct filter* filters, const Ref<const MatrixXb>& masked_indivs, struct geno_block* gblock, vector<variant_block> &all_snps_info, vset& setinfo, vector<snp>& snpinfo, mstream& sout){

  // identify which snps are in each mask
  set_snp_masks(start, bs, all_snps_info, setinfo, snpinfo, sout);
  // identify which snps are in each aaf bin
  set_snp_aafs(start, bs, params->set_aaf, all_snps_info, setinfo, snpinfo, sout);

  colset = keepmask.col(0).array() && keepaaf.col(n_aaf_bins - 1).array(); // take last aaf

  int nkept = colset.count();
  if(nkept == 0) return;
  nmasks_total = nkept;

  Gtmp = MatrixXd::Constant(params->n_samples, nkept + 1, -3); // add full mask
  if(!take_max && !take_comphet) nsites = nkept - 1; // to compute AAF with sum

  // update mask using LOO
#if defined(_OPENMP)
  setNbThreads(1);
  // use MT 
#pragma omp parallel for schedule(dynamic)
#endif
  for(int i = 0; i < colset.size(); i++){
    if(!colset(i)) continue;

    ArrayXb colkeep_loo = colset;
    colkeep_loo(i) = false; // mask snp

    int ix = colkeep_loo.head(i+1).count();
    bool has_non_missing;
    double ds;
    MapArXd maskvec (Gtmp.col(ix).data(), params->n_samples, 1);

    for(size_t k = 0; k < params->n_samples; k++){
      if( !filters->ind_in_analysis(k) ) continue;

      if( !params->strict_mode || (params->strict_mode && masked_indivs(k,0)) ){

        if(take_max){ // max rule to combine variants across sites
          ds = max( maskvec(k), (colkeep_loo).select(gblock->Gmat.row(k).transpose().array(),-3).maxCoeff());
        } else {

          // sum rule (ignore missing)
          ds = 0;
          for(int l = 0; l < colkeep_loo.size(); l++)
            if(colkeep_loo(l) && (gblock->Gmat(k,l) != -3)) has_non_missing = true, ds += gblock->Gmat(k,l);

          if(maskvec(k) != -3) ds += maskvec(k);
          else if( ds == 0 && !has_non_missing ) ds = -3;

        }

        maskvec(k) = ds;
      }

    }
  }

#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

  // compute full mask
  for(int i = 0; i < colset.size(); i++){
    if(!colset(i)) continue;

    int ix = colset.head(i+1).count();
    bool has_non_missing;
    double ds;
    MapArXd maskvec (Gtmp.rightCols(1).data(), params->n_samples, 1); // in last column
    maskvec = Gtmp.col(ix); // start from LOO mask for 1st site

    for(size_t k = 0; k < params->n_samples; k++){
      if( !filters->ind_in_analysis(k) ) continue;

      if( !params->strict_mode || (params->strict_mode && masked_indivs(k,0)) ){

        if(take_max) ds = max( maskvec(k), gblock->Gmat(k, i));
        else {

          // sum rule (ignore missing)
          ds = 0;
          if(gblock->Gmat(k,i) != -3) has_non_missing = true, ds += gblock->Gmat(k,i);
          if(maskvec(k) != -3) ds += maskvec(k);
          else if( ds == 0 && !has_non_missing ) ds = -3;
        }
        maskvec(k) = ds;
      }
    }
   // cerr << endl << Gtmp.col(ix).head(3) << "\n\n" << gblock->Gmat.col(i).head(3) << "\n\n" << maskvec.head(3) << endl;
    break;
  }
  if(!take_max && !take_comphet) {
    nsites.conservativeResize( nkept + 1);
    nsites(nkept) = nkept; // to compute AAF with sum for full mask
  }

}


void GenoMask::tally_masks(struct param* params, struct filter* filters, const Ref<const MatrixXb>& masked_indivs){

#if defined(_OPENMP)
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

        SparseMatrix<double> gv, mv;
        // sparse of current
        mv = maskvec.matrix().sparseView();
        // sparse of previous
        gv = Gtmp.col(index_start-1).sparseView();
        // aggregate to current
        maskvec = MatrixXd(gv.cwiseMax(mv)).array();

      } else {

        SparseVector<double> gv;

        // add previous
        gv = Gtmp.col(index_start-1).sparseView();

        // sum rule (ignore -3)
        for (SparseVector<double>::InnerIterator it(gv); it; ++it) {
          int l = it.index();
          ds = it.value();

          if( !filters->ind_in_analysis(l) || (ds == -3)) continue;
          if( !params->strict_mode || (params->strict_mode && masked_indivs(l,0)) ){
            if( maskvec(l) == -3 ) maskvec(l) = ds;
            else maskvec(l) += ds;
          }
        }

        // for genotype counts, identify when (-3) is 0
        non_missing.col(index_start).array() = non_missing.col(index_start).array() || non_missing.col(index_start-1).array();

      }

    }
  }
#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif
}

void GenoMask::computeMasks(struct param* params, struct filter* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& ymat, struct geno_block* gblock, vector<variant_block> &all_snps_info, vset& setinfo, vector<snp>& snpinfo, mstream& sout){

  // check size of vblock vector
  if(((int) all_snps_info.size()) < nmasks_total) all_snps_info.resize( nmasks_total );

  // tally masks
  tally_masks(params, filters, masked_indivs);

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
  //cerr << "Npass="<< n_mask_pass << "/" << nmasks_total <<endl;
  //cerr << endl << Gtmp.block(0,0,10,5) << endl;

  // reset indices
  setinfo.snp_indices.resize(n_mask_pass); 
  snpinfo.resize(params->n_variants + n_mask_pass); 
  // update Gmat
  gblock->Gmat.resize(gblock->Gmat.rows(), n_mask_pass);
  vector<variant_block> tmp_snp_info; 
  snp tmpsnp;

  // store masks for testing (ignore those that failed filters)
  int k = 0;
  for(size_t i = 0; i < masks.size(); i++){
    for(int j = 0; j < n_aaf_bins; j++){

      int index_start = i * n_aaf_bins + j;

      std::ostringstream buffer;

      // mask + aaf
      buffer <<  masks[i].name << ".";
      if(j==0) buffer << "singleton";
      else buffer << aafs[j-1];

      // save in snpinfo
      tmpsnp.chrom = setinfo.chrom;
      tmpsnp.ID = setinfo.ID + "." + masks[i].region_name + buffer.str();
      tmpsnp.physpos = setinfo.physpos;
      tmpsnp.allele1 = "ref";
      tmpsnp.allele2 = buffer.str();

      if(params->write_masks && in_bed(index_start)) {
        write_genovec(index_start);
        write_genobim(tmpsnp);
        if(write_setlist) append_setlist(index_start, tmpsnp.ID);
      }

      if(!colset(index_start)) continue;

      // update snpinfo
      tmpsnp.offset = params->n_variants + k; // new index in snpinfo vec.
      snpinfo[ params->n_variants + k ] = tmpsnp;

      // save index
      setinfo.snp_indices[k] = tmpsnp.offset;

      // store mask in G
      gblock->Gmat.col(k) = Gtmp.col(index_start);
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
  // use MT 
#pragma omp parallel for schedule(dynamic)
#endif
  for(int i = 0; i <= colset.size(); i++){

    int index_start;

    if(i < colset.size()){
      if(!colset(i)) continue;
      index_start = colset.head(i+1).count() - 1;
    } else index_start = colset.count(); // full mask

    // compute mask
    buildMask(index_start, setinfo.chrom, params, filters, masked_indivs, ymat, &all_snps_info[index_start]);

    if(i < colset.size()) colset(i) = !all_snps_info[index_start].ignored;

  }
#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

  n_mask_pass = colset.count();
  //cerr << "Npass="<< n_mask_pass << "/" << nmasks_total <<endl;

  // reset indices
  setinfo.snp_indices.resize(n_mask_pass+1); 
  snpinfo.resize(params->n_variants + n_mask_pass+1); 
  // update Gmat
  gblock->Gmat.resize(gblock->Gmat.rows(), n_mask_pass+1);
  vector<variant_block> tmp_snp_info; 
  snp tmpsnp;

  // store masks for testing (ignore those that failed filters)
  int k = 0;
  for(int i = 0; i <= colset.size(); i++){
    if(i < colset.size() && !colset(i)) continue;

    int index_start = (i < colset.size() ? in_bed.head(i+1).count() - 1 :  in_bed.count());

    std::ostringstream buffer;
    buffer <<  masks[0].name << "." ;
    if(n_aaf_bins == 1) buffer << "singleton";
    else buffer << aafs[0];

    // update snpinfo
    tmpsnp.chrom = setinfo.chrom;
    if(i < colset.size()) { // lovo mask
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
    tmp_snp_info.push_back(all_snps_info[index_start]);
    k++;

  }

  all_snps_info = tmp_snp_info;

}

void GenoMask::set_snp_masks(const int start, const int bs, vector<variant_block> &all_snps_info, vset& setinfo, vector<snp>& snpinfo, mstream& sout){

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

void GenoMask::set_snp_aafs(const int start, const int bs, const bool aaf_given, vector<variant_block> &all_snps_info, vset& setinfo, vector<snp>& snpinfo, mstream& sout){

  double upper;
  ArrayXb colkeep = ArrayXb::Constant( bs, true );// these will be nested
  keepaaf = MatrixXb::Constant(bs, n_aaf_bins, true);

  // go through each aaf cutoff (also includes singletons)
  for(int i = (n_aaf_bins-1); i >= 0; i--){
    if(i>0) upper = aafs[i-1];

    // get snps who match with mask
    for(int j = 0; j < bs; j++){
      if(all_snps_info[j].ignored || !colkeep(j) ){
        colkeep(j) = false;
        continue;
      }
       
      if( i == 0 ) colkeep(j) = all_snps_info[j].singleton;
      else if(aaf_given) colkeep(j) = (snpinfo[ setinfo.snp_indices[start + j] ].aaf <= upper);
      else colkeep(j) = (all_snps_info[j].af1 <= upper);
      //cerr << snpinfo[ setinfo.snp_indices[start + j] ].aaf  << " " << all_snps_info[j].af1 << endl;
      //if(i==0 && all_snps_info[j].singleton) cerr << snpinfo[ setinfo.snp_indices[start + j] ].ID << endl;
    }

    keepaaf.col(i) = colkeep.matrix();
    //cerr << i+1 << "/" << n_aaf_bins << ": " <<  upper << "--"<< keepaaf.col(i).count()<< endl;

  }

}


void GenoMask::buildMask(const int isnp, const int chrom, struct param* params, struct filter* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& ymat, variant_block* snp_data){

  int ns = 0, hc_val, lval, nmales = 0;
  double ds, total = 0, mac = 0, mval;

  MapArXd maskvec (Gtmp.col(isnp).data(), params->n_samples, 1);
  // reset variant info
  prep_snp_stats(snp_data, params);
  snp_data->fastSPA = params->use_SPA && take_max;

  // if comphet rule, threshold to 2
  if(take_comphet) maskvec = (maskvec >= 2).select(2, maskvec);

  // if dosages were given and writing to PLINK bed, convert dosages to hardcalls
  if(params->dosage_mode && params->write_masks) maskvec = maskvec.round();

  // get counts
  for (int i = 0, index = 0; i < filters->ind_ignore.size(); i++) {

    // skip samples that were ignored from the analysis
    if( filters->ind_ignore(i) ) continue;
    ds = 0;

    if( filters->ind_in_analysis(index) ){
      if( !params->strict_mode || (params->strict_mode && masked_indivs(index,0)) ){

        ds = maskvec(index);
        // distinguish missing from 0 for sum rule
        if(!w_loo && !take_max && (ds == -3) && non_missing(index,isnp)) 
          ds = 0;

        if( ds != -3 ){
          lval = 2, mval = ds;
          if(params->test_mode && (chrom == params->nChrom)) {
            mval = ds * 0.5 * (2 - params->sex[i]);
            lval = params->sex[i];
          }
          total += ds;
          mac += mval;
          nmales += lval;
          ns++;

          // counts by trait
          if(filters->has_missing(index)) update_trait_counts(index, ds, mval, lval, 0, snp_data, masked_indivs);

          // get genotype counts (convert to hardcall)
          if( params->htp_out && (take_max || take_comphet) ) {
            hc_val = (int) (ds + 0.5); // round to nearest integer 0/1/2
            update_genocounts(params->binary_mode, index, hc_val, snp_data->genocounts, masked_indivs, ymat);
          }

        }
      }
    }

    // force masked entries to be 0
    maskvec(index++) = ds;
  }
  //cerr << maskvec.matrix().transpose().array().head(5) << endl << endl;
  if(params->write_masks) make_genovec(isnp, maskvec, filters);

  // check MAC
  if(params->test_mode){
    if(chrom != params->nChrom) mac = total; // use MAC assuming diploid coding

    // only do this when masks is in [0,2]
    if(take_max || take_comphet){

      // get counts by trait 
      snp_data->mac += mac; // aac
      snp_data->ns += ns; // ns
      snp_data->nmales += nmales; // nmales

      if(chrom != params->nChrom) {
        mac = min( mac, 2 * ns - mac );
        snp_data->mac = snp_data->mac.min( 2 * snp_data->ns.cast<double>() - snp_data->mac );
      } else {
        mac = min(mac, 2 * ns - nmales - mac); // males are 0/1
        snp_data->mac = snp_data->mac.min( 2 * snp_data->ns.cast<double>() - snp_data->nmales.cast<double>() - snp_data->mac );
      }
    }

    if(mac < params->min_MAC_mask) { 
      snp_data->ignored = true; return;
    }
    snp_data->ignored_trait = snp_data->mac < params->min_MAC_mask;
  }

  // get counts by trait 
  snp_data->af += total;

  total /= ns;
  snp_data->af1 = total / 2; // all traits
  snp_data->af /= 2 * snp_data->ns.cast<double>(); // single trait

  if(!take_max && !take_comphet) {
    snp_data->af1 /= nsites(isnp); // take average AAF across sites for sum rule
    snp_data->af /= nsites(isnp); 
  }


  if(params->use_SPA) {
    // switch to minor allele
    snp_data->flipped = ((!take_max && !take_comphet) || (params->test_type > 0)) ? false : (total > 1); // skip for DOM/REC test

    if(snp_data->flipped){
      maskvec = ( maskvec != -3).select( 2 -  maskvec, maskvec);
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

    total = ((maskvec != -3) && filters->ind_in_analysis).select(maskvec, 0).sum() / ns;
    if(total < params->numtol) {
      snp_data->ignored = true;
      return;
    }
  }


  // deal with missing data & prep for spa
  for( size_t i = 0; i < params->n_samples; ++i ) {
    ds = maskvec(i);

    // keep track of number of entries filled so avoid using clear
    if( params->use_SPA && (snp_data->fastSPA) && filters->ind_in_analysis(i) && ds > 0 ) 
      update_nnz_spa(i, params->n_samples, snp_data);

    // impute missing
    mean_impute_g(maskvec(i), total, filters->ind_in_analysis(i), masked_indivs(i,0), params->strict_mode);

  }

}



// compute MAF from AAF
void GenoMask::get_mafs(const int bs, ArrayXd& mafvec, vector<variant_block> &all_snps_info){

  for(int j = 0; j < bs; j++){
    mafvec(j) = min( all_snps_info[j].af1, 1 - all_snps_info[j].af1 );
  }

}


void GenoMask::write_info(struct param* params, struct filter* filters, mstream& sout){

  // write fam file
  write_famfile(params, filters, sout);

  // prepare ofstream for bim file
  string fname = gfile_prefix + ".bim";
  outfile_bim.open(fname.c_str(), std::ios::out);

  // write magic number to bed file
  uchar header[3] = {0x6c, 0x1b, 0x01};
  fname = gfile_prefix + ".bed";
  outfile_bed.open(fname.c_str(), std::ios::out | std::ios::binary);
  outfile_bed.write( reinterpret_cast<char*> (&header[0]), sizeof(uchar) * 3);

  // number of bytes [=ceil(N/4.0)]
  gblock_size = (filters->ind_in_analysis.count() + 3) >> 2;
  gvec.resize(nmasks_total);
  for(int i = 0; i < nmasks_total; i++) gvec[i].resize(gblock_size);

}

// write to fam
void GenoMask::write_famfile(struct param* params, struct filter* filters, mstream& sout){

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
        "0\t0\t" << params->sex[i] << "\t-9\n"; // 1 is male and 0 is female/unknown
    }

    index++;
  }

  out.closeFile();

  params->FIDvec.clear();
}

// convert to bits
void GenoMask::make_genovec(const int isnp, Ref<const ArrayXd> mask, struct filter* filters){

  int byte, bit_start, hc;
  setBitsZero(isnp);

  for(int i = 0, index = 0; i < mask.size(); i++){
    if( !filters->ind_in_analysis(i) ) continue;

    byte = index >> 2;
    bit_start = (index & 3) <<1; 

    hc = (int) (mask(i) + 0.5); // round to nearest int
    // ignore 2 since it corresponds to 00 (count number of ref alleles)
    if(hc != 2) set_gvalue(isnp, byte, bit_start, hc);

    index++;
  }

}

void GenoMask::write_genovec(const int isnp){

  outfile_bed.write( reinterpret_cast<char*> (&gvec[isnp][0]), gblock_size);

}

void GenoMask::set_gvalue(const int isnp, const int byte, const int bit_start, const int val){

  if(val == -3) gvec[isnp][byte] |= (1<<bit_start); // 01
  else if(val == 1) gvec[isnp][byte] |= (2<<bit_start); //10
  else if(val == 0) gvec[isnp][byte] |= (3<<bit_start); //11

}

void GenoMask::setBitsZero(const int isnp){
  for(size_t i = 0; i < gblock_size; i++) gvec[isnp][i] &= 0u;
}

void GenoMask::build_map(map<string,vector<int>>& mask_map){

  for(size_t i = 0; i < masks.size(); i++){
    vector<int> myints;
    for(int j = 0; j < n_aaf_bins; j++){
      // collect indices
      int index_start = i * n_aaf_bins + j;
      myints.push_back(index_start);
    }
    // insert in map
    mask_map.insert( std::make_pair( masks[i].name, myints ) );
  }

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
    if( tmp_str_vec.size() < 2 ){
      sout << "ERROR: Line " << lineread << " has too few entries.\n" ;
      exit(EXIT_FAILURE);
    }

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

  if(nfiles < 1) {
    sout << "ERROR : All set list files have unknown masks.\n";
    exit(EXIT_FAILURE);
  }

  sout << " n_files = " << nfiles << endl;
  write_setlist = true;

  // open file for writing
  setfiles.resize(nfiles);
  for(int i = 0; i < nfiles; i++) {
    line = prefix + "_" + suffix[i] + ".setlist";
    setfiles[i] = new Files;
    setfiles[i]->openForWrite( line, sout );
  }
  list_masks.resize(nfiles);

}


void GenoMask::write_genobim(const struct snp tsnp){

  // write mask info to bim file using ref-last
  // CHR ID 0 BP ALT REF 
  outfile_bim << tsnp.chrom << "\t" << tsnp.ID << "\t0\t" << tsnp.physpos << "\t" << tsnp.allele2 << "\t" << tsnp.allele1 << endl;

}

void GenoMask::append_setlist(int imask, string mname){

  // add mask name
  for(size_t i = 0; i < setfiles_index[imask].size(); i++)
    list_masks[ setfiles_index[imask][i] ].push_back(mname);

}

void GenoMask::make_setlist(string sname, int chr, uint32_t pos){

  for(size_t i = 0; i < setfiles.size(); i++){
    std::ostringstream buffer;

    // add set name
    buffer << sname << " " << chr << " " << pos << " " ;

    // add masks
    for(size_t j = 0; j < list_masks[i].size(); j++)
      buffer << list_masks[i][j] << ((j+1) == list_masks[i].size() ? "":",");


    (*setfiles[i]) << buffer.str() << endl;
  }

}

void GenoMask::closeFiles(){
  outfile_bim.close();
  outfile_bed.close();
  if(write_setlist){
    for(size_t i = 0; i < setfiles.size(); i++) {
      setfiles[i]->closeFile();
      delete setfiles[i];
    }
  }
}
