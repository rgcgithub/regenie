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

#include <limits.h> /* for PATH_MAX */
#include <chrono>
#include <boost/filesystem.hpp>
#include "Regenie.hpp"
#include "Geno.hpp"
#include "Joint_Tests.hpp"
#include "Step1_Models.hpp"
#include "Step2_Models.hpp"
#include "Files.hpp"
#include "Pheno.hpp"
#include "Masks.hpp"
#include "Data.hpp"

using namespace std;
using namespace Eigen;
using namespace boost;
namespace fs = boost::filesystem;


template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

using namespace Eigen;
using boost::math::normal;
using boost::math::chi_squared;

Data::Data() { // @suppress("Class members should be properly initialized")
}

Data::~Data() {
  // TODO Auto-generated destructor stub
}


void Data::run() {

  if(params.test_mode) {
    if(params.streamBGEN) check_bgen(files.bgen_file, &params);

    if( params.snp_set ) test_joint();
    else test_snps_fast();
    //else test_snps(); // disabled after v1.0.7


  } else {
    sout << "Fitting null model\n";

    // set number of threads
    setNbThreads(params.threads);
    // set up file for reading
    file_read_initialization();
    // if splitting l0 into many jobs
    if(params.split_l0) set_parallel_l0();
    // read phenotype and covariate files
    read_pheno_and_cov(&files, &params, &in_filters, &pheno_data, &m_ests, sout);
    // adjust for covariates
    prep_run(&files, &params, &pheno_data, &m_ests, sout);
    // set number of blocks and block size and ridge parameters
    set_blocks();
    // some initializations
    setmem();
    // level 0
    level_0_calculations();
    // level 1 ridge
    if(!params.binary_mode) // QT
      if(params.use_loocv) ridge_level_1_loocv(&files, &params, &pheno_data, &l1_ests, sout);
      else ridge_level_1(&files, &params, &l1_ests, sout);
    else // BT
      if(params.use_loocv) ridge_logistic_level_1_loocv(&files, &params, &pheno_data, &m_ests, &l1_ests, sout);
      else ridge_logistic_level_1(&files, &params, &pheno_data, &l1_ests, masked_in_folds, sout);
    // output results
    output();
  }
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          read in files
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::file_read_initialization() {

  // prepare genotype data
  files.chr_counts.assign(params.nChrom, 0.0);

  // for l0 in parallel
  if(params.run_l0_only) prep_parallel_l0();

  try {
    if(params.file_type == "bed") read_bed_bim_fam(&files, &params, &in_filters, snpinfo, chr_map, sout);
    else if(params.file_type == "pgen") read_pgen_pvar_psam(&files, &params, &in_filters, &Gblock, snpinfo, chr_map, sout);
    else prep_bgen(&files, &params, &in_filters, snpinfo, chr_map, Gblock.bgen, sout);
  } catch (bad_alloc& badAlloc)
  {
    cerr << "ERROR: bad_alloc caught, not enough memory (" << badAlloc.what() << ")\n";
    exit(EXIT_FAILURE);
  }

  if( params.setMinINFO && !params.dosage_mode )
    sout << "WARNING: Dosages are not present in the genotype file. Option --minINFO is skipped.\n";

  params.nvs_stored = snpinfo.size();

  if(!params.test_mode && !params.force_run && ((int)params.nvs_stored > params.max_step1_variants)){
    sout << "ERROR: It is not recommened to use more than " << params.max_step1_variants << " variants in step 1 (otherwise use '--force-step1'). " << params.webinfo << endl;
    exit(EXIT_FAILURE);
  }

}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          adjust for covariates in G
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::residualize_genotypes() {

  sout << "   -residualizing and scaling genotypes..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  // mask missing individuals
  Gblock.Gmat.array().rowwise() *= in_filters.ind_in_analysis.matrix().transpose().array().cast<double>();
  if(params.strict_mode) Gblock.Gmat.array().rowwise() *= pheno_data.masked_indivs.col(0).matrix().transpose().array().cast<double>();

  // residuals (centered)
  MatrixXd beta = Gblock.Gmat * pheno_data.new_cov;
  Gblock.Gmat -= beta * pheno_data.new_cov.transpose();

  // scaling (use [N-C] where C=#covariates)
  if(params.strict_mode) scale_G = Gblock.Gmat.rowwise().norm() / sqrt(pheno_data.Neff(0) - params.ncov);
  else scale_G = Gblock.Gmat.rowwise().norm() / sqrt(in_filters.ind_in_analysis.cast<double>().sum() - params.ncov);

  // check sd
  MatrixXd::Index minIndex;
  if(scale_G.array().minCoeff(&minIndex) < params.numtol) {

    if(!params.test_mode) {

      sout << "!! Uh-oh, SNP " << snpinfo[in_filters.step1_snp_count+minIndex].ID << " has low variance (=" << scale_G(minIndex,0) << ").\n";
      exit( EXIT_FAILURE );

    } else {

      Gblock.bad_snps = scale_G.array() < params.numtol;
      for(int i = 0; i < Gblock.Gmat.rows(); i++){
        // make snps polymorphic
        if(Gblock.bad_snps(i)) {
          Gblock.Gmat.row(i).head(20).array() += 1; // +1 to first 20 entries
          scale_G(i) = 1;
          if(params.verbose) sout << "WARNING: Ignoring SNP with low variance.\n";
        }
      }
    }

  }

  Gblock.Gmat.array().colwise() /= scale_G.array();

  // to use MAF dependent prior on effect size [only for step 1]
  // multiply by [p*(1-p)]^(1+alpha)/2
  if( !params.test_mode && (params.alpha_prior != -1) ){
    Gblock.Gmat.array().colwise() *= pow(Gblock.snp_afs.col(0).array() * (1-Gblock.snp_afs.col(0).array()), 0.5 * (params.alpha_prior + 1) );
  }

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;

}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          step 1: prepare for level 0
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::set_parallel_l0(){

  // compute the number of blocks
  set_blocks();

  // Make master file for L0 jobs
  write_l0_master();

  // exit software
  exit_early();
}

void Data::write_l0_master(){

  string fout = files.split_file + ".master";

  sout << " * running level 0 in parallel across " << params.total_n_block << " genotype blocks\n";
  if(params.njobs <= 1){
    sout << "ERROR: Number of jobs must be >1.\n";
    exit(EXIT_FAILURE);
  } else if(params.njobs > params.total_n_block){
    sout << "   -WARNING: Number of jobs cannot be greater than number of blocks.\n";
    params.njobs = params.total_n_block;
  }
  sout << "   -using " << params.njobs << " jobs\n";
  sout << "   -master file written to [" << fout << "]\n";
  sout << "   -variant list files written to [" << files.split_file << "_job*.snplist]\n";

  // open master
  ofstream ofile;
  openStream_write(&ofile, fout, ios::out, sout);

  // header
  ofile << params.nvs_stored << " " << params.block_size << endl;

  // split blocks in chunks of ~B/njobs
  int nall = params.total_n_block / params.njobs;
  int remainder = params.total_n_block - nall * params.njobs;
  int nb = 0, ns = 0, bcount = 0, scount = 0, jcount = 0;
  map<int, vector<int> >::iterator itr;

  int btarget = nall + (bcount < remainder ? 1 : 0);
  for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) {
    int chrom_nsnps = itr->second[0];
    int chrom_nb = ceil(chrom_nsnps * 1.0 / params.block_size);
    if(chrom_nb == 0) continue;

    for(int bb = 0; bb < chrom_nb ; bb++) {

      int bs = params.block_size;
      if((bb + 1) * params.block_size > chrom_nsnps) 
        bs = chrom_nsnps - (bb * params.block_size) ;

      ns+=bs;
      nb++, bcount++;

      if( nb == btarget ){
        string fname = files.split_file + "_job" + to_string( jcount+1 );
        // write in master
        ofile << fname << " " << btarget << " " << ns << endl;
        // write snplist
        write_snplist(fname, scount, ns);

        jcount++;
        scount += ns;
        ns = nb = 0;
        btarget = nall + (jcount < remainder ? 1 : 0);
      }
    }
  }

  if((bcount != params.total_n_block) || (jcount !=params.njobs)){
    sout << "ERROR: could not create master file.\n";
    exit(EXIT_FAILURE);
  }

  ofile.close();

}

void Data::write_snplist(string fname, int start, int ns){

  string fout = fname + ".snplist";
  ofstream ofile;
  openStream_write(&ofile, fout, ios::out, sout);

  for(int i = 0; i < ns; i++) 
    ofile << snpinfo[start+i].ID << endl;

  ofile.close();
}

void Data::set_blocks() {

  params.total_n_block = 0, total_chrs_loco = 0;
  int blocks_left = params.n_block;
  uint32_t n_analyzed = in_filters.ind_in_analysis.cast<int>().sum();
  map<int, vector<int> >::iterator itr;
  map<int, vector<int> > m1;

  // compute number of blocks for each chromosome
  for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) {
    int chrom_nsnps = itr->second[0];
    int nb = ceil(chrom_nsnps * 1.0 / params.block_size);
    if(params.n_block > 0) {
      if(blocks_left > 0) {
        int minb = min(nb, blocks_left);
        //sout << << endl;
        itr->second[1] = minb;
        params.total_n_block += minb;
        blocks_left -= minb;
      }
    } else {
      itr->second[1] = nb;
      params.total_n_block += nb;
    }

    // track how many chromosome will have blups
    if(itr->second[1] > 0) total_chrs_loco++;
    m1.insert(pair<int, vector<int> >(itr->first, itr->second));
  }
  chr_map = m1;
  //sout << "#chrs = "<< chr_map.size() << ";#loco chrs = "<< total_chrs_loco << endl;

  if(params.total_n_block == 0){
    sout << "ERROR: Total number of blocks must be > 0.\n";
    exit( EXIT_FAILURE );
  }

  if(params.split_l0) return;
  else if(params.run_l0_only) {
    if((params.parallel_nBlocks != params.total_n_block) || (params.parallel_nSnps!= (int)params.n_variants)){
    sout << "ERROR: Number of variants/blocks in file don't match with that in master file.\n";
      exit( EXIT_FAILURE );
    }
  } else if(params.run_l1_only) prep_parallel_l1();

  // set ridge params
  for(int i = 0; i < params.n_ridge_l0; i++) params.lambda[i] = (params.run_l0_only ? params.parallel_nGeno : params.n_variants) * (1 - params.lambda[i]) / params.lambda[i];
  for(int i = 0; i < params.n_ridge_l1; i++) {
    params.tau[i] =  (params.total_n_block *  params.n_ridge_l0) * (1 - params.tau[i]) / params.tau[i];
    // Assuming input tau[i] is total SNP heritability on the liability scale= m * 3/pi^2 * (1-h2) / h2
    if(params.binary_mode) params.tau[i] *= 3 / (M_PI * M_PI);
  }

  // for BTs: check if the sample size is lower than 5K (if so, use loocv)
  if( params.binary_mode && ( n_analyzed < 5000) ) {
    if(!params.use_loocv){
      sout << "   -WARNING: Sample size is less than 5,000 so using LOOCV instead of " << params.cv_folds << "-fold CV.\n";
      params.use_loocv = true;
    }
  }

  /*
  // check block size vs sample size
  if(params.use_loocv && params.block_size > n_analyzed){
    sout << "ERROR: Block size must be smaller than the number of samples to perform LOOCV!\n";
    exit( EXIT_FAILURE );
  }
  */
  if(params.use_loocv) params.cv_folds = params.n_samples;

  uint32_t neff_folds = params.use_loocv ? n_analyzed : params.cv_folds;

  // summarize block sizes and ridge params
  sout << left << std::setw(20) << " * # threads" << ": [" << params.threads << "]\n";
  sout << left << std::setw(20) << " * block size" << ": [" << params.block_size << "]\n";
  sout << left << std::setw(20) << " * # blocks" << ": [" << params.total_n_block << "]\n";
  sout << left << std::setw(20) << " * # CV folds" << ": [" << neff_folds << "]\n";
  if(!params.run_l1_only){
    int nv_tot = (params.run_l0_only ? params.parallel_nGeno : params.n_variants);
    sout << left << std::setw(20) << " * ridge data_l0" << ": [" << params.n_ridge_l0 << " : ";
    for(int i = 0; i < params.n_ridge_l0; i++) sout << nv_tot / ( nv_tot + params.lambda[i]) << " ";
    sout << "]\n";
  }
  if(!params.run_l0_only){
    sout << left << std::setw(20) << " * ridge data_l1" << ": [" << params.n_ridge_l1 << " : ";
    for(int i = 0; i < params.n_ridge_l1; i++) {
      if(!params.binary_mode)
        sout << (params.total_n_block *  params.n_ridge_l0) / (params.total_n_block *  params.n_ridge_l0 + params.tau[i] ) << " ";
      else
        sout << (params.total_n_block *  params.n_ridge_l0) / ( (params.total_n_block *  params.n_ridge_l0) + (M_PI * M_PI) * params.tau[i] / 3 ) << " ";
    }
    sout << "]\n";
  }

  // if using maf dependent prior
  if(!params.test_mode && (params.alpha_prior != -1) ) sout << " * applying a MAF dependent prior to the SNP effect sizes in level 0 models (alpha=" << params.alpha_prior << ")\n";

  // print approx. amount of memory needed
  print_usage_info(&params, &files, sout);

  // if within sample predictions are used in level 1
  if (params.within_sample_l0) {
    sout << " * using within-sample predictions from level 0 as features at level 1\n";
  }

}


void Data::set_folds() {
  // set up folds
  params.cv_sizes.resize(params.cv_folds);

  if(params.strict_mode && !params.use_loocv){
    // assign folds accounting for missingness
    uint32_t target_size_folds = floor( pheno_data.Neff(0) / params.cv_folds );
    if( target_size_folds < 1 ){
      sout << "ERROR: Not enough samples are present for " << params.cv_folds<<"-fold CV.\n";
      exit( EXIT_FAILURE );
    }

    uint32_t n_non_miss = 0, cum_size_folds = 0;
    int cur_fold = 0;
    for(size_t i = 0; i < params.n_samples; i++){

      if( pheno_data.masked_indivs(i, 0) ) n_non_miss++;

      if( n_non_miss == target_size_folds){
        params.cv_sizes[cur_fold] = i - cum_size_folds + 1;
        cum_size_folds += params.cv_sizes[cur_fold];
        n_non_miss = 0, cur_fold++;
      } else if( cur_fold == (params.cv_folds - 1) ){
        params.cv_sizes[cur_fold] = params.n_samples - i;
        break;
      }

      //sout << i << " " << cur_fold << " " << n_non_miss << " " << (pheno_data.masked_indivs(i, 0)?"Y":"N") << " "<< target_size_folds << endl;
    }

  } else {

    // assign folds for individuals in analysis
    if( !params.use_loocv ){

      uint32_t target_size_folds = floor( in_filters.ind_in_analysis.cast<double>().sum() / params.cv_folds );
      if( target_size_folds < 1 ){
        sout << "ERROR: Not enough samples are present for " << params.cv_folds<<"-fold CV.\n";
        exit( EXIT_FAILURE );
      }

      uint32_t n_non_miss = 0, cum_size_folds = 0;
      int cur_fold = 0;
      for(size_t i = 0; i < params.n_samples; i++){

        if( in_filters.ind_in_analysis(i) ) n_non_miss++;

        if( n_non_miss == target_size_folds){
          params.cv_sizes[cur_fold] = i - cum_size_folds + 1;
          cum_size_folds += params.cv_sizes[cur_fold];
          n_non_miss = 0, cur_fold++;
        } else if( cur_fold == (params.cv_folds - 1) ){
          params.cv_sizes[cur_fold] = params.n_samples - i;
          break;
        }

        //sout << i << " " << cur_fold << " " << n_non_miss << " " << in_filters.ind_in_analysis(i) << " "<< target_size_folds << endl;
      }

    } else { // loocv
      for(int i = 0; i < params.cv_folds; i++) params.cv_sizes[i] = 1;
    }

  }


  // check sd(Y) in folds
  if(!params.use_loocv && params.binary_mode){

    uint32_t cum_size_folds = 0;
    MatrixXd phenos = ( pheno_data.phenotypes_raw.array() * pheno_data.masked_indivs.array().cast<double>()).matrix();

    for(int i = 0; i < (params.cv_folds - 1); i++) {
      ArrayXd sum = phenos.block(cum_size_folds,0,params.cv_sizes[i],params.n_pheno).colwise().sum();
      ArrayXd n_cv = pheno_data.masked_indivs.block(cum_size_folds,0,params.cv_sizes[i],params.n_pheno).cast<double>().colwise().sum();
      ArrayXd sd_phenos = (sum/n_cv) * (1 - sum/n_cv);

      if( sd_phenos.minCoeff() < params.numtol ){
        sout << "ERROR: One of the folds has only cases/controls! Either use smaller #folds (option --cv) or use LOOCV (option --loocv).\n";
        exit( EXIT_FAILURE );
      }
      cum_size_folds += params.cv_sizes[i];
    }

  }

  // only used for K-fold CV
  if(!params.use_loocv && !params.within_sample_l0){
    l1_ests.X_folds.resize(params.cv_folds);
    l1_ests.XtY.resize(params.cv_folds);	
  }

}


void Data::setmem() {
  sout << " * setting memory..." << flush;

  set_folds();
  l1_ests.cumsum_values.resize(6);
  predictions.resize(1);
  predictions[0] = MatrixXd::Zero(params.n_samples, total_chrs_loco);

  if (params.within_sample_l0) {
    l1_ests.pred_mat.resize(params.n_pheno);
    l1_ests.pred_pheno.resize(params.n_pheno);
  } else if(!params.use_loocv) l1_ests.beta_hat_level_1.resize(params.n_pheno);

  if(!params.use_loocv) {
    l1_ests.test_pheno.resize(params.n_pheno);
    l1_ests.test_mat.resize(params.n_pheno);
  } else l1_ests.test_mat_conc.resize(params.n_pheno);

  if(params.binary_mode){
    if (params.within_sample_l0) {
      l1_ests.pred_pheno_raw.resize(params.n_pheno);
      l1_ests.pred_offset.resize(params.n_pheno);
    }
    l1_ests.test_pheno_raw.resize(params.n_pheno);
    if(!params.use_loocv) l1_ests.test_offset.resize(params.n_pheno);
  }
  masked_in_folds.resize(params.cv_folds);
  if(params.print_block_betas) params.beta_print_out.resize(params.n_pheno);

  for(int i = 0; i < params.n_pheno; ++i ) {

    if (params.within_sample_l0) {
      l1_ests.pred_mat[i].resize(params.cv_folds);
      l1_ests.pred_pheno[i].resize(params.cv_folds);
    } else if(!params.use_loocv) l1_ests.beta_hat_level_1[i].resize(params.cv_folds);

    if(!params.use_loocv) {
      l1_ests.test_pheno[i].resize(params.cv_folds);
      l1_ests.test_mat[i].resize(params.cv_folds);
    } else l1_ests.test_mat_conc[i] = MatrixXd::Zero(params.n_samples, params.n_ridge_l0 * ( params.write_l0_pred ? 1 : params.total_n_block) );

    if(params.binary_mode) {
      if (params.within_sample_l0) {
        l1_ests.pred_pheno_raw[i].resize(params.cv_folds);
        l1_ests.pred_offset[i].resize(params.cv_folds);
      }
      l1_ests.test_pheno_raw[i].resize(params.cv_folds);
      if(!params.use_loocv) l1_ests.test_offset[i].resize(params.cv_folds);
    }

    for(int j = 0; j < params.cv_folds; ++j ) {

      if (params.within_sample_l0) {
        l1_ests.pred_mat[i][j] = MatrixXd::Zero(params.n_samples - params.cv_sizes[j], params.total_n_block * params.n_ridge_l0);
        l1_ests.pred_pheno[i][j] = MatrixXd::Zero(params.n_samples - params.cv_sizes[j], 1);
      } else if(!params.use_loocv) l1_ests.beta_hat_level_1[i][j] = MatrixXd::Zero(params.total_n_block * params.n_ridge_l0, params.n_ridge_l1);

      if(!params.use_loocv) {
        l1_ests.test_pheno[i][j] = MatrixXd::Zero(params.cv_sizes[j], 1);
        l1_ests.test_mat[i][j] = MatrixXd::Zero(params.cv_sizes[j], params.n_ridge_l0 * ( params.write_l0_pred ? 1 : params.total_n_block));
      }

      if(params.binary_mode) {
        if (params.within_sample_l0) {
          l1_ests.pred_pheno_raw[i][j] = MatrixXd::Zero(params.n_samples - params.cv_sizes[j], 1);
          l1_ests.pred_offset[i][j] = MatrixXd::Zero(params.n_samples - params.cv_sizes[j], 1);
        }
        l1_ests.test_pheno_raw[i][j] = MatrixXd::Zero(params.cv_sizes[j], 1);
        if(!params.use_loocv) l1_ests.test_offset[i][j] = MatrixXd::Zero(params.cv_sizes[j], 1);
      }

      if(i == 0) masked_in_folds[j] = MatrixXb::Constant(params.cv_sizes[j], params.n_pheno, false);

    }
  }
  sout << "done\n" << endl;
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          step 1: level 0
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::level_0_calculations() {

  if(params.run_l1_only) {
    set_mem_l1(&files, &params, &in_filters, &m_ests, &Gblock, &pheno_data, &l1_ests, masked_in_folds, sout);
    sout << " (skipping to level 1 models)";
    return;
  }

  int block = 0;
  if(params.print_block_betas) params.print_snpcount = 0;
  ridgel0 l0;

  if(!params.use_loocv){
    l0.G_folds.resize(params.cv_folds);
    l0.GtY.resize(params.cv_folds);
  }

  // open streams to write level 0 predictions
  if(params.write_l0_pred){
    string fout_p;
    files.write_preds_files.resize(params.n_pheno);
    for(int ph = 0; ph < params.n_pheno; ph++){
      files.write_preds_files[ph] = std::make_shared<ofstream>();
      fout_p = files.loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
      openStream_write(files.write_preds_files[ph].get(), fout_p, ios::out | ios::binary, sout);
    }
  }

  // start level 0
  for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {

    int chrom = files.chr_read[itr];
    if( !in_map(chrom, chr_map) ) continue;

    int chrom_nsnps = chr_map[chrom][0];
    int chrom_nb = chr_map[chrom][1];
    if(chrom_nb == 0) continue;

    sout << "Chromosome " << chrom << endl;
    //sout << "Ns="<< chrom_nsnps << endl;

    for(int bb = 0; bb < chrom_nb ; bb++) {

      int bs = params.block_size;
      if((bb + 1) * params.block_size > chrom_nsnps) 
        bs = chrom_nsnps - (bb * params.block_size) ;

      Gblock.Gmat = MatrixXd::Zero(bs, params.n_samples);
      if(params.alpha_prior != -1) Gblock.snp_afs = MatrixXd::Zero(bs, 1);

      get_G(block, bs, chrom, in_filters.step1_snp_count, snpinfo, &params, &files, &Gblock, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, sout);

      // residualize and scale genotypes
      residualize_genotypes();

      // calc working matrices for ridge regressions across folds
      calc_cv_matrices(bs, &l0);

      // calc level 0 ridge regressions
      if(params.use_loocv)
        ridge_level_0_loocv(block, &files, &params, &in_filters, &m_ests, &Gblock, &pheno_data, snpinfo, &l0, &l1_ests, sout);
      else
        ridge_level_0(block, &files, &params, &in_filters, &m_ests, &Gblock, &pheno_data, snpinfo, &l0, &l1_ests, masked_in_folds, sout);

      block++; in_filters.step1_snp_count += bs;
    }

  }

  // close streams
  if(params.write_l0_pred){
    for(int ph = 0; ph < params.n_pheno; ph++)
      files.write_preds_files[ph]->close();
  }

  if(params.early_exit) {
    sout << "\nDone printing out level 0 predictions. There are " <<
      params.n_samples << " rows and " <<
      params.total_n_block * params.n_ridge_l0 << " columns " <<
      "stored in column-major order. Exiting...\n";
    exit_early();
  } else if(params.run_l0_only) exit_early();

  // free up memory not used anymore
  Gblock.Gmat.resize(0,0);
  if(params.write_l0_pred && (params.n_pheno > 1) ) {
    // free level 0 predictions for (P-1) indices in test_mat
    for(int ph = 1; ph < params.n_pheno; ++ph ) {
      if(!params.use_loocv){
        for(int i = 0; i < params.cv_folds; ++i ) l1_ests.test_mat[ph][i].resize(0,0);
        l1_ests.test_mat[ph].resize(0);
      } else {
        l1_ests.test_mat_conc[ph].resize(0,0);
      }
    }
  }

}

void Data::calc_cv_matrices(const int bs, struct ridgel0* l0) {

  sout << "   -calc working matrices..." << flush;
  auto t2 = std::chrono::high_resolution_clock::now();

  if(!params.use_loocv){
    l0->GGt.setZero(bs,bs);
    l0->GTY.setZero(bs,params.n_pheno);
    uint32_t cum_size_folds = 0;

    for( int i = 0; i < params.cv_folds; ++i ) {
      l0->GtY[i] = Gblock.Gmat.block(0, cum_size_folds, bs, params.cv_sizes[i]) * pheno_data.phenotypes.block(cum_size_folds, 0, params.cv_sizes[i], params.n_pheno);
      l0->G_folds[i] = Gblock.Gmat.block(0, cum_size_folds, bs, params.cv_sizes[i]) * Gblock.Gmat.block(0, cum_size_folds, bs, params.cv_sizes[i]).transpose();
      l0->GGt += l0->G_folds[i];
      l0->GTY += l0->GtY[i];
      cum_size_folds += params.cv_sizes[i];
    }
  } else {
    l0->GGt = Gblock.Gmat * Gblock.Gmat.transpose();
    l0->GTY = Gblock.Gmat * pheno_data.phenotypes;
    SelfAdjointEigenSolver<MatrixXd> esG(l0->GGt);
    l0->GGt_eig_vec = esG.eigenvectors();
    l0->GGt_eig_val = esG.eigenvalues();
    l0->Wmat = l0->GGt_eig_vec.transpose() * l0->GTY;
  }

  sout << "done";
  auto t3 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
  sout << " (" << duration.count() << "ms) "<< endl;
}

// identify which block to analyze
void Data::prep_parallel_l0(){

  int tmpi;
  string line;
  string fin = files.split_file; // master file
  std::vector< string > tmp_str_vec ;
  ifstream infile;

  // print info
  sout << " * running jobs in parallel (job #" << params.job_num << ")\n";

  infile.open(fin.c_str(), ios::in);
  if (!infile.is_open()) {
    sout << "ERROR : Cannot read file " << fin  << endl ;
    exit(EXIT_FAILURE);
  }

  // check header
  if(!getline(infile, line)){
    sout << "ERROR: Cannot read header line in master file.\n"; 
    exit(EXIT_FAILURE);
  }
  if( (sscanf( line.c_str(), "%d %d", &params.parallel_nGeno, &tmpi ) != 2) || (tmpi != params.block_size) ){
    sout << "ERROR: Invalid header line.\n"; 
    exit(EXIT_FAILURE);
  }

  // skip to line job_num
  int nskip=1;
  while( (nskip++ < params.job_num) && !infile.eof() )
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  if( (--nskip != params.job_num) || infile.eof() ){
    sout << "ERROR: Could not read line " << params.job_num+1 << " (check number of lines in file).\n";
    exit(EXIT_FAILURE);
  }
  
  // read in line
  getline(infile, line);
  char tmp_chr[MAXFILELEN];
  if( sscanf( line.c_str(), "%s %d %d", tmp_chr, &params.parallel_nBlocks, &params.parallel_nSnps ) != 3 ){
    sout << "ERROR: Could not read line " << params.job_num + 1 << " (check number of lines and format in file).\n"; 
    exit(EXIT_FAILURE);
  }
  files.loco_tmp_prefix = tmp_chr;
  files.file_snps_include = files.loco_tmp_prefix + ".snplist";

  infile.close();
  //cerr << files.loco_tmp_prefix << " " << params.parallel_nBlocks << " " << params.parallel_nSnps << endl;

}


void Data::prep_parallel_l1(){

  int nblocks, lineread, nb ,ns; // make sure all blocks are read
  uint32_t nsnps;
  string line;
  string fin = files.split_file; // master file
  std::vector< string > tmp_str_vec ;
  ifstream infile;

  infile.open(fin.c_str(), ios::in);
  if (!infile.is_open()) {
    sout << "ERROR : Cannot read file " << fin  << endl ;
    exit(EXIT_FAILURE);
  }

  // check header
  if(!getline(infile, line)){
    sout << "ERROR: Cannot read header line in master file.\n"; 
    exit(EXIT_FAILURE);
  }
  if( (sscanf( line.c_str(), "%d %d", &params.parallel_nGeno, &nb ) != 2) || (nb != params.block_size) ){
    sout << "ERROR: Invalid header line.\n"; 
    exit(EXIT_FAILURE);
  }

  nblocks = 0, nsnps = 0, lineread=0;
  while( getline(infile, line) ){

    char tmp_chr[MAXFILELEN];
    if( sscanf( line.c_str(), "%s %d %d", tmp_chr, &nb, &ns ) != 3 ){
      sout << "ERROR: Could not read line " << params.job_num + 1 << " (check number of lines and format in file).\n"; 
      exit(EXIT_FAILURE);
    }

    files.bstart.push_back( nblocks );
    files.btot.push_back( nb );
    files.mprefix.push_back( string(tmp_chr) );

    // check params
    if( (files.bstart[lineread] < 0) || (files.bstart[lineread]>params.total_n_block) || (files.btot[lineread] < 0) ){
      sout << "ERROR: Invalid block information in master file at line " << lineread + 2 << ".\n";
      exit(EXIT_FAILURE);
    }

    nblocks += nb; // update # blocks
    nsnps += ns;
    lineread++;
  }

  if((nblocks != params.total_n_block) || (nsnps != params.n_variants)){
    sout << "ERROR: Number of blocks/variants in master file '" << fin << "' doesn't match that in the analysis.\n";
    exit(EXIT_FAILURE);
  }

  // print info
  params.job_num = lineread;
  sout << " * using results from running " << params.job_num << " parallel jobs at level 0.\n";

  infile.close();

}

void Data::exit_early(){

    runtime.stop();
    sout << "\nElapsed time : " << std::chrono::duration<double>(runtime.end - runtime.begin).count() << "s\n";
    sout << "End time: " << ctime(&runtime.end_time_info) << endl;
    exit( EXIT_SUCCESS );

}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          Evaluate level 1 output
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::output() {

  int min_index;
  double performance_measure, rsq, sse, ll_avg, min_val;
  string pfile, out_blup_list, out_prs_list, loco_filename, prs_filename;
  string fullpath_str, path_prs;
  Files outb, outp;

  sout << "Output\n" << "------\n";

  if(params.make_loco || params.binary_mode){
    out_blup_list = files.out_file + "_pred.list";
    outb.openForWrite(out_blup_list, sout);
  }
  if(params.print_prs){
    out_prs_list = files.out_file + "_prs.list";
    outp.openForWrite(out_prs_list, sout);
  }

  for(int ph = 0; ph < params.n_pheno; ++ph ) {

    sout << "phenotype " << ph+1 << " (" << files.pheno_names[ph] << ") : " ;
    loco_filename = files.out_file + "_" + to_string(ph + 1) + ".loco" + (params.gzOut ? ".gz" : "");
    prs_filename = files.out_file + "_" + to_string(ph + 1) + ".prs" + (params.gzOut ? ".gz" : "");

    if( params.make_loco || params.binary_mode || params.print_prs ) {

      fullpath_str = get_fullpath(loco_filename);
      if(params.print_prs) path_prs = get_fullpath(prs_filename);

      if( !params.binary_mode ) { // for quantitative traits
        outb << files.pheno_names[ph]  << " " <<  fullpath_str << endl;
        if(params.print_prs) outp << files.pheno_names[ph]  << " " <<  path_prs << endl;
      } else { // for binary traits - check level 1 ridge converged
        if( !l1_ests.pheno_l1_not_converged(ph) ) {
          outb << files.pheno_names[ph]  << " " << fullpath_str << endl;
          if(params.print_prs) outp << files.pheno_names[ph]  << " " <<  path_prs << endl;
        } else {
          if(params.write_l0_pred) rm_l0_files(ph); // cleanup level 0 predictions
          sout << "Level 1 logistic did not converge. LOCO predictions calculations are skipped.\n\n";
          continue;
        }
      }

    }
    sout << endl;

    min_index = 0;
    min_val = 1e10;

    // determine optimal parameter by cv using: QT: MSE, BT: -loglik
    for(int j = 0; j < params.n_ridge_l1; ++j ) {
      if(!params.binary_mode)
        performance_measure = l1_ests.cumsum_values[2](ph, j) + l1_ests.cumsum_values[3](ph,j) - 2 * l1_ests.cumsum_values[4](ph,j);
      else
        performance_measure = l1_ests.cumsum_values[5](ph, j);

      performance_measure /= pheno_data.Neff(ph);

      if( performance_measure < min_val) {
        min_index = j;
        min_val = performance_measure;
      }
    }

    for(int j = 0; j < params.n_ridge_l1; ++j ) {
      if(!params.binary_mode)
        sout << "  " << setw(5) << (params.total_n_block *  params.n_ridge_l0)  / (params.total_n_block * params.n_ridge_l0 + params.tau[j] ) ;
      else
        sout << "  " << setw(5) << (params.total_n_block *  params.n_ridge_l0) / ( (params.total_n_block *  params.n_ridge_l0) + (M_PI * M_PI) * params.tau[j] / 3 );

      // output Rsq and MSE
      rsq = l1_ests.cumsum_values[4](ph,j) - l1_ests.cumsum_values[0](ph,j) * l1_ests.cumsum_values[1](ph,j) / pheno_data.Neff(ph); // num = Sxy - SxSy/n
      rsq = (rsq * rsq) / ((l1_ests.cumsum_values[2](ph,j) - l1_ests.cumsum_values[0](ph,j) * l1_ests.cumsum_values[0](ph,j) / pheno_data.Neff(ph)) * (l1_ests.cumsum_values[3](ph,j) - l1_ests.cumsum_values[1](ph,j) * l1_ests.cumsum_values[1](ph,j) / pheno_data.Neff(ph))); // num^2 / ( (Sx2 - Sx^2/n)* (Sy2 - Sy^2/n) )
      sse = l1_ests.cumsum_values[2](ph, j) + l1_ests.cumsum_values[3](ph,j) - 2 * l1_ests.cumsum_values[4](ph,j); // Sx2 + Sy2 - SxSy
      if(params.binary_mode) ll_avg = l1_ests.cumsum_values[5](ph, j) / pheno_data.Neff(ph);

      sout << " : ";
      sout << "Rsq = " << rsq << ", MSE = " << sse/pheno_data.Neff(ph);
      if(params.binary_mode) sout << ", -logLik/N = " << ll_avg;
      if(j == min_index) sout << "<- min value";
      sout << endl;
    }

    if(!params.binary_mode){
      if(params.use_loocv) make_predictions_loocv(ph, min_index);
      else make_predictions(ph, min_index);
    } else if(params.use_loocv) make_predictions_binary_loocv(ph, min_index);
    else make_predictions_binary(ph, min_index);


    // delete file used to store l0 predictions
    if(params.write_l0_pred && params.rm_l0_pred)
      rm_l0_files(ph);

  }

  if(params.make_loco || params.binary_mode){
    outb.closeFile();
    sout << "List of blup files written to: [" << out_blup_list << "]\n";
  }
  if(params.print_prs) {
    outp.closeFile();
    sout << "List of files with whole genome PRS written to: [" << 
      out_prs_list << "]\n";
  }

}

void Data::rm_l0_files(int ph){

  string pfile;

  if(!params.run_l1_only){
    pfile = files.loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
    remove(pfile.c_str());
  } else {
    for(size_t i = 0; i < files.bstart.size(); i++){
      pfile = files.mprefix[i] + "_l0_Y" + to_string(ph+1);
      remove(pfile.c_str());
      if(ph==0){
        pfile = files.mprefix[i] + ".snplist";
        remove(pfile.c_str());
      }
    }
  }

}

// convert filename to full path
std::string get_fullpath(std::string fname){

  string fout;

  try {

    // convert to full path using boost filesystem library
    // this can generate errors due to LC_ALL locale being invalid
    fs::path fullpath;
    fullpath = fs::absolute(fname);
    fout = fullpath.make_preferred().string();

  } catch ( std::runtime_error& ex ) {

    try {

      // use realpath
      char buf[PATH_MAX];
      char *res = realpath(fname.c_str(), buf);
      if(res) fout = string(buf);
      else fout = fname; // if failed to get full path

    } catch ( const std::bad_alloc& ) {
      fout = fname; // if failed to get full path
    }

  }

  return fout;

}

void Data::make_predictions(const int ph, const  int val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  int bs_l1 = params.total_n_block * params.n_ridge_l0;
  int ph_eff = params.write_l0_pred ? 0 : ph;
  string outname, in_pheno;
  ifstream infile;
  ofstream ofile;

  MatrixXd X1, X2, beta_l1, beta_avg;
  MatrixXd ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);

  // read in level 0 predictions from file
  if(params.write_l0_pred)
    read_l0(ph, ph_eff, &files, &params, &l1_ests, sout);


  if(params.within_sample_l0){
    X1 = l1_ests.test_mat[ph_eff][0].transpose() * l1_ests.test_mat[ph_eff][0];
    X2 = l1_ests.test_mat[ph_eff][0].transpose() * l1_ests.test_pheno[ph][0];
    for(int i = 1; i < params.cv_folds; ++i ) {
      X1 += l1_ests.test_mat[ph_eff][i].transpose() * l1_ests.test_mat[ph_eff][i];
      X2 += l1_ests.test_mat[ph_eff][i].transpose() * l1_ests.test_pheno[ph][i];
    }
    beta_l1 = (X1 + params.tau[val] * ident_l1).llt().solve(X2);
  } else if(params.print_block_betas) {
    beta_avg = MatrixXd::Zero(bs_l1, 1);
    for(int i = 0; i < params.cv_folds; ++i ) {
      beta_avg += l1_ests.beta_hat_level_1[ph][i].col(val);
    }
    beta_avg /= params.cv_folds;
  }

  // if specified, write betas to file (open in append mode)
  if(!params.within_sample_l0 && params.print_block_betas) {
    outname = files.out_file + "_level1.betas";
    openStream_write(&ofile, outname, ios::out | ios::app, sout);
    ofile << ph + 1 << " ";
    ofile << beta_avg.transpose() << endl;
    ofile.close();
  }

  // sout << "\nFor tau[" << val <<"] = " << params.tau[val] << endl <<  beta_l1 << endl ;
  int ctr = 0, chr_ctr = 0;
  int nn, cum_size_folds;

  for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {
    int chrom = files.chr_read[itr];
    if( !in_map(chrom, chr_map) ) continue;

    nn = chr_map[chrom][1] * params.n_ridge_l0;
    if(nn > 0) {
      cum_size_folds = 0;
      for(int i = 0; i < params.cv_folds; ++i ) {
        if(!params.within_sample_l0) beta_l1 = l1_ests.beta_hat_level_1[ph][i].col(val);
        predictions[0].block(cum_size_folds, chr_ctr, params.cv_sizes[i], 1) = l1_ests.test_mat[ph_eff][i].block(0, ctr, params.cv_sizes[i], nn) * beta_l1.block(ctr, 0, nn, 1);
        cum_size_folds += params.cv_sizes[i];
      }
      chr_ctr++;
      ctr += nn;
    }
  }

  write_predictions(ph);

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl << endl;
}


void Data::make_predictions_loocv(const int ph, const  int val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  int bs_l1 = params.total_n_block * params.n_ridge_l0;
  int ph_eff = params.write_l0_pred ? 0 : ph;
  string in_pheno;
  ifstream infile;
  MatrixXd Xmat_chunk, Yvec_chunk,  Z1, Z2, b0, xtx;
  VectorXd w1, w2, Vw2, zvec;
  RowVectorXd calFactor;
  ArrayXd dl_inv;

  uint64 max_bytes = params.chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (target_size * bs_l1) matrix ]
  int nchunk = ceil( params.cv_folds * bs_l1 * sizeof(double) * 1.0 / max_bytes );
  if (params.verbose) sout << nchunk << " chunks..." << flush;
  int chunk, size_chunk, target_size = params.cv_folds / nchunk;
  int j_start;


  // read in level 0 predictions from file
  if(params.write_l0_pred)
    read_l0(ph, ph_eff, &files, &params, &l1_ests, sout);

  // fit model on whole data again for optimal ridge param
  xtx = l1_ests.test_mat_conc[ph_eff].transpose() * l1_ests.test_mat_conc[ph_eff];
  SelfAdjointEigenSolver<MatrixXd> eigX(xtx);
  zvec = l1_ests.test_mat_conc[ph_eff].transpose() * pheno_data.phenotypes.col(ph);
  w1 = eigX.eigenvectors().transpose() * zvec;
  dl_inv = (eigX.eigenvalues().array() + params.tau[val]).inverse();
  w2 = (w1.array() * dl_inv).matrix();
  Vw2 = eigX.eigenvectors() * w2;

  for(chunk = 0; chunk < nchunk; ++chunk ) {
    size_chunk = chunk == nchunk - 1? params.cv_folds - target_size * chunk : target_size;
    j_start = chunk * target_size;
    Xmat_chunk = l1_ests.test_mat_conc[ph_eff].block(j_start, 0, size_chunk, bs_l1);
    Yvec_chunk = pheno_data.phenotypes.block(j_start, ph, size_chunk, 1);

    Z1 = (Xmat_chunk * eigX.eigenvectors()).transpose();
    Z2 = dl_inv.matrix().asDiagonal() * Z1;
    calFactor = (Z1.array() * Z2.array()).matrix().colwise().sum();
    b0 = eigX.eigenvectors() * Z2;
    b0.array().rowwise() *= (w2.transpose() * Z1 - Yvec_chunk.transpose()).array() / (1 - calFactor.array());
    b0.colwise() += Vw2;

    int ctr = 0, chr_ctr = 0;
    int nn;

    for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {
      int chrom = files.chr_read[itr];
      if( !in_map(chrom, chr_map) ) continue;

      nn = chr_map[chrom][1] * params.n_ridge_l0;
      if(nn > 0) {
        predictions[0].block(j_start, chr_ctr, size_chunk, 1) = (l1_ests.test_mat_conc[ph_eff].block(j_start, ctr, size_chunk, nn).array() * b0.block(ctr, 0, nn, size_chunk).transpose().array()).rowwise().sum();
        chr_ctr++;
        ctr += nn;
      }
    }
  }

  write_predictions(ph);

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl << endl;
}


// predictions for binary traits
void Data::make_predictions_binary(const int ph, const  int val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  int bs_l1 = params.total_n_block * params.n_ridge_l0;
  int ph_eff = params.write_l0_pred ? 0 : ph;
  string in_pheno;
  ifstream infile;
  ArrayXd etavec, pivec, wvec, zvec, score;
  MatrixXd betaold, betanew, XtW, XtWX, XtWZ, Xconc;
  MatrixXd ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);

  // read in level 0 predictions from file
  if(params.write_l0_pred)
    read_l0(ph, ph_eff, &files, &params, &l1_ests, sout);

  // fit model using out-of-sample level 0 predictions from whole data
  if(params.within_sample_l0){
    betaold = MatrixXd::Zero(bs_l1, 1);

    int niter_cur = 0;
    while(niter_cur++ < params.niter_max_ridge){

      XtWX = MatrixXd::Zero(bs_l1, bs_l1);
      XtWZ = MatrixXd::Zero(bs_l1, 1);

      for(int i = 0; i < params.cv_folds; ++i ) {
        etavec = (l1_ests.test_offset[ph][i] + l1_ests.test_mat[ph_eff][i] * betaold).array();
        pivec = 1 - 1/(etavec.exp() + 1);
        wvec =  pivec * (1 - pivec);
        zvec = (etavec - l1_ests.test_offset[ph][i].array()) + (l1_ests.test_pheno_raw[ph][i].array() - pivec) / wvec;

        XtW = l1_ests.test_mat[ph_eff][i].transpose() * wvec.matrix().asDiagonal();
        XtWX += XtW * l1_ests.test_mat[ph_eff][i];
        XtWZ += XtW * zvec.matrix();
      }
      betanew = (XtWX + params.tau[val] * ident_l1).llt().solve(XtWZ);
      // compute score
      score = ArrayXd::Zero(betanew.rows());
      for(int i = 0; i < params.cv_folds; ++i ) {
        etavec = (l1_ests.test_offset[ph][i] + l1_ests.test_mat[ph_eff][i] * betanew).array();
        pivec = 1 - 1/(etavec.exp() + 1);
        score += (l1_ests.test_mat[ph_eff][i].transpose() * (l1_ests.test_pheno_raw[ph][i].array() - pivec).matrix()).array();
      }
      score -= params.tau[val] * betanew.array();

      // stopping criterion
      if( score.abs().maxCoeff() < params.l1_ridge_eps) break;

      betaold = betanew;
    }
  }

  // compute predictor for each chr
  int ctr = 0, chr_ctr = 0;
  int nn, cum_size_folds;

  for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {
    int chrom = files.chr_read[itr];
    if( !in_map(chrom, chr_map) ) continue;

    nn = chr_map[chrom][1] * params.n_ridge_l0;
    if(nn > 0) {
      cum_size_folds = 0;
      for(int i = 0; i < params.cv_folds; ++i ) {
        if(!params.within_sample_l0) betanew = l1_ests.beta_hat_level_1[ph][i].col(val);
        predictions[0].block(cum_size_folds, chr_ctr, params.cv_sizes[i], 1) = l1_ests.test_mat[ph_eff][i].block(0, ctr, params.cv_sizes[i], nn) * betanew.block(ctr, 0, nn, 1);
        cum_size_folds += params.cv_sizes[i];
      }
      chr_ctr++;
      ctr += nn;
    }
  }

  write_predictions(ph);

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl << endl;
}


void Data::make_predictions_binary_loocv(const int ph, const int val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  int bs_l1 = params.total_n_block * params.n_ridge_l0;
  int ph_eff = params.write_l0_pred ? 0 : ph;
  double v2;
  string in_pheno;
  ifstream infile;

  ArrayXd betaold, etavec, pivec, wvec, zvec, betanew, score;
  MatrixXd XtWX, XtWZ, V1, beta_final;
  LLT<MatrixXd> Hinv;
  MatrixXd ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);

  uint64 max_bytes = params.chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (bs_l1 * target_size) matrix ]
  int nchunk = ceil( params.cv_folds * bs_l1 * sizeof(double) * 1.0 / max_bytes );
  int chunk, size_chunk, target_size = params.cv_folds / nchunk;
  int j_start;

  // read in level 0 predictions from file
  if(params.write_l0_pred)
    read_l0(ph, ph_eff, &files, &params, &l1_ests, sout);

  // fit logistic on whole data again for optimal ridge param
  betaold = ArrayXd::Zero(bs_l1);
  int niter_cur = 0;
  while(niter_cur++ < params.niter_max_ridge){

    get_wvec(ph, etavec, pivec, wvec, betaold, pheno_data.masked_indivs, m_ests.offset_logreg.col(ph), l1_ests.test_mat_conc[ph_eff], params.l1_ridge_eps);
    zvec = (pheno_data.masked_indivs.col(ph).array()).select( (etavec - m_ests.offset_logreg.col(ph).array()) + (pheno_data.phenotypes_raw.col(ph).array() - pivec) / wvec, 0);
    V1 = l1_ests.test_mat_conc[ph_eff].transpose() * wvec.matrix().asDiagonal();
    XtWX = V1 * l1_ests.test_mat_conc[ph_eff];
    XtWZ = V1 * zvec.matrix();
    Hinv.compute( XtWX + params.tau[val] * ident_l1 );
    betanew = (Hinv.solve(XtWZ)).array();

    // get the score
    get_wvec(ph, etavec, pivec, wvec, betanew, pheno_data.masked_indivs, m_ests.offset_logreg.col(ph), l1_ests.test_mat_conc[ph_eff], params.l1_ridge_eps);
    score = ( l1_ests.test_mat_conc[ph_eff].transpose() * (pheno_data.masked_indivs.col(ph).array()).select(pheno_data.phenotypes_raw.col(ph).array() - pivec, 0).matrix()).array() ;
    score -= params.tau[val] * betanew;

    if( score.abs().maxCoeff() < params.l1_ridge_eps) break;

    betaold = betanew;
  }

  // compute Hinv
  //zvec = (etavec - m_ests.offset_logreg.col(ph).array()) + (pheno_data.phenotypes_raw.col(ph).array() - pivec) / wvec;
  V1 = l1_ests.test_mat_conc[ph_eff].transpose() * wvec.matrix().asDiagonal();
  XtWX = V1 * l1_ests.test_mat_conc[ph_eff];
  Hinv.compute( XtWX + params.tau[val] * ident_l1 );

  // loo estimates
  for(chunk = 0; chunk < nchunk; ++chunk ) {
    size_chunk = chunk == nchunk - 1? params.cv_folds - target_size * chunk : target_size;
    j_start = chunk * target_size;
    if( (chunk == 0) || (chunk == nchunk - 1) ) beta_final = MatrixXd::Zero(bs_l1, size_chunk);

    Ref<MatrixXd> Xmat_chunk = l1_ests.test_mat_conc[ph_eff].block(j_start, 0, size_chunk, bs_l1); // n x k
    Ref<MatrixXd> Yvec_chunk = pheno_data.phenotypes_raw.block(j_start, ph, size_chunk, 1);

    V1 = Hinv.solve( Xmat_chunk.transpose() ); // k x n
    for(int i = 0; i < size_chunk; ++i ) {
      v2 = Xmat_chunk.row(i) * V1.col(i);
      v2 *= wvec(j_start + i);
      beta_final.col(i) = (betanew - V1.col(i).array() * (Yvec_chunk(i,0) - pivec(j_start + i)) / (1 - v2)).matrix();
    }

    // compute predictor for each chr
    int ctr = 0, chr_ctr = 0;
    int nn;

    for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {
      int chrom = files.chr_read[itr];
      if( !in_map(chrom, chr_map) ) continue;

      nn = chr_map[chrom][1] * params.n_ridge_l0;

      if(nn > 0) {
        predictions[0].block(j_start, chr_ctr, size_chunk, 1) = ( l1_ests.test_mat_conc[ph_eff].block(j_start, ctr, size_chunk, nn).array() * beta_final.block(ctr, 0, nn, size_chunk).transpose().array() ).matrix().rowwise().sum();
        chr_ctr++;
        ctr += nn;
      }
    }
  }

  write_predictions(ph);

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl << endl;
}


void Data::write_predictions(const int ph){
  // output predictions to file
  string out;
  Files ofile;
  MatrixXd pred, prs;

  // for the per chromosome predictions (not used)
  if(params.write_blups) {

    out = files.out_file + "_" + to_string(ph+1) + (params.gzOut ? ".gz" : "");
    sout << "writing file " << out << "..." << flush;
    ofile.openForWrite(out, sout);

    // enforce all chromosomes are printed
    pred = MatrixXd::Zero(predictions[0].rows(), params.nChrom);

    int chr, nn, chr_ctr = 0;
    for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {
      chr = files.chr_read[itr];
      if( !in_map(chr, chr_map) ) continue;

      nn = chr_map[chr][1];
      if(nn > 0){
        pred.col(chr - 1) = predictions[0].col(chr_ctr);
        ++chr_ctr;
      }
    }

    // header line : FID_IID for all individuals
    ofile << write_ID_header();

    // for each row: print chromosome then blups
    for(chr = 0; chr < params.nChrom; chr++) 
      ofile << write_chr_row(chr+1, ph, pred.col(chr));

    ofile.closeFile();

  }

  if(params.make_loco || params.binary_mode){

    out = files.out_file + "_" + to_string(ph+1) + ".loco" + (params.gzOut ? ".gz" : "");
    sout << "writing LOCO predictions..." << flush;
    ofile.openForWrite(out, sout);

    // output LOCO predictions G_loco * beta_loco for each autosomal chr
    pred.resize(predictions[0].rows(), params.nChrom);
    pred.colwise() = predictions[0].rowwise().sum();

    int chr, nn, chr_ctr = 0;
    for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {
      chr = files.chr_read[itr];
      if( !in_map(chr, chr_map) ) continue;

      nn = chr_map[chr][1];
      if(nn > 0) {
        pred.col(chr - 1) -= predictions[0].col(chr_ctr);
        ++chr_ctr;
      }
    }

    // header line : FID_IID for all individuals
    ofile << write_ID_header();

    // print loco predictions for each chromosome
    for(chr = 0; chr < params.nChrom; chr++) 
      ofile << write_chr_row(chr+1, ph, pred.col(chr));

    ofile.closeFile();

  }


  if(params.print_prs){

    out = files.out_file + "_" + to_string(ph+1) + ".prs" + (params.gzOut ? ".gz" : "");
    sout << "writing whole genome PRS..." << flush;
    ofile.openForWrite(out, sout);

    // output predictions sum(G * beta)
    prs.resize(predictions[0].rows(), 1);
    prs = predictions[0].rowwise().sum();

    // header line : FID_IID for all individuals
    ofile << write_ID_header();

    // print prs (set chr=0)
    ofile << write_chr_row(0, ph, prs.col(0));

    ofile.closeFile();

  }

}

std::string Data::write_ID_header(){

  uint32_t index;
  string out, id_index;
  std::ostringstream buffer;
  map<string, uint32_t >::iterator itr_ind;

  buffer << "FID_IID ";
  for (itr_ind = params.FID_IID_to_ind.begin(); itr_ind != params.FID_IID_to_ind.end(); ++itr_ind) {
    id_index = itr_ind->first;
    index = itr_ind->second;

    // check individual was included in analysis, if not then skip
    if( !in_filters.ind_in_analysis( index ) ) continue;
    buffer << id_index << " ";
  }
  buffer << endl;

  return buffer.str();

}


std::string Data::write_chr_row(const int chr, const int ph, const Eigen::MatrixXd& pred){

  uint32_t index;
  string out, id_index;
  std::ostringstream buffer;
  map<string, uint32_t >::iterator itr_ind;

  buffer << chr << " ";
  for (itr_ind = params.FID_IID_to_ind.begin(); itr_ind != params.FID_IID_to_ind.end(); ++itr_ind) {
    id_index = itr_ind->first;
    index = itr_ind->second;

    // check individual was included in analysis, if not then skip
    if( !in_filters.ind_in_analysis( index ) ) continue;

    // print prs
    if( pheno_data.masked_indivs(index, ph) )
      buffer << pred(index, 0) << " ";
    else
      buffer << "NA ";
  }
  buffer << endl;

  return buffer.str();

}


/*

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Testing mode (approx. single-threaded)
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::test_snps() {

  sout << "Association testing mode using parallelization in Eigen\n";
  chi_squared chisq(1);
  normal nd(0,1);
  std::chrono::high_resolution_clock::time_point t1, t2;

  double chisq_val, bhat, se_b, pval_log, pval_raw, info=0;
  double chisq_thr = quantile(complement(chisq, params.alpha_pvalue));
  double zcrit = quantile(complement(nd, .025));
  uint32_t n_failed_tests = 0;
  uint32_t n_ignored_snps = 0;
  uint32_t n_skipped_snps = 0;
  string out, tmpstr;
  vector < string > out_split;
  // output files
  Files ofile;
  // use pointer to class since it contains non-copyable elements
  vector < Files* > ofile_split;
  MatrixXd WX, GW, sqrt_denum, scaleG_pheno;

  setNbThreads(params.threads); // set threads
  file_read_initialization(); // set up files for reading
  read_pheno_and_cov(&files, &params, &in_filters, &pheno_data, &m_ests, sout);   // read phenotype and covariate files
  prep_run(&files, &params, &pheno_data, &m_ests, sout); // check blup files and adjust for covariates
  set_blocks_for_testing();   // set number of blocks
  print_usage_info(&params, &files, sout);
  print_test_info();
  setup_output(&ofile, out, ofile_split, out_split); // result file

  // set memory for matrices used to store estimates under H0
  m_ests.Y_hat_p = MatrixXd::Zero(params.n_samples, params.n_pheno);
  m_ests.Gamma_sqrt = MatrixXd::Zero(params.n_samples, params.n_pheno);
  m_ests.Xt_Gamma_X_inv.resize(params.n_pheno);
  // set covariates for firth
  if(params.firth){
    firth_est.covs_firth = MatrixXd::Zero(params.n_samples, pheno_data.new_cov.cols() + 1);
    firth_est.covs_firth.leftCols(pheno_data.new_cov.cols()) << pheno_data.new_cov;
  }


  // start analyzing each chromosome
  int block = 0, chrom, chrom_nsnps, chrom_nb, bs;
  uint32_t snp_count = 0;

  for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {
    chrom = files.chr_read[itr];
    if( !in_map(chrom, chr_map) ) continue;

    chrom_nsnps = chr_map[chrom][0];
    chrom_nb = chr_map[chrom][1];
    if(chrom_nb == 0) continue;

    sout << "Chromosome " << chrom << " [" << chrom_nb << " blocks in total]\n";

    // read polygenic effect predictions from step 1
    blup_read_chr(chrom);

    // compute phenotype residual (adjusting for BLUP [and covariates for BTs])
    if(params.binary_mode) compute_res_bin(chrom);
    else compute_res();


    // analyze by blocks of SNPs
    for(int bb = 0; bb < chrom_nb ; bb++) {

      bs = ((bb +1) * params.block_size > chrom_nsnps) ? chrom_nsnps - (bb * params.block_size) : params.block_size;

      // resize elements
      if( (bb == 0) || ((bb +1) * params.block_size > chrom_nsnps) ) {
        Gblock.Gmat = MatrixXd::Zero(bs,params.n_samples);
        stats = MatrixXd::Zero(bs, params.n_pheno);
        Gblock.snp_afs = MatrixXd::Zero(bs, 1);
        Gblock.snp_mac = MatrixXd::Zero(bs, 1);
        if(params.dosage_mode) Gblock.snp_info = MatrixXd::Zero(bs, 1);
        if(params.binary_mode) sqrt_denum = MatrixXd::Zero(bs, params.n_pheno);
        else scaleG_pheno = MatrixXd::Zero(bs, params.n_pheno);
        if(params.use_SPA) {
          spa_est.SPA_pvals = MatrixXd::Zero(bs, params.n_pheno);
          Gblock.snp_flipped.resize(bs);
        }
        if(params.htp_out) {
          Gblock.genocounts.resize(bs);
          for( int j = 0; j < bs; ++j ) Gblock.genocounts[j] = MatrixXd::Zero(6, params.n_pheno);
        }
      }
      Gblock.bad_snps = ArrayXb::Constant(bs, false);

      // get genotype matrix for block (mean impute)
      get_G(block, bs, chrom, snp_count, snpinfo, &params, &files, &Gblock, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, sout);


      if(!params.binary_mode || params.firth_approx){
        // residualize and scale genotypes (and identify monomorphic if present)
        // only do it for firth approx. test with BTs (ests. are unchanged for other tests)
        residualize_genotypes();
      } else {
        scale_G = ArrayXd::Ones(bs); // no scaling is applied to the SNPs
        // ensure that snps which failed MAC filter are all polymorphic
        if( Gblock.bad_snps.cast<int>().sum() > 0 ) {
          for(int i = 0; i < bs ; i++){
            if(Gblock.bad_snps(i)) {
              Gblock.Gmat.row(i).head(10).array() += 1; // +1 to first 10 entries
              if(params.verbose) sout << "WARNING: Ignoring SNP with low variance.\n";
            }
          }
        }
      }
      n_ignored_snps += Gblock.bad_snps.cast<int>().sum();


      // perform assoc. testing
      t1 = std::chrono::high_resolution_clock::now();
      sout << "   -computing and writing association test statistics..." << flush;

      if(params.binary_mode) {

        for( int i = 0; i < params.n_pheno; ++i ) {

          // project out covariates from G
          WX = m_ests.Gamma_sqrt.col(i).asDiagonal() * pheno_data.new_cov;
          GW = Gblock.Gmat * m_ests.Gamma_sqrt.col(i).asDiagonal();
          Gblock.Gmat_tmp = GW - ((GW * WX) * m_ests.Xt_Gamma_X_inv[i]) * WX.transpose();
          Gblock.Gmat_tmp.array().rowwise() *= pheno_data.masked_indivs.col(i).transpose().array().cast<double>();
          denum_tstat = Gblock.Gmat_tmp.rowwise().squaredNorm();
          sqrt_denum.col(i) = denum_tstat.array().sqrt();

          // score test stat for BT
          stats.col(i).array() = (Gblock.Gmat_tmp * res.col(i)).array() / sqrt_denum.col(i).array();

          if(params.use_SPA) run_SPA_test(i);
        }

      } else {

        // score test stat for QT
        if( params.strict_mode )
          stats = (Gblock.Gmat * res) / sqrt( pheno_data.Neff(0) - params.ncov );
        else {
          // compute GtG for each phenotype (different missing patterns)
          for( int i = 0; i < params.n_pheno; ++i ) {
            scaleG_pheno.col(i) = ( Gblock.Gmat.array().rowwise() * pheno_data.masked_indivs.col(i).transpose().array().cast<double>()).square().matrix().rowwise().sum();
          }
          stats = ( (Gblock.Gmat * res).array() / scaleG_pheno.array().sqrt() ).matrix();
        }

      }

      // write stats to file
      for( int i = 0; i < bs; ++i ) {

        // don't print anything for monomorphic snps
        if(Gblock.bad_snps(i)) {
          snp_count++;
          continue;
        }
        if( params.dosage_mode ) info = Gblock.snp_info(i,0);

        if(!params.htp_out) tmpstr = print_sum_stats_head(snp_count, Gblock.snp_afs(i, 0), info, test_string);

        if(!params.split_by_pheno) ofile << tmpstr;

        for( int j = 0; j < params.n_pheno; ++j ) {
          if(params.split_by_pheno) {
            if(!params.htp_out) (*ofile_split[j]) << tmpstr;
            else  (*ofile_split[j]) <<  print_sum_stats_head_htp(snp_count, j, model_type);
          }

          chisq_val = stats(i,j) * stats(i,j);
          // test statistic & pvalue
          if(!params.use_SPA) {
            pval_log = check_pval(chisq_val, chrom, i, j);
          } else {
            pval_log = spa_est.SPA_pvals(i, j);
            pval_converged = (pval_log != params.missing_value_double);
          }

          // summary stats
          if( !params.binary_mode ){

            // estimate & SE for QT
            if( params.strict_mode )
              bhat = stats(i,j) * ( pheno_data.scale_Y(j) * p_sd_yres(j)) / ( sqrt(pheno_data.Neff(j) - params.ncov) * scale_G(i) );
            else
              bhat = stats(i,j) * ( pheno_data.scale_Y(j) * p_sd_yres(j)) / ( sqrt(scaleG_pheno(i,j)) * scale_G(i) );
            se_b = bhat / stats(i,j);

          } else {

            // with Firth, get sum. stats from Firth logistic regression
            if( params.firth && (chisq_val > chisq_thr) && pval_converged ){
              pval_raw = max(params.nl_dbl_dmin, pow(10, -pval_log)); // to prevent overflow
              chisq_val = quantile(complement(chisq, pval_raw));
              bhat = firth_est.bhat_firth;

              // compute SE from beta & pvalue
              if( params.back_correct_se && (chisq_val > 0) )
                se_b = fabs(bhat) / sqrt(chisq_val);
              else
                se_b = firth_est.se_b_firth;

            } else {
              se_b = 1 / sqrt_denum(i, j);
              // with SPA, calculate test stat based on SPA p-value
              if( params.use_SPA && (chisq_val > chisq_thr) && pval_converged ){
                pval_raw = max(params.nl_dbl_dmin, pow(10, -pval_log)); // to prevent overflow
                chisq_val = quantile(complement(chisq, pval_raw));
                bhat = sgn(stats(i,j)) * sqrt(chisq_val);
              } else bhat = stats(i,j);
              bhat *= se_b;
              if( params.use_SPA && Gblock.snp_flipped[i] ) bhat *= -1;
            }

            bhat /= scale_G(i);
            se_b /= scale_G(i);

          }

          if( !pval_converged ) n_failed_tests++;

          if(!params.split_by_pheno) ofile << print_sum_stats(bhat, se_b, chisq_val, pval_log, j, pval_converged);
          else if(!params.htp_out) (*ofile_split[j]) << print_sum_stats(bhat, se_b, chisq_val, pval_log, pval_converged);
          else (*ofile_split[j]) << print_sum_stats_htp(bhat, se_b, chisq_val, pval_log, Gblock.snp_afs(i, 0), info, Gblock.snp_mac(i, 0), Gblock.genocounts[i], zcrit, j, pval_converged);

          if(params.split_by_pheno) (*ofile_split[j]) << endl;
        }
        if(!params.split_by_pheno) ofile << endl;

        snp_count++;
      }

      sout << "done";
      t2 = std::chrono::high_resolution_clock::now();
      auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
      sout << " (" << duration1.count() << "ms) "<< endl;

      block++;
    }

  }

  sout << endl;
  if(!params.split_by_pheno){
    sout << "Association results stored in file : " << out << endl;
    ofile.closeFile();
  } else {
    sout << "Association results stored separately for each trait " << ( params.htp_out ? "(HTPv4 format) " : "" ) << "in files : \n";
    for( int j = 0; j < params.n_pheno; ++j ) {
      ofile_split[j]->closeFile();
      delete ofile_split[j];
      sout << "* [" << out_split[j] << "]\n";
    }
    sout << endl;
  }

  if(params.firth || params.use_SPA) {
    sout << "Number of tests with " << (params.firth ? "Firth " : "SPA ");
    sout << "correction : (" << n_corrected << "/" << (snpinfo.size() - n_skipped_snps - n_ignored_snps) * params.n_pheno << ")" <<  endl;
    sout << "Number of failed tests : (" << n_failed_tests << "/" << n_corrected << ")\n";
  }

  sout << "Number of ignored SNPs due to low MAC " << ( params.setMinINFO ? "or info score " : "");
  sout << ": " << n_ignored_snps << endl;

}
*/

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Functions needed in testing mode
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::setup_output(Files* ofile, string& out, std::vector<Files*>& ofile_split, std::vector< string >& out_split){

  string tmpstr, mask_header;
  if(params.build_mask) mask_header = build_mask_header();

  if(params.getCorMat){// header N,M
    out = files.out_file + ".corr";
    if(params.cor_out_txt){
      sout << " * computing correlation matrix\n  + output to text file ["<<out<<"]\n";
      sout << "  + n_snps = " << params.n_variants <<"\n\n";
      ofile->openForWrite(out, sout);
    } else {
      sout << " * computing correlation matrix (storing R^2 values)\n  + output to binary file ["<<out<<"]\n";
      sout << "  + n_snps = " << params.n_variants <<"\n\n";
      ofile->openBinMode(out, std::ios_base::out | std::ios_base::binary, sout);
      ArrayXi vals(2);
      vals << params.n_samples , params.n_variants;
      //cerr << vals << endl;
      ofile->writeBinMode(vals, sout);
    }
    return;
  }

  if(!params.split_by_pheno){ // single file

    out = files.out_file + ".regenie" + (params.gzOut ? ".gz" : "");
    ofile->openForWrite(out, sout);
    (*ofile) << (params.build_mask ? mask_header : "" ) << print_header_output();

  } else { // split results in separate files for each phenotype

    out_split.resize( params.n_pheno );
    ofile_split.resize( params.n_pheno );

    // header of output file
    if(!params.htp_out) tmpstr = print_header_output_single();
    else tmpstr = print_header_output_htp();

    for(int i = 0; i < params.n_pheno; i++) {
      out_split[i] = files.out_file + "_" + files.pheno_names[i] + ".regenie" + (params.gzOut ? ".gz" : "");
      ofile_split[i] = new Files;
      ofile_split[i]->openForWrite( out_split[i], sout );
      (*ofile_split[i]) << (params.build_mask ? mask_header : "" ) << tmpstr;
    }

  }

}

void Data::print_test_info(){

  if(params.getCorMat) return;

  if(params.write_masks) {
    bm.write_info(&params, &in_filters, sout);
    sout << " * user specified to write masks (in PLINK bed format)\n";
    if(params.dosage_mode) sout << "   +dosages will be converted to hardcalls\n";
    if(params.write_setlist) bm.prep_setlists(files.new_sets, files.out_file, sout);
  }

  sout << " * using minimum MAC of " << (params.build_mask ? params.min_MAC_mask : params.min_MAC) << 
    " (" << (params.build_mask ? "masks" : "variants") << " with lower MAC are ignored)\n";
  if(params.setMinINFO) sout << " * using minimum imputation info score of " << params.min_INFO << " (variants with lower info score are ignored)\n";
  if(params.firth || params.use_SPA) {
    sout << " * using " << (params.firth_approx ? "fast ": "") << (params.firth ? "Firth ": "SPA ");
    sout << "correction for logistic regression p-values less than " << params.alpha_pvalue << endl;
    if(params.back_correct_se) sout << "    - using back-correction to compute Firth SE\n";
    n_corrected = 0;
  }
  // if testing select chromosomes
  if( params.select_chrs ) sout << " * user specified to test only on select chromosomes\n";
  sout << endl;


  if(params.test_type == 0) test_string = "ADD";
  else if(params.test_type == 1) test_string = "DOM";
  else test_string = "REC";


  if(params.htp_out){
    if(params.binary_mode & params.firth) correction_type = "-FIRTH";
    else if(params.binary_mode & params.use_SPA) correction_type = "-SPA";
    else if(params.binary_mode) correction_type = "-LOG";
    else correction_type = "-LR";

    if(params.skip_blups) model_type = test_string + correction_type;
    else model_type = test_string + "-WGR" + correction_type;
  }

  if( params.joint_test ) {
    jt.get_test_info(&params, test_string, sout);
    jt.out_file_prefix = files.out_file;
  }

}

std::string Data::print_header_output(){

  int i;
  std::ostringstream buffer;

  buffer << "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ " << ( params.dosage_mode ? "INFO ":"") << "TEST ";
  for(i = 1; i < params.n_pheno; i++) buffer << "BETA.Y" << i << " SE.Y" << i << " CHISQ.Y" << i << " LOG10P.Y" << i << " ";
  // last phenotype
  buffer << "BETA.Y" << i << " SE.Y" << i <<  " CHISQ.Y" << i << " LOG10P.Y" << i <<  endl;

  return buffer.str();
}


std::string Data::print_header_output_single(){

  std::ostringstream buffer;

  buffer << "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ " << 
    ( !params.build_mask && params.dosage_mode ? "INFO ":"") << 
    "N TEST BETA SE CHISQ LOG10P";
  if(params.joint_test) buffer << " DF";
  buffer << endl;

  return buffer.str();
}

std::string Data::print_header_output_htp(){

  std::ostringstream buffer;

  buffer << "Name" << "\t" << "Chr" << "\t" << "Pos" << "\t" << "Ref" << "\t" << "Alt" << "\t" << "Trait" << "\t" << "Cohort" << "\t" << "Model" << "\t" << "Effect" << "\t" << "LCI_Effect" << "\t" << "UCI_Effect" << "\t" << "Pval" << "\t" << "AAF" << "\t" << "Num_Cases"<< "\t" << "Cases_Ref" << "\t" << "Cases_Het" << "\t" << "Cases_Alt" << "\t" << "Num_Controls" << "\t" << "Controls_Ref" << "\t" << "Controls_Het"<< "\t"<< "Controls_Alt" << "\t" << "Info\n";

  return buffer.str();
}

std::string Data::print_sum_stats_head(const int snp_count){

  std::ostringstream buffer;

  buffer << snpinfo[snp_count].chrom << " " << snpinfo[snp_count].physpos << " "<< snpinfo[snp_count].ID << " "<< snpinfo[snp_count].allele1 << " "<< snpinfo[snp_count].allele2 << " " ;

  return buffer.str();
}

std::string Data::print_sum_stats_head_htp(const int snp_count, const int ph, const string model){

  std::ostringstream buffer;

  buffer << snpinfo[snp_count].ID << "\t"<< snpinfo[snp_count].chrom << "\t" << snpinfo[snp_count].physpos << "\t"<< snpinfo[snp_count].allele1 << "\t"<< snpinfo[snp_count].allele2 << "\t" << files.pheno_names[ph] << "\t" << params.cohort_name << "\t" << model << "\t";

  return buffer.str();
}

// native format - all phenos
std::string Data::print_sum_stats(const double beta, const double se, const double chisq, const double pv, const int ph, const bool test_pass){

  std::ostringstream buffer;
  bool last = ((ph+1)==params.n_pheno);

  buffer << beta << ' ' << se << ' ' << chisq << ' ';
  if(test_pass) buffer << pv;
  else buffer << "NA";
  buffer << (last ? "" : " ");

  return buffer.str();
}

// native format - single pheno
std::string Data::print_sum_stats(const double af, const double info, const int n, const string model, const double beta, const double se, const double chisq, const double pv, const bool test_pass){

  std::ostringstream buffer;

  buffer << af << " " ;
  if(!params.build_mask && params.dosage_mode) buffer << info << " ";
  buffer << n << " " << model << " " << beta << ' ' << se << ' ' << chisq << ' ';
  if(test_pass) buffer << pv;
  else buffer << "NA";

  return buffer.str();
}


std::string Data::print_sum_stats(const double beta, const double se, const double chisq, const double pv, const bool test_pass){

  std::ostringstream buffer;

  buffer << beta << ' ' << se << ' ' << chisq << ' ';
  if(test_pass) buffer << pv;
  else buffer << "NA";

  return buffer.str();
}


std::string Data::print_sum_stats_htp(const double beta, const double se, const double chisq, const double pv, const double af, const double info, const double mac, const MatrixXd& genocounts, const double zcrit, const int ph, const bool test_pass){

  std::ostringstream buffer;
  double outp_val = -1, effect_val, outse_val;
  normal nd(0,1);

  if( test_pass ) {
    outp_val = max(params.nl_dbl_dmin, pow(10, - pv)); // to prevent overflow
    if(outp_val == 1) outp_val = 1 - 1e-7;
  } 

  // Effect / CI bounds / Pvalue columns
  if(!params.binary_mode || ( params.firth & (outp_val >= 0) ) ){

    if(!params.binary_mode) // QT
      buffer << beta << "\t" << (beta - zcrit * se) << "\t" << (beta + zcrit * se) << "\t" << outp_val << "\t";
    else // BT (on OR scale)
      buffer << exp(beta) << "\t" << exp(beta - zcrit * se) << "\t" << exp(beta + zcrit * se) << "\t" << outp_val << "\t";

  } else {
    // for tests that failed, print NA
    if(outp_val<0){
      buffer << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t";
    } else{
      // for spa or uncorrected logistic score test
      // compute allelic OR
      effect_val = (2*genocounts(3,ph)+genocounts(4,ph)+.5)*(2*genocounts(2,ph)+genocounts(1,ph)+.5)/(2*genocounts(5,ph)+genocounts(4,ph)+.5)/(2*genocounts(0,ph)+genocounts(1,ph)+.5);
      // compute SE = log(allelic OR) / zstat
      outse_val = fabs(log(effect_val)) / quantile(complement(nd, outp_val/2 ));

      buffer << effect_val << "\t" << effect_val * exp(- zcrit * outse_val) << "\t" << effect_val * exp(zcrit * outse_val) << "\t" << outp_val << "\t";
    }
  }

  // print out AF
  buffer << af << "\t";

  // print counts in cases
  buffer << (int) genocounts.block(0,ph,3,1).sum() << "\t" << (int) genocounts(0,ph) << "\t" << (int) genocounts(1,ph) << "\t" << (int) genocounts(2,ph) << "\t";

  // print counts in controls
  if(params.binary_mode){
    buffer << (int) genocounts.block(3,ph,3,1).sum() << "\t" << (int) genocounts(3,ph) << "\t" << (int) genocounts(4,ph) << "\t" << (int) genocounts(5,ph);
  } else buffer << "NA\tNA\tNA\tNA";

  // info column
  if(params.binary_mode){

    if(outp_val<0){ // only have NA
      buffer << "\t" << "REGENIE_BETA=NA;" << "REGENIE_SE=NA" << (params.firth ? "" : ";SE=NA");
    } else {
      buffer << "\t" << "REGENIE_BETA=" << beta << ";" << "REGENIE_SE=" << se;
      // SPA/uncorrected logistic => also print SE from allelic OR
      if(!params.firth) buffer << ";SE=" << outse_val;
    }

  } else buffer << "\t" << "REGENIE_SE=" << se; // fot QTs

  // info score
  if(!params.build_mask && params.dosage_mode) buffer << ";INFO=" << info;

  // mac
  buffer << ";MAC=" << mac;

  return buffer.str();
}

void Data::print_cor(Files* ofile){

  int bits = 16; // break [0,1] into 2^bits intervals
  double mult = (1ULL << bits) - 1; // map to 0,...,2^bits-1

  MatrixXd LDmat = (Gblock.Gmat.transpose() * Gblock.Gmat) / (params.n_samples - params.ncov);
  if(params.cor_out_txt){
    IOFormat Fmt(StreamPrecision, 0, " ", "\n", "", "");
    (*ofile) << LDmat;
    ofile->closeFile();
    exit_early();
  }

  ArrayXt vals;
  vals.resize( (Gblock.Gmat.cols() * (Gblock.Gmat.cols() - 1)) / 2 );

  for(int i = 0, k = 0; i < LDmat.rows(); i++){
    for(int j = i+1; j < LDmat.cols(); j++){
      vals(k++) = LDmat(i,j) * LDmat(i,j) * mult + 0.5; // round to nearest integer
    }
  }

  //cerr << LDmat.block(0,0,5,5) << "\n\n" << vals.head(5);

  ofile->writeBinMode(vals, sout);
  ofile->closeFile();

  exit_early();

}




/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    prep for association test
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::blup_read_chr(const int chrom) {
  string line, filename, tmp_pheno;
  std::vector< string > id_strings, tmp_str_vec ;
  double in_blup;
  uint32_t indiv_index;
  Files blupf;

  // skip reading if specified by user or if PRS is given (same for all chromosomes)
  if( params.use_prs || params.skip_blups ) return;

  m_ests.blups = MatrixXd::Zero(params.n_samples, params.n_pheno);

  sout << "   -reading loco predictions for the chromosome..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  // read blup file for each phenotype
  for(int ph = 0; ph < params.n_pheno; ph++) {

    int i_pheno = files.pheno_index[ph];
    ArrayXb read_indiv = ArrayXb::Constant(params.n_samples, false);
    blupf.openForRead(files.blup_files[ph], sout);

    // check header
    blupf.readLine(line);
    id_strings = string_split(line,"\t ");
    if( id_strings[0] != "FID_IID") {
      sout << "ERROR: Header of blup file must start with FID_IID.\n";
      exit(EXIT_FAILURE);
    }

    // skip to chr
    blupf.ignoreLines(chrom-1);

    blupf.readLine(line);
    tmp_str_vec = string_split(line,"\t ");

    // check number of entries is same as in header
    if(tmp_str_vec.size() != id_strings.size()) {
      sout << "ERROR: blup file for phenotype [" << files.pheno_names[i_pheno] << "] has different number of entries on line " << chrom + 1 << " compared to the header (=" << tmp_str_vec.size() << " vs " << id_strings.size() << ").\n";
      exit(EXIT_FAILURE);
    }

    // check starts with chromosome number
    if(chrStrToInt(tmp_str_vec[0], params.nChrom) != chrom) {
      sout << "ERROR: blup file for phenotype [" << files.pheno_names[i_pheno] << "] start with `" << tmp_str_vec[0]<< "`" <<
        "instead of chromosome number=" << chrom << ".\n";
      exit(EXIT_FAILURE);
    }

    // read blup data
    for( size_t filecol = 1; filecol < id_strings.size(); filecol++ ) {

      // ignore sample if it is not in genotype data
      if (!in_map(id_strings[filecol], params.FID_IID_to_ind)) continue;
      indiv_index = params.FID_IID_to_ind[id_strings[filecol]];

      // ignore sample if it is not included in analysis
      if(!in_filters.ind_in_analysis(indiv_index)) continue;

      // check if duplicate
      if( !read_indiv(indiv_index) ){
        read_indiv(indiv_index) = true;
      } else {
        sout << "ERROR: Individual appears more than once in blup file [" << files.blup_files[ph] <<"]: FID_IID=" << id_strings[filecol] << endl;
        exit(EXIT_FAILURE);
      }

      in_blup = convertDouble( tmp_str_vec[filecol], &params, sout);

      // if blup is NA then individual must be ignored in analysis for the phenotype (ie mask = 0)
      if (in_blup == params.missing_value_double){
        if(pheno_data.masked_indivs(indiv_index, i_pheno)){
          sout << "ERROR: Individual (FID_IID=" << id_strings[filecol] << ") has missing blup prediction at chromosome " << chrom <<" for phenotype " << files.pheno_names[i_pheno]<< ". ";
          sout << "Either set their phenotype to `NA`, specify to ignore them using option '--remove', or skip reading predictions with option '--ignore-pred'.\n" << params.err_help ;
          exit(EXIT_FAILURE);
        };
      } else m_ests.blups(indiv_index, i_pheno) = in_blup;
    }

    // force all non-masked samples to have loco predictions
    //   -> this should not occur as masking of absent samples is done in blup_read() function
    if( (pheno_data.masked_indivs.col(i_pheno).array() && read_indiv).cast<int>().sum() < pheno_data.masked_indivs.col(i_pheno).cast<int>().sum() ){
      sout << "ERROR: All samples included in the analysis (for phenotype " <<
        files.pheno_names[i_pheno]<< ") must have LOCO predictions in file : " << files.blup_files[ph] << "\n";
      exit(EXIT_FAILURE);
    }

    blupf.closeFile();
  }

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;

}

void Data::set_blocks_for_testing() {

  params.total_n_block = 0;
  int blocks_left = params.n_block;
  int nchr = 0;

  if(params.getCorMat) params.block_size = params.n_variants;

  map<int, vector<int> >::iterator itr;
  map<int, vector<int> > m1;
  for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) {
    int chrom_nsnps = itr->second[0];
    int nb = ceil((double) chrom_nsnps / params.block_size);

    if(params.n_block > 0) {
      if(blocks_left > 0) {
        int minb = min(nb, blocks_left);
        //sout << << endl;
        itr->second[1] = minb;
        params.total_n_block += minb;
        blocks_left -= minb;
      }
    } else {
      itr->second[1] = nb;
      params.total_n_block += nb;
    }
    if(params.getCorMat && (itr->second[1] > 0)) nchr++;
    m1.insert(pair<int, vector<int> >(itr->first, itr->second));
  }
  chr_map = m1;

  if(params.getCorMat && (nchr > 1 || params.n_variants < 2)){
    sout << "ERROR: can only compute LD matrix for a single chromosome (use --chr/--chrList/--range) and >=2 variants.\n"; exit(EXIT_FAILURE);
  }

  // summarize block sizes
  sout << left << std::setw(20) << " * # threads" << ": [" << params.threads << "]\n";
  sout << left << std::setw(20) << " * block size" << ": [" << params.block_size << "]\n";
  sout << left << std::setw(20) << " * # blocks" << ": [" << params.total_n_block << "]\n";

}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    functions used for assoc. test
/////////////////////////////////////////////////
/////////////////////////////////////////////////
double Data::check_pval(double tstat, int chrom, int snp,  int ph){

  chi_squared chisq(1);
  double chisq_thr = quantile(chisq, 1 - params.alpha_pvalue);
  double pval, logp;
  pval_converged = false;

  // if quantitative trait, or firth isn't used, or Tstat < threshold, no correction & return pvalue
  if(!params.binary_mode || !params.firth || tstat <= chisq_thr){
    pval = cdf(complement(chisq, tstat));

    // check if pvalue from boost 0 (when tstat is very large)
    // use log10p = log10( 2*f(tstat) ) with f=pdf of chisq(1)
    if(pval == 0) logp = log10(2) - 0.5 * log10( 2 * M_PI * tstat ) - 0.5 * tstat * M_LOG10E ;
    else logp = log10(pval);

    pval_converged = true;
    logp *= -1;
    return logp;
  }
  n_corrected++;

  // firth
  return run_firth_correction(chrom, snp, ph);

}

double Data::run_firth_correction(int chrom, int snp,  int ph){

  bool converged_l0, converged_l1;
  double dif_deviance, pval, logp;
  chi_squared chisq(1);

  // add tested SNP to set of covariates
  firth_est.covs_firth.rightCols(1) = Gblock.Gmat.row(snp).transpose();

  // obtain null deviance (set SNP effect to 0 and compute max. pen. LL)
  if(!params.firth_approx){
    converged_l0 = fit_firth_logistic(chrom, ph, true, &params, &pheno_data, &m_ests, &firth_est, sout);
    if(!converged_l0) return params.missing_value_double;
  }

  // fit full model and compute deviance
  converged_l1 = fit_firth_logistic(chrom, ph, false, &params, &pheno_data, &m_ests, &firth_est, sout);
  if(!converged_l1) return params.missing_value_double;

  // compute LR stat
  dif_deviance = firth_est.deviance_logistic;

  // in case of numerical errors
  if( dif_deviance < 0 )
    return params.missing_value_double;

  pval = cdf(complement(chisq, dif_deviance));
  if(pval == 0) logp = log10(2) - 0.5 * log10( 2 * M_PI * dif_deviance ) - 0.5 * dif_deviance * M_LOG10E ;
  else logp = log10(pval);
  logp *= -1;
  pval_converged = true;

  return logp;

}

void Data::run_SPA_test(int ph){

  int bs = Gblock.Gmat_tmp.rows(), index_j, nnz;
  double pval, logp, zstat, zstat_sq, score_num, tval, limK1_low, limK1_high, root_K1;
  chi_squared chisq(1);
  double chisq_thr = quantile(chisq, 1 - params.alpha_pvalue);
  ArrayXd Gmu;

  for(int snp = 0; snp < bs; ++snp) {

    // skip ignored snps
    if(Gblock.bad_snps(snp)){
      spa_est.SPA_pvals(snp, ph) = params.missing_value_double;
      continue;
    }

    // check test statistic relative to threshold
    zstat = stats(snp, ph);
    zstat_sq = zstat * zstat;
    if( zstat_sq <= chisq_thr){
      pval = cdf(complement(chisq, zstat_sq));
      if(pval == 0) logp = log10(2) - 0.5 * log10( 2 * M_PI * zstat_sq ) - 0.5 * zstat_sq * M_LOG10E ;
      else logp = log10(pval);
      logp *= -1;
      spa_est.SPA_pvals(snp, ph) = logp;
      continue;
    }
    n_corrected++;

    // start SPA - check how many non-zero entries there are
    nnz = 0;
    for( std::size_t j = 0; j < Gblock.non_zero_indices_G[snp].size(); ++j ) {
      index_j = Gblock.non_zero_indices_G[snp][j];
      if(!pheno_data.masked_indivs(index_j,ph)) continue;
      nnz++;
    }
    fastSPA = nnz < 0.5 * pheno_data.Neff(ph);

    // compute needed quantities
    spa_est.val_c = sqrt( denum_tstat(snp) );  // sqrt( G'WG )
    score_num = zstat * spa_est.val_c;
    spa_est.Gmod = Gblock.Gmat_tmp.row(snp).transpose().array() / m_ests.Gamma_sqrt.col(ph).array() * pheno_data.masked_indivs.col(ph).array().cast<double>();
    Gmu = spa_est.Gmod * m_ests.Y_hat_p.col(ph).array();
    spa_est.val_a = Gmu.sum();

    if(fastSPA){
      spa_est.val_b = denum_tstat(snp);
      spa_est.val_d = 0;
      for( std::size_t j = 0; j < Gblock.non_zero_indices_G[snp].size(); ++j ) {
        index_j = Gblock.non_zero_indices_G[snp][j];
        if(!pheno_data.masked_indivs(index_j,ph)) continue;
        spa_est.val_b -= Gblock.Gmat_tmp(snp, index_j) * Gblock.Gmat_tmp(snp, index_j);
        spa_est.val_d += Gmu(index_j);
      }
    }

    // check if K'(t)= s can be solved
    limK1_low = (spa_est.Gmod < 0).select(spa_est.Gmod, 0 ).sum() - spa_est.val_a ;
    limK1_high = (spa_est.Gmod > 0).select(spa_est.Gmod, 0 ).sum() - spa_est.val_a ;
    if( score_num < limK1_low || score_num > limK1_high ){
      if(params.verbose) sout << "WARNING: SPA failed (solution to K'(t)=s is infinite)";
      spa_est.SPA_pvals(snp, ph) = params.missing_value_double;
      continue;
    }

    // keep track of whether obs stat is positive
    spa_est.pos_score = zstat  > 0;
    tval = fabs(zstat);

    // solve K'(t)= tval using a mix of Newton-Raphson and bisection method
    root_K1 = solve_K1(tval, fastSPA, denum_tstat(snp), snp, ph, &params, &m_ests, &spa_est, &Gblock, pheno_data.masked_indivs.col(ph), sout);
    if( root_K1 == params.missing_value_double ){
      spa_est.SPA_pvals(snp, ph) = params.missing_value_double;
      continue;
    }

    // compute pvalue
    spa_est.SPA_pvals(snp, ph) = get_SPA_pvalue(root_K1, tval, fastSPA, denum_tstat(snp), snp, ph, &params, &m_ests, &spa_est, &Gblock, pheno_data.masked_indivs.col(ph), sout);
  }

}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Testing mode (multi-threaded with OpenMP)
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::test_snps_fast() {

  sout << "Association testing mode";

  std::chrono::high_resolution_clock::time_point t1, t2;
  normal nd(0,1);
  double zcrit = quantile(complement(nd, .025));
  string out, tmpstr;
  vector < string > out_split;
  // output files
  Files ofile;
  // use pointer to class since it contains non-copyable elements
  vector < Files* > ofile_split;

#if defined(_OPENMP)
  omp_set_num_threads(params.threads); // set threads in OpenMP
  sout << " with " << (params.streamBGEN? "fast " : "") << "multithreading using OpenMP";
#endif
  sout << endl;

  setNbThreads(params.threads);
  file_read_initialization(); // set up files for reading
  read_pheno_and_cov(&files, &params, &in_filters, &pheno_data, &m_ests, sout);   // read phenotype and covariate files
  prep_run(&files, &params, &pheno_data, &m_ests, sout); // check blup files and adjust for covariates
  set_blocks_for_testing();   // set number of blocks
  print_usage_info(&params, &files, sout);
  print_test_info();
  setup_output(&ofile, out, ofile_split, out_split); // result file


  // start analyzing each chromosome
  int block = 0, chrom, chrom_nsnps, chrom_nb, bs;
  tally snp_tally;

  vector< variant_block > block_info;
  m_ests.Y_hat_p = MatrixXd::Zero(params.n_samples, params.n_pheno);
  m_ests.Gamma_sqrt = MatrixXd::Zero(params.n_samples, params.n_pheno);
  m_ests.Xt_Gamma_X_inv.resize(params.n_pheno);
  if(params.firth_approx){ // set covariates for firth
    firth_est.covs_firth = MatrixXd::Zero(params.n_samples, pheno_data.new_cov.cols() + 1);
    firth_est.covs_firth.leftCols(pheno_data.new_cov.cols()) << pheno_data.new_cov;
  }


  for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {

    chrom = files.chr_read[itr];
    if( !in_map(chrom, chr_map) ) continue;

    chrom_nsnps = chr_map[chrom][0];
    chrom_nb = chr_map[chrom][1];
    if(chrom_nb == 0) continue;

    sout << "Chromosome " << chrom << " [" << chrom_nb << " blocks in total]\n";

    if(!params.getCorMat){
      // read polygenic effect predictions from step 1
      blup_read_chr(chrom);

      // compute phenotype residual (adjusting for BLUP [and covariates for BTs])
      if(params.binary_mode) compute_res_bin(chrom);
      else compute_res();
    }

    // analyze by blocks of SNPs
    for(int bb = 0; bb < chrom_nb ; bb++) {

      sout << " block [" << block + 1 << "] : " << flush;

      bs = params.block_size;
      if( ((bb + 1) * params.block_size) > chrom_nsnps)
        bs = chrom_nsnps - (bb * params.block_size);

      Gblock.Gmat.resize(params.n_samples, bs);
      block_info.resize(bs);

      // read SNP, impute missing & compute association test statistic
      analyze_block(chrom, bs, &snp_tally, block_info);

      if(params.getCorMat) print_cor(&ofile);

      // print the results
      for(int isnp = 0; isnp < bs; isnp++) {
        uint32_t snpindex = snp_tally.snp_count + isnp;

        if( block_info[isnp].ignored ) {
          snp_tally.n_ignored_snps++;
          continue;
        }

        if(!params.htp_out) tmpstr = print_sum_stats_head(snpindex);

        if(!params.split_by_pheno) {
          ofile << tmpstr << block_info[isnp].af1 << " "; 
          if(!params.build_mask && params.dosage_mode) ofile << block_info[isnp].info1 << " ";
          ofile << "NA " << test_string;
        }

        for(int j = 0; j < params.n_pheno; ++j) {

          if( block_info[isnp].ignored_trait(j) ) {
            snp_tally.n_ignored_tests++;
            continue;
          }

          if( (params.firth || params.use_SPA) && block_info[isnp].is_corrected[j] ) n_corrected++;
          if( block_info[isnp].test_fail[j] ) snp_tally.n_failed_tests++;

          if(params.split_by_pheno) {
            if(!params.htp_out) (*ofile_split[j]) << tmpstr;
            else  (*ofile_split[j]) <<  print_sum_stats_head_htp(snpindex, j, model_type);
          }


          if(!params.split_by_pheno) ofile << print_sum_stats(block_info[isnp].bhat(j), block_info[isnp].se_b(j), block_info[isnp].chisq_val(j), block_info[isnp].pval_log(j), j, !block_info[isnp].test_fail[j]);
          else if(!params.htp_out) (*ofile_split[j]) << print_sum_stats(block_info[isnp].af(j), block_info[isnp].info(j), block_info[isnp].ns(j),test_string, block_info[isnp].bhat(j), block_info[isnp].se_b(j), block_info[isnp].chisq_val(j), block_info[isnp].pval_log(j), !block_info[isnp].test_fail[j]);
          else (*ofile_split[j]) << print_sum_stats_htp(block_info[isnp].bhat(j), block_info[isnp].se_b(j), block_info[isnp].chisq_val(j), block_info[isnp].pval_log(j), block_info[isnp].af(j), block_info[isnp].info(j), block_info[isnp].mac(j), block_info[isnp].genocounts, zcrit, j, !block_info[isnp].test_fail[j]);

          if(params.split_by_pheno) (*ofile_split[j]) << endl;
        }

        if(!params.split_by_pheno) ofile << endl;

      }

      snp_tally.snp_count += bs;
      block++;
    }

  }

  sout << endl;

  if(!params.split_by_pheno){
    sout << "Association results stored in file : " << out << endl;
    ofile.closeFile();
  } else {
    sout << "Association results stored separately for each trait " << ( params.htp_out ? "(HTPv4 format) " : "" ) << "in files : \n";
    for( int j = 0; j < params.n_pheno; ++j ) {
      ofile_split[j]->closeFile();
      delete ofile_split[j];
      sout << "* [" << out_split[j] << "]\n";
    }
    sout << endl;
  }

  if(params.firth || params.use_SPA) {
    sout << "Number of tests with " << (params.firth ? "Firth " : "SPA ");
    sout << "correction : (" << n_corrected << "/" << (snpinfo.size() - snp_tally.n_skipped_snps - snp_tally.n_ignored_snps) * params.n_pheno << ")" <<  endl;
    sout << "Number of failed tests : (" << snp_tally.n_failed_tests << "/" << n_corrected << ")\n";
  }

  sout << "Number of ignored tests due to low MAC ";
  if( params.setMinINFO ) sout << "or info score ";
  sout << ": " << snp_tally.n_ignored_snps * params.n_pheno + snp_tally.n_ignored_tests << endl;

}

// test SNPs in block
void Data::analyze_block(const int &chrom, const int &n_snps, tally* snp_tally, vector<variant_block> &all_snps_info){
  auto t1 = std::chrono::high_resolution_clock::now();

  const int start = snp_tally->snp_count;
  vector< vector < uchar > > snp_data_blocks;
  vector< uint32_t > insize, outsize;

  if((params.file_type == "bgen") && params.streamBGEN){
    uint64 pos_skip;
    ifstream bfile;
    snp_data_blocks.resize( n_snps );
    insize.resize(n_snps); outsize.resize(n_snps);
    bfile.open( files.bgen_file, ios::in | ios::binary );

    for(int isnp = 0; isnp < n_snps; isnp++) {

      // extract genotype data blocks single-threaded
      pos_skip = snpinfo[start + isnp].offset;
      bfile.seekg( pos_skip, ios_base::beg );

      readChunkFromBGEN(&bfile, &size1, &size2, &(snp_data_blocks[isnp]));
      insize[isnp] = size1;
      outsize[isnp] = size2;

    }
    bfile.close();

  } else if((params.file_type == "bgen") && !params.streamBGEN) readChunkFromBGENFileToG(n_snps, chrom, start, snpinfo, &params, &Gblock, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, all_snps_info, sout);
  else if(params.file_type == "pgen") readChunkFromPGENFileToG(start, n_snps, chrom, &params, &in_filters, &Gblock, pheno_data.masked_indivs, pheno_data.phenotypes_raw, snpinfo, all_snps_info);
  else if(params.file_type == "bed"){
    // read in N/4 bytes from bed file for each snp
    snp_data_blocks.resize( n_snps );
    for(int isnp = 0; isnp < n_snps;) {

      jumpto_bed( snpinfo[start + isnp].offset, &files );
      snp_data_blocks[isnp].resize(files.bed_block_size);
      files.bed_ifstream.read( reinterpret_cast<char *> (&snp_data_blocks[isnp][0]), files.bed_block_size);

      isnp++;
      snp_index_counter++;
    }
  }


  // start openmp for loop
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(int isnp = 0; isnp < n_snps; isnp++) {
    uint32_t snp_index = start + isnp;
    chi_squared chisq(1);
    double pval_raw;

    // to store variant information
    variant_block* block_info = &(all_snps_info[isnp]);

    if((params.file_type == "bgen") && params.streamBGEN){ // uncompress and extract the dosages
      parseSnpfromBGEN(isnp, chrom, &(snp_data_blocks[isnp]), insize[isnp], outsize[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &snpinfo[snp_index], &Gblock, block_info, sout);
    } else if(params.file_type == "bed"){ // extract hardcalls
      parseSnpfromBed(isnp, chrom, snp_data_blocks[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &Gblock, block_info);
    }

    // for QTs (or BTs with firth approx): project out covariates & scale
    residualize_geno(isnp, block_info);

    // skip SNP if fails filters
    if( block_info->ignored || params.getCorMat ) continue;
    
    block_info->pval_log = ArrayXd::Zero(params.n_pheno);
    block_info->bhat = ArrayXd::Zero(params.n_pheno);
    block_info->se_b = ArrayXd::Zero(params.n_pheno);
    block_info->test_fail.assign(params.n_pheno, false);
    block_info->is_corrected.assign(params.n_pheno, true);
    MapArXd Geno (Gblock.Gmat.col(isnp).data(), params.n_samples, 1);


    if(params.binary_mode) {
      MatrixXd tmpG, WX, GW;
      block_info->stats = ArrayXd::Zero(params.n_pheno);
      block_info->denum = ArrayXd::Zero(params.n_pheno);

      for( int i = 0; i < params.n_pheno; ++i ) {

        if( block_info->ignored_trait(i) ) 
          continue;

        // project out covariates from G
        WX = m_ests.Gamma_sqrt.col(i).asDiagonal() * pheno_data.new_cov;
        // apply mask
        GW = (Geno * m_ests.Gamma_sqrt.col(i).array() * pheno_data.masked_indivs.col(i).array().cast<double>()).matrix();
        tmpG = GW - WX * (m_ests.Xt_Gamma_X_inv[i] * (WX.transpose() * GW));
        block_info->denum(i) = tmpG.squaredNorm();

        // score test stat for BT
        block_info->stats(i) = (tmpG.array() * res.col(i).array()).sum() / sqrt( block_info->denum(i) );

        // check this
        if(params.use_SPA) run_SPA_test_snp(block_info, i, tmpG);
      }

    } else {

      // score test stat for QT
      if( params.strict_mode ) {
        Geno *= in_filters.ind_in_analysis.cast<double>();
        block_info->stats = (res.array().colwise() * Geno).matrix().transpose().rowwise().sum() / sqrt( pheno_data.Neff(0) - params.ncov );
      } else {
        // compute GtG for each phenotype (different missing patterns)
        block_info->scale_fac_pheno = (pheno_data.masked_indivs.cast<double>().array().colwise() * Geno).matrix().colwise().squaredNorm();
        block_info->stats = (Geno.matrix().transpose() * res).transpose().array() / block_info->scale_fac_pheno.sqrt();
      }

    }

    block_info->chisq_val = block_info->stats.square();

    for( int i = 0; i < params.n_pheno; ++i ) {

      if( block_info->ignored_trait(i) ) 
        continue;

      // test statistic & pvalue
      if(!params.use_SPA) check_pval_snp(block_info, chrom, i, isnp);

      // summary stats
      if( !params.binary_mode ){
        // estimate & SE for QT
        if( params.strict_mode )
          block_info->bhat(i) = block_info->stats(i) * ( pheno_data.scale_Y(i) * p_sd_yres(i)) / ( sqrt(pheno_data.Neff(i) - params.ncov) * block_info->scale_fac );
        else
          block_info->bhat(i) = block_info->stats(i) * ( pheno_data.scale_Y(i) * p_sd_yres(i) ) / ( sqrt(block_info->scale_fac_pheno(i)) * block_info->scale_fac );
        block_info->se_b(i) = block_info->bhat(i) / block_info->stats(i);
      } else {
        // with Firth, get sum. stats from Firth logistic regression
        if( params.firth && block_info->is_corrected[i] && !block_info->test_fail[i] ){
          pval_raw = max(params.nl_dbl_dmin, pow(10, - block_info->pval_log(i))); // to prevent overflow
          block_info->chisq_val(i) = quantile(complement(chisq, pval_raw));

          // compute SE from beta & pvalue
          if( params.back_correct_se && (block_info->chisq_val(i) > 0) )
            block_info->se_b(i) = fabs(block_info->bhat(i)) / sqrt(block_info->chisq_val(i));

        } else {
          block_info->se_b(i) = 1 / sqrt(block_info->denum(i));
          // with SPA, calculate test stat based on SPA p-value
          if( params.use_SPA && block_info->is_corrected[i] && !block_info->test_fail[i] ){
            pval_raw = max(params.nl_dbl_dmin, pow(10, - block_info->pval_log(i))); // to prevent overflow
            block_info->chisq_val(i) = quantile(complement(chisq, pval_raw));
            block_info->bhat(i) = sgn(block_info->stats(i)) * sqrt(block_info->chisq_val(i));
          } else block_info->bhat(i) = block_info->stats(i);
          block_info->bhat(i) *= block_info->se_b(i);
          if( params.use_SPA && block_info->flipped ) block_info->bhat(i) *= -1;
        }
        block_info->bhat(i) /= block_info->scale_fac;
        block_info->se_b(i) /= block_info->scale_fac;
      }

    }

    //if( (isnp==0) || (isnp == (n_snps-1)) ) cout << "G"<<isnp+1<<" MAF = " <<  block_info.MAF << endl;
  }

#if defined(_OPENMP)
  setNbThreads(params.threads);
#endif

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << "done (" << duration.count() << "ms) "<< endl;
}


void Data::residualize_geno(int isnp, variant_block* snp_data){

  if(snp_data->ignored) return;

  if(!params.binary_mode || params.firth_approx){
    // project out covariates
    MatrixXd beta = pheno_data.new_cov.transpose() * Gblock.Gmat.col(isnp);
    Gblock.Gmat.col(isnp) -= pheno_data.new_cov * beta;

    // scale
    snp_data->scale_fac = Gblock.Gmat.col(isnp).norm();
    snp_data->scale_fac /= sqrt( (params.strict_mode ? pheno_data.Neff(0) : in_filters.ind_in_analysis.cast<float>().sum()) - params.ncov );

    if( snp_data->scale_fac < params.numtol ) {
      snp_data->ignored = true;
      return;
    }
    Gblock.Gmat.col(isnp).array() /= snp_data->scale_fac;

  } else snp_data->scale_fac = 1;

}

void Data::compute_res(){

  res = pheno_data.phenotypes - m_ests.blups;
  res.array() *= pheno_data.masked_indivs.array().cast<double>();

  p_sd_yres = res.colwise().norm();
  p_sd_yres.array() /= sqrt(pheno_data.Neff - params.ncov);
  res.array().rowwise() /= p_sd_yres.array();

}

void Data::compute_res_bin(int chrom){

  fit_null_logistic(chrom, &params, &pheno_data, &m_ests, sout); // for all phenotypes

  res = pheno_data.phenotypes_raw - m_ests.Y_hat_p;
  res.array() /= m_ests.Gamma_sqrt.array();
  res.array() *= pheno_data.masked_indivs.array().cast<double>();

  // if using firth approximation, fit null penalized model with only covariates and store the estimates (to be used as offset when computing LRT in full model)
  if(params.firth_approx) fit_null_firth(chrom, &firth_est, &pheno_data, &m_ests, &files, &params, sout);

}

void Data::check_pval_snp(variant_block* block_info, int chrom, int ph, int isnp){

  chi_squared chisq(1);
  double chisq_thr = quantile(chisq, 1 - params.alpha_pvalue);
  double pval, logp, LRstat;
  double tstat = block_info->chisq_val(ph);

  // if quantitative trait, or firth isn't used, or Tstat < threshold, no correction done
  if(!params.binary_mode || !params.firth || tstat <= chisq_thr){
    pval = cdf(complement(chisq, tstat));

    // check if pvalue from boost 0 (when tstat is very large)
    if(pval == 0) logp = log10(2) - 0.5 * log10( 2 * M_PI * tstat ) - 0.5 * tstat * M_LOG10E ;
    else logp = log10(pval);

    logp *= -1;
    block_info->pval_log(ph) = logp;
    block_info->is_corrected[ph] = false;
    return;
  }

  // firth
  run_firth_correction_snp(chrom, ph, isnp, block_info);
  if(block_info->test_fail[ph]) return;

  LRstat = block_info->dif_deviance;
  pval = cdf(complement(chisq, LRstat));
  if(pval == 0) logp = log10(2) - 0.5 * log10( 2 * M_PI * LRstat ) - 0.5 * LRstat * M_LOG10E ;
  else logp = log10(pval);
  logp *= -1;

  block_info->pval_log(ph) = logp;
  return;
}


void Data::run_firth_correction_snp(int chrom, int ph, int isnp, variant_block* block_info){

  // obtain null deviance (set SNP effect to 0 and compute max. pen. LL)
  if(!params.firth_approx){
    fit_firth_logistic_snp(chrom, ph, isnp, true, &params, &pheno_data, &m_ests, &firth_est, &Gblock, block_info, sout);
    if(block_info->test_fail[ph]) return ;
  }

  // fit full model and compute deviance
  fit_firth_logistic_snp(chrom, ph, isnp, false, &params, &pheno_data, &m_ests, &firth_est, &Gblock, block_info, sout);

}


void Data::run_SPA_test_snp(variant_block* block_info, int ph, const VectorXd& Gtmp){

  int index_j;
  double pval, logp, zstat, zstat_sq, score_num, tval, limK1_low, limK1_high, root_K1;
  chi_squared chisq(1);
  double chisq_thr = quantile(chisq, 1 - params.alpha_pvalue);
  ArrayXd Gmu;

  // check test statistic relative to threshold
  zstat = block_info->stats(ph);
  zstat_sq = zstat * zstat;
  if( zstat_sq <= chisq_thr){
    pval = cdf(complement(chisq, zstat_sq));
    if(pval == 0) logp = log10(2) - 0.5 * log10( 2 * M_PI * zstat_sq ) - 0.5 * zstat_sq * M_LOG10E ;
    else logp = log10(pval);
    logp *= -1;
    block_info->pval_log(ph) = logp;
    block_info->is_corrected[ph] = false;
    return;
  }

  // compute needed quantities
  block_info->val_c = sqrt( block_info->denum(ph) );  // sqrt( G'WG )
  score_num = zstat * block_info->val_c;
  block_info->Gmod = Gtmp.array() / m_ests.Gamma_sqrt.col(ph).array() * pheno_data.masked_indivs.col(ph).array().cast<double>();
  Gmu = block_info->Gmod * m_ests.Y_hat_p.col(ph).array();
  block_info->val_a = Gmu.sum();

  if(block_info->fastSPA){
    block_info->val_b = block_info->denum(ph);
    block_info->val_d = 0;
    for( std::size_t j = 0; j < block_info->n_non_zero; ++j ) {
      index_j = block_info->non_zero_indices[j];
      if(!pheno_data.masked_indivs(index_j,ph)) continue;
      block_info->val_b -= Gtmp(index_j) * Gtmp(index_j);
      block_info->val_d += Gmu(index_j);
    }
  }

  // check if K'(t)= s can be solved
  limK1_low = (block_info->Gmod < 0).select(block_info->Gmod, 0 ).sum() - block_info->val_a ;
  limK1_high = (block_info->Gmod > 0).select(block_info->Gmod, 0 ).sum() - block_info->val_a ;
  if( score_num < limK1_low || score_num > limK1_high ){
    if(params.verbose) sout << "WARNING: SPA failed (solution to K'(t)=s is infinite)";
    block_info->test_fail[ph] = true;
    return;
  }

  // keep track of whether obs stat is positive
  block_info->pos_score = zstat  > 0;
  tval = fabs(zstat);

  // solve K'(t)= tval using a mix of Newton-Raphson and bisection method
  root_K1 = solve_K1_snp(tval, ph, &params, &m_ests, block_info, pheno_data.masked_indivs.col(ph), sout);
  if( root_K1 == params.missing_value_double ){
    block_info->test_fail[ph] = true;
    return;
  }

  // compute pvalue
  block_info->pval_log(ph) = get_SPA_pvalue_snp(root_K1, tval, ph, &params, &m_ests, block_info, pheno_data.masked_indivs.col(ph), sout);

}




/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Testing mode (joint tests)
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::test_joint() {

  
  sout << "Association testing mode (joint tests)";

  std::chrono::high_resolution_clock::time_point t1, t2;
  normal nd(0,1);
  double zcrit = quantile(complement(nd, .025));
  string out, tmpstr;
  vector < string > out_split;
  // output files
  Files ofile;
  // use pointer to class since it contains non-copyable elements
  vector < Files* > ofile_split;

  // set some parameters
  params.split_by_pheno = true; // always split files
  if( params.build_mask ){
    params.min_MAC_mask = params.min_MAC; // for association tests
    params.min_MAC = 0.5; // set this so can retain singletons (0.5 for dosages)
    bm.take_max = params.mask_rule_max;
    bm.take_comphet = params.mask_rule_comphet;
    if(!bm.take_max && !bm.take_comphet) params.htp_out = false; // due to genocounts with sum rule
    if(params.write_masks) bm.gfile_prefix = files.out_file + "_masks";
  }

#if defined(_OPENMP)
  omp_set_num_threads(params.threads); // set threads in OpenMP
  sout << " with " << (params.streamBGEN? "fast " : "") << "multithreading using OpenMP";
#endif
  sout << endl;

  setNbThreads(params.threads);
  file_read_initialization(); // set up files for reading
  read_pheno_and_cov(&files, &params, &in_filters, &pheno_data, &m_ests, sout);   // read phenotype and covariate files
  prep_run(&files, &params, &pheno_data, &m_ests, sout); // check blup files and adjust for covariates
  set_groups_for_testing();   // set groups of snps to test jointly
  print_usage_info(&params, &files, sout);
  print_test_info();
  if(!params.skip_test) setup_output(&ofile, out, ofile_split, out_split); // result file


  // start analyzing each chromosome
  int block = 0, chrom, chrom_nb, bs;
  tally snp_tally;

  m_ests.Y_hat_p = MatrixXd::Zero(params.n_samples, params.n_pheno);
  m_ests.Gamma_sqrt = MatrixXd::Zero(params.n_samples, params.n_pheno);
  m_ests.Xt_Gamma_X_inv.resize(params.n_pheno);
  if(params.firth_approx){ // set covariates for firth
    firth_est.covs_firth = MatrixXd::Zero(params.n_samples, pheno_data.new_cov.cols() + 1);
    firth_est.covs_firth.leftCols(pheno_data.new_cov.cols()) << pheno_data.new_cov;
  }
  if(params.joint_test) jt.scale_denum = (params.strict_mode ? pheno_data.Neff(0) : in_filters.ind_in_analysis.cast<float>().sum()) - jt.ncovars; // for gates


  for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {

    chrom = files.chr_read[itr];
    if( !in_map(chrom, chr_map) ) continue;

    chrom_nb = chr_map[chrom][1];

    // if no sets in chromosome, skip
    if(chrom_nb == 0)  continue;

    sout << "Chromosome " << chrom << " [" << chrom_nb << " sets in total]\n";

    // read polygenic effect predictions from step 1
    blup_read_chr(chrom);

    // compute phenotype residual (adjusting for BLUP [and covariates for BTs])
    if(params.binary_mode) compute_res_bin(chrom);
    else compute_res();


    // analyze by blocks of SNPs
    for(int bb = 0; bb < chrom_nb ; bb++) {

      sout << " set [" << block + 1 << "] : " << flush;

      bs = jt.setinfo[chrom - 1][bb].snp_indices.size();
      sout << bs << " variants..." << flush;

      vector< variant_block > block_info;
      block_info.resize(bs);
      if(params.joint_test && !params.build_mask) Gblock.Gmat.resize(params.n_samples, bs);

      // compute single snp association test statistic
      get_sum_stats(chrom, bb, block_info);

      // update number of variants (if masks were built)
      bs = jt.setinfo[chrom - 1][bb].snp_indices.size();
      jt.nvars = bs;

      if(params.skip_test){ // not performing the assoc tests
        snp_tally.snp_count += bs;
        block++;
        continue;
      }

      // print the results
      for(int isnp = 0; isnp < bs; isnp++) {
        uint32_t snpindex = jt.setinfo[chrom - 1][bb].snp_indices[ isnp ];

        if( block_info[isnp].ignored ) {
          snp_tally.n_ignored_snps++;
          jt.nvars--;
          continue;
        }

        if(!params.p_joint_only ){
          if(!params.htp_out) tmpstr = print_sum_stats_head(snpindex);
          if(!params.split_by_pheno) {
            ofile << tmpstr << block_info[isnp].af1 << " "; 
            if(!params.build_mask && params.dosage_mode) ofile << block_info[isnp].info1 << " ";
            ofile << "NA " << test_string;
          }
        }

        for(int j = 0; j < params.n_pheno; ++j) {

          if( block_info[isnp].ignored_trait(j) ) {
            snp_tally.n_ignored_tests++;
            continue;
          }

          if( (params.firth || params.use_SPA) && block_info[isnp].is_corrected[j] ) n_corrected++;
          if( block_info[isnp].test_fail[j] ) snp_tally.n_failed_tests++;

          if(!params.p_joint_only){

            if(params.split_by_pheno) {
              if(!params.htp_out) (*ofile_split[j]) << tmpstr;
              else  (*ofile_split[j]) <<  print_sum_stats_head_htp(snpindex, j, model_type);
            }

            if(!params.split_by_pheno) ofile << print_sum_stats(block_info[isnp].bhat(j), block_info[isnp].se_b(j), block_info[isnp].chisq_val(j), block_info[isnp].pval_log(j), j, !block_info[isnp].test_fail[j]) << (params.joint_test? " 1" : ""); // DF
            else if(!params.htp_out) (*ofile_split[j]) << print_sum_stats(block_info[isnp].af(j), block_info[isnp].info(j), block_info[isnp].ns(j), test_string, block_info[isnp].bhat(j), block_info[isnp].se_b(j), block_info[isnp].chisq_val(j), block_info[isnp].pval_log(j), !block_info[isnp].test_fail[j]) << (params.joint_test? " 1" : ""); // DF
            else (*ofile_split[j]) << print_sum_stats_htp(block_info[isnp].bhat(j), block_info[isnp].se_b(j), block_info[isnp].chisq_val(j), block_info[isnp].pval_log(j), block_info[isnp].af(j), block_info[isnp].info(j), block_info[isnp].mac(j), block_info[isnp].genocounts, zcrit, j, !block_info[isnp].test_fail[j]) << (params.joint_test? ";DF=1" : "");

            if(params.split_by_pheno) (*ofile_split[j]) << endl;
          }
        }

        if(!params.split_by_pheno && !params.p_joint_only) ofile << endl;

      }

      if( params.joint_test ){

        // compute and print set-based test result
        t1 = std::chrono::high_resolution_clock::now();
        sout << "     -computing joint association tests..." << flush;

        jt.get_variant_names(chrom, bb, snpinfo);
        for(int j = 0; j < params.n_pheno; ++j) 
          (*ofile_split[j]) << jt.apply_joint_test(chrom, bb, j, &pheno_data, res.col(j), &Gblock, block_info, files.pheno_names[j], &params);

        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        sout << "done (" << duration.count() << "ms) "<< endl;
      }

      snp_tally.snp_count += bs;
      block++;
    }


  }

  sout << endl;

  if(!params.skip_test) {
    if(!params.split_by_pheno){
      sout << "Association results stored in file : " << out << endl;
      ofile.closeFile();
    } else {
      sout << "Association results stored separately for each trait " << ( params.htp_out ? "(HTPv4 format) " : "" ) << "in files : \n";
      for( int j = 0; j < params.n_pheno; ++j ) {
        ofile_split[j]->closeFile();
        delete ofile_split[j];
        sout << "* [" << out_split[j] << "]\n";
      }
      sout << endl;
    }

    if(params.firth || params.use_SPA) {
      sout << "Number of tests with " << (params.firth ? "Firth " : "SPA ");
      sout << "correction : " << n_corrected << endl;
      sout << "Number of failed tests : (" << snp_tally.n_failed_tests << "/" << n_corrected << ")\n";
    }
  }

  sout << "Number of ignored tests due to low MAC ";
  if( params.setMinINFO ) sout << "or info score ";
  sout << ": " << snp_tally.n_ignored_snps * params.n_pheno + snp_tally.n_ignored_tests << endl;

  if(params.write_masks){
    bm.closeFiles();
    sout << "\nMasks written to : [" << files.out_file << "_masks.{bed,bim,fam}]\n";
  }

}



void Data::set_groups_for_testing() {

  int blocks_left = params.n_block;
  params.total_n_block = 0;

  // annotate variants by categories
  if(params.build_mask) {
    get_masks_info(&files, &params, &in_filters, bm.annotations, bm.regions, bm.masks, bm.mask_out, bm.all_masks, snpinfo, sout);
    bm.setBins(&params, sout);
  }

  // read list of variant sets to use for joint test
  read_setlist(&files, &params, &in_filters, jt.setinfo, snpinfo, bm.all_masks, bm.max_aaf, sout);

  // delete snpID map
  in_filters.snpID_to_ind.clear();

  // for each chromosome, count number of variant sets
  map<int, vector<int> >::iterator itr;
  map<int, vector<int> > m1;
  for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) {
    int chrom = itr->first;
    int nb = jt.setinfo[ chrom - 1].size();

    if(params.n_block > 0) {
      if(blocks_left > 0) {
        int minb = min(nb, blocks_left);
        //sout << << endl;
        itr->second[1] = minb;
        params.total_n_block += minb;
        blocks_left -= minb;
      }
    } else {
      itr->second[1] = nb;
      params.total_n_block += nb;
    }
    m1.insert(pair<int, vector<int> >(itr->first, itr->second));
  }
  chr_map = m1;


  // summarize block sizes
  sout << left << std::setw(20) << " * # threads" << ": [" << params.threads << "]\n";
  sout << left << std::setw(20) << " * # tested sets" << ": [" << params.total_n_block << "]\n";
  sout << left << std::setw(20) << " * max block size" << ": [" << params.block_size << "]\n";

  if(params.build_mask) sout << " * rule used to build masks : " << params.mask_rule << endl;

}

std::string Data::build_mask_header(){

  std::ostringstream buffer;

  // header = ##MASKS=<Mask1="X,X";Mask2="X,X";...;MaskK="X,X">
  buffer << "##MASKS=<";
  for(size_t i = 0; i < bm.mask_out.size(); i++)
    buffer << bm.mask_out[i][0] << "=\"" << bm.mask_out[i][1] << "\"" << ((i+1) < bm.mask_out.size() ? ";" : "");

  buffer << ">\n";

  return buffer.str();
}


// test SNPs in block
void Data::get_sum_stats(const int chrom, const int varset, vector<variant_block>& all_snps_info){

  auto t1 = std::chrono::high_resolution_clock::now();

  int n_snps = jt.setinfo[chrom - 1][varset].snp_indices.size();
  vector< vector < uchar > > snp_data_blocks;
  vector< uint32_t > insize, outsize;


  // read in markers and if applicable build masks
  if(!params.build_mask) readChunk(chrom, varset, 0, n_snps, snp_data_blocks, insize, outsize, all_snps_info);
  else {
    getMask(chrom, varset, snp_data_blocks, insize, outsize, all_snps_info);
    n_snps = jt.setinfo[chrom - 1][varset].snp_indices.size();
    //cerr << "M=" << n_snps << endl;
    if(params.skip_test || (n_snps == 0)) return;

    // starting association testing with built masks
    t1 = std::chrono::high_resolution_clock::now();
    sout << "     -computing association tests..." << flush;
  }


  setNbThreads(1);
  // start openmp for loop
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
  for(int isnp = 0; isnp < n_snps; isnp++) {
    uint32_t snp_index = jt.setinfo[chrom - 1][varset].snp_indices[isnp];
    chi_squared chisq(1);
    double pval_raw;

    // to store variant information
    variant_block* block_info = &(all_snps_info[isnp]);

    if( !params.build_mask ){
      if((params.file_type == "bgen") && params.streamBGEN){ // uncompress and extract the dosages
        parseSnpfromBGEN(isnp, chrom, &(snp_data_blocks[isnp]), insize[isnp], outsize[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &snpinfo[snp_index], &Gblock, block_info, sout);
      } else if(params.file_type == "bed"){ // extract hardcalls
        parseSnpfromBed(isnp, chrom, snp_data_blocks[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &Gblock, block_info);
      }
    }

    // for QTs (or BTs with firth approx): project out covariates & scale
    residualize_geno(isnp, block_info);

    // skip SNP if fails filters
    if( block_info->ignored ) continue;

    block_info->pval_log = ArrayXd::Zero(params.n_pheno);
    block_info->bhat = ArrayXd::Zero(params.n_pheno);
    block_info->se_b = ArrayXd::Zero(params.n_pheno);
    block_info->test_fail.assign(params.n_pheno, false);
    block_info->is_corrected.assign(params.n_pheno, true);
    MapArXd Geno (Gblock.Gmat.col(isnp).data(), params.n_samples, 1);


    if(params.binary_mode) {
      MatrixXd tmpG, WX, GW;
      block_info->stats = ArrayXd::Zero(params.n_pheno);
      block_info->denum = ArrayXd::Zero(params.n_pheno);

      for( int i = 0; i < params.n_pheno; ++i ) {

        // project out covariates from G
        WX = m_ests.Gamma_sqrt.col(i).asDiagonal() * pheno_data.new_cov;
        // apply mask
        GW = (Geno * m_ests.Gamma_sqrt.col(i).array() * pheno_data.masked_indivs.col(i).array().cast<double>()).matrix();
        tmpG = GW - WX * (m_ests.Xt_Gamma_X_inv[i] * (WX.transpose() * GW));
        block_info->denum(i) = tmpG.squaredNorm();

        // score test stat for BT
        block_info->stats(i) = (tmpG.array() * res.col(i).array()).sum() / sqrt( block_info->denum(i) );

        // check this
        if(params.use_SPA) run_SPA_test_snp(block_info, i, tmpG);
      }

    } else {

      // score test stat for QT
      if( params.strict_mode ) {
        Geno *= in_filters.ind_in_analysis.cast<double>();
        block_info->stats = (res.array().colwise() * Geno).matrix().transpose().rowwise().sum() / sqrt( in_filters.ind_in_analysis.cast<float>().sum() );
      } else {
        // compute GtG for each phenotype (different missing patterns)
        block_info->scale_fac_pheno = (pheno_data.masked_indivs.cast<double>().array().colwise() * Geno).matrix().colwise().squaredNorm();
        block_info->stats = (Geno.matrix().transpose() * res).transpose().array() / block_info->scale_fac_pheno.sqrt();
      }

    }

    block_info->chisq_val = block_info->stats.square();

    for( int i = 0; i < params.n_pheno; ++i ) {

      // test statistic & pvalue
      if(!params.use_SPA) check_pval_snp(block_info, chrom, i, isnp);

      // summary stats
      if( !params.binary_mode ){
        // estimate & SE for QT
        if( params.strict_mode )
          block_info->bhat(i) = block_info->stats(i) * ( pheno_data.scale_Y(i) * p_sd_yres(i)) / ( sqrt(pheno_data.Neff(i)) * block_info->scale_fac );
        else
          block_info->bhat(i) = block_info->stats(i) * ( pheno_data.scale_Y(i) * p_sd_yres(i) ) / ( sqrt(block_info->scale_fac_pheno(i)) * block_info->scale_fac );
        block_info->se_b(i) = block_info->bhat(i) / block_info->stats(i);
      } else {
        // with Firth, get sum. stats from Firth logistic regression
        if( params.firth && block_info->is_corrected[i] && !block_info->test_fail[i] ){
          pval_raw = max(params.nl_dbl_dmin, pow(10, - block_info->pval_log(i))); // to prevent overflow
          block_info->chisq_val(i) = quantile(complement(chisq, pval_raw));

          // compute SE from beta & pvalue
          if( params.back_correct_se && (block_info->chisq_val(i) > 0) )
            block_info->se_b(i) = fabs(block_info->bhat(i)) / sqrt(block_info->chisq_val(i));

        } else {
          block_info->se_b(i) = 1 / sqrt(block_info->denum(i));
          // with SPA, calculate test stat based on SPA p-value
          if( params.use_SPA && block_info->is_corrected[i] && !block_info->test_fail[i] ){
            pval_raw = max(params.nl_dbl_dmin, pow(10, - block_info->pval_log(i))); // to prevent overflow
            block_info->chisq_val(i) = quantile(complement(chisq, pval_raw));
            block_info->bhat(i) = sgn(block_info->stats(i)) * sqrt(block_info->chisq_val(i));
          } else block_info->bhat(i) = block_info->stats(i);
          block_info->bhat(i) *= block_info->se_b(i);
          if( params.use_SPA && block_info->flipped ) block_info->bhat(i) *= -1;
        }
        block_info->bhat(i) /= block_info->scale_fac;
        block_info->se_b(i) /= block_info->scale_fac;
      }

    }

    //if( (isnp==0) || (isnp == (n_snps-1)) ) cout << "G"<<isnp+1<<" MAF = " <<  block_info.MAF << endl;
  }

  setNbThreads(params.threads);

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << "done (" << duration.count() << "ms) "<< endl;

}


void Data::readChunk(const int chrom, const int varset, const int start, const int n_snps, vector< vector < uchar > >& snp_data_blocks, vector<uint32_t>& insize, vector<uint32_t>& outsize, vector<variant_block>& all_snps_info){

  if((params.file_type == "bgen") && params.streamBGEN){
    uint64 pos_skip, index;
    ifstream bfile;
    snp_data_blocks.resize( n_snps );
    insize.resize(n_snps); outsize.resize(n_snps);
    bfile.open( files.bgen_file, ios::in | ios::binary );

    for(int isnp = 0; isnp < n_snps; isnp++) {

      // extract genotype data blocks single-threaded
      index = jt.setinfo[chrom - 1][varset].snp_indices[start + isnp];
      pos_skip = snpinfo[index].offset;
      bfile.seekg( pos_skip );

      readChunkFromBGEN(&bfile, &size1, &size2, &(snp_data_blocks[isnp]));
      insize[isnp] = size1;
      outsize[isnp] = size2;

    }
    bfile.close();

  } else if((params.file_type == "bgen") && !params.streamBGEN) readChunkFromBGENFileToG(n_snps, chrom, start, jt.setinfo[chrom - 1][varset].snp_indices, snpinfo, &params, &Gblock, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, all_snps_info);
  else if(params.file_type == "pgen") readChunkFromPGENFileToG(start, n_snps, jt.setinfo[chrom - 1][varset].snp_indices, chrom, &params, &in_filters, &Gblock, pheno_data.masked_indivs, pheno_data.phenotypes_raw, snpinfo, all_snps_info);
  else {
    // read in N/4 bytes from bed file for each snp
    uint64 index;
    snp_data_blocks.resize( n_snps );

    for(int isnp = 0; isnp < n_snps; isnp++) {

      index = snpinfo[ jt.setinfo[chrom - 1][varset].snp_indices[start + isnp] ].offset;
      files.bed_ifstream.seekg(3 + index * files.bed_block_size, ios_base::beg);

      snp_data_blocks[isnp].resize(files.bed_block_size);
      files.bed_ifstream.read( reinterpret_cast<char *> (&snp_data_blocks[isnp][0]), files.bed_block_size);

    }
  }

}

void Data::getMask(const int chrom, const int varset, vector< vector < uchar > >& snp_data_blocks, vector<uint32_t>& insize, vector<uint32_t>& outsize, vector<variant_block>& all_snps_info){

  auto t1 = std::chrono::high_resolution_clock::now();

  // do it in chunks to reduce memory usage
  int n_snps = jt.setinfo[chrom - 1][varset].snp_indices.size();
  int nchunks = ceil( n_snps * 1.0 / params.block_size );
  int bsize = params.block_size; // default number of SNPs to read at a time
  int nvar_read = 0;
  if(params.mask_loo) {
    bm.nmasks_total = n_snps;
    nchunks = 1;
    bsize = n_snps; 
  }

  if(params.verbose) sout << nchunks << " chunks";
  sout << "\n     -reading in genotypes and building masks..." << flush;

  bm.prepMasks(params.n_samples, jt.setinfo[chrom - 1][varset].ID);  
  Gblock.Gmat.resize(params.n_samples, bsize);


  for(int i = 0; i < nchunks; i++){

    if( i == (nchunks-1) ) {
      bsize = n_snps - i * bsize;// use remainder number of variants
      Gblock.Gmat.resize(params.n_samples, bsize);
    }

    readChunk(chrom, varset, nvar_read, bsize, snp_data_blocks, insize, outsize, all_snps_info);


    // build genotype matrix
#if defined(_OPENMP)
    setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
    for(int isnp = 0; isnp < bsize; isnp++) {

      uint32_t snp_index = jt.setinfo[chrom - 1][varset].snp_indices[nvar_read + isnp];

      variant_block* block_info = &(all_snps_info[isnp]);

      if(params.file_type == "bgen"){ // uncompress and extract the dosages
        parseSnpfromBGEN(isnp, chrom, &(snp_data_blocks[isnp]), insize[isnp], outsize[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &snpinfo[snp_index], &Gblock, block_info, sout);
      } else if(params.file_type == "bed"){ // extract hardcalls
        parseSnpfromBed(isnp, chrom, snp_data_blocks[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &Gblock, block_info);
      }

    }
#if defined(_OPENMP)
    setNbThreads(params.threads);
#endif

    // update mask (taking max/sum)
    if(params.mask_loo)
      bm.updateMasks_loo(nvar_read, bsize, &params, &in_filters, pheno_data.masked_indivs, &Gblock, all_snps_info, jt.setinfo[chrom - 1][varset], snpinfo, sout);
    else
      bm.updateMasks(nvar_read, bsize, &params, &in_filters, pheno_data.masked_indivs, &Gblock, all_snps_info, jt.setinfo[chrom - 1][varset], snpinfo, sout);

    nvar_read += bsize;
  }

  // check mask and store in setinfo & snpinfo
  if(params.mask_loo)
    bm.computeMasks_loo(&params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &Gblock, all_snps_info, jt.setinfo[chrom - 1][varset], snpinfo, sout);
  else
    bm.computeMasks(&params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &Gblock, all_snps_info, jt.setinfo[chrom - 1][varset], snpinfo, sout);


  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;

}
