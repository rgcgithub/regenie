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

#include <limits.h> /* for PATH_MAX */
#include <chrono>

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wparentheses"
#endif
#include <boost/filesystem.hpp>
#include <boost/exception/all.hpp>
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include "Regenie.hpp"
#include "Files.hpp"
#include "Geno.hpp"
#include "Joint_Tests.hpp"
#include "Step1_Models.hpp"
#include "Step2_Models.hpp"
#include "Pheno.hpp"
#include "MultiTrait_Tests.hpp"
#include "Ordinal.hpp"
#include "HLM.hpp"
#include "SKAT.hpp"
#include "Interaction.hpp"
#include "Masks.hpp"
#include "Data.hpp"

#ifdef WITH_HTSLIB
#include "remeta/regenie_ld_matrix_writer.hpp"
#endif

using namespace std;
using namespace Eigen;
using namespace boost;
namespace fs = boost::filesystem;


using boost::math::normal;
using boost::math::chi_squared;

Data::Data() { // @suppress("Class members should be properly initialized")
}

Data::~Data() {
  // TODO Auto-generated destructor stub
}


void Data::run() {

  // set number of threads
  set_threads(&params);

  if(params.streamBGEN) check_bgen(files.bgen_file, params.file_type, params.zlib_compress, params.streamBGEN, params.BGENbits, params.nChrom);

  if(params.test_mode)  // step 2
    run_step2();
  else  // step 1
    run_step1();

}

void Data::run_step1(){

  sout << "Fitting null model\n";

  // set up file for reading
  file_read_initialization();
  // if splitting l0 into many jobs
  if(params.split_l0) set_parallel_l0();
  // read phenotype and covariate files
  read_pheno_and_cov(&files, &params, &in_filters, &pheno_data, &m_ests, &Gblock, sout);
  // adjust for covariates
  prep_run(&files, &in_filters, &params, &pheno_data, &m_ests, sout);
  // set number of blocks and block size and ridge parameters
  set_blocks();
  // some initializations
  setmem();
  // level 0
  level_0_calculations();
  // print y/x/logreg offset used for level 1 
  if(params.debug) write_inputs();
  // prep for level 1 models
  prep_l1_models();
  // level 1 ridge
  if(params.trait_mode == 0){ // QT
    if(params.use_loocv) ridge_level_1_loocv(&files, &params, &pheno_data, &l1_ests, sout);
    else ridge_level_1(&files, &params, &pheno_data, &l1_ests, sout);
  } else if(params.trait_mode == 1){ // BT
    if(params.use_loocv) ridge_logistic_level_1_loocv(&files, &params, &pheno_data, &m_ests, &l1_ests, sout);
    else ridge_logistic_level_1(&files, &params, &pheno_data, &l1_ests, masked_in_folds, sout);
  } else if(params.trait_mode == 2){ // CT
    if(params.use_loocv) ridge_poisson_level_1_loocv(&files, &params, &pheno_data, &m_ests, &l1_ests, sout);
    else ridge_poisson_level_1(&files, &params, &pheno_data, &l1_ests, masked_in_folds, sout);
  }
  // output results
  output();

}

void Data::run_step2(){

  // allocate per thread if using OpenMP
  Gblock.thread_data.resize(params.neff_threads);

  if(params.getCorMat) ld_comp();
  else if( params.snp_set ) test_joint();
  else if (params.trait_set) test_multitrait();
  else if (params.multiphen) test_multiphen();
  else test_snps_fast();

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

  if( params.condition_snps )
    get_conditional_vars(in_filters.condition_snp_names, &files, &params, sout);

  if(params.file_type == "bed") read_bed_bim_fam(&files, &params, &in_filters, snpinfo, chr_map, sout);
  else if(params.file_type == "pgen") read_pgen_pvar_psam(&files, &params, &in_filters, &Gblock, snpinfo, chr_map, sout);
  else prep_bgen(&files, &params, &in_filters, snpinfo, chr_map, Gblock.bgen, sout);

  params.nvs_stored = snpinfo.size();
  //if(params.getCorMat) params.block_size = params.n_variants;

  if(!params.test_mode && !params.force_run && ((int)params.nvs_stored > params.max_step1_variants))
    throw "it is not recommened to use more than " + to_string( params.max_step1_variants ) + 
      " variants in step 1 (otherwise use '--force-step1'). " + params.webinfo ;

  if( params.setMinINFO && !params.dosage_mode )
    sout << "WARNING: Dosages are not present in the genotype file. Option --minINFO is skipped.\n";

}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
////          adjust for covariates in G
/////////////////////////////////////////////////
/////////////////////////////////////////////////

// only for step 1
void Data::residualize_genotypes() {

  sout << "   -residualizing and scaling genotypes..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  // mask missing individuals
  Gblock.Gmat.array().rowwise() *= in_filters.ind_in_analysis.matrix().transpose().array().cast<double>();

  // residuals (centered)
  MatrixXd beta = Gblock.Gmat * pheno_data.new_cov;
  Gblock.Gmat -= beta * pheno_data.new_cov.transpose();

  // scaling (use [N-C] where C=#covariates)
  scale_G = Gblock.Gmat.rowwise().norm() / sqrt(params.n_analyzed - params.ncov);

  // check sd
  MatrixXd::Index minIndex;
  if(scale_G.array().minCoeff(&minIndex) < params.numtol) 
    throw "!! Uh-oh, SNP " + snpinfo[in_filters.step1_snp_count+minIndex].ID + 
      " has low variance (=" + to_string( scale_G(minIndex,0) ) + ").";

  Gblock.Gmat.array().colwise() /= scale_G.array();

  // to use MAF dependent prior on effect size [only for step 1]
  // multiply by [p*(1-p)]^(1+alpha)/2
  if(params.alpha_prior != -1) 
    Gblock.Gmat.array().colwise() *= pow(Gblock.snp_afs.col(0).array() * (1-Gblock.snp_afs.col(0).array()), 0.5 * (params.alpha_prior + 1) );


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

  if(params.njobs <= 1)
    throw "number of jobs must be >1.";
  else if(params.njobs > params.total_n_block){

    sout << "   -WARNING: Number of jobs cannot be greater than number of blocks.\n";
    params.njobs = params.total_n_block;

  }

  sout << "   -using " << params.njobs << " jobs\n";
  sout << "   -master file written to [" << fout << "]\n";
  sout << "   -variant list files written to [" << files.split_file << "_job*.snplist]\n";

  // open master
  ofstream ofile;
  openStream(&ofile, fout, ios::out, sout);

  // header
  ofile << params.nvs_stored << " " << params.block_size << endl;

  // split blocks in chunks of ~B/njobs
  int nall = params.total_n_block / params.njobs;
  int remainder = params.total_n_block - nall * params.njobs;
  int nb = 0, bs, ns = 0, bcount = 0, scount = 0, jcount = 0;
  int btarget = nall + (jcount < remainder ? 1 : 0);
  map<int, vector<int> >::iterator itr;

  for (itr = chr_map.begin(); itr != chr_map.end(); ++itr) {
    int chrom_nsnps = itr->second[0];
    int chrom_nb = ceil(chrom_nsnps * 1.0 / params.block_size);
    if(chrom_nb == 0) continue;

    for(int bb = 0; bb < chrom_nb ; bb++) {

      get_block_size(params.block_size, chrom_nsnps, bb, bs);

      ns+=bs;
      nb++, bcount++;

      if( nb == btarget ){
        string fname = files.split_file + "_job" + to_string( jcount+1 );
        // write in master
        ofile << fname << " " << btarget << " " << ns << endl;
        // write snplist
        writeSnplist(fname, scount, ns, snpinfo, sout);

        jcount++;
        scount += ns;
        ns = nb = 0;
        btarget = nall + (jcount < remainder ? 1 : 0);
      }
    }
  }

  if((bcount != params.total_n_block) || (jcount !=params.njobs))
    throw "could not create master file.";

  ofile.close();

}

void Data::set_blocks() {

  params.total_n_block = 0, total_chrs_loco = 0;
  int blocks_left = params.n_block;
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
    m1[ itr->first ] = itr->second;
  }
  chr_map = m1;
  //sout << "#chrs = "<< chr_map.size() << ";#loco chrs = "<< total_chrs_loco << endl;

  if(params.total_n_block == 0)
    throw "total number of blocks must be > 0.";

  if(params.split_l0) return;
  else if(params.run_l0_only) {
    if((params.parallel_nBlocks != params.total_n_block) || (params.parallel_nSnps!= (int)params.n_variants))
      throw "number of variants/blocks in file (=" + to_string(params.parallel_nSnps) + "/" + to_string(params.total_n_block) +
        ") don't match with that in master file (=" + to_string(params.n_variants) + "/" + to_string(params.parallel_nBlocks) +").";
  } else if(params.run_l1_only) prep_parallel_l1();

  // for BTs: check if the sample size is lower than 5K (if so, force loocv)
  if( (params.trait_mode == 1) && !params.use_loocv && ( params.n_analyzed < 5000) ) {
    sout << "   -WARNING: Sample size is less than 5,000 so using LOOCV instead of " << params.cv_folds << "-fold CV.\n";
    params.use_loocv = true;
  }

  /*
  // check block size vs sample size
  if(params.use_loocv && params.block_size > params.n_analyzed)
    throw "block size must be smaller than the number of samples to perform LOOCV!";
  */
  if(params.use_loocv) params.cv_folds = params.n_samples;

  uint32_t neff_folds = params.use_loocv ? params.n_analyzed : params.cv_folds;

  // summarize block sizes and ridge params
  sout << left << std::setw(20) << " * # threads" << ": [" << params.threads << "]\n";
  sout << left << std::setw(20) << " * block size" << ": [" << params.block_size << "]\n";
  sout << left << std::setw(20) << " * # blocks" << ": [" << params.total_n_block << "] for " << params.nvs_stored << " variants\n";
  sout << left << std::setw(20) << " * # CV folds" << ": [" << neff_folds << "]\n";

  if(!params.run_l1_only){
    IOFormat Fmt(FullPrecision, DontAlignCols, " ", " ", "", "","","");
    sout << left << std::setw(20) << " * ridge data_l0" << ": [" << params.n_ridge_l0 << " : " << params.lambda.format(Fmt) << " ]\n";
  }

  if(!params.run_l0_only){
    IOFormat Fmt(FullPrecision, DontAlignCols, " ", " ", "", "","","");
    sout << left << std::setw(20) << " * ridge data_l1" << ": [" << params.n_ridge_l1 << " : " << params.tau[0].format(Fmt) << " ]\n";
  }

  // if using maf dependent prior
  if(!params.test_mode && (params.alpha_prior != -1) ) 
    sout << " * applying a MAF dependent prior to the SNP effect sizes in level 0 models (alpha=" << params.alpha_prior << ")\n";

  // print approx. amount of memory needed
  print_usage_info(&params, &files, sout);

  // storing null estimates from firth
  if(params.write_null_firth ) 
    sout << " * writing null Firth estimates to file\n";

  // if within sample predictions are used in level 1
  if (params.within_sample_l0) 
    sout << " * using within-sample predictions from level 0 as features at level 1\n";

}


void Data::set_folds() {

  // set up folds
  params.cv_sizes.resize(params.cv_folds, 1);

  // assign folds for individuals in analysis
  if( !params.use_loocv ){

    uint32_t target_size_folds = floor( params.n_analyzed / params.cv_folds );
    if( target_size_folds < 1 )
      throw "not enough samples are present for " + to_string( params.cv_folds ) + "-fold CV.";

    uint32_t n_non_miss = 0, cum_size_folds = 0;
    int cur_fold = 0;
    for(size_t i = 0; i < params.n_samples; i++){

      if( in_filters.ind_in_analysis(i) ) n_non_miss++;

      if( n_non_miss == target_size_folds){
        params.cv_sizes(cur_fold) = i - cum_size_folds + 1;
        cum_size_folds += params.cv_sizes(cur_fold);
        n_non_miss = 0, cur_fold++;
      } else if( cur_fold == (params.cv_folds - 1) ){
        params.cv_sizes(cur_fold) = params.n_samples - i;
        break;
      }

      //sout << i << " " << cur_fold << " " << n_non_miss << " " << in_filters.ind_in_analysis(i) << " "<< target_size_folds << endl;
    }

  } else // loocv
    params.cv_sizes = ArrayXi::Constant(params.cv_folds, 1);;


  // check sd(Y) in folds
  if(!params.use_loocv && params.trait_mode){

    int minIndex;
    uint32_t cum_size_folds = 0;
    ArrayXd sum, n_cv, sd_phenos;

    for(int i = 0; i < params.cv_folds; i++) {

      MatrixXb M = pheno_data.masked_indivs.middleRows(cum_size_folds, params.cv_sizes(i)); // nxp
      MatrixXd Y = (pheno_data.phenotypes_raw.middleRows(cum_size_folds, params.cv_sizes(i)).array() * M.array().cast<double>()).matrix().transpose(); // pxn

      sum = params.pheno_pass.select( Y.array().rowwise().sum() , 10);

      // BTs
      if(params.trait_mode == 1){
        n_cv = params.pheno_pass.select( M.transpose().array().rowwise().count().cast<double>() , 100);
        sd_phenos = (sum/n_cv) * (1 - sum/n_cv);

        if( sd_phenos.minCoeff(&minIndex) < params.numtol )
          throw "one of the folds has only cases/controls for phenotype '" + files.pheno_names[minIndex] 
            + "'. Either use smaller #folds (option --cv) or use LOOCV (option --loocv).";
      } else if(params.trait_mode == 2){

        if( sum.minCoeff(&minIndex) == 0 )
          throw "one of the folds has only zero counts for phenotype '" + files.pheno_names[minIndex] 
            + "'. Either use smaller #folds (option --cv) or use LOOCV (option --loocv).";
      }

      cum_size_folds += params.cv_sizes(i);
    }

  }

  // only used for K-fold CV
  if(!params.use_loocv && !params.within_sample_l0){
    l1_ests.X_folds.resize(params.cv_folds);
    l1_ests.XtY.resize(params.cv_folds);	
  }

}


void Data::setmem() {

  bool is_set = false;
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

  if(params.trait_mode){ // non-QT
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

    if( !params.pheno_pass(i) && (!params.write_l0_pred || (i!=0)) ) continue;

    if (params.within_sample_l0) {
      l1_ests.pred_mat[i].resize(params.cv_folds);
      l1_ests.pred_pheno[i].resize(params.cv_folds);
    } else if(!params.use_loocv) l1_ests.beta_hat_level_1[i].resize(params.cv_folds);

    if(!params.use_loocv) {
      l1_ests.test_pheno[i].resize(params.cv_folds);
      l1_ests.test_mat[i].resize(params.cv_folds);
    } else l1_ests.test_mat_conc[i] = MatrixXd::Zero(params.n_samples, params.n_ridge_l0 * ( params.write_l0_pred ? 1 : params.total_n_block) );

    if(params.trait_mode) {
      if (params.within_sample_l0) {
        l1_ests.pred_pheno_raw[i].resize(params.cv_folds);
        l1_ests.pred_offset[i].resize(params.cv_folds);
      }
      l1_ests.test_pheno_raw[i].resize(params.cv_folds);
      if(!params.use_loocv) l1_ests.test_offset[i].resize(params.cv_folds);
    }

    for(int j = 0; j < params.cv_folds; ++j ) {

      if (params.within_sample_l0) {
        l1_ests.pred_mat[i][j] = MatrixXd::Zero(params.n_samples - params.cv_sizes(j), params.total_n_block * params.n_ridge_l0);
        l1_ests.pred_pheno[i][j] = MatrixXd::Zero(params.n_samples - params.cv_sizes(j), 1);
      } else if(!params.use_loocv) l1_ests.beta_hat_level_1[i][j] = MatrixXd::Zero(params.total_n_block * params.n_ridge_l0, params.n_ridge_l1);

      if(!params.use_loocv) {
        l1_ests.test_pheno[i][j] = MatrixXd::Zero(params.cv_sizes(j), 1);
        l1_ests.test_mat[i][j] = MatrixXd::Zero(params.cv_sizes(j), params.n_ridge_l0 * ( params.write_l0_pred ? 1 : params.total_n_block));
      }

      if(params.trait_mode) {
        if (params.within_sample_l0) {
          l1_ests.pred_pheno_raw[i][j] = MatrixXd::Zero(params.n_samples - params.cv_sizes(j), 1);
          l1_ests.pred_offset[i][j] = MatrixXd::Zero(params.n_samples - params.cv_sizes(j), 1);
        }
        l1_ests.test_pheno_raw[i][j] = MatrixXd::Zero(params.cv_sizes(j), 1);
        if(!params.use_loocv) l1_ests.test_offset[i][j] = MatrixXd::Zero(params.cv_sizes(j), 1);
      }

    }

    if(!is_set) {// only done once
      for(int j = 0; j < params.cv_folds; ++j ) 
        masked_in_folds[j] = MatrixXb::Constant(params.cv_sizes(j), params.n_pheno, false);
      is_set = true;
    }


  }

  sout << "done\n\n";

}

void Data::get_block_size(int const& target, int const& total, int const& block, int& bs){

  if( ((block + 1) * target) > total)
    bs = total - (block * target);
  else
    bs = target;

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

  int block = 0, bs;
  if(params.print_block_betas) params.print_snpcount = 0;
  ridgel0 l0;

  // set ridge params
  params.lambda = (params.run_l0_only ? params.parallel_nGeno : params.n_variants) * (1 - params.lambda) / params.lambda;

  if(!params.use_loocv){
    l0.G_folds.resize(params.cv_folds);
    l0.GtY.resize(params.cv_folds);
  }

  // open streams to write level 0 predictions
  if(params.write_l0_pred){
    string fout_p;
    files.write_preds_files.resize(params.n_pheno);
    for(int ph = 0; ph < params.n_pheno; ph++){
      if( !params.pheno_pass(ph) ) continue;
      files.write_preds_files[ph] = std::make_shared<ofstream>();
      fout_p = files.loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
      openStream(files.write_preds_files[ph].get(), fout_p, ios::out | ios::binary, sout);
    }
  }

  Files ofile_p;
  if(params.test_l0){ // check strength of association for each block
    params.l0_pvals_file = files.out_file + "_p.txt";
    ofile_p.openForWrite(params.l0_pvals_file, sout);
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

      get_block_size(params.block_size, chrom_nsnps, bb, bs);

      Gblock.Gmat = MatrixXd::Zero(bs, params.n_samples);
      if(params.alpha_prior != -1) Gblock.snp_afs = MatrixXd::Zero(bs, 1);

      get_G(block, bs, chrom, in_filters.step1_snp_count, snpinfo, &params, &files, &Gblock, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, sout);

      // residualize and scale genotypes
      residualize_genotypes();

      // calc working matrices for ridge regressions across folds
      calc_cv_matrices(&l0);

      // test association for block
      if(params.test_l0)
        test_assoc_block(chrom, block, l0, ofile_p, params);

      // calc level 0 ridge regressions
      if(params.use_loocv)
        ridge_level_0_loocv(block, &files, &params, &in_filters, &m_ests, &Gblock, &pheno_data, snpinfo, &l0, &l1_ests, sout);
      else
        ridge_level_0(block, &files, &params, &in_filters, &m_ests, &Gblock, &pheno_data, snpinfo, &l0, &l1_ests, masked_in_folds, sout);

      block++; in_filters.step1_snp_count += bs;
    }

  }

  // close streams
  if(params.write_l0_pred)
    for( auto &yfile : files.write_preds_files)
      if(yfile->is_open()) yfile->close();

  if(params.test_l0)
    ofile_p.closeFile();

  if(params.early_exit) {
    sout << "\nDone printing out level 0 predictions. There are " <<
      params.n_samples << " rows and " <<
      params.total_n_block * params.n_ridge_l0 << " columns " <<
      "stored in column-major order. Exiting...\n";
    exit_early();
  } else if(params.run_l0_only) {
    sout << "\nDone writing level 0 predictions to file.\n";
    exit_early();
  }

  // free up memory not used anymore
  Gblock.Gmat.resize(0,0);
  if(params.write_l0_pred && (params.n_pheno > 1) ) 
    // free level 0 predictions for (P-1) indices in test_mat
    for(int ph = 1; ph < params.n_pheno; ++ph ) {
      if( !params.pheno_pass(ph) ) continue;
      if(!params.use_loocv){ // k-fold
        for(int i = 0; i < params.cv_folds; ++i ) 
          l1_ests.test_mat[ph][i].resize(0,0);
        l1_ests.test_mat[ph].resize(0);
      } else l1_ests.test_mat_conc[ph].resize(0,0); // loocv
    }

}

void Data::calc_cv_matrices(struct ridgel0* l0) {

  sout << "   -calc working matrices..." << flush;
  auto t2 = std::chrono::high_resolution_clock::now();
  int bs = Gblock.Gmat.rows();

  if(!params.use_loocv){ // k-fold

    l0->GGt.setZero(bs,bs);
    l0->GTY.setZero(bs,params.n_pheno);
    uint32_t cum_size_folds = 0;

    for( int i = 0; i < params.cv_folds; ++i ) {
      MapMatXd Gmat (&(Gblock.Gmat(0,cum_size_folds)), bs, params.cv_sizes(i));
      l0->GtY[i] = Gmat * pheno_data.phenotypes.middleRows(cum_size_folds, params.cv_sizes(i));
      l0->GTY += l0->GtY[i];
      l0->G_folds[i] = Gmat * Gmat.transpose();
      l0->GGt += l0->G_folds[i];

      cum_size_folds += params.cv_sizes(i);
    }

  } else { // loocv

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

// select which level 0 predictors to use at level 1
void Data::prep_l1_models(){

  int bs_l1 = params.total_n_block * params.n_ridge_l0;
  // arrayxb of which level 0 predictors to keep (default is all)
  l1_ests.l0_colkeep = MatrixXb::Constant(bs_l1, params.n_pheno, true);

  if(params.select_l0){
    // read in pvals for level 0 blocks
    int lineread = 0;
    string line;
    std::vector< string > tmp_str_vec ;
    Files fClass;
    fClass.openForRead(params.l0_pvals_file, sout);

    l1_ests.l0_pv_block.resize(params.total_n_block, params.n_pheno);
    l1_ests.chrom_block.resize(params.total_n_block);

    while( fClass.readLine(line) ){
      tmp_str_vec = string_split(line," ");
      if(lineread >= params.total_n_block)
        throw "number of blocks in file is greater than that analyzed in run.";
      l1_ests.chrom_block(lineread) = atoi(tmp_str_vec[0].c_str());
      if( (int)tmp_str_vec.size() > (params.n_pheno + 2))
        throw "number of phenotypes in file is greater than that analyzed in run.";
      for(int i = 0; i < params.n_pheno; i++)
        l1_ests.l0_pv_block(lineread, i) = convertDouble(tmp_str_vec[i + 2], &params, sout);
      lineread++;
    }
  }

  // set ridge params
  ArrayXd base_tau = params.tau[0];
  params.tau.assign(params.n_pheno, base_tau);

  // for chr map
  l1_ests.chrom_map_ndiff = ArrayXi::Zero(params.nChrom);

}

// identify which block to analyze
void Data::prep_parallel_l0(){

  int tmpi;
  string line, fin = files.split_file; // master file
  std::vector< string > tmp_str_vec ;
  ifstream infile;

  // print info
  sout << " * running jobs in parallel (job #" << params.job_num << ")\n";

  openStream(&infile, fin, ios::in, sout);

  // check header
  if(!getline(infile, line))
    throw "cannot read header line in master file."; 

  if( (sscanf( line.c_str(), "%d %d", &params.parallel_nGeno, &tmpi ) != 2) || (tmpi != params.block_size) )
    throw "invalid header line in master file."; 

  // skip to line job_num
  int nskip=1;
  while( (nskip++ < params.job_num) && !infile.eof() )
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  if( (--nskip != params.job_num) || infile.eof() )
    throw "could not read line " + to_string( params.job_num+1 ) + " (check number of lines in file).";
  
  // read in line
  getline(infile, line);
  char tmp_chr[MAXFILELEN];

  if( sscanf( line.c_str(), "%s %d %d", tmp_chr, &params.parallel_nBlocks, &params.parallel_nSnps ) != 3 )
    throw "could not read line " + to_string( params.job_num + 1 ) + " (check number of lines and format in file)."; 

  files.loco_tmp_prefix = tmp_chr;
  files.file_snps_include.resize(1);
  files.file_snps_include[0] = files.loco_tmp_prefix + ".snplist";

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

  openStream(&infile, fin, ios::in, sout);

  // check header
  if(!getline(infile, line))
    throw "cannot read header line in master file."; 
  if( (sscanf( line.c_str(), "%d %d", &params.parallel_nGeno, &nb ) != 2) || (nb != params.block_size) )
    throw "invalid header line in master file."; 

  nblocks = 0, nsnps = 0, lineread=0;
  while( getline(infile, line) ){

    char tmp_chr[MAXFILELEN];
    if( sscanf( line.c_str(), "%s %d %d", tmp_chr, &nb, &ns ) != 3 )
      throw "could not read line " + to_string( params.job_num + 1 ) + " (check number of lines and format in file)."; 

    files.bstart.push_back( nblocks );
    files.btot.push_back( nb );
    files.mprefix.push_back( string(tmp_chr) );

    // check params
    if( (files.bstart[lineread] < 0) || (files.bstart[lineread]>params.total_n_block) || (files.btot[lineread] < 0) )
      throw "invalid block information in master file at line " + to_string( lineread + 2 ) + ".";

    nblocks += nb; // update # blocks
    nsnps += ns;
    lineread++;
  }

  if((nblocks != params.total_n_block) || (nsnps != params.n_variants))
    throw "number of blocks/variants in master file '" + fin + "' doesn't match that in the analysis.";

  // print info
  params.job_num = lineread;
  sout << " * using results from running " << params.job_num << " parallel jobs at level 0\n";

  infile.close();

}


void Data::write_inputs(){

  // write Y
  IOFormat Fmt(FullPrecision, DontAlignCols, " ", "\n", "", "","","");
  ofstream ofile;
  openStream(&ofile, files.out_file + "_y.txt", ios::out, sout);
  if(params.trait_mode == 0)
    ofile << pheno_data.phenotypes.format(Fmt) << "\n";
  else
    ofile << pheno_data.phenotypes_raw.format(Fmt) << "\n";
  ofile.close();

  // write X
  openStream(&ofile, files.out_file + "_x.txt", ios::out, sout);
  ofile << pheno_data.new_cov.format(Fmt) << "\n";
  ofile.close();

  // write offset
  if(params.trait_mode != 0){
    openStream(&ofile, files.out_file + "_offset.txt", ios::out, sout);
    ofile << m_ests.offset_nullreg.format(Fmt) << "\n";
    ofile.close();
  }
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
  double rate=0, zv;
  string pfile, out_blup_list, out_prs_list, out_firth_list, loco_filename, prs_filename, firth_filename;
  string fullpath_str, path_prs, path_firth;
  Files outb, outp, outf;

  sout << "Output\n------\n";

  if(params.make_loco || params.trait_mode){
    out_blup_list = files.out_file + "_pred.list";
    outb.openForWrite(out_blup_list, sout);
  }

  if(params.print_prs){
    out_prs_list = files.out_file + "_prs.list";
    outp.openForWrite(out_prs_list, sout);
  }

  if(params.write_null_firth){
    out_firth_list = files.out_file + "_firth.list";
    outf.openForWrite(out_firth_list, sout);
  }

  for(int ph = 0; ph < params.n_pheno; ++ph ) {
    if( !params.pheno_pass(ph) ) continue;

    sout << "phenotype " << ph+1 << " (" << files.pheno_names[ph] << ") : " ;
    loco_filename = files.out_file + "_" + to_string(ph + 1) + ".loco" + (params.gzOut ? ".gz" : "");
    prs_filename = files.out_file + "_" + to_string(ph + 1) + ".prs" + (params.gzOut ? ".gz" : "");
    firth_filename = files.out_file + "_" + to_string(ph + 1) + ".firth" + (params.gzOut ? ".gz" : "");

    if( params.make_loco || params.trait_mode || params.print_prs ) {

      fullpath_str = (params.use_rel_path ? loco_filename : get_fullpath(loco_filename));
      if(params.print_prs) path_prs = (params.use_rel_path ? prs_filename : get_fullpath(prs_filename));

      if(params.trait_mode == 0) { // for quantitative traits

        outb << files.pheno_names[ph]  << " " <<  fullpath_str << endl;
        if(params.print_prs) 
          outp << files.pheno_names[ph]  << " " <<  path_prs << endl;

      } else { // check level 1 ridge converged

        if( !l1_ests.pheno_l1_not_converged(ph) ) {
          outb << files.pheno_names[ph]  << " " << fullpath_str << endl;
          if(params.print_prs) 
            outp << files.pheno_names[ph]  << " " <<  path_prs << endl;

        } else { // failed level 1

          if(params.write_l0_pred && params.rm_l0_pred) 
            rm_l0_files(ph); // cleanup level 0 predictions
          sout << "Level 1 model did not converge. LOCO predictions calculations are skipped.\n\n";
          continue;

        }
      }
    }
    sout << endl;

    min_index = 0;
    min_val = 1e10;

    // determine optimal parameter by cv using: QT: MSE, nonQT: -loglik
    for(int j = 0; j < params.n_ridge_l1; ++j ) {
      if(params.trait_mode == 0)
        performance_measure = l1_ests.cumsum_values[2](ph, j) + l1_ests.cumsum_values[3](ph,j) - 2 * l1_ests.cumsum_values[4](ph,j);
      else
        performance_measure = l1_ests.cumsum_values[5](ph, j);

      performance_measure /= pheno_data.Neff(ph);

      if( performance_measure < min_val) {
        min_index = j;
        min_val = performance_measure;
      }
    }

    if(params.trait_mode == 2)
      rate = pheno_data.phenotypes_raw.col(ph).sum() / pheno_data.Neff(ph); // separate for each trait

    for(int j = 0; j < params.n_ridge_l1; ++j ) {

      if(params.trait_mode == 2){
        zv = exp(l1_ests.l0_colkeep.col(ph).count() / params.tau[ph](j)) - 1; 
        sout << "  " << setw(5) << rate * zv / (1 + rate * zv);
      } else sout << "  " << setw(5) << l1_ests.l0_colkeep.col(ph).count() / (l1_ests.l0_colkeep.col(ph).count() + (params.trait_mode == 1? (M_PI * M_PI / 3) : 1) * params.tau[ph](j) );

      // output Rsq and MSE
      rsq = l1_ests.cumsum_values[4](ph,j) - l1_ests.cumsum_values[0](ph,j) * l1_ests.cumsum_values[1](ph,j) / pheno_data.Neff(ph); // num = Sxy - SxSy/n
      rsq = (rsq * rsq) / ((l1_ests.cumsum_values[2](ph,j) - l1_ests.cumsum_values[0](ph,j) * l1_ests.cumsum_values[0](ph,j) / pheno_data.Neff(ph)) * (l1_ests.cumsum_values[3](ph,j) - l1_ests.cumsum_values[1](ph,j) * l1_ests.cumsum_values[1](ph,j) / pheno_data.Neff(ph))); // num^2 / ( (Sx2 - Sx^2/n)* (Sy2 - Sy^2/n) )
      sse = l1_ests.cumsum_values[2](ph, j) + l1_ests.cumsum_values[3](ph,j) - 2 * l1_ests.cumsum_values[4](ph,j); // Sx2 + Sy2 - SxSy
      if(params.trait_mode) ll_avg = l1_ests.cumsum_values[5](ph, j) / pheno_data.Neff(ph);

      sout << " : " 
        << "Rsq = " << rsq;
      if(params.trait_mode!=2) 
        sout  << ", MSE = " << sse/pheno_data.Neff(ph);
      if(params.trait_mode) 
        sout << ", -logLik/N = " << ll_avg;
      if(j == min_index) 
        sout << "<- min value";
      sout << endl;

    }

    if(params.trait_mode == 0){
      if(params.use_loocv) 
        make_predictions_loocv(ph, min_index);
      else 
        make_predictions(ph, min_index);
    } else if(params.trait_mode == 1){
      if(params.l1_full_samples) 
        make_predictions_binary_loocv_full(ph, min_index);
      else if(params.use_loocv) 
        make_predictions_binary_loocv(ph, min_index);
      else 
        make_predictions_binary(ph, min_index);
    } else if(params.trait_mode == 2){
      if(params.use_loocv) 
        make_predictions_count_loocv(ph, min_index);
      else 
        make_predictions_count(ph, min_index);
    }

    // check if firth estimates converged (should have been written to file)
    if(params.write_null_firth && file_exists(firth_filename)){
      path_firth = (params.use_rel_path ? firth_filename : get_fullpath(firth_filename));
      outf << files.pheno_names[ph]  << " " <<  path_firth << endl;
    }

    // delete file used to store l0 predictions
    if(params.write_l0_pred && params.rm_l0_pred)
      rm_l0_files(ph);

  }

  if(params.make_loco || (params.trait_mode!=0)){
    outb.closeFile();
    sout << "List of blup files written to: [" 
      << out_blup_list << "]\n";
  }
  if(params.write_null_firth) {
    outf.closeFile();
    sout << "List of files with null Firth estimates written to: [" 
      << out_firth_list << "]\n";
  }
  if(params.print_prs) {
    outp.closeFile();
    sout << "List of files with whole genome PRS written to: [" 
      << out_prs_list << "]\n";
  }

}

void Data::rm_l0_files(int const& ph){

  string pfile;

  if(!params.run_l1_only){
    pfile = files.loco_tmp_prefix + "_l0_Y" + to_string(ph+1);
    remove(pfile.c_str());
  } else {
    for(auto const& pfx : files.mprefix){
      pfile = pfx + "_l0_Y" + to_string(ph+1);
      if(file_exists(pfile)) remove(pfile.c_str()); // l0 predictions
      pfile = pfx + ".snplist";
      if(file_exists(pfile)) remove(pfile.c_str()); // snplist
    }
  }

}

// convert filename to full path
std::string get_fullpath(std::string fname){

  string fout;
  fs::path fullpath;

  try {

    // convert to full path using boost filesystem library
    // this can generate errors due to LC_ALL locale being invalid
    fullpath = fs::absolute(fname);
    fout = fullpath.make_preferred().string();

  } catch ( std::runtime_error& ex ) {

    // to avoid boost::filesystem error
    setenv("LC_ALL", "C", 1);

    try {

      // try again
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
      } catch ( std::runtime_error& ex ) {
        fout = fname; // if failed to get full path
      }

    }
  }

  return fout;

}

void Data::make_predictions(int const& ph, int const& val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();
  int ph_eff = params.write_l0_pred ? 0 : ph;

  // read in level 0 predictions from file
  if(params.write_l0_pred)
    read_l0(ph, ph_eff, &files, &params, &l1_ests, sout);
  check_l0(ph, ph_eff, &params, &l1_ests, &pheno_data, sout, true);

  int bs_l1 = l1_ests.test_mat[ph_eff][0].cols();
  MatrixXd ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);
  MatrixXd X1, X2, beta_l1, beta_avg;
  string outname;
  ofstream ofile;

  if(params.within_sample_l0){
    X1 = l1_ests.test_mat[ph_eff][0].transpose() * l1_ests.test_mat[ph_eff][0];
    X2 = l1_ests.test_mat[ph_eff][0].transpose() * l1_ests.test_pheno[ph][0];
    for(int i = 1; i < params.cv_folds; ++i ) {
      X1 += l1_ests.test_mat[ph_eff][i].transpose() * l1_ests.test_mat[ph_eff][i];
      X2 += l1_ests.test_mat[ph_eff][i].transpose() * l1_ests.test_pheno[ph][i];
    }
    beta_l1 = (X1 + params.tau[ph](val) * ident_l1).llt().solve(X2);
  } else if(params.print_block_betas) {
    beta_avg = MatrixXd::Zero(bs_l1, 1);
    for(int i = 0; i < params.cv_folds; ++i ) {
      beta_avg += l1_ests.beta_hat_level_1[ph][i].col(val);
    }
    beta_avg /= params.cv_folds;
  }

  // if specified, write betas to file (open in append mode)
  if(params.print_block_betas && !params.within_sample_l0 ) {
    outname = files.out_file + "_level1.betas";
    openStream(&ofile, outname, ios::out | ios::app, sout);
    ofile << ph + 1 << " ";
    ofile << beta_avg.transpose() << endl;
    ofile.close();
  }

  // sout << "\nFor tau[" << val <<"] = " << params.tau[ph](val) << endl <<  beta_l1 << endl ;
  int ctr = 0, chr_ctr = 0;
  int nn, cum_size_folds;

  for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {
    int chrom = files.chr_read[itr];
    if( !in_map(chrom, chr_map) ) continue;

    nn = chr_map[chrom][1] * params.n_ridge_l0 - l1_ests.chrom_map_ndiff(chrom-1);
    if(nn > 0) {
      cum_size_folds = 0;
      for(int i = 0; i < params.cv_folds; ++i ) {
        if(!params.within_sample_l0) beta_l1 = l1_ests.beta_hat_level_1[ph][i].col(val);
        predictions[0].block(cum_size_folds, chr_ctr, params.cv_sizes(i), 1) = l1_ests.test_mat[ph_eff][i].block(0, ctr, params.cv_sizes(i), nn) * beta_l1.block(ctr, 0, nn, 1);
        cum_size_folds += params.cv_sizes(i);
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


void Data::make_predictions_loocv(int const& ph, int const& val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();
  int ph_eff = params.write_l0_pred ? 0 : ph;

  // read in level 0 predictions from file
  if(params.write_l0_pred)
    read_l0(ph, ph_eff, &files, &params, &l1_ests, sout);
  check_l0(ph, ph_eff, &params, &l1_ests, &pheno_data, sout, true);

  int bs_l1 = l1_ests.test_mat_conc[ph_eff].cols();
  MatrixXd b0, xtx, tmpMat, HX_chunk;
  VectorXd zvec, bvec;
  ArrayXd calFactor, yres;

  uint64 max_bytes = params.chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (target_size * bs_l1) matrix ]
  int nchunk = ceil( params.cv_folds * bs_l1 * sizeof(double) * 1.0 / max_bytes );
  if (params.verbose) sout << nchunk << " chunks..." << flush;
  int chunk, size_chunk, target_size = params.cv_folds / nchunk;
  int j_start;

  xtx = l1_ests.test_mat_conc[ph_eff].transpose() * l1_ests.test_mat_conc[ph_eff];
  xtx.diagonal().array() += params.tau[ph](val) * l1_ests.ridge_param_mult;
  zvec = l1_ests.test_mat_conc[ph_eff].transpose() * pheno_data.phenotypes.col(ph);

  // fit model on whole data again for optimal ridge param
  SelfAdjointEigenSolver<MatrixXd> eigMat(xtx);
  tmpMat = eigMat.eigenvectors() * (1/eigMat.eigenvalues().array()).matrix().asDiagonal() * eigMat.eigenvectors().transpose();
  bvec = tmpMat * zvec;
  yres = (pheno_data.phenotypes.col(ph) - l1_ests.test_mat_conc[ph_eff] * bvec).array();

  for(chunk = 0; chunk < nchunk; ++chunk ) {
    size_chunk = chunk == nchunk - 1? params.cv_folds - target_size * chunk : target_size;
    j_start = chunk * target_size;

    HX_chunk = tmpMat * l1_ests.test_mat_conc[ph_eff].middleRows(j_start, size_chunk).transpose(); // k x Nc
    calFactor = (l1_ests.test_mat_conc[ph_eff].middleRows(j_start, size_chunk).array() * HX_chunk.transpose().array()).matrix().rowwise().sum().array();
    b0 = bvec.rowwise().replicate(size_chunk) - HX_chunk * (yres.segment(j_start, size_chunk)/(1-calFactor)).matrix().asDiagonal() ;

    int ctr = 0, chr_ctr = 0;
    int nn;

    for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {
      int chrom = files.chr_read[itr];
      if( !in_map(chrom, chr_map) ) continue;

      nn = chr_map[chrom][1] * params.n_ridge_l0 - l1_ests.chrom_map_ndiff(chrom-1);
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
void Data::make_predictions_binary(int const& ph, int const& val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();
  int ph_eff = params.write_l0_pred ? 0 : ph;

  // read in level 0 predictions from file
  if(params.write_l0_pred)
    read_l0(ph, ph_eff, &files, &params, &l1_ests, sout);
  check_l0(ph, ph_eff, &params, &l1_ests, &pheno_data, sout, true);

  int bs_l1 = l1_ests.test_mat[ph_eff][0].cols();
  ArrayXd etavec, pivec, wvec, zvec, score;
  MatrixXd betaold, betanew, XtW, XtWX, XtWZ;
  MatrixXd ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);

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
      XtWX.diagonal() += (params.tau[ph](val) * l1_ests.ridge_param_mult).matrix();
      betanew = XtWX.llt().solve(XtWZ);
      // compute score
      score = ArrayXd::Zero(betanew.rows());
      for(int i = 0; i < params.cv_folds; ++i ) {
        etavec = (l1_ests.test_offset[ph][i] + l1_ests.test_mat[ph_eff][i] * betanew).array();
        pivec = 1 - 1/(etavec.exp() + 1);
        score += (l1_ests.test_mat[ph_eff][i].transpose() * (l1_ests.test_pheno_raw[ph][i].array() - pivec).matrix()).array();
      }
      score -= params.tau[ph](val) * l1_ests.ridge_param_mult * betanew.array();

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

    nn = chr_map[chrom][1] * params.n_ridge_l0 - l1_ests.chrom_map_ndiff(chrom-1);
    if(nn > 0) {
      cum_size_folds = 0;
      for(int i = 0; i < params.cv_folds; ++i ) {
        if(!params.within_sample_l0) betanew = l1_ests.beta_hat_level_1[ph][i].col(val);
        predictions[0].block(cum_size_folds, chr_ctr, params.cv_sizes(i), 1) = l1_ests.test_mat[ph_eff][i].block(0, ctr, params.cv_sizes(i), nn) * betanew.block(ctr, 0, nn, 1);
        cum_size_folds += params.cv_sizes(i);
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

void Data::make_predictions_binary_loocv_full(int const& ph, int const& val) {

  sout << "  * making predictions (using all samples)..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();
  int ph_eff = params.write_l0_pred ? 0 : ph;

  // read in level 0 predictions from file
  if(params.write_l0_pred)
    read_l0(ph, ph_eff, &files, &params, &l1_ests, sout);
  check_l0(ph, ph_eff, &params, &l1_ests, &pheno_data, sout, true);

  int bs_l1 = l1_ests.test_mat_conc[ph_eff].cols();
  ArrayXd beta, pivec, wvec;
  MatrixXd XtWX, V1;

  uint64 max_bytes = params.chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (bs_l1 * target_size) matrix ]
  int nchunk = ceil( params.cv_folds * bs_l1 * sizeof(double) * 1.0 / max_bytes );
  int target_size = params.cv_folds / nchunk;

  MapArXd Y (pheno_data.phenotypes_raw.col(ph).data(), pheno_data.phenotypes_raw.rows());
  MapMatXd X (l1_ests.test_mat_conc[ph_eff].data(), pheno_data.phenotypes_raw.rows(), bs_l1);
  MapArXd offset (m_ests.offset_nullreg.col(ph).data(), pheno_data.phenotypes_raw.rows());
  MapArXb mask (pheno_data.masked_indivs.col(ph).data(), pheno_data.phenotypes_raw.rows());

  // fit logistic on whole data again for optimal ridge param
  beta = ArrayXd::Zero(bs_l1);
  run_log_ridge_loocv(params.tau[ph](val), l1_ests.ridge_param_mult, target_size, nchunk, beta, pivec, wvec, Y, X, offset, mask, &params, sout);

  // use estimates from this model directly
  // compute predictor for each chr
  int ctr = 0, chr_ctr = 0;
  int nn;

  for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {
    int chrom = files.chr_read[itr];
    if( !in_map(chrom, chr_map) ) continue;

    nn = chr_map[chrom][1] * params.n_ridge_l0 - l1_ests.chrom_map_ndiff(chrom-1);

    if(nn > 0) {
      predictions[0].col(chr_ctr) = l1_ests.test_mat_conc[ph_eff].middleCols(ctr, nn) * beta.segment(ctr, nn).matrix();
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

void Data::make_predictions_binary_loocv(int const& ph, int const& val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();
  int ph_eff = params.write_l0_pred ? 0 : ph;

  // read in level 0 predictions from file
  if(params.write_l0_pred)
    read_l0(ph, ph_eff, &files, &params, &l1_ests, sout);
  check_l0(ph, ph_eff, &params, &l1_ests, &pheno_data, sout, true);

  int bs_l1 = l1_ests.test_mat_conc[ph_eff].cols();
  ArrayXd beta, pivec, wvec, v2;
  MatrixXd XtWX, V1, beta_final;
  LLT<MatrixXd> Hinv;

  uint64 max_bytes = params.chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (bs_l1 * target_size) matrix ]
  int nchunk = ceil( params.cv_folds * bs_l1 * sizeof(double) * 1.0 / max_bytes );
  int chunk, size_chunk, target_size = params.cv_folds / nchunk;
  int j_start;

  MapArXd Y (pheno_data.phenotypes_raw.col(ph).data(), pheno_data.phenotypes_raw.rows());
  MapMatXd X (l1_ests.test_mat_conc[ph_eff].data(), pheno_data.phenotypes_raw.rows(), bs_l1);
  MapArXd offset (m_ests.offset_nullreg.col(ph).data(), pheno_data.phenotypes_raw.rows());
  MapArXb mask (pheno_data.masked_indivs.col(ph).data(), pheno_data.phenotypes_raw.rows());

  // fit logistic on whole data again for optimal ridge param
  beta = ArrayXd::Zero(bs_l1);
  run_log_ridge_loocv(params.tau[ph](val), l1_ests.ridge_param_mult, target_size, nchunk, beta, pivec, wvec, Y, X, offset, mask, &params, sout);

  // compute Hinv
  //zvec = (etavec - m_ests.offset_nullreg.col(ph).array()) + (pheno_data.phenotypes_raw.col(ph).array() - pivec) / wvec;
  XtWX = (params.tau[ph](val) * l1_ests.ridge_param_mult).matrix().asDiagonal(); // compute XtWX in chunks
  for(chunk = 0; chunk < nchunk; ++chunk){
    size_chunk = ( chunk == nchunk - 1 ? params.cv_folds - target_size * chunk : target_size );
    j_start = chunk * target_size;

    Ref<MatrixXd> Xmat_chunk = X.middleRows(j_start, size_chunk); // n x k
    Ref<ArrayXd> w_chunk = wvec.segment(j_start, size_chunk);
    Ref<ArrayXb> mask_chunk = mask.segment(j_start, size_chunk);

    XtWX.noalias() += Xmat_chunk.transpose() * mask_chunk.select(w_chunk,0).matrix().asDiagonal() * Xmat_chunk;
  }
  Hinv.compute( XtWX );

  // loo estimates
  beta_final = MatrixXd::Zero(bs_l1, target_size);
  for(chunk = 0; chunk < nchunk; ++chunk ) {
    size_chunk = chunk == nchunk - 1? params.cv_folds - target_size * chunk : target_size;
    j_start = chunk * target_size;
    if( chunk == (nchunk - 1) ) beta_final = MatrixXd::Zero(bs_l1, size_chunk);

    Ref<MatrixXd> Xmat_chunk = X.middleRows(j_start, size_chunk); // n x k
    Ref<ArrayXd> Yvec_chunk = Y.segment(j_start, size_chunk);
    Ref<ArrayXd> p_chunk = pivec.segment(j_start, size_chunk);
    Ref<ArrayXd> w_chunk = wvec.segment(j_start, size_chunk);

    V1 = Hinv.solve( Xmat_chunk.transpose() ); // k x n
    v2 = (Xmat_chunk.array() * V1.transpose().array()).rowwise().sum() * w_chunk;
    beta_final.array().colwise() = beta;
    beta_final -= V1 * ((Yvec_chunk - p_chunk)/(1-v2)).matrix().asDiagonal();

    // compute predictor for each chr
    int ctr = 0, chr_ctr = 0;
    int nn;

    for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {
      int chrom = files.chr_read[itr];
      if( !in_map(chrom, chr_map) ) continue;

      nn = chr_map[chrom][1] * params.n_ridge_l0 - l1_ests.chrom_map_ndiff(chrom-1);

      if(nn > 0) {
        predictions[0].block(j_start, chr_ctr, size_chunk, 1) = ( X.block(j_start, ctr, size_chunk, nn).array() * beta_final.middleRows(ctr, nn).transpose().array() ).matrix().rowwise().sum();
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


// predictions for count phenotypes
void Data::make_predictions_count(int const& ph, int const& val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();
  int ph_eff = params.write_l0_pred ? 0 : ph;

  // read in level 0 predictions from file
  if(params.write_l0_pred)
    read_l0(ph, ph_eff, &files, &params, &l1_ests, sout);
  check_l0(ph, ph_eff, &params, &l1_ests, &pheno_data, sout, true);

  int bs_l1 = l1_ests.test_mat[ph_eff][0].cols();
  ArrayXd etavec, pivec, zvec, score;
  MatrixXd betaold, betanew, XtW, XtWX, XtWZ;
  MatrixXd ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);

  // fit model using out-of-sample level 0 predictions from whole data
  if(params.within_sample_l0)
    throw "--within is not supported for count phenotypes";

  // compute predictor for each chr
  int ctr = 0, chr_ctr = 0;
  int nn, cum_size_folds;

  for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {
    int chrom = files.chr_read[itr];
    if( !in_map(chrom, chr_map) ) continue;

    nn = chr_map[chrom][1] * params.n_ridge_l0 - l1_ests.chrom_map_ndiff(chrom-1);
    if(nn > 0) {
      cum_size_folds = 0;
      for(int i = 0; i < params.cv_folds; ++i ) {
        betanew = l1_ests.beta_hat_level_1[ph][i].col(val);
        predictions[0].block(cum_size_folds, chr_ctr, params.cv_sizes(i), 1) = l1_ests.test_mat[ph_eff][i].block(0, ctr, params.cv_sizes(i), nn) * betanew.block(ctr, 0, nn, 1);
        cum_size_folds += params.cv_sizes(i);
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


void Data::make_predictions_count_loocv(int const& ph, int const& val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();
  int ph_eff = params.write_l0_pred ? 0 : ph;

  if(params.write_l0_pred)
    read_l0(ph, ph_eff, &files, &params, &l1_ests, sout);
  check_l0(ph, ph_eff, &params, &l1_ests, &pheno_data, sout, true);

  int bs_l1 = l1_ests.test_mat_conc[ph_eff].cols();
  double v2;
  ArrayXd beta, pivec;
  MatrixXd XtWX, V1, beta_final;
  LLT<MatrixXd> Hinv;
  MatrixXd ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);

  uint64 max_bytes = params.chunk_mb * 1e6;
  // amount of RAM used < max_mb [ creating (bs_l1 * target_size) matrix ]
  int nchunk = ceil( params.cv_folds * bs_l1 * sizeof(double) * 1.0 / max_bytes );
  int chunk, size_chunk, target_size = params.cv_folds / nchunk;
  int j_start;

  // read in level 0 predictions from file
  MapArXd Y (pheno_data.phenotypes_raw.col(ph).data(), pheno_data.phenotypes_raw.rows());
  MapMatXd X (l1_ests.test_mat_conc[ph_eff].data(), pheno_data.phenotypes_raw.rows(), bs_l1);
  MapArXd offset (m_ests.offset_nullreg.col(ph).data(), pheno_data.phenotypes_raw.rows());
  MapArXb mask (pheno_data.masked_indivs.col(ph).data(), pheno_data.phenotypes_raw.rows());

  // fit logistic on whole data again for optimal ridge param
  beta = ArrayXd::Zero(bs_l1);
  run_ct_ridge_loocv(params.tau[ph](val), l1_ests.ridge_param_mult, target_size, nchunk, beta, pivec, Y, X, offset, mask, &params, sout);

  // compute Hinv
  //zvec = (etavec - m_ests.offset_nullreg.col(ph).array()) + (pheno_data.phenotypes_raw.col(ph).array() - pivec) / wvec;
  XtWX = (params.tau[ph](val) * l1_ests.ridge_param_mult).matrix().asDiagonal();
  for(chunk = 0; chunk < nchunk; ++chunk){
    size_chunk = ( chunk == nchunk - 1 ? params.cv_folds - target_size * chunk : target_size );
    j_start = chunk * target_size;

    Ref<MatrixXd> Xmat_chunk = X.block(j_start, 0, size_chunk, bs_l1); // n x k
    Ref<MatrixXd> w_chunk = pivec.matrix().block(j_start, 0, size_chunk,1);

    XtWX += Xmat_chunk.transpose() * w_chunk.asDiagonal() * Xmat_chunk;
  }
  Hinv.compute( XtWX );

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
      v2 *= pivec(j_start + i);
      beta_final.col(i) = (beta - V1.col(i).array() * (Yvec_chunk(i,0) - pivec(j_start + i)) / (1 - v2)).matrix();
    }

    // compute predictor for each chr
    int ctr = 0, chr_ctr = 0;
    int nn;

    for (size_t itr = 0; itr < files.chr_read.size(); ++itr) {
      int chrom = files.chr_read[itr];
      if( !in_map(chrom, chr_map) ) continue;

      nn = chr_map[chrom][1] * params.n_ridge_l0 - l1_ests.chrom_map_ndiff(chrom-1);

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


void Data::write_predictions(int const& ph){
  // output predictions to file
  string out, header;
  Files ofile;
  MatrixXd pred, prs;

  // get header line once
  if(params.write_blups || params.make_loco || params.trait_mode || params.print_prs)
    header = write_ID_header();

  // for the per chromosome predictions (not used)
  if(params.write_blups) {

    out = files.out_file + "_" + to_string(ph+1) + (params.gzOut ? ".gz" : "");
    sout << "writing file " << out << "..." << flush;
    ofile.openForWrite(out, sout);

    // enforce all chromosomes are printed
    pred = MatrixXd::Zero(predictions[0].rows(), params.nChrom);

    int nn, chr_ctr = 0;
    for(auto const& chr : files.chr_read){

      if( !in_map(chr, chr_map) ) continue;

      nn = chr_map[chr][1];
      if(nn > 0){
        pred.col(chr - 1) = predictions[0].col(chr_ctr);
        ++chr_ctr;
      }

    }

    // header line : FID_IID for all individuals
    ofile << header;

    // for each row: print chromosome then blups
    for(int chr = 0; chr < params.nChrom; chr++) 
      ofile << write_chr_row(chr+1, ph, pred.col(chr));

    ofile.closeFile();

  }

  // output LOCO predictions G_loco * beta_loco for each autosomal chr
  if(params.make_loco || params.trait_mode){

    out = files.out_file + "_" + to_string(ph+1) + ".loco" + (params.gzOut ? ".gz" : "");
    sout << "writing LOCO predictions..." << flush;
    ofile.openForWrite(out, sout);

    pred.resize(predictions[0].rows(), params.nChrom);
    pred.colwise() = predictions[0].rowwise().sum();

    int nn, chr_ctr = 0;
    for(auto const& chr : files.chr_read){

      if( !in_map(chr, chr_map) ) continue;

      nn = chr_map[chr][1];
      if(nn > 0) {
        pred.col(chr - 1) -= predictions[0].col(chr_ctr);
        ++chr_ctr;
      }

    }

    // header line : FID_IID for all individuals
    ofile << header;

    // print loco predictions for each chromosome
    for(int chr = 0; chr < params.nChrom; chr++) 
      ofile << write_chr_row(chr+1, ph, pred.col(chr));

    ofile.closeFile();

  }

  if(params.write_null_firth){ // store null estimates for Firth

    bool has_converged = true;
    double dev, tstat;
    ArrayXd se, etavec, pivec;
    IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","");

    out = files.out_file + "_" + to_string(ph+1) + ".firth" + (params.gzOut ? ".gz" : "");
    sout << "writing null approximate Firth estimates..." << flush;
    ofile.openForWrite(out, sout);

    ArrayXd bhat = ArrayXd::Zero(params.ncov);
    MapArXd Y (pheno_data.phenotypes_raw.col(ph).data(), pheno_data.phenotypes_raw.rows());
    MapArXb mask (pheno_data.masked_indivs.col(ph).data(), pheno_data.masked_indivs.rows());

    for(int chr = 0; chr < params.nChrom; chr++) {
      // fit null approximate Firth 
      // use warm starts from previous chromosomes
      if(!fit_firth(ph, Y, pheno_data.new_cov, pred.col(chr).array(), mask, pivec, etavec, bhat, se, params.ncov, dev, false, tstat, params.maxstep_null, params.niter_max_firth_null, 50 * params.numtol, &params)){
        has_converged = false;
        break;
      }
      ofile << chr + 1 << " " << bhat.matrix().transpose().format(Fmt) << endl;
    }

    if(!has_converged){ // remove file
      sout << "WARNING: Firth failed to converge";
      remove(out.c_str());
    } else
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
    ofile << header;

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

    // check individual was included in analysis, if not then skip
    index = itr_ind->second;
    if( !in_filters.ind_in_analysis( index ) ) continue;

    id_index = itr_ind->first;
    buffer << id_index << " ";

  }

  buffer << endl;
  return buffer.str();

}


std::string Data::write_chr_row(int const& chr, int const& ph, const Ref<const VectorXd>& pred){

  uint32_t index;
  string out;
  std::ostringstream buffer;
  map<string, uint32_t >::iterator itr_ind;

  buffer << chr << " ";
  for (itr_ind = params.FID_IID_to_ind.begin(); itr_ind != params.FID_IID_to_ind.end(); ++itr_ind) {

    // check individual was included in analysis, if not then skip
    index = itr_ind->second;
    if( !in_filters.ind_in_analysis( index ) ) continue;

    // print prs
    if( pheno_data.masked_indivs(index, ph) )
      buffer << pred(index) << " ";
    else
      buffer << "NA ";
  }

  buffer << endl;
  return buffer.str();

}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Functions needed in testing mode
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::setup_output(Files* ofile, string& out, std::vector<std::shared_ptr<Files>>& ofile_split, std::vector< string >& out_split){

  if(params.getCorMat){ // header N,M
    out = files.out_file + ".corr";
    string runmode = (params.dosage_mode ? "in dosage mode" : "in hard-call mode");
    if(params.cor_out_txt){
      sout << " * computing correlation matrix " + runmode + "\n  + output to text file ["<<out<<"]\n";
      ofile->openForWrite(out, sout);
      if(params.skip_scaleG) (*ofile) << params.extract_vars_order.size() << " " << params.n_samples << "\n";
    } else {
      sout << " * computing correlation matrix " + runmode + " (storing R^2 values)\n  + output to binary file ["<<out<<"]\n";
      ofile->openMode(out, std::ios_base::out | std::ios_base::binary, sout);
      ArrayXi vals(2);
      vals << params.n_samples, params.extract_vars_order.size();
      //cerr << vals << endl;
      ofile->writeBinMode(vals, sout);
    }
    sout << "  + list of snps written to [" << out << ".snplist]\n";
    sout << "  + n_snps = " << params.extract_vars_order.size() <<"\n\n";
    return;
  }

  // info for output file
  string tmpstr = (params.htp_out ? print_header_output_htp() : print_header_output(&params));
  string mask_header = (params.build_mask ? bm.build_header() : "");
  string gz_ext = (params.gzOut ? ".gz" : "");

  if( !params.split_by_pheno ){
    out = files.out_file + ".regenie" + gz_ext;
    ofile->openForWrite( out, sout );
    (*ofile) << mask_header << tmpstr;
    // write dictionary file
    Files dict_file;
    dict_file.openForWrite( files.out_file + ".regenie.Ydict", sout);
    for(int i = 0; i < params.n_pheno; i++) 
      dict_file << "Y" << i+1 << " " << files.pheno_names[i] << endl;
    dict_file.closeFile();
    return;
  }

  out_split.resize( params.n_pheno );
  ofile_split.resize( params.n_pheno );

  for(int i = 0; i < params.n_pheno; i++) {
    out_split[i] = files.out_file + "_" + files.pheno_names[i] + ".regenie" + gz_ext;
    ofile_split[i] = std::make_shared<Files>();
    ofile_split[i]->openForWrite( out_split[i], sout );
    (*ofile_split[i]) << mask_header << tmpstr;
  }

}

void Data::print_test_info(){

  if(params.getCorMat) { params.with_flip = false; return; }

  if(params.write_masks) {

    bm.write_info(&params, &in_filters, sout);
    sout << " * user specified to write masks (in PLINK bed format)\n";
    if(params.dosage_mode) 
      sout << "   +dosages will be converted to hardcalls\n";
    if(params.write_setlist) 
      bm.prep_setlists(files.new_sets, files.out_file, sout);

  }
  if(params.write_mask_snplist) 
    bm.prep_snplist(files.out_file, sout);

  sout << " * using minimum MAC of " << (params.build_mask ? params.min_MAC_mask : params.min_MAC) << 
    " (" << (params.build_mask ? "masks" : "variants") << " with lower MAC are ignored)\n";
  if(params.forced_MAC > 0) sout << "   -using threshold of " << params.forced_MAC << " for subset of specified variants\n";
  if(params.setMinINFO) 
    sout << " * using minimum imputation info score of " << params.min_INFO << " (variants with lower info score are ignored)\n";
  if((params.test_type == 2) && (params.minHOMs > 0))
    sout << " * ignoring variants (masks for gene-based tests) with fewer than " << params.minHOMs << " homALT carriers\n";

  if(params.firth || params.use_SPA) {

    sout << " * using " << (params.firth_approx ? "fast ": "") << (params.firth ? "Firth ": "SPA ");
    sout << "correction for logistic regression p-values less than " << params.alpha_pvalue << endl;
    if(params.back_correct_se) sout << "    - using back-correction to compute Firth SE\n";
    if(params.firth && params.use_adam) sout << "    - using " << (params.adam_mini? "mini-":"") << "batch ADAM to get starting values\n";
    n_corrected = 0;

  }

  // if testing select chromosomes
  if( params.select_chrs ) 
    sout << " * user specified to test only on select chromosomes\n";


  switch(params.test_type){
    case 0:
      test_string = "ADD";
      break;
    case 1:
      test_string = "DOM";
      break;
    case 2:
      test_string = "REC";
      break;
    default:
      throw "unrecognized test value";
  }
  wgr_string = ( params.skip_blups && !params.interaction_prs && !params.blup_cov ?  "" : "-WGR" );


  if(params.htp_out){
    if((params.trait_mode==1) & params.firth) correction_type = "-FIRTH";
    else if((params.trait_mode==1) & params.use_SPA) correction_type = "-SPA";
    else if(params.trait_mode==1) correction_type = "-LOG";
    else if(params.trait_mode==2) correction_type = "-POISSON";
    else correction_type = "-LR";

    model_type = test_string + wgr_string + correction_type;
  }

  if(params.gwas_condtl) // specify main sum stats is conditional gwas
    params.condtl_suff = "-CONDTL";

  params.with_flip = params.with_flip && !params.build_mask && params.trait_mode && (params.test_type == 0);

  if( params.joint_test ) {
    jt.out_file_prefix = files.out_file;
    params.with_flip = jt.get_test_info(&params, test_string, sout) && params.with_flip;
    sout << " * list of joint tests run on burden masks: " << get_test_list(jt.test_list, jt.joint_tests_map) << "\n";
  }

  normal nd(0,1);
  chi_squared chisq(1);
  params.zcrit = quantile(complement(nd, .025));
  params.chisq_thr = quantile(chisq, 1 - params.alpha_pvalue);
  params.z_thr = sqrt(params.chisq_thr);

}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    prep for association test
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::set_blocks_for_testing() {

  params.total_n_block = 0;
  int blocks_left = params.n_block;
  int nchr = 0;

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
    m1[ itr->first ] = itr->second;
  }
  chr_map = m1;

  if(params.getCorMat && (nchr > 1))
    throw "can only compute LD matrix for a single chromosome (use --chr/--chrList/--range).";

  // summarize block sizes
  sout << left << std::setw(20) << " * # threads" << ": [" << params.threads << "]\n";
  sout << left << std::setw(20) << " * block size" << ": [" << params.block_size << "]\n";
  if(!params.getCorMat) {
    sout << left << std::setw(20) << " * # blocks" << ": [" << params.total_n_block << "]\n";
    if(params.start_block > 1) sout << "    + skipping to block #" << params.start_block << endl;
  }

  // storing null estimates from firth
  if(params.use_null_firth) 
    sout << " * reading null Firth estimates using file : [" << files.null_firth_file << "]\n";
  if(params.write_null_firth ) 
    sout << " * writing null Firth estimates to file\n";

}

void Data::set_nullreg_mat(){

  m_ests.Y_hat_p = MatrixXd::Zero(params.n_samples, params.n_pheno);
  m_ests.Gamma_sqrt = MatrixXd::Zero(params.n_samples, params.n_pheno);
  m_ests.X_Gamma.resize(params.n_pheno);

  // for firth  approx
  if(params.firth_approx){
    firth_est.beta_null_firth = MatrixXd::Zero(pheno_data.new_cov.cols() + 1, params.n_pheno);
    if(params.test_mode) firth_est.cov_blup_offset = MatrixXd::Zero(params.n_samples, params.n_pheno);

    // open streams to write firth null estimates
    if(params.write_null_firth){
      string fout_p;
      firth_est.firth_est_files.resize(params.n_pheno);
      for(int ph = 0; ph < params.n_pheno; ph++){
        firth_est.firth_est_files[ph] = std::make_shared<Files>();
        fout_p = files.out_file + "_" + to_string(ph + 1) + ".firth" + (params.gzOut ? ".gz" : "");
        firth_est.firth_est_files[ph]->openForWrite(fout_p, sout);
      }
      if(params.compute_all_chr){
        params.use_null_firth = false; // make sure nothing is read
        files.null_firth_file = get_firth_est_allChr(files, in_filters, m_ests, firth_est, pheno_data, params, sout);;
        params.write_null_firth = false;
        params.use_null_firth = true;
      }
    }


    if(params.use_null_firth) // get files with firth estimates
      check_beta_start_firth(files, params, sout);
  }
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Testing mode (multi-threaded with OpenMP)
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::test_snps_fast() {

  sout << "Association testing mode";

  string out;
  vector < string > out_split;
  // output files
  Files ofile;
  // use pointer to class since it contains non-copyable elements
  vector < std::shared_ptr<Files> > ofile_split;

#if defined(_OPENMP)
  sout << " with " << (params.streamBGEN? "fast " : "") << "multithreading using OpenMP";
#endif
  sout << endl;

  file_read_initialization(); // set up files for reading
  read_pheno_and_cov(&files, &params, &in_filters, &pheno_data, &m_ests, &Gblock, sout);   // read phenotype and covariate files
  prep_run(&files, &in_filters, &params, &pheno_data, &m_ests, sout); // check blup files and adjust for covariates
  set_blocks_for_testing();   // set number of blocks
  print_usage_info(&params, &files, sout);
  print_test_info();
  setup_output(&ofile, out, ofile_split, out_split); // result file
  if(params.w_interaction && (params.trait_mode==0) && !params.no_robust && !params.force_robust) 
    nullHLM.prep_run(&pheno_data, &params);
  if(params.trait_mode) set_nullreg_mat();
  sout << endl;

  // start analyzing each chromosome
  bool block_init_pass = false;
  int block = 0, chrom_nsnps, chrom_nb, bs;
  tally snp_tally;
  vector< variant_block > block_info;
  initialize_thread_data(Gblock.thread_data, params);


  for(auto const& chrom : files.chr_read) {

    if( !in_map(chrom, chr_map) ) continue;

    chrom_nsnps = chr_map[chrom][0];
    chrom_nb = chr_map[chrom][1];
    if(chrom_nb == 0) continue;

    // If specified starting block
    if(!block_init_pass && (params.start_block > (block + chrom_nb)) ) {
      snp_tally.snp_count += chrom_nsnps;
      block += chrom_nb;
      continue;
    }

    sout << "Chromosome " << chrom << " [" << chrom_nb << " blocks in total]\n";

    if(!params.getCorMat){
      // read polygenic effect predictions from step 1
      blup_read_chr(false, chrom, m_ests, files, in_filters, pheno_data, params, sout);

      // compute phenotype residual (adjusting for BLUP [and covariates for non-QTs])
      if(params.trait_mode == 1) compute_res_bin(chrom);
      else if(params.trait_mode == 2) compute_res_count(chrom);
      else compute_res();

      // print y/x/logreg offset used for level 1 
      if(params.debug) write_inputs();
    }

    // analyze by blocks of SNPs
    for(int bb = 0; bb < chrom_nb ; bb++) {

      get_block_size(params.block_size, chrom_nsnps, bb, bs);

      // If specified starting block
      if(!block_init_pass && (params.start_block > (block+1)) ) {
        snp_tally.snp_count += bs;
        block++;
        continue;
      } else if(!block_init_pass) block_init_pass = true;

      sout << " block [" << block + 1 << "/" << params.total_n_block << "] : " << flush;

      allocate_mat(Gblock.Gmat, params.n_samples, bs);
      block_info.resize(bs);

      // read SNP, impute missing & compute association test statistic
      analyze_block(chrom, bs, &snp_tally, block_info);

      // print the results
      for (auto const& snp_data : block_info){

        if( snp_data.ignored ) {
          snp_tally.n_ignored_snps++;
          continue;
        }

        snp_tally.n_ignored_tests += snp_data.ignored_trait.count();
        if(params.firth || params.use_SPA) {
          n_corrected += (!snp_data.ignored_trait && snp_data.is_corrected).count();
          snp_tally.n_failed_tests += (!snp_data.ignored_trait && snp_data.test_fail).count();
          if(params.w_interaction) {
            n_corrected += (2 + params.ncov_interaction) * snp_data.is_corrected_inter.count(); // main, inter & joint
            snp_tally.n_failed_tests += (2 + params.ncov_interaction) * (snp_data.is_corrected_inter && snp_data.test_fail_inter).count(); // main, inter & joint
          }
        }

        for(int j = 0; j < params.n_pheno; ++j) {

          if( !params.pheno_pass(j) || snp_data.ignored_trait(j) ) {
            if(!params.split_by_pheno) // if using single file, print NAs for snp/trait sum stats
              ofile << snp_data.sum_stats[j];

            continue;
          }

          if(params.split_by_pheno)
            (*ofile_split[j]) << snp_data.sum_stats[j]; // add test info
          else
            ofile << snp_data.sum_stats[j]; // add test info
        }

      }

      snp_tally.snp_count += bs;
      block++;
    }

  }

  sout << print_summary(&ofile, out, ofile_split, out_split, n_corrected, snp_tally, files, firth_est, params);

}

// test SNPs in block
void Data::analyze_block(int const& chrom, int const& n_snps, tally* snp_tally, vector<variant_block> &all_snps_info){

  auto t1 = std::chrono::high_resolution_clock::now();
  const int start = snp_tally->snp_count;
  vector< vector < uchar > > snp_data_blocks;
  vector< uint32_t > insize, outsize;

  vector<uint64> indices(n_snps);
  std::iota(indices.begin(), indices.end(), start);

  readChunk(indices, chrom, snp_data_blocks, insize, outsize, all_snps_info);

  // analyze using openmp
  compute_tests_mt(chrom, indices, snp_data_blocks, insize, outsize, all_snps_info);
  //compute_tests_st(chrom, indices, snp_data_blocks, insize, outsize, all_snps_info); // this is slower


  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << "done (" << duration.count() << "ms) "<< endl;
}


void Data::compute_res(){

  if(params.blup_cov) // blup as covariate
    get_lm_resid(res, m_ests.blups, pheno_data.phenotypes); 
  else res = pheno_data.phenotypes - m_ests.blups;
  res.array() *= pheno_data.masked_indivs.array().cast<double>();

  if(params.rerint | params.rerintcov) {
    residualize_res();
  }

  p_sd_yres = res.colwise().norm();
  p_sd_yres.array() /= sqrt(pheno_data.Neff - params.ncov_analyzed); // if blup is cov
  res.array().rowwise() /= p_sd_yres.array();

  if(params.w_interaction && (params.trait_mode==0) && !params.no_robust && !params.force_robust) 
    HLM_fitNull(nullHLM, m_ests, pheno_data, files, params, sout);
}

// two-stage rinting, as described in Sofer et al., 2020, ttps://www.ncbi.nlm.nih.gov/pmc/articles/PMC6416071/
// --apply-rerint = RN-Resid-Unadj (Table 1, Sofer et al., 2020)
// --apply-rerint-cov = RN-Resid-Adj (Table 1, Sofer et al., 2020)
void Data::residualize_res() {
  // for each residual, apply rank-inverse normal transformation
  for(int ph = 0; ph < res.cols(); ph++) {
    rint_pheno(res.col(ph), pheno_data.masked_indivs.col(ph).array());
  }

  // further project covariates out from residuals
  if(params.rerintcov) {
    MatrixXd beta = res.transpose() * pheno_data.new_cov;
    res -= ( (pheno_data.new_cov * beta.transpose()).array() * pheno_data.masked_indivs.array().cast<double>() ).matrix();
  }

  // respect masked individuals (needed here?) 
  res.array() *= pheno_data.masked_indivs.array().cast<double>();

  // performa scaling of reisidualized phenotypes, similar to residualize_phenotypes in Pheno.cpp
  // compute the scale of residuals 
  pheno_data.scale_Y = res.colwise().norm().array() / sqrt(pheno_data.Neff.matrix().transpose().array() - params.ncov_analyzed);
  // set sd for phenotypes which are ignored to 1
  pheno_data.scale_Y = params.pheno_pass.select(pheno_data.scale_Y.transpose().array(), 1).matrix().transpose();
  // check sd is not 0 
  MatrixXd::Index minIndex;
  if(pheno_data.scale_Y.minCoeff(&minIndex) < params.numtol)
    throw "some phenotype residuals has sd=0.";
  // scale residuals
  res.array().rowwise() /= pheno_data.scale_Y.array();
}


void Data::compute_res_bin(int const& chrom){

  fit_null_logistic(false, chrom, &params, &pheno_data, &m_ests, &files, sout); // for all phenotypes

  res = pheno_data.phenotypes_raw - m_ests.Y_hat_p;
  res.array() /= m_ests.Gamma_sqrt.array();
  res.array() *= pheno_data.masked_indivs.array().cast<double>();

  // if using firth approximation, fit null penalized model with only covariates and store the estimates (to be used as offset when computing LRT in full model)
  if(params.firth_approx) fit_null_firth(false, chrom, &firth_est, &pheno_data, &m_ests, &files, &params, sout);

}

void Data::compute_res_count(int const& chrom){

  fit_null_poisson(chrom, &params, &pheno_data, &m_ests, &files, sout); // for all phenotypes

  res = pheno_data.phenotypes_raw - m_ests.Y_hat_p;
  res.array() /= m_ests.Gamma_sqrt.array();
  res.array() *= pheno_data.masked_indivs.array().cast<double>();

}

void Data::compute_tests_mt(int const& chrom, vector<uint64> indices,vector< vector < uchar > >& snp_data_blocks, vector< uint32_t > insize, vector< uint32_t >& outsize, vector<variant_block> &all_snps_info){
  
  size_t const bs = indices.size();
  ArrayXb err_caught = ArrayXb::Constant(bs, false);

  if( !params.build_mask && (((params.file_type == "bgen") && params.streamBGEN) || params.file_type == "bed") ) {
    // start openmp for loop
#if defined(_OPENMP)
    setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
    for(size_t isnp = 0; isnp < bs; isnp++) {
      uint32_t const snp_index = indices[isnp];
      // to store variant information
      variant_block* block_info = &(all_snps_info[isnp]);
      parseSNP(isnp, chrom, &(snp_data_blocks[isnp]), insize[isnp], outsize[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &snpinfo[snp_index], &Gblock, block_info, sout);
    }
#if defined(_OPENMP)
    setNbThreads(params.threads);
#endif
  }

  // for QTs: project out covariates & scale
  // do all at once for all variants
  MatrixXd Gtmp_res;
  int neff = residualize_gmat(false, pheno_data.new_cov, Gblock.Gmat, Gtmp_res, params);

  // start openmp for loop
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(size_t isnp = 0; isnp < bs; isnp++) {
    uint32_t const snp_index = indices[isnp];

    int thread_num = 0;
#if defined(_OPENMP)
    thread_num = omp_get_thread_num();
#endif

    // to store variant information
    variant_block* block_info = &(all_snps_info[isnp]);
    reset_thread(&(Gblock.thread_data[thread_num]), params);

    // check if g is sparse (not for QT without strict mode)
    if(params.trait_mode || params.strict_mode)
      check_sparse_G(isnp, thread_num, &Gblock, params.n_samples, in_filters.ind_in_analysis);

    if(params.w_interaction) {
      if(params.interaction_snp && (snpinfo[snp_index].ID == in_filters.interaction_cov))
        block_info->skip_int = true;
      get_interaction_terms(isnp, thread_num, &pheno_data, &Gblock, block_info, nullHLM, &params, sout);
    }

    // for QTs: project out covariates & scale
    check_res_geno(isnp, thread_num, block_info, false, neff, Gtmp_res, &Gblock, params);

    // skip SNP if fails filters
    if( block_info->ignored || params.getCorMat ) continue;
    
    reset_stats(block_info, params);

    try {
      // if ran vc tests, print out results before mask test
      if((block_info->sum_stats_vc.size() > 0) && !params.p_joint_only)
        print_vc_sumstats(snp_index, "ADD", wgr_string, block_info, snpinfo, files, &params);

      compute_score(isnp, snp_index, chrom, thread_num, test_string + params.condtl_suff, model_type + params.condtl_suff, res, p_sd_yres, params, pheno_data, Gblock, block_info, snpinfo, m_ests, firth_est, files, sout);

      // for joint test, store logp
      if( params.joint_test ) block_info->pval_log = Gblock.thread_data[thread_num].pval_log;

      if(params.w_interaction) 
        apply_interaction_tests(snp_index, isnp, thread_num, res, p_sd_yres, model_type, test_string, &pheno_data, nullHLM, &in_filters, &files, &Gblock, block_info, snpinfo, &m_ests, &firth_est, &params, sout);
    } catch (...) {
      err_caught(isnp) = true;
      block_info->sum_stats[0] = boost::current_exception_diagnostic_information();
      continue;
    }
  }

#if defined(_OPENMP)
  setNbThreads(params.threads);
#endif

  // check no errors
  if(err_caught.any())
    for(int i = 0; i < err_caught.size(); i++)
      if(err_caught(i)) throw all_snps_info[i].sum_stats[0];

}

/*
void Data::compute_tests_st(int const& chrom, vector<uint64> indices,vector< vector < uchar > >& snp_data_blocks, vector< uint32_t > insize, vector< uint32_t >& outsize, vector<variant_block> &all_snps_info){

  size_t const bs = indices.size();

  // start openmp for loop
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(size_t isnp = 0; isnp < bs; isnp++) {
    uint32_t const snp_index = indices[isnp];

    // to store variant information
    variant_block* block_info = &(all_snps_info[isnp]);

    if( !params.build_mask )
      parseSNP(isnp, chrom, &(snp_data_blocks[isnp]), insize[isnp], outsize[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &snpinfo[snp_index], &Gblock, block_info, sout);

    if(params.w_interaction) {
      if(params.interaction_snp && (snpinfo[snp_index].ID == in_filters.interaction_cov))
        block_info->skip_int = true;
      get_interaction_terms(isnp, &pheno_data, &Gblock, block_info, nullHLM, &params, sout);
    }

    // for QTs (or BTs with firth approx): project out covariates & scale
    residualize_geno(isnp, block_info, false, pheno_data.new_cov, &Gblock, &params);

    // skip SNP if fails filters
    if( block_info->ignored || params.getCorMat ) continue;

    reset_stats(block_info, params);
  }

#if defined(_OPENMP)
  setNbThreads(params.threads);
#endif

  int npass = 0;
  for(size_t isnp = 0; (isnp < bs) && (npass == 0); isnp++)
    if(!all_snps_info[isnp].ignored) npass++;

  // skip block if all fails filters
  if( (npass == 0) || params.getCorMat ) return;

  compute_score(indices, chrom, test_string, model_type, res, p_sd_yres, params, pheno_data, Gblock, all_snps_info, snpinfo, m_ests, firth_est, files, sout);

  //if( (isnp==0) || (isnp == (n_snps-1)) ) cout << "G"<<isnp+1<<" MAF = " <<  block_info.MAF << endl;

  if(params.w_interaction) {
    // start openmp for loop
#if defined(_OPENMP)
    setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
    for(size_t isnp = 0; isnp < bs; isnp++) {
      uint32_t const snp_index = indices[isnp];
      variant_block* block_info = &(all_snps_info[isnp]);
      apply_interaction_tests(snp_index, isnp, res, p_sd_yres, model_type, test_string, &pheno_data, nullHLM, &in_filters, &files, &Gblock, block_info, snpinfo, &m_ests, &firth_est, &params, sout);
    }
#if defined(_OPENMP)
    setNbThreads(params.threads);
#endif
  }

}
*/

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Testing mode (joint tests)
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::test_joint() {

  
  sout << "Association testing mode (joint tests)";

  std::chrono::high_resolution_clock::time_point t1, t2;
  string out;
  vector < string > out_split, tmp_str;
  // output files
  Files ofile;
  // use pointer to class since it contains non-copyable elements
  vector < std::shared_ptr<Files> > ofile_split;

  // set some parameters
  if( params.build_mask ) bm.prep_run(params, files);

#if defined(_OPENMP)
  sout << " with " << (params.streamBGEN? "fast " : "") << "multithreading using OpenMP";
#endif
  sout << endl;

  file_read_initialization(); // set up files for reading
  read_pheno_and_cov(&files, &params, &in_filters, &pheno_data, &m_ests, &Gblock, sout);   // read phenotype and covariate files
  prep_run(&files, &in_filters, &params, &pheno_data, &m_ests, sout); // check blup files and adjust for covariates
  set_groups_for_testing();   // set groups of snps to test jointly
  print_usage_info(&params, &files, sout);
  print_test_info();
  if(params.w_interaction && (params.trait_mode==0) && !params.no_robust && !params.force_robust) 
    nullHLM.prep_run(&pheno_data, &params);
  if(!params.skip_test) setup_output(&ofile, out, ofile_split, out_split); // result file
  if(params.trait_mode) set_nullreg_mat();
  sout << endl;

#ifdef WITH_HTSLIB
  if (params.remeta_save_ld) {
    for (size_t i = 0; i < files.pheno_names.size(); ++i) {
      if(params.pheno_pass(i)) {
        remeta_sumstats.skat_matrix_writers.emplace_back(
          RegenieLDMatrixWriter(
            files.out_file + "_" + files.pheno_names[i],
            params.pheno_counts(i, 0)
          )
        );
      } else {
        remeta_sumstats.skat_matrix_writers.emplace_back(
          RegenieLDMatrixWriter()
        );
      }
    }
    remeta_sumstats.sparsity_threshold = params.remeta_ld_spr;
  }
#endif

  // start analyzing each chromosome
  bool block_init_pass = false;
  int block = 0, chrom_nb, bs;
  tally snp_tally;
  vector< variant_block > block_info;
  if(params.joint_test) jt.scale_denum = params.n_analyzed - jt.ncovars; // for gates
  initialize_thread_data(Gblock.thread_data, params);

  for (auto const& chrom : files.chr_read){

    if( !in_map(chrom, chr_map) ) continue;

    chrom_nb = chr_map[chrom][1];

    // if no sets in chromosome, skip
    if(chrom_nb == 0)  continue;

    // If specified starting block
    if(!block_init_pass && (params.start_block > (block + chrom_nb)) ) {
      for(int bb = 0; bb < chrom_nb ; bb++)
        snp_tally.snp_count += jt.setinfo[chrom - 1][bb].snp_indices.size();
      block += chrom_nb;
      continue;
    }

    sout << "Chromosome " << chrom << " [" << chrom_nb << " sets in total]\n";

    // read polygenic effect predictions from step 1
    blup_read_chr(false, chrom, m_ests, files, in_filters, pheno_data, params, sout);

    // compute phenotype residual (adjusting for BLUP [and covariates for BTs])
    if(params.trait_mode == 1) compute_res_bin(chrom);
    else if(params.trait_mode == 2) compute_res_count(chrom);
    else compute_res();


    // analyze by blocks of SNPs
    for(int bb = 0; bb < chrom_nb ; bb++) {

      bs = jt.setinfo[chrom - 1][bb].snp_indices.size();

      // If specified starting block
      if(!block_init_pass && (params.start_block > (block+1)) ) {
        snp_tally.snp_count += bs;
        block++;
        continue;
      } else if(!block_init_pass) block_init_pass = true;

      sout << " set [" << block + 1 << "/" << params.total_n_block << "] : " << jt.setinfo[chrom - 1][bb].ID << " - " << bs << " variants..." << flush;
      if(params.joint_test && !params.build_mask) allocate_mat(Gblock.Gmat, params.n_samples, bs);
      block_info.resize(bs);

      // compute single snp association test statistic
      get_sum_stats(chrom, bb, block_info);

      // update number of variants (if masks were built)
      bs = block_info.size();
      jt.nvars = bs;

      if(params.skip_test) { // skip assoc tests
        snp_tally.snp_count += bs;
        block++;
        continue;
      }

      // tally the results
      for (auto const& snp_data : block_info){

        if( snp_data.ignored ) {
          snp_tally.n_ignored_snps++;
          jt.nvars--;
          continue;
        }

        snp_tally.n_ignored_tests += snp_data.ignored_trait.count();
        if(params.firth || params.use_SPA) {
          n_corrected += (!snp_data.ignored_trait && snp_data.is_corrected).count();
          snp_tally.n_failed_tests += (!snp_data.ignored_trait && snp_data.test_fail).count();
          if(params.w_interaction) {
            n_corrected += (2 + params.ncov_interaction) * snp_data.is_corrected_inter.count(); // main, inter & joint
            snp_tally.n_failed_tests += (2 + params.ncov_interaction) * (snp_data.is_corrected_inter && snp_data.test_fail_inter).count(); // main, inter & joint
          }
        }

        for(int j = 0; j < params.n_pheno; ++j) {

          if( !params.pheno_pass(j) || snp_data.ignored_trait(j) ) {
            if(!params.split_by_pheno) // if using single file, print NAs for snp/trait sum stats
              ofile << snp_data.sum_stats[j];

            continue;
          }

          if(params.split_by_pheno)
            (*ofile_split[j]) << snp_data.sum_stats[j]; // add test info
          else
            ofile << snp_data.sum_stats[j]; // add test info
        }

      }

      if( params.joint_test && block_info.size()){

        // compute and print set-based test result
        t1 = std::chrono::high_resolution_clock::now();
        sout << "     -computing joint association tests..." << flush;

        jt.get_variant_names(chrom, bb, snpinfo);
        tmp_str = jt.apply_joint_test(chrom, bb, &pheno_data, res, &Gblock, block_info, files.pheno_names, &params);

        for(int j = 0; j < params.n_pheno; ++j) {
          if(params.split_by_pheno)
            (*ofile_split[j]) << tmp_str[j];
          else
            ofile << tmp_str[j]; // add test info
        }

        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        sout << "done (" << duration.count() << "ms) "<< endl;
      }

      snp_tally.snp_count += bs;
      block++;
    }

  }

  sout << print_summary(&ofile, out, ofile_split, out_split, n_corrected, snp_tally, files, firth_est, params);

  if(params.write_masks)  bm.closeFiles();

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
    m1[ itr->first ] = itr->second;
  }
  chr_map = m1;


  // summarize block sizes
  sout << left << std::setw(20) << " * # threads" << ": [" << params.threads << "]\n";
  sout << left << std::setw(20) << " * # tested sets" << ": [" << params.total_n_block << "]\n";
  if(params.start_block > params.total_n_block)
    throw "Starting set > number of sets analyzed";
  else if(params.start_block > 1) sout << "    + skipping to set #" << params.start_block << endl;
  sout << left << std::setw(20) << " * max block size" << ": [" << params.block_size << "]\n";

  if(params.build_mask) sout << " * rule used to build masks : " << params.mask_rule << endl;
  if(params.vc_test) {

    if(params.skato_rho.size() == 0){
      params.skato_rho = ArrayXd::Zero(1); // assume rho=0
      sout << " * computing gene-based tests for each set of variants included in a mask\n";
    } else {
      if(CHECK_BIT(params.vc_test,1) && (params.skato_rho(0) != 0) ){
        ArrayXd tmp_rho (params.skato_rho.size()+1); tmp_rho(0) = 0; tmp_rho.tail(params.skato_rho.size()) = params.skato_rho; // insert rho=0 for skat
        params.skato_rho = tmp_rho;
      }
      IOFormat Fmt(StreamPrecision, DontAlignCols, ",", "", "", "","","");
      sout << " * computing gene-based tests for each set of variants included in a mask (rho=[" << params.skato_rho.matrix().transpose().format(Fmt) << "])\n";
    }

    sout << "  -variants with MAC <= " << params.skat_collapse_MAC << " are collapsed into a mask\n";
    if(params.vc_multiply_weights)
      sout << "  -user-provided weights will be multiplied by default weights [from Beta(MAF,"<< params.skat_a1 <<","<< params.skat_a2 <<")] for SKAT/ACAT tests\n";
    else if(params.vc_with_weights)
      sout << "  -user-provided weights will be used for gene-based tests\n";
    else
      sout << "  -weights are obtained from Beta(MAF,"<< params.skat_a1 <<","<< params.skat_a2 <<")\n";
    sout << "  -list of gene-based tests run: " << get_test_list(params.vc_test, params.vc_tests_map) << "\n";

    // set max rho to 0.999 as it will o.w. cause issue with skat-o p-value calculation
    if(params.skato_rho.size() > 1) params.skato_rho = params.skato_rho.min(0.999); 

    // single p per gene
    if(params.apply_gene_pval_strategy)
      sout << " * applying ACAT to output overall gene p-value\n";
  }

  if(params.remeta_save_ld)
    sout << " * saving SKAT LD matrices for REMETA\n";

}

// test SNPs in block
void Data::get_sum_stats(int const& chrom, int const& varset, vector<variant_block>& all_snps_info){

  vector< vector < uchar > > snp_data_blocks;
  vector< uint32_t > insize, outsize;

  if (params.mask_loo) {
    getMask_loo(chrom, varset, snp_data_blocks, insize, outsize, all_snps_info);
  } else {

    auto t1 = std::chrono::high_resolution_clock::now();
    vset* set_info = &(jt.setinfo[chrom - 1][varset]);

    // read in markers and if applicable build masks
    if(!params.build_mask) readChunk(set_info->snp_indices, chrom, snp_data_blocks, insize, outsize, all_snps_info);
    else {
      getMask(chrom, varset, snp_data_blocks, insize, outsize, all_snps_info);
      // update size with new masks
      int n_snps = set_info->snp_indices.size();
      //cerr << "M=" << n_snps << endl;
      if(params.skip_test || (n_snps == 0)) return;

      // starting association testing with built masks
      t1 = std::chrono::high_resolution_clock::now();
      sout << "     -computing association tests..." << flush;
    }

    // analyze using openmp
    compute_tests_mt(chrom, set_info->snp_indices, snp_data_blocks, insize, outsize, all_snps_info);

    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    sout << "done (" << duration.count() << "ms) "<< endl;

  }

}


void Data::readChunk(vector<uint64>& indices, int const& chrom, vector< vector < uchar > >& snp_data_blocks, vector<uint32_t>& insize, vector<uint32_t>& outsize, vector<variant_block>& all_snps_info){


  int const n_snps = indices.size();

  if((params.file_type == "bgen") && params.streamBGEN){
    snp_data_blocks.resize( n_snps );
    insize.resize(n_snps); outsize.resize(n_snps);
    vector<uint64> offsets(n_snps);
    for (int i = 0; i < n_snps; i++) offsets[i] = snpinfo[indices[i]].offset;

    readChunkFromBGEN(&files.geno_ifstream, insize, outsize, snp_data_blocks, offsets);

  } else if((params.file_type == "bgen") && !params.streamBGEN) 
    readChunkFromBGENFileToG(indices, chrom, snpinfo, &params, Gblock.Gmat, Gblock.bgen, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, all_snps_info, sout);
  else if(params.file_type == "pgen") 
    readChunkFromPGENFileToG(indices, chrom, &params, &in_filters, Gblock.Gmat, Gblock.pgr, pheno_data.masked_indivs, pheno_data.phenotypes_raw, snpinfo, all_snps_info);
  else {

    snp_data_blocks.resize( n_snps );
    for(int isnp = 0; isnp < n_snps; isnp++) {

      jumpto_bed( snpinfo[ indices[isnp] ].offset, files.bed_block_size, files.geno_ifstream);
      snp_data_blocks[isnp].resize(files.bed_block_size);
      files.geno_ifstream.read( reinterpret_cast<char *> (&snp_data_blocks[isnp][0]), files.bed_block_size);

    }
  }

}

void Data::getMask(int const& chrom, int const& varset, vector< vector < uchar > >& snp_data_blocks, vector<uint32_t>& insize, vector<uint32_t>& outsize, vector<variant_block>& all_snps_info){

  auto t1 = std::chrono::high_resolution_clock::now();
  vset* set_info = &(jt.setinfo[chrom - 1][varset]);

  // do it in chunks to reduce memory usage
  bool last_chunk = false;
  int n_snps = set_info->snp_indices.size(), nvar_read = 0;
  int nchunks, bsize; 
  SpMat vc_sparse_gmat;
  if(params.use_max_bsize) { // process all variants at once
    nchunks = 1;
    bsize = n_snps; 
  } else { // process variants in blocks
    nchunks = ceil( n_snps * 1.0 / params.block_size );
    bsize = params.block_size; // default number of SNPs to read at a time
  }
  //if(params.mask_loo) bm.nmasks_total = n_snps;

  // custom user weights
  ArrayXd snp_weights = ArrayXd::Constant(n_snps, 1), vc_weights, vc_weights_acat;
  if( params.vc_with_weights )
    if(!get_custom_weights(set_info->ID, snp_weights, snpinfo, set_info->snp_indices)){
      sout << "\n     -WARNING: all variants have 0 weights (set will be skipped)\n";
      set_info->snp_indices.resize(0);
      all_snps_info.resize(0);
      return;
    }

  if(params.verbose) sout << nchunks << " chunks";
  sout << "\n     -reading in genotypes" << ( params.vc_test ? ", computing gene-based tests" : "" ) << " and building masks..." << flush;

  if(params.debug) sout << "(1)" << print_mem() << "..." << flush;
  bm.prepMasks(params.n_samples, set_info->ID);
  allocate_mat(Gblock.Gmat, params.n_samples, bsize);
  if(params.vc_test) {
    set_info->Jmat = MatrixXb::Constant(n_snps + bm.nmasks_total, bm.nmasks_total, false); // MxKm (last S rows are for ultra-rare masks)
    set_info->ultra_rare_ind = ArrayXb::Constant(n_snps, false); // identify which vars are rare
    set_info->vc_rare_mask.resize(params.n_samples, bm.nmasks_total); // identify which vars are rare
    set_info->vc_rare_mask.setZero();
    set_info->vc_rare_mask_non_missing = MatrixXb::Constant(params.n_samples, bm.nmasks_total, false); // distinguish 0 from missing
    vc_sparse_gmat.resize(params.n_samples, n_snps + bm.nmasks_total); // store wG
    vc_sparse_gmat.setZero();
    vc_weights = ArrayXd::Zero(n_snps + bm.nmasks_total, 1);
    vc_weights.head( n_snps ) = snp_weights;
    vc_weights_acat = vc_weights;
  }
  if(params.debug) sout << "(2)" << print_mem() << "..." << flush;


  for(int i = 0; i < nchunks; i++){

    last_chunk = ( i == (nchunks-1) );
    if( last_chunk ) {
      bsize = n_snps - i * bsize;// use remainder number of variants
      allocate_mat(Gblock.Gmat, params.n_samples, bsize);
    }

    vector<uint64> indices (set_info->snp_indices.begin() + nvar_read, set_info->snp_indices.begin() + nvar_read + bsize);
    readChunk(indices, chrom, snp_data_blocks, insize, outsize, all_snps_info);

    // build genotype matrix
    if( ((params.file_type == "bgen") && params.streamBGEN) || params.file_type == "bed") {
#if defined(_OPENMP)
      setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
      for(int isnp = 0; isnp < bsize; isnp++) {

        uint32_t snp_index = indices[isnp];

        variant_block* block_info = &(all_snps_info[isnp]);
        parseSNP(isnp, chrom, &(snp_data_blocks[isnp]), insize[isnp], outsize[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &snpinfo[snp_index], &Gblock, block_info, sout);

      }
#if defined(_OPENMP)
      setNbThreads(params.threads);
#endif
    }

    // update mask (taking max/sum)
    /*if(params.mask_loo)
      bm.updateMasks_loo(nvar_read, bsize, &params, &in_filters, pheno_data.masked_indivs, &Gblock, all_snps_info, *set_info, snpinfo, sout);
    else*/
      bm.updateMasks(nvar_read, bsize, &params, &in_filters, pheno_data.masked_indivs, &Gblock, snp_weights, all_snps_info, *set_info, snpinfo, sout);

    if(params.vc_test) // get G and w
      update_vc_gmat(vc_sparse_gmat, vc_weights, vc_weights_acat, set_info->ultra_rare_ind, nvar_read, bsize, params, in_filters.ind_in_analysis, Gblock.Gmat, all_snps_info, set_info->Jmat);

    /*
    if(params.debug){ 
      cerr << "GG.diag()=\n" << Gblock.Gmat.array().square().colwise().sum() << "\n";
      MatrixXd Gtmp_res;
      residualize_gmat(false, pheno_data.new_cov, Gblock.Gmat, Gtmp_res, params);
      cerr << "GrGr.diag()=\n" << Gtmp_res.array().square().colwise().sum() << "\n";
      cerr << "WGGW.diag()=\n" << vc_weights.head(bsize).square() * Gtmp_res.array().square().colwise().sum().matrix().transpose().array() << "\n";
    }
    */

    nvar_read += bsize;
  }

  // check mask and store in setinfo & snpinfo
  /*if(params.mask_loo)
    bm.computeMasks_loo(&params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &Gblock, all_snps_info, *set_info, snpinfo, sout);
  else*/
    bm.computeMasks(&params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &Gblock, all_snps_info, *set_info, snpinfo, sout);
  if(params.debug) sout << "(3)" << print_mem() << "..." << flush;

  if(params.vc_test) {
    #ifdef WITH_HTSLIB
      remeta_sumstats.skat_snplist = &bm.remeta_snplist;
      remeta_sumstats.gene_name = &set_info->ID;
    #endif
    compute_vc_masks(vc_sparse_gmat, vc_weights, vc_weights_acat, set_info->vc_rare_mask, set_info->vc_rare_mask_non_missing, pheno_data.new_cov, m_ests, firth_est, res, pheno_data.phenotypes_raw, pheno_data.masked_indivs, set_info->Jmat, all_snps_info, in_filters.ind_in_analysis, params, remeta_sumstats); 

    set_info->Jmat.resize(0,0);
    set_info->ultra_rare_ind.resize(0);
    set_info->vc_rare_mask.setZero(); set_info->vc_rare_mask.resize(0,0); set_info->vc_rare_mask.data().squeeze();
    set_info->vc_rare_mask_non_missing.resize(0,0);
  }
  if(params.debug) sout << "(4)" << print_mem() << "..." << flush;

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;

}

void Data::getMask_loo(int const& chrom, int const& varset, vector< vector < uchar > >& snp_data_blocks, vector<uint32_t>& insize, vector<uint32_t>& outsize, vector<variant_block>& all_snps_info){

  MeasureTime mt;
  vset* set_info = &(jt.setinfo[chrom - 1][varset]);
  vector<uint64> orig_indices = set_info->snp_indices;
  vector<variant_block> out_snps_info;
  int n_snps = orig_indices.size();

  // read in variants in chunks storing as sparse matrix
  sout << "\n     -reading in genotypes..." << flush;
  mt.start_ms();
  if(params.debug) sout << "(0)" << print_mem() << "...";
  uint64 snp_index; 
  SpMat rv_mat(params.n_samples, n_snps);
  ArrayXb in_lovo_mask(n_snps), ur_variant(n_snps), var_flip(n_snps);
  in_lovo_mask = false; ur_variant = false; var_flip = false;
  ArrayXd Gvec(params.n_samples), var_mafs(n_snps);

  for(int snp = 0, j = 0; snp < n_snps; snp++){
    snp_index = orig_indices[ snp ];
    read_snp(false, snpinfo[ snp_index ].offset, Gvec, in_filters.ind_in_analysis, in_filters.ind_ignore, &files, Gblock.pgr, &params, false);
    in_lovo_mask(snp) = bm.check_in_lovo_mask(Gvec, in_filters, set_info->ID, snpinfo[ snp_index ], ur_variant(j), var_flip(j), var_mafs(j), chrom, &params);
    if( in_lovo_mask(snp) )
      rv_mat.col(j++) = Gvec.matrix().sparseView();
  }

  int n_snps_lovo = in_lovo_mask.count();
  ArrayXi col_indices_lovo_mask = get_true_indices( in_lovo_mask );
  // remove cols/entries not in any lovo masks
  rv_mat.conservativeResize(params.n_samples, n_snps_lovo);
  rv_mat.makeCompressed();
  ur_variant.conservativeResize(n_snps_lovo, 1);
  var_flip.conservativeResize(n_snps_lovo, 1);
  var_mafs.conservativeResize(n_snps_lovo, 1);
  ArrayXi col_indices_lovo_mask_non_ur = get_true_indices( !ur_variant );

  // custom user weights
  ArrayXd snp_weights = ArrayXd::Constant(n_snps_lovo, 1);
  if( params.vc_with_weights )
    if(!get_custom_weights(set_info->ID, snp_weights, snpinfo, col_indices_lovo_mask, orig_indices)){
      sout << "\n     -WARNING: all variants have 0 weights (set will be skipped)\n";
      return;
    }

  if(params.debug) sout << n_snps_lovo << " variants...(1)" << print_mem() << "...";
  sout << mt.stop_ms() << "\n";

  // generate LOVO masks in chunks
  bool last_chunk = false;
  ArrayXi lovo_masks_indices = check_lovo_snplist(col_indices_lovo_mask, orig_indices, snpinfo, params.masks_loo_snpfile); // if computing a subset of the lovo masks
  int neff_lovo = lovo_masks_indices.size(), nchunks, bsize, nvar_read = 0; // default number of SNPs to read at a time
  bsize = min(neff_lovo, 128);
  nchunks = ceil( neff_lovo * 1.0 / bsize );

  sout << "     -splitting into " <<  nchunks << " chunk" << (nchunks > 1 ? "s" : "") << " of size " << bsize << " (" << neff_lovo << " LOVO masks in total)\n";
  mt.start_ms();

  // For SKAT tests
  ArrayXd vc_weights, vc_weights_acat;
  SpMat vc_sparse_gmat; // contains single variants + ultra-rare masks (incl.for full mask)
  if(params.vc_test) {
    MeasureTime mt_skat;
    if(params.debug) mt_skat.start_ms();
    vc_sparse_gmat.resize(params.n_samples, n_snps_lovo); // store wG
    vc_sparse_gmat.setZero();
    vc_weights = snp_weights;
    vc_weights_acat = vc_weights;
    // store G & SKAT weights
    update_vc_gmat(vc_sparse_gmat, vc_weights, vc_weights_acat, rv_mat, ur_variant, var_flip, var_mafs, in_filters.ind_in_analysis, params);
    if(params.trait_mode == 1) // if need to apply cc correction
      check_cc_correction(vc_sparse_gmat, vc_weights, pheno_data.new_cov, m_ests, firth_est, res, pheno_data.phenotypes_raw, pheno_data.masked_indivs, params); 
    if(params.debug) sout << "(0)skat prep..."<< mt_skat.stop_ms() << "\n";
  }

  for(int i = 0; i < nchunks; i++) {

    last_chunk = ( i == (nchunks-1) );
    if( last_chunk )
      bsize = neff_lovo - i * bsize; // use remainder number of variants
    set_info->snp_indices = orig_indices;
    MeasureTime mt_chunk;

    sout << "      +chunk #" << i+1 << " (" << bsize << " LOVO masks)\n       -" << ( params.vc_test ? "computing gene-based tests and " : "" ) << "building masks..." << flush;
    mt_chunk.start_ms();

    bm.nmasks_total = bsize + last_chunk;
    bm.prepMasks(params.n_samples, set_info->ID);  
    if(!bm.take_max && !bm.take_comphet) {
      bm.nsites = ArrayXi::Constant(bm.nmasks_total, n_snps_lovo - 1); // each loo mask has (n-1) sites included for AAF calculation
      if( last_chunk ) bm.nsites.tail(1) += 1; // last entry is for full mask
    }

    ArrayXi chunk_indices = lovo_masks_indices.segment(nvar_read, bsize);
    ArrayXb in_chunk = ArrayXb::Constant(n_snps_lovo, false);
    in_chunk(chunk_indices) = true;
    if(params.debug) sout << "(1)" << print_mem() << "..." << flush;

    // collapse mask for variants not in chunk
    ArrayXd mask_excl_chunk = ArrayXd::Constant(params.n_samples, -3);
    ArrayXd ur_mask_excl_chunk = ArrayXd::Constant(params.n_samples, -3);
    if((!in_chunk).any())
      bm.collapse_mask_chunk(get_true_indices(!in_chunk), rv_mat, ur_variant, var_flip, snp_weights, mask_excl_chunk, ur_mask_excl_chunk, in_filters.ind_in_analysis);
    if(params.debug) {
      sout << mt_chunk.stop_ms() << "\n";
      mt_chunk.start_ms();
    }

    // update mask (taking max/sum)
    bm.updateMasks_loo(chunk_indices, last_chunk, rv_mat, ur_variant, var_flip, snp_weights, mask_excl_chunk, ur_mask_excl_chunk, in_filters.ind_in_analysis, *set_info, params.threads);
    if(params.debug) {
      sout << "(2)" << print_mem() << "..." << flush;
      sout << mt_chunk.stop_ms() << "\n";
      mt_chunk.start_ms();
    }

    // check mask and store in setinfo & snpinfo
    bm.computeMasks_loo(col_indices_lovo_mask(chunk_indices), last_chunk, &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &Gblock, all_snps_info, *set_info, snpinfo, sout);
    if(params.debug) {
      sout << "(3)" << print_mem() << "..." << flush;
      sout << mt_chunk.stop_ms() << "\n";
      mt_chunk.start_ms();
    }

    if(params.vc_test) {// run skat/acat

      MatrixXb Jmat = MatrixXb::Constant(n_snps_lovo + bm.nmasks_total, bm.nmasks_total, false); // rows: snps + ur masks; cols: lovo sets
      Jmat(col_indices_lovo_mask_non_ur, all).array() = true; // non-ur snps
      for(int j = 0; j < chunk_indices.size(); j++) // apply lovo
        Jmat( chunk_indices(j), j) = false;

      SpMat vc_sparse_gmat_chunk(params.n_samples, Jmat.rows()); // cols: snps + ur masks
      vc_sparse_gmat_chunk.reserve(vc_sparse_gmat.nonZeros());
      vc_sparse_gmat_chunk.leftCols(vc_sparse_gmat.cols()) = vc_sparse_gmat;

      ArrayXd vc_weights_chunk(Jmat.rows()), vc_weights_acat_chunk(Jmat.rows());
      vc_weights_chunk.head(vc_weights.size()) = vc_weights;
      vc_weights_acat_chunk.head(vc_weights_acat.size()) = vc_weights_acat;

      compute_vc_masks(vc_sparse_gmat_chunk, vc_weights_chunk, vc_weights_acat_chunk, set_info->vc_rare_mask, set_info->vc_rare_mask_non_missing, pheno_data.new_cov, m_ests, firth_est, res, pheno_data.phenotypes_raw, pheno_data.masked_indivs, Jmat, all_snps_info, in_filters.ind_in_analysis, params, remeta_sumstats); 

    }

    if(params.debug) sout << "(4)" << print_mem() << "..." << flush;
    nvar_read += bsize;
    sout << mt_chunk.stop_ms() << "\n";

    if(params.skip_test || (set_info->snp_indices.size() == 0)) continue;

    // burden association tests with built masks
    sout << "       -computing association tests..." << flush;
    mt_chunk.start_ms();
    compute_tests_mt(chrom, set_info->snp_indices, snp_data_blocks, insize, outsize, all_snps_info);
    sout << mt_chunk.stop_ms() << "\n";
    out_snps_info.insert(out_snps_info.end(), all_snps_info.begin(), all_snps_info.end());
  }

  if(params.vc_test) {
    set_info->vc_rare_mask.setZero(); set_info->vc_rare_mask.resize(0,0); set_info->vc_rare_mask.data().squeeze();
    set_info->vc_rare_mask_non_missing.resize(0,0);
  }

  sout << "     -> " << mt.stop_ms() << "\n";
  all_snps_info = out_snps_info;

}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Testing mode (multi-trait tests)
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::test_multitrait() 
{
  sout << "Association testing mode (multi-trait tests)";

  std::chrono::high_resolution_clock::time_point t1, t2;
  string out;
  vector<string> out_split, tmp_str;
  // output files
  Files ofile;
  // use pointer to class since it contains non-copyable elements
  vector <std::shared_ptr<Files>> ofile_split;

#if defined(_OPENMP)
  sout << " with " << (params.streamBGEN? "fast " : "") << "multithreading using OpenMP";
#endif
  sout << endl;

  // Set up 
  file_read_initialization(); // set up files for reading
  read_pheno_and_cov(&files, &params, &in_filters, &pheno_data, &m_ests, &Gblock, sout);   // read phenotype and covariate files
  prep_run(&files, &in_filters, &params, &pheno_data, &m_ests, sout); // check blup files and adjust for covariates
  set_blocks_for_testing();   // set number of blocks
  print_usage_info(&params, &files, sout);
  print_test_info();
  setup_output(&ofile, out, ofile_split, out_split); // result files
  sout << endl;

  // Set up mt for all chr: verbose level, masks 
  /* if(params.mt_out_all) mt.verbose = 3; */
  /* mt.verbose = 3; */
  /* if(params.mt_precomp) mt.precomp = true; */
  mt.precomp = true;

  mt.setup_masks(pheno_data.masked_indivs);

  // Loop 1: start analyzing each chromosome
  bool block_init_pass = false;
  int block = 0, chrom_nsnps, chrom_nb, bs;
  tally snp_tally;
  vector< variant_block > block_info;
  /* initialize_thread_data(Gblock.thread_data, params); */

  for(auto const& chrom : files.chr_read) {
    if( !in_map(chrom, chr_map) ) continue;

    chrom_nsnps = chr_map[chrom][0];
    chrom_nb = chr_map[chrom][1];
    if(chrom_nb == 0) continue;

    // If specified starting block
    if(!block_init_pass && (params.start_block > (block + chrom_nb)) ) {
      snp_tally.snp_count += chrom_nsnps;
      block += chrom_nb;
      continue;
    }

    sout << "Chromosome " << chrom << " [" << chrom_nb << " blocks in total]\n";

    // read polygenic effect predictions from step 1
    blup_read_chr(false, chrom, m_ests, files, in_filters, pheno_data, params, sout);

    // compute phenotype residual (adjusting for BLUP)
    if(params.trait_mode == 0) {
      compute_res();
    } else {
      throw std::runtime_error("multi-trait tests only for QTs");
    }

    // Set up mt for each chr: matrix of traits Y
    mt.setup_yres(res);

    /* const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n"); */
    /* string f_cory = files.out_file + ".regenie.Ycor.chr" + to_string(chrom) + ".txt"; */
    /* mt.compute_cory(mt.Yres, mt.Mask); */
    /* ofstream out_cory(f_cory.c_str()); */
    /* out_cory << mt.Ryy.format(CSVFormat); */
    /* out_cory.close(); */

    // analyze by blocks of SNPs
    for(int bb = 0; bb < chrom_nb ; bb++) {
      get_block_size(params.block_size, chrom_nsnps, bb, bs);

      // If specified starting block
      if(!block_init_pass && (params.start_block > (block+1)) ) {
        snp_tally.snp_count += bs;
        block++;
        continue;
      } else {
        if(!block_init_pass) block_init_pass = true;
      }

      sout << " block [" << block + 1 << "/" << params.total_n_block << "] : " << flush;

      allocate_mat(Gblock.Gmat, params.n_samples, bs);
      block_info.resize(bs);

      // read SNP, impute missing & compute association test statistic
      analyze_block_multitrait(chrom, bs, &snp_tally, block_info);

      // print the results
      if(params.split_by_pheno) {
        throw std::runtime_error("test_multitrait: split_by_pheno");
      }

      for (auto const& snp_data : block_info){
        if( snp_data.ignored ) {
          snp_tally.n_ignored_snps++;
          continue;
        }
        snp_tally.n_ignored_tests += snp_data.ignored_trait.count();

        /* size_t n_trait_sets = 1; */
        size_t j = 0;
        ofile << snp_data.sum_stats_mt[j]; // add test info
      }

      snp_tally.snp_count += bs;
      block++;
    }

  }

  sout << print_summary(&ofile, out, ofile_split, out_split, n_corrected, snp_tally, files, firth_est, params);
}

// test SNPs in block for multi-trait tests
void Data::analyze_block_multitrait(int const& chrom, int const& n_snps, tally* snp_tally, vector<variant_block> &all_snps_info){

  auto t1 = std::chrono::high_resolution_clock::now();
  const int start = snp_tally->snp_count;
  vector< vector < uchar > > snp_data_blocks;
  vector< uint32_t > insize, outsize;

  vector<uint64> indices(n_snps);
  std::iota(indices.begin(), indices.end(), start);

  readChunk(indices, chrom, snp_data_blocks, insize, outsize, all_snps_info);

  // analyze using openmp
  /* compute_tests_mt(chrom, indices, snp_data_blocks, insize, outsize, all_snps_info); */
  compute_tests_mt_multitrait(chrom, indices, snp_data_blocks, insize, outsize, all_snps_info);

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << "done (" << duration.count() << "ms) "<< endl;
}

void Data::compute_tests_mt_multitrait(int const& chrom, vector<uint64> indices,vector< vector < uchar > >& snp_data_blocks, vector< uint32_t > insize, vector< uint32_t >& outsize, vector<variant_block> &all_snps_info){
  size_t const bs = indices.size();
  ArrayXb err_caught = ArrayXb::Constant(bs, false);

  // start openmp for loop
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(size_t isnp = 0; isnp < bs; isnp++) {
    uint32_t const snp_index = indices[isnp];

    int thread_num = 0;
#if defined(_OPENMP)
    thread_num = omp_get_thread_num();
#endif

    // to store variant information
    variant_block* block_info = &(all_snps_info[isnp]);
    reset_thread(&(Gblock.thread_data[thread_num]), params);

    parseSNP(isnp, chrom, &(snp_data_blocks[isnp]), insize[isnp], outsize[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &snpinfo[snp_index], &Gblock, block_info, sout);

    // for QTs: project out covariates & scale
    residualize_geno(isnp, thread_num, block_info, false, pheno_data.new_cov, &Gblock, &params);

    // skip SNP if fails filters
    if( block_info->ignored ) continue;
    
    reset_stats(block_info, params);

    try {
      // run multi-trait tests & save summary stats
      // v1: store results in mt 
      // - doesn't work for multi-threaded calculations
      /* mt.apply_tests_snp(isnp, Gblock, res, p_sd_yres, params); */
      /* string tmp_str = mt.print_sumstats(isnp, snp_index, test_string + params.condtl_suff, model_type + params.condtl_suff, block_info, snpinfo, &params); */
      // v2: store results outside mt, i.e. separately for each thread
      /* MTestsResults mt_results_i = mt.run_tests_snp(isnp, Gblock, res, p_sd_yres, params); */
      // v3: pre-load Yres/Y0res
      MTestsResults mt_results_i = mt.run_tests_snp_precomp(isnp, Gblock, params);
      /* string tmp_str = mt.print_sumstats(isnp, snp_index, test_string + params.condtl_suff, model_type + params.condtl_suff, block_info, snpinfo, &params); */
      string tmp_str = mt.print_sumstats(mt_results_i, isnp, snp_index, test_string + params.condtl_suff, model_type + params.condtl_suff, block_info, snpinfo, &params);

      block_info->sum_stats_mt[0].append(tmp_str);
    } catch (...) {
      err_caught(isnp) = true;
      block_info->sum_stats[0] = boost::current_exception_diagnostic_information();
      continue;
    }
  }

#if defined(_OPENMP)
  setNbThreads(params.threads);
#endif

  // check no errors
  if(err_caught.any())
    for(int i = 0; i < err_caught.size(); i++)
      if(err_caught(i)) throw all_snps_info[i].sum_stats[0];

}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    Testing mode (MultiPhen test)
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::test_multiphen() 
{
  sout << "Association testing mode (MultiPhen test)";

  std::chrono::high_resolution_clock::time_point t1, t2;
  string out;
  vector<string> out_split, tmp_str;
  // output files
  Files ofile;
  // use pointer to class since it contains non-copyable elements
  vector <std::shared_ptr<Files>> ofile_split;

#if defined(_OPENMP)
  sout << " with " << (params.streamBGEN? "fast " : "") << "multithreading using OpenMP";
#endif
  sout << endl;

  // Set up 
  file_read_initialization(); // set up files for reading
  read_pheno_and_cov(&files, &params, &in_filters, &pheno_data, &m_ests, &Gblock, sout);   // read phenotype and covariate files
  prep_run(&files, &in_filters, &params, &pheno_data, &m_ests, sout); // check blup files and adjust for covariates
  set_blocks_for_testing();   // set number of blocks
  print_usage_info(&params, &files, sout);
  print_test_info();
  setup_output(&ofile, out, ofile_split, out_split); // result files
  sout << endl;

  // Set up mt
  prep_multiphen();

  // Loop 1: start analyzing each chromosome
  bool block_init_pass = false;
  int block = 0, chrom_nsnps, chrom_nb, bs;
  tally snp_tally;
  vector< variant_block > block_info;
  /* initialize_thread_data(Gblock.thread_data, params); */

  for(auto const& chrom : files.chr_read) {
    if( !in_map(chrom, chr_map) ) continue;

    chrom_nsnps = chr_map[chrom][0];
    chrom_nb = chr_map[chrom][1];
    if(chrom_nb == 0) continue;

    // If specified starting block
    if(!block_init_pass && (params.start_block > (block + chrom_nb)) ) {
      snp_tally.snp_count += chrom_nsnps;
      block += chrom_nb;
      continue;
    }

    sout << "Chromosome " << chrom << " [" << chrom_nb << " blocks in total]\n";

    // read polygenic effect predictions from step 1
    blup_read_chr(false, chrom, m_ests, files, in_filters, pheno_data, params, sout);

    // compute phenotype residual (adjusting for BLUP)
    if(params.trait_mode == 0) {
      compute_res();
      set_multiphen();
    } else {
      throw std::runtime_error("MultiPhen test for QTs only");
    }

    // analyze by blocks of SNPs
    for(int bb = 0; bb < chrom_nb ; bb++) {
      get_block_size(params.block_size, chrom_nsnps, bb, bs);

      // If specified starting block
      if(!block_init_pass && (params.start_block > (block+1)) ) {
        snp_tally.snp_count += bs;
        block++;
        continue;
      } else {
        if(!block_init_pass) block_init_pass = true;
      }

      sout << " block [" << block + 1 << "/" << params.total_n_block << "] : " << flush;

      allocate_mat(Gblock.Gmat, params.n_samples, bs);
      block_info.resize(bs);

      // read SNP, impute missing & compute association test statistic
      analyze_block_multiphen(chrom, bs, &snp_tally, block_info);

      // print the results
      if(params.split_by_pheno) {
        // ignore split_by_pheno for MultiPhen
        // throw std::runtime_error("test_multiphen: split_by_pheno");
      }

      for (auto const& snp_data : block_info){
        if( snp_data.ignored ) {
          snp_tally.n_ignored_snps++;
          continue;
        }
        snp_tally.n_ignored_tests += snp_data.ignored_trait.count();

        /* size_t n_trait_sets = 1; */
        size_t j = 0;
        ofile << snp_data.sum_stats_multiphen[j]; // add test info
      }

      snp_tally.snp_count += bs;
      block++;
    }

  }

  sout << print_summary(&ofile, out, ofile_split, out_split, n_corrected, snp_tally, files, firth_est, params);
}

// test SNPs in block for multi-trait tests
void Data::analyze_block_multiphen(int const& chrom, int const& n_snps, tally* snp_tally, vector<variant_block> &all_snps_info){

  auto t1 = std::chrono::high_resolution_clock::now();
  const int start = snp_tally->snp_count;
  vector< vector < uchar > > snp_data_blocks;
  vector< uint32_t > insize, outsize;

  vector<uint64> indices(n_snps);
  std::iota(indices.begin(), indices.end(), start);

  readChunk(indices, chrom, snp_data_blocks, insize, outsize, all_snps_info);

  // analyze using openmp
  compute_tests_mt_multiphen(chrom, indices, snp_data_blocks, insize, outsize, all_snps_info);

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << "done (" << duration.count() << "ms) "<< endl;
}

void Data::compute_tests_mt_multiphen(int const& chrom, vector<uint64> indices,vector< vector < uchar > >& snp_data_blocks, vector< uint32_t > insize, vector< uint32_t >& outsize, vector<variant_block> &all_snps_info){
  size_t const bs = indices.size();
  ArrayXb err_caught = ArrayXb::Constant(bs, false);

  // start openmp for loop
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(size_t isnp = 0; isnp < bs; isnp++) {
    uint32_t const snp_index = indices[isnp];

    int thread_num = 0;
#if defined(_OPENMP)
    thread_num = omp_get_thread_num();
#endif

    // to store variant information
    variant_block* block_info = &(all_snps_info[isnp]);
    reset_thread(&(Gblock.thread_data[thread_num]), params);

    parseSNP(isnp, chrom, &(snp_data_blocks[isnp]), insize[isnp], outsize[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &snpinfo[snp_index], &Gblock, block_info, sout);

    // for QTs: project out covariates & scale
    /* residualize_geno(isnp, thread_num, block_info, false, pheno_data.new_cov, &Gblock, &params); */

    // skip SNP if fails filters
    if( block_info->ignored ) continue;
    
    reset_stats(block_info, params);

    try {
      // load one SNP into Gmat
      MapMatXd Gmat(Gblock.Gmat.col(isnp).data(), params.n_samples, 1);
      // create a copy of mphen for every SNP (for parallel processing)
      MultiPhen mphen_i = mphen;

      // run MultiPhen test & save summary stats
      /* cout << snpinfo[snp_index].ID << endl; */
      mphen_i.run(Gmat, pheno_data.cov_phenotypes, pheno_data.new_cov.cols() - 1, params.n_pheno); // the last 2 arg.: #cov excluding intercept; #phenotypes

      /* string tmp_str = mphen_i.print_sumstats(isnp, snp_index, test_string + params.condtl_suff, model_type + params.condtl_suff, block_info, snpinfo, &params); */
      /* print_sum_stats_line(int const& snp_index, int const& i, string const& tmpstr, string const& test_string, string const& model_type, variant_block* block_info, data_thread* dt_thr, vector<snp> const& snpinfo, struct in_files const& files, struct param const& params){ */
      std::ostringstream buffer;
      if(params.htp_out) buffer << print_sum_stats_head_htp(snp_index, "MultiPhen", model_type + params.condtl_suff, snpinfo, &params) << mphen_i.print_sum_stats_htp(block_info, &params);
      else buffer << mphen_i.print_sumstats(isnp, snp_index, test_string + params.condtl_suff, model_type + params.condtl_suff, block_info, snpinfo, &params); 
      /* else buffer << print_sum_stats( */
      /*     (params.split_by_pheno ? block_info->af(i) : block_info->af1), */ 
      /*     block_info->af_case(i),block_info->af_control(i), */ 
      /*     (params.split_by_pheno ? block_info->info(i) : block_info->info1), */
      /*     (params.split_by_pheno ? block_info->ns(i) : block_info->ns1), */ 
      /*     block_info->ns_case(i), */ 
      /*     block_info->ns_control(i), */ 
      /*     test_string, */ 
      /*     dt_thr->bhat(i), dt_thr->se_b(i), dt_thr->chisq_val(i), dt_thr->pval_log(i), !block_info->test_fail(i), 1, &params, (i+1)); */ 

    /* else */  
    /*   buffer << (!params.split_by_pheno && (i>0) ? "" : tmpstr) << print_sum_stats((params.split_by_pheno ? block_info->af(i) : block_info->af1), block_info->af_case(i),block_info->af_control(i), (params.split_by_pheno ? block_info->info(i) : block_info->info1), (params.split_by_pheno ? block_info->ns(i) : block_info->ns1), block_info->ns_case(i), block_info->ns_control(i), test_string, dt_thr->bhat(i), dt_thr->se_b(i), dt_thr->chisq_val(i), dt_thr->pval_log(i), !block_info->test_fail(i), 1, &params, (i+1)); */
    /* return buffer.str(); */
    /* if(params.htp_out) */ 
    /*   buffer <<  print_sum_stats_head_htp(snp_index, files.pheno_names[i], model_type, snpinfo, &params) << 
     *   print_sum_stats_htp(dt_thr->bhat(i), dt_thr->se_b(i), dt_thr->chisq_val(i), dt_thr->pval_log(i), 
     *   block_info->af(i), block_info->info(i), block_info->mac(i), block_info->genocounts, 
     *   i, !block_info->test_fail(i), 1, &params, dt_thr->scores(i), dt_thr->cal_factor(i)); */
    /* else */  
    /*   buffer << (!params.split_by_pheno && (i>0) ? "" : tmpstr) << print_sum_stats((params.split_by_pheno ? block_info->af(i) : block_info->af1), block_info->af_case(i),block_info->af_control(i), (params.split_by_pheno ? block_info->info(i) : block_info->info1), (params.split_by_pheno ? block_info->ns(i) : block_info->ns1), block_info->ns_case(i), block_info->ns_control(i), test_string, dt_thr->bhat(i), dt_thr->se_b(i), dt_thr->chisq_val(i), dt_thr->pval_log(i), !block_info->test_fail(i), 1, &params, (i+1)); */
      std::string tmp_str = buffer.str();

      // 1 set of traits
      block_info->sum_stats_multiphen[0].append(tmp_str);
    } catch (...) {
      err_caught(isnp) = true;
      block_info->sum_stats[0] = boost::current_exception_diagnostic_information();
      continue;
    }
  }

#if defined(_OPENMP)
  setNbThreads(params.threads);
#endif

  // check no errors
  if(err_caught.any())
    for(int i = 0; i < err_caught.size(); i++)
      if(err_caught(i)) throw all_snps_info[i].sum_stats[0];

}

void Data::prep_multiphen()
{
  // user parameters 
  mphen.test = params.multiphen_test;
  mphen.optim = params.multiphen_optim;
  mphen.pval_thr = params.multiphen_thr;
  mphen.firth_mult = params.multiphen_firth_mult;
  mphen.firth_binom = (params.multiphen_firth_mult > 0);
  mphen.firth_multinom = (params.multiphen_firth_mult > 0);
  mphen.tol = params.multiphen_tol;
  mphen.trace = params.multiphen_trace;
  mphen.verbose = params.multiphen_verbose;
  // parameters for model fitting
  mphen.maxit = params.multiphen_maxit; mphen.maxit2 = params.multiphen_maxit2;
  mphen.strict = params.multiphen_strict;
  mphen.check_step = (params.multiphen_maxstep > 0);
  mphen.max_step = params.multiphen_maxstep; 
  mphen.reuse_start = true;
  mphen.mac_approx_offset = params.multiphen_approx_offset;
  mphen.pseudo_stophalf = params.multiphen_pseudo_stophalf;
  mphen.reset_start = params.multiphen_reset_start;
  mphen.offset_mode = params.multiphen_offset;

  // prepare new matrix of covariates X + matrix of phenotypes Y
  unsigned int n_samples = pheno_data.new_cov.rows(), n_cov1 = pheno_data.new_cov.cols();
  unsigned int n_cov = n_cov1 - 1;
  unsigned int n_phen = params.n_pheno;

  pheno_data.cov_phenotypes.resize(n_samples, n_cov + 2*n_phen + 2); // +2 intercepts
  // column # 1 = Intercept
  pheno_data.cov_phenotypes.col(0) = ArrayXd::Constant(n_samples, 1.0);
  // next n_cov columns = covariates **without** intercept
  // new_cov has intercept in the last column
  if(n_cov) {
    pheno_data.cov_phenotypes.leftCols(n_cov1).rightCols(n_cov) = pheno_data.new_cov.leftCols(n_cov);
  }
  // next n_phen columns = phenotypes (skipped here & to be filled in for each chr.)
  // next & the last column = Intercept
  pheno_data.cov_phenotypes.rightCols(1) = ArrayXd::Constant(n_samples, 1.0);

  // v2
  if(!params.strict_mode) throw std::runtime_error("--strict mode is required for MultiPhen test");

  VectorXb Mask = pheno_data.masked_indivs.col(0);
  for(unsigned int i = 1; i < pheno_data.masked_indivs.cols(); i++) {
    Mask.col(0).array() = Mask.col(0).array() || pheno_data.masked_indivs.col(i).array();
  }
  pheno_data.cov_phenotypes.array().colwise() *= Mask.array().cast<double>().array();
  mphen.setup_x(Mask, pheno_data.cov_phenotypes, n_cov, n_phen, true, false); // (ignored by MultiPhen) pos_intercept_first = true, pos_phen_first = false
}

void Data::set_multiphen()
{
  unsigned int n_samples = pheno_data.new_cov.rows(), n_cov1 = pheno_data.new_cov.cols();
  unsigned int n_cov = n_cov1 - 1;
  unsigned int n_phen = params.n_pheno;

  if(pheno_data.cov_phenotypes.rows() != n_samples) throw std::runtime_error("#rows in cov_phenotypes");
  if(pheno_data.cov_phenotypes.cols() != n_cov + 2*n_phen + 2) throw std::runtime_error("#rows in cov_phenotypes");

  // v2
  for(unsigned i = n_cov1, k = 0; k < n_phen; i++, k++) {
    pheno_data.cov_phenotypes.col(i) = mphen.Mask.select(res.col(k), 0.0);
  }
  for(unsigned i = n_cov1 + n_phen + 1, k = 0; k < n_phen; i++, k++) {
    pheno_data.cov_phenotypes.col(i) = mphen.Mask.select(res.col(k), 0.0);
  }
  /* pheno_data.cov_phenotypes.rightCols(n_phen + 1).leftCols(n_phen) = res; */
  /* pheno_data.cov_phenotypes.array().colwise() *= mphen.Mask.array().cast<double>().array(); */

  // v1
  /* if(!params.strict_mode) throw std::runtime_error("--strict mode is required for MultiPhen test"); */
  /* mphen.setup_x(pheno_data.masked_indivs.col(0), pheno_data.cov_phenotypes, n_cov, n_phen, true, false); // (ignored by MultiPhen) pos_intercept_first = true, pos_phen_first = false */
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
////    for LD computation
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void Data::ld_comp() {

  sout << "LD computation";

  string out;
  vector < string > out_split, tmp_str;
  // output files
  Files ofile;
  // use pointer to class since it contains non-copyable elements
  vector < std::shared_ptr<Files> > ofile_split;

  // set some parameters
  if( params.build_mask ) bm.prep_run(params, files);

#if defined(_OPENMP)
  sout << " with " << (params.streamBGEN? "fast " : "") << "multithreading using OpenMP";
#endif
  sout << endl;

  file_read_initialization(); // set up files for reading
  read_pheno_and_cov(&files, &params, &in_filters, &pheno_data, &m_ests, &Gblock, sout);   // read phenotype and covariate files
  prep_run(&files, &in_filters, &params, &pheno_data, &m_ests, sout); // check blup files and adjust for covariates
  if(params.build_mask)
    set_groups_for_testing();   // set groups of snps to test jointly
  else
    set_blocks_for_testing();   // set number of blocks
  print_usage_info(&params, &files, sout);
  print_test_info();
  setup_output(&ofile, out, ofile_split, out_split); // result file
  sout << endl;

  // start analyzing each chromosome
  initialize_thread_data(Gblock.thread_data, params);
  params.ld_n = params.extract_vars_order.size();

  if(params.dosage_mode) // with dosages, avoid using sparse matrix
    compute_ld_dosages(&ofile);
  else // hard-calls only so use sparse matrix
    compute_ld_hardcalls(&ofile);

  return;
}

void Data::get_G_indices(ArrayXi& indices_ld, map<string, int>& colnames_Gmat){

  map<string, uint32_t >::iterator itr;
  int i_absent = colnames_Gmat.size();
  for (itr = params.extract_vars_order.begin(); itr != params.extract_vars_order.end(); ++itr) 
    if(in_map(itr->first, colnames_Gmat))
      indices_ld(itr->second) = colnames_Gmat[ itr->first ];
    else
      indices_ld(itr->second) = i_absent++; // cols for absent sv/masks are the same (ie 0 vector)

}

void Data::write_snplist(ArrayXb& is_absent){

  string const out = files.out_file + ".corr.snplist";
  map<string, uint32_t >::iterator itr;
  Files ofile;
  IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","\n");

  Eigen::Array<std::string,Eigen::Dynamic,1> ID_sorted (params.extract_vars_order.size());
  for (itr = params.extract_vars_order.begin(); itr != params.extract_vars_order.end(); ++itr) 
      ID_sorted( itr->second ) = itr->first;
  // write SNP list
  ofile.openForWrite(out, sout);
  ofile << ID_sorted.format(Fmt);
  ofile.closeFile();

  if(is_absent.any()){
    sout << " WARNING: there were variants" << (params.build_mask ? "/masks" : "") << " not found in the data; these were kept in the LD matrix.\n" <<
      "  + list is written to [" << files.out_file << ".corr.forcedIn.snplist]\n";
    ofile.openForWrite(files.out_file + ".corr.forcedIn.snplist", sout);
    ofile << ID_sorted(get_true_indices(is_absent)).format(Fmt);
    ofile.closeFile();
  }

}

void Data::compute_ld_dosages(Files* ofile){

  ArrayXb ld_var_absent = ArrayXb::Constant(params.ld_n, true);
  map<string, int> colnames_ld_mat;// to track id of cols in full_mat
  ArrayXi indices_ld(params.ld_n);

  MatrixXd LD = MatrixXd::Zero(params.ld_n, params.ld_n);

  // LD matrix will have first SVs then the burden masks
  for(size_t isnp = 0; isnp < params.ld_sv_offsets.size(); isnp++) {
    uint32_t snp_index = params.ld_sv_offsets[isnp];
    colnames_ld_mat[ snpinfo[ snp_index ].ID ] = colnames_ld_mat.size();
    ld_var_absent(params.extract_vars_order[snpinfo[ snp_index ].ID]) = false;
  }

  // build and read in masks (use sparse matrix)
  SpMat Gmask;
  MatrixXd Gmask_X;
  if(params.build_mask) {
    get_G_masks(Gmask, ld_var_absent, colnames_ld_mat);
    // project covariates
    Gmask_X = Gmask.transpose() * pheno_data.new_cov; // MxK
  }

  // to set columns of LD mat in right order 
  get_G_indices(indices_ld, colnames_ld_mat);

  // compute LD blocks for SVs with burden masks
  int nblocks_sv = ceil( params.ld_sv_offsets.size() * 1.0 / params.block_size );
  if(params.debug) cout << print_mem() << "\n";
  sout << "** Computing LD matrix " << (params.skip_scaleG ? "(=GtG) " : "") << "**\n";
  if(nblocks_sv > 0) sout << "  -> splitting across " << nblocks_sv << " SV blocks\n";
  MeasureTime mt_chunk;

  for(int snp_row = 0; snp_row < nblocks_sv; snp_row++){
    int row_start = params.block_size * snp_row;
    int row_nsnps = (snp_row == (nblocks_sv - 1)) ? (params.ld_sv_offsets.size() - snp_row *  params.block_size) : params.block_size;
    MatrixXd Grow(params.n_samples, row_nsnps);

    sout << "     - row " << snp_row + 1 << "\n";
    mt_chunk.start_ms();

    // read in G
    get_G_svs(snp_row, row_nsnps);
    Grow = Gblock.Gmat;
    // project covariates
    MatrixXd GtX_row = Grow.transpose() * pheno_data.new_cov; // MxK
    // compute diagonal
    LD.block(row_start, row_start, row_nsnps, row_nsnps).noalias() = -GtX_row * GtX_row.transpose();
    LD.block(row_start, row_start, row_nsnps, row_nsnps) += Grow.transpose() * Grow;
    sout << "       -> LD diagonal block computation..." << (params.debug ? print_mem() : "") << "..." << mt_chunk.stop_ms() << "\n";
    if((nblocks_sv > 1) && (snp_row < (nblocks_sv - 1))) {
      sout << "       -> computing LD with other variants (" << nblocks_sv - snp_row - 1 << " blocks)... " << flush;
      mt_chunk.start_ms();
    }

    for(int snp_col = (snp_row + 1); snp_col < nblocks_sv; snp_col++){
      int col_start = params.block_size * snp_col;
      int col_nsnps = (snp_col == (nblocks_sv - 1)) ? (params.ld_sv_offsets.size() - snp_col *  params.block_size) : params.block_size;
      if(params.debug) sout << snp_col - snp_row << "..." << flush;

      get_G_svs(snp_col, col_nsnps);
      // project covariates
      MatrixXd GtX_col = Gblock.Gmat.transpose() * pheno_data.new_cov; // MxK

      // compute ld block
      LD.block(row_start, col_start, row_nsnps, col_nsnps).noalias() = -GtX_row * GtX_col.transpose();
      LD.block(row_start, col_start, row_nsnps, col_nsnps) += Grow.transpose() * Gblock.Gmat;
    }
    if((nblocks_sv > 1) && (snp_row < (nblocks_sv - 1))) sout << mt_chunk.stop_ms() << "\n";

    // compute LD block for burden mask
    if(Gmask.cols() > 0){
      mt_chunk.start_ms();
      LD.block(row_start, params.ld_sv_offsets.size(), row_nsnps, Gmask.cols()).noalias() = -GtX_row * Gmask_X.transpose();
      LD.block(row_start, params.ld_sv_offsets.size(), row_nsnps, Gmask.cols()) += Grow.transpose() * Gmask;
      sout << "       -> computing LD with burden masks..." << mt_chunk.stop_ms() << "\n";
    }
  }

  // compute LD diagonal block for burden mask
  if(Gmask.cols() > 0){
    mt_chunk.start_ms();
    LD.block(params.ld_sv_offsets.size(), params.ld_sv_offsets.size(), Gmask.cols(), Gmask.cols()).noalias() = -Gmask_X * Gmask_X.transpose();
    LD.block(params.ld_sv_offsets.size(), params.ld_sv_offsets.size(), Gmask.cols(), Gmask.cols()) += Gmask.transpose() * Gmask;
    sout << "     - computing LD between burden masks..." << mt_chunk.stop_ms() << "\n";
  }

  // write out LD matrix
  print_ld(LD, indices_ld, ld_var_absent, ofile);

}


void Data::get_G_masks(SpMat& Gmat, ArrayXb& is_absent, map<string, int>& colnames_Gmat){

  MeasureTime mt;
  int block = 0, chrom_nb, bs;
  vector< variant_block > block_info;
  MatrixXd burden_mat;

  sout << "** Building burden masks **\n";
  mt.start_ms();

  // start analyzing each chromosome
  for (auto const& chrom : files.chr_read){

    if( !in_map(chrom, chr_map) ) continue;
    chrom_nb = chr_map[chrom][1];
    // if no sets in chromosome, skip
    if(chrom_nb == 0)  continue;

    // go through each set
    for(int bb = 0; bb < chrom_nb ; bb++) {

      vector< vector < uchar > > snp_data_blocks;
      vector< uint32_t > insize, outsize;
      vector<int> indices_mask_keep;

      vset* set_info = &(jt.setinfo[chrom - 1][bb]);
      bs = set_info->snp_indices.size();

      sout << " set [" << block + 1 << "/" << params.total_n_block << "] : " << set_info->ID << " - " << bs << " variants..." << flush;

      // build the masks
      block_info.resize(bs);
      getMask(chrom, bb, snp_data_blocks, insize, outsize, block_info);

      // store only the ones used in ld matrix
      for(size_t mask = 0; mask < set_info->snp_indices.size(); mask++)
        if(in_map(snpinfo[ set_info->snp_indices[mask] ].ID, params.extract_vars_order)){
          is_absent(params.extract_vars_order[ snpinfo[ set_info->snp_indices[mask] ].ID ]) = false;
          colnames_Gmat[ snpinfo[ set_info->snp_indices[mask] ].ID ] = colnames_Gmat.size(); // burden are after SV
          indices_mask_keep.push_back(mask);
        }

      // store in Gmat
      if(indices_mask_keep.size() > 0){
        burden_mat.conservativeResize(Gblock.Gmat.rows(), burden_mat.cols() + indices_mask_keep.size());
        burden_mat.rightCols(indices_mask_keep.size()) = Gblock.Gmat(Eigen::placeholders::all, indices_mask_keep);
      }

      block++;
    }

  }

  Gmat = burden_mat.sparseView();

  if(params.debug) cout << print_mem() << "...";
  sout << " -> " << mt.stop_ms() << "\n\n";

}


void Data::get_G_svs(int const& sv_block, int const& bsize){

  int chrom; 
  if( bsize == 0) return;
  int nvar_read = params.block_size * sv_block;

  vector<variant_block> all_snps_info;
  vector< vector < uchar > > snp_data_blocks;
  vector< uint32_t > insize, outsize;

  allocate_mat(Gblock.Gmat, params.n_samples, bsize);
  all_snps_info.resize(bsize);
  chrom = snpinfo[ params.ld_sv_offsets[0] ].chrom;

  vector<uint64> indices (params.ld_sv_offsets.begin() + nvar_read, params.ld_sv_offsets.begin() + nvar_read + bsize);
  readChunk(indices, chrom, snp_data_blocks, insize, outsize, all_snps_info);

#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(int isnp = 0; isnp < bsize; isnp++) {

    uint32_t snp_index = indices[isnp];
    variant_block* block_info = &(all_snps_info[isnp]);

    // build genotype matrix
    if( ((params.file_type == "bgen") && params.streamBGEN) || params.file_type == "bed") 
      parseSNP(isnp, chrom, &(snp_data_blocks[isnp]), insize[isnp], outsize[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &snpinfo[snp_index], &Gblock, block_info, sout);

    // impute missing if present
    if(block_info->ns1 < params.n_analyzed){
      MapArXd Geno (Gblock.Gmat.col(isnp).data(), params.n_samples, 1);
      mean_impute_g(block_info->af1*2, Geno, in_filters.ind_in_analysis);
    }

  }
#if defined(_OPENMP)
  setNbThreads(params.threads);
#endif

}

void Data::print_ld(MatrixXd& LDmat, ArrayXi& indices_ld, ArrayXb& is_absent, Files* ofile){

  MeasureTime mt;
  int bits = 16; // break [0,1] into 2^bits intervals
  double mult = (1ULL << bits) - 1; // map to 0,...,2^bits-1

  // write list of snps to file (corresponding to columns in LD matrix)
  write_snplist(is_absent);

  // only upper tri is loaded
  if(params.debug) cout << "     - raw covariance matrix[1:5,1:5]:\n" << LDmat.block(0,0,min(params.ld_n,5),min(params.ld_n,5)) << "\n" << print_mem() << "\n";

  // check if any of the diagonal entries are negative (but numerically zero -- due to rounding error)
  ArrayXb sd_G_zero = (LDmat.diagonal().array() < 0) && (LDmat.diagonal().array().abs() < params.tol) ;
  if(sd_G_zero.any()) {// set entries in LD matrix to 0
    ArrayXi ind_0 = get_true_indices(sd_G_zero);
    LDmat(ind_0,all).array() = 0; LDmat(all,ind_0).array() = 0;
  }

  if(!params.skip_scaleG) { // get cormat
    ArrayXd sds = (LDmat.diagonal().array() <= 0).select(sqrt(params.numtol), LDmat.diagonal().array().sqrt()); // bug fix for negative but numerically zero diagonal entries
    LDmat.diagonal().array() = sds.square();
  if(params.debug) cout << "     - thresholded covariance matrix[1:5,1:5]:\n" << LDmat.block(0,0,min(params.ld_n,5),min(params.ld_n,5)) << "\n" << print_mem() << "\n";
    LDmat = (1/sds).matrix().asDiagonal() * LDmat * (1/sds).matrix().asDiagonal();
  if(params.debug) cout << "     - correlation matrix[1:5,1:5]:\n" << LDmat.block(0,0,min(params.ld_n,5),min(params.ld_n,5)) << "\n" << print_mem() << "\n";
  } else 
    LDmat.diagonal().array() = LDmat.diagonal().array().max(params.numtol);

  // print corr
  sout << "     - writing to file...";
  mt.start_ms();

  if(params.ld_sparse_thr > 0){ // apply sparse threshold to LD matrix for off diagonal entries

    double out_val;
    // first diagonal entries (single line)
    ArrayXd sds = LDmat.diagonal().array().sqrt();
    IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","\n");
    (*ofile) << sds(indices_ld).matrix().transpose().format(Fmt);
    // off diagonal entries above thr based on corr (fmt = row/col/value [1-based])
    for(int i = 0; i < LDmat.rows(); i++)
      for(int j = i+1; j < LDmat.cols(); j++){
        if( indices_ld(i) < indices_ld(j) )
          out_val = LDmat(indices_ld(i),indices_ld(j)) / sds(indices_ld(i)) / sds(indices_ld(j));
        else
          out_val = LDmat(indices_ld(j),indices_ld(i)) / sds(indices_ld(i)) / sds(indices_ld(j));
        if(fabs(out_val) >= params.ld_sparse_thr)
          (*ofile) << i+1 << " " << j+1 << " " << out_val << "\n";
      }
    ofile->closeFile();

  } else if(params.cor_out_txt){

    IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","");
    MatrixXd full_LDmat = LDmat.selfadjointView<Eigen::Upper>();
    (*ofile) << full_LDmat(indices_ld, indices_ld).format(Fmt);
    ofile->closeFile();

  } else {

    ArrayXt vals;
    vals.resize( (LDmat.rows() * (LDmat.rows() - 1)) / 2 ); // m choose 2

    for(int i = 0, k = 0; i < LDmat.rows(); i++)
      for(int j = i+1; j < LDmat.cols(); j++)
        if( indices_ld(i) < indices_ld(j) )
          vals(k++) = LDmat(indices_ld(i),indices_ld(j)) * LDmat(indices_ld(i),indices_ld(j)) * mult + 0.5; // round to nearest integer
        else
          vals(k++) = LDmat(indices_ld(j),indices_ld(i)) * LDmat(indices_ld(j),indices_ld(i)) * mult + 0.5; // round to nearest integer

    //cerr << "\norig:\n" << LDmat.block(0,0,5,5).array().square().matrix() << "\nbin:\n" << 
     // vals.head(5) << "\n-->" << vals.size() << endl;

    ofile->writeBinMode(vals, sout);
    ofile->closeFile();
  }

  sout << " -> " << mt.stop_ms() << "\n";

  exit_early();

}


void Data::compute_ld_hardcalls(Files* ofile){

  ArrayXb ld_var_absent = ArrayXb::Constant(params.ld_n, true);
  map<string, int> colnames_ld_mat;// to track id of cols in full_mat
  ArrayXi indices_ld(params.ld_n);
	SpMat Gmat(params.n_samples, params.ld_n);

	// read in SVs
	get_G_svs(Gmat, ld_var_absent, colnames_ld_mat);

	// read in masks
	if(params.build_mask) get_G_masks_hc(Gmat, ld_var_absent, colnames_ld_mat);

	// to set columns of LD mat in right order 
	get_G_indices(indices_ld, colnames_ld_mat);

	// compute LD matrix
	print_ld(Gmat, indices_ld, ld_var_absent, ofile);

}

void Data::get_G_svs(SpMat& Gmat, ArrayXb& is_absent, map<string, int>& colnames_Gmat){

  int n_snps = params.ld_sv_offsets.size();
  if( n_snps == 0) return;

  bool last_chunk = false;
  int nchunks, bsize, chrom, nvar_read = 0; 
  vector<variant_block> all_snps_info;
  vector< vector < uchar > > snp_data_blocks;
  vector< uint32_t > insize, outsize;

  // read in variants in chunks storing as sparse matrix
  nchunks = ceil( n_snps * 1.0 / params.block_size );
  sout << "** reading in single variant genotypes **\n  + " << n_snps << " variants in total split across " << nchunks << " blocks\n";
  MeasureTime mt, mt_chunk;
  mt.start_ms();
  if(params.debug) cerr << print_mem() << "...";

  // do it in chunks to reduce memory usage when reading as dense
  bsize = params.block_size; // default number of SNPs to read at a time
  allocate_mat(Gblock.Gmat, params.n_samples, bsize);
  all_snps_info.resize(bsize);
  chrom = snpinfo[ params.ld_sv_offsets[0] ].chrom;
  for(int i = 0; i < nchunks; i++){

    sout << "  block [" << i + 1 << "/" << nchunks << "] : reading in genotypes..." << flush;
    mt_chunk.start_ms();

    last_chunk = ( i == (nchunks-1) );
    if( last_chunk ) {
      bsize = n_snps - i * bsize;// use remainder number of variants
      allocate_mat(Gblock.Gmat, params.n_samples, bsize);
    }

    vector<uint64> indices (params.ld_sv_offsets.begin() + nvar_read, params.ld_sv_offsets.begin() + nvar_read + bsize);
    readChunk(indices, chrom, snp_data_blocks, insize, outsize, all_snps_info);

#if defined(_OPENMP)
    setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
    for(int isnp = 0; isnp < bsize; isnp++) {

      uint32_t snp_index = indices[isnp];
      variant_block* block_info = &(all_snps_info[isnp]);

      // build genotype matrix
      if( ((params.file_type == "bgen") && params.streamBGEN) || params.file_type == "bed") 
        parseSNP(isnp, chrom, &(snp_data_blocks[isnp]), insize[isnp], outsize[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &snpinfo[snp_index], &Gblock, block_info, sout);

      // impute missing if present
      MapArXd Geno (Gblock.Gmat.col(isnp).data(), params.n_samples, 1);
      mean_impute_g(block_info->af1*2, Geno, in_filters.ind_in_analysis);

      // check if in LD matrix 
      is_absent(params.extract_vars_order[snpinfo[ snp_index ].ID]) = false;
    }
#if defined(_OPENMP)
    setNbThreads(params.threads);
#endif
    // convert to sparse
    sout << mt_chunk.stop_ms() << "...converting to sparse..."; mt_chunk.start_ms();
    Gmat.middleCols(nvar_read, bsize) = Gblock.Gmat.sparseView();
    // store ID in Gmat (can't do it multithreaded)
    for(int isnp = 0; isnp < bsize; isnp++) {
      uint32_t snp_index = indices[isnp];
      colnames_Gmat[ snpinfo[ snp_index ].ID ] = colnames_Gmat.size();
    }

    if(params.debug) cout << print_mem() << "...";
    sout << mt_chunk.stop_ms() << "\n";
    nvar_read += bsize;
  }

  if(params.debug) cout << print_mem() << "...";
  sout << " -> " << mt.stop_ms() << "\n";

}

void Data::get_G_masks_hc(SpMat& Gmat, ArrayXb& is_absent, map<string, int>& colnames_Gmat){

  MeasureTime mt;
  int block = 0, chrom_nb, bs;
  int nvar_read = colnames_Gmat.size();
  vector< variant_block > block_info;

  sout << "\n** Building burden masks **\n";
  mt.start_ms();

  // start analyzing each chromosome
  for (auto const& chrom : files.chr_read){

    if( !in_map(chrom, chr_map) ) continue;
    chrom_nb = chr_map[chrom][1];
    // if no sets in chromosome, skip
    if(chrom_nb == 0)  continue;

    // go through each set
    for(int bb = 0; bb < chrom_nb ; bb++) {

      vector< vector < uchar > > snp_data_blocks;
      vector< uint32_t > insize, outsize;
      vector<int> indices_mask_keep;

      vset* set_info = &(jt.setinfo[chrom - 1][bb]);
      bs = set_info->snp_indices.size();

      sout << " set [" << block + 1 << "/" << params.total_n_block << "] : " << set_info->ID << " - " << bs << " variants..." << flush;

      // build the masks
      block_info.resize(bs);
      getMask(chrom, bb, snp_data_blocks, insize, outsize, block_info);

      // store only the ones used in ld matrix
      for(size_t mask = 0; mask < set_info->snp_indices.size(); mask++)
        if(in_map(snpinfo[ set_info->snp_indices[mask] ].ID, params.extract_vars_order)){
          is_absent(params.extract_vars_order[ snpinfo[ set_info->snp_indices[mask] ].ID ]) = false;
          colnames_Gmat[ snpinfo[ set_info->snp_indices[mask] ].ID ] = colnames_Gmat.size();
          indices_mask_keep.push_back(mask);
        }

      // store in Gmat
      if(indices_mask_keep.size() > 0){
        Gmat.middleCols(nvar_read, indices_mask_keep.size()) = Gblock.Gmat(Eigen::placeholders::all, indices_mask_keep).sparseView();
        nvar_read += indices_mask_keep.size();
      }
      block++;
    }

  }

  if(params.debug) cout << print_mem() << "...";
  sout << " -> " << mt.stop_ms() << "\n";

}

void Data::print_ld(SpMat& Gmat, ArrayXi& indices_ld, ArrayXb& is_absent, Files* ofile){

  MeasureTime mt;
  int bits = 16; // break [0,1] into 2^bits intervals
  double mult = (1ULL << bits) - 1; // map to 0,...,2^bits-1

  sout << "\n** computing LD matrix " << (params.skip_scaleG ? "(=GtG) " : "") << "**\n";
  mt.start_ms();

  // write list of snps to file (corresponding to columns in LD matrix)
  write_snplist(is_absent);

  // get LD matrix - first project covariates
  MatrixXd GtX = Gmat.transpose() * pheno_data.new_cov; // MxK
  MatrixXd LDmat = -GtX * GtX.transpose();
  LDmat += Gmat.transpose() * Gmat;
  if(params.debug) cout << "     - raw covariance matrix[1:5,1:5]:\n" << LDmat.block(0,0,min(params.ld_n,5),min(params.ld_n,5)) << "\n" << print_mem() << "\n";

  // check if any of the diagonal entries are negative (but numerically zero -- due to rounding error)
  ArrayXb sd_G_zero = (LDmat.diagonal().array() < 0) && (LDmat.diagonal().array().abs() < params.tol) ;
  if(sd_G_zero.any()) {// set entries in LD matrix to 0
    ArrayXi ind_0 = get_true_indices(sd_G_zero);
    LDmat(ind_0,all).array() = 0; LDmat(all,ind_0).array() = 0;
  }

  if(!params.skip_scaleG) { // get cormat
    ArrayXd sds = (LDmat.diagonal().array() <= 0).select(sqrt(params.numtol), LDmat.diagonal().array().sqrt()); // bug fix for negative but numerically zero diagonal entries
    LDmat.diagonal().array() = sds.square();
  if(params.debug) cout << "     - thresholded covariance matrix[1:5,1:5]:\n" << LDmat.block(0,0,min(params.ld_n,5),min(params.ld_n,5)) << "\n" << print_mem() << "\n";
    LDmat = (1/sds).matrix().asDiagonal() * LDmat * (1/sds).matrix().asDiagonal();
  if(params.debug) cout << "     - correlation matrix[1:5,1:5]:\n" << LDmat.block(0,0,min(params.ld_n,5),min(params.ld_n,5)) << "\n" << print_mem() << "\n";
  } else 
    LDmat.diagonal().array() = LDmat.diagonal().array().max(params.numtol);
  sout << " -> " << mt.stop_ms() << "\n";

  // print corr
  sout << "\n** writing to file **\n";
  mt.start_ms();

  if(params.ld_sparse_thr > 0){ // apply sparse threshold to LD matrix for off diagonal entries

    double out_val;
    // first diagonal entries (single line)
    ArrayXd sds = LDmat.diagonal().array().sqrt();
    IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","\n");
    (*ofile) << sds(indices_ld).matrix().transpose().format(Fmt);
    // off diagonal entries above thr based on corr (fmt = row/col/value [1-based])
    for(int i = 0; i < LDmat.rows(); i++)
      for(int j = i+1; j < LDmat.cols(); j++){
        out_val = LDmat(indices_ld(i),indices_ld(j)) / sds(indices_ld(i)) / sds(indices_ld(j));
        if(fabs(out_val) >= params.ld_sparse_thr)
          (*ofile) << i+1 << " " << j+1 << " " << out_val << "\n";
      }
    ofile->closeFile();

  } else if(params.cor_out_txt){

    IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","");
    (*ofile) << LDmat(indices_ld, indices_ld).format(Fmt);
    ofile->closeFile();

  } else {

    ArrayXt vals;
    vals.resize( (LDmat.rows() * (LDmat.rows() - 1)) / 2 ); // m choose 2

    for(int i = 0, k = 0; i < LDmat.rows(); i++)
      for(int j = i+1; j < LDmat.cols(); j++)
        vals(k++) = LDmat(indices_ld(i),indices_ld(j)) * LDmat(indices_ld(i),indices_ld(j)) * mult + 0.5; // round to nearest integer

    //cerr << "\norig:\n" << LDmat.block(0,0,5,5).array().square().matrix() << "\nbin:\n" << 
     // vals.head(5) << "\n-->" << vals.size() << endl;

    ofile->writeBinMode(vals, sout);
    ofile->closeFile();
  }

  sout << " -> " << mt.stop_ms() << "\n";

  exit_early();

}

