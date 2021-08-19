/*

   This file is part of the regenie software package.

   Copyright (c) 2020-2021 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

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
#include "Files.hpp"
#include "Geno.hpp"
#include "Joint_Tests.hpp"
#include "Step1_Models.hpp"
#include "Step2_Models.hpp"
#include "Pheno.hpp"
#include "HLM.hpp"
#include "Interaction.hpp"
#include "Masks.hpp"
#include "Data.hpp"

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
  // level 1 ridge
  if(params.trait_mode == 0){ // QT
    if(params.use_loocv) ridge_level_1_loocv(&files, &params, &pheno_data, &l1_ests, sout);
    else ridge_level_1(&files, &params, &l1_ests, sout);
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

  if( params.snp_set ) test_joint();
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
  if(params.getCorMat) params.block_size = params.n_variants;

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
      throw "number of variants/blocks in file don't match with that in master file.";
  } else if(params.run_l1_only) prep_parallel_l1();

  // set ridge params
  params.lambda = (params.run_l0_only ? params.parallel_nGeno : params.n_variants) * (1 - params.lambda) / params.lambda;
  if(params.trait_mode == 2){
    ArrayXd base_tau = params.tau[0];
    params.tau.assign(params.n_pheno, base_tau);
    for(int i = 0; i < params.n_pheno; i++){
      double rate = pheno_data.phenotypes_raw.col(i).sum() / pheno_data.Neff(i); // masked entries are 0
      params.tau[i] = (params.total_n_block * params.n_ridge_l0) / (1 + params.tau[i] / (rate * (1 - params.tau[i]))).log();
      //cerr << endl << params.tau[i].matrix().transpose() << endl;
    }
  } else {
    params.tau[0] = (params.total_n_block * params.n_ridge_l0) * (1 - params.tau[0]) / params.tau[0];
    // Assuming input tau is total SNP heritability on the liability scale= m * 3/pi^2 * (1-h2) / h2
    if(params.trait_mode == 1) params.tau[0] *= 3 / (M_PI * M_PI);
  }

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
    int nv_tot = (params.run_l0_only ? params.parallel_nGeno : params.n_variants);
    sout << left << std::setw(20) << " * ridge data_l0" << ": [" << params.n_ridge_l0 << " : ";
    for(int i = 0; i < params.lambda.size(); i++)
      sout << nv_tot / ( nv_tot + params.lambda(i)) << " ";
    sout << "]\n";
  }

  if(!params.run_l0_only){
    sout << left << std::setw(20) << " * ridge data_l1" << ": [" << params.n_ridge_l1 << " : ";
    for(int i = 0; i < params.tau[0].size(); i++)
      if(params.trait_mode == 2){
        double rate = pheno_data.phenotypes_raw.col(0).sum() / pheno_data.Neff(0); // only use trait 1
        double zv = exp(params.total_n_block * params.n_ridge_l0 / params.tau[0](i)) - 1; 
        sout <<  rate * zv / (1 + rate * zv) << " ";
      } else 
        sout << (params.total_n_block * params.n_ridge_l0) / (params.total_n_block * params.n_ridge_l0 + (params.trait_mode == 1 ? (M_PI * M_PI / 3) : 1) * params.tau[0](i)) << " ";
    sout << "]\n";
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
    MatrixXd phenos = ( pheno_data.phenotypes_raw.array() * pheno_data.masked_indivs.array().cast<double>()).matrix();

    for(int i = 0; i < params.cv_folds; i++) {
      sum = phenos.block(cum_size_folds,0,params.cv_sizes(i),params.n_pheno).colwise().sum();

      // BTs
      if(params.trait_mode == 1){
        n_cv = pheno_data.masked_indivs.block(cum_size_folds,0,params.cv_sizes(i),params.n_pheno).cast<double>().colwise().sum();
        sd_phenos = (sum/n_cv) * (1 - sum/n_cv);

        if( sd_phenos.minCoeff(&minIndex) < params.numtol )
          throw "one of the folds has only cases/controls for phenotype '" + files.pheno_names[minIndex] 
            + "'. Either use smaller #folds (option --cv) or use LOOCV (option --loocv).";
      } else if(params.trait_mode == 2){

        if( sum.maxCoeff(&minIndex) <= 0 )
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

      if(i == 0) masked_in_folds[j] = MatrixXb::Constant(params.cv_sizes(j), params.n_pheno, false);

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
      openStream(files.write_preds_files[ph].get(), fout_p, ios::out | ios::binary, sout);
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

      get_block_size(params.block_size, chrom_nsnps, bb, bs);

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
  if(params.write_l0_pred)
    for( auto &yfile : files.write_preds_files)
      yfile->close();

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
  if(params.write_l0_pred && (params.n_pheno > 1) ) {
    // free level 0 predictions for (P-1) indices in test_mat
    for(int ph = 1; ph < params.n_pheno; ++ph ) {
      if(!params.use_loocv){ // k-fold
        for(int i = 0; i < params.cv_folds; ++i ) 
          l1_ests.test_mat[ph][i].resize(0,0);
        l1_ests.test_mat[ph].resize(0);
      } else { // loocv
        l1_ests.test_mat_conc[ph].resize(0,0);
      }
    }
  }

}

void Data::calc_cv_matrices(int const& bs, struct ridgel0* l0) {

  sout << "   -calc working matrices..." << flush;
  auto t2 = std::chrono::high_resolution_clock::now();

  if(!params.use_loocv){ // k-fold

    l0->GGt.setZero(bs,bs);
    l0->GTY.setZero(bs,params.n_pheno);
    uint32_t cum_size_folds = 0;

    for( int i = 0; i < params.cv_folds; ++i ) {
      l0->GtY[i] = Gblock.Gmat.block(0, cum_size_folds, bs, params.cv_sizes(i)) * pheno_data.phenotypes.block(cum_size_folds, 0, params.cv_sizes(i), params.n_pheno);
      l0->GTY += l0->GtY[i];
      l0->G_folds[i] = Gblock.Gmat.block(0, cum_size_folds, bs, params.cv_sizes(i)) * Gblock.Gmat.block(0, cum_size_folds, bs, params.cv_sizes(i)).transpose();
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
    throw "invalid header line."; 

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
    throw "invalid header line."; 

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
  double rate, zv;
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

      } else { // for binary traits - check level 1 ridge converged

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
        zv = exp(params.total_n_block * params.n_ridge_l0 / params.tau[ph](j)) - 1; 
        sout << "  " << setw(5) << rate * zv / (1 + rate * zv);
      } else sout << "  " << setw(5) << (params.total_n_block *  params.n_ridge_l0) / (params.total_n_block * params.n_ridge_l0 + (params.trait_mode == 1? (M_PI * M_PI / 3) : 1) * params.tau[0](j) );

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
      if(params.use_loocv) 
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
      remove(pfile.c_str()); // l0 predictions
      if(ph==0){
        pfile = pfx + ".snplist";
        remove(pfile.c_str()); // snplist
      }
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

  int bs_l1 = params.total_n_block * params.n_ridge_l0;
  int ph_eff = params.write_l0_pred ? 0 : ph;
  string outname, in_pheno;
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
    beta_l1 = (X1 + params.tau[0](val) * ident_l1).llt().solve(X2);
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

  // sout << "\nFor tau[" << val <<"] = " << params.tau[0](val) << endl <<  beta_l1 << endl ;
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

  int bs_l1 = params.total_n_block * params.n_ridge_l0;
  int ph_eff = params.write_l0_pred ? 0 : ph;
  string in_pheno;
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
  dl_inv = (eigX.eigenvalues().array() + params.tau[0](val)).inverse();
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
void Data::make_predictions_binary(int const& ph, int const& val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  int bs_l1 = params.total_n_block * params.n_ridge_l0;
  int ph_eff = params.write_l0_pred ? 0 : ph;
  string in_pheno;
  ArrayXd etavec, pivec, wvec, zvec, score;
  MatrixXd betaold, betanew, XtW, XtWX, XtWZ;
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
      betanew = (XtWX + params.tau[0](val) * ident_l1).llt().solve(XtWZ);
      // compute score
      score = ArrayXd::Zero(betanew.rows());
      for(int i = 0; i < params.cv_folds; ++i ) {
        etavec = (l1_ests.test_offset[ph][i] + l1_ests.test_mat[ph_eff][i] * betanew).array();
        pivec = 1 - 1/(etavec.exp() + 1);
        score += (l1_ests.test_mat[ph_eff][i].transpose() * (l1_ests.test_pheno_raw[ph][i].array() - pivec).matrix()).array();
      }
      score -= params.tau[0](val) * betanew.array();

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


void Data::make_predictions_binary_loocv(int const& ph, int const& val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  int bs_l1 = params.total_n_block * params.n_ridge_l0;
  int ph_eff = params.write_l0_pred ? 0 : ph;
  double v2;
  string in_pheno;

  ArrayXd beta, pivec, wvec;
  MatrixXd XtWX, V1, beta_final;
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

  MapArXd Y (pheno_data.phenotypes_raw.col(ph).data(), pheno_data.phenotypes_raw.rows());
  MapMatXd X (l1_ests.test_mat_conc[ph_eff].data(), pheno_data.phenotypes_raw.rows(), bs_l1);
  MapArXd offset (m_ests.offset_nullreg.col(ph).data(), pheno_data.phenotypes_raw.rows());
  MapArXb mask (pheno_data.masked_indivs.col(ph).data(), pheno_data.phenotypes_raw.rows());

  // fit logistic on whole data again for optimal ridge param
  beta = ArrayXd::Zero(bs_l1);
  run_log_ridge_loocv(params.tau[0](val), target_size, nchunk, beta, pivec, wvec, Y, X, offset, mask, &params, sout);

  // compute Hinv
  //zvec = (etavec - m_ests.offset_nullreg.col(ph).array()) + (pheno_data.phenotypes_raw.col(ph).array() - pivec) / wvec;
  XtWX = MatrixXd::Zero(bs_l1, bs_l1);
  for(chunk = 0; chunk < nchunk; ++chunk){
    size_chunk = ( chunk == nchunk - 1 ? params.cv_folds - target_size * chunk : target_size );
    j_start = chunk * target_size;

    Ref<MatrixXd> Xmat_chunk = X.block(j_start, 0, size_chunk, bs_l1); // n x k
    Ref<MatrixXd> w_chunk = wvec.matrix().block(j_start, 0, size_chunk,1);

    XtWX += Xmat_chunk.transpose() * w_chunk.asDiagonal() * Xmat_chunk;
  }
  Hinv.compute( XtWX + params.tau[0](val) * ident_l1 );

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
      beta_final.col(i) = (beta - V1.col(i).array() * (Yvec_chunk(i,0) - pivec(j_start + i)) / (1 - v2)).matrix();
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


// predictions for count phenotypes
void Data::make_predictions_count(int const& ph, int const& val) {

  sout << "  * making predictions..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  int bs_l1 = params.total_n_block * params.n_ridge_l0;
  int ph_eff = params.write_l0_pred ? 0 : ph;
  string in_pheno;
  ArrayXd etavec, pivec, zvec, score;
  MatrixXd betaold, betanew, XtW, XtWX, XtWZ;
  MatrixXd ident_l1 = MatrixXd::Identity(bs_l1,bs_l1);

  // read in level 0 predictions from file
  if(params.write_l0_pred)
    read_l0(ph, ph_eff, &files, &params, &l1_ests, sout);

  // fit model using out-of-sample level 0 predictions from whole data
  if(params.within_sample_l0)
    throw "--within is not supported for count phenotypes";

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

  int bs_l1 = params.total_n_block * params.n_ridge_l0;
  int ph_eff = params.write_l0_pred ? 0 : ph;
  double v2;
  string in_pheno;

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
  if(params.write_l0_pred)
    read_l0(ph, ph_eff, &files, &params, &l1_ests, sout);

  MapArXd Y (pheno_data.phenotypes_raw.col(ph).data(), pheno_data.phenotypes_raw.rows());
  MapMatXd X (l1_ests.test_mat_conc[ph_eff].data(), pheno_data.phenotypes_raw.rows(), bs_l1);
  MapArXd offset (m_ests.offset_nullreg.col(ph).data(), pheno_data.phenotypes_raw.rows());
  MapArXb mask (pheno_data.masked_indivs.col(ph).data(), pheno_data.phenotypes_raw.rows());

  // fit logistic on whole data again for optimal ridge param
  beta = ArrayXd::Zero(bs_l1);
  run_ct_ridge_loocv(params.tau[ph](val), target_size, nchunk, beta, pivec, Y, X, offset, mask, &params, sout);

  // compute Hinv
  //zvec = (etavec - m_ests.offset_nullreg.col(ph).array()) + (pheno_data.phenotypes_raw.col(ph).array() - pivec) / wvec;
  XtWX = MatrixXd::Zero(bs_l1, bs_l1);
  for(chunk = 0; chunk < nchunk; ++chunk){
    size_chunk = ( chunk == nchunk - 1 ? params.cv_folds - target_size * chunk : target_size );
    j_start = chunk * target_size;

    Ref<MatrixXd> Xmat_chunk = X.block(j_start, 0, size_chunk, bs_l1); // n x k
    Ref<MatrixXd> w_chunk = pivec.matrix().block(j_start, 0, size_chunk,1);

    XtWX += Xmat_chunk.transpose() * w_chunk.asDiagonal() * Xmat_chunk;
  }
  Hinv.compute( XtWX + params.tau[ph](val) * ident_l1 );

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
      if(!fit_firth(ph, Y, pheno_data.new_cov, pred.col(chr).array(), mask, pivec, etavec, bhat, se, params.ncov, dev, false, tstat, params.maxstep_null, params.niter_max_firth_null, params.numtol, &params)){
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


std::string Data::write_chr_row(int const& chr, int const& ph, const Ref<const MatrixXd>& pred){

  uint32_t index;
  string out, id_index;
  std::ostringstream buffer;
  map<string, uint32_t >::iterator itr_ind;

  buffer << chr << " ";
  for (itr_ind = params.FID_IID_to_ind.begin(); itr_ind != params.FID_IID_to_ind.end(); ++itr_ind) {

    // check individual was included in analysis, if not then skip
    index = itr_ind->second;
    if( !in_filters.ind_in_analysis( index ) ) continue;

    // print prs
    id_index = itr_ind->first;
    if( pheno_data.masked_indivs(index, ph) )
      buffer << pred(index, 0) << " ";
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
    if(params.cor_out_txt){
      sout << " * computing correlation matrix\n  + output to text file ["<<out<<"]\n";
      ofile->openForWrite(out, sout);
    } else {
      sout << " * computing correlation matrix (storing R^2 values)\n  + output to binary file ["<<out<<"]\n";
      ofile->openMode(out, std::ios_base::out | std::ios_base::binary, sout);
      ArrayXi vals(2);
      vals << params.n_samples , params.n_variants;
      //cerr << vals << endl;
      ofile->writeBinMode(vals, sout);
    }
    sout << "  + list of snps written to [" << out << ".snplist]\n";
    sout << "  + n_snps = " << params.n_variants <<"\n\n";
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
  if(params.setMinINFO) 
    sout << " * using minimum imputation info score of " << params.min_INFO << " (variants with lower info score are ignored)\n";

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


  if(params.htp_out){
    if((params.trait_mode==1) & params.firth) correction_type = "-FIRTH";
    else if((params.trait_mode==1) & params.use_SPA) correction_type = "-SPA";
    else if(params.trait_mode==1) correction_type = "-LOG";
    else if(params.trait_mode==2) correction_type = "-POISSON";
    else correction_type = "-LR";

    if(params.skip_blups) model_type = test_string + correction_type;
    else model_type = test_string + "-WGR" + correction_type;
  }

  if(params.gwas_condtl) // specify main sum stats is conditional gwas
    params.condtl_suff = "-CONDTL";

  params.with_flip = params.with_flip && !params.build_mask && params.trait_mode && (params.test_type == 0);

  if( params.joint_test ) {
    params.with_flip = jt.get_test_info(&params, test_string, sout) && params.with_flip;
    jt.out_file_prefix = files.out_file;
  }

  normal nd(0,1);
  chi_squared chisq(1);
  params.zcrit = quantile(complement(nd, .025));
  params.chisq_thr = quantile(chisq, 1 - params.alpha_pvalue);
  params.z_thr = sqrt(params.chisq_thr);

}


/// When computing & outputting LD info
// write list of variants used to compute LD
void Data::write_snplist(int const& bs, vector<variant_block> const& all_snps_info){

  bool reject = false;
  string const out = files.out_file + ".corr.snplist";
  string const out_reject = files.out_file + ".corr.exclude";
  std::ostringstream buffer, buffer_reject;
  Files ofile;

  for(int snp = 0; snp < bs; snp++) {
    if(all_snps_info[snp].ignored){
      reject = true;
      buffer_reject << snpinfo[snp].ID << endl;
    } else buffer << snpinfo[snp].ID << endl;
  }

  if(reject){
    // write SNP list to ignore
    ofile.openForWrite(out_reject, sout);
    ofile << buffer_reject.str();
    ofile.closeFile();

    throw "SNPs with low variance are present. Use '--exclude " + out_reject + "' to ignore them.";
  }

  // write SNP list
  ofile.openForWrite(out, sout);
  ofile << buffer.str();
  ofile.closeFile();

}

void Data::print_cor(Files* ofile){

  int bits = 16; // break [0,1] into 2^bits intervals
  double mult = (1ULL << bits) - 1; // map to 0,...,2^bits-1

  // get LD matrix
  sout << "   -computing LD matrix..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  MatrixXd LDmat = (Gblock.Gmat.transpose() * Gblock.Gmat) / (params.n_samples - params.ncov);

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms)\n   -writing to file..." << flush;


  // print corr
  if(params.cor_out_txt){

    IOFormat Fmt(StreamPrecision, DontAlignCols, " ", "\n", "", "","","");
    (*ofile) << LDmat.format(Fmt);
    ofile->closeFile();

  } else {

    ArrayXt vals;
    vals.resize( (Gblock.Gmat.cols() * (Gblock.Gmat.cols() - 1)) / 2 );

    for(int i = 0, k = 0; i < LDmat.rows(); i++)
      for(int j = i+1; j < LDmat.cols(); j++){
        vals(k++) = LDmat(i,j) * LDmat(i,j) * mult + 0.5; // round to nearest integer
      }

    //cerr << "\norig:\n" << LDmat.block(0,0,5,5).array().square().matrix() << "\nbin:\n" << 
     // vals.head(5) << "\n-->" << vals.size() << endl;

    ofile->writeBinMode(vals, sout);
    ofile->closeFile();
  }

  sout << "done";
  t1 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t2);
  sout << " (" << duration.count() << "ms) "<< endl;

  exit_early();

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

  if(params.getCorMat && (nchr > 1 || params.n_variants < 2))
    throw "can only compute LD matrix for a single chromosome (use --chr/--chrList/--range) and >=2 variants.";

  // summarize block sizes
  sout << left << std::setw(20) << " * # threads" << ": [" << params.threads << "]\n";
  sout << left << std::setw(20) << " * block size" << ": [" << params.block_size << "]\n";
  sout << left << std::setw(20) << " * # blocks" << ": [" << params.total_n_block << "]\n";
  if(params.start_block > 1) sout << "    + skipping to block #" << params.start_block << endl;

  // storing null estimates from firth
  if(params.use_null_firth) 
    sout << " * reading null Firth estimates using file : [" << files.null_firth_file << "]\n";
  if(params.write_null_firth ) 
    sout << " * writing null Firth estimates to file\n";

}

void Data::set_nullreg_mat(){

  m_ests.Y_hat_p = MatrixXd::Zero(params.n_samples, params.n_pheno);
  m_ests.Gamma_sqrt = MatrixXd::Zero(params.n_samples, params.n_pheno);
  m_ests.Xt_Gamma_X_inv.resize(params.n_pheno);
  m_ests.X_Gamma.resize(params.n_pheno);

  // for firth  approx
  if(params.firth_approx){ 
    firth_est.covs_firth = MatrixXd::Zero(params.n_samples, pheno_data.new_cov.cols() + 1);
    firth_est.covs_firth.leftCols(pheno_data.new_cov.cols()) = pheno_data.new_cov;
    firth_est.beta_null_firth = MatrixXd::Zero(firth_est.covs_firth.cols(), params.n_pheno);

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

      if(params.getCorMat) {
        write_snplist(bs, block_info);
        print_cor(&ofile);
      }

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

          if( snp_data.ignored_trait(j) ) 
            continue;

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

  if((params.file_type == "bgen") && params.streamBGEN){
    snp_data_blocks.resize( n_snps );
    insize.resize(n_snps); outsize.resize(n_snps);
    vector<uint64> offsets(n_snps);
    for (int i = 0; i < n_snps; i++) offsets[i] = snpinfo[indices[i]].offset;

    readChunkFromBGEN(&files.geno_ifstream, insize, outsize, snp_data_blocks, offsets);

  } else if((params.file_type == "bgen") && !params.streamBGEN) 
    readChunkFromBGENFileToG(indices, chrom, snpinfo, &params, &Gblock, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, all_snps_info, sout);
  else if(params.file_type == "pgen") 
    readChunkFromPGENFileToG(indices, chrom, &params, &in_filters, &Gblock, pheno_data.masked_indivs, pheno_data.phenotypes_raw, snpinfo, all_snps_info);
  else if(params.file_type == "bed"){

    snp_data_blocks.resize( n_snps );
    for(int isnp = 0; isnp < n_snps; isnp++) {

      jumpto_bed( snpinfo[indices[isnp]].offset, files.bed_block_size, files.geno_ifstream);
      snp_data_blocks[isnp].resize(files.bed_block_size);
      files.geno_ifstream.read( reinterpret_cast<char *> (&snp_data_blocks[isnp][0]), files.bed_block_size);

    }

  }

  // analyze using openmp
  compute_tests_mt(chrom, indices, snp_data_blocks, insize, outsize, all_snps_info);
  //compute_tests_st(chrom, indices, snp_data_blocks, insize, outsize, all_snps_info); // this is slower


  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << "done (" << duration.count() << "ms) "<< endl;
}


void Data::compute_res(){

  res = pheno_data.phenotypes - m_ests.blups;
  res.array() *= pheno_data.masked_indivs.array().cast<double>();

  p_sd_yres = res.colwise().norm();
  p_sd_yres.array() /= sqrt(pheno_data.Neff - params.ncov);
  res.array().rowwise() /= p_sd_yres.array();

  if(params.w_interaction && (params.trait_mode==0) && !params.no_robust && !params.force_robust) 
    HLM_fitNull(nullHLM, m_ests, pheno_data, files, params, sout);

}

void Data::compute_res_bin(int const& chrom){

  fit_null_logistic(chrom, &params, &pheno_data, &m_ests, &files, sout); // for all phenotypes

  res = pheno_data.phenotypes_raw - m_ests.Y_hat_p;
  res.array() /= m_ests.Gamma_sqrt.array();
  res.array() *= pheno_data.masked_indivs.array().cast<double>();
  if(params.debug) cerr << endl << res.topRows(4) << endl;

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

    if( !params.build_mask )
      parseSNP(isnp, chrom, &(snp_data_blocks[isnp]), insize[isnp], outsize[isnp], &params, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, &snpinfo[snp_index], &Gblock, block_info, sout);

    // check if g is sparse (not for QT without strict mode)
    if(params.trait_mode || params.strict_mode)
      check_sparse_G(isnp, thread_num, &Gblock, params.n_samples, in_filters.ind_in_analysis);

    if(params.w_interaction) {
      if(params.interaction_snp && (snpinfo[snp_index].ID == in_filters.interaction_cov))
        block_info->skip_int = true;
      get_interaction_terms(isnp, thread_num, &pheno_data, &Gblock, block_info, nullHLM, &params, sout);
    }

    // for QTs: project out covariates & scale
    residualize_geno(isnp, thread_num, block_info, false, pheno_data.new_cov, &Gblock, &params);

    // skip SNP if fails filters
    if( block_info->ignored || params.getCorMat ) continue;
    
    reset_stats(block_info, params);

    compute_score(isnp, snp_index, chrom, thread_num, test_string + params.condtl_suff, model_type + params.condtl_suff, res, p_sd_yres, params, pheno_data, Gblock, block_info, snpinfo, m_ests, firth_est, files, sout);

    // for joint test, store logp
    if( params.joint_test ) block_info->pval_log = Gblock.thread_data[thread_num].pval_log;

    if(params.w_interaction) 
      apply_interaction_tests(snp_index, isnp, thread_num, res, p_sd_yres, model_type, test_string, &pheno_data, nullHLM, &in_filters, &files, &Gblock, block_info, snpinfo, &m_ests, &firth_est, &params, sout);
  }

#if defined(_OPENMP)
  setNbThreads(params.threads);
#endif

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
  vector < string > out_split;
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

      sout << " set [" << block + 1 << "/" << params.total_n_block << "] : " << bs << " variants..." << flush;

      if(params.joint_test && !params.build_mask) allocate_mat(Gblock.Gmat, params.n_samples, bs);
      block_info.resize(bs);

      // compute single snp association test statistic
      get_sum_stats(chrom, bb, block_info);

      // update number of variants (if masks were built)
      bs = jt.setinfo[chrom - 1][bb].snp_indices.size();
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

          if( snp_data.ignored_trait(j) ) 
            continue;

          if(params.split_by_pheno)
            (*ofile_split[j]) << snp_data.sum_stats[j]; // add test info
          else
            ofile << snp_data.sum_stats[j]; // add test info
        }

      }

      if( params.joint_test ){

        // compute and print set-based test result
        t1 = std::chrono::high_resolution_clock::now();
        sout << "     -computing joint association tests..." << flush;

        jt.get_variant_names(chrom, bb, snpinfo);
        for(int j = 0; j < params.n_pheno; ++j) 
          if(params.split_by_pheno)
            (*ofile_split[j]) << jt.apply_joint_test(chrom, bb, j, &pheno_data, res.col(j), &Gblock, block_info, files.pheno_names[j], &params);
          else
            ofile << jt.apply_joint_test(chrom, bb, j, &pheno_data, res.col(j), &Gblock, block_info, files.pheno_names[j], &params); // add test info

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

}

// test SNPs in block
void Data::get_sum_stats(int const& chrom, int const& varset, vector<variant_block>& all_snps_info){

  auto t1 = std::chrono::high_resolution_clock::now();

  vector< vector < uchar > > snp_data_blocks;
  vector< uint32_t > insize, outsize;

  // read in markers and if applicable build masks
  if(!params.build_mask) readChunk(jt.setinfo[chrom - 1][varset].snp_indices, chrom, snp_data_blocks, insize, outsize, all_snps_info);
  else {
    getMask(chrom, varset, snp_data_blocks, insize, outsize, all_snps_info);
    // update size with new masks
    int n_snps = jt.setinfo[chrom - 1][varset].snp_indices.size();
    //cerr << "M=" << n_snps << endl;
    if(params.skip_test || (n_snps == 0)) return;

    // starting association testing with built masks
    t1 = std::chrono::high_resolution_clock::now();
    sout << "     -computing association tests..." << flush;
  }

  // analyze using openmp
  compute_tests_mt(chrom, jt.setinfo[chrom - 1][varset].snp_indices, snp_data_blocks, insize, outsize, all_snps_info);

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << "done (" << duration.count() << "ms) "<< endl;

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
    readChunkFromBGENFileToG(indices, chrom, snpinfo, &params, &Gblock, &in_filters, pheno_data.masked_indivs, pheno_data.phenotypes_raw, all_snps_info, sout);
  else if(params.file_type == "pgen") 
    readChunkFromPGENFileToG(indices, chrom, &params, &in_filters, &Gblock, pheno_data.masked_indivs, pheno_data.phenotypes_raw, snpinfo, all_snps_info);
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
  allocate_mat(Gblock.Gmat, params.n_samples, bsize);


  for(int i = 0; i < nchunks; i++){

    if( i == (nchunks-1) ) {
      bsize = n_snps - i * bsize;// use remainder number of variants
      allocate_mat(Gblock.Gmat, params.n_samples, bsize);
    }

    vector<uint64> indices (jt.setinfo[chrom - 1][varset].snp_indices.begin() + nvar_read, jt.setinfo[chrom - 1][varset].snp_indices.begin() + nvar_read + bsize);
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
