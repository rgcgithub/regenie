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

#include "cxxopts.hpp"
#include <regex>
#include <chrono>
#include <time.h>
#include "Regenie.hpp"
#include "Files.hpp"
#include "Geno.hpp"
#include "Joint_Tests.hpp"
#include "MultiTrait_Tests.hpp"
#include "Ordinal.hpp"
#include "survival_data.hpp"
#include "cox_score.hpp"
#include "Step1_Models.hpp"
#include "Step2_Models.hpp"
#include "Pheno.hpp"
#include "Masks.hpp"
#include "HLM.hpp"
#include "Data.hpp"

#include <boost/exception/all.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace std;
using namespace Eigen;
using namespace boost;


mstream::mstream(){ }
mstream::~mstream(){ }
MeasureTime::MeasureTime(){ }
MeasureTime::~MeasureTime(){ }


int main( int argc, char** argv ) {

  Data data;
  read_params_and_check(argc, argv, &data.params, &data.files, &data.in_filters, &data.runtime, data.sout);

  try {// after opening sout

    // for rng
    std::mt19937_64 rng_rd(data.params.rng_seed);
    data.params.rng_rd = &rng_rd;

    data.run();

  } catch (bad_alloc& badAlloc) {
    data.sout << "ERROR: bad_alloc caught, not enough memory (" << badAlloc.what() << ")\n";
    exit(EXIT_FAILURE);
  } catch (const std::string& msg){ 
    data.sout << "ERROR: " << msg << endl;
    exit(EXIT_FAILURE);
  } catch (const char* msg) {
    std::string str_msg = msg;
    data.sout <<  "ERROR: " <<  str_msg << endl;
    exit(EXIT_FAILURE);
  } catch (boost::exception const& e) {
    data.sout << "ERROR: " << boost::diagnostic_information(e) << endl;
    exit(EXIT_FAILURE);
  } catch (std::exception const&  e) {
    data.sout << "ERROR: " << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch (...) {
    data.sout << boost::current_exception_diagnostic_information() << endl;
    exit(EXIT_FAILURE);
}

  data.runtime.stop();

  data.sout << "\nElapsed time : " << std::chrono::duration<double>(data.runtime.end - data.runtime.begin).count() << "s" << endl;
  data.sout << "End time: " << ctime(&data.runtime.end_time_info) << endl; 

}


void print_header(std::ostream& o){

  std::ostringstream oss;
  string vnumber;
  // adjust spacing for version with Boost Iostream library (`.gz` suffix)
#ifndef HAS_BOOST_IOSTREAM
  oss << "  ";
#endif
  oss << "REGENIE v" << VERSION_NUMBER; 
  vnumber = oss.str();

  int out_width = 6;
  int total_width = vnumber.size() + out_width * 2;

  o << left << std::setw(14) << " " << "|" << std::string(total_width, '=')<< "|" << endl;
  o << left << std::setw(14) << " " << "|" << left << std::setw(out_width) << " " <<
    left << std::setw(total_width - out_width) << vnumber << "|" << endl;
  o << left << std::setw(14) << " " << "|" << std::string(total_width, '=')<< "|\n\n";

  o << "Copyright (c) 2020-2024 Joelle Mbatchou, Andrey Ziyatdinov and Jonathan Marchini." << endl;
  o << "Distributed under the MIT License.\n";
#ifdef HAS_BOOST_IOSTREAM
  o << "Compiled with Boost Iostream library.\n";
#endif

#if defined(WITH_HTSLIB)
  o << "Compiled with HTSlib.\n";
#endif

  // adding BLAS/LAPACK external routines
#if defined(WITH_MKL)
  o << "Using Intel MKL with Eigen.\n";
#elif defined(WITH_OPENBLAS)
  o << "Using BLAS/LAPACK routines from OpenBLAS with Eigen.\n";
#endif

  o << "\n";
}


void read_params_and_check(int& argc, char *argv[], struct param* params, struct in_files* files, struct filter* filters, MeasureTime* mt, mstream& sout) {

  cxxopts::Options AllOptions(argv[0], "");

  AllOptions.add_options()
    ("h,help", "print list of available options")
    ("helpFull", "print list of all available options")
    ;

  // add main options
  AllOptions.add_options("Main")
    ("step", "specify if fitting null model (=1) or association testing (=2)", cxxopts::value<int>(params->run_mode),"INT")
    ("bed", "prefix to PLINK .bed/.bim/.fam files", cxxopts::value<std::string>(files->bed_prefix),"PREFIX")
    ("pgen", "prefix to PLINK2 .pgen/.pvar/.psam files", cxxopts::value<std::string>(files->pgen_prefix),"PREFIX")
    ("bgen", "BGEN file", cxxopts::value<std::string>(files->bgen_file),"FILE")
    ("sample", "sample file corresponding to BGEN file", cxxopts::value<std::string>(files->sample_file),"FILE")
    ("bgi", "index bgi file corresponding to BGEN file", cxxopts::value<std::string>(files->bgi_file),"FILE")
    ("ref-first", "use the first allele as the reference for BGEN or PLINK bed/bim/fam input format [default assumes reference is last]")
    ("keep", "comma-separated list of files listing samples to retain in the analysis (no header; starts with FID IID)", cxxopts::value<std::string>(),"FILE")
    ("remove", "comma-separated list of files listing samples to remove from the analysis (no header; starts with FID IID)", cxxopts::value<std::string>(),"FILE")
    ("extract", "comma-separated list of files with IDs of variants to retain in the analysis", cxxopts::value<std::string>(),"FILE")
    ("exclude", "comma-separated list of files with IDs of variants to remove from the analysis", cxxopts::value<std::string>(),"FILE")
    ("p,phenoFile", "phenotype file (header required starting with FID IID)", cxxopts::value<std::string>(files->pheno_file),"FILE")
    ("phenoCol", "phenotype name in header (use for each phenotype to keep; can use parameter expansion {i:j})", cxxopts::value< std::vector<std::string> >(),"STRING")
    ("phenoColList", "comma separated list of phenotype names to keep (can use parameter expansion {i:j})", cxxopts::value<std::string>(),"STRING,..,STRING")
    ("eventColList", "comma separated list of event status names to keep (can use parameter expansion {i:j})", cxxopts::value<std::string>(),"STRING,..,STRING")
    ("c,covarFile", "covariate file (header required starting with FID IID)", cxxopts::value<std::string>(files->cov_file),"FILE")
    ("covarCol", "covariate name in header (use for each covariate to keep; can use parameter expansion {i:j})", cxxopts::value< std::vector<std::string> >(),"STRING")
    ("covarColList", "comma separated list of covariate names to keep (can use parameter expansion {i:j})", cxxopts::value<std::string>(),"STRING,..,STRING")
    ("catCovarList", "comma separated list of categorical covariates", cxxopts::value<std::string>(),"STRING,..,STRING")
    ("o,out", "prefix for output files", cxxopts::value<std::string>(files->out_file),"PREFIX")
    ("qt", "analyze phenotypes as quantitative")
    ("bt", "analyze phenotypes as binary")
    ("t2e", "analyze phenotypes as time to event")
    ("1,cc12", "use control=1,case=2,missing=NA encoding for binary traits")
    ("b,bsize", "size of genotype blocks", cxxopts::value<int>(params->block_size),"INT")
    ("cv", "number of cross validation (CV) folds", cxxopts::value<int>(params->cv_folds),"INT(=5)")
    ("loocv", "use leave-one out cross validation (LOOCV)")
    ("l0", "number of ridge parameters to use when fitting models within blocks [evenly spaced in (0,1)]", cxxopts::value<int>(params->n_ridge_l0),"INT(=5)")
    ("l1", "number of ridge parameters to use when fitting model across blocks [evenly spaced in (0,1)]", cxxopts::value<int>(params->n_ridge_l1),"INT(=5)")
    ("lowmem", "reduce memory usage by writing level 0 predictions to temporary files")
    ("lowmem-prefix", "prefix where to write the temporary files in step 1 (default is to use prefix from --out)", cxxopts::value<std::string>(files->loco_tmp_prefix),"PREFIX")
    ("split-l0", "split level 0 across N jobs and set prefix of output files", cxxopts::value<std::string>(),"PREFIX,N")
    ("run-l0", "run level 0 for job K in {1..N} using master file created from '--split-l0'", cxxopts::value<std::string>(),"FILE,K")
    ("run-l1", "run level 1 using master file from '--split-l0'", cxxopts::value<std::string>(files->split_file),"FILE")
    ("l1-phenoList", "run level 1 for a subset of the phenotypes (specified as comma-separated list)", cxxopts::value<std::string>(),"STRING,...,STRING")
    ("keep-l0", "avoid deleting the level 0 predictions written on disk after fitting the level 1 models")
    ("strict", "remove all samples with missingness at any of the traits")
    ("print-prs", "also output polygenic predictions without using LOCO (=whole genome PRS)")
    ("gz", "compress output files (gzip format)")
    ("apply-rint", "apply Rank-Inverse Normal Transformation to quantitative traits")
    ("apply-rerint", "apply Rank-Inverse Normal Transformation to residualized quantitative traits in step 2")
    ("apply-rerint-cov", "apply Rank-Inverse Normal Transformation to residualized quantitative traits and project covariates out in step 2")
    ("threads", "number of threads", cxxopts::value<int>(params->threads),"INT")
    ("pred", "file containing the list of predictions files from step 1", cxxopts::value<std::string>(files->blup_list_file),"FILE")
    ("ignore-pred", "skip reading predictions from step 1 (equivalent to linear/logistic regression with only covariates)")
    ("use-prs", "when using whole genome PRS step 1 output in '--pred'")
    ("write-samples", "write IDs of samples included for each trait (only in step 2)")
    ("minMAC", "minimum minor allele count (MAC) for tested variants", cxxopts::value<double>(params->min_MAC),"FLOAT(=5)")
    ("minINFO", "minimum imputation info score (Impute/Mach R^2) for tested variants", cxxopts::value<double>(params->min_INFO),"DOUBLE(=0)")
    ("no-split", "combine asssociation results into a single for all traits")
    ("firth", "use Firth correction for p-values less than threshold")
    ("approx", "use approximation to Firth correction for computational speedup")
    ("spa", "use Saddlepoint approximation (SPA) for p-values less than threshold")
    ("pThresh", "P-value threshold below which to apply Firth/SPA correction", cxxopts::value<double>(params->alpha_pvalue),"FLOAT(=0.05)")
    ("write-null-firth", "store coefficients from null models with approximate Firth for step 2")
    ("compute-all", "store Firth estimates for all chromosomes")
    ("use-null-firth", "use stored coefficients for null model in approximate Firth", cxxopts::value<std::string>(files->null_firth_file),"FILE")
    ("chr", "specify chromosome to test in step 2 (use for each chromosome)", cxxopts::value< std::vector<std::string> >(),"STRING")
    ("chrList", "Comma separated list of chromosomes to test in step 2", cxxopts::value<std::string>(),"STRING,..,STRING")
    ("range", "to specify a physical position window for variants to test in step 2", cxxopts::value<std::string>(),"CHR:MINPOS-MAXPOS")
    ("sex-specific", "for sex-specific analyses (male/female)", cxxopts::value<std::string>(),"STRING")
    ("af-cc", "print effect allele frequencies among cases/controls for step 2")
    ("test", "'additive', 'dominant' or 'recessive' (default is additive test)", cxxopts::value<std::string>(),"STRING")
    ("htp", "output association files in step 2 in HTP format specifying the cohort name)", cxxopts::value<std::string>(params->cohort_name),"STRING")
    ("condition-list", "file with list of variants to include as covariates", cxxopts::value<std::string>(files->condition_snps_list),"FILE")
    ("condition-file", "optional genotype file which contains the variants to include as covariates", cxxopts::value<std::string>(),"FORMAT,FILE")
    ("condition-file-sample", "sample file accompanying BGEN file with the conditional variants", cxxopts::value<std::string>(files->condition_snps_info.sample),"FILE")
    ("interaction", "perform interaction testing with a quantitative/categorical covariate", cxxopts::value<std::string>(filters->interaction_cov),"STRING")
    ("interaction-snp", "perform interaction testing with a variant", cxxopts::value<std::string>(filters->interaction_cov),"STRING")
    ("interaction-file", "optional genotype file which contains the variant for GxG interaction test", cxxopts::value<std::string>(),"FORMAT,FILE")
    ("interaction-file-sample", "sample file accompanying BGEN file with the interacting variant", cxxopts::value<std::string>(files->interaction_snp_info.sample),"FILE")
    ("interaction-file-reffirst", "use the first allele as the reference for the BGEN or PLINK file with the interacting variant [default assumes reference is last]")
    ("interaction-prs", "perform interaction testing with the full PRS from step 1")
    ("force-condtl", "to also condition on interacting SNP in the marginal GWAS test")
    ("no-condtl", "to print out all main effects in GxE interaction test")
    ("rare-mac", "minor allele count (MAC) threshold below which to use HLM for interaction testing with QTs", cxxopts::value<double>(params->rareMAC_inter),"FLOAT(=1000)")
    ("set-list", "file with sets definition", cxxopts::value<std::string>(files->set_file),"FILE")
    ("extract-sets", "comma-separated list of files with IDs of sets to retain in the analysis", cxxopts::value<std::string>(),"FILE")
    ("exclude-sets", "comma-separated list of files with IDs of sets to remove from the analysis", cxxopts::value<std::string>(),"FILE")
    ("extract-setlist", "comma separated list of sets to retain in the analysis", cxxopts::value<std::string>(),"STRING")
    ("exclude-setlist", "comma separated list of sets to remove from the analysis", cxxopts::value<std::string>(),"STRING")
    ("anno-file", "file with variant annotations", cxxopts::value<std::string>(files->anno_file),"FILE")
    ("anno-labels", "file with labels to annotations", cxxopts::value<std::string>(files->anno_labs_file),"FILE")
    ("mask-def", "file with mask definitions", cxxopts::value<std::string>(files->mask_file),"FILE")
    ("aaf-file", "file with AAF to use when building masks", cxxopts::value<std::string>(files->aaf_file),"FILE")
    ("set-singletons", "use 0/1 indicator in third column of AAF file to specify singleton variants")
    ("aaf-bins", "comma separated list of AAF bins cutoffs for building masks", cxxopts::value<std::string>(),"FLOAT,..,FLOAT")
    ("build-mask", "rule to construct masks, can be 'max', 'sum' or 'comphet' (default is max)", cxxopts::value<std::string>(params->mask_rule),"STRING")
    ("vc-tests", "comma separated list of tests to compute for each set of variants included in a mask [skat/skato/skato-acat/acatv/acato]", cxxopts::value<std::string>(),"STRING,..,STRING")
    ("vc-maxAAF", "maximum AAF for variants included in gene-based tests", cxxopts::value<double>(params->vc_maxAAF),"FLOAT(=1)")
    ("weights-col", "column index (1-based) for user-defined weights in annotation file", cxxopts::value<int>(params->vc_weight_col))
    ("joint", "comma spearated list of joint tests to perform", cxxopts::value<std::string>(params->burden),"STRING")
    ("singleton-carrier", "define singletons as variants with a single carrier in the sample")
    ("write-mask", "write masks in PLINK bed/bim/fam format")
    ("mask-lovo", "apply Leave-One-Variant-Out (LOVO) scheme when building masks (<set_name>,<mask_name>,<aaf_cutoff>)", cxxopts::value<std::string>(),"STRING")
    ("mask-lodo", "apply Leave-One-Domain-Out (LODO) scheme when building masks (<set_name>,<mask_name>,<aaf_cutoff>)", cxxopts::value<std::string>(),"STRING")
    ("skip-test", "skip computing association tests after building masks")
    ("check-burden-files", "check annotation file, set list file and mask file for consistency")
    ("strict-check-burden", "to exit early if the annotation, set list and mask definition files don't agree")
    ("force-qt", "force QT run for traits with few unique values")
    ("par-region", "build code to identify PAR region boundaries on chrX", cxxopts::value<std::string>(params->build_code),"STRING(=hg38)")
    ;


  // extended options
  AllOptions.add_options("Additional")
    ("v,verbose", "verbose screen output")
    ("version", "print version number and exit")
    ("minCaseCount", "minimum number of cases per trait", cxxopts::value<int>(params->mcc),"INT=10")
    ("tpheno-file", "transposed phenotype file (each row is a phenotype)", cxxopts::value<std::string>(files->pheno_file),"FILE")
    ("tpheno-indexCol", "index of column which contain phenotype name", cxxopts::value<uint32_t>(filters->tpheno_indexCol),"INT")
    ("tpheno-ignoreCols", "comma separated list of indexes for columns to ignore (can use parameter expansion {i:j})", cxxopts::value<std::string>(),"INT,...,INT")
    ("iid-only", "to specify if header in transposed phenotype file only contains sample IID")
    ("extract-or", "file with IDs of variants to retain in the analysis regardless of MAC", cxxopts::value<std::string>(),"FILE")
    ("exclude-or", "file with IDs of variants to remove from the analysis if MAC falls below threshold", cxxopts::value<std::string>(),"FILE")
    ("setl0", "comma separated list of ridge parameters to use when fitting models within blocks", cxxopts::value<std::string>(), "FLOAT,..,FLOAT")
    ("setl1", "comma separated list of ridge parameters to use when fitting model across blocks", cxxopts::value<std::string>(), "FLOAT,..,FLOAT")
    ("use-relative-path", "use relative paths for Step 1 pred.list file")
    ("phenoExcludeList", "comma separated list of phenotype names to ignore (can use parameter expansion {i:j})", cxxopts::value<std::string>(),"STRING,..,STRING")
    ("covarExcludeList", "comma separated list of covariates to ignore (can use parameter expansion {i:j})", cxxopts::value<std::string>(),"STRING,..,STRING")
    ("nauto", "number of autosomal chromosomes", cxxopts::value<int>(),"INT")
    ("exact-p", "output uncapped p-values in the summary statistic file with HTP format")
    ("skip-dosage-comp", "skip dosage compensation for males in chrX non-PAR regions")
    ("maxCatLevels", "maximum number of levels for categorical covariates", cxxopts::value<int>(params->max_cat_levels),"INT(=10)")
    ("max-condition-vars", "maximum number of variants to include as covariates", cxxopts::value<uint32_t>(params->max_condition_vars),"INT(=10000)")
    ("nb", "number of blocks to use", cxxopts::value<int>(params->n_block),"INT")
    ("starting-block", "start run at a specific block/set number for step 2", cxxopts::value<int>(params->start_block),"INT")
    ("force-step1", "run step 1 for more than 1M variants (not recommended)")
    ("write-mask-snplist", "file with list of variants that went into each mask")
    ("minHOMs", "minimum number of homozygote ALT carriers in recessive test", cxxopts::value<double>(params->minHOMs),"FLOAT(=0)")
    ("skat-params", "a1,a2 values for variant weights computed from Beta(MAF,a1,a2) used in gene-based tests", cxxopts::value<std::string>(),"FLOAT,FLOAT(=1,25)")
    ("skato-rho", "comma-separated list of rho values used for SKATO", cxxopts::value<std::string>(),"FLOAT,..,FLOAT")
    ("vc-MACthr", "MAC threshold below which to collapse variants for gene-based tests", cxxopts::value<int>(params->skat_collapse_MAC),"INT(=10)")
    ("lovo-snplist", "list of variants to generate LOVO masks for", cxxopts::value<std::string>(params->masks_loo_snpfile),"FILE")
    ("joint-only", "only output p-values from joint tests")
    ("force-ltco", "use a Leave-Two-Chromosome-Out (LTCO) scheme by specifying additional chromosome to exclude from step 1 LOCO predictions", cxxopts::value<int>(params->ltco_chr),"INT")
    ("niter", "maximum number of iterations for logistic regression", cxxopts::value<int>(params->niter_max),"INT(=50)")
    ("maxstep-null", "maximum step size in null Firth logistic regression", cxxopts::value<int>(params->maxstep_null),"INT(=25)")
    ("maxiter-null", "maximum number of iterations in null Firth logistic regression", cxxopts::value<int>(params->niter_max_firth_null),"INT(=1000)")
    ("skip-fast-firth", "skip fast implementation of approximate Firth for variants below MAC 50")
    ("force-impute", "keep and impute missing observations when in step 2 (default is to drop missing for each trait)")
    ("firth-se", "Compute SE for Firth based on effect size estimate and LRT p-value")
    ("print-pheno", "Print phenotype name when writing sample IDs to file (only for step 2)")
    ("compute-corr", "compute LD matrix (output R^2 values to binary file)")
    ("output-corr-text", "output matrix of Pearson correlations to text file")
    ("forcein-vars", "retain variants from extract file not present in genetic data file for the LD matrix")
    ("ld-extract", "file with list of variants & masks to compute LD matrix", cxxopts::value<string>(params->ld_list_file),"FILE")
    ("skip-scaleG", "compute LD matrix based on unscaled genotypes")
    ("sparse-thr", "threshold used to sparsify the LD matrix", cxxopts::value<double>(params->ld_sparse_thr),"FLOAT(=0)")
    ("print-vcov", "print variance-covariance matrix for interaction test to file")
    ;

  // extra options
  AllOptions.add_options("Extra")
    ("print", "print estimated effect sizes from level 0 and level 1 models")
    ("within", "use within-sample predictions as input when fitting model across blocks in step 1")
    ("early-exit", "Exit program after fitting level 0 models (avoid deleting temporary prediction files from level 0)")
    ("print-cov-betas", "Print covariate effects to file (assumes no multi-colinearity)")
    ("prior-alpha", "alpha value used when speifying the MAF-dependent prior on SNP effect sizes", cxxopts::value<double>(params->alpha_prior),"FLOAT(=-1)")
    ("prs-cov", "include step 1 predictions as covariate rather than offset")
    ("test-l0", "test association for each level 0 block")
    ("l0-pval-thr", "p-value threshold for identifying top SNPs at level 0", cxxopts::value<double>(params->l0_snp_pval_thr),"FLOAT")
    ("select-l0", "file with p-values for each level 0 block (use as flag if with --test-l0)", cxxopts::value<std::string>(params->l0_pvals_file)->implicit_value(""),"FILE")
    ("rm-l0-pct", "remove least x% significant blocks from level 1 models", cxxopts::value<double>(params->rm_l0_pct),"FLOAT(=0)")
    ("l0-event", "use event status as response in level 0 in time-to-event analysis")
    ("l1-full", "use all samples for final L1 model in Step 1 logistic ridge with LOOCV")
    ("prop-zero-thr", "min. proportion of zeros needed to sparsify the genotype vector", cxxopts::value<double>(params->prop_zero_thr),"FLOAT(=0.5)")
    ("force-robust", "use robust SE instead of HLM for rare variant GxE test with quantitative traits")
    ("force-hc4", "use HC4 instead of HC3 robust SE for rare variant GxE test with quantitative traits")
    ("no-robust", "don't use robust SEs or HLM for GxE test")
    ("write-setlist", "file with list of masks to combine as sets", cxxopts::value<std::string>(files->new_sets),"FILE")
    ("sbat-napprox", "number of random draws to use for approximate SBAT test", cxxopts::value<int>(params->nnls_napprox),"INT(=10)")
    ("sbat-adapt", "use adaptive strategy to compute p-value using fewer weights (k=2)")
    ("sbat-mtw", "re-use SBAT weights across all traits")
    ("sbat-verbose", "To output detailed SBAT test results")
    ("acat-beta", "parameters for Beta(a,b) used for weights in ACAT joint test", cxxopts::value<std::string>(), "a,b(=1,25)")
    ("hlm-novquad", "remove quadratic term for E in variance function of HLM model (only for GxE interaction test)")
    ("rgc-gene-p", "apply optimal strategy to extract single p-value per gene")
    ("rgc-gene-def", "file with list of mask groups to run single p-value strategy", cxxopts::value<std::string>(params->genep_mask_sets_file))
    ("skip-sbat", "skip running SBAT test for --rgc-gene-p")
    ("multiply-weights", "multiply the user defined weights by the default SKAT weights in SKAT/ACAT tests")
    ("skip-cf-burden", "skip computing per-mask calibration factor for SKAT tests")
    ("force-mac-filter", "apply a seperate MAC filter on a subset of the SNPs", cxxopts::value<std::string>(), "snpfile,MAC")
    ("use-adam", "use ADAM to fit penalized logistic models")
    ("adam-mini", "use mini-batch for ADAM")
    ("ct", "analyze phenotypes as counts")
    ("seed", "specify seed for random number generation", cxxopts::value<uint>(params->rng_seed))
    ("debug", "more verbose screen output for debugging purposes")
    ("mt", "run multi-trait tests")
    ("mcc", "apply MCC test for quantitative traits")
    ("mcc-skew", "absolute phenotypic skewness to activate MCC [default value is 0]", cxxopts::value<double>(params->mcc_skew),"FLOAT(=0)")
    ("mcc-thr", "threshold to apply MCC if activated [default value is 0.01]", cxxopts::value<double>(params->mcc_thr),"FLOAT(=0.01)")
    ("remeta-save-ld", "store SKAT matrices for use with remeta")
    ("remeta-ld-spr", "sparsity threshold for SKAT matrices", cxxopts::value<double>(params->remeta_ld_spr),"FLOAT(=0.01)")
    ("multiphen", "run MultiPhen test")
    ("multiphen-thr", "threshold to apply LRT for MultiPhen [default value is 0.01]", cxxopts::value<double>(params->multiphen_thr),"FLOAT(=0.001)")
    ("multiphen-test", "type of MultiPhen test", cxxopts::value<std::string>(params->multiphen_test),"STRING")
    ("multiphen-optim", "type of MultiPhen optimization algorithm", cxxopts::value<std::string>(params->multiphen_optim),"STRING")
    ("multiphen-tol", "toleance level for Firth [default value is 1e-4]", cxxopts::value<double>(params->multiphen_tol),"FLOAT(=0.0001)")
    ("multiphen-trace", "trace model fitting performance for MultiPhen")
    ("multiphen-firth-mult", "Firth penalty multiplier [default value is 1]", cxxopts::value<double>(params->multiphen_firth_mult),"FLOAT(=1.0)")
    ("multiphen-verbose", "MultiPhen verbose level", cxxopts::value<int>(params->multiphen_verbose),"INT(=0)")
    ("multiphen-maxstep", "Maximum step in IRLS for MultiPhen [default value is 100]", cxxopts::value<double>(params->multiphen_maxstep),"FLOAT(=25.0)")
    ("multiphen-approx-offset", "MAC to disable MultiPhen offset approximation", cxxopts::value<int>(params->multiphen_approx_offset),"INT(=-1)")
    ("multiphen-maxit", "MultiPhen maximum number of IRLS iterations", cxxopts::value<int>(params->multiphen_maxit),"INT(=150)")
    ("multiphen-maxit2", "MultiPhen maximum number of step-halving IRLS iterations", cxxopts::value<int>(params->multiphen_maxit2),"INT(=5)")
    ("multiphen-strict", "strict mode for MultiPhen IRLS")
    ("multiphen-pseudo-stophalf", "Threshold to stop step-halving in pseudo algorithm [default value is 0.0]", cxxopts::value<double>(params->multiphen_pseudo_stophalf),"FLOAT(=0.0)")
    ("multiphen-reset-start", "reset start values when failed convergence in MultiPhen")
    ("multiphen-offset", "offset mode for MultiPhen", cxxopts::value<std::string>(params->multiphen_offset),"STRING")
    ("t2e-event-l0", "Use event as reponse in level0 for time-to-event phenotype")
    ("t2e-l1-pi6", "use heritability to get penalty")
    ("coxnofirth", "not using firth in cox model, the test uses likelihood ratio test")
    ("coxscore-exact", "use exact score variance")
    ("nocov-approx", "skip adjusting for covariates in score test")
   ;

  try
  {
    bool acato_use_all_rhos = false;

    //AllOptions.parse_positional({"htp"});
    auto vm = AllOptions.parse(argc, argv);
    auto arguments = vm.arguments();
    map<string, bool> valid_args;
    for(const auto &kv: arguments)
      valid_args[ kv.key() ] = true;

    // help menu
    if (vm.count("help")){
      print_header(std::cout);
      std::cout << AllOptions.help({"", "Main"}) << '\n' << params->webinfo << "\n\n";
      exit(EXIT_SUCCESS);
    } else if (vm.count("helpFull")) {
      print_header(std::cout);
      std::cout << AllOptions.help({"", "Main", "Additional"}) << '\n' << params->webinfo << "\n\n";
      exit(EXIT_SUCCESS);
    } else if(vm.count("version")) {
      std::cout << "v" << VERSION_NUMBER << "\n";
      exit(EXIT_SUCCESS);
    }


    if( vm.unmatched().size() > 0 ) {
      std::cout << "\nERROR: There are unmatched arguments:\n";
      for(auto cn :  vm.unmatched())
        cout << "'" << cn << "' ";
      std::cout << "(Make sure there are no spaces in the options arguments)\n";
      exit(EXIT_FAILURE);
    }
    
    if (!vm.count("out")){
      print_header(std::cout);
      std::cout << "ERROR: You must provide an output prefix using '--out'" << '\n' << params->webinfo << "\n\n";
      exit(EXIT_FAILURE);
    }


    // Print output to file and to stdout
    // print command line arguments
    start_log(files->out_file, mt, sout);
    vector< string > tmp_str_vec;

    if( (vm.count("bgen") + vm.count("bed")  + vm.count("pgen"))  != 1 )
      throw "must use either --bed,--bgen or --pgen.";

    if( vm.count("bgen") ) params->file_type = "bgen";
    if( vm.count("bed") ) params->file_type = "bed";
    if( vm.count("pgen") ) params->file_type = "pgen";
    if( vm.count("sample") ) params->bgenSample = true;
    if( vm.count("ref-first") ) params->ref_first = true;
    if( vm.count("bt") ) params->trait_mode = 1;
    if( vm.count("ct") ) params->trait_mode = 2;
    if( vm.count("t2e") ) params->trait_mode = 3;
    if( vm.count("1") ) params->CC_ZeroOne = false;
    if( vm.count("loocv") ) params->use_loocv = true;
    if( vm.count("apply-rint") && !vm.count("bt")) params->rint = true;
    if( vm.count("apply-rerint") && !vm.count("bt")) params->rerint = true;
    if( vm.count("apply-rerint-cov") && !vm.count("bt")) params->rerintcov = true;
    if( vm.count("strict") ) params->strict_mode = true;
    if( vm.count("print-prs") ) params->print_prs = true;
    if( vm.count("use-relative-path") ) params->use_rel_path = true;
    if( vm.count("ignore-pred") ) params->skip_blups = true;
    if( vm.count("use-prs") ) params->use_prs = true;
    if( vm.count("prs-cov") ) params->blup_cov = true;
    if( vm.count("force-impute") ) params->rm_missing_qt = false;
    if( vm.count("no-split") ) params->split_by_pheno = false;
    if( vm.count("approx") ) params->firth_approx = true;
    if( vm.count("approx") && vm.count("skip-fast-firth") ) params->skip_fast_firth = true;
    if( vm.count("nauto") ) params->nChrom = vm["nauto"].as<int>() + 1;
    if( vm.count("maxstep-null") | vm.count("maxiter-null") ) params->fix_maxstep_null = true;
    if( vm.count("firth") ) params->firth = true;
    if( vm.count("write-null-firth") ) params->write_null_firth = true;
    if( vm.count("use-null-firth") ) params->use_null_firth = true;
    if( vm.count("compute-all") ) params->compute_all_chr = true;
    if( vm.count("spa") ) params->use_SPA = true;
    if( vm.count("minMAC") ) params->setMinMAC = true;
    if( vm.count("minINFO") ) params->setMinINFO = true;
    if( vm.count("htp") ) params->htp_out = params->split_by_pheno = true;
    if( vm.count("exact-p") ) params->uncapped_pvals = true;
    if( vm.count("multiphen") ) params->split_by_pheno = false;
    if( vm.count("af-cc") ) params->af_cc = true;
    if( vm.count("tpheno-file") ) params->transposedPheno = true;
    if( vm.count("v") ) params->verbose = true;
    if( vm.count("debug") ) params->verbose = params->debug = true;
    if( vm.count("range") ) params->set_range = true;
    if( vm.count("print") ) params->print_block_betas = true;
    if( vm.count("print-cov-betas") ) params->print_cov_betas = true;
    if( vm.count("test-l0") ) params->test_l0 = true;
    if( vm.count("l0-event") ) params->l0_event = true;
    if( vm.count("select-l0") ) params->select_l0 = true;
    //if( vm.count("nostream") ) params->streamBGEN = params->fastMode = false;
    //if( vm.count("within") ) params->within_sample_l0 = true;
    if( vm.count("write-samples") ) params->write_samples = true;
    if( vm.count("print-pheno") ) params->print_pheno_name = true;
    if( vm.count("early-exit") ) params->early_exit = true;
    if( vm.count("force-step1") ) params->force_run = true;
    if( (params->run_mode == 1) && vm.count("bt") && vm.count("loocv") && vm.count("l1-full") ) params->l1_full_samples = true;
    if( vm.count("lowmem") ) params->write_l0_pred = true;
    if( vm.count("keep-l0") ) params->rm_l0_pred = false;
    if( vm.count("split-l0") ) params->split_l0 = true;
    if( vm.count("run-l0") ) { params->run_l0_only = params->write_l0_pred = params->keep_snps = true;}
    if( vm.count("run-l1") ) params->run_l1_only = params->write_l0_pred = true;
    if( vm.count("firth") && vm.count("firth-se") ) params->back_correct_se = true;
    if( vm.count("use-adam") ) params->use_adam = true;
    if( vm.count("adam-mini") ) params->adam_mini = true;
    if( vm.count("niter") ) params->niter_max_ridge = params->niter_max;
    if( vm.count("force-ltco") ) params->w_ltco = true;
    if( vm.count("joint") ) params->joint_test = true;
    if( vm.count("joint-only") ) params->p_joint_only = true;
    if( vm.count("nocov-approx") ) params->skip_cov_res = true;
    if( vm.count("mt") ) params->trait_set = true;
    if( vm.count("mcc") ) params->mcc_test = true;
    if( vm.count("multiphen") ) params->multiphen = true;
    if( vm.count("multiphen-trace") ) params->multiphen_trace = true;
    if( vm.count("multiphen-strict") ) params->multiphen_strict = true;
    if( vm.count("multiphen-reset-start") ) params->multiphen_reset_start = true;
    if( vm.count("multiphen-strict") ) params->multiphen_strict = true;
    if( vm.count("aaf-file") ) params->set_aaf = true;
    if( vm.count("aaf-file") && vm.count("set-singletons") ) params->aaf_file_wSingletons = true;
    if( vm.count("singleton-carrier") ) params->singleton_carriers = true;
    if( vm.count("mask-lovo") ) params->mask_loo = true;
    if( vm.count("mask-lodo") ) params->mask_lodo = true;
    if( vm.count("write-mask") ) params->write_masks = true;
    if( vm.count("write-setlist") ) params->write_setlist = true;
    if( vm.count("write-mask-snplist") ) params->write_mask_snplist = true;
    if( vm.count("skip-test") ) params->skip_test = true;
    if( vm.count("check-burden-files") ) params->check_mask_files = true;
    if( vm.count("strict-check-burden") ) params->strict_check_burden = true;
    if( vm.count("force-qt") ) params->force_qt_run = true;
    if( vm.count("weights-col") ) params->vc_with_weights = true;
    if( vm.count("multiply-weights") ) params->vc_multiply_weights = true;
    if( vm.count("skip-dosage-comp") ) params->skip_dosage_comp = true;
    if( vm.count("sbat-verbose") ) params->nnls_out_all = true;
    if( vm.count("sbat-adapt") ) params->nnls_adaptive = true;
    if( vm.count("sbat-mtw") ) params->nnls_mt_weights = true;
    if( vm.count("skip-cf-burden") ) params->skip_cf_burden = true;
    if( vm.count("condition-list") ) { params->condition_snps = true;params->rm_snps = true;}
    if( vm.count("force-robust") ) params->force_robust = true;
    if( vm.count("force-hc4") ) params->force_robust = params->force_hc4 = true;
    if( vm.count("no-robust") ) params->no_robust = true;
    if( vm.count("hlm-novquad") ) params->hlm_vquad = false;
    if( vm.count("print-vcov") ) params->print_vcov = true;
    if( vm.count("compute-corr") || vm.count("output-corr-text") ) {
      params->getCorMat = true;
      params->cormat_force_vars = (vm.count("forcein-vars") && vm.count("extract")) || vm.count("ld-extract");
      params->skip_scaleG = vm.count("skip-scaleG");
      params->run_mode = 2;
      params->skip_blups = params->strict_mode = true;
      params->trait_mode = 0;
      params->min_MAC = 0.5;
      if(vm.count("output-corr-text") || vm.count("skip-scaleG")) params->cor_out_txt = true;
      if(vm.count("exclude")) throw "cannot use --exclude with --compute-corr (use --extract instead)";
      if(vm.count("write-mask")){
        sout << "WARNING: option --write-mask cannot be used when computing LD.\n" ;
        params->write_masks = false; valid_args[ "write-mask" ] = false;
      }
    }
    if( vm.count("gz") ) {
# if defined(HAS_BOOST_IOSTREAM)
      // only works when compiled with boost IO library
      params->gzOut = true;
# else
      sout << "WARNING: REGENIE was not compiled with Boost Iostream library so ignoring option '--gz'.\n";
      valid_args[ "gz" ] = false;
#endif
    }


    if( vm.count("phenoColList") ) {
      params->select_phenos = true;
      tmp_str_vec = string_split(vm["phenoColList"].as<string>(),",");
      for( size_t i = 0; i < tmp_str_vec.size(); i++) {
        for(auto cn : check_name(tmp_str_vec[i], sout))
          filters->pheno_colKeep_names[cn] = true;
      }
    }
    if( vm.count("phenoCol") ) {
      params->select_phenos = true;
      tmp_str_vec = vm["phenoCol"].as<std::vector<string>>();
      for( size_t i = 0; i < tmp_str_vec.size(); i++)
        for(auto cn : check_name(tmp_str_vec[i], sout))
          filters->pheno_colKeep_names[cn] = true;
    }
    if( vm.count("phenoExcludeList") ) {
      params->select_phenos_rm = true;
      tmp_str_vec = string_split(vm["phenoExcludeList"].as<string>(),",");
      for( size_t i = 0; i < tmp_str_vec.size(); i++) {
        for(auto cn : check_name(tmp_str_vec[i], sout))
          filters->pheno_colRm_names[cn] = true;
      }
    }
    if( (params->trait_mode == 3) && vm.count("eventColList") && vm.count("phenoCol") )
      throw "You must specify TTE phenotypes using '--phenoColList' (matching in order with events in '--eventColList').";
    if( (params->trait_mode == 3) && (!vm.count("eventColList") || !vm.count("phenoColList")) ) 
      throw "You must specify both '--phenoColList' and '--eventColList' (same order) for time-to-event analysis.";
    if( vm.count("eventColList") ) { // time-to-event names map
      if( params->trait_mode != 3) 
        throw "Option --eventColList must be used with '--t2e' for time-to-event analysis";
      params->select_phenos = true;
      tmp_str_vec = string_split(vm["eventColList"].as<string>(),",");
      vector< string > tmp_str_vec_time = string_split(vm["phenoColList"].as<string>(),",");
      for( size_t i = 0; i < tmp_str_vec.size(); i++) {
        files->t2e_map[tmp_str_vec_time[i]] = tmp_str_vec[i];
        for(auto cn : check_name(tmp_str_vec[i], sout))
          filters->pheno_colKeep_names[cn] = true;
      }
      params->t2e_event_l0 = vm.count("t2e-event-l0");
      params->t2e_l1_pi6 = vm.count("t2e-l1-pi6");
      params->cox_nofirth = vm.count("coxnofirth");
      params->coxscore_exact = vm.count("coxscore-exact");
    }
    if( vm.count("covarColList") ) {
      params->select_covs = true;
      tmp_str_vec = string_split(vm["covarColList"].as<string>(),",");
      for( size_t i = 0; i < tmp_str_vec.size(); i++){
        for(auto cn : check_name(tmp_str_vec[i], sout))
          filters->cov_colKeep_names[cn] = true;
      }
    }
    if( vm.count("covarCol") ) {
      params->select_covs = true;
      tmp_str_vec = vm["covarCol"].as<std::vector<string>>();
      for( size_t i = 0; i < tmp_str_vec.size(); i++)
        for(auto cn : check_name(tmp_str_vec[i], sout))
          filters->cov_colKeep_names[cn] = true;
    }
    if( vm.count("covarExcludeList") ) {
      params->select_covs_rm = true;
      tmp_str_vec = string_split(vm["covarExcludeList"].as<string>(),",");
      for( size_t i = 0; i < tmp_str_vec.size(); i++) {
        for(auto cn : check_name(tmp_str_vec[i], sout))
          filters->cov_colRm_names[cn] = true;
      }
    }
    if( vm.count("catCovarList") ) {
      tmp_str_vec = string_split(vm["catCovarList"].as<string>(),",");
      for( size_t i = 0; i < tmp_str_vec.size(); i++)
        for(auto cn : check_name(tmp_str_vec[i], sout))
          filters->cov_colKeep_names[cn] = false;
    }
    if( (params->run_mode ==2) && (vm.count("interaction") || vm.count("interaction-snp")) ) {
      params->w_interaction = true;
      if(vm.count("interaction-snp")) params->interaction_snp = params->w_ltco =  true;
      check_inter_var(filters->interaction_cov, filters->interaction_cov_null_level, sout);
      if(!vm.count("interaction-snp") && !in_map(filters->interaction_cov,filters->cov_colKeep_names))
        filters->cov_colKeep_names[filters->interaction_cov] = true; // assume qt
      if(vm.count("no-condtl") || (vm.count("interaction-snp") && !vm.count("force-condtl")) )
        params->gwas_condtl = false;
    }
    if( (params->run_mode ==2) && vm.count("interaction-prs") ) {
      params->w_interaction = true;
      params->interaction_prs = true;
      filters->interaction_cov = "PRS";
      if(vm.count("no-condtl") || (!vm.count("force-condtl")) )
        params->gwas_condtl = false;
    }
    if( vm.count("tpheno-ignoreCols") ) {
      tmp_str_vec = string_split(vm["tpheno-ignoreCols"].as<string>(),",");
      for( size_t i = 0; i < tmp_str_vec.size(); i++) {
        for(auto cn : check_name(tmp_str_vec[i], sout))
          filters->tpheno_colrm[ stoi(cn) ] = true;
      }
    }
    if( vm.count("chrList") ) {
      params->select_chrs = true;
      tmp_str_vec = string_split(vm["chrList"].as<string>(),",");
      for( size_t ichr = 0; ichr < tmp_str_vec.size(); ichr++)
        for(auto cn : check_name(tmp_str_vec[ichr], sout))
          filters->chrKeep_test[ chrStrToInt(cn, params->nChrom) ] = true;
    }
    if( vm.count("chr") ) {
      params->select_chrs = true;
      tmp_str_vec = vm["chr"].as<std::vector<string>>();
      for( size_t ichr = 0; ichr < tmp_str_vec.size(); ichr++)
        filters->chrKeep_test[ chrStrToInt(tmp_str_vec[ichr], params->nChrom) ] = true;
    }
    if( vm.count("keep") ){
      files->file_ind_include = string_split(vm["keep"].as<string>(),",");
      params->keep_indivs = true;
    }
    if( vm.count("remove") ){
      files->file_ind_exclude = string_split(vm["remove"].as<string>(),",");
      params->rm_indivs = true;
    }
    if( vm.count("extract") ){
      files->file_snps_include = string_split(vm["extract"].as<string>(),",");
      params->keep_snps = true;
    }
    if( !vm.count("run-l0") && vm.count("exclude") ){
      files->file_snps_exclude = string_split(vm["exclude"].as<string>(),",");
      params->rm_snps = true;
    }
    if( vm.count("extract-or") ){
      files->file_snps_include_or = string_split(vm["extract-or"].as<string>(),",");
      params->keep_or = true;
    }
    if( vm.count("exclude-or") ){
      files->file_snps_exclude_or = string_split(vm["exclude-or"].as<string>(),",");
      params->rm_or = true;
    }
    if( vm.count("extract-sets") ){
      files->file_sets_include = string_split(vm["extract-sets"].as<string>(),",");
      params->keep_sets = true;
    }
    if( vm.count("exclude-sets") ){
      files->file_sets_exclude = string_split(vm["exclude-sets"].as<string>(),",");
      params->rm_sets = true;
    }
    if( vm.count("extract-setlist") ) {
      params->set_select_list = params->keep_sets = true;
      files->file_sets_include.resize(1);
      files->file_sets_include[0] = vm["extract-setlist"].as<string>();
    }
    if( vm.count("exclude-setlist") ) {
      params->set_select_list = params->rm_sets = true;
      files->file_sets_exclude.resize(1);
      files->file_sets_exclude[0] = vm["exclude-setlist"].as<string>();
    }
    if( vm.count("split-l0") ) { // Format: FILE,INT
      tmp_str_vec = string_split(vm["split-l0"].as<string>(),",");
      if(tmp_str_vec.size() != 2 )
        throw "wrong format for --split-l0 (must be FILE,INT).";
      files->split_file = tmp_str_vec[0];
      params->njobs = atoi( tmp_str_vec[1].c_str() );
    }
    if( vm.count("run-l0") ) { // Format: FILE,INT
      tmp_str_vec = string_split(vm["run-l0"].as<string>(),",");
      if(tmp_str_vec.size() != 2 )
        throw "wrong format for --run-l0 (must be FILE,INT).";
      files->split_file = tmp_str_vec[0];
      params->job_num = atoi( tmp_str_vec[1].c_str() );
      if(params->job_num < 1 )
        throw "invalid job number for --run-l0 (must be >=1).";
    }
    if( vm.count("condition-file") ) {
      tmp_str_vec = string_split(vm["condition-file"].as<string>(),",");
      if(tmp_str_vec.size()<2)
        throw "invalid option input for --condition-file";
      if((tmp_str_vec[0] != "bgen") && (tmp_str_vec[0] != "bed") && (tmp_str_vec[0] != "pgen"))
        throw "invalid file format for --condition-file (either bed/bge/pgen)";
      files->condition_snps_info.format = tmp_str_vec[0];
      files->condition_snps_info.file = tmp_str_vec[1];
      params->condition_file = true;
    }
    if( vm.count("interaction-file") ) {
      tmp_str_vec = string_split(vm["interaction-file"].as<string>(),",");
      if(tmp_str_vec.size()<2)
        throw "invalid option input for --interaction-file";
      if((tmp_str_vec[0] != "bgen") && (tmp_str_vec[0] != "bed") && (tmp_str_vec[0] != "pgen"))
        throw "invalid file format for --interaction-file (either bed/bge/pgen)";
      files->interaction_snp_info.format = tmp_str_vec[0];
      files->interaction_snp_info.file = tmp_str_vec[1];
      files->interaction_snp_info.ref_first = vm.count("interaction-file-reffirst") && (tmp_str_vec[0] != "pgen");
      params->interaction_file = true;
    }
    if( vm.count("test") ) {
      if( vm["test"].as<string>() == "additive") params->test_type = 0; 
      else if( vm["test"].as<string>() == "dominant") params->test_type = 1; 
      else if( vm["test"].as<string>() == "recessive") params->test_type = 2; 
      else throw "unrecognized argument for option --test, must be either 'additive', 'dominant' or 'recessive'.";
    }
    if( vm.count("range") ) { // Format: Chr:min-max
      char tmp_chr[20];
      double p0 = -1, p1 = -1;
      string tmpd = vm["range"].as<string>();

      if(sscanf( tmpd.c_str(), "%[^:]:%lf-%lf", tmp_chr, &p0, &p1 ) != 3
          || (p0 < 0) || (p1 < 0) ) 
        //cerr << tmp_chr << "\t" << p0 << "\t" << p1 << endl;
        throw "wrong format for --range (must be CHR:MINPOS-MAXPOS).";

      tmpd = tmp_chr;
      params->range_chr = chrStrToInt(tmpd, params->nChrom);
      params->range_min = min(p0,p1);
      params->range_max = max(p0,p1);
    }
    if(vm.count("sex-specific")){
      if(vm["sex-specific"].as<string>() == "male")
        params->sex_specific = 1;
      else if(vm["sex-specific"].as<string>() == "female")
        params->sex_specific = 2;
      else throw "unrecognized argument for --sex-specific (should be either male/female)";
    }

    if( vm.count("build-mask") ) {
      if( params->mask_rule == "max") params->mask_rule_max = true; 
      else if( params->mask_rule == "sum") params->mask_rule_max = false; 
      else if( params->mask_rule == "comphet") { 
        params->mask_rule_max = false, params->mask_rule_comphet = true; 
      } else throw "unrecognized argument for option --build-mask (=" + params->mask_rule + ").";
      if((params->mask_rule == "sum") && params->htp_out){
        sout << "WARNING: option --htp cannot be used with '--build-mask sum' and will be ignored.\n";
        params->htp_out = false; valid_args[ "htp" ] = false;
      }
    }
    if( vm.count("acat-beta") ) {
      tmp_str_vec = string_split(vm["acat-beta"].as<string>(),",");
      params->acat_a1 = convertDouble( tmp_str_vec[0], params, sout);
      params->acat_a2 = convertDouble( tmp_str_vec[1], params, sout);
    }
    if( vm.count("vc-tests") ) {
      tmp_str_vec = string_split(vm["vc-tests"].as<string>(),",");
      for( size_t i = 0; i < tmp_str_vec.size(); i++)
        if(in_map(tmp_str_vec[i], params->vc_tests_map)) BIT_SET(params->vc_test, params->vc_tests_map[tmp_str_vec[i]]);
        else if(tmp_str_vec[i] == "acato-full") {acato_use_all_rhos = true; BIT_SET(params->vc_test, params->vc_tests_map["acato"]);}
        else throw "unrecognized VC test: '" + tmp_str_vec[i] + "' (accepted=skat/skato/skato-acat/acatv/acato)";
    }
    if( vm.count("rgc-gene-p") && vm.count("anno-file") && vm.count("mask-def") ) {
      params->apply_gene_pval_strategy = params->joint_test = true;
      if(!vm.count("vc-maxAAF")) params->vc_maxAAF = 0.01;
      if(params->burden != "") params->burden.append(",");
      params->burden.append("acat");
      if(!params->trait_mode && !vm.count("skip-sbat")) params->burden.append(",sbat");
      if(params->test_type == 0){
        BIT_SET(params->vc_test, params->vc_tests_map["acatv"]);
        BIT_SET(params->vc_test, params->vc_tests_map["skato-acat"]);
      } else {
        sout << "WARNING: SKATO/ACATV will be skipped for non-additive tests.\n";
        params->vc_test = 0;
      }
      if(vm.count("rgc-gene-def")) check_file (params->genep_mask_sets_file, "rgc-gene-def");
    } else if(vm.count("rgc-gene-p") || vm.count("rgc-gene-def")) {
      valid_args[ "rgc-gene-p" ] = false; // option is ignored
      valid_args[ "rgc-gene-def" ] = false; // option is ignored
    }

    if( CHECK_BIT(params->vc_test, params->vc_tests_map["acato"]) ) {// acato
      BIT_SET(params->vc_test, params->vc_tests_map["acatv"]); // acatv
      params->skato_rho.resize(2,1); params->skato_rho << 0, 1; // skat & burden
    }
    if( acato_use_all_rhos || ( (params->vc_test>>2)&3 ) ) {// skato/skato-acat or acato with all rhos
      params->skato_rho.resize(8,1); 
      params->skato_rho << 0, 0.1*0.1, 0.2*0.2, 0.3*0.3, 0.4*0.4, 0.5*0.5, 0.5, 1;
    }
    if( params->vc_test && vm.count("skat-params") ) {
      tmp_str_vec = string_split(vm["skat-params"].as<string>(),",");
      params->skat_a1 = convertDouble( tmp_str_vec[0], params, sout);
      params->skat_a2 = convertDouble( tmp_str_vec[1], params, sout);
    }
    if( ((params->vc_test>>1)&15) && vm.count("skato-rho") ) {
      if(acato_use_all_rhos)
        sout << "WARNING: ACATO will use the user-specified rho values for SKATO models.\n" ;
      tmp_str_vec = string_split(vm["skato-rho"].as<string>(),",");
      params->skato_rho = get_unit_params(true, "--skato-rho", tmp_str_vec, params, sout);
      if(params->skato_rho.size() > 1) BIT_SET(params->vc_test, 2);
    }

    if ( params->run_mode == 1 ) params->test_mode = false;
    else if (params->run_mode == 2 ) params->test_mode = true;
    else throw "specify which mode regenie should be running using option --step.";

    if(!params->test_mode) {

      // loocv only used with out-of-sample predictions
      if(params->use_loocv && params->within_sample_l0) {
        sout << "WARNING: option --loocv cannot be used with option --within.\n" ;
        params->use_loocv = false; valid_args[ "loocv" ] = false;
      }

      // writing of level 0 predictions only available when using out-of-sample predictions
      if(params->write_l0_pred && params->within_sample_l0){
        sout << "WARNING: option --lowmem cannot be used with option --within.\n" ;
        params->write_l0_pred = false; valid_args[ "lowmem" ] = valid_args[ "lowmem-prefix" ] = false;
      }

      // user specified ridge parameters to use at l0
      if( vm.count("setl0") ) {
        params->user_ridge_params_l0 = true;
        tmp_str_vec = string_split(vm["setl0"].as<string>(),",");
        params->lambda = get_unit_params(false, "--l0", tmp_str_vec, params, sout);
        params->n_ridge_l0 = params->lambda.size();
      } else set_ridge_params(params->n_ridge_l0, params->lambda, sout);

      // user specified ridge parameters to use at l1
      params->tau.resize(1); // may be assigned for each trait
      if( vm.count("setl1") ) {
        params->user_ridge_params_l1 = true;
        tmp_str_vec = string_split(vm["setl1"].as<string>(),",");
        params->tau[0] = get_unit_params(false, "--l1", tmp_str_vec, params, sout);
        params->n_ridge_l1 = params->tau[0].size();
      } else set_ridge_params(params->n_ridge_l1, params->tau[0], sout);

      if( params->run_l1_only && vm.count("l1-phenoList") ) {
        tmp_str_vec = string_split(vm["l1-phenoList"].as<string>(),",");
        for( size_t i = 0; i < tmp_str_vec.size(); i++) {
          for(auto cn : check_name(tmp_str_vec[i], sout))
            params->select_pheno_l1[cn] = true;
        }
      }

      // firth only done in test mode
      if(params->firth) params->firth = false;
      if(params->use_SPA) params->use_SPA = false;
      valid_args[ "firth" ] = valid_args[ "spa" ] = valid_args[ "approx" ] = false;

      params->test_type = 0;
      if( vm.count("range") ) {
        params->set_range = false; valid_args[ "range" ] =false;
        sout << "WARNING: option --range only works for step 2.\n";
      }
      if(params->rm_or || params->keep_or){
        sout << "WARNING: Options --extract-or/--exclude-or only work in step 2.\n";
        params->rm_or = params->keep_or = false; valid_args[ "extract-or" ] = valid_args[ "exclude-or" ] = false;
      }

    } 
    if(params->firth && (params->trait_mode!=1 && params->trait_mode!=3)) {
      // firth correction is only applied to binary traits and time-to-event traits
      sout << "WARNING: option --firth will not be applied (it is only run with binary traits and time-to-event traits).\n";
      params->firth = false; valid_args[ "firth" ] = valid_args[ "approx" ] = false;
    } 
    if(params->use_SPA && (params->trait_mode!=1)) {
      // SPA is only applied to binary traits
      sout << "WARNING: option --spa will not be applied (it is only run with binary traits).\n";
      params->use_SPA = false; valid_args[ "spa" ] = false;
    }

    if(vm.count("covarExcludeList") && !vm.count("covarFile")) {
      params->select_covs_rm = false; valid_args[ "covarExcludeList" ] = false;
    }

    if(params->test_mode && params->use_loocv) {params->use_loocv = false;valid_args[ "loocv" ] = false;}

    if( (vm.count("write-samples") || vm.count("write-mask")) && vm.count("bgen") && !vm.count("sample") )
      throw "must specify sample file (using --sample) if writing sample IDs to file.";

    if( vm.count("test") && (params->run_mode !=2)) 
      throw "can only use --test in step 2 (association testing).";
    if( (params->test_type > 0) && params->vc_test) 
      throw "cannot use --test with --vc-tests.";
    if(params->skip_dosage_comp && params->test_type)
      throw "cannot use --skip-dosage-comp with --test.";
    if( !params->getCorMat && params->joint_test ){
      if( (params->test_type > 0) && !vm.count("rgc-gene-p")) 
        throw "cannot use --test with --joint.";
      else if ( vm.count("sbat-napprox") && params->nnls_napprox < 1 )
        throw "must pass positive integer for --sbat-napprox.";
      params->snp_set = true;
    }
    if(vm.count("sparse-thr")){
     if(!vm.count("skip-scaleG") )
      throw "cannot use --sparse-thr without --skip-scaleG";
     else if((params->ld_sparse_thr < 0) || (params->ld_sparse_thr >=1))
      throw "invalid value passed in --sparse-thr (must be in [0,1)";
    }
    if(vm.count("ld-extract") && !vm.count("compute-corr"))
      throw "must use --ld-extract with --compute-corr";
    if(vm.count("ld-extract") && (vm.count("extract-sets")+vm.count("exclude-sets")+vm.count("extract-setlist")+vm.count("exclude-setlist")))
      throw "cannot use --ld-extract with --extract-sets/--exclude-sets";
    if( vm.count("write-null-firth") && vm.count("use-prs") )
      throw "cannot use --write-null-firth with --use-prs";

    if( vm.count("anno-file") || vm.count("mask-def") ){

      if( params->getCorMat && !vm.count("ld-extract") )
        throw "must use --ld-extract if building masks in LD matrix.";

      if(vm.count("anno-labels")) params->w_anno_lab = true;

      if( !(vm.count("anno-file") && vm.count("mask-def")) )
        throw "must use --anno-file with --mask-def.";

      if( (params->test_type > 0) && !(params->mask_rule_max || params->mask_rule_comphet) )
        throw "only additive test allowed when using 'sum' in --build-mask.";

      if(params->write_masks && !params->mask_rule_max && !params->mask_rule_comphet )
        throw "cannot write masks when using 'sum' in --build-mask.";

      // store aaf bins if given
      if( vm.count("aaf-bins") ) 
        tmp_str_vec = string_split(vm["aaf-bins"].as<string>(),",");
      else if(vm.count("rgc-gene-p")) tmp_str_vec = std::vector<std::string>({ "0.00001","0.0001","0.001","0.01" }); 
      else tmp_str_vec.resize(0);
      params->mbins = tmp_str_vec;

      if( vm.count("mask-lovo") ) {
        int cstart = 1;
        tmp_str_vec = string_split(vm["mask-lovo"].as<string>(),",");
        if(tmp_str_vec.size() < 3)
          throw "wrong format for option --mask-lovo.";
        else if ( tmp_str_vec.size() == 4 ) {
          params->w_regions = true; cstart++;
        }
        params->mask_loo_set = tmp_str_vec[0];
        if(params->w_regions) params->mask_loo_region = tmp_str_vec[cstart-1];
        params->mask_loo_name = tmp_str_vec[cstart];
        params->mbins.resize(1);
        params->mbins[0] = tmp_str_vec[cstart+1]; // either singleton or AAF cutoff
        if(params->vc_test){
          if(params->mbins[0] == "all") params->vc_maxAAF = 1;
          else if(params->mbins[0] != "singleton") params->vc_maxAAF = convertDouble( params->mbins[0], params, sout);
        }
        if(params->write_masks){
          sout << "WARNING: cannot use --write-mask with --mask-lovo.\n";
          params->write_masks = false; valid_args[ "write-mask" ] = false;
        }
        if(params->joint_test)
          throw "cannot use --joint with --mask-lovo";
        valid_args[ "vc-maxAAF" ] = valid_args[ "aaf-bins" ] = false;
      } else if (vm.count("lovo-snplist"))
        throw "cannot use --lovo-snplist without --mask-lovo";

      if( vm.count("mask-lodo") ) {
        tmp_str_vec = string_split(vm["mask-lodo"].as<string>(),",");
        if(tmp_str_vec.size() != 3)
          throw "wrong format for option --mask-lodo.";
        else if(vm.count("mask-lovo"))
          throw "cannot use --mask-lovo with --mask-lodo.";
        params->w_regions = true;
        params->mask_loo_set = tmp_str_vec[0];
        params->mask_loo_name = tmp_str_vec[1];
        params->mbins.resize(1);
        params->mbins[0] = tmp_str_vec[2]; // either singleton or AAF cutoff
        if(params->vc_test){
          if(params->mbins[0] == "all") params->vc_maxAAF = 1;
          else if(params->mbins[0] != "singleton") params->vc_maxAAF = convertDouble( params->mbins[0], params, sout);
        }
        if(params->write_masks){
          sout << "WARNING: cannot use --write-mask with --mask-lodo.\n";
          params->write_masks = false; valid_args[ "write-mask" ] = false;
        }
        valid_args[ "vc-maxAAF" ] = valid_args[ "aaf-bins" ] = false;
      }

      params->snp_set = true;
      params->build_mask = true;

    }
    if( params->test_mode && vm.count("force-mac-filter") ) {
      tmp_str_vec = string_split(vm["force-mac-filter"].as<string>(),",");
      params->forced_MAC_snpfile = tmp_str_vec[0];
      params->forced_MAC = convertDouble( tmp_str_vec[1], params, sout);
      if(params->forced_MAC < 0.5) throw "MAC must be greater than 0.5 for --force-mac-filter";
      if(params->rm_or || params->keep_or) throw "option --force-mac-filter cannot be used with --extract-or/--exclude-or";
      if(params->build_mask) throw "option --force-mac-filter cannot be used when building masks";
    } else valid_args[ "force-mac-filter" ] = false;

    if(!params->build_mask && params->write_masks) {params->write_masks = false; valid_args[ "write-mask" ] = false;}
    if(!params->build_mask && params->check_mask_files) {params->check_mask_files = false; valid_args[ "check-burden-files" ] = false;}
    if(!params->build_mask && params->strict_check_burden) {params->strict_check_burden = false; valid_args[ "strict-check-burden" ] = false;}
    if(!params->build_mask && params->write_mask_snplist) {params->write_mask_snplist = false; valid_args[ "write-mask-snplist" ] = false;}
    if(!(params->write_masks || params->write_mask_snplist) && params->skip_test) {params->skip_test = false; valid_args[ "skip-test" ] = false;}
    if(!params->w_interaction) params->gwas_condtl = false;
    if(!params->write_masks && params->write_setlist) {
      sout << "WARNING: must use --write-setlist with --write-mask.\n";
      params->write_setlist = false; valid_args[ "write-setlist" ] = false;
    }
    if((vm.count("1") || vm.count("cc12")) && (params->trait_mode != 1 || params->trait_mode != 3)) valid_args[ "1" ] = valid_args[ "cc12" ] = false;
    if( vm.count("write-mask-snplist") && (vm.count("mask-lovo") || vm.count("mask-lodo")) ) {
      sout << "WARNING: cannot use --write-mask-snplist with LOVO/LODO.\n";
      params->write_mask_snplist = false; valid_args[ "write-mask-snplist" ] = false;
    }
    if( vm.count("write-setlist") && (vm.count("mask-lovo") || vm.count("mask-lodo")) ) {
      sout << "WARNING: cannot use --write-setlist with LOVO/LODO.\n";
      params->write_setlist = false; valid_args[ "write-setlist" ] = false;
    }

    if( params->snp_set && !vm.count("set-list") )
      throw "must specify set list (using --set-list).";

    if( params->snp_set && 
        (vm.count("extract-sets")+vm.count("exclude-sets")+vm.count("extract-setlist")+vm.count("exclude-setlist"))>1 
      )
      throw "must use only one of --extract-sets/--exclude-sets/--extract-setlist/--exclude-setlist.";
    if( params->w_interaction && params->vc_test ){
      sout << "WARNING: skipping non-burden gene-based tests for GxG/GxE mode.\n";
      params->vc_test = 0; params->apply_gene_pval_strategy = params->joint_test = false;
      valid_args[ "vc-tests" ] = valid_args[ "joint" ] = valid_args[ "sbat-adapt" ] = valid_args[ "rgc-gene-p" ] = valid_args[ "rgc-gene-def" ] = valid_args[ "vc-maxAAF" ] = valid_args[ "vc-MACthr" ] = false;
    }

    if(!params->test_mode && params->setMinMAC){
      sout << "WARNING: option --minMAC only works in step 2 of REGENIE.\n";
      params->setMinMAC = false; valid_args[ "minMAC" ] = false;
    }
    if(params->test_mode && params->min_MAC < 0.5)
      throw "minimum MAC must be at least 0.5.";
    if(!params->test_mode && params->setMinINFO){
      sout << "WARNING: option --minINFO only works in step 2 of REGENIE.\n";
      params->setMinINFO = false; valid_args[ "minINFO" ] = false;
    }
    if( !params->split_by_pheno && params->w_interaction){
      sout << "WARNING: option --no-split does not work for interaction tests.\n";
      params->split_by_pheno = true; valid_args[ "no-split" ] = false;
    }
    if( !params->split_by_pheno && params->nnls_out_all){
      sout << "WARNING: option --no-split does not work with --sbat-verbose.\n";
      params->split_by_pheno = true; valid_args[ "no-split" ] = false;
    }
    if( vm.count("no-split") && vm.count("htp")){
      sout << "WARNING: option --no-split cannot be used with --htp and will be ignored.\n";
      valid_args[ "no-split" ] = false;
    }
    if(params->uncapped_pvals && !params->htp_out){
      sout << "WARNING: option --exact-p must be used with --htp.\n";
      params->uncapped_pvals = false; valid_args[ "exact-p" ] = false;
    }
    if( (!params->test_mode || (params->trait_mode!=1) || params->htp_out || !params->split_by_pheno) && params->af_cc ) {
      sout << "WARNING: disabling option --af-cc (only for BTs in step 2 in native output format split by trait).\n";
      params->af_cc = false; valid_args[ "af-cc" ] = false;
    }
    if(params->test_mode) check_build_code(params);
    if(params->rm_snps && params->keep_snps )
      sout << "WARNING: only variants which satisfy both extract/exclude options will be kept.\n";

    if(params->test_mode && (params->min_INFO < 0 || params->min_INFO > 1) )
      throw "minimum info score must be in [0,1].";
    if( params->rm_missing_qt && (params->strict_mode || params->trait_mode || !params->test_mode) ) params->rm_missing_qt = false;

    if( !vm.count("bsize") && !params->snp_set ) 
      throw "must specify the block size using '--bsize'.";
    else if(vm.count("bsize") && ( params->block_size < 1 ))
      throw "block size must be at least 1.";
    if(params->set_aaf && !params->build_mask) params->set_aaf = false;
    if(params->run_l0_only && params->test_l0)
      throw "cannot use --test-l0 with --run-l0";
    if(params->test_l0 && params->print_block_betas) 
      throw "cannot use --test-l0 with --print";
    if(params->test_l0 && (params->l0_pvals_file != ""))
      throw "--select-l0 must be specified without an argument";
    if(params->print_cov_betas && (params->w_interaction || params->blup_cov))
      throw "cannot use --print-cov-betas with interaction tests or --prs-cov";
    if(params->print_cov_betas && !params->test_mode)
      throw "can only use --print-cov-betas in step 2";

    // determine number of threads if not specified
    if(params->threads < 1)
      params->threads = std::max(1u, std::thread::hardware_concurrency() - 1); //may return 0 when not able to detect

    // check parallel l0
    if(params->test_mode && 
        (vm.count("split-l0")||vm.count("run-l0")||vm.count("run-l1")) ) {
      sout << "WARNING: options --split-l0/--run-l0/--run-l1 only work in step 1.\n";
      params->split_l0 = params->run_l0_only = params->run_l1_only = false;
      valid_args[ "split-l0" ] = valid_args[ "run-l0" ] = valid_args[ "run-l1" ] = false;
    } else if( vm.count("nb") && 
        (vm.count("split-l0")||vm.count("run-l0")||vm.count("run-l1")) ) {
      sout << "WARNING: options --split-l0/--run-l0/--run-l1 cannot be used with --nb.\n";
      params->split_l0 = params->run_l0_only = params->run_l1_only = false;
      valid_args[ "split-l0" ] = valid_args[ "run-l0" ] = valid_args[ "run-l1" ] = false;
    }
    if(params->test_l0 && 
        (vm.count("split-l0")||vm.count("run-l0")||vm.count("run-l1")) ) 
      throw "cannot use --test-l0 with --split-l0/--run-l0/--run-l1";
    if( vm.count("run-l0") || vm.count("run-l1") ) 
      check_file(files->split_file, "run-l0/l1");

    // set Firth as default if both Firth and SPA are specified
    if(params->use_SPA && params->firth) {
      sout << "WARNING: only one of --firth/--spa can be used. Only Firth will be used.\n";
      params->use_SPA = false; valid_args[ "spa" ] = false;
    }
    params->mk_snp_map = params->rm_snps || params->keep_snps || params->rm_or || params->keep_or || params->snp_set || params->getCorMat || (params->forced_MAC > 0);
    params->keep_snp_map = params->rm_or || params->keep_or || params->snp_set || params->getCorMat || (params->forced_MAC > 0);

    // check fallback pvalue threshold
    if((params->alpha_pvalue < params->nl_dbl_dmin) || (params->alpha_pvalue > 1 - params->numtol) )
      throw "Fallback p-value threshold must be in (0,1).";
    if(params->firth_approx && !params->firth) {
      params->firth_approx = false; valid_args[ "approx" ] = false;
    }
    if( params->skip_cf_burden && !(params->use_SPA || params->firth) ) {
      params->skip_cf_burden = false; valid_args[ "skip-cf-burden" ] = false;
    }

    // check arguments for logistic regression 
    if(params->trait_mode && (params->niter_max < 1))
      throw "invalid argument for --niter (must be positive integer).";
    if(params->firth && (params->maxstep_null < 1))
      throw "invalid argument for --maxstep-null (must be a positive integer).";
    if(params->firth && (params->niter_max_firth_null < 1))
      throw "invalid argument for --maxiter-null (must be a positive integer).";
    if(params->nChrom < 2)
      throw "invalid argument for --nauto (must be > 1).";
    if(params->set_range && (params->range_chr == -1))
      throw "unrecognized chromosome in --range.";
    if(params->rm_indivs && params->keep_indivs )
      throw "cannot use both --keep and --remove.";
    if(params->rm_or && params->keep_or )
      throw "cannot use both --extract-or and --exclude-or.";
    if(params->condition_file && !params->condition_snps )
      throw "must use --condition-list if using --condition-file.";
    if(params->interaction_file && !params->interaction_snp )
      throw "must use --interaction-snp if using --interaction-file.";
    if( !vm.count("covarFile") && vm.count("interaction") && !vm.count("interaction-snp") )
      throw "must use --covarFile if using --interaction.";

    if( params->test_mode && params->select_chrs && in_map(-1, filters->chrKeep_test) )
      throw "invalid chromosome specified by --chr/--chrList.";

    if(params->test_mode && !params->skip_blups && !vm.count("pred")) 
      throw "must specify --pred if using --step 2 (otherwise use --ignore-pred).";
    if(vm.count("interaction") || vm.count("interaction-snp") || vm.count("interaction-prs")){
      if(!vm.count("interaction-snp") && !vm.count("interaction-prs") && (!vm.count("covarFile") || !params->test_mode) )
        throw "can only use --interaction with --covarFile in step 2.";
      if( (vm.count("interaction") + vm.count("interaction-prs") + vm.count("interaction-snp")) > 1 ) 
        throw "must only specify single interacting variable";
      if(params->use_SPA)
        throw "cannot use --interaction with SPA test.";
      if(vm.count("interaction-snp") && vm.count("use-prs"))
        throw "cannot use --interaction-snp with full PRS.";
      if(vm.count("firth") && !vm.count("approx")){
        sout << "WARNING: using approximate Firth for association testing.\n";
        params->firth_approx = true;
      }
    }
    if(params->skip_test && params->vc_test)
      throw "cannot use '--skip-test' with SKAT/SKATO/ACATO\n";
    if(params->vc_test && params->firth && !params->firth_approx){
      sout << "WARNING: Using approximate Firth for association testing.\n";
      params->firth_approx = true;
    }
    if((params->skato_rho.size() > 0) && (params->skato_rho <0 || params->skato_rho >1).any())
      throw "rho values for SKAT-O must be in [0,1]";
    if(params->singleton_carriers && params->aaf_file_wSingletons){
      sout << "WARNING: Ignoring option --singleton-carrier when using --set-singletons.\n";
      params->singleton_carriers = false; valid_args[ "singleton-carrier" ] = false;
    }

    //params->use_max_bsize = params->mask_loo;
    if( (params->trait_mode==2) && params->w_interaction)
      throw "cannot use interaction tests with count phenotypes.";
    if( params->interaction_prs && !(vm.count("use-prs") || vm.count("pred")) )
      throw "must supply step 1 predictions.";

    if(vm.count("force-ltco") && vm.count("use-prs"))
      throw "cannot use LTCO with full PRS.";
    if(vm.count("use-null-firth") && !params->firth_approx) 
      throw "option --use-null-firth only wors with approximate Firth test.";
    if(vm.count("write-null-firth") && 
        ( (params->test_mode && !params->firth_approx) || (!params->test_mode && (params->trait_mode!=1)) ) ) {
      sout << "WARNING: option --write-null-firth only works for BTs with approximate Firth test.\n";
      params->write_null_firth = false; valid_args[ "write-null-firth" ] = false;
    }
    if( (filters->cov_colKeep_names.size() > 0) && !vm.count("covarFile") )
      throw "you specified covariates without specifying a covariate file (using --covarFile).";

    if(params->transposedPheno){
      if(vm.count("phenoFile") ) 
        throw "cannot use both --phenoFile and --tpheno-file.";
      if(!vm.count("tpheno-indexCol") ) 
        throw "must specify --tpheno-indexCol with --tpheno-file.";
      if(vm.count("iid-only"))
        params->tpheno_iid_only = true;
    }
    if(vm.count("starting-block")){
      if(!params->test_mode)
        throw "option --starting-block only works in step 2";
      else if(params->start_block<1)
        throw "starting block must be >=1";
      if(vm.count("nb")) params->n_block += params->start_block - 1;
    }
    if(vm.count("sex-specific") && (params->file_type == "bgen") && !params->bgenSample)
      throw "must specifying sample file using --sample for sex-specific analyses";
    if(params->blup_cov && (!(params->use_prs || !params->skip_blups) || !params->test_mode || (params->firth && !params->firth_approx) )){
      params->blup_cov = false; valid_args[ "prs-cov" ] = false;
    }

    if(params->test_mode && (params->file_type == "pgen") && !params->fastMode)
      throw "cannot use --nostream with PGEN format.";

    // check apply-rint options
    if(params->rerint & params->rerintcov)
      throw "must select one of the two options, --apply-rerint or --apply-rerint-cov";

    // check multi-trait settings
    if(params->trait_set) {
      if(!params->strict_mode) 
        throw "--strict mode is required for multi-trait tests";
      if(params->split_by_pheno) 
        throw "--no-split mode is required for multi-trait tests";
    }

    if(params->mcc_skew < 0) {
        throw "absolute phenotypic skewness must be positive";
    }
    if(params->mcc_skew > 0) {
      if(!params->mcc_test) {
        throw "--mcc must be on when specifying absolute phenotypic skewness";
      }
    }
    if(params->mcc_test) {
      // convert mcc thr. from raw to -log10 scale
      if((params->mcc_thr > 1) && (params->mcc_thr <= 0))
        throw "--mcc-thr must be in (0; 1]";
      if(params->mcc_thr < 1) 
        params->mcc_apply_thr = true;
      params->mcc_thr_nlog10 = -log10(params->mcc_thr); // -log10 transformation
      // debug
      /* cout << "mcc_test = " << params->mcc_test << " | mcc_apply_thr = " << params->mcc_apply_thr << " | mcc_thr  = " << params->mcc_thr << " | mcc_thr_nlog10 = " << params->mcc_thr_nlog10 << " | mcc_skew = " << params->mcc_skew << endl; */
    }

    // check MultiPhen-trait settings
    if(params->multiphen) {
      if(!params->strict_mode) throw "--strict mode is required for MultiPhen test";
      /* if(params->split_by_pheno) throw "--no-split mode is required for MultiPhen test"; */
      if((params->multiphen_thr > 1) && (params->multiphen_thr <= 0)) throw "--multiphen-thr must be in (0; 1]";
      params->n_tests_multitrait = 1; // a single test = MultiPhen
      params->split_by_multitrait = false; // no split of output files
    }

    // check input files
    if(params->file_type == "bgen") {
      check_file (files->bgen_file, "bgen"); 
      if(params->bgenSample) check_file (files->sample_file, "sample"); 
      if(files->bgi_file != "") {
        check_file (files->bgi_file, "bgi");
        params->with_bgi = true;
      } else {
        files->bgi_file = files->bgen_file + ".bgi";
        params->with_bgi = file_exists (files->bgi_file) ;
      }
    }
    if(vm.count("covarFile")) check_file(files->cov_file,"covarFile");
    if(!params->getCorMat) check_file(files->pheno_file,"phenoFile"); 
    if(params->file_type == "bed"){
      vector<string> suffs = {".bed",".bim",".fam"};
      check_file(files->bed_prefix, suffs, "bed");
    }
    if(params->file_type == "pgen"){
      vector<string> suffs = {".pgen",".pvar",".psam"};
      check_file(files->pgen_prefix, suffs, "pgen");
    }
    if(params->keep_indivs)
      for(auto cn : files->file_ind_include)
        check_file(cn, "keep");
    if(params->rm_indivs)
      for(auto cn : files->file_ind_exclude)
        check_file(cn, "remove");
    if(!vm.count("run-l0") && params->keep_snps)
      for(auto cn : files->file_snps_include)
        check_file(cn, "extract");
    if(params->rm_snps)
      for(auto cn : files->file_snps_exclude)
        check_file(cn, "exclude");
    if(params->keep_or)
      for(auto cn : files->file_snps_include_or)
        check_file(cn, "extract-or");
    if(params->rm_or)
      for(auto cn : files->file_snps_exclude_or)
      check_file(cn, "exclude-or");
    if(params->snp_set) {
      check_file(files->set_file, "set-list");
      if(!vm.count("extract-setlist") && params->keep_sets)
        for(auto cn : files->file_sets_include)
          check_file(cn, "extract-sets");
      if(!vm.count("exclude-setlist") && params->rm_sets)
        for(auto cn : files->file_sets_exclude)
          check_file(cn, "exclude-sets");
    }
    if(params->forced_MAC > 0) check_file(params->forced_MAC_snpfile, "force-mac-filter");
    if(params->select_l0 && !params->test_l0)
      check_file(params->l0_pvals_file, "select-l0");
    if(vm.count("ld-extract"))
      check_file(params->ld_list_file, "ld-extract");
    if(params->build_mask){
      check_file(files->anno_file, "anno-file");
      check_file(files->mask_file, "mask-def");
      if(vm.count("anno-labels")) check_file(files->anno_labs_file, "anno-labels");
    }
    if(params->set_aaf) check_file(files->aaf_file, "aaf-file");
    if(params->condition_snps) {
      check_file(files->condition_snps_list, "condition-list");
      if(params->condition_file && (files->condition_snps_info.format == "bgen")) {
        check_file(files->condition_snps_info.file, "condition-file");
        // optional sample file?
        files->condition_snps_info.with_sample = files->condition_snps_info.sample != "";
        if(files->condition_snps_info.with_sample) check_file(files->condition_snps_info.sample, "condition-file-sample");
        files->condition_snps_info.with_bgi = file_exists (files->condition_snps_info.file + ".bgi") ;
      } else if(params->condition_file && (files->condition_snps_info.format == "bed")) {
        vector<string> suffs = {".bed",".bim",".fam"};
        check_file(files->condition_snps_info.file, suffs, "condition-file");
      } else if(params->condition_file && (files->condition_snps_info.format == "pgen")) {
        vector<string> suffs = {".pgen",".pvar",".psam"};
        check_file(files->condition_snps_info.file, suffs, "condition-file");
      }
    }
    if(params->interaction_snp && params->interaction_file) {
      if(files->interaction_snp_info.format == "bgen") {
        check_file(files->interaction_snp_info.file, "interaction-file");
        // optional sample file?
        files->interaction_snp_info.with_sample = files->interaction_snp_info.sample != "";
        if(files->interaction_snp_info.with_sample) check_file(files->interaction_snp_info.sample, "interaction-file-sample");
        files->interaction_snp_info.with_bgi = file_exists (files->interaction_snp_info.file + ".bgi") ;
      } else if(files->interaction_snp_info.format == "bed") {
        vector<string> suffs = {".bed",".bim",".fam"};
        check_file(files->interaction_snp_info.file, suffs, "interaction-file");
      } else if(files->interaction_snp_info.format == "pgen") {
        vector<string> suffs = {".pgen",".pvar",".psam"};
        check_file(files->interaction_snp_info.file, suffs, "interaction-file");
      }
    }
    if(vm.count("lovo-snplist")) check_file(params->masks_loo_snpfile, "lovo-snplist");

    if(vm.count("remeta-save-ld") > 0) {
  #ifndef WITH_HTSLIB
    throw "--remeta-save-ld option requires compilation with HTSlib";
  #else
      params->remeta_save_ld = true;
      if(vm.count("remeta-ld-spr") > 0) {
        params->remeta_ld_spr = vm["remeta-ld-spr"].as<double>();
      }
      if(params->skat_collapse_MAC > 0) {
        throw "--remeta-save-ld option requires --vc-MACthr 0";
      }
  #endif
    }

    check_seed(params->rng_seed, vm.count("seed"));
    print_args(arguments, valid_args, sout);

  } catch (const cxxopts::OptionException& e) {
    if (sout.coss.is_open())
      print_header(sout.coss);
    print_header(cout);
    sout << "ERROR: " << e.what() << endl << params->err_help << "\n";
    exit(EXIT_FAILURE);
  } catch (const std::string& msg) {// after opening sout
    sout <<  "ERROR: " <<  msg << "\n" <<  params->err_help << "\n";
    exit(EXIT_FAILURE);
  } catch (const char* msg) {// after opening sout
    std::string str_msg = msg;
    sout <<  "ERROR: " <<  str_msg << "\n" <<  params->err_help << "\n";
    exit(EXIT_FAILURE);
  }

  return;
}

void check_file(string const& infile, string const& option_name){

  if(infile == "") 
    throw "Invalid argument (=' ') specified for option --" + option_name;
  else if(!file_exists (infile))
    throw infile + " doesn't exist for option --" + option_name;

}

void check_file(string const& infile, vector<string> const& suffixes, string const& option_name){

  if(infile == "") 
    throw "Invalid file argument (=' ') specified for option --" + option_name;
  for(auto suffix : suffixes)
    // allow for gzipped bim/fam/psam/pvar files
    if(!file_exists (infile + suffix) && ((suffix == ".bed") || (suffix == ".pgen") || !file_exists(infile + suffix + ".gz")))
      throw infile + suffix + " doesn't exist for option --" + option_name;

}


void start_log(const string& out_file, MeasureTime* mt, mstream& sout){

  string log_name = out_file + ".log";
  sout.coss.open(log_name.c_str(), ios::out | ios::trunc); 
  if (!sout.coss.is_open()) {
    print_header(cout);
    cout << "ERROR: Cannot write log file '" << log_name << "'\n" ;
    exit(EXIT_FAILURE);
  } 

  mt->init();
  sout << "Start time: " << ctime( &(mt->start_time_info) ) << endl; 
  print_header(sout.coss);
  print_header(cout);
  sout << "Log of output saved in file : " << log_name << endl<< endl;

}

template <typename T> 
void print_args(T arguments, map<string,bool>& amap, mstream& sout){

  // print options
  sout << "Options in effect:\n";

  for(size_t counter = 0; counter < arguments.size(); counter++){	  
    if(!amap[arguments[counter].key()]) continue;

    sout << "  --" << arguments[counter].key();
    if(arguments[counter].value() != "true") sout << " " << arguments[counter].value(); 

    if(counter < (arguments.size() - 1)) sout << " \\";
    sout << "\n";
  }

  sout << "\n";

}

ArrayXd get_unit_params(bool const& incl_bound, string const& opt, vector<string> const& str_vec, struct param const* params, mstream& sout){

  std::vector<double> vals;

  for( size_t val = 0; val < str_vec.size(); val++)
    vals.push_back(convertDouble( str_vec[val], params, sout));
  std::sort(vals.begin(), vals.end());
  vals.erase( unique( vals.begin(), vals.end() ), vals.end() );

  ArrayXd vvals = MapArXd( vals.data(), vals.size() ); 
  // check parameters
  if( incl_bound && ((vvals<0) || (vvals>1)).any() )
    throw "must specify values for " + opt + " in [0,1].";
  else if( !incl_bound && ((vvals<=0) || (vvals>=1)).any() )
    throw "must specify values for " + opt + " in (0,1).";

  return vvals;

}

void set_ridge_params(int const& nparams, ArrayXd& vec, mstream& sout){

  if(nparams < 2)
    throw "number of ridge parameters must be at least 2 (=" + to_string( nparams ) + ")";

  // endpoints are 0.01 and 0.99 
  double step = 1.0 / ( nparams - 1 );
  vec = ArrayXd::LinSpaced(nparams, 0, nparams-1) * step;
  vec.head(1) = 0.01;
  vec.tail(1) = 0.99;

}

void print_usage_info(struct param const* params, struct in_files* files, mstream& sout){

  double total_ram;
  string ram_unit;

  ///// Memory usage
  if(!params->test_mode){
    // Step 1
    // 4P + max( B + PRT, PRT) + #chrs [P:#traits;R=#ridge l0;T=#predictions from l0]
    int t_eff = ( params->write_l0_pred ? 1 : params->total_n_block );
    int p_eff = ( params->write_l0_pred ? 1 : params->n_pheno );
    int b_eff = params->total_n_block;

    total_ram = 4 * params->n_pheno + params->nChrom + params->ncov;
    total_ram += std::max( params->block_size + params->n_pheno * params->n_ridge_l0 * t_eff, p_eff * params->n_ridge_l0 * b_eff );
  } else {
    // Step 2
    // 3P + B
    total_ram = params->n_pheno * 3 + params->block_size + params->ncov * 2; // y, mask, y_resid, g, X, X getbasis projection
    if(params->trait_mode) {
      total_ram += 3 * params->n_pheno + params->block_size + params->n_pheno * params->ncov; // y_raw, gamma_hat, gamma_hat_mask, g_resid
      if(params->use_SPA) total_ram += 0.5 * params->block_size; // non_zero_indices of g (4 bytes)
      if(params->firth_approx) total_ram += params->n_pheno; // cov offset
      if(params->start_block > params->total_n_block)
        throw "Starting block > number of blocks analyzed";
    } else total_ram += params->block_size; // for Gresid
    if((params->file_type == "bed") && params->fastMode) total_ram += params->block_size/4.0/sizeof(double); //for extracting snp_data_block
    if(params->use_max_bsize) total_ram += params->block_size; // loo masks
    // for Hmat (G_E, G, G*E )
    if(params->w_interaction) 
      total_ram += params->threads * ((params->gwas_condtl ? 1 : 2) * params->ncov_interaction + 1); 
  }

  total_ram *= params->n_samples * sizeof(double);
  total_ram += params->nvs_stored * sizeof(struct snp);
  if( params->getCorMat ){ // M^2 (x2 with txt output) + 3NB
      total_ram += (params->cor_out_txt && (params->ld_sparse_thr == 0) ? 2 : 1) * params->extract_vars_order.size() * params->extract_vars_order.size() * sizeof(double);
      total_ram += params->n_samples * params->block_size * sizeof(double);
  }
  if( params->use_loocv ) total_ram += params->chunk_mb * 1e6; // max amount of memory used for LOO computations involved
  if( params->mask_loo ) total_ram += 1e9; // at most 1GB
  if( params->vc_test ) total_ram += 2 * params->max_bsize * params->max_bsize * sizeof(double); // MxM matrices
  total_ram /= 1000.0 * 1000.0; 
  if( total_ram > 1000 ) {
    total_ram /= 1000.0; 
    ram_unit = "GB";
  } else ram_unit = "MB";

  int ram_int = (int) ceil( total_ram );
  sout << " * approximate memory usage : " << ram_int << ram_unit << endl;

  ///// Disk space usage
  if(!params->test_mode && !params->run_l1_only && params->write_l0_pred){
    if(files->loco_tmp_prefix.empty()) files->loco_tmp_prefix = files->out_file;
    sout << " * writing level 0 predictions to disk" << endl;
    sout << "   -" << (params->rm_l0_pred ? "temporary " : "") << "files will have prefix [" << files->loco_tmp_prefix << "_l0_Y]" << endl;
    // N*P*T*R
    int b_eff = params->total_n_block;
    total_ram = params->n_pheno * b_eff * params->n_ridge_l0;
    total_ram *= params->n_samples * sizeof(double);
    total_ram /= 1024.0 * 1024.0; 
    if( total_ram > 1000 ) {
      total_ram /= 1024.0; 
      ram_unit = "GB";
    } else ram_unit = "MB";
    int ram_int = (int) ceil( total_ram );
    sout << "   -approximate disk space needed : " << ram_int << ram_unit << endl;
  }

  if(params->debug)
    sout << " * rng seed : " << params->rng_seed << "\n";
}

int chrStrToInt(const string& chrom, const int& nChrom) {

  // if label is chr1, chr2,...
  string s_chr = std::regex_replace(chrom, std::regex(R"(^chr)"), "");

  if (isdigit(s_chr[0])) {
    int chr = atoi(s_chr.c_str());
    if((chr >= 1) && (chr <= nChrom)) return chr;
  } else if ( (s_chr == "X") || (s_chr == "XY") || (s_chr == "Y") || (s_chr == "PAR1") || (s_chr == "PAR2") ) return nChrom;

  return -1;
}

vector<string> check_name(string const& str, mstream& sout){

  int imin, imax;
  size_t pos_start = 0, pos_end; 
  string name, pref, suf, strerror;
  strerror = "invalid string expansion (=" + str + ").\n";
  vector<string> strout;

  if(str.size() == 0) return strout;

  pos_end = str.find("{"); 
  if(pos_end == std::string::npos) {
    strout.push_back(str); return strout;
  }

  try {
    // prefix if present
    name = str.substr (pos_start, pos_end - pos_start);
    pref = name;

    // find :
    pos_start = pos_end + 1, pos_end = str.find(":"); 
    if(pos_end == std::string::npos) throw strerror;
    name = str.substr (pos_start, pos_end - pos_start);
    imin = stoi( name );

    // find }
    pos_start = pos_end+1, pos_end = str.find("}"); 
    if(pos_end == std::string::npos) throw strerror;
    name = str.substr (pos_start, pos_end - pos_start);
    imax = stoi( name );

  } catch (const std::invalid_argument& ia){ 
    throw strerror ;
  } 

  // suffix is present
  suf = str.substr (pos_end+1, std::string::npos);

  for(int j = imin; j <= imax; j++){
    name = pref + to_string(j) + suf;
    strout.push_back(name);
  }

  return strout;
}

void check_build_code(struct param* params){
  vector<string> valid_codes = { "b36", "b37", "b38", "hg18", "hg19", "hg38"};

  if (std::find(valid_codes.begin(), valid_codes.end(), params->build_code) == valid_codes.end()){ // format: <end_par1>,<start_par2>
    int min_npar, max_npar;
    if((sscanf( params->build_code.c_str(), "%d,%d", &min_npar, &max_npar ) != 2) || (min_npar < 1) || (max_npar < 1) || (max_npar < min_npar)) 
      throw "invalid build code given (valid ones are '" + print_sv(valid_codes, "|") + "' or [start,end] position of the non-par region)"; 
    params->par1_max_bound = min_npar - 1;
    params->par2_min_bound = max_npar + 1;

  } else if((params->build_code == "b36") || (params->build_code == "hg18")){
    params->par1_max_bound = 2709520, params->par2_min_bound = 154584238;
  } else if ((params->build_code == "b37") || (params->build_code == "hg19")){
    params->par1_max_bound = 2699520, params->par2_min_bound = 154931044;
  } else{
    params->par1_max_bound = 2781479, params->par2_min_bound = 155701383;
  }

}

double convertDouble(const string& val, struct param const* params, mstream& sout){

  if(val == params->missing_pheno_str)
    return params->missing_value_double;
  else if( (val == "nan") || (val == "inf") )
    return params->missing_value_double;

  double dval;
  if(sscanf(val.c_str(), "%lf", &dval) != 1)
    throw "could not convert value to double: '" + val + "'";

  return dval;
}

float convertFloat(const string& val, struct param const* params, mstream& sout){

  if(val == params->missing_pheno_str)
    return params->missing_value_float;
  else if( (val == "nan") || (val == "inf") )
    return params->missing_value_float;

  float dval;
  if(sscanf(val.c_str(), "%f", &dval) != 1)
    throw "could not convert value to float: '" + val + "'";

  return dval;
}

string convert_logp_raw(double const& logp, double const& log_dbl_min){

  char pval_str[256];

  if(logp <= 3) {
    sprintf(pval_str, "%f", pow(10, -logp));
  } else if(logp <= log_dbl_min) {
    sprintf(pval_str, "%g", pow(10, -logp));
  } else {
    double thr = log(9.95)/log(10);
    int base = ceil(logp);
    double res = base - logp;
    if(res>=thr) {res = 0; base++;}
    sprintf(pval_str, "%.1fe-%d", pow(10, res), base);
  }

  return( string(pval_str) );
}

// convert to numerical category using map
double convertNumLevel(const string& val, std::map<std::string,int>& cmap, struct param const* params, mstream& sout){

  if(val == params->missing_pheno_str)
    return params->missing_value_double;
  else if( (val == "nan") || (val == "inf") )
    return params->missing_value_double;

  if(in_map(val, cmap)) 
    return cmap[val];

  // add to map
  int newcat = cmap.size();
  cmap[val] = newcat;

  return newcat;
}

// for strings with format: str[lvl]
void check_inter_var(std::string& str, std::string& lvl, mstream& sout){

  string name;
  size_t pos_start = 0, pos_end; 

  // check if contains "["
  pos_end = str.find("["); 
  if(pos_end == std::string::npos) 
    return ;
  name = str.substr (pos_start, pos_end - pos_start);

  // find "]"
  pos_start = pos_end + 1, pos_end = str.find("]"); 
  if(pos_end == std::string::npos) 
    throw "ERROR: Invalid string :" + str ;

  lvl = str.substr (pos_start, pos_end - pos_start);

  str = name;
}

// comma separated strings
std::string print_csv(const vector<string>& vlist){
  return print_sv(vlist, ",");
}

// semi-colon separated strings
std::string print_scsv(const vector<string>& vlist){
  return print_sv(vlist, ";");
}

template <typename T>
std::string print_sv(const std::vector<T>& vlist, const string& delim)
{
  std::ostringstream buffer;
  if(!vlist.empty()) {
    std::copy(std::begin(vlist), std::end(vlist) - 1, std::ostream_iterator<T>(buffer, delim.c_str()));
    buffer << vlist.back(); // last element
  }
  return buffer.str();
}

void removeCarriageReturn(std::string& str) {
  if (!str.empty() && str.back() == '\r') {
    str.pop_back();
  }
}

Eigen::ArrayXi get_true_indices(const Ref<const ArrayXb>&  bool_arr){

  ArrayXi v_indices ( bool_arr.count() );
  for(int i = 0, j = 0; i < bool_arr.size(); i++)
    if(bool_arr(i)) v_indices(j++) = i;

  return v_indices;
}

void get_both_indices(std::vector<Eigen::ArrayXi>& res, const Eigen::Ref<const ArrayXb>& bool_arr){

  res.resize(2);
  int Ntot = bool_arr.size();
  res[0].resize(bool_arr.count()); // true entries
  res[1].resize(Ntot - res[0].size()); // false entries
  for(int i = 0, j_t = 0, j_f = 0; i < Ntot; i++) {
    if(bool_arr(i)) res[0](j_t++) = i;
    else res[1](j_f++) = i;
  }

}
void get_both_indices(std::vector<Eigen::ArrayXi>& res, const Eigen::Ref<const ArrayXb>& bool_arr, const Eigen::Ref<const ArrayXb>& mask){

  res.resize(2);
  int Ntot = mask.count();
  res[0].resize((mask && bool_arr).count()); // true entries
  res[1].resize(Ntot - res[0].size()); // false entries
  for(int i = 0, j_t = 0, j_f = 0; i < bool_arr.size(); i++) {
    if(mask(i)){
      if(bool_arr(i)) res[0](j_t++) = i;
      else res[1](j_f++) = i;
    }
  }

}

// get logp from chisq(1)
void get_logp(double& logp, const double& Tstat){

  boost::math::chi_squared chisq1(1);

  if( (Tstat < 0) && (fabs(Tstat) < 1e-6)){logp = 0; return;} // num err
  else if(Tstat<0) {logp = -1; return;} // fail
  double pv = cdf(complement(chisq1, Tstat));

  if(pv == 0) logp = log10(2) - 0.5 * log10( 2 * M_PI * Tstat ) - 0.5 * Tstat * M_LOG10E ;
  else logp = log10(pv);

  logp *= -1;

}

// get logp & chisq1 from pv
void get_logp(const double& pv, double& logp, double& Tstat, double const& dbl_dmin){

  if((pv < 0) || (pv > 1)) { // fail
    logp = -1; 
    Tstat = 0;
    return;
  }

  boost::math::chi_squared chisq1(1);

  double pval = max(dbl_dmin, pv); // to prevent underflow
  Tstat = quantile(complement(chisq1, pval)); // chisq stat
  logp = -log10(pval); // -log10p

}

// get logp from chisq(k)
void get_logp(double& logp, const double& Tstat, double const& df){

  boost::math::chi_squared chisqK(df);

  if( (Tstat < 0) && (fabs(Tstat) < 1e-6)){logp = 0; return;} // num err
  else if(Tstat<0) {logp = -1; return;} // fail
  double pv = cdf(complement(chisqK, Tstat));

  if(pv == 0) logp = log10(2) - 0.5 * df * log10(2) - boost::math::lgamma(df * 0.5) / log(10) + 0.5 * (df-2) * log10(Tstat) - 0.5 * Tstat * M_LOG10E ;
  else logp = log10(pv);

  logp *= -1;

}

// get chisq1 & pval from logp
void get_chisq_stat_pv(double& pv, double& Tstat, double const& logp, double const& dbl_dmin, double const& log10_dbl_dmin){

  if(logp<0) { // fail
    pv = -1; 
    Tstat = 0;
    return;
  }

  boost::math::chi_squared chisq1(1);

  if(logp > log10_dbl_dmin){
    double val = logp * log(100) + log(2/M_PI);
    Tstat = val - log(val); // approximation for small p-values using Lambert W function
    pv = dbl_dmin; // prevent underflow
  } else {
    pv = pow(10, -logp);
    Tstat = quantile(complement(chisq1, pv)); // chisq stat
  }

}

void allocate_mat(MatrixXd& M, int const& nrows, int const& ncols){
  M.resize(nrows, ncols);
}

std::string print_mat_dims(MatrixXd const& mat){
  std::ostringstream buffer;
  buffer << "#rows=" << mat.rows() << " | #cols=" <<  mat.cols();
  return buffer.str();
}

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

void print_obj(const Ref<const MatrixXd>& mat, string const& fname){
  // write obj to file
  IOFormat Fmt(FullPrecision, DontAlignCols, " ", "\n", "", "","","");
  ofstream ofile;
  ofile.open(fname);
  ofile << mat.format(Fmt) << "\n";
  ofile.close();
}

int get_mem(){ // in MB
    FILE* file = fopen("/proc/self/status", "r");
    double result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line) / 1024.0;
            break;
        }
    }
    fclose(file);
    return result;
}

std::string print_mem(){
  return "memory usage=" + to_string( get_mem() ) + "MB";
}

void set_threads(struct param* params) {

#if defined(_OPENMP)
  omp_set_num_threads(params->threads); // set threads in OpenMP
  params->neff_threads = params->threads;
#endif
#if defined(WITH_MKL)
  mkl_set_num_threads(params->threads);
#endif
  setNbThreads(params->threads);

}

void check_seed(uint& seed, bool const& skip_gen_seed) {
  if(!skip_gen_seed) {
    std::random_device rd;
    seed = rd();
  }
  srand(seed);
}
