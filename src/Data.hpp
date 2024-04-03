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

#ifndef DATA_H
#define DATA_H

class Data {

  public:
    // class elements
    mstream sout;
    MeasureTime runtime;
    param params;
    in_files files;
    filter in_filters;
    std::vector<snp> snpinfo;
    phenodt pheno_data;
    geno_block Gblock;
    std::map<int, std::vector<int>> chr_map; // first=chr; second=[# SNPs analyzed, #blocks, # SNPs in file]
    ests m_ests;
    ridgel1 l1_ests;
    f_ests firth_est;
    // HLM
    HLM nullHLM; // for null model fitting of HLM
    remeta_sumstat_writer remeta_sumstats;

    std::string model_type, correction_type, test_string, wgr_string;

    uint32_t n_corrected = 0; // to keep track of how many SNPs require correction
    bool pval_converged = false; // keep track of whether SPA/Firth converged
    bool fastSPA; // use fast approx. for rare SNPs

    std::vector < MatrixXb > masked_in_folds;
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > predictions;

    uint32_t total_chrs_loco;
    Eigen::MatrixXd blup;
    Eigen::VectorXd denum_tstat;
    Eigen::MatrixXd res, stats, W_hat;
    Eigen::RowVectorXd p_sd_yres;
    Eigen::VectorXd scale_G; // keep track of sd(Y) (1xP) and sd(G) (M*1)
    MultiPhen mphen;

    // function definitions
    void run();
    void run_step1();
    void run_step2();

    void file_read_initialization();
    void residualize_genotypes();
    void scale_genotypes(bool);
    void get_block_size(int const&,int const&,int const&,int&);

    // step 1 
    void set_parallel_l0();
    void write_l0_master();
    void prep_parallel_l0();
    void prep_parallel_l1();
    void set_blocks();
    void set_folds();
    void setmem();
    void calc_cv_matrices(struct ridgel0*);
    void level_0_calculations();
    void prep_l1_models();
    void write_inputs(); 
    void exit_early();
    // output of step 1
    void output();
    void make_predictions(int const&,int const&);
    void make_predictions_loocv(int const&,int const&);
    void make_predictions_binary(int const&,int const&);
    void make_predictions_binary_loocv_full(int const&,int const&);
    void make_predictions_binary_loocv(int const&,int const&);
    void make_predictions_count(int const&,int const&);
    void make_predictions_count_loocv(int const&,int const&);
    void write_predictions(int const&);
    std::string write_ID_header();
    std::string write_chr_row(int const&,int const&,const Eigen::Ref<const Eigen::VectorXd>&);
    void rm_l0_files(int const& ph);

    // step 2 main functions
    void test_snps();
    void set_blocks_for_testing();
    void print_test_info();
    void set_nullreg_mat();
    void compute_res();
    void residualize_res();
    void compute_res_bin(int const&);
    void compute_res_count(int const&);
    void setup_output(Files*,std::string&,std::vector<std::shared_ptr<Files>>&,std::vector<std::string>&);

    // step 2 using multithreading in eigen
    double check_pval(double const&,int const&,int const&,int const&);
    double run_firth_correction(int const&,int const&,int const&);
    void run_SPA_test(int const&);

    // step2 using multithreading in openmp
    void test_snps_fast();
    void analyze_block(int const&,int const&,tally*,std::vector<variant_block>&);
    void compute_tests_mt(int const&,std::vector<uint64>,std::vector<std::vector <uchar>>&,std::vector<uint32_t>,std::vector<uint32_t>&,std::vector<variant_block>&);
    void compute_tests_st(int const&,std::vector<uint64>,std::vector<std::vector <uchar>>&,std::vector<uint32_t>,std::vector<uint32_t>&,std::vector<variant_block>&);

    // step 2 with joint tests
    JTests jt;
    GenoMask bm;
    void test_joint();
    void set_groups_for_testing();
    void get_sum_stats(int const&,int const&,std::vector<variant_block>&);
    void readChunk(std::vector<uint64>&,int const&,std::vector<std::vector<uchar>>&,std::vector<uint32_t>&,std::vector<uint32_t>&,std::vector<variant_block>&);
    void getMask(int const&,int const&,std::vector<std::vector<uchar>>&,std::vector<uint32_t>&,std::vector<uint32_t>&,std::vector<variant_block>&);
    void getMask_loo(int const&,int const&,std::vector<std::vector<uchar>>&,std::vector<uint32_t>&,std::vector<uint32_t>&,std::vector<variant_block>&);

    // step 2 with multi-trait tests
    MTests mt;
    void test_multitrait();
    void analyze_block_multitrait(int const&,int const&,tally*,std::vector<variant_block>&);
    void compute_tests_mt_multitrait(int const&,std::vector<uint64>,std::vector<std::vector <uchar>>&,std::vector<uint32_t>,std::vector<uint32_t>&,std::vector<variant_block>&);
    void prep_multitrait(); 

    // step 2 with MultiPhen test
    /* MTests mt; */
    void test_multiphen();
    void analyze_block_multiphen(int const&,int const&,tally*,std::vector<variant_block>&);
    void compute_tests_mt_multiphen(int const&,std::vector<uint64>,std::vector<std::vector <uchar>>&,std::vector<uint32_t>,std::vector<uint32_t>&,std::vector<variant_block>&);
    void prep_multiphen(); 
    void set_multiphen();

    // for LD computation
    void ld_comp();
    void get_G_indices(Eigen::ArrayXi&,std::map<std::string,int>&);
    void write_snplist(ArrayXb&);
    // dosage-mode
    void compute_ld_dosages(Files*);
    void get_G_masks(SpMat&,ArrayXb&,std::map<std::string,int>&);
    void get_G_svs(int const&,int const&);
    void print_ld(MatrixXd&,Eigen::ArrayXi&,ArrayXb&,Files*);
    // hard-call mode
    void compute_ld_hardcalls(Files*);
    void get_G_svs(SpMat&,ArrayXb&,std::map<std::string,int>&);
    void get_G_masks_hc(SpMat&,ArrayXb&,std::map<std::string,int>&);
    void print_ld(SpMat&,Eigen::ArrayXi&,ArrayXb&,Files*);
    
    Data();
    ~Data();
};

// extra function
std::string get_fullpath(std::string);

#endif
