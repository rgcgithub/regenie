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
    spa_ests spa_est;

    std::string model_type, correction_type, test_string;

    uint n_corrected = 0; // to keep track of how many SNPs require correction
    bool pval_converged = false; // keep track of whether SPA/Firth converged
    bool fastSPA; // use fast approx. for rare SNPs

    std::vector < MatrixXb > masked_in_folds;
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > predictions;

    uint32_t size1, size2, snp_index_counter, total_chrs_loco;
    Eigen::MatrixXd blup;
    Eigen::VectorXd denum_tstat;
    Eigen::MatrixXd res, stats, W_hat;
    Eigen::RowVectorXd p_sd_yres;
    Eigen::VectorXd scale_G; // keep track of sd(Y) (1xP) and sd(G) (M*1)

    // function definitions
    void run();
    void file_read_initialization();
    void residualize_genotypes();
    void scale_genotypes(bool);

    // step 1 
    void set_parallel_l0();
    void write_l0_master();
    void write_snplist(std::string,int,int);
    void prep_parallel_l0();
    void prep_parallel_l1();
    void set_blocks();
    void set_folds();
    void setmem();
    void calc_cv_matrices(const int,struct ridgel0*);
    void level_0_calculations();
    void exit_early();
    // output of step 1
    void output();
    void make_predictions(const int,const int);
    void make_predictions_loocv(const int,const int);
    void make_predictions_binary(const int,const int);
    void make_predictions_binary_loocv(const int,const int);
    void write_predictions(const int);
    std::string write_ID_header();
    std::string write_chr_row(const int,const int,const Eigen::MatrixXd&);
    void rm_l0_files(int ph);

    // step 2 main functions
    void test_snps();
    void set_blocks_for_testing();
    void print_test_info();
    void blup_read_chr(const int);
    void compute_res();
    void compute_res_bin(int);
    void setup_output(Files*,std::string&,std::vector<Files*>&,std::vector<std::string>&);
    std::string print_header_output();
    std::string print_header_output_single();
    std::string print_header_output_htp();
    std::string print_sum_stats_head(const int);
    std::string print_sum_stats_head_htp(const int,const int,const std::string);
    std::string print_sum_stats(const double,const double,const double,const double,const int,const bool);
    std::string print_sum_stats(const double,const double,const int,const std::string,const double,const double,const double,const double,const bool);
    std::string print_sum_stats(const double,const double,const double,const double,const bool);
    std::string print_sum_stats_htp(const double,const double,const double,const double,const double,const double,const double,const Eigen::MatrixXd&,const double,const int,const bool);

    // step 2 using multithreading in eigen
    double check_pval(double,int,int,int);
    double run_firth_correction(int,int,int);
    void run_SPA_test(int);

    // step2 using multithreading in openmp
    void test_snps_fast();
    void analyze_block(const int&,const int&,tally*,std::vector<variant_block>&);
    void residualize_geno(int,variant_block*);
    void check_pval_snp(variant_block*,int,int,int);
    void run_firth_correction_snp(int,int,int,variant_block*);
    void run_SPA_test_snp(variant_block*,int,const Eigen::VectorXd&);

    // step 2 with joint tests
    JTests jt;
    GenoMask bm;
    void test_joint();
    void set_groups_for_testing();
    std::string build_mask_header();
    void get_sum_stats(const int,const int,std::vector<variant_block>&);
    void readChunk(const int,const int,const int,const int,std::vector<std::vector<uchar>>&,std::vector<uint32_t>&,std::vector<uint32_t>&,std::vector<variant_block>&);
    void getMask(const int,const int,std::vector<std::vector<uchar>>&,std::vector<uint32_t>&,std::vector<uint32_t>&,std::vector<variant_block>&);


    Data();
    ~Data();
};

// extra function
std::string get_fullpath(std::string);

#endif
