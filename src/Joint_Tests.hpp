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

#ifndef JOINT_H
#define JOINT_H


class JTests {

  public:
    std::map<std::string, uint16_t> joint_tests_map = { {"minp", 0}, {"ftest", 1}, {"gates", 2}, {"nnls", 3}, {"acat", 4}, {"nnls_pos", 5}, {"nnls_neg", 6}, {"gene_p", 7} };
    std::vector<std::string> test_names = {"MINP","F","GATES","NNLS","ACAT","NNLS_POS","NNLS_NEG", "GENE_P"};
    std::map <std::string, std::map <std::string, bool>> gene_p_tests;
    bool genep_all_masks = true;
    uint16_t test_list, qr_tests = 7ULL<<1; // 01110000 for ftest, gates and nnls

    // store variant set info (for each chr)
    std::vector<std::vector<vset>> setinfo;

    // for testing
    int df_test, nvars, ncovars; // df, number of variants which passed filters,#covariates
    bool nnls_verbose_out = false;
    bool nnls_normalize = true, nnls_strict = false;
    int nnls_napprox, nnls_verbose = 0;
    double acat_a1,acat_a2;
    bool valid_snp_mode, debug_mode;
    double pval, plog, zval, scale_denum = 0;
    double pval_nnls_pos, pval_nnls_neg;
    double tol = 1e-6, qr_tol = 1e-7, nnls_tol = 1e-10; // qr threshold used in R
    double nl_dbl_dmin = 10.0 * std::numeric_limits<double>::min();
    std::string burden_type, burden_str, burden_model;
    std::string out_file_prefix; // prefix of output files
    Eigen::ArrayXi colKeep; // keep track of linearly independent columns in Gmat
    Eigen::MatrixXd Gtmp;
    ArrayXb good_vars;
    Eigen::ArrayXd log10pv; 
    std::vector<int> indices_vars;
    std::vector<std::string> variant_names;
    bool apply_single_p = false;

    // for prep.
    bool get_test_info(const struct param*,const std::string&,mstream&);
    bool set_vars(const int&,const int&,std::vector<variant_block> const&);
    void compute_qr_G(const Eigen::Ref<const MatrixXb>&,struct geno_block const*);

    // assoc. tests
    std::vector<std::string> apply_joint_test(const int&,const int&,struct phenodt const*,const Eigen::Ref<const Eigen::MatrixXd>&,struct geno_block const*,std::vector<variant_block>&,std::vector<std::string> const&,struct param const*);
    void compute_minp();
    void compute_acat(const int&,const int&,const std::vector<variant_block>&);
    void compute_ftest(const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&);
    void compute_nnls(const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,bool const&);
    void compute_gates(const int&,const std::vector<variant_block>&);
    double get_me(const Eigen::Ref<const Eigen::MatrixXd>&);

    // final acat round
    void check_class_genep(std::string const&,std::map<std::string,bool> const&);
    void add_class(std::string const&,std::vector<std::string> const&,std::map<std::string,bool> const&);
    void run_single_p_acat(int const&,const int&,const int&,int const&,const std::string&,std::vector<variant_block>&,std::map<std::string, double>&,struct geno_block const*,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const MatrixXb>&,std::vector<std::vector<std::string>>&,struct param const*);
    std::string print_gene_output(const std::string&,const std::string&,const int&,const int&,const int&,const std::string&,struct param const*);
    std::string print_sum_stats_gene(const std::string&,const std::string&,const int&,const int&,const int&,struct param const*);
    std::string print_sum_stats_htp_gene(const std::string&,const std::string&,const int&,const int&,const std::string&,struct param const*);

    // print results
    std::string print_output(const int&,const int&,const int&,const int&,const std::string&,struct param const*);
    std::string print_output(const int&,const std::string&,const int&,const int&,const int&,const std::string&,struct param const*);
    std::string print_sum_stats(const std::string&,const int&,const int&,const int&,struct param const*);
    std::string print_sum_stats_htp(const std::string&,const int&,const int&,const std::string&,struct param const*);

    void get_variant_names(int const&,int const&,std::vector<snp> const&);
    void reset_vals();
    void get_pv(const double&);

    JTests();
    ~JTests();
};

double get_acat(const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&);
double get_acat(const Eigen::Ref<const Eigen::ArrayXd>&); // uniform weights
bool valid_pval(double const&);
void map_to_vec(int&,std::map<std::string,double>&,Eigen::ArrayXd&);

#endif
