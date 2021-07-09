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

#ifndef JOINT_H
#define JOINT_H


class JTests {

  public:
    std::vector<std::string> test_names = {"MINP","F","GATES","NNLS","ACAT","NNLS_POS","NNLS_NEG"};
    uchar qr_tests = 7<<1; // 01110000 for ftest, gates and nnls

    // store variant set info (for each chr)
    std::vector<std::vector<vset>> setinfo;

    // for testing
    uchar test_list; // 8 tests max
    int df_test, nvars, ncovars; // df, number of variants which passed filters,#covariates
    bool nnls_verbose_out = false;
    bool nnls_normalize = true, nnls_strict = false;
    int nnls_napprox, nnls_verbose = 0;
    double acat_a1,acat_a2, acat_pv_thr = 1e-15;
    bool valid_snp_mode;
    double pval, plog, zval, scale_denum = 0;
    double pval_nnls_pos, pval_nnls_neg;
    double tol = 1e-6, qr_tol = 1e-7, nnls_tol = 1e-10; // qr threshold used in R
    double nl_dbl_dmin = 10.0 * std::numeric_limits<double>::min();
    std::string burden_type, burden_str, burden_model;
    std::string out_file_prefix; // prefix of output files
    Eigen::ArrayXi colKeep; // keep track of linearly independent columns in Gmat
    Eigen::MatrixXd Gtmp;
    std::vector<bool> good_vars;
    std::vector<int> indices_vars;
    std::vector<std::string> variant_names;

    // for prep.
    bool get_test_info(const struct param*,const std::string&,mstream&);
    bool set_vars(const int&,const int&,std::vector<variant_block> const&);
    void compute_qr_G(const Eigen::Ref<const MatrixXb>&,struct geno_block const*);

    // assoc. tests
    std::string apply_joint_test(const int&,const int&,const int&,struct phenodt const*,const Eigen::Ref<const Eigen::MatrixXd>&,struct geno_block const*,std::vector<variant_block>&,const std::string&,struct param const*);
    void compute_minp(const int&,const int&,const std::vector<variant_block>&);
    void compute_acat(const int&,const int&,const std::vector<variant_block>&);
    void compute_ftest(const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&);
    void compute_nnls(const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&);
    void compute_gates(const int&,const std::vector<variant_block>&);
    double get_me(const Eigen::Ref<const Eigen::MatrixXd>&);

    // print results
    std::string print_output(const int&,const int&,const int&,const int&,const std::string&,struct param const*);
    std::string print_sum_stats(const int&,const int&,const int&,const int&,struct param const*);
    std::string print_sum_stats_htp(const int&,const int&,const int&,const std::string&,struct param const*);

    void get_variant_names(int const&,int const&,std::vector<snp> const&);
    void reset_vals();
    void get_pv(const double&);
    int test_ind(const std::string&);

    JTests();
    ~JTests();
};

#endif
