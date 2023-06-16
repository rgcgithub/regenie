/* 

   This file is part of the regenie software package.

   Copyright (c) 2020-2023 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

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

#ifndef MULTITRAIT_H
#define MULTITRAIT_H

using namespace std;
using namespace Eigen;

#define MULTITRAIT_N_TESTS 3

class MTests {

  public:
    // information copied from Data class
    Eigen::ArrayXd Neff; // number of non-missing samples (per trait)
    
    // list of multi-trait tests
    unsigned int n_tests;
    std::map<std::string, uint16_t> mt_tests_map = { {"omnibus", 0} };
    std::vector<std::string> test_names = {"Omnibus"};
    uint16_t test_list;

    // dimensions
    unsigned int N; // number of samples
    unsigned int q; // number of traits                
    unsigned int M; // number of variables (to test for association)
                    
    // parameters for Bayesian test
    double prior_a0;
    double prior_Q0;
    double prior_Mbeta0;
    double prior_Vbeta0;

    // data for traits
    MatrixXd Ryy;
    MatrixXi N_Ryy;

    // test results
    VectorXd stats;
    VectorXd pvals;
    /* vector<string> sum_stats; */
    vector<vector<double>> logp_tests;

    // parameters of association tests
    double nl_dbl_dmin = 10.0 * std::numeric_limits<double>::min();

    //----------
    // Methods
    //----------
    void setup();
    void check_setup();

    // data prep.
    void compute_ncory(const Ref<const MatrixXb>& M);
    void compute_cory(const Ref<const MatrixXd>& Y, const Ref<const MatrixXb>& M);

    // assoc. tests
    void apply_tests(const int&,const int&,struct phenodt const*,const Eigen::Ref<const Eigen::MatrixXd>&,struct geno_block const*,std::vector<variant_block>&,std::vector<std::string> const&,struct param const*);
    void apply_tests_snp(int const&, struct geno_block& , const Ref<const MatrixXd>&, const Ref<const RowVectorXd>&, struct param const&);
    
    void assoc_manova(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G);
    void assoc_omnibus0(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G);
    void assoc_bayes(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G);

    // print results 
    /* string print_multitrait_output(int const& snp_index, std::string const& test_string, std::string const& wgr_string, */ 
    /*     variant_block* block_info, std::vector<snp> const& snpinfo, */
    /*     struct in_files const& files,struct param const* params); */
    string print_sumstats( int const& snp_index, uint32_t const&, string const& test_string, string const& wgr_string, variant_block* block_info, vector<snp> const& snpinfo, struct param const* params);

    MTests();
    ~MTests();
};


#endif
