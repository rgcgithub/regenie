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

#include "NNLS.hpp"

#define MULTITRAIT_N_TESTS 14
// 0. MANOVA -- assoc_manova
// 1. Omnibus 0 --  assoc_omnibus0
// 2. Bayesian test (Baes Factor or BF) -- assoc_bayes
// 3. NNLS 0 = Non-negative least squares 0  --  assoc_nnls0
// 4. SumZ  --  assoc_omnibus0
// 5. NP-MANOVA 0 = Non-Parametric MANOVA 0 -- assoc_manova
// 6. Het. Omnibus 0 --  assoc_omnibus0
// 7. CPC 0  -- assoc_cpc0
// 8. Robust CPC 0 = RCPC0SUMCHI2 -- assoc_rcpc0
// 9. Robust CPC 0 = RCPC0FISHER -- assoc_rcpc0
// 10. Robust CPC 0 = RCPC0ACAT -- assoc_rcpc0
// 11. Adjusted CPC 0 = ACPC0SUMCHI2  --  assoc_cpc0
// 12. Adjusted CPC 0 = ACPC0FISHER   --  assoc_cpc0
// 13. Adjusted CPC 0 = ACPC0ACAT     --  assoc_cpc0
// * 0 = strict on missing

class MTestsResults 
{
  public:
    unsigned int n_tests; // numbre of multi-trait tests
    unsigned int q; // number of traits                
    unsigned int M; // number of variables (to test for association)
                      
    // test results
    std::vector<std::vector<double>> stats_mt; // dim 1 = mt; dim 2 = variants
    std::vector<std::vector<double>> logp_mt; // dim 1 = mt; dim 2 = variants
    std::vector<std::vector<vector<double>>> logp_univar; // dim 1 = mt; dim 2 = variants; dim 3 = traits
    std::vector<std::vector<double>> zscore_univar; // dim 1 = variants; dim 2 = traits
    std::vector<std::vector<double>> zscore_cpc; // dim 1 = variants; dim 2 = traits
    std::vector<std::vector<double>> zscore_acpc; // dim 1 = variants; dim 2 = traits
    std::vector<std::vector<double>> zscore_rcpc; // dim 1 = variants; dim 2 = traits
                                                
    //----------
    // Methods
    //----------
    void setup(unsigned int, unsigned int, unsigned int);

    MTestsResults(unsigned int, unsigned int, unsigned int);
    ~MTestsResults();
};

class MTests {

  public:
    int verbose; // verbose level 
    bool precomp; // precompute products like YtY to speed up computation

    // NNLS object
    NNLS nnls;

    // MCC for ACPC
    double mcc_skew_abs; // skewness thr. to apply MCC
    double mcc_z2; // Z-score squared thr. to apply MCC
    Eigen::ArrayXd skew_PC; // (absolute) skewness of Robust PCs
    unsigned int n_skewed_pc; // number of PCs with abs. skewness > thr.

    // information copied from Data class
    MatrixXb Mask0; // masked samples = samples with missing values on phenotypes; 0 = missing in any phenotype
    double Neff0; // number of non-missing samples (per trait)
    MatrixXb Mask; // masked samples = samples with missing values on phenotypes; 0 = missing in a given phenotype; column = trait
    Eigen::VectorXd Neff; // number of non-missing samples (per trait)

    Eigen::MatrixXd Yres; // matrix of residualized traits (per-trait missing patterns)
    Eigen::MatrixXd Y0res; // matrix of residualized traits (commont missing pattern across traits)
    
    // list of multi-trait tests
    unsigned int n_tests;
    std::map<std::string, uint16_t> mt_tests_map = { {"manova", 0}, {"omnibus0", 1}, {"bayes", 2} };
    std::vector<std::string> test_names = {"MANOVA", "Omnibus0", "Bayes"};
    uint16_t test_list;

    // dimensions
    unsigned int n_traits; // number of traits                
    /* unsigned int N; // number of samples */
    /* unsigned int q; // number of traits */                
    /* unsigned int M; // number of variables (to test for association) */
                    
    // parameters for Bayesian test
    double prior_a0;
    double prior_Q0;
    double prior_Mbeta0;
    double prior_Vbeta0;

    // precompute quantities
    Eigen::MatrixXd precomp0_YtY; // cross-product: Y^T Y
    Eigen::MatrixXd precomp0_Syy; // covariance of Y: cov(Y)
    Eigen::MatrixXd precomp0_Syy_inv; // inverse of covariance of Y
    // MANOVA-specific
    double precomp0_ld0;
    // Bayes-specific
    double precomp0_LL_M0;
    // Hiearachical Omnibus
    Eigen::VectorXd precomp0_lambdas_Syy, precomp0_lambdas_norm_Syy;
    // Robust Omnibus
    Eigen::MatrixXd PC_Y0res;
    Eigen::MatrixXd RPC_Y0res;

    // data for traits
    Eigen::MatrixXd Ryy;
    Eigen::MatrixXi N_Ryy;

    // test results
    /* VectorXd stats; */
    /* VectorXd pvals; */
    /* /1* vector<string> sum_stats; *1/ */
    /* vector<vector<double>> logp_mt; // 1 = mt; 2 = variants */
    /* vector<vector<vector<double>>> logp_univar; // 1 = mt; 2 = variants; 3 = traits */

    // parameters of association tests
    double nl_dbl_dmin = 10.0 * std::numeric_limits<double>::min();

    //----------
    // Methods
    //----------
    void setup();
    void setup_masks(const MatrixXb &);
    void setup_yres(const Eigen::MatrixXd &);
    void check_setup_no_data();
    void check_setup_data();

    // data prep.
    void compute_ncory(const Eigen::Ref<const MatrixXb>& M);
    void compute_cory(const Eigen::Ref<const Eigen::MatrixXd>& Y, const Eigen::Ref<const MatrixXb>& M);
    void compute_skew_pc();

    // assoc. tests
    void apply_tests(const int&,const int&,struct phenodt const*,const Eigen::Ref<const Eigen::MatrixXd>&,struct geno_block const*,std::vector<variant_block>&,std::vector<std::string> const&,struct param const*);
    void apply_tests_snp(int const&, struct geno_block& , const Eigen::Ref<const Eigen::MatrixXd>&, const Eigen::Ref<const Eigen::RowVectorXd>&, struct param const&);
    MTestsResults run_tests_snp(int const&, struct geno_block& , const Eigen::Ref<const Eigen::MatrixXd>&, const Eigen::Ref<const Eigen::RowVectorXd>&, struct param const&);
    MTestsResults run_tests_snp_precomp(int const&, struct geno_block& , struct param const&);
    
    void assoc_manova(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G, MTestsResults&);
    void assoc_omnibus0(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G, MTestsResults&);
    void assoc_cpc0(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G, MTestsResults&);
    void assoc_rcpc0(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G, MTestsResults&);
    void assoc_bayes(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G, MTestsResults&);
    void assoc_nnls0(const Eigen::MatrixXd &Y, const Eigen::MatrixXd& G, MTestsResults&);

    // print results 
    /* string print_multitrait_output(int const& snp_index, std::string const& test_string, std::string const& wgr_string, */ 
    /*     variant_block* block_info, std::vector<snp> const& snpinfo, */
    /*     struct in_files const& files,struct param const* params); */
    string print_sumstats(const MTestsResults&, int const& snp_index, uint32_t const&, string const& test_string, string const& wgr_string, variant_block* block_info, vector<snp> const& snpinfo, struct param const* params);

    // debug
    void dump_data(const Eigen::MatrixXd &, const Eigen::MatrixXd &, const MatrixXb &, const MatrixXb &);

    MTests();
    ~MTests();
};
#endif
