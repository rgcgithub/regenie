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

#ifndef MASK_H
#define MASK_H


class GenoMask {

  public:
    std::map <std::string, anno_name> annotations; // store identifier as 1 byte vector
    std::map <std::string, std::map <std::string, uint64>> regions; // store identifier as 1 byte vector
    std::vector <maskinfo> masks, base_masks;
    std::vector <std::vector <std::string>> mask_out, list_masks;//contains mask info
    Eigen::VectorXd aafs;
    Eigen::ArrayXi nsites;
    ArrayXb colset;
    MatrixXb keepaaf, keepmask, non_missing;
    Eigen::MatrixXd Gtmp; // holds mask
    std::ofstream outfile_bed, outfile_bim;
    std::vector<std::shared_ptr<Files>> setfiles;// for written setlist
    std::vector<std::vector<int>> setfiles_index;// for written setlist
    Files snplist_out;
    std::vector <std::vector <std::string>> list_snps;//contains snplist

    double tol = 1e-6;
    double minAAF = 1e-7, default_aaf = .01;
    int n_aaf_bins, max_aaf_bins = 12, nmasks_total;
    int n_mask_pass = 0; // number of masks generated
    bool write_setlist = false, write_masks = false, write_snplist = false, verbose = false;
    bool w_regions = false, w_loo = false, w_lodo = false, w_vc_tests = false, w_vc_cust_weights = false;
    bool take_max = true, take_comphet = false; // either max comphet or sum
    bool force_singleton = false; // allow user to specify singleton variants
    std::string gfile_prefix;
    uint64 gblock_size; // number of bytes to use for bed file format
    double max_aaf = -1, vc_aaf, vc_collapse_MAC; // maximum AAF to consider
    uint64 all_masks = 0ULL; // keep track of all annotations considered in analysis
    std::vector<std::vector<uchar>> gvec;
    uchar last_byte_correction_factor = 0u;

    bool remeta_save_ld = false;
    std::vector<std::string> remeta_snplist; // list of snps contained in any mask

    // functions
    void prep_run(struct param&,struct in_files const&);
    void setBins(struct param*,mstream&);
    void prepMasks(const int&,const std::string&);
    void set_snp_masks(const int&,const int&,std::vector<variant_block> const &,vset const&,std::vector<snp>&,mstream&);
    void set_snp_aafs(const int&,const int&,const bool&,std::vector<variant_block> const&,vset&,std::vector<snp> const&,mstream&);
    bool check_in_lovo_mask(const Eigen::Ref<const Eigen::ArrayXd>&,struct filter const&,std::string const&,snp&,bool&,bool&,double&,int const&,struct param const*);
    void updateMasks(const int&,const int&,struct param*,struct filter*,const Eigen::Ref<const MatrixXb>&,struct geno_block*,const Eigen::Ref<const Eigen::ArrayXd>&,std::vector<variant_block>&,vset&,std::vector<snp>&,mstream&);
    void apply_rule(SpVec&,SpVec const&,const Eigen::Ref<const ArrayXb>&,bool const&);
    void apply_rule(Eigen::Ref<Eigen::ArrayXd>,SpVec const&,const Eigen::Ref<const ArrayXb>&,bool const&);
    void apply_rule(Eigen::Ref<Eigen::ArrayXd>,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const ArrayXb>&,bool const&);
    void collapse_mask_chunk(const Eigen::Ref<const Eigen::ArrayXi>&,SpMat const&,const Eigen::Ref<const ArrayXb>&,const Eigen::Ref<const ArrayXb>&,const Eigen::Ref<const Eigen::ArrayXd>&,Eigen::Ref<Eigen::ArrayXd>,Eigen::Ref<Eigen::ArrayXd>,const Eigen::Ref<const ArrayXb>&);
    void updateMasks_loo(const Eigen::Ref<const Eigen::ArrayXi>&,bool const&,SpMat const&,const Eigen::Ref<const ArrayXb>&,const Eigen::Ref<const ArrayXb>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const Eigen::ArrayXd>&,const Eigen::Ref<const ArrayXb>&,vset&,int const&);
    void tally_masks(struct param const*,struct filter const*,const Eigen::Ref<const MatrixXb>&,SpMat&,MatrixXb&);
    void computeMasks(struct param*,struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,struct geno_block*,std::vector<variant_block>&,vset&,std::vector<snp>&,mstream&);
    void computeMasks_loo(const Eigen::Ref<const Eigen::ArrayXi>&,bool const&,struct param*,struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,struct geno_block*,std::vector<variant_block>&,vset&,std::vector<snp>&,mstream&);
    void buildMask(const int&,const int&,struct param const*,struct filter const*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,variant_block*);

    void get_mafs(const int&,Eigen::ArrayXd&,std::vector<variant_block> const&);

    void write_info(struct param*,struct filter const*,mstream&);
    void write_famfile(struct param*,struct filter const*,mstream&);
    void reset_gvec();
    void make_genovec(const int&,Eigen::Ref<const Eigen::ArrayXd>,struct filter const*);
    void write_genovec(const int&);
    void set_gvalue(const int&,const int&,const int&,const int&);
    void write_genobim(struct snp const&);
    void setAllBitsZero(const int&);
    void setAllBitsOne(const int&);
    std::string build_header();
    void build_map(std::map<std::string,std::vector<int>>&);
    void prep_snplist(const std::string&,mstream& sout);
    void append_snplist(int const&,ArrayXb const&, int const&,vset const&,std::vector<snp> const&);
    void make_snplist(int const&,std::string const&);
    void prep_setlists(const std::string&,const std::string&,mstream& sout);
    void append_setlist(int const&,std::string const&);
    void make_setlist(std::string const&,int const&,uint32_t const&);
    void closeFiles();


    GenoMask();
    ~GenoMask();

};

Eigen::ArrayXi check_lovo_snplist(const Eigen::Ref<const Eigen::ArrayXi>&,std::vector<uint64> const&,std::vector<snp> const&,std::string const&);
Eigen::ArrayXi get_index_vec_loo(int const&,int const&);

#endif
