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

#ifndef MASK_H
#define MASK_H


class GenoMask {

  public:
    std::map <std::string, anno_name> annotations; // store identifier as 1 byte vector
    std::map <std::string, std::map <std::string, uchar>> regions; // store identifier as 1 byte vector
    std::vector <maskinfo> masks, base_masks;
    std::vector <std::vector <std::string>> mask_out, list_masks;//contains mask info
    std::vector<double> aafs;
    Eigen::ArrayXi nsites;
    ArrayXb colset;
    MatrixXb keepaaf, keepmask, non_missing;
    Eigen::MatrixXd Gtmp; // holds mask
    std::ofstream outfile_bed, outfile_bim;
    std::vector<Files*> setfiles;// for written setlist
    std::vector<std::vector<int>> setfiles_index;// for written setlist

    double tol = 1e-6;
    double minAAF = 1e-7, default_aaf = .01;
    int n_aaf_bins, max_aaf_bins = 12, nmasks_total;
    uint32_t n_mask_pass = 0; // number of masks generated
    bool write_setlist = false, w_regions = false, w_loo = false;
    bool take_max = true, take_comphet = false; // either max comphet or sum
    std::string gfile_prefix;
    uint64 gblock_size; // number of bytes to use for bed file format
    double max_aaf = -1; // maximum AAF to consider
    uint64 all_masks = 0u; // keep track of all annotations considered in analysis
    std::vector<std::vector<uchar>> gvec;

    // functions
    void setBins(struct param*,mstream&);
    void prepMasks(const int,const std::string&);
    void set_snp_masks(const int,const int,std::vector<variant_block> &,vset&,std::vector<snp>&,mstream&);
    void set_snp_aafs(const int,const int,const bool,std::vector<variant_block> &,vset&,std::vector<snp>&,mstream&);
    void updateMasks(const int,const int,struct param*,struct filter*,const Eigen::Ref<const MatrixXb>&,struct geno_block*,std::vector<variant_block>&,vset&,std::vector<snp>&,mstream&);
    void updateMasks_loo(const int,const int,struct param*,struct filter*,const Eigen::Ref<const MatrixXb>&,struct geno_block*,std::vector<variant_block>&,vset&,std::vector<snp>&,mstream&);
    void tally_masks(struct param*,struct filter*,const Eigen::Ref<const MatrixXb>&);
    void computeMasks(struct param*,struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,struct geno_block*,std::vector<variant_block>&,vset&,std::vector<snp>&,mstream&);
    void computeMasks_loo(struct param*,struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,struct geno_block*,std::vector<variant_block>&,vset&,std::vector<snp>&,mstream&);
    void buildMask(const int,const int,struct param*,struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,variant_block*);

    void get_mafs(const int,Eigen::ArrayXd&,std::vector<variant_block>&);

    void write_info(struct param*,struct filter*,mstream&);
    void write_famfile(struct param*,struct filter*,mstream&);
    void make_genovec(const int,Eigen::Ref<const Eigen::ArrayXd>,struct filter*);
    void write_genovec(const int);
    void set_gvalue(const int,const int,const int,const int);
    void write_genobim(const struct snp);
    void setBitsZero(const int);
    void build_map(std::map<std::string,std::vector<int>>&);
    void prep_setlists(const std::string&,const std::string&,mstream& sout);
    void append_setlist(int,std::string);
    void make_setlist(std::string,int,uint32_t);
    void closeFiles();


    GenoMask();
    ~GenoMask();

};

#endif
