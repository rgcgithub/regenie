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

#ifndef GENO_H
#define GENO_H

#include "bgen_to_vcf.hpp"
#include "pgenlibr.h"

#define BIT_SET(a,b) ((a) |= (1ULL<<(b)))
#define CHECK_BIT(a,b) ((a) & (1ULL<<(b)))

struct annoinfo {
  uchar regionid = 0u;
  uint64 id = 0ULL;
};

struct snp {
  int chrom;
  std::string ID;
  //double genpos; 
  uint32_t physpos;
  std::string allele1, allele2;
  double MAF;
  uint64 offset;
  // for masks
  std::map <std::string, annoinfo> anno; // annotation
  float aaf = -1;
} ;

struct tally {
  uint32_t snp_count = 0;
  uint32_t n_failed_tests = 0;
  uint32_t n_ignored_snps = 0;
  uint32_t n_ignored_tests = 0;
  uint32_t n_skipped_snps = 0;
};

struct geno_block {
  BgenParser bgen;
  PgenReader pgr;
  std::vector<double> genobuf;
  Eigen::MatrixXd Gmat;
  Eigen::MatrixXd Gmat_tmp;
  std::vector < std::vector < uint > > non_zero_indices_G;
  std::vector < Eigen::MatrixXd > genocounts;
  Eigen::MatrixXd snp_afs;
  Eigen::MatrixXd snp_mac;
  Eigen::MatrixXd snp_info;
  ArrayXb bad_snps;
  std::vector<bool> snp_flipped;
};

struct variant_block {
  double scale_fac, af1, info1;
  Eigen::ArrayXi ns, nmales;
  Eigen::ArrayXd af, mac, info;
  Eigen::ArrayXd Gmod;
  Eigen::MatrixXd genocounts;
  std::vector<bool> test_fail;
  std::vector<bool> is_corrected; // for firth/spa
  Eigen::ArrayXd scale_fac_pheno;
  Eigen::ArrayXd denum;
  Eigen::ArrayXd stats;
  Eigen::ArrayXd chisq_val;
  Eigen::ArrayXd pval_log;
  Eigen::ArrayXd bhat;
  Eigen::ArrayXd se_b;
  // for spa
  bool flipped, pos_score;
  double val_a, val_b, val_c, val_d; 
  std::vector <uint> non_zero_indices;
  // firth
  double dif_deviance;
  Eigen::MatrixXd beta_null_firth;
  // reset each time
  bool ignored = false;
  ArrayXb ignored_trait;
  bool fastSPA = true;
  uint n_non_zero = 0;
  // for masks
  bool singleton;
};


struct findID {
  uint32_t index;
  bool is_found;
};



void check_bgen(const std::string,struct param*);
void prep_bgen(struct in_files*,struct param*,struct filter*,std::vector<snp>&,std::map<int,std::vector<int>>&,BgenParser&,mstream&);
void read_bgen_sample(const std::string,struct param*,std::vector<std::string> &,mstream&);
void read_bgi_file(BgenParser&,struct in_files*,struct param*,struct filter*,std::vector<snp>&,mstream&);

void read_bed_bim_fam(struct in_files*,struct param*,struct filter*,std::vector<snp>&,std::map<int,std::vector<int>>&,mstream&);
void read_bim(struct in_files*,struct param*,struct filter*,std::vector<snp>&,mstream&);
void read_fam(struct in_files*,struct param*,mstream&);
void prep_bed(const uint32_t&, struct in_files*,mstream&);

void read_pgen_pvar_psam(struct in_files*,struct param*,struct filter*,struct geno_block*,std::vector<snp>&,std::map<int,std::vector<int>>&,mstream&);
uint64 read_pvar(struct in_files*,struct param*,struct filter*,std::vector<snp>&,mstream&);
void read_psam(struct in_files*,struct param*,mstream&);
void prep_pgen(const uint32_t,const uint32_t,struct in_files*,struct filter*,struct geno_block*,mstream&);

void check_snps_include_exclude(struct in_files*,struct param*,struct filter*,std::vector<snp>&,std::map<int,std::vector<int>>&,mstream&);
void set_snps_to_keep(struct in_files*,struct param*,struct filter*,std::vector<snp>&,mstream&);
void set_snps_to_rm(struct in_files*,struct param*,struct filter*,std::vector<snp>&,mstream&);

void check_samples_include_exclude(struct in_files*,struct param*,struct filter*,mstream&);
void set_IDs_to_rm(struct in_files*,struct filter*,struct param*,mstream&);
void set_IDs_to_keep(struct in_files*,struct filter*,struct param*,mstream&);

void get_G(const int,const int,const int,const uint32_t,std::vector<snp>&,struct param*,struct in_files*,struct geno_block*,struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,mstream&);

void readChunkFromBGENFileToG(const int,const int,const uint32_t,std::vector<snp>&,struct param*,struct geno_block*,struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,mstream&);
void readChunkFromBedFileToG(const int,const int,const uint32_t,std::vector<snp>&,struct param*,struct in_files*,struct geno_block*,struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,mstream&);
void readChunkFromPGENFileToG(const int,const uint32_t,std::vector<snp>&,struct param*,struct geno_block*,struct filter*,const Eigen::Ref<const MatrixXb>&,mstream&);

void readChunkFromBGENFileToG(const int,const int,const uint32_t,std::vector<snp>&,struct param*,struct geno_block*,struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,std::vector<variant_block>&,mstream&);
void readChunkFromBGEN(std::istream*,uint32_t*,uint32_t*,std::vector<uchar>*);
void parseSnpfromBGEN(const int,const int&,std::vector<uchar>*,const uint32_t,const uint32_t,struct param*,const struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,const snp*,struct geno_block*,variant_block*,mstream&);
void parseSnpfromBed(const int,const int&,const std::vector<uchar>,struct param*,const struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,struct geno_block*,variant_block*);
void readChunkFromPGENFileToG(const int&,const int&,const int&,struct param*,struct filter*,struct geno_block*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,const std::vector<snp>&,std::vector<variant_block>&);

void skip_snps(uint64,struct param*,struct in_files*,struct geno_block*);
void jumpto_bed(uint64,struct in_files*);
void prep_snp_stats(variant_block*,struct param*);
void update_trait_counts(int,double,double,int,double,variant_block*,const Eigen::Ref<const MatrixXb>&);
void update_genocounts(bool,int,int,Eigen::MatrixXd&,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&);
void compute_mac(bool,double&,double,int,int,variant_block*,struct param*);
void compute_aaf_info(double&,int,double,variant_block*,struct param*);
void update_nnz_spa(uint32_t,uint32_t,variant_block*);
void mean_impute_g(double &,const double,const bool,const bool,const bool);

bool in_chrList(const int,struct filter*);
std::string bgi_chrList(struct filter*);
bool in_range(int,uint32_t,struct param*);
findID getIndivIndex(const std::string&,const std::string&,struct param*,mstream&);


// for snp-set methods
struct vset {
  int chrom = 0;
  uint32_t physpos;
  std::string ID;
  std::vector<uint64> snp_indices;
} ;
struct anno_name {
  std::string name;
  uint64 id;
};
struct maskinfo {
  std::string name, region_name = "";
  uint64 id = 0ULL;
  uchar region = 0u;
};

void read_setlist(const struct in_files*,struct param*,struct filter*,std::vector<std::vector<vset>>&,std::vector<snp>&,const uint64,const double,mstream&);
void check_sets_include_exclude(bool,const struct in_files*,struct param*,struct filter*,std::vector<std::vector<vset>>&,mstream&);
void set_sets_to_keep(int&,const struct in_files*,struct param*,struct filter*,mstream&);
void set_sets_to_rm(int&,const struct in_files*,struct param*,struct filter*,mstream&);

void get_masks_info(const struct in_files*,struct param*,struct filter*,std::map<std::string,anno_name>&,std::vector<maskinfo>&,std::vector<std::vector<std::string>>&,uint64&,std::vector<snp>&,mstream& sout);
void read_anno_cat(const struct in_files*,struct param*,std::map<std::string,anno_name>&,mstream& sout);
void read_anno(struct param*,const struct in_files*,struct filter*,std::map<std::string,anno_name>&,std::map <std::string, int>&,std::vector<snp>&,mstream& sout);
void read_aafs(const double,const struct in_files*,struct filter*,std::vector<snp>&,mstream& sout);
void read_masks(const struct in_files*,struct param*,std::map<std::string,anno_name>&,std::map <std::string, int>,std::vector<maskinfo>&,std::vector<std::vector<std::string>>&,uint64&,mstream& sout);

void readChunkFromBGENFileToG(const int,const int,const uint32_t,std::vector<uint64>&,std::vector<snp>&,struct param*,struct geno_block*,struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,std::vector<variant_block>&);
void readChunkFromPGENFileToG(const int,const int,std::vector<uint64>&,const int&,struct param*,struct filter*,struct geno_block*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,const std::vector<snp>&,std::vector<variant_block>&);

#endif
