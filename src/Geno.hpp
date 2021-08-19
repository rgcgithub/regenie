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

#ifndef GENO_H
#define GENO_H

#include "bgen_to_vcf.hpp"
#include "pgenlibr.h"

#define BIT_SET(a,b) ((a) |= (1ULL<<(b)))
#define BIT_UNSET(a,b) ((a) &= ~(1ULL << (b)))
#define CHECK_BIT(a,b) ((a) & (1ULL<<(b)))

struct annoinfo {
  uint16_t regionid = 0ULL;
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
  bool MAC_fail_if_checked = true; // for extract/exclude OR
} ;

struct tally {
  uint32_t snp_count = 0;
  uint32_t n_failed_tests = 0;
  uint32_t n_ignored_snps = 0;
  uint32_t n_ignored_tests = 0;
};

// for step 2 per thread
struct data_thread {
  SpVec Gsparse;
  Eigen::ArrayXd scale_fac_pheno;
  Eigen::MatrixXd Gres;
  Eigen::ArrayXd Gmod;
  Eigen::ArrayXd denum;
  Eigen::ArrayXd stats;
  Eigen::ArrayXd chisq_val;
  Eigen::ArrayXd pval_log;
  Eigen::ArrayXd bhat;
  Eigen::ArrayXd se_b;
  // for spa
  bool pos_score;
  double val_a, val_b, val_c, val_d; 
  // firth
  double dif_deviance;
  Eigen::MatrixXd beta_null_firth;
  // reset each time
  bool fastSPA = true;
  bool is_sparse = false;
};

struct geno_block {
  uint32_t ns, nv;
  BgenParser bgen;
  PgenReader pgr;
  Eigen::MatrixXd Gmat;
  Eigen::MatrixXd snp_afs;
  std::vector<data_thread> thread_data;
};

struct variant_block {
  bool ignored, flipped;
  double scale_fac, af1, info1, ns1;
  Eigen::ArrayXi ns, ns_case, nmales;
  Eigen::ArrayXd af, af_case, af_control, mac, info;
  Eigen::MatrixXd genocounts;
  ArrayXb ignored_trait;
  ArrayXb test_fail;
  ArrayXb is_corrected; // for firth/spa
  // for masks
  bool singleton;
  // for joint test
  Eigen::ArrayXd pval_log;
  // interaction test
  bool skip_int, fitHLM;
  ArrayXb is_corrected_inter; // for firth
  ArrayXb test_fail_inter; // for firth
  // association test info
  std::vector<std::string> sum_stats;
};

// for conditional analyses
struct cond_geno_info {
  bool dosage_mode, zlib_compress, streamBGEN;
  PgenReader pgr;
  ArrayXb sample_keep; // keep track of samples in analysis
  Eigen::ArrayXi sample_index; // index of samples in analysis
};

struct findID {
  uint32_t index;
  bool is_found;
};



void check_bgen(const std::string&,std::string const&,bool&,bool&,uint&,int const&);
void prep_bgen(struct in_files*,struct param*,struct filter*,std::vector<snp>&,std::map<int,std::vector<int>>&,BgenParser&,mstream&);
void read_bgen_sample(const std::string&,struct param*,std::vector<std::string> &,mstream&);
void read_bgen_sample(const std::string&,std::vector<std::string> &,mstream&);
void read_bgi_file(BgenParser&,struct in_files*,struct param*,struct filter*,std::vector<snp>&,mstream&);
void read_bgi_file(BgenParser&,struct in_files*,struct param*,struct filter*,mstream&);

void read_bed_bim_fam(struct in_files*,struct param*,struct filter*,std::vector<snp>&,std::map<int,std::vector<int>>&,mstream&);
void read_bim(struct in_files*,struct param*,struct filter*,std::vector<snp>&,mstream&);
void read_fam(struct in_files*,struct param*,mstream&);
void prep_bed(const uint32_t&, struct in_files*,mstream&);

void read_pgen_pvar_psam(struct in_files*,struct param*,struct filter*,struct geno_block*,std::vector<snp>&,std::map<int,std::vector<int>>&,mstream&);
uint64 read_pvar(struct in_files*,struct param*,struct filter*,std::vector<snp>&,mstream&);
void read_psam(struct in_files*,struct param*,mstream&);
void prep_pgen(struct in_files const*,struct filter const*,struct geno_block*,struct param*,mstream&);

ArrayXb check_in_map_from_files(std::map<std::string,uint>&,std::vector<std::string> const&,mstream&);
ArrayXb check_in_map_from_files_IDs(std::vector<std::string> const&,struct param*,mstream&);
void check_snps_include_exclude(struct in_files*,struct param*,struct filter*,std::vector<snp>&,std::map<int,std::vector<int>>&,mstream&);
void check_snps_include_exclude_or(struct in_files*,struct param*,struct filter*,std::vector<snp>&,mstream&);
void check_samples_include_exclude(struct in_files const*,struct param*,struct filter*,mstream&);

void get_G(const int&,const int&,const int&,const uint32_t&,std::vector<snp> const&,struct param const*,struct in_files*,struct geno_block*,struct filter const*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,mstream&);

void readChunkFromBGENFileToG(const int&,const int&,const uint32_t&,std::vector<snp> const&,struct param const*,struct geno_block*,struct filter const*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,mstream&);
void readChunkFromBGENFileToG_fast(const int&,const int&,const uint32_t&,std::vector<snp> const&,struct param const*,struct in_files*,struct geno_block*,struct filter const*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,mstream&);
void readChunkFromBedFileToG(const int&,const int&,const uint32_t&,std::vector<snp> const&,struct param const*,struct in_files*,struct geno_block*,struct filter const*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,mstream&);
void readChunkFromPGENFileToG(const int&,const uint32_t&,std::vector<snp> const&,struct param const*,struct geno_block*,struct filter const*,const Eigen::Ref<const MatrixXb>&,mstream&);

void readChunkFromBGENFileToG(std::vector<uint64> const&,const int&,std::vector<snp> const&,struct param const*,struct geno_block*,struct filter const*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,std::vector<variant_block>&,mstream&);
void readChunkFromBGEN(std::istream*,std::vector<uint32_t>&,std::vector<uint32_t>&,std::vector<std::vector<uchar>>&,std::vector<uint64>&);
void parseSNP(const int&,const int&,std::vector<uchar>*,const uint32_t&,const uint32_t&,struct param const*,struct filter const*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,const snp*,struct geno_block*,variant_block*,mstream&);
void parseSnpfromBGEN(const int&,const int&,std::vector<uchar>*,const uint32_t&,const uint32_t&,struct param const*,struct filter const*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,const snp*,struct geno_block*,variant_block*,mstream&);
void parseSnpfromBed(const int&,const int&,const std::vector<uchar>&,struct param const*,struct filter const*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,const snp*,struct geno_block*,variant_block*);
void readChunkFromPGENFileToG(std::vector<uint64> const&,const int&,struct param const*,struct filter const*,struct geno_block*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,std::vector<snp> const&,std::vector<variant_block>&);

void skip_snps(uint64 const&,struct param const*,struct in_files*,struct geno_block*);
void jumpto_bed(uint64 const&,uint64 const&,std::ifstream&);
void buildLookupTable(std::vector<Eigen::ArrayXd>&);
void prep_snp_stats(variant_block*,struct param const*);
void initialize_thread_data(std::vector<data_thread>&,struct param const&);
void reset_thread(data_thread*,struct param const&);
void reset_stats(variant_block*,struct param const&);
void update_trait_counts(int const&,double const&,double const&,int const&,double const&,variant_block*,const Eigen::Ref<const MatrixXb>&);
void update_genocounts(bool const&,int const&,int const&,Eigen::MatrixXd&,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&);
void update_af_cc(int const&,double const&,variant_block*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&);
void compute_mac(bool const&,double&,double const&,int const&,bool const&,variant_block*,struct param const*);
void compute_aaf_info(double&,double const&,variant_block*,struct param const*);
void flip_geno(double&,Eigen::Ref<Eigen::ArrayXd>,variant_block*,struct param const*);
void check_sparse_G(int const&,int const&,struct geno_block*,uint32_t const&,const Eigen::Ref<const ArrayXb>&);
void mean_impute_g(double &,const double&,const bool&);
void mean_impute_g(const double&,Eigen::Ref<Eigen::ArrayXd>,const Eigen::Ref<const ArrayXb>&);
void residualize_geno(int const&,int const&,variant_block*,bool const&,const Eigen::Ref<const Eigen::MatrixXd>&,struct geno_block*,struct param const*);
void writeSnplist(std::string const&,int const&,int const&,std::vector<snp> const&,mstream&);

bool in_chrList(const int&,struct filter const*);
std::string bgi_chrList(struct filter*,const int&);
std::string bgi_chrList(const int&,const int&);
std::string bgi_rsidList(std::map <std::string, uint64>&);
bool in_range(int const&,uint32_t const&,struct param const*);
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
  uint16_t region = 0ULL;
};

void read_setlist(const struct in_files*,struct param*,struct filter*,std::vector<std::vector<vset>>&,std::vector<snp>&,const uint64,const double,mstream&);
void check_sets_include_exclude(bool const&,const struct in_files*,struct param*,struct filter*,std::vector<std::vector<vset>>&,mstream&);
void check_in_map_from_files_sets(bool const&,std::map<std::string,std::vector<int>>&,std::vector<std::string> const&,bool const&,mstream&);

void get_masks_info(const struct in_files*,struct param*,struct filter*,std::map<std::string,anno_name>&,std::map <std::string, std::map <std::string,uint16_t>>&,std::vector<maskinfo>&,std::vector<std::vector<std::string>>&,uint64&,std::vector<snp>&,mstream& sout);
void read_anno_cat(const struct in_files*,struct param*,std::map<std::string,anno_name>&,mstream& sout);
void read_anno(struct param*,const struct in_files*,struct filter*,std::map<std::string,anno_name>&,std::map <std::string, std::map <std::string,uint16_t>>&,std::vector<snp>&,mstream& sout);
void read_aafs(const double,const struct in_files*,struct filter*,std::vector<snp>&,mstream& sout);
void read_masks(const struct in_files*,struct param*,std::map<std::string,anno_name>&,std::vector<maskinfo>&,std::vector<std::vector<std::string>>&,uint64&,mstream& sout);

void read_snp(uint64 const&,Eigen::Ref<Eigen::ArrayXd>,Eigen::Ref<ArrayXb>,struct filter*,struct in_files*,struct geno_block*,struct param*);
void read_snp_bed(uint64 const&,Eigen::Ref<Eigen::ArrayXd>,Eigen::Ref<ArrayXb>,struct filter*,struct in_files*,struct param*);
void read_snp_pgen(uint64 const&,Eigen::Ref<Eigen::ArrayXd>,Eigen::Ref<ArrayXb>,struct geno_block*,struct param*);
void read_snp_bgen(uint64 const&,Eigen::Ref<Eigen::ArrayXd>,Eigen::Ref<ArrayXb>,struct filter*,std::string const&,struct param*,int const&);
void code_snp(Eigen::MatrixXd&,Eigen::Ref<ArrayXb>,uint64 const&,struct filter*,struct in_files*,struct param*,mstream&);

// for conditional analyses
void get_conditional_vars(std::map<std::string,uint64>&,struct in_files*,struct param const*,mstream&);
void get_snps_offset(std::map<std::string,uint64>&,std::map<std::string,uint32_t>&,std::vector<snp> const&,mstream&);
void get_snps_offset(std::map<std::string,uint64>&,std::map<std::string,uint64>&,mstream&);
Eigen::MatrixXd extract_from_genofile(Eigen::Ref<ArrayXb>,struct filter*,struct in_files*,struct param*,mstream&);
void setup_bgen(struct cond_geno_info&,std::map<std::string,uint64>&,Eigen::Ref<ArrayXb>,struct in_files*,struct param*,struct filter*,mstream&);
void read_snps_bgen(std::map<std::string,uint64>&,Eigen::Ref<Eigen::MatrixXd>,struct cond_geno_info&,Eigen::Ref<ArrayXb>,std::string const&,struct param*);
void read_snps_bgen(std::map<std::string,uint64>&,Eigen::Ref<Eigen::MatrixXd>,struct cond_geno_info&,Eigen::Ref<ArrayXb>,std::string const&);
void setup_pgen(struct cond_geno_info&,std::map<std::string,uint64>&,Eigen::Ref<ArrayXb>,struct in_files*,struct param*,mstream&);
uint32_t read_pvar(std::map<std::string,uint64>&,struct in_files*);
uint32_t read_psam(struct cond_geno_info&,Eigen::Ref<ArrayXb>,struct in_files*,struct param*);
void prep_pgen(uint32_t&,uint32_t&,struct cond_geno_info&,struct in_files const*);
void read_snps_pgen(std::map<std::string,uint64>&,Eigen::Ref<Eigen::MatrixXd>,struct cond_geno_info&,Eigen::Ref<ArrayXb>);
void setup_bed(struct cond_geno_info&,std::map<std::string,uint64>&,Eigen::Ref<ArrayXb>,struct in_files*,struct param*,mstream&);
uint32_t read_bim(std::map<std::string,uint64>&,struct in_files*);
uint32_t read_fam(struct cond_geno_info&,Eigen::Ref<ArrayXb>,struct in_files*,struct param*);
void prep_bed(uint32_t&,struct cond_geno_info&,struct in_files const*);
void read_snps_bed(std::map<std::string,uint64>&,Eigen::Ref<Eigen::MatrixXd>,struct cond_geno_info&,Eigen::Ref<ArrayXb>,std::string const&,struct param*,mstream&);

#endif
