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

#include "pgenlibr.h"

struct snp {
  int chrom;
  std::string ID;
  double genpos; 
  uint32_t physpos;
  std::string allele1, allele2;
  double MAF;
  uint64 offset;
  bool mask = false;
} ;

struct tally {
  uint32_t snp_count = 0;
  uint32_t n_failed_tests = 0;
  uint32_t n_ignored_snps = 0;
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
  Eigen::MatrixXd snp_info;
  Eigen::ArrayXd bad_snps;
  std::vector<bool> snp_flipped;
};

struct variant_block {
  double af, info, scale_fac;
  double dif_deviance;
  double val_a, val_b, val_c, val_d; 
  Eigen::ArrayXd Geno;
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
  std::vector <uint> non_zero_indices;
  bool flipped, pos_score;
  // firth
  Eigen::MatrixXd beta_null_firth;
  // reset each time
  bool ignored = false;
  bool fastSPA = true;
  uint n_non_zero = 0;
};


struct findID {
  uint32_t index;
  bool is_found;
};



void check_bgen(const std::string,struct param*);
void prep_bgen(struct in_files*,struct param*,struct filter*,std::vector<snp>&,std::map<int,std::vector<int>>&,BgenParser&,mstream&);
void read_bgen_sample(const std::string,const int,std::vector<std::string> &,mstream&);

void read_bed_bim_fam(struct in_files*,struct param*,struct filter*,std::vector<snp>&,std::map<int,std::vector<int>>&,mstream&);
void read_bim(struct in_files*,struct param*,struct filter*,std::vector<snp>&,std::map<int,std::vector<int>>&,mstream&);
void read_fam(struct in_files*,struct param*,mstream&);
void prep_bed(struct in_files*,struct param*,mstream&);

void read_pgen_pvar_psam(struct in_files*,struct param*,struct filter*,struct geno_block*,std::vector<snp>&,std::map<int,std::vector<int>>&,mstream&);
void read_pvar(uint32_t&,struct in_files*,struct param*,struct filter*,std::vector<snp>&,std::map<int,std::vector<int>>&,mstream&);
void read_psam(struct in_files*,struct param*,mstream&);
void prep_pgen(uint32_t&,struct in_files*,struct param*,struct geno_block*,mstream&);

void set_snps_to_keep(struct in_files*,struct filter*,mstream&);
void set_snps_to_rm(struct in_files*,struct filter*,mstream&);
void set_IDs_to_rm(struct in_files*,struct filter*,struct param*,mstream&);
void set_IDs_to_keep(struct in_files*,struct filter*,struct param*,mstream&);

void get_G(const int,const int,const int,uint32_t&,std::vector<snp>&,struct param*,struct in_files*,struct geno_block*,struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,mstream&);

void readChunkFromBGENFileToG(const int,const int,uint32_t&,std::vector<snp>&,struct param*,struct geno_block*,struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,mstream&);
void readChunkFromBedFileToG(const int,uint32_t&,std::vector<snp>&,struct param*,struct in_files*,struct geno_block*,struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,mstream&);
void readChunkFromPGENFileToG(const int,uint32_t&,std::vector<snp>&,struct param*,struct geno_block*,struct filter*,const Eigen::Ref<const MatrixXb>&,mstream&);

void readChunkFromBGEN(std::istream*,uint32_t*,uint32_t*,std::vector<uchar>*);
void parseSnpfromBGEN(std::vector<uchar>*,const uint32_t,const uint32_t,const struct param*,const struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,const snp*,variant_block*,mstream&);
void parseSnpfromBed(const std::vector<uchar>,const struct param*,const struct filter*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,variant_block*);
void readChunkFromPGENFileToG(const int&,const int&,struct param*,struct filter*,struct geno_block*,const Eigen::Ref<const MatrixXb>&,const Eigen::Ref<const Eigen::MatrixXd>&,std::vector<variant_block>&);

void skip_snps(const int&,struct param*,struct in_files*,struct geno_block*);


findID getIndivIndex(const std::string&,const std::string&,struct param*,mstream&);


#endif
