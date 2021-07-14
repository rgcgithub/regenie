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

#ifndef PHENO_H
#define PHENO_H

struct rank_pair {
  double val;
  size_t index;
};

struct phenodt {

  Eigen::MatrixXd new_cov;
  Eigen::MatrixXd interaction_cov, interaction_cov_res;
  Eigen::ArrayXd scl_inter_X;
  std::vector<Eigen::MatrixXd> Hmat;
  std::vector<Eigen::ArrayXd> scf_i;
  Eigen::MatrixXd phenotypes;
  Eigen::MatrixXd phenotypes_raw;
  MatrixXb masked_indivs;
  Eigen::ArrayXd Neff; // number of non-missing samples (per trait)
  Eigen::RowVectorXd scale_Y;

};


void read_pheno_and_cov(struct in_files*,struct param*,struct filter*,struct phenodt*, struct ests*,struct geno_block*,mstream&);
void pheno_read(struct param*,struct in_files*,struct filter*,struct phenodt*,Eigen::Ref<ArrayXb>,mstream&);
void tpheno_read(struct param*,struct in_files*,struct filter*,struct phenodt*,Eigen::Ref<ArrayXb>,mstream&);
void rm_phenoCols(Eigen::Ref<ArrayXb>,struct in_files*,struct param*,struct phenodt*,mstream&);
void covariate_read(struct param*,struct in_files*,struct filter*,struct phenodt*,Eigen::Ref<ArrayXb>,mstream&);
void setMasks(struct param*,struct filter*,struct phenodt*,mstream&);
void print_cc_info(struct param*,struct in_files*,struct phenodt*,mstream&);
void extract_interaction_snp(struct param*,struct in_files*,struct filter*,struct phenodt*,struct geno_block*,Eigen::Ref<ArrayXb>,mstream&);
int check_categories(std::vector<std::string>&,std::vector<std::map<std::string,int>>&,struct param*,struct filter*,mstream&);
Eigen::MatrixXd get_dummies(const Eigen::Ref<const Eigen::ArrayXd>&);
bool add_square_term(const Eigen::Ref<const Eigen::MatrixXd>&);
void extract_names(std::vector<std::string>&,std::map<std::string,int>&);
int getBasis(Eigen::MatrixXd&,struct param const*);
void QRcheck(Eigen::MatrixXd&,struct param*);
void check_sd(const Eigen::Ref<const Eigen::MatrixXd>&,int const&,double const&);
void residualize_phenotypes(struct param const*,struct phenodt*,const std::vector<std::string>&,mstream&);
bool residualize_matrix(Eigen::MatrixXd&,Eigen::ArrayXd&,struct param const*,struct phenodt*,mstream&);
void apply_QR(Eigen::MatrixXd&,struct param const*,bool const&);
void rescale_mat(Eigen::Ref<Eigen::MatrixXd>,struct param const*);
void prep_run(struct in_files*,struct filter*,struct param*,struct phenodt*,struct ests*,mstream&);
void check_blup(Eigen::Ref<ArrayXb>,struct in_files*,struct param*,mstream&);
void check_blup(std::map<std::string,bool>&,struct in_files*,struct param*,mstream&);
bool has_blup(std::string const&,std::map<std::string,bool> const&,struct param const*,mstream&);
void blup_read(struct in_files*,struct param*,struct phenodt*,struct ests*,mstream&);
void write_ids(struct in_files const*,struct param*,struct phenodt const*,mstream&);
void check_str(std::string&);
void apply_rint(struct phenodt*,struct param const*);
void rint_pheno(Eigen::Ref<Eigen::MatrixXd>,Eigen::Ref<ArrayXb>);
bool cmp_rank_pair(struct rank_pair&,struct rank_pair&);

#endif

