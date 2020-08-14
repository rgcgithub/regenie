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

#ifndef PHENO_H
#define PHENO_H


struct phenodt {

  Eigen::MatrixXd new_cov;
  Eigen::MatrixXd phenotypes;
  Eigen::MatrixXd phenotypes_raw;
  MatrixXb masked_indivs;
  Eigen::ArrayXd Neff; // number of non-missing samples (per trait)
  Eigen::RowVectorXd scale_Y;

};


void read_pheno_and_cov(struct in_files*,struct param*,struct filter*,struct phenodt*, struct ests*,mstream&);
void pheno_read(struct param*,struct in_files*,struct filter*,struct phenodt*,ArrayXb&,mstream&);
void covariate_read(struct param*,struct in_files*,struct filter*,struct phenodt*,ArrayXb&,mstream&);
void getCovBasis(Eigen::MatrixXd&,struct param*);
void residualize_phenotypes(struct param*,struct phenodt*,const std::vector<std::string>&,mstream&);
void prep_run(struct in_files*,struct param*,struct phenodt*,struct ests*,mstream&);
void blup_read(struct in_files*,struct param*,struct phenodt*,struct ests*,mstream&);
double convertDouble(const std::string&,struct param*,mstream&);

#endif

