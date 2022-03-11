/*

   This file is part of the regenie software package.

   Copyright (c) 2020-2022 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

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


#ifndef INTERACTION_H
#define INTERACTION_H

// for interaction testing
void get_interaction_terms(const int&,const int&,struct phenodt*,struct geno_block*,variant_block*,HLM&,struct param const*,mstream&);
void apply_interaction_tests(const int&,const int&,const int&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::RowVectorXd>&,std::string const&,std::string const&,struct phenodt*,HLM&,struct filter const*,struct in_files*,struct geno_block*,variant_block*,std::vector<snp> const&,struct ests*,struct f_ests*,struct param const*,mstream&);
void apply_interaction_tests_qt(const int&,const int&,const int&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::RowVectorXd>&,std::string const&,std::string const&,struct phenodt*,struct filter const*,struct in_files*,struct geno_block*,variant_block*,std::vector<snp> const&,struct param const*,mstream&);
void apply_interaction_tests_HLM(const int&,const int&,const int&,const Eigen::Ref<const Eigen::MatrixXd>&,const Eigen::Ref<const Eigen::RowVectorXd>&,std::string const&,std::string const&,struct phenodt*,HLM&,struct filter const*,struct in_files*,struct geno_block*,variant_block*,std::vector<snp> const&,struct param const*,mstream&);
void apply_interaction_tests_bt(const int&,const int&,const int&,std::string const&,std::string const&,struct phenodt*,struct filter const*,struct in_files*,struct geno_block*,variant_block*,std::vector<snp> const&,struct ests*,struct f_ests*,struct param const*,mstream&);
std::string apply_interaction_tests_firth(const int&,const int&,const int&,const int&,std::string const&,std::string const&,struct phenodt*,struct filter const*,struct in_files const*,struct geno_block*,variant_block*,std::vector<snp> const&,struct ests const*,struct f_ests const*,struct param const*,mstream&);


#endif
