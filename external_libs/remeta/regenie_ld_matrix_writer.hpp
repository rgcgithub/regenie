#ifdef WITH_HTSLIB
#ifndef REGENIE_LD_MATRIX_WRITER_H
#define REGENIE_LD_MATRIX_WRITER_H

#include <string>
#include <vector>
using namespace std;

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include "bgz_writer.hpp"

struct sparse_matrix_entry {
  int32_t i;  // row index
  int32_t j;  // col index
  float data; // value at entry i,j
};

class RegenieLDMatrixWriter : public BgzWriter {
 public:
  RegenieLDMatrixWriter();

  RegenieLDMatrixWriter(string file_prefix, int sample_size);

  void write_matrix_dense(const MatrixXd& ld_mat,
                          const string& gene_name,
                          const vector<string>& variant_ids);

  void write_matrix_sparse(const MatrixXd& ld_mat,
                           const string& gene_name,
                           const vector<string>& variant_ids,
                           const double& sparsity_threshold);

  /* 
    These functions write sparse matrices in pieces. Writing each matrix
    requires 3 steps:
      1. write_sparse_header  : Writes the matrix header to the file and
                                writes its address to the index.
      2. write_sparse_entry   : Writes matrix entry to a file. This can be
                                called multiple times.
      3. write_sparse_footer  : Writes a footer signifying the end of the matrix.
  */
  void write_sparse_header(const string& gene_name,
                           const VectorXd& variances,
                           const vector<string>& variant_ids,
                           const double& sparsity_threshold);

  void write_sparse_entry(const sparse_matrix_entry& entry);

  void write_sparse_footer();

  void open(string file_prefix, int sample_size);

  void close();
 
 private:
  BgzWriter idx;

  void write_idx_entry(const string& gene_name,
                       const vector<string>& variant_ids,
                       const int64_t& addr);
};

#endif
#endif