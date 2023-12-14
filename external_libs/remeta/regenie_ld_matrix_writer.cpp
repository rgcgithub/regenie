#ifdef WITH_HTSLIB
#include "regenie_ld_matrix_writer.hpp"

#include <Eigen/Dense>
using Eigen::VectorXd;
using Eigen::DiagonalMatrix;

void cov_to_corr(VectorXd& variances, MatrixXd& corr, const MatrixXd& cov) {
  variances = cov.diagonal();
  MatrixXd tmp = (variances.array() > 0).select(
                    (variances.array() > 0).select(variances, 1)
                                           .array()
                                           .sqrt()
                                           .inverse()
                                           .matrix(),
                    0).asDiagonal();
  corr = tmp * cov * tmp;
}

RegenieLDMatrixWriter::RegenieLDMatrixWriter()
 : BgzWriter()
 , idx() {}

RegenieLDMatrixWriter::RegenieLDMatrixWriter(string file_prefix, int sample_size)
 : BgzWriter(file_prefix + ".rg.ld", "w")
 , idx(file_prefix + ".rg.ld.idx.gz", "w") {
  this->write((int32_t)sample_size);
  if (sizeof(float) != 4) {
    throw runtime_error("bad float size: sizeof(float) != 4");
  }
}

void RegenieLDMatrixWriter::write_matrix_dense(const MatrixXd& ld_mat,
                                        const string& gene_name,
                                        const vector<string>& variant_ids) {
  if (this->is_closed()) {
    throw runtime_error("operating on a closed file");
  }

  int64_t addr = this->tell();
  size_t nrows = ld_mat.rows();
  size_t ncols = ld_mat.cols();
  this->write('d'); // dense
  this->write((int32_t)nrows);
  this->write((int32_t)0); // symmetry with sparse format

  if (nrows != ncols || nrows != variant_ids.size()) {
    throw runtime_error("dimension mismatch when writing LD matrix");
  } else if ( ((ld_mat - ld_mat.transpose()).array().abs() > 1e-3).any() ) {
    throw runtime_error("LD matrix must be symmetric");
  }

  for (size_t i = 0; i < nrows; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      this->write((float)ld_mat(i, j));
    }
  }

  this->write_idx_entry(gene_name, variant_ids, addr);
}

void RegenieLDMatrixWriter::write_matrix_sparse(const MatrixXd& ld_mat,
                                         const string& gene_name,
                                         const vector<string>& variant_ids,
                                         const double& sparsity_threshold) {
  size_t nrows = ld_mat.rows();
  size_t ncols = ld_mat.cols();
  if (nrows != ncols || nrows != variant_ids.size()) {
    throw runtime_error("dimension mismatch when writing LD matrix");
  } else if ( ((ld_mat - ld_mat.transpose()).array().abs() > 1e-3).any() ) {
    throw runtime_error("LD matrix should be symmetric.");
  } else if ( (ld_mat.diagonal().array() < 0).any()) {
    throw runtime_error("Diagonal elements of LD matrix should be non-negative.");
  }

  VectorXd variances(nrows);
  MatrixXd corr(nrows, ncols);
  cov_to_corr(variances, corr, ld_mat);
  this->write_sparse_header(gene_name, variances, variant_ids, sparsity_threshold);

  for (size_t i = 0; i < nrows; ++i) {
    for (size_t j = 0; j < i; ++j) {
      if (abs(corr(i, j)) > sparsity_threshold) {
        this->write_sparse_entry(
          sparse_matrix_entry {
            (int32_t)i,
            (int32_t)j,
            (float)corr(i, j)
          }
        );
      }
    }
  }
  this->write_sparse_footer();
}

void RegenieLDMatrixWriter::write_sparse_header(const string& gene_name,
                                         const VectorXd& variances,
                                         const vector<string>& variant_ids,
                                         const double& sparsity_threshold) {
  if (this->is_closed()) {
    throw runtime_error("operating on a closed file");
  }
  if (variant_ids.size() == 0) {
    throw runtime_error("writing an empty matrix");
  }
  
  int64_t addr = this->tell();
  size_t nrows = variances.size();
  this->write('s'); // sparse
  this->write((int32_t)nrows);
  this->write((float)sparsity_threshold);
  for (size_t i = 0; i < nrows; ++i) {
    this->write((float)variances[i]);
  }
  this->write_idx_entry(gene_name, variant_ids, addr);
}

void RegenieLDMatrixWriter::write_sparse_entry(const sparse_matrix_entry& entry) {
  this->write(entry);
}

void RegenieLDMatrixWriter::write_sparse_footer() {
  this->write(sparse_matrix_entry {
    (int32_t)-1,
    (int32_t)-1,
    (float)0
  });
}

void RegenieLDMatrixWriter::write_idx_entry(const string& gene_name,
                                     const vector<string>& variant_ids,
                                     const int64_t& addr) {
  string ids = "";
  for (size_t i = 0; i < variant_ids.size() - 1; ++i) {
    ids += variant_ids[i] + ",";
  }
  ids += variant_ids[variant_ids.size() - 1];
  idx.write(gene_name + "\t" + to_string(addr) + "\t" + ids + "\n");                            
}

void RegenieLDMatrixWriter::open(string file_prefix, int sample_size) {
  BgzWriter::open(file_prefix + ".metamat", "w");
  this->idx.open(file_prefix + ".metamat.idx.gz", "w");
  this->write((int32_t)sample_size);
  if (sizeof(float) != 4) {
    throw runtime_error("bad float size: sizeof(float) != 4");
  }
}

void RegenieLDMatrixWriter::close() {
  this->idx.close();
  BgzWriter::close();
}
#endif