/*
 *
 * File derived from pgenlibr R library:
 * https://github.com/chrchang/plink-ng/tree/master/2.0/pgenlibr
 *
 * License info obtained from DESCRIPTION file:
 * https://github.com/chrchang/plink-ng/blob/master/2.0/pgenlibr/DESCRIPTION
 * -----------------------------------------------------
    Package: pgenlibr
    Type: Package
    Title: PLINK 2 Binary (.pgen) Reader
    Version: 0.2
    Date: 2019-07-10
    Author: Christopher Chang
    Maintainer: Christopher Chang <chrchang@alumni.caltech.edu>
    Description: A thin wrapper over PLINK 2's core libraries which provides an R
    interface for reading .pgen files.  A minimal .pvar loader is also included.
    License: LGPL (>= 3)
    Imports: Rcpp (>= 1.0.1)
    LinkingTo: Rcpp
 * -----------------------------------------------------

 *  Modified by Joelle Mbatchou - June 29 2020
 *  - removed functions that were for R
 *  - split file to header (added link to several standard C++ libraries)
 *  - modified remaining functions to be fully C/C++ compatible 
 *
 * This file remains under LGPL v3 license (license is in same directory as this file)
 */


#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <memory>
#include "pvar_ffi_support.h"
#include "pgenlib_ffi_support.h"
#include "include/pgenlib_read.h"


class PgenReader {
public:
  PgenReader();

  void Load(std::string filename, uint32_t cur_sample_ct, std::vector<int> sample_subset_1based, int nthr);

  uint32_t GetRawSampleCt() const;

  uint32_t GetSubsetSize() const;

  uint32_t GetVariantCt() const;

  uint32_t GetAlleleCt(uint32_t variant_idx) const;

  uint32_t GetMaxAlleleCt() const;

  bool HardcallPhasePresent() const;
  
  bool DosagePresent() const;

  void ReadIntHardcalls(std::vector<int>& buf, int variant_idx, int allele_idx);

  void ReadHardcalls(double* buf, size_t const& n, int const& thr, int variant_idx, int allele_idx);

  void Read(double* buf, size_t const& n, int const& thr, int variant_idx, int allele_idx);

  void Close();

  ~PgenReader();

private:
  plink2::PgenFileInfo* _info_ptr;
  uintptr_t* _allele_idx_offsetsp = nullptr;
  //plink2::RefcountedWptr* _allele_idx_offsetsp;
  plink2::RefcountedWptr* _nonref_flagsp;

  // have all below be threads specific
  std::vector<plink2::PgenReader*> _state_ptr;
  std::vector<plink2::PgrSampleSubsetIndex> _subset_index;
  std::vector<std::shared_ptr<plink2::PgenVariant>> _pgv;
  std::vector<uintptr_t*> _subset_include_interleaved_vec;
  std::vector<uint32_t*> _subset_cumulative_popcounts;
  std::vector<uint32_t> _subset_size;
  std::vector<uintptr_t*> _subset_include_vec;

  /*
  // kPglNypTransposeBatch (= 256) variants at a time, and then transpose
  uintptr_t* _multivar_vmaj_geno_buf;
  uintptr_t* _multivar_vmaj_phasepresent_buf;
  uintptr_t* _multivar_vmaj_phaseinfo_buf;
  uintptr_t* _multivar_smaj_geno_batch_buf;
  uintptr_t* _multivar_smaj_phaseinfo_batch_buf;
  uintptr_t* _multivar_smaj_phasepresent_batch_buf;
*/

  void SetSampleSubsetInternal(std::vector<int>& sample_subset_1based, int const& thr);
  void ReadAllelesPhasedInternal(int variant_idx);
};

