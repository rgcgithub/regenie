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

#include "Regenie.hpp"
#include "Files.hpp"
#include "Geno.hpp"
#include "db/sqlite3.hpp"

using namespace std;
using namespace Eigen;
using namespace boost;



void prep_bgen(struct in_files* files, struct param* params, struct filter* filters, vector<snp>& snpinfo, map<int, vector<int>>& chr_map, BgenParser& bgen, mstream& sout){

  uint32_t nOutofOrder = 0;
  uint64 lineread = 0;
  std::string chromosome, rsid;
  uint32_t position ;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;
  std::vector< string > tmp_ids ;
  snp tmp_snp;
  BgenParser bgen_tmp;

  sout << left << std::setw(20) << " * bgen" << ": [" << files->bgen_file << "]" << endl;
  // open file and print file info
  bgen_tmp.open( files->bgen_file ) ;
  bgen_tmp.summarise( sout.coss ) ;
  bgen_tmp.summarise( cerr ) ;

  // get info for variants
  if( params->with_bgi ) read_bgi_file(bgen_tmp, files, params, filters, snpinfo, sout);
  else {
    tmp_snp.offset = bgen_tmp.get_position();
    while(bgen_tmp.read_variant( &chromosome, &position, &rsid, &alleles )) {

      bgen_tmp.ignore_probs();
      assert(alleles.size() == 2) ; // only bi-allelic allowed

      tmp_snp.chrom = chrStrToInt(chromosome, params->nChrom);
      if (tmp_snp.chrom == -1) {
        sout << "ERROR: Unknown chromosome code in bgen file."<< endl;
        exit(EXIT_FAILURE);
      }

      if( files->chr_read.empty() || (tmp_snp.chrom != files->chr_read.back()) ) files->chr_read.push_back(tmp_snp.chrom);

      tmp_snp.physpos = position;
      tmp_snp.ID = rsid;
      if( params->ref_first ) { // reference is first (i.e. allele0)
        tmp_snp.allele1 = alleles[0];
        tmp_snp.allele2 = alleles[1];
      } else {
        tmp_snp.allele1 = alleles[1];
        tmp_snp.allele2 = alleles[0]; // switch so allele0 is ALT
      }

      // check if snps are in order (same chromosome & non-decreasing positions)
      if (!snpinfo.empty() && (tmp_snp.chrom == snpinfo.back().chrom) && ( (tmp_snp.physpos < snpinfo.back().physpos) )) nOutofOrder++;

      lineread++;

      // if specified chrlist/range
      if(
          (params->select_chrs && !in_chrList(tmp_snp.chrom, filters))
          ||
          (params->set_range && !in_range(chromosome, position, params))
        ) {
        // go to next variant (get its offset first)
        tmp_snp.offset = bgen_tmp.get_position();
        continue;
      }

      // make list of variant IDs if inclusion/exclusion file is given
      if(params->rm_snps || params->keep_snps) 
        filters->snpID_to_ind.insert( std::make_pair( tmp_snp.ID, snpinfo.size() ) );

      // keep track of how many included snps per chromosome there are
      files->chr_counts[tmp_snp.chrom-1]++;
      snpinfo.push_back(tmp_snp);

      tmp_snp.offset = bgen_tmp.get_position();
    }

    if (!params->test_mode && (nOutofOrder > 0)) sout << "WARNING: Total number of snps out-of-order in bgen file : " << nOutofOrder << endl;
  }

  // check if should mask snps
  check_snps_include_exclude(files, params, filters, snpinfo, chr_map, sout);


  // get info on samples
  params->n_samples  = bgen_tmp.number_of_samples();

  // get sample IDs (from sample file or directly from bgen file)
  if( params->bgenSample ) {
    read_bgen_sample(files->sample_file, params, tmp_ids, sout);
  } else {
    bgen_tmp.get_sample_ids(
        [&tmp_ids]( std::string const& id ) { tmp_ids.push_back( id ) ; } );
    // assume all females
    params->sex.resize(params->n_samples, false);
  }

  // check duplicates -- if not, store in map
  for(size_t i = 0; i < params->n_samples; i++) {
    if (params->FID_IID_to_ind.find(tmp_ids[i]) != params->FID_IID_to_ind.end()) {
      sout << "ERROR: Duplicate individual in bgen file : FID_IID=" << tmp_ids[i] << endl;
      exit(EXIT_FAILURE);
    }
    params->FID_IID_to_ind.insert( std::make_pair( tmp_ids[i], i ) );
  }

  // check if should mask samples
  check_samples_include_exclude(files, params, filters, sout);

  // prepare file
  if( !params->streamBGEN ) {
    // setup file for reading the genotype probabilities later
    bgen.open( files->bgen_file ) ;
  }

  if (params->test_mode) params->dosage_mode = true;
}


// read .bgi file to get SNP info
void read_bgi_file(BgenParser& bgen, struct in_files* files, struct param* params, struct filter* filters, std::vector<snp>& snpinfo, mstream& sout){

  int nalleles;
  uint64 lineread = 0, variant_bgi_size, variant_bgen_size;
  string bgi_file = files->bgen_file + ".bgi";
  string sql_query = "SELECT * FROM Variant";
  snp tmp_snp;
  sqlite3* db;
  sqlite3_stmt* stmt;

  uint32_t n_variants = bgen.number_of_variants();
  uint32_t position ;
  std::string chromosome, rsid;
  std::vector< std::string > alleles ;

  // edit sql statement if chromosome position range is given
  if( params->set_range ){
    string tmp_q = sql_query + " WHERE chromosome=" + params->range_chr + " AND position>=" + to_string(params->range_min) + " AND position<=" + to_string(params->range_max);
    sql_query = tmp_q;
  } else if( params->select_chrs ){
    
    string tmp_q = sql_query + " WHERE chromosome IN (" + bgi_chrList(filters) + ")";
    sql_query = tmp_q;
  }

  sout << "   -index bgi file [" << bgi_file<< "]" << endl;
  if( sqlite3_open( bgi_file.c_str(), &db ) != SQLITE_OK ) {
    sout <<  "ERROR: " << sqlite3_errmsg(db) << endl;
    exit(EXIT_FAILURE);
  }


  // header: chromosome|position|rsid|number_of_alleles|allele1|allele2|file_start_position|size_in_bytes
  if( sqlite3_prepare_v2( db, sql_query.c_str(), -1, &stmt, NULL ) != SQLITE_OK ){
    sout << "ERROR: " << sqlite3_errmsg(db) << endl;
    exit(EXIT_FAILURE);
  }

  bool done = false;
  uint32_t nOutofOrder = 0;
  while (!done) {
    switch (sqlite3_step(stmt)) {
      case SQLITE_ROW:

        chromosome = std::string( (char *) sqlite3_column_text(stmt, 0) );
        tmp_snp.chrom = chrStrToInt(chromosome, params->nChrom);
        if (tmp_snp.chrom == -1) {
          sout << "ERROR: Unknown chromosome code in bgi file (=" << chromosome << ").\n";
          exit(EXIT_FAILURE);
        }
        if( files->chr_read.empty() || (tmp_snp.chrom != files->chr_read.back()) ) files->chr_read.push_back(tmp_snp.chrom);

        tmp_snp.physpos = strtoul( (char *) sqlite3_column_text(stmt, 1), NULL, 10);
        tmp_snp.ID = std::string( (char *) sqlite3_column_text(stmt, 2) );
        nalleles = atoi( (char *) sqlite3_column_text(stmt, 3) );
        assert(nalleles == 2) ; // only bi-allelic allowed
        if( params->ref_first ){ // reference is first
          tmp_snp.allele1 = std::string( (char *) sqlite3_column_text(stmt, 4) );
          tmp_snp.allele2 = std::string( (char *) sqlite3_column_text(stmt, 5) );
        } else {
          tmp_snp.allele1 = std::string( (char *) sqlite3_column_text(stmt, 5) );
          tmp_snp.allele2 = std::string( (char *) sqlite3_column_text(stmt, 4) ); // switch so allele0 is ALT
        }
        tmp_snp.offset = strtoull( (char *) sqlite3_column_text(stmt, 6), NULL, 10);

        // check if matches with info from bgenparser for first read variant
        if( snpinfo.empty() ){
          bgen.jumpto(tmp_snp.offset);
          bgen.read_variant( &chromosome, &position, &rsid, &alleles );
          bgen.ignore_probs();
          variant_bgen_size = bgen.get_position() - tmp_snp.offset;
          variant_bgi_size = strtoull( (char *) sqlite3_column_text(stmt, 7), NULL, 10);
          assert( tmp_snp.ID == rsid );
          assert( variant_bgi_size == variant_bgen_size );
        }

        // check if snps are in order (same chromosome & non-decreasing positions)
        if (!snpinfo.empty()
            && (tmp_snp.chrom == snpinfo.back().chrom)
            && ( (tmp_snp.physpos < snpinfo.back().physpos) ))
          nOutofOrder++;

        lineread++;

        // make list of variant IDs if inclusion/exclusion file is given
        if(params->rm_snps || params->keep_snps)
          filters->snpID_to_ind.insert( std::make_pair( tmp_snp.ID, snpinfo.size() ) );

        // keep track of how many included snps per chromosome there are
        files->chr_counts[tmp_snp.chrom-1]++;
        snpinfo.push_back(tmp_snp);
        break;

      case SQLITE_DONE:
        done = true;
        break;

      default:
        sout << "ERROR: Failed reading file (" << sqlite3_errmsg(db) << ").\n";
        exit(EXIT_FAILURE);
    }
  }

  sqlite3_finalize(stmt);
  sqlite3_close(db);

  if( !params->set_range && !params->select_chrs) assert( lineread == n_variants );
  if (!params->test_mode && (nOutofOrder > 0)) sout << "WARNING: Total number of snps out-of-order in bgen file : " << nOutofOrder << endl;

}


void read_bgen_sample(const string sample_file, struct param* params, std::vector<string> &ids, mstream& sout){

  int nline = 0;
  string FID, IID, line, tmp_str;
  std::vector<string> IDvec;
  ifstream myfile;
  if( params->write_samples ) IDvec.resize(2);

  sout << "   -sample file: " << sample_file << endl;
  myfile.open (sample_file, ios::in);
  if (!myfile.is_open()) {   
    sout << "ERROR: Cannot open sample file." << endl;
    exit(EXIT_FAILURE);
  }

  // read fid/iid information
  while (getline (myfile,line)) {
    std::istringstream iss(line);

    if( !(iss >> FID >> IID) ){
      sout << "ERROR: Incorrectly formatted sample file at line" << ids.size() + 1 << endl;
      exit(EXIT_FAILURE);
    }

    // check first two lines for correct format
    if(nline == 0){
      if( (FID != "ID_1") || (IID != "ID_2") ) {
        sout << "ERROR: Header of the sample file must start with: ID_1 ID_2" << endl;
        exit(EXIT_FAILURE);
      }
    } else if(nline == 1){
      if( (FID != "0") || (IID != "0") ) {
        sout << "ERROR: Second line of sample file must start with: 0 0" << endl;
        exit(EXIT_FAILURE);
      }
    } else {
      tmp_str = FID + "_" + IID;
      ids.push_back(tmp_str);
      if(params->write_samples) {
        IDvec[0] = FID;
        IDvec[1] = IID;
        params->FIDvec.push_back(IDvec);
      }
    }

    // get sex into IID (if no sex column, set to 0)
    if( !(iss >> FID >> IID) ){
      params->sex.push_back(false);
    } else {
      // save males as 1 (else assume all females=0)
      if( IID == "1" ) params->sex.push_back(true);
      else params->sex.push_back(false);
    }

    nline++;
  }

  if( params->n_samples != ids.size() ){
    sout << "ERROR: Number of samples in BGEN file does not match that in the sample file." << endl;
    exit(EXIT_FAILURE);
  }

  myfile.close();
}


void read_bed_bim_fam(struct in_files* files, struct param* params, struct filter* filters, vector<snp>& snpinfo, map<int,vector<int>>& chr_map, mstream& sout) {

  uint32_t nsamples_bed;
  read_bim(files, params, filters, snpinfo, sout);
  // check if should mask snps
  check_snps_include_exclude(files, params, filters, snpinfo, chr_map, sout);

  read_fam(files, params, sout);
  nsamples_bed = params->n_samples;
  // check if should mask samples
  check_samples_include_exclude(files, params, filters, sout);

  prep_bed(nsamples_bed, files, sout);
}


void read_bim(struct in_files* files, struct param* params, struct filter* filters, vector<snp>& snpinfo, mstream& sout) {

  uint32_t nOutofOrder = 0;
  int minChr_read = 0; // enforce that chromosomes in file are sorted
  uint64 lineread = 0;
  std::vector< string > tmp_str_vec ;
  snp tmp_snp;
  string line, fname;
  ifstream myfile;

  fname = files->bed_prefix + ".bim";
  sout << left << std::setw(20) << " * bim" << ": [" << fname << "] " << flush;
  myfile.open(fname.c_str());
  if (!myfile.is_open()) {   
    sout << "ERROR: Cannot open bim file." << endl;
    exit(EXIT_FAILURE);
  }
  //if(params->set_range) cerr << params->range_chr << "\t" << params->range_min << "\t" << params->range_max<< endl;

  while (getline(myfile, line)) {
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 6 ){
      sout << "ERROR: Incorrectly formatted bim file at line " << snpinfo.size()+1 << endl;
      exit(EXIT_FAILURE);
    }

    tmp_snp.chrom = chrStrToInt(tmp_str_vec[0], params->nChrom);
    tmp_snp.ID = tmp_str_vec[1];
    //tmp_snp.genpos = std::stod( tmp_str_vec[2]);
    tmp_snp.physpos = std::stoul( tmp_str_vec[3],nullptr,0);
    if( params->ref_first ){ // reference is first
      tmp_snp.allele1 = tmp_str_vec[4];
      tmp_snp.allele2 = tmp_str_vec[5];
    } else { // reference is last
      tmp_snp.allele1 = tmp_str_vec[5];
      tmp_snp.allele2 = tmp_str_vec[4];
    }
    tmp_snp.offset = lineread;

    if (tmp_snp.chrom == -1) {
      sout << "ERROR: Unknown chromosome code in bim file at line " << snpinfo.size()+1 << endl;
      exit(EXIT_FAILURE);
    }

    if( files->chr_read.empty() || (tmp_snp.chrom != files->chr_read.back() ) ) {
      files->chr_read.push_back(tmp_snp.chrom);
      if( tmp_snp.chrom <= minChr_read ){
        sout << "ERROR: Chromosomes in bim file are not in ascending order.\n";
        exit(EXIT_FAILURE);
      } else minChr_read = tmp_snp.chrom;
    }

    // check if snps are in order (same chromosome & non-decreasing positions)
    if (!snpinfo.empty() && (tmp_snp.chrom == snpinfo.back().chrom) && ( (tmp_snp.physpos < snpinfo.back().physpos) )) nOutofOrder++;

    lineread++;

    // if specified chrlist/range
    if(
        (params->select_chrs && !in_chrList(tmp_snp.chrom, filters))
        ||
        (params->set_range && !in_range(tmp_str_vec[0], tmp_snp.physpos, params))
      ) continue;

    // make list of variant IDs if inclusion/exclusion file is given
    if(params->rm_snps || params->keep_snps)
      filters->snpID_to_ind.insert( std::make_pair( tmp_snp.ID, snpinfo.size() ) );

    // keep track of how many included snps per chromosome there are
    files->chr_counts[tmp_snp.chrom-1]++;
    snpinfo.push_back(tmp_snp);
  }

  sout << "n_snps = " << lineread << endl;

  if (!params->test_mode && (nOutofOrder > 0)) sout << "WARNING: Total number of snps out-of-order in bim file : " << nOutofOrder << endl;

  myfile.close();

}


void read_fam(struct in_files* files, struct param* params, mstream& sout) {

  int lineread = 0;
  string line, tmp_id, fname;
  std::vector< string > tmp_str_vec, IDvec;
  ifstream myfile;
  if( params->write_samples ) IDvec.resize(2);

  fname = files->bed_prefix + ".fam";
  sout << left << std::setw(20) << " * fam" << ": [" << fname << "] ";
  myfile.open(fname.c_str());
  if (!myfile.is_open()) {   
    sout << "ERROR: Cannot open fam file." << endl;
    exit(EXIT_FAILURE);
  }

  while (getline(myfile, line)) {
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 6 ){
      sout << "ERROR: Incorrectly formatted fam file at line " << lineread + 1 << endl;
      exit(EXIT_FAILURE);
    }

    tmp_id = tmp_str_vec[0] + "_" + tmp_str_vec[1];

    // check duplicates -- if not, store in map
    if (params->FID_IID_to_ind.find(tmp_id) != params->FID_IID_to_ind.end()) {
      sout << "ERROR: Duplicate individual in fam file : FID_IID=" << tmp_id << endl;
      exit(EXIT_FAILURE);
    }
    params->FID_IID_to_ind.insert( std::make_pair( tmp_id, lineread ) );
    if(params->write_samples) {
      IDvec[0] = tmp_str_vec[0];
      IDvec[1] = tmp_str_vec[1];
      params->FIDvec.push_back(IDvec);
    }

    // save males as 1 (else assume all females=0)
    if( tmp_str_vec[4] == "1" ) params->sex.push_back(true);
    else params->sex.push_back(false);

    lineread++;
  }

  myfile.close();
  params->n_samples = lineread;

  sout << "n_samples = " << params->n_samples << endl;
}


void prep_bed(const uint32_t& nsamples, struct in_files* files, mstream& sout) {

  string fname;

  fname = files->bed_prefix + ".bed";
  sout << left << std::setw(20) << " * bed" << ": [" << fname << "]" << endl;
  files->bed_ifstream.open(fname.c_str(), std::ios::in | std::ios::binary);
  if (!files->bed_ifstream.is_open()) {   
    sout << "ERROR: Cannot open bed file." << endl;
    exit(EXIT_FAILURE);
  }

  uchar header[3];
  files->bed_ifstream.read( reinterpret_cast<char *> (&header[0]), 3);
  if ( (header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01) ) {
    sout << "ERROR: Incorrect magic number in bed file.\n";
    exit(EXIT_FAILURE);
  }

  // size of genotype block [(n+3)/4 = ceil(n/4.0)]
  files->bed_block_size = (nsamples+3)>>2;
  files->inbed.resize( files->bed_block_size );
}



void read_pgen_pvar_psam(struct in_files* files, struct param* params, struct filter* filters, struct geno_block* gblock, vector<snp>& snpinfo, map<int,vector<int>>& chr_map, mstream& sout) {

  uint32_t pgen_nvariants, pgen_nsamples;

  pgen_nvariants = read_pvar(files, params, filters, snpinfo, sout);;
  // check if should mask snps
  check_snps_include_exclude(files, params, filters, snpinfo, chr_map, sout);

  read_psam(files, params, sout);
  sout << "n_samples = " << params->n_samples << endl;
  pgen_nsamples = params->n_samples;
  // check if should mask samples
  check_samples_include_exclude(files, params, filters, sout);

  prep_pgen(pgen_nsamples, pgen_nvariants, files, filters, gblock, sout);

  params->dosage_mode = gblock->pgr.DosagePresent();

}


uint64 read_pvar(struct in_files* files, struct param* params, struct filter* filters, vector<snp>& snpinfo, mstream& sout) {

  uint32_t nOutofOrder = 0;
  int minChr_read = 0; // enforce that chromosomes in file are sorted
  uint64 lineread = 0;
  std::vector< string > tmp_str_vec ;
  snp tmp_snp;
  string line, fname;
  ifstream myfile;

  fname = files->pgen_prefix + ".pvar";
  sout << left << std::setw(20) << " * pvar" << ": [" << fname << "] " << flush;
  myfile.open(fname.c_str());
  if (!myfile.is_open()) {
    sout << "ERROR: Cannot open pvar file." << endl;
    exit(EXIT_FAILURE);
  }

  while (getline(myfile, line)) { // skip to main header line
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 1 ){
      sout << "ERROR: No blank lines should be before the header line in pvar file.\n";
      exit(EXIT_FAILURE);
    }

    if( tmp_str_vec[0] == "#CHROM" ) break;
  }

  // check header
  if( (tmp_str_vec.size() < 5) || (tmp_str_vec[1] != "POS") || (tmp_str_vec[2] != "ID") || (tmp_str_vec[3] != "REF") || (tmp_str_vec[4] != "ALT") ){
    cerr << "ERROR: Header of pvar file does not have correct format.\n";
    exit(EXIT_FAILURE);
  }

  while (getline(myfile, line)) {
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 5 ){
      sout << "ERROR: Incorrectly formatted pvar file at line " << snpinfo.size()+1 << endl;
      exit(EXIT_FAILURE);
    }

    tmp_snp.chrom = chrStrToInt(tmp_str_vec[0], params->nChrom);
    tmp_snp.physpos = std::stoul( tmp_str_vec[1],nullptr,0);
    tmp_snp.ID = tmp_str_vec[2];
    tmp_snp.allele1 = tmp_str_vec[3];
    tmp_snp.allele2 = tmp_str_vec[4];
    tmp_snp.offset = lineread; // store index in file

    if (tmp_snp.chrom == -1) {
      sout << "ERROR: Unknown chromosome code in pvar file at line " << snpinfo.size()+1 << endl;
      exit(EXIT_FAILURE);
    }

    if( files->chr_read.empty() || (tmp_snp.chrom != files->chr_read.back() ) ) {
      files->chr_read.push_back(tmp_snp.chrom);
      if( tmp_snp.chrom <= minChr_read ){
        sout << "ERROR: Chromosomes in pvar file are not in ascending order.\n";
        exit(EXIT_FAILURE);
      } else minChr_read = tmp_snp.chrom;
    }

    // check if snps are in order (same chromosome & non-decreasing positions)
    if (!snpinfo.empty() && (tmp_snp.chrom == snpinfo.back().chrom) && ( (tmp_snp.physpos < snpinfo.back().physpos) )) nOutofOrder++;

    lineread++;

    // if specified chrlist/range
    if(
        (params->select_chrs && !in_chrList(tmp_snp.chrom, filters))
        ||
        (params->set_range && !in_range(tmp_str_vec[0], tmp_snp.physpos, params))
      ) continue;

    // make list of variant IDs if inclusion/exclusion file is given
    if(params->rm_snps || params->keep_snps)
      filters->snpID_to_ind.insert( std::make_pair( tmp_snp.ID, snpinfo.size() ) );

    // keep track of how many included snps per chromosome there are
    files->chr_counts[tmp_snp.chrom-1]++;
    snpinfo.push_back(tmp_snp);
  }

  sout << "n_snps = " <<  lineread << endl;

  if (!params->test_mode && (nOutofOrder > 0)) sout << "WARNING: Total number of snps out-of-order in bim file : " << nOutofOrder << endl;

  myfile.close();

  return lineread;
}


void read_psam(struct in_files* files, struct param* params, mstream& sout) {

  int lineread = 0, sex_col = 0;
  bool col_found = false;
  string line, tmp_id, fname;
  std::vector< string > tmp_str_vec, IDvec;
  ifstream myfile;
  if( params->write_samples ) IDvec.resize(2);

  fname = files->pgen_prefix + ".psam";
  sout << left << std::setw(20) << " * psam" << ": [" << fname << "] " << flush;
  myfile.open(fname.c_str());
  if (!myfile.is_open()) {
    sout << "ERROR: Cannot open psam file." << endl;
    exit(EXIT_FAILURE);
  }

  while (getline(myfile, line)) { // skip to main header line
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 1 ){
      cerr << "ERROR: No blank lines should be before the header line in psam file.\n";
      exit(EXIT_FAILURE);
    }

    if( tmp_str_vec[0] == "#FID" ) break;
  }

  // check header
  if( (tmp_str_vec.size() < 2) || (tmp_str_vec[1] != "IID")){
    cerr << "ERROR: Header does not have the correct format.\n";
    exit(EXIT_FAILURE);
  }
  // find if sex column is present
  auto scol = find(tmp_str_vec.begin(), tmp_str_vec.end(), "SEX");
  col_found = scol != tmp_str_vec.end();
  if(col_found) sex_col = std::distance(tmp_str_vec.begin(), scol);

  while (getline(myfile, line)) {
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 3 ){
      sout << "ERROR: Incorrectly formatted psam file at line " << lineread + 1 << endl;
      exit(EXIT_FAILURE);
    }

    tmp_id = tmp_str_vec[0] + "_" + tmp_str_vec[1];

    // check duplicates -- if not, store in map
    if (params->FID_IID_to_ind.find(tmp_id ) != params->FID_IID_to_ind.end()) {
      sout << "ERROR: Duplicate individual in fam file : FID_IID=" << tmp_id << endl;
      exit(EXIT_FAILURE);
    }
    params->FID_IID_to_ind.insert( std::make_pair( tmp_id, lineread ) );
    if(params->write_samples) {
      IDvec[0] = tmp_str_vec[0];
      IDvec[1] = tmp_str_vec[1];
      params->FIDvec.push_back(IDvec);
    }

    // save males as 1 (else assume all females=0)
    if( col_found && tmp_str_vec[sex_col] == "1" ) params->sex.push_back(true);
    else params->sex.push_back(false);

    lineread++;
  }

  myfile.close();
  params->n_samples = lineread;
}


void prep_pgen(const uint32_t pgen_ns, const uint32_t pgen_nv, struct in_files* files, struct filter* filters, struct geno_block* gblock, mstream& sout){

  int pgen_samples, pgen_variants, pgen_ac;
  vector<int> subset_indices_1based;
  string fname;

  fname = files->pgen_prefix + ".pgen";
  sout << left << std::setw(20) << " * pgen" << ": [" << fname << "] " << flush;

  // set subset when samples have been excluded from analysis
  if( filters->ind_in_analysis.size() < pgen_ns ){
    // need to create vector of indices to keep (1-based)
    for( size_t i = 0; i < pgen_ns; i++)
      if(!filters->ind_ignore(i))
        subset_indices_1based.push_back(i+1);
  }

  gblock->pgr.Load(fname, pgen_ns, subset_indices_1based);
  pgen_samples = gblock->pgr.GetRawSampleCt();
  pgen_variants = gblock->pgr.GetVariantCt();
  pgen_ac = gblock->pgr.GetMaxAlleleCt();

  if(pgen_samples != (int) pgen_ns){
    cerr << "ERROR: Number of samples in pgen file and psam file don't match.\n";
    exit(EXIT_FAILURE);
  }
  if(pgen_variants != (int) pgen_nv){
    cerr << "ERROR: Number of variants in pgen file and pvar file don't match.\n";
    exit(EXIT_FAILURE);
  }
  if(pgen_ac != 2){
    cerr << "ERROR: Only bi-allelic variants are accepted.\n";
    exit(EXIT_FAILURE);
  }

  gblock->genobuf.resize(filters->ind_in_analysis.size());
  sout << endl;

}


// determine if snps should be included/excluded for step 1
void check_snps_include_exclude(struct in_files* files, struct param* params, struct filter* filters, vector<snp>& snpinfo, map<int,vector<int>>& chr_map, mstream& sout){

  vector<snp> tmp_snpinfo;
  params->n_variants = snpinfo.size(); // current variants count
  if(params->set_range)
    sout << "   -number of variants after filtering on range = " << params->n_variants << endl;

  // set all masks to false
  filters->geno_mask.assign(params->n_variants, false);

  // if inclusion/exclusion file is given
  if(params->rm_snps || params->keep_snps) {

    // apply masking to snps
    if( params->rm_snps ) set_snps_to_rm(files, params, filters, snpinfo, sout);
    else if( params->keep_snps ) set_snps_to_keep(files, params, filters, snpinfo, sout);

    // delete snpID map
    filters->snpID_to_ind.clear();

    // make snpinfo only with kept elements
    tmp_snpinfo.reserve( params->n_variants );
    for(size_t i = 0; i < filters->geno_mask.size(); i++){
      if(filters->geno_mask[i]) continue;
      //cerr << snpinfo[i].ID << endl;
      tmp_snpinfo.push_back( snpinfo[i] );
    }

    snpinfo.clear();
    std::vector<snp>().swap(snpinfo); // free memory
    snpinfo = tmp_snpinfo;

  }

  // check nonzero
  if(params->n_variants == 0){
    sout << "ERROR: No variant left to include in analysis.\n";
    exit(EXIT_FAILURE);
  }

  // go through each chromosome in order & save number of snps
  // and save how many are actually read
  vector<int> tmp_v;
  tmp_v.resize(2, 0);
  for(size_t j = 0; j < files->chr_read.size(); j++){
    int i = files->chr_read[j];
    tmp_v[0] = files->chr_counts[i-1];
    chr_map.insert(pair<int, vector<int> >(i, tmp_v));
  }


  if(params->rm_snps || params->keep_snps) {
    if(params->keep_snps) {
      sout << "   -keeping only variants specified in [" << files->file_snps_include << "]" << endl;
      if(filters->geno_mask.size() == params->n_variants) params->keep_snps = false;
    } else if(params->rm_snps) {
      sout << "   -removing variants specified in [" << files->file_snps_exclude << "]" << endl;
      if(filters->geno_mask.size() == params->n_variants) params->rm_snps = false;
    }
    sout << "     +number of variants remaining in the analysis = " << params->n_variants << endl;
  }
}


// snps to retain in step 1 analysis
void set_snps_to_keep(struct in_files* files, struct param* params, struct filter* filters, vector<snp>& snpinfo, mstream& sout) {

  uint64 snp_pos;
  string line;
  std::vector< string > tmp_str_vec ;
  Files myfile;

  myfile.openForRead (files->file_snps_include, sout);

  // set all the masks to true
  params->n_variants = 0;
  std::fill(filters->geno_mask.begin(), filters->geno_mask.end(), true);

  // set chr counts to 0
  std::fill(files->chr_counts.begin(), files->chr_counts.end(), 0);

  while( myfile.readLine(line) ){
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 1 ){
      sout << "ERROR: Incorrectly formatted file specified by --extract." << endl;
      exit(EXIT_FAILURE);
    }

    if (filters->snpID_to_ind.find(tmp_str_vec[0]) != filters->snpID_to_ind.end()) {
      snp_pos = filters->snpID_to_ind[ tmp_str_vec[0] ];
      filters->geno_mask[ snp_pos ] = false;
      // adjust counts
      files->chr_counts[ snpinfo[ snp_pos ].chrom - 1 ]++;
      params->n_variants++;
      //cerr << (filters->geno_mask[ snp_pos ] ? "N" : "Y") << " " << tmp_str_vec[0] << " "<< snp_pos << " " << params->n_variants << endl;
    }
  }

  myfile.closeFile();

}


// snps to exclude from step 1 analysis
void set_snps_to_rm(struct in_files* files, struct param* params, struct filter* filters, vector<snp>& snpinfo, mstream& sout) {

  uint64 snp_pos;
  string line;
  std::vector< string > tmp_str_vec ;
  Files myfile;

  myfile.openForRead (files->file_snps_exclude, sout);

  while( myfile.readLine(line) ){
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 1 ){
      sout << "ERROR: Incorrectly formatted file specified by --exclude." << endl;
      exit(EXIT_FAILURE);
    }

    if (filters->snpID_to_ind.find(tmp_str_vec[0]) != filters->snpID_to_ind.end()) {
      snp_pos = filters->snpID_to_ind[ tmp_str_vec[0] ];
      filters->geno_mask[ snp_pos ] = true;
      // adjust counts
      files->chr_counts[ snpinfo[ snp_pos ].chrom - 1 ]--;
      params->n_variants--;
    }
  }

  myfile.closeFile();
}

void check_samples_include_exclude(struct in_files* files, struct param* params, struct filter* filters, mstream& sout){

  uint32_t ind_pos = 0, cum_pos;
  string ind_ID;
  std::map <std::string, uint32_t> new_map;
  std::map <std::string, uint32_t>::iterator itr;
  vector< string > allIDs;
  vector< vector<string> > newFIDs;

  //  keep track of samples to remove
  filters->ind_in_analysis = ArrayXb::Constant(params->n_samples, true);

  if( params->rm_indivs )
    set_IDs_to_rm(files, filters, params, sout);
  else if( params->keep_indivs )
    set_IDs_to_keep(files, filters, params, sout);

  // to keep track of individual to exclude (i.e. not stored in memory)
  // this is used when reading in genotypes
  filters->ind_ignore = !filters->ind_in_analysis;

  if( params->rm_indivs || params->keep_indivs ) {

    // need to re-assign indices
    // retrieve all sample IDs (need to keep same order as in genotype file)
    allIDs.resize( params->n_samples );
    for (itr = params->FID_IID_to_ind.begin(); itr != params->FID_IID_to_ind.end(); ++itr) {
      ind_ID = itr->first;
      ind_pos = itr->second;
      allIDs[ ind_pos ] = ind_ID;
    }

    // create new map
    cum_pos = 0;
    for( size_t j = 0; j < params->n_samples; j++){
      if( !filters->ind_ignore(j) ){
        new_map.insert( std::make_pair( allIDs[j] , cum_pos ) );
        if( params->write_samples ) newFIDs.push_back( params->FIDvec[j] );
        cum_pos++;
      }
    }

    // save map
    params->FID_IID_to_ind = new_map;
    if( params->write_samples ) params->FIDvec = newFIDs;

    // resize ind_in_analysis
    filters->ind_in_analysis = ArrayXb::Constant(cum_pos, true);
  }


  params->n_samples = filters->ind_in_analysis.cast<int>().sum();
}


void set_IDs_to_keep(struct in_files* files, struct filter* filters, struct param* params, mstream& sout) {

  uint32_t n_kept = 0;
  string line;
  std::vector< string > tmp_str_vec ;
  findID person;
  Files myfile;

  filters->ind_in_analysis.fill(false);

  // track individuals to include -> remaining are ignored
  myfile.openForRead (files->file_ind_include, sout);

  while( myfile.readLine(line) ){
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 2 ){
      sout << "ERROR: Incorrectly formatted file specified by --keep." << endl;
      exit(EXIT_FAILURE);
    }

    person = getIndivIndex(tmp_str_vec[0], tmp_str_vec[1], params, sout);
    if(!person.is_found) continue;

    filters->ind_in_analysis(person.index) = true;
    n_kept++;
  }

  sout << "   -keeping only individuals specified in [" << files->file_ind_include<< "]" << endl;

  // check size
  if( n_kept < 1 ) {
    sout << "ERROR: None of the individuals are in the genotype file.\n";
    exit(EXIT_FAILURE);
  }

  sout << "     +number of genotyped individuals to keep in the analysis = " << n_kept << endl;

  myfile.closeFile();
}

void set_IDs_to_rm(struct in_files* files, struct filter* filters, struct param* params, mstream& sout) {

  uint32_t n_rm = 0;
  string line;
  std::vector< string > tmp_str_vec ;
  findID person;
  Files myfile;

  myfile.openForRead (files->file_ind_exclude, sout);

  while( myfile.readLine(line) ){
    boost::algorithm::split(tmp_str_vec, line, is_any_of("\t "));

    if( tmp_str_vec.size() < 2 ){
      sout << "ERROR: Incorrectly formatted file specified by --remove." << endl;
      exit(EXIT_FAILURE);
    }

    person = getIndivIndex(tmp_str_vec[0], tmp_str_vec[1], params, sout);
    if(!person.is_found) continue;

    filters->ind_in_analysis(person.index) = false;
    n_rm++;
  }

  sout << "   -removing individuals specified in [" << files->file_ind_exclude<< "]" << endl;

  if( n_rm == params->n_samples ){
    sout << "ERROR: No individuals remain in the analysis.\n";
    exit(EXIT_FAILURE);
  }

  sout << "     +number of genotyped individuals to exclude from the analysis = " << n_rm << endl;

  myfile.closeFile();
}


void get_G(const int block, const int bs, const int chrom, const uint32_t snpcount, vector<snp>& snpinfo, struct param* params, struct in_files* files, struct geno_block* gblock, struct filter* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, mstream& sout){

  auto t1 = std::chrono::high_resolution_clock::now();
  sout << " block [" << block + 1 << "] : " << flush;

  // prepare vector to store non_zero indices (for SPA)
  if(params->use_SPA) {
    gblock->non_zero_indices_G.resize(bs);
    for( std::size_t i = 0; i < gblock->non_zero_indices_G.size(); ++i )
      gblock->non_zero_indices_G[i].clear();
  }

  // set genotype counts to 0
  if(params->htp_out){
    for( int i = 0; i < bs; ++i ) gblock->genocounts[i].setZero();
  }

  if(params->file_type == "bed")
    readChunkFromBedFileToG(bs, chrom, snpcount, snpinfo, params, files, gblock, filters, masked_indivs, phenotypes_raw, sout);
  else if(params->file_type == "pgen")
    readChunkFromPGENFileToG(bs, snpcount, snpinfo, params, gblock, filters, masked_indivs, sout);
  else
    readChunkFromBGENFileToG(bs, chrom, snpcount, snpinfo, params, gblock, filters, masked_indivs, phenotypes_raw, sout);

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;
}


void readChunkFromBGENFileToG(const int bs, const int chrom, const uint32_t snpcount, vector<snp>& snpinfo, struct param* params, struct geno_block* gblock, struct filter* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, mstream& sout) {

  int ns, hc_val, nmales;
  uint32_t index ;
  double ds, total, mac, info_num;
  bool switch_alleles;
  std::string chromosome, rsid;
  uint32_t position ;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;

  for(int snp = 0; snp < bs; snp++) {

    // set to correct position
    gblock->bgen.jumpto( snpinfo[ snpcount + snp ].offset );
    gblock->bgen.read_variant( &chromosome, &position, &rsid, &alleles );

    //sout << "["<< chrom << "]SNPid stored ("<< snpinfo[snpcount+bs].chrom <<") = " << snpinfo[snpcount+bs].ID<< "/ SNPIDread ("<<chromosome<<")= " << rsid << endl; exit 1;

    assert(chrStrToInt(chromosome, params->nChrom) == chrom);
    gblock->bgen.read_probs( &probs ) ;

    ns = 0, hc_val = 0, index = 0, nmales = 0;
    total = 0, mac = 0, info_num = 0;
    for( std::size_t i = 0; i < probs.size(); ++i ) {

      // skip samples that were ignored from the analysis
      if( filters->ind_ignore(i) ) continue;

      ds = 0;
      for( std::size_t j = 1; j < probs[i].size(); ++j ) ds += probs[i][j] * j;

      if(ds != -3) {
        ds = params->ref_first ? ds : (2 - ds); // if ref-first, no need to switch

        if( filters->ind_in_analysis(index) ){
          if( !params->strict_mode || (params->strict_mode && masked_indivs(index,0)) ){
            total += ds;
            // compute MAC using 0.5*g for males for variants on sex chr (males coded as diploid)
            // sex is 1 for males and 0 o.w.
            if(params->test_mode && (chrom == params->nChrom)) {
              mac +=  ds * 0.5 * (2 - params->sex[i]);
              if(params->sex[i]) nmales++;
            }

            if( params->ref_first )
              info_num += 4 * probs[i][2] + probs[i][1] - ds * ds;
            else
              info_num += 4 * probs[i][0] + probs[i][1] - ds * ds;

            ns++;
          }

          // get genotype counts (convert to hardcall)
          if( params->htp_out ) {
            hc_val = (int) (ds + 0.5); // round to nearest integer (0/1/2)
            update_genocounts(params->binary_mode, index, hc_val, gblock->genocounts[snp], masked_indivs, phenotypes_raw);
          }
        }
      }

      gblock->Gmat(snp, index) = ds;
      index++;
    }

    if( params->test_mode){
      if(chrom != params->nChrom) {
        mac = total; // use MAC assuming diploid coding
        mac = min( mac, 2 * ns - mac );
      } else mac = min(mac, 2 * ns - nmales - mac); // males are 0/1

      if(mac < params->min_MAC) { 
        gblock->bad_snps(snp) = true;
      }
      if( params->htp_out ) gblock->snp_mac(snp,0) = mac;
    }

    //sout << "SNP#" << snp + 1 << "AC=" << mac << " BAD="<< (bad_snps(snp)?"BAD":"GOOD")<< endl;
    total /= ns;
    if( (params->alpha_prior != -1) || params->test_mode) gblock->snp_afs(snp, 0) = total / 2;

    if(params->test_mode) {
      if( (gblock->snp_afs(snp, 0) == 0) || (gblock->snp_afs(snp, 0) == 1) ) gblock->snp_info(snp, 0) = 1;
      else gblock->snp_info(snp, 0) = 1 - info_num / (2 * ns * gblock->snp_afs(snp, 0) * (1 - gblock->snp_afs(snp, 0)));

      if( params->setMinINFO && ( gblock->snp_info(snp, 0) < params->min_INFO) ) {
        gblock->bad_snps(snp) = true;
      }
    }

    if(params->use_SPA) {
      // switch to minor allele
      switch_alleles = (params->test_type > 0) ? false : (total > 1); // skip for DOM/REC test

      if(switch_alleles){
        gblock->Gmat.row(snp).array() = ( gblock->Gmat.row(snp).array() != -3).select( 2 - gblock->Gmat.row(snp).array(), gblock->Gmat.row(snp).array() );
        total = 2 - total;
        gblock->snp_flipped[snp] = true;
      } else gblock->snp_flipped[snp] = false;
    }

    // apply dominant/recessive encoding & recompute mean
    if(params->test_type > 0){
      index = 0;
      for( std::size_t i = 0; i < probs.size(); ++i ) {
        // skip samples that were ignored from the analysis
        if( filters->ind_ignore(i) ) continue;

        if( (gblock->Gmat(snp, index) != -3)  && filters->ind_in_analysis(index) &&
            (!params->strict_mode || (params->strict_mode && masked_indivs(index,0))) ){
          if(params->test_type == 1){ //dominant
            gblock->Gmat(snp, index) = params->ref_first ? (probs[i][1] + probs[i][2]) : (probs[i][0] + probs[i][1]);
          } else if(params->test_type == 2){ //recessive
            gblock->Gmat(snp, index) = params->ref_first ? probs[i][2] : probs[i][0];
          }
        }
        index++;
      }
      total = ((gblock->Gmat.row(snp).transpose().array()!= -3) && filters->ind_in_analysis).select(gblock->Gmat.row(snp).transpose().array(), 0).sum() / ns;
      if(total < params->numtol) gblock->bad_snps(snp) = true;
    }

    // deal with missing data and center SNPs
    for( std::size_t i = 0; i < params->n_samples; ++i ) {
      ds = gblock->Gmat(snp, i);
      if( params->use_SPA && filters->ind_in_analysis(i) && ds > 0 ) gblock->non_zero_indices_G[snp].push_back(i);

      // impute missing
      mean_impute_g(gblock->Gmat(snp, i), total, filters->ind_in_analysis(i), masked_indivs(i,0), params->strict_mode);

    }

  }

  if(!params->verbose) sout << bs << " snps ";
}


void readChunkFromBedFileToG(const int bs, const int chrom, const uint32_t snpcount, vector<snp>& snpinfo, struct param* params, struct in_files* files, struct geno_block* gblock, struct filter* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, mstream& sout) {

  int hc, ns, byte_start, bit_start, nmales;
  uint32_t index ;
  double total, mac;
  bool switch_alleles;
  // mapping matches the switch of alleles done when reading bim
  const int maptogeno[4] = {2, -3, 1, 0};

  // only for step 1
  for(int j = 0; j < bs; j++) {

    ns = 0, total = 0, mac = 0, index = 0, nmales = 0;
    
    // set to correct position
    jumpto_bed(snpinfo[snpcount + j].offset, files);
    files->bed_ifstream.read( reinterpret_cast<char *> (&files->inbed[0]), files->bed_block_size);

    for (int i = 0; i < filters->ind_ignore.size(); i++) {

      // skip samples that were ignored from the analysis
      if( filters->ind_ignore(i) ) continue;

      byte_start = i>>2; // 4 samples per byte
      bit_start = (i&3)<<1; // 2 bits per sample
      hc = maptogeno[ (files->inbed[byte_start] >> bit_start)&3 ];
      if(params->ref_first && (hc != -3)) hc = 2 - hc;
      gblock->Gmat(j, index) = hc;

      if(hc != -3) {
        if( filters->ind_in_analysis(index) ){
          if( !params->strict_mode || (params->strict_mode && masked_indivs(index,0)) ){
            total += hc;
            // compute MAC using 0.5*g for males for variants on sex chr (males coded as diploid)
            // sex is 1 for males and 0 o.w.
            if(params->test_mode && (chrom == params->nChrom)) {
              mac +=  hc * 0.5 * (2 - params->sex[i]);
              if(params->sex[i]) nmales++;
            }
            ns++;
          }

          // get genotype counts
          if( params->htp_out ) 
            update_genocounts(params->binary_mode, index, hc, gblock->genocounts[j], masked_indivs, phenotypes_raw);

        }
      }
      index++;
    }

    if( params->test_mode){
      if(chrom != params->nChrom) {
        mac = total; // use MAC assuming diploid coding
        mac = min( mac, 2 * ns - mac );
      } else mac = min(mac, 2 * ns - nmales - mac); // males are 0/1

      if(mac < params->min_MAC) { 
        gblock->bad_snps(j) = true; 
      }
      if( params->htp_out ) gblock->snp_mac(j,0) = mac;
    }

    total /= ns;
    if((params->alpha_prior != -1) || params->test_mode) gblock->snp_afs(j, 0) = total / 2;

    if(params->use_SPA) {
      // switch to minor allele
      switch_alleles = (params->test_type > 0) ? false : (total > 1); // skip for DOM/REC test

      if(switch_alleles){
        gblock->Gmat.row(j).array() = ( gblock->Gmat.row(j).array() != -3).select( 2 - gblock->Gmat.row(j).array(), gblock->Gmat.row(j).array() );
        total = 2 - total;
        gblock->snp_flipped[j] = true;
      } else gblock->snp_flipped[j] = false;
    }

    // apply dominant/recessive encoding & recompute mean
    if(params->test_type > 0){
      if(params->test_type == 1){ //dominant
        gblock->Gmat.row(j).array() = (gblock->Gmat.row(j).array() == 2).select(1, gblock->Gmat.row(j).array());
      } else if(params->test_type == 2){ //recessive
        gblock->Gmat.row(j).array() = (gblock->Gmat.row(j).array() >= 1).select(gblock->Gmat.row(j).array() - 1, gblock->Gmat.row(j).array());
      }
      total = ((gblock->Gmat.row(j).transpose().array() != -3) && filters->ind_in_analysis).select(gblock->Gmat.row(j).transpose().array(), 0).sum() / ns;
      if(total < params->numtol) gblock->bad_snps(j) = true;
    }

    //if(j<5) sout << "\nj="<< j+1 << ":" <<  gblock->Gmat.row(j).array().head(5);
    // deal with missing data and center SNPs
    for (size_t i = 0; i < params->n_samples; i++) {
      hc = gblock->Gmat(j, i);
      if( params->use_SPA && (hc > 0) ) gblock->non_zero_indices_G[j].push_back(i);

      // impute missing
      mean_impute_g(gblock->Gmat(j, i), total, filters->ind_in_analysis(i), masked_indivs(i,0), params->strict_mode);

    }

  }

  sout << bs << " snps ";

}


// only for step 1
void readChunkFromPGENFileToG(const int bs, const uint32_t snpcount, vector<snp>& snpinfo, struct param* params, struct geno_block* gblock, struct filter* filters, const Ref<const MatrixXb>& masked_indivs, mstream& sout) {

  int ns;
  double total;

  for(int j = 0; j < bs; j++) {

    // read genotype data
    if( params->dosage_mode )
      gblock->pgr.Read(gblock->genobuf, snpinfo[snpcount+j].offset, 1);
    else
      gblock->pgr.ReadHardcalls(gblock->genobuf, snpinfo[snpcount+j].offset, 1);

    ns = 0, total = 0;
    for (size_t i = 0; i < params->n_samples; i++) {

      gblock->Gmat(j, i) = gblock->genobuf[i];

      if(gblock->genobuf[i] != -3.0) {
        if( filters->ind_in_analysis(i) ){
          if( !params->strict_mode || (params->strict_mode && masked_indivs(i,0)) ){
            total += gblock->genobuf[i];
            ns++;
          }
        }
      }
    }

    total /= ns;
    if( params->alpha_prior != -1) gblock->snp_afs(j, 0) = total / 2;

    //if(j<5) sout << "\nj="<< j+1 << ":" <<  gblock->Gmat.row(j).array().head(5);
    // deal with missing data and center SNPs
    for (size_t i = 0; i < params->n_samples; i++) {

      // impute missing
      mean_impute_g(gblock->Gmat(j, i), total, filters->ind_in_analysis(i), masked_indivs(i,0), params->strict_mode);

    }

  }

  sout << bs << " snps ";
}


// check if uses Layout 2 (v1.2/1.3) and compressed using zlib's compress() function & check for first SNP if precision for probabilities is 8 bits
void check_bgen(const string bgen_file, struct param* params){

  // for non-bgen file input, skip check
  if(params->file_type != "bgen") return;

  BgenParser bgen_ck;
  bgen_ck.open( bgen_file ) ;
  bool layoutV2 = bgen_ck.get_layout();
  bool zlib_compress = bgen_ck.get_compression();
  if( !layoutV2 || !zlib_compress ){
    params->streamBGEN = false;
    return;
  }
  uint64 first_snp = bgen_ck.get_position();

  uint minploidy = 0, maxploidy = 0, phasing = 0, bits_prob = 0;
  uint16_t SNPID_size = 0, RSID_size = 0, chromosome_size = 0 , numberOfAlleles = 0 ;
  uint32_t position = 0, allele_size = 0, nindivs = 0;
  string allele, tmp_buffer;

  // check bits only for first snp
  //cout << endl << "Snp1 pos:" << first_snp << endl;
  ifstream bfile;
  bfile.open( bgen_file, ios::in | ios::binary );
  bfile.seekg( first_snp );
  // snpid
  bfile.read( reinterpret_cast<char *> (&SNPID_size), 2 );
  tmp_buffer.resize(SNPID_size);
  bfile.read( reinterpret_cast<char *> (&tmp_buffer[0]), SNPID_size );
  // rsid
  bfile.read( reinterpret_cast<char *> (&RSID_size), 2) ;
  tmp_buffer.resize(RSID_size);
  bfile.read( reinterpret_cast<char *> (&tmp_buffer[0]), RSID_size );
  //cout << "RSID:" << tmp_buffer ;
  // chromosome
  bfile.read( reinterpret_cast<char *> (&chromosome_size), 2 );
  tmp_buffer.resize(chromosome_size);
  bfile.read( reinterpret_cast<char *> (&tmp_buffer[0]), chromosome_size );
  assert( chrStrToInt(tmp_buffer , params->nChrom) > 0 );
  //cout << ",CHR:" << tmp_buffer ;
  // position
  bfile.read( reinterpret_cast<char *> (&position), 4 );
  //cout << ",POS:" << position << endl;
  // number of alleles
  bfile.read( reinterpret_cast<char *> (&numberOfAlleles), 2 );
  assert( numberOfAlleles == 2 ); // only diploid
  //cout << ",Nalleles:" << numberOfAlleles ;
  // alleles
  bfile.read( reinterpret_cast<char *> (&allele_size), 4 );
  tmp_buffer.resize(allele_size);
  bfile.read( reinterpret_cast<char *> (&tmp_buffer[0]), allele_size );
  //cout << ",A0:"<<tmp_buffer ;
  bfile.read( reinterpret_cast<char *> (&allele_size), 4 );
  tmp_buffer.resize(allele_size);
  bfile.read( reinterpret_cast<char *> (&tmp_buffer[0]), allele_size );
  //cout << ",A1:"<<tmp_buffer ;

  // set genotype data block
  vector < uchar > geno_block, geno_block_uncompressed;
  uint32_t size_block = 0, size_block_post_compression = 0;
  bfile.read( reinterpret_cast<char *> (&size_block), 4 );
  bfile.read( reinterpret_cast<char *> (&size_block_post_compression), 4);
  //cout << ",block size:"<<size_block  << ",block size post compress:" << size_block_post_compression << endl;
  geno_block.resize(size_block - 4);
  geno_block_uncompressed.resize(size_block_post_compression);
  bfile.read( reinterpret_cast<char *> (&geno_block[0]), size_block - 4);

  // uncompress the block using zlib
  uLongf dest_size = size_block_post_compression;
  if( (uncompress( &(geno_block_uncompressed[0]), &dest_size, &geno_block[0], size_block - 4) != Z_OK) || (dest_size != size_block_post_compression) ){
    params->streamBGEN = false;
    return;
  }

  // stream to uncompressed block
  uchar *buffer = &geno_block_uncompressed[0];
  // sample size
  std::memcpy(&nindivs, &(buffer[0]), 4);
  //cout << "N:"<< nindivs ;
  assert( ((int) nindivs) == bgen_ck.number_of_samples() );
  buffer += 4;
  // num alleles
  std::memcpy(&numberOfAlleles, &(buffer[0]), 2);
  //cout << ",allele:"<< numberOfAlleles ;
  assert( numberOfAlleles == 2 );
  buffer += 2;
  // ploidy
  std::memcpy(&minploidy, &(buffer[0]), 1);
  //cout << ",minP:"<< minploidy ;
  assert( minploidy == 2 );
  buffer ++;
  std::memcpy(&maxploidy, &(buffer[0]), 1);
  //cout << ",maxP:"<< maxploidy ;
  assert( maxploidy == 2 );
  buffer ++;

  /* //to identify missing when getting dosages
     vector < uchar > ploidy_n;
     ploidy_n.resize( nindivs );
     std::memcpy(&(ploidy_n[0]), &(buffer[0]), nindivs);
     */
  buffer += nindivs;

  // phasing
  std::memcpy(&phasing, &(buffer[0]), 1);
  //cout << ",phasing:"<< phasing ;
  buffer ++;

  // bits per probability
  std::memcpy(&bits_prob, &(buffer[0]), 1);
  //cout << ",bits:"<< bits_prob ;
  buffer ++;
  if( (phasing != 0) || (bits_prob != 8) ){
    params->streamBGEN = false;
    return;
  }

}


// for step 2 (using MT in openmp)
void readChunkFromBGEN(std::istream* bfile, uint32_t* size1, uint32_t* size2, vector<uchar>* geno_block){

  uint16_t SNPID_size = 0, RSID_size = 0, chromosome_size = 0 , numberOfAlleles = 0 ;
  uint32_t position = 0, allele_size = 0;
  string tmp_buffer;

  // snpid
  bfile->read( reinterpret_cast<char *> (&SNPID_size), 2 );
  tmp_buffer.resize(SNPID_size);
  bfile->read( reinterpret_cast<char *> (&tmp_buffer[0]), SNPID_size );
  // rsid
  bfile->read( reinterpret_cast<char *> (&RSID_size), 2) ;
  tmp_buffer.resize(RSID_size);
  bfile->read( reinterpret_cast<char *> (&tmp_buffer[0]), RSID_size );
  // chromosome
  bfile->read( reinterpret_cast<char *> (&chromosome_size), 2 );
  tmp_buffer.resize(chromosome_size);
  bfile->read( reinterpret_cast<char *> (&tmp_buffer[0]), chromosome_size );
  // position
  bfile->read( reinterpret_cast<char *> (&position), 4 );
  // number of alleles
  bfile->read( reinterpret_cast<char *> (&numberOfAlleles), 2 );
  // alleles
  bfile->read( reinterpret_cast<char *> (&allele_size), 4 );
  tmp_buffer.resize(allele_size);
  bfile->read( reinterpret_cast<char *> (&tmp_buffer[0]), allele_size );
  bfile->read( reinterpret_cast<char *> (&allele_size), 4 );
  tmp_buffer.resize(allele_size);
  bfile->read( reinterpret_cast<char *> (&tmp_buffer[0]), allele_size );

  // set genotype data block
  bfile->read( reinterpret_cast<char *> (size1), 4 );
  bfile->read( reinterpret_cast<char *> (size2), 4);
  geno_block->resize(*size1 - 4);
  bfile->read( reinterpret_cast<char *> (&((*geno_block)[0])), *size1 - 4);

  return;
}


void parseSnpfromBGEN(const int isnp, const int &chrom, vector<uchar>* geno_block, const uint32_t insize, const uint32_t outsize, const struct param* params, const struct filter* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, const snp* infosnp, struct geno_block* gblock, variant_block* snp_data, mstream& sout){

  uint minploidy = 0, maxploidy = 0, phasing = 0, bits_prob = 0;
  uint16_t numberOfAlleles = 0 ;
  uint32_t nindivs = 0;
  uint32_t index;
  string tmp_buffer;

  MapArXd Geno (gblock->Gmat.col(isnp).data(), params->n_samples, 1);
  // reset variant info
  snp_data->ignored = false;
  snp_data->fastSPA = params->use_SPA;
  snp_data->n_non_zero = 0;

  // set genotype data block
  vector < uchar > geno_block_uncompressed;
  geno_block_uncompressed.resize(outsize);

  // uncompress the block using zlib
  uLongf dest_size = outsize;
  if( (uncompress( &(geno_block_uncompressed[0]), &dest_size, &((*geno_block)[0]), insize) != Z_OK) || (dest_size != outsize) ){
    sout << "ERROR: Failed to decompress genotype data block for variant: " << infosnp->ID << endl;
    exit(EXIT_FAILURE);
  }

  // stream to uncompressed block
  uchar *buffer = &geno_block_uncompressed[0];
  // sample size in file
  std::memcpy(&nindivs, &(buffer[0]), 4);
  assert( nindivs == filters->ind_ignore.size() );
  buffer += 4;
  // num alleles
  std::memcpy(&numberOfAlleles, &(buffer[0]), 2);
  assert( numberOfAlleles == 2 );
  buffer += 2;
  // ploidy
  std::memcpy(&minploidy, &(buffer[0]), 1);
  assert( minploidy == 2 );
  buffer ++;
  std::memcpy(&maxploidy, &(buffer[0]), 1);
  assert( maxploidy == 2 );
  buffer ++;

  //to identify missing when getting dosages
  vector < uchar > ploidy_n;
  ploidy_n.resize( nindivs );
  std::memcpy(&(ploidy_n[0]), &(buffer[0]), nindivs);
  buffer += nindivs;

  // phasing
  std::memcpy(&phasing, &(buffer[0]), 1);
  buffer++;

  // bits per probability
  std::memcpy(&bits_prob, &(buffer[0]), 1);
  buffer++;

  // get dosages (can compute mean as going along (and identify non-zero entries if SPA is used)
  bool missing;
  int ns = 0, hc_val, nmales = 0;
  double prob0, prob1, prob2, total = 0, mac = 0, ds, info_num = 0;
  Geno = ArrayXd::Zero(params->n_samples);
  snp_data->genocounts = MatrixXd::Zero(6, params->n_pheno);

  // parse genotype probabilities block
  index = 0;
  for(size_t i = 0; i < nindivs; i++) {

    // skip samples that were ignored from the analysis
    if( filters->ind_ignore(i) ) {
      buffer+=2;
      continue;
    }

    missing = ((ploidy_n[i]) & 0x80);
    if(missing) {
      // bug fix (with imputed data this case should not occur)
      Geno(index++) = -3;
      buffer+=2;
      continue;
    }

    prob0 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;
    prob1 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;
    prob2 = std::max( 1 - prob0 - prob1, 0.0);

    Geno(index) = prob1 + 2 * prob2;
    if(!params->ref_first) Geno(index) = 2 - Geno(index); // switch allele0 to ALT

    if( filters->ind_in_analysis(index) ){
      if( !params->strict_mode || (params->strict_mode && masked_indivs(index,0)) ){
        total += Geno(index);
        // compute MAC using 0.5*g for males for variants on sex chr (males coded as diploid)
        // sex is 1 for males and 0 o.w.
        if(params->test_mode && (chrom == params->nChrom)) {
          mac +=  Geno(index) * 0.5 * (2 - params->sex[i]);
          if(params->sex[i]) nmales++;
        }

        if( params->ref_first )
          info_num += 4 * prob2 + prob1 - Geno(index) * Geno(index);
        else
          info_num += 4 * prob0 + prob1 - Geno(index) * Geno(index);

        ns++;
      }

      // get genotype counts (convert to hardcall)
      if( params->htp_out ) {
        hc_val = (int) (Geno(index) + 0.5); // round to nearest integer 0/1/2
        update_genocounts(params->binary_mode, index, hc_val, snp_data->genocounts, masked_indivs, phenotypes_raw);
      }

    }
    index++;
  }

  // check MAC
  if( params->test_mode){
    if(chrom != params->nChrom) mac = total; // use MAC assuming diploid coding


    if(chrom != params->nChrom) mac = min( mac, 2 * ns - mac );
    else mac = min(mac, 2 * ns - nmales - mac); // males are 0/1

    if(mac < params->min_MAC) { 
      snp_data->ignored = true;return;
    }
    if( params->htp_out ) snp_data->mac = mac;
  }

  total /= ns;
  snp_data->af = total / 2;
  // impute info score
  if( snp_data->af == 0 || snp_data->af == 1 )
    snp_data->info = 1;
  else
    snp_data->info = 1 - info_num / (2 * ns * snp_data->af * (1 - snp_data->af));

  // check INFO score
  if( params->setMinINFO && ( snp_data->info < params->min_INFO) ) {
    snp_data->ignored = true;
    return;
  }

  if(params->use_SPA) {
    // switch to minor allele
    snp_data->flipped = (params->test_type > 0) ? false : (total > 1); // skip for DOM/REC test

    if(snp_data->flipped){
      Geno = ( Geno != -3).select( 2 - Geno, Geno );
      total = 2 - total;
    }
  }

  // apply dominant/recessive encoding & recompute mean
  if(params->test_type > 0){
    // go over data block again
    buffer -= 2 * nindivs;
    index = 0;
    for(size_t i = 0; i < nindivs; i++) {

      // skip samples that were ignored from the analysis
      if( filters->ind_ignore(i) ) {
        buffer+=2;
        continue;
      }

      missing = ((ploidy_n[i]) & 0x80);
      if(missing) {
        buffer+=2;
        continue;
      }
      prob0 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;
      prob1 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;
      prob2 = std::max( 1 - prob0 - prob1, 0.0);

      if(filters->ind_in_analysis(index)){
        if(params->test_type == 1){ //dominant
          Geno(index) = params->ref_first ? (prob1 + prob2) : (prob0 + prob1);
        } else if(params->test_type == 2){ //recessive
          Geno(index) = params->ref_first ? prob2 : prob0;
        }
      }
      index++;
    }
    total = ((Geno != -3) && filters->ind_in_analysis).select(Geno, 0).sum() / ns;
    if(total < params->numtol) {
      snp_data->ignored = true;
      return;
    }
  }

  // deal with missing data & prep for spa
  for( size_t i = 0; i < params->n_samples; ++i ) {
    ds = Geno(i);

    // keep track of number of entries filled so avoid using clear
    if( params->use_SPA && (snp_data->fastSPA) && filters->ind_in_analysis(i) && ds > 0 ) 
      update_nnz_spa(i, params->n_samples, snp_data);

    // impute missing
    mean_impute_g(Geno(i), total, filters->ind_in_analysis(i), masked_indivs(i,0), params->strict_mode);
  }

  return;
}


void parseSnpfromBed(const int isnp, const int &chrom, const vector<uchar> geno_block, const struct param* params, const struct filter* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, struct geno_block* gblock, variant_block* snp_data){

  int hc, ns, byte_start, bit_start, nmales;
  uint32_t index ;
  double total, mac;
  // mapping matches the switch of alleles done when reading bim
  const int maptogeno[4] = {2, -3, 1, 0};

  MapArXd Geno (gblock->Gmat.col(isnp).data(), params->n_samples, 1);
  Geno = ArrayXd::Zero(params->n_samples);
  snp_data->genocounts = MatrixXd::Zero(6, params->n_pheno);
  // reset variant info
  snp_data->ignored = false;
  snp_data->fastSPA = params->use_SPA;
  snp_data->n_non_zero = 0;

  ns = 0, total = 0, mac = 0, index = 0, nmales = 0;
  for (int i = 0; i < filters->ind_ignore.size(); i++) {

    // skip samples that were ignored from the analysis
    if( filters->ind_ignore(i) ) continue;

    byte_start = i>>2; // 4 samples per byte
    bit_start = (i&3)<<1; // 2 bits per sample
    hc = maptogeno[ (geno_block[byte_start] >> bit_start)&3 ];
    if(params->ref_first && (hc != -3)) hc = 2 - hc;
    Geno(index) = hc;

    if(hc != -3) {
      if( filters->ind_in_analysis(index) ){
        if( !params->strict_mode || (params->strict_mode && masked_indivs(index,0)) ){
          total += hc;
          // compute MAC using 0.5*g for males for variants on sex chr (males coded as diploid)
          // sex is 1 for males and 0 o.w.
          if(params->test_mode && (chrom == params->nChrom)) {
            mac +=  hc * 0.5 * (2 - params->sex[i]);
            if(params->sex[i]) nmales++;
          }
          ns++;
        }

        // get genotype counts
        if( params->htp_out ) 
          update_genocounts(params->binary_mode, index, hc, snp_data->genocounts, masked_indivs, phenotypes_raw);

      }
    }
    index++;
  }

  // check MAC
  if( params->test_mode){
    if(chrom != params->nChrom) mac = total; // use MAC assuming diploid coding

    if(chrom != params->nChrom) mac = min( mac, 2 * ns - mac );
    else mac = min(mac, 2 * ns - nmales - mac); // males are 0/1

    if(mac < params->min_MAC) { 
      snp_data->ignored = true;
      return;
    }
    if( params->htp_out ) snp_data->mac = mac;
  }

  total /= ns;
  snp_data->af = total / 2;

  if(params->use_SPA) {
    // switch to minor allele
    snp_data->flipped = (params->test_type > 0) ? false : (total > 1); // skip for DOM/REC test

    if(snp_data->flipped){
      Geno = ( Geno != -3).select( 2 - Geno, Geno);
      total = 2 - total;
    }
  }

  // apply dominant/recessive encoding & recompute mean
  if(params->test_type > 0){
    if(params->test_type == 1){ //dominant
      Geno = (Geno == 2).select(1, Geno);
    } else if(params->test_type == 2){ //recessive
      Geno = (Geno >= 1).select(Geno - 1, Geno);
    }
    total = ((Geno != -3) && filters->ind_in_analysis).select(Geno, 0).sum() / ns;
    if(total < params->numtol) {
      snp_data->ignored = true;
      return;
    }
  }


  // deal with missing data & prep for spa
  for( size_t i = 0; i < params->n_samples; ++i ) {
    hc = Geno(i);

    // keep track of number of entries filled so avoid using clear
    if( params->use_SPA && (snp_data->fastSPA) && filters->ind_in_analysis(i) && hc > 0 ) 
      update_nnz_spa(i, params->n_samples, snp_data);

    // impute missing
    mean_impute_g(Geno(i), total, filters->ind_in_analysis(i), masked_indivs(i,0), params->strict_mode);

  }

}


// step 2
void readChunkFromPGENFileToG(const int &start, const int &bs, const int &chrom, struct param* params, struct filter* filters, struct geno_block* gblock, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, const vector<snp>& snpinfo, vector<variant_block> &all_snps_info){

  int hc, ns, index, nmales;
  double ds, total, mac, eij2 = 0;
  int cur_index;

  for(int j = 0; j < bs; j++) {

    variant_block* snp_data = &(all_snps_info[j]);
    MapArXd Geno (gblock->Gmat.col(j).data(), params->n_samples, 1);
    // reset variant info
    Geno = ArrayXd::Zero(params->n_samples);
    snp_data->genocounts = MatrixXd::Zero(6, params->n_pheno);
    snp_data->ignored = false;
    snp_data->fastSPA = params->use_SPA;
    snp_data->n_non_zero = 0;

    ns = 0, total = 0, mac = 0, index = 0, nmales = 0;
    if( params->dosage_mode ) eij2 = 0;

    // read genotype data
    cur_index = snpinfo[ start + j ].offset;
    if( params->dosage_mode )
      gblock->pgr.Read(gblock->genobuf, cur_index, 1);
    else
      gblock->pgr.ReadHardcalls(gblock->genobuf, cur_index, 1);


    for (int i = 0; i < filters->ind_ignore.size(); i++) {

      // skip samples that were ignored from the analysis
      if( filters->ind_ignore(i) ) continue;

      Geno(index) = gblock->genobuf[index];

      if(gblock->genobuf[index] != -3.0) {
        if( filters->ind_in_analysis(index) ){
          if( !params->strict_mode || (params->strict_mode && masked_indivs(index,0)) ){
            total += gblock->genobuf[index];
            // compute MAC using 0.5*g for males for variants on sex chr (males coded as diploid)
            // sex is 1 for males and 0 o.w.
            if(params->test_mode && (chrom == params->nChrom)) {
              mac +=  gblock->genobuf[index] * 0.5 * (2 - params->sex[i]);
              if(params->sex[i]) nmales++;
            }
            if( params->dosage_mode ) eij2 += gblock->genobuf[index] * gblock->genobuf[index];
            ns++;
          }

          // get genotype counts
          if( params->htp_out ) {
            hc = (int) (Geno(index) + 0.5); // round to nearest integer 0/1/2
            update_genocounts(params->binary_mode, index, hc, snp_data->genocounts, masked_indivs, phenotypes_raw);
          }

        }
      }
      index++;
    }

    // check MAC
    if( params->test_mode){
      if(chrom != params->nChrom) {
        mac = total; // use MAC assuming diploid coding
        mac = min( mac, 2 * ns - mac );
      } else mac = min(mac, 2 * ns - nmales - mac); // males are 0/1

      if(mac < params->min_MAC) { 
        snp_data->ignored = true; 
        continue;
      }
      if( params->htp_out ) snp_data->mac = mac;
    }

    total /= ns;
    snp_data->af = total / 2;

    // mach r2 info score
    if( params->dosage_mode ){
      if( snp_data->af == 0 || snp_data->af == 1 )
        snp_data->info = 1;
      else
        snp_data->info = (eij2 / ns - total * total) / (2 * snp_data->af * (1 - snp_data->af));

      // check INFO score
      if( params->setMinINFO && ( snp_data->info < params->min_INFO) ) {
        snp_data->ignored = true;
        continue;
      }

    }

    if(params->use_SPA) {
      // switch to minor allele
      snp_data->flipped = total > 1;
      if( params->test_type > 0) snp_data->flipped = false; // skip for DOM/REC test
      if(snp_data->flipped){
        Geno = ( Geno != -3.0).select( 2 - Geno, Geno);
        total = 2 - total;
      }
    }

    // apply dominant/recessive encoding & recompute mean
    // pgen does not contain genotype probs for dosages so convert to hardcalls
    if(params->test_type > 0){
      for( size_t i = 0; i < params->n_samples; ++i ) {
        if( (Geno(i) == -3.0) || !filters->ind_in_analysis(i) ) continue;
        hc = (int) (Geno(i) + 0.5);

        if(params->test_type == 1){ //dominant
          Geno(i) = (hc == 2 ? 1 : hc);
        } else if(params->test_type == 2){ //recessive
          Geno(i) = (hc >= 1 ? hc - 1 : hc);
        }
      }
      total = ((Geno != -3) && filters->ind_in_analysis).select(Geno, 0).sum() / ns;
      if( params->test_mode && (total < params->numtol) ) {
        snp_data->ignored = true;
        continue;
      }
    }

    // deal with missing data & prep for spa
    for( size_t i = 0; i < params->n_samples; ++i ) {
      ds = Geno(i);

      // keep track of number of entries filled so avoid using clear
      if( params->use_SPA && (snp_data->fastSPA) && filters->ind_in_analysis(i) && ds > 0 ) 
        update_nnz_spa(i, params->n_samples, snp_data);

      // impute missing
      mean_impute_g(Geno(i), total, filters->ind_in_analysis(i), masked_indivs(i,0), params->strict_mode);

    }

  }

}

bool in_chrList(const int snp_chr, struct filter* filters){

  return ( filters->chrKeep_test.find(snp_chr)!= filters->chrKeep_test.end() );
}

string bgi_chrList(struct filter* filters){

  int nchr_kept = filters->chrKeep_test.size();
  std::ostringstream buffer;
  map<int, bool >::iterator itr;

  int ic = 0;
  for (itr = filters->chrKeep_test.begin(); itr != filters->chrKeep_test.end(); ++itr) {
    buffer << "'" << to_string(itr->first) << "'";
    if(++ic < nchr_kept) buffer << ",";
  }

  return buffer.str();
}

bool in_range(string snp_chr, uint32_t snp_pos, struct param* params){

  if( snp_chr != params->range_chr ) return false;
  else if ( snp_pos < params->range_min ) return false;
  else if( snp_pos > params->range_max ) return false;

  return true; 
}


void skip_snps(uint64 offset, struct param* params, struct in_files* files, struct geno_block* gblock){

  // set to new position based on offset
  if(params->file_type == "bed") jumpto_bed(offset, files);
  else if(params->file_type == "bgen") gblock->bgen.jumpto(offset);

}

// jump to given snp index in bed file (+magic number)
void jumpto_bed(uint64 offset, struct in_files* files){
    files->bed_ifstream.seekg( 3 + offset * files->bed_block_size, ios_base::beg);
}

void update_genocounts(bool binary_mode, int ind, int hc, MatrixXd& genocounts, const Ref<const MatrixXb>& mask, const Ref<const MatrixXd>& ymat){

  if( !binary_mode ) {
    genocounts.row(hc) += mask.row(ind).cast<double>();
  } else {
    genocounts.row(hc).array() += mask.row(ind).array().cast<double>() * ymat.row(ind).array();
    genocounts.row(3 + hc).array() += mask.row(ind).array().cast<double>() * (1 - ymat.row(ind).array());
  }

}


// to track non-zero entries
void update_nnz_spa(uint32_t ind, uint32_t nsamples, variant_block* snp_data){

  snp_data->n_non_zero++;
  if(snp_data->n_non_zero > (nsamples * 0.5)) {
    snp_data->fastSPA = false;
  } else {
    if(snp_data->non_zero_indices.size() < snp_data->n_non_zero)
      snp_data->non_zero_indices.push_back(ind);
    else
      snp_data->non_zero_indices[snp_data->n_non_zero - 1] = ind;
  }

}

// mean impute (only individuals who are not masked)
void mean_impute_g(double &geno, const double mu, const bool in_analysis, const bool mask, const bool strict_mode){

  if(geno != -3 && in_analysis && (!strict_mode || (strict_mode && mask) ) ){
    geno -= mu;
  } else {
    geno = 0;
  }

}

findID getIndivIndex(const string &FID, const string &IID, struct param* params, mstream& sout){

  string tmp_str;
  findID indiv;

  // get ID of individual
  tmp_str = FID + "_" + IID;

  // check individual is in genotype data
  indiv.is_found = ( params->FID_IID_to_ind.find(tmp_str) != params->FID_IID_to_ind.end() );

  if(indiv.is_found)
    indiv.index = params->FID_IID_to_ind[tmp_str];

  return indiv;
}

