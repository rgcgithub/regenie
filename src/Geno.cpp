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

#include "Regenie.hpp"
#include "Files.hpp"
#include "Geno.hpp"
#include "db/sqlite3.hpp"

using namespace std;
using namespace Eigen;
using namespace boost;



void prep_bgen(struct in_files* files, struct param* params, struct filter* filters, vector<snp>& snpinfo, map<int, vector<int>>& chr_map, BgenParser& bgen, mstream& sout){

  bool interaction_snp_found = false;
  uint32_t nOutofOrder = 0, lineread = 0;
  std::string chromosome, rsid, msg;
  uint32_t position ;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;
  std::vector< string > tmp_ids ;
  snp tmp_snp;
  BgenParser bgen_tmp;

  sout << left << std::setw(20) << " * bgen" << ": [" << files->bgen_file << "]" << endl;
  // open file and print file info
  bgen_tmp.open( files->bgen_file ) ;
  sout << bgen_tmp.summarise( ) << 
    // also add in the number of bits
    (params->BGENbits > 0 ? " with " + to_string(params->BGENbits) + "-bit encoding" : "") << ".\n";

  // get info for variants
  if( params->with_bgi ) read_bgi_file(bgen_tmp, files, params, filters, snpinfo, sout);
  else {
    tmp_snp.offset = bgen_tmp.get_position();
    while(bgen_tmp.read_variant( &chromosome, &position, &rsid, &alleles )) {

      assert(alleles.size() == 2) ; // only bi-allelic allowed
      // check phasing for first variant
      if(lineread == 0){
        bgen_tmp.read_probs( &probs ) ;

        if( probs[0].size() != 3 ) // unphased only 
          throw "only unphased bgen are supported.";

      } else bgen_tmp.ignore_probs();

      tmp_snp.chrom = chrStrToInt(chromosome, params->nChrom);
      if (tmp_snp.chrom == -1) 
        throw "unknown chromosome code in bgen file.";

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

      // if using GxG interaction test
      if(params->interaction_snp && (tmp_snp.ID == filters->interaction_cov)){
        if(!params->interaction_file) {
          params->interaction_snp_offset = tmp_snp.offset;
          params->ltco_chr = tmp_snp.chrom;
          interaction_snp_found = true;
        }
        // go to next variant (get its offset first)
        tmp_snp.offset = bgen_tmp.get_position();
        continue;
      }

      // if specified chrlist/range
      if(
          (params->select_chrs && !in_chrList(tmp_snp.chrom, filters))
          ||
          (params->set_range && !in_range(tmp_snp.chrom, position, params))
        ) {
        // go to next variant (get its offset first)
        tmp_snp.offset = bgen_tmp.get_position();
        continue;
      }

      // make list of variant IDs if inclusion/exclusion file is given
      if(params->mk_snp_map){
        if (in_map(tmp_snp.ID, filters->snpID_to_ind)) { // ignore duplicate
          tmp_snp.offset = bgen_tmp.get_position();
          continue;
        }
        filters->snpID_to_ind[ tmp_snp.ID ] = snpinfo.size();
      }

      // keep track of how many included snps per chromosome there are
      files->chr_counts[tmp_snp.chrom-1]++;

      snpinfo.push_back(tmp_snp);

      tmp_snp.offset = bgen_tmp.get_position();
    }

    if(params->interaction_snp && !params->interaction_file && !interaction_snp_found)
      throw "SNP specified for GxG interaction test was not found.";

    if (!params->test_mode && (nOutofOrder > 0)) 
      sout << "WARNING: Total number of snps out-of-order in bgen file : " << nOutofOrder << endl;
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
    // set to unknown
    params->sex = ArrayXi::Constant(params->n_samples, 0);
  }

  // check duplicates -- if not, store in map
  for(size_t i = 0; i < params->n_samples; i++) {

    if (in_map(tmp_ids[i], params->FID_IID_to_ind)) 
      throw "duplicate individual in bgen file : FID_IID =" + tmp_ids[i];

    params->FID_IID_to_ind[ tmp_ids[i] ] = i;
  }

  // check if should mask samples
  check_samples_include_exclude(files, params, filters, sout);

  // setup file for reading the genotype probabilities later
  if( !params->streamBGEN ) 
    bgen.open( files->bgen_file ) ;
  else
    openStream(&files->geno_ifstream, files->bgen_file, ios::in | ios::binary, sout);

  if (params->test_mode) params->dosage_mode = true;
}


// read .bgi file to get SNP info
void read_bgi_file(BgenParser& bgen, struct in_files* files, struct param* params, struct filter* filters, std::vector<snp>& snpinfo, mstream& sout){

  bool interaction_snp_found = false;
  int nalleles;
  uint32_t lineread = 0;
  uint64 variant_bgi_size, variant_bgen_size;
  string bgi_file = files->bgen_file + ".bgi";
  string sql_query = "SELECT * FROM Variant", cnd1 = "";
  snp tmp_snp;
  sqlite3* db;
  sqlite3_stmt* stmt;

  uint32_t n_variants = bgen.number_of_variants();
  uint32_t position ;
  std::string chromosome, rsid;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;

  // edit sql statement if chromosome position range is given
  if( params->set_range ){
    cnd1 = " WHERE ( chromosome IN (" + bgi_chrList(params->range_chr, params->nChrom) + ") AND position>=" + to_string(params->range_min) + " AND position<=" + to_string(params->range_max) + ")";
  } else if( params->select_chrs ){
    cnd1 = " WHERE ( chromosome IN (" + bgi_chrList(filters, params->nChrom) + " ) )";
  }
  // with GxG tests
  if(params->interaction_snp && (cnd1.size() > 0)){ // bug fix - only use this if querying on chrs/range
    cnd1.append(" OR ( rsid = '" + filters->interaction_cov + "' )" );
  }
  sql_query.append( cnd1 );

  sout << "   -index bgi file [" << bgi_file<< "]" << endl;
  if( sqlite3_open( bgi_file.c_str(), &db ) != SQLITE_OK ) 
    throw  sqlite3_errmsg(db);


  // header: chromosome|position|rsid|number_of_alleles|allele1|allele2|file_start_position|size_in_bytes
  if( sqlite3_prepare_v2( db, sql_query.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
    throw sqlite3_errmsg(db);

  bool done = false;
  uint32_t nOutofOrder = 0;
  while (!done) {
    switch (sqlite3_step(stmt)) {
      case SQLITE_ROW:

        chromosome = std::string( (char *) sqlite3_column_text(stmt, 0) );
        tmp_snp.chrom = chrStrToInt(chromosome, params->nChrom);
        if (tmp_snp.chrom == -1) 
          throw "unknown chromosome code in bgi file (=" + chromosome + ").";
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
          bgen.read_probs( &probs ) ;
          if( probs[0].size() != 3 ) // unphased only 
            throw "only unphased bgen are supported.";
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

        // if using GxG interaction test
        if(params->interaction_snp && (tmp_snp.ID == filters->interaction_cov)){
          if(!params->interaction_file) {
            params->interaction_snp_offset = tmp_snp.offset;
            params->ltco_chr = tmp_snp.chrom;
            interaction_snp_found = true;
          }
          continue; // don't save it
        }

        // make list of variant IDs if inclusion/exclusion file is given
        if(params->mk_snp_map){
          if (in_map(tmp_snp.ID, filters->snpID_to_ind))
            continue; // don't save it
          filters->snpID_to_ind[ tmp_snp.ID ] = snpinfo.size();
        }

        // keep track of how many included snps per chromosome there are
        files->chr_counts[tmp_snp.chrom-1]++;

        snpinfo.push_back(tmp_snp);
        break;

      case SQLITE_DONE:
        done = true;
        break;

      default:
        throw "failed reading file (" + std::string( sqlite3_errmsg(db) ) + ").";
    }
  }

  sqlite3_finalize(stmt);
  sqlite3_close(db);

  if(params->interaction_snp && !params->interaction_file && !interaction_snp_found)
    throw "SNP specified for GxG interaction test was not found.";

  if( !params->set_range && !params->select_chrs) assert( lineread == n_variants );
  if (!params->test_mode && (nOutofOrder > 0)) sout << "WARNING: Total number of snps out-of-order in bgen file : " << nOutofOrder << endl;

}

void read_bgi_file(string const& setting, BgenParser& bgen, geno_file_info* ext_file_info, map <string, uint64>* variant_names, struct param* params, mstream& sout){

  uint32_t lineread = 0;
  uint64 offset, variant_bgi_size, variant_bgen_size;
  string bgi_file = ext_file_info->file + ".bgi";
  string sql_query, rsid;
  map <string, uint64> tmp_map;
  sqlite3* db;
  sqlite3_stmt* stmt;

  int nalleles, chrom;
  uint32_t position ;
  std::string tmp_str, chr;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;

  // sql statement to pass variants to keep
  sql_query = "SELECT * FROM Variant WHERE rsid IN (" + bgi_rsidList((*variant_names)) + " )";

  sout << "      -index bgi file [" << bgi_file<< "]" << endl;
  if( sqlite3_open( bgi_file.c_str(), &db ) != SQLITE_OK ) 
    throw  sqlite3_errmsg(db);

  // header: chromosome|position|rsid|number_of_alleles|allele1|allele2|file_start_position|size_in_bytes
  if( sqlite3_prepare_v2( db, sql_query.c_str(), -1, &stmt, NULL ) != SQLITE_OK )
    throw sqlite3_errmsg(db);

  bool done = false;
  while (!done) {
    switch (sqlite3_step(stmt)) {
      case SQLITE_ROW:

        nalleles = atoi( (char *) sqlite3_column_text(stmt, 3) );
        assert(nalleles == 2) ; // only bi-allelic allowed

        rsid = std::string( (char *) sqlite3_column_text(stmt, 2) );
        offset = strtoull( (char *) sqlite3_column_text(stmt, 6), NULL, 10);
        // make list of variant IDs
        tmp_map[ rsid ] = offset;

        // check if matches with info from bgenparser for first read variant
        if( lineread++ == 0 ){
          bgen.jumpto(offset);
          bgen.read_variant( &chr, &position, &tmp_str, &alleles );
          bgen.read_probs( &probs ) ;
          if( probs[0].size() != 3 ) // unphased only 
            throw "only unphased bgen are supported.";
          variant_bgen_size = bgen.get_position() - offset;
          variant_bgi_size = strtoull( (char *) sqlite3_column_text(stmt, 7), NULL, 10);
          assert( tmp_str == rsid );
          assert( variant_bgi_size == variant_bgen_size );
          if(setting == "interaction") {
            chrom = chrStrToInt(chr, params->nChrom);
            if (chrom <= 0) 
              throw "unknown chromosome code in bgen file.";
            params->ltco_chr = chrom;
          }
        }

        break;

      case SQLITE_DONE:
        done = true;
        break;

      default:
        throw "failed reading file (" + std::string( sqlite3_errmsg(db) ) + ").";
    }
  }

  sqlite3_finalize(stmt);
  sqlite3_close(db);

  if(tmp_map.size() > params->max_condition_vars) // not relevant for gxg(=1)
    throw "number of variants used for conditional analysis is greater than maximum of " + to_string(params->max_condition_vars) + " (otherwise use --max-condition-vars)";
  else if(tmp_map.size() == 0)
    throw "no variants were found in the BGEN file";

  // replace with new map
  (*variant_names) = tmp_map;

}


void read_bgen_sample(const string& sample_file, struct param* params, std::vector<string> &ids, mstream& sout){

  int nline = 0;
  string FID, IID, line, tmp_str, fname;
  std::vector<int> sex;
  std::vector<string> IDvec;
  Files myfile;
  if( params->write_samples || params->write_masks) IDvec.resize(2);

  fname = sample_file;
  if(!file_exists (fname)) fname.append(".gz");
  sout << "   -sample file: " << fname << endl;
  myfile.openForRead(fname, sout);

  // read fid/iid information
  while (myfile.readLine(line)) {
    std::istringstream iss(line);

    if( !(iss >> FID >> IID) )
      throw "incorrectly formatted sample file at line" + to_string( ids.size() + 1 );

    // check first two lines for correct format
    if(nline == 0){

      if( (FID != "ID_1") || (IID != "ID_2") ) 
        throw "header of the sample file must start with: ID_1 ID_2";

    } else if(nline == 1){

      if( (FID != "0") || (IID != "0") ) 
        throw "second line of sample file must start with: 0 0.";

    } else {

      tmp_str = FID + "_" + IID;
      ids.push_back(tmp_str);
      if(params->write_samples || params->write_masks) {
        IDvec[0] = FID;
        IDvec[1] = IID;
        params->FIDvec.push_back(IDvec);
      }

      // get sex into IID (if no sex column, set to 0)
      if( !(iss >> FID >> IID) ) sex.push_back(0);
      else if( (IID == "0") || (IID == "NA") ) sex.push_back(0);
      else if( IID == "1" ) sex.push_back(1);
      else if( IID == "2" ) sex.push_back(2);
      else throw "unrecognized sex code in file : '" + IID + "'";

    }

    nline++;
  }

  if( params->n_samples != ids.size() )
    throw "number of samples in BGEN file does not match that in the sample file.";

  params->sex = Map<ArrayXi>(sex.data(), params->n_samples, 1);

  myfile.closeFile();
}

void read_bgen_sample(const string& sample_file, std::vector<string> &ids, mstream& sout){

  int nline = 0;
  string FID, IID, line, fname;
  Files myfile;

  fname = sample_file;
  if(!file_exists (fname)) fname.append(".gz");
  sout << "      -sample file: " << fname << endl;
  myfile.openForRead(fname, sout);

  // read fid/iid information
  while (myfile.readLine(line)) {
    std::istringstream iss(line);

    if( !(iss >> FID >> IID) )
      throw "incorrectly formatted sample file at line" + to_string( nline + 1 );

    // check first two lines for correct format
    if(nline == 0){

      if( (FID != "ID_1") || (IID != "ID_2") ) 
        throw "header of the sample file must start with: ID_1 ID_2";

    } else if(nline == 1){

      if( (FID != "0") || (IID != "0") ) 
        throw "second line of sample file must start with: 0 0.";

    } else ids.push_back(FID + "_" + IID);

    nline++;
  }

  myfile.closeFile();

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

  // build lookup table
  buildLookupTable(params->bed_lookup_table);
}


void read_bim(struct in_files* files, struct param* params, struct filter* filters, vector<snp>& snpinfo, mstream& sout) {

  bool interaction_snp_found = false;
  uint32_t nOutofOrder = 0;
  int minChr_read = 0; // enforce that chromosomes in file are sorted
  uint64 lineread = 0;
  std::vector< string > tmp_str_vec ;
  snp tmp_snp;
  string line, fname;
  Files myfile;

  fname = files->bed_prefix + ".bim";
  if(!file_exists (fname)) fname.append(".gz");
  sout << left << std::setw(20) << " * bim" << ": [" << fname << "] " << flush;
  myfile.openForRead(fname, sout);

  //if(params->set_range) cerr << params->range_chr << "\t" << params->range_min << "\t" << params->range_max<< endl;

  while (myfile.readLine(line)) {
    tmp_str_vec = string_split(line,"\t ");

    if( tmp_str_vec.size() < 6 )
      throw "incorrectly formatted bim file at line " + to_string( snpinfo.size()+1 );

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

    if (tmp_snp.chrom == -1) 
      throw "unknown chromosome code in bim file at line " + to_string( snpinfo.size()+1 );

    if( files->chr_read.empty() || (tmp_snp.chrom != files->chr_read.back() ) ) {
      files->chr_read.push_back(tmp_snp.chrom);
      if( tmp_snp.chrom <= minChr_read )
        throw "chromosomes in bim file are not in ascending order.";
      else 
        minChr_read = tmp_snp.chrom;
    }

    // check if snps are in order (same chromosome & non-decreasing positions)
    if (!snpinfo.empty() && (tmp_snp.chrom == snpinfo.back().chrom) && ( (tmp_snp.physpos < snpinfo.back().physpos) )) nOutofOrder++;

    lineread++;

    // if using GxG interaction test
    if(params->interaction_snp && (tmp_snp.ID == filters->interaction_cov)){
      if(!params->interaction_file) {
        params->interaction_snp_offset = tmp_snp.offset;
        params->ltco_chr = tmp_snp.chrom;
        interaction_snp_found = true;
      }
      continue;
    }

    // if specified chrlist/range
    if(
        (params->select_chrs && !in_chrList(tmp_snp.chrom, filters))
        ||
        (params->set_range && !in_range(tmp_snp.chrom, tmp_snp.physpos, params))
      ) continue;

    // make list of variant IDs if inclusion/exclusion file is given
    if(params->mk_snp_map){
      if (in_map(tmp_snp.ID, filters->snpID_to_ind))
        continue; // skip duplicate
      filters->snpID_to_ind[ tmp_snp.ID ] = snpinfo.size();
    }

    // keep track of how many included snps per chromosome there are
    files->chr_counts[tmp_snp.chrom-1]++;

    snpinfo.push_back(tmp_snp);
  }

  sout << "n_snps = " << lineread << endl;

  if(params->interaction_snp && !params->interaction_file && !interaction_snp_found)
    throw "SNP specified for GxG interaction test was not found.";

  if (!params->test_mode && (nOutofOrder > 0)) sout << "WARNING: Total number of snps out-of-order in bim file : " << nOutofOrder << endl;

  myfile.closeFile();

}

uint32_t read_bim(map<string, vector<uint64>>& index_map, geno_file_info* ext_file_info, struct param* params, mstream& sout) {

  int chrom;
  uint64 lineread = 0;
  std::vector< string > tmp_str_vec ;
  std::vector< uint64 > tmp_v = std::vector< uint64 >(2);
  string line, fname;
  Files myfile;

  fname = ext_file_info->file + ".bim";
  if(!file_exists (fname)) fname.append(".gz");
  myfile.openForRead(fname, sout);

  while (myfile.readLine(line)) {
    tmp_str_vec = string_split(line,"\t ");
    if( tmp_str_vec.size() < 6 )
      throw "incorrectly formatted bim file at line " + to_string( lineread+1 );
    chrom = chrStrToInt(tmp_str_vec[0], params->nChrom);
    if (chrom <= 0) 
      throw "unknown chromosome code in bgen file.";
    tmp_v[0] = lineread++, tmp_v[1] = chrom;
    index_map[ tmp_str_vec[1] ] = tmp_v;
  }

  myfile.closeFile();
  return index_map.size();
}


void read_fam(struct in_files* files, struct param* params, mstream& sout) {

  int lineread = 0;
  string line, tmp_id, fname;
  std::vector<int> sex;
  std::vector< string > tmp_str_vec, IDvec;
  Files myfile;
  if( params->write_samples || params->write_masks) IDvec.resize(2);

  fname = files->bed_prefix + ".fam";
  if(!file_exists (fname)) fname.append(".gz");
  sout << left << std::setw(20) << " * fam" << ": [" << fname << "] ";
  myfile.openForRead(fname, sout);

  while (myfile.readLine(line)) {
    tmp_str_vec = string_split(line,"\t ");

    if( tmp_str_vec.size() < 6 )
      throw "incorrectly formatted fam file at line " + to_string( lineread + 1 );

    tmp_id = tmp_str_vec[0] + "_" + tmp_str_vec[1];

    // check duplicates -- if not, store in map
    if (in_map(tmp_id, params->FID_IID_to_ind)) 
      throw "duplicate individual in fam file : FID_IID=" + tmp_id ;

    params->FID_IID_to_ind[ tmp_id ] = lineread;
    if(params->write_samples || params->write_masks) {
      IDvec[0] = tmp_str_vec[0];
      IDvec[1] = tmp_str_vec[1];
      params->FIDvec.push_back(IDvec);
    }

    // store sex
    if( tmp_str_vec[4] == "0" ) sex.push_back(0);
    else if( tmp_str_vec[4] == "1" ) sex.push_back(1);
    else if( tmp_str_vec[4] == "2" ) sex.push_back(2);
    else throw "unrecognized sex code in file : '" + tmp_str_vec[4] + "'";

    lineread++;
  }

  myfile.closeFile();
  params->n_samples = lineread;
  params->sex = Map<ArrayXi>(sex.data(), params->n_samples, 1);

  sout << "n_samples = " << params->n_samples << endl;
}

uint32_t read_fam(struct ext_geno_info& ginfo, geno_file_info* ext_file_info, Ref<ArrayXb> mask, struct param* params, mstream& sout) {

  uint32_t position;
  string line, fname;
  std::vector< string > tmp_str_vec, tmp_ids;
  Files myfile;

  fname = ext_file_info->file + ".fam";
  if(!file_exists (fname)) fname.append(".gz");
  myfile.openForRead(fname, sout);

  while (myfile.readLine(line)) {
    tmp_str_vec = string_split(line,"\t ");
    if( tmp_str_vec.size() < 6 )
      throw "incorrectly formatted fam file at line " + to_string( tmp_ids.size() + 1 );
    tmp_ids.push_back(tmp_str_vec[0] + "_" + tmp_str_vec[1]);
  }

  myfile.closeFile();

  // check if included in the analysis (if yes, store IDs)
  ginfo.sample_keep.resize(tmp_ids.size());
  ginfo.sample_index.resize(tmp_ids.size());
  for(size_t i = 0; i < tmp_ids.size(); i++) {
    ginfo.sample_keep(i) = in_map(tmp_ids[i], params->FID_IID_to_ind); 
    if(ginfo.sample_keep(i)) {
      position = params->FID_IID_to_ind[ tmp_ids[i] ];
      if(mask(position)) // analyzed sample
        ginfo.sample_index(i) = position;
      else
        ginfo.sample_keep(i) = false; 
    }
  }

  if(!ginfo.sample_keep.any())
    throw "none of the analyzed samples are present in the file";

  return tmp_ids.size();
}


void prep_bed(const uint32_t& nsamples, struct in_files* files, mstream& sout) {

  string fname;

  fname = files->bed_prefix + ".bed";
  sout << left << std::setw(20) << " * bed" << ": [" << fname << "]" << endl;
  openStream(&files->geno_ifstream, fname, std::ios::in | std::ios::binary, sout);

  uchar header[3];
  files->geno_ifstream.read( reinterpret_cast<char *> (&header[0]), 3);
  if ( (header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01) ) 
    throw "incorrect magic number in bed file.";

  // size of genotype block [(n+3)/4 = ceil(n/4.0)]
  files->bed_block_size = (nsamples+3)>>2;
  files->inbed.resize( files->bed_block_size );
}


void read_pgen_pvar_psam(struct in_files* files, struct param* params, struct filter* filters, struct geno_block* gblock, vector<snp>& snpinfo, map<int,vector<int>>& chr_map, mstream& sout) {

  gblock->nv = read_pvar(files, params, filters, snpinfo, sout);

  // check if should mask snps
  check_snps_include_exclude(files, params, filters, snpinfo, chr_map, sout);

  read_psam(files, params, sout);
  sout << "n_samples = " << params->n_samples << endl;
  gblock->ns = params->n_samples;
  // check if should mask samples
  check_samples_include_exclude(files, params, filters, sout);

  prep_pgen(files, filters, gblock, params, sout);

}

uint64 read_pvar(struct in_files* files, struct param* params, struct filter* filters, vector<snp>& snpinfo, mstream& sout) {

  bool interaction_snp_found = false;
  uint32_t nOutofOrder = 0;
  int minChr_read = 0; // enforce that chromosomes in file are sorted
  uint64 lineread = 0;
  std::vector< string > tmp_str_vec ;
  snp tmp_snp;
  string line, fname;
  Files myfile;

  fname = files->pgen_prefix + ".pvar";
  if(!file_exists (fname)) fname.append(".gz");
  sout << left << std::setw(20) << " * pvar" << ": [" << fname << "] " << flush;
  myfile.openForRead(fname, sout);

  while (myfile.readLine(line)) { // skip to main header line
    tmp_str_vec = string_split(line,"\t ");

    if( tmp_str_vec.size() < 1 )
      throw "no blank lines should be before the header line in pvar file.";

    if( tmp_str_vec[0] == "#CHROM" ) break;
  }

  // check header
  if( (tmp_str_vec.size() < 5) || 
      (tmp_str_vec[1] != "POS") || 
      (tmp_str_vec[2] != "ID") || 
      (tmp_str_vec[3] != "REF") || 
      (tmp_str_vec[4] != "ALT") )
    throw "header of pvar file does not have correct format.";

  while (myfile.readLine(line)) {
    tmp_str_vec = string_split(line,"\t ");

    if( tmp_str_vec.size() < 5 )
      throw "incorrectly formatted pvar file at line " + to_string( snpinfo.size()+1 );

    tmp_snp.chrom = chrStrToInt(tmp_str_vec[0], params->nChrom);
    tmp_snp.physpos = std::stoul( tmp_str_vec[1],nullptr,0);
    tmp_snp.ID = tmp_str_vec[2];
    tmp_snp.allele1 = tmp_str_vec[3];
    tmp_snp.allele2 = tmp_str_vec[4];
    tmp_snp.offset = lineread; // store index in file

    if (tmp_snp.chrom == -1) 
      throw "unknown chromosome code in pvar file at line " + to_string( snpinfo.size()+1 ) ;

    if( files->chr_read.empty() || (tmp_snp.chrom != files->chr_read.back() ) ) {
      files->chr_read.push_back(tmp_snp.chrom);
      if( tmp_snp.chrom <= minChr_read )
        throw "chromosomes in pvar file are not in ascending order.";
      else 
        minChr_read = tmp_snp.chrom;
    }

    // check if snps are in order (same chromosome & non-decreasing positions)
    if (!snpinfo.empty() && (tmp_snp.chrom == snpinfo.back().chrom) && ( (tmp_snp.physpos < snpinfo.back().physpos) )) nOutofOrder++;

    lineread++;

    // if using GxG interaction test
    if(params->interaction_snp && (tmp_snp.ID == filters->interaction_cov)){
        if(!params->interaction_file) {
          params->interaction_snp_offset = tmp_snp.offset;
          params->ltco_chr = tmp_snp.chrom;
          interaction_snp_found = true;
        }
      continue;
    }

    // if specified chrlist/range
    if((params->select_chrs && !in_chrList(tmp_snp.chrom, filters)) ||
        (params->set_range && !in_range(tmp_snp.chrom, tmp_snp.physpos, params))) 
      continue;

    // make list of variant IDs if inclusion/exclusion file is given
    if(params->mk_snp_map){
      if (in_map(tmp_snp.ID, filters->snpID_to_ind)) 
        continue; // skip duplicate
      filters->snpID_to_ind[ tmp_snp.ID ] = snpinfo.size();
    }

    // keep track of how many included snps per chromosome there are
    files->chr_counts[tmp_snp.chrom-1]++;

    snpinfo.push_back(tmp_snp);
  }

  sout << "n_snps = " <<  lineread << endl;

  if(params->interaction_snp && !params->interaction_file && !interaction_snp_found)
    throw "SNP specified for GxG interaction test was not found.";

  if (!params->test_mode && (nOutofOrder > 0)) sout << "WARNING: Total number of snps out-of-order in bim file : " << nOutofOrder << endl;

  myfile.closeFile();

  return lineread;
}

uint32_t read_pvar(map<string, vector<uint64>>& index_map, geno_file_info* ext_file_info, struct param* params, mstream& sout) {

  int chrom;
  uint64 lineread = 0;
  std::vector< string > tmp_str_vec ;
  std::vector< uint64 > tmp_v = std::vector< uint64 >(2);
  string line, fname;
  Files myfile;

  fname = ext_file_info->file + ".pvar";
  if(!file_exists (fname)) fname.append(".gz");
  myfile.openForRead(fname, sout);

  while (myfile.readLine(line)) { // skip to main header line
    tmp_str_vec = string_split(line,"\t ");
    if( tmp_str_vec.size() < 1 )
      throw "no blank lines should be before the header line in pvar file.";
    if( tmp_str_vec[0] == "#CHROM" ) break;
  }

  // check header
  if( (tmp_str_vec.size() < 5) || 
      (tmp_str_vec[1] != "POS") || 
      (tmp_str_vec[2] != "ID") || 
      (tmp_str_vec[3] != "REF") || 
      (tmp_str_vec[4] != "ALT") )
    throw "header of pvar file does not have correct format.";

  while (myfile.readLine(line)) {
    tmp_str_vec = string_split(line,"\t ");
    if( tmp_str_vec.size() < 5 )
      throw "incorrectly formatted pvar file at line " + to_string( lineread+1 );
    chrom = chrStrToInt(tmp_str_vec[0], params->nChrom);
    if (chrom <= 0) 
      throw "unknown chromosome code in bgen file.";
    tmp_v[0] = lineread++, tmp_v[1] = chrom;
    index_map[ tmp_str_vec[2] ] = tmp_v;
  }

  myfile.closeFile();
  return index_map.size();
}

void read_psam(struct in_files* files, struct param* params, mstream& sout) {

  int lineread = 0, sex_col = 0;
  bool col_found = false;
  string line, tmp_id, fname;
  std::vector<int> sex;
  std::vector< string > tmp_str_vec, IDvec;
  Files myfile;
  if( params->write_samples || params->write_masks) IDvec.resize(2);

  fname = files->pgen_prefix + ".psam";
  if(!file_exists (fname)) fname.append(".gz");
  sout << left << std::setw(20) << " * psam" << ": [" << fname << "] " << flush;
  myfile.openForRead(fname, sout);

  while (myfile.readLine(line)) { // skip to main header line
    tmp_str_vec = string_split(line,"\t ");

    if( tmp_str_vec.size() < 1 )
      throw "no blank lines should be before the header line in psam file.";

    if( tmp_str_vec[0] == "#FID" ) 
      break;
  }

  // check header
  if( (tmp_str_vec.size() < 2) || (tmp_str_vec[1] != "IID"))
    throw "header does not have the correct format.";

  // find if sex column is present
  auto scol = find(tmp_str_vec.begin(), tmp_str_vec.end(), "SEX");
  col_found = scol != tmp_str_vec.end();
  if(col_found) sex_col = std::distance(tmp_str_vec.begin(), scol);

  while (myfile.readLine(line)) {
    tmp_str_vec = string_split(line,"\t ");

    if( tmp_str_vec.size() < 3 )
      throw "incorrectly formatted psam file at line " + to_string( lineread + 1 ) ;

    tmp_id = tmp_str_vec[0] + "_" + tmp_str_vec[1];

    // check duplicates -- if not, store in map
    if (in_map(tmp_id, params->FID_IID_to_ind)) 
      throw "duplicate individual in fam file : FID_IID=" + tmp_id ;

    params->FID_IID_to_ind[ tmp_id ] = lineread;
    if(params->write_samples || params->write_masks) {
      IDvec[0] = tmp_str_vec[0];
      IDvec[1] = tmp_str_vec[1];
      params->FIDvec.push_back(IDvec);
    }

    // store sex
    if( col_found ){
      if(( tmp_str_vec[sex_col] == "0") || (tmp_str_vec[sex_col] == "NA") ) sex.push_back(0);
      else if( tmp_str_vec[sex_col] == "1" ) sex.push_back(1);
      else if( tmp_str_vec[sex_col] == "2" ) sex.push_back(2);
      else throw "unrecognized sex code in file : '" + tmp_str_vec[sex_col] + "'";
    } else sex.push_back(0);

    lineread++;
  }

  myfile.closeFile();
  params->n_samples = lineread;
  params->sex = Map<ArrayXi>(sex.data(), params->n_samples, 1);

}

uint32_t read_psam(struct ext_geno_info& ginfo, geno_file_info* ext_file_info, Ref<ArrayXb> mask, struct param* params, mstream& sout) {

  uint32_t position;
  string line, fname;
  std::vector< string > tmp_str_vec, tmp_ids;
  Files myfile;

  fname = ext_file_info->file + ".psam";
  if(!file_exists (fname)) fname.append(".gz");
  myfile.openForRead(fname, sout);

  while (myfile.readLine(line)) { // skip to main header line
    tmp_str_vec = string_split(line,"\t ");
    if( tmp_str_vec.size() < 1 )
      throw "no blank lines should be before the header line in psam file.";
    if( tmp_str_vec[0] == "#FID" ) 
      break;
  }

  // check header
  if( (tmp_str_vec.size() < 2) || (tmp_str_vec[1] != "IID"))
    throw "header does not have the correct format.";

  while (myfile.readLine(line)) {
    tmp_str_vec = string_split(line,"\t ");
    if( tmp_str_vec.size() < 2 )
      throw "incorrectly formatted psam file at line " + to_string( tmp_ids.size() + 1 ) ;
    tmp_ids.push_back(tmp_str_vec[0] + "_" + tmp_str_vec[1]);
  }

  myfile.closeFile();

  // check if included in the analysis (if yes, store IDs)
  ginfo.sample_keep.resize(tmp_ids.size());
  ginfo.sample_index.resize(tmp_ids.size());
  for(size_t i = 0; i < tmp_ids.size(); i++) {
    ginfo.sample_keep(i) = in_map(tmp_ids[i], params->FID_IID_to_ind); 
    if(ginfo.sample_keep(i)) {
      position = params->FID_IID_to_ind[ tmp_ids[i] ];
      if(mask(position)) // analyzed sample
        ginfo.sample_index(i) = position;
      else
        ginfo.sample_keep(i) = false; 
    }
  }

  if(ginfo.sample_keep.count() == 0)
    throw "none of the analyzed samples are present in the file";

  return tmp_ids.size();
}


void prep_pgen(struct in_files const* files, struct filter const* filters, struct geno_block* gblock, struct param* params, mstream& sout){

  int pgen_samples, pgen_variants, pgen_ac;
  vector<int> subset_indices_1based;
  string fname;

  // need to know maximum block size before loading pgen
  fname = files->pgen_prefix + ".pgen";
  sout << left << std::setw(20) << " * pgen" << ": [" << fname << "] " << endl;

  // set subset when samples have been excluded from analysis
  if( filters->ind_in_analysis.size() < gblock->ns ){
    // need to create vector of indices to keep (1-based)
    for( size_t i = 0; i < gblock->ns; i++)
      if(!filters->ind_ignore(i))
        subset_indices_1based.push_back(i+1);
  }

  gblock->pgr.Load(fname, gblock->ns, subset_indices_1based, params->threads);
  pgen_samples = gblock->pgr.GetRawSampleCt();
  pgen_variants = gblock->pgr.GetVariantCt();
  pgen_ac = gblock->pgr.GetMaxAlleleCt();

  if(pgen_samples != (int) gblock->ns)
    throw "number of samples in pgen file and psam file don't match.";
  if(pgen_variants != (int) gblock->nv)
    throw "number of variants in pgen file and pvar file don't match.";
  if(pgen_ac != 2)
    throw "only bi-allelic variants are accepted.";

  params->dosage_mode = gblock->pgr.DosagePresent();

}

void prep_pgen(uint32_t& nsamples, uint32_t& nvars, struct ext_geno_info& ginfo, geno_file_info* ext_file_info){

  vector<int> subset_indices_1based;
  string fname = ext_file_info->file + ".pgen";

  // set subset when samples have been excluded from analysis
  if( ginfo.sample_keep.count() < nsamples )
    for(size_t i = 0; i < nsamples; i++)
      if(ginfo.sample_keep(i))
        subset_indices_1based.push_back(i+1);

  ginfo.pgr.Load(fname, nsamples, subset_indices_1based, 1);
  if(ginfo.pgr.GetRawSampleCt() != nsamples)
    throw "number of samples in pgen file and psam file don't match.";
  if(ginfo.pgr.GetVariantCt() != nvars)
    throw "number of variants in pgen file and pvar file don't match.";
  if(ginfo.pgr.GetMaxAlleleCt() != 2)
    throw "only bi-allelic variants are accepted.";

  ginfo.dosage_mode = ginfo.pgr.DosagePresent();
}

// determine if snps should be included/excluded for step 1
void check_snps_include_exclude(struct in_files* files, struct param* params, struct filter* filters, vector<snp>& snpinfo, map<int,vector<int>>& chr_map, mstream& sout){

  uint32_t tmppos = 0;
  vector<snp> tmp_snpinfo;
  std::map <std::string, uint32_t> tmp_map;

  params->n_variants = snpinfo.size(); // current variants count

  if( params->condition_snps && !params->condition_file)
    get_snps_offset(filters->condition_snp_names, filters->snpID_to_ind, snpinfo, sout);

  if(params->set_range)
    sout << "   -number of variants after filtering on range = " << params->n_variants << endl;

  // if inclusion/exclusion file is given
  if(params->rm_snps || params->keep_snps) {

    assert( snpinfo.size() == filters->snpID_to_ind.size() ); // should be the same
    ArrayXb geno_mask, geno_mask_rm, geno_mask_keep;// true = keep, false = rm
    geno_mask_keep = geno_mask_rm = ArrayXb::Constant(params->n_variants, true);

    // apply masking to snps
   if( params->keep_snps ) {
      sout << "   -keeping variants specified by --extract\n";
      if(params->cormat_force_vars && (params->ld_list_file == "")) geno_mask_keep = check_in_map_from_files(filters->snpID_to_ind, files->file_snps_include, params, sout);
      else geno_mask_keep = check_in_map_from_files(filters->snpID_to_ind, files->file_snps_include, sout);
    }
    if( params->rm_snps ) {
      sout << "   -removing variants specified by --exclude\n";
      geno_mask_rm = !check_in_map_from_files(filters->snpID_to_ind, files->file_snps_exclude, sout);
    } 
    geno_mask = geno_mask_rm && geno_mask_keep;

    if(geno_mask.all()) // no snps to remove
      params->rm_snps = params->keep_snps = false;
    else {

      // delete snpID map
      filters->snpID_to_ind.clear();

      // set chr counts to 0
      std::fill(files->chr_counts.begin(), files->chr_counts.end(), 0);

      // make snpinfo only with kept elements
      params->n_variants = geno_mask.count();
      tmp_snpinfo.reserve( params->n_variants );
      for(int i = 0; i < geno_mask.size(); i++){

        if(!geno_mask(i)) continue;

        //cerr << snpinfo[i].ID << endl;
        tmp_snpinfo.push_back( snpinfo[i] );
        files->chr_counts[ snpinfo[ i ].chrom - 1 ]++;
        // remake map if needed
        if( params->keep_snp_map ) tmp_map[ snpinfo[i].ID ] = tmppos;

        tmppos++;
      }

      snpinfo.clear();
      std::vector<snp>().swap(snpinfo); // free memory
      snpinfo = tmp_snpinfo;
      if( params->keep_snp_map ) filters->snpID_to_ind = tmp_map;

    }
  }

  // check nonzero
  if(params->n_variants == 0)
    throw "no variant left to include in analysis.";
  if(params->rm_snps || params->keep_snps)
    sout << "   -number of variants remaining in the analysis = " << params->n_variants << endl;

  // go through each chromosome in order & save number of snps
  // and save how many are actually read
  vector<int> tmp_v;
  tmp_v.resize(2, 0);
  for(size_t j = 0; j < files->chr_read.size(); j++){
    int i = files->chr_read[j];
    tmp_v[0] = files->chr_counts[i-1];
    chr_map[ i ] = tmp_v;
  }

  if(params->getCorMat)
    check_ld_list(filters->snpID_to_ind, files, params, sout);

  // with OR
  check_snps_include_exclude_or(files, params, filters, snpinfo, sout);

}

// determine if snps should be included/excluded for step 2 using OR filter with MAC
void check_snps_include_exclude_or(struct in_files* files, struct param* params, struct filter* filters, vector<snp>& snpinfo, mstream& sout){

  if(!(params->rm_or || params->keep_or)) 
    return;

  assert( snpinfo.size() == filters->snpID_to_ind.size() ); // should be the same
  ArrayXb geno_mask;// if true, check MAC

  if( params->rm_or ) {
    sout << "   -removing variants specified by --exclude-or and with MAC below threshold\n";
    geno_mask = check_in_map_from_files(filters->snpID_to_ind, files->file_snps_exclude_or, sout);
  } else if( params->keep_or ) {
    sout << "   -keeping only variants specified by --extract-or or with MAC above threshold\n";
    geno_mask = !check_in_map_from_files(filters->snpID_to_ind, files->file_snps_include_or, sout);
  }

  for(int i = 0; i < geno_mask.size(); i++)
    snpinfo[ i ].MAC_fail_if_checked = geno_mask(i);

  // not needed if not using sets
  if(!params->snp_set)
    filters->snpID_to_ind.clear();

}

void check_samples_include_exclude(struct in_files const* files, struct param* params, struct filter* filters, mstream& sout){

  bool keep_ids = params->write_samples || params->write_masks;
  uint32_t ind_pos = 0, cum_pos;
  string ind_ID;
  std::map <std::string, uint32_t> new_map;
  std::map <std::string, uint32_t>::iterator itr;
  vector< string > allIDs;
  vector< vector<string> > newFIDs;

  // check number of samples
  if( params->n_samples == 0 )
    throw "no samples remaining in the analysis.";

  if( params->rm_indivs ){
    sout << "   -removing individuals specified by --remove\n";
    filters->ind_in_analysis = !check_in_map_from_files_IDs(files->file_ind_exclude, params, sout);
  } else if( params->keep_indivs ){
    sout << "   -keeping only individuals specified by --keep\n";
    filters->ind_in_analysis = check_in_map_from_files_IDs(files->file_ind_include, params, sout);
  } else
    filters->ind_in_analysis = ArrayXb::Constant(params->n_samples, true);

  // for sex-specific analyses
  if(params->sex_specific > 0){
    if(params->sex_specific == 1) // male-only
      filters->ind_in_analysis = filters->ind_in_analysis && (params->sex == 1);
    else if(params->sex_specific == 2) // female-only
      filters->ind_in_analysis = filters->ind_in_analysis && (params->sex == 2);
    sout << "   -keeping only " << ( params->sex_specific == 1 ? "male" : "female" ) << " individuals in the analysis\n";
  }

  // keep track of individual to exclude (i.e. not stored in memory)
  filters->ind_ignore = !filters->ind_in_analysis;

  if( !(filters->ind_in_analysis.all()) ) {

    if( !(filters->ind_in_analysis.any()) )
      throw "no samples remaining in the analysis.";

    // need to re-assign indices
    // retrieve all sample IDs (need to keep same order as in genotype file)
    allIDs.resize( params->n_samples );
    for (itr = params->FID_IID_to_ind.begin(); itr != params->FID_IID_to_ind.end(); ++itr) {
      ind_ID = itr->first;
      ind_pos = itr->second;
      allIDs[ ind_pos ] = ind_ID;
    }

    // create new map
    if( keep_ids ) 
      newFIDs.reserve( filters->ind_in_analysis.count() );
    cum_pos = 0;
    for( size_t j = 0; j < params->n_samples; j++){

      if( filters->ind_ignore(j) ) continue;

        new_map[ allIDs[j] ] = cum_pos;
        if(keep_ids) newFIDs.push_back( params->FIDvec[j] );
        cum_pos++;

    }

    // save map
    params->FID_IID_to_ind = new_map;
    if(keep_ids) params->FIDvec = newFIDs;

    // resize ind_in_analysis
    filters->ind_in_analysis = ArrayXb::Constant(cum_pos, true);
    sout << "   -number of genotyped individuals remaining in the analysis = " << cum_pos << endl;

  }

  params->n_samples = filters->ind_in_analysis.count();
}

ArrayXb check_in_map_from_files(map <string, uint32_t>& map_ID, vector<string> const& file_list, struct param* params, mstream& sout) {

  uint32_t lineread = 0;
  string line;
  std::vector< string > tmp_str_vec ;
  Files myfile;
  ArrayXb mask = ArrayXb::Constant( map_ID.size() , false); 

  // only allow a single extract file
  if(file_list.size() > 1) throw "cannot have multiple extract files";

  for(auto fin : file_list) {

    myfile.openForRead (fin, sout);

    while( myfile.readLine(line) ){
      tmp_str_vec = string_split(line,"\t ");

      if( tmp_str_vec.size() < 1 )
        throw "incorrectly formatted file.";
      if( in_map(tmp_str_vec[0], params->extract_vars_order) ) 
        continue; // ignore duplicates

      if( in_map(tmp_str_vec[0], map_ID) )
        mask( map_ID[ tmp_str_vec[0] ] ) = true;

      params->extract_vars_order[ tmp_str_vec[0] ] = lineread++;
    }

    myfile.closeFile();
  }

  return mask;

}

ArrayXb check_in_map_from_files(map <string, uint32_t>& map_ID, vector<string> const& file_list, mstream& sout) {

  string line;
  std::vector< string > tmp_str_vec ;
  Files myfile;
  ArrayXb mask = ArrayXb::Constant( map_ID.size() , false); 

  for(auto fin : file_list) {

    myfile.openForRead (fin, sout);

    while( myfile.readLine(line) ){
      tmp_str_vec = string_split(line,"\t ");

      if( tmp_str_vec.size() < 1 )
        throw "incorrectly formatted file.";

      if( in_map(tmp_str_vec[0], map_ID) ) 
        mask( map_ID[ tmp_str_vec[0] ] ) = true;
    }

    myfile.closeFile();
  }

  return mask;

}

ArrayXb check_in_map_from_files_IDs(vector<string> const& file_list, struct param* params, mstream& sout) {

  uint32_t nids = params->n_samples;
  findID person;
  string line;
  std::vector< string > tmp_str_vec ;
  Files myfile;
  ArrayXb mask = ArrayXb::Constant(nids, false); 

  for(auto fin : file_list) {

    myfile.openForRead (fin, sout);

    while( myfile.readLine(line) ){
      tmp_str_vec = string_split(line,"\t ");

      if( tmp_str_vec.size() < 2 )
        throw "incorrectly formatted file.";

      person = getIndivIndex(tmp_str_vec[0], tmp_str_vec[1], params, sout);
      if(!person.is_found) continue;
      mask(person.index) = true;
    }

    myfile.closeFile();
  }

  return mask;

}

void check_ld_list(map <string, uint32_t>& map_ID, struct in_files* files, struct param* params, mstream& sout) {

  if(params->ld_list_file == "") {
    if(params->extract_vars_order.size() == 0) // use all genotyped variants
      params->extract_vars_order = map_ID; 
    map<string, uint32_t >::iterator itr;
    for (itr = params->extract_vars_order.begin(); itr != params->extract_vars_order.end(); ++itr) 
      if(in_map(itr->first, map_ID))
        params->ld_sv_offsets.push_back( map_ID[ itr->first ] );
    return;
  }

  string line;
  std::vector< string > tmp_str_vec, set_keep_names;
  Files myfile;

  myfile.openForRead (params->ld_list_file, sout);

  while( myfile.readLine(line) ){
    tmp_str_vec = string_split(line,"\t ");

    if( tmp_str_vec.size() < 2 )
      throw "incorrectly formatted file (fewer than 2 entries)";

    if( in_map(tmp_str_vec[1], params->extract_vars_order) ) 
      continue; // ignore duplicates

    if(tmp_str_vec[0] == "sv"){ // single variant

      // check if in geno file & if so, store index
      if( in_map(tmp_str_vec[1], map_ID) ) 
        params->ld_sv_offsets.push_back( map_ID[ tmp_str_vec[1] ] );

    } else if(tmp_str_vec[0] == "mask"){ // mask

      if( tmp_str_vec.size() < 3 )
        throw "incorrectly formatted file (fewer than 3 entries)";
      // store gene name for extraction
      set_keep_names.push_back( tmp_str_vec[2] );

    } else throw "unrecognized entry in first column (=" + tmp_str_vec[0] + "). Should be sv/mask";

    params->extract_vars_order[ tmp_str_vec[1] ] = params->extract_vars_order.size();
  }

  params->keep_sets = params->set_select_list = set_keep_names.size() > 0;
  if(params->keep_sets)
    files->file_sets_include = {print_csv(set_keep_names)};

  myfile.closeFile();
}

// only used in step 1
void get_G(const int& block, const int& bs, const int& chrom, const uint32_t& snpcount, vector<snp> const& snpinfo, struct param const* params, struct in_files* files, struct geno_block* gblock, struct filter const* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, mstream& sout){

  auto t1 = std::chrono::high_resolution_clock::now();
  sout << " block [" << block + 1 << "] : " << flush;

  if(params->file_type == "bed")
    readChunkFromBedFileToG(bs, chrom, snpcount, snpinfo, params, files, gblock, filters, masked_indivs, phenotypes_raw, sout);
  else if(params->file_type == "pgen")
    readChunkFromPGENFileToG(bs, snpcount, snpinfo, params, gblock, filters, masked_indivs, sout);
  else if(params->streamBGEN)
    readChunkFromBGENFileToG_fast(bs, chrom, snpcount, snpinfo, params, files, gblock, filters, masked_indivs, phenotypes_raw, sout);
  else
    readChunkFromBGENFileToG(bs, chrom, snpcount, snpinfo, params, gblock, filters, masked_indivs, phenotypes_raw, sout);

  sout << bs << " snps ";

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;
}

// step 1 using BGEN library API
void readChunkFromBGENFileToG(const int& bs, const int& chrom, const uint32_t& snpcount, vector<snp> const& snpinfo, struct param const* params, struct geno_block* gblock, struct filter const* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, mstream& sout) {

  int ns;
  uint32_t index ;
  double ds, total;
  std::string chromosome, rsid;
  uint32_t position ;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;

  for(int snp = 0; snp < bs; snp++) {

    // set to correct position
    gblock->bgen.jumpto( snpinfo[ snpcount + snp ].offset );
    gblock->bgen.read_variant( &chromosome, &position, &rsid, &alleles );

    //sout << "["<< chrom << "]SNPid stored ("<< snpinfo[snpcount+snp].chrom <<") = " << snpinfo[snpcount+snp].ID<< "/ SNPIDread ("<<chromosome<<")= " << rsid << endl; exit(EXIT_FAILURE);

    assert(chrStrToInt(chromosome, params->nChrom) == chrom);
    gblock->bgen.read_probs( &probs ) ;

    ns = 0, index = 0, total = 0;
    for( std::size_t i = 0; i < probs.size(); ++i ) {

      // skip samples that were ignored from the analysis
      if( filters->ind_ignore(i) ) continue;

      ds = 0;
      for( std::size_t j = 1; j < probs[i].size(); ++j ) ds += probs[i][j] * j;

      if(ds != -3) {
        ds = params->ref_first ? ds : (2 - ds); // if ref-first, no need to switch

        if( filters->ind_in_analysis(index) ){
            total += ds;
            ns++;
        }
      }
      gblock->Gmat(snp, index) = ds;
      index++;
    }

    total /= ns;
    if( params->alpha_prior != -1 ) gblock->snp_afs(snp, 0) = total / 2;

      // impute missing
    for (size_t i = 0; i < params->n_samples; ++i ) 
      mean_impute_g(gblock->Gmat(snp, i), total, filters->ind_in_analysis(i));

  }

}

// step 1 BGEN faster file reading using OpenMP
void readChunkFromBGENFileToG_fast(const int& bs, const int& chrom, const uint32_t& start, vector<snp> const& snpinfo, struct param const* params, struct in_files* files, struct geno_block* gblock, struct filter const* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, mstream& sout) {

  vector< vector < uchar > > snp_data_blocks;
  vector< uint32_t > insize(bs), outsize(bs);
  vector<uint64> indices(bs);

  snp_data_blocks.resize( bs );
  for (int i = 0; i < bs; i++) indices[i] = snpinfo[start + i].offset;

  readChunkFromBGEN(&files->geno_ifstream, insize, outsize, snp_data_blocks, indices);

  // unpack data for each variant
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(int isnp = 0; isnp < bs; isnp++) {

    uint32_t const snpindex = start + isnp;

    uint minploidy = 0, maxploidy = 0, phasing = 0, bits_prob = 0;
    uint16_t numberOfAlleles = 0 ;
    uint32_t nindivs = 0, index;
    string tmp_buffer;
    vector<uchar>* geno_block = &snp_data_blocks[isnp];

    // set genotype data block
    vector < uchar > geno_block_uncompressed;
    geno_block_uncompressed.resize(outsize[isnp]);

    // uncompress the block
    bool compress_fail;
    if(params->zlib_compress){ // using zlib
      uLongf dest_size = outsize[isnp];
      compress_fail = (uncompress( &(geno_block_uncompressed[0]), &dest_size, &((*geno_block)[0]), insize[isnp] - 4) != Z_OK) || (dest_size != outsize[isnp]);
    } else { // using zstd
      size_t const dest_size = ZSTD_decompress(&(geno_block_uncompressed[0]), outsize[isnp], &((*geno_block)[0]), insize[isnp] - 4) ;
      compress_fail = (dest_size != outsize[isnp]);
    }
    // check it was successful
    if( compress_fail )
      throw "failed to decompress genotype data block for variant: " + snpinfo[ snpindex ].ID;

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
    assert( phasing == 0 );
    buffer++;

    // bits per probability
    std::memcpy(&bits_prob, &(buffer[0]), 1);
    assert( bits_prob == 8 );
    buffer++;

    // get dosages 
    int ns = 0;
    double prob0, prob1, prob2, total = 0;

    // parse genotype probabilities block
    index = 0;
    for(size_t i = 0; i < nindivs; i++) {

      // skip samples that were ignored from the analysis
      if( filters->ind_ignore(i) ) {
        buffer+=2;
        continue;
      }

      if(ploidy_n[i] & 0x80) {
        gblock->Gmat(isnp, index++) = -3;
        buffer+=2;
        continue;
      }

      prob0 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;
      prob1 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;
      prob2 = std::max( 1 - prob0 - prob1, 0.0);

      if(params->ref_first) 
        gblock->Gmat(isnp, index) = prob1 + 2 * prob2;
      else // switch allele0 to ALT
        gblock->Gmat(isnp, index) = prob1 + 2 * prob0;

      if( filters->ind_in_analysis(index) ){
          total += gblock->Gmat(isnp, index);
          ns++;
      }
      index++;
    }
    total /= ns;

    if (params->alpha_prior != -1) gblock->snp_afs(isnp, 0) = total / 2;

    // impute missing
    for (size_t i = 0; i < params->n_samples; ++i ) 
      mean_impute_g(gblock->Gmat(isnp, i), total, filters->ind_in_analysis(i));

  }
#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

}

// only for step 1
void readChunkFromBedFileToG(const int& bs, const int& chrom, const uint32_t& snpcount, vector<snp> const& snpinfo, struct param const* params, struct in_files* files, struct geno_block* gblock, struct filter const* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, mstream& sout) {

  int const nbl = files->bed_data_blocks.size();
  uint32_t const nmax = filters->ind_ignore.size();

  // allocate memory if needed
  if( nbl < bs ){
    files->bed_data_blocks.resize(bs);
    for (int i = nbl; i < bs; i++)
      files->bed_data_blocks[i].resize(files->bed_block_size);
  }
  // read in N/4 bytes from bed file for each snp
  for(int j = 0; j < bs; j++) {
    // set to correct position
    jumpto_bed( snpinfo[snpcount + j].offset, files->bed_block_size, files->geno_ifstream);
    files->geno_ifstream.read( reinterpret_cast<char *> (&files->bed_data_blocks[j][0]), files->bed_block_size);
  }

#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(int j = 0; j < bs; j++) {

    int hc, ns;
    uint32_t i, index ;
    double total;
    ArrayXd geno4; // genotype values for 4 samples at a time

    ns = 0, total = 0, i = 0, index = 0;

    for (size_t byte_start = 0; byte_start < files->bed_block_size; byte_start++) {

      geno4 = params->bed_lookup_table[ files->bed_data_blocks[j][byte_start] ];

      for(int bit_start = 0; bit_start < 4; bit_start++, i++){

        // skip remainder past N samples
        if(i >= nmax) break;

        // skip samples that were ignored from the analysis
        if( filters->ind_ignore(i) ) continue;

        hc = geno4(bit_start);
        if(params->ref_first && (hc != -3)) hc = 2 - hc;
        gblock->Gmat(j, index) = hc;

        if( filters->ind_in_analysis(index) && (hc != -3) ){
          total += hc;
          ns++;
        }
        index++;
      }
    }
    total /= ns;
    if(params->alpha_prior != -1) gblock->snp_afs(j, 0) = total / 2;

    // impute missing
    for (size_t i = 0; i < params->n_samples; i++) 
      mean_impute_g(gblock->Gmat(j, i), total, filters->ind_in_analysis(i));

  }

#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

}


// only for step 1
void readChunkFromPGENFileToG(const int& bs, const uint32_t& snpcount, vector<snp> const& snpinfo, struct param const* params, struct geno_block* gblock, struct filter const* filters, const Ref<const MatrixXb>& masked_indivs, mstream& sout) {

  ArrayXb oob_err = ArrayXb::Constant(bs, false);

#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(int j = 0; j < bs; j++) {

    int thread_num = 0;
#if defined(_OPENMP)
    thread_num = omp_get_thread_num();
#endif
    //cerr << "#" << thread_num << endl;

    double total;
    ArrayXb keep_indices;
    // G is MxN, but need to pass g as column vector
    ArrayXd g (params->n_samples, 1);

    // read genotype data
    if( params->dosage_mode ){
      gblock->pgr.Read(g.data(), params->n_samples, thread_num, snpinfo[snpcount+j].offset, 1);
    } else
      gblock->pgr.ReadHardcalls(g.data(), params->n_samples, thread_num, snpinfo[snpcount+j].offset, 1);

    oob_err(j) = ((g < -3) || (g > 2)).any();
    if(oob_err(j)) continue;

    gblock->Gmat.row(j) = g.matrix().transpose();

    keep_indices = filters->ind_in_analysis && (g != -3.0);
    total = keep_indices.select(g,0).sum() / keep_indices.count();

    if( params->alpha_prior != -1) gblock->snp_afs(j, 0) = total / 2;

    // impute missing
    for (size_t i = 0; i < params->n_samples; i++) 
      mean_impute_g(gblock->Gmat(j, i), total, filters->ind_in_analysis(i));

  }
#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

  if(oob_err.any()) 
    throw "there is a variant in the block that has a value not in [0,2] or missing";

}


// check if uses Layout 2 (v1.2/1.3) & check for first SNP if precision for probabilities is 8 bits
void check_bgen(const string& bgen_file, string const& file_type, bool& zlib_compress, bool& streamBGEN, uint& BGENbits, int const& nChrom){

  // for non-bgen file input, skip check
  if(file_type != "bgen") return;

  BgenParser bgen_ck;
  bgen_ck.open( bgen_file ) ;
  bool layoutV2 = bgen_ck.get_layout();
  zlib_compress = bgen_ck.get_compression();
  if( !layoutV2 ){
    streamBGEN = false;
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
  assert( chrStrToInt(tmp_buffer , nChrom) > 0 );
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

  // uncompress the block
  //cout << "zlib:"<< std::boolalpha << zlib_compress ;
  if(zlib_compress){ // using zlib
    uLongf dest_size = size_block_post_compression;
    if( (uncompress( &(geno_block_uncompressed[0]), &dest_size, &geno_block[0], size_block - 4) != Z_OK) || (dest_size != size_block_post_compression) ){
      streamBGEN = false;
      return;
    }
  } else { // using zstd
    size_t const dest_size = ZSTD_decompress(&(geno_block_uncompressed[0]), size_block_post_compression, &geno_block[0], size_block - 4) ;
    //cerr << size_block_post_compression << " " << dest_size << " " << size_block - 4 <<endl;
    if( dest_size != size_block_post_compression ){
      streamBGEN = false;
      return;
    }
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
  assert( phasing == 0 ); // must be unphased
  buffer ++;

  // bits per probability
  std::memcpy(&bits_prob, &(buffer[0]), 1);
  //cout << ",bits:"<< bits_prob ;
  BGENbits = bits_prob;;
  if( bits_prob != 8 ){
    streamBGEN = false;
    return;
  }

  streamBGEN = true;
  bfile.close();
}


// for step 2 (using BGEN library API)
void readChunkFromBGENFileToG(vector<uint64> const& indices, const int& chrom, vector<snp> const& snpinfo, struct param const* params, Ref<MatrixXd> Gmat, BgenParser& bgen, struct filter const* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, vector<variant_block> &all_snps_info, mstream& sout) {

  int const bs = indices.size();
  int hc_val, lval, ncarriers, nmales;
  uint32_t index ;
  double ds, total, mac, mval, ival, info_num, sum_pos;
  std::string chromosome, rsid;
  uint32_t position ;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;

  for(int snp = 0; snp < bs; snp++) {

    variant_block* snp_data = &(all_snps_info[snp]);
    struct snp const* snp_info = &(snpinfo[indices[snp]]);
    MapArXd Geno (Gmat.col(snp).data(), params->n_samples, 1);

    // reset variant info
    Geno = 0;
    prep_snp_stats(snp_data, params);

    hc_val = 0, index = 0, ncarriers = 0, nmales = 0;
    total = 0, mac = 0, info_num = 0;
    bool non_par = in_non_par(chrom, snp_info->physpos, params);

    // set to correct position
    bgen.jumpto( snp_info->offset );
    bgen.read_variant( &chromosome, &position, &rsid, &alleles );
    bgen.read_probs( &probs ) ;
    //sout << "["<< chrom << "]SNPid stored ("<< snp_info->chrom <<") = " << snp_info->ID<< "/ SNPIDread ("<<chromosome<<")= " << rsid << endl; exit(-1);
    assert( snp_info->ID == rsid );
    //assert(chrStrToInt(chromosome, params->nChrom) == chrom);

    for( std::size_t i = 0; i < probs.size(); ++i ) {

      // skip samples that were ignored from the analysis
      if( filters->ind_ignore(i) ) continue;

      ds = 0;
      for( std::size_t j = 1; j < probs[i].size(); ++j ) ds += probs[i][j] * j;

      if(ds != -3) {
        ds = params->ref_first ? ds : (2 - ds); // if ref-first, no need to switch

        if( filters->ind_in_analysis(index) ){
          // compute MAC using 0.5*g for males for variants on sex chr (males coded as diploid)
          // sex is 1 for males and 0 o.w.
          lval = 0, mval = ds;
          if(params->test_mode && non_par) {
            lval = (params->sex(i) == 1);
            mval = ds * 0.5 * (2 - lval);
          }
          
          if( params->ref_first )
            ival = 4 * probs[i][2] + probs[i][1] - ds * ds;
          else
            ival = 4 * probs[i][0] + probs[i][1] - ds * ds;

          // check if carrier
          if(params->build_mask && params->singleton_carriers && (ds >= 0.5)) ncarriers ++;

          total += ds;
          mac += mval;
          nmales += lval;
          info_num += ival;
          snp_data->ns1++;

          // counts by trait
          if(filters->has_missing(index)) update_trait_counts(index, ds, mval, lval, ival, snp_data, masked_indivs);

          // get genotype counts (convert to hardcall)
          if( params->htp_out ) {
            // counts for males are 0/2
            if(params->test_mode && non_par && (lval>0)) 
              hc_val = (ds < 1 ? 0 : 2);
            else
              hc_val = (int) (ds + 0.5); // round to nearest integer (0/1/2)
            update_genocounts(params->trait_mode==1, index, hc_val, snp_data->genocounts, masked_indivs, phenotypes_raw);
          } else if( params->af_cc )
            update_af_cc(index, ds, snp_data, masked_indivs, phenotypes_raw);
        }
      }

      Geno(index) = ds;
      index++;
    }

    // check MAC
    if( params->test_mode){
      compute_mac(!non_par, mac, total, nmales, ncarriers, snp_info->MAC_fail_if_checked, snp_data, params);
      if(snp_data->ignored) continue;
    }

    //sout << "SNP#" << snp + 1 << "AC=" << mac << endl;
    compute_aaf_info(total, info_num, snp_data, params);

    if(params->test_mode && params->setMinINFO && ( snp_data->info1 < params->min_INFO) ) {
      snp_data->ignored = true; continue;
    }

    // for SPA switch effect allele to minor allele
    flip_geno(total, Geno, snp_data, params);

    // apply dominant/recessive encoding & recompute mean
    if(!params->build_mask && (params->test_type > 0)){
      index = 0;
      for( std::size_t i = 0; i < probs.size(); ++i ) {
        // skip samples that were ignored from the analysis
        if( filters->ind_ignore(i) ) continue;

        if( filters->ind_in_analysis(index) && (Geno(index) != -3) ){
          if(params->test_type == 1){ //dominant
            Geno(index) = params->ref_first ? (probs[i][1] + probs[i][2]) : (probs[i][0] + probs[i][1]);
          } else if(params->test_type == 2){ //recessive
            Geno(index) = params->ref_first ? probs[i][2] : probs[i][0];
          }
        }
        index++;
      }

      sum_pos = ((Geno != -3) && filters->ind_in_analysis).select(Geno, 0).sum();
      if((params->test_type == 2) && (sum_pos < params->minHOMs)) { // filter on homALT carriers
        snp_data->ignored = true;
        continue;
      }

      total = sum_pos / snp_data->ns1;
      if(total < params->numtol) {
        snp_data->ignored = true;
        continue;
      }
    }

    // impute missing
    if(!params->build_mask)
      mean_impute_g(total, Geno, filters->ind_in_analysis);
  }

}

// for step 2 (read in raw data)
void readChunkFromBGEN(std::istream* bfile, vector<uint32_t>& insize, vector<uint32_t>& outsize, vector<vector<uchar>>& snp_data_blocks, vector<uint64>& indices){

  uint16_t SNPID_size = 0, RSID_size = 0, chromosome_size = 0 , numberOfAlleles = 0 ;
  uint32_t position = 0, allele_size = 0;
  int n_snps = indices.size();
  string tmp_buffer;

  // extract genotype data blocks single-threaded
  for(int isnp = 0; isnp < n_snps; isnp++) {
    //if(isnp % 100 == 0) cerr << "At #" << isnp+1 << endl;

    vector<uchar>* geno_block = &(snp_data_blocks[isnp]);
    uint32_t* size1 = &insize[isnp];
    uint32_t* size2 = &outsize[isnp];

    bfile->seekg( indices[isnp] );

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

  }

}

void parseSNP(const int& isnp, const int &chrom, vector<uchar>* geno_block, const uint32_t& insize, const uint32_t& outsize, struct param const* params, struct filter const* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, const snp* infosnp, struct geno_block* gblock, variant_block* snp_data, mstream& sout){

  if( ((params->file_type == "bgen") && !params->streamBGEN) || params->file_type == "pgen")
    return;

  if(params->file_type == "bgen") // uncompress and extract the dosages
    parseSnpfromBGEN(isnp, chrom, geno_block, insize, outsize, params,filters, masked_indivs, phenotypes_raw, infosnp, gblock, snp_data, sout);
  else if(params->file_type == "bed") // extract hardcalls
    parseSnpfromBed(isnp, chrom, *geno_block, params, filters, masked_indivs, phenotypes_raw, infosnp, gblock, snp_data);

}


void parseSnpfromBGEN(const int& isnp, const int &chrom, vector<uchar>* geno_block, const uint32_t& insize, const uint32_t& outsize, struct param const* params, struct filter const* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, const snp* infosnp, struct geno_block* gblock, variant_block* snp_data, mstream& sout){

  uint minploidy = 0, maxploidy = 0, phasing = 0, bits_prob = 0;
  uint16_t numberOfAlleles = 0 ;
  uint32_t nindivs = 0;
  uint32_t index;
  string tmp_buffer;

  MapArXd Geno (gblock->Gmat.col(isnp).data(), params->n_samples, 1);
  Geno = 0;
  // reset variant info
  prep_snp_stats(snp_data, params);

  // set genotype data block
  vector < uchar > geno_block_uncompressed;
  geno_block_uncompressed.resize(outsize);

  // uncompress the block
  bool compress_fail;
  if(params->zlib_compress){ // using zlib
    uLongf dest_size = outsize;
    compress_fail = (uncompress( &(geno_block_uncompressed[0]), &dest_size, &((*geno_block)[0]), insize - 4) != Z_OK) || (dest_size != outsize);
  } else { // using zstd
    size_t const dest_size = ZSTD_decompress(&(geno_block_uncompressed[0]), outsize, &((*geno_block)[0]), insize - 4) ;
    //cerr << outsize << " " << dest_size << " " << insize - 4 << endl;
    compress_fail = (dest_size != outsize);
  }
  // check it was successful
  if( compress_fail )
    throw "failed to decompress genotype data block for variant: " + infosnp->ID;

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
  bool non_par = in_non_par(chrom, infosnp->physpos, params);
  int hc_val, lval, ncarriers = 0, nmales = 0;
  double prob0, prob1, prob2, total = 0, mac = 0, mval, ival, info_num = 0, sum_pos;

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

    if(params->ref_first) 
      Geno(index) = prob1 + 2 * prob2;
    else 
      Geno(index) = prob1 + 2 * prob0; // switch allele0 to ALT

    if( filters->ind_in_analysis(index) ){
      // compute MAC using 0.5*g for males for variants on sex chr (males coded as diploid)
      // sex is 1 for males and 0 o.w.
      lval = 0, mval = Geno(index);
      if(params->test_mode && non_par) {
        lval = (params->sex(i) == 1);
        mval =  Geno(index) * 0.5 * (2 - lval);
      }

      if( params->ref_first )
        ival = 4 * prob2 + prob1 - Geno(index) * Geno(index);
      else
        ival = 4 * prob0 + prob1 - Geno(index) * Geno(index);

      // check if carrier
      if(params->build_mask && params->singleton_carriers) ncarriers += (int) (Geno(index) >= 0.5); // round dosages

      total += Geno(index);
      mac += mval;
      nmales += lval;
      info_num += ival;
      snp_data->ns1++;

      // counts by trait
      if(filters->has_missing(index)) update_trait_counts(index, Geno(index), mval, lval, ival, snp_data, masked_indivs);

      // get genotype counts (convert to hardcall)
      if( params->htp_out ) {
        // counts for males are 0/2
        if(params->test_mode && non_par && (lval>0)) 
          hc_val = (Geno(index) < 1 ? 0 : 2);
        else
          hc_val = (int) (Geno(index) + 0.5); // round to nearest integer 0/1/2
        update_genocounts(params->trait_mode==1, index, hc_val, snp_data->genocounts, masked_indivs, phenotypes_raw);
      } else if( params->af_cc )
        update_af_cc(index, Geno(index), snp_data, masked_indivs, phenotypes_raw);

    }
    index++;
  }

  // check MAC
  if( params->test_mode){
    compute_mac(!non_par, mac, total, nmales, ncarriers, infosnp->MAC_fail_if_checked, snp_data, params);
    if(snp_data->ignored) return;
  }

  compute_aaf_info(total, info_num, snp_data, params);

  // check INFO score
  if( params->setMinINFO && ( snp_data->info1 < params->min_INFO) ) {
    snp_data->ignored = true;
    return;
  }

  // for SPA switch effect allele to minor allele
  flip_geno(total, Geno, snp_data, params);

  // apply dominant/recessive encoding & recompute mean
  if(!params->build_mask && (params->test_type > 0)){
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
        index++; // bug fix
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
    sum_pos = ((Geno != -3) && filters->ind_in_analysis).select(Geno, 0).sum();
    if((params->test_type == 2) && (sum_pos < params->minHOMs)) { // filter on homALT carriers
      snp_data->ignored = true;
      return;
    }
    
    total = sum_pos / snp_data->ns1;
    if(total < params->numtol) {
      snp_data->ignored = true;
      return;
    }
  }

  // impute missing
  if(!params->build_mask)
    mean_impute_g(total, Geno, filters->ind_in_analysis);

  return;
}


void parseSnpfromBed(const int& isnp, const int &chrom, const vector<uchar>& bed_data, struct param const* params, struct filter const* filters, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, const snp* infosnp, struct geno_block* gblock, variant_block* snp_data){

  int hc, lval, ncarriers = 0, nmales;
  bool non_par = in_non_par(chrom, infosnp->physpos, params);
  uint32_t const nmax = filters->ind_ignore.size();
  uint32_t i, index ;
  double total, mac, mval, sum_pos;
  ArrayXd geno4; // genotype values for 4 samples at a time

  MapArXd Geno (gblock->Gmat.col(isnp).data(), params->n_samples, 1);
  Geno = ArrayXd::Zero(params->n_samples);
  // reset variant info
  prep_snp_stats(snp_data, params);

  total = 0, mac = 0, i = 0, index = 0, nmales = 0;
  for (size_t byte_start = 0; byte_start < bed_data.size(); byte_start++) {

    geno4 = params->bed_lookup_table[ bed_data[byte_start] ];

    for(int bit_start = 0; bit_start < 4; bit_start++, i++){

      // skip remainder past N samples
      if(i >= nmax) break;

      // skip samples that were ignored from the analysis
      if( filters->ind_ignore(i) ) continue;

      hc = geno4(bit_start);
      if(params->ref_first && (hc != -3)) hc = 2 - hc;
      Geno(index) = hc;

      if( filters->ind_in_analysis(index) && (hc != -3) ){
        // compute MAC using 0.5*g for males for variants on sex chr (males coded as diploid)
        // sex is 1 for males and 0 o.w.
        lval = 0, mval = hc;
        if(params->test_mode && non_par) {
          lval = (params->sex(i) == 1);
          mval = hc * 0.5 * (2 - lval);
          // check if not 0/2
          if( (lval == 1) && (hc == 1) ) cerr << "WARNING: genotype is 1 for a male on chrX at " << infosnp->ID << " (males should coded as diploid).";
        }

        // check if carrier
        if(params->build_mask && params->singleton_carriers) ncarriers += (int) (hc >= 1); 

        total += hc;
        mac += mval;
        nmales += lval;
        snp_data->ns1++;

        // counts by trait
        if(filters->has_missing(index)) update_trait_counts(index, Geno(index), mval, lval, 0, snp_data, masked_indivs);

        // get genotype counts
        if( params->htp_out ) 
          update_genocounts(params->trait_mode==1, index, hc, snp_data->genocounts, masked_indivs, phenotypes_raw);
        else if( params->af_cc )
          update_af_cc(index, Geno(index), snp_data, masked_indivs, phenotypes_raw);

      }
      index++;
    }
  }

  // check MAC
  if( params->test_mode){
    compute_mac(!non_par, mac, total, nmales, ncarriers, infosnp->MAC_fail_if_checked, snp_data, params);
    if(snp_data->ignored) return;
  }

  compute_aaf_info(total, 0, snp_data, params);

  // for SPA switch effect allele to minor allele
  flip_geno(total, Geno, snp_data, params);

  // apply dominant/recessive encoding & recompute mean
  if(!params->build_mask && (params->test_type > 0)){
    if(params->test_type == 1){ //dominant
      Geno = (Geno == 2).select(1, Geno);
    } else if(params->test_type == 2){ //recessive
      Geno = (Geno >= 1).select(Geno - 1, Geno);
    }

    sum_pos = ((Geno != -3) && filters->ind_in_analysis).select(Geno, 0).sum();
    if((params->test_type == 2) && (sum_pos < params->minHOMs)) { // filter on homALT carriers
      snp_data->ignored = true;
      return;
    }

    total = sum_pos / snp_data->ns1;
    if(total < params->numtol) {
      snp_data->ignored = true;
      return;
    }
  }

  // impute missing
  if(!params->build_mask)
    mean_impute_g(total, Geno, filters->ind_in_analysis);

}


// step 2
void readChunkFromPGENFileToG(vector<uint64> const& indices, const int &chrom, struct param const* params, struct filter const* filters, Ref<MatrixXd> Gmat, PgenReader& pgr, const Ref<const MatrixXb>& masked_indivs, const Ref<const MatrixXd>& phenotypes_raw, vector<snp> const& snpinfo, vector<variant_block> &all_snps_info){

  int const bs = indices.size();
  ArrayXb oob_err = ArrayXb::Constant(bs, false);

#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(int j = 0; j < bs; j++) {

    int thread_num = 0;
#if defined(_OPENMP)
    thread_num = omp_get_thread_num();
#endif

    int hc, cur_index, index, lval, nmales, ncarriers;
    double total, mac, mval, ival, eij2 = 0, sum_pos;
    ArrayXb keep_index;

    variant_block* snp_data = &(all_snps_info[j]);
    struct snp const* snp_info = &(snpinfo[ indices[j] ]);
    MapArXd Geno (Gmat.col(j).data(), params->n_samples, 1);

    // reset variant info
    prep_snp_stats(snp_data, params);

    mac = 0, index = 0, nmales = 0;
    bool non_par = in_non_par(chrom, snp_info->physpos, params);
    if( params->dosage_mode ) eij2 = 0;

    // read genotype data
    cur_index = snp_info->offset;
    if( params->dosage_mode )
      pgr.Read(Geno.data(), Geno.size(), thread_num, cur_index, 1);
    else
      pgr.ReadHardcalls(Geno.data(), Geno.size(), thread_num, cur_index, 1);

    oob_err(j) = ((Geno < -3) || (Geno > 2)).any();
    if(oob_err(j)) continue;

    keep_index = filters->ind_in_analysis && (Geno != -3.0);
    total = keep_index.select(Geno,0).sum();
    snp_data->ns1 = keep_index.count();
    //cerr << "ID: " << snp_info->ID << "\nG bounds: " << 
    //  (Geno * keep_index.cast<double>()).minCoeff() << " - " << (Geno * keep_index.cast<double>()).maxCoeff() << "\n\n";

    for (int i = 0; i < filters->ind_ignore.size(); i++) {

      // skip samples that were ignored from the analysis
      if( filters->ind_ignore(i) ) continue;

      if( keep_index(index) ){
        // compute MAC using 0.5*g for males for variants on sex chr non-PAR (males coded as diploid)
        // sex is 1 for males and 0 o.w.
        ival = 0, lval = 0, mval = Geno(index);
        if(params->test_mode && non_par) {
          lval = (params->sex(i) == 1);
          mval *= 0.5 * (2 - lval);
          // check if not 0/2
          if( !params->dosage_mode && (lval == 1) && (Geno(index) == 1) )
            cerr << "WARNING: genotype is 1 for a male on chrX at " << snp_info->ID << " (males should coded as diploid).";
        }

        if( params->dosage_mode ) ival = Geno(index) * Geno(index);

        mac += mval;
        nmales += lval;
        eij2 += ival;

        // counts by trait
        if(filters->has_missing(index)) update_trait_counts(index, Geno(index), mval, lval, ival, snp_data, masked_indivs);

        // get genotype counts
        if( params->htp_out ) {
          // counts for males are 0/2
          if(params->test_mode && non_par && (lval>0)) 
            hc = (Geno(index) < 1 ? 0 : 2);
          else
            hc = (int) (Geno(index) + 0.5); // round to nearest integer 0/1/2
          update_genocounts(params->trait_mode==1, index, hc, snp_data->genocounts, masked_indivs, phenotypes_raw);
        } else if( params->af_cc )
          update_af_cc(index, Geno(index), snp_data, masked_indivs, phenotypes_raw);

      }
      index++;
    }

    // check MAC
    if( params->test_mode){
      ncarriers = (keep_index && (Geno >= 0.5)).count(); // check carriers
      compute_mac(!non_par, mac, total, nmales, ncarriers, snp_info->MAC_fail_if_checked, snp_data, params);
      if(snp_data->ignored) continue;
    }

    compute_aaf_info(total, eij2, snp_data, params);

    // check INFO score
    if( params->dosage_mode && params->setMinINFO && ( snp_data->info1 < params->min_INFO) ) {
      snp_data->ignored = true; continue;
    }

    // for SPA switch effect allele to minor allele
    flip_geno(total, Geno, snp_data, params);

    // apply dominant/recessive encoding & recompute mean
    // pgen does not contain genotype probs for dosages so convert to hardcalls
    if(!params->build_mask && (params->test_type > 0)){
      for( size_t i = 0; i < params->n_samples; ++i ) {
        if( (Geno(i) == -3.0) || !filters->ind_in_analysis(i) ) continue;
        hc = (int) (Geno(i) + 0.5);

        if(params->test_type == 1){ //dominant
          Geno(i) = (hc == 2 ? 1 : hc);
        } else if(params->test_type == 2){ //recessive
          Geno(i) = (hc >= 1 ? hc - 1 : hc);
        }
      }

      sum_pos = ((Geno != -3) && filters->ind_in_analysis).select(Geno, 0).sum();
      if((params->test_type == 2) && (sum_pos < params->minHOMs)) { // filter on homALT carriers
        snp_data->ignored = true;
        continue;
      }

      total = sum_pos / snp_data->ns1;
      if( params->test_mode && (total < params->numtol) ) {
        snp_data->ignored = true;
        continue;
      }
    }

    // impute missing
    if(!params->build_mask)
      mean_impute_g(total, Geno, filters->ind_in_analysis);

  }
#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

  if(oob_err.any()) 
    throw "there is a variant in the block that has a value not in [0,2] or missing";

}

bool in_chrList(const int& snp_chr, struct filter const* filters){
  return in_map(snp_chr, filters->chrKeep_test);
}

string bgi_chrList(struct filter* filters, const int& nChrom){// for --chr/--chrList

  string fmt;
  vector<string> clist;
  map<int, bool >::iterator itr;

  for (itr = filters->chrKeep_test.begin(); itr != filters->chrKeep_test.end(); ++itr) {
    // add X and chrX format
    fmt = to_string(itr->first);
    clist.push_back( fmt );
    fmt = "'chr" + to_string(itr->first) + "'";
    clist.push_back( fmt );

    if(itr->first < 10){ // add 0X and chr0x format
      fmt = "'0" + to_string(itr->first) + "'";
      clist.push_back( fmt );
      fmt = "'chr0" + to_string(itr->first) + "'";
      clist.push_back( fmt );
    } else if(itr->first == nChrom){ // add XY, X, PARs
      clist.push_back( "'X'" );
      clist.push_back( "'chrX'" );
      clist.push_back( "'XY'" );
      clist.push_back( "'chrXY'" );
      clist.push_back( "'PAR1'" );
      clist.push_back( "'chrPAR1'" );
      clist.push_back( "'PAR2'" );
      clist.push_back( "'chrPAR2'" );
    }
  }

  return print_csv(clist);
}

string bgi_chrList(const int& range_chr, const int& nChrom){// for range

  string fmt = to_string(range_chr);
  vector<string> clist;

  // add X and chrX format
  clist.push_back( fmt );
  fmt = "'chr" + to_string(range_chr) + "'";
  clist.push_back( fmt );

  if(range_chr < 10){ // add 0X and chr0X format
    fmt = "'0" + to_string(range_chr) + "'";
    clist.push_back( fmt );
    fmt = "'chr0" + to_string(range_chr) + "'";
    clist.push_back( fmt );
  } else if(range_chr == nChrom){ // add XY, X, PARs
    clist.push_back( "'X'" );
    clist.push_back( "'chrX'" );
    clist.push_back( "'XY'" );
    clist.push_back( "'chrXY'" );
    clist.push_back( "'PAR1'" );
    clist.push_back( "'chrPAR1'" );
    clist.push_back( "'PAR2'" );
    clist.push_back( "'chrPAR2'" );
  }

  return print_csv(clist);
}

string bgi_rsidList(std::map <std::string, uint64>& rsids){// list of snp names

  std::map <std::string, uint64>::iterator itr;
  vector<string> clist;

  for (itr = rsids.begin(); itr != rsids.end(); ++itr)
    clist.push_back( "'" + itr->first + "'" );

  return print_csv(clist);
}

bool in_range(int const& snp_chr, uint32_t const& snp_pos, struct param const* params){

  if( (snp_chr != params->range_chr) ||
      (snp_pos < params->range_min) || 
      (snp_pos > params->range_max) )
    return false;

  return true; 
}

bool in_non_par(int const& snp_chr, uint32_t const& snp_pos, struct param const* params){

  // if not on chrX, return false
  if(snp_chr != params->nChrom) return false;

  // in par1 or par2
  if( (snp_pos <= params->par1_max_bound) || 
      (snp_pos >= params->par2_min_bound) )
    return false;

  // in non-par chrX
  return true; 
}


void skip_snps(uint64 const& offset, struct param const* params, struct in_files* files, struct geno_block* gblock){

  // set to new position based on offset
  if(params->file_type == "bed") 
    jumpto_bed(offset, files->bed_block_size, files->geno_ifstream);
  else if(params->file_type == "bgen") 
    gblock->bgen.jumpto(offset);

}

// jump to given snp index in bed file (+magic number)
void jumpto_bed(uint64 const& offset, uint64 const& bed_block_size, std::ifstream& bed_ifstream){
  bed_ifstream.seekg( 3 + offset * bed_block_size, ios_base::beg);
}

// create table for all possible values in 1 PLINK byte
void buildLookupTable(vector<ArrayXd>& lookup_table){

  uchar plink_byte;
  int bit_start;
  const int nvals = 256;
  // using 'ref-last':
  //  00 -> hom. alt
  //  10 -> missing
  //  01 -> het
  //  11 -> hom. ref
  const int maptogeno[4] = {2, -3, 1, 0};

  lookup_table.assign(nvals, ArrayXd::Zero(4));

  for(size_t i = 0; i < nvals; i++){
    plink_byte = i;

    for(int j=0; j<4; j++){
      bit_start = j<<1; // 2 bits per sample
      lookup_table[i](j) = maptogeno[ (plink_byte >> bit_start)&3 ];
    }

  }

}

void prep_snp_stats(variant_block* snp_data, struct param const* params){

  // reset variant info
  snp_data->af = ArrayXd::Zero(params->n_pheno);
  snp_data->af_case = ArrayXd::Zero(params->n_pheno);
  snp_data->af_control = ArrayXd::Zero(params->n_pheno);
  snp_data->mac = ArrayXd::Zero(params->n_pheno);
  snp_data->info = ArrayXd::Zero(params->n_pheno);
  snp_data->cf_burden = ArrayXd::Constant(params->n_pheno, -1);
  snp_data->nmales = ArrayXi::Zero(params->n_pheno);
  snp_data->ns = ArrayXi::Zero(params->n_pheno);
  snp_data->ns_case = ArrayXi::Zero(params->n_pheno);
  snp_data->ns_control= ArrayXi::Zero(params->n_pheno);
  snp_data->genocounts = MatrixXd::Zero(6, params->n_pheno);
  snp_data->ignored = false;
  snp_data->skip_int = false;
  snp_data->fitHLM = false;
  snp_data->flipped = false;
  snp_data->ns1 = 0;
  snp_data->ignored_trait = ArrayXb::Constant(params->n_pheno, false);

}

void initialize_thread_data(vector<data_thread>& all_snp_data, struct param const& params){

  for(size_t i = 0; i < all_snp_data.size(); i++){
    data_thread* snp_data = &(all_snp_data[i]);

    snp_data->chisq_val = ArrayXd::Zero(params.n_pheno);
    snp_data->pval_log = ArrayXd::Zero(params.n_pheno);
    snp_data->bhat = ArrayXd::Zero(params.n_pheno);
    snp_data->se_b = ArrayXd::Zero(params.n_pheno);
    snp_data->scores = ArrayXd::Zero(params.n_pheno);
    snp_data->cal_factor = ArrayXd::Zero(params.n_pheno);
    if(params.trait_mode){
      snp_data->stats = ArrayXd::Zero(params.n_pheno);
      snp_data->denum = ArrayXd::Zero(params.n_pheno);
    }
  }
}

void reset_thread(data_thread* snp_data, struct param const& params){

    snp_data->chisq_val = 0;
    snp_data->pval_log = 0;
    snp_data->bhat = 0;
    snp_data->se_b = 0;
    snp_data->scores = params.missing_value_double;
    snp_data->cal_factor = -1;
    if(params.trait_mode){
      snp_data->stats = 0;
      snp_data->denum = 0;
    }
    snp_data->is_sparse = false;
    snp_data->fastSPA = params.use_SPA && (!params.build_mask || (params.mask_rule_max || params.mask_rule_comphet));
}

void reset_stats(variant_block* snp_data, struct param const& params){

    snp_data->test_fail = ArrayXb::Constant(params.n_pheno, false);
    snp_data->is_corrected = ArrayXb::Constant(params.n_pheno, params.firth || params.use_SPA);
    if(params.w_interaction && params.firth) {
      snp_data->is_corrected_inter = ArrayXb::Constant(params.n_pheno, false);
      snp_data->test_fail_inter = ArrayXb::Constant(params.n_pheno, true);
    }
    if( params.joint_test ) snp_data->pval_log = ArrayXd::Zero(params.n_pheno);

    snp_data->sum_stats.resize( params.n_pheno );
    std::fill(snp_data->sum_stats.begin(), snp_data->sum_stats.end(), "");

    // multi-trait test results
    if(params.trait_set) {
      snp_data->sum_stats_mt.resize(1); // current only 1 trait set
      std::fill(snp_data->sum_stats_mt.begin(), snp_data->sum_stats_mt.end(), "");
    }

}

void update_trait_counts(int const& index, double const& genoValue, double const& macValue, int const& sexValue, double const& infoValue, variant_block* snp_data, const Ref<const MatrixXb>& mask){

  ArrayXi imask = 1 - mask.row(index).cast<int>().array(); // get masked samples

  // will subtract from total computed on all analyzed samples (masked & unmasked)
  snp_data->af -= genoValue * imask.cast<double>();
  snp_data->mac -= macValue * imask.cast<double>();
  snp_data->info -= infoValue * imask.cast<double>();
  snp_data->nmales -= imask * sexValue;
  snp_data->ns -= imask;

}

void update_genocounts(bool const& binary_mode, int const& ind, int const& hc, MatrixXd& genocounts, const Ref<const MatrixXb>& mask, const Ref<const MatrixXd>& ymat){

  if( !binary_mode ) {
    genocounts.row(hc) += mask.row(ind).cast<double>();
  } else {
    genocounts.row(hc).array() += mask.row(ind).array().cast<double>() * ymat.row(ind).array();
    genocounts.row(3 + hc).array() += mask.row(ind).array().cast<double>() * (1 - ymat.row(ind).array());
  }

}

void update_af_cc(int const& ind, double const& genoValue, variant_block* snp_data, const Ref<const MatrixXb>& mask, const Ref<const MatrixXd>& ymat){

  // only compute in cases as N-case=control
  snp_data->af_case += genoValue * mask.row(ind).array().cast<double>() * ymat.row(ind).array();
  snp_data->ns_case += mask.row(ind).array().cast<int>() * ymat.row(ind).array().cast<int>();

}

void compute_mac(bool const& auto_chrom, double& mac, double const& total, int const& nmales, int const& ncarriers, bool const& MAC_fail_if_checked, variant_block* snp_data, struct param const* params){

  if(auto_chrom) mac = total; // use MAC assuming diploid coding
  //cerr << snp_data->mac << endl << endl; 
  snp_data->mac1 = mac; // across all traits

  // for masks, identify singletons
  if(params->build_mask && !params->singleton_carriers) snp_data->singleton = ( ((int)(mac+0.5)) == 1 ); // use AAC (round for dosages)
  else if(params->build_mask && params->singleton_carriers) snp_data->singleton = (ncarriers == 1);

  // get counts by trait 
  snp_data->mac += mac; // aac
  snp_data->ns += snp_data->ns1; // ns
  snp_data->nmales += nmales; // nmales

  if(auto_chrom) {
    mac = min( mac, 2 * snp_data->ns1 - mac );
    snp_data->mac = snp_data->mac.min( 2 * snp_data->ns.cast<double>() - snp_data->mac ); // mac for each trait
  } else {
    mac = min(mac, 2 * snp_data->ns1 - nmales - mac); // males are 0/1
    snp_data->mac = snp_data->mac.min( 2 * snp_data->ns.cast<double>() - snp_data->nmales.cast<double>() - snp_data->mac );
  }

  snp_data->ignored_trait = MAC_fail_if_checked;
  snp_data->ignored_trait = snp_data->ignored_trait && (snp_data->mac < params->min_MAC);
  //cerr << snp_data->ignored_trait.cast<double>() << endl << endl; exit(EXIT_FAILURE);
  if((mac < params->min_MAC) && MAC_fail_if_checked) 
    snp_data->ignored = true;

}

void compute_aaf_info(double& total, double const& info_num, variant_block* snp_data, struct param const* params){

  // get counts by trait 
  snp_data->af += total;
  snp_data->info += info_num;

  if(params->af_cc){
    snp_data->af_control = snp_data->af - snp_data->af_case;
    snp_data->af_case /= 2 * snp_data->ns_case.cast<double>();
    snp_data->ns_control = snp_data->ns - snp_data->ns_case;
    snp_data->af_control /= 2 * snp_data->ns_control.cast<double>();
  }

  if(params->vc_test) snp_data->ac1 = total; // for skat
  total /= snp_data->ns1;
  snp_data->af1 = total / 2; // all traits
  snp_data->af /= 2 * snp_data->ns.cast<double>(); // single trait

  if(params->test_mode && params->dosage_mode){

    // all traits
    if( (snp_data->af1 == 0) || (snp_data->af1 == 1) ) snp_data->info1 = 1;
    else if(params->file_type == "bgen") snp_data->info1 = 1 - info_num / (2 * snp_data->ns1 * snp_data->af1 * (1 - snp_data->af1)); // impute
    else snp_data->info1 = (info_num / snp_data->ns1 - total * total) / (2 * snp_data->af1 * (1 - snp_data->af1)); // mach r2 info score

    // single trait
    if(params->file_type == "bgen") snp_data->info = ((snp_data->af == 0) || (snp_data->af == 1)).select(1, 1 - snp_data->info / (2 * snp_data->ns.cast<double>() * snp_data->af * (1 - snp_data->af)) );
    else snp_data->info = ((snp_data->af == 0) || (snp_data->af == 1)).select(1, (snp_data->info / snp_data->ns.cast<double>() - 4 * snp_data->af.square()) / (2 * snp_data->af * (1 - snp_data->af)) );

    if(params->setMinINFO) 
      snp_data->ignored_trait = snp_data->ignored_trait || (snp_data->info < params->min_INFO);

  }

}

void flip_geno(double& total, Ref<ArrayXd> Geno, variant_block* snp_data, struct param const* params){

  if(!params->with_flip) return;

  // switch to minor allele
  snp_data->flipped = (total > 1);

  if(snp_data->flipped){
    Geno = ( Geno != -3.0 ).select( 2 - Geno, Geno);
    total = 2 - total;
  }

}

// for rarer variants, use sparse format
void check_sparse_G(int const& isnp, int const& thread_num, struct geno_block* gblock, uint32_t const& nsamples, const Ref<const ArrayXb>& mask){

  data_thread* snp_data = &(gblock->thread_data[thread_num]);
  MapArXd Geno ( gblock->Gmat.col(isnp).data(), nsamples, 1);

  snp_data->is_sparse = (mask && (Geno != 0)).count() <= (nsamples * 0.5);
  if(snp_data->is_sparse) // get nonzero entries
    snp_data->Gsparse = mask.select(Geno,0).matrix().sparseView();

  // for SPA
  if( snp_data->fastSPA ) snp_data->fastSPA = snp_data->is_sparse;
}

// mean impute (only individuals who are not masked)
void mean_impute_g(double &geno, const double& mu, const bool& in_analysis){
  if (!in_analysis) // zero individuals masked
    geno = 0;
  else if(geno == -3) 
    geno = mu;
}
// impute all at once
void mean_impute_g(const double& mu, Ref<ArrayXd> Geno, const Ref<const ArrayXb>& in_analysis){
  Geno = (!in_analysis).select(0, Geno);
  Geno = (in_analysis && (Geno == -3)).select(mu, Geno);
}

findID getIndivIndex(const string &FID, const string &IID, struct param* params, mstream& sout){

  string tmp_str;
  findID indiv;

  // get ID of individual
  tmp_str = FID + "_" + IID;

  // check individual is in genotype data
  indiv.is_found = in_map(tmp_str, params->FID_IID_to_ind);

  if(indiv.is_found)
    indiv.index = params->FID_IID_to_ind[tmp_str];

  return indiv;
}

void residualize_geno(int const& isnp, int const& thread_num, variant_block* snp_data, bool const& force, const Ref<const MatrixXd>& X, struct geno_block* gblock, struct param const* params){

  if(snp_data->ignored) return;

  if((params->trait_mode==0) || force){
    MatrixXd beta;
    data_thread* dt_thr = &(gblock->thread_data[thread_num]);

    // project out covariates
    if(dt_thr->is_sparse) 
      beta = X.transpose() * dt_thr->Gsparse;
    else
      beta = X.transpose() * gblock->Gmat.col(isnp);

    gblock->Gmat.col(isnp) -= X * beta;

    // scale
    snp_data->scale_fac = gblock->Gmat.col(isnp).norm();
    snp_data->scale_fac /= sqrt( params->n_analyzed - X.cols() );

    if( snp_data->scale_fac < params->numtol ) {
      snp_data->ignored = true;
      return;
    }
    gblock->Gmat.col(isnp).array() /= snp_data->scale_fac;

  } else snp_data->scale_fac = 1;

}

int residualize_gmat(bool const& force, const Ref<const MatrixXd>& X, const Ref<const MatrixXd>& Graw, MatrixXd& Gres, struct param const& params){
  if((params.trait_mode==0) || force){
    MatrixXd beta = X.transpose() * Graw;
    Gres = Graw - X * beta;
  }
  return (params.n_analyzed - X.cols());
}

void check_res_geno(int const& isnp, int const& thread_num, variant_block* snp_data, bool const& force, int const& neff, const Ref<const MatrixXd>& Gres, struct geno_block* gblock, struct param const& params){

  if(snp_data->ignored) return;

  if((params.trait_mode==0) || force){

    // already computed
    gblock->Gmat.col(isnp) = Gres.col(isnp);
    if(params.skip_scaleG) return; // don't scale

    // scale
    snp_data->scale_fac = Gres.col(isnp).norm();
    snp_data->scale_fac /= sqrt( neff );

    if( snp_data->scale_fac < params.numtol ) {
      snp_data->ignored = true;
      return;
    }
    gblock->Gmat.col(isnp).array() /= snp_data->scale_fac;

  } else snp_data->scale_fac = 1;

}

void writeSnplist(string const& fname, int const& start, int const& ns, vector<snp> const& snpinfo, mstream& sout){

  ofstream ofile;
  openStream(&ofile, fname + ".snplist", ios::out, sout);

  for(int i = 0; i < ns; i++) 
    ofile << snpinfo[start+i].ID << endl;

  ofile.close();
}



// joint testing
void read_setlist(const struct in_files* files, struct param* params, struct filter* filters, vector< vector<vset> >& setinfo, vector<snp>& snpinfo, const uint64 all_masks, const double mask_max_aaf, mstream& sout) {

  bool bsize_set, all_in_geno, loo_found = false, all_w_anno, same_chr = true, no_AAF = false;
  int n_sets_incomplete = 0, n_sets_ignored = 0, n_sets_analyzed = 0;
  uint32_t lineread = 0, snp_index;
  std::vector< string > tmp_str_vec, tmp_snp_id, set_problem ;
  std::vector<int> tmpvec(3);
  string line, fname;
  Files myfile;
  ofstream report_file;

  // for snps with no anno for the set
  annoinfo ainfo_null;
  ainfo_null.regionid = 65535; // any region (set all bits to 1)
  BIT_SET(ainfo_null.id, 0);

  sout << left << std::setw(20) << " * set file" << ": [" << files->set_file << "] " << flush;
  myfile.openForRead (files->set_file, sout);
  if(params->check_mask_files) {
    line = files->out_file + "_masks_report.txt";
    openStream(&report_file, line, ios::out | ios::app, sout);
    report_file << "\n## set file: [" << files->set_file << "]\n## list of variants not in annotation or genetic data input files\n";
  }

  setinfo.resize( params->nChrom );

  // check block size
  if(params->use_max_bsize) bsize_set = false;
  else bsize_set = params->block_size >= 2;

  if(!bsize_set) params->block_size = 0;

  // if extract/exclude for sets
  if(params->keep_sets) tmpvec[2]=0;
  else if(params->rm_sets) tmpvec[2]=1;

  while (myfile.readLine(line)) {

    all_in_geno = all_w_anno = true;
    vset tmp_set;
    if(params->check_mask_files) set_problem.resize(0);

    tmp_str_vec = string_split(line,"\t ,");

    // at least 4 columns: set name | set chr | set position | variant list 
    if( tmp_str_vec.size() < 4 )
      throw "incorrectly formatted file at line " + to_string( lineread+1 ) + " (has " + to_string(tmp_str_vec.size()) + " columns)";

    // name of set
    tmp_set.ID = tmp_str_vec[0];

    // check set if using LOO 
    if(params->mask_loo || params->mask_lodo) {
      if (params->mask_loo_set != tmp_set.ID) {
        lineread++;
        continue;
      } else loo_found = true;
    }

    // chr of set
    tmp_set.chrom = chrStrToInt(tmp_str_vec[1], params->nChrom);
    if (tmp_set.chrom == -1) 
      throw "unknown chromosome code in set list file.";
    //// check if it is in chrlist
    if(params->select_chrs && !in_chrList(tmp_set.chrom, filters)) {
      lineread++;
      continue;
    }

    // position of set
    tmp_set.physpos = std::stoul( tmp_str_vec[2],nullptr,0);

    // for each variant in set, get index in genotype file
    for (size_t i = 3; i < tmp_str_vec.size(); i++){

      // check variant is in genotype file
      if (!in_map(tmp_str_vec[i], filters->snpID_to_ind)) {
        if(params->check_mask_files) set_problem.push_back(tmp_str_vec[i]);
        all_in_geno = false; continue;// mark as incomplete
      }

      // get index in geno file
      snp_index = filters->snpID_to_ind[ tmp_str_vec[i] ];
      struct snp* snp_info = &(snpinfo[ snp_index ]);

      // check chromosome
      if( tmp_set.chrom != snp_info->chrom )
        same_chr = false;

      if( params->build_mask ){
        // check annotation for set has been given for variant
        // else, assign to default annotation category 0
        if (!in_map(tmp_set.ID, snp_info->anno)) {
          all_w_anno = false;
          if(params->check_mask_files) set_problem.push_back(tmp_str_vec[i]);
          snp_info->anno[ tmp_set.ID ] = ainfo_null;
        }

        // check that variant has category in at least one of the masks
        if( (snp_info->anno[tmp_set.ID].id & all_masks) == 0 )  
          continue;
      }

      // if AAF is user defined, check it has been given for the variants
      if(params->set_aaf) {
        if(snp_info->aaf < 0) // don't add variant to set
        { no_AAF=true; continue;}
        // check that variant has AAF < max mask AAF (unless singleton)
        else if( (mask_max_aaf > 0) && (snp_info->aaf > mask_max_aaf) ) 
          continue;
      }

      // add index
      tmp_set.snp_indices.push_back(snp_index);
    }

    if(!all_in_geno || !all_w_anno ) {
      if(!all_w_anno && params->strict_check_burden) params->fail_check = true;
      if(params->check_mask_files)
        report_file << tmp_set.ID << " " << print_csv(set_problem) << endl; 
      if( tmp_set.snp_indices.size() > 0 ) n_sets_incomplete++;
      else { n_sets_ignored++; continue; } //ignore set
    }

    // sort and retain unique values
    std::sort(tmp_set.snp_indices.begin(), tmp_set.snp_indices.end());
    tmp_set.snp_indices.erase( unique( tmp_set.snp_indices.begin(), tmp_set.snp_indices.end() ), tmp_set.snp_indices.end() );

    // check how many variants are present
    if( !params->build_mask && (tmp_set.snp_indices.size() > params->max_set_size) ) 
      throw "set '" + tmp_set.ID + "' is larger than maximum allowed (=" + to_string( params->max_set_size ) + ").";

    // if not set, fix block size to maximum number of variants in set
    if(tmp_set.snp_indices.size() > params->max_bsize) params->max_bsize = tmp_set.snp_indices.size();

    // add to map if needed
    if( !(params->mask_loo || params->mask_lodo) && (params->keep_sets || params->rm_sets) ){
      tmpvec[0] = tmp_set.chrom;
      tmpvec[1] = setinfo[tmp_set.chrom - 1].size();
      filters->setID_to_ind[ tmp_set.ID ] = tmpvec;
    }

    // add to list of sets (check unique set names?)
    setinfo[tmp_set.chrom - 1].push_back(tmp_set);
    n_sets_analyzed++; lineread++;

    if(loo_found) break; // stop reading after LOO set 
  }

  myfile.closeFile();

  if(n_sets_analyzed == 0)
    throw "no sets are left to be analyzed.";

  sout << "n_sets = " << n_sets_analyzed << endl;

  if(!same_chr) sout << "WARNING: Detected at least one set where variants are not all in the same chromosome.\n";
  // report
  if(n_sets_incomplete > 0) sout << "WARNING: Detected " << n_sets_incomplete << " sets with variants not in genetic data or annotation files.\n";
  if(n_sets_ignored > 0) sout << "WARNING: Detected " << n_sets_ignored << " sets with only unknown variants (these are ignored).\n";

  if(params->check_mask_files) {
    report_file << "->Detected " << n_sets_incomplete << " sets with variants not in genetic data or annotation files.\n";
    report_file << "->Detected " << n_sets_ignored << " sets with only unknown variants.\n";
    report_file.close();
    sout << "     +report on burden input files written to [" << files->out_file + "_masks_report.txt]\n";
  }
  if(params->strict_check_burden && params->fail_check){
    string msg;
    if(params->check_mask_files) msg = " Check report for details.";
    else msg = " For more details, re-run with '--check-burden-files'.";
    throw "Annotation/Set list/Mask definition files don't agree." + msg;
  }

  if(no_AAF) sout << "WARNING: Variants in the set list file not in the AAF file will be ignored.\n";
  if( !bsize_set ) params->block_size = params->max_bsize;

  if( !(params->mask_loo || params->mask_lodo) && (params->keep_sets || params->rm_sets) ) 
    check_sets_include_exclude(bsize_set, files, params, filters, setinfo, sout);

}


// determine if sets should be included/excluded
void check_sets_include_exclude(bool const& bsize_set, const struct in_files* files, struct param* params, struct filter* filters, vector< vector<vset> >& setinfo, mstream& sout){

  uint32_t nsets = 0;
  unsigned long bsize = 0;
  vector< vector<vset> > tmp_setinfo;
  map<string, vector<int> >::iterator itr;

  //cerr << nsets << endl;

  // apply masking to sets
  if( params->rm_sets ) {
    sout << "   -removing specified sets\n";
    check_in_map_from_files_sets(false, filters->setID_to_ind, files->file_sets_exclude, params->set_select_list, sout);
  } else if( params->keep_sets ) {
    sout << "   -keeping only specified sets\n";
    check_in_map_from_files_sets(true, filters->setID_to_ind, files->file_sets_include, params->set_select_list, sout);
  }

  // re-make setinfo only with kept elements
  tmp_setinfo.resize( setinfo.size() );
  for (itr = filters->setID_to_ind.begin(); itr != filters->setID_to_ind.end(); ++itr) {

    if(itr->second[2] == 0) continue;

    tmp_setinfo[itr->second[0] - 1].push_back( setinfo[itr->second[0] - 1][itr->second[1]] );

    // track max set size
    if(!bsize_set) 
      bsize = max(bsize, setinfo[itr->second[0] - 1][itr->second[1]].snp_indices.size());

    nsets++;
  }

  // check nonzero
  if(nsets == 0)
    throw "no set left to include in analysis.";

  // free memory
  for(size_t i = 0; i < setinfo.size(); i++){
    setinfo[i].clear();
    std::vector<vset>().swap(setinfo[i]);
  }
  setinfo = tmp_setinfo;

  if(!bsize_set) params->block_size = bsize;

  // delete setID map
  filters->setID_to_ind.clear();

  sout << "     +number of sets remaining in the analysis = " << nsets << endl;

}

void check_in_map_from_files_sets(bool const& keep, map <string, vector<int>>& map_ID, vector<string> const& file_list, bool const& csv_list, mstream& sout) {

  int keep_int = (int) keep; // 0 for rm and 1 for keep
  string name;
  Files myfile;

  // user gave a comma-seperated list
  if(csv_list){
    for(auto const& setname : string_split(file_list[0],","))
      if (in_map(setname, map_ID)) 
        map_ID[ setname ][2] = keep_int;
    return;
  }

  // user gave a list of files
  for(auto fin : file_list) {

    myfile.openForRead (fin, sout);

    while( myfile.readLine(name) ){// assume single column with setname
      if (in_map(name, map_ID)) 
        map_ID[ name ][2] = keep_int;
    }

    myfile.closeFile();
  }

}

void get_masks_info(const struct in_files* files, struct param* params, struct filter* filters, map<std::string, anno_name>& anno_map, std::map <std::string, std::map <std::string, uint16_t>>& regions, vector<maskinfo>& mask_map, std::vector <std::vector<string>>& mask_out, uint64& all_masks, vector<snp>& snpinfo, mstream& sout) {

  // read annotation categories if specified
  if(params->w_anno_lab) read_anno_cat(files, params, anno_map, sout);

  // read annotations
  read_anno(params, files, filters, anno_map, regions, snpinfo, sout);

  if(params->set_aaf) read_aafs(params->tol, files, filters, snpinfo, params->aaf_file_wSingletons, sout);

  // read masks
  read_masks(files, params, anno_map, mask_map, mask_out, all_masks, sout);

}

void read_anno_cat(const struct in_files* files, struct param* params, map<string, anno_name>& anno_map, mstream& sout) {

  int lineread = 0, cval;
  uint64 null_cat = 0ULL;
  std::vector< string > tmp_str_vec ;
  string line;
  anno_name new_anno;
  Files myfile;

  sout << left << std::setw(20) << " * annotation labels" << ": [" << files->anno_labs_file << "] " << flush;
  myfile.openForRead (files->anno_labs_file, sout);


  while (myfile.readLine(line)) {

    new_anno.id = null_cat;
    tmp_str_vec = string_split(line,"\t ,");

    if( tmp_str_vec.size() != 2 )
      throw "incorrectly formatted file at line " + to_string( lineread+1 );

    // name of category
    new_anno.name = tmp_str_vec[1];
    cval = atoi( tmp_str_vec[0].c_str() );

    // check value is in 0-max 
    if( (cval < 0) || (cval >= (int)params->max_cat) )
      throw "category must be <= " + to_string( params->max_cat - 1 ) + 
        " on line " + to_string( lineread+1 ) + " (=" + tmp_str_vec[0] +  ").";

    // check category has not been specified
    if (in_map(tmp_str_vec[0], anno_map)) 
      throw "duplicate category on line " + to_string(lineread+1) + " (=" + tmp_str_vec[0] + ").";

    // set bit for category
    BIT_SET(new_anno.id, cval);

    // insert in map
    anno_map[ tmp_str_vec[0] ] = new_anno;

    lineread++;
  }

  // insert category 0 if not already given
  line = "0";
  if (!in_map(line, anno_map)) {
    new_anno.name = "NULL";
    new_anno.id = null_cat;
    BIT_SET(new_anno.id, 0);
    anno_map[ "0" ] = new_anno;
    lineread++; // count in the category
  }
  myfile.closeFile();

  sout << "n_categories = " << lineread << endl;
}

void read_anno(struct param* params, const struct in_files* files, struct filter* filters, map<string, anno_name>& anno_map, std::map <std::string, std::map <std::string, uint16_t>>& regions, vector<snp>& snpinfo, mstream& sout) {

  int lineread = 0, col_cat = 2, nregions = 0;
  uint32_t snp_pos, ncat = 0, n_anno_read = 0;
  uint64 null_id = 0ULL;
  uint16_t null_region = 0ULL;
  double set_weight = 0;
  anno_name new_anno;
  annoinfo ainfo;
  std::vector< string > tmp_str_vec ;
  string line, sname, gname;
  Files myfile;

  if(!params->w_anno_lab) { // add NULL category
    new_anno.name = "NULL";
    new_anno.id = null_id;
    BIT_SET(new_anno.id, ncat++);
    anno_map[ new_anno.name ] = new_anno;
  }

  sout << left << std::setw(20) << " * annotations " << ": [" << files->anno_file << "] " << endl;
  myfile.openForRead (files->anno_file, sout);
  if(params->vc_with_weights && (params->vc_weight_col < 4))
   throw "invalid column index specified for user-defined weights (=" + to_string( params->vc_weight_col );

  while (myfile.readLine(line)) {

    ainfo.id = null_id;
    ainfo.regionid = null_region;

    tmp_str_vec = string_split(line,"\t ,");
    if(lineread == 0) {
      // for LOVO with region
      if((params->mask_loo || params->mask_lodo) && params->w_regions && (tmp_str_vec.size() != 4))
        throw "annotation file has fewer than 4 columns for LOVO.";
      params->w_regions = !params->vc_with_weights && (tmp_str_vec.size() == 4);
      //cerr << std::boolalpha << params->w_regions << endl;
      if(params->w_regions)  col_cat = 3; // set label column
    }

    // variants | set_name | region (optional) | annotation (unique)
    if( (!params->w_regions && !params->vc_with_weights && (tmp_str_vec.size() < 3)) || 
        (params->w_regions && (tmp_str_vec.size() != 4)) || 
        (params->vc_with_weights && ((int)tmp_str_vec.size() < params->vc_weight_col)) 
        ) 
      throw "incorrectly formatted file at line " + to_string(lineread+1);

    // name of variant
    sname = tmp_str_vec[0];
    // check it is in genotype file
    if (!in_map(sname, filters->snpID_to_ind)) {
      lineread++; continue;
    }
    snp_pos = filters->snpID_to_ind[ sname ];
    struct snp* snp_info = &(snpinfo[ snp_pos ]);

    // set name
    gname = tmp_str_vec[1];
    if (!params->w_regions && in_map(gname, snp_info->anno)) 
      throw "duplicate variant annotations at line " + to_string( lineread+1 ) + ".";

    // check if matches with LOVO gene
    if((params->mask_loo || params->mask_lodo) && (gname != params->mask_loo_set)){
      lineread++; continue;
    }

    // get regions
    if(params->w_regions){

      // check if matches with LOVO region
      if(params->mask_loo && (tmp_str_vec[col_cat-1] != params->mask_loo_region)){
        lineread++; continue;
      }

      // check if new set
      if (!in_map(gname, regions)){ // create new map with region for set
        BIT_SET(ainfo.regionid, 0); // set first bit
        std::map <std::string, uint16_t> gene_region_map;
        gene_region_map[tmp_str_vec[col_cat-1]] = ainfo.regionid;
        regions[gname] = gene_region_map;
        nregions++;
      } else if (!in_map(tmp_str_vec[col_cat-1], regions[gname])) { // add region for set

        if(regions[gname].size() >= params->nmax_regions) 
          throw "cannot have more than " + to_string(params->nmax_regions) + " domains per set.";

        BIT_SET(ainfo.regionid, regions[gname].size()); // set bit for new region
        regions[gname][tmp_str_vec[col_cat-1]] = ainfo.regionid;
        nregions++;

      } else ainfo.regionid = regions[gname][tmp_str_vec[col_cat-1]]; 

    }

    // check category is in map
    if (!in_map(tmp_str_vec[col_cat], anno_map)) {

      if(params->w_anno_lab) 
        throw "unknown category at line " + to_string( lineread+1 ) +  " (=" + tmp_str_vec[col_cat] + ".";
      else { 
        // check # categories 
        if( ncat >= params->max_cat) 
          throw "cannot have more than " + to_string( params->max_cat ) + " categories (including NULL category).";

        // add to map
        new_anno.name = tmp_str_vec[col_cat];
        new_anno.id = null_id;
        BIT_SET(new_anno.id, ncat++);
        anno_map[ new_anno.name ] = new_anno;

      }
    }

    // with multiple regions for same variant & gene
    // annotation must be the same
    if (params->w_regions && in_map(gname, snp_info->anno) && 
        (snp_info->anno[gname].id != anno_map[ tmp_str_vec[col_cat] ].id) ) 
      throw "inconsistent variant annotation at line " + to_string( lineread+1 ) +  ".";

    // set bit for category
    ainfo.id |= anno_map[ tmp_str_vec[col_cat] ].id;

    //insert in snpinfo
    if (in_map(gname, snp_info->anno)) 
      snp_info->anno[gname].regionid |= ainfo.regionid;
    else
      snp_info->anno[ gname ] = ainfo;
    //if(lineread <5) cerr << snp_info->ID << "--" << ainfo.id << " " << (int) ainfo.regionid <<  endl; 

    // if using custom weights in VC tests
    if( params->vc_with_weights ){
        set_weight = convertDouble(tmp_str_vec[params->vc_weight_col - 1], params, sout);
        if( set_weight < 0 ) throw "weight = " + tmp_str_vec[params->vc_weight_col - 1] + " for variant " + sname  + " in set " + gname;
        snp_info->set_weight[ gname ] = set_weight;
    }
    
    n_anno_read++, lineread++;
  }

  myfile.closeFile();

  if(n_anno_read == 0)
    throw "annotation information could not be read. Perhaps check variant IDs matches those in the genotype file?";
  if(!params->w_anno_lab) {
    if(ncat == 0)
      throw "there are no annotation categories read from file.";
    sout << "   +number of annotations categories = " << ncat << endl;
  }
  if(params->w_regions) {
    if(nregions == 0)
      throw "there are no domains read from file.";
    sout << "   +number of domains across all sets = " << nregions << endl;
  }

}

void read_aafs(const double tol, const struct in_files* files, struct filter* filters, vector<snp>& snpinfo, bool const& wSingletons, mstream& sout) {

  int lineread = 0, id_col = 0, aaf_col = 1, singleton_col = -1, npass = 0;
  float aaf;
  uint32_t snp_pos, ncols_min = 2;
  std::vector< string > tmp_str_vec ;
  string line, sname;
  Files myfile;

  sout << left << std::setw(20) << " * user-given AAFs " << ": [" << files->aaf_file << "] " << endl;
  myfile.openForRead (files->aaf_file, sout);

  // check if there is a header line
  myfile.readLine(line);
  tmp_str_vec = string_split(line,"\t ,");
  if( tmp_str_vec.size() < ncols_min ) 
    throw "incorrectly formatted file at line " + to_string( lineread+1 );
  if( startswith(tmp_str_vec[0].c_str(), "#") ){
    // find ID column
    id_col = find_col(tmp_str_vec, "ID");
    // find AAF column
    aaf_col = find_col(tmp_str_vec, "ALT_FREQS");
    // check if columns were found
    if( (id_col < 0) || (aaf_col < 0) ) throw "could not find 'ID' or 'ALT_FREQS' in header";
    ncols_min = max(id_col, aaf_col) + 1;
  } else {
    if(wSingletons) {
      if(tmp_str_vec.size() < 3) throw "not enough columns in AAF file in line 1";
      singleton_col = 2;
      sout << left << std::setw(20) << "  -using third column to identify singleton variants\n";
      ncols_min = 3;
    }
    if (in_map(tmp_str_vec[id_col], filters->snpID_to_ind)){ // read in AAF for variant
      snp_pos = filters->snpID_to_ind[ tmp_str_vec[id_col] ];
      aaf = stof( tmp_str_vec[aaf_col] );
      snpinfo[ snp_pos ].aaf = aaf;
      if(wSingletons) snpinfo[ snp_pos ].force_singleton = check_singleton_column( tmp_str_vec[singleton_col] );
      npass++;
    }
  }
  lineread++;


  while (myfile.readLine(line)) {

    tmp_str_vec = string_split(line,"\t ,");

    if( tmp_str_vec.size() < ncols_min ) 
      throw "incorrectly formatted file at line " + to_string( lineread+1 );

    // name of variant
    sname = tmp_str_vec[id_col];

    // check it is in genotype file
    if (!in_map(sname, filters->snpID_to_ind)) {
      lineread++;
      continue;
    }
    snp_pos = filters->snpID_to_ind[ sname ];

    aaf = stof( tmp_str_vec[aaf_col] );

    /* // not necessary (other checks to remove monomorphic masks)
    if( (aaf < tol) || (aaf > (1-tol)) )
      throw "invalid AAF given at line " + to_string( lineread+1 );
    */
    snpinfo[ snp_pos ].aaf = aaf;
    if(wSingletons) snpinfo[ snp_pos ].force_singleton = check_singleton_column( tmp_str_vec[singleton_col] );

    npass++;
    lineread++;
  }

  myfile.closeFile();
  if( !npass ) throw "could not process any variant in the AAF file";

}

bool check_singleton_column(string const& col_str){

  if(col_str == "0") return false;
  else if(col_str == "1") return true;
  else throw "unindentified value in third column ('=" + col_str + "')";

}

void read_masks(const struct in_files* files, struct param* params, map<string, anno_name>& anno_map, vector<maskinfo>& minfo, std::vector <std::vector<string>>& mask_out, uint64& all_masks, mstream& sout) {

  bool valid_mask;
  int lineread = 0, ncat = 0, n_with_missing = 0, n_non_valid = 0;
  uint64 id;
  maskinfo tmp_mask;
  std::vector< string > tmp_str_vec, mask_str, anno_problem;
  mask_str.resize(2);
  string line;
  Files myfile;
  ofstream report_file;

  sout << left << std::setw(20) << " * masks " << ": [" << files->mask_file << "] " << flush;
  myfile.openForRead (files->mask_file, sout);

  if(params->check_mask_files) {
    line = files->out_file + "_masks_report.txt";
    openStream(&report_file, line, ios::out, sout);
    report_file << "## mask file: [" << files->mask_file << "]\n## list of unknown annnotations in mask file\n";
  }

  while (myfile.readLine(line)) {

    valid_mask = true;
    id = 0ULL;
    if(params->check_mask_files) anno_problem.resize(0);

    tmp_str_vec = string_split(line,"\t ,");
    ncat = tmp_str_vec.size() - 1;

    if( ncat < 1 ) 
      throw "incorrectly formatted file at line " + to_string( lineread+1 );

    // mask name
    tmp_mask.name = tmp_str_vec[0];
    mask_str[0] = tmp_mask.name;

    // check if using LOO (then single mask)
    if((params->mask_loo || params->mask_lodo) && (params->mask_loo_name != tmp_mask.name)) {
      lineread++;
      continue;
    }

    // go through each category to define mask
    std::vector< string > s_vec;
    for(int i = 0; i < ncat; i++){

      // check it is in map
      if (!in_map(tmp_str_vec[i+1], anno_map)) {
        if( tmp_str_vec[i+1].size() > 0 ){
          valid_mask = false;
          if(params->strict_check_burden) params->fail_check = true;
          if(params->check_mask_files) anno_problem.push_back(tmp_str_vec[i+1]);
        }
        continue;
      }
      s_vec.push_back( anno_map[ tmp_str_vec[i+1] ].name );

      // set bit for category
      id |= anno_map[ tmp_str_vec[i+1] ].id;
    }

    if(!valid_mask) { // one of the categories is unrecognized
      if(params->check_mask_files)
        report_file << tmp_mask.name << " " << print_csv(anno_problem) << endl;
      if(id == 0) { n_non_valid++; continue; }
      else n_with_missing++;
    }

    tmp_mask.id = id;
    mask_str[1] = print_csv( s_vec );
    //if(lineread<5)cerr << tmp_mask.name << "--" << tmp_mask.id << endl; 

    // save mask
    mask_out.push_back(mask_str);
    minfo.push_back(tmp_mask);
    params->mask_map[ tmp_mask.name ] = true;

    // take union across all categories read
    all_masks |= id;

    lineread++;
  }

  myfile.closeFile();

  sout << "n_masks = " << minfo.size() << endl;

  // report
  if(n_with_missing > 0) sout << "WARNING: Detected " << n_with_missing << " masks with unknown annotations.\n";
  if(n_non_valid > 0) sout << "WARNING: Detected " << n_non_valid << " masks with only unknown annotations (these are ignored).\n";
  if(params->check_mask_files) {
    report_file << "->Detected " << n_with_missing << " masks with unknown annotations.\n";
    report_file << "->Detected " << n_non_valid << " masks with only unknown annotations.\n";
    report_file.close();
  }

  if(minfo.size() == 0)
    throw "no masks are left to be included in the analysis.";
}


// read a single variant
void read_snp(bool const& mean_impute, uint64 const& offset, Ref<ArrayXd> Geno, Ref<ArrayXb> mask, const Eigen::Ref<const ArrayXb>& ind_ignore, struct in_files* files, PgenReader& pgr, struct param* params, bool const& check_miss){

  Geno = 0;
  if(params->file_type == "bed")
    read_snp_bed(offset, Geno, mask, ind_ignore, files, params);
  else if(params->file_type == "pgen")
    read_snp_pgen(offset, Geno, mask, pgr, params->dosage_mode);
  else
    read_snp_bgen(offset, Geno, mask, ind_ignore, files->bgen_file, params->ref_first, 0);

  if(check_miss){
    // mask missing or impute with mean
    if(mean_impute) {
      double meanG = (mask && (Geno != -3)).select(Geno,0).sum() / (mask && (Geno != -3)).count();
      Geno = (mask && (Geno == -3)).select(meanG, Geno); 
    } else  mask = (Geno != -3).select(mask, false);
  }

}

void read_snp_bed(uint64 const& offset, Ref<ArrayXd> Geno, Ref<ArrayXb> mask, const Eigen::Ref<const ArrayXb>& ind_ignore, struct in_files* files, struct param* params){

  int hc;
  uint32_t const nmax = ind_ignore.size();
  uint32_t i = 0, index = 0;
  ArrayXd geno4; // genotype values for 4 samples at a time

  // set to correct position
  jumpto_bed(offset, files->bed_block_size, files->geno_ifstream);
  files->geno_ifstream.read( reinterpret_cast<char *> (&files->inbed[0]), files->bed_block_size);

  for (size_t byte_start = 0; byte_start < files->bed_block_size; byte_start++) {

    geno4 = params->bed_lookup_table[ files->inbed[byte_start] ];

    for(int bit_start = 0; bit_start < 4; bit_start++, i++){

      // skip remainder past N samples
      if(i >= nmax) break;

      // skip samples that were ignored from the analysis
      if( ind_ignore(i) ) continue;

      if(mask(index)){
        hc = geno4(bit_start);
        if(params->ref_first && (hc != -3)) hc = 2 - hc;
        Geno(index) = hc;
      }

      index++;
    }
  }

}

void read_snp_pgen(uint64 const& offset, Ref<ArrayXd> Geno, Ref<ArrayXb> mask, PgenReader& pgr, bool const& dosage_mode){

  // read genotype data
  if( dosage_mode )
    pgr.Read(Geno.data(), Geno.size(), 0, offset, 1);
  else
    pgr.ReadHardcalls(Geno.data(), Geno.size(), 0, offset, 1);

  Geno *= mask.cast<double>();

}

// using bgen library API
// ttype is 0: add, 1:dom, 2:rec
void read_snp_bgen(uint64 const& offset, Ref<ArrayXd> Geno, Ref<ArrayXb> mask, const Eigen::Ref<const ArrayXb>& ind_ignore, string const& bgen_file, bool const& ref_first, int const& ttype){

  uint32_t index = 0;
  double ds;
  std::string chromosome, rsid;
  uint32_t position ;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;

  // open file
  BgenParser bgen_tmp;
  bgen_tmp.open( bgen_file ) ;
  bgen_tmp.jumpto( offset );
  bgen_tmp.read_variant( &chromosome, &position, &rsid, &alleles );
  bgen_tmp.read_probs( &probs ) ;

  for( std::size_t i = 0; i < probs.size(); ++i ) {
    // skip samples that were ignored from the analysis
    if( ind_ignore(i) ) continue;

    if(mask(index)){

      // get dosage from file
      ds = 0;
      for( std::size_t j = 1; j < probs[i].size(); ++j ) ds += probs[i][j] * j;

      // coding
      if(ds!= -3){
        if(ttype == 1) //dominant
          ds = ref_first ? (probs[i][1] + probs[i][2]) : (probs[i][0] + probs[i][1]);
        else if(ttype == 2) //recessive
          ds = ref_first ? probs[i][2] : probs[i][0];
        else // additive
          ds = ref_first ? ds : (2 - ds); // if ref-first, no need to switch
      }

      Geno(index) = ds;
    }

    index++;
  }

}

// check how to code G_E variant: add/dom/rec/cat
void code_snp(MatrixXd& Gcov, Ref<ArrayXb> mask, uint64 const& offset, struct filter* filters, struct in_files* files, struct param* params, mstream& sout){

  string const ttype = filters->interaction_cov_null_level;

  if( (ttype == "add" ) || (ttype == "add-homdev" ) || (ttype.size() == 0) ){ // additive
    params->add_homdev = (ttype == "add-homdev");
    if(params->gwas_condtl && params->add_homdev)
      throw "'add-homdev' coding cannot be used with --force-condtl";
    return;
  }

  MapArXd Geno ( Gcov.col(0).data(), params->n_samples, 1);

  if(ttype == "dom"){ // dominant

    if(params->file_type == "bgen")
      read_snp_bgen(offset, Geno, mask, filters->ind_ignore, files->bgen_file, params->ref_first, 1);
    else {

      if((params->file_type == "pgen") && params->dosage_mode)
        sout <<  "     +converting dosages to hardcalls for the variant\n";
      Geno = Geno.round(); // convert to hardcalls
      Geno = (mask && (Geno >= 1)).cast<double>();

    }

  } else if(ttype == "rec"){ //recessive

    if(params->file_type == "bgen")
      read_snp_bgen(offset, Geno, mask, filters->ind_ignore, files->bgen_file, params->ref_first, 2);
    else {

      if((params->file_type == "pgen") && params->dosage_mode)
        sout <<  "     +converting dosages to hardcalls for the variant\n";
      Geno = Geno.round(); // convert to hardcalls
      Geno = (mask && (Geno == 2)).cast<double>();

    }

  } else if(ttype == "cat"){ // categorical

    // convert to hardcalls
    Geno = Geno.round();

    // create 2 dummy predictors for 1 and 2
    MatrixXd newGeno = MatrixXd::Zero(params->n_samples, 2);
    newGeno.col(0).array() = (mask && (Geno == 1)).cast<double>();
    newGeno.col(1).array() = (mask && (Geno == 2)).cast<double>();

    Gcov = newGeno;
    params->interaction_cat = true;
    vector<string> vecstr{ "1","2" };
    params->interaction_lvl_names = vecstr;

  } else throw "unrecognized coding for GxG variant (can be either add/dom/rec/cat/add-homdev).";

}


void get_conditional_vars(map<string, uint64>& snps, struct in_files* files, struct param const* params, mstream& sout) {

  string line;
  std::vector< string > tmp_str_vec ;
  Files myfile;

  myfile.openForRead (files->condition_snps_list, sout);

  // get list of variants
  while( myfile.readLine(line) ){
    tmp_str_vec = string_split(line,"\t ");

    if( tmp_str_vec.size() < 1 )
      throw "incorrectly formatted file (" + files->condition_snps_list + ")";

    snps[tmp_str_vec[0]] = 0;
  }

  myfile.closeFile();

  if(snps.size() > params->max_condition_vars)
    throw "number of variants used for conditional analysis is greater than maximum of " + to_string(params->max_condition_vars) + " (otherwise use --max-condition-vars)";
  else if(snps.size() == 0)
    throw "no variants for conditional analysis given in file " + files->condition_snps_list;

  // make sure variants will be ignored by adding to exclude list
  files->file_snps_exclude.push_back(files->condition_snps_list);

}

void get_snps_offset(map<string, uint64>& snps, map<string, uint32_t>& index_map, vector<snp> const& snpinfo, mstream& sout){

  uint32_t nstart = snps.size();
  std::map <std::string, uint64> snps_found;
  std::map <std::string, uint64>::iterator itr;

  for (itr = snps.begin(); itr != snps.end(); ++itr) {
    if (in_map(itr->first, index_map))
      snps_found[ itr->first ] = snpinfo[ index_map[ itr->first ] ].offset;
  }

  snps = snps_found;

  if(snps.size() == 0)
    throw "none of the variants were found in the genotype file";
  else if(snps.size() != nstart) // enforce this
    throw to_string( nstart - snps.size() ) + " of the variants could not be found in the genotype file";

}

void get_snps_offset(map<string, uint64>& snps, map<string, vector<uint64>>& index_map, mstream& sout){

  std::map <std::string, uint64> snps_found;
  std::map <std::string, uint64>::iterator itr;

  for (itr = snps.begin(); itr != snps.end(); ++itr) 
    if (in_map(itr->first, index_map))
      snps_found[ itr->first ] = index_map[ itr->first ][0];

  snps = snps_found;

  if(snps.size() == 0)
    throw "none of the conditional variants were found in the genotype file";

}

MatrixXd extract_from_genofile(string const& setting, bool const& mean_impute, Ref<ArrayXb> mask, struct filter* filters, struct in_files* files, struct param* params, mstream& sout){

  ext_geno_info geno_info;
  map <string, vector<uint64>> tmp_map;
  map <string, uint64> tmp_map_inter;
  uint32_t nstart;
  // pointers to various information
  geno_file_info* ext_file_info;
  map <string, uint64>* variant_names;
  
  if(setting == "interaction") {
    nstart = 1;
    ext_file_info = &(files->interaction_snp_info);
    tmp_map_inter[filters->interaction_cov] = 0;
    variant_names = &(tmp_map_inter);
  } else if(setting == "conditional") {
    nstart = filters->condition_snp_names.size();
    ext_file_info = &(files->condition_snps_info);
    variant_names = &(filters->condition_snp_names);
  } else throw "unrecognized input to extract from external genotype file";


  if(ext_file_info->format == "bgen")
    setup_bgen(setting, geno_info, ext_file_info, variant_names, tmp_map, mask, files, params, filters, sout);
  else if(ext_file_info->format == "pgen")
    setup_pgen(geno_info, ext_file_info, tmp_map, mask, params, sout);
  else if(ext_file_info->format == "bed")
    setup_bed(geno_info, ext_file_info, tmp_map, mask, params, sout);

  if(params->debug) cerr << geno_info.sample_keep.count() << " " << tmp_map.size() << "\n\n" << geno_info.sample_index.head(10) << "\n\n";

  
  if(!ext_file_info->with_bgi) { // filter using map
    get_snps_offset((*variant_names), tmp_map, sout);
    if(setting == "interaction") params->ltco_chr = tmp_map[filters->interaction_cov][1];
  }

  // check number of variants
  if(variant_names->size() > params->max_condition_vars)
    throw "number of variants used for conditional analysis is greater than maximum of " + to_string(params->max_condition_vars) + " (otherwise use --max-condition-vars)";
  else if(variant_names->size() != nstart)
    throw to_string( nstart - variant_names->size() ) + " of the variants could not be found in the genotype file";

  if(setting == "conditional") 
    sout <<  "      -n_used = " << variant_names->size() << endl;
  MatrixXd Gmat = MatrixXd::Constant(params->n_samples, variant_names->size(), -3); // set all to missing

  // read in variants & impute if missing
  // note variant with all missing will be captured in intercept
  if((ext_file_info->format == "bgen") && geno_info.streamBGEN)
    read_snps_bgen(mean_impute, (*variant_names), Gmat, geno_info, mask, ext_file_info->file, params);
  else if(ext_file_info->format == "bgen")
    read_snps_bgen(mean_impute, (*variant_names), Gmat, geno_info, mask, ext_file_info->file);
  else if(ext_file_info->format == "pgen")
    read_snps_pgen(mean_impute, (*variant_names), Gmat, geno_info, mask);
  else if(ext_file_info->format == "bed")
    read_snps_bed(mean_impute, (*variant_names), Gmat, geno_info, mask, ext_file_info->file, params, sout);

  // if did not impute missing with mean, mask samples
  if(!mean_impute){
    mask = ((Gmat.array() == -3).rowwise().any()).select(false, mask);
    Gmat.array().colwise() *= mask.cast<double>();
  }

  // if ref-first for GxG (inactive for pgen)
  // bug fix: bgen api already assumes ref-first
  if((setting == "interaction") && files->interaction_snp_info.ref_first && !((ext_file_info->format == "bgen") && !geno_info.streamBGEN))
    Gmat.array() = (2 - Gmat.array()).colwise() * mask.cast<double>();

  return Gmat;
}

// for conditional analyses
void setup_bgen(string const& setting, struct ext_geno_info& ginfo, geno_file_info* ext_file_info, map <string, uint64>* variant_names, map<string, vector<uint64>>& index_map, Ref<ArrayXb> mask, struct in_files* files, struct param* params, struct filter* filters, mstream& sout){

  sout << "      -extracting variants from file [" << ext_file_info->file << "]\n";

  int chrom;
  uint32_t lineread = 0;
  uint BGENbits;
  uint64 offset;
  std::vector< string > tmp_ids ;
  std::vector< uint64 > tmp_v = std::vector< uint64 >(2);
  BgenParser bgen_tmp;

  uint32_t position ;
  std::string chromosome, rsid, msg;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;

  // check if can use faster file stream
  check_bgen(ext_file_info->file, ext_file_info->format, ginfo.zlib_compress, ginfo.streamBGEN, BGENbits, params->nChrom);

  // open file and print file info
  bgen_tmp.open( ext_file_info->file ) ;

  // get info for variants
  if( ext_file_info->with_bgi ) read_bgi_file(setting, bgen_tmp, ext_file_info, variant_names, params, sout);
  else {
    offset = bgen_tmp.get_position();
    while(bgen_tmp.read_variant( &chromosome, &position, &rsid, &alleles )) {

      assert(alleles.size() == 2) ; // only bi-allelic allowed
      // check phasing for first variant
      if(lineread++ == 0){
        bgen_tmp.read_probs( &probs ) ;
        if( probs[0].size() != 3 ) // unphased only 
          throw "only unphased bgen are supported.";
      } else bgen_tmp.ignore_probs();

      chrom = chrStrToInt(chromosome, params->nChrom);
      if (chrom <= 0) 
        throw "unknown chromosome code in bgen file.";

      // make list of variant IDs
      tmp_v[0] = offset, tmp_v[1] = chrom;
      index_map[ rsid ] = tmp_v;

      offset = bgen_tmp.get_position();
    }

  }

  // get sample IDs (from sample file or directly from bgen file)
  if( ext_file_info->with_sample ) {
    read_bgen_sample(ext_file_info->sample, tmp_ids, sout);
  } else {
    bgen_tmp.get_sample_ids(
        [&tmp_ids]( std::string const& id ) { tmp_ids.push_back( id ) ; } );
  }

  // check if included in the analysis (if yes, store IDs)
  ginfo.sample_keep.resize(tmp_ids.size());
  ginfo.sample_index.resize(tmp_ids.size());
  for(size_t i = 0; i < tmp_ids.size(); i++) {
    ginfo.sample_keep(i) = in_map(tmp_ids[i], params->FID_IID_to_ind); 
    if(ginfo.sample_keep(i)) {
      position = params->FID_IID_to_ind[ tmp_ids[i] ];
      if(mask(position)) // analyzed sample
        ginfo.sample_index(i) = position;
      else
        ginfo.sample_keep(i) = false; 
    }
  }

  if(ginfo.sample_keep.count() == 0)
    throw "none of the analyzed samples are present in the file";

}

// fast streaming
void read_snps_bgen(bool const& mean_impute, map<string, uint64>& snp_map, Ref<MatrixXd> Gmat, struct ext_geno_info& ginfo, Ref<ArrayXb> mask, string const& bgen_file, struct param* params){

  int bs = snp_map.size();
  vector< vector < uchar > > snp_data_blocks;
  vector< uint32_t > insize, outsize;
  vector<uint64> indices;
  ArrayXb read_error = ArrayXb::Constant(bs, false);
  std::map <std::string, uint64>::iterator itr;
  std::ifstream bgen_ifstream;

  snp_data_blocks.resize( bs );
  insize.resize( bs );
  outsize.resize( bs );
  indices.reserve( bs );
  for (itr = snp_map.begin(); itr != snp_map.end(); ++itr)
    indices.push_back(itr->second);
  std::sort(indices.begin(), indices.end());// sort indices to read in order

  bgen_ifstream.open( bgen_file, ios::in | ios::binary);
  readChunkFromBGEN(&bgen_ifstream, insize, outsize, snp_data_blocks, indices);


  // unpack data for each variant
#if defined(_OPENMP)
  setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
#endif
  for(int isnp = 0; isnp < bs; isnp++) {

    uint minploidy = 0, maxploidy = 0, phasing = 0, bits_prob = 0;
    uint16_t numberOfAlleles = 0 ;
    uint32_t nindivs = 0, index;
    string tmp_buffer;
    vector<uchar>* geno_block = &snp_data_blocks[isnp];

    // set genotype data block
    vector < uchar > geno_block_uncompressed;
    geno_block_uncompressed.resize(outsize[isnp]);

    // uncompress the block
    bool compress_fail;
    if(ginfo.zlib_compress){ // using zlib
      uLongf dest_size = outsize[isnp];
      compress_fail = (uncompress( &(geno_block_uncompressed[0]), &dest_size, &((*geno_block)[0]), insize[isnp] - 4) != Z_OK) || (dest_size != outsize[isnp]);
    } else { // using zstd
      size_t const dest_size = ZSTD_decompress(&(geno_block_uncompressed[0]), outsize[isnp], &((*geno_block)[0]), insize[isnp] - 4) ;
      compress_fail = (dest_size != outsize[isnp]);
    }
    // check it was successful
    if( compress_fail ){
      read_error(isnp) = true;
      continue; // don't use throw as not thread-safe
    }

    // stream to uncompressed block
    uchar *buffer = &geno_block_uncompressed[0];
    // sample size in file
    std::memcpy(&nindivs, &(buffer[0]), 4);
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
    assert( phasing == 0 );
    buffer++;
    // bits per probability
    std::memcpy(&bits_prob, &(buffer[0]), 1);
    assert( bits_prob == 8 );
    buffer++;

    // get dosages 
    int ns = 0;
    double prob0, prob1, total = 0;
    MapArXd Geno (Gmat.col(isnp).data(), Gmat.rows(), 1);

    // parse genotype probabilities block
    for(size_t i = 0; i < nindivs; i++) {

      // skip samples that were ignored from the analysis
      if( !ginfo.sample_keep(i) ) {
        buffer+=2;
        continue;
      }
      index = ginfo.sample_index(i);

      if(ploidy_n[i] & 0x80) {
        Geno(index) = -3;
        buffer+=2;
        continue;
      }

      prob0 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;
      prob1 = double((*reinterpret_cast< uint8_t const* >( buffer++ ))) / 255.0;
      Geno(index) = prob1 + 2 * prob0;

      if( Geno(index) != -3 ){
          total += Geno(index);
          ns++;
      }
    }

    if(ns==0) {Geno=-3; continue;} // mask all samples
    if(mean_impute) mean_impute_g(total/ns, Geno, mask);

  }
#if defined(_OPENMP)
  setNbThreads(params->threads);
#endif

  bgen_ifstream.close();
  if(read_error.any())
    throw "failed to decompress genotype data block.";

}

void read_snps_bgen(bool const& mean_impute, map<string, uint64>& snp_map, Ref<MatrixXd> Gmat, struct ext_geno_info& ginfo, Ref<ArrayXb> mask, string const& bgen_file){

  int index = 0, count = 0;
  double ds, total, ns;
  std::map <std::string, uint64>::iterator itr;

  std::string chromosome, rsid;
  uint32_t position ;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;

  // open file
  BgenParser bgen_tmp;
  bgen_tmp.open( bgen_file ) ;

  for (itr = snp_map.begin(); itr != snp_map.end(); ++itr, count++) {

    MapArXd Geno (Gmat.col(count).data(), Gmat.rows(), 1);

    bgen_tmp.jumpto( itr->second );
    bgen_tmp.read_variant( &chromosome, &position, &rsid, &alleles );
    bgen_tmp.read_probs( &probs ) ;

    total = 0, ns = 0;
    for( std::size_t i = 0; i < probs.size(); ++i ) {
      // skip samples that were ignored from the analysis
      if( !ginfo.sample_keep(i) ) continue;
      index = ginfo.sample_index(i);

      // get dosage from file
      ds = 0;
      for( std::size_t j = 1; j < probs[i].size(); ++j ) ds += probs[i][j] * j;
      // does not matter if ref-first/ref-last since used as covar

      if(ds != -3) {total += ds; ns++;}
      Geno(index) = ds;
    }

    if(ns==0) {Geno=-3; continue;} // mask all samples

    // impute missing
    if(mean_impute) mean_impute_g(total/ns, Geno, mask);

  }

}

void setup_pgen(struct ext_geno_info& ginfo, geno_file_info* ext_file_info, map<string, vector<uint64>>& index_map, Ref<ArrayXb> mask, struct param* params, mstream& sout) {

 sout << "      -extracting variants using PGEN file prefix [" << ext_file_info->file << "]\n";

  uint32_t nv = read_pvar(index_map, ext_file_info, params, sout);
  uint32_t ns = read_psam(ginfo, ext_file_info, mask, params, sout);
  //cerr << "Nsamples=" << ns << "\tNvariants=" << nv << endl;
  prep_pgen(ns, nv, ginfo, ext_file_info);

}

void read_snps_pgen(bool const& mean_impute, map<string, uint64>& snp_map, Ref<MatrixXd> Gmat, struct ext_geno_info& ginfo, Ref<ArrayXb> mask){

  int index = 0, count = 0;
  double total, ns;
  std::map <std::string, uint64>::iterator itr;
  ArrayXd Gread (ginfo.sample_keep.count()); // some analyzed samples may not be in file

  for (itr = snp_map.begin(); itr != snp_map.end(); ++itr, count++) {

    MapArXd Geno (Gmat.col(count).data(), Gmat.rows(), 1);

    // read genotype data
    if( ginfo.dosage_mode )
      ginfo.pgr.Read(Gread.data(), Gread.size(), 0, itr->second, 1);
    else
      ginfo.pgr.ReadHardcalls(Gread.data(), Gread.size(), 0, itr->second, 1);

    total = 0, ns = 0;
    for(int i_raw = 0, i = 0; i_raw < ginfo.sample_keep.size(); ++i_raw ) {
      // skip samples that were ignored from the analysis
      if( !ginfo.sample_keep(i_raw) ) continue;
      index = ginfo.sample_index(i_raw);

      if(Gread(i) != -3) {total += Gread(i); ns++;}
      Geno(index) = Gread(i++);
    }

    if(ns==0) {Geno=-3; continue;} // mask all samples

    // impute missing
    if(mean_impute) mean_impute_g(total/ns, Geno, mask);

  }

}


void setup_bed(struct ext_geno_info& ginfo, geno_file_info* ext_file_info, map<string, vector<uint64>>& index_map, Ref<ArrayXb> mask, struct param* params, mstream& sout) {

 sout << "      -extracting variants using BED file prefix [" << ext_file_info->file << "]\n";

  uint32_t nv = read_bim(index_map, ext_file_info, params, sout);
  uint32_t ns = read_fam(ginfo, ext_file_info, mask, params, sout);
  if(params->debug) cerr << "Nsamples=" << ns << "\tNvariants=" << nv << endl;

  // check if need to make lookup table
  if(params->bed_lookup_table.size() == 0)
    buildLookupTable(params->bed_lookup_table);
}


void read_snps_bed(bool const& mean_impute, map<string, uint64>& snp_map, Ref<MatrixXd> Gmat, struct ext_geno_info& ginfo, Ref<ArrayXb> mask, string const& bed_prefix, struct param* params, mstream& sout){

  int index = 0, count = 0, hc;
  uint32_t const nmax = ginfo.sample_keep.size();
  uint64 bed_block_size;
  double total, ns;
  std::map <std::string, uint64>::iterator itr;
  ArrayXd geno4; // genotype values for 4 samples at a time
  ifstream bed_ifstream;
  vector<uchar> inbed;

  // open file and check header
  uchar header[3];
  string fname = bed_prefix + ".bed";
  openStream(&bed_ifstream, fname, std::ios::in | std::ios::binary, sout);
  bed_ifstream.read( reinterpret_cast<char *> (&header[0]), 3);
  if ( (header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01) ) 
    throw "incorrect magic number in bed file.";
  // size of genotype block [(n+3)/4 = ceil(n/4.0)]
  bed_block_size = (nmax+3)>>2;
  inbed.resize( bed_block_size );


  for (itr = snp_map.begin(); itr != snp_map.end(); ++itr, count++) {

    MapArXd Geno (Gmat.col(count).data(), Gmat.rows(), 1);

    // set to correct position
    jumpto_bed(itr->second, bed_block_size, bed_ifstream);
    bed_ifstream.read( reinterpret_cast<char *> (&inbed[0]), bed_block_size);
    total = 0, ns = 0;

    for (size_t byte_start = 0, i = 0; byte_start < bed_block_size; byte_start++) {

      geno4 = params->bed_lookup_table[ inbed[byte_start] ];

      for(int bit_start = 0; bit_start < 4; bit_start++, i++){

        // skip remainder past N samples
        if(i >= nmax) break;
        // skip samples that were ignored from the analysis
        if( !ginfo.sample_keep(i) ) continue;
        index = ginfo.sample_index(i);

        hc = geno4(bit_start);
        if(hc != -3) {total += hc; ns++;}
        Geno(index) = hc;
      }
    }

    if(ns==0) {Geno=-3; continue;} // mask all samples

    // impute missing
    if(mean_impute) mean_impute_g(total/ns, Geno, mask);

  }

}
