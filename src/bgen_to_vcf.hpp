
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include "genfile/bgen/bgen.hpp"

// ProbSetter is a callback object appropriate for passing to bgen::read_genotype_data_block() or
// the synonymous method of genfile::bgen::View. See the comment in bgen.hpp above
// bgen::read_genotype_data_block(), or the bgen wiki for a description of the API.
// The purpose of this object is to store genotype probability values in the desired
// data structure (which here is a vector of vectors of doubles).
struct ProbSetter {
  typedef std::vector< std::vector< double > > Data ;
  ProbSetter( Data* result ):
    m_result( result ),
    m_sample_i(0)
  {}

  // Called once allowing us to set storage.
  void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
    m_result->clear() ;
    m_result->resize( number_of_samples ) ;
  }

  // If present with this signature, called once after initialise()
  // to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
  // This enables us to set up storage for the data ahead of time.
  void set_min_max_ploidy( uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries ) {
    for( std::size_t i = 0; i < m_result->size(); ++i ) {
      m_result->at( i ).reserve( max_entries ) ;
    }
  }

  // Called once per sample to determine whether we want data for this sample
  bool set_sample( std::size_t i ) {
    m_sample_i = i ;
    // Yes, here we want info for all samples.
    return true ;
  }

  // Called once per sample to set the number of probabilities that are present.
  void set_number_of_entries(
      std::size_t ploidy,
      std::size_t number_of_entries,
      genfile::OrderType order_type,
      genfile::ValueType value_type
      ) {
    assert( value_type == genfile::eProbability ) ;
    m_result->at( m_sample_i ).resize( number_of_entries ) ;
    m_entry_i = 0 ;
  }

  // Called once for each genotype (or haplotype) probability per sample.
  void set_value( uint32_t, double value ) {
    m_result->at( m_sample_i ).at( m_entry_i++ ) = value ;
  }

  // Ditto, but called if data is missing for this sample.
  void set_value( uint32_t, genfile::MissingValue value ) {
    // Here we encode missing probabilities with -1
    m_result->at( m_sample_i ).at( m_entry_i++ ) = -1 ;
  }

  // If present with this signature, called once after all data has been set.
  void finalise() {
    // nothing to do in this implementation.
  }

  private:
  Data* m_result ;
  std::size_t m_sample_i ;
  std::size_t m_entry_i ;
} ;

// BgenParser is a thin wrapper around the core functions in genfile/bgen/bgen.hpp.
// This class tracks file state and handles passing the right callbacks.
struct BgenParser {


  BgenParser():
    m_state( e_NotOpen ),
    m_have_sample_ids( false )
  {}

  void open( std::string const& filename ) {
    m_filename =filename;
    m_state = e_NotOpen;
    m_have_sample_ids = false;

    // Open the stream
    m_stream.reset(new std::ifstream( filename, std::ifstream::binary )) ;
    if( !*m_stream ) {
      throw std::invalid_argument( filename ) ;
    }
    m_state = e_Open ;

    // Read the offset, header, and sample IDs if present.
    genfile::bgen::read_offset( *m_stream, &m_offset ) ;
    genfile::bgen::read_header_block( *m_stream, &m_context ) ;
    if( m_context.flags & genfile::bgen::e_SampleIdentifiers ) {
      genfile::bgen::read_sample_identifier_block(
          *m_stream, m_context,
          [this]( std::string id ) { m_sample_ids.push_back( id ) ; }
          ) ;
      m_have_sample_ids = true ;
    }

    // Jump to the first variant data block.
    m_stream->seekg( m_offset + 4 ) ;

    // We keep track of state (though it's not really needed for this implementation.)
    m_state = e_ReadyForVariant ;
    //std::cout << m_state << " " << e_ReadyForVariant << std::endl;
  }

  // modified by J Mbatchou (03/19/21)
  std::ostream& summarise( std::ostream& o ) const {

    o << "   -summary : bgen file ("
      << ( m_context.flags & genfile::bgen::e_Layout2 ? "v1.2 layout" : "v1.1 layout" )
      << ", ";
    if((m_context.flags & genfile::bgen::e_CompressedSNPBlocks) == genfile::bgen::e_ZlibCompression)
      o << "zlib ";
    else if((m_context.flags & genfile::bgen::e_CompressedSNPBlocks) == genfile::bgen::e_ZstdCompression)
      o << "zstd ";
    o << ( m_context.flags & genfile::bgen::e_CompressedSNPBlocks ? "compressed" : "uncompressed" ) << ")"
      << " with " 
      << m_context.number_of_samples << " " << ( m_have_sample_ids ? "named" : "anonymous" ) << " samples and "
      << m_context.number_of_variants << " variants" ;
    return o ;
  }

  // Added by JMbatchou (04/04/21)
  std::string summarise( ) {
    std::ostringstream o; 
    o << "   -summary : bgen file ("
      << ( m_context.flags & genfile::bgen::e_Layout2 ? "v1.2 layout" : "v1.1 layout" )
      << ", ";
    if((m_context.flags & genfile::bgen::e_CompressedSNPBlocks) == genfile::bgen::e_ZlibCompression)
      o << "zlib ";
    else if((m_context.flags & genfile::bgen::e_CompressedSNPBlocks) == genfile::bgen::e_ZstdCompression)
      o << "zstd ";
    o << ( m_context.flags & genfile::bgen::e_CompressedSNPBlocks ? "compressed" : "uncompressed" ) << ")"
      << " with " 
      << m_context.number_of_samples << " " << ( m_have_sample_ids ? "named" : "anonymous" ) << " samples and "
      << m_context.number_of_variants << " variants" ;
    return o.str() ;
  }

  bool get_layout () {
    return (m_context.flags & genfile::bgen::e_Layout2);
  }

  bool get_compression () {
    uint32_t const compressionType = (m_context.flags & genfile::bgen::e_CompressedSNPBlocks);
    return (compressionType == genfile::bgen::e_ZlibCompression) ;
  }

  uint64_t get_position() {
    return m_stream->tellg();
  }

  // added by JMbatchou (12/23/2020)
  void jumpto( uint64_t pos ){
    m_stream->seekg( pos, std::ios_base::beg);
    m_state = e_ReadyForVariant ;
  }

  int number_of_samples() const {
    return m_context.number_of_samples ;
  }

  int number_of_variants() const {
    return m_context.number_of_variants ;
  }
  // Report the sample IDs in the file using the given setter object
  // (If there are no sample IDs in the file, we report a dummy identifier).
  template< typename Setter >
    void get_sample_ids( Setter setter ) {
      if( m_have_sample_ids ) {
        for( std::size_t i = 0; i < m_context.number_of_samples; ++i ) {
          setter( m_sample_ids[i] ) ;
        }
      } else {
        for( std::size_t i = 0; i < m_context.number_of_samples; ++i ) {
          setter( "(unknown_sample_" + std::to_string( i+1 ) + ")" ) ;
        }
      }
    }

  // Attempt to read identifying information about a variant from the bgen file, returning
  // it in the given fields.
  // If this method returns true, data was successfully read, and it should be safe to call read_probs()
  // or ignore_probs().
  // If this method returns false, data was not successfully read indicating the end of the file.
  bool read_variant(
      std::string* chromosome,
      uint32_t* position,
      std::string* rsid,
      std::vector< std::string >* alleles
      ) {
    assert( m_state == e_ReadyForVariant ) ;
    std::string SNPID ; // read but ignored in this toy implementation

    if(
        genfile::bgen::read_snp_identifying_data(
          *m_stream, m_context,
          &SNPID, rsid, chromosome, position,
          [&alleles]( std::size_t n ) { alleles->resize( n ) ; },
          [&alleles]( std::size_t i, std::string const& allele ) { alleles->at(i) = allele ; }
          )
      ) {
      m_state = e_ReadyForProbs ;
      return true ;
    } else {
      return false ;
    }
  }

  // Read genotype probability data for the SNP just read using read_variant()
  // After calling this method it should be safe to call read_variant() to fetch
  // the next variant from the file.
  void read_probs( std::vector< std::vector< double > >* probs ) {
    assert( m_state == e_ReadyForProbs ) ;
    ProbSetter setter( probs ) ;
    genfile::bgen::read_and_parse_genotype_data_block< ProbSetter >(
        *m_stream,
        m_context,
        setter,
        &m_buffer1,
        &m_buffer2
        ) ;
    m_state = e_ReadyForVariant ;
  }

  // Ignore genotype probability data for the SNP just read using read_variant()
  // After calling this method it should be safe to call read_variant()
  // to fetch the next variant from the file.
  void ignore_probs() {
    genfile::bgen::ignore_genotype_data_block( *m_stream, m_context ) ;
    m_state = e_ReadyForVariant ;
  }

  private:
  std::string  m_filename ;
  std::unique_ptr< std::istream > m_stream ;

  // bgen::Context object holds information from the header block,
  // including bgen flags
  genfile::bgen::Context m_context ;

  // offset byte from top of bgen file.
  uint32_t m_offset ;

  // We keep track of our state in the file.
  // Not strictly necessary for this implentation but makes it clear that
  // calls must be read_variant() followed by read_probs() (or ignore_probs())
  // repeatedly.
  enum State { e_NotOpen = 0, e_Open = 1, e_ReadyForVariant = 2, e_ReadyForProbs = 3, eComplete = 4 } ;
  State m_state ;

  // If the BGEN file contains samples ids, they will be read here.
  bool m_have_sample_ids ;
  std::vector< std::string > m_sample_ids ;

  // Buffers, these are used as working space by bgen implementation.
  std::vector< genfile::byte_t > m_buffer1, m_buffer2 ;
} ;
