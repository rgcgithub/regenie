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


#include "Regenie.hpp"
#include "Files.hpp"

namespace fs = boost::filesystem;

Files::Files(){
}
Files::~Files(){
}

// Open file (either regular or gzipped)
void Files::openForRead(std::string const& filename, mstream& sout){

  read_mode = true;
  is_gz = isGzipped(filename, true);
  //std::cerr << filename << " - gzip = " << std::boolalpha << is_gz << std::endl;

  // only used if compiled with boost iostream
# if not defined(HAS_BOOST_IOSTREAM)
  if(is_gz) 
    throw "cannot read gzip file if compilation is not done with the Boost Iostream library (i.e. 'make HAS_BOOST_IOSTREAM=1').";
#endif

  std::ios_base::openmode mode = (is_gz ? std::ios_base::in | std::ios_base::binary : std::ios_base::in ); 

  openStream(&infile, filename, mode, sout);

# if defined(HAS_BOOST_IOSTREAM)
  if(is_gz){
    ingzfile.push(boost::iostreams::gzip_decompressor());
    ingzfile.push(infile);
  }
#endif

}

bool Files::readLine(std::string& line){

# if defined(HAS_BOOST_IOSTREAM)
  if(is_gz) 
    return static_cast<bool>( getline(ingzfile, line) );
#endif

  return  static_cast<bool>( getline(infile, line) );
}


void Files::ignoreLines(int const& nlines){

  int linenumber=0;

  if(nlines < 1) return;

  while(linenumber++ < nlines){

    if(is_gz) {
# if defined(HAS_BOOST_IOSTREAM)
      ingzfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
#endif
    } else 
      infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  }
}


// Open file for writing
void Files::openForWrite(std::string const& filename, mstream& sout){

  read_mode = false;
  is_gz = isGzipped(filename, false);

  // only used if compiled with boost iostream
# if not defined(HAS_BOOST_IOSTREAM)
  if(is_gz) 
    throw "cannot write gzip file if compilation is not done with the Boost Iostream library (i.e. 'make HAS_BOOST_IOSTREAM=1').";
#endif

  std::ios_base::openmode mode = (is_gz ? std::ios_base::out | std::ios_base::binary : std::ios_base::out ); 

  openStream(&outfile, filename, mode, sout);

# if defined(HAS_BOOST_IOSTREAM)
  if(is_gz){
    outgzfile.push(boost::iostreams::gzip_compressor());
    outgzfile.push(outfile);
  }
#endif

}

void Files::closeFile(){

  if( read_mode ){

# if defined(HAS_BOOST_IOSTREAM)
    if(is_gz) 
      ingzfile.reset();
#endif
    infile.close();

  } else {

# if defined(HAS_BOOST_IOSTREAM)
    if(is_gz) 
      outgzfile.reset();
#endif
    outfile.close();

  }
}

// Check file extension
bool Files::isGzipped(std::string const& filename, bool const& check_file) {

  // require all gzipped file to end in .gz
  if( fs::extension(filename) != ".gz" )
    return false;

  // open file and check first 2 bytes (should equal 0x1f8b)
  if(check_file){
    infile.open(filename, std::ios_base::in | std::ios_base::binary);
    if (infile.fail()) 
      throw "cannot read file : " + filename ;

    uchar header[2];
    infile.read( reinterpret_cast<char *> (&header[0]), 2);
    infile.close();

    if ( (header[0] != 0x1f) || (header[1] != 0x8b) ) 
      return false;
  }

  return true;
}

void Files::openMode(std::string const& filename, std::ios_base::openmode mode, mstream& sout){

  if(mode & std::ios_base::out){
    read_mode = false;
    openStream(&outfile, filename, mode, sout);
  } else {
    read_mode = true;
    openStream(&infile, filename, mode, sout);
  }

}

// Split string by tokens
std::vector<std::string> string_split(std::string const& s, const char* delims) {

  std::vector<std::string> out;

  if(s.size() == 0) return out;

  const char* p = s.c_str(); //beginning of string
  const char* q = strpbrk(p+1, delims);//to first delimiter

  for( ; q != NULL; q = strpbrk(p, delims)){
    out.push_back( std::string(p,q) );// add to vector using range constructor
    p = q + 1;
  }

  // check string after last delimiter
  if(p && (p[0] != '\0')) out.push_back( std::string(p) );

  return(out);

}


