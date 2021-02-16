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

namespace fs = boost::filesystem;

Files::Files(){
}
Files::~Files(){
}

// Open file (either regular or gzipped)
void Files::openForRead(std::string filename, mstream& sout){

  is_gz = checkFileExtension(filename);

  // only used if compiled with boost iostream
# if not defined(HAS_BOOST_IOSTREAM)
  if(is_gz) {
    sout << "ERROR: Cannot read gzip file if compilation is not done with the Boost Iostream library (i.e. 'make HAS_BOOST_IOSTREAM=1').\n";
    exit(EXIT_FAILURE);
  }
#endif

  std::ios_base::openmode mode = (is_gz ? std::ios_base::in | std::ios_base::binary : std::ios_base::in ); 

  infile.open(filename.c_str(), mode);
  if (!infile) {    
    sout << "ERROR: Cannot open file : " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

# if defined(HAS_BOOST_IOSTREAM)
  if(is_gz){
    ingzfile.push(boost::iostreams::gzip_decompressor());
    ingzfile.push(infile);
  }
#endif

}

bool Files::readLine(std::string& line){

# if defined(HAS_BOOST_IOSTREAM)
  if(is_gz) return static_cast<bool>( getline(ingzfile, line) );
#endif

  return  static_cast<bool>( getline(infile, line) );
}


void Files::ignoreLines(int nlines){

  int linenumber=0;
  if(nlines < 1) return;

  while(linenumber++ < nlines){
# if defined(HAS_BOOST_IOSTREAM)
    if(is_gz) ingzfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    else infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
# else
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
#endif
  }
}


// Open file for writing
void Files::openForWrite(std::string filename, mstream& sout){

  read_mode = false;
  is_gz = checkFileExtension(filename);

  // only used if compiled with boost iostream
# if not defined(HAS_BOOST_IOSTREAM)
  if(is_gz) {
    sout << "ERROR: Cannot write gzip file if compilation is not done with the Boost Iostream library (i.e. 'make HAS_BOOST_IOSTREAM=1').\n";
    exit(EXIT_FAILURE);
  }
#endif

  std::ios_base::openmode mode = (is_gz ? std::ios_base::out | std::ios_base::binary : std::ios_base::out ); 

  outfile.open(filename.c_str(), mode);
  if (!outfile) {    
    sout << "ERROR: Cannot open file : " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

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
    if(is_gz) ingzfile.reset();
#endif
    infile.close();
  } else {
# if defined(HAS_BOOST_IOSTREAM)
    if(is_gz) outgzfile.reset();
#endif
    outfile.close();
  }
}

// Check file extension
bool Files::checkFileExtension(std::string filename) {
  return fs::extension(filename) == ".gz";
}

void Files::openBinMode(std::string filename, std::ios_base::openmode mode, mstream& sout){

  try {
    if(mode & std::ios_base::out){
      read_mode = false;
      outfile.open(filename.c_str(), mode);
      if (!outfile) throw filename;    
    } else {
      infile.open(filename.c_str(), mode);
      if (!infile) throw filename;    
    }
  } catch (const std::string& fname) {
    sout << "ERROR: Cannot open file : " << fname << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Files::writeBinMode(Eigen::ArrayXi& vals, mstream& sout){

  outfile.write( reinterpret_cast<char *> (&vals(0)), vals.size() * sizeof(vals(0)) );
  if (outfile.fail()) {    
    sout << "ERROR: Cannot write values to file.\n";
    exit(EXIT_FAILURE);
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


void openStream_write(std::ofstream* ofs, std::string const& fname, std::ios_base::openmode mode, mstream& sout){

  ofs->open(fname, mode);
  if (ofs->fail()) {    
    sout << "ERROR: Cannot write to file : " << fname << std::endl;
    exit(EXIT_FAILURE);
  }

  return;
}
