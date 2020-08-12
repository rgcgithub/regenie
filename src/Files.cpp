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

  // only used if compiled with boost iostream
# if defined(HAS_BOOST_IOSTREAM)
  is_gz = checkFileExtension(filename);
#else
  is_gz = false;
#endif

  std::ios_base::openmode mode = (is_gz ? std::ios_base::in | std::ios_base::binary : std::ios_base::in ); 

  infile.open(filename.c_str(), mode);
  if (!infile) {    
    sout << "ERROR: Cannot open file : " << filename << std::endl;
    exit(-1);
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
  // only used if compiled with boost iostream
# if defined(HAS_BOOST_IOSTREAM)
  is_gz = checkFileExtension(filename);
#else
  is_gz = false;
#endif

  std::ios_base::openmode mode = (is_gz ? std::ios_base::out | std::ios_base::binary : std::ios_base::out ); 

  outfile.open(filename.c_str(), mode);
  if (!outfile) {    
    sout << "ERROR: Cannot open file : " << filename << std::endl;
    exit(-1);
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
