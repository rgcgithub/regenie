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

// Open file (either regular or gzipped
void Files::open(std::string filename, mstream& sout){

  // only used if compiled with boost iostream
# if defined(HAS_BOOST_IOSTREAM)
  is_gz = checkFileExtension(filename);
#else
  is_gz = false;
#endif

  std::ios_base::openmode mode = (is_gz ? std::ios_base::in | std::ios_base::binary : std::ios_base::in ); 

  myfile.open(filename.c_str(), mode);
  if (!myfile.is_open()) {    
    sout << "ERROR: Cannot open file : " << filename << std::endl;
    exit(-1);
  }

# if defined(HAS_BOOST_IOSTREAM)
  if(is_gz){
    mygzfile.push(boost::iostreams::gzip_decompressor());
    mygzfile.push(myfile);
  }
#endif

}

bool Files::readLine(std::string& line){

# if defined(HAS_BOOST_IOSTREAM)
  if(is_gz) return static_cast<bool>( getline(mygzfile, line) );
#endif

  return  static_cast<bool>( getline(myfile, line) );
}

void Files::closeFile(){
# if defined(HAS_BOOST_IOSTREAM)
  if(is_gz) mygzfile.reset();
#endif
  myfile.close();
}

// Check file extension
bool Files::checkFileExtension(std::string filename) {
  return fs::extension(filename) == ".gz";
}
