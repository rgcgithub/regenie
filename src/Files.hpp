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
#ifndef RFILES_H
#define RFILES_H

#include <fstream>
#include <iostream>
#include <boost/filesystem.hpp>

# if defined(HAS_BOOST_IOSTREAM)
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#endif

class Files {

  public:
    // variables
    bool is_gz;
    bool read_mode = true;

    // for reading
    std::ifstream myfile;
# if defined(HAS_BOOST_IOSTREAM)
    boost::iostreams::filtering_istream mygzfile;
#endif


    // functions
    bool checkFileExtension(std::string filename);
    void openForRead(std::string filename,mstream&);
    bool readLine(std::string& line);
    void ignoreLines(int);
    void closeFile();


    Files();
    ~Files();
};

#endif
