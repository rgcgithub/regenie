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

#include <boost/filesystem.hpp>

# if defined(HAS_BOOST_IOSTREAM)
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#endif

class Files {

  public:
    // variables
    bool is_gz = false;
    bool read_mode = true;

    // for reading
    std::ifstream infile;
    // for writing
    std::ofstream outfile;

# if defined(HAS_BOOST_IOSTREAM)
    boost::iostreams::filtering_istream ingzfile;
    boost::iostreams::filtering_ostream outgzfile;
#endif


    // functions
    bool checkFileExtension(std::string filename);
    void openForRead(std::string filename,mstream&);
    bool readLine(std::string& line);
    void ignoreLines(int);
    void openForWrite(std::string filename,mstream&);
    void closeFile();
    void openBinMode(std::string filename,std::ios_base::openmode,mstream&);
    void writeBinMode(Eigen::ArrayXi&,mstream&);
    void writeBinMode(ArrayXt&,mstream&);

    // to write to file
    template <class S>
      Files& operator<< (const S& val)
      {
# if defined(HAS_BOOST_IOSTREAM)
        if(is_gz) {
          outgzfile << val;
          return *this;
        }
#endif
        outfile << val;
        return *this;
      }

    // for std::endl
    Files& operator<< (std::ostream& (*pfun)(std::ostream&))
    {
# if defined(HAS_BOOST_IOSTREAM)
        if(is_gz) {
          pfun(outgzfile);
          return *this;
        }
#endif
      pfun(outfile);
      return *this;
    };

    Files();
    ~Files();
};

std::vector<std::string> string_split(std::string const&,const char*);
void openStream_write(std::ofstream*,std::string const&,std::ios_base::openmode,mstream&);

#endif
