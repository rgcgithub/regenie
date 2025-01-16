/* bgz_writer.hpp
* Author: Tyler Joseph
* 
* HTSlib wrapper to write bgzip files.
* 
* Example:
*   // Create a Bgz file for writing:
*   BgzWriter writer("myfile.gz", "w");
*
*   // Write to an uncompressed file:
*   BgzWriter writer("myfile", "wu");
* 
*   // Or append to a file:
*   BgzWriter writer("myfile.gz", "a");
* 
*   // Write to the file:
*   writer.write(string);
* 
*   // Close:
*   writer.close()
*/
#ifdef WITH_HTSLIB
/* bgz_writer.hpp
* Author: Tyler Joseph
* 
* HTSlib wrapper to write bgzip files.
* 
* Example:
*   // Create a Bgz file for writing:
*   BgzWriter writer = BgzWriter("myfile.gz", "w");
*
*   // Write to an uncompressed file:
*   BgzWriter writer = BgzWriter("myfile", "wu");
* 
*   // Or append to a file:
*   BgzWriter writer = BgzWriter("myfile.gz", "a");
* 
*   // Write to the file:
*   writer.write(string);
* 
*   // Close:
*   writer.close()
*/

#ifndef BGZ_WRITER_H
#define BGZ_WRITER_H

#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>
#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

class BgzWriter {
 public:
  BgzWriter();
  BgzWriter(std::string filepath, std::string mode);
  ~BgzWriter();

  /*
    We don't want to have duplicates of the same file open for writing, so
    we need to delete the default copy constructor and assignment operator.
  */
  BgzWriter(const BgzWriter& other) = delete;
  BgzWriter& operator=(BgzWriter other) = delete;

  // Allow BgzWriter to be placed in containers using its move constructor
  BgzWriter(BgzWriter&& other);

  void write(std::string s);

  // this gets overridden in ld_matrix_writer.hpp, but
  // it's annoying that we have to make these a virtual function
  virtual void open(std::string filepath, std::string mode);
  virtual void close();

  bool is_closed() { return closed; }

  int64_t tell() { return bgzf_tell(this->bgzf); }

  template<typename T>
  void write(T data) {
    ssize_t size = sizeof(data);
    if (this->closed) {
      throw std::runtime_error("attempted to write closed file " + this->filepath);
    }

    ssize_t bytes_written = bgzf_write(
      this->bgzf,
      static_cast<char *>(static_cast<void *>(&data)),
      size
    );
    if (bytes_written != size) {
      throw std::runtime_error("failed to write " + this->filepath);
    }
  }

 private:
  std::string filepath;
  std::string mode;
  BGZF* bgzf;
  kstring_t buffer;
  bool closed;
};

#endif
#endif