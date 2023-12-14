#ifdef WITH_HTSLIB
#include "bgz_writer.hpp"

using namespace std;

BgzWriter::BgzWriter()
 : filepath("")
 , mode("")
 , bgzf(nullptr)
 , closed(true) {
  this->buffer = KS_INITIALIZE;
}

BgzWriter::BgzWriter(string filepath, string mode)
 : filepath(filepath)
 , mode(mode)
 , closed(false) {
  if (mode != "w" && mode != "a" && mode != "wu" && mode != "au") {
    throw runtime_error("invalid write mode " + mode + " in BgzWriter");
  }
  this->bgzf = bgzf_open(filepath.c_str(), mode.c_str());
  this->buffer = KS_INITIALIZE;
  if (this->bgzf == NULL) {
    throw runtime_error("failed to open "+ filepath);
  }
}

BgzWriter::~BgzWriter() {
  if (this->bgzf != nullptr && !this->is_closed()) {
    bgzf_close(this->bgzf);
  }
  ks_free(&this->buffer);
}

BgzWriter::BgzWriter(BgzWriter&& other)
 : filepath(std::move(other.filepath))
 , mode(std::move(other.mode))
 , bgzf(std::exchange(other.bgzf, nullptr))
 , closed(std::move(other.closed)) {
  this->buffer = KS_INITIALIZE;
}

void BgzWriter::write(string s) {
  if (this->closed) {
    throw runtime_error("attempted to write to closed file " + this->filepath);
  }

  ssize_t bytes_written = bgzf_write(this->bgzf, s.c_str(), s.size());
  if (bytes_written != (int)s.size()) {
    throw runtime_error("failed to write " + this->filepath);
  }
}

void BgzWriter::open(string filepath, string mode) {
  this->close();
  this->filepath = filepath;
  this->mode = mode;
  this->closed = false;
  this->bgzf = bgzf_open(filepath.c_str(), mode.c_str());
  if (this->bgzf == NULL) {
    throw runtime_error("failed to open "+ filepath);
  }
}

void BgzWriter::close() {
  bgzf_close(this->bgzf);
  this->bgzf = NULL;
  this->closed = true;
}
#endif