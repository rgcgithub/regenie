# Thanks to Nathan Weeks for sharing this Dockerfile
# minor edits were made
# Filename: Dockerfile

# make this global 
ARG LIB_INSTALL
ARG LIB_INSTALL2

FROM ubuntu:18.04 AS builder

ARG BOOST_IO
ARG LIB_INSTALL
ARG STATIC
ARG CMAKE_VERSION_MAJOR=3.13
ARG CMAKE_VERSION_MINOR=0
ARG HTSLIB_VERSION=1.18

WORKDIR /src

ADD http://cmake.org/files/v${CMAKE_VERSION_MAJOR}/cmake-${CMAKE_VERSION_MAJOR}.${CMAKE_VERSION_MINOR}-Linux-x86_64.sh cmake_install.sh
ADD http://code.enkre.net/bgen/tarball/release/v1.1.7 v1.1.7.tgz
ADD https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2 htslib-$HTSLIB_VERSION.tar.bz2

# install BGEN and HTSlib libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
      gcc \
      g++ \
      make \
      libz-dev \
      libbz2-dev \
      liblzma-dev \
      libcurl4-openssl-dev \
      libssl-dev \
      python3 \
      gfortran \
      zlib1g-dev \
      $LIB_INSTALL \
      && tar -xf htslib-$HTSLIB_VERSION.tar.bz2 \
      && cd htslib-$HTSLIB_VERSION/ \
      && ./configure \
      && make \
      && make install \
      && cd .. \
      && sh cmake_install.sh --prefix=/usr/local --skip-license --exclude-subdir \
      && rm cmake_install.sh \
      && tar -xzf v1.1.7.tgz \
      && rm v1.1.7.tgz \
      && cd v1.1.7 \
      && python3 waf configure \
      && python3 waf

COPY . /src/regenie

WORKDIR /src/regenie

RUN BGEN_PATH=/src/v1.1.7 HAS_BOOST_IOSTREAM=$BOOST_IO HTSLIB_PATH=/usr/local/lib/ STATIC=$STATIC cmake . \
      && make

FROM ubuntu:18.04
ARG LIB_INSTALL2

RUN apt-get update && apt-get install -y --no-install-recommends \
      libgomp1 gfortran $LIB_INSTALL2 \
      && rm -rf /var/lib/apt/lists/*

COPY --from=builder /src/regenie/regenie /usr/local/bin