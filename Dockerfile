# Thanks to Nathan Weeks for sharing this Dockerfile
# minor edits were made
# Filename: Dockerfile

# make this global 
ARG LIB_INSTALL
ARG LIB_INSTALL2


FROM ubuntu:18.04 AS builder

ARG LIB_INSTALL
ARG SHARED=0
ARG CMAKE_VER=3.17.5
ARG BGEN_PATH=/src/v1.1.7
ARG SRC_DIR=/src/regenie
ARG CPU_COUNT=1

WORKDIR /src

ADD http://code.enkre.net/bgen/tarball/release/v1.1.7 v1.1.7.tgz
ADD https://github.com/Kitware/CMake/releases/download/v${CMAKE_VER}/cmake-${CMAKE_VER}-Linux-x86_64.sh cmake.sh

RUN apt-get update && apt-get install -y --no-install-recommends \
      g++ \
      make \
      python3 \
      zlib1g-dev \
      libeigen3-dev \
      $LIB_INSTALL \
      && bash cmake.sh --prefix=/usr/local --skip-license \
      && rm cmake.sh \
      && tar -xzf v1.1.7.tgz \
      && rm v1.1.7.tgz \
      && cd v1.1.7 \
      && python3 waf configure \
      && python3 waf

COPY . /src/regenie

WORKDIR /src/regenie/build

RUN cmake \
      -DBUILD_SHARED_LIBS:BOOL="${SHARED}" \
      -DCMAKE_PREFIX_PATH:PATH=/usr/local \
      -DCMAKE_INSTALL_PREFIX:PATH=/usr/local \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_BUILD_TYPE="Release" \
      -DBGEN_PATH="${BGEN_PATH}" \
      -DWITH_MKL:BOOL=OFF \
      -DWITH_OPENBLAS:BOOL=OFF \
      -S "${SRC_DIR}" \
    && make VERBOSE=1 -j${CPU_COUNT} regenie \
    && make install

FROM ubuntu:18.04
ARG LIB_INSTALL2

RUN apt-get update && apt-get install -y --no-install-recommends \
      libgomp1 $LIB_INSTALL2 \
      && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/bin/regenie /usr/local/bin

# Avoid this to keep image more for general usage
# ENTRYPOINT ["/usr/local/bin/regenie"]

