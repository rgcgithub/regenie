#!/bin/bash

set -ex

# https://bioconda.github.io/troubleshooting.html#zlib-errors
export CFLAGS="-I$PREFIX/include"
export LDFLAGS="-L$PREFIX/lib"
export CPATH=${PREFIX}/include

export MKL_THREADING_LAYER="GNU"

mkdir -p build

cmake \
  -DBUILD_SHARED_LIBS:BOOL=ON \
  -DCMAKE_PREFIX_PATH:PATH=${PREFIX} \
  -DCMAKE_INSTALL_PREFIX:PATH=${PREFIX} \
  -DCMAKE_BUILD_TYPE="Release" \
  -S "${SRC_DIR}" \
  -B build

make -C build -j${CPU_COUNT:=1} regenie
make -C build install

# bash test/test_conda.sh --path "${SRC_DIR}" --gz
