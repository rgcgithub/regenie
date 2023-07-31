# Makefile for Linux and Mac OSX systems for REGENIE
#
# * User needs to specify BGEN_PATH which is the directory
# 	where the BGEN library is installed
# * If the Boost Iostream library is installed on the system,
# 	user can specify to link to it during compilation by
# 	setting  HAS_BOOST_IOSTREAM to 1
#
###############################################################
#
# Optional: To use external BLAS/LAPACK routines in Eigen	
#
# For Intel MKL, add path of installed library to MKLROOT
# For OpenBLAS, add path of installed library to OPENBLAS_ROOT
#   note: it also requires lapacke library
#  
# For static compilation on Linux systems, set STATIC=1
#   -> this excludes GLIBC
#

BGEN_PATH     =
HAS_BOOST_IOSTREAM := 0
MKLROOT       = 
OPENBLAS_ROOT = 
STATIC       := 0

############

# Use only if not set
CXX          ?= g++
CXXFLAGS      = -O3 -Wall -pedantic -ffast-math -std=c++11 -Wno-unused-local-typedefs -Wno-deprecated-declarations -Wno-long-long -Wno-c11-extensions -fPIC

EFILE         = regenie
CFLAGS       ?=

# check BGEN_PATH is set
ifneq ($(MAKECMDGOALS),clean)
ifeq ($(strip $(BGEN_PATH)),)
 $(warning You did not set the path to BGEN library, i.e. "BGEN_PATH")
endif
endif


# detect OS architecture and add flags
UNAME_S      := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
 INC          = -I${BGEN_PATH}/3rd_party/boost_1_55_0
 CFLAGS      += -fopenmp
 ifeq ($(strip $(STATIC)),1)
  LPATHS      = -static-libgcc -static-libstdc++
  DLIBS       = -Wl,-Bdynamic
 endif
else ifeq ($(UNAME_S),Darwin)
 RGFLAGS     += -stdlib=libc++
endif


RG_VERSION    = $(shell cat VERSION)

## for docker
DFILE         = ./Dockerfile
TEST_SCRIPT   = ./test/test_docker.sh
ifeq ($(strip $(STATIC)),1)
 ifneq ($(strip $(MKLROOT)),)
  DFILE       = ./Dockerfile_mkl # only for static linking
 endif
endif


## for boost iostream
ifeq ($(HAS_BOOST_IOSTREAM),1)
 RG_VERSION  := $(RG_VERSION).gz
 RGFLAGS     += -DHAS_BOOST_IOSTREAM
 ifeq ($(strip $(STATIC)),1)
  SLIBS       = -Wl,-Bstatic -lboost_iostreams
 else
  DLIBS      += -lboost_iostreams
  LIB_BIO2    = libboost-iostreams-dev ## for docker build
 endif
 LIB_BIO      = libboost-iostreams-dev ## for docker build
endif


# Intel MKL or OpenBLAS
ifneq ($(strip $(MKLROOT)),)
 ifeq ($(UNAME_S),Linux)
  RGFLAGS    += -DWITH_MKL -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
  INC        += -I${MKLROOT}/include/
	# static linking
  ifeq ($(strip $(STATIC)),1)
   SLIBS     += -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group
   DLIBS     += -lgomp -lpthread
	# dynamic linking
  else
   LIBMKL     = -L${MKLROOT}/lib/intel64/
   DLIBS     += -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread
  endif
 endif

else ifneq ($(strip $(OPENBLAS_ROOT)),)
 ifeq ($(UNAME_S),Linux)
  RGFLAGS    += -DWITH_OPENBLAS -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
  INC        += -I${OPENBLAS_ROOT}/include/
  # static linking
  ifeq ($(strip $(STATIC)),1)
   SLIBS     += -Wl,-rpath=${OPENBLAS_ROOT}/lib/ -llapack -llapacke -lopenblas
  # dynamic linking
  else
   DLIBS     += -Wl,-rpath=${OPENBLAS_ROOT}/lib/ -llapack -llapacke -lopenblas
  endif
 endif
endif


# pass on version number to software
RGFLAGS      += -DVERSION_NUMBER=\"$(RG_VERSION)\"


OBJECTS       = $(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp))

PGEN_PATH     = ./external_libs/pgenlib/
INC          += -I${PGEN_PATH} -I${PGEN_PATH}/simde/ -I${PGEN_PATH}/include/ -I./external_libs/cxxopts/include/ -I./external_libs/LBFGSpp/include/ -I${BGEN_PATH} -I./external_libs/eigen-3.4.0/ -I${BGEN_PATH}/genfile/include/ -I${BGEN_PATH}/3rd_party/boost_1_55_0/ -I${BGEN_PATH}/3rd_party/zstd-1.1.0/lib -I${BGEN_PATH}/db/include/ -I${BGEN_PATH}/3rd_party/sqlite3 -I./external_libs/

LPATHS       += ${LIBMKL} -L${BGEN_PATH}/build/ -L${BGEN_PATH}/build/3rd_party/zstd-1.1.0/ -L${BGEN_PATH}/build/db/ -L${BGEN_PATH}/build/3rd_party/sqlite3/ -L${BGEN_PATH}/build/3rd_party/boost_1_55_0 -L/usr/lib/

LIBS         += ${SLIBS} -lbgen -lzstd -ldb  -lsqlite3 -lboost
LIBS         += -lz ${DLIBS} -lm -ldl -lgfortran



.PHONY: docker-build docker-test debug clean

all: ${EFILE}

${EFILE}: libMvtnorm libqf libquad pgenlib ${OBJECTS}
	${CXX} ${CXXFLAGS} ${RGFLAGS} ${CFLAGS} -o ${EFILE} ${OBJECTS} ./external_libs/mvtnorm/libMvtnorm.a ./external_libs/qf/qf.a ./external_libs/quadpack/libquad.a ./external_libs/pgenlib/pgenlib.a ${LPATHS} ${LIBS}

%.o: %.cpp
	${CXX} ${CXXFLAGS} ${RGFLAGS} -o $@ -c $< ${INC} ${CFLAGS}

libMvtnorm: 
		(cd ./external_libs/mvtnorm/;$(MAKE))

libqf: 
		(cd ./external_libs/qf/;$(MAKE))

libquad: 
		(cd ./external_libs/quadpack/;$(MAKE))

pgenlib: 
		(cd ./external_libs/pgenlib/;$(MAKE))

#####
## For use with Docker
# create Docker image
docker-build:
	@echo "Building docker image for REGENIE v${RG_VERSION}"
ifeq ($(HAS_BOOST_IOSTREAM),1)
	@echo Compiling with Boost Iostream library
endif
ifeq ($(STATIC),1)
 ifneq ($(strip $(MKLROOT)),)
	@echo Compiling with Intel MKL library
 endif
	@echo Linking = static
 else
	@echo Linking = dynamic
endif

	@docker build --rm -f ${DFILE} \
		--no-cache --pull \
		--build-arg BOOST_IO=${HAS_BOOST_IOSTREAM} \
		--build-arg LIB_INSTALL=${LIB_BIO} \
		--build-arg LIB_INSTALL2=${LIB_BIO2} \
		--build-arg STATIC=${STATIC} \
		-t regenie:v${RG_VERSION} .

docker-test:
	@${TEST_SCRIPT} . "regenie:v${RG_VERSION}" ${HAS_BOOST_IOSTREAM}
####


debug: CXXFLAGS  = -O0 -g -std=c++11 -fPIC
debug: ${EFILE}

clean:
	rm -rf ${EFILE} ./src/*.o ./build/
	(cd ./external_libs/mvtnorm/;$(MAKE) clean)
	(cd ./external_libs/qf/;$(MAKE) clean)
	(cd ./external_libs/quadpack/;$(MAKE) clean)
	(cd ./external_libs/pgenlib/;$(MAKE) clean)
