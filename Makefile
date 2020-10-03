# Makefile for Linux and Mac OSX systems for REGENIE
#
# * User needs to specify BGEN_PATH which is the directory
# 	where the BGEN library is installed
# * If the Boost Iostream library is installed on the system,
# 	user can specify to link to it during compilation by
# 	setting  HAS_BOOST_IOSTREAM to 1
#
#
# Optional: To use external BLAS/LAPACK routines in Eigen	
# * If Intel MKL is installed on system, add path to MKLROOT
# * Else if OpenBLAS is installed on system, add path to OPENBLAS_ROOT
#
#

BGEN_PATH     =
HAS_BOOST_IOSTREAM := 0
MKLROOT       = 
OPENBLAS_ROOT = 

############

CXX           = g++
CXXFLAGS      = -O3 -Wall -ffast-math -std=c++11
EFILE         = regenie
CFLAGS        =


# detect OS architecture and add flags
UNAME_S      := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
  INC         = -I${BGEN_PATH}/3rd_party/boost_1_55_0
  CFLAGS     += -fopenmp
else ifeq ($(UNAME_S),Darwin)
  RGFLAGS    += -stdlib=libc++
endif


RG_VERSION    = $(shell cat VERSION)

## for docker
DFILE         = ./Dockerfile
TEST_SCRIPT   = ./test/test_docker.sh

## for boost iostream
ifeq ($(HAS_BOOST_IOSTREAM),1)
  RG_VERSION := $(RG_VERSION).gz
  RGFLAGS    += -DHAS_BOOST_IOSTREAM
	LBIO        = -lboost_iostreams
  LIB_BIO     = libboost-iostreams-dev ## for docker build
endif


# Intel MKL or OpenBLAS
ifneq ($(strip $(MKLROOT)),)
	ifeq ($(UNAME_S),Linux)
		RGFLAGS  += -DWITH_MKL -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
		INC      += -I${MKLROOT}/include/
		LPATHS    = -L${MKLROOT}/lib/intel64
		LLAPACK   = -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -lblas -llapack -llapacke
	endif
else ifneq ($(strip $(OPENBLAS_ROOT)),)
	ifeq ($(UNAME_S),Linux)
		RGFLAGS  += -DWITH_OPENBLAS -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
		INC      += -I${OPENBLAS_ROOT}/include/
		LLAPACK   = -Wl,-rpath,${OPENBLAS_ROOT}/lib/ -llapack -llapacke -lopenblas -lgfortran
	endif
endif


# pass on version number to software
RGFLAGS      += -DVERSION_NUMBER=\"$(RG_VERSION)\"


PGEN_PATH     = ./external_libs/pgenlib/
PGEN_OBJECTS  = $(patsubst %.cc,%.o,$(wildcard ${PGEN_PATH}include/*.cc)) $(patsubst %.cpp,%.o,$(wildcard ${PGEN_PATH}*.cpp))
OBJECTS       = $(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp)) ${PGEN_OBJECTS}


INC          += -I${PGEN_PATH} -I${PGEN_PATH}/include/ -I${BGEN_PATH} -I${BGEN_PATH}/genfile/include/ -I${BGEN_PATH}/3rd_party/zstd-1.1.0/lib -I${BGEN_PATH}/db/include/ -I${BGEN_PATH}/3rd_party/sqlite3 -I./external_libs/

LPATHS        = -L${BGEN_PATH}/build/ -L${BGEN_PATH}/build/3rd_party/zstd-1.1.0/ -L${BGEN_PATH}/build/db/ -L${BGEN_PATH}/build/3rd_party/sqlite3/ -L${BGEN_PATH}/build/3rd_party/boost_1_55_0 -L/usr/lib/

LIBS          = -lbgen -lzstd -ldb  -lsqlite3 -lboost ${LBIO} ${LLAPACK}
LIBS         += -ldl -lz



.PHONY: docker-build docker-test debug clean

all: ${EFILE}

${EFILE}: ${OBJECTS}
	${CXX} ${CXXFLAGS} ${RGFLAGS} ${CFLAGS} -o ${EFILE} ${OBJECTS} ${LPATHS} ${LIBS}

%.o: %.cpp
	${CXX} ${CXXFLAGS} ${RGFLAGS} -o $@ -c $< ${INC} ${CFLAGS}

%.o: %.cc
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC} ${CFLAGS}


#####
## For use with Docker
# create Docker image
docker-build:
	@echo "Building docker image for REGENIE v${RG_VERSION}"
ifeq ($(HAS_BOOST_IOSTREAM),1)
	@echo Compiling with Boost Iostream library
endif
	@docker build -f ${DFILE} \
		--no-cache --pull \
		--build-arg BOOST_IO=${HAS_BOOST_IOSTREAM} \
		--build-arg LIB_INSTALL=${LIB_BIO} \
		-t regenie:v${RG_VERSION} .

docker-test:
	@${TEST_SCRIPT} . "regenie:v${RG_VERSION}" ${HAS_BOOST_IOSTREAM}
####


debug: CXXFLAGS  = -O0 -g
debug: ${EFILE}

clean:
	rm -f ${EFILE} ./src/*.o ${PGEN_PATH}/*.o ${PGEN_PATH}/include/*.o
