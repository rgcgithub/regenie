CXX = g++
CXXFLAGS = -O3 -w -ffast-math -std=c++11
EFILE = regenie
CFLAGS = 

# specify BGEN library path 
BGEN_PATH = 
# specify if compiling Boost iostream library (set to 1 if yes)
HAS_BOOST_IOSTREAM := 0


# detect OS architecture and add flags
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	INC = -I${BGEN_PATH}/3rd_party/boost_1_55_0
	CFLAGS += -fopenmp
else ifeq ($(UNAME_S),Darwin)
	CXXFLAGS += -stdlib=libc++
endif


PGEN_PATH = ./external_libs/pgenlib/
PGEN_OBJECTS = $(patsubst %.cc,%.o,$(wildcard ${PGEN_PATH}include/*.cc)) $(patsubst %.cpp,%.o,$(wildcard ${PGEN_PATH}*.cpp))   

LPATHS = -L${BGEN_PATH}/build/ -L${BGEN_PATH}/build/3rd_party/zstd-1.1.0/ -L${BGEN_PATH}/build/db/ -L${BGEN_PATH}/build/3rd_party/sqlite3/ -L${BGEN_PATH}/build/3rd_party/boost_1_55_0 -L/usr/lib/
LIBS = -lbgen -lzstd -ldb  -lsqlite3 -lboost -lz
INC += -I${PGEN_PATH} -I${PGEN_PATH}/include/ -I${BGEN_PATH} -I${BGEN_PATH}/genfile/include/ -I${BGEN_PATH}/3rd_party/zstd-1.1.0/lib -I./external_libs/ 


## specify dockerfile 
REGENIE_VERSION := v1.0.5.5
DFILE = ./Dockerfile
TEST_SCRIPT=./test/test.sh

## add boost iostream
ifeq ($(HAS_BOOST_IOSTREAM),1)
	LIBS += -lboost_iostreams
	CXXFLAGS += -DHAS_BOOST_IOSTREAM
	REGENIE_VERSION := $(REGENIE_VERSION).gz
	LIB_BIO = libboost-iostreams-dev ## for docker build
endif

OBJECTS = $(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp)) ${PGEN_OBJECTS} 

.PHONY: all docker-build docker-test debug clean

all: ${EFILE}

${EFILE}: ${OBJECTS}
	${CXX} ${CXXFLAGS} ${CFLAGS} -o ${EFILE} ${OBJECTS} ${LPATHS} ${LIBS}

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC} ${CFLAGS}

%.o: %.cc
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC} ${CFLAGS}

#####
## For use with Docker
# create Docker image
docker-build:	
	@echo "Building docker image for REGENIE ${REGENIE_VERSION}"
ifeq ($(HAS_BOOST_IOSTREAM),1)
	@echo Compiling with Boost Iostream library
endif
	@docker build -f ${DFILE} \
		--no-cache --pull \
		--build-arg BOOST_IO=${HAS_BOOST_IOSTREAM} \
		--build-arg LIB_INSTALL=${LIB_BIO} \
		-t regenie:${REGENIE_VERSION} .

docker-test:	
	@${TEST_SCRIPT} . "regenie:${REGENIE_VERSION}" ${HAS_BOOST_IOSTREAM}
####

debug: CXXFLAGS = -O0 -g
debug: ${EFILE}

clean:
	rm -f ${EFILE} ./src/*.o ${PGEN_PATH}/*.o ${PGEN_PATH}/include/*.o
