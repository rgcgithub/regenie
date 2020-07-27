CXX = g++
CXXFLAGS = -O3 -w -ffast-math -std=c++11
EFILE = regenie
CFLAGS = 

# specify BGEN library path 
BGEN_PATH = 

# detect OS architecture and add flags
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	INC = -I${BGEN_PATH}/3rd_party/boost_1_55_0
	CFLAGS += -fopenmp
endif
ifeq ($(UNAME_S),Darwin)
	CXXFLAGS += -stdlib=libc++
endif


PGEN_PATH = ./external_libs/pgenlib/
PGEN_OBJECTS = $(patsubst %.cc,%.o,$(wildcard ${PGEN_PATH}include/*.cc)) $(patsubst %.cpp,%.o,$(wildcard ${PGEN_PATH}*.cpp))   

LPATHS = -L${BGEN_PATH}/build/ -L${BGEN_PATH}/build/3rd_party/zstd-1.1.0/ -L${BGEN_PATH}/build/db/ -L${BGEN_PATH}/build/3rd_party/sqlite3/ -L${BGEN_PATH}/build/3rd_party/boost_1_55_0 -L/usr/lib/
LIBS = -lbgen -lzstd -ldb  -lsqlite3 -lboost -lz
INC += -I${PGEN_PATH} -I${PGEN_PATH}/include/ -I${BGEN_PATH} -I${BGEN_PATH}/genfile/include/ -I${BGEN_PATH}/3rd_party/zstd-1.1.0/lib -I./external_libs/ 

OBJECTS = $(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp)) ${PGEN_OBJECTS} 


all: ${EFILE}

${EFILE}: ${OBJECTS}
	${CXX} ${CXXFLAGS} ${CFLAGS} -o ${EFILE} ${OBJECTS} ${LPATHS} ${LIBS} 

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC} ${CFLAGS}

%.o: %.cc
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC} ${CFLAGS}

debug: CXXFLAGS = -O0 -g
debug: ${EFILE}

clean:
	rm -f ${EFILE} ./src/*.o ${PGEN_PATH}/*.o ${PGEN_PATH}/include/*.o
