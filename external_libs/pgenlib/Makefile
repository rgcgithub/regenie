#
# Makefile to compile PGEN library
#
AR         ?= ar
CXX        ?= g++
CXXFLAGS    = -O3 -Wall -std=c++11
CFLAGS      = 

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
 CXXFLAGS += -arch x86_64
endif

OBJECTS       = $(patsubst %.cc,%.o,$(wildcard ./include/*.cc)) $(patsubst %.cpp,%.o,$(wildcard *.cpp))  $(patsubst %.cc,%.o,$(wildcard *.cc))
INC           =  -I./simde/ -I./include/

pgenlib.a: ${OBJECTS}
	${AR} rcs $@ $^ 

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

%.o: %.cc
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

clean:
	rm -f *.o *.a ./include/*.o
