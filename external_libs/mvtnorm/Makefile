FLAG  =-O2 -DNDEBUG -fPIC
AR   ?= ar
CXX  ?= g++
CC   ?= gcc
# make defines FC as `f77` itself, so `?=` would never take effect here.
ifeq ($(origin FC),default)
FC   := gfortran
endif

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_S),Darwin)
 ifeq ($(UNAME_M),x86_64)
  FLAG += -arch x86_64
 endif
endif

libMvtnorm.a: mvtnorm.o mvt.o randomF77.o
	${AR} rcs $@ $^ 

mvtnorm.o: mvtnorm.cpp mvtnorm.h
	${CXX} $(FLAG) -c mvtnorm.cpp -o mvtnorm.o

mvt.o: mvt.f
	${FC} $(FLAG) -c mvt.f

randomF77.o: randomF77.c
	${CC} $(FLAG) -c randomF77.c

clean:
	rm -f *.o *.a
