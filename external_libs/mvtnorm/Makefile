FLAG  =-O2 -DNDEBUG -fPIC
AR   ?= ar
CXX  ?= g++
CC   ?= gcc
FC   ?= gfortran

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
