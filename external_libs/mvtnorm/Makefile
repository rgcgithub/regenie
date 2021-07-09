FLAG=-O2 -DNDEBUG -fPIC

libMvtnorm.a: mvtnorm.o mvt.o randomF77.o
	ar rcs $@ $^ 

mvtnorm.o: mvtnorm.cpp mvtnorm.h
	g++ $(FLAG) -c mvtnorm.cpp -o mvtnorm.o

mvt.o: mvt.f
	gfortran $(FLAG) -c mvt.f

randomF77.o: randomF77.c
	gcc $(FLAG) -c randomF77.c

clean:
	rm -f *.o *.a
