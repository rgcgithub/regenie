FLAG=-O2 -DNDEBUG
libMvtnorm.a: mvtnorm.o mvt.o  randomF77.o
	ar rcs $@ $^ 

mvtnorm.o: mvtnorm.cpp mvtnorm.h
	g++ $(FLAG) -c mvtnorm.cpp -o mvtnorm.o

randomF77.o: randomF77.c
	gcc $(FLAG) -c randomF77.c

mvt.o: mvt.f
	gfortran $(FLAG) -c mvt.f

clean:
	rm -f *.o *.a
