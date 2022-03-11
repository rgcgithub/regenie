FLAG=-O2 -DNDEBUG -fPIC -std=legacy

OBJECTS = $(patsubst %.f,%.o,$(wildcard ./*.f))

libquad.a: ${OBJECTS}
	ar rcs $@ $^ 

%.o: %.f
	gfortran $(FLAG) -c $<

clean:
	rm -f *.o *.a
