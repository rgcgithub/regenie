FLAG=-O2 -DNDEBUG -fPIC -std=legacy -arch x86_64
AR   ?= ar
FC   ?= gfortran

OBJECTS = $(patsubst %.f,%.o,$(wildcard ./*.f))

libquad.a: ${OBJECTS}
	${AR} rcs $@ $^ 

%.o: %.f
	${FC} $(FLAG) -c $<

clean:
	rm -f *.o *.a
