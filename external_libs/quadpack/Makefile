FLAG=-O2 -DNDEBUG -fPIC -std=legacy
AR   ?= ar
FC   ?= gfortran

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
 FLAG += -arch x86_64
endif

OBJECTS = $(patsubst %.f,%.o,$(wildcard ./*.f))

libquad.a: ${OBJECTS}
	${AR} rcs $@ $^ 

%.o: %.f
	${FC} $(FLAG) -c $<

clean:
	rm -f *.o *.a
