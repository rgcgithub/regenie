FLAG=-O2 -DNDEBUG -fPIC -std=legacy
AR   ?= ar
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

OBJECTS = $(patsubst %.f,%.o,$(wildcard ./*.f))

libquad.a: ${OBJECTS}
	${AR} rcs $@ $^ 

%.o: %.f
	${FC} $(FLAG) -c $<

clean:
	rm -f *.o *.a
