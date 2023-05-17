AR       ?= ar
CXX      ?= g++
CXXFLAGS  = -O2 -DNDEBUG -fPIC

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
 CXXFLAGS += -arch x86_64
endif

qf.a: qfc.o
	${AR} rcs $@ $^ 

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

clean:
	rm -f *.o *.a
