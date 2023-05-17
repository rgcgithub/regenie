AR       ?= ar
CXX      ?= g++
CXXFLAGS  = -O2 -DNDEBUG -fPIC -arch x86_64

qf.a: qfc.o
	${AR} rcs $@ $^ 

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

clean:
	rm -f *.o *.a
