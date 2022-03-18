AR   		?= ar
CXX     ?= g++
CXXFLAGS = -O2 -DNDEBUG -fPIC

qf.a: qfc.o
	${AR} rcs $@ $^ 

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

clean:
	rm -f *.o *.a
