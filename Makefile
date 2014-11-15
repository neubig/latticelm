# CXX=g++
# CC=g++
CXXFLAGS=-g -O3 -Wall -std=c++0x
LDFLAGS=-g -O3 -lfst -ldl -std=c++0x

all: latticelm

latticelm: latticelm.h pylm.h lexfst.h ${ADDLD}
	${CXX} -o latticelm mainlatticelm.cc ${LDFLAGS} 

clean:
	rm -f latticelm
