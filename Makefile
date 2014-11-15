# CXX=g++
# CC=g++
CXXFLAGS=-g -O3 -Wall
LDFLAGS=-g -O3 -lfst -ldl

all: latticelm

latticelm: latticelm.h pylm.h lexfst.h ${ADDLD}
	${CXX} -o latticelm mainlatticelm.cc ${LDFLAGS} 

clean:
	rm -f latticelm
