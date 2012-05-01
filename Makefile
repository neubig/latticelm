CXX=g++
CC=g++
CXXFLAGS=-g -O3 -Wall
LDFLAGS=-g -O3 -ldl -lfst

all: latticelm

latticelm: latticelm.h pylm.h lexfst.h ${ADDLD}
	${CXX} ${LDFLAGS} -o latticelm mainlatticelm.cc

clean:
	rm -f latticelm
