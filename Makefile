# CXX=g++
# CC=g++
FSTPATH=/Users/neubig/usr
LDFLAGS=-g -O3 -lfst -ldl -std=c++0x -I${FSTPATH}/include -L${FSTPATH}/lib

all: latticelm

latticelm: latticelm.h pylm.h lexfst.h ${ADDLD}
	${CXX} -o latticelm mainlatticelm.cc ${LDFLAGS} 

clean:
	rm -f latticelm
