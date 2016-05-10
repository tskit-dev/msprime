CC=h5cc
# CC=clang
# INCLUDES:=$(shell pkg-config gsl hdf5 libconfig --cflags)
CFLAGS=-std=c99 -pedantic -Werror -Wall -W \
  -Wno-unused-parameter\
  -Wmissing-prototypes -Wstrict-prototypes \
  -Wconversion -Wshadow -Wpointer-arith \
  -Wcast-qual -Wcast-align \
  -Wwrite-strings -Wnested-externs \
  -fshort-enums -fno-common -Dinline= -g -O2\
  -DH5_NO_DEPRECATED_SYMBOLS\
  ${INCLUDES}
# LDFLAGS:=$(shell pkg-config gsl hdf5 libconfig --libs)
LDFLAGS=-lgsl -lgslcblas -lhdf5 -lconfig


SRC=msprime.c fenwick.c fenwick.h msprime.h err.h main.c tree_sequence.c \
    object_heap.c object_heap.h newick.c hapgen.c recomb_map.c mutgen.c \
    environment.c vargen.c
COMPILED=msprime.c fenwick.c main.c tree_sequence.c object_heap.c newick.c \
    hapgen.c recomb_map.c mutgen.c environment.c vargen.c

all: main test_fenwick

avl.o: avl.c
	${CC} -Wall -g -O3 -c avl.c

main: ${SRC} avl.o
	${CC} ${CFLAGS} -o main ${COMPILED} ${LDFLAGS} avl.o
	
test_fenwick: fenwick.c fenwick.h test_fenwick.c
	${CC} ${CFLAGS} -o test_fenwick test_fenwick.c fenwick.c ${LDFLAGS}
	
tags:
	etags *.c *.h 

clean:
	rm -f main test_fenwick avl.o
