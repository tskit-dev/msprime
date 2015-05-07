CFLAGS=-Wall -g -O3
#CFLAGS=-Wall -g -O3 -march=native 
CFLAGS=-std=c99 -pedantic -Werror -Wall -W \
  -Wno-unused-parameter\
  -Wmissing-prototypes -Wstrict-prototypes \
  -Wconversion -Wshadow -Wpointer-arith \
  -Wcast-qual -Wcast-align \
  -Wwrite-strings -Wnested-externs \
  -fshort-enums -fno-common -Dinline= -g -O2\
  -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial
# TODO fix the paths here.
LDFLAGS=-lgsl -lgslcblas -lconfig -lm -lhdf5
SRC=msprime.c fenwick.c fenwick.h msprime.h err.h main.c tree_sequence.c \
    object_heap.c object_heap.h
COMPILED=msprime.c fenwick.c main.c tree_sequence.c object_heap.c

all: main test_fenwick

avl.o: avl.c
	gcc -Wall -g -O3 -c avl.c

main: ${SRC} avl.o
	gcc ${CFLAGS} -o main ${COMPILED} ${LDFLAGS} avl.o
	
test_fenwick: fenwick.c fenwick.h test_fenwick.c
	gcc ${CFLAGS} -o test_fenwick test_fenwick.c fenwick.c ${LDFLAGS}
	
tags:
	etags *.c *.h 

clean:
	rm -f main test_fenwick avl.o
