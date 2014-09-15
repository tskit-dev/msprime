CFLAGS=-Wall -g -O3
#CFLAGS=-Wall -g -O3 -march=native 
LDFLAGS=-lgsl -lgslcblas -lconfig -lavl
SRC=msprime.c util.c util.h fenwick.c fenwick.h msprime.h
COMPILED=msprime.c util.c fenwick.c main.c

all: main test_fenwick

main: ${SRC}
	gcc ${CFLAGS} -o main ${COMPILED} ${LDFLAGS}
	
test_fenwick: fenwick.c fenwick.h test_fenwick.c
	gcc ${CFLAGS} -o test_fenwick test_fenwick.c fenwick.c util.c ${LDFLAGS}
	
tags:
	etags *.c *.h *.py

clean:
	rm -f main test_fenwick
