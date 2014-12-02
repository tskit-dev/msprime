CFLAGS=-Wall -g -O3
#CFLAGS=-Wall -g -O3 -march=native 
LDFLAGS=-lgsl -lgslcblas -lconfig
SRC=msprime.c fenwick.c fenwick.h msprime.h err.h main.c tree_file.c hapgen.c
COMPILED=msprime.c fenwick.c main.c avl.c tree_file.c hapgen.c

all: main test_fenwick

main: ${SRC}
	gcc ${CFLAGS} -o main ${COMPILED} ${LDFLAGS}
	
test_fenwick: fenwick.c fenwick.h test_fenwick.c
	gcc ${CFLAGS} -o test_fenwick test_fenwick.c fenwick.c ${LDFLAGS}
	
tags:
	etags *.c *.h 

clean:
	rm -f main test_fenwick
