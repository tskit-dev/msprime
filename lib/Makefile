CC=h5cc 
CFLAGS=-std=c99 -pedantic -Werror -Wall -W \
  -Wno-unused-parameter\
  -Wmissing-prototypes -Wstrict-prototypes \
  -Wconversion -Wshadow -Wpointer-arith \
  -Wcast-qual -Wcast-align \
  -Wwrite-strings -Wnested-externs \
  -fshort-enums -fno-common -Dinline= -g -O2\
  -DH5_NO_DEPRECATED_SYMBOLS
LDFLAGS=-lgsl -lgslcblas -lhdf5 -lm

SRC=msprime.c fenwick.c fenwick.h msprime.h err.h tree_sequence.c \
    object_heap.c object_heap.h newick.c hapgen.c recomb_map.c mutgen.c \
    environment.c vargen.c vcf.c
COMPILED=msprime.c fenwick.c tree_sequence.c object_heap.c newick.c \
    hapgen.c recomb_map.c mutgen.c environment.c vargen.c vcf.c

all: main tests

avl.o: avl.c
	${CC} -Wall -g -O3 -c avl.c

main: main.c ${SRC} avl.o
	${CC} ${CFLAGS} -o main main.c ${COMPILED} ${LDFLAGS} -lconfig avl.o

tests: tests.c ${SRC} avl.o
	${CC} -coverage ${CFLAGS} -o tests tests.c ${COMPILED} ${LDFLAGS} -lcunit avl.o

tags:
	etags *.c *.h 

clean:
	rm -f main test_fenwick *.o *.gcda *.gcno

travis-tests: CC=gcc 
travis-tests: CFLAGS=-DH5_NO_DEPRECATED_SYMBOLS -coverage
travis-tests: avl.o tests 

coverity-tests: CC=gcc 
coverity-tests: CFLAGS+=-I/usr/include/hdf5/serial
coverity-tests: LDFLAGS+=-L/usr/lib/x86_64-linux-gnu/hdf5/serial
coverity-tests: avl.o tests main

coverity: 
	# We must use CC=gcc for this target or coverity won't 
	# work properly.
	make clean
	/home/jk/admin/software/cov-analysis-linux64-8.5.0/bin/cov-build --dir cov-int make coverity-tests
	tar -zcvf cov-int.tar.gz cov-int
