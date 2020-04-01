CC=gcc
EXTRA_CFLAGS=-std=c99 -pedantic -Werror -Wall -W \
  -Wmissing-prototypes -Wstrict-prototypes \
  -Wconversion -Wshadow -Wpointer-arith \
  -Wcast-align -Wcast-qual \
  -Wwrite-strings -Wnested-externs \
  -fshort-enums -fno-common -Dinline= 
CFLAGS=-g -O2 -Itskit/c -Itskit/c/kastore/c
LDFLAGS=-lgsl -lgslcblas -lm 

HEADERS=msprime.h util.h 
COMPILED=msprime.o fenwick.o object_heap.o \
    recomb_map.o mutgen.o avl.o util.o \
    tskit/c/tsk_core.o\
    tskit/c/tsk_tables.o\
    tskit/c/tsk_trees.o\
    tskit/c/kastore/c/kastore.o

all: main simulation_tests 

# We need a seperate rule for avl.c as it won't pass the strict checks.
avl.o: avl.c
	${CC} -Wall -g -O2 -c avl.c

argtable3.o: argtable3.c
	${CC} -Dlint -Wall -g -O2 -c argtable3.c

%.o : %.c ${HEADERS}
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

main: CFLAGS+=${EXTRA_CFLAGS}
main: main.c ${COMPILED} ${HEADERS} argtable3.o 
	${CC} ${CFLAGS} ${EXTRA_CFLAGS} -o main main.c ${COMPILED} argtable3.o ${LDFLAGS} -lconfig 

simulation_tests: simulation_tests.c ${COMPILED} ${HEADERS}
	${CC} ${CFLAGS} -Wall -o simulation_tests simulation_tests.c ${COMPILED} ${LDFLAGS} -lcunit 

tags:
	ctags *.c *.h

clean:
	rm -f main tests simulation_tests table_example table_example_cpp *.o *.gcda *.gcno

debian-tests: CFLAGS=--coverage 
debian-tests: simulation_tests

coverity-tests: CC=gcc 
coverity-tests: main simulation_tests

coverity: 
	# We must use CC=gcc for this target or coverity won't 
	# work properly.
	make clean
	/home/jk/admin/software/cov-analysis-linux64-8.5.0/bin/cov-build --dir cov-int make coverity-tests
	tar -zcvf cov-int.tar.gz cov-int
