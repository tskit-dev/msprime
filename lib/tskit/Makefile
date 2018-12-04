CC=gcc
CFLAGS=-g -O2 -std=c99 -pedantic -Werror -Wall -W \
  -Wmissing-prototypes -Wstrict-prototypes \
  -Wconversion -Wshadow -Wpointer-arith \
  -Wcast-align -Wcast-qual \
  -Wwrite-strings -Wnested-externs \
  -fshort-enums -fno-common -Dinline= \
  -I../kastore/c
LDFLAGS=-lm
OBJECTS=tsk_core.o tsk_tables.o tsk_trees.o tsk_genotypes.o \
	tsk_convert.o tsk_stats.o kastore.o 
EXECUTABLES=main test_core test_tables test_trees test_genotypes

all: old_tests ${EXECUTABLES} libtskit.a

main: main.c ${OBJECTS} argtable3.o
	${CC} ${CFLAGS} -o $@ $^ ${LDFLAGS} 

libtskit.a: ${OBJECTS}
	ar rcs $@ $^

test_%: test_%.c testlib.o ${OBJECTS}
	${CC} ${CFLAGS} -o $@ $^ ${LDFLAGS} -lcunit

# We can't turn on all the usual compiler checks because of CUnit.
testlib.o: testlib.c
	${CC} -I../kastore/c -std=c99 -Wall -g -c $^ 

argtable3.o: argtable3.c
	${CC} -Dlint -Wall -g -O2 -c argtable3.c

old_tests: old_tests.c ${OBJECTS}
	${CC} -I../kastore/c -std=c99 -Wall -g $^ ${LDFLAGS} -lcunit -lgsl -lgslcblas

kastore.o: ../kastore/c/kastore.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< 

clean:
	rm -f *.o ${EXECUTABLES}

tags:
	ctags *.c *.h
