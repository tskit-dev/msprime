CFLAGS=-Wall -g -O3
#CFLAGS=-Wall -g -O3 -march=native 
LDFLAGS=-lgsl -lgslcblas -lconfig -lavl
SRC=msprime.c util.c util.h bit.c bit.h msprime.h
COMPILED=msprime.c util.c bit.c main.c

all: main test_bit

main: ${SRC}
	gcc ${CFLAGS} -o main ${COMPILED} ${LDFLAGS}
	
test_bit: bit.c bit.h test_bit.c
	gcc ${CFLAGS} -o test_bit test_bit.c bit.c util.c ${LDFLAGS}
	
tags:
	etags *.c *.h *.py

clean:
	rm -f main test_bit
