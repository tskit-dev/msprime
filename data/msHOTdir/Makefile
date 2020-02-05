CC=gcc
CFLAGS=-ggdb -Wno-unused-result -Wno-return-type
SRC=msGCHOT.c streecGCHOT.c rand1.c
all: msHOT msHOT_summary_stats

msHOT_summary_stats: ${SRC}
	${CC} ${CFLAGS} -DSUPPRESS_OUTPUT -DSUMMARY_STATS -o msHOT_summary_stats ${SRC} -lm

msHOT: ${SRC}
	${CC} ${CFLAGS} -o msHOT ${SRC} -lm

clean:
	rm msHOT msHOT_summary_stats
