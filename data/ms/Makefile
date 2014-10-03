CFLAGS=-ggdb -Wno-unused-result
SRC=ms.c streec.c rand1.c
all: ms ms_summary_stats 

ms_summary_stats: ${SRC} 
	gcc ${CFLAGS} -DSUPPRESS_OUTPUT -DSUMMARY_STATS -o ms_summary_stats ${SRC} -lm

ms: ${SRC} 
	gcc ${CFLAGS} -o ms ${SRC} -lm

clean:
	rm -f ms ms_summary_stats
