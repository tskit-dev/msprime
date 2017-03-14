all: scrm ms 

scrm:
	wget https://github.com/scrm/scrm/releases/download/v1.7.2/scrm-src.tar.gz
	tar -xf scrm-src.tar.gz
	cd scrm-1.7.2 && ./configure && make
	cp scrm-1.7.2/scrm ./

ms:	
	make -C msdir
	cp msdir/ms ./
	cp msdir/ms_summary_stats ./
	cp msdir/sample_stats ./

clean:
	rm -fR scrm* ms ms_summary_stats sample_stats
