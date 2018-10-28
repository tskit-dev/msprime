all: scrm ms slim

slim:
	rm -fR SLiM*
	wget http://benhaller.com/slim/SLiM.zip
	unzip SLiM.zip
	mkdir SLiM/build
	# Arbitrarily using 4 threads to build here..
	cd SLiM/build && cmake .. && make -j 4
	cp SLiM/build/slim ./

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
