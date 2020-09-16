all: scrm ms msHOT slim seq-gen msms

seq-gen:
	rm -fR seq-gen Seq-Gen* 1.3.4.tar.gz
	wget https://github.com/rambaut/Seq-Gen/archive/1.3.4.tar.gz
	tar -zxf 1.3.4.tar.gz
	make -C Seq-Gen-1.3.4/source
	cp Seq-Gen-1.3.4/source/seq-gen ./

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

msHOT:
	make -C msHOTdir
	cp msHOTdir/msHOT ./
	cp msHOTdir/msHOT_summary_stats ./

discoal:
	wget https://github.com/kern-lab/discoal/archive/v0.1.5.tar.gz
	tar -xf v0.1.5.tar.gz
	# To allow larger num_sites, a macro in the header file is updated
	sed -i -e 's/^#define MAXSITES.*$$/#define MAXSITES 500000/' \
		discoal-0.1.5/discoal.h 
	cd discoal-0.1.5/ && make
	cp discoal-0.1.5/discoal ./

msms:
	wget https://www.mabs.at/ewing/msms/msms3.2rc-b163.jar
	mv msms3.2rc-b163.jar msms.jar

msprime-0.7.4:
	python3 -m venv msprime-0.7.4
	bash -c "source msprime-0.7.4/bin/activate; pip install msprime==0.7.4"

clean:
	rm -fR scrm* ms ms_summary_stats msHOT msHOT_summary_stats sample_stats \
		discoal* v0.1.4.tar.gz* seq-gen SeqGen* SLiM msms.jar msprime*
