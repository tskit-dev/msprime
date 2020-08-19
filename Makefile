# simple makefile for development.

SRC=msprime/_msprimemodule.c

ext3: ${SRC}
	python3 setup.py build_ext --inplace

allchecks: ${SRC}
	CFLAGS="-std=c99 -Wall -Wextra -Werror -Wno-unused-parameter -Wno-missing-field-initializers -Wno-cast-function-type" \
	python3 setup.py build_ext --inplace

ext3-coverage: ${SRC}
	rm -fR build
	CFLAGS="-coverage" python3 setup.py build_ext --inplace

docs: ext3
	cd docs && make clean && make html
	
tags:
	ctags -f TAGS msprime/*.c lib/*.[c,h] msprime/*.py tests/*.py

clean:
	rm -fR build
	rm -f msprime/*.o msprime/*.so
