SHELL := /bin/bash

# simple makefile for development.
SRC=msprime/_msprimemodule.c

# The default target builds the C module in the simplest way.
cmodule: ${SRC}
	python3 setup.py build_ext --inplace

test:
	gcc \
        -undefined dynamic_lookup \
	-Wl,-rpath,${CONDA_ENV}/lib \
        -Wno-unused-result \
	-Wsign-compare -Wunreachable-code -Wall -Wstrict-prototypes \
        -DNDEBUG -g -fwrapv -O3 -arch x86_64 \
        -Ilib \
        -Igit-submodules/tskit/c \
        -Igit-submodules/tskit/c/tskit \
        -Igit-submodules/tskit/c/subprojects/kastore \
        -I${CONDA_ENV}/include \
        -I${CONDA_ENV}/lib \
        -L${CONDA_ENV}/include \
        -L${CONDA_ENV}/lib \
        -lgsl -lgslcblas \
        lib/forward.c \
        lib/tests/test_forward.c \
        -o test_forward.out \
        -std=c99
		

# allchecks turns on as many checks as make sense when building
# Python-C extensions.
allchecks: ${SRC}
	CFLAGS="-std=c99 -Wall -Wextra -Werror -Wno-unused-parameter" && \
	CFLAGS+=" -Wno-missing-field-initializers -Wno-cast-function-type" && \
	CFLAGS+=" --coverage" && \
	export CFLAGS && python3 setup.py build_ext --inplace

# Turn on coverage builds
coverage: ${SRC}
	rm -fR build
	CFLAGS="-coverage" python3 setup.py build_ext --inplace

# Format the C code ready for a PR
clang-format:
	clang-format -i lib/tests/* lib/*.[c,h]

tags:
	ctags -f TAGS msprime/*.c lib/*.[c,h] msprime/*.py tests/*.py


clean:
	rm -fR build
	rm -f msprime/*.so
