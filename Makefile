SHELL := /bin/bash

# simple makefile for development.
SRC=msprime/_msprimemodule.c

# The default target builds the C module in the simplest way.
# NOTE: delete the "uv run" from the start if you're not using uv
cmodule: ${SRC}
	uv run python setup.py build_ext --inplace

# allchecks turns on as many checks as make sense when building
# Python-C extensions.
allchecks: ${SRC}
	CFLAGS="-std=c99 -Wall -Wextra -Werror -Wno-unused-parameter" && \
	CFLAGS+=" -Wno-missing-field-initializers -Wno-cast-function-type" && \
	CFLAGS+=" --coverage" && \
	export CFLAGS && uv run python3 setup.py build_ext --inplace

# Turn on coverage builds
coverage: ${SRC}
	rm -fR build
	CFLAGS="-coverage" uv run python3 setup.py build_ext --inplace

# Format the C code ready for a PR
clang-format:
	clang-format -i lib/tests/* lib/*.[c,h]

tags:
	ctags -f TAGS msprime/*.c lib/*.[c,h] msprime/*.py tests/*.py algorithms.py


clean:
	rm -fR build
	rm -f msprime/*.so
