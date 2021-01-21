
all:
	# Need to set PYTHONPATH so that we pick up the local msprime
	# TODO add -W and -n options to make build stricter.
	PYTHONPATH=${PWD}/.. jupyter-book build -W --keep-going ./

clean:
	rm -fR _build
