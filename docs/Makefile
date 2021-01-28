
# Need to set PYTHONPATH so that we pick up the local msprime
PYPATH=${PWD}/..
MSP_VERSION:=$(shell PYTHONPATH=${PYPATH} \
   python3 -c 'import msprime; print(msprime.__version__.split("+")[0])')

dev:
	# TODO add -n option to make build stricter.
	PYTHONPATH=${PYPATH} jupyter-book build -W --keep-going ./

dist:
	@echo Building distribution for msprime version ${MSP_VERSION}
	sed -i s/__MSPRIME_VERSION__/${MSP_VERSION}/g _config.yml
	# TODO add -n option to make build stricter.
	PYTHONPATH=${PYPATH} jupyter-book build -W --keep-going ./

clean:
	rm -fR _build
