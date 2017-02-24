=======
msprime
=======

Msprime is a reimplementation of Hudson's classical ms program for modern
datasets. It provides a convenient Python API and HDF5-based storage format, as
well as the ``mspms`` command line interface. This program provides a fully
``ms`` compatible interface, and can be used as a drop-in replacement in
existing workflows.

Msprime can simulate the coalescent with recombination much
faster than programs based on the Sequentially Markov Coalescent
for large sample sizes and has very reasonable memory requirements. Msprime
makes it possible to simulate chromosome sized regions with hundreds of
thousands of samples.

If you use ``msprime`` in your work, please cite the `PLOS Computational
Biology <http://dx.doi.org/10.1371/journal.pcbi.1004842>`_ paper.
See `here <https://msprime.readthedocs.org/en/stable/CITATION.html>`_ for
full citation details.

Please see the `documentation <https://msprime.readthedocs.org/en/stable/>`_
for further details.

Msprime is very portable, and provides a number of installation options.
See `here <https://msprime.readthedocs.org/en/stable/installation.html>`_ for
details.
