
.. _sec_utilities:

=========
Utilities
=========

*********
Rate maps
*********

.. todo:: Add some high-level content here.

.. autoclass:: msprime.RateMap
    :members:

.. autofunction:: msprime.read_hapmap

--------------
Deprecated API
--------------

.. autoclass:: msprime.RecombinationMap
    :members:

*******
Logging
*******

Msprime uses the Python :mod:`logging` infrastructure to help debugging
complex simulations. Messages at the INFO level are high-level information
about the state of the current simulation, such as switches in simulation
model, etc. While it is straightforward to set up logging messages using
the built-in Python methods, the `daiquiri
<https://daiquiri.readthedocs.io/en/latest/>`_ is a bit more convenient.
For example,

.. literalinclude:: examples/logging_info.py

Running this code we would get output something like the following:

.. code-block:: none

    2020-06-24 16:02:20,983 [11188] INFO     msprime.ancestry: Running model {'name': 'dtwf'} until max time: 100.000000
    2020-06-24 16:02:20,984 [11188] INFO     msprime.ancestry: Running model {'name': 'hudson'} until max time: inf
    2020-06-24 16:02:20,984 [11188] INFO     msprime.ancestry: Completed at time=1449.4 nodes=19 edges=18


When running larger simulations and trying to figure out when
they might finish, it can be helpful to use the DEBUG logging output.
For example,

.. literalinclude:: examples/logging_debug.py

which gives us:

.. code-block:: none

    2020-06-24 16:11:33,373 [11729] INFO     msprime.ancestry: Running model {'name': 'hudson'} until max time: inf
    2020-06-24 16:11:33,396 [11729] DEBUG    msprime.ancestry: time=0.0444792 ancestors=90162
    2020-06-24 16:11:33,418 [11729] DEBUG    msprime.ancestry: time=0.0979881 ancestors=80314
    2020-06-24 16:11:33,438 [11729] DEBUG    msprime.ancestry: time=0.167461 ancestors=70518
    2020-06-24 16:11:33,466 [11729] DEBUG    msprime.ancestry: time=0.258912 ancestors=60734
    2020-06-24 16:11:33,501 [11729] DEBUG    msprime.ancestry: time=0.386571 ancestors=51002
    2020-06-24 16:11:33,542 [11729] DEBUG    msprime.ancestry: time=0.575475 ancestors=41292
    2020-06-24 16:11:33,596 [11729] DEBUG    msprime.ancestry: time=0.872121 ancestors=31704
    2020-06-24 16:11:33,668 [11729] DEBUG    msprime.ancestry: time=1.42071 ancestors=22228
    2020-06-24 16:11:33,763 [11729] DEBUG    msprime.ancestry: time=2.75437 ancestors=13042
    2020-06-24 16:11:33,906 [11729] DEBUG    msprime.ancestry: time=8.80469 ancestors=4752
    2020-06-24 16:11:34,018 [11729] DEBUG    msprime.ancestry: time=277.713 ancestors=416
    2020-06-24 16:11:34,034 [11729] DEBUG    msprime.ancestry: time=6352.95 ancestors=142
    2020-06-24 16:11:34,042 [11729] DEBUG    msprime.ancestry: time=20955.5 ancestors=104
    2020-06-24 16:11:34,050 [11729] DEBUG    msprime.ancestry: time=40321.8 ancestors=83
    2020-06-24 16:11:34,056 [11729] DEBUG    msprime.ancestry: time=72735.7 ancestors=74
    2020-06-24 16:11:34,059 [11729] INFO     msprime.ancestry: Completed at time=235063 nodes=207110 edges=233940


In this example we run a reasonably large simulation and turn on
the DEBUG output. This will then periodically (every 10,000 simulation
events) print out the current time in the simulation, and the
number of extant ancestral lineages.

.. warning:: The format of these logging messages is not fixed and may change
    arbitrarily in the future. If you need to obtain the information within
    them, please open an issue on GitHub so that we can provide a documented
    API for this.


