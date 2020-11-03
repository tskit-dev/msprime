
.. _sec_likelihood:

=====================
Computing likelihoods
=====================

.. todo:: Review the content in this page and decide if the structure is
    correct. Arguably this could be moved to the "utilities" page, but
    this does seem like a fundamental feature of msprime like simulating
    ancestry and mutations and not a "utility".

``msprime`` provides the capability to evaluate the sampling probabilities:
that of a stored tree sequence for a given diploid effective population size
:math:`N_e` and per-link, per-generation recombination probability :math:`r`
under the standard ancestral recombination graph; and that of a pattern of
mutations given a tree sequence and per-site, per-generation mutation
probability :math:`\mu` under the infinite sites model.

***
API
***

.. autofunction:: msprime.log_arg_likelihood

.. autofunction:: msprime.unnormalised_log_mutation_likelihood

********
Examples
********

.. todo:: Refactor this to use the jupyter execute and to follow
    similar conventions as the ancestry.rst page.


``msprime`` can be used to evaluate the sampling probability of a tree sequence
for a given effective population size and per-site, per-generation recombination
rate, as well as the probability of a configuration of infinite sites mutations
given a tree sequence and a per-site, per-generation mutation probability. In
both cases, the tree sequence must conform to the ``record_full_arg`` option of
the :func:`.simulate` function. The following example illustrates the evaluation
of these log likelihoods:

.. code-block:: python

    def likelihood_example():
        ts = msprime.simulate(
            sample_size=5, recombination_rate=0.1, mutation_rate=0.1,
            record_full_arg=True, random_seed=42)
        print(msprime.log_arg_likelihood(ts, recombination_rate=0.1, Ne=1))
        print(msprime.log_arg_likelihood(ts, recombination_rate=1, Ne=1))
        print(msprime.log_arg_likelihood(ts, recombination_rate=1, Ne=10))
        print(msprime.unnormalised_log_mutation_likelihood(ts, mu=0))
        print(msprime.unnormalised_log_mutation_likelihood(ts, mu=0.1))
        print(msprime.unnormalised_log_mutation_likelihood(ts, mu=1))

Running this code we get::

    -11.22279995534112
    -14.947399100839986
    -22.154011926066893
    -inf
    -5.665181028073889
    -7.087195080578711
