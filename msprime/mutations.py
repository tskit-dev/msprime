#
# Copyright (C) 2018 University of Oxford
#
# This file is part of msprime.
#
# msprime is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# msprime is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with msprime.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Module responsible for generating mutations on a given tree sequence.
"""
from __future__ import division
from __future__ import print_function

import json

import _msprime
import msprime.simulations as simulations
import msprime.provenance as provenance


# Alphabets for mutations.
BINARY = 0
NUCLEOTIDES = 1


class InfiniteSites(object):
    """
    The "infinitely many sites" mutation model. In this model each mutation
    corresponds to a unique site, which has a floating-point position chose
    uniformly along the sequence. As a result, each site is associated with
    exactly one mutation.

    By default, the ancestral and derived states in this model are
    represented by the characters "0" and "1". Thus, the ancestral state
    at a site is always "0" and the derived state for a mutation is
    always "1". However, by specifying the ``alphabet=NUCLEOTIDES`` we
    can generate mutations from the nucleotide characters ACGT. In this
    case, for each mutation an ancestral state is chosen uniformly from
    these letters. The derived state is then chosen uniformly from the
    *remaining* characters so that the ancestral and derived states
    are always distinct.
    """
    def __init__(self, alphabet=BINARY):
        if alphabet not in [BINARY, NUCLEOTIDES]:
            raise ValueError("Bad alphabet")
        self.alphabet = alphabet


def mutate(tree_sequence, rate=None, random_seed=None, model=None, keep=False):
    """
    Simulates mutations on the specified ancestry and returns the resulting
    :class:`.TreeSequence`. Mutations are generated at the specified rate in
    measured generations. Mutations are generated under the infinite sites
    model, and so the rate of new mutations is per unit of sequence length per
    generation.

    If a random seed is specified, this is used to seed the random number
    generator. If the same seed is specified and all other parameters are equal
    then the same mutations will be generated. If no random seed is specified
    then one is generated automatically.

    If the ``model`` parameter is specified, this determines the model under
    which mutations are generated. Currently only the :class:`.InfiniteSites`
    mutation model is supported. This parameter is useful if you wish to obtain
    sequences with letters from the nucleotide alphabet rather than the default
    0/1 states. By default mutations from the infinite sites model with a binary
    alphabet are generated.

    By default, sites and mutations in the parameter tree sequence are
    discarded. If the ``keep`` parameter is true, however, *additional*
    mutations are simulated. Under the infinite sites mutation model, all new
    mutations generated will occur at distinct positions from each other and
    from any existing mutations (by rejection sampling).

    :param TreeSequence tree_sequence: The tree sequence onto which we
        wish to throw mutations.
    :param float rate: The rate of mutation per generation.
    :param int random_seed: The random seed. If this is `None`, a
        random seed will be automatically generated. Valid random
        seeds must be between 1 and :math:`2^{32} - 1`.
    :param MutationModel model: The mutation model to use when generating
        mutations. If not specified or None, the :class:`.InfiniteSites`
        mutation model is used.
    :return: The :class:`.TreeSequence` object  resulting from overlaying
        mutations on the input tree sequence.
    :rtype: :class:`.TreeSequence`
    """
    try:
        tables = tree_sequence.dump_tables()
    except AttributeError:
        raise ValueError("First argument must be a TreeSequence instance.")
    if random_seed is None:
        random_seed = simulations._get_random_seed()
    rng = _msprime.RandomGenerator(int(random_seed))
    if model is None:
        model = InfiniteSites()
    try:
        alphabet = model.alphabet
    except AttributeError:
        raise TypeError("model must be an InfiniteSites instance")
    mutation_generator = _msprime.MutationGenerator(rng, rate, alphabet=alphabet)
    mutation_generator.generate(tables.ll_tables, keep=keep)
    parameters = {"rate": rate, "random_seed": random_seed}
    provenance_dict = provenance.get_provenance_dict("mutate", parameters)
    tables.provenances.add_row(json.dumps(provenance_dict))
    return tables.tree_sequence()
