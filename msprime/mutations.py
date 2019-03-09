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
import json
import sys
import tskit

import _msprime
import msprime.simulations as simulations
import msprime.provenance as provenance


# Alphabets for mutations.
BINARY = 0
NUCLEOTIDES = 1


class InfiniteSites(object):
    """
    The "infinitely many sites" mutation model. In this model each mutation
    corresponds to a unique site, which has a floating-point position chosen
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


def mutate(
        tree_sequence, rate=None, random_seed=None, model=None, keep=False,
        start_time=None, end_time=None):
    """
    Simulates mutations on the specified ancestry and returns the resulting
    :class:`tskit.TreeSequence`. Mutations are generated at the specified rate in
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

    The time interval over which mutations can occur may be controlled
    using the ``start_time`` and ``end_time`` parameters. The ``start_time``
    defines the lower bound (in time-ago) on this interval and ``max_time``
    the upper bound. Note that we may have mutations associated with
    nodes with time <= ``start_time`` since mutations store the node at the
    bottom (i.e., towards the leaves) of the branch that they occur on.

    :param tskit.TreeSequence tree_sequence: The tree sequence onto which we
        wish to throw mutations.
    :param float rate: The rate of mutation per generation. (Default: 0).
    :param int random_seed: The random seed. If this is `None`, a
        random seed will be automatically generated. Valid random
        seeds must be between 1 and :math:`2^{32} - 1`.
    :param MutationModel model: The mutation model to use when generating
        mutations. If not specified or None, the :class:`.InfiniteSites`
        mutation model is used.
    :param bool keep: Whether to keep existing mutations (default: False).
    :param float start_time: The minimum time at which a mutation can
        occur. (Default: no restriction.)
    :param float end_time: The maximum time at which a mutation can occur
        (Default: no restriction).
    :return: The :class:`tskit.TreeSequence` object  resulting from overlaying
        mutations on the input tree sequence.
    :rtype: :class:`tskit.TreeSequence`
    """
    try:
        tables = tree_sequence.tables
    except AttributeError:
        raise ValueError("First argument must be a TreeSequence instance.")
    if random_seed is None:
        random_seed = simulations._get_random_seed()
    random_seed = int(random_seed)

    rng = _msprime.RandomGenerator(random_seed)
    if model is None:
        model = InfiniteSites()
    try:
        alphabet = model.alphabet
    except AttributeError:
        raise TypeError("model must be an InfiniteSites instance")
    if rate is None:
        rate = 0
    rate = float(rate)
    keep = bool(keep)

    parameters = {
        "command": "mutate", "rate": rate, "random_seed": random_seed, "keep": keep}

    if start_time is None:
        start_time = -sys.float_info.max
    else:
        start_time = float(start_time)
        parameters["start_time"] = start_time

    if end_time is None:
        end_time = sys.float_info.max
    else:
        end_time = float(end_time)
        parameters["end_time"] = end_time
    # TODO Add a JSON representation of the model to the provenance.
    provenance_dict = provenance.get_provenance_dict(parameters)

    if start_time > end_time:
        raise ValueError("start_time must be <= end_time")

    mutation_generator = _msprime.MutationGenerator(
        rng, rate, alphabet=alphabet, start_time=start_time, end_time=end_time)
    lwt = _msprime.LightweightTableCollection()
    lwt.fromdict(tables.asdict())
    mutation_generator.generate(lwt, keep=keep)

    tables = tskit.TableCollection.fromdict(lwt.asdict())
    tables.provenances.add_row(json.dumps(provenance_dict))
    return tables.tree_sequence()
