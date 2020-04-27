#
# Copyright (C) 2018-2020 University of Oxford
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
import inspect
import sys

import numpy as np
import tskit

import _msprime
from . import provenance
from . import utils


class MutationModel(_msprime.MutationModel):
    """
    Superclass of mutation models. Allows you to build your own mutation model.

    TODO document properly.
    """

    def asdict(self):
        # This version of asdict makes sure that we have sufficient parameters
        # to call the contructor and recreate the class. However, this means
        # that subclasses *must* have an instance variable of the same name.
        # This is essential for Provenance round-tripping to work.
        return {
            key: getattr(self, key)
            for key in inspect.signature(self.__init__).parameters.keys()
            if hasattr(self, key)
        }

    def __str__(self):
        alleles = " ".join([x.decode() for x in self.alleles])
        s = "Mutation model with alleles {}\n".format(alleles)
        s += "  root distribution: {}\n".format(
            " ".join(map(str, self.root_distribution))
        )
        s += "  transition matrix:\n"
        for row in self.transition_matrix:
            s += "     {}\n".format(" ".join(map(str, row)))
        return s


class BinaryMutations(MutationModel):
    """
    The simplest mutational model with 0/1 states.

    TODO document properly.
    """

    def __init__(self):
        alleles = [b"0", b"1"]
        root_distribution = [1, 0]
        transition_matrix = [[0, 1], [1, 0]]
        alleles = [b"0", b"1"]
        super().__init__(alleles, root_distribution, transition_matrix)


class JukesCantor(MutationModel):
    """
    The Jukes-Cantor mutation model.

    .. todo: documentation
    """

    def __init__(self):
        alleles = [b"A", b"C", b"T", b"G"]
        root_distribution = [0.25, 0.25, 0.25, 0.25]
        transition_matrix = np.zeros((4, 4))
        transition_matrix[:] = 1 / 3
        np.fill_diagonal(transition_matrix, 0)
        super().__init__(alleles, root_distribution, transition_matrix)


# Pre 1.0, we had these constants to define the alphabet for the
# simple infinite sites mutations. These should be deprecated and removed
# along with the InfiniteSites class.
BINARY = 0
NUCLEOTIDES = 1


class InfiniteSites(MutationModel):
    # This mutation model is defined for backwards compatability, and is a remnant
    # of an earlier design. The class should be formally deprecated and removed at
    # some point.
    def __init__(self, alphabet=BINARY):
        self.alphabet = alphabet
        models = {BINARY: BinaryMutations(), NUCLEOTIDES: JukesCantor()}
        if alphabet not in models:
            raise ValueError("Bad alphabet")
        model = models[alphabet]
        super().__init__(
            model.alleles, model.root_distribution, model.transition_matrix
        )


class MutationMap(object):
    # TODO Get rid of this class and use IntervalMap (or maybe IntervalRateMap) instead.
    # See
    # https://github.com/tskit-dev/msprime/issues/902
    # https://github.com/tskit-dev/msprime/issues/920
    def __init__(self, position, rate):
        self.position = position
        self.rate = rate
        self._ll_map = _msprime.IntervalMap(position=position, value=rate)

    def asdict(self):
        return {"position": self.position, "rate": self.rate}


def _simple_mutation_generator(rate, sequence_length, rng):
    """
    Factory function used to create a low-level mutation generator that
    produces binary infinite sites mutations suitable for use in
    msprime.simulate() and the mspms CLI.
    """
    # Add a ``discrete`` parameter and pass to the MutationGenerator
    # constructor.
    if rate is None:
        return None
    rate_map = MutationMap(position=[0, sequence_length], rate=[rate, 0])
    return _msprime.MutationGenerator(rng, rate_map._ll_map, BinaryMutations())


def mutate(
    tree_sequence,
    rate=None,
    random_seed=None,
    model=None,
    keep=False,
    start_time=None,
    end_time=None,
    discrete=False,
):
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
    :param float rate: The rate of mutation per generation, as either a
        single number (for a uniform rate) or as a
        :class:`.MutationMap`. (Default: 0).
    :param int random_seed: The random seed. If this is `None`, a
        random seed will be automatically generated. Valid random
        seeds must be between 1 and :math:`2^{32} - 1`.
    :param MutationModel model: The mutation model to use when generating
        mutations. If not specified or None, the :class:`.BinaryMutations`
        mutation model is used.
    :param bool keep: Whether to keep existing mutations (default: False).
    :param float start_time: The minimum time ago at which a mutation can
        occur. (Default: no restriction.)
    :param float end_time: The maximum time ago at which a mutation can occur
        (Default: no restriction).
    :param bool discrete: Whether to generate mutations at only integer positions
        along the genome.  Default is False, which produces infinite-sites
        mutations at floating-point positions.
    :return: The :class:`tskit.TreeSequence` object  resulting from overlaying
        mutations on the input tree sequence.
    :rtype: :class:`tskit.TreeSequence`
    """
    try:
        tables = tree_sequence.tables
    except AttributeError:
        raise ValueError("First argument must be a TreeSequence instance.")
    seed = random_seed
    if random_seed is None:
        seed = utils.get_random_seed()
    else:
        seed = int(seed)

    if rate is None:
        rate = 0
    try:
        rate = float(rate)
        rate_map = MutationMap(
            position=[0.0, tree_sequence.sequence_length], rate=[rate, 0.0]
        )
    except TypeError:
        rate_map = rate
    if not isinstance(rate_map, MutationMap):
        raise TypeError("rate must be a float or a MutationMap")

    if start_time is None:
        start_time = -sys.float_info.max
    else:
        start_time = float(start_time)
    if end_time is None:
        end_time = sys.float_info.max
    else:
        end_time = float(end_time)
    if start_time > end_time:
        raise ValueError("start_time must be <= end_time")
    keep = bool(keep)
    discrete = bool(discrete)

    if model is None:
        model = BinaryMutations()
    if not isinstance(model, MutationModel):
        raise TypeError("model must be a MutationModel")

    argspec = inspect.getargvalues(inspect.currentframe())
    parameters = {
        "command": "mutate",
        **{arg: argspec.locals[arg] for arg in argspec.args},
    }
    parameters["random_seed"] = seed
    encoded_provenance = provenance.json_encode_provenance(
        provenance.get_provenance_dict(parameters)
    )

    rng = _msprime.RandomGenerator(seed)
    mutation_generator = _msprime.MutationGenerator(
        random_generator=rng, rate_map=rate_map._ll_map, model=model
    )
    lwt = _msprime.LightweightTableCollection()
    lwt.fromdict(tables.asdict())
    mutation_generator.generate(
        lwt, keep=keep, start_time=start_time, end_time=end_time, discrete=discrete
    )

    tables = tskit.TableCollection.fromdict(lwt.asdict())
    tables.provenances.add_row(encoded_provenance)
    return tables.tree_sequence()
