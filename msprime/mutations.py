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
from . import core
from . import provenance
from _msprime import BaseMutationModel

_ACGT_ALLELES = ["A", "C", "G", "T"]
_AMINO_ACIDS = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
]


# NOTE we're doing this hack of monkey-patching asdict onto the
# MutationModel subclasses because of a bug in sphinx, where it's
# not possible to have muliple inheritance with mocked out classes.
def safe_asdict(self):
    # This version of asdict makes sure that we have sufficient parameters
    # to call the contructor and recreate the class. However, this means
    # that subclasses *must* have an instance variable of the same name.
    # This is essential for Provenance round-tripping to work.
    return {
        key: getattr(self, key)
        for key in inspect.signature(self.__init__).parameters.keys()
        if hasattr(self, key)
    }


# TODO Change this to MatrixMutationModel
class MutationModel(_msprime.MatrixMutationModel):
    """
    Superclass of mutation models. Allows you to build your own mutation model.

    TODO document properly.
    """

    asdict = safe_asdict

    def __str__(self):
        s = f"Mutation model with alleles {self.alleles}\n"
        s += "  root distribution: {}\n".format(
            " ".join(map(str, self.root_distribution))
        )
        s += "  transition matrix:\n"
        for row in self.transition_matrix:
            s += "     {}\n".format(" ".join(map(str, row)))
        return s


class SlimMutationModel(_msprime.SlimMutationModel):
    asdict = safe_asdict


class InfiniteAllelesMutationModel(_msprime.InfiniteAllelesMutationModel):
    asdict = safe_asdict


class BinaryMutations(MutationModel):
    """
    The simplest mutational model with 0/1 states.

    TODO document properly.
    """

    def __init__(self):
        alleles = ["0", "1"]
        root_distribution = [1, 0]
        transition_matrix = [[0, 1], [1, 0]]
        super().__init__(alleles, root_distribution, transition_matrix)


class JukesCantor(MutationModel):
    """
    The Jukes-Cantor mutation model.

    .. todo: documentation
    """

    def __init__(self):
        alleles = _ACGT_ALLELES
        root_distribution = [0.25, 0.25, 0.25, 0.25]
        transition_matrix = np.zeros((4, 4))
        transition_matrix[:] = 1 / 3
        np.fill_diagonal(transition_matrix, 0)
        super().__init__(alleles, root_distribution, transition_matrix)


class HKY(MutationModel):
    """
    The Hasegawa, Kishino and Yano mutation model (Hasegawa et al. 1985).

    .. todo: documentation
    """

    def __init__(self, kappa=1.0, equilibrium_frequencies=None, root_distribution=None):
        alleles = _ACGT_ALLELES
        if equilibrium_frequencies is None:
            equilibrium_frequencies = np.array(
                [0.25, 0.25, 0.25, 0.25], dtype="float64"
            )
        if root_distribution is None:
            root_distribution = equilibrium_frequencies.copy()

        transition_matrix = np.full((4, 4), 1.0)
        np.fill_diagonal(transition_matrix, 0.0)
        # positions in transition matrix for transversions
        transition_pos = ((0, 1, 2, 3), (2, 3, 0, 1))
        transition_matrix[transition_pos] = kappa
        transition_matrix *= equilibrium_frequencies
        row_sums = transition_matrix.sum(axis=1, dtype="float64")
        transition_matrix /= max(row_sums)
        np.fill_diagonal(transition_matrix, 1.0 - row_sums / max(row_sums))

        super().__init__(alleles, root_distribution, transition_matrix)


class F84(MutationModel):
    """
    The F84 mutation model (Felsenstein and Churchill, 1996).
    2 types of events
    I: no change/transition
    II: no change/transition/transversion
    .. todo: documentation
    """

    def __init__(self, kappa=1.0, equilibrium_frequencies=None, root_distribution=None):
        alleles = _ACGT_ALLELES
        if equilibrium_frequencies is None:
            equilibrium_frequencies = [0.25, 0.25, 0.25, 0.25]
        if root_distribution is None:
            root_distribution = equilibrium_frequencies.copy()

        transition_matrix = np.full((4, 4), 1.0, dtype="float64")
        np.fill_diagonal(transition_matrix, 0.0)
        # positions in transition matrix for transversions
        transition_pos_AG = ((0, 2), (2, 0))
        transition_pos_CT = ((1, 3), (3, 1))
        p_AG = equilibrium_frequencies[0] + equilibrium_frequencies[2]
        p_CT = equilibrium_frequencies[1] + equilibrium_frequencies[3]
        gamma = kappa - 1
        transition_matrix[transition_pos_AG] = np.float64(1 + gamma / p_AG)
        transition_matrix[transition_pos_CT] = np.float64(1 + gamma / p_CT)
        transition_matrix *= equilibrium_frequencies
        row_sums = transition_matrix.sum(axis=1)
        transition_matrix = transition_matrix / max(row_sums)
        row_sums = transition_matrix.sum(axis=1, dtype="float64")
        np.fill_diagonal(transition_matrix, 1.0 - row_sums)

        super().__init__(alleles, root_distribution, transition_matrix)


class GTR(MutationModel):
    """
    The  Generalised time-reversible mutation model (Tavaré et al. 1986).
    .. todo: documentation
    """

    def __init__(
        self, relative_rates, equilibrium_frequencies=None, root_distribution=None
    ):
        alleles = _ACGT_ALLELES
        assert len(relative_rates) == 6
        if equilibrium_frequencies is None:
            equilibrium_frequencies = [0.25, 0.25, 0.25, 0.25]
        if root_distribution is None:
            root_distribution = equilibrium_frequencies

        transition_matrix = np.zeros((4, 4))
        # relative_rates: [A->C, A->G,A->T,C->G,C->T,G->T]
        tri_upper = np.triu_indices_from(transition_matrix, k=1)
        transition_matrix[tri_upper] = relative_rates
        transition_matrix += transition_matrix.T
        transition_matrix *= equilibrium_frequencies
        row_sums = transition_matrix.sum(axis=1)
        transition_matrix = transition_matrix / max(row_sums)
        row_sums = transition_matrix.sum(axis=1, dtype="float64")
        np.fill_diagonal(transition_matrix, 1.0 - row_sums)

        super().__init__(alleles, root_distribution, transition_matrix)


class BLOSUM62(MutationModel):
    """
    values copied from Seqgen: http://tree.bio.ed.ac.uk/software/seqgen/
    original paper: Henikoff, S., and J. G. Henikoff. 1992. PNAS USA 89:10915-10919.
    values included in: Yu,Y.-K., Wootton,J.C. and Altschul,S.F. (2003)
    The compositional adjustment of amino acid substitution matrices.
    Proc. Natl Acad. Sci., USA, 100, 15688–15693.
    """

    def __init__(self):
        alleles = _AMINO_ACIDS
        num_alleles = len(alleles)
        root_distribution = [
            0.074,
            0.052,
            0.045,
            0.054,
            0.025,
            0.034,
            0.054,
            0.074,
            0.026,
            0.068,
            0.099,
            0.058,
            0.025,
            0.047,
            0.039,
            0.057,
            0.051,
            0.013,
            0.032,
            0.073,
        ]
        relative_rates = [
            0.735790389698,
            0.485391055466,
            1.297446705134,
            0.543161820899,
            0.500964408555,
            3.180100048216,
            1.459995310470,
            0.227826574209,
            0.397358949897,
            0.240836614802,
            1.199705704602,
            3.020833610064,
            1.839216146992,
            1.190945703396,
            0.329801504630,
            1.170949042800,
            1.360574190420,
            1.240488508640,
            3.761625208368,
            0.140748891814,
            5.528919177928,
            1.955883574960,
            0.418763308518,
            1.355872344485,
            0.798473248968,
            0.418203192284,
            0.609846305383,
            0.423579992176,
            0.716241444998,
            1.456141166336,
            2.414501434208,
            0.778142664022,
            0.354058109831,
            2.435341131140,
            1.626891056982,
            0.539859124954,
            0.605899003687,
            0.232036445142,
            0.283017326278,
            0.418555732462,
            0.774894022794,
            0.236202451204,
            0.186848046932,
            0.189296292376,
            0.252718447885,
            0.800016530518,
            0.622711669692,
            0.211888159615,
            0.218131577594,
            0.831842640142,
            0.580737093181,
            0.372625175087,
            0.217721159236,
            0.348072209797,
            3.890963773304,
            1.295201266783,
            5.411115141489,
            1.593137043457,
            1.032447924952,
            0.285078800906,
            3.945277674515,
            2.802427151679,
            0.752042440303,
            1.022507035889,
            0.406193586642,
            0.445570274261,
            1.253758266664,
            0.983692987457,
            0.648441278787,
            0.222621897958,
            0.767688823480,
            2.494896077113,
            0.555415397470,
            0.459436173579,
            0.984311525359,
            3.364797763104,
            6.030559379572,
            1.073061184332,
            0.492964679748,
            0.371644693209,
            0.354861249223,
            0.281730694207,
            0.441337471187,
            0.144356959750,
            0.291409084165,
            0.368166464453,
            0.714533703928,
            1.517359325954,
            2.064839703237,
            0.266924750511,
            1.773855168830,
            1.173275900924,
            0.448133661718,
            0.494887043702,
            0.730628272998,
            0.356008498769,
            0.858570575674,
            0.926563934846,
            0.504086599527,
            0.527007339151,
            0.388355409206,
            0.374555687471,
            1.047383450722,
            0.454123625103,
            0.233597909629,
            4.325092687057,
            1.122783104210,
            2.904101656456,
            1.582754142065,
            1.197188415094,
            1.934870924596,
            1.769893238937,
            1.509326253224,
            1.117029762910,
            0.357544412460,
            0.352969184527,
            1.752165917819,
            0.918723415746,
            0.540027644824,
            1.169129577716,
            1.729178019485,
            0.914665954563,
            1.898173634533,
            0.934187509431,
            1.119831358516,
            1.277480294596,
            1.071097236007,
            0.641436011405,
            0.585407090225,
            1.179091197260,
            0.915259857694,
            1.303875200799,
            1.488548053722,
            0.488206118793,
            1.005451683149,
            5.151556292270,
            0.465839367725,
            0.426382310122,
            0.191482046247,
            0.145345046279,
            0.527664418872,
            0.758653808642,
            0.407635648938,
            0.508358924638,
            0.301248600780,
            0.341985787540,
            0.691474634600,
            0.332243040634,
            0.888101098152,
            2.074324893497,
            0.252214830027,
            0.387925622098,
            0.513128126891,
            0.718206697586,
            0.720517441216,
            0.538222519037,
            0.261422208965,
            0.470237733696,
            0.958989742850,
            0.596719300346,
            0.308055737035,
            4.218953969389,
            0.674617093228,
            0.811245856323,
            0.717993486900,
            0.951682162246,
            6.747260430801,
            0.369405319355,
            0.796751520761,
            0.801010243199,
            4.054419006558,
            2.187774522005,
            0.438388343772,
            0.312858797993,
            0.258129289418,
            1.116352478606,
            0.530785790125,
            0.524253846338,
            0.253340790190,
            0.201555971750,
            8.311839405458,
            2.231405688913,
            0.498138475304,
            2.575850755315,
            0.838119610178,
            0.496908410676,
            0.561925457442,
            2.253074051176,
            0.266508731426,
            1.000000000000,
        ]
        transition_matrix = np.zeros((num_alleles, num_alleles))
        tril = np.tril_indices(num_alleles, k=-1)
        transition_matrix[tril] = relative_rates
        transition_matrix += np.tril(transition_matrix).T
        transition_matrix *= root_distribution
        row_sums = transition_matrix.sum(axis=1)
        transition_matrix = transition_matrix / max(row_sums)
        row_sums = transition_matrix.sum(axis=1, dtype="float64")
        np.fill_diagonal(transition_matrix, 1.0 - row_sums)
        super().__init__(alleles, root_distribution, transition_matrix)


class PAM(MutationModel):
    """
    Dayhoff DCMut as described in Kosiol, C., and Goldman, N. 2005.
    Different versions of the Dayhoff rate matrix.
    Molecular Biology and Evolution 22:193-199.
    Values copied from PAML (http://abacus.gene.ucl.ac.uk/software/paml.html)
    original paper: Dayhoff, M.O., Schwartz, R.M., Orcutt, B.C. (1978).
    A model of evolutionary change in proteins. Atlas of Protein Sequence Structur.,
    Vol5, Suppl. 3, National Biomedical Research Foundation, Washington DC, pp. 345-352.
    """

    def __init__(self):
        alleles = _AMINO_ACIDS
        num_alleles = len(alleles)
        root_distribution = [
            0.087127,
            0.040904,
            0.040432,
            0.046872,
            0.033474,
            0.038255,
            0.049530,
            0.088612,
            0.033619,
            0.036886,
            0.085357,
            0.080481,
            0.014753,
            0.039772,
            0.050680,
            0.069577,
            0.058542,
            0.010494,
            0.029916,
            0.064717,
        ]
        relative_rates = [
            0.267828,
            0.984474,
            0.327059,
            1.199805,
            0.000000,
            8.931515,
            0.360016,
            0.232374,
            0.000000,
            0.000000,
            0.887753,
            2.439939,
            1.028509,
            1.348551,
            0.000000,
            1.961167,
            0.000000,
            1.493409,
            11.388659,
            0.000000,
            7.086022,
            2.386111,
            0.087791,
            1.385352,
            1.240981,
            0.107278,
            0.281581,
            0.811907,
            0.228116,
            2.383148,
            5.290024,
            0.868241,
            0.282729,
            6.011613,
            0.439469,
            0.106802,
            0.653416,
            0.632629,
            0.768024,
            0.239248,
            0.438074,
            0.180393,
            0.609526,
            0.000000,
            0.076981,
            0.406431,
            0.154924,
            0.341113,
            0.000000,
            0.000000,
            0.730772,
            0.112880,
            0.071514,
            0.443504,
            2.556685,
            0.258635,
            4.610124,
            3.148371,
            0.716913,
            0.000000,
            1.519078,
            0.830078,
            0.267683,
            0.270475,
            0.460857,
            0.180629,
            0.717840,
            0.896321,
            0.000000,
            0.000000,
            0.000000,
            1.127499,
            0.304803,
            0.170372,
            0.000000,
            3.332732,
            5.230115,
            2.411739,
            0.183641,
            0.136906,
            0.138503,
            0.000000,
            0.000000,
            0.000000,
            0.000000,
            0.153478,
            0.475927,
            1.951951,
            1.565160,
            0.000000,
            0.921860,
            2.485920,
            1.028313,
            0.419244,
            0.133940,
            0.187550,
            1.526188,
            0.507003,
            0.347153,
            0.933709,
            0.119152,
            0.316258,
            0.335419,
            0.170205,
            0.110506,
            4.051870,
            1.531590,
            4.885892,
            0.956097,
            1.598356,
            0.561828,
            0.793999,
            2.322243,
            0.353643,
            0.247955,
            0.171432,
            0.954557,
            0.619951,
            0.459901,
            2.427202,
            3.680365,
            0.265745,
            2.271697,
            0.660930,
            0.162366,
            0.525651,
            0.340156,
            0.306662,
            0.226333,
            1.900739,
            0.331090,
            1.350599,
            1.031534,
            0.136655,
            0.782857,
            5.436674,
            0.000000,
            2.001375,
            0.224968,
            0.000000,
            0.000000,
            0.000000,
            0.000000,
            0.000000,
            0.270564,
            0.000000,
            0.461776,
            0.000000,
            0.000000,
            0.762354,
            0.000000,
            0.740819,
            0.000000,
            0.244139,
            0.078012,
            0.946940,
            0.000000,
            0.953164,
            0.000000,
            0.214717,
            0.000000,
            1.265400,
            0.374834,
            0.286572,
            0.132142,
            0.000000,
            6.952629,
            0.000000,
            0.336289,
            0.417839,
            0.608070,
            2.059564,
            0.240368,
            0.158067,
            0.178316,
            0.484678,
            0.346983,
            0.367250,
            0.538165,
            0.438715,
            8.810038,
            1.745156,
            0.103850,
            2.565955,
            0.123606,
            0.485026,
            0.303836,
            1.561997,
            0.000000,
            0.279379,
        ]
        transition_matrix = np.zeros((num_alleles, num_alleles))
        tril = np.tril_indices(num_alleles, k=-1)
        transition_matrix[tril] = relative_rates
        transition_matrix += np.tril(transition_matrix).T
        transition_matrix *= root_distribution
        row_sums = transition_matrix.sum(axis=1)
        transition_matrix = transition_matrix / max(row_sums)
        row_sums = transition_matrix.sum(axis=1, dtype="float64")
        np.fill_diagonal(transition_matrix, 1.0 - row_sums)
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


class MutationMap:
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
        seed = core.get_random_seed()
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
    if not isinstance(model, BaseMutationModel):
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
