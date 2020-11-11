#
# Copyright (C) 2019 University of Oxford
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
Module responsible for computing likelihoods.
"""
import math

from msprime import _msprime


def log_mutation_likelihood(ts, mutation_rate):
    """
    Returns the unnormalised log probability of the stored pattern of mutations
    on the stored tree sequence, assuming infinite sites mutation. In particular,
    each stored site must only contain a single mutation. The omitted normalising
    constant depends on the pattern of mutations, but not on the tree sequence or
    the mutation rate.

    The function first computes the probability of the overall number of mutations
    :math:`M` from the Poisson probability mass function

    .. math::
        e^{-T \\mu / 2} \\frac{(T \\mu / 2)^M}{M!},

    where :math:`T` is the total area of ancestral material in the tree sequence
    stored in units of generations, and :math:`\\mu` is the per-site, per-generation
    mutation probability. Each mutation then contributes an individual factor of
    :math:`l / T`, where :math:`l` is the total branch length on which the mutation
    could have arisen while appearing on all of the required lineages, again stored
    in generations.

    .. warning::
        If a tree at the site of a mutation contains unary nodes, then :math:`l` could
        span more than one edge. In particular, we do not constrain mutations to take
        place on the edge directly above the node on which they have been recorded,
        but rather on any edge which would yield the same configuration of SNPs at the
        leaves of the tree sequence.

    :param tskit.TreeSequence ts: The tree sequence object with mutations.
    :param float mutation_rate: The per-site, per-generation mutation probablity.
        Must be non-negative.
    :return: The unnormalised log probability of the observed SNPs given the tree
        sequence. If the mutation rate is set to zero and the tree sequence contains
        at least one mutation, then returns `-float("inf")`.
    """
    tables = ts.tables
    time = tables.nodes.time
    total_material = 0
    for e in tables.edges:
        total_material += (e.right - e.left) * (time[e.parent] - time[e.child])
    number_of_mutations = len(tables.mutations)
    if mutation_rate == 0:
        if number_of_mutations == 0:
            ret = 0
        else:
            ret = -float("inf")
    else:
        ret = (
            number_of_mutations * math.log(total_material * mutation_rate)
            - total_material * mutation_rate
        )
        for tree in ts.trees():
            for site in tree.sites():
                mutation = site.mutations[0]
                child = mutation.node
                parent = tree.parent(child)
                potential_branch_length = tree.branch_length(child)
                while (
                    tree.parent(parent) is not None and len(tree.children(parent)) == 1
                ):
                    child = parent
                    parent = tree.parent(child)
                    potential_branch_length += tree.branch_length(child)
                child = mutation.node
                while len(tree.children(child)) == 1:
                    child = tree.children(child)[0]
                    potential_branch_length += tree.branch_length(child)
                ret += math.log(potential_branch_length / total_material)
    return ret


def log_arg_likelihood(ts, recombination_rate, Ne=1):
    """
    Returns the log probability of the stored tree sequence under the Hudson ARG.
    An exact expression for this probability is given in equation (1) of
    `Kuhner et al. (2000) <https://www.genetics.org/content/156/3/1393>`_.

    We assume branch lengths stored in generations, resulting in a coalescence
    rate of :math:`1 / (2 N_e)` per pair of lineages.

    .. warning::
        The stored tree sequence must store the full realisation of the ARG,
        including all recombination events and all common ancestor events,
        regardless of whether the recombinations cause a change in the ancestral
        tree or whether the common ancestor events cause coalescence of ancestral
        material. See :ref:`sec_ancestry_full_arg` for details of this
        data structure, and how to generate them using ``msprime``.

    .. warning::
        This method only supports continuous genomes.
        See :ref:`sec_ancestry_discrete_genome` for how these can be specified
        when simulating tree sequences using ``msprime``.

    :param tskit.TreeSequence ts: The tree sequence object.
    :param float recombination_rate: The per-link, per-generation recombination
        probability. Must be non-negative.
    :param float Ne: The diploid effective population size.
    :return: The log probability of the tree sequence under the Hudson ancestral
        recombination graph model. If the recombination rate is zero and the tree
        sequence contains at least one recombination event, then returns
        `-DBL_MAX`.
    """
    # Get the tables into the format we need to interchange with the low-level code.
    lw_tables = _msprime.LightweightTableCollection()
    lw_tables.fromdict(ts.tables.asdict())
    return _msprime.log_likelihood_arg(
        lw_tables, Ne=Ne, recombination_rate=recombination_rate
    )
