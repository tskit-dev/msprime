import random

import numpy as np
import tskit


def wright_fisher(N, T, L=100, random_seed=None):
    """
    Simulate a Wright-Fisher population of N haploid individuals with L
    discrete loci for T generations. Based on Algorithm W from
    https://www.biorxiv.org/content/biorxiv/early/2018/01/16/248500.full.pdf
    """
    random.seed(random_seed)
    tables = tskit.TableCollection(L)
    P = np.arange(N, dtype=int)
    # Mark the initial generation as samples so that we remember these nodes.
    for _ in range(N):
        tables.nodes.add_row(time=T, flags=tskit.NODE_IS_SAMPLE)
    t = T
    while t > 0:
        t -= 1
        Pp = P.copy()
        for j in range(N):
            u = tables.nodes.add_row(time=t, flags=0)
            Pp[j] = u
            a = random.randint(0, N - 1)
            b = random.randint(0, N - 1)
            x = random.randint(1, L - 1)
            tables.edges.add_row(0, x, P[a], u)
            tables.edges.add_row(x, L, P[b], u)
        P = Pp

    # Now do some table manipulations to ensure that the tree sequence
    # that we output has the form that msprime needs to finish the
    # simulation. Much of the complexity here is caused by the tables API
    # not allowing direct access to memory, which will change soon.

    # Mark the extant population as samples also
    flags = tables.nodes.flags
    flags[P] = tskit.NODE_IS_SAMPLE
    tables.nodes.set_columns(flags=flags, time=tables.nodes.time)
    tables.sort()
    # Simplify with respect to the current generation, but ensuring we keep the
    # ancient nodes from the initial population.
    tables.simplify()
    # Unmark the initial generation as samples
    flags = tables.nodes.flags
    time = tables.nodes.time
    flags[:] = 0
    flags[time == 0] = tskit.NODE_IS_SAMPLE
    # The final tables must also have at least one population which
    # the samples are assigned to
    tables.populations.add_row()
    tables.nodes.set_columns(
        flags=flags, time=time, population=np.zeros_like(tables.nodes.population)
    )
    return tables.tree_sequence()
