import msprime
from itertools import count
import random
import numpy as np


def random_breakpoint():
    return min(1.0, max(0.0, 2*random.random()-0.5))


def random_mutations(rate):
    nmuts = np.random.poisson(lam=rate)
    return [random.random() for _ in range(nmuts)]


def random_allele():
    return random.choice(['A', 'C', 'G', 'T'])


def wf_sim(N, ngens, survival=0.0, mutation_rate=0.0, debug=False, seed=None):
    '''
    SIMPLE simulation of a bisexual, haploid Wright-Fisher population of size N
    for ngens generations, in which each individual survives with probability
    survival and only those who die are replaced.  The chromosome is 1.0
    Morgans long, and the mutation rate is in units of
    mutations/Morgan/generation.
    '''
    if seed is not None:
        random.seed(seed)
    # initial population
    init_ts = msprime.simulate(N, recombination_rate=1.0)

    nodes = msprime.NodeTable()
    edgesets = msprime.EdgesetTable()
    sites = msprime.SiteTable()
    mutations = msprime.MutationTable()
    init_ts.dump_tables(nodes=nodes, edgesets=edgesets,
                        sites=sites, mutations=mutations)
    nodes.set_columns(time=nodes.time + ngens + 1,
                      flags=nodes.flags,
                      population=nodes.population)
    # searchable
    mut_positions = {}

    # get ready to record things
    labels = count(1 + max(init_ts.samples()), 1)
    pop = [next(labels) for k in range(N)]

    for t in range(ngens, -1, -1):
        if debug:
            print("t:", t)
            print("pop:", pop)

        dead = [(random.random() > survival) for k in pop]
        j = 0
        if debug:
            print("Replacing", sum(dead), "individuals.")
        while j < N:
            while not dead[j]:
                j += 1
            # this is: offspring ID, lparent, rparent, breakpoint
            offspring = next(labels)
            lparent = random.choice(pop)
            rparent = random.choice(pop)
            bp = random_breakpoint()
            muts = random_mutations(mutation_rate)
            if debug:
                print("--->", offspring, lparent, rparent, bp)
            pop[j] = offspring
            j += 1
            nodes.add_row(time=t)
            if bp > 0.0:
                edgesets.add_row(left=0.0, right=bp,
                                 parent=lparent, children=(offspring,))
            if bp < 1.0:
                edgesets.add_row(left=bp, right=1.0,
                                 parent=rparent, children=(offspring,))
            for mut in muts:
                if mut not in mut_positions:
                    mut_positions[mut] = sites.num_rows
                    sites.add_row(site=mut, ancestral_state=random_allele())
                mutations.add_row(site=mut_positions[mut],
                                  node=offspring, derived_state=random_allele())

    if debug:
        print("Done, now sorting.")

    nodes.set_columns(time=nodes.time,
                      flags=[(msprime.NODE_IS_SAMPLE if u in pop else 0)
                             for u in range(nodes.num_rows)],
                      population=nodes.population)
    # msprime.sort_tables(nodes=nodes, edgesets=edgesets, sites=sites,
    #                     mutations=mutations)

    # ts = msprime.load_tables(nodes=nodes, edgesets=edgesets, sites=sites,
    #                          mutations=mutations)

    if debug:
        print("Done.")

    return msprime.TableTuple(nodes=nodes, edgesets=edgesets, sites=sites,
                              mutations=mutations)
