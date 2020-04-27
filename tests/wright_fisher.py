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
A simple Wright-Fisher simulator.
"""
import random

import numpy as np
import tskit

import msprime


class WrightFisherSimulator:
    """
    SIMPLE simulation of a bisexual, haploid Wright-Fisher population of size N
    for ngens generations, in which each individual survives with probability
    survival and only those who die are replaced.  If num_loci is None,
    the chromosome is 1.0 Morgans long, and the mutation rate is in units of
    mutations/Morgan/generation. If num_loci not None, a discrete recombination
    model is used where breakpoints are chosen uniformly from 1 to num_loci - 1.
    """

    def __init__(
        self,
        N,
        survival=0.0,
        seed=None,
        deep_history=True,
        debug=False,
        initial_generation_samples=False,
        num_loci=None,
    ):
        self.N = N
        self.num_loci = num_loci
        self.survival = survival
        self.deep_history = deep_history
        self.debug = debug
        self.initial_generation_samples = initial_generation_samples
        self.seed = seed
        self.rng = random.Random(seed)

    def random_breakpoint(self):
        if self.num_loci is None:
            return min(1.0, max(0.0, 2 * self.rng.random() - 0.5))
        else:
            return self.rng.randint(1, self.num_loci - 1)

    def run(self, ngens):
        L = 1
        if self.num_loci is not None:
            L = self.num_loci
        tables = tskit.TableCollection(sequence_length=L)
        tables.populations.add_row()
        if self.deep_history:
            # initial population
            init_ts = msprime.simulate(
                self.N, recombination_rate=1.0, length=L, random_seed=self.seed
            )
            init_tables = init_ts.dump_tables()
            flags = init_tables.nodes.flags
            if not self.initial_generation_samples:
                flags = np.zeros_like(init_tables.nodes.flags)
            tables.nodes.set_columns(time=init_tables.nodes.time + ngens, flags=flags)
            tables.edges.set_columns(
                left=init_tables.edges.left,
                right=init_tables.edges.right,
                parent=init_tables.edges.parent,
                child=init_tables.edges.child,
            )
        else:
            flags = 0
            if self.initial_generation_samples:
                flags = tskit.NODE_IS_SAMPLE
            for _ in range(self.N):
                tables.nodes.add_row(flags=flags, time=ngens, population=0)

        pop = list(range(self.N))
        for t in range(ngens - 1, -1, -1):
            if self.debug:
                print("t:", t)
                print("pop:", pop)

            dead = [self.rng.random() > self.survival for k in pop]
            # sample these first so that all parents are from the previous gen
            new_parents = [
                (self.rng.choice(pop), self.rng.choice(pop)) for k in range(sum(dead))
            ]
            k = 0
            if self.debug:
                print("Replacing", sum(dead), "individuals.")
            for j in range(self.N):
                if dead[j]:
                    # this is: offspring ID, lparent, rparent, breakpoint
                    offspring = len(tables.nodes)
                    tables.nodes.add_row(time=t, population=0)
                    lparent, rparent = new_parents[k]
                    k += 1
                    bp = self.random_breakpoint()
                    if self.debug:
                        print("--->", offspring, lparent, rparent, bp)
                    pop[j] = offspring
                    if bp > 0.0:
                        tables.edges.add_row(
                            left=0.0, right=bp, parent=lparent, child=offspring
                        )
                    if bp < L:
                        tables.edges.add_row(
                            left=bp, right=L, parent=rparent, child=offspring
                        )

        if self.debug:
            print("Done! Final pop:")
            print(pop)
        flags = tables.nodes.flags
        flags[pop] = tskit.NODE_IS_SAMPLE
        tables.nodes.set_columns(
            flags=flags, time=tables.nodes.time, population=tables.nodes.population
        )
        return tables


def wf_sim(
    N,
    ngens,
    survival=0.0,
    deep_history=True,
    debug=False,
    seed=None,
    initial_generation_samples=False,
    num_loci=None,
):
    sim = WrightFisherSimulator(
        N,
        survival=survival,
        deep_history=deep_history,
        debug=debug,
        seed=seed,
        initial_generation_samples=initial_generation_samples,
        num_loci=num_loci,
    )
    return sim.run(ngens)
