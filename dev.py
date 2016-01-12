"""
Simple client code for development purposes.
"""

from __future__ import print_function
from __future__ import division

import math

import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot

import msprime

def mutations():
    recomb_rates = [(10, 0.05), (20, 0.1), (30, 0), (40, 0.05)]
    for x, rate in recomb_rates:
        print(x, rate)
    max_rate = max(rate for _, rate in recomb_rates)
    print("max_ rate = ", max_rate)

    tree_sequence = msprime.simulate(10, 100, max_rate, random_seed=1)
    for tree in tree_sequence.trees():
        print(tree.get_interval())


def physical_to_genetic(x, recomb_rates):
    s = 0
    last_phys_x = 0
    j = 0
    while j < len(recomb_rates) and x > recomb_rates[j][0]:
        phys_x, recomb_rate = recomb_rates[j]
        s += (phys_x - last_phys_x) * recomb_rate
        j += 1
        last_phys_x = phys_x
    if x != last_phys_x:
        _, recomb_rate = recomb_rates[j]
        s += (x - last_phys_x) * recomb_rate
    return s



def plot_distance_maps(recomb_rates):
    # Plot the piecewise map of physical distance to recombination rate
    x = np.zeros(2 * len(recomb_rates))
    y = np.copy(x)
    last_phys_x = 0
    j = 0
    for phys_x, recomb_rate in recomb_rates:
        x[j] = last_phys_x
        y[j] = recomb_rate
        j += 1
        x[j] = phys_x
        y[j] = recomb_rate
        last_phys_x = phys_x
        j += 1
    pyplot.plot(x, y)
    pyplot.ylim(-0.01, 1.01)
    pyplot.savefig("phys_recomb_rate.png")

    pyplot.clf()

    x = np.zeros(1 + len(recomb_rates))
    y = np.copy(x)
    j = 1
    s = 0
    last_phys_x = 0
    for phys_x, recomb_rate in recomb_rates:
        s += (phys_x - last_phys_x) * recomb_rate
        y[j] = s
        x[j] = phys_x
        j += 1
        last_phys_x = phys_x
    pyplot.plot(x, y)
    physical_dist = 25.6
    genetic_dist = physical_to_genetic(physical_dist, recomb_rates)
    pyplot.axvline(x=physical_dist, color="green")
    pyplot.axhline(y=genetic_dist, color="green")
    pyplot.savefig("phys_genetic_distance.png")


if __name__ == "__main__":
    # mutations()

    plot_distance_maps(
        [(10, 0.1), (11, 1), (20, 0.1), (21, 1), (30, 0.1)]
    )

