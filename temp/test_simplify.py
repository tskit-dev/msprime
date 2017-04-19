import unittest
import itertools
import six

import msprime

from simplify import *

import math
import random

def assert_paintings_equal(a,b):
    assert(len(a) == len(b))
    for k in range(len(a)):
        assert(len(a[k]) == 3)
        assert(len(b[k]) == 3)
        for j in range(3):
            assert(a[k][j] == b[k][j])
    return True

# test add_paintings

f = [[1.5, 2.0, 1]]
g = [[0.0, 2.0, 1], [5.0, 7.0, 2]]
h = []

def this_add_fun(x,left,right,parent,time):
    if all([(z is -1) for z in x]):
        return -1
    else:
        return sum([z for z in x if z is not -1])

fgh = add_paintings([f, g, h], parent=0, time=0.0, 
                    add_fun=this_add_fun)
true_fgh = [[0.0, 1.5, 1], [1.5, 2.0, 2], [5.0, 7.0, 2]]
assert_paintings_equal(fgh, true_fgh)

# test remove_paint
painting = [(0.0, 2.0, 1), (2.5, 4.5, 2), (5.0, 7.0, 2), (10, 12, 1)]
removed = remove_paint(painting, 1.0, 7.0)
assert_paintings_equal(removed, [(1.0, 2.0, 1), (2.5, 4.5, 2), (5.0, 7.0, 2)])
assert_paintings_equal(painting, [(0.0, 1.0, 1), (10, 12, 1)])

painting = [(0.0, 2.0, 1), (2.5, 4.5, 2), (5.0, 7.0, 2), (10, 12, 1)]
removed = remove_paint(painting, 2.5, 6.0)
assert_paintings_equal(removed, [(2.5, 4.5, 2), (5.0, 6.0, 2)])
assert_paintings_equal(painting, [(0.0, 2.0, 1), (6.0, 7.0, 2), (10, 12, 1)])

painting = [(0.0, 2.0, 1), (2.5, 4.5, 2), (5.0, 7.0, 2), (10, 12, 1)]
removed = remove_paint(painting, 2.0, 9.0)
assert_paintings_equal(removed, [(2.5, 4.5, 2), (5.0, 7.0, 2)])
assert_paintings_equal(painting, [(0.0, 2.0, 1), (10, 12, 1)])


random.seed(23)

ts = msprime.simulate(5, recombination_rate=5)

simp_ts = simplify(ts, [0,1,2])

ts_simp = ts.simplify([0,1,2])

st = simp_ts.trees()
ts = ts_simp.trees()

print("----------------")

for t in simp_ts.trees():
    print(t.interval, t)

print("----------------")

for t in ts_simp.trees():
    print(t.interval, t)
