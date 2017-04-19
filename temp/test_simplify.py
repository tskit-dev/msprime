import msprime
from simplify import *

import math
import random

def stepfun_equal(a,b):
    assert(len(a) == 2)
    assert(len(b) == 2)
    assert(len(a[0]) == len(a[1])+1)
    assert(len(b[0]) == len(b[1])+1)
    assert(len(a[0]) == len(b[0]))
    assert(len(a[1]) == len(b[1]))
    for j in [0,1]:
        for k in range(len(a[j])):
            assert(a[j][k] == b[j][k])
    return True

def intervals_equal(a,b):
    assert(len(a) == len(b))
    for k in range(len(a)):
        assert(len(a[k]) == 3)
        assert(len(b[k]) == 3)
        for j in range(3):
            assert(a[k][j] == b[k][j])
    return True

# test add_stepfuns

f = [[0.0, 1.5, 2.0, math.inf], [-1, 1, -1]]
g = [[0.0, 2.0, 5.0, 7.0, math.inf], [1, -1, 2, -1]]
h = [[0.0, math.inf], [-1]]
fgh = add_stepfuns([f, g, h], parent=0, time=0.0, add_fun=lambda x,l,r,p,t: sum([z for z in x if z is not -1]))
true_fgh = [[0.0, 1.5, 2.0, 5.0, 7.0, math.inf],
            [1, 2, 0, 2, 0]]
stepfun_equal(fgh, true_fgh)

# test intervals_to_stepfun
#  and stepfun_to_intervals
fi = [(1.5, 2.0, 1)]
stepfun_equal(f, intervals_to_stepfun(fi))
intervals_equal(fi, stepfun_to_intervals(f))
gi = [(0.0, 2.0, 1), (5.0, 7.0, 2)]
stepfun_equal(g, intervals_to_stepfun(gi))
intervals_equal(gi, stepfun_to_intervals(g))
hi = []
stepfun_equal(h, intervals_to_stepfun(hi))
intervals_equal(hi, stepfun_to_intervals(h))

# test remove_paint
painting = [(0.0, 2.0, 1), (2.5, 4.5, 2), (5.0, 7.0, 2), (10, 12, 1)]
removed = remove_paint(painting, 1.0, 7.0)
intervals_equal(removed, [(1.0, 2.0, 1), (2.5, 4.5, 2), (5.0, 7.0, 2)])
intervals_equal(painting, [(0.0, 1.0, 1), (10, 12, 1)])

painting = [(0.0, 2.0, 1), (2.5, 4.5, 2), (5.0, 7.0, 2), (10, 12, 1)]
removed = remove_paint(painting, 2.5, 6.0)
intervals_equal(removed, [(2.5, 4.5, 2), (5.0, 6.0, 2)])
intervals_equal(painting, [(0.0, 2.0, 1), (6.0, 7.0, 2), (10, 12, 1)])

painting = [(0.0, 2.0, 1), (2.5, 4.5, 2), (5.0, 7.0, 2), (10, 12, 1)]
removed = remove_paint(painting, 2.0, 9.0)
intervals_equal(removed, [(2.5, 4.5, 2), (5.0, 7.0, 2)])
intervals_equal(painting, [(0.0, 2.0, 1), (10, 12, 1)])


random.seed(23)

ts = msprime.simulate(5, recombination_rate=10)

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
