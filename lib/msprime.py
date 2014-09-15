"""
Python version of the msprime multilocus coalescent algorithm.
"""
from __future__ import print_function
from __future__ import division

import sys
import math
import random

import bintrees

import fenwick


class Segment(object):
    """
    A class representing a single segment. Each segment has a left
    and right, denoting the loci over which it spans, a value and a
    next, giving the next in the chain.
    """
    def __init__(self, index):
        self.left = None
        self.right = None
        self.value = None
        self.prev = None
        self.next = None
        self.index = index

    def __str__(self):
        s = "({0}:{1}-{2}->{3}: prev={4} next={5})".format(self.index,
                self.left, self.right, self.value, repr(self.prev),
                repr(self.next))
        return s


class ReferenceSimulator(object):
    """
    A reference implementation of the multi locus simulation algorithm.
    """
    def __init__(self, n, m, r, max_segments=100):
        self.n = n
        self.m = m
        self.r = r
        self.max_segments = max_segments
        self.L = fenwick.FenwickTree(self.max_segments)
        self.segment_stack = []
        self.segments = [None for j in range(self.max_segments + 1)]
        for j in range(self.max_segments):
            s = Segment(j + 1)
            self.segments[j + 1] = s
            self.segment_stack.append(s)
        self.P = [None for j in range(n)]
        for j in range(n):
            x = self.alloc_segment(1, m, j + 1)
            self.L.set_value(x.index, m - 1)
            self.P[j] = x
        self.B = bintrees.AVLTree()
        self.B[1] = n + 1
        self.B[m + 1] = 0
        self.C = []
        self.t = 0
        self.num_ca_events = 0
        self.num_re_events = 0
        # for the asymptote output
        self.picture_num = 0

    def alloc_segment(self, left, right, value, prev=None, next=None):
        """
        Pops a new segment off the stack and sets its properties.
        """
        s = self.segment_stack.pop()
        s.left = left
        s.right = right
        s.value = value
        s.next = next
        s.prev = prev
        return s

    def free_segment(self, u):
        """
        Frees the specified segment making it ready for reuse and
        setting its weight to zero.
        """
        self.L.set_value(u.index, 0)
        self.segment_stack.append(u)

    def simulate(self):
        """
        Simulates the algorithm until all loci have coalesced.
        """
        while len(self.P) != 0:
            #self.population_size.append(len(self.P))
            #self.num_links.append(self.L.get_total())
            #print(len(self.P), "\t", self.L.get_total())
            #self.print_state()
            #self.output_asy_state()
            self.verify()
            lambda_r = self.r * self.L.get_total()
            if lambda_r == 0:
                t_recomb = sys.float_info.max
            else:
                t_recomb = random.expovariate(lambda_r)
            lambda_c = len(self.P) * (len(self.P) - 1)
            t_coancestry = random.expovariate(lambda_c)
            if t_recomb < t_coancestry:
                self.t += t_recomb
                self.recombination_event()
            else:
                self.t += t_coancestry
                self.coancestry_event()


    def recombination_event(self):
        """
        Implements a recombination event.
        """
        self.num_re_events += 1
        l = random.randint(1, self.L.get_total())
        # Get the segment containing the l'th link
        j = self.L.find(l)
        g = self.L.get_cumulative_frequency(j) - l
        y = self.segments[j]
        k = y.right - g - 1
        x = y.prev
        if y.left <= k:
            # Make new segment
            z = self.alloc_segment(k + 1, y.right, y.value, None, y.next)
            if y.next is not None:
                y.next.prev = z
            y.next = None
            y.right = k
            self.L.increment(y.index, k - z.right)
            if k + 1 not in self.B:
                self.B[k + 1] = self.B[self.B.floor_key(k)]
        else:
            # split the link between x and y.
            x.next = None
            y.prev = None
            z = y
        self.L.set_value(z.index, z.right - z.left)
        self.P.append(z)

    def coancestry_event(self):
        """
        Implements a coancestry event.
        """
        self.num_ca_events += 1
        # Choose two ancestors uniformly.
        j = random.randint(0, len(self.P) - 1)
        x = self.P[j]
        del self.P[j]
        j = random.randint(0, len(self.P) - 1)
        y = self.P[j]
        del self.P[j]
        z = None
        while x is not None or y is not None:
            alpha = None
            if x is None or y is None:
                if x is not None:
                    alpha = x
                    x = None
                if y is not None:
                    alpha = y
                    y = None
            else:
                if y.left < x.left:
                    beta = x
                    x = y
                    y = beta
                if x.right < y.left:
                    alpha = x
                    x = x.next
                    alpha.next = None
                elif x.left != y.left:
                    alpha = self.alloc_segment(x.left, y.left - 1, x.value)
                    x.left = y.left
                else:
                    l = x.left
                    eta = self.B[l]
                    r_max = min(x.right, y.right)
                    self.B[l] = eta + 1
                    r = self.B.succ_key(l)
                    while self.B[r] == eta and r < r_max:
                        self.B[r] = eta + 1
                        r = self.B.succ_key(r)
                    r -= 1
                    self.C.append((l, r, x.value, y.value, eta, self.t))
                    if eta < 2 * self.n - 1:
                        alpha = self.alloc_segment(l, r, eta)
                    if x.right == r:
                        self.free_segment(x)
                        x = x.next
                    else:
                        x.left = r + 1
                    if y.right == r:
                        self.free_segment(y)
                        y = y.next
                    else:
                        y.left = r + 1
            if alpha is not None:
                l = alpha.left
                if z == None:
                    self.P.append(alpha)
                else:
                    z.next = alpha
                    l = z.right
                alpha.prev = z
                z = alpha
                self.L.set_value(alpha.index, alpha.right - l)


    def fill_tree(self, k):
        """
        Returns the tree for the specified locus filled out from the
        extant coalescence records.
        """
        pi = [0 for j in range(2 * self.n)]
        tau = [0 for j in range(2 * self.n)]
        for l, r, c1, c2, p, t in self.C:
            if l <= k <= r:
                pi[c1] = p
                pi[c2] = p
                tau[p] = t
        return pi, tau

    def output_asy_state(self):
        """
        Writes the state of the algorithm out as an asymptote picture.
        """
        self.picture_num += 1
        filename = "output/state_{0}.asy".format(self.picture_num)
        with open(filename, "w") as f:
            print("size(16cm,0);", file=f)
            b_table = "\\begin{array}{cc}"
            for j, eta in self.B.items():
                b_table += "{0}&{1}\\\\ ".format(j, eta)
            b_table += "\\end{array}"
            print("string b_table = \"${0}$\";".format(b_table), file=f)
            print("label(\"$B$\", (0,0), N);", file=f)
            print("label(b_table, (0,0), S);", file=f)
            l_table = "\\begin{array}{ccc}"
            for j in range(1, self.max_segments + 1):
                s = self.L.get_frequency(j)
                c = self.L.get_cumulative_frequency(j)
                l_table += "s_{{{0}}}&{1}&{2}\\\\ ".format(j, s, c)
            l_table += "\\end{array}"
            print("string l_table = \"${0}$\";".format(l_table), file=f)
            print("label(\"$L$\", (20,0), N);", file=f)
            print("label(l_table, (20,0), S);", file=f)
            # TODO add a seperate function for this and come up with a pictorial
            # representaion.
            for j, u in enumerate(self.P):
                s = "$"
                while u is not None:
                    s += "({0},{1},{2},s_{{{3}}})".format(u.left, u.right, u.value,
                            u.index)
                    u = u.next
                s += "$"
                print("label(\"{0}\", (-100,{1}), E);".format(s, j * -5), file=f)
            c_table = "\\begin{array}{cccccc}\\\\"
            c_table += "\\ell& r& c_1& c_2& p& t\\\\"
            for l, r, c1, c2, p, t in self.C:
                c_table += "{0}&{1}&{2}&{3}&{4}&{5:.3f}\\\\ ".format(l, r, c1, c2, p, t)
            c_table += "\\end{array}"
            print("string c_table = \"${0}$\";".format(c_table), file=f)
            print("label(c_table, (-100,{1}), E);".format(c_table, -50), file=f)

    def print_state(self):
        print("State @ time ", self.t)
        print("Links = ", self.L.get_total())
        print("Population:", len(self.P))
        q = 0
        for u in self.P:
            s = ""
            while u is not None:
                s += "({0}-{1}->{2}({3}))".format(u.left, u.right, u.value,
                        u.index)
                u = u.next
            print("\t" + s)
        print("Breakpoints", len(self.B))
        for j, eta in self.B.items():
            pi, tau = self.fill_tree(j)
            if j <= self.m:
                print("\t", j, ":", eta, ":", pi, ":", tau)
        print("Fenwick tree:", self.L.get_total())
        for j in range(1, self.max_segments + 1):
            s = self.L.get_frequency(j)
            if s != 0:
                print("\t", j, "->", s,
                        self.L.get_cumulative_frequency(j))
        print("Coalescence records: ")
        for rec in self.C:
            print("\t", rec)
        self.verify()

    def verify(self):
        """
        Checks that the state of the simulator is consistent.
        """
        q = 0
        for u in self.P:
            assert u.prev == None
            left = u.left
            right = u.left
            while u is not None:
                assert u.left <= u.right
                s = u.right - right
                assert s == self.L.get_frequency(u.index)
                right = u.right
                v = u.next
                if v is not None:
                    #print(u, v)
                    assert v.prev == u
                u = v
            q += right - left
        assert q == self.L.get_total()
        assert self.B[self.m + 1] == 0

    def verify_end(self):
        """
        Verify the state of the simulation at the end.
        """
        # Verify that the coalescence records correspond
        # to the breakpoints.
        B = set([self.m + 1])
        for l, r, c1, c2, p, t in self.C:
            B.add(l)
        assert B == set(self.B.keys())

class PopulationModel(object):
    """
    Superclass of PopulationModels.
    """
    def __init__(self):
        self.start_time = 0.0
        self.initial_size = 1.0

    def set_start_time(self, start_time):
        """
        Sets the time that this PopulationModel begins at to the
        specified value.
        """
        self.start_time = start_time

    def set_initial_size(self, initial_size):
        """
        Sets the size of the population at the start of the period
        defined by this population model to the specified value.
        """
        self.initial_size = initial_size

class ConstantPopulation(PopulationModel):
    """
    Class representing a population of constant size. This size
    is expressed as a fraction of the initial population size.
    """
    def __init__(self, size):
        self.size = size

    def get_size(self, t):
        """
        Returns the size of the population at time t as a fraction
        of the population size at time 0.
        """
        return self.size

    def get_waiting_time(self, lambda_c, t):
        """
        Returns the waiting time until the common ancestor event
        when the unscaled rate of coancestry events is lambda_c.
        """
        t_wait = random.expovariate(lambda_c) * self.size
        return t_wait

class ExponentialPopulation(PopulationModel):
    """
    Class representing an exponentially growing population in which
    the population size after start_time is given by exp(alpha * (t -
    start_time)).
    """
    def __init__(self, alpha):
        self.start_time = 0.0
        self.alpha = alpha

    def get_size(self, t):
        """
        Returns the size of the population at time t as a fraction
        of the population size at time 0. This is given by
        exp(-alpha * (t - start_time))
        """
        return math.exp(-self.alpha * (t - self.start_time))

    def get_waiting_time(self, lambda_c, t):
        dt = t - self.start_time
        z = self.alpha * self.initial_size * math.exp(-self.alpha * dt)
        z = 1  + z * random.expovariate(lambda_c)
        ret = sys.float_info.max
        # if z is <= 0 no coancestry can occur
        if z > 0:
            ret =math.log(z) / self.alpha
        return ret



class DemographicSimulator(ReferenceSimulator):
    """
    Adds demographic events to the reference simulator.
    """
    def __init__(self, n, m, r, max_segments=100):
        super(DemographicSimulator, self).__init__(n, m , r, max_segments)
        self.population_models = []
        self.population_models.append((-1, ConstantPopulation(1.0)))
        self.population_models.append((sys.float_info.max, None))

    def add_demographic_event(self, time, pop_model):
        """
        Adds a demographic event at the specified time in the population's
        history. After this time, the size of the population is driven
        by the specified population model.
        """
        pop_model.set_start_time(time)
        self.population_models.append((time, pop_model))

    def simulate(self):
        """
        Simulates the algorithm until all loci have coalesced.
        """
        self.population_models.sort(key=lambda x: x[0])
        pop_size = 1.0
        start_time, pop_model = self.population_models.pop(0)
        while len(self.P) != 0:
            self.verify()
            lambda_r = self.r * self.L.get_total()
            if lambda_r == 0:
                t_recomb = sys.float_info.max
            else:
                t_recomb = random.expovariate(lambda_r)
            lambda_c = len(self.P) * (len(self.P) - 1)
            t_coancestry = pop_model.get_waiting_time(lambda_c, self.t)
            if t_recomb < t_coancestry:
                waiting_time = t_recomb
                event_method = self.recombination_event
            else:
                waiting_time = t_coancestry
                event_method = self.coancestry_event
            if self.t + waiting_time >= self.population_models[0][0]:
                old_model = pop_model
                start_time, pop_model = self.population_models.pop(0)
                pop_model.set_initial_size(old_model.get_size(self.t))
                self.t = start_time
            else:
                self.t += waiting_time
                event_method()



def harmonic_number(n):
    """
    Returns the nth Harmonic number.
    """
    return sum(1 / k for k in range(1, n + 1))

def run_replicates():
    random.seed(1)
    n = 15
    m = 20
    rho = 0 # 0.1
    for j in range(1000):
        s = DemographicSimulator(n, m, rho, 100)
        s.add_demographic_event(0.5, ExponentialPopulation(0.1))
        s.add_demographic_event(0.7, ExponentialPopulation(2.0))
        s.simulate()
        print(s.t, "\t", s.num_re_events, s.num_ca_events)



def main():
    import numpy as np
    from matplotlib import pyplot
    random.seed(1)
    n = 15
    m = 1
    rho = 1e-19
    #s = ReferenceSimulator(n, m, rho, 100)
    s = DemographicSimulator(n, m, rho, 100)
    s.add_demographic_event(0.0, ExponentialPopulation(0.1))
    s.add_demographic_event(0.1, ExponentialPopulation(0.2))
    s.add_demographic_event(0.4, ConstantPopulation(2))
    s.simulate()
    s.print_state()
    s.verify_end()

if __name__ == "__main__":
    #main()
    run_replicates()
