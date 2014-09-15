"""
An implementation of the Fenwick tree (also known as a Binary Indexed Tree)
based on "A new data structure for cumulative frequency tables", Software 
Practice and Experience, Vol 24, No 3, pp 327 336 Mar 1994. This 
implementation supports any non-negative frequencies, and the search 
procedure always returns the smallest index such that its cumulative 
frequency <= f. This search procedure is a slightly modified version 
of that presented in Tech Report 110, "A new data structure for cumulative 
frequency tables: an improved frequency-to-symbol algorithm." available 
at https://www.cs.auckland.ac.nz/~peter-f/FTPfiles/TechRep110.ps
"""
from __future__ import print_function
from __future__ import division

import random

class FenwickTree(object):
    """
    A Fenwick Tree to represent cumulative frequency tables over
    integers. Each index from 1 to max_index initially has a 
    zero frequency.
    """
    def __init__(self, max_index):
        assert max_index > 0
        self.__max_index = max_index
        self.__tree = [0 for j in range(max_index + 1)]
        # Compute the binary logarithm of max_index
        u = self.__max_index
        while u != 0:
            self.__log_max_index = u
            u -= (u & -u)

    def get_total(self):
        """
        Returns the total cumulative frequency over all indexes.
        """
        #return self.__tree[self.__max_index]
        return self.get_cumulative_frequency(self.__max_index)

    def increment(self, index, v):
        """
        Increments the frequency of the specified index by the specified 
        value.
        """
        assert 0 < index <= self.__max_index 
        j = index
        while j <= self.__max_index:
            self.__tree[j] += v
            j += (j & -j)
    
    def set_value(self, index, v):
        """
        Sets the frequency at the specified index to the specified value.
        """
        f = self.get_frequency(index)
        self.increment(index, v - f)

    def get_cumulative_frequency(self, index):
        """
        Returns the cumulative frequency of the specified index. 
        """
        assert 0 < index <= self.__max_index 
        j = index
        s = 0
        while j > 0:
            s += self.__tree[j]
            j -= (j & -j)
        return s 

    def get_frequency(self, index):
        """
        Returns the frequency of the specified index. 
        """
        assert 0 < index <= self.__max_index 
        j = index
        v = self.__tree[j] 
        p = j & (j - 1)
        j -= 1
        while p != j:
            v -= self.__tree[j]
            j = j & (j - 1)
        return v
    
    def find(self, v):
        """
        Returns the smallest index with cumulative sum >= v.
        """
        j = 0
        s = v 
        half = self.__log_max_index
        while half > 0:
            # Skip non-existant entries
            while j + half > self.__max_index:
                half >>= 1
            k = j + half
            if s > self.__tree[k]:
                j = k
                s -= self.__tree[j]
            half >>= 1
        return j + 1


class NaiveFenwickTree(object):
    """
    A class providing the same interface as the Fenwick tree, but implemented
    in a naive fashion for testing.
    """
    def __init__(self, max_index):
        self.__max_index = max_index
        self.__v = [0 for j in range(max_index + 1)]
        self.__c = [0 for j in range(max_index + 1)]

    def increment(self, index, v):
        self.__v[index] += v
        for j in range(index, self.__max_index + 1):
            self.__c[j] += v

    def get_cumulative_frequency(self, index):
        return self.__c[index]

    def get_frequency(self, index):
        return self.__v[index] 
    
    def find(self, v):
        """
        Returns the smallest index with cumulative sum >= v.
        """
        j = 1 
        while j < self.__max_index and self.__c[j] < v:
            j += 1
        return j


def test_simple_cases():
    for n in range(1, 20):
        v = range(n + 1) 
        f = FenwickTree(n)
        t = NaiveFenwickTree(n)
        for j in range(1, n + 1):
            f.increment(j, v[j])
            t.increment(j, v[j])
            c = sum(v[:j + 1])
            assert f.get_frequency(j) == v[j]
            assert t.get_frequency(j) == v[j]
            assert f.get_cumulative_frequency(j) == c
            assert t.get_cumulative_frequency(j) == c
            assert f.find(c) == j
            assert t.find(c) == j


def test_list(values):
    """
    Testst the specified list of values.
    """
    n = len(values)
    v = [0] + values
    s = sum(v)
    f = FenwickTree(n)
    g = NaiveFenwickTree(n)
    for j in range(1, n + 1):
        f.increment(j, v[j])
        g.increment(j, v[j])
    for j in range(1, n + 1):
        assert f.get_frequency(j) == g.get_frequency(j)
        c = f.get_cumulative_frequency(j) 
        assert c == g.get_cumulative_frequency(j)
        if v[j] != 0:
            assert f.find(c) == j
            assert g.find(c) == j
        # Generate a random value and find it
        if s != 0:
            c = random.randint(1, s)
            assert f.find(c) == g.find(c)
    # test extreme values
    assert f.find(0) == 1
    assert g.find(0) == 1
    assert f.find(v[0] - 1) == 1
    assert g.find(v[0] - 1) == 1
    # k is the largest non-zero index
    if s == 0:
        k = 1
    else:
        k = n
        while v[k] == 0 and k > 0:
            k -= 1
    for val in [s, s + 1, 2 * s]:
        assert f.find(s) == k
        assert g.find(s) == k


def test_random_cases(N):
    
    for n in range(1, N):
        v = [random.randint(1, 10) for j in range(n)]
        test_list(v)
        # test more sparse values
        v = [random.randint(0, 1) for j in range(n)]
        test_list(v)
        v = range(n)
        test_list(v)
        v = list(reversed(range(n)))
        test_list(v)
        random.shuffle(v)
        test_list(v)


def main():
    test_simple_cases() 
    test_random_cases(200)

if __name__ == "__main__":
    main()
