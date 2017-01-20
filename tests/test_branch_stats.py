
import unittest
import random

import msprime
import _msprime

##
# This tests implementation of the algorithm described in branch-lengths-methods.md
##


def build_tree_sequence(records, mutations=[]):
    ts = _msprime.TreeSequence()
    ts.load_records(records)
    ts.set_mutations(mutations)
    return msprime.TreeSequence(ts)


class BranchStatsTestCase(unittest.TestCase):
    """
    Tests of branch statistic computation.
    """
    random_seed = 123456

    def assertListAlmostEqual(self, x, y):
        for a, b in zip(x, y):
            self.assertAlmostEqual(a, b)

    def check_vectorization(self, ts):
        samples = random.sample(ts.samples(), 3)
        A = [[samples[0]], [samples[1]], [samples[2]]]

        def f(x):
            return [float((x[0] > 0) != (x[1] > 0)),
                    float((x[0] > 0) != (x[2] > 0)),
                    float((x[1] > 0) != (x[2] > 0))]

        self.assertListAlmostEqual(
                msprime.branch_stats_vector_node_iter(ts, A, f, method='mutations'),
                [ts.pairwise_diversity(samples=[samples[0], samples[1]]),
                 ts.pairwise_diversity(samples=[samples[0], samples[2]]),
                 ts.pairwise_diversity(samples=[samples[1], samples[2]])])
        self.assertListAlmostEqual(
                msprime.branch_stats_vector_node_iter(ts, A, f, method='length'),
                [msprime.branch_length_diversity(ts, A[0], A[1]),
                 msprime.branch_length_diversity(ts, A[0], A[2]),
                 msprime.branch_length_diversity(ts, A[1], A[2])])
        self.assertListAlmostEqual(
                msprime.branch_stats_vector(ts, A, f),
                [msprime.branch_length_diversity(ts, A[0], A[1]),
                 msprime.branch_length_diversity(ts, A[0], A[2]),
                 msprime.branch_length_diversity(ts, A[1], A[2])])

    def check_pairwise_diversity(self, ts):
        samples = random.sample(ts.samples(), 2)
        A_one = [[samples[0]], [samples[1]]]
        A_many = [random.sample(ts.samples(), 2),
                  random.sample(ts.samples(), 2)]
        for A in (A_one, A_many):
            n = [len(a) for a in A]

            def f(x):
                return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/float(n[0]*n[1])

            self.assertAlmostEqual(
                    msprime.branch_stats_node_iter(ts, A, f, method='length'),
                    msprime.branch_length_diversity(ts, A[0], A[1]))
            self.assertAlmostEqual(
                    msprime.branch_stats(ts, A, f),
                    msprime.branch_length_diversity(ts, A[0], A[1]))

    def check_pairwise_diversity_mutations(self, ts):
        samples = random.sample(ts.samples(), 2)
        A = [[samples[0]], [samples[1]]]
        n = [len(a) for a in A]

        def f(x):
            return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/float(n[0]*n[1])

        self.assertAlmostEqual(
                msprime.branch_stats_node_iter(ts, A, f, method='mutations'),
                ts.pairwise_diversity(samples=samples))

    def check_Y_stat(self, ts):
        samples = random.sample(ts.samples(), 3)
        A = [[samples[0]], samples[1:3]]

        def f(x):
            return float(((x[0] == 1) and (x[1] == 0)) or ((x[0] == 0) and (x[1] == 2)))

        self.assertAlmostEqual(
                msprime.branch_stats(ts, A, f),
                msprime.branch_length_Y(ts, A[0][0], A[1][0], A[1][1]))
        self.assertAlmostEqual(
                msprime.branch_stats_node_iter(ts, A, f, method='length'),
                msprime.branch_length_Y(ts, A[0][0], A[1][0], A[1][1]))

    def check_f4_stat(self, ts):
        samples = random.sample(ts.samples(), 4)
        A_zero = [[samples[0]], [samples[0]], [samples[1]], [samples[1]]]
        A_f1 = [[samples[0]], [samples[1]], [samples[0]], [samples[1]]]
        A_one = [[samples[0]], [samples[1]], [samples[2]], [samples[3]]]
        A_many = [random.sample(ts.samples(), 3),
                  random.sample(ts.samples(), 3),
                  random.sample(ts.samples(), 3),
                  random.sample(ts.samples(), 3)]
        for A in (A_zero, A_f1, A_one, A_many):

            def f(x):
                return ((float(x[0])/len(A[0])-float(x[1])/len(A[1]))
                        * (float(x[2])/len(A[2])-float(x[3])/len(A[3])))

            self.assertAlmostEqual(
                    msprime.branch_stats(ts, A, f),
                    msprime.branch_length_f4(ts, A[0], A[1], A[2], A[3]))
            self.assertAlmostEqual(
                    msprime.branch_stats_node_iter(ts, A, f, method='length'),
                    msprime.branch_length_f4(ts, A[0], A[1], A[2], A[3]))

    def test_pairwise_diversity(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, recombination_rate=100)
        self.check_pairwise_diversity(ts)
        self.check_pairwise_diversity_mutations(ts)

    def test_Y_stat(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, recombination_rate=100)
        self.check_Y_stat(ts)

    def test_f4(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, recombination_rate=100)
        self.check_f4_stat(ts)

    def test_vectorization(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, recombination_rate=100)
        self.check_vectorization(ts)

    def test_case_1(self):
        # With mutations:
        #
        # 1.0          6
        # 0.7         / \                                    5
        #            /   X                                  / \
        # 0.5       X     4                4               /   4
        #          /     / \              / \             /   X X
        # 0.4     X     X   \            X   3           X   /   \
        #        /     /     X          /   / X         /   /     \
        # 0.0   0     1       2        1   0   2       0   1       2
        #          (0.0, 0.2),        (0.2, 0.8),       (0.8, 1.0)
        #

        true_diversity_01 = 2*(1 * (0.2-0) + 0.5 * (0.8-0.2) + 0.7 * (1.0-0.8))
        true_diversity_02 = 2*(1 * (0.2-0) + 0.4 * (0.8-0.2) + 0.7 * (1.0-0.8))
        true_diversity_12 = 2*(0.5 * (0.2-0) + 0.5 * (0.8-0.2) + 0.5 * (1.0-0.8))

        records = [
            msprime.CoalescenceRecord(
                left=0.2, right=0.8, node=3, children=(0, 2),
                time=0.4, population=0),
            msprime.CoalescenceRecord(
                left=0.0, right=0.2, node=4, children=(1, 2),
                time=0.5, population=0),
            msprime.CoalescenceRecord(
                left=0.2, right=0.8, node=4, children=(1, 3),
                time=0.5, population=0),
            msprime.CoalescenceRecord(
                left=0.8, right=1.0, node=4, children=(1, 2),
                time=0.5, population=0),
            msprime.CoalescenceRecord(
                left=0.8, right=1.0, node=5, children=(0, 4),
                time=0.7, population=0),
            msprime.CoalescenceRecord(
                left=0.0, right=0.2, node=6, children=(0, 4),
                time=1.0, population=0)]

        ll_ts = _msprime.TreeSequence()
        ll_ts.load_records(records)
        ts = msprime.TreeSequence(ll_ts)

        ts.set_mutations([
            (0.1, 0),
            (0.15, 0),
            (0.15, 1),
            (0.05, 4),
            (0.1, 2),
            (0.3, 1),
            (0.6, 2),
            (0.9, 0),
            (0.95, 1),
            (0.95, 2)])

        self.check_pairwise_diversity(ts)
        self.check_pairwise_diversity_mutations(ts)
        self.check_Y_stat(ts)
        self.check_vectorization(ts)

        # diversity between 0 and 1
        A = [[0], [1]]

        def f(x):
            return float((x[0] > 0) != (x[1] > 0))

        # branch lengths:
        self.assertAlmostEqual(msprime.branch_length_diversity(ts, [0], [1]),
                               true_diversity_01)
        self.assertAlmostEqual(msprime.branch_stats(ts, A, f),
                               true_diversity_01)
        self.assertAlmostEqual(msprime.branch_stats_node_iter(ts, A, f),
                               true_diversity_01)

        # mean diversity between [0, 1] and [0, 2]:
        true_mean_diversity = (0 + true_diversity_02
                               + true_diversity_01 + true_diversity_12)/4
        A = [[0, 1], [0, 2]]
        n = [len(a) for a in A]

        def f(x):
            return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/4.0

        # branch lengths:
        self.assertAlmostEqual(msprime.branch_length_diversity(ts, A[0], A[1]),
                               true_mean_diversity)
        self.assertAlmostEqual(msprime.branch_stats(ts, A, f),
                               true_mean_diversity)
        self.assertAlmostEqual(msprime.branch_stats_node_iter(ts, A, f),
                               true_mean_diversity)

        # Y-statistic for (0/12)
        A = [[0], [1, 2]]

        def f(x):
            return ((x[0] == 1) and (x[1] == 0)) or ((x[0] == 0) and (x[1] == 2))

        # branch lengths:
        true_Y = 0.2*(1 + 0.5) + 0.6*(0.4) + 0.2*(0.7+0.2)
        self.assertAlmostEqual(msprime.branch_length_Y(ts, 0, 1, 2), true_Y)
        self.assertAlmostEqual(msprime.branch_stats(ts, A, f), true_Y)
        self.assertAlmostEqual(msprime.branch_stats_node_iter(ts, A, f), true_Y)

    def test_case_2(self):
        # Here are the trees:
        # t                  |              |              |             |
        #
        # 0       --3--      |     --3--    |     --3--    |    --3--    |    --3--
        #        /  |  \     |    /  |  \   |    /     \   |   /     \   |   /     \
        # 1     4   |   5    |   4   |   5  |   4       5  |  4       5  |  4       5
        #       |\ / \ /|    |   |\   \     |   |\     /   |  |\     /   |  |\     /|
        # 2     | 6   7 |    |   | 6   7    |   | 6   7    |  | 6   7    |  | 6   7 |
        #       | |\ /| |    |   |  \  |    |   |  \  |    |  |  \       |  |  \    | ...
        # 3     | | 8 | |    |   |   8 |    |   |   8 |    |  |   8      |  |   8   |
        #       | |/ \| |    |   |  /  |    |   |  /  |    |  |  / \     |  |  / \  |
        # 4     | 9  10 |    |   | 9  10    |   | 9  10    |  | 9  10    |  | 9  10 |
        #       |/ \ / \|    |   |  \   \   |   |  \   \   |  |  \   \   |  |  \    |
        # 5     0   1   2    |   0   1   2  |   0   1   2  |  0   1   2  |  0   1   2
        #
        #                    |   0.0 - 0.1  |   0.1 - 0.2  |  0.2 - 0.4  |  0.4 - 0.5
        # ... continued:
        # t                  |             |             |             |
        #
        # 0         --3--    |    --3--    |    --3--    |    --3--    |    --3--
        #          /     \   |   /     \   |   /     \   |   /     \   |   /  |  \
        # 1       4       5  |  4       5  |  4       5  |  4       5  |  4   |   5
        #         |\     /|  |   \     /|  |   \     /|  |   \     /|  |     /   /|
        # 2       | 6   7 |  |    6   7 |  |    6   7 |  |    6   7 |  |    6   7 |
        #         |  \    |  |     \    |  |       /  |  |    |  /  |  |    |  /  |
        # 3  ...  |   8   |  |      8   |  |      8   |  |    | 8   |  |    | 8   |
        #         |  / \  |  |     / \  |  |     / \  |  |    |  \  |  |    |  \  |
        # 4       | 9  10 |  |    9  10 |  |    9  10 |  |    9  10 |  |    9  10 |
        #         |    /  |  |   /   /  |  |   /   /  |  |   /   /  |  |   /   /  |
        # 5       0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2
        #
        #         0.5 - 0.6  |  0.6 - 0.7  |  0.7 - 0.8  |  0.8 - 0.9  |  0.9 - 1.0

        # divergence betw 0 and 1
        true_diversity_01 = 2*(0.6*4 + 0.2*2 + 0.2*5)
        # divergence betw 1 and 2
        true_diversity_12 = 2*(0.2*5 + 0.2*2 + 0.3*5 + 0.3*4)
        # divergence betw 0 and 2
        true_diversity_02 = 2*(0.2*5 + 0.2*4 + 0.3*5 + 0.1*4 + 0.2*5)
        # mean divergence between 0, 1 and 0, 2
        true_mean_diversity = (0 + true_diversity_02
                               + true_diversity_01 + true_diversity_12)/4
        # Y(0;1, 2)
        true_Y = 0.2*4 + 0.2*(4+2) + 0.2*4 + 0.2*2 + 0.2*(5+1)

        records = [
                msprime.CoalescenceRecord(
                    left=0.5, right=1.0, node=10, children=(1,),
                    time=5.0-4.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.0, right=0.4, node=10, children=(2,),
                    time=5.0-4.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.6, right=1.0, node=9, children=(0,),
                    time=5.0-4.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.0, right=0.5, node=9, children=(1,),
                    time=5.0-4.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.8, right=1.0, node=8, children=(10,),
                    time=5.0-3.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.2, right=0.8, node=8, children=(9, 10),
                    time=5.0-3.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.0, right=0.2, node=8, children=(9,),
                    time=5.0-3.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.7, right=1.0, node=7, children=(8,),
                    time=5.0-2.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.0, right=0.2, node=7, children=(10,),
                    time=5.0-2.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.8, right=1.0, node=6, children=(9,),
                    time=5.0-2.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.0, right=0.7, node=6, children=(8,),
                    time=5.0-2.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.4, right=1.0, node=5, children=(2, 7),
                    time=5.0-1.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.1, right=0.4, node=5, children=(7,),
                    time=5.0-1.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.6, right=0.9, node=4, children=(6,),
                    time=5.0-1.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.0, right=0.6, node=4, children=(0, 6),
                    time=5.0-1.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.9, right=1.0, node=3, children=(4, 5, 6),
                    time=5.0-0.0, population=0),
                msprime.CoalescenceRecord(
                    left=0.1, right=0.9, node=3, children=(4, 5),
                    time=5.0-0.0, population=0),
                msprime.CoalescenceRecord(
                        left=0.0, right=0.1, node=3, children=(4, 5, 7),
                        time=5.0-0.0, population=0),
               ]

        ll_ts = _msprime.TreeSequence()
        ll_ts.load_records(records)
        ts = msprime.TreeSequence(ll_ts)

        self.check_pairwise_diversity(ts)
        self.check_pairwise_diversity_mutations(ts)
        self.check_Y_stat(ts)
        self.check_vectorization(ts)

        # divergence between 0 and 1
        A = [[0], [1]]

        def f(x):
            return (x[0] > 0) != (x[1] > 0)

        # branch lengths:
        self.assertAlmostEqual(msprime.branch_length_diversity(ts, [0], [1]),
                               true_diversity_01)
        self.assertAlmostEqual(msprime.branch_stats(ts, A, f),
                               true_diversity_01)
        self.assertAlmostEqual(msprime.branch_stats_node_iter(ts, A, f),
                               true_diversity_01)

        # mean divergence between 0, 1 and 0, 2
        A = [[0, 1], [0, 2]]
        n = [len(a) for a in A]

        def f(x):
            return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/4.0

        # branch lengths:
        self.assertAlmostEqual(msprime.branch_length_diversity(ts, A[0], A[1]),
                               true_mean_diversity)
        self.assertAlmostEqual(msprime.branch_stats(ts, A, f),
                               true_mean_diversity)
        self.assertAlmostEqual(msprime.branch_stats_node_iter(ts, A, f),
                               true_mean_diversity)

        # Y-statistic for (0/12)
        A = [[0], [1, 2]]

        def f(x):
            return ((x[0] == 1) and (x[1] == 0)) or ((x[0] == 0) and (x[1] == 2))

        # branch lengths:
        self.assertAlmostEqual(msprime.branch_length_Y(ts, 0, 1, 2), true_Y)
        self.assertAlmostEqual(msprime.branch_stats(ts, A, f), true_Y)
        self.assertAlmostEqual(msprime.branch_stats_node_iter(ts, A, f), true_Y)
