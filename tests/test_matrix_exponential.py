import textwrap

import demes
import numpy as np

import msprime


def _matrix_exponential_eigen(A):
    """
    Returns the matrix exponential of A using eigendecomposition.
    https://en.wikipedia.org/wiki/Matrix_exponential
    """
    d, Y = np.linalg.eig(A)
    Yinv = np.linalg.pinv(Y)
    D = np.diag(np.exp(d))
    B = np.matmul(Y, np.matmul(D, Yinv))
    return np.real_if_close(B, tol=1000)


class TestMatrixExponential:
    """
    Test cases for the matrix exponential function.
    """

    def verify(self, A, E=None):
        assert np.max(np.diag(A)) <= 0
        assert np.min(A - np.diag(np.diag(A))) >= 0
        assert np.max(np.sum(A, 1)) <= 1e-10
        E1 = _matrix_exponential_eigen(A)
        E2 = msprime.demography._matrix_exponential(A)
        assert np.min(E2) >= 0
        assert np.allclose(np.sum(E2, 1), np.ones((np.shape(A)[0],)))
        assert E1.shape == E2.shape
        assert np.allclose(E1, E2)
        if E is not None:
            assert E.shape == E2.shape
            assert np.allclose(E, E2)

    def test_zeros(self):
        for j in range(1, 10):
            A = np.zeros((j, j))
            self.verify(A, np.eye(j))

    def test_ones_minus_diagonal(self):
        # If we got to larger values we start getting complex number results.
        # (k x k) matrices of ones, but with (-k) on the diagonal, for k >= 2.
        for j in range(2, 5):
            A = np.ones((j, j))
            A = A - (j * np.eye(j))
            E = np.exp(-j) * np.eye(j) + (1 - np.exp(-j)) * np.ones((j, j)) / j
            self.verify(A, E)

    def test_identity_exp(self):
        # (-1) * np.eye(k), compared to exp(-1) * np.eye(k)
        for k in range(2, 5):
            A = (-1) * np.eye(k)
            for i in range(k):
                A[i, (i + 1) % k] = 1.0
            self.verify(A)

    def test_diagonal_minus_row_sum(self):
        for k in range(2, 5):
            A = np.arange(k**2).reshape(k, k)
            np.fill_diagonal(A, -1 * np.sum(A, 1) + np.diag(A))
            self.verify(A)

        rng = np.random.default_rng(1234)
        for k in [2, 5, 30]:
            for _ in range(20):
                A = rng.uniform(size=(k, k))
                np.fill_diagonal(A, -1 * np.sum(A, 1) + np.diag(A))
                self.verify(A)


class TestDemographyTrajectoriesEdgeCase:
    def test_edge_case(self):
        # This test exposes unstable behaviour in some matrix exponential
        # functions, namely scipy.linalg.expm() and _matrix_exponential_eigen()
        # above. The instability leads to coalescence rate trajectories that
        # have negative values and a probability vector that is non monotonic.
        # See issue #1775.
        model_str = textwrap.dedent(
            """\
            time_units: generations
            defaults:
              epoch: {start_size: 1000}
            demes:
            - name: A
            - name: B
              ancestors: [A]
              start_time: 6000
            - name: C
              ancestors: [B]
              start_time: 2000
            - name: D
              ancestors: [C]
              start_time: 1000
            migrations:
            - demes: [A, D]
              rate: 1e-5
            """
        )
        model = demes.loads(model_str)
        ddb = msprime.Demography.from_demes(model).debug()
        last_N = max(ddb.population_size_history[:, ddb.num_epochs - 1])
        last_epoch = ddb.epoch_start_time[-1]
        T = np.linspace(0, last_epoch + 12 * last_N, 101)
        rates, probs = ddb.coalescence_rate_trajectory(
            T, {"A": 1, "C": 1}, double_step_validation=False
        )
        assert np.all(rates >= 0)
        assert np.all(probs >= 0)
        assert np.all(probs <= 1)
        assert np.all(np.diff(probs) <= 0)
