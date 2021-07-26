import textwrap

import demes
import numpy as np
import pytest
import scipy.linalg

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
        assert np.max(np.sum(A, 1)) <= 0
        E1 = scipy.linalg.expm(A)
        E2 = msprime.demography._matrix_exponential(A)
        E3 = _matrix_exponential_eigen(A)
        assert np.min(E2) >= 0
        assert np.allclose(np.sum(E2, 1), np.ones((np.shape(A)[0],)))
        assert E1.shape == E2.shape == E3.shape
        assert np.allclose(E1, E2)
        assert np.allclose(E2, E3)
        if E is not None:
            assert E.shape == E2.shape == E3.shape
            assert np.allclose(E, E2)
            assert np.allclose(E, E3)

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


class TestDemographyTrajectoriesEdgeCase:
    @pytest.mark.slow
    def test_edge_case(self, monkeypatch):
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
        rates1, P1 = ddb.coalescence_rate_trajectory(
            T, {"A": 1, "C": 1}, double_step_validation=False
        )
        assert np.all(rates1 >= 0)
        assert np.all(P1 >= 0)
        assert np.all(P1 <= 1)
        assert np.all(np.diff(P1) <= 0)

        # scipy.linalg.expm() is unstable for this example.
        # Although results are similar to above, we get probabilities > 1.
        monkeypatch.setattr(
            msprime.demography, "_matrix_exponential", scipy.linalg.expm
        )
        rates2, P2 = ddb.coalescence_rate_trajectory(
            T, {"A": 1, "C": 1}, double_step_validation=False
        )
        assert np.allclose(rates1, rates2)
        assert np.allclose(P1, P2, rtol=1e-2)
        assert np.all(rates2 >= 0)
        assert np.all(P2 >= 0)
        # Don't confirm failure here in case scipy improves in the future.
        # with pytest.raises(AssertionError):
        #    assert np.all(P2 <= 1)
        # with pytest.raises(AssertionError):
        #    assert np.all(np.diff(P2) <= 0)

        # _matrix_exponential_eigen() is unstable for this example, and results
        # are not even in the right ball park. The trajectory produces rates
        # below zero and a probability vector that isn't close to monotonic.
        monkeypatch.setattr(
            msprime.demography, "_matrix_exponential", _matrix_exponential_eigen
        )
        rates3, P3 = ddb.coalescence_rate_trajectory(
            T, {"A": 1, "C": 1}, double_step_validation=False
        )
        # Fails on Linux with openblas, but not on Github Action's MacOS blas.
        # with pytest.raises(AssertionError):
        #    assert np.all(rates3 >= 0)
        with pytest.raises(AssertionError):
            assert np.all(np.diff(P3) <= 0)
