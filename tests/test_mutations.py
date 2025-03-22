#
# Copyright (C) 2018-2021 University of Oxford
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
Test cases for mutation generation.
"""
from __future__ import annotations

import numpy as np
import pytest

import msprime


class TestMatrixMutationModel:

    @pytest.mark.parametrize(
        "p, m, lo, hi",
        [
            (0.9, 0.9, 2, 50),
            (0.5, 0.25, 50, 500),
        ],
    )
    def test_TPM(self, p, m, lo, hi):
        model = msprime.TPM(p=p, m=m, lo=lo, hi=hi)
        self.validate_model(model)
        self.validate_stationary(model)

        # test for ratio of TPM trans mat vals
        for i in range(lo, hi + 1):
            if lo < i < hi:
                ii = i - lo
                exp = (p + (m / (1 - (1 - m) ** (hi - i))) * (1 - p)) / (
                    p + (m / (1 - (1 - m) ** (i - lo)) * (1 - p))
                )
                obs = (
                    model.transition_matrix[ii, ii + 1]
                    / model.transition_matrix[ii, ii - 1]
                )
                assert np.isclose(exp, obs)
