#
# Copyright (C) 2020 University of Oxford
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
Test that the documentation examples do what they are supposed to.
"""
import contextlib
import io
import sys
from unittest import mock

import stdpopsim

import msprime
from docs import examples


@contextlib.contextmanager
def capture_stdout():
    new_out = io.StringIO()
    old_out = sys.stdout
    try:
        sys.stdout = new_out
        yield sys.stdout
    finally:
        sys.stdout = old_out


class TestDocumentationExamples:
    def test_ooa_model(self):
        correct_model = stdpopsim.get_species("HomSap").get_demographic_model(
            "OutOfAfrica_3G09"
        )
        ooa_docs = examples.out_of_africa()
        pops = []
        for pop_config in ooa_docs["population_configurations"]:
            pops.append(stdpopsim.Population(id=None, description=None))
            pop_config.sample_size = None

        local_model = stdpopsim.DemographicModel(
            id=None,
            description=None,
            long_description=None,
            generation_time=None,
            populations=pops,
            population_configurations=ooa_docs["population_configurations"],
            migration_matrix=ooa_docs["migration_matrix"],
            demographic_events=ooa_docs["demographic_events"],
        )
        correct_model.verify_equal(local_model)

    def test_segregating_sites(self):
        with capture_stdout() as stdout:
            examples.segregating_sites(10, 5, 10)
        output = stdout.getvalue().splitlines()
        assert len(output) == 3
        assert output[1].startswith("Observed")
        assert output[2].startswith("Analytical")

    def test_migration(self):
        with capture_stdout() as stdout:
            examples.migration_example(10)
        output = stdout.getvalue().splitlines()
        assert len(output) == 2
        assert output[0].startswith("Observed")
        assert output[1].startswith("Predicted")

    def test_logging_info(self):
        with mock.patch("daiquiri.setup"):
            examples.logging_info_example()

    def test_logging_debug(self):
        with mock.patch("daiquiri.setup"):
            examples.logging_debug_example()


def ooa_model():
    """
    Returns the Gutenkunst et al three population OOA model.

    THIS IS A DRAFT.

    Basic conversion of the Demes model. Needs some review, but
    not much point in putting in effort before getting the
    population split event in place.
    """
    populations = [
        msprime.Population(
            name="ancestral",
            description="Equilibrium/root population",
            initial_size=7300,
        ),
        msprime.Population(
            name="AMH", description="Anatomically modern humans", initial_size=12300
        ),
        msprime.Population(
            name="OOA",
            description="Bottleneck out-of-Africa population",
            initial_size=2100,
        ),
        msprime.Population(
            name="YRI",
            description="Yoruba in Ibadan, Nigeria",
            initial_size=12300,
        ),
        msprime.Population(
            name="CEU",
            description=(
                "Utah Residents (CEPH) with Northern and Western European Ancestry"
            ),
            initial_size=1000,
            growth_rate=0.004,
        ),
        msprime.Population(
            name="CHB",
            description="Han Chinese in Beijing, China",
            initial_size=510,
            growth_rate=0.0055,
        ),
    ]

    demography = msprime.Demography(populations)
    # Set the migration rates between extant populations
    demography.set_symmetric_migration_rate(["CEU", "CHB"], 9.6e-5)
    demography.set_symmetric_migration_rate(["YRI", "CHB"], 1.9e-5)
    demography.set_symmetric_migration_rate(["YRI", "CEU"], 3e-5)

    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_OOA = 21.2e3 / generation_time
    T_AMH = 140e3 / generation_time
    T_ANC = 220e3 / generation_time
    yri_ooa_migration_rate = 25e-5

    demography.events = [
        # CEU and CHB merge into the OOA population
        # Zero out the migration matrix first.
        msprime.MigrationRateChange(T_OOA, rate=0),
        msprime.MassMigration(time=T_OOA, source="CEU", dest="OOA", proportion=1),
        msprime.MassMigration(time=T_OOA, source="CHB", dest="OOA", proportion=1),
        # Set the migration rate between OOA and YRI
        msprime.MigrationRateChange(
            time=T_OOA, source="YRI", dest="OOA", rate=yri_ooa_migration_rate
        ),
        msprime.MigrationRateChange(
            time=T_OOA, source="OOA", dest="YRI", rate=yri_ooa_migration_rate
        ),
        # OOA and YRI merge into AMH population
        msprime.MigrationRateChange(T_AMH, rate=0),
        msprime.MassMigration(time=T_AMH, source="OOA", dest="AMH", proportion=1),
        msprime.MassMigration(time=T_AMH, source="YRI", dest="AMH", proportion=1),
        # The AMH population becomes the ancestral population
        msprime.MassMigration(time=T_ANC, source="YRI", dest="ancestral", proportion=1),
    ]
    return demography


def test_ooa():
    demography = ooa_model()
    assert demography.num_populations == 6
