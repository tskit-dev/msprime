#
# Copyright (C) 2018-2020 University of Oxford
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
Configuration and fixtures for pytest. Only put test-suite wide fixtures in here. Module
specific fixtures should live in their modules.

To use a fixture in a test simply refer to it by name as an argument. This is called
dependancy injection. Note that all fixtures should have the suffix "_fixture" to make
it clear in test code.

For example to use the `hapmap_strings_fixture` fixture in a test:

class Something:
    def test_something(self, hapmap_strings_fixture):
        msprime.read_hapmap(hapmap_strings_fixture['original'])

Fixtures can be parameterised etc. see https://docs.pytest.org/en/stable/fixture.html

Note that fixtures have a "scope" for example `ts` below is only created once per
test session and re-used for subsequent tests.
"""
import pytest
from pytest import fixture


def pytest_addoption(parser):
    """
    Add an option to skip tests marked with `@pytest.mark.slow`
    """
    parser.addoption(
        "--skip-slow", action="store_true", default=False, help="Skip slow tests"
    )


def pytest_configure(config):
    """
    Add docs on the "slow" marker
    """
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--skip-slow"):
        skip_slow = pytest.mark.skip(reason="--skip-slow specified")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)


@fixture(scope="session")
def hapmap_strings_fixture():
    """
    A few curated real-world examples of hapmap files
    """
    return {
        "original": """\
            chr pos        rate                    cM
            1   4283592    3.79115663174456        0
            1   4361401    0.0664276817058413      0.294986106359414
            1   7979763   10.9082897515584         0.535345505591925
            1   8007051    0.0976780648822495      0.833010916332456
            1   8762788    0.0899929572085616      0.906829844052373
            1   9477943    0.0864382908650907      0.971188757364862
            1   9696341    4.76495005895746        0.990066707213216
            1   9752154    0.0864316558730679      1.25601286485381
            1   9881751    0.0                     1.26721414815999""",
        "zero_start": """\
            Chromosome  Position(bp)    Rate(cM/Mb)     Map(cM)
            chr1        55550           2.981822        0.000000
            chr1        82571           2.082414        0.080572
            chr1        88169           2.081358        0.092229
            chr1        254996          3.354927        0.439456
            chr1        564598          0.665287        1.478148
            chr1        182973428       2.512769        122.832331
            chr1        183630013       0.000000        124.482178""",
        "nonzero_start": """\
            chrom   pos     recomb_rate     pos_cm
            chr10   48232   0.1614  0.002664
            chr10   48486   0.1589  0.002705
            chr10   50009   0.159   0.002947
            chr10   52147   0.1574  0.003287
            chr10   52541   0.1592  0.003349
            chr10   64718   0.1611  0.005287
            chr10   66015   0.1631  0.005496
            chr10   67284   0.1648  0.005703
            chr10   67994   0.1658  0.00582
            chr10   68368   0.1677  0.005882
            chr10   68839   0.1688  0.005961
            chr10   73953   0.1697  0.006824
            chr10   76169   0  0.0072""",
    }
