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
Utilities for working with intervals and interval maps.
"""
import bisect
import gzip
import warnings

import numpy as np

import _msprime


class IntervalMap:
    """
    A class mapping a numeric value to a set of adjacent intervals along the
    genome.
    """

    def __init__(self, position, value=None):
        size = len(self.position) - 1
        if size < 1:
            raise ValueError("Must have at least two positions to define an interval")
        # We take copies to make sure there are no unintended consequences
        # later when a user modifies the arrays in this map.
        self.position = np.array(position, copy=True)
        if value is None:
            self.value = np.zeros(size)
        else:
            self.value = np.array(value, copy=True)
        if self.value.shape[0] != size:
            raise ValueError(
                "Value array must have one less entry than the position array."
            )
        if self.position[0] != 0:
            raise ValueError("First position must be zero")
        if np.any(self.position[:1] >= self.position[1:]):
            raise ValueError("Position values must be in increasing order")
        # TODO continue.


# TODO refactor to use the new IntervalMap API but without breaking the old
# exported API.


class RecombinationMap:
    """
    A RecombinationMap represents the changing rates of recombination
    along a chromosome. This is defined via two lists of numbers:
    ``positions`` and ``rates``, which must be of the same length.
    Given an index j in these lists, the rate of recombination
    per base per generation is ``rates[j]`` over the interval
    ``positions[j]`` to ``positions[j + 1]``. Consequently, the first
    position must be zero, and by convention the last rate value
    is also required to be zero (although it is not used).

    .. warning::
        The ``num_loci`` parameter is deprecated.

    :param list positions: The positions (in bases) denoting the
        distinct intervals where recombination rates change. These can
        be floating point values.
    :param list rates: The list of rates corresponding to the supplied
        ``positions``. Recombination rates are specified per base,
        per generation.
    :param int num_loci: **This parameter is deprecated**.
        The maximum number of non-recombining loci
        in the underlying simulation. By default this is set to
        the largest possible value, allowing the maximum resolution
        in the recombination process. However, for a finite sites
        model this can be set to smaller values.
    """

    def __init__(self, positions, rates, num_loci=None, map_start=0):
        if num_loci is not None:
            if num_loci == positions[-1]:
                warnings.warn("num_loci is no longer supported and should not be used.")
            else:
                raise ValueError(
                    "num_loci does not match sequence length. "
                    "To set a discrete number of recombination sites, "
                    "scale positions to span the desired number of loci "
                    "and set discrete=True"
                )
        self._ll_recombination_map = _msprime.RecombinationMap(positions, rates)
        self.map_start = map_start

    @classmethod
    def uniform_map(cls, length, rate, num_loci=None):
        """
        Returns a :class:`.RecombinationMap` instance in which the recombination
        rate is constant over a chromosome of the specified length.
        The legacy ``num_loci`` option is no longer supported and should not be used.

        :param float length: The length of the chromosome.
        :param float rate: The rate of recombination per unit of sequence length
            along this chromosome.
        :param int num_loci: This parameter is no longer supported.
        """
        return cls([0, length], [rate, 0], num_loci=num_loci)

    @classmethod
    def read_hapmap(cls, filename):
        """
        Parses the specified file in HapMap format. These files must be
        white-space-delimited, and contain a single header line (which is
        ignored), and then each subsequent line contains the starting position
        and recombination rate for the segment from that position (inclusive)
        to the starting position on the next line (exclusive). Starting
        positions of each segment are given in units of bases, and
        recombination rates in centimorgans/Megabase. The first column in this
        file is ignored, as are additional columns after the third (Position is
        assumed to be the second column, and Rate is assumed to be the third).
        If the first starting position is not equal to zero, then a
        zero-recombination region is inserted at the start of the chromosome.

        A sample of this format is as follows::

            Chromosome	Position(bp)	Rate(cM/Mb)	Map(cM)
            chr1	55550	        2.981822	0.000000
            chr1	82571	        2.082414	0.080572
            chr1	88169	        2.081358	0.092229
            chr1	254996	        3.354927	0.439456
            chr1	564598	        2.887498	1.478148
            ...
            chr1	182973428	2.512769	122.832331
            chr1	183630013	0.000000	124.482178

        :param str filename: The name of the file to be parsed. This may be
            in plain text or gzipped plain text.
        """
        positions = []
        rates = []
        if filename.endswith(".gz"):
            f = gzip.open(filename)
        else:
            f = open(filename)
        try:
            # Skip the header line
            f.readline()
            for j, line in enumerate(f):
                pos, rate, = map(float, line.split()[1:3])
                if j == 0:
                    map_start = pos
                    if pos != 0:
                        positions.append(0)
                        rates.append(0)
                positions.append(pos)
                # Rate is expressed in centimorgans per megabase, which
                # we convert to per-base rates
                rates.append(rate * 1e-8)
            if rate != 0:
                raise ValueError(
                    "The last rate provided in the recombination map must be zero"
                )
        finally:
            f.close()
        return cls(positions, rates, map_start=map_start)

    @property
    def mean_recombination_rate(self):
        """
        Return the weighted mean recombination rate
        across all windows of the entire recombination map.
        """
        chrom_length = self._ll_recombination_map.get_sequence_length()

        positions = self._ll_recombination_map.get_positions()
        positions_diff = self._ll_recombination_map.get_positions()[1:]
        positions_diff = np.append(positions_diff, chrom_length)
        window_sizes = positions_diff - positions

        weights = window_sizes / chrom_length
        if self.map_start != 0:
            weights[0] = 0
        rates = self._ll_recombination_map.get_rates()

        return np.average(rates, weights=weights)

    def slice(self, start=None, end=None, trim=False):  # noqa: A003
        """
        Returns a subset of this recombination map between the specified end
        points. If start is None, it defaults to 0. If end is None, it defaults
        to the end of the map. If trim is True, remove the flanking
        zero recombination rate regions such that the sequence length of the
        new recombination map is end - start.
        """
        positions = self.get_positions()
        rates = self.get_rates()

        if start is None:
            i = 0
            start = 0
        if end is None:
            end = positions[-1]
            j = len(positions)

        if (
            start < 0
            or end < 0
            or start > positions[-1]
            or end > positions[-1]
            or start > end
        ):
            raise IndexError(f"Invalid subset: start={start}, end={end}")

        if start != 0:
            i = bisect.bisect_left(positions, start)
            if start < positions[i]:
                i -= 1
        if end != positions[-1]:
            j = bisect.bisect_right(positions, end, lo=i)

        new_positions = list(positions[i:j])
        new_rates = list(rates[i:j])
        new_positions[0] = start
        if end > new_positions[-1]:
            new_positions.append(end)
            new_rates.append(0)
        else:
            new_rates[-1] = 0
        if trim:
            new_positions = [pos - start for pos in new_positions]
        else:
            if new_positions[0] != 0:
                if new_rates[0] == 0:
                    new_positions[0] = 0
                else:
                    new_positions.insert(0, 0)
                    new_rates.insert(0, 0.0)
            if new_positions[-1] != positions[-1]:
                new_positions.append(positions[-1])
                new_rates.append(0)
        return self.__class__(new_positions, new_rates)

    def __getitem__(self, key):
        """
        Use slice syntax for obtaining a recombination map subset. E.g.
            >>> recomb_map_4m_to_5m = recomb_map[4e6:5e6]
        """
        if not isinstance(key, slice) or key.step is not None:
            raise TypeError("Only interval slicing is supported")
        start, end = key.start, key.stop
        if start is not None and start < 0:
            start += self.get_sequence_length()
        if end is not None and end < 0:
            end += self.get_sequence_length()
        return self.slice(start=start, end=end, trim=True)

    def get_ll_recombination_map(self):
        return self._ll_recombination_map

    def physical_to_genetic(self, physical_x):
        return self._ll_recombination_map.position_to_mass(physical_x)

    def physical_to_discrete_genetic(self, physical_x):
        raise ValueError("Discrete genetic space is no longer supported")

    def genetic_to_physical(self, genetic_x):
        return self._ll_recombination_map.mass_to_position(genetic_x)

    def get_total_recombination_rate(self):
        return self._ll_recombination_map.get_total_recombination_rate()

    def get_per_locus_recombination_rate(self):
        raise ValueError("Genetic loci are no longer supported")

    def get_size(self):
        return self._ll_recombination_map.get_size()

    def get_num_loci(self):
        raise ValueError("num_loci is no longer supported")

    def get_positions(self):
        # For compatability with existing code we convert to a list
        return list(self._ll_recombination_map.get_positions())

    def get_rates(self):
        # For compatability with existing code we convert to a list
        return list(self._ll_recombination_map.get_rates())

    def get_sequence_length(self):
        return self._ll_recombination_map.get_sequence_length()

    def get_length(self):
        # Deprecated: use sequence_length instead
        return self.get_sequence_length()

    def asdict(self):
        return {
            "positions": self.get_positions(),
            "rates": self.get_rates(),
            "map_start": self.map_start,
        }
