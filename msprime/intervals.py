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
import warnings

import numpy as np


# TODO this is a minimal implementation of the functionality we need to get
# varying-rates-of-things-along the genome into msprime. It's deliberately
# lightweight, and unencumbered by any backing classes from the low-level
# model. The intention would be to add a bunch of useful functionality here
# using numpy APIs.
class RateMap:
    """
    A class mapping a numeric value to a set of adjacent intervals along the
    genome.
    """

    def __init__(self, position, rate, map_start=0):
        # We take copies to make sure there are no unintended consequences
        # later when a user modifies the arrays in this map.
        self.position = np.array(position, dtype=float, copy=True)
        self.rate = np.array(rate, dtype=float, copy=True)
        self.map_start = map_start
        size = len(self.rate)
        if size < 1:
            raise ValueError("Must have at least two positions")
        if self.position.shape[0] != size + 1:
            raise ValueError(
                "Rate array must have one less entry than the position array."
            )
        if self.position[0] != 0:
            raise ValueError("First position must be zero")
        if np.any(self.position[:1] >= self.position[1:]):
            raise ValueError("Position values must be in increasing order")
        if np.any(self.rate < 0):
            raise ValueError("Rates must be non-negative")

        self.cumulative = np.insert(np.cumsum(np.diff(self.position) * self.rate), 0, 0)

    def __len__(self):
        return len(self.rate)

    @staticmethod
    def uniform(sequence_length, rate):
        return RateMap([0, sequence_length], [rate])

    def asdict(self):
        return {"position": self.position, "rate": self.rate}

    @property
    def sequence_length(self):
        return self.position[-1]

    @property
    def total_mass(self):
        return self.cumulative[-1]

    @property
    def size(self):
        return self.rate.shape[0]

    @property
    def mean_rate(self):
        """
        Return the weighted mean across all windows of the entire map.
        """
        window_sizes = self.position[1:] - self.position[:-1]
        weights = window_sizes / self.sequence_length
        if self.map_start != 0:
            weights[0] = 0
        return np.average(self.rate, weights=weights)

    def slice(self, start=None, end=None, trim=False):  # noqa: A003
        """
        Returns a subset of this rate map between the specified end
        points. If start is None, it defaults to 0. If end is None, it defaults
        to the end of the map. If trim is True, remove the flanking
        zero rate regions such that the sequence length of the
        new rate map is end - start.
        """
        if start is None:
            i = 0
            start = 0
        if end is None:
            end = self.position[-1]
            j = len(self.position)

        if (
            start < 0
            or end < 0
            or start > self.position[-1]
            or end > self.position[-1]
            or start > end
        ):
            raise IndexError(f"Invalid subset: start={start}, end={end}")

        if start != 0:
            i = np.searchsorted(self.position, start, side="left")
            if start < self.position[i]:
                i -= 1
        if end != self.position[-1]:
            j = i + np.searchsorted(self.position[i:], end, side="right")
            if end > self.position[j - 1]:
                j += 1

        position = self.position[i:j].copy()
        rate = self.rate[i : j - 1].copy()
        position[0] = start
        position[-1] = end
        map_start = 0

        if trim:
            position -= start
        else:
            # Prepend or extend zero-rate region at start of map.
            if position[0] != 0:
                map_start = position[0]  # TODO: is this what we want here?
                if rate[0] == 0:
                    position[0] = 0
                else:
                    position = np.insert(position, 0, 0)
                    rate = np.insert(rate, 0, 0)
            # Append or extend zero-rate region at end of map.
            if position[-1] != self.position[-1]:
                if rate[-1] == 0:
                    position[-1] = self.position[-1]
                else:
                    position = np.append(position, self.position[-1])
                    rate = np.append(rate, 0)

        return self.__class__(position, rate, map_start=map_start)

    def __getitem__(self, key):
        """
        Use slice syntax for obtaining a rate map subset. E.g.
            >>> rate_map_4m_to_5m = rate_map[4e6:5e6]
        """
        if not isinstance(key, slice) or key.step is not None:
            raise TypeError("Only interval slicing is supported")
        start, end = key.start, key.stop
        if start is not None and start < 0:
            start += self.sequence_length
        if end is not None and end < 0:
            end += self.sequence_length
        return self.slice(start=start, end=end, trim=True)


# The RecombinationMap class is deprecated since 1.0. We maintain the
# functionality where it is possible to do so.

# TODO update the documentation to make it clear that this is a legacy
# interface and is deprecated.


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
        self.map = RateMap(positions, rates[:-1])
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
        Parses the specified file in HapMap format.

        .. warning::
            This method is deprecated, use the module-level
            :func:`read_hapmap` function instead.

        :param str filename: The name of the file to be parsed. This may be
            in plain text or gzipped plain text.
        :return: A RecombinationMap object.
        """
        warnings.warn(
            "RecombinationMap.read_hapmap() is deprecated. "
            "Use msprime.read_hapmap() instead.",
            FutureWarning,
        )
        rate_map = read_hapmap(filename)
        rate = np.append(rate_map.rate, 0)
        return cls(rate_map.position, rate, map_start=rate_map.map_start)

    @property
    def mean_recombination_rate(self):
        """
        Return the weighted mean recombination rate
        across all windows of the entire recombination map.
        """
        return self.map.mean_rate

    def get_total_recombination_rate(self):
        """
        Returns the effective recombination rate for this genetic map.
        This is the weighted mean of the rates across all intervals.
        """
        return self.map.cumulative[-1]

    def physical_to_genetic(self, x):
        return np.interp(x, self.map.position, self.map.cumulative)

    def genetic_to_physical(self, genetic_x):
        if self.map.cumulative[-1] == 0:
            # If we have a zero recombination rate throughout then everything
            # except L maps to 0.
            return self.get_sequence_length() if genetic_x > 0 else 0
        if genetic_x == 0:
            return self.map.position[0]
        index = np.searchsorted(self.map.cumulative, genetic_x) - 1
        y = (
            self.map.position[index]
            + (genetic_x - self.map.cumulative[index]) / self.map.rate[index]
        )
        return y

    def physical_to_discrete_genetic(self, physical_x):
        raise ValueError("Discrete genetic space is no longer supported")

    def get_per_locus_recombination_rate(self):
        raise ValueError("Genetic loci are no longer supported")

    def get_num_loci(self):
        raise ValueError("num_loci is no longer supported")

    def get_size(self):
        return self.map.size + 1

    def get_positions(self):
        return list(self.map.position)

    def get_rates(self):
        return list(self.map.rate) + [0]

    def get_sequence_length(self):
        return self.map.sequence_length

    def get_length(self):
        # Deprecated: use get_sequence_length() instead
        return self.get_sequence_length()

    def asdict(self):
        return self.map.asdict()


def read_hapmap(filename):
    # Black barfs with an INTERNAL_ERROR trying to reformat this docstring,
    # so we explicitly disable reformatting here.
    # fmt: off
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
    :return: A RateMap object.
    """
    # fmt: on
    hapmap = np.loadtxt(filename, skiprows=1, usecols=(1, 2))
    position = hapmap[:, 0]
    # Rate is expressed in centimorgans per megabase, which
    # we convert to per-base rates
    rate = 1e-8 * hapmap[:, 1]

    map_start = position[0]
    if map_start != 0:
        position = np.insert(position, 0, 0)
        rate = np.insert(rate, 0, 0)
    if rate[-1] != 0:
        raise ValueError("The last rate provided in the recombination map must be zero")
    return RateMap(position, rate[:-1], map_start=map_start)
