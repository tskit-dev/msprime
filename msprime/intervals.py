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

    :param list position: A list of :math:`n+1` positions, starting at 0, and ending
        in the sequence length over which the RateMap will apply.
    :param list rate: A list of :math:`n` positive rates that apply between each
        position.
    :param float start_position: When returning a mean_rate over the genome, ignore the
        part of the map from zero up to ``start_position``, as long as the rate is 0 in
        that region. If ``None`` (default), the ``start_position`` is 0.
    :param float end_position: When returning a mean_rate over the genome, ignore the
        part of the map from ``end_position`` onwards, as long as the rate is 0 in that
        region. If ``None`` (default), the ``end_position`` is the end of the map.
    """

    def __init__(
        self,
        position,
        rate,
        *,
        start_position=None,
        end_position=None,
    ):
        # Making the arrays read-only guarantees rate and cumulative mass stay in sync
        self._position = np.array(position, dtype=float)
        self._position.flags.writeable = False
        size = len(self._position)
        if size < 2:
            raise ValueError("Must have at least two positions")
        if self._position[0] != 0:
            raise ValueError("First position must be zero")
        if np.any(np.diff(self._position) <= 0):
            bad_pos = np.where(np.diff(self._position) <= 0)[0] + 1
            raise ValueError(
                f"Position values not strictly increasing at indexes {bad_pos}"
            )

        self._rate = np.array(rate, dtype=float)
        self._rate.flags.writeable = False
        if len(self._rate) != size - 1:
            raise ValueError(
                "Rate array must have one less entry than the position array"
            )
        if np.any(self._rate < 0):
            bad_rates = np.where(self._rate < 0)[0]
            raise ValueError(f"Rate values negative at indexes {bad_rates}")

        self._cumulative_mass = np.insert(
            np.cumsum(np.diff(self._position) * self._rate), 0, 0
        )
        self._cumulative_mass.flags.writeable = False

        if start_position is None:
            start_position = 0
        elif self.get_cumulative_mass(start_position) > 0:
            raise ValueError("Rates before the start_position must all be zero")
        self._start_position = start_position

        if end_position is None:
            end_position = self.sequence_length
        elif self.get_cumulative_mass(end_position) < self.total_mass:
            raise ValueError("Rates after the end_position must all be zero")
        self._end_position = end_position

        # Make all of the internal arrays read-only, so they can't get out of sync

    def __len__(self):  # TODO - remove this, as the length of a ratemap is not obvious
        return len(self.rate)

    @property
    def position(self):
        """
        The positions between which a constant rate applies. There will be
        one more position than the number of rates, with the first position being 0
        and the last being the :attr:`sequence_length`.
        """
        return self._position

    @property
    def rate(self):
        """
        The array of constant rates that apply between each :attr:`position`.
        """
        return self._rate

    @property
    def cumulative_mass(self):
        """
        The integral of the rates along the genome, one value for each
        :attr:`position`. The first value will be 0 and the last will be the
        :attr:`total_mass`.
        """
        return self._cumulative_mass

    @property
    def start_position(self):
        """
        The physical start of the map. This position can be greater than 0 for maps
        which have no recorded rate at the start of a chromosome (often the telomeric
        region). It primarily affects the :attr:`mean_rate` returned for this rate map.
        """
        return self._start_position

    @start_position.setter
    def start_position(self, val):
        raise NotImplementedError(
            "Cannot change a RateMap start_position: use ratemap.slice(... trim=False)"
        )

    @property
    def end_position(self):
        """
        The physical end of the map. This position can be less than the
        :attr:`sequence_length` for maps which have no recorded rate at the end of a
        chromosome (often the telomeric region). It primarily affects the
        :attr:`mean_rate` returned for this rate map.
        """
        return self._end_position

    @end_position.setter
    def end_position(self, val):
        raise NotImplementedError(
            "Cannot change a RateMap end_position: use ratemap.slice(... trim=False)"
        )

    @staticmethod
    def uniform(sequence_length, rate):
        """
        Create a uniform rate map
        """
        return RateMap([0, sequence_length], [rate])

    def asdict(self):
        return {"position": self.position, "rate": self.rate}

    @property
    def sequence_length(self):
        """
        The sequence length covered by this map
        """
        return self.position[-1]

    @property
    def total_mass(self):
        """
        The cumulative total over the entire map
        """
        # Since we insist that the rate up to start_position and past end_position is 0
        # we can just return the cumulative total.
        return self.cumulative_mass[-1]

    @property
    def size(self):  # TODO - should we remove this, for the same reason as for __len__
        return self.rate.shape[0]

    @property
    def mean_rate(self):
        """
        The mean rate over this map, from :attr:`start_position` to
        :attr:`end_position`, weighted by the span covered by each rate.
        """
        span = self.end_position - self.start_position
        return self.total_mass / span

    def get_rate(self, x):
        """
        Return the rate at a list of user-specified positions along the map, i.e.
        interpolated between the fixed :attr:`position` values. Any positions
        beyond the end  of the map will return the final rate value.

        :param numpy.ndarray x: The positions for which to return values.

        :return: An array of rates, the same length as ``x``.
        :rtype: numpy.ndarray(dtype=np.float64)
        """
        loc = np.searchsorted(self.position, x, side="right") - 1
        loc = np.minimum(loc, len(self.rate) - 1)  # Return final rate if off the end
        if np.any(loc < 0):
            raise ValueError("No rate exists for positions below the minimum map pos")
        return self.rate[loc]

    def get_cumulative_mass(self, x):
        """
        Return the cumulative rate at a list of user-specified positions along the map,
        i.e. interpolated between the fixed :attr:`position` values.

        :param numpy.ndarray x: The positions for which to return values.

        :return: An array of cumulative rates, the same length as ``x``
        :rtype: numpy.ndarray(dtype=np.float64)
        """
        if np.any(np.array(x) < 0) or np.any(np.array(x) > self.sequence_length):
            raise ValueError(
                f"Cannot have physical positions < 0 or > {self.sequence_length}"
            )
        return np.interp(x, self.position, self.cumulative_mass)

    def slice(self, start=None, end=None, *, trim=False):  # noqa:A003
        """
        Returns a subset of this rate map between the specified end points.

        :param float start: The start of the region to keep. If ``None``, it defaults to
            0.
        :param float end: The end of the region to keep. If end is ``None``, it defaults
            to the end of the map.
        :param bool trim: If True, remove the flanking regions such that the
            sequence length of the new rate map is ``end`` - ``start``. If ``False``
            (default), do not change the coordinate system, but instead replace each
            flanking region with 0, and set the ``start_position`` and ``end_position``
            of the returned map to ``start`` and ``end``, inserting these positions if
            they do not exist.

        :return: A new RateMap instance
        :rtype: RateMap
        """
        if start is None:
            i = 0
            start = 0
        if end is None:
            end = self.sequence_length
            j = len(self.position)

        if (
            start < 0
            or end < 0
            or start > self.sequence_length
            or end > self.sequence_length
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

        if trim:
            # Return trimmed map with changed coords
            return self.__class__(position - start, rate)

        # Need to check regions before & after sliced region are filled out:
        if start != 0:
            if rate[0] == 0:
                position[0] = 0  # Extend
            else:
                rate = np.insert(rate, 0, 0)  # Prepend
                position = np.insert(position, 0, 0)
        if end != self.position[-1]:
            if rate[-1] == 0:
                position[-1] = self.sequence_length  # Extend
            else:
                rate = np.append(rate, 0)  # Append
                position = np.append(position, self.position[-1])
        return self.__class__(position, rate, start_position=start, end_position=end)

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
        return cls(rate_map.position, rate, map_start=rate_map.start_position)

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
        return self.map.cumulative_mass[-1]

    def physical_to_genetic(self, x):
        return self.map.get_cumulative_mass(x)

    def genetic_to_physical(self, genetic_x):
        if self.map.cumulative_mass[-1] == 0:
            # If we have a zero recombination rate throughout then everything
            # except L maps to 0.
            return self.get_sequence_length() if genetic_x > 0 else 0
        if genetic_x == 0:
            return self.map.position[0]
        index = np.searchsorted(self.map.cumulative_mass, genetic_x) - 1
        y = (
            self.map.position[index]
            + (genetic_x - self.map.cumulative_mass[index]) / self.map.rate[index]
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
    return RateMap(position, rate[:-1], start_position=map_start)
