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

    @staticmethod
    def read_hapmap(
        fileobj,
        sequence_length=None,
        *,
        has_header=True,
        position_col=None,
        rate_col=None,
        map_col=None,
    ):
        # Black barfs with an INTERNAL_ERROR trying to reformat this docstring,
        # so we explicitly disable reformatting here.
        # fmt: off
        """
        Parses the specified file in HapMap format and returns a RateMap.
        HapMap files must white-space-delimited, and by default are assumed to
        contain a single header line (which is ignored). Each subsequent line
        then contains a physical position (in base pairs) and either a genetic
        map position (in centiMorgans) or a recombination rate (in centiMorgans
        per megabase); the rate between the current physical position
        (inclusive) and the physical position on the next line (exclusive) is
        taken as constant. By default, the second column of the file is taken
        as the physical position and the fourth column is taken as the genetic
        position, as seen in the following sample of the format::

            Chromosome	Position(bp)  Rate(cM/Mb)  Map(cM)
            chr10       48232         0.1614       0.002664
            chr10       48486         0.1589       0.002705
            chr10       50009         0.159        0.002947
            chr10       52147         0.1574       0.003287
            ...
            chr10	133762002     3.358        181.129345
            chr10	133766368     0.000        181.144008

        .. note::
            The rows are all assumed to come from the same contig, and the
            first column is currently ignored. Therefore if you have a single
            file containing several contigs or chromosomes, you must must split
            it up into multiple files, and pass each one separately to this
            function.

        :param str fileobj: Filename or file to read. This is passed directly
            to :func:`numpy.loadtxt`, so if the filename extension is .gz or .bz2,
            the file is decompressed first
        :param float sequence_length: The total length of the map. If ``None``,
            then assume it is the last physical position listed in the file.
            Otherwise it must be greater then or equal to the last physical
            position in the file, and the region between the last physical position
            and the sequence_length is padded with a rate of zero.
        :param bool has_header: If True (default), assume the file has a header row
            and ignore the first line of the file.
        :param int position_col: The zero-based index of the column in the file
            specifying the physical position in base pairs. If ``None`` (default)
            assume an index of 1 (i.e. the second column).
        :param int rate_col: The zero-based index of the column in the file
            specifying the rate in cM/Mb. If ``None`` (default) do not use the rate
            column, but calculate rates using the genetic map positions, as
            specified in ``map_col``. If the rate column is used, the
            :attr:`RateMap.start_position` of the returned map is set to the first
            physical position in the file, and the last value in the rate column
            must be zero.
        :param int map_col: The zero-based index of the column in the file
            specifying the genetic map position in centiMorgans. If ``None``
            (default), assume an index of 3 (i.e. the fourth column). If the first
            genetic position is 0, set the :attr:`RateMap.start_position` of the
            returned map to the first physical position in the file. Otherwise, act
            as if an additional row, specifying physical position 0 and genetic
            position 0, exists at the start of the file.
        :return: A RateMap object.
        :rtype: RateMap
        """
        # fmt: on
        column_defs = {}  # column definitions passed to np.loadtxt
        if rate_col is None and map_col is None:
            # Default to map_col
            map_col = 3
        elif rate_col is not None and map_col is not None:
            raise ValueError("Cannot specify both rate_col and map_col")
        if map_col is not None:
            column_defs[map_col] = ("map", float)
        else:
            column_defs[rate_col] = ("rate", float)
        position_col = 1 if position_col is None else position_col
        if position_col in column_defs:
            raise ValueError(
                "Cannot specify the same columns for position_col and "
                "rate_col or map_col"
            )
        column_defs[position_col] = ("pos", int)

        column_names = [c[0] for c in column_defs.values()]
        column_data = np.loadtxt(
            fileobj,
            skiprows=1 if has_header else 0,
            dtype=list(column_defs.values()),
            usecols=list(column_defs.keys()),
            unpack=True,
        )
        data = dict(zip(column_names, column_data))

        if "map" not in data:
            assert "rate" in data
            if data["rate"][-1] != 0:
                raise ValueError("The last entry in the 'rate' column must be zero")
            pos_Mb = data["pos"] / 1e6
            map_pos = np.cumsum(data["rate"][:-1] * np.diff(pos_Mb))
            data["map"] = np.insert(map_pos, 0, 0) / 100
        else:
            data["map"] /= 100  # Convert centiMorgans to Morgans
        if len(data["map"]) == 0:
            raise ValueError("Empty hapmap file")

        # TO DO: read in chrom name from col 0 and poss set as .name
        # attribute on the RateMap

        physical_positions = data["pos"]
        genetic_positions = data["map"]
        start = physical_positions[0]
        end = physical_positions[-1]

        if genetic_positions[0] > 0 and start == 0:
            raise ValueError(
                "The map distance at the start of the chromosome must be zero"
            )
        if start > 0:
            physical_positions = np.insert(physical_positions, 0, 0)
            if genetic_positions[0] > 0:
                # Exception for a map that starts > 0cM: include the start rate
                # in the mean
                start = 0
            genetic_positions = np.insert(genetic_positions, 0, 0)

        if sequence_length is not None:
            if sequence_length < end:
                raise ValueError(
                    "The sequence_length cannot be less that the last physical position "
                    f" ({physical_positions[-1]})"
                )
            if sequence_length > end:
                physical_positions = np.append(physical_positions, sequence_length)
                genetic_positions = np.append(genetic_positions, genetic_positions[-1])

        assert genetic_positions[0] == 0
        return RateMap(
            physical_positions,
            rate=np.diff(genetic_positions) / np.diff(physical_positions),
            start_position=start,
            end_position=end,
        )


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
        # Used as an internal flag for the 0.x simulate() function. This allows
        # us to emulate the discrete-sites behaviour of 0.x code.
        self._is_discrete = num_loci == positions[-1]
        if num_loci is not None and num_loci != positions[-1]:
            raise ValueError(
                "The RecombinationMap interface is deprecated and only "
                "partially supported. If you wish to simulate a number of "
                "discrete loci, you must set num_loci == the sequence length. "
                "If you wish to simulate recombination process on as fine "
                "a map as possible, please omit the num_loci parameter (or set "
                "to None). Otherwise, num_loci is no longer supported and "
                "the behaviour of msprime 0.x cannot be emulated. Please "
                "consider upgrading your code to the version 1.x APIs."
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
            "Use RateMap.read_hapmap() instead.",
            FutureWarning,
        )
        rate_map = RateMap.read_hapmap(filename, position_col=1, rate_col=2)
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
        return len(self.map.position)

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
