#
# Copyright (C) 2020-2021 University of Oxford
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
from __future__ import annotations

import collections.abc
import numbers
import warnings

import numpy as np

from msprime import core


class RateMap(collections.abc.Mapping):
    """
    A class mapping a non-negative rate value to a set of adjacent intervals
    along the genome.

    :param list position: A list of :math:`n+1` positions, starting at 0, and ending
        in the sequence length over which the RateMap will apply.
    :param list rate: A list of :math:`n` positive rates that apply between each
        position.
    """

    # The args are marked keyword only to give us some flexibility in how we
    # create class this in the future.
    def __init__(
        self,
        *,
        position,
        rate,
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
        # TODO make these public and think about the names
        self._unknown = np.isnan(self.rate)
        self._known = ~self._unknown
        self._left = self._position[:-1]
        self._right = self.position[1:]
        mass = self.span[self._known] * self.rate[self._known]
        self._cumulative_mass = np.insert(np.cumsum(mass), 0, 0)
        self._cumulative_mass.flags.writeable = False

    @property
    def position(self):
        """
        The positions between which a constant rate applies. There will be
        one more position than the number of rates, with the first position being 0
        and the last being the :attr:`sequence_length`.
        """
        return self._position

    @property
    def left(self):
        """
        The left position of each interval (inclusive) the map is defined
        over.
        """
        return self._left

    @property
    def right(self):
        """
        The right position of each interval (exclusive) the map is defined
        over.
        """
        return self._right

    @property
    def mid(self):
        """
        Returns the midpoint of each of the intervals.
        """
        mid = self.left + self.span / 2
        mid.flags.writeable = False
        return mid

    @property
    def span(self):
        """
        Returns the span (i.e., ``right - left``) of each of the intervals.
        """
        span = self.right - self.left
        span.flags.writeable = False
        return span

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

    @staticmethod
    def uniform(sequence_length, rate) -> RateMap:
        """
        Create a uniform rate map
        """
        return RateMap(position=[0, sequence_length], rate=[rate])

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
        The cumulative total over the entire map.
        """
        return self.cumulative_mass[-1]

    @property
    def mean_rate(self):
        """
        The mean rate over this map weighted by the span covered by each rate.
        Unknown intervals are excluded.
        """
        total_span = np.sum(self.span[self._known])
        return self.total_mass / total_span

    def get_rate(self, x):
        """
        Return the rate at the specified list of positions.

        :param numpy.ndarray x: The positions for which to return values.
        :return: An array of rates, the same length as ``x``.
        :rtype: numpy.ndarray
        """
        loc = np.searchsorted(self.position, x, side="right") - 1
        if np.any(loc < 0) or np.any(loc >= len(self.rate)):
            raise ValueError("position out of bounds")
        return self.rate[loc]

    def get_cumulative_mass(self, x):
        """
        Return the cumulative mass at a list of positions along the map.

        :param numpy.ndarray x: The positions for which to return values.

        :return: An array of cumulative rates, the same length as ``x``
        :rtype: numpy.ndarray
        """
        x = np.array(x)
        if np.any(x < 0) or np.any(x > self.sequence_length):
            raise ValueError(f"Cannot have positions < 0 or > {self.sequence_length}")
        return np.interp(x, self.position, self.cumulative_mass)

    def find_index(self, x: float) -> int:
        """
        Returns the index of the interval that the specified position falls within,
        such that ``rate_map.left[index] <= x < self.rate_map.right[index]``.

        :param float x: The position to search.
        :return: The index of the interval containing this point.
        :rtype: int
        :raises: KeyError if the position is not contained in any of the intervals.
        """
        if x < 0 or x >= self.sequence_length:
            raise KeyError(f"Position {x} out of bounds")
        index = np.searchsorted(self.position, x, side="left")
        if x < self.position[index]:
            index -= 1
        assert self.left[index] <= x < self.right[index]
        return index

    def copy(self) -> RateMap:
        """
        Returns a deep copy of this RateMap.
        """
        return RateMap(position=self.position.copy(), rate=self.rate.copy())

    def slice(self, left=None, right=None, *, trim=False) -> RateMap:  # noqa: A003
        """
        Returns a subset of this rate map in the specified interval.

        :param float left: The left coordinate (inclusive) of the region to keep.
            If ``None``, defaults to 0.
        :param float right: The right coordinate (exclusive) of the region to keep.
            If ``None``, defaults to the sequence length.
        :param bool trim: If True, remove the flanking regions such that the
            sequence length of the new rate map is ``right`` - ``left``. If ``False``
            (default), do not change the coordinate system and mark the flanking
            regions as "unknown".
        :return: A new RateMap instance
        :rtype: RateMap
        """
        left = 0 if left is None else left
        right = self.sequence_length if right is None else right
        if not (0 <= left < right <= self.sequence_length):
            raise KeyError(f"Invalid slice: left={left}, right={right}")

        i = self.find_index(left)
        j = i + np.searchsorted(self.position[i:], right, side="right")
        if right > self.position[j - 1]:
            j += 1

        position = self.position[i:j].copy()
        rate = self.rate[i : j - 1].copy()
        position[0] = left
        position[-1] = right

        if trim:
            # Return trimmed map with changed coords
            return RateMap(position=position - left, rate=rate)

        # Need to check regions before & after sliced region are filled out:
        if left != 0:
            if np.isnan(rate[0]):
                position[0] = 0  # Extend
            else:
                rate = np.insert(rate, 0, np.nan)  # Prepend
                position = np.insert(position, 0, 0)
        if right != self.position[-1]:
            if np.isnan(rate[-1]):
                position[-1] = self.sequence_length  # Extend
            else:
                rate = np.append(rate, np.nan)  # Append
                position = np.append(position, self.position[-1])
        return RateMap(position=position, rate=rate)

    def unknown_intervals(self):
        num_unknown = np.sum(self._unknown)
        out = np.empty((num_unknown, 2))
        out[:, 0] = self.left[self._unknown]
        out[:, 1] = self.right[self._unknown]
        return out

    def __iter__(self):
        # The clinching argument for using mid here is that if we used
        # left instead we would have
        #   RateMap([0, 1], [0.1]) == RateMap([0, 100], [0.1])
        # by the inherited definition of equality since the dictionary items
        # would be equal.
        # Similarly, we only return the midpoints of known intervals
        # because NaN values are not equal, and we would need to do
        # something to work around this. It seems reasonable that
        # this high-level operation returns the *known* values only
        # anyway.
        yield from self.mid[self._known]

    def __len__(self):
        return np.sum(self._known)

    def __getitem__(self, key):
        if isinstance(key, slice):
            if key.step is not None:
                raise TypeError("Only interval slicing is supported")
            return self.slice(key.start, key.stop)
        if isinstance(key, numbers.Number):
            index = self.find_index(key)
            return self.rate[index]
        raise KeyError("Key {key} not in map")

    def _display_table(self):
        data = [
            [
                f"{left:.10g}",
                f"{right:.10g}",
                f"{mid:.10g}",
                f"{span:.10g}",
                f"{rate:.2f}",
            ]
            for left, right, mid, span, rate in zip(
                self.left, self.right, self.mid, self.span, self.rate
            )
        ]
        return ["left", "right", "mid", "span", "rate"], data

    def __str__(self):
        titles, data = self._display_table()
        data = [[[item] for item in row] for row in data]
        table = core.text_table(
            caption="",
            column_titles=[[title] for title in titles],
            column_alignments="<<>>>",
            data=data,
        )
        return table

    def _repr_html_(self):
        col_titles, data = self._display_table()
        return core.html_table("", col_titles, data)

    def __repr__(self):
        return f"RateMap(position={repr(self.position)}, rate={repr(self.rate)})"

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
            interval from 0 to first physical position in the file is marked as
            unknown, and the last value in the rate column must be zero.
        :param int map_col: The zero-based index of the column in the file
            specifying the genetic map position in centiMorgans. If ``None``
            (default), assume an index of 3 (i.e. the fourth column). If the first
            genetic position is 0 the interval from position 0 to the first
            physical position in the file is marked as unknown. Otherwise, act
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
        rate = np.diff(genetic_positions) / np.diff(physical_positions)
        if start != 0:
            rate[0] = np.nan
        if end != physical_positions[-1]:
            rate[-1] = np.nan
        return RateMap(position=physical_positions, rate=rate)


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

    .. important::
        This class is deprecated (but supported indefinitely);
        please use the :class:`.RateMap` class in new code.
        In particular, note that when specifying ``rates`` in the
        the :class:`.RateMap` class we now require an array
        of length :math:`n - 1` (this class requires an array
        of length :math:`n` in which the last entry is zero).

    :param list positions: The positions (in bases) denoting the
        distinct intervals where recombination rates change. These can
        be floating point values.
    :param list rates: The list of rates corresponding to the supplied
        ``positions``. Recombination rates are specified per base,
        per generation.
    :param int num_loci: **This parameter is no longer supported.**
        Must be either None (meaning a continuous genome of the
        finest possible resolution) or be equal to ``positions[-1]``
        (meaning a discrete genome). Any other value will result in
        an error. Please see the :ref:`sec_legacy_0x_genome_discretisation`
        section for more information.
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
        self.map = RateMap(position=positions, rate=rates[:-1])

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
            This method is deprecated, use the :meth:`.RateMap.read_hapmap`
            method instead.

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
        # Mark anything unknown as 0 for backwards compatibility. This will
        # ensure that simulate() never trims parts of the tree sequence.
        rate = rate_map.rate.copy()
        rate[rate_map._unknown] = 0
        return cls(rate_map.position, np.append(rate, 0))

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
