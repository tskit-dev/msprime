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

    :param list position: A list of :math:`n` positions, starting at 0, and ending in the
        sequence length over which the RateMap will apply.
    :param list rate: A list of :math:`n-1` positive rates that apply between each
        position.
    :param list cumulative_rate: A alternative method to specify rates between each
        position: a list of :math:`n` increasing values starting at 0 providing the
        cumulative rate.
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
        rate=None,
        cumulative_rate=None,
        *,
        start_position=None,
        end_position=None,
    ):
        # Taking copies ensures there are no unintended consequences later if a user
        # modifies the input arrays.

        # Save positions into a read-only array
        self._position = np.array(position, dtype=float, copy=True)
        size = len(self.position)
        if size < 2:
            raise ValueError("Must have at least two positions")
        if self.position[0] != 0:
            raise ValueError("First position must be zero")
        if np.any(self.position[:1] >= self.position[1:]):
            raise ValueError("Position values must be in increasing order")
        self._position.flags.writeable = False

        # Save rates / cumulative rates into 2 read-only arrays
        if rate is None and cumulative_rate is None:
            raise TypeError("One of 'rate' or 'cumulative_rate' must be specified")
        if rate is not None and cumulative_rate is not None:
            raise TypeError("Cannot specify both 'rate' and 'cumulative_rate'")

        if rate is not None:
            self._rate = np.array(rate, dtype=float, copy=True)
            if self.rate.shape[0] != size - 1:
                raise ValueError(
                    "Rate array must have one less entry than the position array"
                )
            if np.any(self.rate < 0):
                raise ValueError("Rates must be non-negative")
            self._cumulative = np.insert(
                np.cumsum(np.diff(self.position) * self.rate), 0, 0
            )
        else:
            self._cumulative = np.array(cumulative_rate, dtype=float, copy=True)
            if self.cumulative.shape[0] != size:
                raise ValueError(
                    "Cumulative rate array must be the same size as the position array"
                )
            if self.cumulative[0] != 0:
                raise ValueError("Cumulative rate array must start with 0")
            self._rate = np.diff(self.cumulative) / np.diff(self.position)
            if np.any(self.rate < 0):
                raise ValueError("Cumulative rates must be strictly increasing")

        if start_position is None:
            start_position = 0
        elif self.cumulative_at(start_position) > 0:
            raise ValueError("Rates before the start_position must all be zero")
        self.start_position = start_position

        if end_position is None:
            end_position = self.sequence_length
        elif self.cumulative_at(end_position) < self.total_mass:
            raise ValueError("Rates after the end_position must all be zero")
        self.end_position = end_position

        # Make all of the internal arrays read-only, so they can't get out of sync
        self._rate.flags.writeable = False
        self._cumulative.flags.writeable = False

    @property
    def position(self):
        return self._position

    @property
    def rate(self):
        return self._rate

    @property
    def cumulative(self):
        return self._cumulative

    def __len__(self):
        return len(self._rate)

    @staticmethod
    def uniform(sequence_length, rate):
        return RateMap([0, sequence_length], [rate])

    def asdict(self):
        return {"position": self._position, "rate": self._rate}

    @property
    def sequence_length(self):
        return self.position[-1]

    @property
    def total_mass(self):
        """
        The cumulative total over the entire map
        """
        # Since we insist that the rate up to start_position and past end_position is 0
        # we can just return the cumulative total.
        return self.cumulative[-1]

    @property
    def size(self):
        return self.rate.shape[0]

    @property
    def mean_rate(self):
        """
        The mean rate over this map, from map.start_position to map.end_position.
        """
        span = self.end_position - self.start_position
        return self.total_mass / span

    def at(self, x):
        """
        Return the rate at a list of specified positions along the map. Any positions
        beyond the end of the map will return the final rate value.

        :param numpy.ndarray x: The positions for which to return values.

        :return: An array of rates, the same length as ``x``.
        :rtype: numpy.ndarray(dtype=np.float64)
        """
        loc = np.searchsorted(self.position, x, side="right") - 1
        if np.any(loc < 0):
            raise ValueError("No rate exists for positions below the minimum map pos")
        return self.rate[np.minimum(loc, len(self.rate) - 1)]

    def cumulative_at(self, x):
        """
        Return the cumulative rate at a list of specified positions along the map.

        :param numpy.ndarray x: The positions for which to return values.

        :return: An array of cumulative rates, the same length as ``x``
        :rtype: numpy.ndarray(dtype=np.float64)
        """
        if np.any(np.array(x) < 0) or np.any(np.array(x) > self.sequence_length):
            raise ValueError(
                f"Cannot have physical positions < 0 or > {self.sequence_length}"
            )
        return np.interp(x, self.position, self.cumulative)

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
            of the returned map to ``start`` and ``end``.

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
    :param str name: A name to use for this map, and stored in the ``.name``
        attribute. This will only be used if the map is written out to a file
    """

    def __init__(self, positions, rates, num_loci=None, name=None, **kwargs):
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
        self.map = RateMap(positions, rates[:-1], **kwargs)
        self.name = name

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
    def read_hapmap(cls, filename, sequence_length=None):
        """
        Parses the specified file in HapMap format. These files must be
        white-space-delimited, and contain a single header line (which is
        ignored) and a set of rows that pertain to a single chromosome.
        Each line after the header is expect to contains a name (e.g. chr1),
        and the starting position and recombination rate for the segment from
        that position (inclusive) to the starting position on the next line
        (exclusive). Starting positions of each segment are given in units of
        bases, and recombination rates in centimorgans/Megabase. A warning is
        given if the name in the first column differs between rows. Otherwise
        only the second (Position) and third (Rate) columns are read: any
        further columns are ignored. If the first starting position is not equal
        to zero, then a zero-recombination region is inserted at the start of
        the chromosome.

        .. warning::
            This method is deprecated, use the module-level
            :func:`read_hapmap` function instead, which reads the genetic position
            in centiMorgans from the file, rather than using the rate column.

        :param str filename: The name of the file to be parsed. This may be
            in plain text or gzipped plain text. It may also be a file object.
        :param str sequence_length: The total length of the map. If ``None``, then
            assume it is the final physical position in the file. Otherwise pad the
            map with a zero rate region between the final position in the file and
            the value specified here.
        :return: A RecombinationMap object.
        """
        warnings.warn(
            "RecombinationMap.read_hapmap() is deprecated. "
            "Use msprime.read_hapmap() instead.",
            FutureWarning,
        )
        rownames, position, rate = np.loadtxt(
            filename,
            skiprows=1,
            dtype=[("chrom", object), ("position", int), ("rate", float)],
            usecols=(0, 1, 2),
            unpack=True,
        )
        if len(rownames) == 0:
            raise ValueError("Empty hapmap file")

        if len(np.unique(rownames)) > 1:
            warnings.warn(
                f"Different chromosome names found in the first column of {filename}. "
                "Check that the file contains data for a single chromosome."
            )
        name = rownames[0]
        # Rate is expressed in centimorgans per megabase: convert to per-base rates
        rate *= 1e-8

        start = position[0]
        end = position[-1]
        if start != 0:
            position = np.insert(position, 0, 0)
            rate = np.insert(rate, 0, 0)
        if rate[-1] != 0:
            raise ValueError(
                "The last rate provided in the recombination map must be zero"
            )
        if sequence_length is not None:
            if sequence_length < end:
                raise ValueError(
                    f"Cannot specify a sequence_length < {end} (the last position)"
                )
            elif sequence_length > end:
                position = np.append(position, sequence_length)
                rate = np.append(rate, 0)
        return cls(position, rate, name=name, start_position=start, end_position=end)

    def write_hapmap(self, filename):
        """
        Write out the recombination map in a format with genetic distances in cM,
        so that it can be read in by the non-deprecated read-hapmap function.
        Omit the rate up to the start_position and after the end_position.
        """
        headers = ["Chromosome", "Position(bp)", "Rate(cM/Mb)", "Map(cM)"]
        name = "Unknown" if self.name is None else self.name
        with open(filename, "wt") as file:
            file.write("\t".join(headers) + "\n")
            for pos, rate, cumul in zip(
                self.position, self.get_rates(), self.cumulative
            ):
                if self.map.start_position <= pos <= self.map.end_position:
                    file.write(
                        "{}\t{:.12g}\t{:g}\t{:g}\n".format(
                            name, pos, rate * 1e8, cumul * 100
                        )
                    )

    @property
    def position(self):
        # For backwards compatibility
        return self.map.position

    @property
    def rate(self):
        # For backwards compatibility
        return self.map.rate

    @property
    def cumulative(self):
        return self.map.cumulative

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

    def at(self, x):
        """
        Return the recombination rates at a set of positions along the chromosome

        :param numpy.ndarray x: The physical positions along the chromosome.

        :return: An array of rates, the same length as ``x``
        :rtype: numpy.ndarray(dtype=np.float64)
        """
        return self.map.at(x)

    def physical_to_genetic(self, x):
        """
        Convert physical positions into genetic coordinates

        :param numpy.ndarray x: The physical positions along the chromosome.

        :return: An array of genetic positions, the same length as ``x``
        :rtype: numpy.ndarray(dtype=np.float64)
        """
        return self.map.cumulative_at(x)

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

    def __len__(self):
        return len(self.map.rate)

    def asdict(self):
        return self.map.asdict()


def read_hapmap(
    filename, *, sequence_length=None, has_header=True, bp_column=1, cM_column=3
):
    # Black barfs with an INTERNAL_ERROR trying to reformat this docstring,
    # so we explicitly disable reformatting here.
    # fmt: off
    """
    Parses the specified file in HapMap format. These files must be
    white-space-delimited, and by default are assumed to contain a single header line
    (which is ignored). Each subsequent line then contains a physical position (in base
    pairs) and a corresponding genetic map position (in centiMorgans). The rate between
    the current physical position (inclusive) and the physical position on the next line
    (exclusive) is taken as constant. Only these two values are used from each line:
    by default they are assumed to occupy the second and fourth columns of the file
    respectively, as exemplified in the following sample of the format::

        Chromosome	Position(bp)  Rate(cM/Mb)  Map(cM)
        chr10       48232         0.1614       0.002664
        chr10       48486         0.1589       0.002705
        chr10       50009         0.159        0.002947
        chr10       52147         0.1574       0.003287
        ...
        chr10	    133762002     3.358        181.129345
        chr10	    133766368     0.000        181.144008

    .. note::
        Positions in a RateMap are zero-based, meaning that rate values in the returned
        RateMap object are defined from position 0 onwards. Thus in the example above,
        a rate of 0.002664/48232 cM/bp (i.e. 0.0552 cM/Mb) is assumed to apply from
        position 0 (inclusive) to 48232 (exclusive). If the map positions use a
        one-based, system, with the first position being 1 rather than 0, this may
        result in an miniscule difference in the first rate value calculated.

    :param str filename: The name of the file to be parsed. This may be
        in plain text or gzipped plain text.
    :param str sequence_length: The total length of the map. If ``None``, then assume
        it is the final physical position in the file. Otherwise pad the map with a zero
        rate region between the final position in the file and the value specified here.
    :param str has_header: If True, assume the file has a header row and ignore the
        first line of the file.
    :param str bp_column: The zero-based index of the column in the file specifying
        the physical position in base pairs (default: 1, i.e. the second column)
    :param str cM_column: The zero-based index of the column in the file specifying
        the genetic map position in centiMorgans (default: 3, i.e. the fourth column)
    :return: A RateMap object.
    """
    # fmt: on
    hapmap = np.loadtxt(
        filename, skiprows=1 if has_header else 0, usecols=(bp_column, cM_column)
    )
    physical_positions = hapmap[:, 0]
    genetic_positions = hapmap[:, 1] / 100  # Convert centiMorgans to Morgans

    if np.any(physical_positions < 0) or np.any(genetic_positions < 0):
        bad_lines = np.unique(
            np.concatenate(
                (
                    np.where(physical_positions < 0)[0],
                    np.where(genetic_positions < 0)[0],
                )
            )
        )
        raise ValueError(f"The following lines have positions < 0: {bad_lines}")

    start = physical_positions[0]
    end = physical_positions[-1]

    if genetic_positions[0] > 0 and start == 0:
        raise ValueError("The map distance at the start of the chromosome must be zero")
    if start > 0:
        physical_positions = np.insert(physical_positions, 0, 0)
        if genetic_positions[0] > 0:
            # Exception for a map that starts > 0cM: include the start rate in the mean
            start = 0
        genetic_positions = np.insert(genetic_positions, 0, 0)

    if sequence_length is not None:
        if sequence_length < end:
            raise ValueError(
                "The sequence_length cannot be less that the last physical position "
                f" ({physical_positions[-1]} cM)"
            )
        if sequence_length > end:
            physical_positions = np.append(physical_positions, sequence_length)
            genetic_positions = np.append(genetic_positions, genetic_positions[-1])

    return RateMap(
        physical_positions,
        cumulative_rate=genetic_positions,
        start_position=start,
        end_position=end,
    )
