#
# Copyright (C) 2016 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
Module responsible for converting tree sequence files from older
formats.
"""
from __future__ import division
from __future__ import print_function


try:
    import h5py
    _h5py_imported = True
except ImportError:
    _h5py_imported = False


import msprime


class Hdf5FileReader(object):
    """
    Class to read msprime tree sequence HDF5 files.
    """
    def __init__(self, filename):
        if not _h5py_imported:
            raise RuntimeError("h5py is required for converting HDF5 files.")
        self._filename = filename
        self._root = h5py.File(self._filename, "r")
        if 'format_version' not in self._root.attrs:
            raise ValueError("Cannot read HDF5 file {}".format(self._filename))
        format_version = self._root.attrs['format_version']
        if format_version[0] != 2:
            raise ValueError("Only version 2.x are supported")

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        self.close()

    def close(self):
        """
        Closes this Hdf5FileReader.
        """
        self._root.close()
        self._root = None

    def records(self):
        """
        Returns an iterator over the coalescence records in the file.
        """
        trees = self._root["trees"]
        left = trees["left"]
        right = trees["right"]
        node = trees["node"]
        children = trees["children"]
        population = trees["population"]
        time = trees["time"]
        num_records = len(left)
        for j in range(num_records):
            yield msprime.CoalescenceRecord(
                left=left[j], right=right[j], node=node[j],
                children=tuple(children[j]), time=time[j],
                population=population[j])

    def mutations(self):
        """
        Returns an iterator over the mutations in the file.
        """
        if "mutations" in self._root:
            mutations = self._root["mutations"]
            position = mutations["position"]
            node = mutations["node"]
            num_mutations = len(node)
            for j in range(num_mutations):
                yield msprime.Mutation(position=position[j], node=node[j])

    def samples(self):
        if "samples" in self._root:
            samples = self._root["samples"]
            population = samples["population"]
            time = None
            if "time" in samples:
                time = samples["time"]
            num_samples = len(samples)
            for j in range(num_samples):
                t = 0
                if time is not None:
                    t = time[j]
                yield msprime.Sample(population=population[j], time=t)
