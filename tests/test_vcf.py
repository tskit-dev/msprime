#
# Copyright (C) 2016 University of Oxford
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
Test cases for VCF output in msprime.
"""
from __future__ import print_function
from __future__ import division

import collections
import math
import os
import tempfile
import unittest

import msprime

import vcf
# Pysam is not available on windows, so we don't make it mandatory here.
_pysam_imported = False
try:
    import pysam
    _pysam_imported = True
except ImportError:
    pass


test_data = []


def setUp():
    Datum = collections.namedtuple(
        "Datum",
        ["tree_sequence", "ploidy", "contig_id", "vcf_file", "sample_names"])
    L = 100
    for ploidy in [1, 2, 3, 5]:
        for contig_id in ["1", "x" * 8]:
            for n in [2, 10]:
                for rho in [0, 0.5]:
                    for mu in [0, 1.0]:
                        ts = msprime.simulate(
                            n * ploidy, length=L, recombination_rate=rho,
                            mutation_rate=mu)
                        fd, file_name = tempfile.mkstemp(prefix="msprime_vcf_")
                        os.close(fd)
                        with open(file_name, "w") as f:
                            ts.write_vcf(f, ploidy, contig_id)
                        sample_names = ["msp_{}".format(j) for j in range(n)]
                        test_data.append(
                            Datum(ts, ploidy, contig_id, file_name, sample_names))


def tearDown():
    for datum in test_data:
        os.unlink(datum.vcf_file)


def write_vcf(tree_sequence, output, ploidy, contig_id):
    """
    Writes a VCF using the sample algorithm as the low level code.
    """
    if tree_sequence.get_sample_size() % ploidy != 0:
        raise ValueError("Sample size must a multiple of ploidy")
    n = tree_sequence.get_sample_size() // ploidy
    sample_names = ["msp_{}".format(j) for j in range(n)]
    last_pos = 0
    positions = []
    for variant in tree_sequence.variants():
        pos = int(round(variant.position))
        if pos <= last_pos:
            pos = last_pos + 1
        positions.append(pos)
        last_pos = pos
    contig_length = int(math.ceil(tree_sequence.get_sequence_length()))
    if len(positions) > 0:
        contig_length = max(positions[-1], contig_length)
    print("##fileformat=VCFv4.2", file=output)
    print("##source=msprime {}".format(msprime.__version__), file=output)
    print(
        '##FILTER=<ID=PASS,Description="All filters passed">',
        file=output)
    print("##contig=<ID={},length={}>".format(contig_id, contig_length), file=output)
    print(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        file=output)
    print(
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "FORMAT", sep="\t", end="", file=output)
    for sample_name in sample_names:
        print("\t", sample_name, sep="", end="", file=output)
    print(file=output)
    for variant in tree_sequence.variants():
        pos = positions[variant.index]
        print(
            contig_id, pos, ".", "A", "T", ".", "PASS", ".", "GT",
            sep="\t", end="", file=output)
        for j in range(n):
            genotype = "|".join(
                str(g) for g in
                variant.genotypes[j * ploidy: j * ploidy + ploidy])
            print("\t", genotype, end="", sep="", file=output)
        print(file=output)


class TestEquality(unittest.TestCase):
    """
    Tests if the VCF file produced by the low level code is the
    same as one we generate here.
    """

    def test_equal(self):
        for datum in test_data:
            with tempfile.TemporaryFile("w+") as f:
                write_vcf(datum.tree_sequence, f, datum.ploidy, datum.contig_id)
                f.seek(0)
                vcf1 = f.read()
            with open(datum.vcf_file) as f:
                vcf2 = f.read()
            self.assertEqual(vcf1, vcf2)


class TestHeaderParsers(unittest.TestCase):
    """
    Tests if we can parse the headers with various tools.
    """
    def test_pyvcf(self):
        for datum in test_data:
            reader = vcf.Reader(filename=datum.vcf_file)
            self.assertEqual(len(reader.contigs), 1)
            contig = reader.contigs[datum.contig_id]
            self.assertEqual(contig.id, datum.contig_id)
            self.assertGreater(contig.length, 0)
            self.assertEqual(len(reader.alts), 0)
            self.assertEqual(len(reader.filters), 1)
            p = reader.filters["PASS"]
            self.assertEqual(p.id, "PASS")
            self.assertEqual(len(reader.formats), 1)
            f = reader.formats["GT"]
            self.assertEqual(f.id, "GT")
            self.assertEqual(len(reader.infos), 0)

    @unittest.skipIf(not _pysam_imported, "pysam not available")
    def test_pysam(self):
        for datum in test_data:
            bcf_file = pysam.VariantFile(datum.vcf_file)
            self.assertEqual(bcf_file.format, "VCF")
            self.assertEqual(bcf_file.version, (4, 2))
            header = bcf_file.header
            self.assertEqual(len(header.contigs), 1)
            contig = header.contigs[0]
            self.assertEqual(contig.name, datum.contig_id)
            self.assertGreater(contig.length, 0)
            self.assertEqual(len(header.filters), 1)
            p = header.filters["PASS"]
            self.assertEqual(p.name, "PASS")
            self.assertEqual(p.description, "All filters passed")
            self.assertEqual(len(header.info), 0)
            self.assertEqual(len(header.formats), 1)
            fmt = header.formats["GT"]
            self.assertEqual(fmt.name, "GT")
            self.assertEqual(fmt.number, 1)
            self.assertEqual(fmt.type, "String")
            self.assertEqual(fmt.description, "Genotype")
            self.assertEqual(len(header.samples), len(datum.sample_names))
            for s1, s2 in zip(header.samples, datum.sample_names):
                self.assertEqual(s1, s2)
            bcf_file.close()


@unittest.skipIf(not _pysam_imported, "pysam not available")
class TestRecordsEqual(unittest.TestCase):
    """
    Tests where we parse the input using PyVCF and Pysam
    """
    def verify_records(self, datum, pyvcf_records, pysam_records):
        self.assertEqual(len(pyvcf_records), len(pysam_records))
        for pyvcf_record, pysam_record in zip(pyvcf_records, pysam_records):
            self.assertEqual(pyvcf_record.CHROM, pysam_record.chrom)
            self.assertEqual(pyvcf_record.POS, pysam_record.pos)
            self.assertEqual(pyvcf_record.ID, pysam_record.id)
            self.assertEqual(pyvcf_record.ALT, list(pysam_record.alts))
            self.assertEqual(pyvcf_record.REF, pysam_record.ref)
            self.assertEqual(pysam_record.filter[0].name, "PASS")
            self.assertEqual(pyvcf_record.FORMAT, "GT")
            self.assertEqual(
                datum.sample_names, list(pysam_record.samples.keys()))
            for value in pysam_record.samples.values():
                self.assertEqual(len(value.alleles), datum.ploidy)
            for j, sample in enumerate(pyvcf_record.samples):
                self.assertEqual(sample.sample, datum.sample_names[j])
                if datum.ploidy > 1:
                    self.assertTrue(sample.phased)
                for call in sample.data.GT.split("|"):
                    self.assertIn(call, ["0", "1"])

    def test_all_records(self):
        for datum in test_data:
            vcf_reader = vcf.Reader(filename=datum.vcf_file)
            bcf_file = pysam.VariantFile(datum.vcf_file)
            pyvcf_records = list(vcf_reader)
            pysam_records = list(bcf_file)
            self.verify_records(datum, pyvcf_records, pysam_records)
            bcf_file.close()


class TestContigLengths(unittest.TestCase):
    """
    Tests that we create sensible contig lengths under a variety of conditions.
    """
    def setUp(self):
        fd, self.temp_file = tempfile.mkstemp(prefix="msprime_vcf_")
        os.close(fd)

    def tearDown(self):
        os.unlink(self.temp_file)

    def get_contig_length(self, ts):
        with open(self.temp_file, "w") as f:
            ts.write_vcf(f)
        reader = vcf.Reader(filename=self.temp_file)
        contig = reader.contigs["1"]
        return contig.length

    def test_no_mutations(self):
        ts = msprime.simulate(10, length=1)
        self.assertEqual(ts.num_mutations, 0)
        contig_length = self.get_contig_length(ts)
        self.assertEqual(contig_length, 1)

    def test_long_sequence(self):
        # Nominal case where we expect the positions to map within the original
        # sequence length
        ts = msprime.simulate(10, length=100, mutation_rate=0.01)
        self.assertGreater(ts.num_mutations, 0)
        contig_length = self.get_contig_length(ts)
        self.assertEqual(contig_length, 100)

    def test_short_sequence(self):
        # Degenerate case where the positions cannot map into the sequence length
        ts = msprime.simulate(10, length=1, mutation_rate=10)
        self.assertGreater(ts.num_mutations, 1)
        contig_length = self.get_contig_length(ts)
        self.assertEqual(contig_length, ts.num_mutations)
