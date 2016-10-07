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
Test cases for VCF output in msprime.
"""
from __future__ import print_function
from __future__ import division

import collections
import os
import tempfile
import unittest

import msprime

import vcf
import pysam

test_data = []


def setUp():
    Datum = collections.namedtuple(
        "Datum",
        ["tree_sequence", "ploidy", "vcf_file", "sample_names"])
    L = 100
    for ploidy in [1, 2, 3, 5]:
        for n in [2, 10]:
            for rho in [0, 0.5]:
                for mu in [0, 1.0]:
                    ts = msprime.simulate(
                        n * ploidy, length=L, recombination_rate=rho,
                        mutation_rate=mu)
                    fd, file_name = tempfile.mkstemp(prefix="msprime_vcf")
                    with open(file_name, "w") as f:
                        ts.write_vcf(f, ploidy)
                    sample_names = [
                        "msp_{}".format(j) for j in range(n)]
                    test_data.append(
                        Datum(ts, ploidy, file_name, sample_names))


def tearDown():
    for datum in test_data:
        os.unlink(datum.vcf_file)


def write_vcf(tree_sequence, output, ploidy):
    """
    Writes a VCF using the sample algorithm as the low level code.
    """
    if tree_sequence.get_sample_size() % ploidy != 0:
        raise ValueError("Sample size must a multiple of ploidy")
    n = tree_sequence.get_sample_size() // ploidy
    sample_names = ["msp_{}".format(j) for j in range(n)]
    print("##fileformat=VCFv4.2", file=output)
    print("##source=msprime {}".format(msprime.__version__), file=output)
    print(
        '##FILTER=<ID=PASS,Description="All filters passed">',
        file=output)
    print("##contig=<ID=1>", file=output)
    print(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        file=output)
    print(
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "FORMAT", sep="\t", end="", file=output)
    for sample_name in sample_names:
        print("\t", sample_name, sep="", end="", file=output)
    print(file=output)
    last_pos = 0
    for variant in tree_sequence.variants():
        pos = int(round(variant.position))
        if pos <= last_pos:
            pos = last_pos + 1
        print(
            "1", pos, ".", "A", "T", ".", "PASS", ".", "GT",
            sep="\t", end="", file=output)
        for j in range(n):
            genotype = "|".join(
                str(g) for g in
                variant.genotypes[j * ploidy: j * ploidy + ploidy])
            print("\t", genotype, end="", sep="", file=output)
        print(file=output)
        last_pos = pos


class TestEquality(unittest.TestCase):
    """
    Tests if the VCF file produced by the low level code is the
    same as one we generate here.
    """

    def test_equal(self):
        for datum in test_data:
            with tempfile.TemporaryFile("w+") as f:
                write_vcf(datum.tree_sequence, f, datum.ploidy)
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
            contig = reader.contigs["1"]
            self.assertEqual(contig.id, "1")
            self.assertEqual(len(reader.alts), 0)
            self.assertEqual(len(reader.filters), 1)
            p = reader.filters["PASS"]
            self.assertEqual(p.id, "PASS")
            self.assertEqual(len(reader.formats), 1)
            f = reader.formats["GT"]
            self.assertEqual(f.id, "GT")
            self.assertEqual(len(reader.infos), 0)

    def test_pysam(self):
        for datum in test_data:
            bcf_file = pysam.VariantFile(datum.vcf_file)
            self.assertEqual(bcf_file.format, "VCF")
            self.assertEqual(bcf_file.version, (4, 2))
            header = bcf_file.header
            self.assertEqual(len(header.contigs), 1)
            contig = header.contigs[0]
            self.assertEqual(contig.name, "1")
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
