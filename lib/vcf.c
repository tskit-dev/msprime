/*
** Copyright (C) 2016 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
**
** This file is part of msprime.
**
** msprime is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** msprime is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with msprime.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "err.h"
#include "msprime.h"

void
vcf_converter_print_state(vcf_converter_t *self, FILE* out)
{
    fprintf(out, "VCF converter state\n");
    fprintf(out, "ploidy = %d\n", self->ploidy);
    fprintf(out, "sample_size = %d\n", self->sample_size);
    fprintf(out, "num_vcf_samples = %d\n", self->num_vcf_samples);
    fprintf(out, "header = %d bytes\n", (int) strlen(self->header));
    fprintf(out, "genotypes = %d bytes: %s", (int) self->genotypes_size,
            self->genotypes);
    fprintf(out, "record = %d bytes\n", (int) self->record_size);
}

static int WARN_UNUSED
vcf_converter_make_header(vcf_converter_t *self)
{
    int ret = MSP_ERR_GENERIC;
    const char *header_prefix =
        "##fileformat=VCFv4.2\n"
        "##source=msprime " MSP_LIBRARY_VERSION_STR "\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##contig=<ID=1>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    const char *sample_pattern = "\tmsp_%d";
    size_t buffer_size, offset;
    uint32_t j;
    int written;

    buffer_size = strlen(header_prefix);
    offset = buffer_size;
    for (j = 0; j < self->num_vcf_samples; j++) {
        written = snprintf(NULL, 0, sample_pattern, j);
        if (written < 0) {
            ret = MSP_ERR_IO;
            goto out;
        }
        buffer_size += (size_t) written;
    }
    buffer_size += 2; /* make room for \n\0 */
    self->header = malloc(buffer_size);
    if (self->header == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memcpy(self->header, header_prefix, offset);
    for (j = 0; j < self->num_vcf_samples; j++) {
        written = snprintf(self->header + offset, buffer_size - offset,
                sample_pattern, j);
        if (written < 0) {
            ret = MSP_ERR_IO;
            goto out;
        }
        offset += (size_t) written;
        assert(offset < buffer_size);
    }
    self->header[buffer_size - 2] = '\n';
    self->header[buffer_size - 1] = '\0';
    ret = 0;
out:
    return ret;
}

static int WARN_UNUSED
vcf_converter_make_record(vcf_converter_t *self)
{
    int ret = MSP_ERR_GENERIC;
    unsigned int ploidy = self->ploidy;
    uint32_t n = self->num_vcf_samples;
    uint32_t j, k;

    self->genotypes_size = 2 * self->sample_size + 1;
    /* it's not worth working out exactly what size the record prefix
     * will be. 1K is plenty for us */
    self->record_size = 1024 + self->genotypes_size;
    self->record = malloc(self->record_size);
    self->genotypes = malloc(self->genotypes_size);
    if (self->record == NULL || self->genotypes == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    /* Set up the genotypes string. We don't want to have to put
     * in tabs and |s for every row so we insert them at the start.
     */
    for (j = 0; j < n; j++) {
        for (k = 0; k < ploidy; k++) {
            self->genotypes[2 * ploidy * j + 2 * k] = '0';
            self->genotypes[2 * ploidy * j + 2 * k + 1] = '|';
        }
        self->genotypes[2 * ploidy * (j + 1) - 1] = '\t';
    }
    self->genotypes[self->genotypes_size - 2] = '\n';
    self->genotypes[self->genotypes_size - 1] = '\0';
    ret = 0;
out:
    return ret;
}

static int WARN_UNUSED
vcf_converter_write_record(vcf_converter_t *self, unsigned long pos,
        char *variant)
{
    int ret = MSP_ERR_GENERIC;
    int written;
    uint32_t j, k;
    size_t offset;
    unsigned int p = self->ploidy;

    written = snprintf(self->record, self->record_size,
            "1\t%lu\t.\tA\tT\t.\tPASS\t.\tGT\t", pos);
    if (written < 0) {
        ret = MSP_ERR_IO;
        goto out;
    }
    offset = (size_t) written;

    for (j = 0; j < self->num_vcf_samples; j++) {
        for (k = 0; k < p; k++) {
            self->genotypes[2 * p * j + 2 * k] = variant[j * p + k];
        }
    }
    assert(offset + self->genotypes_size < self->record_size);
    memcpy(self->record + offset, self->genotypes, self->genotypes_size);
    ret = 0;
out:
    return ret;
}

int WARN_UNUSED
vcf_converter_get_header(vcf_converter_t *self, char **header)
{
    *header = self->header;
    return 0;
}

int WARN_UNUSED
vcf_converter_next(vcf_converter_t *self, char **record)
{
    int ret = -1;
    int err;
    char *variant;
    unsigned long pos;
    double x;

    ret = vargen_next(self->vargen, &x, &variant);
    if (ret < 0) {
        goto out;
    }
    if (ret == 1) {
        /* update pos. We use a simple algorithm to ensure positions
         * are unique. */
        pos = (unsigned long) round(x);
        if (pos <= self->last_position) {
            pos = self->last_position + 1;
        }
        err = vcf_converter_write_record(self, pos, variant);
        if (err != 0) {
            ret = err;
            goto out;
        }
        *record = self->record;
        self->last_position = pos;
    }

out:
    return ret;
}

int WARN_UNUSED
vcf_converter_alloc(vcf_converter_t *self,
        tree_sequence_t *tree_sequence, unsigned int ploidy)
{
    int ret = -1;

    memset(self, 0, sizeof(vcf_converter_t));
    self->ploidy = ploidy;
    self->sample_size = tree_sequence_get_sample_size(tree_sequence);
    if (ploidy < 1 || self->sample_size % ploidy != 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->num_vcf_samples = self->sample_size / self->ploidy;
    self->vargen = malloc(sizeof(vargen_t));
    if (self->vargen == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = vargen_alloc(self->vargen, tree_sequence);
    if (ret != 0) {
        goto out;
    }
    ret = vcf_converter_make_header(self);
    if (ret != 0) {
        goto out;
    }
    ret = vcf_converter_make_record(self);
    if (ret != 0) {
        goto out;
    }
    /* VCF is 1-based, so we must make sure we never have a 0 coordinate */
    self->last_position = 0;
out:
    return ret;
}

int
vcf_converter_free(vcf_converter_t *self)
{

    if (self->header != NULL) {
        free(self->header);
    }
    if (self->genotypes != NULL) {
        free(self->genotypes);
    }
    if (self->record != NULL) {
        free(self->record);
    }
    if (self->vargen != NULL) {
        vargen_free(self->vargen);
        free(self->vargen);
    }
    return 0;
}
