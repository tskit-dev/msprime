/*
** Copyright (C) 2015-2016 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <hdf5.h>

#include <gsl/gsl_math.h>

#include "err.h"
#include "msprime.h"


typedef struct {
    double value;
    uint32_t index;
    double time;
} index_sort_t;

static int
cmp_double(const void *a, const void *b) {
    const double *ia = (const double *) a;
    const double *ib = (const double *) b;
    return (*ia > *ib) - (*ia < *ib);
}

static int
cmp_mutation(const void *a, const void *b) {
    const mutation_t *ia = (const mutation_t *) a;
    const mutation_t *ib = (const mutation_t *) b;
    return (ia->position > ib->position) - (ia->position < ib->position);
}

static int
cmp_mutation_pointer(const void *a, const void *b) {
    mutation_t *const*ia = (mutation_t *const*) a;
    mutation_t *const*ib = (mutation_t *const*) b;
    return cmp_mutation(*ia, *ib);
}

static int
cmp_index_sort(const void *a, const void *b) {
    const index_sort_t *ca = (const index_sort_t *) a;
    const index_sort_t *cb = (const index_sort_t *) b;
    int ret = (ca->value > cb->value) - (ca->value < cb->value);
    /* When comparing equal values, we sort by time */
    if (ret == 0) {
        ret = (ca->time > cb->time) - (ca->time < cb->time);
    }
    return ret;
}

void
tree_sequence_print_state(tree_sequence_t *self, FILE *out)
{
    size_t j, k;

    fprintf(out, "tree_sequence state\n");
    fprintf(out, "sample_size = %d\n", self->sample_size);
    fprintf(out, "provenance = (%d)\n", (int) self->num_provenance_strings);
    for (j = 0; j < self->num_provenance_strings; j++) {
        fprintf(out, "\t'%s'\n", self->provenance_strings[j]);
    }
    fprintf(out, "sequence_length = %f\n", self->sequence_length);
    fprintf(out, "nodes (%d)\n", (int) self->num_nodes);
    for (j = 0; j < self->num_nodes; j++) {
        fprintf(out, "\t%d\t%d\t%f\n", (int) j,
                (int) self->trees.nodes.population[j],
                self->trees.nodes.time[j]);
    }
    fprintf(out, "breakpoints (%d)\n", (int) self->trees.num_breakpoints);
    for (j = 0; j < self->trees.num_breakpoints; j++) {
        fprintf(out, "\t%d\t%f\n", (int) j,
                self->trees.breakpoints[j]);
    }

    fprintf(out, "trees.records = (%d records)\n", (int) self->num_records);
    for (j = 0; j < self->num_records; j++) {
        fprintf(out, "\t%d\t%d\t%d\t%d\t(",
                (int) j,
                self->trees.records.left[j],
                self->trees.records.right[j],
                (int) self->trees.records.node[j]);
        for (k = 0; k < self->trees.records.num_children[j]; k++) {
            fprintf(out, "%d", self->trees.records.children[j][k]);
            if (k < self->trees.records.num_children[j] - 1) {
                fprintf(out, ", ");
            }
        }
        fprintf(out, ")\t|\t%d\t%d\n",
                (int) self->trees.indexes.insertion_order[j],
                (int) self->trees.indexes.removal_order[j]);
    }
    fprintf(out, "mutations = (%d records)\n", (int) self->num_mutations);
    for (j = 0; j < self->num_mutations; j++) {
        fprintf(out, "\t%d\t%f\t%d\n", (int) j, self->mutations.position[j],
               (int) self->mutations.node[j]);
    }
}

/* Allocates the memory required for arrays of values. Assumes that
 * the num_records and num_mutations have been set.
 */
static int
tree_sequence_alloc(tree_sequence_t *self)
{
    int ret = MSP_ERR_NO_MEMORY;
    uint32_t j;

    self->trees.nodes.time = malloc(self->num_nodes * sizeof(double));
    self->trees.nodes.population = malloc(self->num_nodes * sizeof(uint8_t));
    if (self->trees.nodes.time == NULL || self->trees.nodes.population == NULL) {
        goto out;
    }
    self->trees.breakpoints = malloc(
            self->trees.num_breakpoints * sizeof(double));
    if (self->trees.breakpoints == NULL) {
        goto out;
    }
    self->trees.records.left = malloc(self->num_records * sizeof(double));
    self->trees.records.right = malloc(self->num_records * sizeof(double));
    self->trees.records.num_children = malloc(self->num_records * sizeof(uint32_t));
    self->trees.records.children = malloc(self->num_records * sizeof(uint32_t *));
    self->trees.records.node = malloc(self->num_records * sizeof(uint32_t));
    self->trees.records.children_mem = malloc(self->num_child_nodes * sizeof(uint32_t));
    if (self->trees.records.left == NULL
            || self->trees.records.right == NULL
            || self->trees.records.children == NULL
            || self->trees.records.node == NULL
            || self->trees.records.num_children == NULL
            || self->trees.records.children_mem == NULL) {
        goto out;
    }
    self->trees.indexes.insertion_order = malloc(self->num_records * sizeof(uint32_t));
    self->trees.indexes.removal_order = malloc(self->num_records * sizeof(uint32_t));
    if (self->trees.indexes.insertion_order == NULL
            || self->trees.indexes.removal_order == NULL) {
        goto out;
    }
    /* Set the optional fields to their unset values. */
    for (j = 0; j < self->num_nodes; j++) {
        self->trees.nodes.population[j] = MSP_NULL_POPULATION_ID;
        self->trees.nodes.time[j] = 0.0;
    }
    if (self->num_mutations > 0) {
        self->mutations.node = malloc(self->num_mutations * sizeof(uint32_t));
        self->mutations.position = malloc(
                self->num_mutations * sizeof(double));
        if (self->mutations.node == NULL || self->mutations.position == NULL) {
            goto out;
        }
    }
    /* Avoid the potential portability issues with malloc(0) here */
    self->provenance_strings = malloc((1 + self->num_provenance_strings)
            * sizeof(char *));
    if (self->provenance_strings == NULL) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_free(tree_sequence_t *self)
{
    size_t j;

    if (self->provenance_strings != NULL) {
        for (j = 0; j < self->num_provenance_strings; j++) {
            free(self->provenance_strings[j]);
        }
        free(self->provenance_strings);
    }
    if (self->trees.nodes.population != NULL) {
        free(self->trees.nodes.population);
    }
    if (self->trees.nodes.time != NULL) {
        free(self->trees.nodes.time);
    }
    if (self->trees.breakpoints != NULL) {
        free(self->trees.breakpoints);
    }
    if (self->trees.records.left != NULL) {
        free(self->trees.records.left);
    }
    if (self->trees.records.right != NULL) {
        free(self->trees.records.right);
    }
    if (self->trees.records.children != NULL) {
        free(self->trees.records.children);
    }
    if (self->trees.records.num_children != NULL) {
        free(self->trees.records.num_children);
    }
    if (self->trees.records.children_mem != NULL) {
        free(self->trees.records.children_mem);
    }
    if (self->trees.records.node != NULL) {
        free(self->trees.records.node);
    }
    if (self->trees.indexes.insertion_order != NULL) {
        free(self->trees.indexes.insertion_order);
    }
    if (self->trees.indexes.removal_order != NULL) {
        free(self->trees.indexes.removal_order);
    }
    if (self->mutations.node != NULL) {
        free(self->mutations.node);
    }
    if (self->mutations.position != NULL) {
        free(self->mutations.position);
    }
    return 0;
}

int WARN_UNUSED
tree_sequence_add_provenance_string(tree_sequence_t *self,
        const char *provenance_string)
{
    int ret = MSP_ERR_GENERIC;
    char **p, *s;
    size_t size;

    if (provenance_string == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    size = strlen(provenance_string);
    if (size == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->num_provenance_strings++;
    p = realloc(self->provenance_strings,
            self->num_provenance_strings * sizeof(char *));
    if (p == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->provenance_strings = p;
    size++; /* allow for '/0' */
    s = malloc((size) * sizeof(char));
    if (s == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    strncpy(s, provenance_string, size);
    self->provenance_strings[self->num_provenance_strings - 1] = s;
    ret = 0;
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_get_provenance_strings(tree_sequence_t *self,
        size_t *num_provenance_strings, char ***provenance_strings)
{
    *num_provenance_strings = self->num_provenance_strings;
    *provenance_strings = self->provenance_strings;
    return 0;
}

static int
tree_sequence_check(tree_sequence_t *self)
{
    int ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
    uint32_t j, k, left;

    left = UINT32_MAX;
    for (j = 0; j < self->num_records; j++) {
        if (j > 0) {
            /* Input data must be time sorted. */
            if (self->trees.nodes.time[self->trees.records.node[j]]
                    < self->trees.nodes.time[self->trees.records.node[j - 1]]) {
                goto out;
            }
        }
        if (self->trees.records.num_children[j] < 2) {
            goto out;
        }
        left = GSL_MIN(left, self->trees.records.left[j]);
        /* Ensure that children are non-null and in ascending order */
        for (k = 0; k < self->trees.records.num_children[j]; k++) {
            if (k < self->trees.records.num_children[j] - 1) {
                if (self->trees.records.children[j][k]
                        >= self->trees.records.children[j][k + 1]) {
                    goto out;
                }
            }
            if (self->trees.records.children[j][k] == MSP_NULL_NODE) {
                goto out;
            }
        }
        if (self->trees.records.left[j] >= self->trees.records.right[j]) {
            goto out;
        }
    }
    if (left != 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_init_from_records(tree_sequence_t *self,
      size_t num_records, coalescence_record_t *records)
{

    int ret = MSP_ERR_GENERIC;
    uint32_t node;
    size_t j, k, offset;
    double last_breakpoint;
    double *left = NULL;
    index_sort_t *sort_buff = NULL;

    memset(self, 0, sizeof(tree_sequence_t));
    if (num_records == 0) {
        ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
        goto out;
    }
    left = malloc((num_records + 1) * sizeof(double));
    if (left == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    /* Make the first pass through the records to see how many
     * child nodes we have in total as well finding the sample
     * size. Also do some basic error checking.
     */
    self->sample_size = UINT32_MAX;
    self->num_mutations = 0;
    self->num_child_nodes = 0;
    self->sequence_length = 0.0;
    self->num_records = num_records;
    self->num_nodes = 0;
    for (j = 0; j < self->num_records; j++) {
        self->num_child_nodes += records[j].num_children;
        if (records[j].num_children < 2) {
            ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
            goto out;
        }
        self->sample_size = GSL_MIN(self->sample_size, records[j].node);
        self->num_nodes = GSL_MAX(self->num_nodes, records[j].node);
        self->sequence_length = GSL_MAX(self->sequence_length,
                records[j].right);
        left[j] = records[j].left;
    }

    if (self->sample_size < 2 || self->sample_size == MSP_NULL_NODE) {
        ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
        goto out;
    }
    if (self->sequence_length <= 0) {
        ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
        goto out;
    }
    self->num_nodes++;
    left[num_records] = self->sequence_length;
    qsort(left, num_records + 1, sizeof(double), cmp_double);
    self->trees.num_breakpoints = 0;
    last_breakpoint = -1.0;
    for (j = 0; j < num_records + 1; j++) {
        if (left[j] != last_breakpoint) {
            self->trees.num_breakpoints++;
            last_breakpoint = left[j];
        }
    }
    /* Now alloc the memory and store the records. */
    ret = tree_sequence_alloc(self);
    if (ret != 0) {
        goto out;
    }
    /* First store the breakpoints */
    last_breakpoint = -1.0;
    k = 0;
    for (j = 0; j < num_records + 1; j++) {
        if (left[j] != last_breakpoint) {
            self->trees.breakpoints[k] = left[j];
            k++;
            last_breakpoint = left[j];
        }
    }
    /* Set up the nodes and the children pointers */
    offset = 0;
    for (j = 0; j < self->num_records; j++) {
        node = records[j].node;
        if (self->trees.nodes.time[node] == 0.0) {
            self->trees.nodes.time[node] = records[j].time;
        } else if (self->trees.nodes.time[node] != records[j].time) {
            ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
            goto out;
        }
        if (self->trees.nodes.population[node] == MSP_NULL_POPULATION_ID) {
            self->trees.nodes.population[node] = records[j].population_id;
        } else if (self->trees.nodes.population[node] != records[j].population_id) {
            ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
            goto out;
        }
        self->trees.records.node[j] = records[j].node;
        self->trees.records.num_children[j] = records[j].num_children;
        self->trees.records.children[j] = &self->trees.records.children_mem[offset];
        offset += records[j].num_children;
        for (k = 0; k < records[j].num_children; k++) {
            self->trees.records.children[j][k] = records[j].children[k];
        }
    }
    assert(offset == self->num_child_nodes);

    /* Now sort create the indexes and set the breakpoint indexes in
     * left and right. */
    sort_buff = malloc(self->num_records * sizeof(index_sort_t));
    if (sort_buff == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    /* sort by left and increasing time to give us the order in which
     * records should be inserted */
    for (j = 0; j < self->num_records; j++) {
        sort_buff[j].index = (uint32_t ) j;
        sort_buff[j].value = records[j].left;
        sort_buff[j].time = records[j].time;
    }
    qsort(sort_buff, self->num_records, sizeof(index_sort_t), cmp_index_sort);
    k = 0;
    for (j = 0; j < self->num_records; j++) {
        self->trees.indexes.insertion_order[j] = sort_buff[j].index;
        while (self->trees.breakpoints[k] < sort_buff[j].value) {
            k++;
        }
        assert(k < self->trees.num_breakpoints);
        self->trees.records.left[sort_buff[j].index] = (uint32_t) k;
    }
    /* sort by right and decreasing time to give us the order in which
     * records should be removed. */
    for (j = 0; j < self->num_records; j++) {
        sort_buff[j].index = (uint32_t ) j;
        sort_buff[j].value = records[j].right;
        sort_buff[j].time = -records[j].time;
    }
    qsort(sort_buff, self->num_records, sizeof(index_sort_t), cmp_index_sort);
    k = 0;
    for (j = 0; j < self->num_records; j++) {
        self->trees.indexes.removal_order[j] = sort_buff[j].index;
        while (self->trees.breakpoints[k] < sort_buff[j].value) {
            k++;
        }
        assert(k < self->trees.num_breakpoints);
        /* If we can't find the value in breakpoints, it means that
         * we have right coordinates not matching to a left coord */
        if (self->trees.breakpoints[k] != sort_buff[j].value) {
            ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
            goto out;
        }
        self->trees.records.right[sort_buff[j].index] = (uint32_t) k;
    }
    ret = tree_sequence_check(self);
out:
    if (left != NULL) {
        free(left);
    }
    if (sort_buff != NULL) {
        free(sort_buff);
    }
    return ret;
}

int WARN_UNUSED
tree_sequence_load_records(tree_sequence_t *self,
      size_t num_records, coalescence_record_t *records)
{
    int ret = MSP_ERR_GENERIC;

    ret = tree_sequence_init_from_records(self, num_records, records);
    if (ret != 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_create(tree_sequence_t *self, msp_t *sim,
        recomb_map_t *recomb_map, double Ne)
{
    int ret = MSP_ERR_GENERIC;
    size_t j, num_records;
    coalescence_record_t *records = NULL;
    sample_t *samples = NULL;

    ret = msp_get_coalescence_records(sim, &records);
    if (ret != 0) {
        goto out;
    }
    num_records = msp_get_num_coalescence_records(sim);
    ret = tree_sequence_init_from_records(self, num_records, records);
    if (ret != 0) {
        goto out;
    }
    assert(self->sample_size == msp_get_sample_size(sim));
    assert(self->sequence_length == (double) msp_get_num_loci(sim));
    assert(self->num_records == msp_get_num_coalescence_records(sim));

    ret = msp_get_samples(sim, &samples);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_set_samples(self, self->sample_size, samples);
    if (ret != 0) {
        goto out;
    }
    /* Rescale node times into generations */
    for (j = 0; j < self->num_nodes; j++) {
        self->trees.nodes.time[j] *= 4 * Ne;
    }
    self->sequence_length = recomb_map_get_sequence_length(recomb_map);
    ret = recomb_map_genetic_to_phys_bulk(
        recomb_map, self->trees.breakpoints, self->trees.num_breakpoints);
out:
    return ret;
}

/* Reads the metadata for the overall file and updates the basic
 * information in the tree_sequence.
 */
static int
tree_sequence_read_hdf5_metadata(tree_sequence_t *self, hid_t file_id)
{
    int ret = MSP_ERR_HDF5;
    hid_t attr_id, dataspace_id;
    herr_t status;
    int rank;
    hsize_t dims;
    int32_t version[2];
    struct _hdf5_metadata_read {
        const char *prefix;
        const char *name;
        hid_t memory_type;
        size_t size;
        void *dest;
    };
    struct _hdf5_metadata_read fields[] = {
        {"/", "format_version", H5T_NATIVE_UINT32, 2, NULL},
        {"/", "sample_size", H5T_NATIVE_UINT32, 0, NULL},
        {"/", "sequence_length", H5T_NATIVE_DOUBLE, 0, NULL},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _hdf5_metadata_read);
    size_t j;

    fields[0].dest = version;
    fields[1].dest = &self->sample_size;
    fields[2].dest = &self->sequence_length;
    for (j = 0; j < num_fields; j++) {
        attr_id = H5Aopen_by_name(file_id, fields[j].prefix, fields[j].name,
                H5P_DEFAULT, H5P_DEFAULT);
        if (attr_id < 0) {
            goto out;
        }
        dataspace_id = H5Aget_space(attr_id);
        if (dataspace_id < 0) {
            goto out;
        }
        rank = H5Sget_simple_extent_ndims(dataspace_id);
        if (fields[j].size == 0) {
            /* SCALAR's have rank 0 */
            if (rank != 0) {
                ret = MSP_ERR_FILE_FORMAT;
                goto out;
            }
        } else {
            if (rank != 1) {
                ret = MSP_ERR_FILE_FORMAT;
                goto out;
            }
            status = H5Sget_simple_extent_dims(dataspace_id, &dims, NULL);
            if (status < 0) {
                goto out;
            }
            if (dims != fields[j].size) {
                ret = MSP_ERR_FILE_FORMAT;
                goto out;
            }
        }
        status = H5Aread(attr_id, fields[j].memory_type, fields[j].dest);
        if (status < 0) {
            goto out;
        }
        status = H5Sclose(dataspace_id);
        if (status < 0) {
            goto out;
        }
        status = H5Aclose(attr_id);
        if (status < 0) {
            goto out;
        }
    }
    /* Sanity check */
    if (version[0] != MSP_FILE_FORMAT_VERSION_MAJOR) {
        ret = MSP_ERR_UNSUPPORTED_FILE_VERSION;
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_check_hdf5_dimensions(tree_sequence_t *self, hid_t file_id)
{
    int ret = MSP_ERR_HDF5;
    hid_t dataset_id, dataspace_id;
    herr_t status;
    int rank;
    hsize_t dims[2];
    struct _dimension_check {
        const char *name;
        int check_size;
        size_t size;
        int required;
    };
    struct _dimension_check fields[] = {
        {"/mutations/node", 1, self->num_mutations, 1},
        {"/mutations/position", 1, self->num_mutations, 1},
        {"/trees/nodes/population", 1, self->num_nodes, 0},
        {"/trees/nodes/time", 1, self->num_nodes, 0},
        {"/trees/breakpoints", 1, self->trees.num_breakpoints, 0},
        {"/trees/records/left", 1, self->num_records, 1},
        {"/trees/records/right", 1, self->num_records, 1},
        {"/trees/records/node", 1, self->num_records, 1},
        {"/trees/records/num_children", 1, self->num_records, 1},
        {"/trees/records/children", 0, self->num_child_nodes, 1},
        {"/trees/indexes/insertion_order", 1, self->num_records, 1},
        {"/trees/indexes/removal_order", 1, self->num_records, 1},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _dimension_check);
    size_t j;

    for (j = 0; j < 2; j++) {
        fields[j].size = self->num_mutations;
        fields[j].required = self->num_mutations > 0;
    }
    for (j = 0; j < num_fields; j++) {
        if (fields[j].required) {
            dataset_id = H5Dopen(file_id, fields[j].name, H5P_DEFAULT);
            if (dataset_id < 0) {
                goto out;
            }
            dataspace_id = H5Dget_space(dataset_id);
            if (dataspace_id < 0) {
                goto out;
            }
            rank = H5Sget_simple_extent_ndims(dataspace_id);
            if (rank != 1) {
                ret = MSP_ERR_FILE_FORMAT;
                goto out;
            }
            status = H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
            if (status < 0) {
                goto out;
            }
            if (fields[j].check_size && dims[0] != fields[j].size) {
                ret = MSP_ERR_FILE_FORMAT;
                goto out;
            }
            status = H5Sclose(dataspace_id);
            if (status < 0) {
                goto out;
            }
            status = H5Dclose(dataset_id);
            if (status < 0) {
                goto out;
            }
        }
    }
    ret = 0;
out:
    return ret;
}

/* Reads the dimensions for the records and mutations and mallocs
 * space.
 */
static int
tree_sequence_read_hdf5_dimensions(tree_sequence_t *self, hid_t file_id)
{
    int ret = MSP_ERR_HDF5;
    hid_t dataset_id, dataspace_id;
    herr_t status;
    htri_t exists;
    int rank;
    hsize_t dims[2];
    struct _dimension_read {
        const char *name;
        size_t *dest;
        int included;
    };
    struct _dimension_read fields[] = {
        {"/mutations/node", &self->num_mutations, 0},
        {"/provenance", &self->num_provenance_strings, 0},
        {"/trees/breakpoints", &self->trees.num_breakpoints, 1},
        {"/trees/nodes/time", &self->num_nodes, 1},
        {"/trees/records/left", &self->num_records, 1},
        {"/trees/records/children", &self->num_child_nodes, 1},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _dimension_read);
    size_t j;
    /* check if the mutations group exists. This seems a bit awkward, but it's
     * an error to call H5Lexists on /mutations/node if /mutations doesn't
     * exist */
    exists = H5Lexists(file_id, "/mutations", H5P_DEFAULT);
    if (exists < 0) {
        goto out;
    }
    self->num_mutations = 0;
    if (exists) {
        fields[0].included = 1;
    }
    /* check if provenance exists */
    exists = H5Lexists(file_id, "/provenance", H5P_DEFAULT);
    if (exists < 0) {
        goto out;
    }
    self->num_provenance_strings = 0;
    if (exists) {
        fields[1].included = 1;
    }
    for (j = 0; j < num_fields; j++) {
        if (fields[j].included) {
            dataset_id = H5Dopen(file_id, fields[j].name, H5P_DEFAULT);
            if (dataset_id < 0) {
                goto out;
            }
            dataspace_id = H5Dget_space(dataset_id);
            if (dataspace_id < 0) {
                goto out;
            }
            rank = H5Sget_simple_extent_ndims(dataspace_id);
            if (rank != 1) {
                ret = MSP_ERR_FILE_FORMAT;
                goto out;
            }
            status = H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
            if (status < 0) {
                goto out;
            }
            *fields[j].dest = (size_t) dims[0];
            status = H5Sclose(dataspace_id);
            if (status < 0) {
                goto out;
            }
            status = H5Dclose(dataset_id);
            if (status < 0) {
                goto out;
            }
        }
    }
    ret = tree_sequence_check_hdf5_dimensions(self, file_id);
    if (ret != 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_read_hdf5_data(tree_sequence_t *self, hid_t file_id)
{
    herr_t status;
    int ret = MSP_ERR_HDF5;
    hid_t dataset_id;
    struct _hdf5_field_read {
        const char *name;
        hid_t type;
        int empty;
        int required;
        void *dest;
    };
    struct _hdf5_field_read fields[] = {
        {"/provenance", 0, 0, 0, self->provenance_strings},
        {"/mutations/node", H5T_NATIVE_UINT32, 0, 1,
            self->mutations.node},
        {"/mutations/position", H5T_NATIVE_DOUBLE, 0, 1,
            self->mutations.position},
        {"/trees/nodes/population", H5T_NATIVE_UINT8, 0, 0,
            self->trees.nodes.population},
        {"/trees/nodes/time", H5T_NATIVE_DOUBLE, 0, 1,
            self->trees.nodes.time},
        {"/trees/breakpoints", H5T_NATIVE_DOUBLE, 0, 1,
            self->trees.breakpoints},
        {"/trees/records/left", H5T_NATIVE_UINT32, 0, 1,
            self->trees.records.left},
        {"/trees/records/right", H5T_NATIVE_UINT32, 0, 1,
            self->trees.records.right},
        {"/trees/records/node", H5T_NATIVE_UINT32, 0, 1,
            self->trees.records.node},
        {"/trees/records/num_children", H5T_NATIVE_UINT32, 0, 1,
            self->trees.records.num_children},
        {"/trees/records/children", H5T_NATIVE_UINT32, 0, 1,
            self->trees.records.children_mem},
        {"/trees/indexes/insertion_order", H5T_NATIVE_UINT32, 0, 1,
            self->trees.indexes.insertion_order},
        {"/trees/indexes/removal_order", H5T_NATIVE_UINT32, 0, 1,
            self->trees.indexes.removal_order},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _hdf5_field_read);
    size_t j, offset;
    hid_t vlen_str;

    vlen_str = H5Tcopy(H5T_C_S1);
    if (vlen_str < 0) {
        goto out;
    }
    status = H5Tset_size(vlen_str, H5T_VARIABLE);
    if (status < 0) {
        goto out;
    }
    fields[0].type = vlen_str;

    /* TODO We're sort of doing the same thing twice here as
     * the mutations _group_ is optional. However, we can't just
     * mark mutations/node and mutations/position as optional as we
     * would then allow one or the other. This would be an error.
     * However, we should improve this logic as it's a bit messy.
     */
    if (self->num_mutations == 0) {
        fields[1].empty = 1;
        fields[2].empty = 1;
    }
    for (j = 0; j < num_fields; j++) {
        /* Skip any non-required fields that are missing. */
        if (!fields[j].required
                && H5Lexists(file_id, fields[j].name, H5P_DEFAULT) <= 0) {
            continue;
        }
        /* Skip any fields that are marked empty */
        if (fields[j].empty) {
            continue;
        }
        /* If we got this far, read in the field. */
        dataset_id = H5Dopen(file_id, fields[j].name, H5P_DEFAULT);
        if (dataset_id < 0) {
            goto out;
        }
        status = H5Dread(dataset_id, fields[j].type, H5S_ALL,
                H5S_ALL, H5P_DEFAULT, fields[j].dest);
        if (status < 0) {
            goto out;
        }
        status = H5Dclose(dataset_id);
        if (status < 0) {
            goto out;
        }
    }
    /* Now update the children vectors */
    offset = 0;
    for (j = 0; j < self->num_records; j++) {
        assert(offset < self->num_child_nodes);
        self->trees.records.children[j] = &self->trees.records.children_mem[offset];
        offset += self->trees.records.num_children[j];
    }
    assert(offset == self->num_child_nodes);
    status = H5Tclose(vlen_str);
    if (status < 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_load(tree_sequence_t *self, const char *filename, int flags)
{
    int ret = MSP_ERR_GENERIC;
    herr_t status;
    hid_t file_id;

    memset(self, 0, sizeof(tree_sequence_t));
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        ret = MSP_ERR_HDF5;
        goto out;
    }
    ret = tree_sequence_read_hdf5_metadata(self, file_id);
    if (ret < 0) {
        goto out;
    }
    ret = tree_sequence_read_hdf5_dimensions(self, file_id);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_alloc(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_read_hdf5_data(self, file_id);
    if (ret != 0) {
        goto out;
    }
    status = H5Fclose(file_id);
    if (status < 0) {
        ret = MSP_ERR_HDF5;
        goto out;
    }
    if (!(flags & MSP_SKIP_H5CLOSE)) {
        status = H5close();
        if (status < 0) {
            goto out;
        }
    }
    ret = tree_sequence_check(self);
out:
    return ret;
}

static int
tree_sequence_write_hdf5_data(tree_sequence_t *self, hid_t file_id, int flags)
{
    herr_t ret = -1;
    herr_t status;
    hid_t group_id, dataset_id, dataspace_id, plist_id;
    hsize_t dims[1];
    struct _hdf5_field_write {
        const char *name;
        hid_t storage_type;
        hid_t memory_type;
        size_t size;
        void *source;
    };
    struct _hdf5_field_write fields[] = {
        {"/provenance",
            0, 0, /* We must set this afterwards */
            self->num_provenance_strings, self->provenance_strings},
        {"/trees/nodes/population",
            H5T_STD_U8LE, H5T_NATIVE_UINT8,
            self->num_nodes, self->trees.nodes.population},
        {"/trees/nodes/time",
            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
            self->num_nodes, self->trees.nodes.time},
        {"/trees/records/left",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->num_records, self->trees.records.left},
        {"/trees/records/right",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->num_records, self->trees.records.right},
        {"/trees/records/node",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->num_records, self->trees.records.node},
        {"/trees/records/num_children",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->num_records, self->trees.records.num_children},
        {"/trees/records/children",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->num_child_nodes, self->trees.records.children_mem},
        {"/trees/indexes/insertion_order",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->num_records, self->trees.indexes.insertion_order},
        {"/trees/indexes/removal_order",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->num_records, self->trees.indexes.removal_order},
        {"/trees/breakpoints",
            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
            self->trees.num_breakpoints, self->trees.breakpoints},
        {"/mutations/node",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->num_mutations, self->mutations.node},
        {"/mutations/position",
            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
            self->num_mutations, self->mutations.position},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _hdf5_field_write);
    struct _hdf5_group_write {
        const char *name;
        int included;
    };
    struct _hdf5_group_write groups[] = {
        {"/mutations", 1},
        {"/trees", 1},
        {"/trees/nodes", 1},
        {"/trees/records", 1},
        {"/trees/indexes", 1},
    };
    size_t num_groups = sizeof(groups) / sizeof(struct _hdf5_group_write);
    size_t j;
    hid_t vlen_str;

    vlen_str = H5Tcopy(H5T_C_S1);
    if (vlen_str < 0) {
        goto out;
    }
    status = H5Tset_size(vlen_str, H5T_VARIABLE);
    if (status < 0) {
        goto out;
    }
    fields[0].storage_type = vlen_str;
    fields[0].memory_type = vlen_str;

    /* We only create the mutations group if it's non-empty */
    if (self->num_mutations == 0) {
        groups[0].included = 0;
    }
    /* Create the groups */
    for (j = 0; j < num_groups; j++) {
        if (groups[j].included) {
            group_id = H5Gcreate(file_id, groups[j].name, H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
            if (group_id < 0) {
                goto out;
            }
            status = H5Gclose(group_id);
            if (status < 0) {
                goto out;
            }
        }
    }
    /* now write the datasets */
    for (j = 0; j < num_fields; j++) {
        if (fields[j].size > 0) {
            dims[0] = fields[j].size;
            dataspace_id = H5Screate_simple(1, dims, NULL);
            if (dataspace_id < 0) {
                goto out;
            }
            plist_id = H5Pcreate(H5P_DATASET_CREATE);
            if (plist_id < 0) {
                goto out;
            }
            /* Set the chunk size to the full size of the dataset since we
             * always read the full thing.
             */
            status = H5Pset_chunk(plist_id, 1, dims);
            if (status < 0) {
                goto out;
            }
            if (fields[j].memory_type != H5T_NATIVE_DOUBLE &&
                    fields[j].memory_type != vlen_str) {
                /* For integer types, use the scale offset compression */
                status = H5Pset_scaleoffset(plist_id, H5Z_SO_INT,
                         H5Z_SO_INT_MINBITS_DEFAULT);
                if (status < 0) {
                    goto out;
                }
            }
            if (flags & MSP_ZLIB_COMPRESSION) {
                /* Turn on byte shuffling to improve compression */
                status = H5Pset_shuffle(plist_id);
                if (status < 0) {
                    goto out;
                }
                /* Set zlib compression at level 9 (best compression) */
                status = H5Pset_deflate(plist_id, 9);
                if (status < 0) {
                    goto out;
                }
            }
            /* Turn on Fletcher32 checksums for integrity checks */
            status = H5Pset_fletcher32(plist_id);
            if (status < 0) {
                goto out;
            }
            dataset_id = H5Dcreate2(file_id, fields[j].name,
                    fields[j].storage_type, dataspace_id, H5P_DEFAULT,
                    plist_id, H5P_DEFAULT);
            if (dataset_id < 0) {
                goto out;
            }
            if (fields[j].size > 0) {
                /* Don't write zero sized datasets to work-around problems
                 * with older versions of hdf5. */
                status = H5Dwrite(dataset_id, fields[j].memory_type, H5S_ALL,
                        H5S_ALL, H5P_DEFAULT, fields[j].source);
                if (status < 0) {
                    goto out;
                }
            }
            status = H5Dclose(dataset_id);
            if (status < 0) {
                goto out;
            }
            status = H5Pclose(plist_id);
            if (status < 0) {
                goto out;
            }
        }
    }
    status = H5Tclose(vlen_str);
    if (status < 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_write_hdf5_metadata(tree_sequence_t *self, hid_t file_id)
{
    herr_t status = -1;
    hid_t attr_id, dataspace_id;
    hsize_t dims = 1;
    uint32_t version[2] = {
        MSP_FILE_FORMAT_VERSION_MAJOR, MSP_FILE_FORMAT_VERSION_MINOR};

    struct _hdf5_metadata_write {
        const char *name;
        hid_t parent;
        hid_t storage_type;
        hid_t memory_type;
        size_t size;
        void *source;
    };
    struct _hdf5_metadata_write fields[] = {
        {"format_version", 0, H5T_STD_U32LE, H5T_NATIVE_UINT32, 2, NULL},
        {"sample_size", 0, H5T_STD_U32LE, H5T_NATIVE_UINT32, 0, NULL},
        {"sequence_length", 0, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, 0, NULL},
    };
    /* TODO random_seed, population_models, etc. */
    size_t num_fields = sizeof(fields) / sizeof(struct _hdf5_metadata_write);
    size_t j;

    fields[0].source = version;
    fields[1].source = &self->sample_size;
    fields[2].source = &self->sequence_length;

    for (j = 0; j < num_fields; j++) {
        if (fields[j].size == 0) {
            dataspace_id = H5Screate(H5S_SCALAR);
        } else {
            dims = fields[j].size;
            dataspace_id = H5Screate_simple(1, &dims, NULL);
        }
        if (dataspace_id < 0) {
            status = dataspace_id;
            goto out;
        }
        attr_id = H5Acreate(file_id, fields[j].name,
                fields[j].storage_type, dataspace_id, H5P_DEFAULT,
                H5P_DEFAULT);
        if (attr_id < 0) {
            goto out;
        }
        status = H5Awrite(attr_id, fields[j].memory_type, fields[j].source);
        if (status < 0) {
            goto out;
        }
        status = H5Aclose(attr_id);
        if (status < 0) {
            goto out;
        }
        status = H5Sclose(dataspace_id);
        if (status < 0) {
            goto out;
        }
    }
 out:
    return status;
}

int WARN_UNUSED
tree_sequence_dump(tree_sequence_t *self, const char *filename, int flags)
{
    int ret = MSP_ERR_HDF5;
    herr_t status;
    hid_t file_id;

    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        goto out;
    }
    status = tree_sequence_write_hdf5_metadata(self, file_id);
    if (status < 0) {
        goto out;
    }
    status = tree_sequence_write_hdf5_data(self, file_id, flags);
    if (status < 0) {
        goto out;
    }
    status = H5Fclose(file_id);
    if (status < 0) {
        goto out;
    }
    if (!(flags & MSP_SKIP_H5CLOSE)) {
        status = H5close();
        if (status < 0) {
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

double
tree_sequence_get_sequence_length(tree_sequence_t *self)
{
    return self->sequence_length;
}

uint32_t
tree_sequence_get_sample_size(tree_sequence_t *self)
{
    return self->sample_size;
}

uint32_t
tree_sequence_get_num_nodes(tree_sequence_t *self)
{
    return (uint32_t) self->num_nodes;
}

int WARN_UNUSED
tree_sequence_get_sample(tree_sequence_t *self, uint32_t u, sample_t *sample)
{
    int ret = 0;

    if (u >= self->sample_size) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    sample->population_id = self->trees.nodes.population[u];
    sample->time = self->trees.nodes.time[u];
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_get_pairwise_diversity(tree_sequence_t *self,
    uint32_t *samples, uint32_t num_samples, double *pi)
{
    int ret = 0;
    uint32_t j, node;
    sparse_tree_t *tree = NULL;
    sparse_tree_iterator_t *tree_iter = NULL;
    double result, denom, count;

    tree = malloc(sizeof(sparse_tree_t));
    if (tree == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = tree_sequence_alloc_sparse_tree(self, tree,
        samples, num_samples, MSP_COUNT_LEAVES);
    if (ret != 0) {
        goto out;
    }
    tree_iter = malloc(sizeof(sparse_tree_iterator_t));
    if (tree_iter == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = sparse_tree_iterator_alloc(tree_iter, self, tree);
    if (ret != 0) {
        goto out;
    }
    /* Allocation done; move onto main algorithm. */
    denom = (num_samples * ((double) num_samples - 1)) / 2.0;
    result = 0.0;
    while ((ret = sparse_tree_iterator_next(tree_iter)) == 1) {
        for (j = 0; j < tree->num_mutations; j++) {
            node = tree->mutations[j].node;
            count = (double) tree->num_tracked_leaves[node];
            result += count * (num_samples - count) / denom;
        }
    }
    if (ret != 0) {
        goto out;
    }
    *pi = result;
out:
    if (tree != NULL) {
        sparse_tree_free(tree);
        free(tree);
    }
    if (tree_iter != NULL) {
        sparse_tree_iterator_free(tree_iter);
        free(tree_iter);
    }
    return ret;
}

size_t
tree_sequence_get_num_coalescence_records(tree_sequence_t *self)
{
    return self->num_records;
}

size_t
tree_sequence_get_num_mutations(tree_sequence_t *self)
{
    return self->num_mutations;
}

int WARN_UNUSED
tree_sequence_get_record(tree_sequence_t *self, size_t index,
        coalescence_record_t **record, int order)
{
    int ret = 0;
    size_t j;

    if (index >= self->num_records) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    switch (order) {
        case MSP_ORDER_TIME:
            j = index;
            break;
        case MSP_ORDER_LEFT:
            j = self->trees.indexes.insertion_order[index];
            break;
        case MSP_ORDER_RIGHT:
            j = self->trees.indexes.removal_order[index];
            break;
        default:
            ret = MSP_ERR_BAD_ORDERING;
            goto out;
    }
    self->returned_record.left = self->trees.breakpoints[
        self->trees.records.left[j]];
    self->returned_record.right = self->trees.breakpoints[
        self->trees.records.right[j]];
    self->returned_record.node = self->trees.records.node[j];
    self->returned_record.num_children = self->trees.records.num_children[j];
    self->returned_record.children = self->trees.records.children[j];
    self->returned_record.time = self->trees.nodes.time[self->trees.records.node[j]];
    self->returned_record.population_id = self->trees.nodes.population[
        self->trees.records.node[j]];
    *record = &self->returned_record;
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_get_mutations(tree_sequence_t *self, mutation_t *mutations)
{
    int ret = 0;
    size_t j;

    assert(mutations != NULL);
    for (j = 0; j < self->num_mutations; j++) {
        mutations[j].position = self->mutations.position[j];
        mutations[j].node = self->mutations.node[j];
    }
    return ret;
}

/*
 * This is a convenience short cut for sparse tree alloc in the common
 * case where we're allocating it for a given sequence.
 */
int WARN_UNUSED
tree_sequence_alloc_sparse_tree(tree_sequence_t *self, sparse_tree_t *tree,
        uint32_t *tracked_leaves, uint32_t num_tracked_leaves, int flags)
{
    return sparse_tree_alloc(tree, self->sample_size,
            (uint32_t) self->num_nodes, self->num_mutations, tracked_leaves,
            num_tracked_leaves, flags);
}

int WARN_UNUSED
tree_sequence_set_samples(tree_sequence_t *self, size_t sample_size,
        sample_t *samples)
{
    int ret = MSP_ERR_BAD_SAMPLES;
    uint32_t j;

    if (sample_size != self->sample_size) {
        goto out;
    }
    for (j = 0; j < self->sample_size; j++) {
        self->trees.nodes.population[j] = samples[j].population_id;
        if (samples[j].time < 0) {
            goto out;
        }
        self->trees.nodes.time[j] = samples[j].time;
    }
    ret = 0;
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_set_mutations(tree_sequence_t *self, size_t num_mutations,
        mutation_t *mutations)
{
    int ret = -1;
    size_t j;
    mutation_t **mutation_ptrs = NULL;

    if (self->num_mutations > 0) {
        /* any mutations that were there previously are overwritten. */
        if (self->mutations.node != NULL) {
            free(self->mutations.node);
        }
        if (self->mutations.position != NULL) {
            free(self->mutations.position);
        }
    }
    self->num_mutations = 0;
    self->mutations.position = NULL;
    self->mutations.node = NULL;
    if (num_mutations > 0) {
        /* Allocate the storage we need to keep the mutations. */
        mutation_ptrs = malloc(num_mutations * sizeof(mutation_t *));
        self->mutations.node = malloc(num_mutations * sizeof(uint32_t));
        self->mutations.position = malloc(num_mutations * sizeof(double));
        if (mutation_ptrs == NULL || self->mutations.node == NULL
                || self->mutations.position == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        for (j = 0; j < num_mutations; j++) {
            mutation_ptrs[j] = mutations + j;
            if (mutations[j].position < 0
                    || mutations[j].position > self->sequence_length
                    || mutations[j].node == MSP_NULL_NODE
                    || mutations[j].node >= self->num_nodes) {
                ret = MSP_ERR_BAD_MUTATION;
                goto out;
            }
        }
        /* Mutations are required to be sorted in position order. */
        qsort(mutation_ptrs, num_mutations, sizeof(mutation_t *),
                cmp_mutation_pointer);
        self->num_mutations = num_mutations;
        for (j = 0; j < num_mutations; j++) {
            self->mutations.node[j] = mutation_ptrs[j]->node;
            self->mutations.position[j] = mutation_ptrs[j]->position;
        }
    }
    ret = 0;
out:
    if (mutation_ptrs != NULL) {
        free(mutation_ptrs);
    }
    return ret;
}

/* ======================================================== *
 * Tree diff iterator.
 * ======================================================== */

int WARN_UNUSED
tree_diff_iterator_alloc(tree_diff_iterator_t *self,
        tree_sequence_t *tree_sequence)
{
    int ret = 0;

    assert(tree_sequence != NULL);
    memset(self, 0, sizeof(tree_diff_iterator_t));
    self->sample_size = tree_sequence_get_sample_size(tree_sequence);
    self->num_nodes = tree_sequence_get_num_nodes(tree_sequence);
    self->num_records = tree_sequence_get_num_coalescence_records(
            tree_sequence);
    self->tree_sequence = tree_sequence;
    self->insertion_index = 0;
    self->removal_index = 0;
    self->tree_left = 0;
    /* The maximum number of records is to remove and insert all n - 1
     * records */
    self->node_records = malloc(2 * self->sample_size * sizeof(node_record_t));
    if (self->node_records == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
out:
    return ret;
}

int WARN_UNUSED
tree_diff_iterator_free(tree_diff_iterator_t *self)
{
    int ret = 0;
    if (self->node_records != NULL) {
        free(self->node_records);
    }
    return ret;
}

void
tree_diff_iterator_print_state(tree_diff_iterator_t *self, FILE *out)
{
    fprintf(out, "tree_diff_iterator state\n");
    fprintf(out, "num_records = %d\n", (int) self->num_records);
    fprintf(out, "insertion_index = %d\n", (int) self->insertion_index);
    fprintf(out, "removal_index = %d\n", (int) self->removal_index);
    fprintf(out, "tree_left = %d\n", self->tree_left);
}

int WARN_UNUSED
tree_diff_iterator_next(tree_diff_iterator_t *self, double *length,
        node_record_t **nodes_out, node_record_t **nodes_in)
{
    int ret = 0;
    uint32_t k;
    uint32_t last_left = self->tree_left;
    size_t next_node_record = 0;
    tree_sequence_t *s = self->tree_sequence;
    node_record_t *out_head = NULL;
    node_record_t *out_tail = NULL;
    node_record_t *in_head = NULL;
    node_record_t *in_tail = NULL;
    node_record_t *w = NULL;
    int first_tree = self->insertion_index == 0;
    size_t in_count = 0;
    size_t out_count = 0;

    assert(s != NULL);
    if (self->insertion_index < self->num_records) {
        /* First we remove the stale records */
        out_count = 0;
        while (s->trees.records.right[
                s->trees.indexes.removal_order[self->removal_index]]
                    == self->tree_left) {
            k = s->trees.indexes.removal_order[self->removal_index];
            assert(next_node_record < 2 * self->sample_size);
            w = &self->node_records[next_node_record];
            next_node_record++;
            w->node = s->trees.records.node[k];
            w->time = s->trees.nodes.time[w->node];
            w->num_children = s->trees.records.num_children[k];
            w->children = s->trees.records.children[k];
            w->next = NULL;
            if (out_head == NULL) {
                out_head = w;
                out_tail = w;
            } else {
                out_tail->next = w;
                out_tail = w;
            }
            self->removal_index++;
            out_count += w->num_children - 1;
        }
        /* Now insert the new records */
        in_count = 0;
        while (self->insertion_index < self->num_records &&
                s->trees.records.left[
                    s->trees.indexes.insertion_order[self->insertion_index]]
                        == self->tree_left) {
            k = s->trees.indexes.insertion_order[self->insertion_index];
            assert(next_node_record < 2 * self->sample_size);
            w = &self->node_records[next_node_record];
            next_node_record++;
            w->node = s->trees.records.node[k];
            w->time = s->trees.nodes.time[w->node];
            w->num_children = s->trees.records.num_children[k];
            w->children = s->trees.records.children[k];
            w->next = NULL;
            if (in_head == NULL) {
                in_head = w;
                in_tail = w;
            } else {
                in_tail->next = w;
                in_tail = w;
            }
            self->insertion_index++;
            in_count += w->num_children - 1;
        }
        /* Update the left coordinate */
        self->tree_left = s->trees.records.right[
            s->trees.indexes.removal_order[self->removal_index]];
        ret = 1;
    }
    if (first_tree) {
        if (in_count != self->sample_size - 1 || out_count != 0) {
            ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
            goto out;
        }
    } else {
        if (in_count != out_count) {
            ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
            goto out;
        }
    }
    *nodes_out = out_head;
    *nodes_in = in_head;
    *length = s->trees.breakpoints[self->tree_left]
        - s->trees.breakpoints[last_left];
out:
    return ret;
}

/* ======================================================== *
 * sparse tree
 * ======================================================== */

int WARN_UNUSED
sparse_tree_alloc(sparse_tree_t *self, uint32_t sample_size, uint32_t num_nodes,
        size_t max_mutations, uint32_t *tracked_leaves,
        uint32_t num_tracked_leaves, int flags)
{
    int ret = MSP_ERR_NO_MEMORY;
    uint32_t j, u;
    leaf_list_node_t *w;

    memset(self, 0, sizeof(sparse_tree_t));
    if (num_nodes == 0 || sample_size == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->num_nodes = num_nodes;
    self->sample_size = sample_size;
    self->flags = flags;
    self->parent = malloc(self->num_nodes * sizeof(uint32_t));
    self->population = malloc(self->num_nodes * sizeof(uint8_t));
    self->time = malloc(self->num_nodes * sizeof(double));
    self->num_children = malloc(self->num_nodes * sizeof(uint32_t));
    self->children = malloc(self->num_nodes * sizeof(uint32_t *));
    if (self->time == NULL || self->parent == NULL || self->children == NULL
            || self->num_children == NULL || self->population == NULL) {
        goto out;
    }
    /* the maximum possible height of the tree is n + 1, including
     * the null value. */
    self->stack1 = malloc((self->sample_size + 1) * sizeof(uint32_t));
    self->stack2 = malloc((self->sample_size + 1) * sizeof(uint32_t));
    if (self->stack1 == NULL || self->stack2 == NULL) {
        goto out;
    }
    self->max_mutations = max_mutations;
    self->num_mutations = 0;
    self->mutations = malloc(max_mutations * sizeof(mutation_t));
    if (self->mutations == NULL) {
        goto out;
    }
    if (self->flags & MSP_COUNT_LEAVES) {
        self->num_leaves = calloc(self->num_nodes, sizeof(uint32_t));
        self->num_tracked_leaves = calloc(self->num_nodes, sizeof(uint32_t));
        if (self->num_leaves == NULL || self->num_tracked_leaves == NULL) {
            goto out;
        }
        self->leaf_list_head = calloc(self->num_nodes,
                sizeof(leaf_list_node_t *));
        self->leaf_list_tail = calloc(self->num_nodes,
                sizeof(leaf_list_node_t *));
        self->leaf_list_node_mem = calloc(self->sample_size,
                sizeof(leaf_list_node_t));
        if (self->leaf_list_head == NULL || self->leaf_list_tail == NULL
                || self->leaf_list_node_mem == NULL) {
            goto out;
        }
        for (j = 0; j < self->sample_size; j++) {
            self->num_leaves[j] = 1;
            w = &self->leaf_list_node_mem[j];
            w->next = NULL;
            w->node = j;
            self->leaf_list_head[j] = w;
            self->leaf_list_tail[j] = w;
        }
        for (j = 0; j < num_tracked_leaves; j++) {
            u = tracked_leaves[j];
            if (u >= self->sample_size) {
                ret = MSP_ERR_BAD_PARAM_VALUE;
                goto out;
            }
            self->num_tracked_leaves[u] = 1;
        }
    }
    ret = 0;
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_free(sparse_tree_t *self)
{
    if (self->parent != NULL) {
        free(self->parent);
    }
    if (self->population != NULL) {
        free(self->population);
    }
    if (self->time != NULL) {
        free(self->time);
    }
    if (self->children != NULL) {
        free(self->children);
    }
    if (self->num_children != NULL) {
        free(self->num_children);
    }
    if (self->stack1 != NULL) {
        free(self->stack1);
    }
    if (self->stack2 != NULL) {
        free(self->stack2);
    }
    if (self->mutations != NULL) {
        free(self->mutations);
    }
    if (self->num_leaves != NULL) {
        free(self->num_leaves);
    }
    if (self->num_tracked_leaves != NULL) {
        free(self->num_tracked_leaves);
    }
    if (self->leaf_list_head != NULL) {
        free(self->leaf_list_head);
    }
    if (self->leaf_list_tail != NULL) {
        free(self->leaf_list_tail);
    }
    if (self->leaf_list_node_mem != NULL) {
        free(self->leaf_list_node_mem);
    }
    return 0;
}

int WARN_UNUSED
sparse_tree_clear(sparse_tree_t *self)
{
    int ret = 0;
    size_t N = self->num_nodes;
    size_t n = self->sample_size;

    self->left = 0;
    self->right = 0;
    self->root = 0;
    self->index = UINT32_MAX;
    memset(self->parent, (int) MSP_NULL_NODE, N * sizeof(uint32_t));
    memset(self->population, (int) MSP_NULL_POPULATION_ID, N * sizeof(uint8_t));
    memset(self->time, 0, N * sizeof(double));
    memset(self->num_children, 0, N * sizeof(uint32_t));
    memset(self->children, 0, N * sizeof(uint32_t *));
    if (self->flags & MSP_COUNT_LEAVES) {
        memset(self->num_leaves + n, 0, (N - n) * sizeof(uint32_t));
        memset(self->num_tracked_leaves + n, 0, (N - n) * sizeof(uint32_t));
        memset(self->leaf_list_head + n, 0,
                (N - n) * sizeof(leaf_list_node_t *));
        memset(self->leaf_list_tail + n, 0,
                (N - n) * sizeof(leaf_list_node_t *));
    }
    return ret;
}

int WARN_UNUSED
sparse_tree_get_mrca(sparse_tree_t *self, uint32_t u, uint32_t v,
        uint32_t *mrca)
{
    int ret = 0;
    uint32_t w = 0;
    uint32_t *s1 = self->stack1;
    uint32_t *s2 = self->stack2;
    uint32_t j;
    int l1, l2;

    if (u >= self->num_nodes || v >= self->num_nodes) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    j = u;
    l1 = 0;
    while (j != MSP_NULL_NODE) {
        assert(l1 < (int) self->sample_size);
        s1[l1] = j;
        l1++;
        j = self->parent[j];
    }
    s1[l1] = MSP_NULL_NODE;
    j = v;
    l2 = 0;
    while (j != MSP_NULL_NODE) {
        assert(l2 < (int) self->sample_size);
        s2[l2] = j;
        l2++;
        j = self->parent[j];
    }
    s2[l2] = MSP_NULL_NODE;
    do {
        w = s1[l1];
        l1--;
        l2--;
    } while (l1 >= 0 && l2 >= 0 && s1[l1] == s2[l2]);
    *mrca = w;
    ret = 0;
out:
    return ret;
}

static int
sparse_tree_check_node(sparse_tree_t *self, uint32_t u)
{
    int ret = 0;
    if (u > self->num_nodes) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
    }
    return ret;
}

static int
sparse_tree_get_num_leaves_by_traversal(sparse_tree_t *self, uint32_t u,
        uint32_t *num_leaves)
{
    int ret = 0;
    uint32_t *stack = self->stack1;
    uint32_t v, c;
    uint32_t count = 0;
    int stack_top = 0;

    stack[0] = u;
    while (stack_top >= 0) {
        v = stack[stack_top];
        stack_top--;
        if (v < self->sample_size) {
            count++;
        }
        for (c = 0; c < self->num_children[v]; c++) {
            stack_top++;
            stack[stack_top] = self->children[v][c];
        }
    }
    *num_leaves = count;
    return ret;
}

int WARN_UNUSED
sparse_tree_get_num_leaves(sparse_tree_t *self, uint32_t u,
        uint32_t *num_leaves)
{
    int ret = 0;

    ret = sparse_tree_check_node(self, u);
    if (ret != 0) {
        goto out;
    }

    if (self->flags & MSP_COUNT_LEAVES) {
        *num_leaves = self->num_leaves[u];
    } else {
        ret = sparse_tree_get_num_leaves_by_traversal(self, u, num_leaves);
    }
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_get_num_tracked_leaves(sparse_tree_t *self, uint32_t u,
        uint32_t *num_tracked_leaves)
{
    int ret = 0;

    ret = sparse_tree_check_node(self, u);
    if (ret != 0) {
        goto out;
    }
    if (! (self->flags & MSP_COUNT_LEAVES)) {
        ret = MSP_ERR_UNSUPPORTED_OPERATION;
        goto out;
    }
    *num_tracked_leaves = self->num_tracked_leaves[u];
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_get_leaf_list(sparse_tree_t *self, uint32_t u,
        leaf_list_node_t **head, leaf_list_node_t **tail)
{
    int ret = 0;

    ret = sparse_tree_check_node(self, u);
    if (ret != 0) {
        goto out;
    }
    if (! (self->flags & MSP_COUNT_LEAVES)) {
        ret = MSP_ERR_UNSUPPORTED_OPERATION;
        goto out;
    }
    if (self->leaf_list_head[u] == NULL ||
            self->leaf_list_tail[u] == NULL) {
        /* This sigifies that we're trying to get the leaf list
         * for a node that is not in the current tree.
         */
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    *head = self->leaf_list_head[u];
    *tail = self->leaf_list_tail[u];
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_get_root(sparse_tree_t *self, uint32_t *root)
{
    *root = self->root;
    return 0;
}


int WARN_UNUSED
sparse_tree_get_parent(sparse_tree_t *self, uint32_t u, uint32_t *parent)
{
    int ret = 0;

    ret = sparse_tree_check_node(self, u);
    if (ret != 0) {
        goto out;
    }
    *parent = self->parent[u];
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_get_time(sparse_tree_t *self, uint32_t u, double *t)
{
    int ret = 0;

    ret = sparse_tree_check_node(self, u);
    if (ret != 0) {
        goto out;
    }
    *t = self->time[u];
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_get_children(sparse_tree_t *self, uint32_t u,
        uint32_t *num_children, uint32_t **children)
{
    int ret = 0;

    ret = sparse_tree_check_node(self, u);
    if (ret != 0) {
        goto out;
    }
    *num_children = self->num_children[u];
    *children = self->children[u];
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_get_mutations(sparse_tree_t *self, size_t *num_mutations,
        mutation_t **mutations)
{
    *mutations = self->mutations;
    *num_mutations = self->num_mutations;
    return 0;
}

/* ======================================================== *
 * sparse tree iterator
 * ======================================================== */

int WARN_UNUSED
sparse_tree_iterator_alloc(sparse_tree_iterator_t *self,
        tree_sequence_t *tree_sequence, sparse_tree_t *tree)
{
    int ret = MSP_ERR_NO_MEMORY;
    uint32_t j;

    assert(tree_sequence != NULL);
    assert(tree != NULL);
    assert(tree->time != NULL && tree->parent != NULL
            && tree->children != NULL);
    if (tree_sequence_get_num_nodes(tree_sequence) != tree->num_nodes ||
            tree_sequence_get_sample_size(tree_sequence)
                != tree->sample_size ||
            tree_sequence_get_num_mutations(tree_sequence)
                != tree->max_mutations) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    memset(self, 0, sizeof(sparse_tree_iterator_t));
    self->sample_size = tree_sequence_get_sample_size(tree_sequence);
    self->num_nodes = tree_sequence_get_num_nodes(tree_sequence);
    self->num_records = tree_sequence_get_num_coalescence_records(
            tree_sequence);
    self->tree_sequence = tree_sequence;
    self->tree = tree;
    self->tree->sample_size = self->sample_size;
    self->insertion_index = 0;
    self->removal_index = 0;
    self->mutation_index = 0;
    ret = sparse_tree_clear(self->tree);
    if (ret != 0) {
        goto out;
    }
    /* Set the sample attributes */
    for (j = 0; j < self->sample_size; j++) {
        self->tree->population[j] = self->tree_sequence->trees.nodes.population[j];
        self->tree->time[j] = self->tree_sequence->trees.nodes.time[j];
    }
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_iterator_free(sparse_tree_iterator_t *self)
{
    int ret = 0;
    return ret;
}

static void
sparse_tree_iterator_check_state(sparse_tree_iterator_t *self)
{
    uint32_t u, v, j, k, num_leaves;
    int err, found;

    assert(self->tree->num_nodes == self->num_nodes);
    for (j = 0; j < self->sample_size; j++) {
        u = j;
        assert(self->tree->time[u] >= 0.0);
        assert(self->tree->num_children[j] == 0);
        while (self->tree->parent[u] != MSP_NULL_NODE) {
            v = self->tree->parent[u];
            found = 0;
            for (k = 0; k < self->tree->num_children[v]; k++) {
                if (self->tree->children[v][k] == u) {
                    found = 1;
                }
            }
            assert(found);
            u = v;
            assert(self->tree->time[u] > 0.0);
        }
        assert(u == self->tree->root);
    }
    if (self->tree->flags & MSP_COUNT_LEAVES) {
        for (j = 0; j < self->num_nodes; j++) {
            err = sparse_tree_get_num_leaves_by_traversal(self->tree, j,
                    &num_leaves);
            assert(err == 0);
            assert(num_leaves == self->tree->num_leaves[j]);
        }
    }
}

void
sparse_tree_iterator_print_state(sparse_tree_iterator_t *self, FILE *out)
{
    size_t j, k;
    leaf_list_node_t *u;

    fprintf(out, "sparse_tree_iterator state\n");
    fprintf(out, "insertion_index = %d\n", (int) self->insertion_index);
    fprintf(out, "removal_index = %d\n", (int) self->removal_index);
    fprintf(out, "mutation_index = %d\n", (int) self->mutation_index);
    fprintf(out, "num_records = %d\n", (int) self->num_records);
    fprintf(out, "tree.flags = %d\n", self->tree->flags);
    fprintf(out, "tree.left = %f\n", self->tree->left);
    fprintf(out, "tree.right = %f\n", self->tree->right);
    fprintf(out, "tree.root = %d\n", self->tree->root);
    fprintf(out, "tree.index = %d\n", self->tree->index);
    for (j = 0; j < self->tree->num_nodes; j++) {
        fprintf(out, "\t%d\t%d\t%f\t%d\t(", (int) j, self->tree->parent[j],
            self->tree->time[j], self->tree->population[j]);
        for (k = 0; k < self->tree->num_children[j]; k++) {
            fprintf(out, "%d", self->tree->children[j][k]);
            if (k < self->tree->num_children[j] - 1) {
                fprintf(out, ", ");
            }
        }
        fprintf(out, ")");
        if (self->tree->flags & MSP_COUNT_LEAVES) {
            fprintf(out, "\t%d\t%d", self->tree->num_leaves[j],
                    self->tree->num_tracked_leaves[j]);
            fprintf(out, "\t[");
            u = self->tree->leaf_list_head[j];
            if (u != NULL) {
                while (1) {
                    fprintf(out, "%d ", u->node);
                    if (u == self->tree->leaf_list_tail[j]) {
                        break;
                    }
                    u = u->next;
                }
            } else {
                assert(self->tree->leaf_list_tail[j] == NULL);
            }

            fprintf(out, "]");
        }
        fprintf(out, "\n");
    }
    fprintf(out, "mutations = \n");
    for (j = 0; j < self->tree->num_mutations; j++) {
        fprintf(out, "\t%d @ %f\n", self->tree->mutations[j].node,
                self->tree->mutations[j].position);
    }
    sparse_tree_iterator_check_state(self);
}

static inline void
sparse_tree_iterator_propagate_leaf_loss(sparse_tree_iterator_t *self, uint32_t u)
{
    sparse_tree_t *t = self->tree;
    uint32_t all_leaves_diff = t->num_leaves[u];
    uint32_t tracked_leaves_diff = t->num_tracked_leaves[u];
    uint32_t v = u;

    t->leaf_list_head[v] = NULL;
    t->leaf_list_tail[v] = NULL;

    /* propogate this loss up as far as we can */
    while (v != MSP_NULL_NODE) {
        t->num_leaves[v] -= all_leaves_diff;
        t->num_tracked_leaves[v] -= tracked_leaves_diff;
        v = t->parent[v];
    }
}

/* Returns the index of child u in the parent v */
static inline uint32_t
sparse_tree_iterator_get_child_index(sparse_tree_iterator_t *self, uint32_t v,
        uint32_t u)
{
    uint32_t j;
    sparse_tree_t *t = self->tree;

    for (j = 0; j < t->num_children[v] && t->children[v][j] != u; j++);
    assert(j != t->num_children[v]);
    return j;
}

static inline void
sparse_tree_iterator_propagate_leaf_gain(sparse_tree_iterator_t *self, uint32_t u)
{
    sparse_tree_t *t = self->tree;
    uint32_t j, k, v, w, *c;
    uint32_t all_leaves_diff = 0;
    uint32_t tracked_leaves_diff = 0;

    c = t->children[u];
    k = t->num_children[u];
    for (j = 0; j < k; j++) {
        all_leaves_diff += t->num_leaves[c[j]];
        tracked_leaves_diff += t->num_tracked_leaves[c[j]];
    }
    /* propogate this gain up as far as we can */
    v = u;
    while (v != MSP_NULL_NODE) {
        t->num_leaves[v] += all_leaves_diff;
        t->num_tracked_leaves[v] += tracked_leaves_diff;
        v = t->parent[v];
    }

    /* LEAF_LISTS */
    for (j = 1; j < k; j++) {
        assert(t->leaf_list_tail[c[j]] != NULL);
        t->leaf_list_tail[c[j - 1]]->next = t->leaf_list_head[c[j]];
    }
    t->leaf_list_head[u] = t->leaf_list_head[c[0]];
    t->leaf_list_tail[u] = t->leaf_list_tail[c[k - 1]];

    v = u;
    w = t->parent[v];
    while (w != MSP_NULL_NODE) {
        j = sparse_tree_iterator_get_child_index(self, w, v);
        if (j != 0) {
            break;
        }
        t->leaf_list_head[w] = t->leaf_list_head[u];
        v = w;
        w = t->parent[w];
    }

    v = u;
    w = t->parent[v];
    while (w != MSP_NULL_NODE) {
        j = sparse_tree_iterator_get_child_index(self, w, v);
        if (j != t->num_children[w] - 1) {
            break;
        }
        t->leaf_list_tail[w] = t->leaf_list_tail[u];
        v = w;
        w = t->parent[w];
    }
}

static inline void
sparse_tree_iterator_post_propagate_leaf_gain(sparse_tree_iterator_t *self, uint32_t u)
{
    uint32_t v, w, j;
    sparse_tree_t *t = self->tree;

    v = u;
    w = t->parent[v];
    while (w != MSP_NULL_NODE) {
        j = sparse_tree_iterator_get_child_index(self, w, v);
        if (j > 0) {
            t->leaf_list_tail[t->children[w][j - 1]]->next = t->leaf_list_head[v];
        }
        if (j < t->num_children[w] - 1) {
            t->leaf_list_tail[v]->next = t->leaf_list_head[t->children[w][j + 1]];
        }
        v = w;
        w = t->parent[w];
    }
}

int WARN_UNUSED
sparse_tree_iterator_next(sparse_tree_iterator_t *self)
{
    int ret = 0;
    uint32_t j, k, u, v, oldest_child;
    size_t h;
    double oldest_child_time;
    tree_sequence_t *s = self->tree_sequence;
    sparse_tree_t *t = self->tree;
    int first_tree = self->insertion_index == 0;
    /* To detect errors we look for situations when the trees are not
     * completed. The propery is that we should have
     * n - 1 = sum(k - 1) where k is the arity of the node over the
     * whole tree.
     */
    size_t in_count = 0;
    size_t out_count = 0;

    assert(t != NULL && s != NULL);
    if (self->insertion_index < self->num_records) {
        /* First we remove the stale records */
        while (s->trees.breakpoints[s->trees.records.right[
                s->trees.indexes.removal_order[self->removal_index]]]
                    == t->right) {
            k = s->trees.indexes.removal_order[self->removal_index];
            u = s->trees.records.node[k];
            out_count += t->num_children[u] - 1;
            oldest_child_time = -1;
            oldest_child = 0;
            for (j = 0; j < t->num_children[u]; j++) {
                t->parent[t->children[u][j]] = MSP_NULL_NODE;
                if (t->time[t->children[u][j]] > oldest_child_time) {
                    oldest_child = t->children[u][j];
                    oldest_child_time = t->time[t->children[u][j]];
                }
            }
            t->num_children[u] = 0;
            t->children[u] = NULL;
            t->time[u] = 0;
            t->population[u] = MSP_NULL_POPULATION_ID;
            if (u == t->root) {
                t->root = oldest_child;
            }
            self->removal_index++;
            if (t->flags & MSP_COUNT_LEAVES) {
                sparse_tree_iterator_propagate_leaf_loss(self, u);
            }
        }
        /* Update the interval */
        t->left = t->right;
        t->right = s->trees.breakpoints[s->trees.records.right[
            s->trees.indexes.removal_order[self->removal_index]]];

        /* Now insert the new records */
        h = self->insertion_index;
        while (self->insertion_index < self->num_records &&
                s->trees.breakpoints[s->trees.records.left[
                    s->trees.indexes.insertion_order[self->insertion_index]]]
                        == t->left) {
            k = s->trees.indexes.insertion_order[self->insertion_index];
            u = s->trees.records.node[k];
            for (j = 0; j < s->trees.records.num_children[k]; j++) {
                t->parent[s->trees.records.children[k][j]] = u;
            }
            t->num_children[u] = s->trees.records.num_children[k];
            t->children[u] = s->trees.records.children[k];
            in_count += t->num_children[u] - 1;
            t->time[u] = s->trees.nodes.time[u];
            t->population[u] = s->trees.nodes.population[u];
            if (t->time[u] > t->time[t->root]) {
                t->root = u;
            }
            self->insertion_index++;
            if (t->flags & MSP_COUNT_LEAVES) {
                sparse_tree_iterator_propagate_leaf_gain(self, u);
            }
        }
        if (t->flags & MSP_COUNT_LEAVES) {
            /* TODO this should only be done on LEAF_LISTS */
            while (h < self->num_records &&
                    s->trees.breakpoints[s->trees.records.left[
                        s->trees.indexes.insertion_order[h]]] == t->left) {
                k = s->trees.indexes.insertion_order[h];
                u = s->trees.records.node[k];
                sparse_tree_iterator_post_propagate_leaf_gain(self, u);
                h++;
            }
        }
        /* Check for errors. */
        if (first_tree) {
            if (out_count != 0 || in_count != self->sample_size - 1) {
                ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
                goto out;
            }
        } else {
            if (in_count != out_count) {
                ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
                goto out;
            }
        }
        /* In very rare situations, we have to traverse upwards to find the
         * new root.
         */
        while (t->parent[t->root] != MSP_NULL_NODE) {
            t->root = t->parent[t->root];
        }
        /* now update the mutations */
        t->num_mutations = 0;
        while (self->mutation_index < s->num_mutations
                && s->mutations.position[self->mutation_index] < t->right) {
            assert(t->num_mutations < t->max_mutations);
            /* Throw an error if the mutation is for a node not in the
             * tree or a root.
             * */
            v = s->mutations.node[self->mutation_index];
            if (t->parent[v] == MSP_NULL_NODE) {
                ret = MSP_ERR_BAD_MUTATION;
                goto out;
            }
            t->mutations[t->num_mutations].node = v;
            t->mutations[t->num_mutations].position =
                s->mutations.position[self->mutation_index];
            self->mutation_index++;
            t->num_mutations++;
        }
        /* Finally, update the tree index and indicate we have a valid tree */
        t->index++;
        ret = 1;
    }
out:
    return ret;
}
