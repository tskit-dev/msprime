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

#define MSP_DIR_FORWARD 1
#define MSP_DIR_REVERSE -1

typedef struct {
    double value;
    uint32_t index;
    int64_t time;
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
    if (ret == 0) {
        ret = (ca->time > cb->time) - (ca->time < cb->time);
    }
    return ret;
}

static void
tree_sequence_check_state(tree_sequence_t *self)
{
    size_t j;

    for (j = 0; j < self->num_records; j++) {
        assert(self->trees.records.num_children[j] >= 2);
    }
}

void
tree_sequence_print_state(tree_sequence_t *self, FILE *out)
{
    size_t j, k;

    fprintf(out, "tree_sequence state\n");
    fprintf(out, "refcount = %d\n", self->refcount);
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
    if (self->num_mutations > 0) {
        fprintf(out, "tree_mutations\n");
        for (j = 0; j < self->trees.num_breakpoints; j++) {
            fprintf(out, "\ttree %d\t%f\n", (int) j, self->trees.breakpoints[j]);
            for (k = 0; k < self->mutations.num_tree_mutations[j]; k++) {
                fprintf(out, "\t\t%f\t%d\n",
                        self->mutations.tree_mutations[j][k].position,
                        self->mutations.tree_mutations[j][k].node);
            }
        }
    }
    tree_sequence_check_state(self);
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
    if (self->mutations.tree_mutations_mem != NULL) {
        free(self->mutations.tree_mutations_mem);
    }
    if (self->mutations.tree_mutations != NULL) {
        free(self->mutations.tree_mutations);
    }
    if (self->mutations.num_tree_mutations != NULL) {
        free(self->mutations.num_tree_mutations);
    }
    return 0;
}

int
tree_sequence_increment_refcount(tree_sequence_t *self)
{
    /* TODO make this threadsafe. */
    self->refcount++;
    return 0;
}

int
tree_sequence_decrement_refcount(tree_sequence_t *self)
{
    /* TODO make this threadsafe. */
    self->refcount--;
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
    p = realloc(self->provenance_strings,
            (self->num_provenance_strings + 1) * sizeof(char *));
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
    self->provenance_strings[self->num_provenance_strings] = s;
    self->num_provenance_strings++;
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
        /* When comparing equal left values, we sort by time. Since we require
         * that records are provided in sorted order, the index can be
         * taken as a proxy for time. This avoids issues unstable sort
         * algorithms when multiple events occur at the same time. We are
         * actually making the stronger requirement that records must be
         * provided *in the order they happened*, not just in increasing
         * time. See also the removal order index below.
         */
        sort_buff[j].time = (int64_t ) j;
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
        sort_buff[j].time = -1 * (int64_t ) j;
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

/* Sets up the memory for the mutations associated with each tree.
 */
static int
tree_sequence_init_tree_mutations(tree_sequence_t *self)
{
    int ret = MSP_ERR_GENERIC;
    size_t j, tree_index;
    mutation_t *mut;

    self->mutations.tree_mutations_mem = malloc(
            self->num_mutations * sizeof(mutation_t));
    self->mutations.tree_mutations = calloc(self->trees.num_breakpoints,
            sizeof(mutation_t *));
    self->mutations.num_tree_mutations = calloc(self->trees.num_breakpoints,
            sizeof(size_t));
    if (self->mutations.tree_mutations_mem == NULL ||
            self->mutations.tree_mutations == NULL ||
            self->mutations.num_tree_mutations == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    tree_index = 0;
    self->mutations.tree_mutations[0] = self->mutations.tree_mutations_mem;
    for (j = 0; j < self->num_mutations; j++) {
        mut = &self->mutations.tree_mutations_mem[j];
        mut->index = j;
        mut->position = self->mutations.position[j];
        mut->node = self->mutations.node[j];
        assert(tree_index < self->trees.num_breakpoints - 1);
        while (mut->position >= self->trees.breakpoints[tree_index + 1]) {
            tree_index++;
            self->mutations.tree_mutations[tree_index] = mut;
        }
        self->mutations.num_tree_mutations[tree_index]++;
    }
    ret = 0;
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
    uint32_t version[2];

    attr_id = H5Aopen_by_name(file_id, "/", "format_version",
            H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id < 0) {
        goto out;
    }
    dataspace_id = H5Aget_space(attr_id);
    if (dataspace_id < 0) {
        goto out;
    }
    rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 1) {
        ret = MSP_ERR_FILE_FORMAT;
        goto out;
    }
    status = H5Sget_simple_extent_dims(dataspace_id, &dims, NULL);
    if (status < 0) {
        goto out;
    }
    if (dims != 2) {
        ret = MSP_ERR_FILE_FORMAT;
        goto out;
    }
    status = H5Aread(attr_id, H5T_NATIVE_UINT32, version);
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

    /* Sanity check */
    if (version[0] < MSP_FILE_FORMAT_VERSION_MAJOR) {
        ret = MSP_ERR_FILE_VERSION_TOO_OLD;
        goto out;
    }
    if (version[0] > MSP_FILE_FORMAT_VERSION_MAJOR) {
        ret = MSP_ERR_FILE_VERSION_TOO_NEW;
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
        {"/trees/nodes/population", 1, self->num_nodes, 1},
        {"/trees/nodes/time", 1, self->num_nodes, 1},
        {"/trees/breakpoints", 1, self->trees.num_breakpoints, 1},
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
        {"/trees/nodes/population", H5T_NATIVE_UINT8, 0, 1,
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
    }
    status = H5Tclose(vlen_str);
    if (status < 0) {
        goto out;
    }
    /* Now update the children vectors and find sample size */
    self->sample_size = UINT32_MAX;
    offset = 0;
    for (j = 0; j < self->num_records; j++) {
        assert(offset < self->num_child_nodes);
        self->trees.records.children[j] =
            &self->trees.records.children_mem[offset];
        offset += self->trees.records.num_children[j];
        self->sample_size = GSL_MIN(self->sample_size,
                self->trees.records.node[j]);
    }
    self->sequence_length = self->trees.breakpoints[
        self->trees.num_breakpoints - 1];
    ret = tree_sequence_init_tree_mutations(self);
    if (ret != 0) {
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
    hid_t file_id = -1;

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
    ret = tree_sequence_check(self);
out:
    if (file_id >= 0) {
        status = H5Fclose(file_id);
        if (status < 0) {
            ret = MSP_ERR_HDF5;
        }
    }
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
    /* We need to use separate types for storage and memory here because
     * we seem to get a memory leak in HDF5 otherwise.*/
    hid_t filetype_str = -1;
    hid_t memtype_str = -1;

    filetype_str = H5Tcopy(H5T_C_S1);
    if (filetype_str < 0) {
        goto out;
    }
    status = H5Tset_size(filetype_str, H5T_VARIABLE);
    if (status < 0) {
        goto out;
    }
    memtype_str = H5Tcopy(H5T_C_S1);
    if (memtype_str < 0) {
        goto out;
    }
    status = H5Tset_size(memtype_str, H5T_VARIABLE);
    if (status < 0) {
        goto out;
    }
    fields[0].storage_type = filetype_str;
    fields[0].memory_type = memtype_str;

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
                    fields[j].memory_type != memtype_str) {
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
            status = H5Sclose(dataspace_id);
            if (status < 0) {
                goto out;
            }
        }
    }
    ret = 0;
out:
    if (filetype_str != -1) {
        status = H5Tclose(filetype_str);
        if (status < 0) {
            ret = MSP_ERR_HDF5;
        }
    }
    if (memtype_str != -1) {
        status = H5Tclose(memtype_str);
        if (status < 0) {
            ret = MSP_ERR_HDF5;
        }
    }
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
    uint32_t unused_value = 0;

    struct _hdf5_metadata_write {
        const char *name;
        hid_t parent;
        hid_t storage_type;
        hid_t memory_type;
        size_t size;
        void *source;
    };
    struct _hdf5_metadata_write fields[] = {
        {"format_version", 0, H5T_STD_U32LE, H5T_NATIVE_UINT32, 2, version},
        /* These two attributes are vestigial, and are only included to allow
         * older versions of msprime give a better error condition when confronted
         * with a newer file format. Due to a bug in the way that these attributes
         * we loaded, versions of msprime pre 0.4.0 would complain about a missing
         * attribute rather than giving a File format error. These attributes
         * should be removed in a later version of the file format once we can be
         * fairly sure that these old versions of msprime are no longer around.
         */
        {"sample_size", 0, H5T_STD_U32LE, H5T_NATIVE_UINT32, 1, &unused_value},
        {"sequence_length", 0, H5T_IEEE_F64LE, H5T_NATIVE_UINT32, 1, &unused_value},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _hdf5_metadata_write);
    size_t j;

    for (j = 0; j < num_fields; j++) {
        dims = fields[j].size;
        dataspace_id = H5Screate_simple(1, &dims, NULL);
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
    hid_t file_id = -1;

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
    ret = 0;
out:
    if (file_id > 0) {
        status = H5Fclose(file_id);
        if (status < 0) {
            ret = MSP_ERR_HDF5;
        }
    }
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
    double result, denom, count;

    tree = malloc(sizeof(sparse_tree_t));
    if (tree == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = sparse_tree_alloc(tree, self, MSP_LEAF_COUNTS);
    if (ret != 0) {
        goto out;
    }
    ret = sparse_tree_set_tracked_leaves(tree, num_samples, samples);
    if (ret != 0) {
        goto out;
    }
    /* Allocation done; move onto main algorithm. */
    denom = (num_samples * ((double) num_samples - 1)) / 2.0;
    result = 0.0;
    for (ret = sparse_tree_first(tree); ret == 1;
            ret = sparse_tree_next(tree)) {
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

size_t
tree_sequence_get_num_trees(tree_sequence_t *self)
{
    return self->trees.num_breakpoints - 1;
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
tree_sequence_get_mutations(tree_sequence_t *self, mutation_t **mutations)
{
    *mutations = self->mutations.tree_mutations_mem;
    return 0;
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

    /* TODO make this threadsafe! */
    if (self->refcount != 0) {
        ret = MSP_ERR_REFCOUNT_NONZERO;
        goto out;
    }
    if (self->num_mutations > 0) {
        /* any mutations that were there previously are overwritten. */
        if (self->mutations.node != NULL) {
            free(self->mutations.node);
            self->mutations.node = NULL;
        }
        if (self->mutations.position != NULL) {
            free(self->mutations.position);
            self->mutations.position = NULL;
        }
        if (self->mutations.tree_mutations_mem != NULL) {
            free(self->mutations.tree_mutations_mem);
            self->mutations.tree_mutations_mem = NULL;
        }
        if (self->mutations.tree_mutations != NULL) {
            free(self->mutations.tree_mutations);
            self->mutations.tree_mutations = NULL;
        }
        if (self->mutations.num_tree_mutations != NULL) {
            free(self->mutations.num_tree_mutations);
            self->mutations.num_tree_mutations = NULL;
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
        ret = tree_sequence_init_tree_mutations(self);
        if (ret != 0) {
            goto out;
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

static int WARN_UNUSED
sparse_tree_clear(sparse_tree_t *self)
{
    int ret = 0;
    size_t N = self->num_nodes;
    size_t n = self->sample_size;

    self->left = 0;
    self->left_breakpoint = 0;
    self->right = 0;
    self->right_breakpoint = 0;
    self->root = 0;
    self->index = (size_t) -1;
    memset(self->parent, (int) MSP_NULL_NODE, N * sizeof(uint32_t));
    memset(self->population + n, (int) MSP_NULL_POPULATION_ID,
            (N - n) * sizeof(uint8_t));
    memset(self->time + n, 0, (N - n) * sizeof(double));
    memset(self->num_children + n, 0, (N - n) * sizeof(uint32_t));
    memset(self->children + n, 0, (N - n) * sizeof(uint32_t *));
    if (self->flags & MSP_LEAF_COUNTS) {
        memset(self->num_leaves + n, 0, (N - n) * sizeof(uint32_t));
        memset(self->num_tracked_leaves + n, 0, (N - n) * sizeof(uint32_t));
        memset(self->marked, 0, N * sizeof(uint8_t));
    }
    if (self->flags & MSP_LEAF_LISTS) {
        memset(self->leaf_list_head + n, 0,
                (N - n) * sizeof(leaf_list_node_t *));
        memset(self->leaf_list_tail + n, 0,
                (N - n) * sizeof(leaf_list_node_t *));
    }
    return ret;
}

int WARN_UNUSED
sparse_tree_alloc(sparse_tree_t *self, tree_sequence_t *tree_sequence, int flags)
{
    int ret = MSP_ERR_NO_MEMORY;
    uint32_t j, sample_size;
    size_t num_nodes;
    leaf_list_node_t *w;

    memset(self, 0, sizeof(sparse_tree_t));
    if (tree_sequence == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    num_nodes = tree_sequence->num_nodes;
    sample_size = tree_sequence->sample_size;
    self->num_nodes = (uint32_t) num_nodes;
    self->sample_size = sample_size;
    self->tree_sequence = tree_sequence;
    ret = tree_sequence_increment_refcount(self->tree_sequence);
    if (ret != 0) {
        goto out;
    }
    self->flags = flags;
    self->parent = malloc(num_nodes * sizeof(uint32_t));
    self->population = malloc(num_nodes * sizeof(uint8_t));
    self->time = malloc(num_nodes * sizeof(double));
    self->num_children = malloc(num_nodes * sizeof(uint32_t));
    self->children = malloc(num_nodes * sizeof(uint32_t *));
    if (self->time == NULL || self->parent == NULL || self->children == NULL
            || self->num_children == NULL || self->population == NULL) {
        goto out;
    }
    /* the maximum possible height of the tree is n + 1, including
     * the null value. */
    self->stack1 = malloc((sample_size + 1) * sizeof(uint32_t));
    self->stack2 = malloc((sample_size + 1) * sizeof(uint32_t));
    if (self->stack1 == NULL || self->stack2 == NULL) {
        goto out;
    }
    if (self->flags & MSP_LEAF_COUNTS) {
        self->num_leaves = calloc(num_nodes, sizeof(uint32_t));
        self->num_tracked_leaves = calloc(num_nodes, sizeof(uint32_t));
        self->marked = calloc(num_nodes, sizeof(uint8_t));
        if (self->num_leaves == NULL || self->num_tracked_leaves == NULL
                || self->marked == NULL) {
            goto out;
        }
        for (j = 0; j < sample_size; j++) {
            self->num_leaves[j] = 1;
        }
    }
    if (self->flags & MSP_LEAF_LISTS) {
        self->leaf_list_head = calloc(num_nodes, sizeof(leaf_list_node_t *));
        self->leaf_list_tail = calloc(num_nodes, sizeof(leaf_list_node_t *));
        self->leaf_list_node_mem = calloc(sample_size,
                sizeof(leaf_list_node_t));
        if (self->leaf_list_head == NULL || self->leaf_list_tail == NULL
                || self->leaf_list_node_mem == NULL) {
            goto out;
        }
        for (j = 0; j < sample_size; j++) {
            w = &self->leaf_list_node_mem[j];
            w->next = NULL;
            w->node = j;
            self->leaf_list_head[j] = w;
            self->leaf_list_tail[j] = w;
        }
    }
    /* Set the sample attributes */
    for (j = 0; j < self->sample_size; j++) {
        self->population[j] = self->tree_sequence->trees.nodes.population[j];
        self->time[j] = self->tree_sequence->trees.nodes.time[j];
        self->children[j] = NULL;
        self->num_children[j] = 0;
    }
    ret = sparse_tree_clear(self);
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_free(sparse_tree_t *self)
{
    if (self->tree_sequence != NULL) {
        tree_sequence_decrement_refcount(self->tree_sequence);
    }
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
    if (self->num_leaves != NULL) {
        free(self->num_leaves);
    }
    if (self->num_tracked_leaves != NULL) {
        free(self->num_tracked_leaves);
    }
    if (self->marked != NULL) {
        free(self->marked);
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

static int WARN_UNUSED
sparse_tree_reset_tracked_leaves(sparse_tree_t *self)
{
    int ret = 0;

    if (!(self->flags & MSP_LEAF_COUNTS)) {
        ret = MSP_ERR_UNSUPPORTED_OPERATION;
        goto out;
    }
    memset(self->num_tracked_leaves, 0, self->num_nodes * sizeof(uint32_t));
out:
    return ret;
}


int WARN_UNUSED
sparse_tree_set_tracked_leaves(sparse_tree_t *self, uint32_t num_tracked_leaves,
        uint32_t *tracked_leaves)
{
    int ret = MSP_ERR_GENERIC;
    uint32_t j, u;

    /* TODO This is not needed when the sparse tree is new. We should use the
     * state machine to check and only reset the tracked leaves when needed.
     */
    ret = sparse_tree_reset_tracked_leaves(self);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < num_tracked_leaves; j++) {
        u = tracked_leaves[j];
        if (u >= self->sample_size) {
            ret = MSP_ERR_OUT_OF_BOUNDS;
            goto out;
        }
        if (self->num_tracked_leaves[u] != 0) {
            ret = MSP_ERR_DUPLICATE_TRACKED_LEAF;
            goto out;
        }
        /* Propagate this upwards */
        while (u != MSP_NULL_NODE) {
            self->num_tracked_leaves[u] += 1;
            u = self->parent[u];
        }
    }
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_set_tracked_leaves_from_leaf_list(sparse_tree_t *self,
        leaf_list_node_t *head, leaf_list_node_t *tail)
{
    int ret = MSP_ERR_GENERIC;
    leaf_list_node_t *list_node = head;
    uint32_t u;
    int not_done;

    if (head == NULL || tail == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    /* TODO This is not needed when the sparse tree is new. We should use the
     * state machine to check and only reset the tracked leaves when needed.
     */
    ret = sparse_tree_reset_tracked_leaves(self);
    if (ret != 0) {
        goto out;
    }
    not_done = 1;
    while (not_done) {
        u = list_node->node;
        /* Propagate this upwards */
        assert(self->num_tracked_leaves[u] == 0);
        while (u != MSP_NULL_NODE) {
            self->num_tracked_leaves[u] += 1;
            u = self->parent[u];
        }
        not_done = list_node != tail;
        list_node = list_node->next;
    }
out:
    return ret;
}


int WARN_UNUSED
sparse_tree_copy(sparse_tree_t *self, sparse_tree_t *source)
{
    int ret = MSP_ERR_GENERIC;
    size_t N = self->num_nodes;
    size_t n = self->sample_size;

    if (self == source) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (self->tree_sequence != source->tree_sequence) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->left = source->left;
    self->left_breakpoint = source->left_breakpoint;
    self->right = source->right;
    self->right_breakpoint = source->right_breakpoint;
    self->root = source->root;
    self->index = source->index;
    self->num_mutations = source->num_mutations;
    self->mutations = source->mutations;

    memcpy(self->parent, source->parent, N * sizeof(uint32_t));
    memcpy(self->population, source->population, N * sizeof(uint8_t));
    memcpy(self->time, source->time, N * sizeof(double));
    memcpy(self->num_children, source->num_children, N * sizeof(uint32_t));
    memcpy(self->children, source->children, N * sizeof(uint32_t *));
    if (self->flags & MSP_LEAF_COUNTS) {
        if (! (source->flags & MSP_LEAF_COUNTS)) {
            ret = MSP_ERR_UNSUPPORTED_OPERATION;
            goto out;
        }
        memcpy(self->num_leaves + n, source->num_leaves + n,
                (N - n) * sizeof(uint32_t));
    }
    if (self->flags & MSP_LEAF_LISTS) {
        ret = MSP_ERR_UNSUPPORTED_OPERATION;
        goto out;
    }
    ret = 0;
out:
    return ret;
}

/* Returns 0 if the specified sparse trees are equal, 1 if they are
 * not equal, and < 0 if an error occurs.
 *
 * We only consider topological properties of the tree. Optional
 * counts and leaf lists are not considered for equality.
 */
int WARN_UNUSED
sparse_tree_equal(sparse_tree_t *self, sparse_tree_t *other)
{
    int ret = 1;
    int condition;
    size_t N = self->num_nodes;

    if (self->tree_sequence != other->tree_sequence) {
        /* It is an error to compare trees from different tree sequences. */
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    condition = self->index == other->index
        && self->left_breakpoint == other->left_breakpoint
        && self->right_breakpoint == other->right_breakpoint
        && self->root == other->root
        && self->num_mutations == other->num_mutations
        && self->mutations == other->mutations
        && memcmp(self->parent, other->parent, N * sizeof(uint32_t)) == 0
        && memcmp(self->population, other->population,
                N * sizeof(uint8_t)) == 0
        && memcmp(self->time, other->time, N * sizeof(double)) ==  0
        && memcmp(self->num_children, other->num_children,
                N * sizeof(uint32_t)) == 0
        && memcmp(self->children, other->children,
                N * sizeof(uint32_t *)) == 0;
    if (condition) {
        ret = 0;
    }
out:
    return ret;
}

static int
sparse_tree_check_node(sparse_tree_t *self, uint32_t u)
{
    int ret = 0;
    if (u >= self->num_nodes) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
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
        ret = MSP_ERR_OUT_OF_BOUNDS;
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

    if (self->flags & MSP_LEAF_COUNTS) {
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
    if (! (self->flags & MSP_LEAF_COUNTS)) {
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
    if (! (self->flags & MSP_LEAF_LISTS)) {
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

static void
sparse_tree_check_state(sparse_tree_t *self)
{
    uint32_t u, v, j, k, num_leaves;
    int err, found;

    for (j = 0; j < self->sample_size; j++) {
        u = j;
        assert(self->time[u] >= 0.0);
        assert(self->num_children[j] == 0);
        while (self->parent[u] != MSP_NULL_NODE) {
            v = self->parent[u];
            found = 0;
            for (k = 0; k < self->num_children[v]; k++) {
                if (self->children[v][k] == u) {
                    found = 1;
                }
            }
            assert(found);
            u = v;
            assert(self->time[u] > 0.0);
        }
        assert(u == self->root);
    }
    if (self->flags & MSP_LEAF_COUNTS) {
        assert(self->num_leaves != NULL);
        assert(self->num_tracked_leaves != NULL);
        for (j = 0; j < self->num_nodes; j++) {
            err = sparse_tree_get_num_leaves_by_traversal(self, j,
                    &num_leaves);
            assert(err == 0);
            assert(num_leaves == self->num_leaves[j]);
        }
    } else {
        assert(self->num_leaves == NULL);
        assert(self->num_tracked_leaves == NULL);
    }
    if (self->flags & MSP_LEAF_LISTS) {
        assert(self->leaf_list_tail != NULL);
        assert(self->leaf_list_head != NULL);
        assert(self->leaf_list_node_mem != NULL);
    } else {
        assert(self->leaf_list_tail == NULL);
        assert(self->leaf_list_head == NULL);
        assert(self->leaf_list_node_mem == NULL);
    }
}

void
sparse_tree_print_state(sparse_tree_t *self, FILE *out)
{
    size_t j, k;
    leaf_list_node_t *u;

    fprintf(out, "Sparse tree state:\n");
    fprintf(out, "flags = %d\n", self->flags);
    fprintf(out, "left = %f\n", self->left);
    fprintf(out, "left_breakpoint = %d\n", self->left_breakpoint);
    fprintf(out, "right = %f\n", self->right);
    fprintf(out, "right_breakpoint = %d\n", self->right_breakpoint);
    fprintf(out, "root = %d\n", self->root);
    fprintf(out, "index = %d\n", (int) self->index);
    for (j = 0; j < self->num_nodes; j++) {
        fprintf(out, "\t%d\t%d\t%f\t%d\t(", (int) j, self->parent[j],
            self->time[j], self->population[j]);
        for (k = 0; k < self->num_children[j]; k++) {
            fprintf(out, "%d", self->children[j][k]);
            if (k < self->num_children[j] - 1) {
                fprintf(out, ", ");
            }
        }
        fprintf(out, ")");
        if (self->flags & MSP_LEAF_COUNTS) {
            fprintf(out, "\t%d\t%d\t%d", self->num_leaves[j],
                    self->num_tracked_leaves[j], self->marked[j]);
        }
        if (self->flags & MSP_LEAF_LISTS) {
            fprintf(out, "\t[");
            u = self->leaf_list_head[j];
            if (u != NULL) {
                while (1) {
                    fprintf(out, "%d ", u->node);
                    if (u == self->leaf_list_tail[j]) {
                        break;
                    }
                    u = u->next;
                }
            } else {
                assert(self->leaf_list_tail[j] == NULL);
            }

            fprintf(out, "]");
        }
        fprintf(out, "\n");
    }
    fprintf(out, "mutations = \n");
    for (j = 0; j < self->num_mutations; j++) {
        fprintf(out, "\t%d @ %f\n", self->mutations[j].node,
                self->mutations[j].position);
    }
    sparse_tree_check_state(self);
}

/* Methods for positioning the tree along the sequence */

static inline void
sparse_tree_propagate_leaf_count_loss(sparse_tree_t *self,
        uint32_t u)
{
    const uint32_t all_leaves_diff = self->num_leaves[u];
    const uint32_t tracked_leaves_diff = self->num_tracked_leaves[u];
    const uint8_t mark = self->mark;
    uint32_t v = u;

    /* propogate this loss up as far as we can */
    while (v != MSP_NULL_NODE) {
        self->num_leaves[v] -= all_leaves_diff;
        self->num_tracked_leaves[v] -= tracked_leaves_diff;
        self->marked[v] = mark;
        v = self->parent[v];
    }
}

/* Returns the index of child u in the parent v */
static inline uint32_t
sparse_tree_get_child_index(sparse_tree_t *self, uint32_t v,
        uint32_t u)
{
    uint32_t j;

    for (j = 0; j < self->num_children[v] && self->children[v][j] != u; j++);
    assert(j != self->num_children[v]);
    return j;
}

static inline void
sparse_tree_propagate_leaf_count_gain(sparse_tree_t *self, uint32_t u)
{
    uint32_t j, k, v, *c;
    uint32_t all_leaves_diff = 0;
    uint32_t tracked_leaves_diff = 0;
    const uint8_t mark = self->mark;

    c = self->children[u];
    k = self->num_children[u];
    for (j = 0; j < k; j++) {
        all_leaves_diff += self->num_leaves[c[j]];
        tracked_leaves_diff += self->num_tracked_leaves[c[j]];
    }
    /* propogate this gain up as far as we can */
    v = u;
    while (v != MSP_NULL_NODE) {
        self->num_leaves[v] += all_leaves_diff;
        self->num_tracked_leaves[v] += tracked_leaves_diff;
        self->marked[v] = mark;
        v = self->parent[v];
    }
}

static inline void
sparse_tree_propagate_leaf_list_gain(sparse_tree_t *self,
        uint32_t u)
{
    uint32_t j, k, v, w, *c;

    c = self->children[u];
    k = self->num_children[u];
    for (j = 1; j < k; j++) {
        assert(self->leaf_list_tail[c[j]] != NULL);
        self->leaf_list_tail[c[j - 1]]->next = self->leaf_list_head[c[j]];
    }
    self->leaf_list_head[u] = self->leaf_list_head[c[0]];
    self->leaf_list_tail[u] = self->leaf_list_tail[c[k - 1]];

    v = u;
    w = self->parent[v];
    while (w != MSP_NULL_NODE) {
        j = sparse_tree_get_child_index(self, w, v);
        if (j != 0) {
            break;
        }
        self->leaf_list_head[w] = self->leaf_list_head[u];
        v = w;
        w = self->parent[w];
    }

    v = u;
    w = self->parent[v];
    while (w != MSP_NULL_NODE) {
        j = sparse_tree_get_child_index(self, w, v);
        if (j != self->num_children[w] - 1) {
            break;
        }
        self->leaf_list_tail[w] = self->leaf_list_tail[u];
        v = w;
        w = self->parent[w];
    }
}

static inline void
sparse_tree_post_propagate_leaf_list_gain(sparse_tree_t *self, uint32_t u)
{
    uint32_t v, w, j;

    v = u;
    w = self->parent[v];
    while (w != MSP_NULL_NODE) {
        j = sparse_tree_get_child_index(self, w, v);
        if (j > 0) {
            self->leaf_list_tail[self->children[w][j - 1]]->next = self->leaf_list_head[v];
        }
        if (j < self->num_children[w] - 1) {
            self->leaf_list_tail[v]->next = self->leaf_list_head[self->children[w][j + 1]];
        }
        v = w;
        w = self->parent[w];
    }
}

static int
sparse_tree_advance(sparse_tree_t *self, int direction,
        uint32_t *out_breakpoints, uint32_t *out_order, size_t *out_index,
        uint32_t *in_breakpoints, uint32_t *in_order, size_t *in_index,
        int first_tree)
{
    int ret = 0;
    int direction_change = direction * (direction != self->direction);
    ssize_t in = (ssize_t) *in_index + direction_change;
    ssize_t out = (ssize_t) *out_index + direction_change;
    ssize_t h;
    uint32_t j, k, u, oldest_child;
    uint32_t x = in_breakpoints[in_order[in]];
    double oldest_child_time;
    tree_sequence_t *s = self->tree_sequence;
    ssize_t R = (ssize_t) s->num_records;
    size_t in_count = 0;
    size_t out_count = 0;

    while (out_breakpoints[out_order[out]] == x) {
        k = out_order[out];
        u = s->trees.records.node[k];
        out_count += self->num_children[u] - 1;
        oldest_child_time = -1;
        oldest_child = 0;
        for (j = 0; j < self->num_children[u]; j++) {
            self->parent[self->children[u][j]] = MSP_NULL_NODE;
            if (self->time[self->children[u][j]] > oldest_child_time) {
                oldest_child = self->children[u][j];
                oldest_child_time = self->time[self->children[u][j]];
            }
        }
        self->num_children[u] = 0;
        self->children[u] = NULL;
        self->time[u] = 0;
        self->population[u] = MSP_NULL_POPULATION_ID;
        if (u == self->root) {
            self->root = oldest_child;
        }
        if (self->flags & MSP_LEAF_COUNTS) {
            sparse_tree_propagate_leaf_count_loss(self, u);
        }
        if (self->flags & MSP_LEAF_LISTS) {
            self->leaf_list_head[u] = NULL;
            self->leaf_list_tail[u] = NULL;
        }
        out += direction;
    }

    h = in;
    while (in >= 0 && in < R && in_breakpoints[in_order[in]] == x) {
        k = in_order[in];
        u = s->trees.records.node[k];
        for (j = 0; j < s->trees.records.num_children[k]; j++) {
            self->parent[s->trees.records.children[k][j]] = u;
        }
        self->num_children[u] = s->trees.records.num_children[k];
        in_count += self->num_children[u] - 1;
        self->children[u] = s->trees.records.children[k];
        self->time[u] = s->trees.nodes.time[u];
        self->population[u] = s->trees.nodes.population[u];
        if (self->time[u] > self->time[self->root]) {
            self->root = u;
        }
        if (self->flags & MSP_LEAF_COUNTS) {
            sparse_tree_propagate_leaf_count_gain(self, u);
        }
        if (self->flags & MSP_LEAF_LISTS) {
            sparse_tree_propagate_leaf_list_gain(self, u);
        }
        in += direction;
    }
    if (self->flags & MSP_LEAF_LISTS) {
        while (h >= 0 && h < R && in_breakpoints[in_order[h]] == x) {
            k = in_order[h];
            u = s->trees.records.node[k];
            sparse_tree_post_propagate_leaf_list_gain(self, u);
            h += direction;
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
    while (self->parent[self->root] != MSP_NULL_NODE) {
        self->root = self->parent[self->root];
    }

    self->direction = direction;
    self->index = (uint32_t) ((int) self->index + direction);
    *out_index = (size_t) out;
    *in_index = (size_t) in;
    if (s->num_mutations > 0) {
        self->mutations = s->mutations.tree_mutations[self->index];
        self->num_mutations = s->mutations.num_tree_mutations[self->index];
    }

    /* TODO These are all redundant and can be derived from the tree index */
    self->left_breakpoint = (uint32_t) self->index;
    self->right_breakpoint = (uint32_t) self->index + 1;
    self->left = s->trees.breakpoints[self->left_breakpoint];
    self->right = s->trees.breakpoints[self->right_breakpoint];
    ret = 1;
out:
    return ret;
}


int WARN_UNUSED
sparse_tree_first(sparse_tree_t *self)
{
    int ret = 0;
    tree_sequence_t *s = self->tree_sequence;

    /* TODO this is redundant if this is the first usage of the tree. We
     * should add a state machine here so we know what state the tree is
     * in and can take the appropriate actions.
     */
    ret = sparse_tree_clear(self);
    if (ret != 0) {
        goto out;
    }
    self->left_index = 0;
    self->right_index = 0;
    self->direction = MSP_DIR_FORWARD;

    ret = sparse_tree_advance(self, MSP_DIR_FORWARD,
            s->trees.records.right, s->trees.indexes.removal_order,
            &self->right_index, s->trees.records.left,
            s->trees.indexes.insertion_order, &self->left_index, 1);
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_last(sparse_tree_t *self)
{
    int ret = 0;
    tree_sequence_t *s = self->tree_sequence;

    /* TODO this is redundant if this is the first usage of the tree. We
     * should add a state machine here so we know what state the tree is
     * in and can take the appropriate actions.
     */
    ret = sparse_tree_clear(self);
    if (ret != 0) {
        goto out;
    }
    self->left_index = s->num_records - 1;
    self->right_index = s->num_records - 1;
    self->direction = MSP_DIR_REVERSE;
    self->index = tree_sequence_get_num_trees(s);

    ret = sparse_tree_advance(self, MSP_DIR_REVERSE,
            s->trees.records.left, s->trees.indexes.insertion_order,
            &self->left_index, s->trees.records.right,
            s->trees.indexes.removal_order, &self->right_index, 1);
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_next(sparse_tree_t *self)
{
    int ret = 0;
    tree_sequence_t *s = self->tree_sequence;
    size_t num_trees = tree_sequence_get_num_trees(s);

    if (self->index < num_trees - 1) {
        ret = sparse_tree_advance(self, MSP_DIR_FORWARD,
                s->trees.records.right, s->trees.indexes.removal_order,
                &self->right_index, s->trees.records.left,
                s->trees.indexes.insertion_order, &self->left_index, 0);
    }
    return ret;
}

int WARN_UNUSED
sparse_tree_prev(sparse_tree_t *self)
{
    int ret = 0;
    tree_sequence_t *s = self->tree_sequence;

    if (self->index > 0) {
        ret = sparse_tree_advance(self, MSP_DIR_REVERSE,
                s->trees.records.left, s->trees.indexes.insertion_order,
                &self->left_index, s->trees.records.right,
                s->trees.indexes.removal_order, &self->right_index, 0);
    }
    return ret;
}
