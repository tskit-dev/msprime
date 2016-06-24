/*
** Copyright (C) 2015 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
tree_sequence_print_state(tree_sequence_t *self)
{
    size_t j;

    printf("tree_sequence state\n");
    printf("sample_size = %d\n", self->sample_size);
    printf("sequence_length = %f\n", self->sequence_length);
    printf("samples\n");
    for (j = 0; j < self->sample_size; j++) {
        printf("\t%d\t%d\n", (int) j, (int) self->samples.population[j]);
    }
    printf("trees = (%d records)\n", (int) self->num_records);
    printf("\tparameters = '%s'\n", self->trees.parameters);
    printf("\tenvironment = '%s'\n", self->trees.environment);
    for (j = 0; j < self->num_records; j++) {
        printf("\t%d\t%f\t%f\t%d\t%d\t%d\t%f\t%d\t|\t%d\t%d\n",
                (int) j,
                self->trees.left[j],
                self->trees.right[j],
                (int) self->trees.node[j],
                (int) self->trees.children[2 * j],
                (int) self->trees.children[2 * j + 1],
                self->trees.time[j],
                (int) self->trees.population[j],
                (int) self->trees.insertion_order[j],
                (int) self->trees.removal_order[j]);
    }
    printf("mutations = (%d records)\n", (int) self->num_mutations);
    printf("\tparameters = '%s'\n", self->mutations.parameters);
    printf("\tenvironment = '%s'\n", self->mutations.environment);
    for (j = 0; j < self->num_mutations; j++) {
        printf("\t%d\t%f\n", (int) self->mutations.node[j],
                self->mutations.position[j]);
    }

}

/* Allocates the memory required for arrays of values. Assumes that
 * the num_records and num_mutations have been set.
 */
static int
tree_sequence_alloc(tree_sequence_t *self)
{
    int ret = MSP_ERR_NO_MEMORY;

    self->samples.population = malloc(self->sample_size * sizeof(uint8_t));
    if (self->samples.population == NULL) {
        goto out;
    }
    self->trees.left = malloc(self->num_records * sizeof(double));
    self->trees.right = malloc(self->num_records * sizeof(double));
    self->trees.children = malloc(2 * self->num_records * sizeof(uint32_t));
    self->trees.node = malloc(self->num_records * sizeof(uint32_t));
    self->trees.population = malloc(self->num_records * sizeof(uint8_t));
    self->trees.time = malloc(self->num_records * sizeof(double));
    self->trees.insertion_order = malloc(self->num_records * sizeof(uint32_t));
    self->trees.removal_order = malloc(self->num_records * sizeof(uint32_t));
    if (self->trees.left == NULL || self->trees.right == NULL
            || self->trees.children == NULL || self->trees.node == NULL
            || self->trees.time == NULL || self->trees.population == NULL
            || self->trees.insertion_order == NULL
            || self->trees.removal_order == NULL) {
        goto out;
    }
    /* Set the optional fields to their unset values. */
    memset(self->samples.population, MSP_NULL_POPULATION_ID, self->sample_size);
    memset(self->trees.population, MSP_NULL_POPULATION_ID, self->num_records);
    if (self->num_mutations > 0) {
        self->mutations.node = malloc(self->num_mutations * sizeof(uint32_t));
        self->mutations.position = malloc(
                self->num_mutations * sizeof(double));
        if (self->mutations.node == NULL || self->mutations.position == NULL) {
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

int
tree_sequence_free(tree_sequence_t *self)
{
    if (self->samples.population != NULL) {
        free(self->samples.population);
    }
    if (self->trees.left != NULL) {
        free(self->trees.left);
    }
    if (self->trees.right != NULL) {
        free(self->trees.right);
    }
    if (self->trees.children != NULL) {
        free(self->trees.children);
    }
    if (self->trees.node != NULL) {
        free(self->trees.node);
    }
    if (self->trees.population != NULL) {
        free(self->trees.population);
    }
    if (self->trees.time != NULL) {
        free(self->trees.time);
    }
    if (self->trees.insertion_order != NULL) {
        free(self->trees.insertion_order);
    }
    if (self->trees.removal_order != NULL) {
        free(self->trees.removal_order);
    }
    if (self->trees.parameters != NULL) {
        free(self->trees.parameters);
    }
    if (self->trees.environment != NULL) {
        free(self->trees.environment);
    }
    if (self->mutations.node != NULL) {
        free(self->mutations.node);
    }
    if (self->mutations.position != NULL) {
        free(self->mutations.position);
    }
    if (self->mutations.parameters != NULL) {
        free(self->mutations.parameters);
    }
    if (self->mutations.environment != NULL) {
        free(self->mutations.environment);
    }
    return 0;
}

static int
tree_sequence_make_indexes(tree_sequence_t *self)
{
    int ret = 0;
    uint32_t j;
    index_sort_t *sort_buff = NULL;

    sort_buff = malloc(self->num_records * sizeof(index_sort_t));
    if (sort_buff == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    /* sort by left and increasing time to give us the order in which
     * records should be inserted */
    for (j = 0; j < self->num_records; j++) {
        sort_buff[j].index = j;
        sort_buff[j].value = self->trees.left[j];
        sort_buff[j].time = self->trees.time[j];
    }
    qsort(sort_buff, self->num_records, sizeof(index_sort_t), cmp_index_sort);
    for (j = 0; j < self->num_records; j++) {
        self->trees.insertion_order[j] = sort_buff[j].index;
    }
    /* sort by right and decreasing time to give us the order in which
     * records should be removed. */
    for (j = 0; j < self->num_records; j++) {
        sort_buff[j].index = j;
        sort_buff[j].value = self->trees.right[j];
        sort_buff[j].time = -self->trees.time[j];
    }
    qsort(sort_buff, self->num_records, sizeof(index_sort_t), cmp_index_sort);
    for (j = 0; j < self->num_records; j++) {
        self->trees.removal_order[j] = sort_buff[j].index;
    }
    /* set the num_nodes value */
    self->num_nodes = self->trees.node[self->num_records - 1] + 1;
out:
    if (sort_buff != NULL) {
        free(sort_buff);
    }
    return ret;
}

static int
tree_sequence_remap_coordinates(tree_sequence_t *self, recomb_map_t *recomb_map)
{
    int ret = 0;
    size_t j;
    double *phys_left = NULL;
    double *phys_right = NULL;

    phys_left = malloc(self->num_records * sizeof(double));
    phys_right = malloc(self->num_records * sizeof(double));
    if (phys_left == NULL || phys_right == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    /* we know that these orders ensure that the left and right coordinates
     * are sorted, and we use this to make bulk remapping to physical coordinates
     * efficient.
     */
    for (j = 0; j < self->num_records; j++) {
        phys_left[j] = self->trees.left[self->trees.insertion_order[j]];
        phys_right[j] = self->trees.right[self->trees.removal_order[j]];
    }
    ret = recomb_map_genetic_to_phys_bulk(
        recomb_map, phys_left, self->num_records);
    if (ret != 0) {
        goto out;
    }
    ret = recomb_map_genetic_to_phys_bulk(
        recomb_map, phys_right, self->num_records);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < self->num_records; j++) {
        self->trees.left[self->trees.insertion_order[j]] = phys_left[j];
        self->trees.right[self->trees.removal_order[j]] = phys_right[j];
    }
out:
    if (phys_left != NULL) {
        free(phys_left);
    }
    if (phys_right != NULL) {
        free(phys_right);
    }
    return ret;
}

int
tree_sequence_create(tree_sequence_t *self, msp_t *sim,
        recomb_map_t *recomb_map, double Ne)
{
    int ret = -1;
    uint32_t j;
    coalescence_record_t *records = NULL;
    sample_t *samples = NULL;
    char *parameters;

    memset(self, 0, sizeof(tree_sequence_t));
    self->num_records = msp_get_num_coalescence_records(sim);
    assert(self->num_records > 0);
    self->sample_size = sim->sample_size;
    self->sequence_length = recomb_map_get_sequence_length(recomb_map);
    self->num_mutations = 0;
    ret = tree_sequence_alloc(self);
    if (ret != 0) {
        goto out;
    }
    records = malloc(self->num_records * sizeof(coalescence_record_t));
    samples = malloc(self->sample_size * sizeof(sample_t));
    if (records == NULL || samples == NULL) {
        goto out;
    }
    ret = msp_get_coalescence_records(sim, records);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < self->num_records; j++) {
        self->trees.left[j] = records[j].left;
        self->trees.right[j] = records[j].right;
        assert(self->trees.left[j] <= sim->num_loci);
        assert(self->trees.right[j] <= sim->num_loci);
        self->trees.node[j] = records[j].node;
        self->trees.population[j] = records[j].population_id;
        self->trees.children[2 * j] = records[j].children[0];
        self->trees.children[2 * j + 1] = records[j].children[1];
        /* Rescale time into generations. */
        self->trees.time[j] = records[j].time * 4 * Ne;
    }
    ret = msp_get_samples(sim, samples);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < self->sample_size; j++) {
        self->samples.population[j] = samples[j].population_id;
    }
    ret = tree_sequence_make_indexes(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_remap_coordinates(self, recomb_map);
    if (ret != 0) {
        goto out;
    }
    parameters = msp_get_configuration_json(sim);
    assert(parameters != NULL);
    self->trees.parameters = malloc(strlen(parameters) + 1);
    if (self->trees.parameters == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    strcpy(self->trees.parameters, parameters);
    ret = msp_encode_environment(&self->trees.environment);
out:
    if (records != NULL) {
        free(records);
    }
    if (samples != NULL) {
        free(samples);
    }
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
        int dimensions;
        size_t size;
        int required;
    };
    struct _dimension_check fields[] = {
        {"/trees/left", 1, 0, 1},
        {"/trees/right", 1, 0, 1},
        {"/trees/node", 1, 0, 1},
        {"/trees/population", 1, 0, 0},
        {"/trees/children", 2, 0, 1},
        {"/trees/time", 1, 0, 1},
        {"/mutations/node", 1, 0, 1},
        {"/mutations/position", 1, 0, 1},
        {"/samples/population", 1, 0, 0},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _dimension_check);
    size_t j;

    for (j = 0; j < 6; j++) {
        fields[j].size = self->num_records;
    }
    for (j = 6; j < 8; j++) {
        fields[j].size = self->num_mutations;
        fields[j].required = self->num_mutations > 0;
    }
    fields[8].size = self->sample_size;
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
            if (rank != fields[j].dimensions) {
                ret = MSP_ERR_FILE_FORMAT;
                goto out;
            }
            status = H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
            if (status < 0) {
                goto out;
            }
            if (dims[0] != fields[j].size) {
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
        {"/trees/left", NULL, 1},
        {"/mutations/node", NULL, 0},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _dimension_read);
    size_t j;

    fields[0].dest = &self->num_records;
    fields[1].dest = &self->num_mutations;
    /* check if the mutations group exists */
    exists = H5Lexists(file_id, "/mutations", H5P_DEFAULT);
    if (exists < 0) {
        goto out;
    }
    self->num_mutations = 0;
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
        {"/trees/left", H5T_NATIVE_DOUBLE, 0, 1, NULL},
        {"/trees/right", H5T_NATIVE_DOUBLE, 0, 1, NULL},
        {"/trees/node", H5T_NATIVE_UINT32, 0, 1, NULL},
        {"/trees/population", H5T_NATIVE_UINT8, 0, 0, NULL},
        {"/trees/children", H5T_NATIVE_UINT32, 0, 1, NULL},
        {"/trees/time", H5T_NATIVE_DOUBLE, 0, 1, NULL},
        {"/mutations/node", H5T_NATIVE_UINT32, 0, 1, NULL},
        {"/mutations/position", H5T_NATIVE_DOUBLE, 0, 1, NULL},
        {"/samples/population", H5T_NATIVE_UINT8, 0, 0, NULL},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _hdf5_field_read);
    size_t j;

    fields[0].dest = self->trees.left;
    fields[1].dest = self->trees.right;
    fields[2].dest = self->trees.node;
    fields[3].dest = self->trees.population;
    fields[4].dest = self->trees.children;
    fields[5].dest = self->trees.time;
    fields[6].dest = self->mutations.node;
    fields[7].dest = self->mutations.position;
    fields[8].dest = self->samples.population;
    /* TODO We're sort of doing the same thing twice here as
     * the mutations _group_ is optional. However, we can't just
     * mark mutations/node and mutations/position as optional as we
     * would then allow one or the other. This would be an error.
     * However, we should improve this logic as it's a bit messy.
     */
    if (self->num_mutations == 0) {
        fields[6].empty = 1;
        fields[7].empty = 1;
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
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_read_hdf5_provenance(tree_sequence_t *self, hid_t file_id)
{
    int ret = MSP_ERR_HDF5;
    hid_t attr_id, atype, type_class, atype_mem;
    herr_t status;
    size_t size;
    struct _hdf5_string_read {
        const char *prefix;
        const char *name;
        char **dest;
        int included;
    };
    struct _hdf5_string_read fields[] = {
        {"trees", "environment", NULL, 1},
        {"trees", "parameters", NULL, 1},
        {"mutations", "environment", NULL, 0},
        {"mutations", "parameters", NULL, 0},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _hdf5_string_read);
    size_t j;
    char *str;

    fields[0].dest = &self->trees.environment;
    fields[1].dest = &self->trees.parameters;
    if (self->num_mutations > 0) {
        fields[2].included = 1;
        fields[3].included = 1;
        fields[2].dest = &self->mutations.environment;
        fields[3].dest = &self->mutations.parameters;
    }

    for (j = 0; j < num_fields; j++) {
        if (fields[j].included) {
            attr_id = H5Aopen_by_name(file_id, fields[j].prefix, fields[j].name,
                    H5P_DEFAULT, H5P_DEFAULT);
            if (attr_id < 0) {
                goto out;
            }
            atype = H5Aget_type(attr_id);
            if (atype < 0) {
                goto out;
            }
            type_class = H5Tget_class(atype);
            if (type_class < 0) {
                goto out;
            }
            if (type_class != H5T_STRING) {
                ret = MSP_ERR_FILE_FORMAT;
                goto out;
            }
            atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
            if (atype_mem < 0) {
                goto out;
            }
            size = H5Tget_size(atype_mem);
            str = malloc(size + 1);
            if (str == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            status = H5Aread(attr_id, atype_mem, str);
            if (status < 0) {
                goto out;
            }
            str[size] = '\0';
            *fields[j].dest = str;
            status = H5Tclose(atype);
            if (status < 0) {
                goto out;
            }
            status = H5Tclose(atype_mem);
            if (status < 0) {
                goto out;
            }
            status = H5Aclose(attr_id);
            if (status < 0) {
                goto out;
            }
        }
    }
    ret = 0;
out:
    return ret;
}

int
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
    ret = tree_sequence_read_hdf5_provenance(self, file_id);
    if (ret < 0) {
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
    ret = tree_sequence_make_indexes(self);
    if (ret != 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_write_hdf5_data(tree_sequence_t *self, hid_t file_id, int flags)
{
    herr_t ret = -1;
    herr_t status;
    hid_t group_id, dataset_id, dataspace_id, plist_id;
    hsize_t dims[2];
    struct _hdf5_field_write {
        const char *name;
        hid_t storage_type;
        hid_t memory_type;
        int dimensions;
        size_t size;
        void *source;
    };
    struct _hdf5_field_write fields[] = {
        {"/trees/left", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, 1, 0, NULL},
        {"/trees/right", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, 1, 0, NULL},
        {"/trees/node", H5T_STD_U32LE, H5T_NATIVE_UINT32, 1, 0, NULL},
        {"/trees/population", H5T_STD_U8LE, H5T_NATIVE_UINT8, 1, 0, NULL},
        {"/trees/children", H5T_STD_U32LE, H5T_NATIVE_UINT32, 2, 0, NULL},
        {"/trees/time", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, 1, 0, NULL},
        {"/mutations/node", H5T_STD_U32LE, H5T_NATIVE_UINT32, 1, 0, NULL},
        {"/mutations/position", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, 1, 0, NULL},
        {"/samples/population", H5T_STD_U8LE, H5T_NATIVE_UINT8, 1, 0, NULL},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _hdf5_field_write);
    struct _hdf5_group_write {
        const char *name;
        int included;
    };
    struct _hdf5_group_write groups[] = {
        {"/trees", 1},
        {"/mutations", 1},
        {"/samples", 1},
    };
    size_t num_groups = sizeof(groups) / sizeof(struct _hdf5_group_write);
    size_t j;

    fields[0].source = self->trees.left;
    fields[1].source = self->trees.right;
    fields[2].source = self->trees.node;
    fields[3].source = self->trees.population;
    fields[4].source = self->trees.children;
    fields[5].source = self->trees.time;
    fields[6].source = self->mutations.node;
    fields[7].source = self->mutations.position;
    fields[8].source = self->samples.population;
    for (j = 0; j < 6; j++) {
        fields[j].size = self->num_records;
    }
    for (j = 6; j < 8; j++) {
        fields[j].size = self->num_mutations;
    }
    fields[8].size = self->sample_size;
    if (self->num_mutations == 0) {
        groups[1].included = 0;
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
            dims[1] = 2; /* unused except for children */
            dataspace_id = H5Screate_simple(fields[j].dimensions, dims, NULL);
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
            status = H5Pset_chunk(plist_id, fields[j].dimensions, dims);
            if (status < 0) {
                goto out;
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
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_write_hdf5_provenance(tree_sequence_t *self, hid_t file_id)
{
    herr_t ret = -1;
    herr_t status;
    hid_t group_id, dataspace_id, attr_id, type_id;
    struct _hdf5_string_write {
        const char *group;
        const char *name;
        const char *value;
        int included;
    };
    struct _hdf5_string_write fields[] = {
        {"trees", "environment", NULL, 1},
        {"trees", "parameters", NULL, 1},
        {"mutations", "environment", NULL, 0},
        {"mutations", "parameters", NULL, 0},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _hdf5_string_write);
    size_t j;

    fields[0].value = self->trees.environment;
    fields[1].value = self->trees.parameters;
    if (self->num_mutations > 0) {
        fields[2].included = 1;
        fields[2].value = self->mutations.environment;
        fields[3].included = 1;
        fields[3].value = self->mutations.parameters;
    }
    for (j = 0; j < num_fields; j++) {
        if (fields[j].included) {
            assert(fields[j].value != NULL);
            group_id = H5Gopen(file_id, fields[j].group, H5P_DEFAULT);
            if (group_id < 0) {
                goto out;
            }
            dataspace_id = H5Screate(H5S_SCALAR);
            if (dataspace_id < 0) {
                goto out;
            }
            type_id = H5Tcopy(H5T_C_S1);
            if (type_id < 0) {
                goto out;
            }
            status = H5Tset_size(type_id, strlen(fields[j].value));
            if (status < 0) {
                goto out;
            }
            attr_id = H5Acreate(group_id, fields[j].name, type_id,
                    dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
            if (attr_id < 0) {
                goto out;
            }
            status = H5Awrite(attr_id, type_id, fields[j].value);
            if (status < 0) {
                goto out;
            }
            status = H5Aclose(attr_id);
            if (status < 0) {
                goto out;
            }
            status = H5Tclose(type_id);
            if (status < 0) {
                goto out;
            }
            status = H5Sclose(dataspace_id);
            if (status < 0) {
                goto out;
            }
            status = H5Gclose(group_id);
            if (status < 0) {
                goto out;
            }
        }
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

int
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
    status = tree_sequence_write_hdf5_provenance(self, file_id);
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
    return self->num_nodes;
}

int
tree_sequence_get_population(tree_sequence_t *self, uint32_t u,
        uint32_t *population_id)
{
    int ret = 0;

    if (u >= self->sample_size) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    *population_id = self->samples.population[u];
out:
    return ret;
}

int
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

/* Returns the parameters for the trees encoded as JSON. This string
 * should NOT be freed by client code.
 */
char *
tree_sequence_get_simulation_parameters(tree_sequence_t *self)
{
    return self->trees.parameters;
}

/* Returns the parameters for the mutations encoded as JSON. This string
 * should NOT be freed by client code. This is NULL if mutations have
 * not been generated.
 */
char *
tree_sequence_get_mutation_parameters(tree_sequence_t *self)
{
    return self->mutations.parameters;
}

int
tree_sequence_get_record(tree_sequence_t *self, size_t index,
        coalescence_record_t *record, int order)
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
            j = self->trees.insertion_order[index];
            break;
        case MSP_ORDER_RIGHT:
            j = self->trees.removal_order[index];
            break;
        default:
            ret = MSP_ERR_BAD_ORDERING;
            goto out;
    }
    record->left = self->trees.left[j];
    record->right = self->trees.right[j];
    record->node = self->trees.node[j];
    record->children[0] = self->trees.children[2 * j];
    record->children[1] = self->trees.children[2 * j + 1];
    record->time = self->trees.time[j];
    record->population_id = self->trees.population[j];
out:
    return ret;
}

int
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
int
tree_sequence_alloc_sparse_tree(tree_sequence_t *self, sparse_tree_t *tree,
        uint32_t *tracked_leaves, uint32_t num_tracked_leaves, int flags)
{
    return sparse_tree_alloc(tree, self->sample_size, self->num_nodes,
            self->num_mutations, tracked_leaves, num_tracked_leaves, flags);
}

int
tree_sequence_set_mutations(tree_sequence_t *self, size_t num_mutations,
        mutation_t *mutations, const char *parameters, const char* environment)
{
    int ret = -1;
    size_t j, len;
    mutation_t **mutation_ptrs = NULL;

    if (self->num_mutations > 0) {
        /* any mutations that were there previously are overwritten. */
        if (self->mutations.node != NULL) {
            free(self->mutations.node);
        }
        if (self->mutations.position != NULL) {
            free(self->mutations.position);
        }
        if (self->mutations.parameters != NULL) {
            free(self->mutations.parameters);
        }
        if (self->mutations.environment != NULL) {
            free(self->mutations.environment);
        }
    }
    self->num_mutations = 0;
    self->mutations.position = NULL;
    self->mutations.node = NULL;
    self->mutations.parameters = NULL;
    self->mutations.environment = NULL;
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
        /* Make copies of the environment and parameters strings. */
        assert(parameters != NULL);
        len = strlen(parameters) + 1;
        self->mutations.parameters = malloc(len);
        if (self->mutations.parameters == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        strncpy(self->mutations.parameters, parameters, len);
        assert(environment != NULL);
        len = strlen(environment) + 1;
        self->mutations.environment = malloc(len);
        if (self->mutations.environment == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        strncpy(self->mutations.environment, environment, len);
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

int
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

int
tree_diff_iterator_free(tree_diff_iterator_t *self)
{
    int ret = 0;
    if (self->node_records != NULL) {
        free(self->node_records);
    }
    return ret;
}

void
tree_diff_iterator_print_state(tree_diff_iterator_t *self)
{
    printf("tree_diff_iterator state\n");
    printf("num_records = %d\n", (int) self->num_records);
    printf("insertion_index = %d\n", (int) self->insertion_index);
    printf("removal_index = %d\n", (int) self->removal_index);
    printf("tree_left = %f\n", self->tree_left);
}

int
tree_diff_iterator_next(tree_diff_iterator_t *self, double *length,
        node_record_t **nodes_out, node_record_t **nodes_in)
{
    int ret = 0;
    uint32_t k;
    double last_left = self->tree_left;
    size_t next_node_record = 0;
    tree_sequence_t *s = self->tree_sequence;
    node_record_t *out_head = NULL;
    node_record_t *out_tail = NULL;
    node_record_t *in_head = NULL;
    node_record_t *in_tail = NULL;
    node_record_t *w = NULL;

    assert(s != NULL);
    if (self->insertion_index < self->num_records) {
        /* First we remove the stale records */
        while (s->trees.right[s->trees.removal_order[self->removal_index]]
                == self->tree_left) {
            k = s->trees.removal_order[self->removal_index];
            assert(next_node_record < 2 * self->sample_size);
            w = &self->node_records[next_node_record];
            next_node_record++;
            w->time = s->trees.time[k];
            w->node = s->trees.node[k];
            w->children[0] = s->trees.children[2 * k];
            w->children[1] = s->trees.children[2 * k + 1];
            w->next = NULL;
            if (out_head == NULL) {
                out_head = w;
                out_tail = w;
            } else {
                out_tail->next = w;
                out_tail = w;
            }
            self->removal_index++;
        }
        /* Now insert the new records */
        while (self->insertion_index < self->num_records &&
                s->trees.left[s->trees.insertion_order[self->insertion_index]]
                == self->tree_left) {
            k = s->trees.insertion_order[self->insertion_index];
            assert(next_node_record < 2 * self->sample_size);
            w = &self->node_records[next_node_record];
            next_node_record++;
            w->time = s->trees.time[k];
            w->node = s->trees.node[k];
            w->children[0] = s->trees.children[2 * k];
            w->children[1] = s->trees.children[2 * k + 1];
            w->next = NULL;
            if (in_head == NULL) {
                in_head = w;
                in_tail = w;
            } else {
                in_tail->next = w;
                in_tail = w;
            }
            self->insertion_index++;
        }
        /* Update the left coordinate */
        self->tree_left = s->trees.right[s->trees.removal_order[
            self->removal_index]];
        ret = 1;
    }
    *nodes_out = out_head;
    *nodes_in = in_head;
    *length = self->tree_left - last_left;
    return ret;
}

/* ======================================================== *
 * sparse tree
 * ======================================================== */

int
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
    self->children = malloc(2 * self->num_nodes * sizeof(uint32_t));
    if (self->time == NULL || self->parent == NULL || self->children == NULL
            || self->population == NULL) {
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

int
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

int
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
    memset(self->children, (int) MSP_NULL_NODE, 2 * N * sizeof(uint32_t));
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

int
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
        } else if (self->children[2 * v] != MSP_NULL_NODE) {
            for (c = 0; c < 2; c++) {
                stack_top++;
                stack[stack_top] = self->children[2 * v + c];
            }
        }
    }
    *num_leaves = count;
    return ret;
}

int
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

int
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

int
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

int
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

/* ======================================================== *
 * sparse tree iterator
 * ======================================================== */

int
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
        self->tree->population[j] = self->tree_sequence->samples.population[j];
    }
out:
    return ret;
}

int
sparse_tree_iterator_free(sparse_tree_iterator_t *self)
{
    int ret = 0;
    return ret;
}

static void
sparse_tree_iterator_check_state(sparse_tree_iterator_t *self)
{
    uint32_t u, v, j, num_leaves;
    int err;

    assert(self->tree->num_nodes == self->num_nodes);
    for (j = 0; j < self->sample_size; j++) {
        u = j;
        assert(self->tree->time[u] == 0.0);
        assert(self->tree->children[2 * j] == MSP_NULL_NODE);
        assert(self->tree->children[2 * j + 1] == MSP_NULL_NODE);
        while (self->tree->parent[u] != MSP_NULL_NODE) {
            assert(self->tree->population[u] != MSP_NULL_POPULATION_ID);
            v = self->tree->parent[u];
            assert(self->tree->children[2 * v] == u
                    || self->tree->children[2 * v + 1] == u);
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
sparse_tree_iterator_print_state(sparse_tree_iterator_t *self)
{
    size_t j;
    uint32_t u;

    printf("sparse_tree_iterator state\n");
    printf("insertion_index = %d\n", (int) self->insertion_index);
    printf("removal_index = %d\n", (int) self->removal_index);
    printf("mutation_index = %d\n", (int) self->mutation_index);
    printf("num_records = %d\n", (int) self->num_records);
    printf("tree.flags = %d\n", self->tree->flags);
    printf("tree.left = %f\n", self->tree->left);
    printf("tree.right = %f\n", self->tree->right);
    printf("tree.root = %d\n", self->tree->root);
    printf("tree.index = %d\n", self->tree->index);
    for (j = 0; j < self->tree->num_nodes; j++) {
        printf("\t%d\t%d\t%d\t%d\t%f\t%d", (int) j, self->tree->parent[j],
                self->tree->children[2 * j], self->tree->children[2 * j + 1],
                self->tree->time[j], self->tree->population[j]);
        if (self->tree->flags & MSP_COUNT_LEAVES) {
            printf("\t%d\t%d", self->tree->num_leaves[j],
                    self->tree->num_tracked_leaves[j]);
            u = 0;
            if (self->tree->leaf_list_head[j] != NULL) {
                u = self->tree->leaf_list_head[j]->node;
            }
            printf("\t%d", u);
            u = 0;
            if (self->tree->leaf_list_tail[j] != NULL) {
                u = self->tree->leaf_list_tail[j]->node;
            }
            printf("\t%d", u);
        }
        printf("\n");
    }
    printf("mutations = \n");
    for (j = 0; j < self->tree->num_mutations; j++) {
        printf("\t%d @ %f\n", self->tree->mutations[j].node,
                self->tree->mutations[j].position);
    }
    sparse_tree_iterator_check_state(self);
}

int
sparse_tree_iterator_next(sparse_tree_iterator_t *self)
{
    int ret = 0;
    uint32_t j, k, u, v, c[2], all_leaves_diff, tracked_leaves_diff;
    tree_sequence_t *s = self->tree_sequence;
    sparse_tree_t *t = self->tree;

    assert(t != NULL && s != NULL);
    if (self->insertion_index < self->num_records) {
        /* First we remove the stale records */
        while (s->trees.right[s->trees.removal_order[self->removal_index]]
                == t->right) {
            k = s->trees.removal_order[self->removal_index];
            u = s->trees.node[k];
            c[0] = s->trees.children[2 * k];
            c[1] = s->trees.children[2 * k + 1];
            for (j = 0; j < 2; j++) {
                t->parent[c[j]] = MSP_NULL_NODE;
                t->children[2 * u + j] = MSP_NULL_NODE;
            }
            t->time[u] = 0;
            t->population[u] = MSP_NULL_POPULATION_ID;
            if (u == t->root) {
                t->root = GSL_MAX(c[0], c[1]);
            }
            self->removal_index++;
            if (t->flags & MSP_COUNT_LEAVES) {
                all_leaves_diff = t->num_leaves[u];
                tracked_leaves_diff = t->num_tracked_leaves[u];
                /* propogate this loss up as far as we can */
                v = u;
                while (v != MSP_NULL_NODE) {
                    t->num_leaves[v] -= all_leaves_diff;
                    t->num_tracked_leaves[v] -= tracked_leaves_diff;
                    t->leaf_list_head[v] = NULL;
                    t->leaf_list_tail[v] = NULL;
                    v = t->parent[v];
                }
            }
        }
        /* Update the interval */
        t->left = t->right;
        t->right = s->trees.right[s->trees.removal_order[self->removal_index]];
        /* Now insert the new records */
        while (self->insertion_index < self->num_records &&
                s->trees.left[s->trees.insertion_order[self->insertion_index]]
                == t->left) {
            k = s->trees.insertion_order[self->insertion_index];
            u = s->trees.node[k];
            c[0] = s->trees.children[2 * k];
            c[1] = s->trees.children[2 * k + 1];
            for (j = 0; j < 2; j++) {
                t->parent[c[j]] = u;
                t->children[2 * u + j] = c[j];
            }
            t->time[u] = s->trees.time[k];
            t->population[u] = s->trees.population[k];
            if (u >t->root) {
                t->root = u;
            }
            self->insertion_index++;
            if (t->flags & MSP_COUNT_LEAVES) {
                all_leaves_diff = t->num_leaves[c[0]] + t->num_leaves[c[1]];
                tracked_leaves_diff = t->num_tracked_leaves[c[0]]
                    + t->num_tracked_leaves[c[1]];
                /* propogate this gain up as far as we can */
                v = u;
                while (v != MSP_NULL_NODE) {
                    t->num_leaves[v] += all_leaves_diff;
                    t->num_tracked_leaves[v] += tracked_leaves_diff;
                    c[0] = t->children[2 * v];
                    c[1] = t->children[2 * v + 1];
                    if (t->leaf_list_head[c[0]] == NULL) {
                        t->leaf_list_head[v] = t->leaf_list_head[c[1]];
                        t->leaf_list_tail[v] = t->leaf_list_tail[c[1]];
                    } else if (t->leaf_list_head[c[1]] == NULL) {
                        t->leaf_list_head[v] = t->leaf_list_head[c[0]];
                        t->leaf_list_tail[v] = t->leaf_list_tail[c[0]];
                    } else {
                        t->leaf_list_head[v] = t->leaf_list_head[c[0]];
                        t->leaf_list_tail[v] = t->leaf_list_tail[c[1]];
                        assert(t->leaf_list_tail[c[0]] != NULL);
                        t->leaf_list_tail[c[0]]->next =
                            t->leaf_list_head[c[1]];
                    }
                    v = t->parent[v];
                }
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
            t->mutations[t->num_mutations].position =
                s->mutations.position[self->mutation_index];
            t->mutations[t->num_mutations].node =
                s->mutations.node[self->mutation_index];
            self->mutation_index++;
            t->num_mutations++;
        }
        /* Finally, update the tree index and indicate we have a valid tree */
        t->index++;
        ret = 1;
    }
    return ret;
}
