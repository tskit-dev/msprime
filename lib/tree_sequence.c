/*
** Copyright (C) 2014 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

#include "err.h"
#include "msprime.h"

/*
 * Comparator for coalescence records. Sort by left.
 */
static int
cmp_coalescence_record(const void *a, const void *b) {
    const coalescence_record_t *ca = (const coalescence_record_t *) a;
    const coalescence_record_t *cb = (const coalescence_record_t *) b;
    return (ca->left > cb->left) - (ca->left < cb->left);
}

int
tree_sequence_create(tree_sequence_t *self, msp_t *sim)
{
    int ret = -1;
    uint32_t j;
    coalescence_record_t *records = NULL;

    printf("Creating tree_sequence\n");
    self->num_records = msp_get_num_coalescence_records(sim);
    self->left = malloc(self->num_records * sizeof(uint32_t));
    self->right = malloc(self->num_records * sizeof(uint32_t));
    self->children = malloc(2 * self->num_records * sizeof(uint32_t));
    self->parent = malloc(self->num_records * sizeof(uint32_t));
    self->time = malloc(self->num_records * sizeof(double));
    if (self->left == NULL || self->right == NULL || self->children == NULL
            || self->parent == NULL || self->time == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    records = malloc(self->num_records * sizeof(coalescence_record_t));
    if (records == NULL) {
        goto out;
    }
    ret = msp_get_coalescence_records(sim, records);
    if (ret != 0) {
        goto out;
    }
    /* Sort the records by the left coordinate*/
    qsort(records, self->num_records, sizeof(coalescence_record_t),
            cmp_coalescence_record);
    for (j = 0; j < self->num_records; j++) {
        self->left[j] = records[j].left;
        self->right[j] = records[j].right;
        self->parent[j] = records[j].parent;
        self->children[j] = records[j].children[0];
        self->children[2 * j] = records[j].children[1];
        self->time[j] = records[j].time;
    }
    ret = 0;
out:
    if (records != NULL) {
        free(records);
    }
    return ret;
}

int
tree_sequence_load(tree_sequence_t *self, const char *filename)
{
    return 0;
}

int
tree_sequence_dump(tree_sequence_t *self, const char *filename)
{
    int ret = -1;
    hid_t file_id, group_id, dataset_id, dataspace_id;
    herr_t status;
    hsize_t dims = self->num_records;

    printf("Dumping tree_sequence to %s\n", filename);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        printf("FIXME\n");
        goto out;
    }
    group_id = H5Gcreate(file_id, "/records", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    printf("group created: %d\n", (int) group_id);
    if (group_id < 0) {
        printf("FIXME\n");
        goto out;
    }

    dataspace_id = H5Screate_simple(1, &dims, NULL);
    if (dataspace_id < 0) {
        printf("FIXME\n");
        goto out;
    }
    /* left */
    dataset_id = H5Dcreate2(file_id, "/records/left", H5T_STD_U32LE,
            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        printf("FIXME\n");
        goto out;
    }
    status = H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, self->left);
    if (status < 0) {
        printf("FIXME\n");
        goto out;
    }
    status = H5Dclose(dataset_id);
    if (status < 0) {
        printf("FIXME\n");
        goto out;
    }
    /* right */
    dataset_id = H5Dcreate2(file_id, "/records/right", H5T_STD_U32LE,
            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        printf("FIXME\n");
        goto out;
    }
    status = H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, self->right);
    if (status < 0) {
        printf("FIXME\n");
        goto out;
    }
    status = H5Dclose(dataset_id);
    if (status < 0) {
        printf("FIXME\n");
        goto out;
    }
    /* parent */
    dataset_id = H5Dcreate2(file_id, "/records/parent", H5T_STD_U32LE,
            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        printf("FIXME\n");
        goto out;
    }
    status = H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, self->parent);
    if (status < 0) {
        printf("FIXME\n");
        goto out;
    }
    status = H5Dclose(dataset_id);
    if (status < 0) {
        printf("FIXME\n");
        goto out;
    }
    /* time */
    dataset_id = H5Dcreate2(file_id, "/records/time", H5T_IEEE_F64LE,
            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        printf("FIXME\n");
        goto out;
    }
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, self->time);
    if (status < 0) {
        printf("FIXME\n");
        goto out;
    }
    status = H5Dclose(dataset_id);
    if (status < 0) {
        printf("FIXME\n");
        goto out;
    }




    status = H5Gclose(group_id);
    if (status < 0) {
        printf("FIXME\n");
        goto out;
    }
    status = H5Sclose(dataspace_id);
    if (status < 0) {
        printf("FIXME\n");
        goto out;
    }
    status = H5Fclose(file_id);
    if (status < 0) {
        printf("FIXME\n");
        goto out;
    }
    ret = 0;
out:
    return ret;
}

int
tree_sequence_free(tree_sequence_t *self)
{
    if (self->left != NULL) {
        free(self->left);
    }
    if (self->right != NULL) {
        free(self->right);
    }
    if (self->children != NULL) {
        free(self->children);
    }
    if (self->parent != NULL) {
        free(self->parent);
    }
    if (self->time != NULL) {
        free(self->time);
    }
    return 0;
}
