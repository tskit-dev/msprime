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
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>

#include <gsl/gsl_math.h>

#include "err.h"
#include "msprime.h"

#define DEFAULT_MAX_ROWS_INCREMENT 1024

static int
cmp_double(const void *a, const void *b) {
    const double *ia = (const double *) a;
    const double *ib = (const double *) b;
    return (*ia > *ib) - (*ia < *ib);
}

static int
sort_unique(size_t *num_rows, double *coordinates)
{
    size_t j, k;

    /* TODO need to write some test cases for this, exploring the edge cases */
    qsort(coordinates, *num_rows, sizeof(double), cmp_double);
    /* Now go through the positions and keep the unique values. */
    k = 0;
    for (j = 1; j < *num_rows; j++) {
        if (coordinates[j] != coordinates[k]) {
            k++;
            coordinates[k] = coordinates[j];
        }
    }
    if (*num_rows > 0) {
        k++;
        coordinates[k] = coordinates[*num_rows - 1];
    }
    assert(k <= *num_rows);
    *num_rows = k;
    return 0;
}

static int
get_index(size_t num_rows, double *coordinates, double x, uint32_t *index)
{
    int ret = MSP_ERR_GENERIC;
    double *result;
    ptrdiff_t diff;

    result = bsearch(&x, coordinates, num_rows, sizeof(double), cmp_double);
    if (result == NULL) {
        ret = MSP_ERR_COORDINATE_NOT_FOUND;
        goto out;
    }
    diff = result - coordinates;
    assert(coordinates[diff] == x);
    *index = (uint32_t) diff;
    ret = 0;
out:
    return ret;
}


static int
expand_column(void **column, size_t new_max_rows, size_t element_size)
{
    int ret = 0;
    void *tmp;

    tmp = realloc((void **) *column, new_max_rows * element_size);
    if (tmp == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    *column = tmp;
out:
    return ret;
}

/*************************
 * mutation table
 *************************/

int
mutation_table_alloc(mutation_table_t *self, size_t max_rows_increment,
        size_t max_total_nodes_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(mutation_table_t));
    self->max_rows_increment = max_rows_increment;
    self->max_total_nodes_increment = max_total_nodes_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    self->max_total_nodes = 0;
    self->total_nodes = 0;
    return ret;
}

int
mutation_table_add_row(mutation_table_t *self, double position,
        uint32_t num_nodes, uint32_t *nodes)
{
    int ret = 0;
    size_t new_size;

    if (self->num_rows == self->max_rows) {
        new_size = self->max_rows + self->max_rows_increment;
        ret = expand_column((void **) &self->position, new_size, sizeof(double));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->num_nodes, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        self->max_rows = new_size;
    }
    if (self->total_nodes + num_nodes >= self->max_total_nodes) {
        new_size = self->max_total_nodes + self->max_total_nodes_increment;
        ret = expand_column((void **) &self->nodes, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        self->max_total_nodes = new_size;
    }
    self->position[self->num_rows] = position;
    self->num_nodes[self->num_rows] = num_nodes;
    memcpy(self->nodes + self->total_nodes, nodes, num_nodes * sizeof(uint32_t));
    self->total_nodes += num_nodes;
    self->num_rows++;
out:
    return ret;
}

int
mutation_table_reset(mutation_table_t *self)
{
    self->num_rows = 0;
    self->total_nodes = 0;
    return 0;
}

int
mutation_table_free(mutation_table_t *self)
{
    msp_safe_free(self->position);
    msp_safe_free(self->num_nodes);
    msp_safe_free(self->nodes);
    return 0;
}

void
mutation_table_print_state(mutation_table_t *self, FILE *out)
{
    size_t j, k, offset;

    printf("mutation_table: %p:%d\t%d\t%d\n", (void *) self,
            (int) self->num_rows, (int) self->max_rows, (int) self->max_rows_increment);
    printf("\tindex\tposition\tnodes\n");
    offset = 0;
    for (j = 0; j < self->num_rows; j++) {
        printf("\t%d\t%f\t", (int) j, self->position[j]);
        for (k = 0; k < self->num_nodes[j]; k++) {
            printf("%d", self->nodes[offset]);
            offset++;
            if (k < self->num_nodes[j] - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}

/*************************
 * node table
 *************************/

int
node_table_alloc(node_table_t *self, size_t max_rows_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(node_table_t));
    if (max_rows_increment == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->max_rows_increment = max_rows_increment;
    self->max_rows = 0;
    self->num_rows = 0;
out:
    return ret;
}

static int
node_table_expand(node_table_t *self, size_t new_size)
{
    int ret = 0;

    if (new_size > self->max_rows) {
        ret = expand_column((void **) &self->flags, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->time, new_size, sizeof(double));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->population, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        self->max_rows = new_size;
    }
out:
    return ret;
}

int
node_table_set_columns(node_table_t *self, size_t num_rows, uint32_t *flags, double *time,
        uint32_t *population)
{
    int ret;

    ret = node_table_expand(self, num_rows);
    if (ret != 0) {
        goto out;
    }
    if (flags == NULL || time == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    memcpy(self->flags, flags, num_rows * sizeof(uint32_t));
    memcpy(self->time, time, num_rows * sizeof(double));
    if (population == NULL) {
        memset(self->population, 0xff, num_rows * sizeof(uint32_t));
    } else {
        memcpy(self->population, population, num_rows * sizeof(uint32_t));
    }
    self->num_rows = num_rows;
out:
    return ret;
}

int
node_table_add_row(node_table_t *self, uint32_t flags, double time, uint32_t population)
{
    int ret = 0;
    size_t new_size;

    if (self->num_rows == self->max_rows) {
        new_size = self->max_rows + self->max_rows_increment;
        ret = node_table_expand(self, new_size);
        if (ret != 0) {
            goto out;
        }
    }
    self->flags[self->num_rows] = flags;
    self->time[self->num_rows] = time;
    self->population[self->num_rows] = population;
    self->num_rows++;
out:
    return ret;
}

int
node_table_reset(node_table_t *self)
{
    self->max_rows = 0;
    return 0;
}

int
node_table_free(node_table_t *self)
{
    msp_safe_free(self->flags);
    msp_safe_free(self->time);
    msp_safe_free(self->population);
    return 0;
}

void
node_table_print_state(node_table_t *self, FILE *out)
{
    size_t j;

    printf("node_table: %p:%d\t%d\t%d\n", (void *) self,
            (int) self->num_rows, (int) self->max_rows, (int) self->max_rows_increment);
    printf("\tindex\tflags\ttime\tpopulation\n");
    for (j = 0; j < self->num_rows; j++) {
        printf("\t%d\t%d\t%f\t%d\n", (int) j, self->flags[j], self->time[j],
                self->population[j]);
    }
}

/*************************
 * edgeset table
 *************************/

int
edgeset_table_alloc(edgeset_table_t *self, size_t max_rows_increment,
        size_t max_total_children_increment, size_t max_coordinates_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(edgeset_table_t));
    self->max_rows_increment = max_rows_increment;
    self->max_total_children_increment = max_total_children_increment;
    self->max_coordinates_increment = max_coordinates_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    self->max_total_children = 0;
    self->total_children = 0;
    self->max_coordinates = 0;
    self->num_coordinates = 0;
    return ret;
}

static int
edgeset_table_expand_main_columns(edgeset_table_t *self, size_t new_size)
{
    int ret = 0;

    if (new_size > self->max_rows) {
        ret = expand_column((void **) &self->left, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->right, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->parent, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->num_children, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        self->max_rows = new_size;
    }
out:
    return ret;
}

static int
edgeset_table_expand_children(edgeset_table_t *self, size_t new_size)
{
    int ret = 0;

    if (new_size > self->max_total_children) {
        ret = expand_column((void **) &self->children, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        self->max_total_children = new_size;
    }
out:
    return ret;
}

static int
edgeset_table_expand_coordinates(edgeset_table_t *self, size_t new_size)
{
    int ret = 0;

    if (new_size > self->max_coordinates) {
        ret = expand_column((void **) &self->coordinates, new_size, sizeof(double));
        if (ret != 0) {
            goto out;
        }
        self->max_coordinates = new_size;
    }
out:
    return ret;
}

int
edgeset_table_add_coordinate(edgeset_table_t *self, double coordinate)
{
    int ret = 0;

    if (self->num_coordinates == self->max_coordinates) {
        ret = edgeset_table_expand_coordinates(self,
                self->max_coordinates + self->max_coordinates_increment);
        if (ret != 0) {
            goto out;
        }
    }
    self->coordinates[self->num_coordinates] = coordinate;
    self->num_coordinates++;
out:
    return ret;
}

int
edgeset_table_finalise_coordinates(edgeset_table_t *self)
{
    int ret = 0;
    ret = sort_unique(&self->num_coordinates, self->coordinates);
    return ret;
}

int
edgeset_table_add_row(edgeset_table_t *self, double left_coordinate,
        double right_coordinate, uint32_t parent, uint32_t num_children,
        uint32_t *children)
{
    int ret = 0;
    uint32_t left, right;

    if (self->num_rows == self->max_rows) {
        ret = edgeset_table_expand_main_columns(self,
                self->max_rows + self->max_rows_increment);
        if (ret != 0) {
            goto out;
        }
    }
    /* Need the loop here in case we have a very large number of children */
    while (self->total_children + num_children >= self->max_total_children) {
        ret = edgeset_table_expand_children(self,
            self->max_total_children + self->max_total_children_increment);
        if (ret != 0) {
            goto out;
        }
    }
    ret = get_index(self->num_coordinates, self->coordinates, left_coordinate, &left);
    if (ret != 0) {
        goto out;
    }
    ret = get_index(self->num_coordinates, self->coordinates, right_coordinate, &right);
    if (ret != 0) {
        goto out;
    }
    assert(left < right);
    self->left[self->num_rows] = left;
    self->right[self->num_rows] = right;
    self->parent[self->num_rows] = parent;
    self->num_children[self->num_rows] = num_children;
    memcpy(self->children + self->total_children, children,
            num_children * sizeof(uint32_t));
    self->total_children += num_children;
    self->num_rows++;
out:
    return ret;
}

int
edgeset_table_set_columns(edgeset_table_t *self,
        size_t num_rows, uint32_t *left, uint32_t *right, uint32_t *parent,
        uint32_t *num_children, size_t total_children, uint32_t *children,
        size_t num_coordinates, double *coordinates)
{
    int ret;

    if (left == NULL || right == NULL || parent == NULL || num_children == NULL
            || children == NULL || coordinates == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    ret = edgeset_table_expand_main_columns(self, num_rows);
    if (ret != 0) {
        goto out;
    }
    ret = edgeset_table_expand_children(self, total_children);
    if (ret != 0) {
        goto out;
    }
    ret = edgeset_table_expand_coordinates(self, num_coordinates);
    if (ret != 0) {
        goto out;
    }
    memcpy(self->left, left, num_rows * sizeof(uint32_t));
    memcpy(self->right, right, num_rows * sizeof(uint32_t));
    memcpy(self->parent, parent, num_rows * sizeof(uint32_t));
    memcpy(self->num_children, num_children, num_rows * sizeof(uint32_t));
    memcpy(self->children, children, total_children * sizeof(uint32_t));
    memcpy(self->coordinates, coordinates, num_coordinates * sizeof(double));
    self->num_rows = num_rows;
    self->total_children = total_children;
    self->num_coordinates = num_coordinates;
out:
    return ret;
}

int
edgeset_table_reset(edgeset_table_t *self)
{
    self->num_rows = 0;
    self->total_children = 0;
    return 0;
}

int
edgeset_table_free(edgeset_table_t *self)
{
    msp_safe_free(self->left);
    msp_safe_free(self->right);
    msp_safe_free(self->parent);
    msp_safe_free(self->children);
    msp_safe_free(self->num_children);
    msp_safe_free(self->coordinates);
    return 0;
}

void
edgeset_table_print_state(edgeset_table_t *self, FILE *out)
{
    size_t j, k, offset;

    printf(
        "edgeset_table: %p:num_rows=%d,max_rows=%d,max_rows_increment%d,"
        "total_children=%d,max_total_children=%d,max_total_children_increment=%d\n",
            (void *) self, (int) self->num_rows, (int) self->max_rows,
            (int) self->max_rows_increment, (int) self->total_children,
            (int) self->max_total_children, (int) self->max_total_children_increment);
    printf("\tindex\tleft\tright\tparent\tchildren\n");
    offset = 0;
    for (j = 0; j < self->num_rows; j++) {
        printf("\t%d\t%d\t%d\t%d\t", (int) j, self->left[j], self->right[j],
                self->parent[j]);
        for (k = 0; k < self->num_children[j]; k++) {
            printf("%d", self->children[offset]);
            if (k < self->num_children[j] - 1) {
                printf(",");
            }
            offset++;
        }
        printf("\n");
    }
    printf("\tcoordinates\n");
    for (j = 0; j < self->num_coordinates; j++) {
        printf("\t\t%d\t%f\n", (int) j, self->coordinates[j]);
    }
}
