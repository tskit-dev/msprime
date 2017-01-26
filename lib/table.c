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
#include <stdbool.h>

#include "err.h"
#include "msprime.h"


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
 * Coordinate table
 *************************/

int
coordinate_table_alloc(coordinate_table_t *self, size_t max_rows_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(coordinate_table_t));
    self->max_rows_increment = max_rows_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    return ret;
}

int
coordinate_table_add_row(coordinate_table_t *self, double position)
{
    int ret = 0;
    size_t new_size;

    if (self->num_rows == self->max_rows) {
        new_size = self->max_rows + self->max_rows_increment;
        ret = expand_column((void **) &self->position, new_size, sizeof(double));
        if (ret != 0) {
            goto out;
        }
        self->max_rows = new_size;
    }
    self->position[self->num_rows] = position;
    self->num_rows++;
out:
    return ret;
}

int
coordinate_table_reset(coordinate_table_t *self)
{
    self->num_rows = 0;
    return 0;
}

int
coordinate_table_free(coordinate_table_t *self)
{
    msp_safe_free(self->position);
    return 0;
}

void
coordinate_table_print_state(coordinate_table_t *self, FILE *out)
{
    size_t j;

    printf("coordinate_table: %p:%d\t%d\t%d\n", (void *) self,
            (int) self->num_rows, (int) self->max_rows, (int) self->max_rows_increment);
    printf("\tindex\tposition\n");
    for (j = 0; j < self->num_rows; j++) {
        printf("\t%d\t%f\n", (int) j, self->position[j]);
    }

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
    self->max_rows_increment = max_rows_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    return ret;
}

int
node_table_add_row(node_table_t *self, uint32_t flags, double time, uint32_t population)
{
    int ret = 0;
    size_t new_size;

    if (self->num_rows == self->max_rows) {
        new_size = self->max_rows + self->max_rows_increment;
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
        size_t max_total_children_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(edgeset_table_t));
    self->max_rows_increment = max_rows_increment;
    self->max_total_children_increment = max_total_children_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    self->max_total_children = 0;
    self->total_children = 0;
    return ret;
}

int
edgeset_table_add_row(edgeset_table_t *self, uint32_t left, uint32_t right,
        uint32_t parent, uint32_t num_children, uint32_t *children)
{
    int ret = 0;
    size_t new_size;

    if (self->num_rows == self->max_rows) {
        new_size = self->max_rows + self->max_rows_increment;
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
    if (self->total_children + num_children >= self->max_total_children) {
        new_size = self->max_total_children + self->max_total_children_increment;
        ret = expand_column((void **) &self->children, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        self->max_total_children = new_size;
    }
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
}
