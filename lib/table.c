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
 * node table
 *************************/

int
node_table_alloc(node_table_t *self, size_t max_rows_increment,
        size_t max_name_length_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(node_table_t));
    if (max_rows_increment == 0 || max_name_length_increment == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->max_rows_increment = max_rows_increment;
    self->max_name_length_increment = max_name_length_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    self->max_name_length = 0;
    self->name_length = 0;
out:
    return ret;
}

static int
node_table_expand_fixed_columns(node_table_t *self, size_t new_size)
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

static int
node_table_expand_name(node_table_t *self, size_t new_size)
{
    int ret = 0;

    if (new_size > self->max_name_length) {
        ret = expand_column((void **) &self->name, new_size, sizeof(char *));
        if (ret != 0) {
            goto out;
        }
        self->max_name_length = new_size;
    }
out:
    return ret;
}

int
node_table_set_columns(node_table_t *self, size_t num_rows, uint32_t *flags, double *time,
        uint32_t *population, size_t name_length, char *name)
{
    int ret;

    if (flags == NULL || time == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (name == NULL && name_length > 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    ret = node_table_expand_fixed_columns(self, num_rows);
    if (ret != 0) {
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
    if (name != NULL) {
        ret = node_table_expand_name(self, name_length);
        if (ret != 0) {
            goto out;
        }
        memcpy(self->name, name, name_length * sizeof(char));
    }
    self->name_length = name_length;

out:
    return ret;
}

int
node_table_add_row(node_table_t *self, uint32_t flags, double time,
        uint32_t population, const char *name)
{
    int ret = 0;
    size_t new_size, name_length;

    if (name == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (self->num_rows == self->max_rows) {
        new_size = self->max_rows + self->max_rows_increment;
        ret = node_table_expand_fixed_columns(self, new_size);
        if (ret != 0) {
            goto out;
        }
    }
    /* include the sentinel \0 */
    name_length = 1 + strlen(name);
    while (self->name_length + name_length
            >= self->max_name_length) {
        new_size = self->max_name_length + self->max_name_length_increment;
        ret = node_table_expand_name(self, new_size);
        if (ret != 0) {
            goto out;
        }
    }
    memcpy(self->name + self->name_length, name, name_length);
    self->name_length += name_length;

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
    self->num_rows = 0;
    self->name_length = 0;
    return 0;
}

int
node_table_free(node_table_t *self)
{
    msp_safe_free(self->flags);
    msp_safe_free(self->time);
    msp_safe_free(self->population);
    msp_safe_free(self->name);
    return 0;
}

void
node_table_print_state(node_table_t *self, FILE *out)
{
    size_t j, k;
    const char *name;

    fprintf(out, "node_table: %p:%d\t%d\t%d\n", (void *) self,
            (int) self->num_rows, (int) self->max_rows, (int) self->max_rows_increment);
    fprintf(out, "\tindex\tflags\ttime\tpopulation\tname\n");
    assert(self->name_length == 0 || self->name_length >= self->num_rows);

    k = 0;
    for (j = 0; j < self->num_rows; j++) {
        fprintf(out, "\t%d\t%d\t%f\t%d\t", (int) j, self->flags[j], self->time[j],
                self->population[j]);
        if (self->name_length == 0) {
            name = "NULL";
        } else {
            assert(k < self->name_length);
            name = self->name + k;
            k += 1 + strlen(name);
        }
        fprintf(out, "%s\n", name);
    }
}

/*************************
 * edgeset table
 *************************/

int
edgeset_table_alloc(edgeset_table_t *self, size_t max_rows_increment,
        size_t max_children_length_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(edgeset_table_t));
    if (max_rows_increment == 0 || max_children_length_increment == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->max_rows_increment = max_rows_increment;
    self->max_children_length_increment = max_children_length_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    self->max_children_length = 0;
    self->children_length = 0;
out:
    return ret;
}

static int
edgeset_table_expand_main_columns(edgeset_table_t *self, size_t new_size)
{
    int ret = 0;

    if (new_size > self->max_rows) {
        ret = expand_column((void **) &self->left, new_size, sizeof(double));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->right, new_size, sizeof(double));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->parent, new_size, sizeof(uint32_t));
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

    if (new_size > self->max_children_length) {
        ret = expand_column((void **) &self->children, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        self->max_children_length = new_size;
    }
out:
    return ret;
}

int
edgeset_table_add_row(edgeset_table_t *self, double left,
        double right, node_id_t parent, uint32_t num_children,
        node_id_t *children)
{
    int ret = 0;

    if (num_children == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (self->num_rows == self->max_rows) {
        ret = edgeset_table_expand_main_columns(self,
                self->max_rows + self->max_rows_increment);
        if (ret != 0) {
            goto out;
        }
    }
    /* Need the loop here in case we have a very large number of children */
    while (1 + self->children_length + num_children >= self->max_children_length) {
        ret = edgeset_table_expand_children(self,
            self->max_children_length + self->max_children_length_increment);
        if (ret != 0) {
            goto out;
        }
    }
    self->left[self->num_rows] = left;
    self->right[self->num_rows] = right;
    self->parent[self->num_rows] = parent;
    memcpy(self->children + self->children_length, children,
            num_children * sizeof(uint32_t));
    self->children[self->children_length + num_children] = MSP_NULL_NODE;
    self->children_length += 1 + num_children;
    self->num_rows++;
out:
    return ret;
}

int
edgeset_table_set_columns(edgeset_table_t *self,
        size_t num_rows, double *left, double *right, node_id_t *parent,
        size_t children_length, node_id_t *children)
{
    int ret;

    if (left == NULL || right == NULL || parent == NULL || children == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    ret = edgeset_table_expand_main_columns(self, num_rows);
    if (ret != 0) {
        goto out;
    }
    ret = edgeset_table_expand_children(self, children_length);
    if (ret != 0) {
        goto out;
    }
    memcpy(self->left, left, num_rows * sizeof(double));
    memcpy(self->right, right, num_rows * sizeof(double));
    memcpy(self->parent, parent, num_rows * sizeof(uint32_t));
    memcpy(self->children, children, children_length * sizeof(uint32_t));
    self->num_rows = num_rows;
    self->children_length = children_length;
out:
    return ret;
}

int
edgeset_table_reset(edgeset_table_t *self)
{
    self->num_rows = 0;
    self->children_length = 0;
    return 0;
}

int
edgeset_table_free(edgeset_table_t *self)
{
    msp_safe_free(self->left);
    msp_safe_free(self->right);
    msp_safe_free(self->parent);
    msp_safe_free(self->children);
    return 0;
}

void
edgeset_table_print_state(edgeset_table_t *self, FILE *out)
{
    size_t j, offset;

    fprintf(out,
        "edgeset_table: %p:num_rows=%d,max_rows=%d,max_rows_increment%d,"
        "children_length=%d,max_children_length=%d,max_children_length_increment=%d\n",
            (void *) self, (int) self->num_rows, (int) self->max_rows,
            (int) self->max_rows_increment, (int) self->children_length,
            (int) self->max_children_length, (int) self->max_children_length_increment);
    fprintf(out, "\tindex\tleft\tright\tparent\tchildren\n");
    offset = 0;
    for (j = 0; j < self->num_rows; j++) {
        fprintf(out, "\t%d\t%.3f\t%.3f\t%d\t", (int) j, self->left[j], self->right[j],
                self->parent[j]);
        while (offset < self->children_length
                && self->children[offset] != MSP_NULL_NODE) {
            fprintf(out, "%d", self->children[offset]);
            assert(offset < self->children_length - 1);
            if (self->children[offset + 1] != MSP_NULL_NODE) {
                fprintf(out, ",");
            }
            offset++;
        }
        offset++;
        fprintf(out, "\n");
    }
    assert(offset == self->children_length);
}

/*************************
 * mutation table
 *************************/

int
mutation_table_alloc(mutation_table_t *self, size_t max_rows_increment,
        size_t max_nodes_length_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(mutation_table_t));
    if (max_rows_increment == 0 || max_nodes_length_increment == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->max_rows_increment = max_rows_increment;
    self->max_nodes_length_increment = max_nodes_length_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    self->max_nodes_length = 0;
    self->nodes_length = 0;
out:
    return ret;
}

static int
mutation_table_expand_main_columns(mutation_table_t *self, size_t new_size)
{
    int ret = 0;

    if (new_size > self->max_rows) {
        ret = expand_column((void **) &self->position, new_size, sizeof(double));
        if (ret != 0) {
            goto out;
        }
        self->max_rows = new_size;
    }
out:
    return ret;
}

static int
mutation_table_expand_nodes(mutation_table_t *self, size_t new_size)
{
    int ret = 0;

    if (new_size > self->max_nodes_length) {
        ret = expand_column((void **) &self->nodes, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        self->max_nodes_length = new_size;
    }
out:
    return ret;
}

int
mutation_table_add_row(mutation_table_t *self, double position,
        uint32_t num_nodes, node_id_t *nodes)
{
    int ret = 0;
    size_t new_size;

    if (num_nodes == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (self->num_rows == self->max_rows) {
        new_size = self->max_rows + self->max_rows_increment;
        ret = mutation_table_expand_main_columns(self, new_size);
        if (ret != 0) {
            goto out;
        }
    }
    while (1 + self->nodes_length + num_nodes >= self->max_nodes_length) {
        new_size = 1 + self->max_nodes_length + self->max_nodes_length_increment;
        ret = mutation_table_expand_nodes(self, new_size);
        if (ret != 0) {
            goto out;
        }
    }
    self->position[self->num_rows] = position;
    memcpy(self->nodes + self->nodes_length, nodes, num_nodes * sizeof(uint32_t));
    self->nodes[self->nodes_length + num_nodes] = MSP_NULL_NODE;
    self->nodes_length += 1 + num_nodes;
    self->num_rows++;
out:
    return ret;
}

int
mutation_table_set_columns(mutation_table_t *self, size_t num_rows, double *position,
        size_t nodes_length, node_id_t *nodes)
{
    int ret = 0;

    if (position == NULL || nodes == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    ret = mutation_table_expand_main_columns(self, num_rows);
    if (ret != 0) {
        goto out;
    }
    ret = mutation_table_expand_nodes(self, nodes_length);
    if (ret != 0) {
        goto out;
    }
    memcpy(self->position, position, num_rows * sizeof(double));
    memcpy(self->nodes, nodes, nodes_length * sizeof(uint32_t));
    self->num_rows = num_rows;
    self->nodes_length = nodes_length;
out:
    return ret;
}

int
mutation_table_reset(mutation_table_t *self)
{
    self->num_rows = 0;
    self->nodes_length = 0;
    return 0;
}

int
mutation_table_free(mutation_table_t *self)
{
    msp_safe_free(self->position);
    msp_safe_free(self->nodes);
    return 0;
}

void
mutation_table_print_state(mutation_table_t *self, FILE *out)
{
    size_t j, offset;

    fprintf(out, "mutation_table: %p:%d\t%d\t%d\n", (void *) self,
            (int) self->num_rows, (int) self->max_rows, (int) self->max_rows_increment);
    fprintf(out, "\tindex\tposition\tnodes\n");
    offset = 0;
    for (j = 0; j < self->num_rows; j++) {
        fprintf(out, "\t%d\t%f\t", (int) j, self->position[j]);
        while (self->nodes[offset] != MSP_NULL_NODE) {
            fprintf(out, "%d", self->nodes[offset]);
            offset++;
            if (self->nodes[offset] != MSP_NULL_NODE) {
                fprintf(out, ",");
            }
        }
        offset++;
        fprintf(out, "\n");
    }
}

/*************************
 * migration table
 *************************/

int
migration_table_alloc(migration_table_t *self, size_t max_rows_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(migration_table_t));
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
migration_table_expand(migration_table_t *self, size_t new_size)
{
    int ret = 0;

    if (new_size > self->max_rows) {
        ret = expand_column((void **) &self->left, new_size, sizeof(double));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->right, new_size, sizeof(double));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->node, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->source, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->dest, new_size, sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->time, new_size, sizeof(double));
        if (ret != 0) {
            goto out;
        }
        self->max_rows = new_size;
    }
out:
    return ret;
}

int
migration_table_set_columns(migration_table_t *self, size_t num_rows, double *left,
        double *right, node_id_t *node, uint32_t *source, uint32_t *dest, double *time)
{
    int ret;

    ret = migration_table_expand(self, num_rows);
    if (ret != 0) {
        goto out;
    }
    if (left == NULL || right == NULL || node == NULL || source == NULL
            || dest == NULL || time == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    memcpy(self->left, left, num_rows * sizeof(double));
    memcpy(self->right, right, num_rows * sizeof(double));
    memcpy(self->node, node, num_rows * sizeof(uint32_t));
    memcpy(self->source, source, num_rows * sizeof(uint32_t));
    memcpy(self->dest, dest, num_rows * sizeof(uint32_t));
    memcpy(self->time, time, num_rows * sizeof(double));
    self->num_rows = num_rows;
out:
    return ret;
}

int
migration_table_add_row(migration_table_t *self, double left, double right,
        node_id_t node, uint32_t source, uint32_t dest, double time)
{
    int ret = 0;
    size_t new_size;

    if (self->num_rows == self->max_rows) {
        new_size = self->max_rows + self->max_rows_increment;
        ret = migration_table_expand(self, new_size);
        if (ret != 0) {
            goto out;
        }
    }
    self->left[self->num_rows] = left;
    self->right[self->num_rows] = right;
    self->node[self->num_rows] = node;
    self->source[self->num_rows] = source;
    self->dest[self->num_rows] = dest;
    self->time[self->num_rows] = time;
    self->num_rows++;
out:
    return ret;
}

int
migration_table_reset(migration_table_t *self)
{
    self->num_rows = 0;
    return 0;
}

int
migration_table_free(migration_table_t *self)
{
    msp_safe_free(self->left);
    msp_safe_free(self->right);
    msp_safe_free(self->node);
    msp_safe_free(self->source);
    msp_safe_free(self->dest);
    msp_safe_free(self->time);
    return 0;
}

void
migration_table_print_state(migration_table_t *self, FILE *out)
{
    size_t j;

    fprintf(out, "migration_table: %p:%d\t%d\t%d\n", (void *) self,
            (int) self->num_rows, (int) self->max_rows, (int) self->max_rows_increment);
    fprintf(out, "\tindex\tleft\tright\tnode\tsource\tdest\tpopulation\n");
    for (j = 0; j < self->num_rows; j++) {
        fprintf(out, "\t%d\t%.3f\t%.3f\t%d\t%d\t%d\t%f\n", (int) j, self->left[j],
                self->right[j], self->node[j], self->source[j], self->dest[j],
                self->time[j]);
    }
}
