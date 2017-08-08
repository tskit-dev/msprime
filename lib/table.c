/*
** Copyright (C) 2015-2016 University of Oxford
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
#include "object_heap.h"
#include "msprime.h"

#define DEFAULT_MAX_ROWS_INCREMENT 1024

#define TABLE_SEP "-----------------------------------------\n"


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
        size_t max_total_name_length_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(node_table_t));
    if (max_rows_increment == 0 || max_total_name_length_increment == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->max_rows_increment = max_rows_increment;
    self->max_total_name_length_increment = max_total_name_length_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    self->max_total_name_length = 0;
    self->total_name_length = 0;
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
        ret = expand_column((void **) &self->population, new_size,
                sizeof(population_id_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->name_length, new_size, sizeof(uint32_t));
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

    if (new_size > self->max_total_name_length) {
        ret = expand_column((void **) &self->name, new_size, sizeof(char *));
        if (ret != 0) {
            goto out;
        }
        self->max_total_name_length = new_size;
    }
out:
    return ret;
}

int
node_table_set_columns(node_table_t *self, size_t num_rows, uint32_t *flags, double *time,
        population_id_t *population, char *name, uint32_t *name_length)
{
    int ret;
    size_t j, total_name_length;

    if (flags == NULL || time == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if ((name == NULL) != (name_length == NULL)) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    ret = node_table_expand_fixed_columns(self, num_rows);
    if (ret != 0) {
        goto out;
    }
    memcpy(self->flags, flags, num_rows * sizeof(uint32_t));
    memcpy(self->time, time, num_rows * sizeof(double));
    if (name == NULL) {
        self->total_name_length = 0;
        memset(self->name_length, 0, num_rows * sizeof(uint32_t));
    } else {
        memcpy(self->name_length, name_length, num_rows * sizeof(uint32_t));
        total_name_length = 0;
        for (j = 0; j < num_rows; j++) {
            total_name_length += name_length[j];
        }
        ret = node_table_expand_name(self, total_name_length);
        if (ret != 0) {
            goto out;
        }
        memcpy(self->name, name, total_name_length * sizeof(char));
        self->total_name_length = total_name_length;
    }
    if (population == NULL) {
        memset(self->population, 0xff, num_rows * sizeof(population_id_t));
    } else {
        memcpy(self->population, population, num_rows * sizeof(population_id_t));
    }
    self->num_rows = num_rows;
out:
    return ret;
}

int
node_table_add_row(node_table_t *self, uint32_t flags, double time,
        population_id_t population, const char *name)
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
    name_length = strlen(name);
    while (self->total_name_length + name_length >= self->max_total_name_length) {
        new_size = self->max_total_name_length + self->max_total_name_length_increment;
        ret = node_table_expand_name(self, new_size);
        if (ret != 0) {
            goto out;
        }
    }
    memcpy(self->name + self->total_name_length, name, name_length);
    self->total_name_length += name_length;
    self->flags[self->num_rows] = flags;
    self->time[self->num_rows] = time;
    self->population[self->num_rows] = population;
    self->name_length[self->num_rows] = (uint32_t) name_length;
    self->num_rows++;
out:
    return ret;
}

int
node_table_reset(node_table_t *self)
{
    self->num_rows = 0;
    self->total_name_length = 0;
    return 0;
}

int
node_table_free(node_table_t *self)
{
    msp_safe_free(self->flags);
    msp_safe_free(self->time);
    msp_safe_free(self->population);
    msp_safe_free(self->name);
    msp_safe_free(self->name_length);
    return 0;
}

void
node_table_print_state(node_table_t *self, FILE *out)
{
    size_t j, k, offset;

    fprintf(out, TABLE_SEP);
    fprintf(out, "node_table: %p:\n", (void *) self);
    fprintf(out, "num_rows          = %d\tmax= %d\tincrement = %d)\n",
            (int) self->num_rows, (int) self->max_rows, (int) self->max_rows_increment);
    fprintf(out, "total_name_length = %d\tmax= %d\tincrement = %d)\n",
            (int) self->total_name_length,
            (int) self->max_total_name_length,
            (int) self->max_total_name_length_increment);
    fprintf(out, TABLE_SEP);
    fprintf(out, "index\tflags\ttime\tpopulation\tname_length\tname\n");
    offset = 0;
    for (j = 0; j < self->num_rows; j++) {
        fprintf(out, "%d\t%d\t%f\t%d\t%d\t", (int) j, self->flags[j], self->time[j],
                (int) self->population[j], self->name_length[j]);
        for (k = 0; k < self->name_length[j]; k++) {
            assert(offset < self->total_name_length);
            fprintf(out, "%c", self->name[offset]);
            offset++;
        }
        fprintf(out, "\n");
    }
}

/*************************
 * edgeset table
 *************************/

int
edgeset_table_alloc(edgeset_table_t *self, size_t max_rows_increment,
        size_t max_total_children_length_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(edgeset_table_t));
    if (max_rows_increment == 0 || max_total_children_length_increment == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->max_rows_increment = max_rows_increment;
    self->max_total_children_length_increment = max_total_children_length_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    self->max_total_children_length = 0;
    self->total_children_length = 0;
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
        ret = expand_column((void **) &self->parent, new_size, sizeof(node_id_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->children_length, new_size,
                sizeof(list_len_t));
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

    if (new_size > self->max_total_children_length) {
        ret = expand_column((void **) &self->children, new_size, sizeof(node_id_t));
        if (ret != 0) {
            goto out;
        }
        self->max_total_children_length = new_size;
    }
out:
    return ret;
}

int
edgeset_table_add_row(edgeset_table_t *self, double left, double right,
        node_id_t parent, node_id_t *children, list_len_t children_length)
{
    int ret = 0;

    if (children_length == 0) {
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
    while (self->total_children_length + children_length
            >= self->max_total_children_length) {
        ret = edgeset_table_expand_children(self,
            self->max_total_children_length + self->max_total_children_length_increment);
        if (ret != 0) {
            goto out;
        }
    }
    self->left[self->num_rows] = left;
    self->right[self->num_rows] = right;
    self->parent[self->num_rows] = parent;
    memcpy(self->children + self->total_children_length, children,
            children_length * sizeof(node_id_t));
    self->children_length[self->num_rows] = children_length;
    self->total_children_length += children_length;
    self->num_rows++;
out:
    return ret;
}

int
edgeset_table_set_columns(edgeset_table_t *self,
        size_t num_rows, double *left, double *right, node_id_t *parent,
        node_id_t *children, list_len_t *children_length)
{
    int ret;
    list_len_t j;
    size_t total_children_length = 0;

    if (left == NULL || right == NULL || parent == NULL || children == NULL
            || children_length == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    for (j = 0; j < num_rows; j++) {
        total_children_length += children_length[j];
    }
    ret = edgeset_table_expand_main_columns(self, num_rows);
    if (ret != 0) {
        goto out;
    }
    ret = edgeset_table_expand_children(self, total_children_length);
    if (ret != 0) {
        goto out;
    }
    memcpy(self->left, left, num_rows * sizeof(double));
    memcpy(self->right, right, num_rows * sizeof(double));
    memcpy(self->parent, parent, num_rows * sizeof(node_id_t));
    memcpy(self->children, children, total_children_length * sizeof(node_id_t));
    memcpy(self->children_length, children_length, num_rows * sizeof(list_len_t));
    self->num_rows = num_rows;
    self->total_children_length += total_children_length;
out:
    return ret;
}

int
edgeset_table_reset(edgeset_table_t *self)
{
    self->num_rows = 0;
    self->total_children_length = 0;
    return 0;
}

int
edgeset_table_free(edgeset_table_t *self)
{
    msp_safe_free(self->left);
    msp_safe_free(self->right);
    msp_safe_free(self->parent);
    msp_safe_free(self->children);
    msp_safe_free(self->children_length);
    return 0;
}

void
edgeset_table_print_state(edgeset_table_t *self, FILE *out)
{
    size_t j, offset;
    list_len_t k;

    fprintf(out, TABLE_SEP);
    fprintf(out, "edgeset_table: %p:\n", (void *) self);
    fprintf(out, "num_rows          = %d\tmax= %d\tincrement = %d)\n",
            (int) self->num_rows, (int) self->max_rows, (int) self->max_rows_increment);
    fprintf(out, "total_children_length   = %d\tmax= %d\tincrement = %d)\n",
            (int) self->total_children_length,
            (int) self->max_total_children_length,
            (int) self->max_total_children_length_increment);
    fprintf(out, TABLE_SEP);
    fprintf(out, "index\tleft\tright\tparent\tchildren_length\tchildren\n");
    offset = 0;
    for (j = 0; j < self->num_rows; j++) {
        fprintf(out, "%d\t%.3f\t%.3f\t%d\t%d\t", (int) j, self->left[j], self->right[j],
                (int) self->parent[j], self->children_length[j]);
        for (k = 0; k < self->children_length[j]; k++) {
            assert(offset < self->total_children_length);
            fprintf(out, "%d", (int) self->children[offset]);
            offset++;
            if (k < self->children_length[j] - 1) {
                fprintf(out, ",");
            }
        }
        fprintf(out, "\n");
    }
    assert(offset == self->total_children_length);
}

/*************************
 * site table
 *************************/

int
site_table_alloc(site_table_t *self, size_t max_rows_increment,
        size_t max_total_ancestral_state_length_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(site_table_t));
    if (max_rows_increment == 0 || max_total_ancestral_state_length_increment == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->max_rows_increment = max_rows_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    self->max_total_ancestral_state_length_increment =
        max_total_ancestral_state_length_increment;
    self->max_total_ancestral_state_length = 0;
    self->total_ancestral_state_length = 0;
out:
    return ret;
}

static int
site_table_expand_main_columns(site_table_t *self, size_t new_size)
{
    int ret = 0;

    if (new_size > self->max_rows) {
        ret = expand_column((void **) &self->position, new_size, sizeof(double));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->ancestral_state_length, new_size,
                sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        self->max_rows = new_size;
    }
out:
    return ret;
}

static int
site_table_expand_ancestral_state(site_table_t *self, size_t new_size)
{
    int ret = 0;

    if (new_size > self->max_total_ancestral_state_length) {
        ret = expand_column((void **) &self->ancestral_state, new_size, sizeof(char));
        if (ret != 0) {
            goto out;
        }
        self->max_total_ancestral_state_length = new_size;
    }
out:
    return ret;
}

int
site_table_add_row(site_table_t *self, double position, const char *ancestral_state,
        list_len_t ancestral_state_length)
{
    int ret = 0;
    size_t new_size;

    if (self->num_rows == self->max_rows) {
        new_size = self->max_rows + self->max_rows_increment;
        ret = site_table_expand_main_columns(self, new_size);
        if (ret != 0) {
            goto out;
        }
    }
    while (self->total_ancestral_state_length + ancestral_state_length >=
            self->max_total_ancestral_state_length) {
        ret = site_table_expand_ancestral_state(self,
                self->max_total_ancestral_state_length +
                self->max_total_ancestral_state_length_increment);
        if (ret != 0) {
            goto out;
        }
    }
    self->position[self->num_rows] = position;
    self->ancestral_state_length[self->num_rows] = (uint32_t) ancestral_state_length;
    memcpy(self->ancestral_state + self->total_ancestral_state_length,
            ancestral_state, ancestral_state_length);
    self->total_ancestral_state_length += ancestral_state_length;
    self->num_rows++;
out:
    return ret;
}

int
site_table_set_columns(site_table_t *self, size_t num_rows, double *position,
        const char *ancestral_state, list_len_t *ancestral_state_length)
{
    int ret = 0;
    size_t total_ancestral_state_length = 0;
    size_t j;

    if (position == NULL || ancestral_state == NULL || ancestral_state_length == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }

    for (j = 0; j < num_rows; j++) {
        total_ancestral_state_length += ancestral_state_length[j];
    }
    ret = site_table_expand_main_columns(self, num_rows);
    if (ret != 0) {
        goto out;
    }
    ret = site_table_expand_ancestral_state(self, total_ancestral_state_length);
    if (ret != 0) {
        goto out;
    }
    memcpy(self->position, position, num_rows * sizeof(double));
    memcpy(self->ancestral_state, ancestral_state,
            total_ancestral_state_length * sizeof(char));
    memcpy(self->ancestral_state_length, ancestral_state_length,
            num_rows * sizeof(uint32_t));
    self->num_rows = num_rows;
    self->total_ancestral_state_length = total_ancestral_state_length;
out:
    return ret;
}

bool
site_table_equal(site_table_t *self, site_table_t *other)
{
    bool ret = false;
    if (self->num_rows == other->num_rows
            && self->total_ancestral_state_length == other->total_ancestral_state_length) {
        ret = memcmp(self->position, other->position,
                self->num_rows * sizeof(double)) == 0
            && memcmp(self->ancestral_state_length, other->ancestral_state_length,
                    self->num_rows * sizeof(list_len_t)) == 0
            && memcmp(self->ancestral_state, other->ancestral_state,
                    self->total_ancestral_state_length * sizeof(char)) == 0;
    }
    return ret;
}

int
site_table_reset(site_table_t *self)
{
    self->num_rows = 0;
    self->total_ancestral_state_length = 0;
    return 0;
}

int
site_table_free(site_table_t *self)
{
    msp_safe_free(self->position);
    msp_safe_free(self->ancestral_state);
    msp_safe_free(self->ancestral_state_length);
    return 0;
}

void
site_table_print_state(site_table_t *self, FILE *out)
{
    size_t j, k, ancestral_state_offset;

    fprintf(out, TABLE_SEP);
    fprintf(out, "site_table: %p:\n", (void *) self);
    fprintf(out, "num_rows = %d\tmax= %d\tincrement = %d)\n",
            (int) self->num_rows, (int) self->max_rows, (int) self->max_rows_increment);
    fprintf(out, "total_ancestral_state_length = %d\tmax= %d\tincrement = %d)\n",
            (int) self->total_ancestral_state_length,
            (int) self->max_total_ancestral_state_length,
            (int) self->max_total_ancestral_state_length_increment);
    fprintf(out, TABLE_SEP);
    fprintf(out, "index\tposition\tancestral_state_length\tancestral_state\n");
    ancestral_state_offset = 0;
    for (j = 0; j < self->num_rows; j++) {
        fprintf(out, "%d\t%f\t%d\t", (int) j, self->position[j],
                self->ancestral_state_length[j]);
        for (k = 0; k < self->ancestral_state_length[j]; k++) {
            fprintf(out, "%c", self->ancestral_state[ancestral_state_offset]);
            ancestral_state_offset++;
        }
        fprintf(out, "\n");
    }
}

/*************************
 * mutation table
 *************************/

int
mutation_table_alloc(mutation_table_t *self, size_t max_rows_increment,
        size_t max_total_derived_state_length_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(mutation_table_t));
    if (max_rows_increment == 0 || max_total_derived_state_length_increment == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->max_rows_increment = max_rows_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    self->max_total_derived_state_length_increment =
        max_total_derived_state_length_increment;
    self->max_total_derived_state_length = 0;
    self->total_derived_state_length = 0;
out:
    return ret;
}

static int
mutation_table_expand_main_columns(mutation_table_t *self, size_t new_size)
{
    int ret = 0;

    if (new_size > self->max_rows) {
        ret = expand_column((void **) &self->site, new_size, sizeof(site_id_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->node, new_size, sizeof(node_id_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->derived_state_length, new_size,
                sizeof(uint32_t));
        if (ret != 0) {
            goto out;
        }
        self->max_rows = new_size;
    }
out:
    return ret;
}

static int
mutation_table_expand_derived_state(mutation_table_t *self, size_t new_size)
{
    int ret = 0;

    if (new_size > self->max_total_derived_state_length) {
        ret = expand_column((void **) &self->derived_state, new_size, sizeof(char));
        if (ret != 0) {
            goto out;
        }
        self->max_total_derived_state_length = new_size;
    }
out:
    return ret;
}

int
mutation_table_add_row(mutation_table_t *self, site_id_t site, node_id_t node,
        const char *derived_state, list_len_t derived_state_length)
{
    int ret = 0;
    size_t new_size;

    if (self->num_rows == self->max_rows) {
        new_size = self->max_rows + self->max_rows_increment;
        ret = mutation_table_expand_main_columns(self, new_size);
        if (ret != 0) {
            goto out;
        }
    }
    while (self->total_derived_state_length + derived_state_length >=
            self->max_total_derived_state_length) {
        ret = mutation_table_expand_derived_state(self,
                self->max_total_derived_state_length +
                self->max_total_derived_state_length_increment);
        if (ret != 0) {
            goto out;
        }
    }
    self->site[self->num_rows] = site;
    self->node[self->num_rows] = node;
    self->derived_state_length[self->num_rows] = (list_len_t) derived_state_length;
    memcpy(self->derived_state + self->total_derived_state_length, derived_state,
            derived_state_length * sizeof(char));
    self->total_derived_state_length += derived_state_length;
    self->num_rows++;
out:
    return ret;
}

int
mutation_table_set_columns(mutation_table_t *self, size_t num_rows, site_id_t *site,
        node_id_t *node, const char *derived_state, uint32_t *derived_state_length)
{
    int ret = 0;
    size_t total_derived_state_length = 0;
    size_t j;

    if (site == NULL || node == NULL || derived_state == NULL
            || derived_state_length == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    ret = mutation_table_expand_main_columns(self, num_rows);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < num_rows; j++) {
        total_derived_state_length += (size_t) derived_state_length[j];
    }
    ret = mutation_table_expand_derived_state(self, total_derived_state_length);
    if (ret != 0) {
        goto out;
    }
    memcpy(self->site, site, num_rows * sizeof(site_id_t));
    memcpy(self->node, node, num_rows * sizeof(node_id_t));
    memcpy(self->derived_state_length, derived_state_length, num_rows * sizeof(node_id_t));
    memcpy(self->derived_state, derived_state, total_derived_state_length * sizeof(char));
    self->num_rows = num_rows;
    self->total_derived_state_length = total_derived_state_length;
out:
    return ret;
}

bool
mutation_table_equal(mutation_table_t *self, mutation_table_t *other)
{
    bool ret = false;
    if (self->num_rows == other->num_rows
            && self->total_derived_state_length == other->total_derived_state_length) {
        ret = memcmp(self->site, other->site, self->num_rows * sizeof(site_id_t)) == 0
            && memcmp(self->node, other->node, self->num_rows * sizeof(node_id_t)) == 0
            && memcmp(self->derived_state_length, other->derived_state_length,
                    self->num_rows * sizeof(list_len_t)) == 0
            && memcmp(self->derived_state, other->derived_state,
                    self->total_derived_state_length * sizeof(char)) == 0;
    }
    return ret;
}

int
mutation_table_reset(mutation_table_t *self)
{
    self->num_rows = 0;
    self->total_derived_state_length = 0;
    return 0;
}

int
mutation_table_free(mutation_table_t *self)
{
    msp_safe_free(self->node);
    msp_safe_free(self->site);
    msp_safe_free(self->derived_state);
    msp_safe_free(self->derived_state_length);
    return 0;
}

void
mutation_table_print_state(mutation_table_t *self, FILE *out)
{
    size_t j, k, offset;

    fprintf(out, TABLE_SEP);
    fprintf(out, "mutation_table: %p:\n", (void *) self);
    fprintf(out, "num_rows = %d\tmax= %d\tincrement = %d)\n",
            (int) self->num_rows, (int) self->max_rows, (int) self->max_rows_increment);
    fprintf(out, "derived_state_length = %d\tmax= %d\tincrement = %d)\n",
            (int) self->total_derived_state_length,
            (int) self->max_total_derived_state_length,
            (int) self->max_total_derived_state_length_increment);
    fprintf(out, TABLE_SEP);
    fprintf(out, "index\tsite\tnode\tderived_state_length\tderived_state\n");
    offset = 0;
    for (j = 0; j < self->num_rows; j++) {
        fprintf(out, "%d\t%d\t%d\t%d\t", (int) j, self->site[j], self->node[j],
                self->derived_state_length[j]);
        for (k = 0; k < self->derived_state_length[j]; k++) {
            fprintf(out, "%c", self->derived_state[offset]);
            offset++;
        }
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
        ret = expand_column((void **) &self->node, new_size, sizeof(node_id_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->source, new_size, sizeof(population_id_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->dest, new_size, sizeof(population_id_t));
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
        double *right, node_id_t *node, population_id_t *source, population_id_t *dest,
        double *time)
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
    memcpy(self->node, node, num_rows * sizeof(node_id_t));
    memcpy(self->source, source, num_rows * sizeof(population_id_t));
    memcpy(self->dest, dest, num_rows * sizeof(population_id_t));
    memcpy(self->time, time, num_rows * sizeof(double));
    self->num_rows = num_rows;
out:
    return ret;
}

int
migration_table_add_row(migration_table_t *self, double left, double right,
        node_id_t node, population_id_t source, population_id_t dest, double time)
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

    fprintf(out, TABLE_SEP);
    fprintf(out, "migration_table: %p:\n", (void *) self);
    fprintf(out, "num_rows = %d\tmax= %d\tincrement = %d)\n",
            (int) self->num_rows, (int) self->max_rows, (int) self->max_rows_increment);
    fprintf(out, TABLE_SEP);
    fprintf(out, "index\tleft\tright\tnode\tsource\tdest\ttime\n");
    for (j = 0; j < self->num_rows; j++) {
        fprintf(out, "%d\t%.3f\t%.3f\t%d\t%d\t%d\t%f\n", (int) j, self->left[j],
                self->right[j], (int) self->node[j], (int) self->source[j],
                (int) self->dest[j], self->time[j]);
    }
}

/*************************
 * simplifier
 *************************/

static int
cmp_node_id(const void *a, const void *b) {
    const node_id_t *ia = (const node_id_t *) a;
    const node_id_t *ib = (const node_id_t *) b;
    return (*ia > *ib) - (*ia < *ib);
}

static int
cmp_overlap_count(const void *a, const void *b) {
    const overlap_count_t *ia = (const overlap_count_t *) a;
    const overlap_count_t *ib = (const overlap_count_t *) b;
    int ret = (ia->start > ib->start) - (ia->start < ib->start);
    return ret;
}

/* For the segment priority queue we want to sort on the left
 * coordinate and to break ties we use the node */
static int
cmp_segment_queue(const void *a, const void *b) {
    const simplify_segment_t *ia = (const simplify_segment_t *) a;
    const simplify_segment_t *ib = (const simplify_segment_t *) b;
    int ret = (ia->left > ib->left) - (ia->left < ib->left);
    if (ret == 0)  {
        ret = (ia->node > ib->node) - (ia->node < ib->node);
    }
    return ret;
}

static void
simplifier_check_state(simplifier_t *self)
{
    size_t j;
    size_t total_segments = 0;
    avl_node_t *avl_node;
    simplify_segment_t *u;

    for (j = 0; j < self->input_nodes.num_rows; j++) {
        if (self->ancestor_map[j] != NULL) {
            for (u = self->ancestor_map[j]; u != NULL; u = u->next) {
                assert(u->left < u->right);
                if (u->next != NULL) {
                    if (u->right > u->next->left) {
                        printf("ERROR!! %p %p\n", (void*) u, (void*) u->next);
                    }
                    assert(u->right <= u->next->left);
                }
                total_segments++;
            }
        }
    }

    for (avl_node = self->merge_queue.head; avl_node != NULL; avl_node = avl_node->next) {
        for (u = (simplify_segment_t *) avl_node->item; u != NULL; u = u->next) {
            assert(u->left < u->right);
            if (u->next != NULL) {
                if (u->right > u->next->left) {
                    printf("ERROR Q!! %p %p\n", (void*) u, (void*) u->next);
                }
                assert(u->right <= u->next->left);
            }
            total_segments++;
        }
    }
    /* if (total_segments != object_heap_get_num_allocated(&self->segment_heap)) { */
    /*     printf("MISMATCH IN SEGMENTS: present = %d, alloced = %d\n", (int) total_segments, */
    /*             (int)object_heap_get_num_allocated(&self->segment_heap)); */
    /* } */
    /* assert(total_segments == object_heap_get_num_allocated(&self->segment_heap)); */
}

static void
print_segment_chain(simplify_segment_t *head)
{
    simplify_segment_t *u;

    for (u = head; u != NULL; u = u->next) {
        printf("(%f,%f->%d)", u->left, u->right, u->node);
    }
}

void
simplifier_print_state(simplifier_t *self)
{
    size_t j;
    avl_node_t *avl_node;
    simplify_segment_t *u;

    printf("--simplifier state--\n");
    printf("===\nInput nodes\n==\n");
    node_table_print_state(&self->input_nodes, stdout);
    printf("===\nOutput tables\n==\n");
    node_table_print_state(self->nodes, stdout);
    edgeset_table_print_state(self->edgesets, stdout);
    site_table_print_state(self->sites, stdout);
    mutation_table_print_state(self->mutations, stdout);
    printf("===\nmemory heaps\n==\n");
    printf("segment_heap:\n");
    object_heap_print_state(&self->segment_heap, stdout);
    printf("avl_node_heap:\n");
    object_heap_print_state(&self->avl_node_heap, stdout);
    printf("overlap_count_heap:\n");
    object_heap_print_state(&self->overlap_count_heap, stdout);
    printf("===\nancestors\n==\n");
    for (j = 0; j < self->input_nodes.num_rows; j++) {
        if (self->ancestor_map[j] != NULL) {
            printf("%d:\t", (int) j);
            print_segment_chain(self->ancestor_map[j]);
            printf("\n");
        }
    }
    printf("===\nmerge queue\n==\n");
    for (avl_node = self->merge_queue.head; avl_node != NULL; avl_node = avl_node->next) {
        u = (simplify_segment_t *) avl_node->item;
        print_segment_chain(u);
        printf("\n");
    }
}

static simplify_segment_t * WARN_UNUSED
simplifier_alloc_segment(simplifier_t *self, double left, double right, node_id_t node,
        simplify_segment_t *prev, simplify_segment_t *next)
{
    simplify_segment_t *seg = NULL;

    if (object_heap_empty(&self->segment_heap)) {
        if (object_heap_expand(&self->segment_heap) != 0) {
            goto out;
        }
    }
    seg = (simplify_segment_t *) object_heap_alloc_object(&self->segment_heap);
    if (seg == NULL) {
        goto out;
    }
    seg->prev = prev;
    seg->next = next;
    seg->left = left;
    seg->right = right;
    seg->node = node;
out:
    return seg;
}

static inline void
simplifier_free_segment(simplifier_t *self, simplify_segment_t *seg)
{
    object_heap_free_object(&self->segment_heap, seg);
}

static inline avl_node_t * WARN_UNUSED
simplifier_alloc_avl_node(simplifier_t *self)
{
    avl_node_t *ret = NULL;

    if (object_heap_empty(&self->avl_node_heap)) {
        if (object_heap_expand(&self->avl_node_heap) != 0) {
            goto out;
        }
    }
    ret = (avl_node_t *) object_heap_alloc_object(&self->avl_node_heap);
out:
    return ret;
}

static inline void
simplifier_free_avl_node(simplifier_t *self, avl_node_t *node)
{
    object_heap_free_object(&self->avl_node_heap, node);
}

static inline overlap_count_t *
simplifier_alloc_overlap_count(simplifier_t *self)
{
    overlap_count_t *ret = NULL;

    if (object_heap_empty(&self->overlap_count_heap)) {
        if (object_heap_expand(&self->overlap_count_heap) != 0) {
            goto out;
        }
    }
    ret = (overlap_count_t *) object_heap_alloc_object(&self->overlap_count_heap);
out:
    return ret;
}

/*
 * Inserts a new overlap_count at the specified position, mapping to the
 * specified number of overlapping segments v.
 */
static int WARN_UNUSED
simplifier_insert_overlap_count(simplifier_t *self, double x, uint32_t v)
{
    int ret = 0;
    avl_node_t *node = simplifier_alloc_avl_node(self);
    overlap_count_t *m = simplifier_alloc_overlap_count(self);

    if (node == NULL || m == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    m->start = x;
    m->count = v;
    avl_init_node(node, m);
    node = avl_insert_node(&self->overlap_counts, node);
    assert(node != NULL);
out:
    return ret;
}

/*
 * Inserts a new overlap_count at the specified position, and copies its
 * count from the containing overlap count.
 */
static int WARN_UNUSED
simplifier_copy_overlap_count(simplifier_t *self, double x)
{
    int ret;
    overlap_count_t search, *nm;
    avl_node_t *node;

    search.start = x;
    avl_search_closest(&self->overlap_counts, &search, &node);
    assert(node != NULL);
    nm = (overlap_count_t *) node->item;
    if (nm->start > x) {
        node = node->prev;
        assert(node != NULL);
        nm = (overlap_count_t *) node->item;
    }
    ret = simplifier_insert_overlap_count(self, x, nm->count);
    return ret;
}



/* Add a new node to the output node table corresponding to the specified input id */
static int
simplifier_record_node(simplifier_t *self, node_id_t input_id)
{
    int ret = 0;

    /* FIXME!! This does not transfer the name properly. */
    ret = node_table_add_row(self->nodes, self->input_nodes.flags[input_id],
            self->input_nodes.time[input_id], self->input_nodes.population[input_id],
            "");
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

/* Records the specified edgeset in the output table */
static int
simplifier_record_edgeset(simplifier_t *self, double left, double right,
        node_id_t parent, node_id_t *children, list_len_t children_length)
{
    int ret = 0;
    bool squash;

    qsort(children, children_length, sizeof(node_id_t), cmp_node_id);
    if (self->last_edgeset.children_length == 0) {
        /* First edgeset */
        self->last_edgeset.left = left;
        self->last_edgeset.right = right;
        self->last_edgeset.parent = parent;
        self->last_edgeset.children_length = children_length;
        memcpy(self->last_edgeset.children, children,
                children_length * sizeof(node_id_t));
    } else {
        squash = false;
        if (self->last_edgeset.children_length == children_length) {
            squash =
                left == self->last_edgeset.right &&
                parent == self->last_edgeset.parent &&
                memcmp(children, self->last_edgeset.children,
                        children_length * sizeof(node_id_t)) == 0;
        }
        if (squash) {
            self->last_edgeset.right = right;
        } else {
            /* Flush the last edgeset */
            ret = edgeset_table_add_row(self->edgesets,
                    self->last_edgeset.left,
                    self->last_edgeset.right,
                    self->last_edgeset.parent,
                    self->last_edgeset.children,
                    self->last_edgeset.children_length);
            if (ret != 0) {
                goto out;
            }
            self->last_edgeset.left = left;
            self->last_edgeset.right = right;
            self->last_edgeset.parent = parent;
            self->last_edgeset.children_length = children_length;
            memcpy(self->last_edgeset.children, children,
                    children_length * sizeof(node_id_t));
        }
    }
out:
    return ret;
}

int
simplifier_alloc(simplifier_t *self,
        node_table_t *nodes, edgeset_table_t *edgesets, migration_table_t *migrations,
        site_table_t *sites, mutation_table_t *mutations,
        node_id_t *samples, size_t num_samples,
        double sequence_length, int flags)
{
    int ret = 0;
    size_t j;
    node_id_t input_node;

    memset(self, 0, sizeof(simplifier_t));
    self->samples = samples;
    self->num_samples = num_samples;
    self->flags = flags;
    self->sequence_length = sequence_length;
    self->nodes = nodes;
    self->edgesets = edgesets;
    self->sites = sites;
    self->mutations = mutations;

    if (num_samples < 2 || nodes->num_rows == 0 || edgesets->num_rows == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    /* Make a copy of the input nodes and clear the table ready for output */
    ret = node_table_alloc(&self->input_nodes, nodes->num_rows,
            nodes->total_name_length + 1);
    if (ret != 0) {
        goto out;
    }
    ret = node_table_set_columns(&self->input_nodes, nodes->num_rows,
            nodes->flags, nodes->time, nodes->population, nodes->name, nodes->name_length);
    if (ret != 0) {
        goto out;
    }
    ret = node_table_reset(self->nodes);
    if (ret != 0) {
        goto out;
    }
    /* Allocate the heaps used for small objects. */
    /* TODO assuming that the number of edgesets is a good guess here. */
    ret = object_heap_init(&self->segment_heap, sizeof(simplify_segment_t),
            edgesets->num_rows, NULL);
    if (ret != 0) {
        goto out;
    }
    ret = object_heap_init(&self->avl_node_heap, sizeof(avl_node_t),
            edgesets->num_rows, NULL);
    if (ret != 0) {
        goto out;
    }
    ret = object_heap_init(&self->overlap_count_heap, sizeof(overlap_count_t),
            edgesets->num_rows, NULL);
    if (ret != 0) {
        goto out;
    }
    /* Make the maps and set the intial state */
    self->ancestor_map = calloc(self->input_nodes.num_rows,
            sizeof(simplify_segment_t *));
    if (self->ancestor_map == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->nodes->num_rows = 0;
    for (j = 0; j < self->num_samples; j++) {
        input_node = samples[j];
        if (input_node < 0 || input_node >= (node_id_t) self->input_nodes.num_rows) {
            ret = MSP_ERR_OUT_OF_BOUNDS;
            goto out;
        }
        if (!(self->input_nodes.flags[input_node] & MSP_NODE_IS_SAMPLE)) {
            ret = MSP_ERR_BAD_SAMPLES;
            goto out;
        }
        if (self->ancestor_map[samples[j]] != NULL) {
            ret = MSP_ERR_DUPLICATE_SAMPLE;
            goto out;
        }
        self->ancestor_map[samples[j]] = simplifier_alloc_segment(self, 0,
                self->sequence_length, (node_id_t) self->nodes->num_rows, NULL, NULL);
        if (self->ancestor_map[samples[j]] == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        simplifier_record_node(self, input_node);
    }
    avl_init_tree(&self->overlap_counts, cmp_overlap_count, NULL);
    avl_init_tree(&self->merge_queue, cmp_segment_queue, NULL);

    /* Allocate the children and segment buffers */
    self->children_buffer_size = 2;
    self->children_buffer = malloc(self->children_buffer_size * sizeof(node_id_t));
    self->last_edgeset.children = malloc(self->children_buffer_size * sizeof(node_id_t));
    self->segment_buffer_size = 64;
    self->segment_buffer = malloc(self->segment_buffer_size * sizeof(simplify_segment_t *));
    if (self->children_buffer == NULL || self->last_edgeset.children == NULL
            || self->segment_buffer == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    /* Set the initial overlap counts */
    ret = simplifier_insert_overlap_count(self, 0, (uint32_t) num_samples);
    if (ret != 0) {
        goto out;
    }
    ret = simplifier_insert_overlap_count(self, self->sequence_length,
            (uint32_t) num_samples + 1);
    if (ret != 0) {
        goto out;
    }
    /* Set up the sites */

    ret = site_table_reset(self->sites);
    if (ret != 0) {
        goto out;
    }
    ret = mutation_table_reset(self->mutations);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int
simplifier_free(simplifier_t *self)
{
    node_table_free(&self->input_nodes);
    object_heap_free(&self->segment_heap);
    object_heap_free(&self->avl_node_heap);
    object_heap_free(&self->overlap_count_heap);
    msp_safe_free(self->ancestor_map);
    msp_safe_free(self->children_buffer);
    msp_safe_free(self->last_edgeset.children);
    msp_safe_free(self->segment_buffer);
    return 0;
}

static node_id_t *
simplifier_get_children_buffer(simplifier_t *self, size_t num_children)
{
    node_id_t *ret = NULL;
    node_id_t *tmp;

    if (self->children_buffer_size < num_children) {
        self->children_buffer_size = num_children;
        msp_safe_free(self->children_buffer);
        self->children_buffer = malloc(self->children_buffer_size * sizeof(node_id_t));
        if (self->children_buffer == NULL) {
            goto out;
        }
        /* Also realloc the buffer for last children */
        tmp = realloc(self->last_edgeset.children, num_children * sizeof(node_id_t));
        if (tmp == NULL) {
            goto out;
        }
        self->last_edgeset.children = tmp;
    }
    ret = self->children_buffer;
out:
    return ret;
}

static simplify_segment_t **
simplifier_get_segment_buffer(simplifier_t *self, size_t size)
{
    simplify_segment_t **ret = NULL;

    if (self->segment_buffer_size < size) {
        self->segment_buffer_size = size;
        msp_safe_free(self->segment_buffer);
        self->segment_buffer = malloc(self->segment_buffer_size
                * sizeof(simplify_segment_t *));
        if (self->segment_buffer == NULL) {
            goto out;
        }
    }
    ret = self->segment_buffer;
out:
    return ret;
}

static int WARN_UNUSED
simplifier_priority_queue_insert(simplifier_t *self, avl_tree_t *Q, simplify_segment_t *u)
{
    int ret = 0;
    avl_node_t *node;

    assert(u != NULL);
    node = simplifier_alloc_avl_node(self);
    if (node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    avl_init_node(node, u);
    node = avl_insert_node(Q, node);
    assert(node != NULL);
out:
    return ret;
}

static int
simplifier_remove_ancestry(simplifier_t *self, double left, double right,
        node_id_t input_id)
{
    int ret = 0;
    simplify_segment_t *x, *y, *z, *w, *next_w;
    double seg_right, out_left, out_right;
    bool overhang_left, overhang_right;

    x = self->ancestor_map[input_id];
    y = NULL;
    z = NULL;
    w = NULL;
    /* TODO this can be clarified a bit I think
     * 1) Put the y-loop explicitly at the start, so that we just pop off the segments
     *    that are before left.
     * 2) Don't try to reuse x as an output segment. Just free the segment, and
     *    deal with w in a separate block. Performance penalty is negligible and
     *    might make the flow clearer.
     * 3) We're definitely leaking segments somewhere. Put back in the assertion check
     *    for this.
     * 4) Also check for AVL node leakage. We should be using and freeing a lot of these.
     */
    while (x != NULL && right > x->left) {
        if (left < x->right && right > x->left) {
            seg_right = x->right;
            out_left = GSL_MAX(left, x->left);
            out_right = GSL_MIN(right, x->right);
            overhang_left = x->left < out_left;
            overhang_right = x->right > out_right;

            if (overhang_left) {
                y = x;
                x->right = out_left;
                next_w = simplifier_alloc_segment(self, out_left, out_right, x->node, w, NULL);
                if (next_w == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
            } else {
                x->prev = w;
                x->right = out_right;
                next_w = x;
            }
            if (w == NULL) {
                ret = simplifier_priority_queue_insert(self, &self->merge_queue, next_w);
                if (ret != 0) {
                    goto out;
                }
            } else {
                w->next = next_w;
            }
            w = next_w;
            if (overhang_right) {
                z = simplifier_alloc_segment(self, out_right, seg_right, x->node, y, x->next);
                if (z == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
                if (x->next != NULL) {
                    x->next->prev = z;
                }
                break;
            }
        } else {
            y = x;
        }
        x = x->next;
    }
    if (w != NULL) {
        w->next = NULL;

        if (!overhang_right) {
            z = x;
        }
        if (y != NULL) {
            y->next = z;
        }
        if (z != NULL) {
            z->prev = y;
        }
        if (y == NULL) {
            /* TODO I think we might be leaking a segment here when we overwrite this */
            self->ancestor_map[input_id] = z;
        }
    }
out:
    return ret;
}

static int WARN_UNUSED
simplifier_merge_ancestors(simplifier_t *self, node_id_t input_id)
{
    int ret = MSP_ERR_GENERIC;
    bool coalescence = 0;
    /* bool defrag_required = 0; */
    node_id_t v, *children;
    double l, r, r_max, next_l;
    uint32_t j, h;
    avl_node_t *node;
    overlap_count_t *nm, search;
    simplify_segment_t *x, *z, *alpha;
    simplify_segment_t **H = NULL;
    avl_tree_t *Q = &self->merge_queue;

    H = simplifier_get_segment_buffer(self, avl_count(Q));
    if (H == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    /* l_min = 0; */
    z = NULL;
    while (avl_count(Q) > 0) {
        h = 0;
        node = Q->head;
        l = ((simplify_segment_t *) node->item)->left;
        r_max = self->sequence_length;
        while (node != NULL && ((simplify_segment_t *) node->item)->left == l) {
            H[h] = (simplify_segment_t *) node->item;
            r_max = GSL_MIN(r_max, H[h]->right);
            h++;
            simplifier_free_avl_node(self, node);
            avl_unlink_node(Q, node);
            node = node->next;
        }
        next_l = 0;
        if (node != NULL) {
            next_l = ((simplify_segment_t *) node->item)->left;
            r_max = GSL_MIN(r_max, next_l);
        }
        alpha = NULL;
        if (h == 1) {
            x = H[0];
            if (node != NULL && next_l < x->right) {
                alpha = simplifier_alloc_segment(self, x->left, next_l, x->node,
                        NULL, NULL);
                if (alpha == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
                x->left = next_l;
            } else {
                alpha = x;
                x = x->next;
                alpha->next = NULL;
            }
            if (x != NULL) {
                ret = simplifier_priority_queue_insert(self, Q, x);
                if (ret != 0) {
                    goto out;
                }
            }
        } else {
            if (!coalescence) {
                coalescence = 1;
                /* l_min = l; */
                ret = simplifier_record_node(self, input_id);
                if (ret != 0) {
                    goto out;
                }
            }
            v = (node_id_t) self->nodes->num_rows - 1;
            /* Insert overlap counts for bounds, if necessary */
            search.start = l;
            node = avl_search(&self->overlap_counts, &search);
            if (node == NULL) {
                ret = simplifier_copy_overlap_count(self, l);
                if (ret < 0) {
                    goto out;
                }
            }
            search.start = r_max;
            node = avl_search(&self->overlap_counts, &search);
            if (node == NULL) {
                ret = simplifier_copy_overlap_count(self, r_max);
                if (ret < 0) {
                    goto out;
                }
            }
            /* Update the extant segments and allocate alpha if the interval
             * has not coalesced. */
            search.start = l;
            node = avl_search(&self->overlap_counts, &search);
            assert(node != NULL);
            nm = (overlap_count_t *) node->item;
            if (nm->count == h) {
                nm->count = 0;
                node = node->next;
                assert(node != NULL);
                nm = (overlap_count_t *) node->item;
                r = nm->start;
            } else {
                r = l;
                while (nm->count != h && r < r_max) {
                    nm->count -= h - 1;
                    node = node->next;
                    assert(node != NULL);
                    nm = (overlap_count_t *) node->item;
                    r = nm->start;
                }
                alpha = simplifier_alloc_segment(self, l, r, v, NULL, NULL);
                if (alpha == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
            }
            /* Create the record and update the priority queue */
            children = simplifier_get_children_buffer(self, h);
            if (children == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            for (j = 0; j < h; j++) {
                x = H[j];
                children[j] = x->node;
                if (x->right == r) {
                    simplifier_free_segment(self, x);
                    x = x->next;
                } else if (x->right > r) {
                    x->left = r;
                }
                if (x != NULL) {
                    ret = simplifier_priority_queue_insert(self, Q, x);
                    if (ret != 0) {
                        goto out;
                    }
                }
            }
            ret = simplifier_record_edgeset(self, l, r, v, children, h);
            if (ret != 0) {
                goto out;
            }
        }
        /* Loop tail; integrate alpha into the global state */
        if (alpha != NULL) {
            if (z == NULL) {
                self->ancestor_map[input_id] = alpha;
            } else {
                /* defrag_required |= */
                /*     z->right == alpha->left && z->value == alpha->value; */
                z->next = alpha;
                /* fenwick_set_value(&self->links, alpha->id, */
                /*         alpha->right - z->right); */
            }
            alpha->prev = z;
            z = alpha;
        }
    }
    /*
    if (defrag_required) {
        ret = simplifier_defrag_simplify_segment_chain(self, z);
        if (ret != 0) {
            goto out;
        }
    }
    if (coalescence) {
        ret = simplifier_conditional_compress_overlap_counts(self, l_min, r_max);
        if (ret != 0) {
            goto out;
        }

    }
    */
    ret = 0;
out:
    return ret;
}


int WARN_UNUSED
simplifier_run(simplifier_t *self)
{
    int ret = 0;
    size_t j;
    node_id_t parent, current_parent, *children;
    list_len_t k, children_length;
    double left, right;
    size_t num_input_edgesets = self->edgesets->num_rows;
    size_t children_offset = 0;

    /* We modify the edgeset table in place, adding records to the
     * start of the table after we have read in and queued the edgesets
     * for each parent. */
    ret = edgeset_table_reset(self->edgesets);
    if (ret != 0) {
        goto out;
    }

    /* printf("START\n"); */
    /* simplifier_print_state(self); */

    if (num_input_edgesets > 0) {
        current_parent = self->edgesets->parent[0];

        for (j = 0; j < num_input_edgesets; j++) {
            assert(j >= self->edgesets->num_rows);
            parent = self->edgesets->parent[j];
            left = self->edgesets->left[j];
            right = self->edgesets->right[j];
            children_length = self->edgesets->children_length[j];
            children = self->edgesets->children + children_offset;
            children_offset += children_length;

            if (parent != current_parent) {
                /* printf("FLUSH: %d: %d segments\n", current_parent, */
                /*         avl_count(&self->merge_queue)); */
                simplifier_check_state(self);
                /* simplifier_print_state(self); */
                ret = simplifier_merge_ancestors(self, current_parent);
                if (ret != 0) {
                    goto out;
                }
                assert(avl_count(&self->merge_queue) == 0);
                /* printf("Done merge\n"); */
                /* simplifier_check_state(self); */
                if (self->input_nodes.time[current_parent]
                        > self->input_nodes.time[parent]) {
                    ret = MSP_ERR_RECORDS_NOT_TIME_SORTED;
                    goto out;
                }
                current_parent = parent;
            }
            /* printf("EDGESET: %d %d: %f-%f %d: %p\n", (int) j, */
            /*         parent, left, right, children_length, (void *) children); */
            for (k = 0; k < children_length; k++) {
                if (self->ancestor_map[children[k]] != NULL) {
                    /* printf("BEFORE REMOVE: %d: Q size = %d\n", children[k], avl_count(&self->merge_queue)); */
                    /* simplifier_print_state(self); */
                    simplifier_check_state(self);
                    ret = simplifier_remove_ancestry(self, left, right, children[k]);
                    if (ret != 0) {
                        goto out;
                    }
                    /* printf("AFTER REMOVE: %d: Q size = %d\n", children[k], avl_count(&self->merge_queue)); */
                    simplifier_check_state(self);
                    /* printf("DONE AFTER REMOVE CHECK\n"); */
                }

            }

        }
        /* printf("FLUSH: %d\n", current_parent); */
        ret = simplifier_merge_ancestors(self, current_parent);
        if (ret != 0) {
            goto out;
        }
        assert(avl_count(&self->merge_queue) == 0);
        simplifier_check_state(self);
    }
    /* Flush the last edgeset */
    ret = edgeset_table_add_row(self->edgesets,
            self->last_edgeset.left,
            self->last_edgeset.right,
            self->last_edgeset.parent,
            self->last_edgeset.children,
            self->last_edgeset.children_length);
    if (ret != 0) {
        goto out;
    }
    /* printf("DONE\n"); */
    /* simplifier_print_state(self); */
out:

    return ret;
}
