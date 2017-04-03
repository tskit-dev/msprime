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
    fprintf(out, "index\tleft\tright\tnode\tsource\tdest\tpopulation\n");
    for (j = 0; j < self->num_rows; j++) {
        fprintf(out, "%d\t%.3f\t%.3f\t%d\t%d\t%d\t%f\n", (int) j, self->left[j],
                self->right[j], (int) self->node[j], (int) self->source[j],
                (int) self->dest[j], self->time[j]);
    }
}
