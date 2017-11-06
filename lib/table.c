/*
** Copyright (C) 2017 University of Oxford
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
#include "object_heap.h"


#define DEFAULT_SIZE_INCREMENT 1024

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

static int
node_table_expand_main_columns(node_table_t *self, size_t additional_rows)
{
    int ret = 0;
    size_t increment = GSL_MAX(additional_rows, self->max_rows_increment);
    size_t new_size = self->max_rows + increment;

    if ((self->num_rows + additional_rows) > self->max_rows) {
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
node_table_expand_name(node_table_t *self, size_t additional_length)
{
    int ret = 0;
    size_t increment = GSL_MAX(additional_length,
            self->max_total_name_length_increment);
    size_t new_size = self->max_total_name_length + increment;

    if ((self->total_name_length + additional_length) > self->max_total_name_length) {
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
node_table_alloc(node_table_t *self, size_t max_rows_increment,
        size_t max_total_name_length_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(node_table_t));
    if (max_rows_increment == 0) {
       max_rows_increment = DEFAULT_SIZE_INCREMENT;
    }
    if (max_total_name_length_increment == 0) {
        max_total_name_length_increment = DEFAULT_SIZE_INCREMENT;
    }
    self->max_rows_increment = max_rows_increment;
    self->max_total_name_length_increment = max_total_name_length_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    self->max_total_name_length = 0;
    self->total_name_length = 0;
    ret = node_table_expand_main_columns(self, 1);
    if (ret != 0) {
        goto out;
    }
    ret = node_table_expand_name(self, 1);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int
node_table_set_columns(node_table_t *self, size_t num_rows, uint32_t *flags, double *time,
        population_id_t *population, char *name, uint32_t *name_length)
{
    int ret;

    ret = node_table_reset(self);
    if (ret != 0) {
        goto out;
    }
    ret = node_table_append_columns(self, num_rows, flags, time, population, name,
            name_length);
out:
    return ret;
}

int
node_table_append_columns(node_table_t *self, size_t num_rows, uint32_t *flags, double *time,
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
    ret = node_table_expand_main_columns(self, num_rows);
    if (ret != 0) {
        goto out;
    }
    memcpy(self->flags + self->num_rows, flags, num_rows * sizeof(uint32_t));
    memcpy(self->time + self->num_rows, time, num_rows * sizeof(double));
    if (name == NULL) {
        self->total_name_length = 0;
        memset(self->name_length + self->num_rows, 0, num_rows * sizeof(uint32_t));
    } else {
        memcpy(self->name_length + self->num_rows, name_length,
                num_rows * sizeof(uint32_t));
        total_name_length = 0;
        for (j = 0; j < num_rows; j++) {
            total_name_length += name_length[j];
        }
        ret = node_table_expand_name(self, total_name_length);
        if (ret != 0) {
            goto out;
        }
        memcpy(self->name + self->total_name_length, name,
                total_name_length * sizeof(char));
        self->total_name_length += total_name_length;
    }
    if (population == NULL) {
        /* Set population to NULL_POPULATION (-1) if not specified */
        memset(self->population + self->num_rows, 0xff,
                num_rows * sizeof(population_id_t));
    } else {
        memcpy(self->population + self->num_rows, population,
                num_rows * sizeof(population_id_t));
    }
    self->num_rows += num_rows;
out:
    return ret;
}

static int
node_table_add_row_internal(node_table_t *self, uint32_t flags, double time,
        population_id_t population, size_t name_length, const char *name)
{
    assert(self->num_rows < self->max_rows);
    assert(self->total_name_length + name_length < self->max_total_name_length);
    memcpy(self->name + self->total_name_length, name, name_length);
    self->total_name_length += name_length;
    self->flags[self->num_rows] = flags;
    self->time[self->num_rows] = time;
    self->population[self->num_rows] = population;
    self->name_length[self->num_rows] = (uint32_t) name_length;
    self->num_rows++;
    return 0;
}

/* TODO this is a bad API: we should include size_t name_length here. */
int
node_table_add_row(node_table_t *self, uint32_t flags, double time,
        population_id_t population, const char *name)
{
    int ret = 0;
    size_t name_length;

    if (name == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    ret = node_table_expand_main_columns(self, 1);
    if (ret != 0) {
        goto out;
    }
    name_length = strlen(name);
    ret = node_table_expand_name(self, name_length);
    if (ret != 0) {
        goto out;
    }
    ret = node_table_add_row_internal(self, flags, time, population, name_length, name);
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
 * edge table
 *************************/

static int
edge_table_expand_columns(edge_table_t *self, size_t additional_rows)
{
    int ret = 0;
    size_t increment = GSL_MAX(additional_rows, self->max_rows_increment);
    size_t new_size = self->max_rows + increment;

    if ((self->num_rows + additional_rows) > self->max_rows) {
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
        ret = expand_column((void **) &self->child, new_size, sizeof(node_id_t));
        if (ret != 0) {
            goto out;
        }
        self->max_rows = new_size;
    }
out:
    return ret;
}

int
edge_table_alloc(edge_table_t *self, size_t max_rows_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(edge_table_t));
    if (max_rows_increment == 0) {
        max_rows_increment = DEFAULT_SIZE_INCREMENT;
    }
    self->max_rows_increment = max_rows_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    ret = edge_table_expand_columns(self, 1);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int
edge_table_add_row(edge_table_t *self, double left, double right, node_id_t parent,
        node_id_t child)
{
    int ret = 0;

    ret = edge_table_expand_columns(self, 1);
    if (ret != 0) {
        goto out;
    }
    self->left[self->num_rows] = left;
    self->right[self->num_rows] = right;
    self->parent[self->num_rows] = parent;
    self->child[self->num_rows] = child;
    self->num_rows++;
out:
    return ret;
}

int
edge_table_set_columns(edge_table_t *self,
        size_t num_rows, double *left, double *right, node_id_t *parent, node_id_t *child)
{
    int ret = 0;

    ret = edge_table_reset(self);
    if (ret != 0) {
        goto out;
    }
    ret = edge_table_append_columns(self, num_rows, left, right, parent, child);
out:
    return ret;
}

int
edge_table_append_columns(edge_table_t *self,
        size_t num_rows, double *left, double *right, node_id_t *parent, node_id_t *child)
{
    int ret;

    if (left == NULL || right == NULL || parent == NULL || child == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    ret = edge_table_expand_columns(self, num_rows);
    if (ret != 0) {
        goto out;
    }
    memcpy(self->left + self->num_rows, left, num_rows * sizeof(double));
    memcpy(self->right + self->num_rows, right, num_rows * sizeof(double));
    memcpy(self->parent + self->num_rows, parent, num_rows * sizeof(node_id_t));
    memcpy(self->child + self->num_rows, child, num_rows * sizeof(node_id_t));
    self->num_rows += num_rows;
out:
    return ret;
}

int
edge_table_reset(edge_table_t *self)
{
    self->num_rows = 0;
    return 0;
}

int
edge_table_free(edge_table_t *self)
{
    msp_safe_free(self->left);
    msp_safe_free(self->right);
    msp_safe_free(self->parent);
    msp_safe_free(self->child);
    return 0;
}

void
edge_table_print_state(edge_table_t *self, FILE *out)
{
    size_t j;

    fprintf(out, TABLE_SEP);
    fprintf(out, "edge_table: %p:\n", (void *) self);
    fprintf(out, "num_rows          = %d\tmax= %d\tincrement = %d)\n",
            (int) self->num_rows, (int) self->max_rows, (int) self->max_rows_increment);
    fprintf(out, TABLE_SEP);
    fprintf(out, "index\tleft\tright\tparent\tchild\n");
    for (j = 0; j < self->num_rows; j++) {
        fprintf(out, "%d\t%.3f\t%.3f\t%d\t%d\t", (int) j, self->left[j], self->right[j],
                self->parent[j], self->child[j]);
        fprintf(out, "\n");
    }
}

/*************************
 * site table
 *************************/

static int
site_table_expand_main_columns(site_table_t *self, size_t additional_rows)
{
    int ret = 0;
    size_t increment = GSL_MAX(additional_rows, self->max_rows_increment);
    size_t new_size = self->max_rows + increment;

    if ((self->num_rows + additional_rows) > self->max_rows) {
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
site_table_expand_ancestral_state(site_table_t *self, size_t additional_length)
{
    int ret = 0;
    size_t increment = GSL_MAX(additional_length,
            self->max_total_ancestral_state_length_increment);
    size_t new_size = self->max_total_ancestral_state_length + increment;

    if ((self->total_ancestral_state_length + additional_length)
            > self->max_total_ancestral_state_length) {
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
site_table_alloc(site_table_t *self, size_t max_rows_increment,
        size_t max_total_ancestral_state_length_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(site_table_t));
    if (max_rows_increment == 0) {
        max_rows_increment = DEFAULT_SIZE_INCREMENT;
    }
    if (max_total_ancestral_state_length_increment == 0) {
        max_total_ancestral_state_length_increment = DEFAULT_SIZE_INCREMENT;
    }
    self->max_rows_increment = max_rows_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    self->max_total_ancestral_state_length_increment =
        max_total_ancestral_state_length_increment;
    self->max_total_ancestral_state_length = 0;
    self->total_ancestral_state_length = 0;
    ret = site_table_expand_main_columns(self, 1);
    if (ret != 0) {
        goto out;
    }
    ret = site_table_expand_ancestral_state(self, 1);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int
site_table_add_row(site_table_t *self, double position, const char *ancestral_state,
        list_len_t ancestral_state_length)
{
    int ret = 0;

    ret = site_table_expand_main_columns(self, 1);
    if (ret != 0) {
        goto out;
    }
    ret = site_table_expand_ancestral_state(self, ancestral_state_length);
    if (ret != 0) {
        goto out;
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
site_table_append_columns(site_table_t *self, size_t num_rows, double *position,
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
    memcpy(self->position + self->num_rows, position, num_rows * sizeof(double));
    memcpy(self->ancestral_state + self->total_ancestral_state_length, ancestral_state,
            total_ancestral_state_length * sizeof(char));
    memcpy(self->ancestral_state_length + self->num_rows, ancestral_state_length,
            num_rows * sizeof(uint32_t));
    self->num_rows += num_rows;
    self->total_ancestral_state_length += total_ancestral_state_length;
out:
    return ret;
}

int
site_table_set_columns(site_table_t *self, size_t num_rows, double *position,
        const char *ancestral_state, list_len_t *ancestral_state_length)
{
    int ret = 0;

    ret = site_table_reset(self);
    if (ret != 0) {
        goto out;
    }
    ret = site_table_append_columns(self, num_rows, position, ancestral_state,
            ancestral_state_length);
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

static int
mutation_table_expand_main_columns(mutation_table_t *self, size_t additional_rows)
{
    int ret = 0;
    size_t increment = GSL_MAX(additional_rows, self->max_rows_increment);
    size_t new_size = self->max_rows + increment;

    if ((self->num_rows + additional_rows) > self->max_rows) {
        ret = expand_column((void **) &self->site, new_size, sizeof(site_id_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->node, new_size, sizeof(node_id_t));
        if (ret != 0) {
            goto out;
        }
        ret = expand_column((void **) &self->parent, new_size, sizeof(mutation_id_t));
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
mutation_table_expand_derived_state(mutation_table_t *self, size_t additional_length)
{
    int ret = 0;
    size_t increment = GSL_MAX(additional_length,
            self->max_total_derived_state_length_increment);
    size_t new_size = self->max_total_derived_state_length + increment;

    if ((self->total_derived_state_length + additional_length)
            > self->max_total_derived_state_length) {
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
mutation_table_alloc(mutation_table_t *self, size_t max_rows_increment,
        size_t max_total_derived_state_length_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(mutation_table_t));
    if (max_rows_increment == 0) {
        max_rows_increment = DEFAULT_SIZE_INCREMENT;
    }
    if (max_total_derived_state_length_increment == 0) {
        max_total_derived_state_length_increment = DEFAULT_SIZE_INCREMENT;
    }
    self->max_rows_increment = max_rows_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    self->max_total_derived_state_length_increment =
        max_total_derived_state_length_increment;
    self->max_total_derived_state_length = 0;
    self->total_derived_state_length = 0;
    ret = mutation_table_expand_main_columns(self, 1);
    if (ret != 0) {
        goto out;
    }
    ret = mutation_table_expand_derived_state(self, 1);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int
mutation_table_add_row(mutation_table_t *self, site_id_t site, node_id_t node,
        mutation_id_t parent, const char *derived_state, list_len_t derived_state_length)
{
    int ret = 0;

    ret = mutation_table_expand_main_columns(self, 1);
    if (ret != 0) {
        goto out;
    }
    ret = mutation_table_expand_derived_state(self, derived_state_length);
    if (ret != 0) {
        goto out;
    }
    self->site[self->num_rows] = site;
    self->node[self->num_rows] = node;
    self->parent[self->num_rows] = parent;
    self->derived_state_length[self->num_rows] = (list_len_t) derived_state_length;
    memcpy(self->derived_state + self->total_derived_state_length, derived_state,
            derived_state_length * sizeof(char));
    self->total_derived_state_length += derived_state_length;
    self->num_rows++;
out:
    return ret;
}

int
mutation_table_append_columns(mutation_table_t *self, size_t num_rows, site_id_t *site,
        node_id_t *node, mutation_id_t *parent,
        const char *derived_state, uint32_t *derived_state_length)
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
    memcpy(self->site + self->num_rows, site, num_rows * sizeof(site_id_t));
    memcpy(self->node + self->num_rows, node, num_rows * sizeof(node_id_t));
    memcpy(self->derived_state_length + self->num_rows, derived_state_length,
            num_rows * sizeof(node_id_t));
    memcpy(self->derived_state + self->total_derived_state_length, derived_state,
            total_derived_state_length * sizeof(char));
    if (parent == NULL) {
        /* If parent is NULL, set all parents to the null mutation */
        memset(self->parent + self->num_rows, 0xff, num_rows * sizeof(mutation_id_t));
    } else {
        memcpy(self->parent + self->num_rows, parent, num_rows * sizeof(mutation_id_t));
    }
    self->num_rows += num_rows;
    self->total_derived_state_length += total_derived_state_length;
out:
    return ret;
}

int
mutation_table_set_columns(mutation_table_t *self, size_t num_rows, site_id_t *site,
        node_id_t *node, mutation_id_t *parent,
        const char *derived_state, uint32_t *derived_state_length)
{
    int ret = 0;

    ret = mutation_table_reset(self);
    if (ret != 0) {
        goto out;
    }
    ret = mutation_table_append_columns(self, num_rows ,site, node, parent,
            derived_state, derived_state_length);
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
            && memcmp(self->parent, other->parent,
                    self->num_rows * sizeof(mutation_id_t)) == 0
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
    msp_safe_free(self->parent);
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
    fprintf(out, "index\tsite\tnode\tparent\tderived_state_length\tderived_state\n");
    offset = 0;
    for (j = 0; j < self->num_rows; j++) {
        fprintf(out, "%d\t%d\t%d\t%d\t%d\t", (int) j, self->site[j], self->node[j],
                self->parent[j], self->derived_state_length[j]);
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

static int
migration_table_expand(migration_table_t *self, size_t additional_rows)
{
    int ret = 0;
    size_t increment = GSL_MAX(additional_rows, self->max_rows_increment);
    size_t new_size = self->max_rows + increment;

    if ((self->num_rows + additional_rows) > self->max_rows) {
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
migration_table_alloc(migration_table_t *self, size_t max_rows_increment)
{
    int ret = 0;

    memset(self, 0, sizeof(migration_table_t));
    if (max_rows_increment == 0) {
        max_rows_increment = DEFAULT_SIZE_INCREMENT;
    }
    self->max_rows_increment = max_rows_increment;
    self->max_rows = 0;
    self->num_rows = 0;
    ret = migration_table_expand(self, 1);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int
migration_table_append_columns(migration_table_t *self, size_t num_rows, double *left,
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
    memcpy(self->left + self->num_rows, left, num_rows * sizeof(double));
    memcpy(self->right + self->num_rows, right, num_rows * sizeof(double));
    memcpy(self->node + self->num_rows, node, num_rows * sizeof(node_id_t));
    memcpy(self->source + self->num_rows, source, num_rows * sizeof(population_id_t));
    memcpy(self->dest + self->num_rows, dest, num_rows * sizeof(population_id_t));
    memcpy(self->time + self->num_rows, time, num_rows * sizeof(double));
    self->num_rows += num_rows;
out:
    return ret;
}

int
migration_table_set_columns(migration_table_t *self, size_t num_rows, double *left,
        double *right, node_id_t *node, population_id_t *source, population_id_t *dest,
        double *time)
{
    int ret;

    ret = migration_table_reset(self);
    if (ret != 0) {
        goto out;
    }
    ret = migration_table_append_columns(self, num_rows, left, right, node, source,
            dest, time);
out:
    return ret;
}

int
migration_table_add_row(migration_table_t *self, double left, double right,
        node_id_t node, population_id_t source, population_id_t dest, double time)
{
    int ret = 0;

    ret = migration_table_expand(self, 1);
    if (ret != 0) {
        goto out;
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
 * sort_tables
 *************************/

typedef struct {
    double left;
    double right;
    node_id_t parent;
    node_id_t child;
    double time;
} edge_sort_t;

typedef struct {
    /* Input tables. */
    node_table_t *nodes;
    edge_table_t *edges;
    site_table_t *sites;
    mutation_table_t *mutations;
    migration_table_t *migrations;
    /* sorting edges */
    edge_sort_t *sorted_edges;
    /* sorting sites */
    site_id_t *site_id_map;
    site_t *sorted_sites;
    mutation_t *sorted_mutations;
    char *ancestral_state_mem;
    char *derived_state_mem;
} table_sorter_t;

static int
cmp_site(const void *a, const void *b) {
    const site_t *ia = (const site_t *) a;
    const site_t *ib = (const site_t *) b;
    /* Compare sites by position */
    return (ia->position > ib->position) - (ia->position < ib->position);
}

static int
cmp_mutation(const void *a, const void *b) {
    const mutation_t *ia = (const mutation_t *) a;
    const mutation_t *ib = (const mutation_t *) b;
    /* Compare mutations by site */
    int ret = (ia->site > ib->site) - (ia->site < ib->site);
    if (ret == 0) {
        /* Within a particular site sort by ID. This ensures that relative ordering
         * within a site is maintained */
        ret = (ia->id > ib->id) - (ia->id < ib->id);
    }
    return ret;
}

static int
cmp_edge(const void *a, const void *b) {
    const edge_sort_t *ca = (const edge_sort_t *) a;
    const edge_sort_t *cb = (const edge_sort_t *) b;

    int ret = (ca->time > cb->time) - (ca->time < cb->time);
    /* If time values are equal, sort by the parent node */
    if (ret == 0) {
        ret = (ca->parent > cb->parent) - (ca->parent < cb->parent);
        /* If the parent nodes are equal, sort by the child ID. */
        if (ret == 0) {
            ret = (ca->child > cb->child) - (ca->child < cb->child);
            /* If the child nodes are equal, sort by the left coordinate. */
            if (ret == 0) {
                ret = (ca->left > cb->left) - (ca->left < cb->left);
            }
        }
    }
    return ret;
}

static int
table_sorter_alloc(table_sorter_t *self, node_table_t *nodes, edge_table_t *edges,
        site_table_t *sites, mutation_table_t *mutations, migration_table_t *migrations)
{
    int ret = 0;

    memset(self, 0, sizeof(table_sorter_t));
    if (nodes == NULL || edges == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->nodes = nodes;
    self->edges = edges;
    self->mutations = mutations;
    self->sites = sites;
    self->migrations = migrations;
    self->sorted_edges = malloc(edges->num_rows * sizeof(edge_sort_t));
    if (self->sorted_edges == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    if (self->sites != NULL) {
        /* If you provide a site table, you must provide a mutation table (even if it is
         * empty */
        if (self->mutations == NULL) {
            ret = MSP_ERR_BAD_PARAM_VALUE;
            goto out;
        }
        self->sorted_sites = malloc(sites->num_rows * sizeof(site_t));
        self->ancestral_state_mem = malloc(sites->total_ancestral_state_length * sizeof(char));
        self->site_id_map = malloc(sites->num_rows * sizeof(site_id_t));
        if (self->sorted_sites == NULL || self->ancestral_state_mem == NULL
                || self->site_id_map == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        self->sorted_mutations = malloc(mutations->num_rows * sizeof(mutation_t));
        self->derived_state_mem = malloc(mutations->total_derived_state_length * sizeof(char));
        if (self->sorted_mutations == NULL || self->derived_state_mem == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
out:
    return ret;
}

static int
table_sorter_sort_edges(table_sorter_t *self, size_t start)
{
    int ret = 0;
    edge_sort_t *e;
    size_t j, k;
    size_t n = self->edges->num_rows - start;

    for (j = 0; j < n; j++) {
        e = self->sorted_edges + j;
        k = start + j;
        e->left = self->edges->left[k];
        e->right = self->edges->right[k];
        e->parent = self->edges->parent[k];
        e->child = self->edges->child[k];
        if (e->parent >= (node_id_t) self->nodes->num_rows) {
            ret = MSP_ERR_OUT_OF_BOUNDS;
            goto out;
        }
        e->time = self->nodes->time[e->parent];
    }
    qsort(self->sorted_edges, n, sizeof(edge_sort_t), cmp_edge);
    /* Copy the edges back into the table. */
    for (j = 0; j < n; j++) {
        e = self->sorted_edges + j;
        k = start + j;
        self->edges->left[k] = e->left;
        self->edges->right[k] = e->right;
        self->edges->parent[k] = e->parent;
        self->edges->child[k] = e->child;
    }
out:
    return ret;
}

static int
table_sorter_sort_sites(table_sorter_t *self)
{
    int ret = 0;
    size_t j, ancestral_state_offset;

    memcpy(self->ancestral_state_mem, self->sites->ancestral_state,
            self->sites->total_ancestral_state_length * sizeof(char));
    ancestral_state_offset = 0;
    for (j = 0; j < self->sites->num_rows; j++) {
        self->sorted_sites[j].id = (site_id_t) j;
        self->sorted_sites[j].position = self->sites->position[j];
        self->sorted_sites[j].ancestral_state_length = self->sites->ancestral_state_length[j];
        self->sorted_sites[j].ancestral_state = self->ancestral_state_mem
            + ancestral_state_offset;
        ancestral_state_offset += self->sites->ancestral_state_length[j];
    }
    /* Sort the sites by position */
    qsort(self->sorted_sites, self->sites->num_rows, sizeof(site_t), cmp_site);
    /* Build the mapping from old site IDs to new site IDs and copy back into the table */
    ancestral_state_offset = 0;
    for (j = 0; j < self->sites->num_rows; j++) {
        self->site_id_map[self->sorted_sites[j].id] = (site_id_t) j;
        self->sites->position[j] = self->sorted_sites[j].position;
        self->sites->ancestral_state_length[j] = self->sorted_sites[j].ancestral_state_length;
        memcpy(self->sites->ancestral_state + ancestral_state_offset,
            self->sorted_sites[j].ancestral_state, self->sorted_sites[j].ancestral_state_length);
        ancestral_state_offset += self->sorted_sites[j].ancestral_state_length;
    }
    return ret;
}

static int
table_sorter_sort_mutations(table_sorter_t *self)
{
    int ret = 0;
    size_t j, derived_state_offset;
    site_id_t site;
    node_id_t node;
    mutation_id_t parent, mapped_parent;
    mutation_id_t *mutation_id_map = malloc(self->mutations->num_rows
            * sizeof(mutation_id_t));

    memcpy(self->derived_state_mem, self->mutations->derived_state,
            self->mutations->total_derived_state_length * sizeof(char));
    derived_state_offset = 0;
    for (j = 0; j < self->mutations->num_rows; j++) {
        site = self->mutations->site[j];
        if (site >= (site_id_t) self->sites->num_rows) {
            ret = MSP_ERR_OUT_OF_BOUNDS;
            goto out;
        }
        node = self->mutations->node[j];
        if (node >= (node_id_t) self->nodes->num_rows) {
            ret = MSP_ERR_OUT_OF_BOUNDS;
            goto out;
        }
        parent = self->mutations->parent[j];
        if (parent != MSP_NULL_MUTATION) {
            if (parent < 0 || parent >= (mutation_id_t) self->mutations->num_rows) {
                ret = MSP_ERR_MUTATION_OUT_OF_BOUNDS;
                goto out;
            }
        }
        self->sorted_mutations[j].id = (mutation_id_t) j;
        self->sorted_mutations[j].site = self->site_id_map[site];
        self->sorted_mutations[j].node = node;
        self->sorted_mutations[j].parent = self->mutations->parent[j];
        self->sorted_mutations[j].derived_state_length =
            self->mutations->derived_state_length[j];
        self->sorted_mutations[j].derived_state = self->derived_state_mem
            + derived_state_offset;
        derived_state_offset += self->mutations->derived_state_length[j];
    }

    qsort(self->sorted_mutations, self->mutations->num_rows, sizeof(mutation_t),
        cmp_mutation);
    /* Make a first pass through the sorted mutations to build the ID map. */
    for (j = 0; j < self->mutations->num_rows; j++) {
        mutation_id_map[self->sorted_mutations[j].id] = (mutation_id_t) j;
    }

    /* Copy the sorted mutations back into the table */
    derived_state_offset = 0;
    for (j = 0; j < self->mutations->num_rows; j++) {
        self->mutations->site[j] = self->sorted_mutations[j].site;
        self->mutations->node[j] = self->sorted_mutations[j].node;
        mapped_parent = MSP_NULL_MUTATION;
        parent = self->sorted_mutations[j].parent;
        if (parent != MSP_NULL_MUTATION) {
            mapped_parent = mutation_id_map[parent];
        }
        self->mutations->parent[j] = mapped_parent;
        self->mutations->derived_state_length[j] =
            self->sorted_mutations[j].derived_state_length;
        memcpy(self->mutations->derived_state + derived_state_offset,
            self->sorted_mutations[j].derived_state,
            self->sorted_mutations[j].derived_state_length * sizeof(char));
        derived_state_offset += self->sorted_mutations[j].derived_state_length;
    }
out:
    msp_safe_free(mutation_id_map);
    return ret;
}

static int
table_sorter_run(table_sorter_t *self, size_t edge_start)
{
    int ret = 0;

    if (edge_start > self->edges->num_rows) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    ret = table_sorter_sort_edges(self, edge_start);
    if (ret != 0) {
        goto out;
    }
    if (self->sites != NULL) {
        ret = table_sorter_sort_sites(self);
        if (ret != 0) {
            goto out;
        }
        ret = table_sorter_sort_mutations(self);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

static void
table_sorter_free(table_sorter_t *self)
{
    msp_safe_free(self->sorted_edges);
    msp_safe_free(self->sorted_sites);
    msp_safe_free(self->site_id_map);
    msp_safe_free(self->sorted_mutations);
    msp_safe_free(self->ancestral_state_mem);
    msp_safe_free(self->derived_state_mem);
}

int
sort_tables(node_table_t *nodes, edge_table_t *edges, migration_table_t *migrations,
        site_table_t *sites, mutation_table_t *mutations, size_t edge_start)
{
    int ret = 0;
    table_sorter_t *sorter = NULL;

    sorter = malloc(sizeof(table_sorter_t));
    if (sorter == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = table_sorter_alloc(sorter, nodes, edges, sites, mutations, migrations);
    if (ret != 0) {
        goto out;
    }
    ret = table_sorter_run(sorter, edge_start);
out:
    if (sorter != NULL) {
        table_sorter_free(sorter);
        free(sorter);
    }
    return ret;
}

/*************************
 * simplifier
 *************************/

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

static int
cmp_mutation_position_map(const void *a, const void *b) {
    const mutation_position_map_t *ia = (const mutation_position_map_t *) a;
    const mutation_position_map_t *ib = (const mutation_position_map_t *) b;
    int ret = (ia->position > ib->position) - (ia->position < ib->position);
    return ret;
}

static void
simplifier_check_state(simplifier_t *self)
{
    size_t j;
    size_t total_segments = 0;
    size_t total_avl_nodes = 0;
    site_id_t site;
    avl_node_t *a;
    simplify_segment_t *u;
    mutation_node_list_t *mnl;
    mutation_position_map_t *mpm;

    for (j = 0; j < self->input_nodes.num_rows; j++) {
        for (u = self->ancestor_map[j]; u != NULL; u = u->next) {
            assert(u->left < u->right);
            if (u->next != NULL) {
                assert(u->right <= u->next->left);
            }
            total_segments++;
        }
        for (a = self->mutation_position_map[j].head; a != NULL; a = a->next) {
            mpm = (mutation_position_map_t *) a->item;
            assert(mpm->head != NULL);
            for (mnl = mpm->head; mnl != NULL; mnl = mnl->next) {
                site = self->input_mutations.site[mnl->mutation_id];
                assert(self->input_mutations.node[mnl->mutation_id] == (node_id_t) j);
                assert(self->input_sites.position[site] == mpm->position);
            }
        }
        total_avl_nodes += avl_count(&self->mutation_position_map[j]);
    }
    for (a = self->merge_queue.head; a != NULL; a = a->next) {
        total_avl_nodes++;
        for (u = (simplify_segment_t *) a->item; u != NULL; u = u->next) {
            assert(u->left < u->right);
            if (u->next != NULL) {
                assert(u->right <= u->next->left);
            }
            total_segments++;
        }
    }
    assert(total_segments == object_heap_get_num_allocated(&self->segment_heap));
    assert(total_avl_nodes == object_heap_get_num_allocated(&self->avl_node_heap));
}

static void
print_segment_chain(simplify_segment_t *head, FILE *out)
{
    simplify_segment_t *u;

    for (u = head; u != NULL; u = u->next) {
        fprintf(out, "(%f,%f->%d)", u->left, u->right, u->node);
    }
}

void
simplifier_print_state(simplifier_t *self, FILE *out)
{
    size_t j;
    avl_node_t *avl_node;
    simplify_segment_t *u;
    mutation_node_list_t *mnl;
    mutation_position_map_t *mpm;

    fprintf(out, "--simplifier state--\n");
    fprintf(out, "===\nInput nodes\n==\n");
    node_table_print_state(&self->input_nodes, out);
    fprintf(out, "===\nInput edges\n==\n");
    edge_table_print_state(&self->input_edges, out);
    fprintf(out, "===\nInput sites\n==\n");
    site_table_print_state(&self->input_sites, out);
    fprintf(out, "===\nInput mutations\n==\n");
    mutation_table_print_state(&self->input_mutations, out);
    fprintf(out, "===\nOutput tables\n==\n");
    node_table_print_state(self->nodes, out);
    edge_table_print_state(self->edges, out);
    site_table_print_state(self->sites, out);
    mutation_table_print_state(self->mutations, out);
    fprintf(out, "===\nmemory heaps\n==\n");
    fprintf(out, "segment_heap:\n");
    object_heap_print_state(&self->segment_heap, out);
    fprintf(out, "avl_node_heap:\n");
    object_heap_print_state(&self->avl_node_heap, out);
    fprintf(out, "===\nancestors\n==\n");
    for (j = 0; j < self->input_nodes.num_rows; j++) {
        if (self->ancestor_map[j] != NULL) {
            fprintf(out, "%d:\t", (int) j);
            print_segment_chain(self->ancestor_map[j], out);
            fprintf(out, "\n");
        }
    }
    fprintf(out, "===\nnode_id map (input->output)\n==\n");
    for (j = 0; j < self->input_nodes.num_rows; j++) {
        if (self->node_id_map[j] != MSP_NULL_NODE) {
            fprintf(out, "%d->%d\n", (int) j, self->node_id_map[j]);
        }
    }
    fprintf(out, "===\nmerge queue\n==\n");
    for (avl_node = self->merge_queue.head; avl_node != NULL; avl_node = avl_node->next) {
        u = (simplify_segment_t *) avl_node->item;
        print_segment_chain(u, out);
        fprintf(out, "\n");
    }
    fprintf(out, "===\nbuffered edges\n==\n");
    for (j = 0; j < self->num_buffered_edges; j++) {
        fprintf(out, "%f\t%f\t%d\t%d\n", self->edge_buffer[j].left,
                self->edge_buffer[j].right, self->edge_buffer[j].parent,
                self->edge_buffer[j].child);
    }
    fprintf(out, "===\nmutation position map\n==\n");
    for (j = 0; j < self->input_nodes.num_rows; j++) {
        if (avl_count(&self->mutation_position_map[j]) > 0) {
            fprintf(out, "node %d:\n", (int) j);
            for (avl_node = self->mutation_position_map[j].head;
                    avl_node != NULL; avl_node = avl_node->next) {
                mpm = (mutation_position_map_t *) avl_node->item;
                fprintf(out, "\t%f -> ", mpm->position);
                for (mnl = mpm->head; mnl != NULL; mnl = mnl->next) {
                    fprintf(out, "%d, ", mnl->mutation_id);
                }
                fprintf(out, "\n");
            }
            fprintf(out, "\n");
        }
    }
    fprintf(out, "===\nmutation node map\n==\n");
    for (j = 0; j < self->input_mutations.num_rows; j++) {
        fprintf(out, "%d\t-> %d\n", (int) j, self->mutation_node_map[j]);
    }
    simplifier_check_state(self);
}

static simplify_segment_t * WARN_UNUSED
simplifier_alloc_segment(simplifier_t *self, double left, double right, node_id_t node,
        simplify_segment_t *next)
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

/* Add a new node to the output node table corresponding to the specified input id */
static int WARN_UNUSED
simplifier_record_node(simplifier_t *self, node_id_t input_id, bool is_sample)
{
    int ret = 0;
    const char *name = self->input_nodes.name + self->node_name_offset[input_id];
    uint32_t flags = self->input_nodes.flags[input_id];

    /* Zero out the sample bit */
    flags &= (uint32_t) ~MSP_NODE_IS_SAMPLE;
    if (is_sample) {
        flags |= MSP_NODE_IS_SAMPLE;
    }
    self->node_id_map[input_id] = (node_id_t) self->nodes->num_rows;
    ret = node_table_add_row_internal(self->nodes, flags,
            self->input_nodes.time[input_id], self->input_nodes.population[input_id],
            self->input_nodes.name_length[input_id], name);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

static int
simplifier_flush_edges(simplifier_t *self)
{
    int ret = 0;
    size_t j, num_output_edges;
    edge_t e;

    if (self->num_buffered_edges > 0) {
        ret = squash_edges(self->edge_buffer, self->num_buffered_edges, &num_output_edges);
        if (ret != 0) {
            goto out;
        }
        /* Flush these edges to the table */
        for (j = 0; j < num_output_edges; j++) {
            e = self->edge_buffer[j];
            ret = edge_table_add_row(self->edges, e.left, e.right, e.parent, e.child);
            if (ret != 0) {
                goto out;
            }
        }
    }
    self->num_buffered_edges = 0;
out:
    return ret;
}

/* Records the specified edge in the output table */
static int
simplifier_record_edge(simplifier_t *self, double left, double right, node_id_t parent,
        node_id_t child)
{
    int ret = 0;
    edge_t *e;

    if (self->num_buffered_edges == self->max_buffered_edges - 1) {
        /* Grow the array. Use a doubling strategy here as we expect this
         * to usually be reasonably small */
        self->max_buffered_edges *= 2;
        e = realloc(self->edge_buffer, self->max_buffered_edges * sizeof(edge_t));
        if (e == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        self->edge_buffer = e;
    }
    assert(self->num_buffered_edges < self->max_buffered_edges);
    e = self->edge_buffer + self->num_buffered_edges;
    e->left = left;
    e->right = right;
    e->parent = parent;
    e->child = child;
    self->num_buffered_edges++;
out:
    return ret;
}

static int
simplifier_init_sites(simplifier_t *self)
{
    int ret = 0;
    size_t j, next_mutation_position_map;
    mutation_position_map_t *mpm, search;
    mutation_node_list_t *mnl;
    node_id_t node;
    site_id_t site;
    avl_node_t *avl_node;
    avl_tree_t *avl_tree;

    self->mutation_id_map = calloc(self->input_mutations.num_rows, sizeof(mutation_id_t));
    self->mutation_node_map = calloc(self->input_mutations.num_rows, sizeof(node_id_t));
    self->mutation_position_map = calloc(self->input_nodes.num_rows, sizeof(avl_tree_t));
    self->mutation_node_list_mem = calloc(self->input_mutations.num_rows,
            sizeof(mutation_node_list_t));
    self->mutation_position_map_mem = calloc(self->input_mutations.num_rows,
            sizeof(mutation_position_map_t));
    if (self->mutation_id_map == NULL || self->mutation_node_map == NULL
            || self->mutation_position_map == NULL
            || self->mutation_node_list_mem == NULL
            || self->mutation_position_map_mem == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memset(self->mutation_id_map, 0xff,
            self->input_mutations.num_rows * sizeof(mutation_id_t));
    memset(self->mutation_node_map, 0xff,
            self->input_mutations.num_rows * sizeof(node_id_t));

    /* The mutation_position_map stores an AVL tree for each input node. Each of
     * these AVL trees maps the sites for mutations carrying the node in question
     * to their positions and list of mutations IDs */
    /* TODO: we should simplify this and remove the AVL trees. See the note above
     * in simplifier_record_mutations */
    for (j = 0; j < self->input_nodes.num_rows; j++) {
        avl_init_tree(self->mutation_position_map + j, cmp_mutation_position_map, NULL);
    }
    next_mutation_position_map = 0;
    for (j = 0; j < self->input_mutations.num_rows; j++) {
        site = self->input_mutations.site[j];
        node = self->input_mutations.node[j];

        avl_tree = self->mutation_position_map + node;
        /* If this position has already been inserted into the avl tree for this
         * node, add this mutation_id to its list. Otherwise alloc and insert a
         * new position into the list */
        search.position = self->input_sites.position[site];

        avl_node = avl_search(avl_tree, &search);
        if (avl_node == NULL) {
            assert(next_mutation_position_map < self->input_mutations.num_rows);
            mpm = self->mutation_position_map_mem + next_mutation_position_map;
            mpm->position = search.position;
            next_mutation_position_map++;
            avl_node = simplifier_alloc_avl_node(self);
            if (avl_node == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            avl_init_node(avl_node, mpm);
            avl_node = avl_insert_node(avl_tree, avl_node);
            assert(avl_node != NULL);
        } else {
            mpm = (mutation_position_map_t *) avl_node->item;
        }
        /* Insert this mutation at the head of the list for this position */
        mnl = self->mutation_node_list_mem + j;
        mnl->mutation_id = (mutation_id_t) j;
        mnl->next = mpm->head;
        mpm->head = mnl;
    }
out:
    return ret;
}

static int
simplifier_check_input(simplifier_t *self)
{
    int ret = MSP_ERR_GENERIC;
    node_id_t num_nodes = (node_id_t) self->nodes->num_rows;
    site_id_t num_sites;
    double *time = self->nodes->time;
    char *node_seen = NULL;
    node_id_t last_parent, parent, child;
    size_t j;

    node_seen = calloc((size_t) num_nodes, sizeof(char));
    /* Check the edges */
    last_parent = self->edges->parent[0];
    for (j = 0; j < self->edges->num_rows; j++) {
        if (self->edges->left[j] >= self->edges->right[j]) {
            ret = MSP_ERR_BAD_EDGE_INTERVAL;
            goto out;
        }
        parent = self->edges->parent[j];
        if (parent < 0 || parent >= num_nodes) {
            ret = MSP_ERR_NODE_OUT_OF_BOUNDS;
            goto out;
        }
        if (parent != last_parent) {
            node_seen[last_parent] = 1;
            if (node_seen[parent] != 0) {
                ret = MSP_ERR_EDGES_NOT_SORTED_PARENT_TIME;
                goto out;
            }
            if (time[last_parent] > time[parent]) {
                ret = MSP_ERR_EDGES_NOT_SORTED_PARENT_TIME;
                goto out;
            }
            last_parent = parent;
        }
        /* Check the children */
        child = self->edges->child[j];
        if (child < 0 || child >= num_nodes) {
            ret = MSP_ERR_NODE_OUT_OF_BOUNDS;
            goto out;
        }
        /* time[child] must be < time[parent] */
        if (time[child] >= time[parent]) {
            ret = MSP_ERR_BAD_NODE_TIME_ORDERING;
            goto out;
        }
    }
    /* Check the samples */
    for (j = 0; j < self->num_samples; j++) {
        if (self->samples[j] < 0 || self->samples[j] >= num_nodes) {
            ret = MSP_ERR_NODE_OUT_OF_BOUNDS;
            goto out;
        }
        if (!(self->nodes->flags[self->samples[j]] & MSP_NODE_IS_SAMPLE)) {
            ret = MSP_ERR_BAD_SAMPLES;
            goto out;
        }

    }
    /* Check the sites */
    for (j = 0; j < self->sites->num_rows; j++) {
        if (self->sites->position[j] < 0
                || self->sites->position[j] >= self->sequence_length) {
            ret = MSP_ERR_BAD_SITE_POSITION;
            goto out;
        }
    }
    /* Check the mutations */
    num_sites = (site_id_t) self->sites->num_rows;
    for (j = 0; j < self->mutations->num_rows; j++) {
        if (self->mutations->site[j] < 0 || self->mutations->site[j] >= num_sites) {
            ret = MSP_ERR_SITE_OUT_OF_BOUNDS;
            goto out;
        }
        if (self->mutations->node[j] < 0 || self->mutations->node[j] >= num_nodes) {
            ret = MSP_ERR_NODE_OUT_OF_BOUNDS;
            goto out;
        }
    }
    ret = 0;
out:
    msp_safe_free(node_seen);
    return ret;
}

static int WARN_UNUSED
simplifier_insert_sample(simplifier_t *self, node_id_t sample_id)
{
    int ret = 0;

    assert(self->ancestor_map[sample_id] == NULL);
    self->ancestor_map[sample_id] = simplifier_alloc_segment(self, 0,
            self->sequence_length, (node_id_t) self->nodes->num_rows, NULL);
    if (self->ancestor_map[sample_id] == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = simplifier_record_node(self, sample_id, true);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

static int
simplifier_init_samples(simplifier_t *self, node_id_t *samples)
{
    int ret = 0;
    size_t j;

    /* Go through the samples to check for errors. */
    for (j = 0; j < self->num_samples; j++) {
        if (samples[j] < 0 || samples[j] > (node_id_t) self->input_nodes.num_rows) {
            ret = MSP_ERR_OUT_OF_BOUNDS;
            goto out;
        }
        if (self->is_sample[samples[j]]) {
            ret = MSP_ERR_DUPLICATE_SAMPLE;
            goto out;
        }
        self->is_sample[samples[j]] = true;
        ret = simplifier_insert_sample(self, samples[j]);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

int
simplifier_alloc(simplifier_t *self, double sequence_length,
        node_id_t *samples, size_t num_samples,
        node_table_t *nodes, edge_table_t *edges, migration_table_t *migrations,
        site_table_t *sites, mutation_table_t *mutations,
        size_t max_buffered_edges, int flags)
{
    int ret = 0;
    size_t j, offset, max_alloc_block, num_nodes_alloc, num_edges_alloc;

    memset(self, 0, sizeof(simplifier_t));
    self->samples = samples;
    self->num_samples = num_samples;
    self->flags = flags;
    self->nodes = nodes;
    self->edges = edges;
    self->sites = sites;
    self->mutations = mutations;


    if (nodes == NULL || edges == NULL || samples == NULL
            || sites == NULL || mutations == NULL || migrations == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (sequence_length == 0) {
        /* infer sequence length from the edges */
        sequence_length = 0.0;
        for (j = 0; j < edges->num_rows; j++) {
            sequence_length = GSL_MAX(sequence_length, edges->right[j]);
        }
        if (sequence_length <= 0.0) {
            ret = MSP_ERR_BAD_SEQUENCE_LENGTH;
            goto out;
        }
    }
    self->sequence_length = sequence_length;

    /* If we have more then 256K blocks or edges just allocate this much */
    max_alloc_block = 256 * 1024;
    /* Need to avoid malloc(0) so make sure we have at least 1. */
    num_nodes_alloc = GSL_MAX(max_alloc_block, 1 + nodes->num_rows);
    num_edges_alloc = GSL_MAX(max_alloc_block, 1 + edges->num_rows);

    /* TODO we can add a flag to skip these checks for when we know they are
     * unnecessary */
    ret = simplifier_check_input(self);
    if (ret != 0) {
        goto out;
    }

    /* Make a copy of the input nodes and clear the table ready for output */
    ret = node_table_alloc(&self->input_nodes, nodes->num_rows, nodes->total_name_length);
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
    /* Build the offset table so we can map node names */
    self->node_name_offset = malloc(num_nodes_alloc * sizeof(size_t));
    if (self->node_name_offset == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    offset = 0;
    for (j = 0; j < self->input_nodes.num_rows; j++) {
        self->node_name_offset[j] = offset;
        offset += self->input_nodes.name_length[j];
    }

    /* Make a copy of the input edges and clear the input table, ready for output. */
    ret = edge_table_alloc(&self->input_edges, edges->num_rows);
    if (ret != 0) {
        goto out;
    }
    ret = edge_table_set_columns(&self->input_edges, edges->num_rows,
            edges->left, edges->right, edges->parent, edges->child);
    if (ret != 0) {
        goto out;
    }
    ret = edge_table_reset(self->edges);
    if (ret != 0) {
        goto out;
    }
    if (max_buffered_edges == 0) {
        max_buffered_edges = 1024;
    }
    self->max_buffered_edges = max_buffered_edges;
    self->num_buffered_edges = 0;
    self->edge_buffer = malloc(self->max_buffered_edges * sizeof(edge_t));
    if (self->edge_buffer == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    /* Make a copy of the input sites and clear the input table, ready for output. */
    ret = site_table_alloc(&self->input_sites, sites->num_rows,
            sites->total_ancestral_state_length);
    if (ret != 0) {
        goto out;
    }
    ret = site_table_set_columns(&self->input_sites, sites->num_rows,
            sites->position, sites->ancestral_state, sites->ancestral_state_length);
    if (ret != 0) {
        goto out;
    }
    ret = site_table_reset(self->sites);
    if (ret != 0) {
        goto out;
    }

    /* Make a copy of the input mutations and clear the input table, ready for output. */
    ret = mutation_table_alloc(&self->input_mutations, mutations->num_rows,
            mutations->total_derived_state_length);
    if (ret != 0) {
        goto out;
    }
    ret = mutation_table_set_columns(&self->input_mutations, mutations->num_rows,
            mutations->site, mutations->node, mutations->parent,
            mutations->derived_state, mutations->derived_state_length);
    if (ret != 0) {
        goto out;
    }
    ret = mutation_table_reset(self->mutations);
    if (ret != 0) {
        goto out;
    }

    /* Allocate the heaps used for small objects. */
    /* TODO assuming that the number of edges is a good guess here. */
    ret = object_heap_init(&self->segment_heap, sizeof(simplify_segment_t),
            num_edges_alloc, NULL);
    if (ret != 0) {
        goto out;
    }
    ret = object_heap_init(&self->avl_node_heap, sizeof(avl_node_t),
            num_edges_alloc, NULL);
    if (ret != 0) {
        goto out;
    }
    /* Make the maps and set the intial state */
    self->ancestor_map = calloc(num_nodes_alloc, sizeof(simplify_segment_t *));
    self->node_id_map = malloc(num_nodes_alloc * sizeof(node_id_t));
    self->is_sample = calloc(num_nodes_alloc, sizeof(bool));
    if (self->ancestor_map == NULL || self->node_id_map == NULL
            || self->is_sample == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memset(self->node_id_map, 0xff, self->input_nodes.num_rows * sizeof(node_id_t));
    self->nodes->num_rows = 0;
    avl_init_tree(&self->merge_queue, cmp_segment_queue, NULL);

    /* Allocate the segment buffers */
    self->segment_buffer_size = 64;
    self->segment_buffer = malloc(self->segment_buffer_size * sizeof(simplify_segment_t *));
    if (self->segment_buffer == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = simplifier_init_sites(self);
    if (ret != 0) {
        goto out;
    }
    ret = simplifier_init_samples(self, samples);
out:
    return ret;
}

int
simplifier_free(simplifier_t *self)
{
    node_table_free(&self->input_nodes);
    edge_table_free(&self->input_edges);
    site_table_free(&self->input_sites);
    mutation_table_free(&self->input_mutations);
    object_heap_free(&self->segment_heap);
    object_heap_free(&self->avl_node_heap);
    msp_safe_free(self->node_name_offset);
    msp_safe_free(self->ancestor_map);
    msp_safe_free(self->node_id_map);
    msp_safe_free(self->segment_buffer);
    msp_safe_free(self->edge_buffer);
    msp_safe_free(self->is_sample);
    msp_safe_free(self->mutation_id_map);
    msp_safe_free(self->mutation_node_map);
    msp_safe_free(self->mutation_position_map);
    msp_safe_free(self->mutation_node_list_mem);
    msp_safe_free(self->mutation_position_map_mem);
    return 0;
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
simplifier_defrag_segment_chain(simplifier_t *self, node_id_t input_id)
{
    simplify_segment_t *x_prev, *x;

    x_prev = self->ancestor_map[input_id];
    assert(x_prev != NULL);
    x = x_prev->next;
    while (x != NULL) {
        if (x_prev->right == x->left && x_prev->node == x->node) {
            x_prev->right = x->right;
            x_prev->next = x->next;
            simplifier_free_segment(self, x);
            x = x_prev->next;
        } else {
            x_prev = x;
            x = x->next;
        }
    }
    return 0;
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
simplifier_record_mutations(simplifier_t *self, node_id_t input_id)
{
    int ret = 0;
    simplify_segment_t *seg = self->ancestor_map[input_id];
    avl_tree_t *avl_tree = &self->mutation_position_map[input_id];
    avl_node_t *a;
    mutation_position_map_t *mpm, search;
    mutation_node_list_t *mnl;
    /* The AVL tree is overkill here really. All we want is a sorted list of the
     * position -> [mutation_ids] for each input ID. We can then co-iterate over
     * the ancestry segments and the list of positions. This will be more efficient
     * unless we have huge numbers of mutations per node */

    if (avl_count(avl_tree) > 0) {
        while (seg != NULL) {
            search.position = seg->left;
            avl_search_closest(avl_tree, &search, &a);
            assert(a != NULL);
            mpm = (mutation_position_map_t *) a->item;
            if (mpm->position < seg->left) {
                a = a->next;
            }
            while (a != NULL
                    && ((mutation_position_map_t *) a->item)->position < seg->right) {
                mpm = (mutation_position_map_t *) a->item;
                assert(mpm->position >= seg->left && mpm->position < seg->right);
                for (mnl = mpm->head; mnl != NULL; mnl = mnl->next) {
                    /* Set the output node for this mutation to the segment's node */
                    self->mutation_node_map[mnl->mutation_id] = seg->node;
                }
                a = a->next;
            }
            seg = seg->next;
        }
    }
    return ret;
}

static int
simplifier_remove_ancestry(simplifier_t *self, double left, double right,
        node_id_t input_id)
{
    int ret = 0;
    simplify_segment_t *x, *y, *head, *last, *x_prev;

    head = self->ancestor_map[input_id];
    x = head;
    last = NULL;
    x_prev = NULL;  /* Keep the compiler happy */
    /* Skip the leading segments before left */
    while (x != NULL && x->right <= left) {
        last = x;
        x = x->next;
    }
    if (x != NULL && x->left < left) {
        /* The left edge of x overhangs. Insert a new segment for the excess. */
        y = simplifier_alloc_segment(self, x->left, left, x->node, NULL);
        if (y == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        x->left = left;
        if (last != NULL) {
            last->next = y;
        }
        last = y;
        if (x == head) {
            head = last;
        }
    }
    if (x != NULL && x->left < right) {
        /* x is the first segment within the target interval, so add it to the
         * output queue */
        ret = simplifier_priority_queue_insert(self, &self->merge_queue, x);
        if (ret != 0) {
            goto out;
        }
        /* Skip over segments strictly within the target interval */
        while (x != NULL && x->right <= right) {
            x_prev = x;
            x = x->next;
        }
        if (x != NULL && x->left < right) {
            /* We have an overhang on the right hand side. Create a new
             * segment for the overhang and terminate the output chain. */
            y = simplifier_alloc_segment(self, right, x->right, x->node, x->next);
            if (y == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            x->right = right;
            x->next = NULL;
            x = y;
        } else if (x_prev != NULL) {
            x_prev->next = NULL;
        }
    }
    /* x is the first segment in the new chain starting after right. */
    if (last == NULL) {
        head = x;
    } else {
        last->next = x;
    }
    self->ancestor_map[input_id] = head;
out:
    return ret;
}

static int WARN_UNUSED
simplifier_merge_ancestors(simplifier_t *self, node_id_t input_id)
{
    int ret = MSP_ERR_GENERIC;
    bool coalescence = false;
    bool defrag_required = 0;
    node_id_t v;
    double l, r, next_l;
    uint32_t j, h;
    avl_node_t *node;
    simplify_segment_t head_sentinel;
    simplify_segment_t *x, *z, *alpha, *head;
    simplify_segment_t **H = NULL;
    avl_tree_t *Q = &self->merge_queue;

    H = simplifier_get_segment_buffer(self, avl_count(Q));
    if (H == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    head_sentinel.left = 0.0;
    head_sentinel.right = 0.0;
    head_sentinel.node = MSP_NULL_NODE;
    head_sentinel.next = NULL;
    head = &head_sentinel;
    z = head;
    while (avl_count(Q) > 0) {
        h = 0;
        node = Q->head;
        l = ((simplify_segment_t *) node->item)->left;
        r = self->sequence_length;
        while (node != NULL && ((simplify_segment_t *) node->item)->left == l) {
            H[h] = (simplify_segment_t *) node->item;
            r = GSL_MIN(r, H[h]->right);
            h++;
            simplifier_free_avl_node(self, node);
            avl_unlink_node(Q, node);
            node = node->next;
        }
        next_l = 0;
        if (node != NULL) {
            next_l = ((simplify_segment_t *) node->item)->left;
            r = GSL_MIN(r, next_l);
        }
        if (h == 1) {
            x = H[0];
            if (node != NULL && next_l < x->right) {
                alpha = simplifier_alloc_segment(self, x->left, next_l, x->node, NULL);
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
            v = self->node_id_map[input_id];
            if (v >= 0 && self->is_sample[input_id]) {
                /* If we have a mapped node over an interval with no other segments,
                 * then we must record this edge as it joins internal samples */
                ret = simplifier_record_edge(self, alpha->left, alpha->right, v,
                        alpha->node);
                if (ret != 0) {
                    goto out;
                }
                /* The node on this segment must be remapped to the parent so that we
                 * have the correct nodes further up in the tree. */
                alpha->node = v;
            }
        } else {
            if (!coalescence) {
                coalescence = true;
                if (self->node_id_map[input_id] == MSP_NULL_NODE) {
                    ret = simplifier_record_node(self, input_id, false);
                    if (ret != 0) {
                        goto out;
                    }
                }
            }
            v = self->node_id_map[input_id];
            alpha = simplifier_alloc_segment(self, l, r, v, NULL);
            if (alpha == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            /* Create the record and update the priority queue */
            for (j = 0; j < h; j++) {
                x = H[j];
                ret = simplifier_record_edge(self, l, r, v, x->node);
                if (ret != 0) {
                    goto out;
                }
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
        }

        /* Loop tail; integrate alpha into the global state */
        assert(alpha != NULL);
        if (z == NULL) {
            self->ancestor_map[input_id] = alpha;
        } else {
        }
        defrag_required |= z->right == alpha->left && z->node == alpha->node;
        z->next = alpha;
        z = alpha;
    }
    if (self->ancestor_map[input_id] == NULL) {
        self->ancestor_map[input_id] = head->next;
    } else {
        /* There is already ancestral material present for this node, which
         * is a sample. We can free the segments computed here because the
         * existing mapping covers the full region */
        x = head->next;
        while (x != NULL) {
            assert(x->node == self->node_id_map[input_id]);
            simplifier_free_segment(self ,x);
            x = x->next;
        }
        x = self->ancestor_map[input_id];
        assert(self->is_sample[input_id]);
        assert(x->left == 0.0);
        assert(x->right == self->sequence_length);
        assert(x->node == self->node_id_map[input_id]);
    }
    if (defrag_required) {
        ret = simplifier_defrag_segment_chain(self, input_id);
        if (ret != 0) {
            goto out;
        }
    }
    ret = simplifier_flush_edges(self);
out:
    return ret;
}

static int WARN_UNUSED
simplifier_output_sites(simplifier_t *self)
{
    int ret = 0;
    site_id_t input_site;
    mutation_id_t input_mutation, mapped_parent ,site_start, site_end;
    site_id_t num_input_sites = (site_id_t) self->input_sites.num_rows;
    mutation_id_t num_input_mutations = (mutation_id_t) self->input_mutations.num_rows;
    mutation_id_t input_parent, num_output_mutations, num_output_site_mutations;
    node_id_t mapped_node;
    bool keep_mutation, keep_site;
    bool filter_zero_mutation_sites = (self->flags & MSP_FILTER_ZERO_MUTATION_SITES);

    /* TODO Implement the checks below for ancestral state and derived state properly when
     * we've added the _offset columns to the tables. */
    assert(self->input_sites.total_ancestral_state_length
            == self->input_sites.num_rows);
    assert(self->input_mutations.total_derived_state_length
            == self->input_mutations.num_rows);

    input_mutation = 0;
    num_output_mutations = 0;
    for (input_site = 0; input_site < num_input_sites; input_site++) {
        site_start = input_mutation;
        num_output_site_mutations = 0;
        while (input_mutation < num_input_mutations
                && self->input_mutations.site[input_mutation] == input_site) {
            mapped_node = self->mutation_node_map[input_mutation];
            if (mapped_node != MSP_NULL_NODE) {
                input_parent = self->input_mutations.parent[input_mutation];
                mapped_parent = MSP_NULL_MUTATION;
                if (input_parent != MSP_NULL_MUTATION) {
                    mapped_parent = self->mutation_id_map[input_parent];
                }
                keep_mutation = true;
                if (mapped_parent == MSP_NULL_MUTATION) {
                    /* If there is no parent and the ancestral state is equal to the
                     * derived state, then we remove this mutation.
                     */
                    if (self->input_mutations.derived_state[input_mutation]
                            == self->input_sites.ancestral_state[input_site]) {
                        keep_mutation = false;
                    }
                }
                if (keep_mutation) {
                    self->mutation_id_map[input_mutation] = num_output_mutations;
                    num_output_mutations++;
                    num_output_site_mutations++;
                }
            }
            input_mutation++;
        }
        site_end = input_mutation;

        keep_site = true;
        if (filter_zero_mutation_sites && num_output_site_mutations == 0) {
            keep_site = false;
        }
        if (keep_site) {
            for (input_mutation = site_start; input_mutation < site_end; input_mutation++) {
                if (self->mutation_id_map[input_mutation] != MSP_NULL_MUTATION) {
                    assert(self->mutations->num_rows
                            == (size_t) self->mutation_id_map[input_mutation]);
                    mapped_node = self->mutation_node_map[input_mutation];
                    assert(mapped_node != MSP_NULL_NODE);
                    mapped_parent = self->input_mutations.parent[input_mutation];
                    if (mapped_parent != MSP_NULL_MUTATION) {
                        mapped_parent = self->mutation_id_map[mapped_parent];
                    }
                    ret = mutation_table_add_row(self->mutations,
                            (site_id_t) self->sites->num_rows,
                            mapped_node, mapped_parent,
                            /* FIXME do this properly when derived_state_offset is
                             * implemented */
                            self->input_mutations.derived_state + input_mutation,
                            self->input_mutations.derived_state_length[input_mutation]);
                    if (ret != 0) {
                        goto out;
                    }
                }
            }
            ret = site_table_add_row(self->sites,
                    self->input_sites.position[input_site],
                    /* FIXME do this properly when ancestral_state_offset is implemented */
                    self->input_sites.ancestral_state + input_site,
                    self->input_sites.ancestral_state_length[input_site]);
            if (ret != 0) {
                goto out;
            }

        }
        assert(num_output_mutations == (mutation_id_t) self->mutations->num_rows);
        input_mutation = site_end;
    }
    assert(input_mutation == num_input_mutations);
out:
    return ret;
}

static int WARN_UNUSED
simplifier_process_parent_edges(simplifier_t *self, node_id_t parent, size_t start,
        size_t end)
{
    int ret = 0;
    size_t j;
    node_id_t child;
    double left, right;

    /* Go through the edges and remove ancestry */
    for (j = start; j < end; j++) {
        assert(parent == self->input_edges.parent[j]);
        child = self->input_edges.child[j];
        if (self->ancestor_map[child] != NULL) {
            ret = simplifier_record_mutations(self, child);
            if (ret != 0) {
                goto out;
            }
            left = self->input_edges.left[j];
            right = self->input_edges.right[j];
            ret = simplifier_remove_ancestry(self, left, right, child);
            if (ret != 0) {
                goto out;
            }
        }
    }
    /* We can now merge the ancestral segments for the parent */
    ret = simplifier_merge_ancestors(self, parent);
    if (ret != 0) {
        goto out;
    }
    assert(avl_count(&self->merge_queue) == 0);
out:
    return ret;
}

int WARN_UNUSED
simplifier_run(simplifier_t *self, node_id_t *node_map)
{
    int ret = 0;
    size_t j, start;
    node_id_t parent, current_parent;
    size_t num_edges = self->input_edges.num_rows;

    if (num_edges > 0) {
        start = 0;
        current_parent = self->input_edges.parent[0];
        for (j = 0; j < num_edges; j++) {
            parent = self->input_edges.parent[j];
            if (parent != current_parent) {
                ret = simplifier_process_parent_edges(self, current_parent, start, j);
                if (ret != 0) {
                    goto out;
                }
                current_parent = parent;
                start = j;
            }
        }
        ret = simplifier_process_parent_edges(self, current_parent, start, num_edges);
        if (ret != 0) {
            goto out;
        }
    }
    /* Record any remaining mutations over the roots */
    for (j = 0; j < self->input_nodes.num_rows; j++) {
        if (self->ancestor_map[j] != NULL) {
            ret = simplifier_record_mutations(self, (node_id_t) j);
            if (ret != 0) {
                goto out;
            }
        }
    }
    ret = simplifier_output_sites(self);
    if (ret != 0) {
        goto out;
    }
    if (node_map != NULL) {
        /* Finally, output the new IDs for the nodes, if required. */
        memcpy(node_map, self->node_id_map, self->input_nodes.num_rows * sizeof(node_id_t));
    }
out:
    return ret;
}
