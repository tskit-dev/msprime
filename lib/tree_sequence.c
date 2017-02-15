/*
** Copyright (C) 2015-2017 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

#include <hdf5.h>

#include <gsl/gsl_math.h>

#include "err.h"
#include "msprime.h"

#define MSP_DIR_FORWARD 1
#define MSP_DIR_REVERSE -1

typedef struct {
    double value;
    node_id_t index;
    int64_t time;
} index_sort_t;

static int
cmp_node_id_t(const void *a, const void *b) {
    const node_id_t *ia = (const node_id_t *) a;
    const node_id_t *ib = (const node_id_t *) b;
    return (*ia > *ib) - (*ia < *ib);
}

static int
cmp_double(const void *a, const void *b) {
    const double *ia = (const double *) a;
    const double *ib = (const double *) b;
    return (*ia > *ib) - (*ia < *ib);
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

static int
cmp_record_time_left(const void *a, const void *b) {
    const coalescence_record_t *ca = (const coalescence_record_t *) a;
    const coalescence_record_t *cb = (const coalescence_record_t *) b;
    int ret = (ca->time > cb->time) - (ca->time < cb->time);
    if (ret == 0) {
        ret = (ca->left > cb->left) - (ca->left < cb->left);
    }
    return ret;
}

static void
tree_sequence_check_state(tree_sequence_t *self)
{
    size_t j;

    for (j = 0; j < self->edgesets.num_records; j++) {
        assert(self->edgesets.num_children[j] >= 1);
    }
}

void
tree_sequence_print_state(tree_sequence_t *self, FILE *out)
{
    size_t j, num_mutations;
    int k;
    int ret;
    mutation_t *mutations;
    sparse_tree_t tree;

    fprintf(out, "tree_sequence state\n");
    fprintf(out, "sample_size = %d\n", (int) self->sample_size);
    fprintf(out, "provenance = (%d)\n", (int) self->num_provenance_strings);
    for (j = 0; j < self->num_provenance_strings; j++) {
        fprintf(out, "\t'%s'\n", self->provenance_strings[j]);
    }
    fprintf(out, "sequence_length = %f\n", self->sequence_length);
    fprintf(out, "nodes (%d)\n", (int) self->nodes.num_records);
    for (j = 0; j < self->nodes.num_records; j++) {
        fprintf(out, "\t%d\t%d\t%f\t%s\n", (int) j,
                (int) self->nodes.population[j],
                self->nodes.time[j],
                self->nodes.name[j]);
    }
    fprintf(out, "edgesets = (%d records)\n", (int) self->edgesets.num_records);
    for (j = 0; j < self->edgesets.num_records; j++) {
        fprintf(out, "\t%d\t%f\t%f\t%d\t(",
                (int) j,
                self->edgesets.left[j],
                self->edgesets.right[j],
                (int) self->edgesets.parent[j]);
        for (k = 0; k < self->edgesets.num_children[j]; k++) {
            fprintf(out, "%d", (int) self->edgesets.children[j][k]);
            if (k < self->edgesets.num_children[j] - 1) {
                fprintf(out, ", ");
            }
        }
        fprintf(out, ")\t|\t%d\t%d\n",
                (int) self->edgesets.indexes.insertion_order[j],
                (int) self->edgesets.indexes.removal_order[j]);
    }
    fprintf(out, "mutations = (%d records)\n", (int) self->mutations.num_records);
    fprintf(out, "\ttotal_nodes = %d\n", (int) self->mutations.total_nodes);
    for (j = 0; j < self->mutations.num_records; j++) {
        fprintf(out, "\t%d\t%f\t", (int) j, self->mutations.position[j]);
        for (k = 0; k < self->mutations.num_nodes[j]; k++) {
            fprintf(out, "\t%d,", (int) self->mutations.nodes[j][k]);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "migrations.records = (%d records)\n",
            (int) self->migrations.num_records);
    for (j = 0; j < self->migrations.num_records; j++) {
        fprintf(out, "\t%d\t%f\t%f\t%d\t%d\t%d\t%f\n", (int) j,
                self->migrations.left[j],
                self->migrations.right[j],
                (int) self->migrations.node[j],
                (int) self->migrations.source[j],
                (int) self->migrations.dest[j],
                self->migrations.time[j]);
    }
    fprintf(out, "memory\n");
    fprintf(out, "\tnodes.num_records = %d\n", (int) self->nodes.num_records);
    fprintf(out, "\tnodes.max_num_records = %d\n", (int) self->nodes.max_num_records);
    fprintf(out, "\tedgesets.num_records = %d\n", (int) self->edgesets.num_records);
    fprintf(out, "\tedgesets.max_num_records = %d\n", (int) self->edgesets.max_num_records);
    fprintf(out, "\tedgesets.total_children = %d\n", (int) self->edgesets.total_children);
    fprintf(out, "\tedgesets.max_total_children = %d\n",
            (int) self->edgesets.max_total_children);
    fprintf(out, "\tmutations.num_records = %d\n", (int) self->mutations.num_records);
    fprintf(out, "\tmutations.max_num_records = %d\n", (int) self->mutations.max_num_records);
    fprintf(out, "\tmutations.total_nodes = %d\n", (int) self->mutations.total_nodes);
    fprintf(out, "\tmutations.max_total_nodes = %d\n", (int) self->mutations.max_total_nodes);
    fprintf(out, "\tmigrations.num_records = %d\n", (int) self->migrations.num_records);
    fprintf(out, "\tmigrations.max_num_records = %d\n", (int) self->migrations.max_num_records);

    fprintf(out, "edgesets.\n");
    ret = sparse_tree_alloc(&tree, self, 0);
    assert(ret == 0);
    for (ret = sparse_tree_first(&tree); ret == 1; ret = sparse_tree_next(&tree)) {
        fprintf(out, "Tree : %d: %f-%f\n", (int) tree.index, tree.left, tree.right);
        ret = sparse_tree_get_mutations(&tree, &num_mutations, &mutations);
        assert(ret == 0);
        for (j = 0; j < num_mutations; j++) {
            fprintf(out, "\t\t%d\t%f\t", (int) mutations[j].index, mutations[j].position);
            for (k = 0; k < (int) mutations[j].num_nodes; k++) {
                fprintf(out, "%d,", (int) mutations[j].nodes[k]);

            }
            fprintf(out, "\n");
        }
    }
    assert(ret == 0);
    sparse_tree_free(&tree);
    tree_sequence_check_state(self);
}

static int
tree_sequence_alloc_mutations(tree_sequence_t *self)
{
    int ret = MSP_ERR_GENERIC;
    size_t size;

    if (self->mutations.total_nodes > self->mutations.max_total_nodes) {
        self->mutations.max_total_nodes = self->mutations.total_nodes;
        size = self->mutations.max_total_nodes;
        msp_safe_free(self->mutations.nodes_mem);
        self->mutations.nodes_mem = malloc(size * sizeof(node_id_t));
        if (self->mutations.nodes_mem == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    if (self->mutations.num_records > self->mutations.max_num_records) {
        self->mutations.max_num_records = self->mutations.num_records;
        size = self->mutations.max_num_records;
        msp_safe_free(self->mutations.nodes);
        msp_safe_free(self->mutations.num_nodes);
        msp_safe_free(self->mutations.position);
        msp_safe_free(self->mutations.ancestral_state);
        msp_safe_free(self->mutations.derived_state);
        msp_safe_free(self->mutations.tree_mutations_mem);
        self->mutations.nodes = malloc(size * sizeof(node_id_t *));
        self->mutations.num_nodes = malloc(size * sizeof(size_t));
        self->mutations.position = malloc(size * sizeof(double));
        self->mutations.ancestral_state = malloc(size * sizeof(char));
        self->mutations.derived_state = malloc(size * sizeof(char));
        self->mutations.tree_mutations_mem = malloc(size * sizeof(mutation_t));
        if (self->mutations.num_nodes == NULL
                || self->mutations.position == NULL
                || self->mutations.ancestral_state == NULL
                || self->mutations.derived_state == NULL
                || self->mutations.nodes == NULL
                || self->mutations.tree_mutations_mem == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_alloc_trees(tree_sequence_t *self)
{
    int ret = MSP_ERR_GENERIC;
    size_t size;

    if (self->nodes.total_name_length > self->nodes.max_total_name_length) {
        size = self->nodes.total_name_length;
        self->nodes.max_total_name_length = size;
        msp_safe_free(self->nodes.name_mem);
        self->nodes.name_mem = malloc(size * sizeof(char));
        if (self->nodes.name_mem == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    if (self->nodes.num_records > self->nodes.max_num_records) {
        size = self->nodes.num_records;
        self->nodes.max_num_records = size;
        msp_safe_free(self->nodes.time);
        msp_safe_free(self->nodes.population);
        msp_safe_free(self->nodes.flags);
        msp_safe_free(self->nodes.name);
        msp_safe_free(self->nodes.name_length);
        self->nodes.flags = malloc(size * sizeof(uint32_t));
        self->nodes.time = malloc(size * sizeof(double));
        self->nodes.population = malloc(size * sizeof(population_id_t));
        self->nodes.name = malloc(size * sizeof(char *));
        self->nodes.name_length = malloc(size * sizeof(size_t));
        if (self->nodes.flags == NULL
                || self->nodes.time == NULL
                || self->nodes.population == NULL
                || self->nodes.name == NULL
                || self->nodes.name_length == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    if (self->edgesets.num_records > self->edgesets.max_num_records) {
        size = self->edgesets.num_records;
        self->edgesets.max_num_records = size;
        msp_safe_free(self->edgesets.left);
        msp_safe_free(self->edgesets.right);
        msp_safe_free(self->edgesets.num_children);
        msp_safe_free(self->edgesets.children);
        msp_safe_free(self->edgesets.parent);
        msp_safe_free(self->edgesets.indexes.insertion_order);
        msp_safe_free(self->edgesets.indexes.removal_order);
        self->edgesets.left = malloc(size * sizeof(double));
        self->edgesets.right = malloc(size * sizeof(double));
        self->edgesets.num_children = malloc(size * sizeof(size_t));
        self->edgesets.children = malloc(size * sizeof(node_id_t *));
        self->edgesets.parent = malloc(size * sizeof(node_id_t));
        self->edgesets.indexes.insertion_order = malloc(size * sizeof(node_id_t));
        self->edgesets.indexes.removal_order = malloc(size * sizeof(node_id_t));
        if (self->edgesets.left == NULL
                || self->edgesets.right == NULL
                || self->edgesets.children == NULL
                || self->edgesets.parent == NULL
                || self->edgesets.num_children == NULL
                || self->edgesets.indexes.insertion_order == NULL
                || self->edgesets.indexes.removal_order == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    if (self->edgesets.total_children > self->edgesets.max_total_children) {
        size = self->edgesets.total_children;
        self->edgesets.max_total_children = size;
        msp_safe_free(self->edgesets.children_mem);
        self->edgesets.children_mem = malloc(size * sizeof(node_id_t));
        if (self->edgesets.children_mem == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_alloc_migrations(tree_sequence_t *self)
{
    int ret = MSP_ERR_GENERIC;
    size_t size;

    if (self->migrations.num_records > self->migrations.max_num_records) {
        size = self->migrations.num_records;
        self->migrations.max_num_records = size;
        msp_safe_free(self->migrations.node);
        msp_safe_free(self->migrations.source);
        msp_safe_free(self->migrations.dest);
        msp_safe_free(self->migrations.left);
        msp_safe_free(self->migrations.right);
        msp_safe_free(self->migrations.time);
        self->migrations.node = malloc(size * sizeof(node_id_t));
        self->migrations.source = malloc(size * sizeof(population_id_t));
        self->migrations.dest = malloc(size * sizeof(population_id_t));
        self->migrations.left = malloc(size * sizeof(double));
        self->migrations.right = malloc(size * sizeof(double));
        self->migrations.time = malloc(size * sizeof(double));
        if (self->migrations.node == NULL
                || self->migrations.source == NULL
                || self->migrations.dest == NULL
                || self->migrations.left == NULL
                || self->migrations.right == NULL
                || self->migrations.time == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_alloc_provenance(tree_sequence_t *self)
{
    int ret = MSP_ERR_GENERIC;

    if (self->num_provenance_strings > 0) {
        msp_safe_free(self->provenance_strings);
        self->provenance_strings = calloc(self->num_provenance_strings, sizeof(char *));
        if (self->provenance_strings == NULL) {
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

/* Allocates the memory required for arrays of values. Assumes that
 * the num_records and num_mutations have been set.
 */
static int
tree_sequence_alloc(tree_sequence_t *self)
{
    int ret = MSP_ERR_NO_MEMORY;

    ret = tree_sequence_alloc_trees(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_alloc_mutations(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_alloc_migrations(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_alloc_provenance(self);
    if (ret != 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_initialise(tree_sequence_t *self)
{
    memset(self, 0, sizeof(tree_sequence_t));
    self->initialised_magic = MSP_INITIALISED_MAGIC;
    return 0;
}

int
tree_sequence_free(tree_sequence_t *self)
{
    size_t j;

    if (self->provenance_strings != NULL) {
        for (j = 0; j < self->num_provenance_strings; j++) {
            free(self->provenance_strings[j]);
        }
        free(self->provenance_strings);
    }
    msp_safe_free(self->nodes.flags);
    msp_safe_free(self->nodes.population);
    msp_safe_free(self->nodes.time);
    msp_safe_free(self->nodes.name);
    msp_safe_free(self->nodes.name_mem);
    msp_safe_free(self->nodes.name_length);
    msp_safe_free(self->edgesets.left);
    msp_safe_free(self->edgesets.right);
    msp_safe_free(self->edgesets.children);
    msp_safe_free(self->edgesets.num_children);
    msp_safe_free(self->edgesets.children_mem);
    msp_safe_free(self->edgesets.parent);
    msp_safe_free(self->edgesets.indexes.insertion_order);
    msp_safe_free(self->edgesets.indexes.removal_order);
    msp_safe_free(self->mutations.nodes);
    msp_safe_free(self->mutations.num_nodes);
    msp_safe_free(self->mutations.position);
    msp_safe_free(self->mutations.ancestral_state);
    msp_safe_free(self->mutations.derived_state);
    msp_safe_free(self->mutations.nodes_mem);
    msp_safe_free(self->mutations.tree_mutations_mem);
    msp_safe_free(self->mutations.tree_mutations);
    msp_safe_free(self->mutations.num_tree_mutations);
    msp_safe_free(self->migrations.node);
    msp_safe_free(self->migrations.source);
    msp_safe_free(self->migrations.dest);
    msp_safe_free(self->migrations.left);
    msp_safe_free(self->migrations.right);
    msp_safe_free(self->migrations.time);
    return 0;
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
    node_id_t u, k, child, node;
    uint32_t j;
    size_t num_coordinates = self->edgesets.num_records + 1;
    double left, *result;
    double *coordinates = malloc(num_coordinates * sizeof(double));

    if (coordinates == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (j = 0; j < self->edgesets.num_records; j++) {
        coordinates[j] = self->edgesets.left[j];
    }
    coordinates[self->edgesets.num_records] = self->sequence_length;
    qsort(coordinates, num_coordinates, sizeof(double), cmp_double);

    if (coordinates[0] != 0.0) {
        /* TODO specific error for this */
        ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
        goto out;
    }

    left = DBL_MAX;
    for (j = 0; j < self->edgesets.num_records; j++) {
        node = self->edgesets.parent[j];
        if (node == MSP_NULL_NODE) {
            ret = MSP_ERR_NULL_NODE_IN_RECORD;
            goto out;
        }
        if (self->edgesets.num_children[j] < 1) {
            ret = MSP_ERR_ZERO_CHILDREN;
            goto out;
        }
        if (j > 0) {
            /* Input data must be time sorted. */
            if (self->nodes.time[node]
                    < self->nodes.time[self->edgesets.parent[j - 1]]) {
                ret = MSP_ERR_RECORDS_NOT_TIME_SORTED;
                goto out;
            }
        }
        left = GSL_MIN(left, self->edgesets.left[j]);
        for (k = 0; k < self->edgesets.num_children[j]; k++) {
            child = self->edgesets.children[j][k];
            assert(node != MSP_NULL_NODE);
            /* Children must be in ascending order */
            if (k < self->edgesets.num_children[j] - 1) {
                if (child >= self->edgesets.children[j][k + 1]) {
                    ret = MSP_ERR_UNSORTED_CHILDREN;
                    goto out;
                }
            }
            /* time[child] must be < time[parent] */
            if (self->nodes.time[child] >= self->nodes.time[node]) {
                ret = MSP_ERR_BAD_NODE_TIME_ORDERING;
                goto out;
            }
        }
        if (self->edgesets.left[j] >= self->edgesets.right[j]) {
            ret = MSP_ERR_BAD_RECORD_INTERVAL;
            goto out;
        }
        result = bsearch(self->edgesets.right + j, coordinates, num_coordinates,
                sizeof(double), cmp_double);
        if (result == NULL) {
            ret = MSP_ERR_BAD_COALESCENCE_RECORD_NONMATCHING_RIGHT;
            goto out;
        }
    }
    if (self->edgesets.num_records > 0 && left != 0) {
        ret = MSP_ERR_BAD_COALESCENCE_RECORDS;
        goto out;
    }

    /* Check the mutations */
    for (j = 0; j < self->mutations.num_records; j++) {
        if (self->mutations.position[j] < 0
                || self->mutations.position[j] >= self->sequence_length
                || self->mutations.num_nodes[j] < 1) {
            ret = MSP_ERR_BAD_MUTATION;
            goto out;
        }
        for (k = 0; k < self->mutations.num_nodes[j]; k++) {
            u = self->mutations.nodes[j][k];
            if (u < 0 || u >= (node_id_t) self->nodes.num_records) {
                ret = MSP_ERR_BAD_MUTATION;
                goto out;
            }
            if (k > 0) {
                if (u < self->mutations.nodes[j][k - 1]) {
                    ret = MSP_ERR_UNSORTED_MUTATION_NODES;
                    goto out;
                }
                if (u == self->mutations.nodes[j][k - 1]) {
                    ret = MSP_ERR_DUPLICATE_MUTATION_NODES;
                    goto out;
                }
            }
        }
        if (j > 0) {
            if (self->mutations.position[j] < self->mutations.position[j - 1]) {
                ret = MSP_ERR_MUTATIONS_NOT_POSITION_SORTED;
                goto out;
            }
        }
    }
    ret = 0;
out:
    if (coordinates != NULL) {
        free(coordinates);
    }
    return ret;
}

static int
tree_sequence_init_nodes(tree_sequence_t *self)
{
    size_t j, k;
    int ret = 0;

    self->sample_size = 0;
    k = 0;
    for (j = 0; j < self->nodes.num_records; j++) {
        if (self->nodes.flags[j] & MSP_NODE_SAMPLE) {
            self->sample_size++;
        }
        self->nodes.name_length[j] = 0;
        self->nodes.name[j] = self->nodes.name_mem + k;
        while (k < self->nodes.total_name_length && self->nodes.name_mem[k] != '\0') {
            self->nodes.name_length[j]++;
            k++;
        }
        k++;
    }
    if (self->nodes.num_records > 0 && self->sample_size == 0) {
        ret = MSP_ERR_BAD_COALESCENCE_RECORDS_SAMPLE_SIZE;
        goto out;
    }
out:
    return ret;
}

static int
tree_sequence_init_edgesets(tree_sequence_t *self)
{
    size_t j, offset;

    offset = 0;
    self->sequence_length = 0.0;
    for (j = 0; j < self->edgesets.num_records; j++) {
        self->sequence_length = GSL_MAX(self->sequence_length, self->edgesets.right[j]);
        assert(offset < self->edgesets.total_children);
        self->edgesets.children[j] = self->edgesets.children_mem + offset;
        offset += (size_t) self->edgesets.num_children[j];
    }
    return 0;
}

/* Initialises memory associated with mutations.
 */
static int
tree_sequence_init_mutations(tree_sequence_t *self)
{
    int ret = 0;
    size_t j, offset;

    offset = 0;
    for (j = 0; j < self->mutations.num_records; j++) {
        assert(offset < self->mutations.total_nodes);
        self->mutations.nodes[j] = self->mutations.nodes_mem + offset;
        offset += (size_t) self->mutations.num_nodes[j];
    }

    /* TODO remove this when we support storing mutation states */
    memset(self->mutations.ancestral_state, '0',
            self->mutations.num_records * sizeof(char));
    memset(self->mutations.derived_state, '1',
            self->mutations.num_records * sizeof(char));
    return ret;
}

/* Initialiases memory associated with the trees.
 */
static int
tree_sequence_init_trees(tree_sequence_t *self)
{
    int ret = MSP_ERR_GENERIC;
    size_t j, k, tree_index;
    mutation_t *mut;
    double last_x = -1;
    double x;
    node_id_t *I = self->edgesets.indexes.insertion_order;

    self->num_trees = 0;
    for (j = 0; j < self->edgesets.num_records; j++) {
        x = self->edgesets.left[I[j]];
        if (x != last_x) {
            self->num_trees++;
            last_x = x;
        }
    }
    if (self->num_trees > 0) {
        /* TODO this is an ugly departure from the other patterns of
         * mallocing and using high-water mark memory semantics. Do we really need
         * to have these?
         */
        msp_safe_free(self->mutations.num_tree_mutations);
        msp_safe_free(self->mutations.tree_mutations);
        self->mutations.num_tree_mutations = malloc(self->num_trees * sizeof(size_t));
        self->mutations.tree_mutations = malloc(self->num_trees * sizeof(mutation_t *));
        if (self->mutations.num_tree_mutations == NULL
                || self->mutations.tree_mutations == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        memset(self->mutations.num_tree_mutations, 0, self->num_trees * sizeof(size_t));
        memset(self->mutations.tree_mutations, 0, self->num_trees * sizeof(mutation_t *));
        if (self->edgesets.num_records > 0) {
            self->mutations.tree_mutations[0] = self->mutations.tree_mutations_mem;
        }
        /* Initialise the tree mutation records */
        for (j = 0; j < self->mutations.num_records; j++) {
            mut = &self->mutations.tree_mutations_mem[j];
            mut->index = j;
            mut->position = self->mutations.position[j];
            mut->ancestral_state = self->mutations.ancestral_state[j];
            mut->derived_state = self->mutations.derived_state[j];
            mut->num_nodes = (size_t) self->mutations.num_nodes[j];
            mut->nodes = self->mutations.nodes[j];
        }
        /* Now update the pointers so we have direct access within the trees */
        tree_index = 0;
        last_x = 0;
        k = 0;
        for (j = 0; j < self->edgesets.num_records; j++) {
            x = self->edgesets.left[I[j]];
            if (x != last_x) {
                self->mutations.tree_mutations[tree_index] =
                    self->mutations.tree_mutations_mem + k;
                last_x = x;
                while (k < self->mutations.num_records
                        && self->mutations.position[k] < x) {
                    self->mutations.num_tree_mutations[tree_index]++;
                    k++;
                }
                tree_index++;
            }
        }
        self->mutations.tree_mutations[tree_index] = self->mutations.tree_mutations_mem + k;
        while (k < self->mutations.num_records &&
                self->mutations.position[k] < self->sequence_length) {
            self->mutations.num_tree_mutations[tree_index]++;
            k++;
        }

    }
    ret = 0;
out:
    return ret;
}

static int WARN_UNUSED
tree_sequence_store_provenance_strings(tree_sequence_t *self,
        size_t num_provenance_strings, char**provenance_strings)
{
    int ret = MSP_ERR_GENERIC;
    char *s;
    size_t j, size;

    ret = tree_sequence_alloc_provenance(self);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < num_provenance_strings; j++) {
        if (provenance_strings[j] == NULL) {
            ret = MSP_ERR_BAD_PARAM_VALUE;
            goto out;
        }
        size = strlen(provenance_strings[j]);
        size++; /* allow for '/0' */
        s = malloc((size) * sizeof(char));
        if (s == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        strncpy(s, provenance_strings[j], size);
        self->provenance_strings[j] = s;
    }
    ret = 0;
out:
    return ret;
}

/* Temporary interface used to translate into load_tables for the simplify
 * function. Remove once simplify has been updated.
 */
static int WARN_UNUSED
tree_sequence_load_records(tree_sequence_t *self,
        size_t num_samples, sample_t *samples,
        size_t num_coalescence_records, coalescence_record_t *records,
        size_t num_mutations, mutation_t *mutations)
{
    int ret = MSP_ERR_GENERIC;
    node_table_t *node_table = NULL;
    edgeset_table_t *edgeset_table = NULL;
    mutation_table_t *mutation_table = NULL;
    size_t j;
    node_id_t last_node;
    coalescence_record_t *cr;

    node_table = malloc(sizeof(node_table_t));
    if (node_table == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = node_table_alloc(node_table, num_samples + num_coalescence_records, 1);
    if (ret != 0) {
        goto out;
    }
    edgeset_table = malloc(sizeof(edgeset_table_t));
    if (edgeset_table == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = edgeset_table_alloc(edgeset_table, num_coalescence_records,
            2 * num_coalescence_records);
    if (ret != 0) {
        goto out;
    }
    mutation_table = malloc(sizeof(mutation_table_t));
    if (mutation_table == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = mutation_table_alloc(mutation_table, num_mutations + 1, num_mutations + 1);
    if (ret != 0) {
        goto out;
    }

    for (j = 0; j < num_samples; j++) {
        ret = node_table_add_row(node_table, MSP_NODE_SAMPLE,
                samples[j].time, samples[j].population_id, "");
        if (ret != 0) {
            goto out;
        }
    }
    last_node = 0;
    for (j = 0; j < num_coalescence_records; j++) {
        cr = &records[j];
        if (cr->node != last_node) {
            assert(cr->node > last_node);
            last_node = cr->node;
            ret = node_table_add_row(node_table, 0, cr->time, cr->population_id, "");
            if (ret != 0) {
                goto out;
            }
        }
        ret = edgeset_table_add_row(edgeset_table, cr->left, cr->right,
                cr->node, cr->num_children, cr->children);
        if (ret != 0) {
            goto out;
        }
    }
    for (j = 0; j < num_mutations; j++) {
        ret = mutation_table_add_row(mutation_table, mutations[j].position,
                mutations[j].num_nodes, mutations[j].nodes);
        if (ret != 0) {
            goto out;
        }
    }
    ret = tree_sequence_load_tables_tmp(self, node_table, edgeset_table,
            NULL, mutation_table, 0, NULL);

out:
    if (node_table != NULL) {
        node_table_free(node_table);
        free(node_table);
    }
    if (edgeset_table != NULL) {
        edgeset_table_free(edgeset_table);
        free(edgeset_table);
    }
    if (mutation_table != NULL) {
        mutation_table_free(mutation_table);
        free(mutation_table);
    }
    return ret;
}

static int WARN_UNUSED
tree_sequence_build_indexes(tree_sequence_t *self)
{
    int ret = MSP_ERR_GENERIC;
    size_t j;
    double x;
    index_sort_t *sort_buff = NULL;

    sort_buff = malloc(self->edgesets.num_records * sizeof(index_sort_t));
    if (sort_buff == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    /* sort by left and increasing time to give us the order in which
     * records should be inserted */
    for (j = 0; j < self->edgesets.num_records; j++) {
        sort_buff[j].index = (node_id_t ) j;
        x = self->edgesets.left[j];
        sort_buff[j].value = x;
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
    qsort(sort_buff, self->edgesets.num_records, sizeof(index_sort_t), cmp_index_sort);
    for (j = 0; j < self->edgesets.num_records; j++) {
        self->edgesets.indexes.insertion_order[j] = sort_buff[j].index;
    }
    /* sort by right and decreasing time to give us the order in which
     * records should be removed. */
    for (j = 0; j < self->edgesets.num_records; j++) {
        sort_buff[j].index = (node_id_t ) j;
        x = self->edgesets.right[j];
        sort_buff[j].value = x;
        sort_buff[j].time = -1 * (int64_t ) j;
    }
    qsort(sort_buff, self->edgesets.num_records, sizeof(index_sort_t), cmp_index_sort);
    for (j = 0; j < self->edgesets.num_records; j++) {
        self->edgesets.indexes.removal_order[j] = sort_buff[j].index;
    }
    ret = 0;
out:
    if (sort_buff != NULL) {
        free(sort_buff);
    }
    return ret;
}

int WARN_UNUSED
tree_sequence_load_tables_tmp(tree_sequence_t *self,
    node_table_t *nodes, edgeset_table_t *edgesets, migration_table_t *migrations,
    mutation_table_t *mutations, size_t num_provenance_strings,
    char **provenance_strings)
{
    int ret = 0;
    size_t j, k, offset;

    /* TODO need to do a lot of input validation here. What do we allow to be
     * null? What are the size restrictions on the tables? */
    /* Do we allow zero nodes and edgesets?? */
    if (nodes == NULL || edgesets == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (edgesets->children_length <= edgesets->num_rows) {
        ret = MSP_ERR_BAD_CHILDREN_ARRAY;
        goto out;
    }

    self->num_provenance_strings = num_provenance_strings;
    self->nodes.num_records = nodes->num_rows;
    self->nodes.total_name_length = nodes->name_length;
    self->edgesets.total_children = edgesets->children_length - edgesets->num_rows;
    self->edgesets.num_records = edgesets->num_rows;
    self->mutations.num_records = 0;
    self->mutations.total_nodes = 0;
    if (mutations != NULL) {
        if (mutations->num_rows > 0 && mutations->nodes_length <= mutations->num_rows) {
            ret = MSP_ERR_BAD_NODES_ARRAY;
            goto out;
        }
        self->mutations.num_records = mutations->num_rows;
        self->mutations.total_nodes = mutations->nodes_length - mutations->num_rows;
    }
    self->migrations.num_records = 0;
    if (migrations != NULL) {
        self->migrations.num_records = migrations->num_rows;
    }
    ret = tree_sequence_alloc(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_store_provenance_strings(self, num_provenance_strings,
            provenance_strings);
    if (ret != 0) {
        goto out;
    }
    memcpy(self->nodes.time, nodes->time, nodes->num_rows * sizeof(double));
    memcpy(self->nodes.flags, nodes->flags, nodes->num_rows * sizeof(uint32_t));
    memcpy(self->nodes.population, nodes->population,
            nodes->num_rows * sizeof(population_id_t));
    memcpy(self->nodes.name_mem, nodes->name, nodes->name_length * sizeof(char));
    ret = tree_sequence_init_nodes(self);
    if (ret != 0) {
        goto out;
    }

    /* Setup the edgesets */
    memcpy(self->edgesets.left, edgesets->left, edgesets->num_rows * sizeof(double));
    memcpy(self->edgesets.right, edgesets->right, edgesets->num_rows * sizeof(double));
    memcpy(self->edgesets.parent, edgesets->parent, edgesets->num_rows * sizeof(node_id_t));
    k = 0;
    offset = 0;
    for (j = 0; j < edgesets->num_rows; j++) {
        self->edgesets.num_children[j] = 0;
        while (k < edgesets->children_length && edgesets->children[k] != MSP_NULL_NODE) {
            self->edgesets.num_children[j]++;
            self->edgesets.children_mem[offset] = edgesets->children[k];
            offset++;
            k++;
        }
        k++;
    }
    if (k != edgesets->children_length) {
        ret = MSP_ERR_BAD_CHILDREN_ARRAY;
        goto out;
    }
    ret = tree_sequence_init_edgesets(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_build_indexes(self);
    if (ret != 0) {
        goto out;
    }

    if (mutations != NULL) {
        memcpy(self->mutations.position, mutations->position,
                mutations->num_rows * sizeof(double));
        offset = 0;
        k = 0;
        for (j = 0; j < mutations->num_rows; j++) {
            self->mutations.num_nodes[j] = 0;
            while (k < mutations->nodes_length && mutations->nodes[k] != MSP_NULL_NODE) {
                self->mutations.num_nodes[j]++;
                self->mutations.nodes_mem[offset] = mutations->nodes[k];
                offset++;
                k++;
            }
            k++;
        }
        if (k != mutations->nodes_length) {
            ret = MSP_ERR_BAD_NODES_ARRAY;
            goto out;
        }
        ret = tree_sequence_init_mutations(self);
        if (ret != 0) {
            goto out;
        }
    }
    if (migrations != NULL) {
        /* Set up the migrations */
        memcpy(self->migrations.left, migrations->left, migrations->num_rows * sizeof(double));
        memcpy(self->migrations.right, migrations->right, migrations->num_rows * sizeof(double));
        memcpy(self->migrations.node, migrations->node, migrations->num_rows * sizeof(node_id_t));
        memcpy(self->migrations.source, migrations->source,
                migrations->num_rows * sizeof(population_id_t));
        memcpy(self->migrations.dest, migrations->dest,
                migrations->num_rows * sizeof(population_id_t));
        memcpy(self->migrations.time, migrations->time, migrations->num_rows * sizeof(double));
    }
    ret = tree_sequence_check(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_init_trees(self);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}


int WARN_UNUSED
tree_sequence_dump_tables_tmp(tree_sequence_t *self,
    node_table_t *nodes, edgeset_table_t *edgesets, migration_table_t *migrations,
    mutation_table_t *mutations, size_t *num_provenance_strings,
    char ***provenance_strings)
{
    int ret = 0;
    uint32_t flags;
    size_t j;
    double left, right;

    if (nodes == NULL || edgesets == NULL
            || num_provenance_strings == NULL || provenance_strings == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }

    ret = node_table_reset(nodes);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < self->nodes.num_records; j++) {
        flags = j < self->sample_size? MSP_NODE_SAMPLE: 0;
        ret = node_table_add_row(nodes, flags,
                self->nodes.time[j], self->nodes.population[j],
                self->nodes.name[j]);
        if (ret != 0) {
            goto out;
        }
    }

    ret = edgeset_table_reset(edgesets);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < self->edgesets.num_records; j++) {
        left = self->edgesets.left[j];
        right = self->edgesets.right[j];
        ret = edgeset_table_add_row(edgesets, left, right,
                self->edgesets.parent[j], (size_t) self->edgesets.num_children[j],
                self->edgesets.children[j]);
        if (ret != 0) {
            goto out;
        }
    }

    if (migrations != NULL) {
        ret = migration_table_reset(migrations);
        if (ret != 0) {
            goto out;
        }
        for (j = 0; j < self->migrations.num_records; j++) {
            ret = migration_table_add_row(migrations,
                    self->migrations.left[j],
                    self->migrations.right[j],
                    self->migrations.node[j],
                    self->migrations.source[j],
                    self->migrations.dest[j],
                    self->migrations.time[j]);
            if (ret != 0) {
                goto out;
            }
        }
    }

    if (mutations != NULL) {
        ret = mutation_table_reset(mutations);
        if (ret != 0) {
            goto out;
        }
        for (j = 0; j < self->mutations.num_records; j++) {
            ret = mutation_table_add_row(mutations,
                    self->mutations.position[j], (size_t) self->mutations.num_nodes[j],
                    self->mutations.nodes[j]);
            if (ret != 0) {
                goto out;
            }
        }
    }

    *num_provenance_strings = self->num_provenance_strings;
    *provenance_strings = self->provenance_strings;

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
    htri_t exists;
    struct _dimension_check {
        const char *name;
        int check_size;
        size_t size;
    };
    struct _dimension_check fields[] = {
        {"/mutations/position", 1, self->mutations.num_records},
        {"/mutations/num_nodes", 1, self->mutations.num_records},
        {"/mutations/nodes", 1, self->mutations.total_nodes},
        {"/nodes/flags", 1, self->nodes.num_records},
        {"/nodes/population", 1, self->nodes.num_records},
        {"/nodes/name_length", 1, self->nodes.num_records},
        {"/nodes/time", 1, self->nodes.num_records},
        {"/edgesets/left", 1, self->edgesets.num_records},
        {"/edgesets/right", 1, self->edgesets.num_records},
        {"/edgesets/parent", 1, self->edgesets.num_records},
        {"/edgesets/num_children", 1, self->edgesets.num_records},
        {"/edgesets/children", 0, self->edgesets.total_children},
        {"/edgesets/indexes/insertion_order", 1, self->edgesets.num_records},
        {"/edgesets/indexes/removal_order", 1, self->edgesets.num_records},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _dimension_check);
    size_t j;

    /* First make sure that the root number make sense */
    if (self->edgesets.num_records > 0) {
        if (self->nodes.num_records == 0) {
            ret = MSP_ERR_FILE_FORMAT;
            goto out;
        }
        if (self->edgesets.total_children == 0) {
            ret = MSP_ERR_FILE_FORMAT;
            goto out;
        }
    }
    if (self->mutations.num_records > 0 &&
            self->mutations.total_nodes < self->mutations.num_records) {
        ret = MSP_ERR_FILE_FORMAT;
        goto out;
    }
    /* Now go though the rest of the fields and make sure they have the
     * right sizes
     */
    for (j = 0; j < num_fields; j++) {
        exists = H5Lexists(file_id, fields[j].name, H5P_DEFAULT);
        if (exists < 0) {
            goto out;
        }
        dims[0] = 0;
        if (exists) {
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
            status = H5Sclose(dataspace_id);
            if (status < 0) {
                goto out;
            }
            status = H5Dclose(dataset_id);
            if (status < 0) {
                goto out;
            }
        }
        if (fields[j].check_size && dims[0] != fields[j].size) {
            ret = MSP_ERR_FILE_FORMAT;
            goto out;
        }
    }
    /* Make sure that the total number of nodes makes sense */
    if (self->mutations.total_nodes < self->mutations.num_records
            || self->edgesets.total_children < self->edgesets.num_records) {
        ret = MSP_ERR_FILE_FORMAT;
        goto out;
    }
    ret = 0;
out:
    return ret;
}

/* Reads the groups within the HDF5 file to ensure that they exist.
 */
static int
tree_sequence_read_hdf5_groups(tree_sequence_t *self, hid_t file_id)
{
    int ret = MSP_ERR_HDF5;
    htri_t exists;
    const char* groups[] = {
        "/edgesets/indexes",
        "/nodes",
        "/edgesets",
        "/mutations"
    };
    size_t num_groups = sizeof(groups) / sizeof(const char *);
    size_t j;

    for (j = 0; j < num_groups; j++) {
        exists = H5Lexists(file_id, groups[j], H5P_DEFAULT);
        if (exists < 0) {
            goto out;
        }
        if (! exists) {
            ret = MSP_ERR_FILE_FORMAT;
            goto out;
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
    hsize_t dims;
    struct _dimension_read {
        const char *name;
        size_t *dest;
    };
    size_t flattened_name_length;
    struct _dimension_read fields[] = {
        {"/mutations/position", &self->mutations.num_records},
        {"/mutations/nodes", &self->mutations.total_nodes},
        {"/provenance", &self->num_provenance_strings},
        {"/nodes/time", &self->nodes.num_records},
        {"/nodes/name", &flattened_name_length},
        {"/edgesets/left", &self->edgesets.num_records},
        {"/edgesets/children", &self->edgesets.total_children},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _dimension_read);
    size_t j;

    for (j = 0; j < num_fields; j++) {
        *fields[j].dest = 0;
        exists = H5Lexists(file_id, fields[j].name, H5P_DEFAULT);
        if (exists < 0) {
            goto out;
        }
        if (exists) {
            dataset_id = H5Dopen(file_id, fields[j].name, H5P_DEFAULT);
            if (dataset_id < 0) {
                ret = MSP_ERR_FILE_FORMAT;
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
            status = H5Sget_simple_extent_dims(dataspace_id, &dims, NULL);
            if (status < 0) {
                goto out;
            }
            *fields[j].dest = (size_t) dims;
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
    self->nodes.total_name_length = flattened_name_length + self->nodes.num_records;
    ret = tree_sequence_check_hdf5_dimensions(self, file_id);
    if (ret != 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_init_node_names(tree_sequence_t *self, size_t flattened_name_length,
        char *flattened_name)
{
    int ret = 0;
    size_t input_offset = 0;
    size_t output_offset = 0;
    size_t j, k;

    for (j = 0; j < self->nodes.num_records; j++) {
        for (k = 0; k < self->nodes.name_length[j]; k++) {
            assert(input_offset < flattened_name_length);
            assert(output_offset < self->nodes.total_name_length);
            self->nodes.name_mem[output_offset] = flattened_name[input_offset];
            input_offset++;
            output_offset++;
        }
        assert(output_offset < self->nodes.total_name_length);
        self->nodes.name_mem[output_offset] = '\0';
        output_offset++;
    }
    return ret;
}

static int
tree_sequence_read_hdf5_data(tree_sequence_t *self, hid_t file_id)
{
    herr_t status;
    int ret = MSP_ERR_HDF5;
    hid_t dataset_id;
    htri_t exists;
    struct _hdf5_field_read {
        const char *name;
        hid_t type;
        void *dest;
    };
    size_t flattened_name_length;
    char *flattened_name = NULL;
    struct _hdf5_field_read fields[] = {
        {"/provenance", 0, self->provenance_strings},
        {"/nodes/name", H5T_NATIVE_CHAR, NULL},
        {"/nodes/name_length", H5T_NATIVE_UINT32, self->nodes.name_length},
        {"/nodes/flags", H5T_NATIVE_UINT32, self->nodes.flags},
        {"/nodes/population", H5T_NATIVE_INT32, self->nodes.population},
        {"/nodes/time", H5T_NATIVE_DOUBLE, self->nodes.time},
        {"/mutations/nodes", H5T_NATIVE_INT32, self->mutations.nodes_mem},
        {"/mutations/num_nodes", H5T_NATIVE_INT32, self->mutations.num_nodes},
        {"/mutations/position", H5T_NATIVE_DOUBLE, self->mutations.position},
        {"/edgesets/left", H5T_NATIVE_DOUBLE, self->edgesets.left},
        {"/edgesets/right", H5T_NATIVE_DOUBLE, self->edgesets.right},
        {"/edgesets/parent", H5T_NATIVE_INT32, self->edgesets.parent},
        {"/edgesets/num_children", H5T_NATIVE_INT32, self->edgesets.num_children},
        {"/edgesets/children", H5T_NATIVE_INT32, self->edgesets.children_mem},
        {"/edgesets/indexes/insertion_order", H5T_NATIVE_INT32,
            self->edgesets.indexes.insertion_order},
        {"/edgesets/indexes/removal_order", H5T_NATIVE_INT32,
            self->edgesets.indexes.removal_order},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _hdf5_field_read);
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
    fields[0].type = vlen_str;

    flattened_name_length = self->nodes.total_name_length - self->nodes.num_records;
    flattened_name = malloc(flattened_name_length * sizeof(char));
    if (flattened_name == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    fields[1].dest = flattened_name;

    for (j = 0; j < num_fields; j++) {
        exists = H5Lexists(file_id, fields[j].name, H5P_DEFAULT);
        if (exists < 0) {
            goto out;
        }
        if (exists) {
            dataset_id = H5Dopen(file_id, fields[j].name, H5P_DEFAULT);
            if (dataset_id < 0) {
                goto out;
            }
            status = H5Dread(dataset_id, fields[j].type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    fields[j].dest);
            if (status < 0) {
                goto out;
            }
            status = H5Dclose(dataset_id);
            if (status < 0) {
                goto out;
            }
        }
    }
    status = H5Tclose(vlen_str);
    if (status < 0) {
        goto out;
    }
    ret = tree_sequence_init_node_names(self, flattened_name_length, flattened_name);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_init_nodes(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_init_edgesets(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_init_mutations(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_init_trees(self);
    if (ret != 0) {
        goto out;
    }
    ret = 0;
out:
    if (flattened_name != NULL) {
        free(flattened_name);
    }
    return ret;
}

int WARN_UNUSED
tree_sequence_load(tree_sequence_t *self, const char *filename, int flags)
{
    int ret = MSP_ERR_GENERIC;
    herr_t status;
    hid_t file_id = -1;

    if (self->initialised_magic != MSP_INITIALISED_MAGIC) {
        ret = MSP_ERR_NOT_INITIALISED;
        goto out;
    }
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        ret = MSP_ERR_HDF5;
        goto out;
    }
    ret = tree_sequence_read_hdf5_metadata(self, file_id);
    if (ret < 0) {
        goto out;
    }
    ret = tree_sequence_read_hdf5_groups(self, file_id);
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
    hsize_t dim, chunk_size;
    char *flattened_name = NULL;
    size_t flattened_name_length;
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
        {"/nodes/name",
            H5T_STD_I8LE, H5T_NATIVE_CHAR, 0, NULL},
        {"/nodes/name_length",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->nodes.num_records, self->nodes.name_length},
        {"/nodes/flags",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->nodes.num_records, self->nodes.flags},
        {"/nodes/population",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->nodes.num_records, self->nodes.population},
        {"/nodes/time",
            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
            self->nodes.num_records, self->nodes.time},
        {"/edgesets/left",
            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
            self->edgesets.num_records, self->edgesets.left},
        {"/edgesets/right",
            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
            self->edgesets.num_records, self->edgesets.right},
        {"/edgesets/parent",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->edgesets.num_records, self->edgesets.parent},
        {"/edgesets/num_children",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->edgesets.num_records, self->edgesets.num_children},
        {"/edgesets/children",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->edgesets.total_children, self->edgesets.children_mem},
        {"/edgesets/indexes/insertion_order",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->edgesets.num_records, self->edgesets.indexes.insertion_order},
        {"/edgesets/indexes/removal_order",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->edgesets.num_records, self->edgesets.indexes.removal_order},
        {"/mutations/nodes",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->mutations.total_nodes, self->mutations.nodes_mem},
        {"/mutations/num_nodes",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->mutations.num_records, self->mutations.num_nodes},
        {"/mutations/position",
            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
            self->mutations.num_records, self->mutations.position},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _hdf5_field_write);
    struct _hdf5_group_write {
        const char *name;
    };
    struct _hdf5_group_write groups[] = {
        {"/mutations"},
        {"/nodes"},
        {"/edgesets"},
        {"/edgesets/indexes"},
    };
    size_t num_groups = sizeof(groups) / sizeof(struct _hdf5_group_write);
    size_t j, k;

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

    assert(self->nodes.total_name_length >= self->nodes.num_records);
    /* Make the array to hold the flattened string */
    flattened_name_length = self->nodes.total_name_length - self->nodes.num_records;
    if (flattened_name_length != 0) {
        flattened_name = malloc(flattened_name_length * sizeof(char));
        if (flattened_name == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        /* fill in the array */
        k = 0;
        for (j = 0; j < self->nodes.total_name_length; j++) {
            if (self->nodes.name_mem[j] != '\0') {
                flattened_name[k] = self->nodes.name_mem[j];
                k++;
            }
        }
        assert(k == flattened_name_length);
        fields[1].size = flattened_name_length;
        fields[1].source = flattened_name;
    }

    /* Create the groups */
    for (j = 0; j < num_groups; j++) {
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
    /* now write the datasets */
    for (j = 0; j < num_fields; j++) {
        dim = fields[j].size;
        /* Never create any 0-sized datasets. This causes all sorts of problems in older
         * versions of HDF5, and so we adopt the protocol of omitting the dataset if it
         * is of zero size.
         */
        if (dim > 0) {
            dataspace_id = H5Screate_simple(1, &dim, &dim);
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
            chunk_size = GSL_MAX(1, fields[j].size);
            status = H5Pset_chunk(plist_id, 1, &chunk_size);
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
    if (flattened_name != NULL) {
        free(flattened_name);
    }
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
    ret = tree_sequence_write_hdf5_data(self, file_id, flags);
    if (ret < 0) {
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

/* Simple attribute getters */

double
tree_sequence_get_sequence_length(tree_sequence_t *self)
{
    return self->sequence_length;
}

size_t
tree_sequence_get_sample_size(tree_sequence_t *self)
{
    return self->sample_size;
}

size_t
tree_sequence_get_num_nodes(tree_sequence_t *self)
{
    return self->nodes.num_records;
}

size_t
tree_sequence_get_num_edgesets(tree_sequence_t *self)
{
    return self->edgesets.num_records;
}

size_t
tree_sequence_get_num_migrations(tree_sequence_t *self)
{
    return self->migrations.num_records;
}

size_t
tree_sequence_get_num_mutations(tree_sequence_t *self)
{
    return self->mutations.num_records;
}

size_t
tree_sequence_get_num_trees(tree_sequence_t *self)
{
    return self->num_trees;
}

/* Accessors for records */

int WARN_UNUSED
tree_sequence_get_pairwise_diversity(tree_sequence_t *self,
    node_id_t *samples, size_t num_samples, double *pi)
{
    int ret = 0;
    size_t j, k;
    node_id_t node;
    sparse_tree_t *tree = NULL;
    double result, denom, count, n;
    mutation_t *mut;

    if (num_samples < 2 || num_samples > self->sample_size) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    n = (double) num_samples;
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
    result = 0.0;
    for (ret = sparse_tree_first(tree); ret == 1; ret = sparse_tree_next(tree)) {
        for (j = 0; j < tree->num_mutations; j++) {
            mut = &tree->mutations[j];
            for (k = 0; k < mut->num_nodes; k++) {
                node = mut->nodes[k];
                count = (double) tree->num_tracked_leaves[node];
                result += count * (n - count);
            }
        }
    }
    if (ret != 0) {
        goto out;
    }
    denom = (n * (n - 1)) / 2.0;
    *pi = result / denom;
out:
    if (tree != NULL) {
        sparse_tree_free(tree);
        free(tree);
    }
    return ret;
}

int WARN_UNUSED
tree_sequence_get_node(tree_sequence_t *self, node_id_t index, node_t *node)
{
    int ret = 0;

    if (index < 0 || index >= (node_id_t) self->nodes.num_records) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    node->time = self->nodes.time[index];
    node->population = self->nodes.population[index];
    node->flags = self->nodes.flags[index];
    node->name = self->nodes.name[index];
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_get_edgeset(tree_sequence_t *self, size_t index, edgeset_t *edgeset)
{
    int ret = 0;

    if (index >= self->edgesets.num_records) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    edgeset->left = self->edgesets.left[index];
    edgeset->right = self->edgesets.right[index];
    edgeset->parent = self->edgesets.parent[index];
    edgeset->num_children = (size_t) self->edgesets.num_children[index];
    edgeset->children = self->edgesets.children[index];
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_get_migration(tree_sequence_t *self, size_t index, migration_t *record)
{
    int ret = 0;

    if (index >= self->migrations.num_records) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    record->node = self->migrations.node[index];
    record->source = self->migrations.source[index];
    record->dest = self->migrations.dest[index];
    record->left = self->migrations.left[index];
    record->right = self->migrations.right[index];
    record->time = self->migrations.time[index];
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_get_mutation(tree_sequence_t *self, size_t index, mutation_t *record)
{
    int ret = 0;

    if (index >= self->mutations.num_records) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    /* TODO is this redundant?? */
    record->index = index;
    record->position = self->mutations.position[index];
    record->num_nodes = (size_t) self->mutations.num_nodes[index];
    record->nodes = self->mutations.nodes[index];
    record->ancestral_state = self->mutations.ancestral_state[index];
    record->derived_state = self->mutations.derived_state[index];
out:
    return ret;
}


/* Compress the node space in the specified set of records and mutations.
 */
static int WARN_UNUSED
tree_sequence_compress_nodes(tree_sequence_t *self, node_id_t *samples, size_t num_samples,
        coalescence_record_t *records, size_t num_records, mutation_t *mutations,
        size_t num_mutations)
{
    int ret = MSP_ERR_GENERIC;
    node_id_t *node_map = NULL;
    node_id_t next_node;
    size_t c, j, k;
    coalescence_record_t *cr;

    node_map = malloc(self->nodes.num_records * sizeof(node_id_t));
    if (node_map == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (j = 0; j < self->nodes.num_records; j++) {
        node_map[j] = MSP_NULL_NODE;
    }
    for (j = 0; j < num_samples; j++) {
        node_map[samples[j]] = (node_id_t) j;
    }
    next_node = (node_id_t) num_samples;
    for (j = 0; j < num_records; j++) {
        cr = &records[j];
        if (node_map[cr->node] == MSP_NULL_NODE) {
            node_map[cr->node] = next_node;
            next_node++;
        }
        cr->node = node_map[cr->node];
        for (c = 0; c < cr->num_children; c++) {
            cr->children[c] = node_map[cr->children[c]];
        }
        qsort(cr->children, cr->num_children, sizeof(node_id_t), cmp_node_id_t);
    }
    for (j = 0; j < num_mutations; j++) {
        for (k = 0; k < mutations[j].num_nodes; k++) {
            mutations[j].nodes[k] = node_map[mutations[j].nodes[k]];
            assert(mutations[j].nodes[k] != MSP_NULL_NODE);
        }
        qsort(mutations[j].nodes, mutations[j].num_nodes, sizeof(node_id_t), cmp_node_id_t);
    }
    ret = 0;
out:
    if (node_map != NULL) {
        free(node_map);
    }
    return ret;
}

/* TODO this needs to be updated to use the new tables/edgesets API. We currently
 * use coalescence_records because it makes it simpler for sorting records by time.
 * This should really be spun into its own class, as this function is far too long.
 */
int WARN_UNUSED
tree_sequence_simplify(tree_sequence_t *self, node_id_t *samples,
        size_t num_samples, int flags, tree_sequence_t *output)
{
    typedef struct {
        bool active;
        double left;
        node_id_t *mapped_children;
        uint32_t num_mapped_children;
    } active_record_t;

    int ret = MSP_ERR_GENERIC;
    node_id_t *parent = NULL;
    node_id_t *num_children = NULL;
    node_id_t **children = NULL;
    node_id_t *mapping = NULL;
    node_id_t *mapped_children = NULL;
    node_id_t *mapped_children_mem = NULL;
    node_id_t *output_mutations_nodes_mem = NULL;
    sample_t *sample_objects = NULL;
    active_record_t *active_records = NULL;
    coalescence_record_t *output_records = NULL;
    mutation_t *output_mutations = NULL;
    node_id_t *I = self->edgesets.indexes.insertion_order;
    node_id_t *O = self->edgesets.indexes.removal_order;
    size_t M = self->edgesets.num_records;
    node_id_t h;
    size_t j, k, next_avl_node, mapped_children_mem_offset, num_output_records,
           num_output_mutations, max_num_child_nodes, max_num_records,
           output_mutations_mem_offset;
    node_id_t u, v, w, c;
    size_t l, node_index, num_mapped_children;
    avl_tree_t visited_nodes;
    avl_node_t *avl_node_mem = NULL;
    node_id_t *avl_node_value_mem = NULL;
    avl_node_t *avl_node;
    active_record_t *ar;
    coalescence_record_t *cr;
    mutation_t *mut;
    bool equal, activate_record, keep;
    double right, x;
    bool filter_root_mutations = flags & MSP_FILTER_ROOT_MUTATIONS;

    if (num_samples < 2) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    parent = malloc(self->nodes.num_records * sizeof(node_id_t));
    children = malloc(self->nodes.num_records * sizeof(node_id_t *));
    num_children = malloc(self->nodes.num_records * sizeof(node_id_t));
    mapping = malloc(self->nodes.num_records * sizeof(node_id_t));
    sample_objects = malloc(num_samples * sizeof(sample_t));
    avl_node_mem = malloc(self->nodes.num_records * sizeof(avl_node_t));
    avl_node_value_mem = malloc(self->nodes.num_records * sizeof(node_id_t));
    active_records = malloc(self->nodes.num_records * sizeof(active_record_t));
    mapped_children = malloc(self->nodes.num_records * sizeof(node_id_t));
    /* TODO work out better bounds for these values */
    max_num_child_nodes = 2 * self->edgesets.total_children;
    max_num_records = 2 * self->edgesets.num_records;
    mapped_children_mem = malloc(max_num_child_nodes * sizeof(node_id_t));
    output_records = malloc(max_num_records * sizeof(coalescence_record_t));
    output_mutations = malloc(self->mutations.num_records * sizeof(mutation_t));
    output_mutations_nodes_mem = malloc(self->mutations.total_nodes *
            sizeof(node_id_t));
    if (parent == NULL || children == NULL || num_children == NULL
            || mapping == NULL || sample_objects == NULL
            || avl_node_mem == NULL || avl_node_value_mem == NULL
            || mapped_children == NULL || active_records == NULL
            || mapped_children_mem == NULL || output_records == NULL
            || output_mutations == NULL || output_mutations_nodes_mem == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    /* Initialise the mapping and tree structures */
    for (u = 0; u < (node_id_t) self->nodes.num_records; u++) {
        parent[u] = MSP_NULL_NODE;
        children[u] = NULL;
        num_children[u] = 0;
        mapping[u] = MSP_NULL_NODE;
        avl_node_mem[u].item = avl_node_value_mem + u;
        active_records[u].active = false;
    }
    for (j = 0; j < num_samples; j++) {
        u = samples[j];
        if (u < 0 || u >= (node_id_t) self->sample_size) {
            ret = MSP_ERR_BAD_SAMPLES;
            goto out;
        }
        if (mapping[u] != MSP_NULL_NODE) {
            ret = MSP_ERR_DUPLICATE_SAMPLE;
            goto out;
        }
        mapping[u] = u;
        sample_objects[j].population_id = self->nodes.population[u];
        sample_objects[j].time = self->nodes.time[u];
    }
    avl_init_tree(&visited_nodes, cmp_node_id_t, NULL);
    mapped_children_mem_offset = 0;
    output_mutations_mem_offset = 0;
    num_output_records = 0;
    num_output_mutations = 0;

    j = 0;
    k = 0;
    l = 0;
    while (j < M) {
        x = self->edgesets.left[I[j]];
        next_avl_node = 0;
        avl_clear_tree(&visited_nodes);

        /* Records out */
        while (self->edgesets.right[O[k]] == x) {
            h = O[k];
            k++;
            u = self->edgesets.parent[h];
            for (c = 0; c < num_children[u]; c++) {
                parent[children[u][c]] = MSP_NULL_NODE;
            }
            num_children[u] = 0;
            children[u] = NULL;
            /* Propagate up to the root and save visited nodes */
            while (u != MSP_NULL_NODE) {
                if (avl_search(&visited_nodes, &u) == NULL) {
                    assert(next_avl_node < self->nodes.num_records);
                    avl_node = &avl_node_mem[next_avl_node];
                    next_avl_node++;
                    *((node_id_t *) avl_node->item) = u;
                    avl_node = avl_insert_node(&visited_nodes, avl_node);
                    assert(avl_node != NULL);
                }

                w = MSP_NULL_NODE;
                for (c = 0; c < num_children[u]; c++) {
                    v = children[u][c];
                    if (mapping[v] != MSP_NULL_NODE) {
                        w = w == MSP_NULL_NODE ? mapping[v]: u;
                    }
                }
                mapping[u] = w;
                u = parent[u];
            }
        }

        /* Records in */
        while (j < M && self->edgesets.left[I[j]] == x) {
            h = I[j];
            j++;
            u = self->edgesets.parent[h];
            num_children[u] = self->edgesets.num_children[h];
            children[u] = self->edgesets.children[h];
            for (c = 0; c < num_children[u]; c++) {
                v = children[u][c];
                parent[v] = u;
            }
            /* Propagate up to the root and save visited nodes */
            while (u != MSP_NULL_NODE) {
                if (avl_search(&visited_nodes, &u) == NULL) {
                    assert(next_avl_node < self->nodes.num_records);
                    avl_node = &avl_node_mem[next_avl_node];
                    next_avl_node++;
                    *((node_id_t *) avl_node->item) = u;
                    avl_node = avl_insert_node(&visited_nodes, avl_node);
                    assert(avl_node != NULL);
                }

                w = MSP_NULL_NODE;
                for (c = 0; c < num_children[u]; c++) {
                    v = children[u][c];
                    if (mapping[v] != MSP_NULL_NODE) {
                        w = w == MSP_NULL_NODE ? mapping[v]: u;
                    }
                }
                mapping[u] = w;
                u = parent[u];
            }
        }

        /* Examine the visited nodes and update the active records */
        for (avl_node = visited_nodes.head; avl_node != NULL;
                avl_node = avl_node->next) {
            u = *((node_id_t *) avl_node->item);
            ar = &active_records[u];
            activate_record = false;
            if (ar->active) {
                /* Compare the mapped children at this node to the record. */
                num_mapped_children = 0;
                for (c = 0; c < num_children[u]; c++) {
                    v = children[u][c];
                    if (mapping[v] != MSP_NULL_NODE) {
                        assert(num_mapped_children < self->nodes.num_records);
                        mapped_children[num_mapped_children] = mapping[v];
                        num_mapped_children++;
                    }
                }
                equal = false;
                if (num_mapped_children == ar->num_mapped_children) {
                    qsort(mapped_children, num_mapped_children, sizeof(node_id_t),
                            cmp_node_id_t);
                    equal = memcmp(ar->mapped_children, mapped_children,
                            num_mapped_children * sizeof(node_id_t)) == 0;
                }
                if (!equal) {
                    ar->active = false;
                    assert(num_output_records < max_num_records);
                    cr = &output_records[num_output_records];
                    num_output_records++;
                    cr->left = ar->left;
                    cr->right = x;
                    cr->node = u;
                    cr->num_children = ar->num_mapped_children;
                    cr->children = ar->mapped_children;
                    cr->time = self->nodes.time[u];
                    cr->population_id = self->nodes.population[u];
                    if (u == mapping[u]) {
                        activate_record = true;
                    }
                }
            } else {
                if (u == mapping[u]) {
                    activate_record = true;
                }
            }
            if (activate_record) {
                ar->active = true;
                ar->left = x;
                ar->num_mapped_children = 0;
                ar->mapped_children = mapped_children_mem + mapped_children_mem_offset;
                for (c = 0; c < num_children[u]; c++) {
                    v = children[u][c];
                    if (mapping[v] != MSP_NULL_NODE) {
                        assert(mapped_children_mem_offset < max_num_child_nodes);
                        mapped_children_mem_offset++;
                        ar->mapped_children[ar->num_mapped_children] = mapping[v];
                        ar->num_mapped_children++;
                    }
                }
                qsort(ar->mapped_children, ar->num_mapped_children, sizeof(node_id_t),
                        cmp_node_id_t);
            }
        }

        /* Update the mutations for this tree */
        right = self->edgesets.right[O[k]];
        while (l < self->mutations.num_records && self->mutations.position[l] < right) {
            node_index = 0;
            mut = NULL;
            for (h = 0; h < self->mutations.num_nodes[l]; h++) {
                u = self->mutations.nodes[l][h];
                if (mapping[u] != MSP_NULL_NODE) {
                    keep = true;
                    if (filter_root_mutations) {
                        /* Traverse up the tree until we find either another node in
                         * the subset tree or the root */
                        v = parent[u];
                        while (v != MSP_NULL_NODE && mapping[v] != v) {
                            v = parent[v];
                        }
                        keep = v != MSP_NULL_NODE;
                    }
                    if (keep) {
                        if (node_index == 0) {
                            assert(num_output_mutations < self->mutations.num_records);
                            mut = &output_mutations[num_output_mutations];
                            num_output_mutations++;
                            assert(output_mutations_mem_offset < self->mutations.total_nodes);
                            mut->nodes = output_mutations_nodes_mem + output_mutations_mem_offset;
                        }
                        mut->nodes[node_index] = mapping[u];
                        node_index++;
                        output_mutations_mem_offset++;
                        mut->num_nodes = node_index;
                        mut->position = self->mutations.position[l];
                        mut->ancestral_state = self->mutations.ancestral_state[l];
                        mut->derived_state = self->mutations.derived_state[l];
                    }
                }
            }
            l++;
        }
    }

    /* After the main loop has completed, find all the records that have not
     * been finished and terminate them.
     */
    x = self->sequence_length;
    for (u = 0; u < (node_id_t) self->nodes.num_records; u++) {
        ar = &active_records[u];
        if (ar->active) {
            assert(num_output_records < max_num_records);
            cr = &output_records[num_output_records];
            num_output_records++;
            cr->left = ar->left;
            cr->right = x;
            cr->node = u;
            cr->time = self->nodes.time[u];
            cr->population_id = self->nodes.population[u];
            cr->num_children = (uint32_t) ar->num_mapped_children;
            cr->children = ar->mapped_children;
        }
    }

    if (num_output_records == 0) {
        ret = MSP_ERR_CANNOT_SIMPLIFY;
        goto out;
    }
    /* Sort the records by time and left coordinate */
    qsort(output_records, num_output_records, sizeof(coalescence_record_t),
            cmp_record_time_left);
    ret = tree_sequence_compress_nodes(self, samples, num_samples,
            output_records, num_output_records, output_mutations,
            num_output_mutations);
    if (ret != 0) {
        goto out;
    }
    /* Alloc a new tree sequence for these records. */
    ret = tree_sequence_load_records(output, num_samples, sample_objects,
            num_output_records, output_records,
            num_output_mutations, output_mutations);
    if (ret != 0) {
        tree_sequence_free(output);
        goto out;
    }
out:
    if (parent != NULL) {
        free(parent);
    }
    if (children != NULL) {
        free(children);
    }
    if (num_children != NULL) {
        free(num_children);
    }
    if (mapping != NULL) {
        free(mapping);
    }
    if (sample_objects != NULL) {
        free(sample_objects);
    }
    if (avl_node_value_mem != NULL) {
        free(avl_node_value_mem);
    }
    if (avl_node_mem != NULL) {
        free(avl_node_mem);
    }
    if (active_records != NULL) {
        free(active_records);
    }
    if (mapped_children != NULL) {
        free(mapped_children);
    }
    if (mapped_children_mem != NULL) {
        free(mapped_children_mem);
    }
    if (output_records != NULL) {
        free(output_records);
    }
    if (output_mutations != NULL) {
        free(output_mutations);
    }
    if (output_mutations_nodes_mem != NULL) {
        free(output_mutations_nodes_mem);
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
    self->num_records = tree_sequence_get_num_edgesets(tree_sequence);
    self->tree_sequence = tree_sequence;
    self->insertion_index = 0;
    self->removal_index = 0;
    self->tree_left = 0;
    self->tree_index = (size_t) -1;
    self->node_records = malloc(self->num_nodes * sizeof(node_record_t));
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
    fprintf(out, "tree_left = %f\n", self->tree_left);
    fprintf(out, "tree_index = %d\n", (int) self->tree_index);
}

int WARN_UNUSED
tree_diff_iterator_next(tree_diff_iterator_t *self, double *length,
        node_record_t **nodes_out, node_record_t **nodes_in)
{
    int ret = 0;
    node_id_t k;
    double last_left = self->tree_left;
    size_t next_node_record = 0;
    tree_sequence_t *s = self->tree_sequence;
    node_record_t *out_head = NULL;
    node_record_t *out_tail = NULL;
    node_record_t *in_head = NULL;
    node_record_t *in_tail = NULL;
    node_record_t *w = NULL;
    size_t num_trees = tree_sequence_get_num_trees(s);

    assert(s != NULL);

    if (self->tree_index + 1 < num_trees) {
        /* First we remove the stale records */
        while (s->edgesets.right[
                s->edgesets.indexes.removal_order[self->removal_index]]
                    == self->tree_left) {
            k = s->edgesets.indexes.removal_order[self->removal_index];
            assert(next_node_record < self->num_nodes);
            w = &self->node_records[next_node_record];
            next_node_record++;
            w->node = s->edgesets.parent[k];
            w->time = s->nodes.time[w->node];
            w->num_children = (size_t) s->edgesets.num_children[k];
            w->children = s->edgesets.children[k];
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
                s->edgesets.left[
                    s->edgesets.indexes.insertion_order[self->insertion_index]]
                        == self->tree_left) {
            k = s->edgesets.indexes.insertion_order[self->insertion_index];
            assert(next_node_record < self->num_nodes);
            w = &self->node_records[next_node_record];
            next_node_record++;
            w->node = s->edgesets.parent[k];
            w->time = s->nodes.time[w->node];
            w->num_children = (size_t) s->edgesets.num_children[k];
            w->children = s->edgesets.children[k];
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
        self->tree_left = s->edgesets.right[
            s->edgesets.indexes.removal_order[self->removal_index]];
        self->tree_index++;
        ret = 1;
    }
    *nodes_out = out_head;
    *nodes_in = in_head;
    *length = 0;
    if (num_trees > 0) {
        *length = self->tree_left - last_left;
    }
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
    self->right = 0;
    self->root = 0;
    self->index = (size_t) -1;
    memset(self->parent, 0xff, N * sizeof(node_id_t));
    memset(self->population + n, 0xff, (N - n) * sizeof(population_id_t));
    memset(self->time + n, 0, (N - n) * sizeof(double));
    memset(self->num_children + n, 0, (N - n) * sizeof(node_id_t));
    memset(self->children + n, 0, (N - n) * sizeof(node_id_t *));
    if (self->flags & MSP_LEAF_COUNTS) {
        memset(self->num_leaves + n, 0, (N - n) * sizeof(node_id_t));
        memset(self->num_tracked_leaves + n, 0, (N - n) * sizeof(node_id_t));
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
    size_t j, sample_size;
    size_t num_nodes;
    leaf_list_node_t *w;

    memset(self, 0, sizeof(sparse_tree_t));
    if (tree_sequence == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    num_nodes = tree_sequence->nodes.num_records;
    sample_size = tree_sequence->sample_size;
    self->num_nodes = num_nodes;
    self->sample_size = sample_size;
    self->tree_sequence = tree_sequence;
    self->flags = flags;
    self->parent = malloc(num_nodes * sizeof(node_id_t));
    self->population = malloc(num_nodes * sizeof(population_id_t));
    self->time = malloc(num_nodes * sizeof(double));
    self->num_children = malloc(num_nodes * sizeof(node_id_t));
    self->children = malloc(num_nodes * sizeof(node_id_t *));
    if (self->time == NULL || self->parent == NULL || self->children == NULL
            || self->num_children == NULL || self->population == NULL) {
        goto out;
    }
    /* the maximum possible height of the tree is num_nodes + 1, including
     * the null value. */
    self->stack1 = malloc((num_nodes + 1) * sizeof(node_id_t));
    self->stack2 = malloc((num_nodes + 1) * sizeof(node_id_t));
    if (self->stack1 == NULL || self->stack2 == NULL) {
        goto out;
    }
    if (self->flags & MSP_LEAF_COUNTS) {
        self->num_leaves = calloc(num_nodes, sizeof(node_id_t));
        self->num_tracked_leaves = calloc(num_nodes, sizeof(node_id_t));
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
            w->node = (node_id_t) j;
            self->leaf_list_head[j] = w;
            self->leaf_list_tail[j] = w;
        }
    }
    /* Set the sample attributes */
    for (j = 0; j < self->sample_size; j++) {
        self->population[j] = self->tree_sequence->nodes.population[j];
        self->time[j] = self->tree_sequence->nodes.time[j];
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
    memset(self->num_tracked_leaves, 0, self->num_nodes * sizeof(node_id_t));
out:
    return ret;
}


int WARN_UNUSED
sparse_tree_set_tracked_leaves(sparse_tree_t *self, size_t num_tracked_leaves,
        node_id_t *tracked_leaves)
{
    int ret = MSP_ERR_GENERIC;
    size_t j;
    node_id_t u;

    /* TODO This is not needed when the sparse tree is new. We should use the
     * state machine to check and only reset the tracked leaves when needed.
     */
    ret = sparse_tree_reset_tracked_leaves(self);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < num_tracked_leaves; j++) {
        u = tracked_leaves[j];
        if (u < 0 || u >= (node_id_t) self->sample_size) {
            ret = MSP_ERR_OUT_OF_BOUNDS;
            goto out;
        }
        if (self->num_tracked_leaves[u] != 0) {
            ret = MSP_ERR_DUPLICATE_SAMPLE;
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
    node_id_t u;
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
    self->right = source->right;
    self->root = source->root;
    self->index = source->index;
    self->num_mutations= source->num_mutations;
    self->mutations = source->mutations;

    memcpy(self->parent, source->parent, N * sizeof(node_id_t));
    memcpy(self->population, source->population, N * sizeof(population_id_t));
    memcpy(self->time, source->time, N * sizeof(double));
    memcpy(self->num_children, source->num_children, N * sizeof(node_id_t));
    memcpy(self->children, source->children, N * sizeof(node_id_t *));
    if (self->flags & MSP_LEAF_COUNTS) {
        if (! (source->flags & MSP_LEAF_COUNTS)) {
            ret = MSP_ERR_UNSUPPORTED_OPERATION;
            goto out;
        }
        memcpy(self->num_leaves + n, source->num_leaves + n,
                (N - n) * sizeof(node_id_t));
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
        && self->left == other->left
        && self->right == other->right
        && self->root == other->root
        && self->num_mutations == other->num_mutations
        && self->mutations == other->mutations
        && memcmp(self->parent, other->parent, N * sizeof(node_id_t)) == 0
        && memcmp(self->population, other->population,
                N * sizeof(population_id_t)) == 0
        && memcmp(self->time, other->time, N * sizeof(double)) ==  0
        && memcmp(self->num_children, other->num_children,
                N * sizeof(node_id_t)) == 0
        && memcmp(self->children, other->children,
                N * sizeof(node_id_t *)) == 0;
    if (condition) {
        ret = 0;
    }
out:
    return ret;
}

static int
sparse_tree_check_node(sparse_tree_t *self, node_id_t u)
{
    int ret = 0;
    if (u < 0 || u >= (node_id_t) self->num_nodes) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
    }
    return ret;
}


int WARN_UNUSED
sparse_tree_get_mrca(sparse_tree_t *self, node_id_t u, node_id_t v,
        node_id_t *mrca)
{
    int ret = 0;
    node_id_t w = 0;
    node_id_t *s1 = self->stack1;
    node_id_t *s2 = self->stack2;
    node_id_t j;
    int l1, l2;

    ret = sparse_tree_check_node(self, u);
    if (ret != 0) {
        goto out;
    }
    ret = sparse_tree_check_node(self, v);
    if (ret != 0) {
        goto out;
    }
    j = u;
    l1 = 0;
    while (j != MSP_NULL_NODE) {
        assert(l1 < (int) self->num_nodes);
        s1[l1] = j;
        l1++;
        j = self->parent[j];
    }
    s1[l1] = MSP_NULL_NODE;
    j = v;
    l2 = 0;
    while (j != MSP_NULL_NODE) {
        assert(l2 < (int) self->num_nodes);
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
sparse_tree_get_num_leaves_by_traversal(sparse_tree_t *self, node_id_t u,
        size_t *num_leaves)
{
    int ret = 0;
    node_id_t *stack = self->stack1;
    node_id_t v, c;
    size_t count = 0;
    int stack_top = 0;

    stack[0] = u;
    while (stack_top >= 0) {
        v = stack[stack_top];
        stack_top--;
        if (v < (node_id_t) self->sample_size) {
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
sparse_tree_get_num_leaves(sparse_tree_t *self, node_id_t u, size_t *num_leaves)
{
    int ret = 0;

    ret = sparse_tree_check_node(self, u);
    if (ret != 0) {
        goto out;
    }

    if (self->flags & MSP_LEAF_COUNTS) {
        *num_leaves = (size_t) self->num_leaves[u];
    } else {
        ret = sparse_tree_get_num_leaves_by_traversal(self, u, num_leaves);
    }
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_get_num_tracked_leaves(sparse_tree_t *self, node_id_t u,
        size_t *num_tracked_leaves)
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
    *num_tracked_leaves = (size_t) self->num_tracked_leaves[u];
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_get_leaf_list(sparse_tree_t *self, node_id_t u,
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
    *head = self->leaf_list_head[u];
    *tail = self->leaf_list_tail[u];
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_get_root(sparse_tree_t *self, node_id_t *root)
{
    *root = self->root;
    return 0;
}


int WARN_UNUSED
sparse_tree_get_parent(sparse_tree_t *self, node_id_t u, node_id_t *parent)
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
sparse_tree_get_time(sparse_tree_t *self, node_id_t u, double *t)
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
sparse_tree_get_children(sparse_tree_t *self, node_id_t u,
        size_t *num_children, node_id_t **children)
{
    int ret = 0;

    ret = sparse_tree_check_node(self, u);
    if (ret != 0) {
        goto out;
    }
    *num_children = (size_t) self->num_children[u];
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
    node_id_t u, v;
    size_t j, k, num_leaves;
    int err, found;

    for (j = 0; j < self->sample_size; j++) {
        u = (node_id_t) j;
        assert(self->time[u] >= 0.0);
        assert(self->num_children[j] == 0);
        while (self->parent[u] != MSP_NULL_NODE) {
            v = self->parent[u];
            found = 0;
            for (k = 0; k < (size_t) self->num_children[v]; k++) {
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
        for (u = 0; u < (node_id_t) self->num_nodes; u++) {
            err = sparse_tree_get_num_leaves_by_traversal(self, u, &num_leaves);
            assert(err == 0);
            assert(num_leaves == (size_t) self->num_leaves[u]);
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
    fprintf(out, "right = %f\n", self->right);
    fprintf(out, "root = %d\n", (int) self->root);
    fprintf(out, "index = %d\n", (int) self->index);
    for (j = 0; j < self->num_nodes; j++) {
        fprintf(out, "\t%d\t%d\t%f\t%d\t(", (int) j, (int) self->parent[j],
            self->time[j], (int) self->population[j]);
        for (k = 0; k < (size_t) self->num_children[j]; k++) {
            fprintf(out, "%d", (int) self->children[j][k]);
            if (k < (size_t) self->num_children[j] - 1) {
                fprintf(out, ", ");
            }
        }
        fprintf(out, ")");
        if (self->flags & MSP_LEAF_COUNTS) {
            fprintf(out, "\t%d\t%d\t%d", (int) self->num_leaves[j],
                    (int) self->num_tracked_leaves[j], self->marked[j]);
        }
        if (self->flags & MSP_LEAF_LISTS) {
            fprintf(out, "\t[");
            u = self->leaf_list_head[j];
            if (u != NULL) {
                while (1) {
                    fprintf(out, "%d ", (int) u->node);
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
        fprintf(out, "\t%f\t", self->mutations[j].position);
        for (k = 0; k < self->mutations[j].num_nodes; k++) {
            fprintf(out, "%d,", (int) self->mutations[j].nodes[k]);
        }
        fprintf(out, "\n");
    }
    sparse_tree_check_state(self);
}

/* Methods for positioning the tree along the sequence */

static inline void
sparse_tree_propagate_leaf_count_loss(sparse_tree_t *self, node_id_t u)
{
    const node_id_t all_leaves_diff = self->num_leaves[u];
    const node_id_t tracked_leaves_diff = self->num_tracked_leaves[u];
    const uint8_t mark = self->mark;
    node_id_t v = u;

    /* propagate this loss up as far as we can */
    while (v != MSP_NULL_NODE) {
        self->num_leaves[v] -= all_leaves_diff;
        self->num_tracked_leaves[v] -= tracked_leaves_diff;
        self->marked[v] = mark;
        v = self->parent[v];
    }
}

static inline void
sparse_tree_propagate_leaf_count_gain(sparse_tree_t *self, node_id_t u)
{
    node_id_t j, k, v, *c;
    node_id_t all_leaves_diff = 0;
    node_id_t tracked_leaves_diff = 0;
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
sparse_tree_update_leaf_lists(sparse_tree_t *self, node_id_t node)
{
    node_id_t u, v, c;
    leaf_list_node_t **head = self->leaf_list_head;
    leaf_list_node_t **tail = self->leaf_list_tail;

    u = node;
    while (u != MSP_NULL_NODE) {
        head[u] = NULL;
        tail[u] = NULL;
        for (c = 0; c < self->num_children[u]; c++) {
            v = self->children[u][c];
            if (head[v] != NULL) {
                assert(tail[v] != NULL);
                if (head[u] == NULL) {
                    head[u] = head[v];
                    tail[u] = tail[v];
                } else {
                    tail[u]->next = head[v];
                    tail[u] = tail[v];
                }
            }
        }
        u = self->parent[u];
    }
}

static int
sparse_tree_advance(sparse_tree_t *self, int direction,
        double *out_breakpoints, node_id_t *out_order, node_id_t *out_index,
        double *in_breakpoints, node_id_t *in_order, node_id_t *in_index,
        int first_tree)
{
    int ret = 0;
    int direction_change = direction * (direction != self->direction);
    node_id_t in = *in_index + direction_change;
    node_id_t out = *out_index + direction_change;
    node_id_t j, k, u, oldest_child;
    double x = in_breakpoints[in_order[in]];
    double oldest_child_time;
    tree_sequence_t *s = self->tree_sequence;
    node_id_t R = (node_id_t) s->edgesets.num_records;

    while (out_breakpoints[out_order[out]] == x) {
        k = out_order[out];
        u = s->edgesets.parent[k];
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
            sparse_tree_update_leaf_lists(self, u);
        }
        out += direction;
    }

    while (in >= 0 && in < R && in_breakpoints[in_order[in]] == x) {
        k = in_order[in];
        u = s->edgesets.parent[k];
        for (j = 0; j < s->edgesets.num_children[k]; j++) {
            self->parent[s->edgesets.children[k][j]] = u;
        }
        self->num_children[u] = s->edgesets.num_children[k];
        self->children[u] = s->edgesets.children[k];
        self->time[u] = s->nodes.time[u];
        self->population[u] = s->nodes.population[u];
        if (self->time[u] > self->time[self->root]) {
            self->root = u;
        }
        if (self->flags & MSP_LEAF_COUNTS) {
            sparse_tree_propagate_leaf_count_gain(self, u);
        }
        if (self->flags & MSP_LEAF_LISTS) {
            sparse_tree_update_leaf_lists(self, u);
        }
        in += direction;
    }
    /* In very rare situations, we have to traverse upwards to find the
     * new root.
     */
    while (self->parent[self->root] != MSP_NULL_NODE) {
        self->root = self->parent[self->root];
    }

    if (direction == MSP_DIR_FORWARD) {
        self->left = x;
        self->right = out_breakpoints[out_order[out]];
    } else {
        self->left = out_breakpoints[out_order[out]];
        self->right = x;
    }
    self->direction = direction;
    self->index = (size_t) ((int64_t) self->index + direction);
    *out_index = out;
    *in_index = in;
    if (s->mutations.num_records > 0) {
        self->mutations = s->mutations.tree_mutations[self->index];
        self->num_mutations = (size_t) s->mutations.num_tree_mutations[self->index];
    }

    ret = 1;
    return ret;
}


int WARN_UNUSED
sparse_tree_first(sparse_tree_t *self)
{
    int ret = 0;
    tree_sequence_t *s = self->tree_sequence;

    if (s->edgesets.num_records > 0) {
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
                s->edgesets.right, s->edgesets.indexes.removal_order,
                &self->right_index, s->edgesets.left,
                s->edgesets.indexes.insertion_order, &self->left_index, 1);
    }
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_last(sparse_tree_t *self)
{
    int ret = 0;
    tree_sequence_t *s = self->tree_sequence;

    if (s->edgesets.num_records > 0) {
        /* TODO this is redundant if this is the first usage of the tree. We
         * should add a state machine here so we know what state the tree is
         * in and can take the appropriate actions.
         */
        ret = sparse_tree_clear(self);
        if (ret != 0) {
            goto out;
        }
        self->left_index = (node_id_t) s->edgesets.num_records - 1;
        self->right_index = (node_id_t) s->edgesets.num_records - 1;
        self->direction = MSP_DIR_REVERSE;
        self->index = tree_sequence_get_num_trees(s);

        ret = sparse_tree_advance(self, MSP_DIR_REVERSE,
                s->edgesets.left, s->edgesets.indexes.insertion_order,
                &self->left_index, s->edgesets.right,
                s->edgesets.indexes.removal_order, &self->right_index, 1);
    }
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
                s->edgesets.right, s->edgesets.indexes.removal_order,
                &self->right_index, s->edgesets.left,
                s->edgesets.indexes.insertion_order, &self->left_index, 0);
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
                s->edgesets.left, s->edgesets.indexes.insertion_order,
                &self->left_index, s->edgesets.right,
                s->edgesets.indexes.removal_order, &self->right_index, 0);
    }
    return ret;
}
