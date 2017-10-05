/*
** Copyright (C) 2015-2017 University of Oxford
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
#include "object_heap.h"
#include "msprime.h"


#define MSP_DIR_FORWARD 1
#define MSP_DIR_REVERSE -1

typedef struct {
    node_id_t index;
    /* These are the sort keys in order */
    double first;
    double second;
    node_id_t third;
    node_id_t fourth;
} index_sort_t;

static int
cmp_index_sort(const void *a, const void *b) {
    const index_sort_t *ca = (const index_sort_t *) a;
    const index_sort_t *cb = (const index_sort_t *) b;
    int ret = (ca->first > cb->first) - (ca->first < cb->first);
    if (ret == 0) {
        ret = (ca->second > cb->second) - (ca->second < cb->second);
        if (ret == 0) {
            ret = (ca->third > cb->third) - (ca->third < cb->third);
            if (ret == 0) {
                ret = (ca->fourth > cb->fourth) - (ca->fourth < cb->fourth);
            }
        }
    }
    return ret;
}

/* Initialise the strings represented by the flattened source array
 * and list of lengths into the specified memory and pointers.
 */
static int WARN_UNUSED
init_string_column(size_t num_rows, char *source, uint32_t *length,
        char **pointers, char *mem)
{
    int ret = 0;
    size_t j, k, mem_offset, source_offset;

    mem_offset = 0;
    source_offset = 0;
    for (j = 0; j < num_rows; j++) {
        pointers[j] = mem + mem_offset;
        for (k = 0; k < length[j]; k++) {
            mem[mem_offset] = source[source_offset];
            source_offset++;
            mem_offset++;
        }
        mem[mem_offset] = '\0';
        mem_offset++;
    }
    return ret;
}

static int WARN_UNUSED
flatten_string_column(size_t total_length, char *mem, char *flattened)
{
    int ret = 0;
    size_t j, k;

    /* fill in the array */
    k = 0;
    for (j = 0; j < total_length; j++) {
        if (mem[j] != '\0') {
            flattened[k] = mem[j];
            k++;
        }
    }
    return ret;
}

static int WARN_UNUSED
validate_length(size_t num_rows, uint32_t *length, size_t total_length)
{
    int ret = 0;
    size_t j;
    size_t sum = 0;

    for (j = 0; j < num_rows; j++) {
        sum += length[j];
    }
    if (sum != total_length) {
        ret = MSP_ERR_LENGTH_MISMATCH;
    }
    return ret;

}

/* ======================================================== *
 * tree sequence
 * ======================================================== */

static void
tree_sequence_check_state(tree_sequence_t *self)
{
}

void
tree_sequence_print_state(tree_sequence_t *self, FILE *out)
{
    size_t j;
    list_len_t k, l;
    site_t site;

    fprintf(out, "tree_sequence state\n");
    fprintf(out, "num_trees = %d\n", (int) self->num_trees);
    fprintf(out, "alphabet = %d\n", (int) self->alphabet);
    fprintf(out, "sequence_length = %f\n", self->sequence_length);
    fprintf(out, "samples = (%d)\n", (int) self->sample_size);
    for (j = 0; j < self->sample_size; j++) {
        fprintf(out, "\t%d\n", (int) self->samples[j]);
    }
    fprintf(out, "provenance = (%d)\n", (int) self->num_provenance_strings);
    for (j = 0; j < self->num_provenance_strings; j++) {
        fprintf(out, "\t'%s'\n", self->provenance_strings[j]);
    }
    fprintf(out, "nodes (%d)\n", (int) self->nodes.num_records);
    for (j = 0; j < self->nodes.num_records; j++) {
        fprintf(out, "\t%d\t%d\t%d\t%f\t'%s'\t%d\n", (int) j,
                self->nodes.flags[j],
                (int) self->nodes.population[j],
                self->nodes.time[j],
                self->nodes.name[j],
                self->nodes.sample_index_map[j]);
    }
    fprintf(out, "edges = (%d records)\n", (int) self->edges.num_records);
    for (j = 0; j < self->edges.num_records; j++) {
        fprintf(out, "\t%d\t%f\t%f\t%d\t%d",
                (int) j,
                self->edges.left[j],
                self->edges.right[j],
                self->edges.parent[j],
                self->edges.child[j]);
        fprintf(out, "\t|\t%d\t%d\n",
                (int) self->edges.indexes.insertion_order[j],
                (int) self->edges.indexes.removal_order[j]);
    }
    fprintf(out, "sites = (%d records)\n", (int) self->sites.num_records);
    for (j = 0; j < self->sites.num_records; j++) {
        fprintf(out, "\t%d\t%f\t%s\n", (int) j, self->sites.position[j],
                self->sites.ancestral_state[j]);
    }
    fprintf(out, "mutations = (%d records)\n", (int) self->mutations.num_records);
    for (j = 0; j < self->mutations.num_records; j++) {
        fprintf(out, "\t%d\t%d\t%d\t%s\n", (int) j, self->mutations.site[j],
                self->mutations.node[j], self->mutations.derived_state[j]);
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
    fprintf(out, "tree_sites = \n");
    for (j = 0; j < self->num_trees; j++) {
        fprintf(out, "tree %d\t%d sites\n", (int) j, self->sites.tree_sites_length[j]);
        for (k = 0; k < self->sites.tree_sites_length[j]; k++) {
            site = self->sites.tree_sites[j][k];
            fprintf(out, "\tsite %d pos = %f ancestral state = %s, %d mutations\n",
                    site.id, site.position, site.ancestral_state, site.mutations_length);
            for (l = 0; l < site.mutations_length; l++) {
                fprintf(out, "\t\tmutation %d node = %d derived_state = %s\n",
                        site.mutations[l].id, site.mutations[l].node,
                        site.mutations[l].derived_state);
            }
        }
    }

    fprintf(out, "memory\n");
    fprintf(out, "\tsample_size = %d\n", (int) self->sample_size);
    fprintf(out, "\tmax_sample_size = %d\n", (int) self->max_sample_size);
    fprintf(out, "\tnodes.num_records = %d\n", (int) self->nodes.num_records);
    fprintf(out, "\tnodes.max_num_records = %d\n", (int) self->nodes.max_num_records);
    fprintf(out, "\tedges.num_records = %d\n", (int) self->edges.num_records);
    fprintf(out, "\tedges.max_num_records = %d\n", (int) self->edges.max_num_records);
    fprintf(out, "\tmutations.num_records = %d\n", (int) self->mutations.num_records);
    fprintf(out, "\tmutations.max_num_records = %d\n", (int) self->mutations.max_num_records);
    fprintf(out, "\tmutations.total_derived_state_length = %d\n",
            (int) self->mutations.total_derived_state_length);
    fprintf(out, "\tmutations.max_total_derived_state_length = %d\n",
            (int) self->mutations.max_total_derived_state_length);
    fprintf(out, "\tmigrations.num_records = %d\n", (int) self->migrations.num_records);
    fprintf(out, "\tmigrations.max_num_records = %d\n", (int) self->migrations.max_num_records);

    tree_sequence_check_state(self);
}

static int
tree_sequence_alloc_mutations(tree_sequence_t *self)
{
    int ret = MSP_ERR_GENERIC;
    size_t size;

    if (self->sites.num_records > self->sites.max_num_records) {
        self->sites.max_num_records = self->sites.num_records;
        size = self->sites.max_num_records;
        msp_safe_free(self->sites.ancestral_state_length);
        msp_safe_free(self->sites.ancestral_state);
        msp_safe_free(self->sites.position);
        msp_safe_free(self->sites.site_mutations);
        msp_safe_free(self->sites.site_mutations_length);
        msp_safe_free(self->sites.tree_sites_mem);
        self->sites.ancestral_state = malloc(size * sizeof(char *));
        self->sites.ancestral_state_length = malloc(size * sizeof(list_len_t));
        self->sites.position = malloc(size * sizeof(double));
        self->sites.site_mutations_length = malloc(size * sizeof(list_len_t));
        self->sites.site_mutations = malloc(size * sizeof(mutation_t *));
        self->sites.tree_sites_mem = malloc(size * sizeof(site_t));
        if (self->sites.ancestral_state == NULL
                || self->sites.ancestral_state_length == NULL
                || self->sites.position == NULL
                || self->sites.site_mutations == NULL
                || self->sites.site_mutations_length == NULL
                || self->sites.tree_sites_mem == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    if (self->sites.total_ancestral_state_length >
            self->sites.max_total_ancestral_state_length) {
        self->sites.max_total_ancestral_state_length =
            self->sites.total_ancestral_state_length;
        size = self->sites.total_ancestral_state_length;
        msp_safe_free(self->sites.ancestral_state_mem);
        self->sites.ancestral_state_mem = malloc(size * sizeof(char));
        if (self->sites.ancestral_state_mem == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    if (self->mutations.total_derived_state_length >
            self->mutations.max_total_derived_state_length) {
        self->mutations.max_total_derived_state_length =
            self->mutations.total_derived_state_length;
        size = self->mutations.total_derived_state_length;
        msp_safe_free(self->mutations.derived_state_mem);
        self->mutations.derived_state_mem = malloc(size * sizeof(char));
        if (self->mutations.derived_state_mem == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    if (self->mutations.num_records > self->mutations.max_num_records) {
        self->mutations.max_num_records = self->mutations.num_records;
        size = self->mutations.max_num_records;
        msp_safe_free(self->mutations.node);
        msp_safe_free(self->mutations.site);
        msp_safe_free(self->mutations.derived_state);
        msp_safe_free(self->mutations.derived_state_length);
        msp_safe_free(self->sites.site_mutations_mem);
        self->mutations.node = malloc(size * sizeof(node_id_t));
        self->mutations.derived_state = malloc(size * sizeof(char *));
        self->mutations.derived_state_length = malloc(size * sizeof(list_len_t));
        self->mutations.site = malloc(size * sizeof(site_id_t));
        self->sites.site_mutations_mem = malloc(size * sizeof(mutation_t));
        if (self->mutations.site == NULL
                || self->mutations.node == NULL
                || self->mutations.derived_state == NULL
                || self->mutations.derived_state_length == NULL
                || self->sites.site_mutations_mem == NULL) {
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
        msp_safe_free(self->nodes.sample_index_map);
        self->nodes.flags = malloc(size * sizeof(uint32_t));
        self->nodes.time = malloc(size * sizeof(double));
        self->nodes.population = malloc(size * sizeof(population_id_t));
        self->nodes.name = malloc(size * sizeof(char *));
        self->nodes.name_length = malloc(size * sizeof(size_t));
        self->nodes.sample_index_map = malloc(size * sizeof(node_id_t));
        if (self->nodes.flags == NULL
                || self->nodes.time == NULL
                || self->nodes.population == NULL
                || self->nodes.name == NULL
                || self->nodes.name_length == NULL
                || self->nodes.sample_index_map == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    if (self->edges.num_records > self->edges.max_num_records) {
        size = self->edges.num_records;
        self->edges.max_num_records = size;
        msp_safe_free(self->edges.left);
        msp_safe_free(self->edges.right);
        msp_safe_free(self->edges.parent);
        msp_safe_free(self->edges.child);
        msp_safe_free(self->edges.indexes.insertion_order);
        msp_safe_free(self->edges.indexes.removal_order);
        self->edges.left = malloc(size * sizeof(double));
        self->edges.right = malloc(size * sizeof(double));
        self->edges.child = malloc(size * sizeof(node_id_t *));
        self->edges.parent = malloc(size * sizeof(node_id_t));
        self->edges.indexes.insertion_order = malloc(size * sizeof(node_id_t));
        self->edges.indexes.removal_order = malloc(size * sizeof(node_id_t));
        if (self->edges.left == NULL
                || self->edges.right == NULL
                || self->edges.child == NULL
                || self->edges.parent == NULL
                || self->edges.indexes.insertion_order == NULL
                || self->edges.indexes.removal_order == NULL) {
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
        msp_safe_free(self->provenance_strings);
    }
    msp_safe_free(self->samples);
    msp_safe_free(self->nodes.flags);
    msp_safe_free(self->nodes.population);
    msp_safe_free(self->nodes.time);
    msp_safe_free(self->nodes.name);
    msp_safe_free(self->nodes.name_mem);
    msp_safe_free(self->nodes.name_length);
    msp_safe_free(self->nodes.sample_index_map);
    msp_safe_free(self->edges.left);
    msp_safe_free(self->edges.right);
    msp_safe_free(self->edges.child);
    msp_safe_free(self->edges.parent);
    msp_safe_free(self->edges.indexes.insertion_order);
    msp_safe_free(self->edges.indexes.removal_order);
    msp_safe_free(self->sites.ancestral_state);
    msp_safe_free(self->sites.ancestral_state_length);
    msp_safe_free(self->sites.ancestral_state_mem);
    msp_safe_free(self->sites.position);
    msp_safe_free(self->sites.tree_sites_mem);
    msp_safe_free(self->sites.tree_sites);
    msp_safe_free(self->sites.tree_sites_length);
    msp_safe_free(self->sites.site_mutations_mem);
    msp_safe_free(self->sites.site_mutations_length);
    msp_safe_free(self->sites.site_mutations);
    msp_safe_free(self->mutations.node);
    msp_safe_free(self->mutations.site);
    msp_safe_free(self->mutations.derived_state);
    msp_safe_free(self->mutations.derived_state_length);
    msp_safe_free(self->mutations.derived_state_mem);
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
    int ret = MSP_ERR_GENERIC;
    node_id_t child, parent, last_parent, last_child;
    size_t j;
    double left, last_left;

    /* TODO there is no deep reason why we have to have at least two samples,
     * but parts of the code do assume this at the moment. If we relax this
     * assumption we must be careful to test thoroughly. */
    if (self->sample_size < 2) {
        ret = MSP_ERR_INSUFFICIENT_SAMPLES;
        goto out;
    }
    for (j = 0; j < self->edges.num_records; j++) {
        parent = self->edges.parent[j];
        child = self->edges.child[j];
        left = self->edges.left[j];
        if (parent == MSP_NULL_NODE) {
            ret = MSP_ERR_NULL_PARENT;
            goto out;
        }
        if (parent < 0 || parent >= (node_id_t) self->nodes.num_records) {
            ret = MSP_ERR_NODE_OUT_OF_BOUNDS;
            goto out;
        }
        if (j > 0) {
            last_parent = self->edges.parent[j - 1];
            last_child = self->edges.child[j - 1];
            last_left = self->edges.left[j - 1];
            /* Input data must sorted by (time[parent], parent, child, left). */
            if (self->nodes.time[parent] < self->nodes.time[last_parent]) {
                ret = MSP_ERR_EDGES_NOT_SORTED_PARENT_TIME;
                goto out;
            }
            if (self->nodes.time[parent] == self->nodes.time[last_parent]) {
                if (parent < last_parent) {
                    ret = MSP_ERR_EDGES_NOT_SORTED_PARENT;
                    goto out;
                }
                if (parent == last_parent) {
                    if (child < last_child) {
                        ret = MSP_ERR_EDGES_NOT_SORTED_CHILD;
                        goto out;
                    }
                    if (child == last_child) {
                        if (left == last_left) {
                            ret = MSP_ERR_DUPLICATE_EDGES;
                            goto out;
                        } else if (left < last_left) {
                            ret = MSP_ERR_EDGES_NOT_SORTED_LEFT;
                            goto out;
                        }
                    }
                }
            }
        }
        if (child == MSP_NULL_NODE) {
            ret = MSP_ERR_NULL_CHILD;
            goto out;
        }
        if (child < 0 || child >= (node_id_t) self->nodes.num_records) {
            ret = MSP_ERR_NODE_OUT_OF_BOUNDS;
            goto out;
        }
        /* time[child] must be < time[parent] */
        if (self->nodes.time[child] >= self->nodes.time[parent]) {
            ret = MSP_ERR_BAD_NODE_TIME_ORDERING;
            goto out;
        }
        if (self->edges.left[j] >= self->edges.right[j]) {
            ret = MSP_ERR_BAD_EDGE_INTERVAL;
            goto out;
        }
    }
    /* Check the sites */
    for (j = 0; j < self->sites.num_records; j++) {
        if (self->sites.ancestral_state_length[j] != 1) {
            ret = MSP_ERR_BAD_ALPHABET;
            goto out;
        }
        if (self->sites.position[j] < 0
                || self->sites.position[j] >= self->sequence_length) {
            ret = MSP_ERR_BAD_SITE_POSITION;
            goto out;
        }
        if (j > 0) {
            if (self->sites.position[j - 1] >= self->sites.position[j]) {
                ret = MSP_ERR_UNSORTED_SITES;
                goto out;
            }
        }
    }
    for (j = 0; j < self->mutations.num_records; j++) {
        if (self->mutations.site[j] < 0
                || self->mutations.site[j] >= (mutation_id_t) self->sites.num_records) {
            ret = MSP_ERR_SITE_OUT_OF_BOUNDS;
            goto out;
        }
        if (self->mutations.node[j] < 0
                || self->mutations.node[j] >= (node_id_t) self->nodes.num_records) {
            ret = MSP_ERR_NODE_OUT_OF_BOUNDS;
            goto out;
        }
        if (j > 0) {
            if (self->mutations.site[j - 1] > self->mutations.site[j]) {
                ret = MSP_ERR_UNSORTED_MUTATIONS;
                goto out;
            }
            if (self->mutations.site[j - 1] == self->mutations.site[j]) {
                /* Also relaxing this assumption because it's too difficult to
                 * enforce after simplify() has been called.
                 */
                /* t1 = self->nodes.time[self->mutations.parent[j - 1]]; */
                /* t2 = self->nodes.time[self->mutations.parent[j]]; */
                /* if (t1 < t2) { */
                /*     ret = MSP_ERR_UNSORTED_MUTATION_NODES; */
                /*     goto out; */
                /* } */
                /* We are relaxing this condition for now, but might want to
                 * reinstate it later. The issue arises in simplify, where we
                 * can't easily derive the correct mutations when there are
                 * lots of them along a path. */
                /* Within a site, nodes must be unique */
                /* if (self->mutations.parent[j - 1] == self->mutations.parent[j]) { */
                /*     ret = MSP_ERR_DUPLICATE_MUTATION_NODES; */
                /*     goto out; */
                /* } */
            }
        }
    }
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_init_nodes(tree_sequence_t *self)
{
    size_t j, k, size;
    int ret = 0;

    /* Determine the sample size */
    self->sample_size = 0;
    for (j = 0; j < self->nodes.num_records; j++) {
        if (self->nodes.flags[j] & MSP_NODE_IS_SAMPLE) {
            self->sample_size++;
        }
    }
    if (self->nodes.num_records > 0 && self->sample_size < 2) {
        ret = MSP_ERR_INSUFFICIENT_SAMPLES;
        goto out;
    }
    /* We alloc the samples list here because it is a special case; we don't know
     * how big it is until we've read in the data.
     */
    if (self->sample_size > self->max_sample_size) {
        size = self->sample_size;
        msp_safe_free(self->samples);
        self->samples = malloc(size * sizeof(node_id_t));
        if (self->samples == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        self->max_sample_size = size;
    }
    k = 0;
    for (j = 0; j < self->nodes.num_records; j++) {
        self->nodes.sample_index_map[j] = -1;
        if (self->nodes.flags[j] & MSP_NODE_IS_SAMPLE) {
            self->samples[k] = (node_id_t) j;
            self->nodes.sample_index_map[j] = (node_id_t) k;
            k++;
        }
    }
    assert(k == self->sample_size);
out:
    return ret;
}

static int
tree_sequence_init_edges(tree_sequence_t *self)
{
    int ret = 0;
    size_t j;

    self->sequence_length = 0.0;
    for (j = 0; j < self->edges.num_records; j++) {
        self->sequence_length = GSL_MAX(self->sequence_length, self->edges.right[j]);
    }
    return ret;
}

static int
tree_sequence_init_sites(tree_sequence_t *self)
{
    site_id_t j;
    mutation_id_t k;
    int ret = 0;
    bool binary = true;
    size_t offset = 0;

    self->alphabet = MSP_ALPHABET_ASCII;
    for (k = 0; k < (mutation_id_t) self->mutations.num_records; k++) {
        ret = tree_sequence_get_mutation(self, k, self->sites.site_mutations_mem + k);
        if (ret != 0) {
            goto out;
        }
        if (self->mutations.derived_state_length[k] != 1) {
            ret = MSP_ERR_BAD_ALPHABET;
            goto out;
        }
        if (!(self->mutations.derived_state[k][0] == '0' ||
               self->mutations.derived_state[k][0] == '1')) {
            binary = false;
        }
    }
    k = 0;
    for (j = 0; j < (site_id_t) self->sites.num_records; j++) {
        if (self->sites.ancestral_state_length[j] != 1) {
            ret = MSP_ERR_BAD_ALPHABET;
            goto out;
        }
        if (self->sites.ancestral_state[j][0] != '0') {
            binary = false;
        }
        if (self->sites.position[j] < 0
                || self->sites.position[j] >= self->sequence_length) {
            ret = MSP_ERR_BAD_SITE_POSITION;
            goto out;
        }
        if (j > 1) {
            if (self->sites.position[j - 1] >= self->sites.position[j]) {
                ret = MSP_ERR_UNSORTED_SITES;
                goto out;
            }
        }
        self->sites.site_mutations[j] = self->sites.site_mutations_mem + offset;
        self->sites.site_mutations_length[j] = 0;
        /* Go through all mutations for this site */
        while (k < (mutation_id_t) self->mutations.num_records
                && self->mutations.site[k] == j) {
            self->sites.site_mutations_length[j]++;
            offset++;
            k++;
        }
        ret = tree_sequence_get_site(self, j, self->sites.tree_sites_mem + j);
        if (ret != 0) {
            goto out;
        }
    }
    if (binary) {
        self->alphabet = MSP_ALPHABET_BINARY;
    }
out:
    return ret;
}

/* Initialises memory associated with the trees.
 */
static int
tree_sequence_init_trees(tree_sequence_t *self)
{
    int ret = MSP_ERR_GENERIC;
    size_t j, k, tree_index;
    site_id_t site;
    double tree_left, tree_right;
    node_id_t *I = self->edges.indexes.insertion_order;
    node_id_t *O = self->edges.indexes.removal_order;

    tree_left = 0;
    self->num_trees = 0;
    j = 0;
    k = 0;
    while (j < self->edges.num_records) {
        while (self->edges.right[O[k]] == tree_left) {
            k++;
        }
        while (j < self->edges.num_records && self->edges.left[I[j]] == tree_left) {
            j++;
        }
        tree_right = self->sequence_length;
        if (j < self->edges.num_records) {
            tree_right = GSL_MIN(self->edges.left[I[j]], self->edges.right[O[k]]);
        }
        tree_left = tree_right;
        self->num_trees++;
    }

    if (self->num_trees > 0) {
        /* TODO this is an ugly departure from the other patterns of
         * mallocing and using high-water mark memory semantics. Do we really need
         * to have these?
         */
        msp_safe_free(self->sites.tree_sites);
        msp_safe_free(self->sites.tree_sites_length);
        self->sites.tree_sites_length = malloc(self->num_trees * sizeof(list_len_t));
        self->sites.tree_sites = malloc(self->num_trees * sizeof(site_t *));
        if (self->sites.tree_sites == NULL || self->sites.tree_sites_length == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        memset(self->sites.tree_sites_length, 0, self->num_trees * sizeof(list_len_t));
        memset(self->sites.tree_sites, 0, self->num_trees * sizeof(site_t *));

        tree_left = 0;
        tree_index = 0;
        site = 0;
        j = 0;
        k = 0;
        while (j < self->edges.num_records) {
            while (self->edges.right[O[k]] == tree_left) {
                k++;
            }
            while (j < self->edges.num_records && self->edges.left[I[j]] == tree_left) {
                j++;
            }
            tree_right = self->sequence_length;
            if (j < self->edges.num_records) {
                tree_right = GSL_MIN(self->edges.left[I[j]], self->edges.right[O[k]]);
            }
            self->sites.tree_sites[tree_index] = self->sites.tree_sites_mem + site;
            while (site < (site_id_t) self->sites.num_records
                    && self->sites.position[site] < tree_right) {
                self->sites.tree_sites_length[tree_index]++;
                site++;
            }
            tree_left = tree_right;
            tree_index++;
        }
        assert(site == (site_id_t) self->sites.num_records);
        assert(tree_index == self->num_trees);
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

static int WARN_UNUSED
tree_sequence_build_indexes(tree_sequence_t *self)
{
    int ret = MSP_ERR_GENERIC;
    size_t j;
    double *time = self->nodes.time;
    index_sort_t *sort_buff = NULL;

    sort_buff = malloc(self->edges.num_records * sizeof(index_sort_t));
    if (sort_buff == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    /* sort by left and increasing time to give us the order in which
     * records should be inserted */
    for (j = 0; j < self->edges.num_records; j++) {
        sort_buff[j].index = (node_id_t ) j;
        sort_buff[j].first = self->edges.left[j];
        sort_buff[j].second = time[self->edges.parent[j]];
        sort_buff[j].third = self->edges.parent[j];
        sort_buff[j].fourth = self->edges.child[j];
    }
    qsort(sort_buff, self->edges.num_records, sizeof(index_sort_t), cmp_index_sort);
    for (j = 0; j < self->edges.num_records; j++) {
        self->edges.indexes.insertion_order[j] = sort_buff[j].index;
    }
    /* sort by right and decreasing parent time to give us the order in which
     * records should be removed. */
    for (j = 0; j < self->edges.num_records; j++) {
        sort_buff[j].index = (node_id_t ) j;
        sort_buff[j].first = self->edges.right[j];
        sort_buff[j].second = -time[self->edges.parent[j]];
        sort_buff[j].third = -self->edges.parent[j];
        sort_buff[j].fourth = -self->edges.child[j];
    }
    qsort(sort_buff, self->edges.num_records, sizeof(index_sort_t), cmp_index_sort);
    for (j = 0; j < self->edges.num_records; j++) {
        self->edges.indexes.removal_order[j] = sort_buff[j].index;
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
    node_table_t *nodes, edge_table_t *edges, migration_table_t *migrations,
    site_table_t *sites, mutation_table_t *mutations,
    size_t num_provenance_strings, char **provenance_strings)
{
    int ret = 0;

    /* TODO need to do a lot of input validation here. What do we allow to be
     * null? What are the size restrictions on the tables? */
    /* Do we allow zero nodes and edges?? */
    if (nodes == NULL || edges == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }

    self->num_provenance_strings = num_provenance_strings;
    self->nodes.num_records = nodes->num_rows;
    /* name_mem here contains terminal NULLs, so we need more space */
    self->nodes.total_name_length = nodes->total_name_length + nodes->num_rows;
    self->edges.num_records = edges->num_rows;
    self->sites.num_records = 0;
    self->sites.total_ancestral_state_length = 0;
    self->mutations.num_records = 0;
    self->mutations.total_derived_state_length = 0;
    if (sites != NULL) {
        /* We need to allow space for NULLs in string columns */
        self->sites.num_records = sites->num_rows;
        self->sites.total_ancestral_state_length =
            sites->total_ancestral_state_length + sites->num_rows;
    }
    if (mutations != NULL) {
        if (sites == NULL) {
            ret = MSP_ERR_BAD_PARAM_VALUE;
            goto out;
        }
        self->mutations.num_records = mutations->num_rows;
        self->mutations.total_derived_state_length =
            mutations->total_derived_state_length + mutations->num_rows;
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
    memcpy(self->nodes.name_length, nodes->name_length, nodes->num_rows * sizeof(uint32_t));
    ret = init_string_column(nodes->num_rows, nodes->name, nodes->name_length,
            self->nodes.name, self->nodes.name_mem);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_init_nodes(self);
    if (ret != 0) {
        goto out;
    }

    /* Setup the edges */
    memcpy(self->edges.left, edges->left, edges->num_rows * sizeof(double));
    memcpy(self->edges.right, edges->right, edges->num_rows * sizeof(double));
    memcpy(self->edges.parent, edges->parent, edges->num_rows * sizeof(node_id_t));
    memcpy(self->edges.child, edges->child, edges->num_rows * sizeof(node_id_t));
    ret = tree_sequence_init_edges(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_build_indexes(self);
    if (ret != 0) {
        goto out;
    }
    if (sites != NULL) {
        memcpy(self->sites.position, sites->position, sites->num_rows * sizeof(double));
        memcpy(self->sites.ancestral_state_length,
                sites->ancestral_state_length,
                sites->num_rows * sizeof(list_len_t));
        ret = init_string_column(sites->num_rows, sites->ancestral_state,
                sites->ancestral_state_length,
                self->sites.ancestral_state,
                self->sites.ancestral_state_mem);
        if (ret != 0) {
            goto out;
        }
    }
    if (mutations != NULL) {
        memcpy(self->mutations.site, mutations->site, mutations->num_rows * sizeof(site_id_t));
        memcpy(self->mutations.node, mutations->node, mutations->num_rows * sizeof(node_id_t));
        memcpy(self->mutations.derived_state_length,
                mutations->derived_state_length,
                mutations->num_rows * sizeof(uint32_t));
        ret = init_string_column(mutations->num_rows, mutations->derived_state,
                mutations->derived_state_length,
                self->mutations.derived_state, self->mutations.derived_state_mem);
        if (ret != 0) {
            goto out;
        }
    }
    ret = tree_sequence_init_sites(self);
    if (ret != 0) {
        goto out;
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
    node_table_t *nodes, edge_table_t *edges, migration_table_t *migrations,
    site_table_t *sites, mutation_table_t *mutations,
    size_t *num_provenance_strings, char ***provenance_strings)
{
    int ret = -1;
    size_t j;
    double left, right;

    if (nodes == NULL || edges == NULL
            || num_provenance_strings == NULL || provenance_strings == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    /* mutation types and mutations must be specified together */
    if ((sites != NULL) != (mutations != NULL)) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    ret = node_table_reset(nodes);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < self->nodes.num_records; j++) {
        ret = node_table_add_row(nodes, self->nodes.flags[j],
                self->nodes.time[j], self->nodes.population[j],
                self->nodes.name[j]);
        if (ret != 0) {
            goto out;
        }
    }

    ret = edge_table_reset(edges);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < self->edges.num_records; j++) {
        left = self->edges.left[j];
        right = self->edges.right[j];
        ret = edge_table_add_row(edges, left, right, self->edges.parent[j],
                self->edges.child[j]);
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
    if (sites != NULL) {
        ret = site_table_reset(sites);
        if (ret != 0) {
            goto out;
        }
        for (j = 0; j < self->sites.num_records; j++) {
            ret = site_table_add_row(sites, self->sites.position[j],
                    self->sites.ancestral_state[j], self->sites.ancestral_state_length[j]);
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
                    self->mutations.site[j], self->mutations.node[j],
                    self->mutations.derived_state[j], self->mutations.derived_state_length[j]);
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
        {"/sites/position", 1, self->sites.num_records},
        {"/sites/ancestral_state_length", 1, self->sites.num_records},
        {"/mutations/site", 1, self->mutations.num_records},
        {"/mutations/node", 1, self->mutations.num_records},
        {"/mutations/derived_state_length", 1, self->mutations.num_records},
        {"/nodes/flags", 1, self->nodes.num_records},
        {"/nodes/population", 1, self->nodes.num_records},
        {"/nodes/name_length", 1, self->nodes.num_records},
        {"/nodes/time", 1, self->nodes.num_records},
        {"/edges/left", 1, self->edges.num_records},
        {"/edges/right", 1, self->edges.num_records},
        {"/edges/parent", 1, self->edges.num_records},
        {"/edges/child", 1, self->edges.num_records},
        {"/edges/indexes/insertion_order", 1, self->edges.num_records},
        {"/edges/indexes/removal_order", 1, self->edges.num_records},
        {"/migrations/left", 1, self->migrations.num_records},
        {"/migrations/right", 1, self->migrations.num_records},
        {"/migrations/node", 1, self->migrations.num_records},
        {"/migrations/source", 1, self->migrations.num_records},
        {"/migrations/dest", 1, self->migrations.num_records},
        {"/migrations/time", 1, self->migrations.num_records},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _dimension_check);
    size_t j;

    /* First make sure that the root number make sense */
    if (self->edges.num_records > 0) {
        if (self->nodes.num_records == 0) {
            ret = MSP_ERR_FILE_FORMAT;
            goto out;
        }
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
        "/edges/indexes",
        "/nodes",
        "/edges",
        "/sites",
        "/mutations",
        "/migrations",
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
    size_t name_length, ancestral_state_length, derived_state_length;
    struct _dimension_read fields[] = {
        {"/sites/position", &self->sites.num_records},
        {"/sites/ancestral_state", &ancestral_state_length},
        {"/mutations/site", &self->mutations.num_records},
        {"/mutations/derived_state", &derived_state_length},
        {"/provenance", &self->num_provenance_strings},
        {"/nodes/time", &self->nodes.num_records},
        {"/nodes/name", &name_length},
        {"/edges/left", &self->edges.num_records},
        {"/migrations/left", &self->migrations.num_records},
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
    self->nodes.total_name_length = name_length + self->nodes.num_records;
    self->sites.total_ancestral_state_length = ancestral_state_length
        + self->sites.num_records;
    self->mutations.total_derived_state_length = derived_state_length
        + self->mutations.num_records;
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
    htri_t exists;
    struct _hdf5_field_read {
        const char *name;
        hid_t type;
        void *dest;
    };
    size_t name_length, ancestral_state_length, derived_state_length;
    char *name = NULL;
    char *ancestral_state = NULL;
    char *derived_state = NULL;
    struct _hdf5_field_read fields[] = {
        {"/provenance", 0, self->provenance_strings},
        {"/nodes/name", H5T_NATIVE_CHAR, NULL},
        {"/sites/ancestral_state", H5T_NATIVE_CHAR, NULL},
        {"/mutations/derived_state", H5T_NATIVE_CHAR, NULL},
        {"/nodes/name_length", H5T_NATIVE_UINT32, self->nodes.name_length},
        {"/nodes/flags", H5T_NATIVE_UINT32, self->nodes.flags},
        {"/nodes/population", H5T_NATIVE_INT32, self->nodes.population},
        {"/nodes/time", H5T_NATIVE_DOUBLE, self->nodes.time},
        {"/sites/position", H5T_NATIVE_DOUBLE, self->sites.position},
        {"/sites/ancestral_state_length", H5T_NATIVE_UINT32,
            self->sites.ancestral_state_length},
        {"/mutations/site", H5T_NATIVE_INT32, self->mutations.site},
        {"/mutations/node", H5T_NATIVE_INT32, self->mutations.node},
        {"/mutations/derived_state_length", H5T_NATIVE_UINT32,
            self->mutations.derived_state_length},
        {"/edges/left", H5T_NATIVE_DOUBLE, self->edges.left},
        {"/edges/right", H5T_NATIVE_DOUBLE, self->edges.right},
        {"/edges/parent", H5T_NATIVE_INT32, self->edges.parent},
        {"/edges/child", H5T_NATIVE_INT32, self->edges.child},
        {"/edges/indexes/insertion_order", H5T_NATIVE_INT32,
            self->edges.indexes.insertion_order},
        {"/edges/indexes/removal_order", H5T_NATIVE_INT32,
            self->edges.indexes.removal_order},
        {"/migrations/left", H5T_NATIVE_DOUBLE, self->migrations.left},
        {"/migrations/right", H5T_NATIVE_DOUBLE, self->migrations.right},
        {"/migrations/node", H5T_NATIVE_INT32, self->migrations.node},
        {"/migrations/source", H5T_NATIVE_INT32, self->migrations.source},
        {"/migrations/dest", H5T_NATIVE_INT32, self->migrations.dest},
        {"/migrations/time", H5T_NATIVE_DOUBLE, self->migrations.time},
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

    name_length = self->nodes.total_name_length - self->nodes.num_records;
    name = malloc(name_length * sizeof(char));
    if (name == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    fields[1].dest = name;
    ancestral_state_length = self->sites.total_ancestral_state_length
        - self->sites.num_records;
    ancestral_state = malloc(ancestral_state_length * sizeof(char));
    if (ancestral_state == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    fields[2].dest = ancestral_state;
    derived_state_length = self->mutations.total_derived_state_length
        - self->mutations.num_records;
    derived_state = malloc(derived_state_length * sizeof(char));
    if (derived_state == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    fields[3].dest = derived_state;

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
    /* Initialise the node name column */
    ret = validate_length(self->nodes.num_records, self->nodes.name_length,
            name_length);
    if (ret != 0) {
        goto out;
    }
    ret = init_string_column(self->nodes.num_records, name, self->nodes.name_length,
            self->nodes.name, self->nodes.name_mem);
    if (ret != 0) {
        goto out;
    }
    /* Initialise the ancestral_state column */
    ret = validate_length(self->sites.num_records,
            self->sites.ancestral_state_length, ancestral_state_length);
    if (ret != 0) {
        goto out;
    }
    ret = init_string_column(self->sites.num_records,
            ancestral_state, self->sites.ancestral_state_length,
            self->sites.ancestral_state, self->sites.ancestral_state_mem);
    if (ret != 0) {
        goto out;
    }
    /* Initialise the derived_state column */
    ret = validate_length(self->mutations.num_records,
            self->mutations.derived_state_length, derived_state_length);
    if (ret != 0) {
        goto out;
    }
    ret = init_string_column(self->mutations.num_records,
            derived_state, self->mutations.derived_state_length,
            self->mutations.derived_state, self->mutations.derived_state_mem);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_init_nodes(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_init_edges(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_init_sites(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_init_trees(self);
    if (ret != 0) {
        goto out;
    }
    ret = 0;
out:
    if (name != NULL) {
        free(name);
    }
    if (ancestral_state != NULL) {
        free(ancestral_state);
    }
    if (derived_state != NULL) {
        free(derived_state);
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
    ret = tree_sequence_check_hdf5_dimensions(self, file_id);
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
    if (flags & MSP_LOAD_EXTENDED_CHECKS) {
        ret = tree_sequence_check(self);
    }
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
    char *flattened_ancestral_state = NULL;
    size_t flattened_ancestral_state_length;
    char *flattened_derived_state = NULL;
    size_t flattened_derived_state_length;
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
        {"/sites/ancestral_state",
            H5T_STD_I8LE, H5T_NATIVE_CHAR, 0, NULL},
        {"/mutations/derived_state",
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
        {"/edges/left",
            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
            self->edges.num_records, self->edges.left},
        {"/edges/right",
            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
            self->edges.num_records, self->edges.right},
        {"/edges/parent",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->edges.num_records, self->edges.parent},
        {"/edges/child",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->edges.num_records, self->edges.child},
        {"/edges/indexes/insertion_order",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->edges.num_records, self->edges.indexes.insertion_order},
        {"/edges/indexes/removal_order",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->edges.num_records, self->edges.indexes.removal_order},
        {"/sites/position",
            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
            self->sites.num_records, self->sites.position},
        {"/sites/ancestral_state_length",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->sites.num_records, self->sites.ancestral_state_length},
        {"/mutations/site",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->mutations.num_records, self->mutations.site},
        {"/mutations/node",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->mutations.num_records, self->mutations.node},
        {"/mutations/derived_state_length",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->mutations.num_records, self->mutations.derived_state_length},
        {"/migrations/left",
            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
            self->migrations.num_records, self->migrations.left},
        {"/migrations/right",
            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
            self->migrations.num_records, self->migrations.right},
        {"/migrations/time",
            H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
            self->migrations.num_records, self->migrations.time},
        {"/migrations/node",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->migrations.num_records, self->migrations.node},
        {"/migrations/source",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->migrations.num_records, self->migrations.source},
        {"/migrations/dest",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->migrations.num_records, self->migrations.dest},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _hdf5_field_write);
    struct _hdf5_group_write {
        const char *name;
    };
    struct _hdf5_group_write groups[] = {
        {"/sites"},
        {"/mutations"},
        {"/nodes"},
        {"/edges"},
        {"/edges/indexes"},
        {"/migrations"},
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

    assert(self->nodes.total_name_length >= self->nodes.num_records);

    /* Make the arrays to hold the flattened strings */
    flattened_name_length = self->nodes.total_name_length - self->nodes.num_records;
    if (flattened_name_length != 0) {
        flattened_name = malloc(flattened_name_length * sizeof(char));
        if (flattened_name == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        ret = flatten_string_column(self->nodes.total_name_length, self->nodes.name_mem,
                flattened_name);
        fields[1].size = flattened_name_length;
        fields[1].source = flattened_name;
    }

    assert(self->sites.total_ancestral_state_length
            >= self->sites.num_records);
    flattened_ancestral_state_length = self->sites.total_ancestral_state_length
        - self->sites.num_records;
    if (flattened_ancestral_state_length != 0) {
        flattened_ancestral_state = malloc(flattened_ancestral_state_length * sizeof(char));
        if (flattened_ancestral_state == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        ret = flatten_string_column(self->sites.total_ancestral_state_length,
                self->sites.ancestral_state_mem, flattened_ancestral_state);
        fields[2].size = flattened_ancestral_state_length;
        fields[2].source = flattened_ancestral_state;
    }

    assert(self->mutations.total_derived_state_length >= self->mutations.num_records);
    flattened_derived_state_length = self->mutations.total_derived_state_length
        - self->mutations.num_records;
    if (flattened_derived_state_length != 0) {
        flattened_derived_state = malloc(flattened_derived_state_length * sizeof(char));
        if (flattened_derived_state == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        ret = flatten_string_column(self->mutations.total_derived_state_length,
                self->mutations.derived_state_mem, flattened_derived_state);
        fields[3].size = flattened_derived_state_length;
        fields[3].source = flattened_derived_state;
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
            if (flags & MSP_DUMP_ZLIB_COMPRESSION) {
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
    if (flattened_ancestral_state != NULL) {
        free(flattened_ancestral_state);
    }
    if (flattened_derived_state != NULL) {
        free(flattened_derived_state);
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
            status = (herr_t) dataspace_id;
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

int
tree_sequence_get_alphabet(tree_sequence_t *self)
{
    return self->alphabet;
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
tree_sequence_get_num_edges(tree_sequence_t *self)
{
    return self->edges.num_records;
}

size_t
tree_sequence_get_num_migrations(tree_sequence_t *self)
{
    return self->migrations.num_records;
}

size_t
tree_sequence_get_num_sites(tree_sequence_t *self)
{
    return self->sites.num_records;
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

bool
tree_sequence_is_sample(tree_sequence_t *self, node_id_t u)
{
    bool ret = false;

    if (u >= 0 && u < (node_id_t) self->nodes.num_records) {
        ret = self->nodes.flags[u] & MSP_NODE_IS_SAMPLE;
    }
    return ret;
}

/* Accessors for records */

int WARN_UNUSED
tree_sequence_get_pairwise_diversity(tree_sequence_t *self,
    node_id_t *samples, size_t num_samples, double *pi)
{
    int ret = 0;
    sparse_tree_t *tree = NULL;
    double result, denom, n, count;
    site_t *sites;
    list_len_t j, k, num_sites;

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
    ret = sparse_tree_alloc(tree, self, MSP_SAMPLE_COUNTS);
    if (ret != 0) {
        goto out;
    }
    ret = sparse_tree_set_tracked_samples(tree, num_samples, samples);
    if (ret != 0) {
        goto out;
    }
    /* Allocation done; move onto main algorithm. */
    result = 0.0;
    for (ret = sparse_tree_first(tree); ret == 1; ret = sparse_tree_next(tree)) {
        ret = sparse_tree_get_sites(tree, &sites, &num_sites);
        if (ret != 0) {
            goto out;
        }
        for (j = 0; j < num_sites; j++) {
            if (sites[j].mutations_length != 1) {
                ret = MSP_ERR_UNSUPPORTED_OPERATION;
                goto out;
            }
            for (k = 0; k < sites[j].mutations_length; k++) {
                count = (double) tree->num_tracked_samples[sites[j].mutations[k].node];
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
tree_sequence_get_edge(tree_sequence_t *self, size_t index, edge_t *edge)
{
    int ret = 0;

    if (index >= self->edges.num_records) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    edge->left = self->edges.left[index];
    edge->right = self->edges.right[index];
    edge->parent = self->edges.parent[index];
    edge->child = self->edges.child[index];
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
tree_sequence_get_mutation(tree_sequence_t *self, mutation_id_t id, mutation_t *record)
{
    int ret = 0;

    if (id < 0 || id >= (mutation_id_t) self->mutations.num_records) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    record->id = id;
    record->index = (size_t) id; // TODO what is this for?
    record->site = self->mutations.site[id];
    record->node = self->mutations.node[id];
    record->derived_state = self->mutations.derived_state[id];
    record->derived_state_length = self->mutations.derived_state_length[id];
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_get_site(tree_sequence_t *self, site_id_t id, site_t *record)
{
    int ret = 0;

    if (id < 0 || id >= (site_id_t) self->sites.num_records) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    record->id = id;
    record->ancestral_state = self->sites.ancestral_state[id];
    record->ancestral_state_length = self->sites.ancestral_state_length[id];
    record->position = self->sites.position[id];
    record->mutations = self->sites.site_mutations[id];
    record->mutations_length = self->sites.site_mutations_length[id];
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_get_samples(tree_sequence_t *self, node_id_t **samples)
{
    *samples = self->samples;
    return 0;
}

int WARN_UNUSED
tree_sequence_get_sample_index_map(tree_sequence_t *self, node_id_t **sample_index_map)
{
    *sample_index_map = self->nodes.sample_index_map;
    return 0;
}

int WARN_UNUSED
tree_sequence_simplify(tree_sequence_t *self, node_id_t *samples, size_t num_samples,
        int flags, tree_sequence_t *output, node_id_t *sample_map)
{
    int ret = 0;
    simplifier_t *simplifier = NULL;
    node_table_t *nodes = NULL;
    edge_table_t *edges = NULL;
    migration_table_t *migrations = NULL;
    site_table_t *sites = NULL;
    mutation_table_t *mutations = NULL;
    size_t num_provenance_strings;
    char **provenance_strings;

    /* Allocate the tables. */
    nodes = malloc(sizeof(node_table_t));
    if (nodes == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = node_table_alloc(nodes, self->nodes.num_records, self->nodes.total_name_length);
    if (ret != 0) {
        goto out;
    }
    edges = malloc(sizeof(node_table_t));
    if (edges == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = edge_table_alloc(edges, self->edges.num_records);
    if (ret != 0) {
        goto out;
    }
    migrations = malloc(sizeof(node_table_t));
    if (migrations == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = migration_table_alloc(migrations, self->migrations.num_records);
    if (ret != 0) {
        goto out;
    }
    sites = malloc(sizeof(node_table_t));
    if (sites == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = site_table_alloc(sites, self->sites.num_records,
            self->sites.total_ancestral_state_length);
    if (ret != 0) {
        goto out;
    }
    mutations = malloc(sizeof(node_table_t));
    if (mutations == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = mutation_table_alloc(mutations, self->mutations.num_records,
            self->mutations.total_derived_state_length);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_dump_tables_tmp(self, nodes, edges, migrations,
            sites, mutations, &num_provenance_strings, &provenance_strings);
    if (ret != 0) {
        goto out;
    }
    simplifier = malloc(sizeof(simplifier_t));
    if (simplifier == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = simplifier_alloc(simplifier, samples, num_samples,
            nodes, edges, migrations, sites, mutations, 0, flags);
    if (ret != 0) {
        goto out;
    }
    ret = simplifier_run(simplifier, sample_map);
    if (ret != 0) {
        goto out;
    }
    /* We are done with the simplifier object, so free it to save some memory */
    simplifier_free(simplifier);
    free(simplifier);
    simplifier = NULL;

    /* TODO fix provenance strings here...*/
    ret = tree_sequence_load_tables_tmp(output, nodes, edges, migrations,
            sites, mutations, 0, NULL);
out:
    if (nodes != NULL) {
        node_table_free(nodes);
        free(nodes);
    }
    if (edges != NULL) {
        edge_table_free(edges);
        free(edges);
    }
    if (migrations != NULL) {
        migration_table_free(migrations);
        free(migrations);
    }
    if (sites != NULL) {
        site_table_free(sites);
        free(sites);
    }
    if (mutations != NULL) {
        mutation_table_free(mutations);
        free(mutations);
    }
    if (simplifier != NULL) {
        simplifier_free(simplifier);
        free(simplifier);
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
    self->num_nodes = tree_sequence_get_num_nodes(tree_sequence);
    self->num_edges = tree_sequence_get_num_edges(tree_sequence);
    self->tree_sequence = tree_sequence;
    self->insertion_index = 0;
    self->removal_index = 0;
    self->tree_left = 0;
    self->tree_index = (size_t) -1;
    self->edge_list_nodes = malloc(self->num_edges * sizeof(edge_list_t));
    if (self->edge_list_nodes == NULL) {
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
    msp_safe_free(self->edge_list_nodes);
    return ret;
}

void
tree_diff_iterator_print_state(tree_diff_iterator_t *self, FILE *out)
{
    fprintf(out, "tree_diff_iterator state\n");
    fprintf(out, "num_edges = %d\n", (int) self->num_edges);
    fprintf(out, "insertion_index = %d\n", (int) self->insertion_index);
    fprintf(out, "removal_index = %d\n", (int) self->removal_index);
    fprintf(out, "tree_left = %f\n", self->tree_left);
    fprintf(out, "tree_index = %d\n", (int) self->tree_index);
}

int WARN_UNUSED
tree_diff_iterator_next(tree_diff_iterator_t *self, double *ret_left, double *ret_right,
        edge_list_t **edges_out, edge_list_t **edges_in)
{
    int ret = 0;
    node_id_t k;
    double left = self->tree_left;
    double right = self->tree_sequence->sequence_length;
    size_t next_edge_list_node = 0;
    tree_sequence_t *s = self->tree_sequence;
    edge_list_t *out_head = NULL;
    edge_list_t *out_tail = NULL;
    edge_list_t *in_head = NULL;
    edge_list_t *in_tail = NULL;
    edge_list_t *w = NULL;
    size_t num_trees = tree_sequence_get_num_trees(s);

    assert(s != NULL);

    if (self->tree_index + 1 < num_trees) {
        /* First we remove the stale records */
        while (left == s->edges.right[
                s->edges.indexes.removal_order[self->removal_index]]) {
            k = s->edges.indexes.removal_order[self->removal_index];
            assert(next_edge_list_node < self->num_edges);
            w = &self->edge_list_nodes[next_edge_list_node];
            next_edge_list_node++;
            w->edge.left = s->edges.left[k];
            w->edge.right = s->edges.right[k];
            w->edge.parent = s->edges.parent[k];
            w->edge.child = s->edges.child[k];
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
        while (self->insertion_index < self->num_edges &&
                left == s->edges.left[
                    s->edges.indexes.insertion_order[self->insertion_index]]) {
            k = s->edges.indexes.insertion_order[self->insertion_index];
            assert(next_edge_list_node < self->num_edges);
            w = &self->edge_list_nodes[next_edge_list_node];
            next_edge_list_node++;
            w->edge.left = s->edges.left[k];
            w->edge.right = s->edges.right[k];
            w->edge.parent = s->edges.parent[k];
            w->edge.child = s->edges.child[k];
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
        if (self->insertion_index < self->num_edges) {
            right = GSL_MIN(
                s->edges.left[s->edges.indexes.insertion_order[self->insertion_index]],
                s->edges.right[s->edges.indexes.removal_order[self->removal_index]]);
        }
        self->tree_index++;
        ret = 1;
    }
    *edges_out = out_head;
    *edges_in = in_head;
    *ret_left = left;
    *ret_right = right;
    /* Set the left coordinate for the next tree */
    self->tree_left = right;
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
    size_t j;
    node_id_t u;
    node_list_t *w;

    self->left = 0;
    self->right = 0;
    self->index = (size_t) -1;
    /* TODO we should profile this method to see if just doing a single loop over
     * the nodes would be more efficient than multiple memsets.
     */
    memset(self->parent, 0xff, N * sizeof(node_id_t));
    memset(self->left_child, 0xff, N * sizeof(node_id_t));
    memset(self->right_child, 0xff, N * sizeof(node_id_t));
    memset(self->left_sib, 0xff, N * sizeof(node_id_t));
    memset(self->right_sib, 0xff, N * sizeof(node_id_t));
    if (self->flags & MSP_SAMPLE_COUNTS) {
        memset(self->num_samples, 0, N * sizeof(node_id_t));
        memset(self->marked, 0, N * sizeof(uint8_t));
        /* We can't reset the tracked samples via memset because we don't
         * know where the tracked samples are.
         */
        for (j = 0; j < self->num_nodes; j++) {
            if (! tree_sequence_is_sample(self->tree_sequence, (node_id_t) j)) {
                self->num_tracked_samples[j] = 0;
            }
        }
    }
    if (self->flags & MSP_SAMPLE_LISTS) {
        memset(self->sample_list_head, 0, N * sizeof(node_list_t *));
        memset(self->sample_list_tail, 0, N * sizeof(node_list_t *));
    }
    /* Set the sample attributes */
    if (self->sample_size > 0) {
        self->left_root = self->samples[0];
    }
    for (j = 0; j < self->sample_size; j++) {
        u = self->samples[j];
        if (self->flags & MSP_SAMPLE_COUNTS) {
            self->num_samples[u] = 1;
        }
        if (self->flags & MSP_SAMPLE_LISTS) {
            w = &self->sample_list_node_mem[j];
            w->next = NULL;
            w->node = (node_id_t) u;
            self->sample_list_head[u] = w;
            self->sample_list_tail[u] = w;
        }
        /* Set initial roots */
        if (j < self->sample_size - 1) {
            self->right_sib[self->samples[j]] = self->samples[j + 1];
        }
        if (j > 0) {
            self->left_sib[self->samples[j]] = self->samples[j - 1];
        }
    }
    return ret;
}

int WARN_UNUSED
sparse_tree_alloc(sparse_tree_t *self, tree_sequence_t *tree_sequence, int flags)
{
    int ret = MSP_ERR_NO_MEMORY;
    size_t sample_size;
    size_t num_nodes;

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
    self->samples = tree_sequence->samples;
    self->flags = flags;
    self->parent = malloc(num_nodes * sizeof(node_id_t));
    self->left_child = malloc(num_nodes * sizeof(node_id_t));
    self->right_child = malloc(num_nodes * sizeof(node_id_t));
    self->left_sib = malloc(num_nodes * sizeof(node_id_t));
    self->right_sib = malloc(num_nodes * sizeof(node_id_t));
    if (self->parent == NULL || self->left_child == NULL || self->right_child == NULL
            || self->left_sib == NULL || self->right_sib == NULL) {
        goto out;
    }
    /* the maximum possible height of the tree is num_nodes + 1, including
     * the null value. */
    self->stack1 = malloc((num_nodes + 1) * sizeof(node_id_t));
    self->stack2 = malloc((num_nodes + 1) * sizeof(node_id_t));
    if (self->stack1 == NULL || self->stack2 == NULL) {
        goto out;
    }
    if (self->flags & MSP_SAMPLE_COUNTS) {
        self->num_samples = calloc(num_nodes, sizeof(node_id_t));
        self->num_tracked_samples = calloc(num_nodes, sizeof(node_id_t));
        self->marked = calloc(num_nodes, sizeof(uint8_t));
        if (self->num_samples == NULL || self->num_tracked_samples == NULL
                || self->marked == NULL) {
            goto out;
        }
    }
    if (self->flags & MSP_SAMPLE_LISTS) {
        self->sample_list_head = calloc(num_nodes, sizeof(node_list_t *));
        self->sample_list_tail = calloc(num_nodes, sizeof(node_list_t *));
        self->sample_list_node_mem = calloc(sample_size, sizeof(node_list_t));
        if (self->sample_list_head == NULL || self->sample_list_tail == NULL
                || self->sample_list_node_mem == NULL) {
            goto out;
        }
    }
    ret = sparse_tree_clear(self);
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_free(sparse_tree_t *self)
{
    msp_safe_free(self->parent);
    msp_safe_free(self->left_child);
    msp_safe_free(self->right_child);
    msp_safe_free(self->left_sib);
    msp_safe_free(self->right_sib);
    msp_safe_free(self->stack1);
    msp_safe_free(self->stack2);
    msp_safe_free(self->num_samples);
    msp_safe_free(self->num_tracked_samples);
    msp_safe_free(self->marked);
    msp_safe_free(self->sample_list_head);
    msp_safe_free(self->sample_list_tail);
    msp_safe_free(self->sample_list_node_mem);
    return 0;
}

static int WARN_UNUSED
sparse_tree_reset_tracked_samples(sparse_tree_t *self)
{
    int ret = 0;

    if (!(self->flags & MSP_SAMPLE_COUNTS)) {
        ret = MSP_ERR_UNSUPPORTED_OPERATION;
        goto out;
    }
    memset(self->num_tracked_samples, 0, self->num_nodes * sizeof(node_id_t));
out:
    return ret;
}


int WARN_UNUSED
sparse_tree_set_tracked_samples(sparse_tree_t *self, size_t num_tracked_samples,
        node_id_t *tracked_samples)
{
    int ret = MSP_ERR_GENERIC;
    size_t j;
    node_id_t u;

    /* TODO This is not needed when the sparse tree is new. We should use the
     * state machine to check and only reset the tracked samples when needed.
     */
    ret = sparse_tree_reset_tracked_samples(self);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < num_tracked_samples; j++) {
        u = tracked_samples[j];
        if (u < 0 || u >= (node_id_t) self->num_nodes) {
            ret = MSP_ERR_OUT_OF_BOUNDS;
            goto out;
        }
        if (! tree_sequence_is_sample(self->tree_sequence, u)) {
            ret = MSP_ERR_BAD_SAMPLES;
            goto out;
        }
        if (self->num_tracked_samples[u] != 0) {
            ret = MSP_ERR_DUPLICATE_SAMPLE;
            goto out;
        }
        /* Propagate this upwards */
        while (u != MSP_NULL_NODE) {
            self->num_tracked_samples[u] += 1;
            u = self->parent[u];
        }
    }
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_set_tracked_samples_from_sample_list(sparse_tree_t *self,
        node_list_t *head, node_list_t *tail)
{
    int ret = MSP_ERR_GENERIC;
    node_list_t *list_node = head;
    node_id_t u;
    int not_done;

    if (head == NULL || tail == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    /* TODO This is not needed when the sparse tree is new. We should use the
     * state machine to check and only reset the tracked samples when needed.
     */
    ret = sparse_tree_reset_tracked_samples(self);
    if (ret != 0) {
        goto out;
    }
    not_done = 1;
    while (not_done) {
        u = list_node->node;
        /* Propagate this upwards */
        assert(self->num_tracked_samples[u] == 0);
        while (u != MSP_NULL_NODE) {
            self->num_tracked_samples[u] += 1;
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
    self->left_root = source->left_root;
    self->index = source->index;
    self->sites = source->sites;
    self->sites_length = source->sites_length;

    memcpy(self->parent, source->parent, N * sizeof(node_id_t));
    memcpy(self->left_child, source->left_child, N * sizeof(node_id_t));
    memcpy(self->right_child, source->right_child, N * sizeof(node_id_t));
    memcpy(self->left_sib, source->left_sib, N * sizeof(node_id_t));
    memcpy(self->right_sib, source->right_sib, N * sizeof(node_id_t));
    if (self->flags & MSP_SAMPLE_COUNTS) {
        if (! (source->flags & MSP_SAMPLE_COUNTS)) {
            ret = MSP_ERR_UNSUPPORTED_OPERATION;
            goto out;
        }
        memcpy(self->num_samples, source->num_samples, N * sizeof(node_id_t));
    }
    if (self->flags & MSP_SAMPLE_LISTS) {
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
 * counts and sample lists are not considered for equality.
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
        && self->sites_length == other->sites_length
        && self->sites == other->sites
        && memcmp(self->parent, other->parent, N * sizeof(node_id_t)) == 0;
    /* We do not check the children for equality here because
     * the ordering of the children within a parent are essentially irrelevant
     * in terms of topology. Depending on the way in which we approach a given
     * tree we can get different orderings within the children, and so the
     * same tree would not be equal to itself. */
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
sparse_tree_get_num_samples_by_traversal(sparse_tree_t *self, node_id_t u,
        size_t *num_samples)
{
    int ret = 0;
    node_id_t *stack = self->stack1;
    node_id_t v;
    size_t count = 0;
    int stack_top = 0;

    stack[0] = u;
    while (stack_top >= 0) {
        v = stack[stack_top];
        stack_top--;
        if (tree_sequence_is_sample(self->tree_sequence, v)) {
            count++;
        }
        v = self->left_child[v];
        while (v != MSP_NULL_NODE) {
            stack_top++;
            stack[stack_top] = v;
            v = self->right_sib[v];
        }
    }
    *num_samples = count;
    return ret;
}

int WARN_UNUSED
sparse_tree_get_num_samples(sparse_tree_t *self, node_id_t u, size_t *num_samples)
{
    int ret = 0;

    ret = sparse_tree_check_node(self, u);
    if (ret != 0) {
        goto out;
    }

    if (self->flags & MSP_SAMPLE_COUNTS) {
        *num_samples = (size_t) self->num_samples[u];
    } else {
        ret = sparse_tree_get_num_samples_by_traversal(self, u, num_samples);
    }
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_get_num_tracked_samples(sparse_tree_t *self, node_id_t u,
        size_t *num_tracked_samples)
{
    int ret = 0;

    ret = sparse_tree_check_node(self, u);
    if (ret != 0) {
        goto out;
    }
    if (! (self->flags & MSP_SAMPLE_COUNTS)) {
        ret = MSP_ERR_UNSUPPORTED_OPERATION;
        goto out;
    }
    *num_tracked_samples = (size_t) self->num_tracked_samples[u];
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_get_sample_list(sparse_tree_t *self, node_id_t u,
        node_list_t **head, node_list_t **tail)
{
    int ret = 0;

    ret = sparse_tree_check_node(self, u);
    if (ret != 0) {
        goto out;
    }
    if (! (self->flags & MSP_SAMPLE_LISTS)) {
        ret = MSP_ERR_UNSUPPORTED_OPERATION;
        goto out;
    }
    *head = self->sample_list_head[u];
    *tail = self->sample_list_tail[u];
out:
    return ret;
}

bool
sparse_tree_is_sample(sparse_tree_t *self, node_id_t u)
{
    return tree_sequence_is_sample(self->tree_sequence, u);
}

size_t
sparse_tree_get_num_roots(sparse_tree_t *self)
{
    size_t num_roots = 0;
    node_id_t u = self->left_root;

    while (u != MSP_NULL_NODE) {
        u = self->right_sib[u];
        num_roots++;
    }
    return num_roots;
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
    *t = self->tree_sequence->nodes.time[u];
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_get_sites(sparse_tree_t *self, site_t **sites, list_len_t *sites_length)
{
    *sites = self->sites;
    *sites_length = self->sites_length;
    return 0;
}

int WARN_UNUSED
sparse_tree_get_newick(sparse_tree_t *self, size_t precision, double time_scale,
        int flags, size_t buffer_size, char *newick_buffer)
{
    int ret = 0;
    newick_converter_t nc;

    ret = newick_converter_alloc(&nc, self, precision, time_scale, flags);
    if (ret != 0) {
        goto out;
    }
    ret = newick_converter_run(&nc, buffer_size, newick_buffer);
out:
    newick_converter_free(&nc);
    return ret;
}

static void
sparse_tree_check_state(sparse_tree_t *self)
{
    node_id_t u, v;
    size_t j, num_samples;
    int err, c;
    site_t site;
    node_id_t *children = malloc(self->num_nodes * sizeof(node_id_t));

    assert(children != NULL);

    for (j = 0; j < self->sample_size; j++) {
        u = self->samples[j];
        while (self->parent[u] != MSP_NULL_NODE) {
            u = self->parent[u];
        }
        /* assert(u == self->root); */
    }
    for (u = 0; u < (node_id_t) self->num_nodes; u++) {
        c = 0;
        for (v = self->left_child[u]; v != MSP_NULL_NODE; v = self->right_sib[v]) {
            assert(self->parent[v] == u);
            children[c] = v;
            c++;
        }
        for (v = self->right_child[u]; v != MSP_NULL_NODE; v = self->left_sib[v]) {
            assert(c > 0);
            c--;
            assert(v == children[c]);
        }
    }
    for (j = 0; j < self->sites_length; j++) {
        site = self->sites[j];
        assert(self->left <= site.position);
        assert(site.position < self->right);
    }

    if (self->flags & MSP_SAMPLE_COUNTS) {
        assert(self->num_samples != NULL);
        assert(self->num_tracked_samples != NULL);
        for (u = 0; u < (node_id_t) self->num_nodes; u++) {
            err = sparse_tree_get_num_samples_by_traversal(self, u, &num_samples);
            assert(err == 0);
            assert(num_samples == (size_t) self->num_samples[u]);
        }
    } else {
        assert(self->num_samples == NULL);
        assert(self->num_tracked_samples == NULL);
    }
    if (self->flags & MSP_SAMPLE_LISTS) {
        assert(self->sample_list_tail != NULL);
        assert(self->sample_list_head != NULL);
        assert(self->sample_list_node_mem != NULL);
    } else {
        assert(self->sample_list_tail == NULL);
        assert(self->sample_list_head == NULL);
        assert(self->sample_list_node_mem == NULL);
    }

    free(children);
}

void
sparse_tree_print_state(sparse_tree_t *self, FILE *out)
{
    size_t j;
    node_list_t *u;
    site_t site;

    fprintf(out, "Sparse tree state:\n");
    fprintf(out, "flags = %d\n", self->flags);
    fprintf(out, "left = %f\n", self->left);
    fprintf(out, "right = %f\n", self->right);
    fprintf(out, "left_root = %d\n", (int) self->left_root);
    fprintf(out, "index = %d\n", (int) self->index);
    fprintf(out, "node\tparent\tlchild\trchild\tlsib\trsib\n");

    for (j = 0; j < self->num_nodes; j++) {
        fprintf(out, "%d\t%d\t%d\t%d\t%d\t%d", (int) j, self->parent[j], self->left_child[j],
                self->right_child[j], self->left_sib[j], self->right_sib[j]);
        if (self->flags & MSP_SAMPLE_COUNTS) {
            fprintf(out, "\t%d\t%d\t%d", (int) self->num_samples[j],
                    (int) self->num_tracked_samples[j], self->marked[j]);
        }
        if (self->flags & MSP_SAMPLE_LISTS) {
            fprintf(out, "\t[");
            u = self->sample_list_head[j];
            if (u != NULL) {
                while (1) {
                    fprintf(out, "%d ", (int) u->node);
                    if (u == self->sample_list_tail[j]) {
                        break;
                    }
                    u = u->next;
                }
            } else {
                assert(self->sample_list_tail[j] == NULL);
            }

            fprintf(out, "]");
        }
        fprintf(out, "\n");
    }
    fprintf(out, "sites = \n");
    for (j = 0; j < self->sites_length; j++) {
        site = self->sites[j];
        fprintf(out, "\t%d\t%f\t%s\n", site.id, site.position, site.ancestral_state);
    }
    sparse_tree_check_state(self);
}

/* Methods for positioning the tree along the sequence */

static inline void
sparse_tree_propagate_sample_count_loss(sparse_tree_t *self, node_id_t parent,
        node_id_t child)
{
    const node_id_t all_samples_diff = self->num_samples[child];
    const node_id_t tracked_samples_diff = self->num_tracked_samples[child];
    const uint8_t mark = self->mark;
    node_id_t v = parent;

    /* propagate this loss up as far as we can */
    while (v != MSP_NULL_NODE) {
        self->num_samples[v] -= all_samples_diff;
        self->num_tracked_samples[v] -= tracked_samples_diff;
        self->marked[v] = mark;
        v = self->parent[v];
    }
}

static inline void
sparse_tree_propagate_sample_count_gain(sparse_tree_t *self, node_id_t parent,
        node_id_t child)
{
    node_id_t v;
    const node_id_t all_samples_diff = self->num_samples[child];
    const node_id_t tracked_samples_diff = self->num_tracked_samples[child];
    const uint8_t mark = self->mark;

    /* propogate this gain up as far as we can */
    v = parent;
    while (v != MSP_NULL_NODE) {
        self->num_samples[v] += all_samples_diff;
        self->num_tracked_samples[v] += tracked_samples_diff;
        self->marked[v] = mark;
        v = self->parent[v];
    }
}

static inline void
sparse_tree_update_sample_lists(sparse_tree_t *self, node_id_t node)
{
    node_id_t u, v;
    node_list_t **head = self->sample_list_head;
    node_list_t **tail = self->sample_list_tail;

    u = node;
    while (u != MSP_NULL_NODE) {
        if (self->tree_sequence->nodes.flags[u] & MSP_NODE_IS_SAMPLE) {
            tail[u] = head[u];
        } else {
            head[u] = NULL;
            tail[u] = NULL;
        }
        v = self->left_child[u];
        while (v != MSP_NULL_NODE) {
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
            v = self->right_sib[v];
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
    node_id_t k, p, c, u, lsib, rsib, lroot, rroot;
    tree_sequence_t *s = self->tree_sequence;
    node_id_t R = (node_id_t) s->edges.num_records;
    double x;
    bool is_detached;

    if (direction == MSP_DIR_FORWARD) {
        x = self->right;
    } else {
        x = self->left;
    }
    while (out >= 0 && out < R && out_breakpoints[out_order[out]] == x) {
        assert(out < R);
        k = out_order[out];
        p = s->edges.parent[k];
        c = s->edges.child[k];
        lsib = self->left_sib[c];
        rsib = self->right_sib[c];
        if (lsib == MSP_NULL_NODE) {
            self->left_child[p] = rsib;
        } else {
            self->right_sib[lsib] = rsib;
        }
        if (rsib == MSP_NULL_NODE) {
            self->right_child[p] = lsib;
        } else {
            self->left_sib[rsib] = lsib;
        }
        self->parent[c] = MSP_NULL_NODE;
        self->left_sib[c] = MSP_NULL_NODE;
        self->right_sib[c] = MSP_NULL_NODE;
        if (self->flags & MSP_SAMPLE_COUNTS) {
            sparse_tree_propagate_sample_count_loss(self, p, c);
        }
        if (self->flags & MSP_SAMPLE_LISTS) {
            sparse_tree_update_sample_lists(self, p);
        }
        /* Update roots */
        if (self->parent[p] == MSP_NULL_NODE && self->left_child[p] == MSP_NULL_NODE) {
            /* p is a detached root no longer in the tree. Remove it from the root
             * list. */
            lroot = self->left_sib[p];
            rroot = self->right_sib[p];
            if (lroot != MSP_NULL_NODE) {
                self->right_sib[lroot] = rroot;
                self->left_root = lroot;
            }
            if (rroot != MSP_NULL_NODE) {
                self->left_sib[rroot] = lroot;
                self->left_root = rroot;
            }
            self->left_sib[p] = MSP_NULL_NODE;
            self->right_sib[p] = MSP_NULL_NODE;
        }
        /*  Add c to the root list */
        lroot = self->left_sib[self->left_root];
        if (lroot != MSP_NULL_NODE) {
            self->right_sib[lroot] = c;
        }
        self->left_sib[c] = lroot;
        self->left_sib[self->left_root] = c;
        self->right_sib[c] = self->left_root;
        self->left_root = c;
        out += direction;
    }

    while (in >= 0 && in < R && in_breakpoints[in_order[in]] == x) {
        k = in_order[in];
        p = s->edges.parent[k];
        c = s->edges.child[k];
        if (self->parent[c] != MSP_NULL_NODE) {
            ret = MSP_ERR_BAD_EDGESET_CONTRADICTORY_CHILDREN;
            goto out;
        }
        self->parent[c] = p;
        u = self->right_child[p];
        is_detached = self->parent[p] == MSP_NULL_NODE
            && self->left_child[p] == MSP_NULL_NODE;
        lsib = self->left_sib[c];
        rsib = self->right_sib[c];
        if (u == MSP_NULL_NODE) {
            self->left_child[p] = c;
            self->left_sib[c] = MSP_NULL_NODE;
            self->right_sib[c] = MSP_NULL_NODE;
        } else {
            self->right_sib[u] = c;
            self->left_sib[c] = u;
            self->right_sib[c] = MSP_NULL_NODE;
        }
        self->right_child[p] = c;
        if (self->flags & MSP_SAMPLE_COUNTS) {
            sparse_tree_propagate_sample_count_gain(self, p, c);
        }
        if (self->flags & MSP_SAMPLE_LISTS) {
            sparse_tree_update_sample_lists(self, p);
        }
        /* Update roots */
        if (is_detached) {
            /* Replace c with p in root list */
            if (lsib != MSP_NULL_NODE) {
                self->right_sib[lsib] = p;
            }
            if (rsib != MSP_NULL_NODE) {
                self->left_sib[rsib] = p;
            }
            self->left_sib[p] = lsib;
            self->right_sib[p] = rsib;
            self->left_root = p;
        } else {
            /* Remove c from root list */
            if (lsib != MSP_NULL_NODE) {
                self->right_sib[lsib] = rsib;
                self->left_root = lsib;
            }
            if (rsib != MSP_NULL_NODE) {
                self->left_sib[rsib] = lsib;
                self->left_root = rsib;
            }
        }
        in += direction;
    }

    /* Ensure that left_root is the left-most root */
    while (self->left_sib[self->left_root] != MSP_NULL_NODE) {
        self->left_root = self->left_sib[self->left_root];
    }

    self->direction = direction;
    self->index = (size_t) ((int64_t) self->index + direction);
    if (direction == MSP_DIR_FORWARD) {
        self->left = x;
        self->right = s->sequence_length;
        if (out >= 0 && out < R) {
            self->right = GSL_MIN(self->right, out_breakpoints[out_order[out]]);
        }
        if (in >= 0 && in < R) {
            self->right = GSL_MIN(self->right, in_breakpoints[in_order[in]]);
        }
    } else {
        self->right = x;
        self->left = 0;
        if (out >= 0 && out < R) {
            self->left = GSL_MAX(self->left, out_breakpoints[out_order[out]]);
        }
        if (in >= 0 && in < R) {
            self->left = GSL_MAX(self->left, in_breakpoints[in_order[in]]);
        }
    }
    assert(self->left < self->right);
    *out_index = out;
    *in_index = in;
    if (s->sites.num_records > 0) {
        self->sites = s->sites.tree_sites[self->index];
        self->sites_length = s->sites.tree_sites_length[self->index];
    }
    ret = 1;
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_first(sparse_tree_t *self)
{
    int ret = 0;
    tree_sequence_t *s = self->tree_sequence;

    if (s->edges.num_records > 0) {
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
        self->right = 0;

        ret = sparse_tree_advance(self, MSP_DIR_FORWARD,
                s->edges.right, s->edges.indexes.removal_order,
                &self->right_index, s->edges.left,
                s->edges.indexes.insertion_order, &self->left_index, 1);
    }
out:
    return ret;
}

int WARN_UNUSED
sparse_tree_last(sparse_tree_t *self)
{
    int ret = 0;
    tree_sequence_t *s = self->tree_sequence;

    if (s->edges.num_records > 0) {
        /* TODO this is redundant if this is the first usage of the tree. We
         * should add a state machine here so we know what state the tree is
         * in and can take the appropriate actions.
         */
        ret = sparse_tree_clear(self);
        if (ret != 0) {
            goto out;
        }
        self->left_index = (node_id_t) s->edges.num_records - 1;
        self->right_index = (node_id_t) s->edges.num_records - 1;
        self->direction = MSP_DIR_REVERSE;
        self->index = tree_sequence_get_num_trees(s);
        self->left = s->sequence_length;

        ret = sparse_tree_advance(self, MSP_DIR_REVERSE,
                s->edges.left, s->edges.indexes.insertion_order,
                &self->left_index, s->edges.right,
                s->edges.indexes.removal_order, &self->right_index, 1);
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
                s->edges.right, s->edges.indexes.removal_order,
                &self->right_index, s->edges.left,
                s->edges.indexes.insertion_order, &self->left_index, 0);
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
                s->edges.left, s->edges.indexes.insertion_order,
                &self->left_index, s->edges.right,
                s->edges.indexes.removal_order, &self->right_index, 0);
    }
    return ret;
}
