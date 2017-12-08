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

/* ======================================================== *
 * tree sequence
 * ======================================================== */

static void
tree_sequence_check_state(tree_sequence_t *self)
{
    size_t j;
    list_len_t k, l;
    site_t site;
    site_id_t site_id = 0;

    for (j = 0; j < self->num_trees; j++) {
        for (k = 0; k < self->sites.tree_sites_length[j]; k++) {
            site = self->sites.tree_sites[j][k];
            assert(site.id == site_id);
            site_id++;
            for (l = 0; l < site.mutations_length; l++) {
                assert(site.mutations[l].site == site.id);
            }
        }
    }
    assert(self->nodes.metadata_offset[0] == 0);
    assert(self->sites.ancestral_state_offset[0] == 0);
    assert(self->mutations.derived_state_offset[0] == 0);
}

void
tree_sequence_print_state(tree_sequence_t *self, FILE *out)
{
    size_t j;
    list_len_t k, l, m;
    site_t site;

    fprintf(out, "tree_sequence state\n");
    fprintf(out, "num_trees = %d\n", (int) self->num_trees);
    fprintf(out, "alphabet = %d\n", (int) self->alphabet);
    fprintf(out, "sequence_length = %f\n", self->sequence_length);
    fprintf(out, "samples = (%d)\n", (int) self->num_samples);
    for (j = 0; j < self->num_samples; j++) {
        fprintf(out, "\t%d\n", (int) self->samples[j]);
    }
    fprintf(out, "provenance = (%d)\n", (int) self->provenance.num_records);
    for (j = 0; j < self->provenance.num_records; j++) {
        fprintf(out, "\t");
        for (k = self->provenance.timestamp_offset[j];
                k < self->provenance.timestamp_offset[j + 1]; k++) {
            fprintf(out, "%c", self->provenance.timestamp[k]);
        }
        fprintf(out, "\t");
        for (k = self->provenance.record_offset[j];
                k < self->provenance.record_offset[j + 1]; k++) {
            fprintf(out, "%c", self->provenance.record[k]);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "nodes (%d)\n", (int) self->nodes.num_records);
    for (j = 0; j < self->nodes.num_records; j++) {
        fprintf(out, "\t%d\t%d\t%d\t%f\t%d\t", (int) j,
                self->nodes.flags[j],
                (int) self->nodes.population[j],
                self->nodes.time[j],
                self->nodes.sample_index_map[j]);
        for (k = self->nodes.metadata_offset[j]; k < self->nodes.metadata_offset[j + 1]; k++) {
            fprintf(out, "%c", self->nodes.metadata[k]);
        }
        fprintf(out, "\n");
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
        fprintf(out, "\t%d\t%f\t", (int) j, self->sites.position[j]);
        for (k = self->sites.ancestral_state_offset[j];
                k < self->sites.ancestral_state_offset[j + 1]; k++) {
            fprintf(out, "%c", self->sites.ancestral_state[k]);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "mutations = (%d records)\n", (int) self->mutations.num_records);
    for (j = 0; j < self->mutations.num_records; j++) {
        fprintf(out, "\t%d\t%d\t%d\t%d\t", (int) j, self->mutations.site[j],
                self->mutations.node[j], self->mutations.parent[j]);
        for (k = self->mutations.derived_state_offset[j];
                k < self->mutations.derived_state_offset[j + 1]; k++) {
            fprintf(out, "%c", self->mutations.derived_state[k]);
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
    fprintf(out, "tree_sites = \n");
    for (j = 0; j < self->num_trees; j++) {
        fprintf(out, "tree %d\t%d sites\n", (int) j, self->sites.tree_sites_length[j]);
        for (k = 0; k < self->sites.tree_sites_length[j]; k++) {
            site = self->sites.tree_sites[j][k];
            fprintf(out, "\tsite %d pos = %f ancestral state = ", site.id, site.position);
            for (l = 0; l < site.ancestral_state_length; l++) {
                fprintf(out, "%c", site.ancestral_state[l]);
            }
            fprintf(out, " %d mutations\n", site.mutations_length);
            for (l = 0; l < site.mutations_length; l++) {
                fprintf(out, "\t\tmutation %d node = %d derived_state = ",
                        site.mutations[l].id, site.mutations[l].node);
                for (m = 0; m < site.mutations[l].derived_state_length; m++) {
                    fprintf(out, "%c", site.mutations[l].derived_state[m]);
                }
                fprintf(out, "\n");
            }
        }
    }

    fprintf(out, "memory\n");
    fprintf(out, "\tnum_samples = %d\n", (int) self->num_samples);
    fprintf(out, "\tmax_num_samples = %d\n", (int) self->max_num_samples);
    fprintf(out, "\tnodes.num_records = %d\n", (int) self->nodes.num_records);
    fprintf(out, "\tnodes.max_num_records = %d\n", (int) self->nodes.max_num_records);
    fprintf(out, "\tedges.num_records = %d\n", (int) self->edges.num_records);
    fprintf(out, "\tedges.max_num_records = %d\n", (int) self->edges.max_num_records);
    fprintf(out, "\tmutations.num_records = %d\n", (int) self->mutations.num_records);
    fprintf(out, "\tmutations.max_num_records = %d\n", (int) self->mutations.max_num_records);
    fprintf(out, "\tmutations.derived_state_length = %d\n",
            (int) self->mutations.derived_state_length);
    fprintf(out, "\tmutations.max_derived_state_length = %d\n",
            (int) self->mutations.max_derived_state_length);
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
        msp_safe_free(self->sites.ancestral_state_offset);
        msp_safe_free(self->sites.position);
        msp_safe_free(self->sites.site_mutations);
        msp_safe_free(self->sites.site_mutations_length);
        msp_safe_free(self->sites.tree_sites_mem);
        self->sites.ancestral_state_offset = malloc((size + 1) * sizeof(list_len_t));
        self->sites.position = malloc(size * sizeof(double));
        self->sites.site_mutations_length = malloc(size * sizeof(list_len_t));
        self->sites.site_mutations = malloc(size * sizeof(mutation_t *));
        self->sites.tree_sites_mem = malloc(size * sizeof(site_t));
        if (self->sites.ancestral_state_offset == NULL
                || self->sites.position == NULL
                || self->sites.site_mutations == NULL
                || self->sites.site_mutations_length == NULL
                || self->sites.tree_sites_mem == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    if (self->sites.ancestral_state_length > self->sites.max_ancestral_state_length) {
        self->sites.max_ancestral_state_length = self->sites.ancestral_state_length;
        size = self->sites.ancestral_state_length;
        msp_safe_free(self->sites.ancestral_state);
        self->sites.ancestral_state = malloc(size * sizeof(char));
        if (self->sites.ancestral_state == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    if (self->mutations.derived_state_length > self->mutations.max_derived_state_length) {
        self->mutations.max_derived_state_length = self->mutations.derived_state_length;
        size = self->mutations.derived_state_length;
        msp_safe_free(self->mutations.derived_state);
        self->mutations.derived_state= malloc(size * sizeof(char));
        if (self->mutations.derived_state == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    if (self->mutations.num_records > self->mutations.max_num_records) {
        self->mutations.max_num_records = self->mutations.num_records;
        size = self->mutations.max_num_records;
        msp_safe_free(self->mutations.node);
        msp_safe_free(self->mutations.parent);
        msp_safe_free(self->mutations.site);
        msp_safe_free(self->mutations.derived_state_offset);
        msp_safe_free(self->sites.site_mutations_mem);
        self->mutations.node = malloc(size * sizeof(node_id_t));
        self->mutations.parent = malloc(size * sizeof(mutation_id_t));
        self->mutations.derived_state_offset = malloc((size + 1) * sizeof(list_len_t));
        self->mutations.site = malloc(size * sizeof(site_id_t));
        self->sites.site_mutations_mem = malloc(size * sizeof(mutation_t));
        if (self->mutations.site == NULL
                || self->mutations.node == NULL
                || self->mutations.parent == NULL
                || self->mutations.derived_state_offset == NULL
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

    if (self->nodes.metadata_length > self->nodes.max_metadata_length) {
        size = self->nodes.metadata_length;
        self->nodes.max_metadata_length = size;
        msp_safe_free(self->nodes.metadata);
        self->nodes.metadata = malloc(size * sizeof(char));
        if (self->nodes.metadata == NULL) {
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
        msp_safe_free(self->nodes.metadata_offset);
        msp_safe_free(self->nodes.sample_index_map);
        self->nodes.flags = malloc(size * sizeof(uint32_t));
        self->nodes.time = malloc(size * sizeof(double));
        self->nodes.population = malloc(size * sizeof(population_id_t));
        self->nodes.metadata_offset = malloc((size + 1) * sizeof(list_len_t));
        self->nodes.sample_index_map = malloc(size * sizeof(node_id_t));
        if (self->nodes.flags == NULL
                || self->nodes.time == NULL
                || self->nodes.population == NULL
                || self->nodes.metadata_offset == NULL
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
    size_t size;

    if (self->provenance.timestamp_length > self->provenance.max_timestamp_length) {
        size = self->provenance.timestamp_length;
        self->provenance.max_timestamp_length = size;
        self->provenance.timestamp = malloc(size * sizeof(char));
        if (self->provenance.timestamp == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    if (self->provenance.record_length > self->provenance.max_record_length) {
        size = self->provenance.record_length;
        self->provenance.max_record_length = size;
        msp_safe_free(self->provenance.record);
        self->provenance.record = malloc(size * sizeof(char));
        if (self->provenance.record == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
    if (self->provenance.num_records > self->provenance.max_num_records) {
        size = self->provenance.num_records;
        self->provenance.max_num_records = size;
        msp_safe_free(self->provenance.timestamp_offset);
        msp_safe_free(self->provenance.record_offset);
        self->provenance.timestamp_offset = malloc((size + 1) * sizeof(list_len_t));
        self->provenance.record_offset = malloc((size + 1) * sizeof(list_len_t));
        if (self->provenance.timestamp_offset == NULL
                || self->provenance.record_offset == NULL) {
            ret = MSP_ERR_NO_MEMORY;
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
    int ret = MSP_ERR_GENERIC;
    size_t num_sites = self->sites.num_records;
    size_t num_mutations = self->mutations.num_records;
    size_t num_nodes = self->nodes.num_records;
    size_t num_provenance_records = self->provenance.num_records;

    /* Force an allocation of at least one node, site and mutation record because of the
     * one-extra we always have for the offsets. This is an ugly hack,
     * as we should be allocing more than zero for everything else as well. This
     * workaround should go away when we start using tables native in here like we
     * should. */
    self->sites.num_records = GSL_MAX(1, num_sites);
    self->mutations.num_records = GSL_MAX(1, num_mutations);
    self->nodes.num_records = GSL_MAX(1, num_nodes);
    self->provenance.num_records = GSL_MAX(1, num_provenance_records);

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
    /* Ensure that the offsets are set even when we have an empty tree sequence. */
    self->nodes.metadata_offset[0] = 0;
    self->sites.ancestral_state_offset[0] = 0;
    self->mutations.derived_state_offset[0] = 0;
    self->provenance.timestamp_offset[0] = 0;
    self->provenance.record_offset[0] = 0;
    ret = 0;
out:
    /* Reset the size values. See above for rationale. */
    self->sites.num_records = num_sites;
    self->mutations.num_records = num_mutations;
    self->nodes.num_records = num_nodes;
    self->provenance.num_records = num_provenance_records;
    return ret;
}

/* TODO remove this method and make load_tables zero out the struct instead.
 * This initialisation logic is only needed for the long-lived tree sequence
 * structure being reused for different simulation replicates, which should
 * be removed.
 */
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
    msp_safe_free(self->samples);
    msp_safe_free(self->nodes.flags);
    msp_safe_free(self->nodes.population);
    msp_safe_free(self->nodes.time);
    msp_safe_free(self->nodes.metadata);
    msp_safe_free(self->nodes.metadata_offset);
    msp_safe_free(self->nodes.sample_index_map);
    msp_safe_free(self->edges.left);
    msp_safe_free(self->edges.right);
    msp_safe_free(self->edges.child);
    msp_safe_free(self->edges.parent);
    msp_safe_free(self->edges.indexes.insertion_order);
    msp_safe_free(self->edges.indexes.removal_order);
    msp_safe_free(self->sites.ancestral_state);
    msp_safe_free(self->sites.ancestral_state_offset);
    msp_safe_free(self->sites.position);
    msp_safe_free(self->sites.tree_sites_mem);
    msp_safe_free(self->sites.tree_sites);
    msp_safe_free(self->sites.tree_sites_length);
    msp_safe_free(self->sites.site_mutations_mem);
    msp_safe_free(self->sites.site_mutations_length);
    msp_safe_free(self->sites.site_mutations);
    msp_safe_free(self->mutations.node);
    msp_safe_free(self->mutations.site);
    msp_safe_free(self->mutations.parent);
    msp_safe_free(self->mutations.derived_state);
    msp_safe_free(self->mutations.derived_state_offset);
    msp_safe_free(self->migrations.node);
    msp_safe_free(self->migrations.source);
    msp_safe_free(self->migrations.dest);
    msp_safe_free(self->migrations.left);
    msp_safe_free(self->migrations.right);
    msp_safe_free(self->migrations.time);
    msp_safe_free(self->provenance.timestamp);
    msp_safe_free(self->provenance.record);
    msp_safe_free(self->provenance.timestamp_offset);
    msp_safe_free(self->provenance.record_offset);
    return 0;
}

static int
check_offset_array(size_t num_rows, size_t total_length, list_len_t *offset)
{
    int ret = MSP_ERR_BAD_OFFSET;
    int j;

    if (offset[0] != 0) {
        goto out;
    }
    for (j = 0; j < (int) num_rows; j++) {
        if (offset[j] > offset[j + 1]) {
            goto out;
        }
    }
    if (offset[num_rows] != total_length) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_check_offsets(tree_sequence_t *self)
{
    int ret = 0;

    ret = check_offset_array(self->nodes.num_records, self->nodes.metadata_length,
            self->nodes.metadata_offset);
    if (ret != 0) {
        goto out;
    }
    ret = check_offset_array(self->sites.num_records, self->sites.ancestral_state_length,
            self->sites.ancestral_state_offset);
    if (ret != 0) {
        goto out;
    }
    ret = check_offset_array(self->mutations.num_records,
            self->mutations.derived_state_length, self->mutations.derived_state_offset);
    if (ret != 0) {
        goto out;
    }
    ret = check_offset_array(self->provenance.num_records,
            self->provenance.timestamp_length, self->provenance.timestamp_offset);
    if (ret != 0) {
        goto out;
    }
    ret = check_offset_array(self->provenance.num_records,
            self->provenance.record_length, self->provenance.record_offset);
    if (ret != 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static int
tree_sequence_check(tree_sequence_t *self)
{
    int ret = MSP_ERR_GENERIC;
    node_id_t child, parent, last_parent, last_child;
    mutation_id_t parent_mut;
    list_len_t length;
    size_t j;
    double left, last_left;
    double *time = self->nodes.time;
    bool *parent_seen = calloc(self->nodes.num_records, sizeof(bool));

    if (parent_seen == NULL) {
        ret = MSP_ERR_NO_MEMORY;
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
        if (parent_seen[parent]) {
            ret = MSP_ERR_EDGES_NONCONTIGUOUS_PARENTS;
            goto out;
        }
        if (j > 0) {
            last_parent = self->edges.parent[j - 1];
            last_child = self->edges.child[j - 1];
            last_left = self->edges.left[j - 1];
            /* Input data must sorted by (time[parent], parent, child, left). */
            if (time[parent] < time[last_parent]) {
                ret = MSP_ERR_EDGES_NOT_SORTED_PARENT_TIME;
                goto out;
            }
            if (time[parent] == time[last_parent]) {
                /* if (parent < last_parent) { */
                /*     ret = MSP_ERR_EDGES_NOT_SORTED_PARENT; */
                /*     goto out; */
                /* } */
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
                } else {
                    parent_seen[last_parent] = true;
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
        if (time[child] >= time[parent]) {
            ret = MSP_ERR_BAD_NODE_TIME_ORDERING;
            goto out;
        }
        if (self->edges.left[j] >= self->edges.right[j]) {
            ret = MSP_ERR_BAD_EDGE_INTERVAL;
            goto out;
        }
        if (self->edges.right[j] > self->sequence_length) {
            ret = MSP_ERR_RIGHT_GREATER_SEQ_LENGTH;
            goto out;
        }
    }

    /* Check the sites */
    for (j = 0; j < self->sites.num_records; j++) {
        length = self->sites.ancestral_state_offset[j + 1]
            - self->sites.ancestral_state_offset[j];
        if (length != 1) {
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
        parent_mut = self->mutations.parent[j];
        if (parent_mut < MSP_NULL_MUTATION
                || parent_mut >= (mutation_id_t) self->mutations.num_records) {
            ret = MSP_ERR_MUTATION_OUT_OF_BOUNDS;
            goto out;
        }
        if (parent_mut == (mutation_id_t) j) {
            ret = MSP_ERR_MUTATION_PARENT_EQUAL;
            goto out;
        }
        if (parent_mut != MSP_NULL_MUTATION) {
            /* Parents must be listed before their children */
            if (parent_mut > (mutation_id_t) j) {
                ret = MSP_ERR_MUTATION_PARENT_AFTER_CHILD;
                goto out;
            }
            if (self->mutations.site[parent_mut] != self->mutations.site[j]) {
                ret = MSP_ERR_MUTATION_PARENT_DIFFERENT_SITE;
                goto out;
            }
        }
        if (j > 0) {
            if (self->mutations.site[j - 1] > self->mutations.site[j]) {
                ret = MSP_ERR_UNSORTED_MUTATIONS;
                goto out;
            }
        }
    }
    ret = 0;
out:
    msp_safe_free(parent_seen);
    return ret;
}

static int
tree_sequence_init_nodes(tree_sequence_t *self)
{
    size_t j, k, size;
    int ret = 0;

    /* Determine the sample size */
    self->num_samples = 0;
    for (j = 0; j < self->nodes.num_records; j++) {
        if (self->nodes.flags[j] & MSP_NODE_IS_SAMPLE) {
            self->num_samples++;
        }
    }
    /* We alloc the samples list here because it is a special case; we don't know
     * how big it is until we've read in the data.
     */
    if (self->num_samples > self->max_num_samples) {
        size = self->num_samples;
        msp_safe_free(self->samples);
        self->samples = malloc(size * sizeof(node_id_t));
        if (self->samples == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        self->max_num_samples = size;
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
    assert(k == self->num_samples);
out:
    return ret;
}

static int
tree_sequence_init_sites(tree_sequence_t *self)
{
    site_id_t j;
    list_len_t k;
    int ret = 0;
    bool binary = true;
    size_t offset = 0;

    self->alphabet = MSP_ALPHABET_ASCII;
    for (k = 0; k < (list_len_t) self->mutations.num_records; k++) {
        ret = tree_sequence_get_mutation(self, (mutation_id_t) k,
                self->sites.site_mutations_mem + k);
        if (ret != 0) {
            goto out;
        }
        if (self->mutations.derived_state_offset[k] != k) {
            ret = MSP_ERR_BAD_ALPHABET;
        }
        if (!(self->mutations.derived_state[k] == '0' ||
               self->mutations.derived_state[k] == '1')) {
            binary = false;
        }
    }
    k = 0;
    for (j = 0; j < (site_id_t) self->sites.num_records; j++) {
        if (self->sites.ancestral_state_offset[j] >= self->sites.ancestral_state_length) {
            ret = MSP_ERR_BAD_OFFSET;
            goto out;
        }
        if (self->sites.ancestral_state[self->sites.ancestral_state_offset[j]] != '0') {
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
        while (k < (list_len_t) self->mutations.num_records
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
    tree_right = self->sequence_length;
    self->num_trees = 0;
    j = 0;
    k = 0;
    while (j < self->edges.num_records || tree_left < self->sequence_length) {
        while (k < self->edges.num_records && self->edges.right[O[k]] == tree_left) {
            k++;
        }
        while (j < self->edges.num_records && self->edges.left[I[j]] == tree_left) {
            j++;
        }
        tree_right = self->sequence_length;
        if (j < self->edges.num_records) {
            tree_right = GSL_MIN(tree_right, self->edges.left[I[j]]);
        }
        if (k < self->edges.num_records) {
             tree_right = GSL_MIN(tree_right, self->edges.right[O[k]]);
        }
        tree_left = tree_right;
        self->num_trees++;
    }
    assert(self->num_trees > 0);

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
    tree_right = self->sequence_length;
    tree_index = 0;
    site = 0;
    j = 0;
    k = 0;
    while (j < self->edges.num_records || tree_left < self->sequence_length) {
        while (k < self->edges.num_records && self->edges.right[O[k]] == tree_left) {
            k++;
        }
        while (j < self->edges.num_records && self->edges.left[I[j]] == tree_left) {
            j++;
        }
        tree_right = self->sequence_length;
        if (j < self->edges.num_records) {
            tree_right = GSL_MIN(tree_right, self->edges.left[I[j]]);
        }
        if (k < self->edges.num_records) {
             tree_right = GSL_MIN(tree_right, self->edges.right[O[k]]);
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
tree_sequence_load_tables(tree_sequence_t *self, double sequence_length,
    node_table_t *nodes, edge_table_t *edges, migration_table_t *migrations,
    site_table_t *sites, mutation_table_t *mutations,
    provenance_table_t *provenance, int flags)
{
    int ret = 0;
    size_t j;

    if (nodes == NULL || edges == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (sequence_length < 0) {
        ret = MSP_ERR_BAD_SEQUENCE_LENGTH;
        goto out;
    }
    if (sequence_length == 0) {
        /* If sequence_length is 0 we infer it to be equal to the largest right
         * value in the edges. If there are no edges, then this is an error. */
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
    self->nodes.num_records = nodes->num_rows;
    self->nodes.metadata_length = nodes->metadata_length;
    self->edges.num_records = edges->num_rows;
    self->sites.num_records = 0;
    self->sites.ancestral_state_length = 0;
    self->mutations.num_records = 0;
    self->mutations.derived_state_length = 0;
    if (sites != NULL) {
        self->sites.num_records = sites->num_rows;
        self->sites.ancestral_state_length = sites->ancestral_state_length;
    }
    if (mutations != NULL) {
        if (sites == NULL) {
            ret = MSP_ERR_BAD_PARAM_VALUE;
            goto out;
        }
        self->mutations.num_records = mutations->num_rows;
        self->mutations.derived_state_length = mutations->derived_state_length;
    }
    self->migrations.num_records = 0;
    if (migrations != NULL) {
        self->migrations.num_records = migrations->num_rows;
    }
    self->provenance.num_records = 0;
    if (provenance != NULL) {
        self->provenance.num_records = provenance->num_rows;
        self->provenance.timestamp_length = provenance->timestamp_length;
        self->provenance.record_length = provenance->record_length;
    }
    ret = tree_sequence_alloc(self);
    if (ret != 0) {
        goto out;
    }
    memcpy(self->nodes.time, nodes->time, nodes->num_rows * sizeof(double));
    memcpy(self->nodes.flags, nodes->flags, nodes->num_rows * sizeof(uint32_t));
    memcpy(self->nodes.population, nodes->population,
            nodes->num_rows * sizeof(population_id_t));
    memcpy(self->nodes.metadata, nodes->metadata, nodes->metadata_length * sizeof(char));
    memcpy(self->nodes.metadata_offset, nodes->metadata_offset,
            (nodes->num_rows + 1) * sizeof(list_len_t));
    ret = tree_sequence_init_nodes(self);
    if (ret != 0) {
        goto out;
    }

    /* Setup the edges */
    memcpy(self->edges.left, edges->left, edges->num_rows * sizeof(double));
    memcpy(self->edges.right, edges->right, edges->num_rows * sizeof(double));
    memcpy(self->edges.parent, edges->parent, edges->num_rows * sizeof(node_id_t));
    memcpy(self->edges.child, edges->child, edges->num_rows * sizeof(node_id_t));
    if (sites != NULL) {
        memcpy(self->sites.position, sites->position, sites->num_rows * sizeof(double));
        memcpy(self->sites.ancestral_state_offset,
                sites->ancestral_state_offset,
                (sites->num_rows + 1) * sizeof(list_len_t));
        memcpy(self->sites.ancestral_state,
                sites->ancestral_state, sites->ancestral_state_length * sizeof(char));
    }
    if (mutations != NULL) {
        memcpy(self->mutations.site, mutations->site, mutations->num_rows * sizeof(site_id_t));
        memcpy(self->mutations.node, mutations->node, mutations->num_rows * sizeof(node_id_t));
        memcpy(self->mutations.parent, mutations->parent,
                mutations->num_rows * sizeof(mutation_id_t));
        memcpy(self->mutations.derived_state_offset,
                mutations->derived_state_offset,
                (mutations->num_rows + 1) * sizeof(list_len_t));
        memcpy(self->mutations.derived_state,
                mutations->derived_state,
                mutations->derived_state_length * sizeof(char));
        if (ret != 0) {
            goto out;
        }
    }
    if (migrations != NULL) {
        memcpy(self->migrations.left, migrations->left, migrations->num_rows * sizeof(double));
        memcpy(self->migrations.right, migrations->right, migrations->num_rows * sizeof(double));
        memcpy(self->migrations.node, migrations->node, migrations->num_rows * sizeof(node_id_t));
        memcpy(self->migrations.source, migrations->source,
                migrations->num_rows * sizeof(population_id_t));
        memcpy(self->migrations.dest, migrations->dest,
                migrations->num_rows * sizeof(population_id_t));
        memcpy(self->migrations.time, migrations->time, migrations->num_rows * sizeof(double));
    }
    if (provenance != NULL) {
        memcpy(self->provenance.timestamp_offset, provenance->timestamp_offset,
                (provenance->num_rows + 1) * sizeof(list_len_t));
        memcpy(self->provenance.timestamp, provenance->timestamp,
                provenance->timestamp_length * sizeof(char));
        memcpy(self->provenance.record_offset, provenance->record_offset,
                (provenance->num_rows + 1) * sizeof(list_len_t));
        memcpy(self->provenance.record, provenance->record,
                provenance->record_length * sizeof(char));
    }
    ret = tree_sequence_check_offsets(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_init_sites(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_check(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_build_indexes(self);
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
tree_sequence_dump_tables(tree_sequence_t *self,
    node_table_t *nodes, edge_table_t *edges, migration_table_t *migrations,
    site_table_t *sites, mutation_table_t *mutations,
    provenance_table_t *provenance, int flags)
{
    int ret = -1;
    size_t j;
    double left, right;
    list_len_t offset, length, timestamp_offset, timestamp_length,
               record_offset, record_length;

    if (nodes == NULL || edges == NULL) {
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
        offset = self->nodes.metadata_offset[j];
        length = self->nodes.metadata_offset[j + 1] - offset;
        ret = node_table_add_row(nodes, self->nodes.flags[j],
                self->nodes.time[j], self->nodes.population[j],
                self->nodes.metadata + offset, length);
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
            offset = self->sites.ancestral_state_offset[j];
            length = self->sites.ancestral_state_offset[j + 1] - offset;
            ret = site_table_add_row(sites, self->sites.position[j],
                    self->sites.ancestral_state + offset, length);
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
            offset = self->mutations.derived_state_offset[j];
            length = self->mutations.derived_state_offset[j + 1] - offset;
            ret = mutation_table_add_row(mutations,
                    self->mutations.site[j], self->mutations.node[j],
                    self->mutations.parent[j], self->mutations.derived_state + offset,
                    length);
            if (ret != 0) {
                goto out;
            }
        }
    }

    if (provenance != NULL) {
        ret = provenance_table_reset(provenance);
        if (ret != 0) {
            goto out;
        }
        for (j = 0; j < self->provenance.num_records; j++) {
            timestamp_offset = self->provenance.timestamp_offset[j];
            timestamp_length = self->provenance.timestamp_offset[j + 1] - timestamp_offset;
            record_offset = self->provenance.record_offset[j];
            record_length = self->provenance.record_offset[j + 1] - record_offset;
            ret = provenance_table_add_row(provenance,
                    self->provenance.timestamp + timestamp_offset, timestamp_length,
                    self->provenance.record + record_offset, record_length);
            if (ret != 0) {
                goto out;
            }
        }
    }

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

    attr_id = H5Aopen_by_name(file_id, "/", "format_version", H5P_DEFAULT, H5P_DEFAULT);
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

    attr_id = H5Aopen_by_name(file_id, "/", "sequence_length", H5P_DEFAULT, H5P_DEFAULT);
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
    if (dims != 1) {
        ret = MSP_ERR_FILE_FORMAT;
        goto out;
    }
    status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &self->sequence_length);
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
    if (self->sequence_length <= 0.0) {
        ret = MSP_ERR_BAD_SEQUENCE_LENGTH;
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
        size_t size;
    };
    struct _dimension_check fields[] = {
        {"/sites/position", self->sites.num_records},
        {"/sites/ancestral_state", self->sites.ancestral_state_length},
        {"/sites/ancestral_state_offset", self->sites.num_records + 1},

        {"/mutations/site", self->mutations.num_records},
        {"/mutations/node", self->mutations.num_records},
        {"/mutations/parent", self->mutations.num_records},
        {"/mutations/derived_state", self->mutations.derived_state_length},
        {"/mutations/derived_state_offset", self->mutations.num_records + 1},

        {"/nodes/flags", self->nodes.num_records},
        {"/nodes/population", self->nodes.num_records},
        {"/nodes/metadata", self->nodes.metadata_length},
        {"/nodes/metadata_offset", self->nodes.num_records + 1},
        {"/nodes/time", self->nodes.num_records},

        {"/edges/left", self->edges.num_records},
        {"/edges/right", self->edges.num_records},
        {"/edges/parent", self->edges.num_records},
        {"/edges/child", self->edges.num_records},
        {"/edges/indexes/insertion_order", self->edges.num_records},
        {"/edges/indexes/removal_order", self->edges.num_records},

        {"/migrations/left", self->migrations.num_records},
        {"/migrations/right", self->migrations.num_records},
        {"/migrations/node", self->migrations.num_records},
        {"/migrations/source", self->migrations.num_records},
        {"/migrations/dest", self->migrations.num_records},
        {"/migrations/time", self->migrations.num_records},

        {"/provenance/timestamp", self->provenance.timestamp_length},
        {"/provenance/timestamp_offset", self->provenance.num_records + 1},
        {"/provenance/record", self->provenance.record_length},
        {"/provenance/record_offset", self->provenance.num_records + 1},
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
        if (dims[0] != fields[j].size) {
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
        "/edges",
        "/edges/indexes",
        "/nodes",
        "/sites",
        "/mutations",
        "/migrations",
        "/provenance",
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
    struct _dimension_read fields[] = {
        {"/sites/position", &self->sites.num_records},
        {"/sites/ancestral_state", &self->sites.ancestral_state_length},
        {"/mutations/site", &self->mutations.num_records},
        {"/mutations/derived_state", &self->mutations.derived_state_length},
        {"/nodes/time", &self->nodes.num_records},
        {"/nodes/metadata", &self->nodes.metadata_length},
        {"/edges/left", &self->edges.num_records},
        {"/migrations/left", &self->migrations.num_records},
        {"/provenance/timestamp_offset", &self->provenance.num_records},
        {"/provenance/timestamp", &self->provenance.timestamp_length},
        {"/provenance/record", &self->provenance.record_length},
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
    /* provenance is a special case because we have no simple columns. We must
     * have at least one rown in the offsets col or we have an error. */
    if (self->provenance.num_records == 0) {
        goto out;
    }
    self->provenance.num_records -= 1;
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
    struct _hdf5_field_read fields[] = {
        {"/nodes/metadata", H5T_NATIVE_CHAR, self->nodes.metadata},
        {"/nodes/metadata_offset", H5T_NATIVE_UINT32, self->nodes.metadata_offset},
        {"/nodes/flags", H5T_NATIVE_UINT32, self->nodes.flags},
        {"/nodes/population", H5T_NATIVE_INT32, self->nodes.population},
        {"/nodes/time", H5T_NATIVE_DOUBLE, self->nodes.time},
        {"/sites/position", H5T_NATIVE_DOUBLE, self->sites.position},
        {"/sites/ancestral_state", H5T_NATIVE_CHAR, self->sites.ancestral_state},
        {"/sites/ancestral_state_offset", H5T_NATIVE_UINT32,
            self->sites.ancestral_state_offset},
        {"/mutations/site", H5T_NATIVE_INT32, self->mutations.site},
        {"/mutations/node", H5T_NATIVE_INT32, self->mutations.node},
        {"/mutations/parent", H5T_NATIVE_INT32, self->mutations.parent},
        {"/mutations/derived_state", H5T_NATIVE_CHAR, self->mutations.derived_state},
        {"/mutations/derived_state_offset", H5T_NATIVE_UINT32,
            self->mutations.derived_state_offset},
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
        {"/provenance/timestamp", H5T_NATIVE_CHAR, self->provenance.timestamp},
        {"/provenance/timestamp_offset", H5T_NATIVE_UINT32,
            self->provenance.timestamp_offset},
        {"/provenance/record", H5T_NATIVE_CHAR, self->provenance.record},
        {"/provenance/record_offset", H5T_NATIVE_UINT32,
            self->provenance.record_offset},
    };
    size_t num_fields = sizeof(fields) / sizeof(struct _hdf5_field_read);
    size_t j;

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
    ret = tree_sequence_check_offsets(self);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_init_nodes(self);
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
    struct _hdf5_field_write {
        const char *name;
        hid_t storage_type;
        hid_t memory_type;
        size_t size;
        void *source;
    };
    struct _hdf5_field_write fields[] = {
        {"/nodes/metadata",
            H5T_STD_I8LE, H5T_NATIVE_CHAR,
            self->nodes.metadata_length, self->nodes.metadata},
        {"/nodes/metadata_offset",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->nodes.num_records + 1, self->nodes.metadata_offset},
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
        {"/sites/ancestral_state",
            H5T_STD_I8LE, H5T_NATIVE_CHAR,
            self->sites.ancestral_state_length, self->sites.ancestral_state},
        {"/sites/ancestral_state_offset",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->sites.num_records + 1, self->sites.ancestral_state_offset},
        {"/mutations/site",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->mutations.num_records, self->mutations.site},
        {"/mutations/node",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->mutations.num_records, self->mutations.node},
        {"/mutations/parent",
            H5T_STD_I32LE, H5T_NATIVE_INT32,
            self->mutations.num_records, self->mutations.parent},
        {"/mutations/derived_state",
            H5T_STD_I8LE, H5T_NATIVE_CHAR,
            self->mutations.derived_state_length, self->mutations.derived_state},
        {"/mutations/derived_state_offset",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->mutations.num_records + 1, self->mutations.derived_state_offset},
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
        {"/provenance/timestamp",
            H5T_STD_I8LE, H5T_NATIVE_CHAR,
            self->provenance.timestamp_length, self->provenance.timestamp},
        {"/provenance/timestamp_offset",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->provenance.num_records + 1, self->provenance.timestamp_offset},
        {"/provenance/record",
            H5T_STD_I8LE, H5T_NATIVE_CHAR,
            self->provenance.record_length, self->provenance.record},
        {"/provenance/record_offset",
            H5T_STD_U32LE, H5T_NATIVE_UINT32,
            self->provenance.num_records + 1, self->provenance.record_offset},
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
        {"/provenance"},
    };
    size_t num_groups = sizeof(groups) / sizeof(struct _hdf5_group_write);
    size_t j;

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
            if (fields[j].memory_type != H5T_NATIVE_DOUBLE) {
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
        {"sequence_length", 0, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, 1,
            &self->sequence_length},
        /* The sample size attribute is vestigial, and only included because
         * older versions of msprime give a better error condition when confronted
         * with a newer file format. Due to a bug in the way that this attribute
         * was loaded, versions of msprime pre 0.4.0 would complain about a missing
         * attribute rather than giving a File format error. */
        {"sample_size", 0, H5T_STD_U32LE, H5T_NATIVE_UINT32, 1, &unused_value},
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
tree_sequence_get_num_samples(tree_sequence_t *self)
{
    return self->num_samples;
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
tree_sequence_get_num_provenances(tree_sequence_t *self)
{
    return self->provenance.num_records;
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

    if (num_samples < 2 || num_samples > self->num_samples) {
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
    list_len_t offset, length;

    if (index < 0 || index >= (node_id_t) self->nodes.num_records) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    node->time = self->nodes.time[index];
    node->population = self->nodes.population[index];
    node->flags = self->nodes.flags[index];
    offset = self->nodes.metadata_offset[index];
    length = self->nodes.metadata_offset[index + 1] - offset;
    node->metadata = self->nodes.metadata+ offset;
    node->metadata_length = length;
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
    list_len_t offset, length;

    if (id < 0 || id >= (mutation_id_t) self->mutations.num_records) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    offset = self->mutations.derived_state_offset[id];
    length = self->mutations.derived_state_offset[id + 1] - offset;
    record->id = id;
    record->index = (size_t) id; // TODO what is this for?
    record->site = self->mutations.site[id];
    record->node = self->mutations.node[id];
    record->parent = self->mutations.parent[id];
    record->derived_state = self->mutations.derived_state + offset;
    record->derived_state_length = length;
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_get_site(tree_sequence_t *self, site_id_t id, site_t *record)
{
    int ret = 0;
    list_len_t offset, length;

    if (id < 0 || id >= (site_id_t) self->sites.num_records) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    offset = self->sites.ancestral_state_offset[id];
    length = self->sites.ancestral_state_offset[id + 1] - offset;
    record->id = id;
    record->ancestral_state = self->sites.ancestral_state + offset;
    record->ancestral_state_length = length;
    record->position = self->sites.position[id];
    record->mutations = self->sites.site_mutations[id];
    record->mutations_length = self->sites.site_mutations_length[id];
out:
    return ret;
}

int WARN_UNUSED
tree_sequence_get_provenance(tree_sequence_t *self, size_t index, provenance_t *provenance)
{
    int ret = 0;
    list_len_t offset, length;

    if (index >= self->provenance.num_records) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    offset = self->provenance.timestamp_offset[index];
    length = self->provenance.timestamp_offset[index + 1] - offset;
    provenance->timestamp = self->provenance.timestamp + offset;
    provenance->timestamp_length = length;
    offset = self->provenance.record_offset[index];
    length = self->provenance.record_offset[index + 1] - offset;
    provenance->record = self->provenance.record + offset;
    provenance->record_length = length;
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
        int flags, tree_sequence_t *output, node_id_t *node_map)
{
    int ret = 0;
    simplifier_t *simplifier = NULL;
    node_table_t *nodes = NULL;
    edge_table_t *edges = NULL;
    migration_table_t *migrations = NULL;
    site_table_t *sites = NULL;
    mutation_table_t *mutations = NULL;
    provenance_table_t *provenance = NULL;

    /* Allocate the tables. */
    nodes = malloc(sizeof(*nodes));
    edges = malloc(sizeof(*edges));
    migrations = malloc(sizeof(*migrations));
    sites = malloc(sizeof(*sites));
    mutations = malloc(sizeof(*mutations));
    provenance = malloc(sizeof(*provenance));
    if (nodes == NULL || edges == NULL || migrations == NULL || mutations == NULL
            || sites == NULL || provenance == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = node_table_alloc(nodes, self->nodes.num_records, self->nodes.metadata_length);
    if (ret != 0) {
        goto out;
    }
    ret = edge_table_alloc(edges, self->edges.num_records);
    if (ret != 0) {
        goto out;
    }
    ret = migration_table_alloc(migrations, self->migrations.num_records);
    if (ret != 0) {
        goto out;
    }
    ret = site_table_alloc(sites, self->sites.num_records,
            self->sites.ancestral_state_length);
    if (ret != 0) {
        goto out;
    }
    ret = mutation_table_alloc(mutations, self->mutations.num_records,
            self->mutations.derived_state_length);
    if (ret != 0) {
        goto out;
    }
    /* Stick with the defaults sizes for provenance_table_alloc */
    ret = provenance_table_alloc(provenance, 0, 0, 0);
    if (ret != 0) {
        goto out;
    }

    ret = tree_sequence_dump_tables(self, nodes, edges, migrations,
            sites, mutations, provenance, 0);
    if (ret != 0) {
        goto out;
    }
    simplifier = malloc(sizeof(*simplifier));
    if (simplifier == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = simplifier_alloc(simplifier, self->sequence_length,
            samples, num_samples, nodes, edges, migrations, sites, mutations,
            0, flags);
    if (ret != 0) {
        goto out;
    }
    ret = simplifier_run(simplifier, node_map);
    if (ret != 0) {
        goto out;
    }
    /* We are done with the simplifier object, so free it to save some memory */
    simplifier_free(simplifier);
    free(simplifier);
    simplifier = NULL;

    ret = tree_sequence_load_tables(output, self->sequence_length, nodes, edges,
            migrations, sites, mutations, provenance, 0);
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
    if (provenance != NULL) {
        provenance_table_free(provenance);
        free(provenance);
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
        while (self->removal_index < self->num_edges &&
                left == s->edges.right[
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
        right = s->sequence_length;
        if (self->insertion_index < self->num_edges) {
            right = GSL_MIN(right, s->edges.left[
                    s->edges.indexes.insertion_order[self->insertion_index]]);
        }
        if (self->removal_index < self->num_edges) {
            right = GSL_MIN(right, s->edges.right[
                    s->edges.indexes.removal_order[self->removal_index]]);
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
    size_t num_samples = self->tree_sequence->num_samples;
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
    memset(self->above_sample, 0, N * sizeof(bool));
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
    self->left_root = MSP_NULL_NODE;
    if (num_samples > 0) {
        self->left_root = self->samples[0];
    }
    for (j = 0; j < num_samples; j++) {
        u = self->samples[j];
        self->above_sample[u] = true;
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
        if (j < num_samples - 1) {
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
    size_t num_samples;
    size_t num_nodes;

    memset(self, 0, sizeof(sparse_tree_t));
    if (tree_sequence == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    num_nodes = tree_sequence->nodes.num_records;
    num_samples = tree_sequence->num_samples;
    self->num_nodes = num_nodes;
    self->tree_sequence = tree_sequence;
    self->samples = tree_sequence->samples;
    self->flags = flags;
    self->parent = malloc(num_nodes * sizeof(node_id_t));
    self->left_child = malloc(num_nodes * sizeof(node_id_t));
    self->right_child = malloc(num_nodes * sizeof(node_id_t));
    self->left_sib = malloc(num_nodes * sizeof(node_id_t));
    self->right_sib = malloc(num_nodes * sizeof(node_id_t));
    self->above_sample = malloc(num_nodes * sizeof(bool));
    if (self->parent == NULL || self->left_child == NULL || self->right_child == NULL
            || self->left_sib == NULL || self->right_sib == NULL
            || self->above_sample == NULL) {
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
        self->sample_list_node_mem = calloc(num_samples, sizeof(node_list_t));
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
    msp_safe_free(self->above_sample);
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
    bool *is_root = calloc(self->num_nodes, sizeof(bool));

    assert(children != NULL);

    for (j = 0; j < self->tree_sequence->num_samples; j++) {
        u = self->samples[j];
        while (self->parent[u] != MSP_NULL_NODE) {
            u = self->parent[u];
        }
        is_root[u] = true;
    }
    if (self->tree_sequence->num_samples == 0) {
        assert(self->left_root == MSP_NULL_NODE);
    } else {
        assert(self->left_sib[self->left_root] == MSP_NULL_NODE);
    }
    /* Iterate over the roots and make sure they are set */
    for (u = self->left_root; u != MSP_NULL_NODE; u = self->right_sib[u]) {
        assert(is_root[u]);
        is_root[u] = false;
    }
    for (u = 0; u < (node_id_t) self->num_nodes; u++) {
        assert(!is_root[u]);
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
    free(is_root);
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
        fprintf(out, "\t%d\t%f\n", site.id, site.position);
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
    node_id_t k, p, c, u, v, root, lsib, rsib, lroot, rroot;
    tree_sequence_t *s = self->tree_sequence;
    node_id_t R = (node_id_t) s->edges.num_records;
    double x;
    bool above_sample;

    if (direction == MSP_DIR_FORWARD) {
        x = self->right;
    } else {
        x = self->left;
    }
    while (out >= 0 && out < R && out_breakpoints[out_order[out]] == x) {
        assert(out < R);
        k = out_order[out];
        out += direction;
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

        /* Update the roots. If c is not above a sample then we have nothing to do
         * as we cannot affect the status of any roots. */
        if (self->above_sample[c]) {
            /* Compute the new above sample status for the nodes from p up to root. */
            v = p;
            root = v;
            above_sample = false;
            while (v != MSP_NULL_NODE && !above_sample) {
                above_sample = s->nodes.flags[v] & MSP_NODE_IS_SAMPLE;
                u = self->left_child[v];
                while (u != MSP_NULL_NODE) {
                    above_sample = above_sample || self->above_sample[u];
                    u = self->right_sib[u];
                }
                self->above_sample[v] = above_sample;
                root = v;
                v = self->parent[v];
            }
            if (!above_sample) {
                /* root is no longer above samples. Remove it from the root list */
                lroot = self->left_sib[root];
                rroot = self->right_sib[root];
                self->left_root = MSP_NULL_NODE;
                if (lroot != MSP_NULL_NODE) {
                    self->right_sib[lroot] = rroot;
                    self->left_root = lroot;
                }
                if (rroot != MSP_NULL_NODE) {
                    self->left_sib[rroot] = lroot;
                    self->left_root = rroot;
                }
                self->left_sib[root] = MSP_NULL_NODE;
                self->right_sib[root] = MSP_NULL_NODE;
            }
            /* Add c to the root list */
            if (self->left_root != MSP_NULL_NODE) {
                lroot = self->left_sib[self->left_root];
                if (lroot != MSP_NULL_NODE) {
                    self->right_sib[lroot] = c;
                }
                self->left_sib[c] = lroot;
                self->left_sib[self->left_root] = c;
            }
            self->right_sib[c] = self->left_root;
            self->left_root = c;
        }
    }

    while (in >= 0 && in < R && in_breakpoints[in_order[in]] == x) {
        k = in_order[in];
        in += direction;
        p = s->edges.parent[k];
        c = s->edges.child[k];
        if (self->parent[c] != MSP_NULL_NODE) {
            ret = MSP_ERR_BAD_EDGESET_CONTRADICTORY_CHILDREN;
            goto out;
        }
        self->parent[c] = p;
        u = self->right_child[p];
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

        /* Update the roots. */
        if (self->above_sample[c]) {
            v = p;
            root = v;
            above_sample = false;
            while (v != MSP_NULL_NODE && !above_sample) {
                above_sample = self->above_sample[v];
                self->above_sample[v] = self->above_sample[v] || self->above_sample[c];
                root = v;
                v = self->parent[v];
            }
            if (! above_sample) {
                /* Replace c with root in root list */
                if (lsib != MSP_NULL_NODE) {
                    self->right_sib[lsib] = root;
                }
                if (rsib != MSP_NULL_NODE) {
                    self->left_sib[rsib] = root;
                }
                self->left_sib[root] = lsib;
                self->right_sib[root] = rsib;
                self->left_root = root;
            } else {
                /* Remove c from root list */
                /* self->left_root = MSP_NULL_NODE; */
                if (lsib != MSP_NULL_NODE) {
                    self->right_sib[lsib] = rsib;
                    self->left_root = lsib;
                }
                if (rsib != MSP_NULL_NODE) {
                    self->left_sib[rsib] = lsib;
                    self->left_root = rsib;
                }
            }
        }
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
    int ret = 1;
    tree_sequence_t *s = self->tree_sequence;

    self->left = 0;
    self->index = 0;
    self->right = s->sequence_length;
    self->sites = s->sites.tree_sites[0];
    self->sites_length = s->sites.tree_sites_length[0];

    if (s->edges.num_records > 0) {
        /* TODO this is redundant if this is the first usage of the tree. We
         * should add a state machine here so we know what state the tree is
         * in and can take the appropriate actions.
         */
        ret = sparse_tree_clear(self);
        if (ret != 0) {
            goto out;
        }
        self->index = (size_t) -1;
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
    int ret = 1;
    tree_sequence_t *s = self->tree_sequence;

    self->left = 0;
    self->right = s->sequence_length;
    self->index = 0;
    self->sites = s->sites.tree_sites[0];
    self->sites_length = s->sites.tree_sites_length[0];

    if (s->edges.num_records > 0) {
        /* TODO this is redundant if this is the first usage of the tree. We
         * should add a state machine here so we know what state the tree is
         * in and can take the appropriate actions.
         */
        ret = sparse_tree_clear(self);
        if (ret != 0) {
            goto out;
        }
        self->index = tree_sequence_get_num_trees(s);
        self->left_index = (node_id_t) s->edges.num_records - 1;
        self->right_index = (node_id_t) s->edges.num_records - 1;
        self->direction = MSP_DIR_REVERSE;
        self->left = s->sequence_length;
        self->right = 0;

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
