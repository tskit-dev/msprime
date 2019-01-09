/*
** Copyright (C) 2015-2018 University of Oxford
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
#include <float.h>

#include <gsl/gsl_randist.h>

#include "msprime.h"

typedef struct {
    const char *ancestral_state;
    const char *derived_state;
} mutation_type_t;

static const mutation_type_t binary_mutation_types[] = {
    {"0", "1"}
};

static const mutation_type_t acgt_mutation_types[] = {
    {"A", "C"},
    {"A", "G"},
    {"A", "T"},
    {"C", "A"},
    {"C", "G"},
    {"C", "T"},
    {"G", "A"},
    {"G", "C"},
    {"G", "T"},
    {"T", "A"},
    {"T", "C"},
    {"T", "G"},
};

static int
cmp_site(const void *a, const void *b) {
    const tsk_site_t *ia = (const tsk_site_t *) a;
    const tsk_site_t *ib = (const tsk_site_t *) b;
    return (ia->position > ib->position) - (ia->position < ib->position);
}

void
mutgen_print_state(mutgen_t *self, FILE *out)
{
    avl_node_t *a;
    tsk_site_t *site;
    tsk_mutation_t *mutation;
    tsk_tbl_size_t j;

    fprintf(out, "Mutgen state\n");
    fprintf(out, "\tmutation_rate = %f\n", self->mutation_rate);
    fprintf(out, "\tstart_time = %f\n", self->start_time);
    fprintf(out, "\tend_time = %f\n", self->end_time);
    tsk_blkalloc_print_state(&self->allocator, out);

    for (a = self->sites.head; a != NULL; a = a->next) {
        site = (tsk_site_t *) a->item;
        fprintf(out, "%f\t%.*s\t%.*s\n", site->position,
                (int) site->ancestral_state_length, site->ancestral_state,
                (int) site->metadata_length, site->metadata);
        for (j = 0; j < site->mutations_length; j++) {
            mutation = site->mutations + j;
            fprintf(out, "\t%d\t%d\t%.*s\t%.*s\n", mutation->node, mutation->parent,
                    (int) mutation->derived_state_length, mutation->derived_state,
                    (int) mutation->metadata_length, mutation->metadata);
        }
    }
}

int MSP_WARN_UNUSED
mutgen_alloc(mutgen_t *self, double mutation_rate, gsl_rng *rng, int alphabet,
        size_t block_size)
{
    int ret = 0;

    assert(rng != NULL);
    memset(self, 0, sizeof(mutgen_t));
    if (! (alphabet == MSP_ALPHABET_BINARY || alphabet == MSP_ALPHABET_NUCLEOTIDE)) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->alphabet = alphabet;
    self->mutation_rate = mutation_rate;
    self->rng = rng;
    self->start_time = -DBL_MAX;
    self->end_time = DBL_MAX;

    avl_init_tree(&self->sites, cmp_site, NULL);
    if (block_size == 0) {
        block_size = 8192;
    }
    /* In practice this is the minimum we can support */
    block_size = GSL_MAX(block_size, 128);
    ret = tsk_blkalloc_alloc(&self->allocator, block_size);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
out:
    return ret;
}

int
mutgen_free(mutgen_t *self)
{
    tsk_blkalloc_free(&self->allocator);
    return 0;
}

int
mutgen_set_time_interval(mutgen_t *self, double start_time, double end_time)
{
    int ret = 0;

    if (end_time < start_time) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->start_time = start_time;
    self->end_time = end_time;
out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_add_mutation(mutgen_t *self, node_id_t node, double position,
        const char *ancestral_state, const char *derived_state)
{
    int ret = 0;
    tsk_site_t *site = tsk_blkalloc_get(&self->allocator, sizeof(*site));
    tsk_mutation_t *mutation = tsk_blkalloc_get(&self->allocator, sizeof(*mutation));
    avl_node_t* avl_node = tsk_blkalloc_get(&self->allocator, sizeof(*avl_node));

    if (site == NULL || mutation == NULL || avl_node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memset(site, 0, sizeof(*site));
    memset(mutation, 0, sizeof(*mutation));
    site->position = position;
    site->ancestral_state = ancestral_state;
    site->ancestral_state_length = 1;
    mutation->derived_state = derived_state;
    mutation->derived_state_length = 1;
    mutation->node = node;
    mutation->parent = TSK_NULL;
    site->mutations = mutation;
    site->mutations_length = 1;

    avl_init_node(avl_node, site);
    avl_node = avl_insert_node(&self->sites, avl_node);
    assert(avl_node != NULL);
out:
    return ret;
}

static int MSP_WARN_UNUSED
mutget_initialise_sites(mutgen_t *self, tsk_tbl_collection_t *tables)
{
    int ret = 0;
    tsk_site_tbl_t *sites = tables->sites;
    tsk_mutation_tbl_t *mutations = tables->mutations;
    mutation_id_t mutation_id;
    site_id_t site_id;
    tsk_site_t *site;
    avl_node_t *avl_node;
    tsk_mutation_t *site_mutations;
    tsk_tbl_size_t j, num_mutations, length;
    char *buff;

    mutation_id = 0;
    for (site_id = 0; site_id < (site_id_t) sites->num_rows; site_id++) {
        j = (tsk_tbl_size_t) mutation_id;
        while (j < mutations->num_rows && mutations->site[j] == site_id) {
            j++;
        }
        num_mutations = j - (tsk_tbl_size_t) mutation_id;

        site = tsk_blkalloc_get(&self->allocator, sizeof(*site));
        avl_node = tsk_blkalloc_get(&self->allocator, sizeof(*avl_node));
        site_mutations = tsk_blkalloc_get(&self->allocator,
                num_mutations * sizeof(*mutations));
        if (site == NULL || avl_node == NULL || site_mutations == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        site->position = sites->position[site_id];
        site->mutations = site_mutations;
        site->mutations_length = num_mutations;
        /* ancestral state */
        length = sites->ancestral_state_offset[site_id + 1]
            - sites->ancestral_state_offset[site_id];
        buff = tsk_blkalloc_get(&self->allocator, length);
        if (buff == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        memcpy(buff, sites->ancestral_state +
            sites->ancestral_state_offset[site_id], length);
        site->ancestral_state = buff;
        site->ancestral_state_length = length;

        /* metadata */
        length = sites->metadata_offset[site_id + 1]
            - sites->metadata_offset[site_id];
        buff = tsk_blkalloc_get(&self->allocator, length);
        if (buff == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        memcpy(buff, sites->metadata +
            sites->metadata_offset[site_id], length);
        site->metadata = buff;
        site->metadata_length = length;

        for (j = 0; j < num_mutations; j++) {
            /* We don't use the mutation ID directly here, so we can use it as a
             * flag to indicate when a mutation was imported from outside or not.
             * This is used when computing updated mutation parents when we
             * re-export to the tables.. */
            site_mutations[j].id = 1;
            site_mutations[j].node = mutations->node[mutation_id];
            site_mutations[j].parent = mutations->parent[mutation_id];

            /* ancestral state */
            length = mutations->derived_state_offset[mutation_id + 1]
                - mutations->derived_state_offset[mutation_id];
            buff = tsk_blkalloc_get(&self->allocator, length);
            if (buff == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            memcpy(buff, mutations->derived_state +
                mutations->derived_state_offset[mutation_id], length);
            site_mutations[j].derived_state = buff;
            site_mutations[j].derived_state_length = length;

            /* metadata */
            length = mutations->metadata_offset[mutation_id + 1]
                - mutations->metadata_offset[mutation_id];
            buff = tsk_blkalloc_get(&self->allocator, length);
            if (buff == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            memcpy(buff, mutations->metadata +
                mutations->metadata_offset[mutation_id], length);
            site_mutations[j].metadata = buff;
            site_mutations[j].metadata_length = length;

            mutation_id++;
        }
        avl_init_node(avl_node, site);
        avl_node = avl_insert_node(&self->sites, avl_node);
        if (avl_node == NULL) {
            ret = MSP_ERR_DUPLICATE_SITE_POSITION;
            goto out;
        }
    }
out:
    return ret;
}

static int
mutgen_populate_tables(mutgen_t *self, tsk_site_tbl_t *sites, tsk_mutation_tbl_t *mutations)
{
    int ret = 0;
    size_t j;
    site_id_t site_id;
    avl_node_t *a;
    tsk_site_t *site;
    tsk_mutation_t *mutation;
    mutation_id_t new_mutations = 0;
    mutation_id_t parent;

    for (a = self->sites.head; a != NULL; a = a->next) {
        site = (tsk_site_t *) a->item;
        site_id = tsk_site_tbl_add_row(sites, site->position, site->ancestral_state,
                site->ancestral_state_length, site->metadata, site->metadata_length);
        if (site_id < 0) {
            ret = msp_set_tsk_error(site_id);
            goto out;
        }
        for (j = 0; j < site->mutations_length; j++) {
            mutation = site->mutations + j;
            parent = mutation->parent;
            if (parent != TSK_NULL) {
                parent += new_mutations;
            }
            ret = tsk_mutation_tbl_add_row(mutations, site_id,
                    mutation->node, parent,
                    mutation->derived_state, mutation->derived_state_length,
                    mutation->metadata, mutation->metadata_length);
            if (ret < 0) {
                ret = msp_set_tsk_error(ret);
                goto out;
            }
            if (mutation->id == 0) {
                /* We use this flag to track the number of extra mutations we have
                 * generated, which enables us to compute the new mutation parent
                 * ID efficiently. See above for more notes on this */
                new_mutations++;
            }
        }
    }
    ret = 0;
out:
    return ret;
}

int MSP_WARN_UNUSED
mutgen_generate(mutgen_t *self, tsk_tbl_collection_t *tables, int flags)
{
    int ret = 0;
    tsk_node_tbl_t *nodes = tables->nodes;
    tsk_edge_tbl_t *edges = tables->edges;
    size_t j, l, branch_mutations;
    double left, right, branch_length, distance, mu, position;
    node_id_t parent, child;
    const mutation_type_t *mutation_types;
    unsigned long num_mutation_types;
    unsigned long type;
    double start_time = self->start_time;
    double end_time = self->end_time;
    double branch_start, branch_end;
    avl_node_t *avl_node;
    tsk_site_t search;

    avl_clear_tree(&self->sites);
    tsk_blkalloc_reset(&self->allocator);

    if (flags & MSP_KEEP_SITES) {
        ret = mutget_initialise_sites(self, tables);
        if (ret != 0) {
            goto out;
        }
    }
    ret = tsk_site_tbl_clear(tables->sites);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
    ret = tsk_mutation_tbl_clear(tables->mutations);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }

    if (self->alphabet == 0) {
        mutation_types = binary_mutation_types;
        num_mutation_types = 1;
    } else {
        mutation_types = acgt_mutation_types;
        num_mutation_types = 12;
    }

    for (j = 0; j < edges->num_rows; j++) {
        left = edges->left[j];
        right = edges->right[j];
        distance = right - left;
        parent = edges->parent[j];
        child = edges->child[j];
        assert(child >= 0 && child < (node_id_t) nodes->num_rows);
        branch_start = GSL_MAX(start_time, nodes->time[child]);
        branch_end = GSL_MIN(end_time, nodes->time[parent]);
        branch_length = branch_end - branch_start;
        mu = branch_length * distance * self->mutation_rate;
        branch_mutations = gsl_ran_poisson(self->rng, mu);
        for (l = 0; l < branch_mutations; l++) {
            /* Rejection sample positions until we get one we haven't seen before. */
            /* TODO add a maximum number of rejections here */
            do {
                position = gsl_ran_flat(self->rng, left, right);
                search.position = position;
                avl_node = avl_search(&self->sites, &search);
            } while (avl_node != NULL);
            assert(left <= position && position < right);
            type = gsl_rng_uniform_int(self->rng, num_mutation_types);
            ret = mutgen_add_mutation(self, child, position,
                    mutation_types[type].ancestral_state,
                    mutation_types[type].derived_state);
            if (ret != 0) {
                goto out;
            }
        }
    }
    ret = mutgen_populate_tables(self, tables->sites, tables->mutations);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}
