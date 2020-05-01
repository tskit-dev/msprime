/*
** Copyright (C) 2015-2020 University of Oxford
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
#include <gsl/gsl_minmax.h>

#include "msprime.h"

typedef struct mutation_t {
    tsk_id_t id;
    tsk_id_t site;
    tsk_id_t node;
    char *derived_state;
    tsk_size_t derived_state_length;
    char *metadata;
    tsk_size_t metadata_length;
    // additions to tsk_mutation_t:
    double time;
    struct mutation_t *parent;
    struct mutation_t *next;
    bool new;
    bool keep;
} mutation_t;

typedef struct {
    double position;
    char *ancestral_state;
    tsk_size_t ancestral_state_length;
    char *metadata;
    tsk_size_t metadata_length;
    // modifications to tsk_site_t:
    mutation_t *mutations;
    size_t mutations_length;
    bool new;
} site_t;

static int
cmp_site(const void *a, const void *b)
{
    const site_t *ia = (const site_t *) a;
    const site_t *ib = (const site_t *) b;
    return (ia->position > ib->position) - (ia->position < ib->position);
}

static int
cmp_mutationp(const void *a, const void *b)
{
    int out;
    /* the extra cast is to avoid a gcc bug in -Werror=cast-qual:
     * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=81631 */
    const mutation_t *ia = *(const mutation_t **) (uintptr_t) a;
    const mutation_t *ib = *(const mutation_t **) (uintptr_t) b;
    out = (ia->new || ib->new) ? 0 : (ia->id - ib->id) - (ib->id - ia->id);
    if (out == 0) {
        out = (ib->time > ia->time) - (ib->time < ia->time);
    }
    return out;
}

static void
insert_mutation(site_t *site, mutation_t *new)
{
    mutation_t *u;
    u = site->mutations;
    new->next = u;
    site->mutations = new;
    site->mutations_length++;
}

static int MSP_WARN_UNUSED
sort_mutations(site_t *site)
{
    int ret = 0;
    size_t k;
    mutation_t *m;
    size_t num_mutations = site->mutations_length;
    mutation_t **p = NULL;
    if (num_mutations > 0) {
        p = malloc(num_mutations * sizeof(*p));
        if (p == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        k = 0;
        for (m = site->mutations; m != NULL; m = m->next) {
            p[k] = m;
            k++;
        }
        assert(k == num_mutations);
        qsort(p, (size_t) num_mutations, sizeof(mutation_t *), &cmp_mutationp);
        site->mutations = p[0];
        for (k = 0; k < num_mutations; k++) {
            m = p[k];
            if (k == num_mutations - 1) {
                m->next = NULL;
            } else {
                m->next = p[k + 1];
            }
        }
    }
out:
    msp_safe_free(p);
    return ret;
}

/*******************
 * Mutation model */

void
mutation_model_print_state(mutation_model_t *self, FILE *out)
{
    size_t j, k;
    double *mutation_row;

    fprintf(out, "mutation_model (%p):: (%d)", (void *) self, (int) self->num_alleles);
    fprintf(out, "\nroot_distribution =");
    for (j = 0; j < self->num_alleles; j++) {
        fprintf(out, " %s (%0.4f),", self->alleles[j], self->root_distribution[j]);
    }
    fprintf(out, "\n\t------------------------------\n");
    for (j = 0; j < self->num_alleles; j++) {
        mutation_row = self->transition_matrix + j * self->num_alleles;
        fprintf(out, "\t");
        for (k = 0; k < self->num_alleles; k++) {
            fprintf(out, " %0.4f", mutation_row[k]);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\n");
}

static int MSP_WARN_UNUSED
mutation_model_check_validity(mutation_model_t *self)
{
    int ret = 0;
    size_t j, k;
    double prob, min_prob;
    double *mutation_row;

    /* Check probabilities sum to one */
    prob = 0;
    min_prob = 0;
    for (j = 0; j < self->num_alleles; j++) {
        prob += self->root_distribution[j];
        min_prob = GSL_MIN(min_prob, self->root_distribution[j]);
    }
    if (!doubles_almost_equal(prob, 1.0, 1e-12) || min_prob < 0.0) {
        ret = MSP_ERR_BAD_ROOT_PROBABILITIES;
        goto out;
    }
    for (j = 0; j < self->num_alleles; j++) {
        prob = 0;
        min_prob = 0;
        mutation_row = self->transition_matrix + j * self->num_alleles;
        for (k = 0; k < self->num_alleles; k++) {
            prob += mutation_row[k];
            min_prob = GSL_MIN(min_prob, mutation_row[k]);
        }
        if (!doubles_almost_equal(prob, 1.0, 1e-12) || min_prob < 0.0) {
            ret = MSP_ERR_BAD_TRANSITION_MATRIX;
            goto out;
        }
    }

out:
    return ret;
}

static int MSP_WARN_UNUSED
mutation_model_copy_alleles(mutation_model_t *self, char **alleles)
{
    int ret = 0;
    size_t j, len;

    for (j = 0; j < self->num_alleles; j++) {
        len = strlen(alleles[j]);
        self->alleles[j] = malloc(len + 1);
        if (self->alleles[j] == NULL) {
            goto out;
        }
        strcpy(self->alleles[j], alleles[j]);
    }
out:
    return ret;
}

int MSP_WARN_UNUSED
mutation_model_alloc(mutation_model_t *self, size_t num_alleles, char **alleles,
    double *root_distribution, double *transition_matrix)
{
    int ret = 0;

    memset(self, 0, sizeof(mutation_model_t));
    if (num_alleles < 2) {
        ret = MSP_ERR_INSUFFICIENT_ALLELES;
        goto out;
    }
    self->alleles = calloc(num_alleles, sizeof(*self->alleles));
    self->root_distribution = malloc(num_alleles * sizeof(*self->root_distribution));
    self->transition_matrix
        = malloc(num_alleles * num_alleles * sizeof(*self->transition_matrix));
    if (self->alleles == NULL || self->root_distribution == NULL
        || self->transition_matrix == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->num_alleles = num_alleles;
    memcpy(self->root_distribution, root_distribution,
        num_alleles * sizeof(*root_distribution));
    memcpy(self->transition_matrix, transition_matrix,
        num_alleles * num_alleles * sizeof(*transition_matrix));
    ret = mutation_model_copy_alleles(self, alleles);
    if (ret != 0) {
        goto out;
    }
    ret = mutation_model_check_validity(self);
out:
    return ret;
}

int
mutation_model_free(mutation_model_t *self)
{
    size_t j;

    if (self->alleles != NULL) {
        for (j = 0; j < self->num_alleles; j++) {
            msp_safe_free(self->alleles[j]);
        }
    }
    msp_safe_free(self->alleles);
    msp_safe_free(self->root_distribution);
    msp_safe_free(self->transition_matrix);
    return 0;
}

size_t
mutation_model_get_num_alleles(mutation_model_t *self)
{
    return self->num_alleles;
}

static tsk_id_t
mutation_model_allele_index(mutation_model_t *self, const char *allele, size_t length)
{
    tsk_id_t ret = -1;
    tsk_size_t j;

    for (j = 0; j < self->num_alleles; j++) {
        if (length == strlen(self->alleles[j])
            && memcmp(allele, self->alleles[j], length) == 0) {
            ret = (tsk_id_t) j;
            break;
        }
    }
    return ret;
}

static const char *binary_alleles[] = { "0", "1" };
static double binary_model_root_distribution[] = { 1.0, 0.0 };
static double binary_model_transition_matrix[] = { 0.0, 1.0, 1.0, 0.0 };

#define ONE_THIRD (1.0 / 3.0)
static const char *acgt_alleles[] = { "A", "C", "G", "T" };
static double jukes_cantor_model_root_distribution[] = { 0.25, 0.25, 0.25, 0.25 };
static double jukes_cantor_model_transition_matrix[]
    = { 0.0, ONE_THIRD, ONE_THIRD, ONE_THIRD, ONE_THIRD, 0.0, ONE_THIRD, ONE_THIRD,
          ONE_THIRD, ONE_THIRD, 0.0, ONE_THIRD, ONE_THIRD, ONE_THIRD, ONE_THIRD, 0.0 };

/* Populates this mutation model instance with simple mutation model.
 *
 * 0 -> Symmetric 0/1 mutations
 * 1 -> Jukes Cantor ACGT mutations.
 *
 * This is only provided for testing purposes to make it easy to
 * obtain a populated model.
 */
int
mutation_model_factory(mutation_model_t *self, int model)
{
    int ret = MSP_ERR_GENERIC;

    /* We need to cast to uintptr_t * first to work around the annoying pedantry
     * about discarding const qualifiers. */
    if (model == 0) {
        ret = mutation_model_alloc(self, 2, (char **) (uintptr_t *) binary_alleles,
            binary_model_root_distribution, binary_model_transition_matrix);
    } else if (model == 1) {
        ret = mutation_model_alloc(self, 4, (char **) (uintptr_t *) acgt_alleles,
            jukes_cantor_model_root_distribution, jukes_cantor_model_transition_matrix);
    }
    return ret;
}

/***********************
 * Mutation generator */

static void
mutgen_check_state(mutgen_t *self)
{
    size_t j;
    avl_node_t *a;
    site_t *s;
    mutation_t *m;

    for (a = self->sites.head; a != NULL; a = a->next) {
        s = (site_t *) a->item;
        m = s->mutations;
        for (j = 0; j < s->mutations_length; j++) {
            assert(m != NULL);
            assert(m->id >= -1);
            assert(m->node >= 0);
            if (j == s->mutations_length - 1) {
                assert(m->next == NULL);
            }
            m = m->next;
        }
        assert(m == NULL);
    }
}

void
mutgen_print_state(mutgen_t *self, FILE *out)
{
    avl_node_t *a;
    site_t *s;
    mutation_t *m;
    tsk_id_t parent_id;

    fprintf(out, "Mutgen state\n");
    fprintf(out, "\trate_map:\n");
    interval_map_print_state(self->rate_map, out);
    fprintf(out, "\tstart_time = %f\n", self->start_time);
    fprintf(out, "\tend_time = %f\n", self->end_time);
    fprintf(out, "\tmodel:\n");
    mutation_model_print_state(self->model, out);
    tsk_blkalloc_print_state(&self->allocator, out);

    for (a = self->sites.head; a != NULL; a = a->next) {
        s = (site_t *) a->item;
        fprintf(out, "%f\t%.*s\t%.*s\t(%d)\t%d\n", s->position,
            (int) s->ancestral_state_length, s->ancestral_state,
            (int) s->metadata_length, s->metadata, s->new, (int) s->mutations_length);
        for (m = s->mutations; m != NULL; m = m->next) {
            parent_id = m->parent == NULL ? TSK_NULL : m->parent->id;
            fprintf(out, "\t(%d)\t%f\t%d\t%d\t%.*s\t%.*s\t(%d)\t%d\n", m->id, m->time,
                m->node, parent_id, (int) m->derived_state_length, m->derived_state,
                (int) m->metadata_length, m->metadata, m->new, m->keep);
        }
    }
    mutgen_check_state(self);
}

int MSP_WARN_UNUSED
mutgen_alloc(mutgen_t *self, gsl_rng *rng, interval_map_t *rate_map,
    mutation_model_t *model, size_t block_size)
{
    int ret = 0;
    size_t j;

    assert(rng != NULL);
    memset(self, 0, sizeof(mutgen_t));
    self->rng = rng;
    self->rate_map = rate_map;
    self->model = model;
    self->start_time = -DBL_MAX;
    self->end_time = DBL_MAX;
    self->block_size = block_size;

    avl_init_tree(&self->sites, cmp_site, NULL);
    if (block_size == 0) {
        block_size = 8192;
    }
    /* In practice this is the minimum we can support */
    block_size = GSL_MAX(block_size, 128);
    ret = tsk_blkalloc_init(&self->allocator, block_size);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
    for (j = 0; j < rate_map->size - 1; j++) {
        if (rate_map->value[j] < 0) {
            ret = MSP_ERR_BAD_MUTATION_MAP_RATE;
            goto out;
        }
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

int MSP_WARN_UNUSED
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
mutgen_init_allocator(mutgen_t *self, tsk_table_collection_t *tables)
{
    int ret = -1;

    tsk_blkalloc_free(&self->allocator);
    if (self->block_size == 0) {
        /* Default */
        self->block_size = 8192;
    }
    /* This is the effective minimum */
    self->block_size = GSL_MAX(self->block_size, 128);
    /* Need to make sure we have enough space to store sites and mutations. We
     * allocate ancestral and derived states, as well as a list of mutations
     * for each site. This ensures that we can always allocate the required amount.
     * We need to add one because the assert trips when the block size is equal
     * to chunk size (probably wrongly).
     */
    self->block_size
        = GSL_MAX(self->block_size, 1 + tables->sites.ancestral_state_length);
    self->block_size = GSL_MAX(self->block_size, 1 + tables->sites.metadata_length);
    self->block_size
        = GSL_MAX(self->block_size, 1 + tables->mutations.derived_state_length);
    self->block_size = GSL_MAX(self->block_size, 1 + tables->mutations.metadata_length);
    self->block_size = GSL_MAX(
        self->block_size, (1 + tables->mutations.num_rows) * sizeof(tsk_mutation_t));
    ret = tsk_blkalloc_init(&self->allocator, self->block_size);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_add_site(mutgen_t *self, double position, bool new, const char *ancestral_state,
    tsk_size_t ancestral_state_length, const char *metadata, tsk_size_t metadata_length,
    avl_node_t **avl_nodep)
{
    int ret = 0;
    char *buff;
    avl_node_t *avl_node;
    site_t *site;

    avl_node = tsk_blkalloc_get(&self->allocator, sizeof(*avl_node));
    site = tsk_blkalloc_get(&self->allocator, sizeof(*site));
    if (site == NULL || avl_node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memset(site, 0, sizeof(*site));
    site->position = position;
    site->new = new;

    /* ancestral state */
    buff = tsk_blkalloc_get(&self->allocator, ancestral_state_length);
    if (buff == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memcpy(buff, ancestral_state, ancestral_state_length);
    site->ancestral_state = buff;
    site->ancestral_state_length = ancestral_state_length;

    /* metadata */
    buff = tsk_blkalloc_get(&self->allocator, metadata_length);
    if (buff == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memcpy(buff, metadata, metadata_length);
    site->metadata = buff;
    site->metadata_length = metadata_length;

    avl_init_node(avl_node, site);
    avl_node = avl_insert_node(&self->sites, avl_node);
    if (avl_node == NULL) {
        ret = MSP_ERR_DUPLICATE_SITE_POSITION;
        goto out;
    }
    *avl_nodep = avl_node;

out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_add_mutation(mutgen_t *self, site_t *site, tsk_id_t id, node_id_t node,
    double time, bool new, const char *derived_state, tsk_size_t derived_state_length,
    const char *metadata, tsk_size_t metadata_length)
{
    int ret = 0;
    char *buff;

    mutation_t *mutation = tsk_blkalloc_get(&self->allocator, sizeof(*mutation));
    if (mutation == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memset(mutation, 0, sizeof(*mutation));
    mutation->id = id, mutation->node = node;
    mutation->parent = NULL;
    mutation->time = time;
    mutation->new = new;
    mutation->keep = true;

    /* derived state */
    buff = tsk_blkalloc_get(&self->allocator, derived_state_length);
    if (buff == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memcpy(buff, derived_state, derived_state_length);
    mutation->derived_state = buff;
    mutation->derived_state_length = derived_state_length;

    /* metadata */
    buff = tsk_blkalloc_get(&self->allocator, metadata_length);
    if (buff == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memcpy(buff, metadata, metadata_length);
    mutation->metadata = buff;
    mutation->metadata_length = metadata_length;

    insert_mutation(site, mutation);
out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_initialise_sites(mutgen_t *self, tsk_table_collection_t *tables)
{
    int ret = 0;
    tsk_site_table_t *sites = &tables->sites;
    tsk_mutation_table_t *mutations = &tables->mutations;
    tsk_node_table_t *nodes = &tables->nodes;
    site_id_t site_id;
    avl_node_t *avl_node;
    const char *state, *metadata;
    tsk_size_t j, length, metadata_length;

    j = 0;
    for (site_id = 0; site_id < (site_id_t) sites->num_rows; site_id++) {
        avl_node = tsk_blkalloc_get(&self->allocator, sizeof(*avl_node));
        if (avl_node == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        state = sites->ancestral_state + sites->ancestral_state_offset[site_id];
        length = sites->ancestral_state_offset[site_id + 1]
                 - sites->ancestral_state_offset[site_id];
        metadata = sites->metadata + sites->metadata_offset[site_id];
        metadata_length
            = sites->metadata_offset[site_id + 1] - sites->metadata_offset[site_id];
        ret = mutgen_add_site(self, sites->position[site_id], false, state, length,
            metadata, metadata_length, &avl_node);
        if (ret != 0) {
            goto out;
        }

        while (j < mutations->num_rows && mutations->site[j] == site_id) {
            assert(j < mutations->num_rows);
            state = mutations->derived_state + mutations->derived_state_offset[j];
            length = mutations->derived_state_offset[j + 1]
                     - mutations->derived_state_offset[j];
            metadata = mutations->metadata + mutations->metadata_offset[j];
            metadata_length
                = mutations->metadata_offset[j + 1] - mutations->metadata_offset[j];
            ret = mutgen_add_mutation(self, (site_t *) avl_node->item, (int) j,
                mutations->node[j], nodes->time[mutations->node[j]], false, state,
                length, metadata, metadata_length);
            if (ret != 0) {
                goto out;
            }
            j++;
        }
    }

out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_populate_tables(
    mutgen_t *self, tsk_site_table_t *sites, tsk_mutation_table_t *mutations)
{
    int ret = 0;
    tsk_id_t site_id, mutation_id, parent_id;
    avl_node_t *a;
    site_t *site;
    mutation_t *m;
    size_t num_mutations;

    site_id = 0;
    for (a = self->sites.head; a != NULL; a = a->next) {
        site = (site_t *) a->item;
        num_mutations = 0;
        for (m = site->mutations; m != NULL; m = m->next) {
            if (m->keep) {
                if (m->parent == NULL) {
                    parent_id = TSK_NULL;
                } else {
                    parent_id = m->parent->id;
                    assert(parent_id != TSK_NULL);
                }
                mutation_id = tsk_mutation_table_add_row(mutations, site_id, m->node,
                    parent_id, m->derived_state, m->derived_state_length, m->metadata,
                    m->metadata_length);
                if (mutation_id < 0) {
                    ret = msp_set_tsk_error(mutation_id);
                    goto out;
                }
                assert(mutation_id > parent_id);
                m->id = mutation_id;
                num_mutations++;
            }
        }
        /* Omit any new sites that have no mutations */
        if ((!site->new) || num_mutations > 0) {
            ret = tsk_site_table_add_row(sites, site->position, site->ancestral_state,
                site->ancestral_state_length, site->metadata, site->metadata_length);
            if (ret < 0) {
                goto out;
            }
            site_id++;
        }
    }
    ret = 0;
out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_place_mutations(
    mutgen_t *self, tsk_table_collection_t *tables, bool discrete_sites)
{
    /* The mutation model for discrete sites is that there is
     * a unit of "mutation mass" on each integer, so that
     * in a segment [left, right) there can be mutations at the
     * integers {ceil(left), ceil(left) + 1, ..., ceil(right) - 1},
     * and the total mutation "length" is ceil(right) - ceil(left). */
    int ret = 0;
    const double *map_position = self->rate_map->position;
    const double *map_rate = self->rate_map->value;
    size_t branch_mutations, map_index;
    size_t j, k;
    tsk_node_table_t *nodes = &tables->nodes;
    tsk_edge_table_t *edges = &tables->edges;
    double left, right, site_left, site_right, edge_right;
    double time, mu, position;
    double branch_start, branch_end, branch_length;
    node_id_t parent, child;
    avl_node_t *avl_node;
    double start_time = self->start_time;
    double end_time = self->end_time;
    site_t search;

    for (j = 0; j < edges->num_rows; j++) {
        left = edges->left[j];
        edge_right = edges->right[j];
        parent = edges->parent[j];
        child = edges->child[j];
        assert(child >= 0 && child < (node_id_t) nodes->num_rows);
        branch_start = GSL_MAX(start_time, nodes->time[child]);
        branch_end = GSL_MIN(end_time, nodes->time[parent]);
        branch_length = branch_end - branch_start;

        map_index = interval_map_get_index(self->rate_map, left);
        right = 0;
        while (right != edge_right) {
            right = GSL_MIN(edge_right, map_position[map_index + 1]);
            site_left = discrete_sites ? ceil(left) : left;
            site_right = discrete_sites ? ceil(right) : right;
            mu = branch_length * (site_right - site_left) * map_rate[map_index];
            branch_mutations = gsl_ran_poisson(self->rng, mu);
            for (k = 0; k < branch_mutations; k++) {
                /* Rejection sample positions until we get one we haven't seen before,
                 * unless we are doing discrete sites. Note that in principle this
                 * could lead to an infinite loop here, but in practise we'd need to
                 * use up all of the doubles before it could happen and so we'd
                 * certainly run out of memory first. */
                do {
                    position = gsl_ran_flat(self->rng, site_left, site_right);
                    if (discrete_sites) {
                        position = floor(position);
                    }
                    search.position = position;
                    avl_node = avl_search(&self->sites, &search);
                } while (avl_node != NULL && !discrete_sites);

                time = gsl_ran_flat(self->rng, branch_start, branch_end);
                assert(site_left <= position && position < site_right);
                assert(branch_start <= time && time < branch_end);
                if (avl_node == NULL) {
                    ret = mutgen_add_site(
                        self, position, true, NULL, 0, NULL, 0, &avl_node);
                    if (ret != 0) {
                        goto out;
                    }
                }
                ret = mutgen_add_mutation(self, (site_t *) avl_node->item, TSK_NULL,
                    child, time, true, NULL, 0, NULL, 0);
                if (ret != 0) {
                    goto out;
                }
            }
            map_index++;
        }
    }
out:
    return ret;
}

static int
mutation_model_pick_allele(mutation_model_t *self, gsl_rng *rng)
{
    double u = gsl_ran_flat(rng, 0.0, 1.0);
    int j = 0;
    while (u > self->root_distribution[j]) {
        u -= self->root_distribution[j];
        j++;
        assert(j < (int) self->num_alleles);
    }
    return j;
}

static int
mutation_model_transition_allele(mutation_model_t *self, int parent_allele, gsl_rng *rng)
{
    double *probs = self->transition_matrix + parent_allele * (int) self->num_alleles;
    double u = gsl_ran_flat(rng, 0.0, 1.0);
    int j = 0;

    while (u > probs[j]) {
        u -= probs[j];
        j++;
        assert(j < (int) self->num_alleles);
    }
    return j;
}

static int
mutgen_set_ancestral_state(mutgen_t *self, site_t *site, int allele)
{
    int ret = 0;
    size_t len = strlen(self->model->alleles[allele]);
    char *buff = tsk_blkalloc_get(&self->allocator, len);

    if (buff == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memcpy(buff, self->model->alleles[allele], len);
    site->ancestral_state = buff;
    site->ancestral_state_length = (tsk_size_t) len;
out:
    return ret;
}

static int
mutgen_set_derived_state(mutgen_t *self, mutation_t *mut, int allele)
{
    int ret = 0;
    size_t len = strlen(self->model->alleles[allele]);
    char *buff = tsk_blkalloc_get(&self->allocator, len);

    if (buff == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memcpy(buff, self->model->alleles[allele], len);
    mut->derived_state = buff;
    mut->derived_state_length = (tsk_size_t) len;
out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_choose_alleles(mutgen_t *self, tsk_id_t *parent, mutation_t **bottom_mutation,
    tsk_size_t num_nodes, site_t *site)
{
    int ret = 0;
    int ai, pa, da;
    mutation_t *mut, *parent_mut;
    tsk_id_t u;

    ret = sort_mutations(site);
    if (ret != 0) {
        goto out;
    }

    if (site->new) {
        ai = mutation_model_pick_allele(self->model, self->rng);
        ret = mutgen_set_ancestral_state(self, site, ai);
        if (ret != 0) {
            goto out;
        }
    } else {
        /* we don't error here if the allele is not found (and so ai < 0)
         * because it's ok if we don't try to add a mutation to the site,
         * e.g., if it is kept and we're doing infinite sites, or if we're
         * mutating a different part of the genome  */
        ai = mutation_model_allele_index(
            self->model, site->ancestral_state, site->ancestral_state_length);
    }

    /* Create a mapping from mutations to nodes in bottom_mutation. If we see
     * more than one mutation at a node, the previously seen one must be the
     * parent of the current one since we assume they are in order. */
    for (mut = site->mutations; mut != NULL; mut = mut->next) {
        u = mut->node;
        assert((tsk_size_t) u < num_nodes);
        while (u != TSK_NULL && bottom_mutation[u] == NULL) {
            u = parent[u];
        }
        if (u == TSK_NULL) {
            pa = ai;
            assert(mut->parent == NULL);
        } else {
            parent_mut = bottom_mutation[u];
            mut->parent = parent_mut;
            if (mut->time > parent_mut->time || (parent_mut->new && !mut->new)) {
                ret = MSP_ERR_MUTATION_GENERATION_OUT_OF_ORDER;
                goto out;
            }
            if (mut->new) {
                pa = mutation_model_allele_index(self->model, parent_mut->derived_state,
                    parent_mut->derived_state_length);
            }
        }
        if (mut->new) {
            if (pa < 0) {
                /* only error if we are actually trying to mutate an unknown allele */
                ret = MSP_ERR_UNKNOWN_ALLELE;
                goto out;
            }
            da = mutation_model_transition_allele(self->model, pa, self->rng);
            if (da == pa) {
                // mark mut for removal
                mut->keep = false;
            } else {
                ret = mutgen_set_derived_state(self, mut, da);
                if (ret != 0) {
                    goto out;
                }
            }
        }
        /* mut->keep defaults to true */
        if (mut->keep) {
            bottom_mutation[mut->node] = mut;
        }
    }
    /* Reset the mapping for the next site */
    for (mut = site->mutations; mut != NULL; mut = mut->next) {
        bottom_mutation[mut->node] = NULL;
    }

out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_apply_mutations(mutgen_t *self, tsk_table_collection_t *tables)
{
    int ret = 0;
    const tsk_id_t *I, *O;
    const tsk_edge_table_t edges = tables->edges;
    const tsk_node_table_t nodes = tables->nodes;
    const tsk_id_t M = (tsk_id_t) edges.num_rows;
    tsk_id_t tj, tk;
    tsk_id_t *parent = NULL;
    mutation_t **bottom_mutation = NULL;
    double left, right;
    avl_node_t *avl_node;
    site_t *site;

    parent = malloc(nodes.num_rows * sizeof(*parent));
    bottom_mutation = malloc(nodes.num_rows * sizeof(*bottom_mutation));
    if (parent == NULL || bottom_mutation == NULL) {
        ret = TSK_ERR_NO_MEMORY;
        goto out;
    }
    memset(parent, 0xff, nodes.num_rows * sizeof(*parent));
    memset(bottom_mutation, 0, nodes.num_rows * sizeof(*bottom_mutation));

    if (!tsk_table_collection_has_index(tables, 0)) {
        ret = tsk_table_collection_build_index(tables, 0);
        if (ret != 0) {
            goto out;
        }
    }

    I = tables->indexes.edge_insertion_order;
    O = tables->indexes.edge_removal_order;
    tj = 0;
    tk = 0;
    left = 0;
    avl_node = self->sites.head;
    while (tj < M || left < tables->sequence_length) {
        while (tk < M && edges.right[O[tk]] == left) {
            parent[edges.child[O[tk]]] = TSK_NULL;
            tk++;
        }
        while (tj < M && edges.left[I[tj]] == left) {
            parent[edges.child[I[tj]]] = edges.parent[I[tj]];
            tj++;
        }
        right = tables->sequence_length;
        if (tj < M) {
            right = TSK_MIN(right, edges.left[I[tj]]);
        }
        if (tk < M) {
            right = TSK_MIN(right, edges.right[O[tk]]);
        }

        /* Tree is now ready. We look at each site on this tree in turn */
        while (avl_node != NULL) {
            site = (site_t *) avl_node->item;
            if (site->position >= right) {
                break;
            }
            ret = mutgen_choose_alleles(
                self, parent, bottom_mutation, nodes.num_rows, site);
            if (ret != 0) {
                goto out;
            }
            avl_node = avl_node->next;
        }
        /* Move on to the next tree */
        left = right;
    }

out:
    msp_safe_free(parent);
    msp_safe_free(bottom_mutation);
    return ret;
}

int MSP_WARN_UNUSED
mutgen_generate(mutgen_t *self, tsk_table_collection_t *tables, int flags)
{
    int ret = 0;
    bool discrete_sites = flags & MSP_DISCRETE_SITES;

    avl_clear_tree(&self->sites);

    ret = mutgen_init_allocator(self, tables);
    if (ret != 0) {
        goto out;
    }
    ret = tsk_table_collection_check_integrity(tables, TSK_CHECK_OFFSETS);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
    if (interval_map_get_sequence_length(self->rate_map) != tables->sequence_length) {
        ret = MSP_ERR_INCOMPATIBLE_MUTATION_MAP;
        goto out;
    }
    if (flags & MSP_KEEP_SITES) {
        ret = mutgen_initialise_sites(self, tables);
        if (ret != 0) {
            goto out;
        }
    }

    ret = tsk_site_table_clear(&tables->sites);
    if (ret != 0) {
        goto out;
    }
    ret = tsk_mutation_table_clear(&tables->mutations);
    if (ret != 0) {
        goto out;
    }
    ret = mutgen_place_mutations(self, tables, discrete_sites);
    if (ret != 0) {
        goto out;
    }
    ret = mutgen_apply_mutations(self, tables);
    if (ret != 0) {
        goto out;
    }
    ret = mutgen_populate_tables(self, &tables->sites, &tables->mutations);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}
