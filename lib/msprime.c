/*
** Copyright (C) 2015-2021 University of Oxford
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
#include <float.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_sf.h>

#include "msprime.h"

/* State machine for the simulator object. */
#define MSP_STATE_NEW 0
#define MSP_STATE_INITIALISED 1
#define MSP_STATE_SIMULATING 2
#define MSP_STATE_DEBUGGING 3

/* Draw a random variable from a truncated Beta(a, b) distribution,
 * by rejecting draws above the truncation point x.
 */
static double
ran_inc_beta_rej(gsl_rng *r, double a, double b, double x)
{
    double ret;
    do {
        ret = gsl_ran_beta(r, a, b);
    } while (ret > x);
    return ret;
}

/* Draw a random variable from a truncated Beta(a, b) distribution,
 * using the inverse transform sampling method. A uniform variable
 * is drawn and then transformed into a Beta(a, b) distribution using
 * the inverse CDF. Truncation is obtained by defining the upper bound
 * of the uniform variable using the incomplete beta function. I.e.
 *      upper_bound = gsl_sf_beta_inc(a, b, x),
 * where x is the truncation point.
 */
static double
ran_inc_beta_its(gsl_rng *r, double a, double b, double upper_bound)
{
    double u = gsl_ran_flat(r, 0, upper_bound);
    return gsl_cdf_beta_Pinv(u, a, b);
}

/* Draw a random variable from a truncated Beta(a, b) distribution,
 * using either the rejection method or the inverse transform sampling method.
 * The rejection method is quicker for a moderate number of rejections, but if
 * the acceptance probability is low, we instead choose inverse transform
 * sampling. The acceptance probability for ran_inc_beta_rej() is equivalent to
 * the upper bound used with ran_inc_beta_its().
 */
static double
ran_inc_beta(gsl_rng *r, double a, double b, double x)
{
    double ub = gsl_sf_beta_inc(a, b, x);
    if (ub < 0.1) {
        return ran_inc_beta_its(r, a, b, ub);
    } else {
        return ran_inc_beta_rej(r, a, b, x);
    }
}

/* Returns the size of the specified population at the specified time */
static double
get_population_size(population_t *pop, double t)
{
    double ret = 0;
    double alpha = pop->growth_rate;
    double dt;

    if (alpha == 0.0) {
        ret = pop->initial_size;
    } else {
        dt = t - pop->start_time;
        ret = pop->initial_size * exp(-alpha * dt);
    }
    return ret;
}

static int
cmp_individual(const void *a, const void *b)
{
    const segment_t *ia = (const segment_t *) a;
    const segment_t *ib = (const segment_t *) b;
    return (ia->id > ib->id) - (ia->id < ib->id);
}

/* For the segment priority queue we want to sort on the left
 * coordinate and to break ties we arbitrarily use the ID */
static int
cmp_segment_queue(const void *a, const void *b)
{
    const segment_t *ia = (const segment_t *) a;
    const segment_t *ib = (const segment_t *) b;
    int ret = (ia->left > ib->left) - (ia->left < ib->left);
    if (ret == 0) {
        ret = (ia->id > ib->id) - (ia->id < ib->id);
    }
    return ret;
}

static int
cmp_node_mapping(const void *a, const void *b)
{
    const node_mapping_t *ia = (const node_mapping_t *) a;
    const node_mapping_t *ib = (const node_mapping_t *) b;
    return (ia->position > ib->position) - (ia->position < ib->position);
}

static int
cmp_sampling_event(const void *a, const void *b)
{
    const sampling_event_t *ia = (const sampling_event_t *) a;
    const sampling_event_t *ib = (const sampling_event_t *) b;
    return (ia->time > ib->time) - (ia->time < ib->time);
}

static int
cmp_pointer(const void *a, const void *b)
{
    const intptr_t ia = (const intptr_t) a;
    const intptr_t ib = (const intptr_t) b;
    return (ia > ib) - (ia < ib);
}

static int
cmp_pedigree_individual_p(const void *a, const void *b)
{
    /* the extra cast is to avoid a gcc bug in -Werror=cast-qual:
     * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=81631 */
    const individual_t *ia = *(const individual_t **) (uintptr_t) a;
    const individual_t *ib = *(const individual_t **) (uintptr_t) b;
    int ret = (ia->time > ib->time) - (ia->time < ib->time);

    if (ret == 0) {
        ret = (ia->id > ib->id) - (ia->id < ib->id);
    }
    return ret;
}

static void
segment_init(void **obj, size_t id)
{
    segment_t *seg = (segment_t *) obj;
    seg->id = id + 1;
}

size_t
msp_get_num_avl_node_blocks(msp_t *self)
{
    return self->avl_node_heap.num_blocks;
}

size_t
msp_get_num_node_mapping_blocks(msp_t *self)
{
    return self->node_mapping_heap.num_blocks;
}

size_t
msp_get_num_segment_blocks(msp_t *self)
{
    uint32_t j;
    size_t total = 0;
    for (j = 0; j < self->num_labels; j++) {
        total += self->segment_heap[j].num_blocks;
    }
    return total;
}

size_t
msp_get_num_common_ancestor_events(msp_t *self)
{
    return self->num_ca_events;
}

size_t
msp_get_num_rejected_common_ancestor_events(msp_t *self)
{
    return self->num_rejected_ca_events;
}

size_t
msp_get_num_recombination_events(msp_t *self)
{
    return self->num_re_events;
}

size_t
msp_get_num_gene_conversion_events(msp_t *self)
{
    return self->num_gc_events;
}

size_t
msp_get_num_internal_gene_conversion_events(msp_t *self)
{
    return self->num_internal_gc_events;
}

size_t
msp_get_num_noneffective_gene_conversion_events(msp_t *self)
{
    return self->num_noneffective_gc_events;
}

double
msp_get_sum_internal_gc_tract_lengths(msp_t *self)
{
    return self->sum_internal_gc_tract_lengths;
}

int
msp_set_start_time(msp_t *self, double start_time)
{
    int ret = 0;

    if (!isfinite(start_time)) {
        ret = MSP_ERR_BAD_START_TIME;
        goto out;
    }
    self->start_time = start_time;
out:
    return ret;
}

int
msp_set_store_migrations(msp_t *self, bool store_migrations)
{
    self->store_migrations = store_migrations;
    return 0;
}

int
msp_set_store_full_arg(msp_t *self, bool store_full_arg)
{
    self->store_full_arg = store_full_arg;
    return 0;
}

int
msp_set_ploidy(msp_t *self, int ploidy)
{
    int ret = 0;
    if (ploidy < 1) {
        ret = MSP_ERR_BAD_PLOIDY;
        goto out;
    }
    self->ploidy = (uint32_t) ploidy;
out:
    return ret;
}

int
msp_set_discrete_genome(msp_t *self, bool is_discrete)
{
    self->discrete_genome = is_discrete;
    return 0;
}

int
msp_set_recombination_map(msp_t *self, size_t size, double *position, double *rate)
{
    int ret = 0;

    rate_map_free(&self->recomb_map);

    ret = rate_map_alloc(&self->recomb_map, size, position, rate);
    if (ret != 0) {
        goto out;
    }
    if (rate_map_get_sequence_length(&self->recomb_map) != self->sequence_length) {
        ret = MSP_ERR_BAD_RATE_MAP;
        goto out;
    }
out:
    return ret;
}

/* Short-cut for msp_set_recombination_map can be used in testing. */
int
msp_set_recombination_rate(msp_t *self, double rate)
{
    double position[] = { 0, self->sequence_length };
    return msp_set_recombination_map(self, 1, position, &rate);
}

int
msp_set_gene_conversion_tract_length(msp_t *self, double tract_length)
{
    int ret = 0;

    if ((self->discrete_genome && tract_length < 1) || tract_length <= 0
        || tract_length > self->sequence_length) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->gc_tract_length = tract_length;

out:
    return ret;
}

int
msp_set_gene_conversion_map(msp_t *self, size_t size, double *position, double *rate)
{
    int ret = 0;

    rate_map_free(&self->gc_map);

    ret = rate_map_alloc(&self->gc_map, size, position, rate);
    if (ret != 0) {
        goto out;
    }
    if (rate_map_get_sequence_length(&self->gc_map) != self->sequence_length) {
        ret = MSP_ERR_BAD_RATE_MAP;
        goto out;
    }
out:
    return ret;
}

/* Short-cut for msp_set_gene_conversion_map can be used in testing. */
int
msp_set_gene_conversion_rate(msp_t *self, double rate)
{
    double position[] = { 0, self->sequence_length };
    return msp_set_gene_conversion_map(self, 1, position, &rate);
}

static inline double
msp_get_recomb_left_bound(msp_t *self, segment_t *seg)
{
    double left_bound;
    if (seg->prev == NULL) {
        left_bound = self->discrete_genome ? seg->left + 1 : seg->left;
    } else {
        left_bound = seg->prev->right;
    }
    return left_bound;
}

static inline double
msp_get_gc_left_bound(msp_t *self, segment_t *seg)
{
    return msp_get_recomb_left_bound(self, seg);
}

/* Set the mass of the specified segment to that between the segment's right endpoint
 * and the right endpoint of the left tail segment.
 */
static void
msp_set_segment_mass(msp_t *self, segment_t *seg)
{
    double left_bound, mass;

    if (self->recomb_mass_index != NULL) {
        left_bound = msp_get_recomb_left_bound(self, seg);
        mass = rate_map_mass_between(&self->recomb_map, left_bound, seg->right);
        fenwick_set_value(&self->recomb_mass_index[seg->label], seg->id, mass);
    }
    if (self->gc_mass_index != NULL) {
        /* NOTE: it looks like the gc_left_bound doesn't actually give us the
         * right distribution of gc events, so we'll probably get rid of this
         * and use the same left bound for both. */
        left_bound = msp_get_gc_left_bound(self, seg);
        mass = rate_map_mass_between(&self->gc_map, left_bound, seg->right);
        fenwick_set_value(&self->gc_mass_index[seg->label], seg->id, mass);
    }
}

/* Add all extant segments into the indexes. */
static void
msp_reindex_segments(msp_t *self)
{
    avl_node_t *node;
    avl_tree_t *population_ancestors;
    segment_t *seg;
    size_t j;
    label_id_t label;

    for (j = 0; j < self->num_populations; j++) {
        for (label = 0; label < (label_id_t) self->num_labels; label++) {
            population_ancestors = &self->populations[j].ancestors[label];
            for (node = population_ancestors->head; node != NULL; node = node->next) {
                for (seg = (segment_t *) node->item; seg != NULL; seg = seg->next) {
                    msp_set_segment_mass(self, seg);
                }
            }
        }
    }
}

/* Setup the mass indexes either after a simulation model change
 * or during msp_initialise */
static int
msp_setup_mass_indexes(msp_t *self)
{
    int ret = 0;
    label_id_t label;
    size_t num_segments;
    bool build_recomb_mass_index, build_gc_mass_index;

    /* For simplicity, we always drop the mass indexes even though
     * sometimes we'll be dropping it just to rebuild */
    if (self->recomb_mass_index != NULL) {
        for (label = 0; label < (label_id_t) self->num_labels; label++) {
            fenwick_free(&self->recomb_mass_index[label]);
        }
        msp_safe_free(self->recomb_mass_index);
        self->recomb_mass_index = NULL;
    }
    if (self->gc_mass_index != NULL) {
        for (label = 0; label < (label_id_t) self->num_labels; label++) {
            fenwick_free(&self->gc_mass_index[label]);
        }
        msp_safe_free(self->gc_mass_index);
        self->gc_mass_index = NULL;
    }

    /* We never build indexes for the DTWF and Pedigree models */
    if (self->model.type == MSP_MODEL_DTWF || self->model.type == MSP_MODEL_WF_PED) {
        build_recomb_mass_index = false;
        build_gc_mass_index = false;
    } else {
        /* For all the other models, we maintain an index only if the total rate
         * > 0. */
        build_recomb_mass_index = rate_map_get_total_mass(&self->recomb_map) > 0;
        build_gc_mass_index = rate_map_get_total_mass(&self->gc_map) > 0;
    }

    num_segments = self->segment_heap->size;
    if (build_recomb_mass_index) {
        self->recomb_mass_index
            = calloc(self->num_labels, sizeof(*self->recomb_mass_index));
        if (self->recomb_mass_index == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        for (label = 0; label < (label_id_t) self->num_labels; label++) {
            ret = fenwick_alloc(&self->recomb_mass_index[label], num_segments);
            if (ret != 0) {
                goto out;
            }
        }
    }
    if (build_gc_mass_index) {
        self->gc_mass_index = calloc(self->num_labels, sizeof(*self->gc_mass_index));
        if (self->gc_mass_index == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        for (label = 0; label < (label_id_t) self->num_labels; label++) {
            ret = fenwick_alloc(&self->gc_mass_index[label], num_segments);
            if (ret != 0) {
                goto out;
            }
        }
    }

    msp_reindex_segments(self);
out:
    return ret;
}

static int
msp_alloc_populations(msp_t *self)
{
    int ret = 0;
    size_t j;
    size_t N = self->num_populations * self->num_populations;

    self->initial_migration_matrix = calloc(N, sizeof(*self->initial_migration_matrix));
    self->migration_matrix = calloc(N, sizeof(*self->migration_matrix));
    self->num_migration_events = calloc(N, sizeof(*self->num_migration_events));
    self->initial_populations
        = calloc(self->num_populations, sizeof(*self->initial_populations));
    self->populations = calloc(self->num_populations, sizeof(*self->populations));

    if (self->migration_matrix == NULL || self->initial_migration_matrix == NULL
        || self->num_migration_events == NULL || self->initial_populations == NULL
        || self->populations == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (j = 0; j < self->num_populations; j++) {
        self->populations[j].potential_destinations
            = malloc(self->num_populations
                     * sizeof(*self->populations[j].potential_destinations));
        if (self->populations[j].potential_destinations == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        /* Set the default sizes and growth rates. */
        self->initial_populations[j].growth_rate = 0.0;
        self->initial_populations[j].initial_size = 1.0;
        self->initial_populations[j].start_time = 0.0;
        self->initial_populations[j].state = MSP_POP_STATE_ACTIVE;
    }
out:
    return ret;
}

int
msp_set_num_labels(msp_t *self, size_t num_labels)
{
    int ret = 0;
    size_t j, k;

    if (num_labels < 1 || num_labels > UINT32_MAX) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }

    /* Free any memory, if it has been allocated */
    for (j = 0; j < self->num_populations; j++) {
        msp_safe_free(self->populations[j].ancestors);
    }
    msp_safe_free(self->segment_heap);

    self->num_labels = (uint32_t) num_labels;
    self->segment_heap = calloc(self->num_labels, sizeof(*self->segment_heap));
    if (self->segment_heap == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    for (j = 0; j < self->num_populations; j++) {
        self->populations[j].ancestors
            = malloc(self->num_labels * sizeof(*self->populations[j].ancestors));
        if (self->populations[j].ancestors == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        for (k = 0; k < num_labels; k++) {
            avl_init_tree(&self->populations[j].ancestors[k], cmp_individual, NULL);
        }
    }
out:
    return ret;
}

int
msp_set_population_configuration(msp_t *self, int population_id, double initial_size,
    double growth_rate, bool initially_active)
{
    int ret = MSP_ERR_BAD_POPULATION_CONFIGURATION;

    if (population_id < 0 || population_id >= (int) self->num_populations) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    if (initial_size < 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->initial_populations[population_id].initial_size = initial_size;
    self->initial_populations[population_id].growth_rate = growth_rate;
    self->initial_populations[population_id].state
        = initially_active ? MSP_POP_STATE_ACTIVE : MSP_POP_STATE_INACTIVE;
    ret = 0;
out:
    return ret;
}

int
msp_set_migration_matrix(msp_t *self, size_t size, double *migration_matrix)
{
    int ret = MSP_ERR_BAD_MIGRATION_MATRIX;
    size_t j, k;
    size_t N = self->num_populations;

    if (N * N != size) {
        goto out;
    }
    /* Check values */
    for (j = 0; j < N; j++) {
        for (k = 0; k < N; k++) {
            if (j == k) {
                if (migration_matrix[j * N + k] != 0.0) {
                    goto out;
                }
            } else {
                if (migration_matrix[j * N + k] < 0.0) {
                    goto out;
                }
            }
        }
    }
    for (j = 0; j < N * N; j++) {
        self->initial_migration_matrix[j] = migration_matrix[j];
    }
    ret = 0;
out:
    return ret;
}

int
msp_set_node_mapping_block_size(msp_t *self, size_t block_size)
{
    int ret = 0;

    if (block_size < 1) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->node_mapping_block_size = block_size;
out:
    return ret;
}

int
msp_set_segment_block_size(msp_t *self, size_t block_size)
{
    int ret = 0;

    if (block_size < 1) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->segment_block_size = block_size;
out:
    return ret;
}

int
msp_set_avl_node_block_size(msp_t *self, size_t block_size)
{
    int ret = 0;

    if (block_size < 1) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->avl_node_block_size = block_size;
out:
    return ret;
}

static segment_t *MSP_WARN_UNUSED
msp_alloc_segment(msp_t *self, double left, double right, tsk_id_t value,
    population_id_t population, label_id_t label, segment_t *prev, segment_t *next)
{
    segment_t *seg = NULL;

    if (object_heap_empty(&self->segment_heap[label])) {
        if (object_heap_expand(&self->segment_heap[label]) != 0) {
            goto out;
        }
        if (self->recomb_mass_index != NULL) {
            if (fenwick_expand(&self->recomb_mass_index[label], self->segment_block_size)
                != 0) {
                goto out;
            }
        }
        if (self->gc_mass_index != NULL) {
            if (fenwick_expand(&self->gc_mass_index[label], self->segment_block_size)
                != 0) {
                goto out;
            }
        }
    }
    seg = (segment_t *) object_heap_alloc_object(&self->segment_heap[label]);
    if (seg == NULL) {
        goto out;
    }
    tsk_bug_assert(left < right);
    if (self->recomb_mass_index != NULL) {
        tsk_bug_assert(fenwick_get_value(&self->recomb_mass_index[label], seg->id) == 0);
    }
    if (self->gc_mass_index != NULL) {
        tsk_bug_assert(fenwick_get_value(&self->gc_mass_index[label], seg->id) == 0);
    }
    seg->prev = prev;
    seg->next = next;
    seg->left = left;
    seg->right = right;
    seg->value = value;
    seg->population = population;
    seg->label = label;
out:
    return seg;
}

static segment_t *MSP_WARN_UNUSED
msp_copy_segment(msp_t *self, const segment_t *seg)
{
    return msp_alloc_segment(self, seg->left, seg->right, seg->value, seg->population,
        seg->label, seg->prev, seg->next);
}

/* Top level allocators and initialisation */

int
msp_alloc(msp_t *self, tsk_table_collection_t *tables, gsl_rng *rng)
{
    int ret = -1;

    memset(self, 0, sizeof(msp_t));
    if (rng == NULL || tables == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }

    /* Use the standard coalescent by default. */
    self->model.type = -1;
    ret = msp_set_simulation_model_hudson(self);
    tsk_bug_assert(ret == 0);
    self->rng = rng;
    self->discrete_genome = true;

    self->tables = tables;
    self->sequence_length = tables->sequence_length;
    if (self->sequence_length <= 0) {
        ret = MSP_ERR_BAD_SEQUENCE_LENGTH;
        goto out;
    }
    self->num_populations = (uint32_t) self->tables->populations.num_rows;
    if (self->num_populations == 0) {
        ret = MSP_ERR_ZERO_POPULATIONS;
        goto out;
    }
    ret = msp_set_recombination_rate(self, 0.0);
    if (ret != 0) {
        goto out;
    }
    ret = msp_set_gene_conversion_rate(self, 0.0);
    if (ret != 0) {
        goto out;
    }
    ret = msp_alloc_populations(self);
    if (ret != 0) {
        goto out;
    }
    ret = msp_set_num_labels(self, 1);
    if (ret != 0) {
        goto out;
    }

    /* If the start_time is not set, we default to the minimum root time */
    self->start_time = -DBL_MAX;
    /* Set the memory defaults */
    self->store_migrations = false;
    self->store_full_arg = false;
    self->avl_node_block_size = 1024;
    self->node_mapping_block_size = 1024;
    self->segment_block_size = 1024;
    /* set up the AVL trees */
    avl_init_tree(&self->breakpoints, cmp_node_mapping, NULL);
    avl_init_tree(&self->overlap_counts, cmp_node_mapping, NULL);
    avl_init_tree(&self->non_empty_populations, cmp_pointer, NULL);
    /* Set up the demographic events */
    self->demographic_events_head = NULL;
    self->demographic_events_tail = NULL;
    self->next_demographic_event = NULL;
    self->state = MSP_STATE_NEW;
    /* Set default to diploid */
    self->ploidy = 2;
out:
    return ret;
}

static int
msp_alloc_memory_blocks(msp_t *self)
{
    int ret = 0;
    uint32_t j;

    /* Allocate the memory heaps */
    ret = object_heap_init(
        &self->avl_node_heap, sizeof(avl_node_t), self->avl_node_block_size, NULL);
    if (ret != 0) {
        goto out;
    }
    ret = object_heap_init(&self->node_mapping_heap, sizeof(node_mapping_t),
        self->node_mapping_block_size, NULL);
    if (ret != 0) {
        goto out;
    }
    /* allocate the segments */
    for (j = 0; j < self->num_labels; j++) {
        ret = object_heap_init(&self->segment_heap[j], sizeof(segment_t),
            self->segment_block_size, segment_init);
        if (ret != 0) {
            goto out;
        }
    }
    /* Allocate the edge records */
    self->num_buffered_edges = 0;
    self->max_buffered_edges = 128;
    self->buffered_edges = malloc(self->max_buffered_edges * sizeof(tsk_edge_t));
    if (self->buffered_edges == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = 0;
out:
    return ret;
}

/*
 * Returns true if the simulation has completed.
 */
int
msp_is_completed(msp_t *self)
{
    size_t n = msp_get_num_ancestors(self);

    return self->state == MSP_STATE_SIMULATING && n == 0;
}

int
msp_free(msp_t *self)
{
    int ret = -1;
    uint32_t j;

    demographic_event_t *de = self->demographic_events_head;
    demographic_event_t *tmp;

    while (de != NULL) {
        tmp = de->next;
        free(de);
        de = tmp;
    }
    for (j = 0; j < self->num_labels; j++) {
        if (self->recomb_mass_index != NULL) {
            fenwick_free(&self->recomb_mass_index[j]);
        }
        if (self->gc_mass_index != NULL) {
            fenwick_free(&self->gc_mass_index[j]);
        }
        if (self->segment_heap != NULL) {
            object_heap_free(&self->segment_heap[j]);
        }
    }
    for (j = 0; j < self->num_populations; j++) {
        msp_safe_free(self->populations[j].ancestors);
        msp_safe_free(self->populations[j].potential_destinations);
    }
    msp_safe_free(self->recomb_mass_index);
    msp_safe_free(self->gc_mass_index);
    msp_safe_free(self->segment_heap);
    msp_safe_free(self->initial_migration_matrix);
    msp_safe_free(self->migration_matrix);
    msp_safe_free(self->num_migration_events);
    msp_safe_free(self->initial_populations);
    msp_safe_free(self->populations);
    msp_safe_free(self->sampling_events);
    msp_safe_free(self->buffered_edges);
    msp_safe_free(self->root_segments);
    msp_safe_free(self->initial_overlaps);
    msp_safe_free(self->pedigree.individuals);
    msp_safe_free(self->pedigree.visit_order);
    /* free the object heaps */
    object_heap_free(&self->avl_node_heap);
    object_heap_free(&self->node_mapping_heap);
    rate_map_free(&self->recomb_map);
    rate_map_free(&self->gc_map);
    if (self->model.free != NULL) {
        self->model.free(&self->model);
    }
    ret = 0;
    return ret;
}

static inline avl_node_t *MSP_WARN_UNUSED
msp_alloc_avl_node(msp_t *self)
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
msp_free_avl_node(msp_t *self, avl_node_t *node)
{
    object_heap_free_object(&self->avl_node_heap, node);
}

static inline node_mapping_t *
msp_alloc_node_mapping(msp_t *self)
{
    node_mapping_t *ret = NULL;

    if (object_heap_empty(&self->node_mapping_heap)) {
        if (object_heap_expand(&self->node_mapping_heap) != 0) {
            goto out;
        }
    }
    ret = (node_mapping_t *) object_heap_alloc_object(&self->node_mapping_heap);
out:
    return ret;
}

static void
msp_free_node_mapping(msp_t *self, node_mapping_t *nm)
{
    object_heap_free_object(&self->node_mapping_heap, nm);
}

/*
 * Returns the segment with the specified id.
 */
static segment_t *
msp_get_segment(msp_t *self, size_t id, label_id_t label)
{
    segment_t *u = object_heap_get_object(&self->segment_heap[label], id - 1);

    tsk_bug_assert(u != NULL);
    tsk_bug_assert(u->id == id);
    return u;
}

static void
msp_free_segment(msp_t *self, segment_t *seg)
{
    object_heap_free_object(&self->segment_heap[seg->label], seg);
    if (self->recomb_mass_index != NULL) {
        fenwick_set_value(&self->recomb_mass_index[seg->label], seg->id, 0);
    }
    if (self->gc_mass_index != NULL) {
        fenwick_set_value(&self->gc_mass_index[seg->label], seg->id, 0);
    }
}

static inline avl_tree_t *
msp_get_segment_population(msp_t *self, segment_t *u)
{
    return &self->populations[u->population].ancestors[u->label];
}

static inline int MSP_WARN_UNUSED
msp_insert_individual(msp_t *self, segment_t *u)
{
    int ret = 0;
    avl_node_t *node;

    tsk_bug_assert(u != NULL);
    node = msp_alloc_avl_node(self);
    if (node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    avl_init_node(node, u);
    node = avl_insert_node(msp_get_segment_population(self, u), node);
    tsk_bug_assert(node != NULL);
out:
    return ret;
}

static inline void
msp_remove_individual(msp_t *self, segment_t *u)
{
    avl_node_t *node;
    avl_tree_t *pop = msp_get_segment_population(self, u);

    tsk_bug_assert(u != NULL);
    node = avl_search(pop, u);
    tsk_bug_assert(node != NULL);
    avl_unlink_node(pop, node);
    msp_free_avl_node(self, node);
}

static void
msp_remove_individuals_from_population(msp_t *self, avl_tree_t *Q)
{
    avl_node_t *node;
    for (node = Q->head; node != NULL; node = node->next) {
        msp_remove_individual(self, (segment_t *) node->item);
    }
}

/* Returns true if the specified breakpoint exists */
static bool
msp_has_breakpoint(msp_t *self, double x)
{
    node_mapping_t search;
    search.position = x;

    return avl_search(&self->breakpoints, &search) != NULL;
}

/*
 * Inserts a new breakpoint at the specified locus left.
 */
static int MSP_WARN_UNUSED
msp_insert_breakpoint(msp_t *self, double left)
{
    int ret = 0;
    avl_node_t *node = msp_alloc_avl_node(self);
    node_mapping_t *m = msp_alloc_node_mapping(self);

    if (node == NULL || m == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    m->position = left;
    m->value = 0;
    avl_init_node(node, m);
    node = avl_insert_node(&self->breakpoints, node);
    tsk_bug_assert(node != NULL);
out:
    return ret;
}

static void
msp_print_segment_chain(msp_t *MSP_UNUSED(self), segment_t *head, FILE *out)
{
    segment_t *s = head;

    fprintf(out, "[pop=%d,label=%d]", s->population, s->label);
    while (s != NULL) {
        fprintf(out, "[(%.14g,%.14g) %d] ", s->left, s->right, (int) s->value);
        s = s->next;
    }
    fprintf(out, "\n");
}

/* TODO remove the left_at_zero option, it's for old GC version that didn't work. */
static void
msp_verify_segment_index(
    msp_t *self, fenwick_t *mass_index_array, rate_map_t *rate_map, bool left_at_zero)
{

    double left, right, left_bound;
    double s, ss, total_mass, alt_total_mass;
    size_t j, k;
    const double epsilon = 1e-10;
    avl_node_t *node;
    segment_t *u;

    for (k = 0; k < self->num_labels; k++) {
        total_mass = 0;
        alt_total_mass = 0;
        for (j = 0; j < self->num_populations; j++) {
            node = (&self->populations[j].ancestors[k])->head;
            while (node != NULL) {
                u = (segment_t *) node->item;
                left = u->left;
                while (u != NULL) {
                    if (u->prev != NULL) {
                        s = rate_map_mass_between(rate_map, u->prev->right, u->right);
                    } else {
                        if (left_at_zero) {
                            left_bound = self->discrete_genome ? 1 : 0;
                        } else {
                            left_bound = self->discrete_genome ? u->left + 1 : u->left;
                        }
                        tsk_bug_assert(left_bound <= u->right);
                        s = rate_map_mass_between(rate_map, left_bound, u->right);
                    }
                    tsk_bug_assert(s >= 0);
                    ss = fenwick_get_value(&mass_index_array[k], u->id);
                    tsk_bug_assert(doubles_almost_equal(s, ss, epsilon));
                    total_mass += ss;
                    right = u->right;
                    u = u->next;
                }
                if (left_at_zero) {
                    left_bound = self->discrete_genome ? 1 : 0;
                } else {
                    left_bound = self->discrete_genome ? left + 1 : left;
                }
                s = rate_map_mass_between(rate_map, left_bound, right);
                alt_total_mass += s;
                node = node->next;
            }
        }
        tsk_bug_assert(doubles_almost_equal(
            total_mass, fenwick_get_total(&mass_index_array[k]), epsilon));
        tsk_bug_assert(doubles_almost_equal(total_mass, alt_total_mass, epsilon));
    }
}

static void
msp_verify_segments(msp_t *self, bool verify_breakpoints)
{
    size_t j, k;
    size_t label_segments = 0;
    size_t total_avl_nodes = 0;
    size_t num_root_segments = 0;
    size_t pedigree_avl_nodes = 0;
    avl_node_t *node;
    segment_t *u;
    individual_t *ind;

    for (j = 0; j < self->input_position.nodes; j++) {
        for (u = self->root_segments[j]; u != NULL; u = u->next) {
            num_root_segments++;
        }
    }

    for (k = 0; k < self->num_labels; k++) {
        label_segments = 0;
        if (k == 0) {
            label_segments += num_root_segments;
        }
        for (j = 0; j < self->num_populations; j++) {
            node = (&self->populations[j].ancestors[k])->head;
            while (node != NULL) {
                u = (segment_t *) node->item;
                tsk_bug_assert(u->prev == NULL);
                while (u != NULL) {
                    label_segments++;
                    tsk_bug_assert(u->population == (population_id_t) j);
                    tsk_bug_assert(u->label == (label_id_t) k);
                    tsk_bug_assert(u->left < u->right);
                    tsk_bug_assert(u->right <= self->sequence_length);
                    if (u->prev != NULL) {
                        tsk_bug_assert(u->prev->next == u);
                    }
                    if (verify_breakpoints && u->left != 0) {
                        tsk_bug_assert(msp_has_breakpoint(self, u->left));
                    }
                    if (self->discrete_genome) {
                        tsk_bug_assert(floor(u->left) == u->left);
                    }
                    u = u->next;
                }
                node = node->next;
            }
        }
        tsk_bug_assert(
            label_segments == object_heap_get_num_allocated(&self->segment_heap[k]));
    }
    total_avl_nodes = msp_get_num_ancestors(self) + avl_count(&self->breakpoints)
                      + avl_count(&self->overlap_counts)
                      + avl_count(&self->non_empty_populations);
    for (j = 0; j < self->pedigree.num_individuals; j++) {
        ind = &self->pedigree.individuals[j];
        for (k = 0; k < self->ploidy; k++) {
            pedigree_avl_nodes += avl_count(&ind->common_ancestors[k]);
        }
    }
    tsk_bug_assert(total_avl_nodes + pedigree_avl_nodes
                   == object_heap_get_num_allocated(&self->avl_node_heap));
    tsk_bug_assert(total_avl_nodes - msp_get_num_ancestors(self)
                       - avl_count(&self->non_empty_populations)
                   == object_heap_get_num_allocated(&self->node_mapping_heap));
    if (self->recomb_mass_index != NULL) {
        msp_verify_segment_index(
            self, self->recomb_mass_index, &self->recomb_map, false);
    }
    if (self->gc_mass_index != NULL) {
        msp_verify_segment_index(self, self->gc_mass_index, &self->gc_map, false);
    }
    /* Check that the mass indexes are set appropriately */
    if (self->model.type == MSP_MODEL_DTWF || self->model.type == MSP_MODEL_WF_PED) {
        tsk_bug_assert(self->recomb_mass_index == NULL);
        tsk_bug_assert(self->gc_mass_index == NULL);
    } else {
        tsk_bug_assert((self->recomb_mass_index != NULL)
                       == (rate_map_get_total_mass(&self->recomb_map) > 0));
        tsk_bug_assert((self->gc_mass_index != NULL)
                       == (rate_map_get_total_mass(&self->gc_map) > 0));
    }
}

typedef struct {
    double seq_length;
    segment_t *overlaps;
} overlap_counter_t;

static int
overlap_counter_alloc(overlap_counter_t *self, double seq_length, int initial_count)
{
    int ret = 0;
    memset(self, 0, sizeof(overlap_counter_t));

    segment_t *overlaps = malloc(sizeof(segment_t));
    if (overlaps == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    overlaps->prev = NULL;
    overlaps->next = NULL;
    overlaps->left = 0;
    overlaps->right = seq_length;
    overlaps->value = initial_count;
    overlaps->population = 0;
    overlaps->label = 0;

    self->seq_length = seq_length;
    self->overlaps = overlaps;

out:
    return ret;
}

static void
overlap_counter_free(overlap_counter_t *self)
{
    segment_t *curr_overlap, *next_overlap;

    tsk_bug_assert(self->overlaps->prev == NULL);
    curr_overlap = self->overlaps;
    while (curr_overlap != NULL) {
        next_overlap = curr_overlap->next;
        free(curr_overlap);
        curr_overlap = next_overlap;
    }
}

/* Find the number of segments that overlap at the given position */
static uint32_t
overlap_counter_overlaps_at(overlap_counter_t *self, double pos)
{
    tsk_bug_assert(pos >= 0 && pos < self->seq_length);
    segment_t *curr_overlap = self->overlaps;
    while (curr_overlap->next != NULL) {
        if (curr_overlap->left <= pos && pos < curr_overlap->right) {
            break;
        }
        curr_overlap = curr_overlap->next;
    }

    return (uint32_t) curr_overlap->value;
}

/* Split the segment at breakpoint and add in another segment
 * from breakpoint to seg.right. Set the original segment's
 * right endpoint to breakpoint.
 */
static void
overlap_counter_split_segment(segment_t *seg, double breakpoint)
{
    segment_t *right_seg = malloc(sizeof(segment_t));
    right_seg->prev = NULL;
    right_seg->next = NULL;
    right_seg->left = breakpoint;
    right_seg->right = seg->right;
    right_seg->value = seg->value;
    right_seg->population = 0;
    right_seg->label = 0;

    if (seg->next != NULL) {
        right_seg->next = seg->next;
        seg->next->prev = right_seg;
    }
    right_seg->prev = seg;
    seg->next = right_seg;
    seg->right = breakpoint;
}

/* Increment the number of segments that span
 * [left, right), creating additional intervals if necessary.
 */
static void
overlap_counter_increment_interval(overlap_counter_t *self, double left, double right)
{
    segment_t *curr_interval = self->overlaps;
    while (left < right) {
        if (curr_interval->left == left) {
            if (curr_interval->right <= right) {
                curr_interval->value++;
                left = curr_interval->right;
                curr_interval = curr_interval->next;
            } else {
                overlap_counter_split_segment(curr_interval, right);
                curr_interval->value++;
                break;
            }
        } else {
            if (curr_interval->right < left) {
                curr_interval = curr_interval->next;
            } else {
                overlap_counter_split_segment(curr_interval, left);
                curr_interval = curr_interval->next;
            }
        }
    }
}

static void
msp_verify_overlaps(msp_t *self)
{
    avl_node_t *node;
    node_mapping_t *nm;
    sampling_event_t se;
    segment_t *u;
    size_t j;
    uint32_t label, count;
    overlap_counter_t counter;

    int ok = overlap_counter_alloc(&counter, self->sequence_length, 0);
    tsk_bug_assert(ok == 0);

    /* add in the overlaps for ancient samples */
    for (j = self->next_sampling_event; j < self->num_sampling_events; j++) {
        se = self->sampling_events[j];
        for (u = self->root_segments[se.sample]; u != NULL; u = u->next) {
            overlap_counter_increment_interval(&counter, u->left, u->right);
        }
    }

    for (label = 0; label < self->num_labels; label++) {
        for (j = 0; j < self->num_populations; j++) {
            for (node = (&self->populations[j].ancestors[label])->head; node != NULL;
                 node = node->next) {
                for (u = (segment_t *) node->item; u != NULL; u = u->next) {
                    overlap_counter_increment_interval(&counter, u->left, u->right);
                }
            }
        }
    }
    for (node = self->overlap_counts.head; node->next != NULL; node = node->next) {
        nm = (node_mapping_t *) node->item;
        count = overlap_counter_overlaps_at(&counter, nm->position);
        tsk_bug_assert(nm->value == count);
    }

    overlap_counter_free(&counter);
}

static void
msp_verify_non_empty_populations(msp_t *self)
{
    tsk_id_t j;
    void *search;
    avl_node_t *avl_node;

    for (avl_node = self->non_empty_populations.head; avl_node != NULL;
         avl_node = avl_node->next) {
        j = (tsk_id_t)(intptr_t) avl_node->item;
        tsk_bug_assert(msp_get_num_population_ancestors(self, j) > 0);
    }

    for (j = 0; j < (tsk_id_t) self->num_populations; j++) {
        search = (void *) (intptr_t) j;
        avl_node = avl_search(&self->non_empty_populations, search);
        if (msp_get_num_population_ancestors(self, j) == 0) {
            tsk_bug_assert(avl_node == NULL);
        } else {
            tsk_bug_assert(avl_node != NULL);
        }
    }
}

static void
msp_verify_migration_destinations(msp_t *self)
{
    tsk_id_t j, k, i;
    tsk_id_t N = (tsk_id_t) self->num_populations;
    double *M = self->migration_matrix;
    population_t *pop;
    bool found;

    for (j = 0; j < N; j++) {
        pop = &self->populations[j];
        for (k = 0; k < (tsk_id_t) pop->num_potential_destinations; k++) {
            tsk_bug_assert(M[j * N + pop->potential_destinations[k]] > 0);
        }
    }
    for (j = 0; j < N; j++) {
        pop = &self->populations[j];
        for (k = 0; k < N; k++) {
            found = false;
            for (i = 0; i < (tsk_id_t) pop->num_potential_destinations; i++) {
                if (pop->potential_destinations[i] == k) {
                    found = true;
                    break;
                }
            }
            tsk_bug_assert(found == (M[j * N + k] != 0));
        }
    }
}

static void
msp_verify_initial_state(msp_t *self)
{
    overlap_count_t *overlap;
    double last_overlap_left = -1;
    tsk_size_t j;
    segment_t *head, *seg, *prev;

    for (overlap = self->initial_overlaps; overlap->left < self->sequence_length;
         overlap++) {
        tsk_bug_assert(overlap->left > last_overlap_left);
        last_overlap_left = overlap->left;
    }
    /* Last overlap should be a sentinal */
    overlap->left = self->sequence_length;
    overlap->count = UINT32_MAX;

    /* First overlap should be 0 */
    tsk_bug_assert(self->initial_overlaps->left == 0);

    /* Check the root segments */
    for (j = 0; j < self->input_position.nodes; j++) {
        head = self->root_segments[j];
        if (head != NULL) {
            prev = NULL;
            for (seg = head; seg != NULL; seg = seg->next) {
                if (prev != NULL) {
                    tsk_bug_assert(prev->next == seg);
                    tsk_bug_assert(seg->prev == prev);
                    tsk_bug_assert(prev->right <= seg->left);
                }
                tsk_bug_assert(seg->left < seg->right);
                tsk_bug_assert(seg->value == (tsk_id_t) j);
                prev = seg;
            }
        }
    }
}

static void
msp_verify_pedigree(msp_t *self)
{
    pedigree_t *ped = &self->pedigree;
    individual_t *ind;
    tsk_size_t j;

    /* TODO can we check more about the pedigree here? */
    for (j = 0; j < ped->num_individuals; j++) {
        ind = &ped->individuals[j];
        tsk_bug_assert(ind->id == (tsk_id_t) j);
    }
}

void
msp_verify(msp_t *self, int options)
{
    msp_verify_initial_state(self);
    msp_verify_segments(self, options & MSP_VERIFY_BREAKPOINTS);
    msp_verify_overlaps(self);
    if (self->model.type == MSP_MODEL_HUDSON && self->state == MSP_STATE_SIMULATING) {
        msp_verify_non_empty_populations(self);
        msp_verify_migration_destinations(self);
    }
    if (self->model.type == MSP_MODEL_WF_PED) {
        msp_verify_pedigree(self);
    }
}

static void
msp_print_individual(msp_t *self, individual_t *ind, FILE *out)
{
    size_t j;

    fprintf(out, "\tID: %d; time: %f, population: %d, Parents: [", ind->id, ind->time,
        ind->population);
    for (j = 0; j < self->ploidy; j++) {
        fprintf(out, " %d", ind->parents[j]);
    }
    fprintf(out, "] common_ancestors: [");
    for (j = 0; j < self->ploidy; j++) {
        fprintf(out, " %d", avl_count(&ind->common_ancestors[j]));
    }
    fprintf(out, " ]\n");
}

static void
msp_pedigree_print_state(msp_t *self, FILE *out)
{
    individual_t *ind;
    pedigree_t *ped = &self->pedigree;
    tsk_size_t j;

    fprintf(out, "---------\n");
    fprintf(out, "Pedigree:\n");
    fprintf(out, "next_individual = %d\n", ped->next_individual);
    for (j = 0; j < ped->num_individuals; j++) {
        ind = &ped->individuals[j];
        msp_print_individual(self, ind, out);
    }
    fprintf(out, "visit order \n");
    for (j = 0; j < ped->num_individuals; j++) {
        ind = ped->visit_order[j];
        msp_print_individual(self, ind, out);
    }
    fprintf(out, "---------\n");
}

static void
msp_print_root_segments(msp_t *self, FILE *out)
{
    segment_t *seg, *head;
    tsk_size_t j;

    fprintf(out, "Root segments\n");
    for (j = 0; j < self->input_position.nodes; j++) {
        head = self->root_segments[j];
        if (head != NULL) {
            fprintf(out, "\t%d", (int) j);
            for (seg = head; seg != NULL; seg = seg->next) {
                fprintf(out, "(%f, %f)", seg->left, seg->right);
            }
            fprintf(out, "\n");
        }
    }
}

static void
msp_print_initial_overlaps(msp_t *self, FILE *out)
{
    overlap_count_t *overlap;

    fprintf(out, "Initial overlaps\n");

    for (overlap = self->initial_overlaps; overlap->left < self->sequence_length;
         overlap++) {
        fprintf(out, "\t%f -> %d\n", overlap->left, (int) overlap->count);
    }
    tsk_bug_assert(overlap->left == self->sequence_length);
    fprintf(out, "\t%f -> %d\n", overlap->left, (int) overlap->count);
}

int
msp_print_state(msp_t *self, FILE *out)
{
    int ret = 0;
    avl_node_t *a;
    node_mapping_t *nm;
    segment_t *u;
    tsk_edge_t *edge;
    demographic_event_t *de;
    sampling_event_t *se;
    double v;
    uint32_t j, k;
    segment_t **ancestors = malloc(msp_get_num_ancestors(self) * sizeof(segment_t *));

    if (ancestors == NULL && msp_get_num_ancestors(self) != 0) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = msp_get_ancestors(self, ancestors);
    if (ret != 0) {
        goto out;
    }
    fprintf(out, "simulation model      = '%s'\n", msp_get_model_name(self));
    if (self->model.type == MSP_MODEL_BETA) {
        fprintf(out, "\tbeta coalescent parameters: alpha = %f, truncation_point = %f\n",
            self->model.params.beta_coalescent.alpha,
            self->model.params.beta_coalescent.truncation_point);
    } else if (self->model.type == MSP_MODEL_DIRAC) {
        fprintf(out, "\tdirac coalescent parameters: psi = %f, c = %f\n",
            self->model.params.dirac_coalescent.psi,
            self->model.params.dirac_coalescent.c);
    } else if (self->model.type == MSP_MODEL_SWEEP) {
        fprintf(out, "\tsweep @ locus = %f\n", self->model.params.sweep.position);
        self->model.params.sweep.print_state(&self->model.params.sweep, out);
    }
    fprintf(out, "L = %.14g\n", self->sequence_length);
    fprintf(out, "discrete_genome = %d\n", self->discrete_genome);
    fprintf(out, "start_time = %f\n", self->start_time);
    fprintf(out, "recombination map:\n");
    rate_map_print_state(&self->recomb_map, out);
    fprintf(out, "gene_conversion_tract_length = %f\n", self->gc_tract_length);
    fprintf(out, "gene conversion map:\n");
    rate_map_print_state(&self->gc_map, out);
    msp_pedigree_print_state(self, out);
    msp_print_root_segments(self, out);
    msp_print_initial_overlaps(self, out);
    fprintf(out, "Sampling events:\n");
    for (j = 0; j < self->num_sampling_events; j++) {
        if (j == self->next_sampling_event) {
            fprintf(out, "  ***");
        }
        se = &self->sampling_events[j];
        fprintf(out, "\t");
        fprintf(out, "%d @ %f in deme %d\n", (int) se->sample, se->time,
            (int) se->population);
    }
    fprintf(out, "Demographic events:\n");
    for (de = self->demographic_events_head; de != NULL; de = de->next) {
        if (de == self->next_demographic_event) {
            fprintf(out, "  ***");
        }
        fprintf(out, "\t");
        de->print_state(self, de, out);
    }
    fprintf(out, "Migration matrix\n");
    for (j = 0; j < self->num_populations; j++) {
        fprintf(out, "\t");
        for (k = 0; k < self->num_populations; k++) {
            fprintf(
                out, "%0.3f ", self->migration_matrix[j * self->num_populations + k]);
        }
        fprintf(out, "\n");
    }

    fprintf(out, "Population sizes (total ancestors=%d)\n",
        (int) msp_get_num_ancestors(self));
    for (j = 0; j < self->num_labels; j++) {
        fprintf(out, "label %d\n", j);
        fprintf(out, "\trecomb_mass = %.14g\n",
            self->recomb_mass_index == NULL
                ? 0
                : fenwick_get_total(&self->recomb_mass_index[j]));
        fprintf(out, "\tgc_mass = %.14g\n",
            self->gc_mass_index == NULL ? 0
                                        : fenwick_get_total(&self->gc_mass_index[j]));
        for (k = 0; k < self->num_populations; k++) {
            fprintf(out, "\tpop_size[%d] = %d\n", k,
                avl_count(&self->populations[k].ancestors[j]));
        }
    }
    fprintf(out, "non_empty_populations = [");
    for (a = self->non_empty_populations.head; a != NULL; a = a->next) {
        j = (uint32_t)(intptr_t) a->item;
        fprintf(out, "%d,", j);
    }
    fprintf(out, "]\n");
    for (j = 0; j < self->num_populations; j++) {
        fprintf(out, "pop[%d]:\n", (int) j);
        fprintf(out, "\tstate        = %d\n", self->populations[j].state);
        fprintf(out, "\tstart_time   = %.14g\n", self->populations[j].start_time);
        fprintf(out, "\tinitial_size = %.14g\n", self->populations[j].initial_size);
        fprintf(out, "\tgrowth_rate  = %.14g\n", self->populations[j].growth_rate);
        fprintf(out, "\tcurrent_size = %.14g\n",
            get_population_size(&self->populations[j], self->time));
        fprintf(out, "\tpotential_destinations = [");
        for (k = 0; k < self->populations[j].num_potential_destinations; k++) {
            fprintf(out, "%d,", self->populations[j].potential_destinations[k]);
        }
        fprintf(out, "]\n");
    }
    fprintf(out, "Time = %f\n", self->time);
    for (j = 0; j < msp_get_num_ancestors(self); j++) {
        fprintf(out, "\t");
        msp_print_segment_chain(self, ancestors[j], out);
    }
    fprintf(out, "Fenwick trees\n");
    for (k = 0; k < self->num_labels; k++) {
        fprintf(out, "=====\nLabel %d\n=====\n", k);
        if (self->recomb_mass_index != NULL) {
            fprintf(out, "**Recomb mass**\n");
            fprintf(out, "size=%d numerical drift = %.17g\n",
                (int) fenwick_get_size(&self->recomb_mass_index[k]),
                fenwick_get_numerical_drift(&self->recomb_mass_index[k]));
            for (j = 1; j <= (uint32_t) fenwick_get_size(&self->recomb_mass_index[k]);
                 j++) {
                u = msp_get_segment(self, j, (label_id_t) k);
                v = fenwick_get_value(&self->recomb_mass_index[k], j);
                if (v != 0) {
                    fprintf(out, "\t%.14f\ti=%d l=%.14g r=%.14g v=%d prev=%p next=%p\n",
                        v, (int) u->id, u->left, u->right, (int) u->value,
                        (void *) u->prev, (void *) u->next);
                }
            }
        }
        if (self->gc_mass_index != NULL) {
            fprintf(out, "**GC mass**\n");
            fprintf(out, "numerical drift = %.17g\n",
                fenwick_get_numerical_drift(&self->gc_mass_index[k]));
            for (j = 1; j <= (uint32_t) fenwick_get_size(&self->gc_mass_index[k]); j++) {
                u = msp_get_segment(self, j, (label_id_t) k);
                v = fenwick_get_value(&self->gc_mass_index[k], j);
                if (v != 0) {
                    fprintf(out, "\t%.14f\ti=%d l=%.14g r=%.14g v=%d prev=%p next=%p\n",
                        v, (int) u->id, u->left, u->right, (int) u->value,
                        (void *) u->prev, (void *) u->next);
                }
            }
        }
    }
    fprintf(out, "Breakpoints = %d\n", avl_count(&self->breakpoints));
    for (a = self->breakpoints.head; a != NULL; a = a->next) {
        nm = (node_mapping_t *) a->item;
        fprintf(out, "\t%.14g -> %d\n", nm->position, (int) nm->value);
    }
    fprintf(out, "Overlap count = %d\n", avl_count(&self->overlap_counts));
    for (a = self->overlap_counts.head; a != NULL; a = a->next) {
        nm = (node_mapping_t *) a->item;
        fprintf(out, "\t%.14g -> %d\n", nm->position, (int) nm->value);
    }
    fprintf(out, "Tables = \n");
    tsk_table_collection_print_state(self->tables, out);

    fprintf(out, "Buffered Edges = %ld\n", (long) self->num_buffered_edges);
    for (j = 0; j < self->num_buffered_edges; j++) {
        edge = &self->buffered_edges[j];
        fprintf(out, "\t%f\t%f\t%d\t%d\n", edge->left, edge->right, edge->parent,
            edge->child);
    }
    fprintf(out, "Memory heaps\n");
    for (j = 0; j < self->num_labels; j++) {
        fprintf(out, "segment_heap[%d]:", j);
        object_heap_print_state(&self->segment_heap[j], out);
    }
    fprintf(out, "avl_node_heap:");
    object_heap_print_state(&self->avl_node_heap, out);
    fprintf(out, "node_mapping_heap:");
    object_heap_print_state(&self->node_mapping_heap, out);
    fflush(out);
    msp_verify(self, 0);
out:
    if (ancestors != NULL) {
        free(ancestors);
    }
    return ret;
}

static int MSP_WARN_UNUSED
msp_record_migration(msp_t *self, double left, double right, tsk_id_t node,
    population_id_t source_pop, population_id_t dest_pop)
{
    int ret = 0;

    ret = tsk_migration_table_add_row(&self->tables->migrations, left, right, node,
        source_pop, dest_pop, self->time, NULL, 0);
    if (ret < 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_flush_edges(msp_t *self)
{
    int ret = 0;
    tsk_size_t j, num_edges;
    tsk_edge_t edge;

    if (self->num_buffered_edges > 0) {
        ret = tsk_squash_edges(
            self->buffered_edges, self->num_buffered_edges, &num_edges);
        if (ret != 0) {
            ret = msp_set_tsk_error(ret);
            goto out;
        }
        for (j = 0; j < num_edges; j++) {
            edge = self->buffered_edges[j];
            ret = tsk_edge_table_add_row(&self->tables->edges, edge.left, edge.right,
                edge.parent, edge.child, NULL, 0);
            if (ret < 0) {
                ret = msp_set_tsk_error(ret);
                goto out;
            }
        }
        self->num_buffered_edges = 0;
    }
    ret = 0;
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_store_node(msp_t *self, uint32_t flags, double time, population_id_t population_id,
    tsk_id_t individual)
{
    int ret = 0;

    ret = msp_flush_edges(self);
    if (ret != 0) {
        goto out;
    }
    ret = tsk_node_table_add_row(
        &self->tables->nodes, flags, time, population_id, individual, NULL, 0);
    if (ret < 0) {
        goto out;
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_store_edge(msp_t *self, double left, double right, tsk_id_t parent, tsk_id_t child)
{
    int ret = 0;
    tsk_edge_t *edge;
    const double *node_time = self->tables->nodes.time;

    tsk_bug_assert(parent < (tsk_id_t) self->tables->nodes.num_rows);
    if (self->num_buffered_edges == self->max_buffered_edges - 1) {
        /* Grow the array */
        self->max_buffered_edges *= 2;
        edge = realloc(
            self->buffered_edges, self->max_buffered_edges * sizeof(tsk_edge_t));
        if (edge == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        self->buffered_edges = edge;
    }
    if (node_time[child] >= node_time[parent]) {
        ret = MSP_ERR_TIME_TRAVEL;
        goto out;
    }
    edge = self->buffered_edges + self->num_buffered_edges;
    edge->left = left;
    edge->right = right;
    edge->parent = parent;
    edge->child = child;
    edge->metadata = NULL;
    edge->metadata_length = 0;
    self->num_buffered_edges++;
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_store_arg_edges(msp_t *self, segment_t *z)
{
    int ret = 0;
    tsk_id_t u = (tsk_id_t) msp_get_num_nodes(self) - 1;
    segment_t *x;

    /* Store edges to the left */
    x = z;
    while (x != NULL) {
        if (x->value != u) {
            ret = msp_store_edge(self, x->left, x->right, u, x->value);
            if (ret != 0) {
                goto out;
            }
            x->value = u;
        }
        x = x->prev;
    }

    /* Store edges to the right */
    x = z;
    while (x != NULL) {
        if (x->value != u) {
            ret = msp_store_edge(self, x->left, x->right, u, x->value);
            if (ret != 0) {
                goto out;
            }
            x->value = u;
        }
        x = x->next;
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_move_individual(msp_t *self, avl_node_t *node, avl_tree_t *source,
    population_id_t dest_pop, label_id_t dest_label)
{
    int ret = 0;
    segment_t *ind, *x, *y, *new_ind;
    double recomb_mass, gc_mass;

    if (self->populations[dest_pop].state != MSP_POP_STATE_ACTIVE) {
        ret = MSP_ERR_POPULATION_INACTIVE_MOVE;
        goto out;
    }

    ind = (segment_t *) node->item;
    avl_unlink_node(source, node);
    msp_free_avl_node(self, node);

    if (self->store_full_arg) {
        ret = msp_store_node(
            self, MSP_NODE_IS_MIG_EVENT, self->time, dest_pop, TSK_NULL);
        if (ret < 0) {
            goto out;
        }
        ret = msp_store_arg_edges(self, ind);
        if (ret != 0) {
            goto out;
        }
    }
    if (ind->label == dest_label) {
        /* Need to set the population and label for each segment. */
        new_ind = ind;
        for (x = ind; x != NULL; x = x->next) {
            if (self->store_migrations) {
                ret = msp_record_migration(
                    self, x->left, x->right, x->value, x->population, dest_pop);
                if (ret != 0) {
                    goto out;
                }
            }
            x->population = dest_pop;
        }
    } else {
        /* Because we are changing to a different Fenwick tree we must allocate
         * new segments each time. */
        new_ind = NULL;
        y = NULL;
        for (x = ind; x != NULL; x = x->next) {
            y = msp_alloc_segment(
                self, x->left, x->right, x->value, x->population, dest_label, y, NULL);
            if (new_ind == NULL) {
                new_ind = y;
            } else {
                y->prev->next = y;
            }
            if (self->recomb_mass_index != NULL) {
                recomb_mass
                    = fenwick_get_value(&self->recomb_mass_index[x->label], x->id);
                fenwick_set_value(
                    &self->recomb_mass_index[y->label], y->id, recomb_mass);
            }
            if (self->gc_mass_index != NULL) {
                gc_mass = fenwick_get_value(&self->gc_mass_index[x->label], x->id);
                fenwick_set_value(&self->gc_mass_index[y->label], y->id, gc_mass);
            }
            msp_free_segment(self, x);
        }
    }
    ret = msp_insert_individual(self, new_ind);
out:
    return ret;
}

/*
 * Inserts a population ID into the set of non-empty populations.
 */
static int MSP_WARN_UNUSED
msp_insert_non_empty_population(msp_t *self, tsk_id_t population)
{
    int ret = 0;
    void *value = (void *) (intptr_t) population;
    avl_node_t *node = msp_alloc_avl_node(self);

    if (node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    avl_init_node(node, value);
    if (avl_insert_node(&self->non_empty_populations, node) == NULL) {
        msp_free_avl_node(self, node);
    }
out:
    return ret;
}

/*
 * Removes a population from the non_empty_population set
 */
static int MSP_WARN_UNUSED
msp_remove_non_empty_population(msp_t *self, tsk_id_t population)
{
    int ret = 0;
    avl_node_t *node;
    void *value = (void *) (intptr_t) population;

    node = avl_search(&self->non_empty_populations, value);
    tsk_bug_assert(node != NULL);
    avl_unlink_node(&self->non_empty_populations, node);
    msp_free_avl_node(self, node);
    return ret;
}

/*
 * Inserts a new overlap_count at the specified locus left, mapping to the
 * specified number of overlapping segments b.
 */
static int MSP_WARN_UNUSED
msp_insert_overlap_count(msp_t *self, double left, uint32_t count)
{
    int ret = 0;
    avl_node_t *node = msp_alloc_avl_node(self);
    node_mapping_t *m = msp_alloc_node_mapping(self);

    if (node == NULL || m == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    m->position = left;
    m->value = count;
    avl_init_node(node, m);
    node = avl_insert_node(&self->overlap_counts, node);
    tsk_bug_assert(node != NULL);
out:
    return ret;
}

/*
 * Inserts a new overlap_count at the specified locus, and copies its
 * node mapping from the containing overlap_count.
 */
static int MSP_WARN_UNUSED
msp_copy_overlap_count(msp_t *self, double k)
{
    int ret;
    node_mapping_t search, *nm;
    avl_node_t *node;

    search.position = k;
    avl_search_closest(&self->overlap_counts, &search, &node);
    tsk_bug_assert(node != NULL);
    nm = (node_mapping_t *) node->item;
    if (nm->position > k) {
        node = node->prev;
        tsk_bug_assert(node != NULL);
        nm = (node_mapping_t *) node->item;
    }
    ret = msp_insert_overlap_count(self, k, nm->value);
    return ret;
}

static int
msp_compress_overlap_counts(msp_t *self, double l, double r)
{
    int ret = 0;
    avl_node_t *node1, *node2;
    node_mapping_t search, *nm1, *nm2;

    search.position = l;
    node1 = avl_search(&self->overlap_counts, &search);
    tsk_bug_assert(node1 != NULL);
    if (node1->prev != NULL) {
        node1 = node1->prev;
    }
    node2 = node1->next;
    do {
        nm1 = (node_mapping_t *) node1->item;
        nm2 = (node_mapping_t *) node2->item;
        if (nm1->value == nm2->value) {
            avl_unlink_node(&self->overlap_counts, node2);
            msp_free_avl_node(self, node2);
            msp_free_node_mapping(self, nm2);
            node2 = node1->next;
        } else {
            node1 = node2;
            node2 = node2->next;
        }
    } while (node2 != NULL && nm2->position <= r);
    return ret;
}

static int MSP_WARN_UNUSED
msp_conditional_compress_overlap_counts(msp_t *self, double l, double r)
{
    int ret = 0;
    double covered_fraction = (r - l) / self->sequence_length;

    /* This is a heuristic to prevent us spending a lot of time pointlessly
     * trying to defragment during the early stages of the simulation.
     * 5% of the overall length seems like a good value and leads to
     * a ~15% time reduction when doing large simulations.
     */
    if (covered_fraction < 0.05) {
        ret = msp_compress_overlap_counts(self, l, r);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

/* Defragment the segment chain ending in z by squashing any redundant
 * segments together */
static int MSP_WARN_UNUSED
msp_defrag_segment_chain(msp_t *self, segment_t *z)
{
    segment_t *y, *x;

    y = z;
    while (y->prev != NULL) {
        x = y->prev;
        if (x->right == y->left && x->value == y->value) {
            x->right = y->right;
            x->next = y->next;
            if (y->next != NULL) {
                y->next->prev = x;
            }
            /* msp_add_segment_mass(self, x, y); */
            msp_set_segment_mass(self, x);
            msp_free_segment(self, y);
        }
        y = x;
    }
    return 0;
}

static double
msp_dtwf_generate_breakpoint(msp_t *self, double start)
{
    double left_bound, mass_to_next_recomb, breakpoint;

    left_bound = self->discrete_genome ? start + 1 : start;
    do {
        mass_to_next_recomb = gsl_ran_exponential(self->rng, 1.0);
    } while (mass_to_next_recomb == 0.0);

    breakpoint
        = rate_map_shift_by_mass(&self->recomb_map, left_bound, mass_to_next_recomb);
    return self->discrete_genome ? floor(breakpoint) : breakpoint;
}

#define PEDIGREE_UNINITIALISED (-1)

static int MSP_WARN_UNUSED
msp_reset_pedigree(msp_t *self)
{
    int ret = 0;
    tsk_size_t j, k;
    pedigree_t *pedigree = &self->pedigree;
    individual_t *ind;
    avl_node_t *a;

    for (j = 0; j < pedigree->num_individuals; j++) {
        ind = &pedigree->individuals[j];
        for (k = 0; k < self->ploidy; k++) {
            for (a = ind->common_ancestors[k].head; a != NULL; a = a->next) {
                avl_unlink_node(&ind->common_ancestors[k], a);
                msp_free_avl_node(self, a);
            }
        }
    }
    pedigree->next_individual = PEDIGREE_UNINITIALISED;
    return ret;
}

static int MSP_WARN_UNUSED
msp_pedigree_add_individual_common_ancestor(
    msp_t *self, tsk_id_t individual_id, segment_t *ancestor, tsk_size_t ploid)
{
    int ret;
    individual_t *ind = &self->pedigree.individuals[individual_id];
    avl_node_t *node = msp_alloc_avl_node(self);

    if (node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    tsk_bug_assert(ind->common_ancestors != NULL);
    tsk_bug_assert(ploid < self->ploidy);

    avl_init_node(node, ancestor);
    node = avl_insert_node(&ind->common_ancestors[ploid], node);
    tsk_bug_assert(node != NULL);

    ret = 0;
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_pedigree_add_sample_ancestry(msp_t *self, segment_t *segment)
{
    int ret = 0;
    tsk_size_t ploid;
    tsk_id_t node_id = segment->value;
    tsk_id_t individual_id;
    individual_t *ind;

    tsk_bug_assert(node_id < (tsk_id_t) self->tables->nodes.num_rows);
    individual_id = self->tables->nodes.individual[node_id];
    tsk_bug_assert(individual_id != TSK_NULL);

    ind = &self->pedigree.individuals[individual_id];

    for (ploid = 0; ploid < ind->ploidy; ploid++) {
        if (node_id == ind->nodes[ploid]) {
            break;
        }
    }
    tsk_bug_assert(ploid < ind->ploidy);
    if (avl_count(&ind->common_ancestors[ploid]) > 0) {
        /* This is where we'd deal with the internal samples. What we'll probably
         * need to do is to go through the ancestry after we've processed the
         * individual and then pad out any gaps in the ancestry segments. */
        ret = MSP_ERR_PEDIGREE_INTERNAL_SAMPLE;
        goto out;
    }
    ret = msp_pedigree_add_individual_common_ancestor(self, ind->id, segment, ploid);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_pedigree_initialise(msp_t *self)
{
    int ret = 0;
    population_t *pop;
    segment_t *segment;
    avl_node_t *a;
    label_id_t label = 0;
    tsk_size_t j;

    if (self->next_demographic_event != NULL || self->store_full_arg) {
        ret = MSP_ERR_UNSUPPORTED_OPERATION;
        goto out;
    }
    if (self->ploidy != 2) {
        ret = MSP_ERR_DTWF_DIPLOID_ONLY;
        goto out;
    }
    if (rate_map_get_total_mass(&self->gc_map) != 0.0) {
        /* Could be, we just haven't implemented it */
        ret = MSP_ERR_DTWF_GC_NOT_SUPPORTED;
        goto out;
    }
    tsk_bug_assert(self->num_labels == 1);
    tsk_bug_assert(self->ploidy == 2);

    /* We're not using the mass indexes, so no point in maintaining them */
    tsk_bug_assert(self->recomb_mass_index == NULL);
    tsk_bug_assert(self->gc_mass_index == NULL);

    for (j = 0; j < self->num_populations; j++) {
        pop = &self->populations[j];
        for (a = pop->ancestors[label].head; a != NULL; a = a->next) {
            segment = (segment_t *) a->item;
            ret = msp_pedigree_add_sample_ancestry(self, segment);
            if (ret != 0) {
                goto out;
            }
        }
    }
    self->pedigree.next_individual = 0;
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_dtwf_recombine(msp_t *self, segment_t *x, segment_t **u, segment_t **v)
{
    int ret = 0;
    int ix;
    double k;
    segment_t *y, *z, *tail;
    segment_t s1, s2;
    segment_t *seg_tails[] = { &s1, &s2 };

    k = msp_dtwf_generate_breakpoint(self, x->left);
    s1.next = NULL;
    s2.next = NULL;
    ix = (int) gsl_rng_uniform_int(self->rng, 2);
    seg_tails[ix]->next = x;
    tsk_bug_assert(x->prev == NULL);

    while (x != NULL) {
        seg_tails[ix] = x;
        y = x->next;

        if (x->right > k) {
            // Make new segment
            tsk_bug_assert(x->left < k);
            self->num_re_events++;
            ix = (ix + 1) % 2;

            if (seg_tails[ix] == &s1 || seg_tails[ix] == &s2) {
                tail = NULL;
            } else {
                tail = seg_tails[ix];
            }
            z = msp_alloc_segment(
                self, k, x->right, x->value, x->population, x->label, tail, x->next);
            if (z == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            msp_set_segment_mass(self, z);
            tsk_bug_assert(z->left < z->right);
            if (x->next != NULL) {
                x->next->prev = z;
            }
            seg_tails[ix]->next = z;
            seg_tails[ix] = z;
            x->next = NULL;
            x->right = k;
            msp_set_segment_mass(self, x);
            tsk_bug_assert(x->left < x->right);
            x = z;
            k = msp_dtwf_generate_breakpoint(self, k);
        } else if (x->right <= k && y != NULL && y->left >= k) {
            // Recombine in gap between segment and the next
            x->next = NULL;
            y->prev = NULL;
            while (y->left >= k) {
                self->num_re_events++;
                ix = (ix + 1) % 2;
                k = msp_dtwf_generate_breakpoint(self, k);
            }
            seg_tails[ix]->next = y;
            if (seg_tails[ix] == &s1 || seg_tails[ix] == &s2) {
                tail = NULL;
            } else {
                tail = seg_tails[ix];
            }
            y->prev = tail;
            msp_set_segment_mass(self, y);
            seg_tails[ix] = y;
            x = y;
        } else {
            // Breakpoint in later segment
            x = y;
        }
    }
    // Remove sentinal segments
    *u = s1.next;
    *v = s2.next;
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_store_arg_recombination(msp_t *self, segment_t *lhs_tail, segment_t *rhs)
{
    int ret = 0;

    /* Store the edges for the LHS */
    ret = msp_store_node(
        self, MSP_NODE_IS_RE_EVENT, self->time, lhs_tail->population, TSK_NULL);
    if (ret < 0) {
        goto out;
    }
    ret = msp_store_arg_edges(self, lhs_tail);
    if (ret != 0) {
        goto out;
    }
    /* Store the edges for the RHS */
    ret = msp_store_node(
        self, MSP_NODE_IS_RE_EVENT, self->time, rhs->population, TSK_NULL);
    if (ret < 0) {
        goto out;
    }
    ret = msp_store_arg_edges(self, rhs);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_store_arg_gene_conversion(
    msp_t *self, segment_t *tail, segment_t *alpha, segment_t *head)
{
    int ret = 0;

    /* If head and tail both null then all material moves in the gene conversion so no
     * split to record */
    if (tail != NULL || head != NULL) {
        tsk_bug_assert(alpha != NULL);
        /* Store the edges for tail & head */
        ret = msp_store_node(
            self, MSP_NODE_IS_GC_EVENT, self->time, alpha->population, TSK_NULL);
        if (ret < 0) {
            goto out;
        }
        ret = msp_store_arg_edges(self, tail);
        if (ret != 0) {
            goto out;
        }
        ret = msp_store_arg_edges(self, head);
        if (ret != 0) {
            goto out;
        }
        /* Store the edges for the alpha section */
        ret = msp_store_node(
            self, MSP_NODE_IS_GC_EVENT, self->time, alpha->population, TSK_NULL);
        if (ret < 0) {
            goto out;
        }
        ret = msp_store_arg_edges(self, alpha);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_choose_uniform_breakpoint(msp_t *self, int label, rate_map_t *rate_map,
    fenwick_t *mass_index_array, bool left_at_zero, double *ret_breakpoint,
    segment_t **ret_seg)
{

    int ret = 0;
    double breakpoint, breakpoint_mass, random_mass, y_cumulative_mass, y_right_mass,
        left_bound;
    segment_t *x, *y;
    fenwick_t *tree = &mass_index_array[label];
    int num_breakpoint_resamplings = 0;
    size_t segment_id;

    do {
        /* Choose a recombination mass uniformly from the total and find the
         * segment y that is associated with this *cumulative* value. */
        random_mass = gsl_ran_flat(self->rng, 0, fenwick_get_total(tree));
        segment_id = fenwick_find(tree, random_mass);
        y = msp_get_segment(self, segment_id, label);
        tsk_bug_assert(fenwick_get_value(tree, y->id) > 0);
        x = y->prev;
        y_cumulative_mass = fenwick_get_cumulative_sum(tree, y->id);
        y_right_mass = rate_map_position_to_mass(rate_map, y->right);
        breakpoint_mass = y_right_mass - (y_cumulative_mass - random_mass);
        breakpoint = rate_map_mass_to_position(rate_map, breakpoint_mass);
        if (self->discrete_genome) {
            breakpoint = floor(breakpoint);
        }

        /* Deal with various quirks that can happen with numerical
         * imprecision from going back and forth between rate mass
         * and physical positions. We try to make this robust by making
         * resampling the default case and only break out of the loop
         * when the conditions we need are explicitly met. */
        if (x == NULL) {
            left_bound = left_at_zero ? 0 : y->left;
            /* if there is no previous segment we cannot have breakpoint
             * <= y->left (or zero, if the left limit is zero) */
            if (left_bound < breakpoint && breakpoint < y->right) {
                break;
            }
        } else {
            tsk_bug_assert(x->right <= y->left);
            if (x->right <= breakpoint && breakpoint < y->right) {
                break;
            }
        }
        num_breakpoint_resamplings++;
        /* Arbitrary limit - if we hit this many resamplings then there's
         * definitely something wrong. */
        if (num_breakpoint_resamplings == 10) {
            ret = MSP_ERR_BREAKPOINT_RESAMPLE_OVERFLOW;
            goto out;
        }
    } while (true);

    *ret_breakpoint = breakpoint;
    *ret_seg = y;
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_recombination_event(msp_t *self, label_id_t label, segment_t **lhs, segment_t **rhs)
{
    int ret = 0;
    double breakpoint;
    segment_t *x, *y, *alpha, *lhs_tail;

    self->num_re_events++;
    tsk_bug_assert(self->recomb_mass_index != NULL);

    ret = msp_choose_uniform_breakpoint(
        self, label, &self->recomb_map, self->recomb_mass_index, false, &breakpoint, &y);
    if (ret != 0) {
        goto out;
    }
    x = y->prev;

    if (y->left < breakpoint) {
        tsk_bug_assert(breakpoint < y->right);
        alpha = msp_alloc_segment(self, breakpoint, y->right, y->value, y->population,
            y->label, NULL, y->next);
        if (alpha == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        if (y->next != NULL) {
            y->next->prev = alpha;
        }
        y->next = NULL;
        y->right = breakpoint;
        msp_set_segment_mass(self, y);
        if (msp_has_breakpoint(self, breakpoint)) {
            self->num_multiple_re_events++;
        } else {
            ret = msp_insert_breakpoint(self, breakpoint);
            if (ret != 0) {
                goto out;
            }
        }
        lhs_tail = y;
        tsk_bug_assert(y->left < y->right);
    } else {
        tsk_bug_assert(x != NULL);
        x->next = NULL;
        y->prev = NULL;
        alpha = y;
        self->num_trapped_re_events++;
        lhs_tail = x;
    }
    tsk_bug_assert(alpha->left < alpha->right);
    msp_set_segment_mass(self, alpha);
    ret = msp_insert_individual(self, alpha);
    if (ret != 0) {
        goto out;
    }
    if (self->store_full_arg) {
        ret = msp_store_arg_recombination(self, lhs_tail, alpha);
        if (ret != 0) {
            goto out;
        }
    }
    if (lhs != NULL) {
        x = lhs_tail;
        /* Seek back to the head of the x chain */
        while (x->prev != NULL) {
            x = x->prev;
        }
        *lhs = x;
        *rhs = alpha;
    }
out:
    return ret;
}

static double
msp_generate_gc_tract_length(msp_t *self)
{
    double tl;
    int num_resamplings = 0;

    /* generate tract length */
    if (self->discrete_genome) {
        tl = gsl_ran_geometric(self->rng, 1 / self->gc_tract_length);
    } else {
        do {
            if (num_resamplings == 10) {
                tl = -1;
                goto out;
            }
            tl = gsl_ran_exponential(self->rng, self->gc_tract_length);
            num_resamplings++;
        } while (tl <= 0);
    }
out:
    return tl;
}

static int MSP_WARN_UNUSED
msp_gene_conversion_event(msp_t *self, label_id_t label)
{
    int ret = 0;
    segment_t *x, *y, *alpha, *head, *tail, *z, *new_individual_head;
    double left_breakpoint, right_breakpoint, tl;
    bool insert_alpha;

    tsk_bug_assert(self->gc_mass_index != NULL);
    self->num_gc_events++;
    self->num_internal_gc_events++;
    ret = msp_choose_uniform_breakpoint(
        self, label, &self->gc_map, self->gc_mass_index, true, &left_breakpoint, &y);
    if (ret != 0) {
        goto out;
    }

    x = y->prev;

    /* generate tract length */
    tl = msp_generate_gc_tract_length(self);
    if (tl == -1) {
        ret = MSP_ERR_TRACTLEN_RESAMPLE_OVERFLOW;
        goto out;
    }
    tsk_bug_assert(tl > 0);
    self->sum_internal_gc_tract_lengths += tl;
    right_breakpoint = left_breakpoint + tl;

    if (y->left >= right_breakpoint) {
        //                  y
        // ...  |   |   ========== ...
        //     lbp rbp
        self->num_noneffective_gc_events++;
        return 0;
    }

    /* Process left break */
    insert_alpha = true;
    if (left_breakpoint <= y->left) {
        //  x             y
        // =====  |  ==========
        //       lbp
        //
        // becomes
        //  x
        // =====         α
        //           ==========
        if (x == NULL) {
            // In this case we *don't* insert alpha because it is already
            // the head of a segment chain
            insert_alpha = false;
        } else {
            x->next = NULL;
        }
        y->prev = NULL;
        alpha = y;
        tail = x;
    } else {
        //  x             y
        // =====     ====|=====
        //              lbp
        //
        // becomes
        //  x         y
        // =====     ====   α
        //               ======
        /* alpha = self->copy_segment(y) */
        alpha = msp_copy_segment(self, y);
        if (alpha == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        alpha->left = left_breakpoint;
        alpha->prev = NULL;
        if (y->next != NULL) {
            y->next->prev = alpha;
        }
        y->next = NULL;
        y->right = left_breakpoint;
        msp_set_segment_mass(self, y);
        tail = y;

        if (!msp_has_breakpoint(self, left_breakpoint)) {
            ret = msp_insert_breakpoint(self, left_breakpoint);
            if (ret != 0) {
                goto out;
            }
        }
    }
    msp_set_segment_mass(self, alpha);

    // Find the segment z that the right breakpoint falls in
    z = alpha;
    while (z != NULL && right_breakpoint >= z->right) {
        z = z->next;
    }

    head = NULL;
    // Process the right break
    if (z != NULL) {
        if (z->left < right_breakpoint) {
            //   tail             z
            // ======
            //       ...  ===|==========
            //              rbp
            //
            // becomes
            //  tail              head
            // =====         ===========
            //      ...   ===
            //             z
            head = msp_copy_segment(self, z);
            if (head == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            head->left = right_breakpoint;
            if (z->next != NULL) {
                z->next->prev = head;
            }
            z->right = right_breakpoint;
            z->next = NULL;
            msp_set_segment_mass(self, z);

            if (!msp_has_breakpoint(self, right_breakpoint)) {
                ret = msp_insert_breakpoint(self, right_breakpoint);
                if (ret != 0) {
                    goto out;
                }
            }
        } else {
            //   tail             z
            // ======
            //   ...   |   =============
            //        rbp
            //
            // becomes
            //  tail             z
            // ======      =============
            //  ...
            if (z->prev != NULL) {
                z->prev->next = NULL;
            }
            head = z;
        }
        if (tail != NULL) {
            tail->next = head;
        }
        head->prev = tail;
        msp_set_segment_mass(self, head);
    }

    //        y            z
    //  |  ========== ... ===== |
    // lbp                     rbp
    // When y and z are the head and tail of the segment chains, then
    // this GC event does nothing. This logic takes care of this situation.
    new_individual_head = NULL;
    if (insert_alpha) {
        new_individual_head = alpha;
    } else if (head != NULL) {
        new_individual_head = head;
    }
    if (new_individual_head != NULL) {
        ret = msp_insert_individual(self, new_individual_head);
    } else {
        self->num_noneffective_gc_events++;
    }
    if (self->store_full_arg) {
        ret = msp_store_arg_gene_conversion(self, tail, alpha, head);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

/* If we're implementing the SMC or SMC', discard this CA event if
 * there aren't any overlapping segments.
 */
static int MSP_WARN_UNUSED
msp_reject_ca_event(msp_t *self, segment_t *a, segment_t *b)
{
    int ret = 0;
    segment_t *x = a;
    segment_t *y = b;
    segment_t *beta;
    double overlap;
    int model = self->model.type;

    if (model == MSP_MODEL_SMC || model == MSP_MODEL_SMC_PRIME) {
        ret = 1;
        while (x != NULL && y != NULL) {
            if (y->left < x->left) {
                beta = x;
                x = y;
                y = beta;
            }
            overlap = x->right - y->left;
            /* For the SMC' overlap must be >= 0, but SMC it must
             * be strictly greater than zero, as it doesn't allow
             * directly adjacent segments to coalesce */
            if ((model == MSP_MODEL_SMC_PRIME && overlap >= 0)
                || (overlap > 0)) { /* SMC */
                ret = 0;
                break;
            }
            x = x->next;
        }
    }
    return ret;
}

static int MSP_WARN_UNUSED
msp_merge_two_ancestors(msp_t *self, population_id_t population_id, label_id_t label,
    segment_t *a, segment_t *b, tsk_id_t new_node_id, segment_t **ret_merged_head)
{
    int ret = 0;
    bool coalescence = false;
    bool defrag_required = false;
    tsk_id_t v;
    double l, r, l_min, r_max;
    avl_node_t *node;
    node_mapping_t *nm, search;
    segment_t *x, *y, *z, *alpha, *beta, *merged_head;

    x = a;
    y = b;
    merged_head = NULL;
    /* Keep GCC happy */
    l_min = 0;
    r_max = 0;

    /* update recomb mass and get ready for loop */
    z = NULL;
    while (x != NULL || y != NULL) {
        alpha = NULL;
        if (x == NULL || y == NULL) {
            if (x != NULL) {
                alpha = x;
                x = NULL;
            }
            if (y != NULL) {
                alpha = y;
                y = NULL;
            }
        } else {
            if (y->left < x->left) {
                beta = x;
                x = y;
                y = beta;
            }
            if (x->right <= y->left) {
                alpha = x;
                x = x->next;
                alpha->next = NULL;
            } else if (x->left != y->left) {
                alpha = msp_alloc_segment(self, x->left, y->left, x->value,
                    x->population, x->label, NULL, NULL);
                if (alpha == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
                x->left = y->left;
            } else {
                l = x->left;
                r_max = GSL_MIN(x->right, y->right);
                if (!coalescence) {
                    coalescence = true;
                    l_min = l;
                    if (new_node_id == TSK_NULL) {
                        new_node_id = msp_store_node(
                            self, 0, self->time, population_id, TSK_NULL);
                        if (new_node_id < 0) {
                            ret = (int) new_node_id;
                            goto out;
                        }
                    }
                }
                v = new_node_id;
                /* Insert overlap counts for bounds, if necessary */
                search.position = l;
                node = avl_search(&self->overlap_counts, &search);
                if (node == NULL) {
                    ret = msp_copy_overlap_count(self, l);
                    if (ret < 0) {
                        goto out;
                    }
                }
                search.position = r_max;
                node = avl_search(&self->overlap_counts, &search);
                if (node == NULL) {
                    ret = msp_copy_overlap_count(self, r_max);
                    if (ret < 0) {
                        goto out;
                    }
                }
                /* Now get overlap count at the left */
                search.position = l;
                node = avl_search(&self->overlap_counts, &search);
                tsk_bug_assert(node != NULL);
                nm = (node_mapping_t *) node->item;
                if (nm->value == 2) {
                    nm->value = 0;
                    node = node->next;
                    tsk_bug_assert(node != NULL);
                    nm = (node_mapping_t *) node->item;
                    r = nm->position;
                } else {
                    r = l;
                    while (nm->value != 2 && r < r_max) {
                        nm->value--;
                        node = node->next;
                        tsk_bug_assert(node != NULL);
                        nm = (node_mapping_t *) node->item;
                        r = nm->position;
                    }
                    alpha = msp_alloc_segment(
                        self, l, r, v, population_id, label, NULL, NULL);
                    if (alpha == NULL) {
                        ret = MSP_ERR_NO_MEMORY;
                        goto out;
                    }
                }
                tsk_bug_assert(v != x->value);
                ret = msp_store_edge(self, l, r, v, x->value);
                if (ret != 0) {
                    goto out;
                }
                ret = msp_store_edge(self, l, r, v, y->value);
                if (ret != 0) {
                    goto out;
                }
                /* Trim the ends of x and y, and prepare for next iteration. */
                if (x->right == r) {
                    beta = x;
                    x = x->next;
                    msp_free_segment(self, beta);
                } else {
                    x->left = r;
                }
                if (y->right == r) {
                    beta = y;
                    y = y->next;
                    msp_free_segment(self, beta);
                } else {
                    y->left = r;
                }
            }
        }
        if (alpha != NULL) {
            if (z == NULL) {
                ret = msp_insert_individual(self, alpha);
                if (ret != 0) {
                    goto out;
                }
                merged_head = alpha;
            } else {
                if (self->store_full_arg) {
                    // we pre-empt the fact that values will be set equal later
                    defrag_required |= z->right == alpha->left;
                } else {
                    defrag_required
                        |= z->right == alpha->left && z->value == alpha->value;
                }
                tsk_bug_assert(z->right <= alpha->left);
                z->next = alpha;
            }
            alpha->prev = z;
            msp_set_segment_mass(self, alpha);
            z = alpha;
        }
    }
    if (self->store_full_arg) {
        if (!coalescence) {
            ret = msp_store_node(
                self, MSP_NODE_IS_CA_EVENT, self->time, population_id, TSK_NULL);
            if (ret < 0) {
                goto out;
            }
        }
        ret = msp_store_arg_edges(self, z);
        if (ret != 0) {
            goto out;
        }
    }
    if (defrag_required) {
        ret = msp_defrag_segment_chain(self, z);
        if (ret != 0) {
            goto out;
        }
    }
    if (coalescence) {
        ret = msp_conditional_compress_overlap_counts(self, l_min, r_max);
        if (ret != 0) {
            goto out;
        }
    }
    if (ret_merged_head != NULL) {
        *ret_merged_head = merged_head;
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_priority_queue_insert(msp_t *self, avl_tree_t *Q, segment_t *u)
{
    int ret = 0;
    avl_node_t *node;

    tsk_bug_assert(u != NULL);
    node = msp_alloc_avl_node(self);
    if (node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    avl_init_node(node, u);
    node = avl_insert_node(Q, node);
    tsk_bug_assert(node != NULL);
out:
    return ret;
}

static segment_t *
msp_priority_queue_pop(msp_t *self, avl_tree_t *Q)
{
    avl_node_t *node = Q->head;
    segment_t *seg = (segment_t *) node->item;
    msp_free_avl_node(self, node);
    avl_unlink_node(Q, node);

    return seg;
}

/* Merge the specified set of ancestors into a single ancestor. This is a
 * generalisation of the msp_common_ancestor_event method where we allow
 * any number of ancestors to merge. The AVL tree is a priority queue in
 * sorted by left coordinate.
 */
static int MSP_WARN_UNUSED
msp_merge_ancestors(msp_t *self, avl_tree_t *Q, population_id_t population_id,
    label_id_t label, tsk_id_t new_node_id, segment_t **ret_merged_head)

{
    int ret = MSP_ERR_GENERIC;
    bool coalescence = false;
    bool defrag_required = false;
    uint32_t j, h;
    double l, r, r_max, next_l, l_min;
    avl_node_t *node;
    node_mapping_t *nm, search;
    segment_t *x, *z, *alpha;
    segment_t **H = NULL;
    segment_t *merged_head = NULL;
    tsk_id_t individual = TSK_NULL;

    H = malloc(avl_count(Q) * sizeof(segment_t *));
    if (H == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    r_max = 0; /* keep compiler happy */
    l_min = 0;
    z = NULL;
    merged_head = NULL;
    while (avl_count(Q) > 0) {
        h = 0;
        node = Q->head;
        l = ((segment_t *) node->item)->left;
        r_max = self->sequence_length;
        while (node != NULL && ((segment_t *) node->item)->left == l) {
            H[h] = (segment_t *) node->item;
            r_max = GSL_MIN(r_max, H[h]->right);
            h++;
            msp_free_avl_node(self, node);
            avl_unlink_node(Q, node);
            node = node->next;
        }
        next_l = 0;
        if (node != NULL) {
            next_l = ((segment_t *) node->item)->left;
            r_max = GSL_MIN(r_max, next_l);
        }
        alpha = NULL;
        if (h == 1) {
            x = H[0];
            if (node != NULL && next_l < x->right) {
                alpha = msp_alloc_segment(self, x->left, next_l, x->value, x->population,
                    x->label, NULL, NULL);
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
                ret = msp_priority_queue_insert(self, Q, x);
                if (ret != 0) {
                    goto out;
                }
            }
        } else {
            coalescence = true;
            if (new_node_id == TSK_NULL) {
                l_min = l;
                new_node_id
                    = msp_store_node(self, 0, self->time, population_id, individual);
                if (new_node_id < 0) {
                    ret = (int) new_node_id;
                    goto out;
                }
            }
            /* Insert overlap counts for bounds, if necessary */
            search.position = l;
            node = avl_search(&self->overlap_counts, &search);
            if (node == NULL) {
                ret = msp_copy_overlap_count(self, l);
                if (ret < 0) {
                    goto out;
                }
            }
            search.position = r_max;
            node = avl_search(&self->overlap_counts, &search);
            if (node == NULL) {
                ret = msp_copy_overlap_count(self, r_max);
                if (ret < 0) {
                    goto out;
                }
            }
            /* Update the extant segments and allocate alpha if the interval
             * has not coalesced. */
            search.position = l;
            node = avl_search(&self->overlap_counts, &search);
            tsk_bug_assert(node != NULL);
            nm = (node_mapping_t *) node->item;
            if (nm->value == h) {
                nm->value = 0;
                node = node->next;
                tsk_bug_assert(node != NULL);
                nm = (node_mapping_t *) node->item;
                r = nm->position;
            } else {
                r = l;
                while (nm->value != h && r < r_max) {
                    nm->value -= h - 1;
                    node = node->next;
                    tsk_bug_assert(node != NULL);
                    nm = (node_mapping_t *) node->item;
                    r = nm->position;
                }
                alpha = msp_alloc_segment(
                    self, l, r, new_node_id, population_id, label, NULL, NULL);
                if (alpha == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
            }
            /* Store the edges and update the priority queue */
            for (j = 0; j < h; j++) {
                x = H[j];
                tsk_bug_assert(new_node_id != x->value);
                ret = msp_store_edge(self, l, r, new_node_id, x->value);
                if (ret != 0) {
                    goto out;
                }
                if (x->right == r) {
                    msp_free_segment(self, x);
                    x = x->next;
                } else if (x->right > r) {
                    x->left = r;
                }
                if (x != NULL) {
                    ret = msp_priority_queue_insert(self, Q, x);
                    if (ret != 0) {
                        goto out;
                    }
                }
            }
        }
        /* Loop tail; integrate alpha into the global state */
        if (alpha != NULL) {
            if (z == NULL) {
                merged_head = alpha;
                ret = msp_insert_individual(self, alpha);
                if (ret != 0) {
                    goto out;
                }
            } else {
                if (self->store_full_arg) {
                    // we pre-empt the fact that values will be set equal later
                    defrag_required |= z->right == alpha->left;
                } else {
                    defrag_required
                        |= z->right == alpha->left && z->value == alpha->value;
                }
                z->next = alpha;
            }
            alpha->prev = z;
            msp_set_segment_mass(self, alpha);
            z = alpha;
        }
    }
    if (self->store_full_arg) {
        if (!coalescence) {
            ret = msp_store_node(
                self, MSP_NODE_IS_CA_EVENT, self->time, population_id, individual);
            if (ret < 0) {
                goto out;
            }
        }
        ret = msp_store_arg_edges(self, z);
        if (ret != 0) {
            goto out;
        }
    }
    if (defrag_required) {
        ret = msp_defrag_segment_chain(self, z);
        if (ret != 0) {
            goto out;
        }
    }
    if (coalescence) {
        ret = msp_conditional_compress_overlap_counts(self, l_min, r_max);
        if (ret != 0) {
            goto out;
        }
    }
    if (ret_merged_head != NULL) {
        *ret_merged_head = merged_head;
    }
    ret = 0;
out:
    if (H != NULL) {
        free(H);
    }
    return ret;
}

/* Merge n ancestors as efficiently as we can by using specialised
 * methods where possible and falling back on msp_merge_ancestors
 * where necessary. */
static int MSP_WARN_UNUSED
msp_merge_n_ancestors(msp_t *self, avl_tree_t *Q, population_id_t population_id,
    label_id_t label, tsk_id_t new_node_id, segment_t **ret_merged_head)
{
    int ret = 0;
    const size_t num_common_ancestors = avl_count(Q);
    population_t *current_pop;
    avl_node_t *a, *avl_node;
    segment_t *merged_head = NULL;
    segment_t *u, *v;

    /* Migrate any of the child segments to this population, if necessary */
    for (a = Q->head; a != NULL; a = a->next) {
        u = (segment_t *) a->item;
        if (u->population != population_id) {
            current_pop = &self->populations[u->population];
            avl_node = avl_search(&current_pop->ancestors[label], u);
            tsk_bug_assert(avl_node != NULL);
            ret = msp_move_individual(
                self, avl_node, &current_pop->ancestors[label], population_id, label);
            if (ret != 0) {
                goto out;
            }
        }
    }

    if (num_common_ancestors == 1) {
        merged_head = msp_priority_queue_pop(self, Q);
    } else if (num_common_ancestors >= 2) {
        msp_remove_individuals_from_population(self, Q);
        if (num_common_ancestors == 2) {
            u = msp_priority_queue_pop(self, Q);
            v = msp_priority_queue_pop(self, Q);
            ret = msp_merge_two_ancestors(
                self, population_id, label, u, v, new_node_id, &merged_head);
        } else {
            ret = msp_merge_ancestors(
                self, Q, population_id, label, new_node_id, &merged_head);
        }
    }
    if (ret_merged_head != NULL) {
        *ret_merged_head = merged_head;
    }
    if (merged_head != NULL) {
        tsk_bug_assert(merged_head->population == population_id);
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_migration_event(msp_t *self, population_id_t source_pop, population_id_t dest_pop)
{
    int ret = 0;
    uint32_t j;
    avl_node_t *node;
    label_id_t label = 0; /* For now only support label 0 */
    avl_tree_t *source = &self->populations[source_pop].ancestors[label];
    size_t index = ((size_t) source_pop) * self->num_populations + (size_t) dest_pop;

    self->num_migration_events[index]++;
    j = (uint32_t) gsl_rng_uniform_int(self->rng, avl_count(source));
    node = avl_at(source, j);
    tsk_bug_assert(node != NULL);
    ret = msp_move_individual(self, node, source, dest_pop, label);
    return ret;
}

static int MSP_WARN_UNUSED
msp_reset_memory_state(msp_t *self)
{
    int ret = 0;
    avl_node_t *node;
    node_mapping_t *nm;
    population_t *pop;
    segment_t *u, *v;
    label_id_t label;
    size_t j;

    for (j = 0; j < self->num_populations; j++) {
        pop = &self->populations[j];
        for (label = 0; label < (label_id_t) self->num_labels; label++) {
            for (node = pop->ancestors[label].head; node != NULL; node = node->next) {
                u = (segment_t *) node->item;
                while (u != NULL) {
                    v = u->next;
                    msp_free_segment(self, u);
                    u = v;
                }
                avl_unlink_node(&pop->ancestors[label], node);
                msp_free_avl_node(self, node);
            }
        }
    }
    for (node = self->breakpoints.head; node != NULL; node = node->next) {
        nm = (node_mapping_t *) node->item;
        avl_unlink_node(&self->breakpoints, node);
        msp_free_avl_node(self, node);
        msp_free_node_mapping(self, nm);
    }
    for (node = self->overlap_counts.head; node != NULL; node = node->next) {
        nm = (node_mapping_t *) node->item;
        avl_unlink_node(&self->overlap_counts, node);
        msp_free_avl_node(self, node);
        msp_free_node_mapping(self, nm);
    }
    return ret;
}

static int
msp_insert_root_segments(msp_t *self, const segment_t *head, segment_t **new_head)
{
    int ret = 0;
    segment_t *copy, *prev;
    const segment_t *seg;
    double breakpoints[2];
    int j;

    prev = NULL;
    for (seg = head; seg != NULL; seg = seg->next) {
        /* Insert breakpoints, if we need to */
        breakpoints[0] = seg->left;
        breakpoints[1] = seg->right;
        for (j = 0; j < 2; j++) {
            if (breakpoints[j] != 0 && breakpoints[j] != self->sequence_length
                && !msp_has_breakpoint(self, breakpoints[j])) {
                ret = msp_insert_breakpoint(self, breakpoints[j]);
                if (ret != 0) {
                    goto out;
                }
            }
        }
        /* Copy the segment and insert into the global state */
        copy = msp_copy_segment(self, seg);
        if (copy == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        if (seg == head && new_head != NULL) {
            *new_head = copy;
        }
        copy->prev = prev;
        if (prev == NULL) {
            ret = msp_insert_individual(self, copy);
            if (ret != 0) {
                goto out;
            }

        } else {
            prev->next = copy;
        }
        msp_set_segment_mass(self, copy);
        prev = copy;
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_insert_sample(msp_t *self, tsk_id_t node)
{
    int ret = 0;
    segment_t *root_seg;
    population_t pop;

    root_seg = self->root_segments[node];
    pop = self->populations[root_seg->population];
    if (pop.state != MSP_POP_STATE_ACTIVE) {
        ret = MSP_ERR_POPULATION_INACTIVE_SAMPLE;
        goto out;
    }
    if (pop.initial_size == 0) {
        ret = MSP_ERR_POP_SIZE_ZERO_SAMPLE;
        goto out;
    }
    ret = msp_insert_root_segments(self, root_seg, NULL);
out:
    return ret;
}

static inline int
msp_allocate_root_segments(msp_t *self, tsk_tree_t *tree, double left, double right,
    segment_t *restrict *root_segments_head, segment_t *restrict *root_segments_tail)
{
    int ret = 0;
    tsk_id_t root;
    segment_t *seg, *tail;
    population_id_t population;
    const population_id_t *restrict node_population = self->tables->nodes.population;
    label_id_t label = 0; /* For now only support label 0 */

    for (root = tsk_tree_get_left_root(tree); root != TSK_NULL;
         root = tree->right_sib[root]) {
        population = node_population[root];
        /* tskit will make sure that population references are good, but
         * we can still have NULL refs. */
        if (population == TSK_NULL) {
            ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
            goto out;
        }
        if (root_segments_head[root] == NULL) {
            seg = msp_alloc_segment(
                self, left, right, root, population, label, NULL, NULL);
            if (seg == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            root_segments_head[root] = seg;
            root_segments_tail[root] = seg;
        } else {
            tail = root_segments_tail[root];
            if (tail->right == left) {
                tail->right = right;
            } else {
                seg = msp_alloc_segment(
                    self, left, right, root, population, label, tail, NULL);
                if (seg == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
                tail->next = seg;
                root_segments_tail[root] = seg;
            }
        }
    }
out:
    return ret;
}

static int
msp_process_input_trees(msp_t *self)
{
    int ret = 0;
    int t_iter;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    uint32_t overlap_count, last_overlap_count;
    tsk_size_t num_trees, num_roots;
    const size_t num_nodes = self->tables->nodes.num_rows;
    overlap_count_t *overlap;
    segment_t **root_segments_tail = NULL;

    /* Initialise the memory for the tree and tree sequence so we can
     * safely free them in all cases */
    memset(&ts, 0, sizeof(ts));
    memset(&tree, 0, sizeof(tree));

    ret = tsk_treeseq_init(&ts, self->tables, TSK_TS_INIT_BUILD_INDEXES);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
    num_trees = tsk_treeseq_get_num_trees(&ts);

    root_segments_tail = calloc(num_nodes + 1, sizeof(*root_segments_tail));
    self->root_segments = calloc(num_nodes + 1, sizeof(*self->root_segments));
    /* We can't have more than num_trees intervals, and allow for one sentinel */
    self->initial_overlaps = calloc(num_trees + 1, sizeof(*self->initial_overlaps));

    if (self->root_segments == NULL || root_segments_tail == NULL
        || self->initial_overlaps == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    ret = tsk_tree_init(&tree, &ts, 0);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }

    overlap = self->initial_overlaps;
    last_overlap_count = UINT32_MAX;
    for (t_iter = tsk_tree_first(&tree); t_iter == 1; t_iter = tsk_tree_next(&tree)) {
        num_roots = tsk_tree_get_num_roots(&tree);
        overlap_count = 0;
        if (num_roots > 1) {
            overlap_count = (uint32_t) num_roots;
            ret = msp_allocate_root_segments(self, &tree, tree.interval.left,
                tree.interval.right, self->root_segments, root_segments_tail);
            if (ret != 0) {
                goto out;
            }
        }
        if (overlap_count != last_overlap_count) {
            overlap->left = tree.interval.left;
            overlap->count = overlap_count;
            overlap++;
            last_overlap_count = overlap_count;
        }
    }
    if (t_iter != 0) {
        ret = msp_set_tsk_error(t_iter);
        goto out;
    }
    overlap->left = self->sequence_length;
    overlap->count = UINT32_MAX;

out:
    tsk_treeseq_free(&ts);
    tsk_tree_free(&tree);
    msp_safe_free(root_segments_tail);
    return ret;
}

static int
msp_reset_population_state(msp_t *self)
{
    int ret = 0;
    overlap_count_t *overlap = self->initial_overlaps;
    while (true) {
        ret = msp_insert_overlap_count(self, overlap->left, overlap->count);
        if (ret != 0) {
            goto out;
        }
        if (overlap->left == self->sequence_length) {
            break;
        }
        overlap++;
    }
out:
    return ret;
}

/* Apply all demographic events at the specified time
 */
static int MSP_WARN_UNUSED
msp_apply_demographic_events(msp_t *self, double time)
{
    int ret = 0;
    demographic_event_t *event;

    while (self->next_demographic_event != NULL
           && self->next_demographic_event->time == time) {
        event = self->next_demographic_event;
        self->time = time;
        tsk_bug_assert(event->change_state != NULL);
        ret = event->change_state(self, event);
        if (ret != 0) {
            goto out;
        }
        self->next_demographic_event = event->next;
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_apply_sampling_events(msp_t *self, double time)
{
    int ret = 0;
    sampling_event_t *se;

    /* Add in all samples with this time */
    while (self->next_sampling_event < self->num_sampling_events
           && self->sampling_events[self->next_sampling_event].time == time) {
        se = &self->sampling_events[self->next_sampling_event];
        self->time = time;
        ret = msp_insert_sample(self, se->sample);
        if (ret != 0) {
            goto out;
        }
        self->next_sampling_event++;
        /* To keep things simple, just insert the population
         * unconditionally for each sample. */
        ret = msp_insert_non_empty_population(self, se->population);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_pedigree_insert_ancient_samples(msp_t *self)
{
    int ret = 0;
    sampling_event_t *se;
    segment_t *root_seg, *new_head;

    while (self->next_sampling_event < self->num_sampling_events
           && self->sampling_events[self->next_sampling_event].time <= self->time) {
        se = self->sampling_events + self->next_sampling_event;
        root_seg = self->root_segments[se->sample];
        ret = msp_insert_root_segments(self, root_seg, &new_head);
        if (ret != 0) {
            goto out;
        }
        ret = msp_pedigree_add_sample_ancestry(self, new_head);
        if (ret != 0) {
            goto out;
        }
        self->next_sampling_event++;
    }
out:
    return ret;
}

static double
msp_get_next_fixed_event_time(msp_t *self)
{
    double sampling_event_time = DBL_MAX;
    double demographic_event_time = DBL_MAX;

    if (self->next_sampling_event < self->num_sampling_events) {
        sampling_event_time = self->sampling_events[self->next_sampling_event].time;
    }
    if (self->next_demographic_event != NULL) {
        demographic_event_time = self->next_demographic_event->time;
    }
    return GSL_MIN(sampling_event_time, demographic_event_time);
}

static int MSP_WARN_UNUSED
msp_apply_fixed_events(msp_t *self, double time)
{
    int ret = 0;

    ret = msp_apply_demographic_events(self, time);
    if (ret != 0) {
        goto out;
    }
    ret = msp_apply_sampling_events(self, time);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int
msp_reset(msp_t *self)
{
    int ret = 0;
    size_t N = self->num_populations;
    double event_time;
    population_id_t population_id;
    population_t *pop, *initial_pop;

    memcpy(&self->model, &self->initial_model, sizeof(self->model));
    ret = msp_reset_pedigree(self);
    if (ret != 0) {
        goto out;
    }
    ret = msp_reset_memory_state(self);
    if (ret != 0) {
        goto out;
    }
    /* Set up the initial segments and algorithm state */
    for (population_id = 0; population_id < (population_id_t) N; population_id++) {
        pop = self->populations + population_id;
        /* Set the initial population parameters */
        initial_pop = &self->initial_populations[population_id];
        pop->growth_rate = initial_pop->growth_rate;
        pop->initial_size = initial_pop->initial_size;
        pop->start_time = 0;
        pop->state = initial_pop->state;
    }
    /* Reset the tables to their correct position for replication */
    ret = tsk_table_collection_truncate(self->tables, &self->input_position);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
    tsk_bug_assert(self->tables->populations.num_rows == self->num_populations);

    ret = msp_reset_population_state(self);
    if (ret != 0) {
        goto out;
    }

    self->next_demographic_event = self->demographic_events_head;
    memcpy(
        self->migration_matrix, self->initial_migration_matrix, N * N * sizeof(double));
    self->next_sampling_event = 0;
    self->num_re_events = 0;
    self->num_gc_events = 0;
    self->num_internal_gc_events = 0;
    self->sum_internal_gc_tract_lengths = 0;
    self->num_noneffective_gc_events = 0;
    self->num_ca_events = 0;
    self->num_rejected_ca_events = 0;
    self->num_trapped_re_events = 0;
    self->num_multiple_re_events = 0;
    memset(self->num_migration_events, 0, N * N * sizeof(size_t));

    if (self->start_time < DBL_MAX) {
        while ((event_time = msp_get_next_fixed_event_time(self)) <= self->start_time) {
            ret = msp_apply_fixed_events(self, event_time);
            if (ret != 0) {
                goto out;
            }
        }
    }
    self->time = self->start_time;
    self->state = MSP_STATE_INITIALISED;
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_initialise_simulation_state(msp_t *self)
{
    int ret = 0;
    tsk_size_t num_samples;
    size_t j;
    double min_root_time;
    segment_t *head;
    tsk_id_t root;
    const double *restrict node_time = self->tables->nodes.time;
    const tsk_id_t *restrict node_population = self->tables->nodes.population;
    tsk_id_t *samples = malloc(self->tables->nodes.num_rows * sizeof(*samples));

    if (samples == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    ret = msp_process_input_trees(self);
    if (ret != 0) {
        goto out;
    }

    /* Process the root segments */
    min_root_time = DBL_MAX;
    for (j = 0; j < self->input_position.nodes; j++) {
        head = self->root_segments[j];
        if (head != NULL) {
            root = head->value;
            min_root_time = GSL_MIN(node_time[root], min_root_time);
        }
    }
    if (self->input_position.nodes == 0) {
        /* When we're initialising for debugging we provide no samples. */
        self->start_time = 0;
    } else {
        /* Make sure that the start-time is no less than the time of the
         * youngest root */
        self->start_time = GSL_MAX(self->start_time, min_root_time);
    }

    num_samples = 0;
    for (j = 0; j < self->input_position.nodes; j++) {
        head = self->root_segments[j];
        if (head != NULL) {
            root = head->value;
            samples[num_samples] = root;
            num_samples++;
        }
    }

    /* Set up the sampling events */
    self->num_sampling_events = num_samples;
    self->sampling_events = NULL;
    if (self->num_sampling_events > 0) {
        self->sampling_events
            = malloc(self->num_sampling_events * sizeof(sampling_event_t));
        if (self->sampling_events == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        for (j = 0; j < num_samples; j++) {
            root = samples[j];
            self->sampling_events[j].sample = root;
            self->sampling_events[j].time = node_time[root];
            self->sampling_events[j].population = node_population[root];
        }
        /* Sort the sampling events by time. */
        qsort(self->sampling_events, self->num_sampling_events, sizeof(sampling_event_t),
            cmp_sampling_event);
    }

    /* ret = msp_compress_overlap_counts(self, 0, self->sequence_length); */
out:
    msp_safe_free(samples);
    return ret;
}

/*
 * Sets up the memory heaps.
 */
int MSP_WARN_UNUSED
msp_initialise(msp_t *self)
{
    int ret = -1;

    /* Bookmark the tables so that we know where to reset once for each
     * simulation. */
    ret = tsk_table_collection_record_num_rows(self->tables, &self->input_position);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }

    ret = msp_alloc_memory_blocks(self);
    if (ret != 0) {
        goto out;
    }
    ret = msp_setup_mass_indexes(self);
    if (ret != 0) {
        goto out;
    }
    ret = msp_initialise_simulation_state(self);
    if (ret != 0) {
        goto out;
    }
    /* Copy the state of the simulation model into the initial model */
    memcpy(&self->initial_model, &self->model, sizeof(self->model));
    ret = msp_reset(self);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

/* In the exceedingly rare cases where gsl_ran_exponential returns
 * 0, we return the smallest representable value > the current time
 * to avoid returning a tree sequence with zero length branches.
 * Note that we can still return 0 from this function if population
 * sizes are extremely small. This is intentional, as it is almost
 * certainly an error to have simulations at such extreme values
 * and the user should be alerted to this. */
static double
handle_zero_waiting_time(double t)
{
    double ret = nextafter(t, DBL_MAX) - t;
    tsk_bug_assert(ret != 0);
    return ret;
}

/* Given the specified rate, return the waiting time until the next common ancestor
 * event for the specified population */
static double
msp_get_common_ancestor_waiting_time_from_rate(
    msp_t *self, population_t *pop, double lambda)
{
    double ret = DBL_MAX;
    double alpha = pop->growth_rate;
    double t = self->time;
    double u, dt, z;

    if (lambda > 0.0) {
        u = gsl_ran_exponential(self->rng, 1.0 / lambda);
        if (alpha == 0.0) {
            ret = self->ploidy * pop->initial_size * u;
        } else {
            dt = t - pop->start_time;
            z = 1 + self->ploidy * alpha * pop->initial_size * exp(-alpha * dt) * u;
            /* if z is <= 0 no coancestry can occur */
            if (z > 0) {
                ret = log(z) / alpha;
            }
        }
        if (u == 0) {
            ret = handle_zero_waiting_time(t);
        }
    }
    return ret;
}

/* Computes the set of non empty populations and the set
 * of populations reachable from each population. */
static int MSP_WARN_UNUSED
msp_compute_population_indexes(msp_t *self)
{
    int ret = 0;
    const tsk_id_t N = (tsk_id_t) self->num_populations;
    tsk_id_t j, k;
    double migration_rate;
    population_t *pop;
    avl_node_t *avl_node;

    /* Set up the possible destinations for each population */
    for (j = 0; j < N; j++) {
        pop = &self->populations[j];
        pop->num_potential_destinations = 0;
        for (k = 0; k < N; k++) {
            migration_rate = self->migration_matrix[j * N + k];
            if (migration_rate > 0) {
                pop->potential_destinations[pop->num_potential_destinations] = k;
                pop->num_potential_destinations++;
            }
        }
    }

    /* Set up the non_empty_populations */
    /* First clear out any existing structures */
    for (avl_node = self->non_empty_populations.head; avl_node != NULL;
         avl_node = avl_node->next) {
        avl_unlink_node(&self->non_empty_populations, avl_node);
        msp_free_avl_node(self, avl_node);
    }
    tsk_bug_assert(avl_count(&self->non_empty_populations) == 0);
    for (j = 0; j < N; j++) {
        if (msp_get_num_population_ancestors(self, j) > 0) {
            ret = msp_insert_non_empty_population(self, j);
            if (ret != 0) {
                goto out;
            }
        }
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_sample_waiting_time(
    msp_t *self, fenwick_t *mass_indexes, label_id_t label, double *ret_t_wait)
{
    int ret = 0;
    double total_mass, t_wait, lambda;
    fenwick_t *mass_index;

    if (mass_indexes == NULL) {
        /* When the mass index is NULL this sigifies a total rate of zero
         * and so we have an infinite waiting time. */
        *ret_t_wait = DBL_MAX;
    } else {
        mass_index = &mass_indexes[label];
        /* In very large simulations, the fenwick tree used as an indexing
         * structure for genomic segments will experience some numerical
         * drift, where the indexed values diverge from the true values
         * associated with segments. We ensure that this drift does not
         * become too large by rebuilding the indexing structure every
         * now and again. */

        if (fenwick_rebuild_required(mass_index)) {
            fenwick_rebuild(mass_index);
            self->num_fenwick_rebuilds++;
        }

        total_mass = fenwick_get_total(mass_index);
        if (!isfinite(total_mass)) {
            ret = MSP_ERR_BREAKPOINT_MASS_NON_FINITE;
            goto out;
        }
        lambda = total_mass;
        t_wait = DBL_MAX;
        if (lambda > 0.0) {
            t_wait = gsl_ran_exponential(self->rng, 1.0 / lambda);
        }

        *ret_t_wait = t_wait;
    }
out:
    return ret;
}

static double
msp_get_total_gc_left(msp_t *self)
{
    double total = 0;

    size_t num_ancestors = msp_get_num_ancestors(self);
    double mean_gc_rate = rate_map_get_total_mass(&self->gc_map) / self->sequence_length;
    total = (double) num_ancestors * mean_gc_rate * self->gc_tract_length;
    return total;
}

static segment_t *
msp_find_gc_left_individual(msp_t *self, label_id_t label, double value)
{
    size_t j, num_ancestors, individual_index;
    avl_tree_t *ancestors;
    avl_node_t *node;
    segment_t *ind;

    double mean_gc_rate = rate_map_get_total_mass(&self->gc_map) / self->sequence_length;
    individual_index = (size_t) floor(value / (mean_gc_rate * self->gc_tract_length));
    for (j = 0; j < self->num_populations; j++) {
        num_ancestors = msp_get_num_population_ancestors(self, (tsk_id_t) j);
        if (individual_index < num_ancestors) {
            ancestors = &self->populations[j].ancestors[label];
            /* Choose the correct individual */
            node = avl_at(ancestors, (unsigned int) individual_index);
            assert(node != NULL);
            ind = (segment_t *) node->item;
            return ind;
        } else {
            individual_index -= num_ancestors;
        }
    }
    return NULL;
}

static double
msp_get_total_gc_left_rate(msp_t *self)
{
    double mean_gc_rate = rate_map_get_total_mass(&self->gc_map) / self->sequence_length;
    double ret = 0;
    double total_gc_left;

    if (mean_gc_rate > 0) {
        total_gc_left = msp_get_total_gc_left(self);
        ret = total_gc_left;
    }
    return ret;
}

static int MSP_WARN_UNUSED
msp_sample_gc_left_waiting_time(msp_t *self, double *ret_t_wait)
{
    int ret = 0;
    double lambda = msp_get_total_gc_left_rate(self);
    double t_wait = DBL_MAX;

    if (lambda > 0.0) {
        t_wait = gsl_ran_exponential(self->rng, 1.0 / lambda);
    }
    *ret_t_wait = t_wait;
    return ret;
}

static int MSP_WARN_UNUSED
msp_gene_conversion_left_event(msp_t *self, label_id_t label)
{
    int ret = 0;
    const double gc_left_total = msp_get_total_gc_left(self);
    double h = gsl_rng_uniform(self->rng) * gc_left_total;
    double tl, bp;
    segment_t *y, *x, *alpha;

    y = msp_find_gc_left_individual(self, label, h);
    assert(y != NULL);

    /* generate tract length */
    tl = msp_generate_gc_tract_length(self);
    if (tl == -1) {
        ret = MSP_ERR_TRACTLEN_RESAMPLE_OVERFLOW;
        goto out;
    }
    tsk_bug_assert(tl > 0);
    self->num_internal_gc_events++;
    self->sum_internal_gc_tract_lengths += tl;

    bp = y->left + tl;

    while (y != NULL && y->right <= bp) {
        y = y->next;
    }

    if (y == NULL) {
        //   last segment
        // ... ==========   |
        //                  bp
        self->num_noneffective_gc_events++;
        return 0;
    }
    tsk_bug_assert(y != NULL);
    self->num_gc_events++;
    x = y->prev;

    if (y->left < bp) {
        //  x          y
        // =====   =====|====
        //              bp
        // becomes
        //  x         y
        // =====   =====
        //              =====
        //                α
        alpha = msp_copy_segment(self, y);
        if (alpha == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        alpha->left = bp;
        alpha->prev = NULL;
        if (alpha->next != NULL) {
            alpha->next->prev = alpha;
        }
        y->next = NULL;
        y->right = bp;
        msp_set_segment_mass(self, y);
        if (!msp_has_breakpoint(self, bp)) {
            ret = msp_insert_breakpoint(self, bp);
            if (ret != 0) {
                goto out;
            }
        }
    } else {
        //  x          y
        // ===== |  =========
        //       bp
        // becomes
        //  x
        // =====
        //          =========
        //              α
        x->next = NULL;
        y->prev = NULL;
        alpha = y;
        // Ensure y points to the last segment left of the break for full ARG recording
        y = x;
    }
    msp_set_segment_mass(self, alpha);
    tsk_bug_assert(alpha->prev == NULL);
    ret = msp_insert_individual(self, alpha);
    if (self->store_full_arg) {
        ret = msp_store_arg_gene_conversion(self, NULL, y, alpha);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

/* The main event loop for continuous time coalescent models. Runs until either
 * coalescence; or the time of a simulated event would have exceeded the
 * specified max_time; or for a specified number of events. The num_events
 * parameter is provided so that higher-level code can run the simulation for
 * smaller time chunks. This is used in the Python interface to check for
 * interrupts, so that long-running simulations can be killed using CTRL-C.
 *
 * Returns:
 * MSP_EXIT_COALESCENCE if the simulation completed to coalescence
 * MSP_EXIT_MAX_EVENTS if the simulation stopped because the maximum number
 *    of events was reached.
 * MSP_EXIT_MAX_TIME if the simulation stopped because the maximum time would
 *    have been exceeded by an event.
 * A negative value if an error occured.
 */
static int MSP_WARN_UNUSED
msp_run_coalescent(msp_t *self, double max_time, unsigned long max_events)
{
    int ret = 0;
    double lambda, t_temp, t_wait, ca_t_wait, re_t_wait, gc_t_wait, gc_left_t_wait,
        mig_t_wait, random_event_time, fixed_event_time;
    uint32_t n;
    tsk_id_t i, pop_id, pop_id_j, pop_id_k, ca_pop_id, mig_source_pop, mig_dest_pop;
    const tsk_id_t N = (tsk_id_t) self->num_populations;
    population_t *pop;
    unsigned long events = 0;
    avl_node_t *avl_node;
    /* Only support a single label for now. */
    label_id_t label = 0;

    ret = msp_compute_population_indexes(self);
    if (ret != 0) {
        goto out;
    }

    while (msp_get_num_ancestors(self) > 0) {
        if (events == max_events) {
            ret = MSP_EXIT_MAX_EVENTS;
            break;
        }
        events++;

        /* Recombination */
        ret = msp_sample_waiting_time(self, self->recomb_mass_index, label, &re_t_wait);
        if (ret != 0) {
            goto out;
        }

        /* Gene conversion */
        gc_t_wait = DBL_MAX;
        ret = msp_sample_waiting_time(self, self->gc_mass_index, label, &gc_t_wait);
        if (ret != 0) {
            goto out;
        }
        gc_left_t_wait = DBL_MAX;
        ret = msp_sample_gc_left_waiting_time(self, &gc_left_t_wait);
        if (ret != 0) {
            goto out;
        }

        /* Common ancestors */
        ca_t_wait = DBL_MAX;
        ca_pop_id = 0;
        for (avl_node = self->non_empty_populations.head; avl_node != NULL;
             avl_node = avl_node->next) {
            pop_id = (tsk_id_t)(intptr_t) avl_node->item;
            t_temp = self->get_common_ancestor_waiting_time(self, pop_id, label);
            if (t_temp < ca_t_wait) {
                ca_t_wait = t_temp;
                ca_pop_id = pop_id;
            }
        }

        /* Migration */
        mig_t_wait = DBL_MAX;
        mig_source_pop = 0;
        mig_dest_pop = 0;
        for (avl_node = self->non_empty_populations.head; avl_node != NULL;
             avl_node = avl_node->next) {
            pop_id_j = (tsk_id_t)(intptr_t) avl_node->item;
            pop = &self->populations[pop_id_j];
            n = avl_count(&pop->ancestors[label]);
            tsk_bug_assert(n > 0);
            for (i = 0; i < (tsk_id_t) pop->num_potential_destinations; i++) {
                pop_id_k = pop->potential_destinations[i];
                lambda = n * self->migration_matrix[pop_id_j * N + pop_id_k];
                tsk_bug_assert(lambda > 0);
                t_temp = gsl_ran_exponential(self->rng, 1.0 / lambda);
                if (t_temp < mig_t_wait) {
                    mig_t_wait = t_temp;
                    /* m[j, k] is the rate at which migrants move from
                     * population k to j forwards in time. Backwards
                     * in time, we move the individual from from
                     * population j into population k.
                     */
                    mig_source_pop = pop_id_j;
                    mig_dest_pop = pop_id_k;
                }
            }
        }

        fixed_event_time = msp_get_next_fixed_event_time(self);
        t_wait = GSL_MIN(mig_t_wait,
            GSL_MIN(gc_t_wait, GSL_MIN(gc_left_t_wait, GSL_MIN(re_t_wait, ca_t_wait))));

        if (fixed_event_time == DBL_MAX && t_wait == DBL_MAX) {
            ret = MSP_ERR_INFINITE_WAITING_TIME;
            goto out;
        }
        random_event_time = self->time + t_wait;

        /* The simulation state is can only changed from this point on. If
         * any of the events would cause the time to be >= max_time, we exit
         */
        if (fixed_event_time < random_event_time) {
            if (fixed_event_time > max_time) {
                ret = MSP_EXIT_MAX_TIME;
                break;
            }
            ret = msp_apply_fixed_events(self, fixed_event_time);
            if (ret != 0) {
                goto out;
            }
            /* Rather than try to reason about the changes that have occured
             * during any demographic events, just recompute the indexes
             * used to track nonempty populations and migration destinations
             */
            ret = msp_compute_population_indexes(self);
            if (ret != 0) {
                goto out;
            }
        } else {
            if (random_event_time > max_time) {
                ret = MSP_EXIT_MAX_TIME;
                break;
            }
            self->time = random_event_time;
            if (re_t_wait == t_wait) {
                ret = msp_recombination_event(self, label, NULL, NULL);
            } else if (gc_t_wait == t_wait) {
                ret = msp_gene_conversion_event(self, label);
            } else if (gc_left_t_wait == t_wait) {
                ret = msp_gene_conversion_left_event(self, label);
            } else if (ca_t_wait == t_wait) {
                ret = self->common_ancestor_event(self, ca_pop_id, label);
                if (ret == 1) {
                    /* The CA event has signalled that this event should be rejected */
                    self->time -= t_wait;
                    ret = 0;
                }
                if (ret != 0) {
                    goto out;
                }
                if (msp_get_num_population_ancestors(self, ca_pop_id) == 0) {
                    ret = msp_remove_non_empty_population(self, ca_pop_id);
                }
            } else {
                ret = msp_migration_event(self, mig_source_pop, mig_dest_pop);
                if (ret != 0) {
                    goto out;
                }
                if (msp_get_num_population_ancestors(self, mig_source_pop) == 0) {
                    ret = msp_remove_non_empty_population(self, mig_source_pop);
                }
                ret = msp_insert_non_empty_population(self, mig_dest_pop);
            }
            if (ret != 0) {
                goto out;
            }
        }
    }
out:
    return ret;
}

static int
msp_pedigree_process_common_ancestors(msp_t *self, individual_t *ind, tsk_size_t ploid)
{
    int ret = 0;
    tsk_id_t node = ind->nodes[ploid];
    tsk_id_t parent = ind->parents[ploid];
    avl_tree_t *common_ancestors = &ind->common_ancestors[ploid];
    segment_t *genome, *parent_ancestry[MSP_MAX_PED_PLOIDY], *seg;
    const tsk_size_t ploidy = self->ploidy;
    tsk_size_t j;

    /* FIXME - assuming 1 label here */
    ret = msp_merge_n_ancestors(
        self, common_ancestors, ind->population, 0, node, &genome);
    if (ret != 0) {
        goto out;
    }
    if (genome != NULL) {
        tsk_bug_assert(genome->prev == NULL);

        if (parent == TSK_NULL) {
            /* If parent is NULL, we are at a pedigree founder and create
             * unary edges to mark where the ancestral material left the
             * pedigree and can be picked up again by later simulations.
             *
             * NOTE arguably we should remove ``genome`` from the population
             * and decrement overlap counts so that we can detect early
             * coalescence. However, this makes it difficult to tell the
             * difference between a true full coalescence, and the ending
             * of the model. It's simpler to make finalise_tables a no-op
             * for the pedigree simulation to avoid adding the unary edges
             * up to the final simulation time.
             */
            for (seg = genome; seg != NULL; seg = seg->next) {
                if (seg->value != node) {
                    ret = msp_store_edge(self, seg->left, seg->right, node, seg->value);
                    if (ret != 0) {
                        goto out;
                    }
                }
            }
        } else {
            if (rate_map_get_total_mass(&self->recomb_map) > 0) {
                ret = msp_dtwf_recombine(
                    self, genome, &parent_ancestry[0], &parent_ancestry[1]);
                if (ret != 0) {
                    goto out;
                }
            } else {
                parent_ancestry[0] = NULL;
                parent_ancestry[1] = NULL;
                j = (tsk_size_t) gsl_rng_uniform_int(self->rng, 2);
                parent_ancestry[j] = genome;
            }
            for (j = 0; j < ploidy; j++) {
                seg = parent_ancestry[j];
                if (seg != NULL) {
                    tsk_bug_assert(seg->prev == NULL);
                    ret = msp_pedigree_add_individual_common_ancestor(
                        self, parent, seg, j);
                    if (ret != 0) {
                        goto out;
                    }
                    if (seg != genome) {
                        ret = msp_insert_individual(self, seg);
                        if (ret != 0) {
                            goto out;
                        }
                    }
                }
            }
        }
    }
    ret = msp_flush_edges(self);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

/* Main loop for simulating through a pedigree.
 *
 * Returns:
 * MSP_EXIT_COALESCENCE if the simulation completed to coalescence
 * MSP_EXIT_MAX_EVENTS if the simulation stopped because the maximum number
 *    of events was reached.
 * MSP_EXIT_MAX_TIME if the simulation stopped because the maximum time would
 *    have been exceeded by an event.
 * A negative value if an error occured.
 */
static int
msp_run_pedigree(msp_t *self, double max_time, unsigned long max_events)
{
    int ret = 0;
    tsk_size_t ploid;
    individual_t *ind = NULL;
    pedigree_t *pedigree = &self->pedigree;
    const tsk_id_t num_individuals = (tsk_id_t) pedigree->num_individuals;
    unsigned long num_events = 0;

    tsk_bug_assert(pedigree != NULL);
    if (pedigree->next_individual == PEDIGREE_UNINITIALISED) {
        ret = msp_pedigree_initialise(self);
        if (ret != 0) {
            goto out;
        }
    }
    /* Visit the individuals in time-increasing order. */
    while (pedigree->next_individual < num_individuals) {
        ind = pedigree->visit_order[pedigree->next_individual];
        /* To keep the semantics as close as possible to DTWF we simulate until
         * the *end* of generation max_time. So, we stop just before processing
         * the first individual in generation max_time + 1. */
        if (ind->time > max_time) {
            ret = MSP_EXIT_MAX_TIME;
            break;
        }
        if (num_events >= max_events) {
            ret = MSP_EXIT_MAX_EVENTS;
            break;
        }
        tsk_bug_assert(ind->time >= self->time);
        self->time = ind->time;
        ret = msp_pedigree_insert_ancient_samples(self);
        if (ret != 0) {
            goto out;
        }
        for (ploid = 0; ploid < ind->ploidy; ploid++) {
            ret = msp_pedigree_process_common_ancestors(self, ind, ploid);
            if (ret != 0) {
                goto out;
            }
        }
        pedigree->next_individual++;
        num_events++;
    }
    if (msp_get_num_ancestors(self) == 0) {
        ret = MSP_EXIT_COALESCENCE;
    } else if (pedigree->next_individual == num_individuals) {
        ret = MSP_EXIT_MODEL_COMPLETE;
    } else {
        tsk_bug_assert(ret == MSP_EXIT_MAX_TIME || ret == MSP_EXIT_MAX_EVENTS);
    }
out:
    return ret;
}

/* List structure for collecting segments by parent */
typedef struct _segment_list_t {
    avl_node_t *node;
    struct _segment_list_t *next;
} segment_list_t;

/* Performs a single generation under the Wright Fisher model */
static int MSP_WARN_UNUSED
msp_dtwf_generation(msp_t *self)
{
    int ret = 0;
    int ix;
    uint32_t N, i, j, k, p;
    size_t segment_mem_offset;
    population_t *pop;
    segment_t *x, *u[2];
    segment_list_t **parents = NULL;
    segment_list_t *segment_mem = NULL;
    segment_list_t *s;
    avl_node_t *a, *node;
    avl_tree_t Q[2];
    /* Only support single structured coalescent label for now. */
    label_id_t label = 0;

    for (i = 0; i < 2; i++) {
        avl_init_tree(&Q[i], cmp_segment_queue, NULL);
    }

    for (j = 0; j < self->num_populations; j++) {

        pop = &self->populations[j];
        if (avl_count(&pop->ancestors[label]) == 0) {
            continue;
        }
        /* For the DTWF, N for each population is the reference population size
         * from the model multiplied by the current population size, rounded to
         * the nearest integer. Thus, the population's size is always relative
         * to the reference model population size (which is also true for the
         * coalescent models. */
        N = (uint32_t) round(get_population_size(pop, self->time));
        if (N == 0) {
            ret = MSP_ERR_DTWF_ZERO_POPULATION_SIZE;
            goto out;
        }

        // Allocate memory for linked list of offspring per parent
        parents = calloc(N, sizeof(segment_list_t *));
        segment_mem = malloc(msp_get_num_ancestors(self) * sizeof(segment_list_t));
        if (parents == NULL || segment_mem == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        // Iterate through ancestors and draw parents
        segment_mem_offset = 0;
        for (a = pop->ancestors[label].head; a != NULL; a = a->next) {
            s = segment_mem + segment_mem_offset;
            segment_mem_offset++;
            p = (uint32_t) gsl_rng_uniform_int(self->rng, N);
            if (parents[p] != NULL) {
                self->num_ca_events++;
            }
            s->next = parents[p];
            s->node = a;
            parents[p] = s;
        }

        // Iterate through offspring of parent k, adding to avl_tree
        for (k = 0; k < N; k++) {
            for (s = parents[k]; s != NULL; s = s->next) {
                node = s->node;
                x = (segment_t *) node->item;
                // Recombine ancestor
                // TODO Should this be the recombination rate going foward from x.left?
                if (rate_map_get_total_mass(&self->recomb_map) > 0) {
                    ret = msp_dtwf_recombine(self, x, &u[0], &u[1]);
                    if (ret != 0) {
                        goto out;
                    }
                    for (i = 0; i < 2; i++) {
                        if (u[i] != NULL && u[i] != x) {
                            ret = msp_insert_individual(self, u[i]);
                            if (ret != 0) {
                                goto out;
                            }
                        }
                    }
                } else {
                    ix = (int) gsl_rng_uniform_int(self->rng, 2);
                    u[0] = NULL;
                    u[1] = NULL;
                    u[ix] = x;
                }
                // Add to AVLTree for each parental chromosome
                for (i = 0; i < 2; i++) {
                    if (u[i] != NULL) {
                        ret = msp_priority_queue_insert(self, &Q[i], u[i]);
                        if (ret != 0) {
                            goto out;
                        }
                    }
                }
            }
            // Merge segments in each parental chromosome
            for (i = 0; i < 2; i++) {
                ret = msp_merge_n_ancestors(
                    self, &Q[i], (population_id_t) j, label, TSK_NULL, NULL);
                if (ret != 0) {
                    goto out;
                }
            }
        }
        free(parents);
        free(segment_mem);
        segment_mem = NULL;
        parents = NULL;
    }
out:
    msp_safe_free(parents);
    msp_safe_free(segment_mem);
    return ret;
}

static int MSP_WARN_UNUSED
msp_store_simultaneous_migration_events(
    msp_t *self, avl_tree_t *nodes, population_id_t source_pop, label_id_t label)
{
    int ret = 0;
    uint32_t j;
    avl_node_t *node;
    avl_tree_t *source;

    source = &self->populations[source_pop].ancestors[label];

    // Choose node to migrate
    j = (uint32_t) gsl_rng_uniform_int(self->rng, avl_count(source));
    node = avl_at(source, j);
    tsk_bug_assert(node != NULL);

    avl_unlink_node(source, node);
    node = avl_insert_node(nodes, node);
    tsk_bug_assert(node != NULL);

    return ret;
}

static int MSP_WARN_UNUSED
msp_simultaneous_migration_event(
    msp_t *self, avl_tree_t *nodes, population_id_t source_pop, population_id_t dest_pop)
{
    int ret = 0;
    avl_node_t *node;
    size_t index;
    label_id_t label = 0; /* For now only support label 0 */

    index = ((size_t) source_pop) * self->num_populations + (size_t) dest_pop;
    self->num_migration_events[index]++;

    // Iterate through nodes in tree and move to new pop
    // -- removed from "nodes" in msp_move_individual()
    while (avl_count(nodes) > 0) {
        node = nodes->head;
        ret = msp_move_individual(self, node, nodes, dest_pop, label);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

/* The main event loop for the Wright Fisher model.
 *
 * Returns:
 * MSP_EXIT_COALESCENCE if the simulation completed to coalescence
 * MSP_EXIT_MAX_EVENTS if the simulation stopped because the maximum number
 *    of events was reached.
 * MSP_EXIT_MAX_TIME if the simulation stopped because the maximum time would
 *    have been exceeded by an event.
 * A negative value if an error occured.
 */
static int MSP_WARN_UNUSED
msp_run_dtwf(msp_t *self, double max_time, unsigned long max_events)
{
    int ret = 0;
    unsigned long events = 0;
    int mig_source_pop, mig_dest_pop;
    sampling_event_t *se;
    uint32_t j, k, i, N;
    unsigned int *n = NULL;
    double *mig_tmp = NULL;
    double sum, cur_time;
    avl_tree_t *node_trees = NULL;
    avl_tree_t *nodes;
    /* Only support a single structured coalescent label at the moment */
    label_id_t label = 0;

    tsk_bug_assert(self->recomb_mass_index == NULL);
    tsk_bug_assert(self->gc_mass_index == NULL);
    if (rate_map_get_total_mass(&self->gc_map) != 0.0) {
        /* Could be, we just haven't implemented it */
        ret = MSP_ERR_DTWF_GC_NOT_SUPPORTED;
        goto out;
    }
    if (self->ploidy != 2) {
        ret = MSP_ERR_DTWF_DIPLOID_ONLY;
        goto out;
    }

    n = malloc(self->num_populations * sizeof(int));
    mig_tmp = malloc(self->num_populations * sizeof(double));
    if (n == NULL || mig_tmp == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    while (msp_get_num_ancestors(self) > 0) {
        if (events == max_events) {
            ret = MSP_EXIT_MAX_EVENTS;
            break;
        }
        events++;
        if (self->time + 1 > max_time) {
            ret = MSP_EXIT_MAX_TIME;
            goto out;
        }
        self->time++;

        /* Following SLiM, we perform migrations prior to selecting
         * parents for the current generation */
        node_trees
            = malloc(self->num_populations * self->num_populations * sizeof(avl_tree_t));
        if (node_trees == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }

        mig_source_pop = 0;
        mig_dest_pop = 0;
        for (j = 0; j < self->num_populations; j++) {
            // For proper sampling, we need to calculate the proportion
            // of non-migrants as well
            sum = 0;
            for (k = 0; k < self->num_populations; k++) {
                mig_tmp[k] = self->migration_matrix[j * self->num_populations + k];
                sum += mig_tmp[k];
            }
            tsk_bug_assert(mig_tmp[j] == 0);

            // Must check that row sums of migration matrix are <=1 in the main
            // loop, as multiple indices can change in the same generation
            if (sum > 1) {
                ret = MSP_ERR_DTWF_MIGRATION_MATRIX_NOT_STOCHASTIC;
                goto out;
            }

            mig_tmp[j] = 1 - sum;
            N = avl_count(&self->populations[j].ancestors[label]);
            gsl_ran_multinomial(self->rng, self->num_populations, N, mig_tmp, n);

            for (k = 0; k < self->num_populations; k++) {
                if (k == j) {
                    continue;
                }
                // Initialize an avl tree for this pair of populations
                nodes = &node_trees[j * self->num_populations + k];
                avl_init_tree(nodes, cmp_individual, NULL);

                /* m[j, k] is the rate at which migrants move from
                 * population k to j forwards in time. Backwards
                 * in time, we move the individual from from
                 * population j into population k.
                 */
                mig_source_pop = (population_id_t) j;
                mig_dest_pop = (population_id_t) k;

                for (i = 0; i < n[k]; i++) {
                    ret = msp_store_simultaneous_migration_events(
                        self, nodes, mig_source_pop, label);
                    if (ret != 0) {
                        goto out;
                    }
                }
            }
        }
        for (j = 0; j < self->num_populations; j++) {
            for (k = 0; k < self->num_populations; k++) {
                if (k == j) {
                    continue;
                }
                nodes = &node_trees[j * self->num_populations + k];
                mig_source_pop = (population_id_t) j;
                mig_dest_pop = (population_id_t) k;
                ret = msp_simultaneous_migration_event(
                    self, nodes, mig_source_pop, mig_dest_pop);
                if (ret != 0) {
                    goto out;
                }
            }
        }
        free(node_trees);
        node_trees = NULL;

        /* Demographic events set the simulation time to the time of the event.
         * In the DTWF, this would prevent more than one event occurring per
         * generation, and throw off the time between generations. We avoid
         * this by saving the current time and returning to it. */
        cur_time = self->time;
        while (self->next_demographic_event != NULL
               && self->next_demographic_event->time <= cur_time) {
            if (self->next_demographic_event->time > max_time) {
                ret = MSP_EXIT_MAX_TIME;
                goto out;
            }
            ret = msp_apply_demographic_events(self, self->next_demographic_event->time);
            if (ret != 0) {
                goto out;
            }
        }
        self->time = cur_time;
        ret = msp_dtwf_generation(self);
        if (ret != 0) {
            goto out;
        }

        while (self->next_sampling_event < self->num_sampling_events
               && self->sampling_events[self->next_sampling_event].time <= self->time) {
            se = self->sampling_events + self->next_sampling_event;
            /* The sampling event doesn't modify the tables, so we don't need to
             * catch it here */
            ret = msp_insert_sample(self, se->sample);
            if (ret != 0) {
                goto out;
            }
            self->next_sampling_event++;
        }
    }
out:
    msp_safe_free(node_trees);
    msp_safe_free(n);
    msp_safe_free(mig_tmp);
    return ret;
}

/* Set up the intial populations for a sweep by moving individuals from
 * label 0 to label 1 with the specified probability.
 */
static int
msp_sweep_initialise(msp_t *self, double switch_proba)
{
    int ret = 0;
    uint32_t j;
    avl_node_t *node, *next;
    avl_tree_t *pop;

    /* We only support one population and two labels for now */
    if (self->num_populations != 1 || self->num_labels != 2) {
        ret = MSP_ERR_UNSUPPORTED_OPERATION;
        goto out;
    }

    /* Move ancestors to new labels. */
    for (j = 0; j < self->num_populations; j++) {
        tsk_bug_assert(avl_count(&self->populations[j].ancestors[1]) == 0);
        pop = &self->populations[j].ancestors[0];
        node = pop->head;
        while (node != NULL) {
            next = node->next;
            if (gsl_rng_uniform(self->rng) < switch_proba) {
                ret = msp_move_individual(self, node, pop, (population_id_t) j, 1);
                if (ret != 0) {
                    goto out;
                }
            }
            node = next;
        }
    }
out:
    return ret;
}

/* Finalise the sweep by moving all lineages back to label 0.
 */
static int
msp_sweep_finalise(msp_t *self)
{
    int ret = 0;
    uint32_t j;
    avl_node_t *node, *next;
    avl_tree_t *pop;

    /* Move ancestors to new labels. */
    for (j = 0; j < self->num_populations; j++) {
        pop = &self->populations[j].ancestors[1];
        node = pop->head;
        while (node != NULL) {
            next = node->next;
            ret = msp_move_individual(self, node, pop, (population_id_t) j, 0);
            if (ret != 0) {
                goto out;
            }
            node = next;
        }
    }
out:
    return ret;
}

/* Migrate the specified individual to the specified new label.
 */
static int
msp_change_label(msp_t *self, segment_t *ind, label_id_t label)
{
    int ret = 0;
    avl_tree_t *pop = &self->populations[ind->population].ancestors[ind->label];
    avl_node_t *node;

    /* Find the this individual in the AVL tree. */
    node = avl_search(pop, ind);
    tsk_bug_assert(node != NULL);
    ret = msp_move_individual(self, node, pop, ind->population, label);
    return ret;
}

static int
msp_sweep_recombination_event(
    msp_t *self, label_id_t label, double sweep_locus, double population_frequency)
{
    int ret = 0;
    segment_t *lhs, *rhs;
    label_id_t new_label;
    double r;

    ret = msp_recombination_event(self, label, &lhs, &rhs);
    if (ret != 0) {
        goto out;
    }
    /* NOTE: we can look at rhs->left when we compare to the sweep site. */
    r = gsl_rng_uniform(self->rng);
    if (sweep_locus < rhs->left) {
        if (r < 1.0 - population_frequency) {
            /* move rhs to other population */
            new_label = (label + 1) % 2;
            ret = msp_change_label(self, rhs, new_label);
            if (ret != 0) {
                goto out;
            }
        }
    } else {
        if (r < 1.0 - population_frequency) {
            /* move lhs to other population */
            new_label = (label + 1) % 2;
            ret = msp_change_label(self, lhs, new_label);
            if (ret != 0) {
                goto out;
            }
        }
    }
out:
    return ret;
}

static int
msp_run_sweep(msp_t *self)
{
    int ret = 0;
    simulation_model_t *model = &self->model;
    size_t curr_step = 1;
    size_t num_steps;
    double *allele_frequency = NULL;
    double *time = NULL;
    double sweep_locus = model->params.sweep.position;
    double sweep_dt;
    size_t j = 0;
    double recomb_mass;
    unsigned long events = 0;
    label_id_t label;
    double rec_rates[] = { 0.0, 0.0 };
    double sweep_pop_sizes[] = { 0.0, 0.0 };
    double event_prob, event_rand, tmp_rand, e_sum, pop_size;
    double p_coal_b, p_coal_B, total_rate, sweep_pop_tot_rate;
    double p_rec_b, p_rec_B;
    double t_start, t_unscaled;
    bool sweep_over;

    if (rate_map_get_total_mass(&self->gc_map) != 0.0) {
        /* Could be, we just haven't implemented it */
        ret = MSP_ERR_SWEEPS_GC_NOT_SUPPORTED;
        goto out;
    }

    /* Keep the compiler happy */
    sweep_pop_tot_rate = 0;
    p_coal_b = 0;
    p_coal_B = 0;
    p_rec_b = 0;
    p_rec_B = 0;
    /* JK: I've removed the time and event limits on this function to simplify
     * things as function is currently 'non-rentrant'; we can't stop in the middle
     * of the sweep and pick it up again from where we left off later. This is
     * different to all the other model runners. Probably the simplest thing to
     * do here is to make the trajectory times absolute wrt to the simulation, and
     * to then increment curr_step accordingly. We can then probably reason about
     * when the call msp_sweep_initialise and msp_sweep_finalise
     * depending on the value of curr_step, and hopefully reintroduce the max_time
     * and max_steps parameters. */

    ret = model->params.sweep.generate_trajectory(
        &model->params.sweep, self, &num_steps, &time, &allele_frequency);

    t_start = self->time;
    sweep_dt = model->params.sweep.trajectory_params.genic_selection_trajectory.dt;
    tsk_bug_assert(sweep_dt > 0);
    if (ret != 0) {
        goto out;
    }
    ret = msp_sweep_initialise(self, allele_frequency[0]);
    if (ret != 0) {
        goto out;
    }

    curr_step = 1;
    while (msp_get_num_ancestors(self) > 0 && curr_step < num_steps) {
        events++;
        /* Set pop sizes & rec_rates */
        for (j = 0; j < self->num_labels; j++) {
            label = (label_id_t) j;
            recomb_mass = self->recomb_mass_index == NULL
                              ? 0
                              : fenwick_get_total(&self->recomb_mass_index[label]);
            sweep_pop_sizes[j] = avl_count(&self->populations[0].ancestors[label]);
            /* We can get small negative rates by numerical jitter which causes
             * problems in later calculations and leads to assertion trips.
             * See https://github.com/tskit-dev/msprime/issues/1966
             */
            rec_rates[j] = TSK_MAX(0, recomb_mass);
        }

        event_prob = 1.0;
        event_rand = gsl_rng_uniform(self->rng);
        sweep_over = false;
        while (event_prob > event_rand && curr_step < num_steps && !sweep_over) {
            pop_size = get_population_size(&self->populations[0], self->time);
            p_coal_B = 0;
            if (avl_count(&self->populations[0].ancestors[1]) > 1) {
                p_coal_B = ((sweep_pop_sizes[1] * (sweep_pop_sizes[1] - 1)) * 0.5)
                           / allele_frequency[curr_step] * sweep_dt;
            }
            p_coal_b = 0;
            if (avl_count(&self->populations[0].ancestors[0]) > 1) {
                p_coal_b = ((sweep_pop_sizes[0] * (sweep_pop_sizes[0] - 1)) * 0.5)
                           / (1.0 - allele_frequency[curr_step]) * sweep_dt;
            }
            p_rec_b = rec_rates[0] * pop_size * self->ploidy * sweep_dt;
            p_rec_B = rec_rates[1] * pop_size * self->ploidy * sweep_dt;
            sweep_pop_tot_rate = p_coal_b + p_coal_B + p_rec_b + p_rec_B;
            /* doing this to build in generality if we want >1 pop */

            total_rate = sweep_pop_tot_rate;
            /* DEBUG
            printf("pop_size: %g sweep_dt: %g sweep_pop_sizes[1]: %g "
                   "sweep_pop_sizes[1]: %g\n",
                pop_size, sweep_dt, sweep_pop_sizes[0], sweep_pop_sizes[1]);
            printf("x: %g p_rec_b: %g p_rec_B: %g p_coal_b: %g p_coal_B: %g\n",
                allele_frequency[curr_step], p_rec_b, p_rec_B, p_coal_b, p_coal_B);
            printf("rec_rates[0]: %g rec_rates[1]: %g ploidy: %d\n", rec_rates[0],
                rec_rates[1], self->ploidy);
            printf("event_prob: %g rand: %g\n", event_prob, event_rand);
            */
            event_prob *= 1.0 - total_rate;
            curr_step++;

            sweep_over = total_rate == 0;
        }
        if (sweep_over) {
            break;
        }

        tmp_rand = gsl_rng_uniform(self->rng);

        e_sum = p_coal_b;
        /* convert time scale */
        pop_size = get_population_size(&self->populations[0], self->time);
        t_unscaled = time[curr_step - 1] * self->ploidy * pop_size;
        tsk_bug_assert(t_unscaled > 0);
        self->time = t_start + t_unscaled;
        /* printf("event time: %g\n", self->time); */
        if (tmp_rand < e_sum / sweep_pop_tot_rate) {
            /* coalescent in b background */
            ret = self->common_ancestor_event(self, 0, 0);
        } else {
            e_sum += p_coal_B;
            if (tmp_rand < e_sum / sweep_pop_tot_rate) {
                /* coalescent in B background */
                ret = self->common_ancestor_event(self, 0, 1);
            } else {
                e_sum += p_rec_b;
                if (tmp_rand < e_sum / sweep_pop_tot_rate) {
                    /* recomb in b background */
                    ret = msp_sweep_recombination_event(
                        self, 0, sweep_locus, (1.0 - allele_frequency[curr_step - 1]));
                } else {
                    /* recomb in B background */
                    ret = msp_sweep_recombination_event(
                        self, 1, sweep_locus, allele_frequency[curr_step - 1]);
                }
            }
        }
        if (ret != 0) {
            goto out;
        }
        /* msp_print_state(self, stdout); */
    }

    /* TODO we should probably support fixed events here using
     * msp_apply_fixed_events() like we do in the coalescent models.
     * The point below about computing population sizes should be easily
     * worked around using msp_compute_population_size(). */

    /* Check if any demographic events should have happened during the
     * event and raise an error if so. This is to keep computing population
     * sizes simple */
    if (self->next_demographic_event != NULL
        && self->next_demographic_event->time <= self->time) {
        ret = MSP_ERR_EVENTS_DURING_SWEEP;
        goto out;
    }
    /* Rule out sampling events for now, but there's no reason we can't
     * support them. */
    if (self->next_sampling_event < self->num_sampling_events
        && self->sampling_events[self->next_sampling_event].time <= self->time) {
        ret = MSP_ERR_EVENTS_DURING_SWEEP;
        goto out;
    }
    ret = msp_sweep_finalise(self);
    if (ret != 0) {
        goto out;
    }
    ret = MSP_EXIT_MODEL_COMPLETE;
out:
    msp_safe_free(time);
    msp_safe_free(allele_frequency);
    return ret;
}

/* Runs the simulation backwards in time until either the sample has coalesced,
 * or specified maximum simulation time has been reached or the specified maximum
 * number of events has been reached.
 */
int MSP_WARN_UNUSED
msp_run(msp_t *self, double max_time, unsigned long max_events)
{
    int ret = 0;
    int err;

    if (self->state == MSP_STATE_INITIALISED) {
        self->state = MSP_STATE_SIMULATING;
    }
    if (self->state != MSP_STATE_SIMULATING) {
        ret = MSP_ERR_BAD_STATE;
        goto out;
    }
    if (self->store_full_arg
        && !(self->model.type == MSP_MODEL_HUDSON || self->model.type == MSP_MODEL_SMC
               || self->model.type == MSP_MODEL_SMC_PRIME
               || self->model.type == MSP_MODEL_BETA
               || self->model.type == MSP_MODEL_DIRAC)) {
        /* We currently only support the full ARG recording on the standard
         * coalescent, SMC, or multiple merger coalescents. */
        ret = MSP_ERR_UNSUPPORTED_OPERATION;
        goto out;
    }

    if (msp_is_completed(self)) {
        /* If the simulation is completed, run() is a no-op for
         * all models. */
        ret = 0;
    } else if (self->model.type == MSP_MODEL_DTWF) {
        ret = msp_run_dtwf(self, max_time, max_events);
    } else if (self->model.type == MSP_MODEL_WF_PED) {
        ret = msp_run_pedigree(self, max_time, max_events);
    } else if (self->model.type == MSP_MODEL_SWEEP) {
        /* FIXME making sweep atomic for now as it's non-rentrant */
        ret = msp_run_sweep(self);
    } else {
        ret = msp_run_coalescent(self, max_time, max_events);
    }

    if (ret < 0) {
        goto out;
    }
    if (ret == MSP_EXIT_MAX_TIME) {
        /* Set the time to the max_time specified. If the tables are finalised
         * after this we will get unary edges on the end of each extant node
         * to this point so that the simulation can be resumed accurately.
         */
        self->time = max_time;
    }
    err = msp_flush_edges(self);
    if (err != 0) {
        ret = err;
        goto out;
    }
out:
    return ret;
}

/* Add in nodes and edges for the remaining segments to the output table. */
static int MSP_WARN_UNUSED
msp_insert_uncoalesced_edges(msp_t *self)
{
    int ret = 0;
    population_id_t pop;
    label_id_t label;
    avl_node_t *a;
    segment_t *seg;
    tsk_id_t node;
    int64_t edge_start;
    tsk_node_table_t *nodes = &self->tables->nodes;
    const double current_time = self->time;
    tsk_bookmark_t bookmark;

    for (pop = 0; pop < (population_id_t) self->num_populations; pop++) {
        for (label = 0; label < (label_id_t) self->num_labels; label++) {
            for (a = self->populations[pop].ancestors[label].head; a != NULL;
                 a = a->next) {
                /* If there are any nodes in the segment chain with the current time,
                 * then we don't make any unary edges for them. This is because (a)
                 * we'd end up edges with the same parent and child time (if we didn't
                 * hack an extra epsilon onto the parent time) and (b), this node
                 * could only have arisen as the result of a coalescence and so this
                 * node really does represent the current ancestor */
                node = TSK_NULL;
                for (seg = (segment_t *) a->item; seg != NULL; seg = seg->next) {
                    if (nodes->time[seg->value] == current_time) {
                        node = seg->value;
                        break;
                    }
                }
                if (node == TSK_NULL) {
                    /* Add a node for this ancestor */
                    node = tsk_node_table_add_row(
                        nodes, 0, current_time, pop, TSK_NULL, NULL, 0);
                    if (node < 0) {
                        ret = msp_set_tsk_error(node);
                        goto out;
                    }
                }

                /* For every segment add an edge pointing to this new node */
                for (seg = (segment_t *) a->item; seg != NULL; seg = seg->next) {
                    if (seg->value != node) {
                        tsk_bug_assert(nodes->time[node] > nodes->time[seg->value]);
                        ret = tsk_edge_table_add_row(&self->tables->edges, seg->left,
                            seg->right, node, seg->value, NULL, 0);
                        if (ret < 0) {
                            ret = msp_set_tsk_error(ret);
                            goto out;
                        }
                    }
                }
            }
        }
    }

    /* Find the first edge with parent == current time */
    edge_start = ((int64_t) self->tables->edges.num_rows) - 1;
    while (edge_start >= 0
           && nodes->time[self->tables->edges.parent[edge_start]] == current_time) {
        edge_start--;
    }
    memset(&bookmark, 0, sizeof(bookmark));
    if (edge_start > 0) {
        bookmark.edges = (tsk_size_t) edge_start;
    }
    /* We know migrations are sorted already */
    bookmark.migrations = self->tables->migrations.num_rows;
    /* Don't sort individuals */
    bookmark.individuals = self->tables->individuals.num_rows;
    ret = tsk_table_collection_sort(self->tables, &bookmark, 0);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
out:
    return ret;
}

int MSP_WARN_UNUSED
msp_finalise_tables(msp_t *self)
{
    int ret = 0;
    tsk_bookmark_t bookmark;

    /* We don't want to add unary edges for the pedigree simulation model */
    if (!msp_is_completed(self) && self->model.type != MSP_MODEL_WF_PED) {
        ret = msp_insert_uncoalesced_edges(self);
        if (ret != 0) {
            goto out;
        }
    }
    ret = tsk_table_collection_build_index(self->tables, 0);
    if (ret == TSK_ERR_EDGES_NOT_SORTED_PARENT_TIME) {
        /* There are rare cases when we are simulating from awkward initial
         * states where the tables we produce in the simulation are not
         * correctly sorted. The simplest course of action here to just let it
         * fail and sort. https://github.com/tskit-dev/msprime/issues/1606
         */
        memset(&bookmark, 0, sizeof(bookmark));
        /* Don't sort migrations or individuals */
        bookmark.migrations = self->tables->migrations.num_rows;
        bookmark.individuals = self->tables->individuals.num_rows;
        ret = tsk_table_collection_sort(self->tables, &bookmark, 0);
        if (ret != 0) {
            goto out;
        }
        ret = tsk_table_collection_build_index(self->tables, 0);
    }
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int
msp_debug_demography(msp_t *self, double *end_time)
{
    int ret = 0;
    double t = GSL_POSINF;
    int first_call = 0;
    demographic_event_t *de;
    sampling_event_t *se;

    if (self->state == MSP_STATE_INITIALISED) {
        self->state = MSP_STATE_DEBUGGING;
        first_call = 1;
    }
    if (self->state != MSP_STATE_DEBUGGING) {
        ret = MSP_ERR_BAD_STATE;
        goto out;
    }
    if (!first_call && self->next_demographic_event != NULL) {
        de = self->next_demographic_event;

        /* Add in historical samples more recent than next demographic event */
        while (self->next_sampling_event < self->num_sampling_events
               && self->sampling_events[self->next_sampling_event].time <= de->time) {
            se = self->sampling_events + self->next_sampling_event;
            ret = msp_insert_sample(self, se->sample);
            if (ret != 0) {
                goto out;
            }
            self->next_sampling_event++;
        }

        ret = msp_apply_demographic_events(self, de->time);
        if (ret != 0) {
            goto out;
        }
    }
    if (self->next_demographic_event != NULL) {
        t = self->next_demographic_event->time;
    }
    *end_time = t;
out:
    return ret;
}

/* Used for high-level debugging. */
int MSP_WARN_UNUSED
msp_compute_population_size(
    msp_t *self, size_t population_id, double time, double *pop_size)
{
    int ret = 0;
    population_t *pop;
    double dt;

    if (population_id > self->num_populations) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    pop = &self->populations[population_id];
    if (pop->growth_rate == 0.0) {
        *pop_size = pop->initial_size;
    } else {
        dt = time - pop->start_time;
        *pop_size = pop->initial_size * exp(-pop->growth_rate * dt);
    }
out:
    return ret;
}

simulation_model_t *
msp_get_model(msp_t *self)
{
    return &self->model;
}

const char *
msp_get_model_name(msp_t *self)
{
    const char *ret;

    switch (self->model.type) {
        case MSP_MODEL_HUDSON:
            ret = "hudson";
            break;
        case MSP_MODEL_SMC:
            ret = "smc";
            break;
        case MSP_MODEL_SMC_PRIME:
            ret = "smc_prime";
            break;
        case MSP_MODEL_DIRAC:
            ret = "dirac";
            break;
        case MSP_MODEL_BETA:
            ret = "beta";
            break;
        case MSP_MODEL_DTWF:
            ret = "dtwf";
            break;
        case MSP_MODEL_WF_PED:
            ret = "fixed_pedigree";
            break;
        case MSP_MODEL_SWEEP:
            ret = "single-sweep";
            break;
        default:
            ret = "BUG: bad model in simulator!";
            break;
    }
    return ret;
}

bool
msp_get_store_migrations(msp_t *self)
{
    return self->store_migrations;
}

size_t
msp_get_num_populations(msp_t *self)
{
    return (size_t) self->num_populations;
}

size_t
msp_get_num_labels(msp_t *self)
{
    return (size_t) self->num_labels;
}

size_t
msp_get_num_population_ancestors(msp_t *self, tsk_id_t population)
{
    tsk_id_t label;
    const population_t *pop = &self->populations[population];
    size_t n = 0;

    for (label = 0; label < (tsk_id_t) self->num_labels; label++) {
        n += avl_count(&pop->ancestors[label]);
    }
    return n;
}

size_t
msp_get_num_ancestors(msp_t *self)
{
    size_t n = 0;
    tsk_id_t j;

    for (j = 0; j < (tsk_id_t) self->num_populations; j++) {
        n += msp_get_num_population_ancestors(self, j);
    }
    return n;
}

size_t
msp_get_num_breakpoints(msp_t *self)
{
    return avl_count(&self->breakpoints);
}

size_t
msp_get_num_nodes(msp_t *self)
{
    return (size_t) self->tables->nodes.num_rows;
}

size_t
msp_get_num_edges(msp_t *self)
{
    return (size_t) self->tables->edges.num_rows;
}

size_t
msp_get_num_migrations(msp_t *self)
{
    return (size_t) self->tables->migrations.num_rows;
}

int MSP_WARN_UNUSED
msp_get_ancestors(msp_t *self, segment_t **ancestors)
{
    int ret = -1;
    avl_node_t *node;
    avl_tree_t *population_ancestors;
    size_t j;
    label_id_t label;
    size_t k = 0;

    for (j = 0; j < self->num_populations; j++) {
        for (label = 0; label < (label_id_t) self->num_labels; label++) {
            population_ancestors = &self->populations[j].ancestors[label];
            for (node = population_ancestors->head; node != NULL; node = node->next) {
                ancestors[k] = (segment_t *) node->item;
                k++;
            }
        }
    }
    ret = 0;
    return ret;
}

int MSP_WARN_UNUSED
msp_get_breakpoints(msp_t *self, size_t *breakpoints)
{
    int ret = -1;
    avl_node_t *node;
    node_mapping_t *nm;
    size_t j = 0;

    for (node = (&self->breakpoints)->head; node != NULL; node = node->next) {
        nm = (node_mapping_t *) node->item;
        breakpoints[j] = (size_t) nm->position;
        j++;
    }
    ret = 0;
    return ret;
}

int MSP_WARN_UNUSED
msp_get_migration_matrix(msp_t *self, double *migration_matrix)
{
    size_t N = self->num_populations;
    size_t j;

    for (j = 0; j < N * N; j++) {
        migration_matrix[j] = self->migration_matrix[j];
    }
    return 0;
}

int MSP_WARN_UNUSED
msp_get_num_migration_events(msp_t *self, size_t *num_migration_events)
{
    size_t N = self->num_populations;

    memcpy(num_migration_events, self->num_migration_events, N * N * sizeof(size_t));
    return 0;
}

int MSP_WARN_UNUSED
msp_get_population_configuration(msp_t *self, size_t population_id, double *initial_size,
    double *growth_rate, int *state)
{
    int ret = 0;
    population_t *pop;

    if (population_id > self->num_populations) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    pop = &self->populations[population_id];
    *initial_size = pop->initial_size;
    *growth_rate = pop->growth_rate;
    *state = pop->state;
out:
    return ret;
}

double
msp_get_time(msp_t *self)
{
    return self->time;
}

/* Demographic events. All times and input parameters are specified in units
 * of generations. When we store these values, we must rescale them into
 * model time, as appropriate. */
static int MSP_WARN_UNUSED
msp_add_demographic_event(msp_t *self, double time, demographic_event_t **event)
{
    int ret = MSP_ERR_GENERIC;
    demographic_event_t *ret_event = NULL;

    if (time < 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }

    if (self->demographic_events_tail != NULL) {
        if (time < self->demographic_events_tail->time) {
            ret = MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS;
            goto out;
        }
    }
    ret_event = calloc(1, sizeof(demographic_event_t));
    if (ret_event == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret_event->time = time;
    /* now insert this event at the end of the chain. */
    if (self->demographic_events_head == NULL) {
        self->demographic_events_head = ret_event;
        self->demographic_events_tail = ret_event;
    } else {
        tsk_bug_assert(self->demographic_events_tail != NULL);
        self->demographic_events_tail->next = ret_event;
        self->demographic_events_tail = ret_event;
    }
    *event = ret_event;
    ret = 0;
out:
    return ret;
}

/* Population parameter change */

static int
msp_change_single_population_parameters(msp_t *self, size_t population_id, double time,
    double initial_size, double growth_rate)
{
    int ret = 0;
    double dt;
    population_t *pop;

    if (population_id >= self->num_populations) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    pop = &self->populations[population_id];
    /* If initial_size is not specified, calculate the initial_size of the
     * population over the coming time period based on the growth rate over
     * the preceding period.
     */
    if (gsl_isnan(initial_size)) {
        dt = time - pop->start_time;
        pop->initial_size = pop->initial_size * exp(-pop->growth_rate * dt);
    } else {
        pop->initial_size = initial_size;
    }
    /* Do not change the growth_rate unless it is specified */
    if (!gsl_isnan(growth_rate)) {
        pop->growth_rate = growth_rate;
    }
    pop->start_time = time;
out:
    return ret;
}

static int
msp_change_population_parameters(msp_t *self, demographic_event_t *event)
{
    int ret = 0;
    population_id_t pid = event->params.population_parameters_change.population;
    double initial_size = event->params.population_parameters_change.initial_size;
    double growth_rate = event->params.population_parameters_change.growth_rate;

    if (pid == -1) {
        for (pid = 0; pid < (int) self->num_populations; pid++) {
            ret = msp_change_single_population_parameters(
                self, (size_t) pid, event->time, initial_size, growth_rate);
            if (ret != 0) {
                goto out;
            }
        }
    } else {
        ret = msp_change_single_population_parameters(
            self, (size_t) pid, event->time, initial_size, growth_rate);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

static void
msp_print_population_parameters_change(
    msp_t *MSP_UNUSED(self), demographic_event_t *event, FILE *out)
{
    fprintf(out,
        "%f\tpopulation_parameters_change: %d -> initial_size=%f, growth_rate=%f\n",
        event->time, (int) event->params.population_parameters_change.population,
        event->params.population_parameters_change.initial_size,
        event->params.population_parameters_change.growth_rate);
}

/* Adds a population parameter change event. Time and growth_rate are measured in
 * units of generations, and the initial size is an absolute value. */
int
msp_add_population_parameters_change(
    msp_t *self, double time, int population_id, double initial_size, double growth_rate)
{
    int ret = -1;
    demographic_event_t *de;
    int N = (int) self->num_populations;

    if (population_id < -1 || population_id >= N) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    if (initial_size < 0) {
        tsk_bug_assert(!gsl_isnan(initial_size));
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (gsl_isnan(initial_size) && gsl_isnan(growth_rate)) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    ret = msp_add_demographic_event(self, time, &de);
    if (ret != 0) {
        goto out;
    }
    de->params.population_parameters_change.population = population_id;
    /* Note we don't rescale the size and rates until we apply the event,
     * because we don't know which model will apply until then */
    de->params.population_parameters_change.initial_size = initial_size;
    de->params.population_parameters_change.growth_rate = growth_rate;
    de->change_state = msp_change_population_parameters;
    de->print_state = msp_print_population_parameters_change;
    ret = 0;
out:
    return ret;
}

/* Migration rate change */

static int MSP_WARN_UNUSED
msp_change_migration_matrix_entry(msp_t *self, size_t index, double rate)
{
    int ret = 0;
    size_t N = self->num_populations;

    if (index >= N * N) {
        ret = MSP_ERR_BAD_MIGRATION_MATRIX_INDEX;
        goto out;
    }
    if (index % (N + 1) == 0) {
        ret = MSP_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX;
        goto out;
    }
    self->migration_matrix[index] = rate;
out:
    return ret;
}

static int
msp_change_migration_rate(msp_t *self, demographic_event_t *event)
{
    int ret = 0;
    int index = event->params.migration_rate_change.matrix_index;
    int N = (int) self->num_populations;
    double rate = event->params.migration_rate_change.migration_rate;

    if (index == -1) {
        for (index = 0; index < N * N; index++) {
            if (index % (N + 1) != 0) {
                ret = msp_change_migration_matrix_entry(self, (size_t) index, rate);
                if (ret != 0) {
                    goto out;
                }
            }
        }
    } else {
        ret = msp_change_migration_matrix_entry(self, (size_t) index, rate);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

static void
msp_print_migration_rate_change(
    msp_t *MSP_UNUSED(self), demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tmigration_rate_change: %d -> %f\n", event->time,
        event->params.migration_rate_change.matrix_index,
        event->params.migration_rate_change.migration_rate);
}

/* Add a migration rate change event. Time and migration rate are measured in
 * units of generations. */
int MSP_WARN_UNUSED
msp_add_migration_rate_change(
    msp_t *self, double time, int source, int dest, double migration_rate)
{
    int ret = -1;
    demographic_event_t *de;
    int N = (int) self->num_populations;
    int matrix_index;

    if (source == -1 && dest == -1) {
        matrix_index = -1;
    } else {
        if (source < 0 || source >= N || dest < 0 || dest >= N) {
            ret = MSP_ERR_BAD_MIGRATION_MATRIX_INDEX;
            goto out;
        }
        if (source == dest) {
            ret = MSP_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX;
            goto out;
        }
        matrix_index = source * N + dest;
    }
    if (migration_rate < 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    ret = msp_add_demographic_event(self, time, &de);
    if (ret != 0) {
        goto out;
    }
    de->params.migration_rate_change.migration_rate = migration_rate;
    de->params.migration_rate_change.matrix_index = matrix_index;
    /* Wait until the event happens to rescale the rate */
    de->change_state = msp_change_migration_rate;
    de->print_state = msp_print_migration_rate_change;
    ret = 0;
out:
    return ret;
}

/* Mass migration */

static int
msp_mass_migration(msp_t *self, demographic_event_t *event)
{
    int ret = 0;
    population_id_t source = event->params.mass_migration.source;
    population_id_t dest = event->params.mass_migration.destination;
    double p = event->params.mass_migration.proportion;
    population_id_t N = (population_id_t) self->num_populations;
    avl_node_t *node, *next;
    avl_tree_t *pop;
    label_id_t label = 0; /* For now only support label 0 */

    /* This should have been caught on adding the event */
    if (source < 0 || source > N || dest < 0 || dest > N) {
        ret = MSP_ERR_ASSERTION_FAILED;
        goto out;
    }
    /*
     * Move lineages from source to dest with probability p.
     */
    pop = &self->populations[source].ancestors[label];
    node = pop->head;
    while (node != NULL) {
        next = node->next;
        if (gsl_rng_uniform(self->rng) < p) {
            ret = msp_move_individual(self, node, pop, dest, label);
            if (ret != 0) {
                goto out;
            }
        }
        node = next;
    }
out:
    return ret;
}

static void
msp_print_mass_migration(msp_t *MSP_UNUSED(self), demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tmass_migration: %d -> %d p = %f\n", event->time,
        (int) event->params.mass_migration.source,
        (int) event->params.mass_migration.destination,
        event->params.mass_migration.proportion);
}

/* Adds a mass migration event. Time is measured in units of generations */
int MSP_WARN_UNUSED
msp_add_mass_migration(
    msp_t *self, double time, int source, int destination, double proportion)
{
    int ret = 0;
    demographic_event_t *de;
    int N = (int) self->num_populations;

    if (source < 0 || source >= N || destination < 0 || destination >= N) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    if (source == destination) {
        ret = MSP_ERR_SOURCE_DEST_EQUAL;
        goto out;
    }
    if (proportion < 0.0 || proportion > 1.0) {
        ret = MSP_ERR_BAD_PROPORTION;
        goto out;
    }
    ret = msp_add_demographic_event(self, time, &de);
    if (ret != 0) {
        goto out;
    }
    de->params.mass_migration.source = source;
    de->params.mass_migration.destination = destination;
    de->params.mass_migration.proportion = proportion;
    de->change_state = msp_mass_migration;
    de->print_state = msp_print_mass_migration;
    ret = 0;
out:
    return ret;
}

static int
msp_deactivate_population(msp_t *self, int population_id)
{
    int ret = 0;
    population_t *pop = &self->populations[population_id];

    if (pop->state == MSP_POP_STATE_INACTIVE) {
        ret = MSP_ERR_DEACTIVATE_INACTIVE_POPULATION;
        goto out;
    }

    if (pop->state == MSP_POP_STATE_PREVIOUSLY_ACTIVE) {
        ret = MSP_ERR_DEACTIVATE_PREVIOUSLY_ACTIVE_POPULATION;
        goto out;
    }

    tsk_bug_assert(pop->state == MSP_POP_STATE_ACTIVE);
    tsk_bug_assert(msp_get_num_population_ancestors(self, population_id) == 0);

    pop->state = MSP_POP_STATE_PREVIOUSLY_ACTIVE;
    /* Set these to zero for tidyness sake */
    pop->initial_size = 0;
    pop->growth_rate = 0;
out:
    return ret;
}

static int
msp_activate_population(msp_t *self, int population_id)
{
    int ret = 0;
    population_t *pop = &self->populations[population_id];

    if (pop->state == MSP_POP_STATE_PREVIOUSLY_ACTIVE) {
        ret = MSP_ERR_POPULATION_PREVIOUSLY_ACTIVE;
        goto out;
    }
    if (pop->state == MSP_POP_STATE_ACTIVE) {
        ret = MSP_ERR_POPULATION_CURRENTLY_ACTIVE;
        goto out;
    }
    pop->state = MSP_POP_STATE_ACTIVE;
out:
    return ret;
}

/* Population split */

static int
msp_population_split(msp_t *self, demographic_event_t *event)
{
    int ret = 0;
    population_id_t *derived = event->params.population_split.derived;
    population_id_t ancestral = event->params.population_split.ancestral;
    size_t num_populations = event->params.population_split.num_derived;
    demographic_event_t mass_migration;
    population_t *pop = &self->populations[ancestral];
    size_t j, k;
    size_t N = self->num_populations;

    if (pop->state != MSP_POP_STATE_ACTIVE) {
        ret = msp_activate_population(self, ancestral);
        if (ret != 0) {
            goto out;
        }
    }

    /* The mass migration moves lineages into the ancestral population */
    mass_migration.params.mass_migration.destination = ancestral;
    mass_migration.params.mass_migration.proportion = 1.0;

    for (j = 0; j < num_populations; j++) {
        pop = &self->populations[derived[j]];
        if (pop->state != MSP_POP_STATE_ACTIVE) {
            ret = MSP_ERR_SPLIT_DERIVED_NOT_ACTIVE;
            goto out;
        }

        /* Turn off all migration to and from derived[j] */
        for (k = 0; k < self->num_populations; k++) {
            self->migration_matrix[((size_t) derived[j] * N) + k] = 0;
            self->migration_matrix[k * N + (size_t) derived[j]] = 0;
        }
        /* Move all lineages out of derived and into ancestral */
        mass_migration.params.mass_migration.source = derived[j];
        ret = msp_mass_migration(self, &mass_migration);
        if (ret != 0) {
            goto out;
        }
        msp_deactivate_population(self, derived[j]);
    }
out:
    return ret;
}

static void
msp_print_population_split(
    msp_t *MSP_UNUSED(self), demographic_event_t *event, FILE *out)
{
    size_t num_populations = event->params.population_split.num_derived;
    size_t j;

    fprintf(out, "%f\tpopulation_split: %d [", event->time, (int) num_populations);
    for (j = 0; j < num_populations; j++) {
        fprintf(out, "%d", event->params.population_split.derived[j]);
        if (j < num_populations - 1) {
            fprintf(out, ", ");
        }
    }
    fprintf(out, "] -> %d \n", event->params.population_split.ancestral);
}

static int MSP_WARN_UNUSED
msp_check_event_populations(
    msp_t *self, size_t num_populations, int32_t *populations, int different_population)
{
    int ret = 0;
    int N = (int) self->num_populations;
    size_t j;
    bool *population_used = calloc(self->num_populations, sizeof(*population_used));

    if (population_used == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    if (num_populations >= MSP_MAX_EVENT_POPULATIONS) {
        ret = MSP_ERR_TOO_MANY_EVENT_POPULATIONS;
        goto out;
    }
    if (different_population < 0 || different_population >= N) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    for (j = 0; j < num_populations; j++) {
        if (populations[j] < 0 || populations[j] >= N) {
            ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
            goto out;
        }
        if (populations[j] == different_population) {
            ret = MSP_ERR_SOURCE_DEST_EQUAL;
            goto out;
        }
        if (population_used[populations[j]]) {
            ret = MSP_ERR_DUPLICATE_POPULATION;
            goto out;
        }
        population_used[populations[j]] = true;
    }
out:
    msp_safe_free(population_used);
    return ret;
}

/* Adds a population split event.
 * Note: we use the int32_t type here so we can be sure we're passing in
 * arrays of the correct type from numpy where we have to specify the size.
 */
int MSP_WARN_UNUSED
msp_add_population_split(
    msp_t *self, double time, size_t num_derived, int32_t *derived, int ancestral)
{
    int ret = 0;
    size_t j;
    demographic_event_t *de;

    ret = msp_check_event_populations(self, num_derived, derived, ancestral);
    if (ret != 0) {
        goto out;
    }

    ret = msp_add_demographic_event(self, time, &de);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < num_derived; j++) {
        de->params.population_split.derived[j] = derived[j];
    }
    de->params.population_split.ancestral = ancestral;
    de->params.population_split.num_derived = num_derived;
    de->change_state = msp_population_split;
    de->print_state = msp_print_population_split;
    ret = 0;
out:
    return ret;
}

/* Admixture */

static int
msp_admixture(msp_t *self, demographic_event_t *event)
{
    int ret = 0;
    population_id_t *ancestral = event->params.admixture.ancestral;
    double *proportion = event->params.admixture.proportion;
    population_id_t derived = event->params.admixture.derived;
    population_t *pop;
    size_t num_ancestral = event->params.admixture.num_ancestral;
    size_t index;
    avl_tree_t *source;
    avl_node_t *node, *next;
    double u;
    label_id_t label = 0; /* For now only support label 0 */

    pop = &self->populations[derived];
    if (pop->state != MSP_POP_STATE_ACTIVE) {
        ret = MSP_ERR_ADMIX_DERIVED_NOT_ACTIVE;
        goto out;
    }
    for (index = 0; index < num_ancestral; index++) {
        pop = &self->populations[ancestral[index]];
        if (pop->state != MSP_POP_STATE_ACTIVE) {
            ret = msp_activate_population(self, ancestral[index]);
            if (ret != 0) {
                goto out;
            }
        }
    }

    /*
     * Move lineages from derived to ancestral[j] with probability
     * proportion[j].
     */
    source = &self->populations[derived].ancestors[label];
    node = source->head;
    while (node != NULL) {
        next = node->next;
        u = gsl_rng_uniform(self->rng);
        index = probability_list_select(u, num_ancestral, proportion);
        ret = msp_move_individual(self, node, source, ancestral[index], label);
        if (ret != 0) {
            goto out;
        }
        node = next;
    }
    msp_deactivate_population(self, derived);
out:
    return ret;
}

static void
msp_print_admixture(msp_t *MSP_UNUSED(self), demographic_event_t *event, FILE *out)
{
    size_t num_populations = event->params.admixture.num_ancestral;
    size_t j;

    fprintf(out, "%f\tadmixture: %d -> [", event->time, event->params.admixture.derived);
    for (j = 0; j < num_populations; j++) {
        fprintf(out, "(%d, p=%f)", event->params.admixture.ancestral[j],
            event->params.admixture.proportion[j]);
        if (j < num_populations - 1) {
            fprintf(out, ", ");
        }
    }
    fprintf(out, "]\n");
}

/* Adds an admixture event.
 * Note: we use the int32_t type here so we can be sure we're passing in
 * arrays of the correct type from numpy where we have to specify the size.
 */
int MSP_WARN_UNUSED
msp_add_admixture(msp_t *self, double time, int derived, size_t num_ancestral,
    int32_t *ancestral, double *proportion)
{
    int ret = 0;
    size_t j;
    demographic_event_t *de;

    ret = msp_check_event_populations(self, num_ancestral, ancestral, derived);
    if (ret != 0) {
        goto out;
    }

    for (j = 0; j < num_ancestral; j++) {
        /* We don't bother checking if the proportions sum to 1 as this
         * is better done at higher levels, and the code is fairly robust
         * to this type of misspecification. */
        if (proportion[j] < 0 || proportion[j] > 1.0) {
            ret = MSP_ERR_BAD_PROPORTION;
            goto out;
        }
    }

    ret = msp_add_demographic_event(self, time, &de);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < num_ancestral; j++) {
        de->params.admixture.ancestral[j] = ancestral[j];
        de->params.admixture.proportion[j] = proportion[j];
    }
    de->params.admixture.derived = derived;
    de->params.admixture.num_ancestral = num_ancestral;
    de->change_state = msp_admixture;
    de->print_state = msp_print_admixture;
out:
    return ret;
}

/* Simple bottlenecks. At some time we coalesce a fraction of the
 * extant lineages into a single ancestor. */

static int
msp_simple_bottleneck(msp_t *self, demographic_event_t *event)
{
    int ret = 0;
    population_id_t population_id = event->params.simple_bottleneck.population;
    double p = event->params.simple_bottleneck.proportion;
    population_id_t N = (population_id_t) self->num_populations;
    avl_node_t *node, *next, *q_node;
    avl_tree_t *pop, Q;
    segment_t *u;
    label_id_t label = 0; /* For now only support label 0 */

    /* This should have been caught on adding the event */
    if (population_id < 0 || population_id > N) {
        ret = MSP_ERR_ASSERTION_FAILED;
        goto out;
    }
    if (self->model.type == MSP_MODEL_DTWF) {
        ret = MSP_ERR_DTWF_UNSUPPORTED_BOTTLENECK;
        goto out;
    }
    avl_init_tree(&Q, cmp_segment_queue, NULL);
    /*
     * Find the individuals that descend from the common ancestor
     * during this simple_bottleneck.
     */
    pop = &self->populations[population_id].ancestors[label];
    node = pop->head;
    while (node != NULL) {
        next = node->next;
        if (gsl_rng_uniform(self->rng) < p) {
            u = (segment_t *) node->item;
            avl_unlink_node(pop, node);
            msp_free_avl_node(self, node);
            q_node = msp_alloc_avl_node(self);
            if (q_node == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            avl_init_node(q_node, u);
            q_node = avl_insert_node(&Q, q_node);
            tsk_bug_assert(q_node != NULL);
        }
        node = next;
    }
    ret = msp_merge_ancestors(self, &Q, population_id, label, TSK_NULL, NULL);
out:
    return ret;
}

static void
msp_print_simple_bottleneck(
    msp_t *MSP_UNUSED(self), demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tsimple_bottleneck: %d I = %f\n", event->time,
        (int) event->params.simple_bottleneck.population,
        event->params.simple_bottleneck.proportion);
}

/* Add a simple bottleneck event. Time is measured in generations */
int MSP_WARN_UNUSED
msp_add_simple_bottleneck(msp_t *self, double time, int population_id, double proportion)
{
    int ret = 0;
    demographic_event_t *de;
    int N = (int) self->num_populations;

    if (population_id < 0 || population_id >= N) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    if (proportion < 0.0 || proportion > 1.0) {
        ret = MSP_ERR_BAD_PROPORTION;
        goto out;
    }
    ret = msp_add_demographic_event(self, time, &de);
    if (ret != 0) {
        goto out;
    }
    de->params.simple_bottleneck.population = population_id;
    de->params.simple_bottleneck.proportion = proportion;
    de->change_state = msp_simple_bottleneck;
    de->print_state = msp_print_simple_bottleneck;
    ret = 0;
out:
    return ret;
}

/* Instantaneous bottlenecks. At a time T1 we have a burst of coalescence
 * equivalent to what would happen in time T2.
 */

static int
msp_instantaneous_bottleneck(msp_t *self, demographic_event_t *event)
{
    int ret = 0;
    population_id_t population_id = event->params.instantaneous_bottleneck.population;
    double T2 = event->params.instantaneous_bottleneck.strength;
    population_id_t N = (population_id_t) self->num_populations;
    tsk_id_t *lineages = NULL;
    tsk_id_t *pi = NULL;
    avl_node_t **avl_nodes = NULL;
    avl_tree_t *sets = NULL;
    tsk_id_t u, parent;
    uint32_t j, k, n, num_roots;
    double rate, t;
    avl_tree_t *pop;
    avl_node_t *node, *set_node;
    segment_t *individual;
    label_id_t label = 0; /* For now only support label 0 */

    if (self->model.type == MSP_MODEL_DTWF) {
        ret = MSP_ERR_DTWF_UNSUPPORTED_BOTTLENECK;
        goto out;
    }
    tsk_bug_assert(population_id >= 0 && population_id < N);
    pop = &self->populations[population_id].ancestors[label];
    n = avl_count(pop);
    lineages = malloc(n * sizeof(tsk_id_t));
    avl_nodes = malloc(n * sizeof(avl_node_t *));
    pi = malloc(2 * n * sizeof(tsk_id_t));
    sets = malloc(2 * n * sizeof(avl_tree_t));
    if (lineages == NULL || avl_nodes == NULL || pi == NULL || sets == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (u = 0; u < (tsk_id_t) n; u++) {
        lineages[u] = u;
    }
    for (u = 0; u < (tsk_id_t)(2 * n); u++) {
        pi[u] = TSK_NULL;
    }
    j = 0;
    for (node = pop->head; node != NULL; node = node->next) {
        avl_nodes[j] = node;
        j++;
    }

    /* Now we implement the Kingman coalescent for these lineages until we have
     * exceeded T2. This is based on the algorithm from Hudson 1990.
     */
    j = n - 1;
    t = 0.0;
    parent = (tsk_id_t) n;
    while (j > 0) {
        rate = (j + 1) * j / 2;
        t += msp_get_common_ancestor_waiting_time_from_rate(
            self, &self->populations[population_id], rate);
        if (t >= T2) {
            break;
        }
        /* Note: there might be issues here if we have very large sample
         * sizes as the uniform_int has a limited range.
         */
        k = (uint32_t) gsl_rng_uniform_int(self->rng, j);
        pi[lineages[k]] = parent;
        lineages[k] = lineages[j];
        j--;
        k = j > 0 ? (uint32_t) gsl_rng_uniform_int(self->rng, j) : 0;
        pi[lineages[k]] = parent;
        lineages[k] = parent;
        parent++;
    }
    num_roots = j + 1;
    for (j = 0; j < num_roots; j++) {
        if (lineages[j] >= (tsk_id_t) n) {
            avl_init_tree(&sets[lineages[j]], cmp_segment_queue, NULL);
        }
    }

    /* Assign each lineage to the set corresponding to a given root.
     * For any root < n, this lineages has not been affected, so we
     * leave it alone.
     */
    for (j = 0; j < n; j++) {
        u = (tsk_id_t) j;
        while (pi[u] != TSK_NULL) {
            u = pi[u];
        }
        if (u >= (tsk_id_t) n) {
            /* Remove this node from the population, and add it into the
             * set for the root at u */
            individual = (segment_t *) avl_nodes[j]->item;
            avl_unlink_node(pop, avl_nodes[j]);
            msp_free_avl_node(self, avl_nodes[j]);
            set_node = msp_alloc_avl_node(self);
            if (set_node == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            avl_init_node(set_node, individual);
            set_node = avl_insert_node(&sets[u], set_node);
            tsk_bug_assert(set_node != NULL);
        }
    }
    for (j = 0; j < num_roots; j++) {
        if (lineages[j] >= (tsk_id_t) n) {
            ret = msp_merge_ancestors(
                self, &sets[lineages[j]], population_id, label, TSK_NULL, NULL);
            if (ret != 0) {
                goto out;
            }
        }
    }
out:
    if (lineages != NULL) {
        free(lineages);
    }
    if (pi != NULL) {
        free(pi);
    }
    if (sets != NULL) {
        free(sets);
    }
    if (avl_nodes != NULL) {
        free(avl_nodes);
    }
    return ret;
}

static void
msp_print_instantaneous_bottleneck(
    msp_t *MSP_UNUSED(self), demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tinstantaneous_bottleneck: %d T2 = %f\n", event->time,
        (int) event->params.instantaneous_bottleneck.population,
        event->params.instantaneous_bottleneck.strength);
}

/* Add an instantaneous bottleneck event. Time and strength are measured in generations
 */
int MSP_WARN_UNUSED
msp_add_instantaneous_bottleneck(
    msp_t *self, double time, int population_id, double strength)
{
    int ret = 0;
    demographic_event_t *de;
    int N = (int) self->num_populations;

    if (population_id < 0 || population_id >= N) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    if (strength < 0.0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    ret = msp_add_demographic_event(self, time, &de);
    if (ret != 0) {
        goto out;
    }
    de->params.instantaneous_bottleneck.population = population_id;
    de->params.instantaneous_bottleneck.strength = strength;
    de->change_state = msp_instantaneous_bottleneck;
    de->print_state = msp_print_instantaneous_bottleneck;
    ret = 0;
out:
    return ret;
}

/* Add a census event at a specified time (given in generations). */
static int
msp_census_event(msp_t *self, demographic_event_t *event)
{
    int ret = 0;
    avl_tree_t *ancestors;
    avl_node_t *node;
    segment_t *seg;
    tsk_id_t i, j;
    tsk_id_t u;

    for (i = 0; i < (int) self->num_populations; i++) {
        for (j = 0; j < (int) self->num_labels; j++) {

            // Get segment from an ancestor in a population.
            ancestors = &self->populations[i].ancestors[j];
            node = ancestors->head;

            while (node != NULL) {
                seg = (segment_t *) node->item;

                while (seg != NULL) {
                    // Add an edge to the edge table.
                    ret = msp_flush_edges(self);
                    if (ret != 0) {
                        goto out;
                    }
                    ret = tsk_node_table_add_row(&self->tables->nodes,
                        MSP_NODE_IS_CEN_EVENT, event->time, i, TSK_NULL, NULL, 0);
                    if (ret < 0) {
                        goto out;
                    }
                    u = (tsk_id_t) ret;
                    // Add an edge joining the segment to the new node.
                    ret = msp_store_edge(self, seg->left, seg->right, u, seg->value);
                    if (ret != 0) {
                        goto out;
                    }
                    // Modify segment node id.
                    seg->value = u;
                    seg = seg->next;
                }
                node = node->next;
            }
        }
    }
    ret = 0;
out:
    return ret;
}

static void
msp_print_census_event(msp_t *MSP_UNUSED(self), demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tcensus_event:\n", event->time);
}

int MSP_WARN_UNUSED
msp_add_census_event(msp_t *self, double time)
{
    int ret = 0;
    demographic_event_t *de;

    ret = msp_add_demographic_event(self, time, &de);
    if (ret != 0) {
        goto out;
    }

    de->change_state = msp_census_event;
    de->print_state = msp_print_census_event;
    ret = 0;

out:
    return ret;
}

/* Add an activate_population_event at a given time */

static int
msp_activate_population_event(msp_t *self, demographic_event_t *event)
{
    population_id_t population = event->params.activate_population.population;

    return msp_activate_population(self, population);
}

static void
msp_print_activate_population_event(
    msp_t *MSP_UNUSED(self), demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tactivate_population_event:\n", event->time);
}

int MSP_WARN_UNUSED
msp_add_activate_population_event(msp_t *self, double time, int population_id)
{
    int ret = 0;
    demographic_event_t *de;
    int N = (int) self->num_populations;

    if (population_id < 0 || population_id >= N) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    ret = msp_add_demographic_event(self, time, &de);
    if (ret != 0) {
        goto out;
    }
    de->params.activate_population.population = population_id;
    de->change_state = msp_activate_population_event;
    de->print_state = msp_print_activate_population_event;
    ret = 0;
out:
    return ret;
}

/* Add a deactivate_population_event at a given time */

static int
msp_deactivate_population_event(msp_t *self, demographic_event_t *event)
{
    population_id_t population = event->params.deactivate_population.population;

    return msp_deactivate_population(self, population);
}

static void
msp_print_deactivate_population_event(
    msp_t *MSP_UNUSED(self), demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tdeactivate_population_event:\n", event->time);
}

int MSP_WARN_UNUSED
msp_add_deactivate_population_event(msp_t *self, double time, int population_id)
{
    int ret = 0;
    demographic_event_t *de;
    int N = (int) self->num_populations;

    if (population_id < 0 || population_id >= N) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    ret = msp_add_demographic_event(self, time, &de);
    if (ret != 0) {
        goto out;
    }
    de->params.deactivate_population.population = population_id;
    de->change_state = msp_deactivate_population_event;
    de->print_state = msp_print_deactivate_population_event;
    ret = 0;
out:
    return ret;
}

/*
 * Model specific implementations.
 *
 * Each model must define a way scaling time from generations to model time
 * and vice versa, along with providing methods for computing the waiting time
 * until the next common ancestor event, and the implementation of one of these
 * events.
 */

/**************************************************************
 * Standard coalescent, including SMC/SMC' variants
 **************************************************************/

static double
msp_std_get_common_ancestor_waiting_time(
    msp_t *self, population_id_t pop_id, label_id_t label)
{
    population_t *pop = &self->populations[pop_id];
    double n = (double) avl_count(&pop->ancestors[label]);
    double lambda = n * (n - 1.0) / 2.0;

    return msp_get_common_ancestor_waiting_time_from_rate(self, pop, lambda);
}

static int MSP_WARN_UNUSED
msp_std_common_ancestor_event(
    msp_t *self, population_id_t population_id, label_id_t label)
{
    int ret = 0;
    uint32_t j, n;
    avl_tree_t *ancestors;
    avl_node_t *x_node, *y_node, *node;
    segment_t *x, *y;

    ancestors = &self->populations[population_id].ancestors[label];
    /* Choose x and y */
    n = avl_count(ancestors);
    j = (uint32_t) gsl_rng_uniform_int(self->rng, n);
    x_node = avl_at(ancestors, j);
    tsk_bug_assert(x_node != NULL);
    x = (segment_t *) x_node->item;
    avl_unlink_node(ancestors, x_node);
    j = (uint32_t) gsl_rng_uniform_int(self->rng, n - 1);
    y_node = avl_at(ancestors, j);
    tsk_bug_assert(y_node != NULL);
    y = (segment_t *) y_node->item;
    avl_unlink_node(ancestors, y_node);

    /* For SMC and SMC' models we reject some events to get the required
     * distribution. */
    if (msp_reject_ca_event(self, x, y)) {
        self->num_rejected_ca_events++;
        /* insert x and y back into the population */
        tsk_bug_assert(x_node->item == x);
        node = avl_insert_node(ancestors, x_node);
        tsk_bug_assert(node != NULL);
        tsk_bug_assert(y_node->item == y);
        node = avl_insert_node(ancestors, y_node);
        tsk_bug_assert(node != NULL);
    } else {
        self->num_ca_events++;
        msp_free_avl_node(self, x_node);
        msp_free_avl_node(self, y_node);
        ret = msp_merge_two_ancestors(self, population_id, label, x, y, TSK_NULL, NULL);
    }
    return ret;
}

/**************************************************************
 * Dirac coalescent
 **************************************************************/

/* Given the specified rate, return the waiting time until the next common ancestor
 * event for the specified population */
static double
msp_dirac_get_common_ancestor_waiting_time_from_rate(
    msp_t *self, population_t *pop, double lambda)
{
    double ret = DBL_MAX;
    double alpha = pop->growth_rate;
    double t = self->time;
    double u, dt, z;

    if (lambda > 0.0) {
        u = gsl_ran_exponential(self->rng, 1.0 / lambda);
        if (alpha == 0.0) {
            if (self->ploidy == 1) {
                ret = pop->initial_size * pop->initial_size * u;
            } else {
                /* For ploidy > 1 we assume N/2 two-parent families, so that the rate
                 * with which 2 lineages belong to a common family is (2/N)^2 */
                ret = pop->initial_size * pop->initial_size * u / 4.0;
            }
        } else {
            dt = t - pop->start_time;
            if (self->ploidy == 1) {
                z = 1
                    + alpha * pop->initial_size * pop->initial_size * exp(-alpha * dt)
                          * u;
            } else {
                /* For ploidy > 1 we assume N/2 two-parent families, so that the rate
                 * with which 2 lineages belong to a common family is (2/N)^2 */
                z = 1
                    + alpha * pop->initial_size * pop->initial_size * exp(-alpha * dt)
                          * u / 4.0;
            }
            /* if z is <= 0 no coancestry can occur */
            if (z > 0) {
                ret = log(z) / alpha;
            }
        }
        if (u == 0) {
            ret = handle_zero_waiting_time(t);
        }
    }
    return ret;
}

static double
msp_dirac_get_common_ancestor_waiting_time(
    msp_t *self, population_id_t pop_id, label_id_t label)
{
    population_t *pop = &self->populations[pop_id];
    unsigned int n = (unsigned int) avl_count(&pop->ancestors[label]);
    double c = self->model.params.dirac_coalescent.c;
    double lambda = n * (n - 1.0) / 2.0;
    if (self->ploidy == 1) {
        lambda += c;
    } else {
        lambda += c / (2.0 * self->ploidy);
    }

    return msp_dirac_get_common_ancestor_waiting_time_from_rate(self, pop, lambda);
}

static int MSP_WARN_UNUSED
msp_dirac_common_ancestor_event(msp_t *self, population_id_t pop_id, label_id_t label)
{
    int ret = 0;
    uint32_t j, n, num_participants, num_parental_copies;
    avl_tree_t *ancestors, Q[4]; /* MSVC won't let us use num_pots here */
    avl_node_t *x_node, *y_node;
    segment_t *x, *y;
    double nC2, p;
    double psi = self->model.params.dirac_coalescent.psi;

    /* We assume haploid reproduction is single-parent, while all other ploidies
     * are two-parent */
    if (self->ploidy == 1) {
        num_parental_copies = 1;
    } else {
        num_parental_copies = 2 * self->ploidy;
    }

    ancestors = &self->populations[pop_id].ancestors[label];
    n = avl_count(ancestors);
    nC2 = gsl_sf_choose(n, 2);
    if (self->ploidy == 1) {
        p = (nC2 / (nC2 + self->model.params.dirac_coalescent.c));
    } else {
        p = (nC2 / (nC2 + self->model.params.dirac_coalescent.c / (2.0 * self->ploidy)));
    }
    if (gsl_rng_uniform(self->rng) < p) {
        /* When 2 * ploidy parental chromosomes are available, Mendelian segregation
         * results in a merger only 1 / (2 * ploidy) of the time. */
        if (self->ploidy == 1
            || gsl_rng_uniform(self->rng) < 1.0 / (2.0 * self->ploidy)) {
            /* Choose x and y */
            n = avl_count(ancestors);
            j = (uint32_t) gsl_rng_uniform_int(self->rng, n);
            x_node = avl_at(ancestors, j);
            tsk_bug_assert(x_node != NULL);
            x = (segment_t *) x_node->item;
            avl_unlink_node(ancestors, x_node);
            j = (uint32_t) gsl_rng_uniform_int(self->rng, n - 1);
            y_node = avl_at(ancestors, j);
            tsk_bug_assert(y_node != NULL);
            y = (segment_t *) y_node->item;
            avl_unlink_node(ancestors, y_node);
            self->num_ca_events++;
            msp_free_avl_node(self, x_node);
            msp_free_avl_node(self, y_node);
            ret = msp_merge_two_ancestors(self, pop_id, label, x, y, TSK_NULL, NULL);
        }
    } else {
        for (j = 0; j < num_parental_copies; j++) {
            avl_init_tree(&Q[j], cmp_segment_queue, NULL);
        }
        num_participants = gsl_ran_binomial(self->rng, psi, n);
        ret = msp_multi_merger_common_ancestor_event(
            self, ancestors, Q, num_participants, num_parental_copies);
        if (ret < 0) {
            goto out;
        }
        /* All the lineages that have been assigned to the particular pots can now be
         * merged.
         */
        for (j = 0; j < num_parental_copies; j++) {
            ret = msp_merge_ancestors(self, &Q[j], pop_id, label, TSK_NULL, NULL);
            if (ret < 0) {
                goto out;
            }
        }
    }
out:
    return ret;
}

/**************************************************************
 * Beta coalescent
 **************************************************************/

static double
beta_compute_juvenile_mean(msp_t *self)
{
    double alpha = self->model.params.beta_coalescent.alpha;
    double m;
    if (self->ploidy > 1) {
        m = 2 + exp(alpha * log(2) + (1 - alpha) * log(3) - log(alpha - 1));
    } else {
        m = 1 + exp((1 - alpha) * log(2) - log(alpha - 1));
    }
    return m;
}

static double
beta_compute_truncation(msp_t *self)
{
    double truncation_point = self->model.params.beta_coalescent.truncation_point;
    double m;
    if (truncation_point < DBL_MAX) {
        m = beta_compute_juvenile_mean(self);
        truncation_point /= (truncation_point + m);
    } else {
        truncation_point = 1.0;
    }
    return truncation_point;
}

static double
beta_compute_timescale(msp_t *self, population_t *pop)
{
    double alpha = self->model.params.beta_coalescent.alpha;
    double truncation_point = beta_compute_truncation(self);
    double m = beta_compute_juvenile_mean(self);
    double pop_size = pop->initial_size;
    double timescale;
    /* For ploidy > 1 we assume N/2 two-parent families, so that the rate
     * with which 2 lineages belong to a common family is based on "population size"
     * N/2 */
    if (self->ploidy > 1) {
        pop_size /= 2.0;
    }
    timescale = exp(alpha * log(m) + (alpha - 1) * log(pop_size) - log(alpha)
                    - gsl_sf_lnbeta(2 - alpha, alpha)
                    - log(gsl_sf_beta_inc(2 - alpha, alpha, truncation_point)));
    return timescale;
}

/* Given the specified rate, return the waiting time until the next common ancestor
 * event for the specified population */
static double
msp_beta_get_common_ancestor_waiting_time_from_rate(
    msp_t *self, population_t *pop, double lambda)
{
    double ret = DBL_MAX;
    double alpha = self->model.params.beta_coalescent.alpha;
    double gamma = pop->growth_rate * (alpha - 1);
    double t = self->time;
    double u, dt, z;

    if (lambda > 0.0) {
        u = gsl_ran_exponential(self->rng, 1.0 / lambda);
        if (gamma == 0.0) {
            ret = beta_compute_timescale(self, pop) * u;
        } else {
            dt = t - pop->start_time;
            z = 1 + gamma * beta_compute_timescale(self, pop) * exp(-gamma * dt) * u;
            /* if z is <= 0 no coancestry can occur */
            if (z > 0) {
                ret = log(z) / gamma;
            }
        }
        if (u == 0) {
            ret = handle_zero_waiting_time(t);
        }
    }
    return ret;
}

static double
msp_beta_get_common_ancestor_waiting_time(
    msp_t *self, population_id_t pop_id, label_id_t label)
{
    population_t *pop = &self->populations[pop_id];
    unsigned int n = (unsigned int) avl_count(&pop->ancestors[label]);
    double lambda = n * (n - 1.0) / 2.0;
    double result
        = msp_beta_get_common_ancestor_waiting_time_from_rate(self, pop, lambda);
    return result;
}

int MSP_WARN_UNUSED
msp_multi_merger_common_ancestor_event(
    msp_t *self, avl_tree_t *ancestors, avl_tree_t *Q, uint32_t k, uint32_t num_pots)
{
    int ret = 0;
    uint32_t j, i, l;
    avl_node_t *node, *q_node;
    segment_t *u;
    uint32_t pot_size;
    uint32_t cumul_pot_size = 0;

    /* In the multiple merger regime we have four different 'pots' that
     * lineages get assigned to, where all lineages in a given pot are merged into
     * a common ancestor.
     */
    for (i = 0; i < num_pots; i++) {
        pot_size = gsl_ran_binomial(self->rng, 1.0 / (num_pots - i), k - cumul_pot_size);
        cumul_pot_size += pot_size;
        if (pot_size > 1) {
            for (l = 0; l < pot_size; l++) {
                j = (uint32_t) gsl_rng_uniform_int(self->rng, avl_count(ancestors));
                node = avl_at(ancestors, j);
                tsk_bug_assert(node != NULL);

                u = (segment_t *) node->item;
                avl_unlink_node(ancestors, node);
                msp_free_avl_node(self, node);

                q_node = msp_alloc_avl_node(self);
                if (q_node == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
                avl_init_node(q_node, u);
                q_node = avl_insert_node(&Q[i], q_node);
                tsk_bug_assert(q_node != NULL);
            }
        }
    }

out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_beta_common_ancestor_event(msp_t *self, population_id_t pop_id, label_id_t label)
{
    int ret = 0;
    uint32_t j, n, num_participants, num_parental_copies;
    avl_tree_t *ancestors, Q[4]; /* MSVC won't let us use num_pots here */
    double alpha = self->model.params.beta_coalescent.alpha;
    double truncation_point = beta_compute_truncation(self);
    double beta_x, u, increment;

    /* We assume haploid reproduction is single-parent, while all other ploidies
     * are two-parent */
    if (self->ploidy == 1) {
        num_parental_copies = 1;
    } else {
        num_parental_copies = 2 * self->ploidy;
    }

    for (j = 0; j < num_parental_copies; j++) {
        avl_init_tree(&Q[j], cmp_segment_queue, NULL);
    }
    ancestors = &self->populations[pop_id].ancestors[label];
    n = avl_count(ancestors);
    beta_x = ran_inc_beta(self->rng, 2.0 - alpha, alpha, truncation_point);

    /* We calculate the probability of accepting the event */
    if (beta_x > 1e-9) {
        u = (n - 1) * log(1 - beta_x) + log(1 + (n - 1) * beta_x);
        u = exp(log(1 - exp(u)) - 2 * log(beta_x) - gsl_sf_lnchoose(n, 2));
    } else {
        /* For very small values of beta_x we need a polynomial expansion
         * for numerical stability */
        u = 0;
        for (j = 2; j <= n; j += 2) {
            increment = (j - 1) * exp(gsl_sf_lnchoose(n, j) + (j - 2) * log(beta_x));
            if (increment / u < 1e-12) {
                /* We truncate the expansion adaptively once the increment
                 * becomes negligible. */
                break;
            }
            u += increment;
        }
        for (j = 3; j <= n; j += 2) {
            increment = (j - 1) * exp(gsl_sf_lnchoose(n, j) + (j - 2) * log(beta_x));
            if (increment / u < 1e-12) {
                break;
            }
            u -= increment;
        }
        u /= gsl_sf_choose(n, 2);
    }

    if (gsl_rng_uniform(self->rng) < u) {
        do {
            /* Rejection sampling for the number of participants */
            num_participants = 2 + gsl_ran_binomial(self->rng, beta_x, n - 2);
        } while (gsl_rng_uniform(self->rng) > 1 / gsl_sf_choose(num_participants, 2));

        ret = msp_multi_merger_common_ancestor_event(
            self, ancestors, Q, num_participants, num_parental_copies);
        if (ret < 0) {
            goto out;
        }

        /* All the lineages that have been assigned to the particular pots can now be
         * merged.
         */
        for (j = 0; j < num_parental_copies; j++) {
            ret = msp_merge_ancestors(self, &Q[j], pop_id, label, TSK_NULL, NULL);
            if (ret < 0) {
                goto out;
            }
        }
    }

out:
    return ret;
}

/**************************************************************
 * Allele frequency trajectory simulation for genic selection
 *
 **************************************************************/

static double
genic_selection_stochastic_forwards(double dt, double freq, double alpha, double u)
{
    /* this is scaled following Ewens chapter 5 e.g.,
     * w_11=1+s; w_12=1+s/2; w_22=1; that is h=0.5 */
    double ux = ((alpha / 2.0) * freq * (1 - freq)) / tanh((alpha / 2.0) * freq);
    int sign = u < 0.5 ? 1 : -1;
    return freq + (ux * dt) + sign * sqrt(freq * (1.0 - freq) * dt);
}

static int
genic_selection_generate_trajectory(sweep_t *self, msp_t *simulator,
    size_t *ret_num_steps, double **ret_time, double **ret_allele_frequency)
{
    int ret = 0;
    genic_selection_trajectory_t trajectory
        = self->trajectory_params.genic_selection_trajectory;
    gsl_rng *rng = simulator->rng;
    size_t max_steps = 64;
    double *time = malloc(max_steps * sizeof(*time));
    double *allele_frequency = malloc(max_steps * sizeof(*allele_frequency));
    double x, t, *tmp, pop_size, sim_time, alpha;
    size_t num_steps;

    if (time == NULL || allele_frequency == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    /* TODO Wrap this in a rejection sample loop and get the population size
     * from the simulator. We can use
     * pop_size = get_population_size(sim->populations[0], time);
     * to do this because we assume there are no demographic events
     * during a sweep */

    x = trajectory.end_frequency;
    num_steps = 0;
    t = 0;
    sim_time = simulator->time; /*time in generations*/
    time[num_steps] = t;
    allele_frequency[num_steps] = x;
    num_steps++;
    while (x > trajectory.start_frequency) {
        if (num_steps + 1 >= max_steps) {
            max_steps *= 2;
            tmp = realloc(time, max_steps * sizeof(*time));
            if (tmp == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            time = tmp;
            tmp = realloc(allele_frequency, max_steps * sizeof(*allele_frequency));
            if (tmp == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            allele_frequency = tmp;
        }
        pop_size = get_population_size(&simulator->populations[0], sim_time);
        alpha = 2 * pop_size * trajectory.s;
        x = 1.0
            - genic_selection_stochastic_forwards(
                  trajectory.dt, 1.0 - x, alpha, gsl_rng_uniform(rng));
        /* need our recored traj to stay in bounds */
        t += trajectory.dt;
        sim_time += trajectory.dt * pop_size * simulator->ploidy;
        if (x > trajectory.start_frequency) {
            allele_frequency[num_steps] = x;
            time[num_steps] = t;
            num_steps++;
        }
    }
    tsk_bug_assert(num_steps < max_steps); /* num_steps + 1 above guarantees this */
    time[num_steps] = t;
    allele_frequency[num_steps] = trajectory.start_frequency;
    num_steps++;

    *ret_num_steps = num_steps;
    *ret_time = time;
    *ret_allele_frequency = allele_frequency;
    time = NULL;
    allele_frequency = NULL;
out:
    msp_safe_free(time);
    msp_safe_free(allele_frequency);
    return ret;
}

static void
genic_selection_print_state(sweep_t *self, FILE *out)
{
    genic_selection_trajectory_t *trajectory
        = &self->trajectory_params.genic_selection_trajectory;

    fprintf(out, "\tGenic selection trajectory\n");
    fprintf(out, "\t\tstart_frequency = %f\n", trajectory->start_frequency);
    fprintf(out, "\t\tend_frequency = %f\n", trajectory->end_frequency);
    fprintf(out, "\t\ts = %f\n", trajectory->s);
    fprintf(out, "\t\tdt = %f\n", trajectory->dt);
}

/**************************************************************
 * Public API for setting simulation models.
 **************************************************************/

static int
msp_set_simulation_model(msp_t *self, int model)
{
    int ret = 0;

    if (model != MSP_MODEL_HUDSON && model != MSP_MODEL_SMC
        && model != MSP_MODEL_SMC_PRIME && model != MSP_MODEL_DIRAC
        && model != MSP_MODEL_BETA && model != MSP_MODEL_DTWF
        && model != MSP_MODEL_WF_PED && model != MSP_MODEL_SWEEP) {
        ret = MSP_ERR_BAD_MODEL;
        goto out;
    }
    if (self->model.type == MSP_MODEL_WF_PED) {
        ret = MSP_ERR_OTHER_MODELS_WITH_PED;
        goto out;
    }
    if (self->model.free != NULL) {
        self->model.free(&self->model);
    }
    self->model.type = model;
    self->get_common_ancestor_waiting_time = msp_std_get_common_ancestor_waiting_time;
    self->common_ancestor_event = msp_std_common_ancestor_event;
    if (self->state != MSP_STATE_NEW) {
        /* We only need to setup the mass indexes if we are already simulating
         * another model */
        ret = msp_setup_mass_indexes(self);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

int
msp_set_simulation_model_hudson(msp_t *self)
{
    return msp_set_simulation_model(self, MSP_MODEL_HUDSON);
}

int
msp_set_simulation_model_smc(msp_t *self)
{
    return msp_set_simulation_model(self, MSP_MODEL_SMC);
}

int
msp_set_simulation_model_smc_prime(msp_t *self)
{
    return msp_set_simulation_model(self, MSP_MODEL_SMC_PRIME);
}

int
msp_set_simulation_model_dtwf(msp_t *self)
{
    return msp_set_simulation_model(self, MSP_MODEL_DTWF);
}

int
msp_set_simulation_model_fixed_pedigree(msp_t *self)
{
    int ret;
    tsk_size_t j, k;
    tsk_treeseq_t ts;
    tsk_node_t tsk_node;
    tsk_individual_t tsk_ind;
    individual_t *ind, *parent;
    tsk_size_t num_samples = 0;
    const tsk_size_t num_individuals = self->tables->individuals.num_rows;
    const tsk_size_t num_nodes = self->tables->nodes.num_rows;
    pedigree_t *pedigree = &self->pedigree;

    ret = tsk_treeseq_init(&ts, self->tables, TSK_TS_INIT_BUILD_INDEXES);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }

    if (self->state == MSP_STATE_SIMULATING) {
        ret = MSP_ERR_OTHER_MODELS_WITH_PED;
        goto out;
    }
    if (num_individuals == 0) {
        ret = MSP_ERR_EMPTY_PEDIGREE;
        goto out;
    }

    pedigree->individuals = calloc(num_individuals, sizeof(*pedigree->individuals));
    pedigree->visit_order = calloc(num_individuals, sizeof(*pedigree->visit_order));
    if (pedigree->individuals == NULL || pedigree->visit_order == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    pedigree->num_individuals = num_individuals;
    pedigree->num_nodes = num_nodes;
    pedigree->next_individual = PEDIGREE_UNINITIALISED;

    for (j = 0; j < num_individuals; j++) {
        ret = tsk_treeseq_get_individual(&ts, (tsk_id_t) j, &tsk_ind);
        tsk_bug_assert(ret == 0);
        ind = &pedigree->individuals[j];
        ind->ploidy = self->ploidy;
        ind->id = (tsk_id_t) j;
        if (tsk_ind.parents_length != self->ploidy) {
            ret = MSP_ERR_PEDIGREE_IND_NOT_TWO_PARENTS;
            goto out;
        }
        if (tsk_ind.nodes_length != self->ploidy) {
            ret = MSP_ERR_PEDIGREE_IND_NOT_DIPLOID;
            goto out;
        }
        for (k = 0; k < self->ploidy; k++) {
            avl_init_tree(&ind->common_ancestors[k], cmp_segment_queue, NULL);
            ind->parents[k] = tsk_ind.parents[k];
        }
        /* Tskit doesn't make any guarantees about how the nodes and pedigree
         * individuals match up, so we have to do some extra checking here. */
        ind->time = DBL_MAX;
        ind->population = -1;
        for (k = 0; k < self->ploidy; k++) {
            ret = tsk_treeseq_get_node(&ts, tsk_ind.nodes[k], &tsk_node);
            assert(ret == 0);
            assert(tsk_node.individual == (tsk_id_t) j);
            if (ind->time == DBL_MAX) {
                ind->time = tsk_node.time;
            } else if (ind->time != tsk_node.time) {
                ret = MSP_ERR_PEDIGREE_IND_NODE_TIME_DISAGREE;
                goto out;
            }
            if (ind->population == -1) {
                ind->population = tsk_node.population;
            } else if (ind->population != tsk_node.population) {
                ret = MSP_ERR_PEDIGREE_IND_NODE_POPULATION_DISAGREE;
                goto out;
            }
            ind->nodes[k] = tsk_node.id;
            if (tsk_node.flags & TSK_NODE_IS_SAMPLE) {
                num_samples++;
            }
        }
        pedigree->visit_order[j] = ind;
    }
    if (num_samples == 0) {
        ret = MSP_ERR_INSUFFICIENT_SAMPLES;
        goto out;
    }
    for (j = 0; j < num_individuals; j++) {
        ind = &pedigree->individuals[j];
        for (k = 0; k < self->ploidy; k++) {
            if (ind->parents[k] != TSK_NULL) {
                parent = &pedigree->individuals[ind->parents[k]];
                if (parent->time <= ind->time) {
                    ret = MSP_ERR_PEDIGREE_TIME_TRAVEL;
                    goto out;
                }
            }
        }
    }
    qsort(pedigree->visit_order, num_individuals, sizeof(individual_t *),
        cmp_pedigree_individual_p);
    ret = msp_set_simulation_model(self, MSP_MODEL_WF_PED);
out:
    tsk_treeseq_free(&ts);
    return ret;
}

int
msp_set_simulation_model_dirac(msp_t *self, double psi, double c)
{
    int ret = 0;

    /* We assume to be in the limit where N is infinite*/
    if (psi <= 0 || psi > 1.0) {
        ret = MSP_ERR_BAD_PSI;
        goto out;
    }
    if (c < 0.0) {
        ret = MSP_ERR_BAD_C;
        goto out;
    }
    ret = msp_set_simulation_model(self, MSP_MODEL_DIRAC);
    if (ret != 0) {
        goto out;
    }
    self->model.params.dirac_coalescent.psi = psi;
    self->model.params.dirac_coalescent.c = c;
    self->get_common_ancestor_waiting_time = msp_dirac_get_common_ancestor_waiting_time;
    self->common_ancestor_event = msp_dirac_common_ancestor_event;
out:
    return ret;
}

int
msp_set_simulation_model_beta(msp_t *self, double alpha, double truncation_point)
{
    int ret = 0;

    if (alpha <= 1.0 || alpha >= 2.0) {
        ret = MSP_ERR_BAD_BETA_MODEL_ALPHA;
        goto out;
    }

    if (truncation_point <= 0.0 || (!isfinite(truncation_point))) {
        ret = MSP_ERR_BAD_TRUNCATION_POINT;
        goto out;
    }

    ret = msp_set_simulation_model(self, MSP_MODEL_BETA);
    if (ret != 0) {
        goto out;
    }

    self->model.params.beta_coalescent.alpha = alpha;
    self->model.params.beta_coalescent.truncation_point = truncation_point;
    self->get_common_ancestor_waiting_time = msp_beta_get_common_ancestor_waiting_time;
    self->common_ancestor_event = msp_beta_common_ancestor_event;
out:
    return ret;
}

int
msp_set_simulation_model_sweep_genic_selection(msp_t *self, double position,
    double start_frequency, double end_frequency, double s, double dt)
{
    int ret = 0;
    simulation_model_t *model = &self->model;
    genic_selection_trajectory_t *trajectory
        = &model->params.sweep.trajectory_params.genic_selection_trajectory;
    double L = self->sequence_length;

    /* Check the inputs to make sure they make sense */
    if (position < 0 || position >= L) {
        ret = MSP_ERR_BAD_SWEEP_POSITION;
        goto out;
    }
    if (start_frequency <= 0.0 || start_frequency >= 1.0 || end_frequency <= 0.0
        || end_frequency >= 1.0) {
        ret = MSP_ERR_BAD_ALLELE_FREQUENCY;
        goto out;
    }
    if (start_frequency >= end_frequency) {
        ret = MSP_ERR_BAD_TRAJECTORY_START_END;
        goto out;
    }
    if (dt <= 0) {
        ret = MSP_ERR_BAD_TIME_DELTA;
        goto out;
    }
    if (s <= 0) {
        ret = MSP_ERR_BAD_SWEEP_GENIC_SELECTION_S;
        goto out;
    }

    ret = msp_set_simulation_model(self, MSP_MODEL_SWEEP);
    if (ret != 0) {
        goto out;
    }
    model->params.sweep.position = position;
    model->params.sweep.generate_trajectory = genic_selection_generate_trajectory;
    model->params.sweep.print_state = genic_selection_print_state;
    trajectory->start_frequency = start_frequency;
    trajectory->end_frequency = end_frequency;
    trajectory->s = s;
    trajectory->dt = dt;
out:
    return ret;
}
