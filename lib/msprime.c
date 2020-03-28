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
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_sf.h>

#include "util.h"
#include "avl.h"
#include "object_heap.h"
#include "fenwick.h"
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
    } while(ret > x);
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

static int
cmp_individual(const void *a, const void *b) {
    const segment_t *ia = (const segment_t *) a;
    const segment_t *ib = (const segment_t *) b;
    return (ia->id > ib->id) - (ia->id < ib->id);
}

/* For pedigree individuals we sort on time and to break ties
 * we arbitrarily use the ID */
static int
cmp_pedigree_individual(const void *a, const void *b) {
    const individual_t *ia = (const individual_t *) a;
    const individual_t *ib = (const individual_t *) b;
    int ret =  (ia->time > ib->time) - (ia->time < ib->time);
    if (ret == 0)  {
        ret = (ia->id > ib->id) - (ia->id < ib->id);
    }
    return ret;
}

/* For the segment priority queue we want to sort on the left
 * coordinate and to break ties we arbitrarily use the ID */
static int
cmp_segment_queue(const void *a, const void *b) {
    const segment_t *ia = (const segment_t *) a;
    const segment_t *ib = (const segment_t *) b;
    int ret = (ia->left > ib->left) - (ia->left < ib->left);
    if (ret == 0)  {
        ret = (ia->id > ib->id) - (ia->id < ib->id);
    }
    return ret;
}

static int
cmp_node_mapping(const void *a, const void *b) {
    const node_mapping_t *ia = (const node_mapping_t *) a;
    const node_mapping_t *ib = (const node_mapping_t *) b;
    return (ia->left > ib->left) - (ia->left < ib->left);
}

static int
cmp_sampling_event(const void *a, const void *b) {
    const sampling_event_t *ia = (const sampling_event_t *) a;
    const sampling_event_t *ib = (const sampling_event_t *) b;
    return (ia->time > ib->time) - (ia->time < ib->time);
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

int
msp_set_start_time(msp_t *self, double start_time)
{
    int ret = 0;

    if (start_time < 0.0) {
        ret = MSP_ERR_BAD_START_TIME;
        goto out;
    }
    self->start_time = self->model.generations_to_model_time(&self->model, start_time);
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
msp_set_dimensions(msp_t *self, size_t num_populations, size_t num_labels)
{
    int ret = 0;
    size_t j, k;

    if (num_populations < 1 || num_populations > UINT32_MAX) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (num_labels < 1 || num_labels > UINT32_MAX) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }

    /* Free any memory, if it has been allocated */
    for (j = 0; j < self->num_populations; j++) {
        msp_safe_free(self->populations[j].ancestors);
    }
    msp_safe_free(self->populations);
    msp_safe_free(self->initial_populations);
    msp_safe_free(self->initial_migration_matrix);
    msp_safe_free(self->migration_matrix);
    msp_safe_free(self->num_migration_events);
    msp_safe_free(self->links);
    msp_safe_free(self->segment_heap);

    self->num_populations = (uint32_t) num_populations;
    self->num_labels = (uint32_t) num_labels;
    self->initial_migration_matrix = calloc(num_populations * num_populations,
            sizeof(double));
    self->migration_matrix = calloc(num_populations * num_populations,
            sizeof(double));
    self->num_migration_events = calloc(num_populations * num_populations,
            sizeof(size_t));
    self->initial_populations = calloc(num_populations, sizeof(population_t));
    self->populations = calloc(num_populations, sizeof(population_t));
    self->links = calloc(self->num_labels, sizeof(fenwick_t));
    self->segment_heap = calloc(self->num_labels, sizeof(object_heap_t));
    if (self->migration_matrix == NULL
            || self->initial_migration_matrix == NULL
            || self->num_migration_events == NULL
            || self->initial_populations == NULL
            || self->populations == NULL
            || self->links == NULL
            || self->segment_heap == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (j = 0; j < num_populations; j++) {
        self->populations[j].ancestors = malloc(self->num_labels * sizeof(avl_tree_t));
        if (self->populations[j].ancestors == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        for (k = 0; k < num_labels; k++) {
            avl_init_tree(&self->populations[j].ancestors[k], cmp_individual, NULL);
        }
        /* Set the default sizes and growth rates. */
        self->initial_populations[j].growth_rate = 0.0;
        self->initial_populations[j].initial_size = 1.0;
        self->initial_populations[j].start_time = 0.0;
    }
out:
    return ret;
}

int
msp_set_gene_conversion_rate(msp_t *self, double rate, double track_length)
{
    int ret = 0;
    double sequence_length = recomb_map_get_sequence_length(&self->recomb_map);

    /* if the rate is zero, we ignore the track length */
    if (rate > 0) {
        if (track_length < 0 || track_length > sequence_length) {
            ret = MSP_ERR_BAD_PARAM_VALUE;
            goto out;
        }
    }
    if (rate < 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    /* NB: the gene_conversion_rate gets rescaled when we call msp_initialise,
     * so that it's in model units. */
    self->gene_conversion_rate = rate;
    self->gene_conversion_track_length = track_length;
out:
    return ret;
}

int
msp_set_num_populations(msp_t *self, size_t num_populations)
{
    return msp_set_dimensions(self, num_populations, 1);
}

int
msp_set_population_configuration(msp_t *self, int population_id, double initial_size,
        double growth_rate)
{
    int ret = MSP_ERR_BAD_POPULATION_CONFIGURATION;
    simulation_model_t *model = &self->model;

    if (population_id < 0 || population_id > (int) self->num_populations) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    if (initial_size <= 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->initial_populations[population_id].initial_size =
        initial_size / model->reference_size;
    self->initial_populations[population_id].growth_rate =
        model->generation_rate_to_model_rate(model, growth_rate);
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
    simulation_model_t *model = &self->model;

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
        self->initial_migration_matrix[j] = model->generation_rate_to_model_rate(
                model, migration_matrix[j]);
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

static segment_t * MSP_WARN_UNUSED
msp_alloc_segment(msp_t *self, double left, double right,
        double left_mass, double right_mass, node_id_t value,
        population_id_t population, label_id_t label,
        segment_t *prev, segment_t *next)
{
    segment_t *seg = NULL;

    if (object_heap_empty(&self->segment_heap[label])) {
        if (object_heap_expand(&self->segment_heap[label]) != 0) {
            goto out;
        }
        if (fenwick_expand(&self->links[label], self->segment_block_size) != 0) {
            goto out;
        }
    }
    seg = (segment_t *) object_heap_alloc_object(&self->segment_heap[label]);
    if (seg == NULL) {
        goto out;
    }
    seg->prev = prev;
    seg->next = next;
    seg->left = left;
    seg->right = right;
    seg->left_mass = left_mass;
    seg->right_mass = right_mass;
    seg->value = value;
    seg->population_id = population;
    seg->label = label;
out:
    return seg;
}

/* Top level allocators and initialisation */

int
msp_alloc(msp_t *self,
        size_t num_samples, sample_t *samples,
        recomb_map_t *recomb_map, tsk_table_collection_t *tables, gsl_rng *rng) {
    int ret = -1;
    size_t j;

    memset(self, 0, sizeof(msp_t));
    if (rng == NULL || recomb_map == NULL || tables == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    ret = recomb_map_copy(&self->recomb_map, recomb_map);
    if (ret != 0) {
        goto out;
    }

    /* Use the standard coalescent with coalescent time units by default. */
    self->model.type = -1;
    ret = msp_set_simulation_model_hudson(self, 0.25);
    assert(ret == 0);
    self->rng = rng;

    self->tables = tables;
    self->sequence_length = recomb_map_get_sequence_length(&self->recomb_map);

    /* And rescale the rates */
    self->gene_conversion_rate = self->model.generation_rate_to_model_rate(
            &self->model, self->gene_conversion_rate);

    if (num_samples > 0) {
        if (num_samples < 2 || samples == NULL || self->tables->nodes.num_rows > 0) {
            ret = MSP_ERR_BAD_PARAM_VALUE;
            goto out;
        }
        self->num_samples = (uint32_t) num_samples;
        self->samples = malloc(num_samples * sizeof(sample_t));
        if (self->samples == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        for (j = 0; j < num_samples; j++) {
            self->samples[j].population_id = samples[j].population_id;
            self->samples[j].time = samples[j].time;
            if (self->samples[j].time < 0) {
                ret = MSP_ERR_BAD_PARAM_VALUE;
                goto out;
            }
        }
    } else {
        if (self->tables->nodes.num_rows == 0) {
            ret = MSP_ERR_BAD_PARAM_VALUE;
            goto out;
        }
        self->from_ts = malloc(sizeof(*self->from_ts));
        if (self->from_ts == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        ret = tsk_treeseq_init(self->from_ts, self->tables, TSK_BUILD_INDEXES);
        if (ret != 0) {
            ret = msp_set_tsk_error(ret);
            goto out;
        }
        ret = tsk_table_collection_record_num_rows(self->tables, &self->from_position);
        if (ret != 0) {
            ret = msp_set_tsk_error(ret);
            goto out;
        }
        if (self->sequence_length != self->tables->sequence_length) {
            ret = MSP_ERR_INCOMPATIBLE_FROM_TS;
            goto out;
        }
    }

    self->start_time = -1;
    /* We have one population and one label by default */
    ret = msp_set_dimensions(self, 1, 1);
    if (ret != 0) {
        goto out;
    }
    /* Set sensible defaults for the sample_config and migration matrix */
    self->initial_migration_matrix[0] = 0.0;
    /* Set the memory defaults */
    self->store_migrations = false;
    self->store_full_arg = false;
    self->avl_node_block_size = 1024;
    self->node_mapping_block_size = 1024;
    self->segment_block_size = 1024;
    /* set up the AVL trees */
    avl_init_tree(&self->breakpoints, cmp_node_mapping, NULL);
    avl_init_tree(&self->overlap_counts, cmp_node_mapping, NULL);
    /* Set up the demographic events */
    self->demographic_events_head = NULL;
    self->demographic_events_tail = NULL;
    self->next_demographic_event = NULL;
    self->state = MSP_STATE_NEW;
    /* Set up pedigree */
    self->pedigree = NULL;
out:
    return ret;
}

static int
msp_alloc_memory_blocks(msp_t *self)
{
    int ret = 0;
    uint32_t j;

    /* Allocate the memory heaps */
    ret = object_heap_init(&self->avl_node_heap, sizeof(avl_node_t),
           self->avl_node_block_size, NULL);
    if (ret != 0) {
        goto out;
    }
    ret = object_heap_init(&self->node_mapping_heap, sizeof(node_mapping_t),
           self->node_mapping_block_size, NULL);
    if (ret != 0) {
        goto out;
    }
    /* allocate the segments and Fenwick trees */
    for (j = 0; j < self->num_labels; j++) {
        ret = object_heap_init(&self->segment_heap[j], sizeof(segment_t),
               self->segment_block_size, segment_init);
        if (ret != 0) {
            goto out;
        }
        ret = fenwick_alloc(&self->links[j], self->segment_block_size);
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
        if (self->links != NULL) {
            fenwick_free(&self->links[j]);
        }
        if (self->segment_heap != NULL) {
            object_heap_free(&self->segment_heap[j]);
        }
    }
    for (j = 0; j < self->num_populations; j++) {
        msp_safe_free(self->populations[j].ancestors);
    }
    msp_safe_free(self->links);
    msp_safe_free(self->segment_heap);
    msp_safe_free(self->initial_migration_matrix);
    msp_safe_free(self->migration_matrix);
    msp_safe_free(self->num_migration_events);
    msp_safe_free(self->initial_populations);
    msp_safe_free(self->populations);
    msp_safe_free(self->samples);
    msp_safe_free(self->sampling_events);
    msp_safe_free(self->buffered_edges);
    /* free the object heaps */
    object_heap_free(&self->avl_node_heap);
    object_heap_free(&self->node_mapping_heap);
    recomb_map_free(&self->recomb_map);
    if (self->from_ts != NULL) {
        tsk_treeseq_free(self->from_ts);
        free(self->from_ts);
    }
    if (self->model.free != NULL) {
        self->model.free(&self->model);
    }
    if (self->pedigree != NULL) {
        msp_free_pedigree(self);
    }
    ret = 0;
    return ret;
}

static inline avl_node_t * MSP_WARN_UNUSED
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

    assert(u != NULL);
    assert(u->id == id);
    return u;
}

static void
msp_free_segment(msp_t *self, segment_t *seg)
{
    object_heap_free_object(&self->segment_heap[seg->label], seg);
    fenwick_set_value(&self->links[seg->label], seg->id, 0);
}

static inline avl_tree_t *
msp_get_segment_population(msp_t *self, segment_t *u)
{
    return &self->populations[u->population_id].ancestors[u->label];
}

static inline int MSP_WARN_UNUSED
msp_insert_individual(msp_t *self, segment_t *u)
{
    int ret = 0;
    avl_node_t *node;

    assert(u != NULL);
    node = msp_alloc_avl_node(self);
    if (node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    avl_init_node(node, u);
    node = avl_insert_node(msp_get_segment_population(self, u), node);
    assert(node != NULL);
out:
    return ret;
}

static inline void
msp_remove_individual(msp_t *self, segment_t *u)
{
    avl_node_t *node;
    avl_tree_t *pop = msp_get_segment_population(self, u);

    assert(u != NULL);
    node = avl_search(pop, u);
    assert(node != NULL);
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
    search.left = x;

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
    m->left = left;
    m->value = 0;
    avl_init_node(node, m);
    node = avl_insert_node(&self->breakpoints, node);
    assert(node != NULL);
out:
    return ret;
}

static void
msp_print_segment_chain(msp_t * MSP_UNUSED(self), segment_t *head, FILE *out)
{
    segment_t *s = head;

    fprintf(out, "[pop=%d,label=%d]", s->population_id, s->label);
    while (s != NULL) {
        fprintf(out, "[(%f-%f) %d] ", s->left, s->right, (int) s->value);
        s = s->next;
    }
    fprintf(out, "\n");
}

static void
msp_verify_segments(msp_t *self, bool verify_breakpoints)
{
    double left, right, l_mass, r_mass;
    double s, ss, total_mass, alt_total_mass;
    size_t j, k;
    size_t label_segments = 0;
    size_t total_avl_nodes = 0;
    avl_node_t *node;
    segment_t *u;

    for (k = 0; k < self->num_labels; k++) {
        total_mass = 0;
        alt_total_mass = 0;
        label_segments = 0;
        for (j = 0; j < self->num_populations; j++) {
            node = (&self->populations[j].ancestors[k])->head;
            while (node != NULL) {
                u = (segment_t *) node->item;
                assert(u->prev == NULL);
                left = u->left;
                while (u != NULL) {
                    label_segments++;
                    assert(u->population_id == (population_id_t) j);
                    assert(u->label == (label_id_t) k);
                    assert(u->left < u->right);
                    assert(u->right <= self->sequence_length);

                    l_mass = recomb_map_position_to_mass(&self->recomb_map, u->left);
                    r_mass = recomb_map_position_to_mass(&self->recomb_map, u->right);
                    assert(u->left_mass == l_mass);
                    assert(u->right_mass == r_mass);
                    if (u->prev != NULL) {
                        s = recomb_map_mass_between(
                                &self->recomb_map, u->prev->right, u->right);
                    } else {
                        s = recomb_map_mass_between_left_exclusive(
                                &self->recomb_map, u->left, u->right);
                    }
                    ss = fenwick_get_value(&self->links[k], u->id);
                    total_mass += ss;
                    assert(doubles_almost_equal(s, ss, 1e-6));
                    if (s == ss) {
                        /* do nothing; just to keep compiler happy - see below also */
                    }
                    if (verify_breakpoints && u->left != 0) {
                        assert(msp_has_breakpoint(self, u->left));
                    }
                    right = u->right;
                    u = u->next;
                }
                s = recomb_map_mass_between_left_exclusive(
                        &self->recomb_map, left, right);
                alt_total_mass += s;
                node = node->next;
            }
        }
        assert(doubles_almost_equal(
                    total_mass, fenwick_get_total(&self->links[k]), 1e-6));
        assert(doubles_almost_equal(total_mass, alt_total_mass, 1e-6));
        assert(label_segments == object_heap_get_num_allocated(&self->segment_heap[k]));
    }
    total_avl_nodes = msp_get_num_ancestors(self)
            + avl_count(&self->breakpoints)
            + avl_count(&self->overlap_counts);
    assert(total_avl_nodes == object_heap_get_num_allocated(
                &self->avl_node_heap));
    assert(total_avl_nodes - msp_get_num_ancestors(self)
            == object_heap_get_num_allocated(&self->node_mapping_heap));
    if (total_avl_nodes == label_segments) {
        /* do nothing - this is just to keep the compiler happy when
         * asserts are turned off.
         */
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
    overlaps->population_id = 0;
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

    assert(self->overlaps->prev == NULL);
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
    assert(pos >= 0 && pos < self->seq_length);
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
    right_seg->population_id = 0;
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
    segment_t *u;
    uint32_t j, label, count;
    overlap_counter_t counter;
    int remaining_samples = (int) (self->num_sampling_events - self->next_sampling_event);

    int ok = overlap_counter_alloc(&counter, self->sequence_length, remaining_samples);
    assert(ok == 0);

    for (label = 0; label < self->num_labels; label++) {
        for (j = 0; j < self->num_populations; j++) {
            for (node = (&self->populations[j].ancestors[label])->head;
                    node != NULL; node = node->next) {
                u = (segment_t *) node->item;
                while (u != NULL) {
                    overlap_counter_increment_interval(&counter, u->left, u->right);
                    u = u->next;
                }
            }
        }
    }
    for (node = self->overlap_counts.head; node->next != NULL; node = node->next) {
        nm = (node_mapping_t *) node->item;
        count = overlap_counter_overlaps_at(&counter, nm->left);
        assert(nm->value == count);
    }

    overlap_counter_free(&counter);
}

void
msp_verify(msp_t *self, int options)
{
    msp_verify_segments(self, options & MSP_VERIFY_BREAKPOINTS);
    msp_verify_overlaps(self);
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
    segment_t **ancestors = malloc(msp_get_num_ancestors(self)
            * sizeof(segment_t *));

    if (ancestors == NULL && msp_get_num_ancestors(self) != 0) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = msp_get_ancestors(self, ancestors);
    if (ret != 0) {
        goto out;
    }
    fprintf(out, "simulation model      = '%s'\n", msp_get_model_name(self));
    fprintf(out, "model reference_size = %f\n", self->model.reference_size);
    if (self->model.type == MSP_MODEL_BETA) {
        fprintf(out, "\tbeta coalescent parameters: alpha = %f, truncation_point = %f\n",
                self->model.params.beta_coalescent.alpha,
                self->model.params.beta_coalescent.truncation_point);
    } else if (self->model.type == MSP_MODEL_DIRAC) {
        fprintf(out, "\tdirac coalescent parameters: psi = %f, c = %f\n",
                self->model.params.dirac_coalescent.psi,
                self->model.params.dirac_coalescent.c);
    } else if (self->model.type == MSP_MODEL_SWEEP) {
        fprintf(out, "\tsweep @ locus = %f\n", self->model.params.sweep.locus);
        self->model.params.sweep.print_state(&self->model.params.sweep, out);
    }
    fprintf(out, "n = %d\n", self->num_samples);
    fprintf(out, "m = %f\n", self->sequence_length);
    fprintf(out, "gene_conversion_rate         = %f\n", self->gene_conversion_rate);
    fprintf(out, "gene_conversion_track_length = %f\n", self->gene_conversion_track_length);
    fprintf(out, "from_ts    = %p\n", (void *) self->from_ts);
    fprintf(out, "start_time = %f\n", self->start_time);
    fprintf(out, "Samples    = \n");
    for (j = 0; j < self->num_samples; j++) {
        fprintf(out, "\t%d\tpopulation=%d\ttime=%f\n", j, (int) self->samples[j].population_id,
                self->samples[j].time);
    }
    fprintf(out, "Sampling events:\n");
    for (j = 0; j < self->num_sampling_events; j++) {
        if (j == self->next_sampling_event) {
            fprintf(out, "  ***");
        }
        se = &self->sampling_events[j];
        fprintf(out, "\t");
        fprintf(out, "%d @ %f in deme %d\n", (int) se->sample, se->time, (int) se->population_id);
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
            fprintf(out, "%0.3f ",
                self->migration_matrix[j * self->num_populations + k]);
        }
        fprintf(out, "\n");
    }

    fprintf(out, "Population sizes\n");
    for (j = 0; j < self->num_labels; j++) {
        fprintf(out, "label %d\n", j);
        fprintf(out, "\trecomb_mass = %f\n",
                fenwick_get_total(&self->links[j]));
        for (k = 0; k < self->num_populations; k++) {
            fprintf(out, "\tpop_size[%d] = %d\n", k,
                avl_count(&self->populations[k].ancestors[j]));
        }
    }
    for (j = 0; j < self->num_populations; j++) {
        fprintf(out, "pop[%d]:\n", (int) j);
        fprintf(out, "\tstart_time = %f\n", self->populations[j].start_time);
        fprintf(out, "\tinitial_size = %f\n", self->populations[j].initial_size);
        fprintf(out, "\tgrowth_rate = %f\n", self->populations[j].growth_rate);
    }
    fprintf(out, "Time = %f\n", self->time);
    for (j = 0; j < msp_get_num_ancestors(self); j++) {
        fprintf(out, "\t");
        msp_print_segment_chain(self, ancestors[j], out);
    }
    fprintf(out, "Fenwick trees\n");
    for (k = 0; k < self->num_labels; k++) {
        fprintf(out, "=====\nLabel %d\n=====\n", k);
        for (j = 1; j <= (uint32_t) fenwick_get_size(&self->links[k]); j++) {
            u = msp_get_segment(self, j, (label_id_t) k);
            v = fenwick_get_value(&self->links[k], j);
            if (v != 0) {
                fprintf(out, "\t%ld\ti=%d l=%f r=%f v=%d prev=%p next=%p\n", (long) v,
                        (int) u->id, u->left, u->right,
                        (int) u->value, (void *) u->prev, (void *) u->next);
            }
        }
    }
    fprintf(out, "Breakpoints = %d\n", avl_count(&self->breakpoints));
    for (a = self->breakpoints.head; a != NULL; a = a->next) {
        nm = (node_mapping_t *) a->item;
        fprintf(out, "\t%f -> %d\n", nm->left, nm->value);
    }
    fprintf(out, "Overlap count = %d\n", avl_count(&self->overlap_counts));
    for (a = self->overlap_counts.head; a != NULL; a = a->next) {
        nm = (node_mapping_t *) a->item;
        fprintf(out, "\t%f -> %d\n", nm->left, nm->value);
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
msp_record_migration(msp_t *self, double left, double right,
        node_id_t node, population_id_t source_pop, population_id_t dest_pop)
{
    int ret = 0;
    double scaled_time = self->model.model_time_to_generations(&self->model, self->time);

    ret = tsk_migration_table_add_row(&self->tables->migrations,
            left, right, node, source_pop, dest_pop, scaled_time);
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
    size_t j, num_edges;
    tsk_edge_t edge;

    if (self->num_buffered_edges > 0) {
        ret = tsk_squash_edges(self->buffered_edges, self->num_buffered_edges, &num_edges);
        if (ret != 0) {
            ret = msp_set_tsk_error(ret);
            goto out;
        }
        for (j = 0; j < num_edges; j++) {
            edge = self->buffered_edges[j];
            ret = tsk_edge_table_add_row(&self->tables->edges,
                    edge.left, edge.right, edge.parent, edge.child);
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
    double scaled_time = self->model.model_time_to_generations(&self->model, time);

    ret = msp_flush_edges(self);
    if (ret != 0) {
        goto out;
    }
    ret = tsk_node_table_add_row(&self->tables->nodes, flags, scaled_time, population_id,
            individual, NULL, 0);
    if (ret < 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_store_edge(msp_t *self, double left, double right, node_id_t parent, node_id_t child)
{
    int ret = 0;
    tsk_edge_t *edge;
    const double *node_time = self->tables->nodes.time;

    assert(parent > child);
    assert(parent < (node_id_t) self->tables->nodes.num_rows);
    if (self->num_buffered_edges == self->max_buffered_edges - 1) {
        /* Grow the array */
        self->max_buffered_edges *= 2;
        edge = realloc(self->buffered_edges, self->max_buffered_edges * sizeof(tsk_edge_t));
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
    self->num_buffered_edges++;
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_store_arg_edges(msp_t *self, segment_t *z)
{
    int ret = 0;
    node_id_t u = (node_id_t) msp_get_num_nodes(self) - 1;
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
    double recomb_mass;

    ind = (segment_t *) node->item;
    avl_unlink_node(source, node);
    msp_free_avl_node(self, node);

    if (self->store_full_arg) {
        ret = msp_store_node(self, MSP_NODE_IS_MIG_EVENT, self->time, dest_pop, TSK_NULL);
        if (ret != 0) {
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
                ret = msp_record_migration(self, x->left, x->right, x->value,
                        x->population_id, dest_pop);
                if (ret != 0) {
                    goto out;
                }
            }
            x->population_id = dest_pop;
        }
    } else {
        /* Because we are changing to a different Fenwick tree we must allocate
         * new segments each time. */
        new_ind = NULL;
        y = NULL;
        for (x = ind; x != NULL; x = x->next) {
            y = msp_alloc_segment(self, x->left, x->right,
                    x->left_mass, x->right_mass, x->value,
                    x->population_id, dest_label, y, NULL);
            if (new_ind == NULL) {
                new_ind = y;
            } else {
                y->prev->next = y;
            }
            recomb_mass = fenwick_get_value(&self->links[x->label], x->id);
            fenwick_increment(&self->links[y->label], y->id, recomb_mass);
            msp_free_segment(self, x);
        }
    }
    ret = msp_insert_individual(self, new_ind);
out:
    return ret;
}

/*
 * Inserts a new overlap_count at the specified locus left, mapping to the
 * specified number of overlapping segments b.
 */
static int MSP_WARN_UNUSED
msp_insert_overlap_count(msp_t *self, double left, uint32_t v)
{
    int ret = 0;
    avl_node_t *node = msp_alloc_avl_node(self);
    node_mapping_t *m = msp_alloc_node_mapping(self);

    if (node == NULL || m == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    m->left = left;
    m->value = v;
    avl_init_node(node, m);
    node = avl_insert_node(&self->overlap_counts, node);
    assert(node != NULL);
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

    search.left = k;
    avl_search_closest(&self->overlap_counts, &search, &node);
    assert(node != NULL);
    nm = (node_mapping_t *) node->item;
    if (nm->left > k) {
        node = node->prev;
        assert(node != NULL);
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

    search.left = l;
    node1 = avl_search(&self->overlap_counts, &search);
    assert(node1 != NULL);
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
    } while (node2 != NULL && nm2->left <= r);
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

/* Updates the mass on the specified segment to account for the additional
 * mass incurred between the left and right masses.
 */
static void
msp_add_segment_mass_between(msp_t *self, segment_t *seg, double l_mass, double r_mass)
{
    fenwick_increment(&self->links[seg->label], seg->id, r_mass - l_mass);
}

/* Add the mass subtended by the endpoints of seg2 to that of seg1
 */
static void
msp_add_segment_mass(msp_t *self, segment_t *seg1, segment_t *seg2)
{
    double mass = seg2->right_mass - seg2->left_mass;
    fenwick_increment(&self->links[seg1->label], seg1->id, mass);
}

/* Subtract the mass subtended by the endpoints of seg2 to that of seg1
 */
static void
msp_subtract_segment_mass(msp_t *self, segment_t *seg1, segment_t *seg2)
{
    double mass = seg2->left_mass - seg2->right_mass;
    fenwick_increment(&self->links[seg1->label], seg1->id, mass);
}

/* Set the mass of the specified segment to that between the segment's right endpoint
 * and the right endpoint of the left tail segment.
 */
static void
msp_set_segment_mass(msp_t *self, segment_t *seg, segment_t *tail_seg)
{
    double mass = seg->right_mass - tail_seg->right_mass;
    fenwick_set_value(&self->links[seg->label], seg->id, mass);
}

/* Set the mass of a specified segment that is not part of a
 * chain to the mass spanned by the segment, exlusive of its
 * endpoints.
 */
static void
msp_set_single_segment_mass(msp_t *self, segment_t *seg)
{
    double mass;
    if (recomb_map_get_discrete(&self->recomb_map)) {
        /* Exclude the left endpoint because breakpoints can't happen there */
        mass = recomb_map_mass_between(&self->recomb_map, seg->left + 1, seg->right);
    } else {
        mass = seg->right_mass - seg->left_mass;
    }

    fenwick_set_value(&self->links[seg->label], seg->id, mass);
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
            x->right_mass = y->right_mass;
            x->next = y->next;
            if (y->next != NULL) {
                y->next->prev = x;
            }
            msp_add_segment_mass(self, x, y);
            msp_free_segment(self, y);
        }
        y = x;
    }
    return 0;
}

/* TODO: is floor the correct operation here or should we round? */
static double
msp_dtwf_generate_breakpoint(msp_t *self, double start)
{
    double k = recomb_map_sample_poisson(&self->recomb_map, self->rng, start);
    assert(k > start);
    return k;
}

int MSP_WARN_UNUSED
alloc_individual(individual_t *ind, size_t ploidy)
{
    int ret;
    size_t i;

    ind->id = -1;
    ind->tsk_id = TSK_NULL;
    ind->sex = -1;
    ind->time = -1;
    ind->queued = false;
    ind->merged = false;

    // Better to allocate these as a block?
    ind->parents = malloc(ploidy * sizeof(individual_t *));
    ind->segments = malloc(ploidy * sizeof(avl_tree_t));
    if (ind->parents == NULL || ind->segments == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (i = 0; i < ploidy; i++) {
        avl_init_tree(&ind->segments[i], cmp_segment_queue, NULL);
        ind->parents[i] = NULL;
    }
    ret = 0;
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_reset_individual(msp_t *self, individual_t *ind)
{
    int ret = 0;
    /* avl_node_t *node; */
    size_t i;

    ind->tsk_id = TSK_NULL;
    ind->queued = false;
    ind->merged = false;

    for (i = 0; i < self->pedigree->ploidy; i++) {
        /* TODO: We don't yet support terminating pedigree simulations before
           reaching the pedigree founders, which means all segments are moved
           back into the population pool before a reset is possible. Might need
           more here when we support early termination. */
        assert(avl_count(&ind->segments[i]) == 0);
    }
    return ret;
}

static int MSP_WARN_UNUSED
msp_reset_pedigree(msp_t *self)
{
    int ret;
    size_t i;
    /* avl_node_t *node; */
    individual_t *ind;

    ind = self->pedigree->inds;
    for (i = 0; i < self->pedigree->num_inds; i++) {
        ret = msp_reset_individual(self, ind);
        if (ret != 0) {
            goto out;
        }
        ind++;
    }
    /* TODO: We don't yet support terminating pedigree simulations before
       reaching the pedigree founders, which means no individuals will remain in
       the pedigree heap when a reset is possibile. Might need more here when we
       support early termination. */
    assert(avl_count(&self->pedigree->ind_heap) == 0);

    self->pedigree->state = MSP_PED_STATE_UNCLIMBED;

    ret = 0;
out:
    return ret;
}

int MSP_WARN_UNUSED
msp_alloc_pedigree(msp_t *self, size_t num_inds, size_t ploidy)
{
    int ret;
    size_t i, num_samples;
    individual_t *ind;

    num_samples = self->num_samples / ploidy;

    self->pedigree = malloc(sizeof(pedigree_t));
    if (self->pedigree == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->pedigree->inds = calloc(num_inds, sizeof(individual_t));
    if (self->pedigree->inds == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ind = self->pedigree->inds;
    for (i = 0; i < num_inds; i++) {
        ret = alloc_individual(ind, ploidy);
        if (ret != 0) {
            goto out;
        }
        ind++;
    }
    self->pedigree->samples = malloc(num_samples * sizeof(individual_t *));
    if (self->pedigree->samples == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    avl_init_tree(&self->pedigree->ind_heap, cmp_pedigree_individual, NULL);

    self->pedigree->num_inds = num_inds;
    self->pedigree->ploidy = ploidy;
    self->pedigree->num_samples = num_samples;
    self->pedigree->state = MSP_PED_STATE_UNCLIMBED;

    ret = 0;
out:
    return ret;
}

int MSP_WARN_UNUSED
msp_free_pedigree(msp_t *self)
{
    int ret;
    individual_t *ind = NULL;
    size_t i;

    ind = self->pedigree->inds;
    if (ind != NULL) {
        assert(self->pedigree->num_inds > 0);
        for (i = 0; i < self->pedigree->num_inds; i++) {
            msp_safe_free(ind->parents);
            msp_safe_free(ind->segments);
            ind++;
        }
    }
    msp_safe_free(self->pedigree->inds);
    msp_safe_free(self->pedigree->samples);
    msp_safe_free(self->pedigree);

    ret = 0;
    return ret;
}

int MSP_WARN_UNUSED
msp_set_pedigree(msp_t *self, size_t num_rows, tsk_id_t *inds, tsk_id_t *parents,
        double *times, tsk_flags_t *is_sample)
{
    int ret;
    size_t i, j;
    tsk_id_t parent_ix;
    tsk_flags_t sample_flag;
    size_t sample_num;
    individual_t *ind = NULL;

    assert(self->pedigree != NULL);
    if (num_rows != self->pedigree->num_inds) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }

    ind = self->pedigree->inds;
    sample_num = 0;
    for (i = 0; i < self->pedigree->num_inds; i++) {
        ind->id = inds[i];

        if (ind->id <= 0) {
            ret = MSP_ERR_BAD_PEDIGREE_ID;
            goto out;
        }

        ind->time = times[i];

        // Link individuals to parents
        for (j = 0; j < self->pedigree->ploidy; j++) {
            parent_ix = parents[i * self->pedigree->ploidy + j];
            if (parent_ix >= 0) {
                *(ind->parents + j) = self->pedigree->inds + parent_ix;
            }
        }

        // Set samples
        sample_flag = is_sample[i];
        if (sample_flag != 0) {
            assert(sample_flag == 1);
            self->pedigree->samples[sample_num] = ind;
            sample_num++;
        }
        ind++;
    }
    if (sample_num != self->pedigree->num_samples) {
        /* TODO: Should be something like MSP_ERR_PEDIGREE_BAD_NUM_SAMPLES */
        ret = MSP_ERR_BAD_PEDIGREE_ID;
        goto out;
    }

    ret = 0;
out:
    return ret;
}

void
msp_check_samples(msp_t *self)
{
    // Samples should have a single segment for each copy of their genome
    size_t i, j;
    individual_t *sample;

    for (i = 0; i < self->pedigree->num_samples; i++) {
        sample = self->pedigree->samples[i];
        for (j = 0; j < self->pedigree->ploidy; j++) {
            assert(avl_count(&sample->segments[j]) == 1);
        }
    }
}

int MSP_WARN_UNUSED
msp_pedigree_load_pop(msp_t *self)
{
    int ret;
    size_t i, sample_ix, parent_ix, ploidy;
    population_t *pop;
    individual_t *sample_ind;
    segment_t *segment;
    avl_node_t *node;
    label_id_t label = 0;

    assert(self->num_populations == 1); // Only support single pop for now
    assert(self->pedigree->ploidy > 0);

    pop = &self->populations[0];
    ploidy = self->pedigree->ploidy;
    if (avl_count(&pop->ancestors[label]) != self->pedigree->num_samples * ploidy) {
        ret = MSP_ERR_BAD_PEDIGREE_NUM_SAMPLES;
        goto out;
    }

    // Move segments from population into pedigree samples
    i = 0;
    while(avl_count(&pop->ancestors[label]) > 0) {
        node = pop->ancestors[label].head;
        sample_ix = i / ploidy;
        sample_ind = self->pedigree->samples[sample_ix];
        parent_ix = i % ploidy;
        segment = node->item;
        avl_unlink_node(&pop->ancestors[label], node);
        msp_free_avl_node(self, node);

        ret = msp_pedigree_add_individual_segment(self, sample_ind, segment, parent_ix);
        if (ret != 0) {
            goto out;
        }
        i++;
    }
    msp_check_samples(self);
    ret = 0;
out:
    return ret;
}

int MSP_WARN_UNUSED
msp_pedigree_add_individual_segment(msp_t *self, individual_t *ind,
        segment_t *segment, size_t parent_ix)
{
    int ret;
    avl_node_t *node;

    assert(ind->segments != NULL);
    assert(parent_ix < self->pedigree->ploidy);

    node = msp_alloc_avl_node(self);
    if (node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    avl_init_node(node, segment);
    node = avl_insert_node(&ind->segments[parent_ix], node);
    assert(node != NULL);

    ret = 0;
out:
    return ret;
}

int MSP_WARN_UNUSED
msp_pedigree_build_ind_queue(msp_t *self)
{
    int ret;
    size_t i;
    individual_t *ind;

    assert(self->pedigree->num_samples > 0);
    assert(self->pedigree->samples != NULL);

    for (i = 0; i < self->pedigree->num_samples; i++) {
        ind = self->pedigree->samples[i];
        ret = msp_pedigree_push_ind(self, ind);
        if (ret != 0) {
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

int MSP_WARN_UNUSED
msp_pedigree_push_ind(msp_t *self, individual_t *ind)
{
    int ret;
    avl_node_t *node;

    assert(ind->queued == false);

    node = msp_alloc_avl_node(self);
    if (node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    avl_init_node(node, ind);
    node = avl_insert_node(&self->pedigree->ind_heap, node);
    assert(node != NULL);
    ind->queued = true;

    ret = 0;
out:
    return ret;
}

int MSP_WARN_UNUSED
msp_pedigree_pop_ind(msp_t *self, individual_t **ind)
{
    int ret;
    avl_node_t *node;

    assert(avl_count(&self->pedigree->ind_heap) > 0);

    node = self->pedigree->ind_heap.head;
    assert(node != NULL);
    *ind = node->item;
    assert((*ind)->queued);
    (*ind)->queued = false;
    msp_free_avl_node(self, node);
    avl_unlink_node(&self->pedigree->ind_heap, node);

    ret = 0;
    return ret;
}

static void
msp_print_individual(msp_t *self, individual_t ind, FILE *out)
{
    size_t j;

    fprintf(out, "ID: %d, TSK_ID %u - Time: %f, Parents: [", ind.id, ind.tsk_id, ind.time);

    for (j = 0; j < self->pedigree->ploidy; j++) {
        if (ind.parents[j] != NULL) {
            fprintf(out, " %d", ind.parents[j]->id);
        } else {
            fprintf(out, " None");
        }
    }
    fprintf(out, " ]\n");
}

void
msp_print_pedigree_inds(msp_t *self, FILE *out)
{
    individual_t ind;
    size_t i;

    assert(self->pedigree != NULL);
    assert(self->pedigree->inds != NULL);
    assert(self->pedigree->num_inds > 0);

    for (i = 0; i < self->pedigree->num_inds; i++) {
        ind = self->pedigree->inds[i];
        assert(ind.id > 0);
        msp_print_individual(self, ind, out);
    }
}

static int MSP_WARN_UNUSED
msp_dtwf_recombine(msp_t *self, segment_t *x, segment_t **u, segment_t **v)
{
    int ret = 0;
    int ix;
    double k;
    double k_mass;
    segment_t *y, *z, *tail;
    segment_t s1, s2;
    segment_t *seg_tails[] = {&s1, &s2};

    k = msp_dtwf_generate_breakpoint(self, x->left);
    s1.next = NULL;
    s2.next = NULL;
    ix = (int) gsl_rng_uniform_int(self->rng, 2);
    seg_tails[ix]->next = x;
    assert(x->prev == NULL);

    while (x != NULL) {
        seg_tails[ix] = x;
        y = x->next;

        if (x->right > k) {
            // Make new segment
            k_mass = recomb_map_position_to_mass(&self->recomb_map, k);
            assert(x->left < k);
            self->num_re_events++;
            ix = (ix + 1) % 2;

            if (seg_tails[ix] == &s1 || seg_tails[ix] == &s2) {
                tail = NULL;
            } else {
                tail = seg_tails[ix];
            }
            z = msp_alloc_segment(self, k, x->right,
                    k_mass, x->right_mass, x->value,
                    x->population_id, x->label, tail, x->next);
            if (z == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            if (z->prev == NULL) {
                msp_set_single_segment_mass(self, z);
            } else {
                msp_set_segment_mass(self, z, z->prev);
            }
            assert(z->left < z->right);
            if (x->next != NULL) {
                x->next->prev = z;
            }
            seg_tails[ix]->next = z;
            seg_tails[ix] = z;
            x->next = NULL;
            x->right = k;
            x->right_mass = k_mass;
            msp_subtract_segment_mass(self, x, z);
            assert(x->left < x->right);
            x = z;
            k = msp_dtwf_generate_breakpoint(self, k);
        }
        else if (x->right <= k && y != NULL && y->left >= k) {
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
            if (y->prev == NULL) {
                msp_set_single_segment_mass(self, y);
            } else {
                msp_set_segment_mass(self, y, y->prev);
            }
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


static double
msp_init_segments_and_compute_breakpoint(msp_t *self, label_id_t label, segment_t **x_ret, segment_t **y_ret)
{
    double h, t, k;
    segment_t *x, *y;
    fenwick_t *tree = &self->links[label];

    h = gsl_rng_uniform(self->rng) * fenwick_get_total(tree);
    y = msp_get_segment(self, fenwick_find(tree, h), label);
    t = fenwick_get_cumulative_sum(tree, y->id);
    x = y->prev;

    do {
        k = recomb_map_mass_to_position(&self->recomb_map, y->right_mass - (t - h));
    } while (y->left >= k && y->prev == NULL);

    *x_ret = x;
    *y_ret = y;

    return k;
}

static int MSP_WARN_UNUSED
msp_recombination_event(msp_t *self, label_id_t label, segment_t **lhs, segment_t **rhs)
{
    int ret = 0;
    double k, k_mass;
    segment_t *x, *y, *z, *lhs_tail;

    self->num_re_events++;
    k = msp_init_segments_and_compute_breakpoint(self, label, &x, &y);
    k_mass = recomb_map_position_to_mass(&self->recomb_map, k);
    if (y->left < k) {
        z = msp_alloc_segment(self, k, y->right,
                k_mass, y->right_mass, y->value,
                y->population_id, y->label, NULL, y->next);
        if (z == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        if (y->next != NULL) {
            y->next->prev = z;
        }
        y->next = NULL;
        y->right = k;
        y->right_mass = k_mass;
        msp_subtract_segment_mass(self, y, z);
        if (msp_has_breakpoint(self, k)) {
            self->num_multiple_re_events++;
        } else {
            ret = msp_insert_breakpoint(self, k);
            if (ret != 0) {
                goto out;
            }
        }
        lhs_tail = y;
    } else {
        assert(x != NULL);
        x->next = NULL;
        y->prev = NULL;
        z = y;
        self->num_trapped_re_events++;
        lhs_tail = x;
    }
    msp_set_single_segment_mass(self, z);
    ret = msp_insert_individual(self, z);
    if (ret != 0) {
        goto out;
    }
    if (self->store_full_arg) {
        /* Store the edges for the LHS */
        ret = msp_store_node(self, MSP_NODE_IS_RE_EVENT, self->time, lhs_tail->population_id, TSK_NULL);
        if (ret != 0) {
            goto out;
        }
        ret = msp_store_arg_edges(self, lhs_tail);
        if (ret != 0) {
            goto out;
        }
        /* Store the edges for the RHS */
        ret = msp_store_node(self, MSP_NODE_IS_RE_EVENT, self->time, z->population_id, TSK_NULL);
        if (ret != 0) {
            goto out;
        }
        ret = msp_store_arg_edges(self, z);
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
        *rhs = z;
    }
out:
    return ret;
}

/* Helper function for doing GC */
static int MSP_WARN_UNUSED
msp_cut_right_break(msp_t *self, segment_t *lhs_tail, segment_t *y, segment_t *new_segment,
        double track_end)
{
    int ret = 0;
    double track_end_mass = recomb_map_position_to_mass(&self->recomb_map, track_end);
    assert(lhs_tail != NULL);
    lhs_tail->next = new_segment;
    msp_set_segment_mass(self, new_segment, lhs_tail);
    if (y->next != NULL){
        y->next->prev = new_segment;
    }
    y->next = NULL;
    y->right = track_end;
    y->right_mass = track_end_mass;
    msp_add_segment_mass_between(self, y, new_segment->right_mass, track_end_mass);
    if (!msp_has_breakpoint(self, track_end)) {
        ret = msp_insert_breakpoint(self, track_end);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

/* Processes a gene conversion event that starts within or between segments.
 */
static int MSP_WARN_UNUSED
msp_gene_conversion_within_event(msp_t *self, label_id_t label)
{
    int ret = 0;
    double h, t;
    double k, k_mass, tl, k_plus_tl_mass;
    size_t segment_id;
    segment_t *x, *y, *y2, *z, *z2, *lhs_tail;
    double recomb_mass = fenwick_get_total(&self->links[label]);

    h = gsl_rng_uniform(self->rng) * recomb_mass;
    assert(h > 0 && h <= recomb_mass);
    /* generate track length */
    tl = gsl_ran_geometric(self->rng, 1.0 / self->gene_conversion_track_length);
    assert(tl > 0);
    segment_id = fenwick_find(&self->links[label], h);
    y = msp_get_segment(self, segment_id, label);

    t = fenwick_get_cumulative_sum(&self->links[label], segment_id);
    k = recomb_map_shift_by_mass(&self->recomb_map, y->right, h - t);
    k_mass = recomb_map_position_to_mass(&self->recomb_map, k);
    k_plus_tl_mass = recomb_map_position_to_mass(&self->recomb_map, k + tl);
    assert(k >= 0 && k < self->sequence_length);
    /* Check if the gene conversion falls between segments and hence has no effect */
    if (y->left >= k + tl){
        self->num_gc_events++;
        self->num_noneffective_gc_events++;
        return 0;
    }

    self->num_gc_events++;
    x = y->prev;

    if (k + tl < y->right) {
        /* Both breaks are within the same segment */
        if (k <= y->left) {
            y->prev = NULL;
            z2 = msp_alloc_segment(self, k + tl, y->right,
                    k_plus_tl_mass, y->right_mass, y->value,
                    y->population_id, y->label, x, y->next);
            if (z2 == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            lhs_tail = x;
            ret = msp_cut_right_break(self, lhs_tail, y, z2, k + tl);
            if (ret != 0){
                goto out;
            }
            z = y;
        } else {
            z = msp_alloc_segment(self, k, k + tl,
                    k_mass, k_plus_tl_mass, y->value,
                    y->population_id, y->label, NULL, NULL);
            z2 = msp_alloc_segment(self, k + tl, y->right,
                    k_plus_tl_mass, y->right_mass, y->value,
                    y->population_id, y->label, y, y->next);
            if (z == NULL || z2 == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            if (y->next != NULL){
                y->next->prev = z2;
            }
            y->next = z2;
            y->right = k;
            y->right_mass = k_mass;
            msp_set_segment_mass(self, z2, y);
            msp_subtract_segment_mass(self, y, z);
            msp_subtract_segment_mass(self, y, z2);
            if (!msp_has_breakpoint(self, k)) {
                ret = msp_insert_breakpoint(self, k);
                if (ret != 0) {
                    goto out;
                }
            }

            if (!msp_has_breakpoint(self, k + tl)) {
                ret = msp_insert_breakpoint(self, k + tl);
                if (ret != 0) {
                    goto out;
                }
            }
            lhs_tail = y;
        }
    } else{
        /* Breaks are in separate segments */

        /* Get the segment y2 containing the end of the conversion tract*/
        y2 = y;
        while(y2 != NULL && k + tl >= y2->right){
            y2 = y2->next;
        }
        /* Process left break */
        if (k <= y->left) {
            if (x != NULL) {
                x->next = NULL;
            }
            y->prev = NULL;
            z = y;
            lhs_tail = x;
        } else {
            z = msp_alloc_segment(self, k, y->right,
                    k_mass, y->right_mass, y->value,
                    y->population_id, y->label, NULL, y->next);
            if (z == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            msp_set_single_segment_mass(self, z);
            if (y->next != NULL) {
                y->next->prev = z;
            }
            y->next = NULL;
            y->right = k;
            y->right_mass = k_mass;
            msp_subtract_segment_mass(self, y, z);
            if (!msp_has_breakpoint(self, k)) {
                ret = msp_insert_breakpoint(self, k);
                if (ret != 0) {
                    goto out;
                }
            }
            lhs_tail = y;
        }

        /* Process right break */
        if (y2 != NULL) {
            if (y2->left < k + tl) {
                z2 = msp_alloc_segment(self, k + tl, y2->right,
                        recomb_map_position_to_mass(&self->recomb_map, k + tl),
                        y2->right_mass, y2->value,
                        y2->population_id, y2->label, lhs_tail, y2->next);
                if (z2 == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
                ret = msp_cut_right_break(self, lhs_tail, y2, z2, k + tl);
                if (ret != 0){
                    goto out;
                }
                if (z2->prev == NULL) {
                    z = z2;
                }
            } else {
                lhs_tail->next = y2;
                y2->prev->next = NULL;
                y2->prev = lhs_tail;
                msp_set_segment_mass(self, y2, lhs_tail);
            }
        }
    }

    /* Update population */
    z->label = label;
    msp_set_single_segment_mass(self, z);
    ret = msp_insert_individual(self, z);
out:
    return ret;
}

/* This is an inefficient function used until we figure out how to use a
 * Fenwick tree to implement it without iterating over the full population */
/* TODO: what does 'cleft' mean? I'm not finding it very enlightening. We should
 * change this to something more meaningful here and in algorithms.py */
static double
msp_get_cleft_total(msp_t *self)
{
    double ret = 0;
    avl_node_t *node;
    avl_tree_t *population_ancestors;
    label_id_t label;
    size_t j;
    segment_t *u;
    double dist, left, right;
    const double track_length = self->gene_conversion_track_length;
    const double x = (track_length - 1) / track_length;

    for (j = 0; j < self->num_populations; j++) {
        for (label = 0; label < (label_id_t) self->num_labels; label++) {
            population_ancestors = &self->populations[j].ancestors[label];
            for (node = population_ancestors->head; node != NULL; node = node->next) {
                u = (segment_t *) node->item;
                left = u->left;
                while (u->next != NULL) {
                    u = u->next;
                }
                right = u -> right;
                dist = right - left;
                ret += 1 - pow(x, dist - 1);
            }
        }
    }
    /* printf("CLEFT_TOTAL = %f\n", ret); */
    return ret;
}

static void
msp_find_cleft_individual(msp_t *self, double rvalue, size_t *segment_id, double *dist)
{
    avl_node_t *node;
    avl_tree_t *population_ancestors;
    label_id_t label;
    size_t j;
    size_t head_id = 0;
    segment_t *u;
    double left, right, distance;
    const double track_length = self->gene_conversion_track_length;
    const double x = (track_length - 1) / track_length;

    distance = 0;
    for (j = 0; j < self->num_populations; j++) {
        for (label = 0; label < (label_id_t) self->num_labels; label++) {
            population_ancestors = &self->populations[j].ancestors[label];
            for (node = population_ancestors->head; node != NULL; node = node->next) {
                if (rvalue > 0){
                    u = (segment_t *) node->item;
                    left = u->left;
                    head_id = u->id;
                    while (u->next != NULL) {
                        u = u->next;
                    }
                    right = u -> right;
                    distance = right - left;
                    rvalue -= 1 - pow(x, distance - 1);
                }
            }
        }
    }
    *segment_id = head_id;
    *dist = distance;
}

/* Processes a gene conversion event that started left of a first
 * segment and does not span the whole segment chain*/
static int MSP_WARN_UNUSED
msp_gene_conversion_left_event(msp_t *self, label_id_t label)
{
    int ret = 0;
    double h, length, p, logp, u;
    segment_t *x, *y, *z;
    double k, tl, k_mass;
    size_t segment_id;
    const double track_length = self->gene_conversion_track_length;

    self->num_gc_events++;
    h = gsl_rng_uniform(self->rng) * msp_get_cleft_total(self);
    /* Get the segment where gc starts from left and the length of the segment chain */
    msp_find_cleft_individual(self, h, &segment_id, &length);
    y = msp_get_segment(self, segment_id, label);
    /* Generate conditional track length */
    assert(length > 0);

    tl = 1.0;
    if (track_length > 1.0) {
        /* p is the proba of continuing the track */
        p = (track_length - 1.0) / track_length;
        logp = log(1.0 - 1.0 / track_length);
        u = gsl_rng_uniform(self->rng);
        tl = floor(1.0 + log(1.0 - u * (1.0 - pow(p, length - 1.0))) / logp);
    }
    k = y->left + tl;
    k_mass = recomb_map_position_to_mass(&self->recomb_map, k);

    while (y->right <= k){
        y = y->next;
    }
    x = y->prev;
    if (y->left < k){
        /*make new segment*/
        z = msp_alloc_segment(self, k, y->right,
                k_mass, y->right_mass, y->value,
                y->population_id, y->label, NULL, y->next);
        if (z == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        if (y->next != NULL){
            y->next->prev = z;
        }
        y->next = NULL;
        y->right = k;
        y->right_mass = k_mass;
        msp_subtract_segment_mass(self, y, z);
        if (!msp_has_breakpoint(self, k)) {
            ret = msp_insert_breakpoint(self, k);
            if (ret != 0) {
                goto out;
            }
        }
    } else {
        /*split the link between x and y*/
        x->next = NULL;
        y->prev = NULL;
        z = y;
    }
    z->label = label;
    msp_set_single_segment_mass(self, z);
    ret = msp_insert_individual(self, z);
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

static void
msp_set_segment_left_endpoint(msp_t *self, segment_t *seg, double left)
{
    seg->left = left;
    seg->left_mass = recomb_map_position_to_mass(&self->recomb_map, left);
}

static int MSP_WARN_UNUSED
msp_merge_two_ancestors(msp_t *self, population_id_t population_id, label_id_t label,
        segment_t *a, segment_t *b)
{
    int ret = 0;
    bool coalescence = false;
    bool defrag_required = false;
    node_id_t v;
    double l, r, l_min, r_max;
    avl_node_t *node;
    node_mapping_t *nm, search;
    segment_t *x, *y, *z, *alpha, *beta;

    x = a;
    y = b;
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
                alpha = msp_alloc_segment(self, x->left, y->left,
                        x->left_mass, y->left_mass, x->value,
                        x->population_id, x->label, NULL, NULL);
                if (alpha == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
                x->left = y->left;
                x->left_mass = y->left_mass;
            } else {
                l = x->left;
                r_max = GSL_MIN(x->right, y->right);
                if (!coalescence) {
                    coalescence = true;
                    l_min = l;
                    ret = msp_store_node(self, 0, self->time, population_id, TSK_NULL);
                    if (ret != 0) {
                        goto out;
                    }
                }
                v = (node_id_t) msp_get_num_nodes(self) - 1;
                /* Insert overlap counts for bounds, if necessary */
                search.left = l;
                node = avl_search(&self->overlap_counts, &search);
                if (node == NULL) {
                    ret = msp_copy_overlap_count(self, l);
                    if (ret < 0) {
                        goto out;
                    }
                }
                search.left = r_max;
                node = avl_search(&self->overlap_counts, &search);
                if (node == NULL) {
                    ret = msp_copy_overlap_count(self, r_max);
                    if (ret < 0) {
                        goto out;
                    }
                }
                /* Now get overlap count at the left */
                search.left = l;
                node = avl_search(&self->overlap_counts, &search);
                assert(node != NULL);
                nm = (node_mapping_t *) node->item;
                if (nm->value == 2) {
                    nm->value = 0;
                    node = node->next;
                    assert(node != NULL);
                    nm = (node_mapping_t *) node->item;
                    r = nm->left;
                } else {
                    r = l;
                    while (nm->value != 2 && r < r_max) {
                        nm->value--;
                        node = node->next;
                        assert(node != NULL);
                        nm = (node_mapping_t *) node->item;
                        r = nm->left;
                    }
                    alpha = msp_alloc_segment(self, l, r,
                            recomb_map_position_to_mass(&self->recomb_map, l),
                            recomb_map_position_to_mass(&self->recomb_map, r),
                            v, population_id, label,
                            NULL, NULL);
                    if (alpha == NULL) {
                        ret = MSP_ERR_NO_MEMORY;
                        goto out;
                    }
                }
                assert(v != x->value);
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
                    msp_set_segment_left_endpoint(self, x, r);
                }
                if (y->right == r) {
                    beta = y;
                    y = y->next;
                    msp_free_segment(self, beta);
                } else {
                    msp_set_segment_left_endpoint(self, y, r);
                }
            }
        }
        if (alpha != NULL) {
            if (z == NULL) {
                ret = msp_insert_individual(self, alpha);
                if (ret != 0) {
                    goto out;
                }
                msp_set_single_segment_mass(self, alpha);
            } else {
                if (self->store_full_arg) {
                    // we pre-empt the fact that values will be set equal later
                    defrag_required |= z->right == alpha->left;
                } else {
                    defrag_required |= z->right == alpha->left && z->value == alpha->value;
                }
                z->next = alpha;
                msp_set_segment_mass(self, alpha, z);
            }
            alpha->prev = z;
            z = alpha;
        }
    }
    if (self->store_full_arg) {
        if (! coalescence) {
            ret = msp_store_node(self, MSP_NODE_IS_CA_EVENT, self->time, population_id, TSK_NULL);
            if (ret != 0) {
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
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_priority_queue_insert(msp_t *self, avl_tree_t *Q, segment_t *u)
{
    int ret = 0;
    avl_node_t *node;

    assert(u != NULL);
    node = msp_alloc_avl_node(self);
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
        label_id_t label, segment_t **merged_segment, tsk_id_t individual)
{
    int ret = MSP_ERR_GENERIC;
    bool coalescence = false;
    bool defrag_required = false;
    bool set_merged = false;
    node_id_t v;
    uint32_t j, h;
    double l, r, r_max, next_l, next_l_mass, l_min;
    avl_node_t *node;
    node_mapping_t *nm, search;
    segment_t *x, *z, *alpha;
    segment_t **H = NULL;

    H = malloc(avl_count(Q) * sizeof(segment_t *));
    if (H == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    r_max = 0; /* keep compiler happy */
    l_min = 0;
    z = NULL;
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
            next_l_mass = ((segment_t *) node->item)->left_mass;
            r_max = GSL_MIN(r_max, next_l);
        }
        alpha = NULL;
        if (h == 1) {
            x = H[0];
            if (node != NULL && next_l < x->right) {
                alpha = msp_alloc_segment(self, x->left, next_l,
                        x->left_mass, next_l_mass,
                        x->value,
                        x->population_id, x->label, NULL, NULL);
                if (alpha == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
                x->left = next_l;
                x->left_mass = next_l_mass;
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
            if (!coalescence) {
                coalescence = true;
                l_min = l;
                ret = msp_store_node(self, 0, self->time, population_id, individual);
                if (ret != 0) {
                    goto out;
                }
            }
            v = (node_id_t) msp_get_num_nodes(self) - 1;
            /* Insert overlap counts for bounds, if necessary */
            search.left = l;
            node = avl_search(&self->overlap_counts, &search);
            if (node == NULL) {
                ret = msp_copy_overlap_count(self, l);
                if (ret < 0) {
                    goto out;
                }
            }
            search.left = r_max;
            node = avl_search(&self->overlap_counts, &search);
            if (node == NULL) {
                ret = msp_copy_overlap_count(self, r_max);
                if (ret < 0) {
                    goto out;
                }
            }
            /* Update the extant segments and allocate alpha if the interval
             * has not coalesced. */
            search.left = l;
            node = avl_search(&self->overlap_counts, &search);
            assert(node != NULL);
            nm = (node_mapping_t *) node->item;
            if (nm->value == h) {
                nm->value = 0;
                node = node->next;
                assert(node != NULL);
                nm = (node_mapping_t *) node->item;
                r = nm->left;
            } else {
                r = l;
                while (nm->value != h && r < r_max) {
                    nm->value -= h - 1;
                    node = node->next;
                    assert(node != NULL);
                    nm = (node_mapping_t *) node->item;
                    r = nm->left;
                }
                alpha = msp_alloc_segment(self, l, r,
                        recomb_map_position_to_mass(&self->recomb_map, l),
                        recomb_map_position_to_mass(&self->recomb_map, r),
                        v, population_id,
                        label, NULL, NULL);
                if (alpha == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
            }
            /* Store the edges and update the priority queue */
            for (j = 0; j < h; j++) {
                x = H[j];
                assert(v != x->value);
                ret = msp_store_edge(self, l, r, v, x->value);
                if (ret != 0) {
                    goto out;
                }
                if (x->right == r) {
                    msp_free_segment(self, x);
                    x = x->next;
                } else if (x->right > r) {
                    msp_set_segment_left_endpoint(self, x, r);
                }
                if (x != NULL) {
                    ret = msp_priority_queue_insert(self, Q, x);
                    if (ret != 0) {
                        goto out;
                    }
                } else {
                    /* If we've fully coalesced, we stop tracking the segment in
                       the pedigree. */
                    if (self->pedigree != NULL &&
                            self->pedigree->state == MSP_PED_STATE_CLIMBING) {
                        if (!set_merged) {
                            *merged_segment = NULL;
                        }
                    }
                }
            }
        }
        /* Loop tail; integrate alpha into the global state */
        if (alpha != NULL) {
            if (z == NULL) {
                msp_set_single_segment_mass(self, alpha);
                /* Pedigree doesn't currently track lineages in Populations, so
                   keep reference to merged segments instead */
                if (self->pedigree != NULL &&
                        self->pedigree->state == MSP_PED_STATE_CLIMBING) {
                    assert(merged_segment != NULL);
                    set_merged = true; //TODO: Must be better way of checking this
                    *merged_segment = alpha;
                } else {
                    ret = msp_insert_individual(self, alpha);
                    if (ret != 0) {
                        goto out;
                    }
                }
            } else {
                if (self->store_full_arg) {
                    // we pre-empt the fact that values will be set equal later
                    defrag_required |= z->right == alpha->left;
                } else {
                    defrag_required |=
                        z->right == alpha->left && z->value == alpha->value;
                }
                z->next = alpha;
                msp_set_segment_mass(self, alpha, z);
            }
            alpha->prev = z;
            z = alpha;
        }
    }
    if (self->store_full_arg) {
        if (! coalescence) {
            ret = msp_store_node(self, MSP_NODE_IS_CA_EVENT, self->time, population_id, individual);
            if (ret != 0) {
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
    ret = 0;
out:
    if (H != NULL) {
        free(H);
    }
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
    assert(node != NULL);
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

static int MSP_WARN_UNUSED
msp_insert_sample(msp_t *self, node_id_t sample, population_id_t population)
{
    int ret = MSP_ERR_GENERIC;
    double seq_len = self->sequence_length;
    segment_t *u;

    u = msp_alloc_segment(self, 0, seq_len,
            0, recomb_map_position_to_mass(&self->recomb_map, seq_len),
            sample, population, 0, NULL, NULL);
    if (u == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = msp_insert_individual(self, u);
    if (ret != 0) {
        goto out;
    }
    msp_set_single_segment_mass(self, u);
out:
    return ret;
}

static inline int
msp_allocate_root_segments(msp_t *self, tsk_tree_t *tree,
        double left, double right,
        segment_t * restrict *root_segments_head,
        segment_t * restrict *root_segments_tail)
{
    int ret = 0;
    tsk_table_collection_t *tables = self->from_ts->tables;
    const population_id_t *node_population = tables->nodes.population;
    population_id_t population;
    node_id_t root;
    segment_t *seg, *tail;
    label_id_t label = 0; /* For now only support label 0 */

    for (root = tree->left_root; root != TSK_NULL; root = tree->right_sib[root]) {
        population = node_population[root];
        /* Reference integrity has alreay been checked, but the null population
         * is still possibile */
        if (population == TSK_NULL) {
            ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
            goto out;
        }
        if (root_segments_head[root] == NULL) {
            seg = msp_alloc_segment(self, left, right,
                    recomb_map_position_to_mass(&self->recomb_map, left),
                    recomb_map_position_to_mass(&self->recomb_map, right),
                    root, population, label,
                    NULL, NULL);
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
                tail->right_mass = recomb_map_position_to_mass(&self->recomb_map, right);
            } else {
                seg = msp_alloc_segment(self, left, right,
                        recomb_map_position_to_mass(&self->recomb_map, left),
                        recomb_map_position_to_mass(&self->recomb_map, right),
                        root, population, label,
                        tail, NULL);
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
msp_reset_from_ts(msp_t *self)
{
    int ret = 0;
    int t_iter;
    tsk_tree_t t;
    node_id_t root;
    segment_t *seg;
    uint32_t num_roots, overlap, last_overlap;
    size_t num_nodes = self->tables->nodes.num_rows;
    segment_t **root_segments_head = calloc(num_nodes, sizeof(*root_segments_head));
    segment_t **root_segments_tail = calloc(num_nodes, sizeof(*root_segments_tail));

    ret = tsk_tree_init(&t, self->from_ts, 0);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
    if (root_segments_head == NULL || root_segments_tail == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    /* Reset the tables to their correct position for replication */
    ret = tsk_table_collection_truncate(self->tables, &self->from_position);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }

    last_overlap = UINT32_MAX;
    for (t_iter = tsk_tree_first(&t); t_iter == 1; t_iter = tsk_tree_next(&t)) {
        num_roots = (uint32_t) tsk_tree_get_num_roots(&t);
        overlap = 0;
        if (num_roots > 1) {
            overlap = num_roots;
            ret = msp_allocate_root_segments(self, &t, t.left, t.right,
                    root_segments_head, root_segments_tail);
            if (ret != 0) {
                goto out;
            }
        }
        if (overlap != last_overlap) {
            ret = msp_insert_overlap_count(self, t.left, overlap);
            if (ret != 0) {
                goto out;
            }
        }
    }
    if (t_iter != 0) {
        ret = msp_set_tsk_error(t_iter);
        goto out;
    }

    ret = msp_insert_overlap_count(self, self->sequence_length, UINT32_MAX);
    if (ret != 0) {
        goto out;
    }

    /* Insert the segment chains into the algorithm state */
    for (root = 0; root < (node_id_t) num_nodes; root++) {
        seg = root_segments_head[root];
        if (seg != NULL) {
            ret = msp_insert_individual(self, seg);
            if (ret != 0) {
                goto out;
            }
            msp_set_single_segment_mass(self, seg);
            for (seg = seg->next; seg != NULL; seg = seg->next) {
                msp_set_segment_mass(self, seg, seg->prev);
            }
        }
    }

out:
    tsk_tree_free(&t);
    msp_safe_free(root_segments_head);
    msp_safe_free(root_segments_tail);
    return ret;
}

static int
msp_reset_from_samples(msp_t *self)
{
    int ret = 0;
    char id_str[100];
    tsk_size_t id_str_len;
    size_t sample_idx, j;
    individual_t *ind = NULL;
    node_id_t u;
    tsk_id_t tsk_ind;

    tsk_table_collection_clear(self->tables);

    self->tables->sequence_length = recomb_map_get_sequence_length(&self->recomb_map);
    for (j = 0; j < self->num_populations; j++) {
        ret = tsk_population_table_add_row(&self->tables->populations, NULL, 0);
        if (ret < 0) {
            ret = msp_set_tsk_error(ret);
            goto out;
        }
    }
    /* Set up the sample */
    for (u = 0; u < (node_id_t) self->num_samples; u++) {
        if (self->samples[u].time <= self->start_time) {
            ret = msp_insert_sample(self, u, self->samples[u].population_id);
            if (ret != 0) {
                goto out;
            }
        }
        tsk_ind = TSK_NULL;
        if (self->pedigree != NULL) {
            /* If we're doing pedigree simulations, assign 'ploidy' nodes
             * per individual. */
            assert(self->pedigree->num_samples * self->pedigree->ploidy ==
                    self->num_samples);
            sample_idx = (size_t) u / self->pedigree->ploidy;

            // TODO: When pedigrees and populations are properly sorted out,
            //       add population to individual here
            ind = self->pedigree->samples[sample_idx];
            id_str_len = (tsk_size_t) ceil(log10(ind->id + 1));
            sprintf(id_str, "%d", ind->id);
            if (ind->tsk_id == TSK_NULL) {
                ret = tsk_individual_table_add_row(&self->tables->individuals, 0,
                        NULL, 0, id_str, id_str_len);
                if (ret < 0) {
                    goto out;
                }
                ind->tsk_id = ret;
            }
            tsk_ind = ind->tsk_id;
        }
        ret = msp_store_node(self, TSK_NODE_IS_SAMPLE,
                self->samples[u].time,
                self->samples[u].population_id, tsk_ind);
        if (ret != 0) {
            goto out;
        }
    }
    ret = msp_insert_overlap_count(self, 0, self->num_samples);
    if (ret != 0) {
        goto out;
    }
    ret = msp_insert_overlap_count(self, self->sequence_length, self->num_samples + 1);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_apply_demographic_events(msp_t *self)
{
    int ret = 0;
    demographic_event_t *event;

    assert(self->next_demographic_event != NULL);
    /* Process all events with equal time in one block. */
    self->time = self->next_demographic_event->time;
    while (self->next_demographic_event != NULL
            && self->next_demographic_event->time == self->time) {
        /* We skip ahead to the start time for the next demographic
         * event, and use its change_state method to update the
         * state of the simulation.
         */
        event = self->next_demographic_event;
        assert(event->change_state != NULL);
        ret = event->change_state(self, event);
        if (ret != 0) {
            goto out;
        }
        self->next_demographic_event = event->next;
    }
out:
    return ret;
}

int
msp_reset(msp_t *self)
{
    // TODO: This will need to be updated to allow replicates within pedigree sims
    int ret = 0;
    size_t N = self->num_populations;
    population_id_t population_id;
    population_t *pop, *initial_pop;

    memcpy(&self->model, &self->initial_model, sizeof(self->model));
    if (self->pedigree != NULL) {
        ret = msp_reset_pedigree(self);
        if (ret != 0) {
            goto out;
        }
    }
    ret = msp_reset_memory_state(self);
    if (ret != 0) {
        goto out;
    }
    /* Set up the initial segments and algorithm state */
    self->time = self->start_time;
    assert(self->time >= 0);
    for (population_id = 0; population_id < (population_id_t) N; population_id++) {
        pop = self->populations + population_id;
        /* Set the initial population parameters */
        initial_pop = self->initial_populations + population_id;
        pop->growth_rate = initial_pop->growth_rate;
        pop->initial_size = initial_pop->initial_size;
        pop->start_time = self->time;
    }
    if (self->from_ts == NULL) {
        ret = msp_reset_from_samples(self);
    } else {
        ret = msp_reset_from_ts(self);
    }
    if (ret != 0) {
        goto out;
    }
    self->next_demographic_event = self->demographic_events_head;
    memcpy(self->migration_matrix, self->initial_migration_matrix,
            N * N * sizeof(double));
    self->next_sampling_event = 0;
    self->num_re_events = 0;
    self->num_gc_events = 0;
    self->num_ca_events = 0;
    self->num_rejected_ca_events = 0;
    self->num_trapped_re_events = 0;
    self->num_multiple_re_events = 0;
    memset(self->num_migration_events, 0, N * N * sizeof(size_t));
    self->state = MSP_STATE_INITIALISED;
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_initialise_from_samples(msp_t *self)
{
    int ret = 0;
    size_t j, k, initial_samples;
    population_id_t pop;

    if (self->start_time < 0) {
        self->start_time = 0;
    }
    initial_samples = 0;
    for (j = 0; j < self->num_samples; j++) {
        /* Check that the sample configuration makes sense */
        pop = self->samples[j].population_id;
        if (pop < 0 || pop >= (population_id_t) self->num_populations) {
            ret = MSP_ERR_BAD_SAMPLES;
            goto out;
        }
        if (self->samples[j].time <= self->start_time) {
            initial_samples++;
        }
    }
    /* Set up the historical sampling events */
    self->num_sampling_events = self->num_samples - initial_samples;
    self->sampling_events = NULL;
    if (self->num_sampling_events > 0) {
        self->sampling_events = malloc(self->num_sampling_events *
                sizeof(sampling_event_t));
        if (self->sampling_events == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        k = 0;
        for (j = 0; j < self->num_samples; j++) {
            if (self->samples[j].time > self->start_time) {
                self->sampling_events[k].sample = (node_id_t) j;
                self->sampling_events[k].time = self->samples[j].time;
                self->sampling_events[k].population_id =
                    self->samples[j].population_id;
                k++;
            }
        }
        assert(k == self->num_sampling_events);
        /* Now we must sort the sampling events by time. */
        qsort(self->sampling_events, self->num_sampling_events,
                sizeof(sampling_event_t), cmp_sampling_event);
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
msp_initialise_from_ts(msp_t *self)
{
    int ret = 0;
    double model_time, root_time;
    uint32_t num_samples;
    size_t num_nodes = self->tables->nodes.num_rows;
    population_id_t pop;
    size_t j;

    if (self->num_populations != self->tables->populations.num_rows) {
        ret = MSP_ERR_INCOMPATIBLE_FROM_TS;
        goto out;
    }
    /* Find the maximum time among the existing nodes */
    num_samples = 0;
    root_time = 0.0;
    for (j = 0; j < num_nodes; j++) {
        model_time = self->model.generations_to_model_time(
                &self->model, self->tables->nodes.time[j]);
        root_time = GSL_MAX(model_time, root_time);
        /* TODO we can catch ancient samples here and insert them as sampling
         * events, if we wish to support this. */
        if (self->tables->nodes.flags[j] & TSK_NODE_IS_SAMPLE) {
            num_samples++;
        }
        pop = self->tables->nodes.population[j];
        if (pop < 0 || pop >= (population_id_t) self->num_populations) {
            ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
            goto out;
        }
    }
    if (self->start_time < 0) {
        self->start_time = root_time;
    } else {
        if (root_time > self->start_time) {
            ret = MSP_ERR_BAD_START_TIME_FROM_TS;
            goto out;
        }
    }
    if (num_samples < 2) {
        ret = MSP_ERR_INSUFFICIENT_SAMPLES;
        goto out;
    }
out:
    return ret;
}

/*
 * Sets up the memory heaps and rescales times and rates into simulation units.
 */
int MSP_WARN_UNUSED
msp_initialise(msp_t *self)
{
    int ret = -1;

    /* These should really be proper checks with a return value */
    assert(self->num_populations >= 1);

    ret = msp_alloc_memory_blocks(self);
    if (ret != 0) {
        goto out;
    }
    if (self->from_ts == NULL) {
        ret = msp_initialise_from_samples(self);
    } else {
        ret = msp_initialise_from_ts(self);
    }
    if (ret != 0) {
        goto out;
    }

    /* Copy the state of the simulation model into the initial model */
    memcpy(&self->initial_model, &self->model, sizeof(self->model));
    /* If any demographic events have time < than the start_time then
     * raise an error */
    if (self->demographic_events_head != NULL) {
        if (self->demographic_events_head->time < self->start_time) {
            ret = MSP_ERR_BAD_DEMOGRAPHIC_EVENT_TIME;
            goto out;
        }
    }
    ret = msp_reset(self);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
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

/* Given the specified rate, return the waiting time until the next common ancestor
 * event for the specified population */
static double
msp_get_common_ancestor_waiting_time_from_rate(msp_t *self, population_t *pop, double lambda)
{
    double ret = DBL_MAX;
    double alpha = pop->growth_rate;
    double t = self->time;
    double u, dt, z;

    if (lambda > 0.0) {
        u = gsl_ran_exponential(self->rng, 1.0 / lambda);
        if (alpha == 0.0) {
            ret = pop->initial_size * u;
        } else {
            dt = t - pop->start_time;
            z = 1 + alpha * pop->initial_size * exp(-alpha * dt) * u;
            /* if z is <= 0 no coancestry can occur */
            if (z > 0) {
                ret = log(z) / alpha;
            }
        }
        if (u == 0) {
            /* In the exceedingly rare cases where gsl_ran_exponential returns
             * 0, we return the smallest representable value > the current time
             * to avoid returning a tree sequence with zero length branches.
             * Note that we can still return 0 from this function if population
             * sizes are extremely small. This is intentional, as it is almost
             * certainly an error to have simulations at such extreme values
             * and the user should be alerted to this. */
            ret = nextafter(t, DBL_MAX) - t;
            assert(ret != 0);
        }
    }
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
 * 0 if the simulation completed to coalescence
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
    double lambda, t_temp, t_wait, ca_t_wait, re_t_wait,
           gc_in_t_wait, gc_left_t_wait, mig_t_wait,
           sampling_event_time, demographic_event_time;
    double recomb_mass, total_recomb_rate;
    uint32_t j, k, n;
    population_id_t ca_pop_id, mig_source_pop, mig_dest_pop;
    unsigned long events = 0;
    sampling_event_t *se;
    /* Only support a single label for now. */
    label_id_t label = 0;

    gc_in_t_wait = DBL_MAX;
    gc_left_t_wait = DBL_MAX;

    while (msp_get_num_ancestors(self) > 0) {
        if (events == max_events) {
            ret = MSP_EXIT_MAX_EVENTS;
            break;
        }
        events++;

        recomb_mass = fenwick_get_total(&self->links[label]);
        /* Recombination */
        lambda = recomb_mass;
        re_t_wait = DBL_MAX;
        if (lambda > 0.0) { /* fenwick_get_total sometimes returns -0.0 */
            re_t_wait = gsl_ran_exponential(self->rng, 1.0 / lambda);
        }

        /* For now, don't compute the GC rates if 0 to avoid slowing down other
         * simulations */
        if (self->gene_conversion_rate > 0) {
            /* Gene conversion within segments */
            if (recomb_mass > 0.0) {
                lambda = recomb_map_mass_to_position(
                    &self->recomb_map, recomb_mass) * self->gene_conversion_rate;
            } else {
                total_recomb_rate = recomb_map_get_total_recombination_rate(
                                        &self->recomb_map);
                if (total_recomb_rate == 0.0) {
                    printf("recombination rate zero and gene conversion rate > 0 currently not supported\n");
                    assert(total_recomb_rate > 0.0);
                } else {
                    lambda = 0.0;
                }
            }
            gc_in_t_wait = DBL_MAX;
            if (lambda != 0.0) {
                gc_in_t_wait = gsl_ran_exponential(self->rng, 1.0 / lambda);
            }

            /* Gene conversion to the left of initial segments */
            gc_left_t_wait = DBL_MAX;
            lambda = msp_get_cleft_total(self) * self->gene_conversion_rate
                * self->gene_conversion_track_length;
            if (lambda != 0.0) {
                gc_left_t_wait = gsl_ran_exponential(self->rng, 1.0 / lambda);
            }
        }

        /* Common ancestors */
        ca_t_wait = DBL_MAX;
        ca_pop_id = 0;
        for (j = 0; j < self->num_populations; j++) {
            t_temp = self->get_common_ancestor_waiting_time(self, (population_id_t) j, label);
            if (t_temp < ca_t_wait) {
                ca_t_wait = t_temp;
                ca_pop_id = (population_id_t) j;
            }
        }

        /* Migration */
        mig_t_wait = DBL_MAX;
        mig_source_pop = 0;
        mig_dest_pop = 0;
        for (j = 0; j < self->num_populations; j++) {
            n = avl_count(&self->populations[j].ancestors[label]);
            for (k = 0; k < self->num_populations; k++) {
                lambda = n * self->migration_matrix[j * self->num_populations + k];
                if (lambda != 0.0) {
                    t_temp = gsl_ran_exponential(self->rng, 1.0 / lambda);
                    if (t_temp < mig_t_wait) {
                        mig_t_wait = t_temp;
                        /* m[j, k] is the rate at which migrants move from
                         * population k to j forwards in time. Backwards
                         * in time, we move the individual from from
                         * population j into population k.
                         */
                        mig_source_pop = (population_id_t) j;
                        mig_dest_pop = (population_id_t) k;
                    }
                }
            }
        }
        t_wait = GSL_MIN(mig_t_wait,
            GSL_MIN(gc_in_t_wait,
                GSL_MIN(gc_left_t_wait,
                    GSL_MIN(re_t_wait, ca_t_wait))));
        if (self->next_demographic_event == NULL
                && self->next_sampling_event == self->num_sampling_events
                && t_wait == DBL_MAX) {
            ret = MSP_ERR_INFINITE_WAITING_TIME;
            goto out;
        }
        t_temp = self->time + t_wait;

        sampling_event_time = DBL_MAX;
        if (self->next_sampling_event < self->num_sampling_events) {
            sampling_event_time = self->sampling_events[self->next_sampling_event].time;
        }
        demographic_event_time = DBL_MAX;
        if (self->next_demographic_event != NULL) {
            demographic_event_time = self->next_demographic_event->time;
        }
        /* The simulation state is can only changed from this point on. If
         * any of the events would cause the time to be >= max_time, we exit
         */
        if (sampling_event_time < t_temp
                && sampling_event_time < demographic_event_time) {
            se = &self->sampling_events[self->next_sampling_event];
            if (se->time >= max_time) {
                ret = MSP_EXIT_MAX_TIME;
                break;
            }
            self->time = se->time;
            /* Add in all samples with this time */
            while (self->next_sampling_event < self->num_sampling_events &&
                    self->sampling_events[self->next_sampling_event].time
                        == sampling_event_time) {
                se = self->sampling_events + self->next_sampling_event;
                ret = msp_insert_sample(self, se->sample, se->population_id);
                if (ret != 0) {
                    goto out;
                }
                self->next_sampling_event++;
            }
        } else if (demographic_event_time < t_temp) {
            if (demographic_event_time >= max_time) {
                ret = MSP_EXIT_MAX_TIME;
                break;
            }
            ret = msp_apply_demographic_events(self);
            if (ret != 0) {
                goto out;
            }
        } else {
            if (t_temp >= max_time) {
                ret = MSP_EXIT_MAX_TIME;
                break;
            }
            self->time = t_temp;
            if (re_t_wait == t_wait) {
                ret = msp_recombination_event(self, label, NULL, NULL);
            } else if (gc_in_t_wait == t_wait) {
                ret = msp_gene_conversion_within_event(self, label);
            } else if (gc_left_t_wait == t_wait) {
                ret = msp_gene_conversion_left_event(self, label);
            } else if (ca_t_wait == t_wait) {
                ret = self->common_ancestor_event(self, ca_pop_id, label);
                if (ret == 1) {
                    /* The CA event has signalled that this event should be rejected */
                    self->time -= t_wait;
                    ret = 0;
                }
            } else {
                ret = msp_migration_event(self, mig_source_pop, mig_dest_pop);
            }
            if (ret != 0) {
                goto out;
            }
        }
    }
out:
    return ret;
}

int MSP_WARN_UNUSED
msp_pedigree_climb(msp_t *self)
{
    int ret, ix;
    char id_str[100];
    size_t i, j;
    tsk_size_t id_str_len;
    tsk_id_t node_tsk_id = TSK_NULL;
    individual_t *ind = NULL;
    individual_t *parent = NULL;
    segment_t *merged_segment = NULL;
    segment_t *u[2]; // Will need to update for different ploidy
    avl_tree_t *segments = NULL;
    /* avl_node_t *node; */

    assert(self->num_populations == 1);
    assert(avl_count(&self->pedigree->ind_heap) > 0);
    assert(self->pedigree->state == MSP_PED_STATE_UNCLIMBED);

    self->pedigree->state = MSP_PED_STATE_CLIMBING;

    while (avl_count(&self->pedigree->ind_heap) > 0) {
        /* NOTE: We don't yet support early termination - need to properly
         handle moving segments back into population (or possibly keep them
         there in the first place) before we can handle that */
        ret = msp_pedigree_pop_ind(self, &ind);
        if (ret != 0) {
            goto out;
        }
        assert(ind->time >= self->time);
        self->time = ind->time;

        for (i = 0; i < self->pedigree->ploidy; i++) {
            parent = ind->parents[i];
            if (parent != NULL && ind->time >= parent->time) {
                ret = MSP_ERR_TIME_TRAVEL;
                goto out;
            }
            segments = ind->segments + i;

            /* This parent may not have contributed any ancestral material
             * to the samples */
            if (avl_count(segments) == 0) {
                continue;
            }

            /* If the parent did contribute, we add them to the individual table */
            // TODO: This adds the parents of all individuals who are reached by
            // climbing - wasteful, since few visited individuals become nodes
            // through CA events
            if (parent != NULL && parent->tsk_id == TSK_NULL) {
                sprintf(id_str, "%d", parent->id);
                id_str_len = (tsk_size_t) ceil(log10(parent->id + 1));
                assert(id_str_len > 0);
                ret = tsk_individual_table_add_row(&self->tables->individuals, 0,
                        NULL, 0, id_str, id_str_len);
                if (ret < 0) {
                    goto out;
                }
                parent->tsk_id = ret;
            }
            node_tsk_id = TSK_NULL;
            if (parent != NULL) {
                node_tsk_id = parent->tsk_id;
            }

            /* Merge segments inherited from this ind and recombine */
            // TODO: Make sure population gets properly set when more than one
            ret = msp_merge_ancestors(self, segments, 0, 0, &merged_segment,
                    node_tsk_id);
            if (ret != 0) {
                goto out;
            }
            if (merged_segment == NULL) {
                // This lineage has coalesced
                continue;
            }
            assert(avl_count(segments) == 0);
            assert(merged_segment->prev == NULL);

            /* If parent is NULL, we are at a pedigree founder and we add the
             * lineage back to its original population */
            if (parent == NULL) {
                ret = msp_insert_individual(self, merged_segment);
                if (ret != 0) {
                    goto out;
                }
                continue;
            }

            /* Recombine and climb to segments to the parents */
            if (recomb_map_get_total_recombination_rate(&self->recomb_map) > 0) {
                ret = msp_dtwf_recombine(self, merged_segment, &u[0], &u[1]);
                if (ret != 0) {
                    goto out;
                }
            } else {
                ix = (int) gsl_rng_uniform_int(self->rng, 2);
                u[0] = NULL;
                u[1] = NULL;
                u[ix] = merged_segment;
            }
            for (j = 0; j < self->pedigree->ploidy; j++) {
                if (u[j] == NULL) {
                    continue;
                }
                /* assert(u[j]->prev == NULL); */
                ret = msp_pedigree_add_individual_segment(self, parent, u[j], j);
                if (ret != 0) {
                    goto out;
                }
            }
            if (parent->queued == false) {
                ret = msp_pedigree_push_ind(self, parent);
                if (ret != 0) {
                    goto out;
                }
            }
        }
        ind->merged = true;
    }
    self->pedigree->state = MSP_PED_STATE_CLIMB_COMPLETE;

    ret = 0;
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
    unsigned int segments_to_merge;
    uint32_t N, i, j, k, p;
    size_t segment_mem_offset;
    population_t *pop;
    segment_t *x, *ind1, *ind2;
    segment_t *u[2];
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
        N = (uint32_t) round(
            get_population_size(pop, self->time) * self->model.reference_size);
        if (N == 0) {
            ret = MSP_ERR_DTWF_ZERO_POPULATION_SIZE;
            goto out;
        }

        // Allocate memory for linked list of offspring per parent
        parents = calloc(N, sizeof(segment_list_t *));
        segment_mem = malloc(msp_get_num_ancestors(self) * sizeof(segment_list_t));
        if (parents == NULL || segment_mem == NULL){
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
                if (recomb_map_get_total_recombination_rate(&self->recomb_map) > 0) {
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
            for (i = 0; i < 2; i ++) {
                segments_to_merge = avl_count(&Q[i]);
                if (segments_to_merge == 1) {
                    msp_priority_queue_pop(self, &Q[i]);
                } else if (segments_to_merge >= 2) {
                    msp_remove_individuals_from_population(self, &Q[i]);
                    if (segments_to_merge == 2) {
                        ind1 = msp_priority_queue_pop(self, &Q[i]);
                        ind2 = msp_priority_queue_pop(self, &Q[i]);
                        ret = msp_merge_two_ancestors(self,
                                (population_id_t) j, label, ind1, ind2);
                    } else {
                        ret = msp_merge_ancestors(self,
                                &Q[i], (population_id_t) j, label, NULL, TSK_NULL);
                    }
                }
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
msp_store_simultaneous_migration_events(msp_t *self, avl_tree_t *nodes,
       population_id_t source_pop, label_id_t label) {
    int ret = 0;
    uint32_t j;
    avl_node_t *node;
    avl_tree_t *source;

    source = &self->populations[source_pop].ancestors[label];

    // Choose node to migrate
    j = (uint32_t) gsl_rng_uniform_int(self->rng, avl_count(source));
    node = avl_at(source, j);
    assert(node != NULL);

    avl_unlink_node(source, node);
    node = avl_insert_node(nodes, node);
    assert(node != NULL);

    return ret;
}

static int MSP_WARN_UNUSED
msp_simultaneous_migration_event(msp_t *self, avl_tree_t *nodes,
        population_id_t source_pop, population_id_t dest_pop) {
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
 * 0 if the simulation completed to coalescence
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

    n = malloc(self->num_populations * sizeof(int));
    mig_tmp = malloc(self->num_populations * sizeof(double));
    if (n == NULL || mig_tmp == NULL) {
        ret  = MSP_ERR_NO_MEMORY;
        goto out;
    }

    while (msp_get_num_ancestors(self) > 0) {
        if (events == max_events) {
            ret = MSP_EXIT_MAX_EVENTS;
            break;
        }
        events++;
        if (self->time + 1 >= max_time) {
            ret = MSP_EXIT_MAX_TIME;
            goto out;
        }
        self->time++;

        /* Following SLiM, we perform migrations prior to selecting
         * parents for the current generation */
        node_trees = malloc(self->num_populations * self->num_populations
                * sizeof(avl_tree_t));
        if (node_trees == NULL){
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
            assert(mig_tmp[j] == 0);

            mig_tmp[j] = 1 - sum;
            N = avl_count(&self->populations[j].ancestors[label]);
            gsl_ran_multinomial(
                    self->rng, self->num_populations, N, mig_tmp, n);

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
        while (self->next_demographic_event != NULL &&
                self->next_demographic_event->time <= cur_time) {
            if (self->next_demographic_event->time >= max_time) {
                ret = MSP_EXIT_MAX_TIME;
                goto out;
            }
            ret = msp_apply_demographic_events(self);
            if (ret != 0) {
                goto out;
            }
        }
        self->time = cur_time;
        ret = msp_dtwf_generation(self);
        if (ret != 0) {
            goto out;
        }

        while (self->next_sampling_event < self->num_sampling_events &&
                self->sampling_events[self->next_sampling_event].time <= self->time) {
            se = self->sampling_events + self->next_sampling_event;
            /* The sampling event doesn't modify the tables, so we don't need to
             * catch it here */
            ret = msp_insert_sample(self, se->sample, se->population_id);
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
        assert(avl_count(&self->populations[j].ancestors[1]) == 0);
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
    avl_tree_t *pop = &self->populations[ind->population_id].ancestors[ind->label];
    avl_node_t *node;

    /* Find the this individual in the AVL tree. */
    node = avl_search(pop, ind);
    assert(node != NULL);
    ret = msp_move_individual(self, node, pop, ind->population_id, label);
    return ret;
}

static int
msp_sweep_recombination_event(msp_t *self, label_id_t label,
        double sweep_locus, double population_frequency)
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
    size_t curr_step = 0;
    size_t num_steps;
    double *allele_frequency = NULL;
    double *time = NULL;
    double sweep_locus = model->params.sweep.locus;
    double sweep_dt;
    size_t j = 0;
    double recomb_mass;
    unsigned long events = 0;
    label_id_t label;
    double rec_rates[] = {0.0, 0.0};
    double sweep_pop_sizes[] = {0.0, 0.0};
    double event_prob, event_rand, tmp_rand, e_sum, pop_size;
    double p_coal_b, p_coal_B, total_rate, sweep_pop_tot_rate;
    double p_rec_b, p_rec_B;
    bool sweep_over;

    /* Keep the compiler happy */
    sweep_pop_tot_rate = 0;
    p_coal_b = 0;
    p_coal_B = 0;

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
    if (ret != 0) {
        goto out;
    }
    ret = msp_sweep_initialise(self, allele_frequency[0]);
    if (ret != 0) {
        goto out;
    }

    while (msp_get_num_ancestors(self) > 0 && curr_step < num_steps) {
        events++;
        /* Set pop sizes & rec_rates */
        for (j = 0; j < self->num_labels; j++) {
            label = (label_id_t) j;
            recomb_mass = fenwick_get_total(&self->links[label]);
            sweep_pop_sizes[j] = avl_count(&self->populations[0].ancestors[label]);
            rec_rates[j] = recomb_mass;
        }

        curr_step++;
        event_prob = 1.0;
        event_rand = gsl_rng_uniform(self->rng);
        sweep_over = false;
        while (event_prob > event_rand && curr_step < num_steps && !sweep_over) {
            sweep_dt = time[curr_step] - time[curr_step - 1];
            /* using pop sizes grabbed from get_population_size */
            pop_size = get_population_size(&self->populations[0], time[curr_step]);
            p_coal_B = ((sweep_pop_sizes[1] * (sweep_pop_sizes[1] - 1) ) * 0.5)
                / allele_frequency[curr_step] * sweep_dt / pop_size;
            p_coal_b = ((sweep_pop_sizes[0] * (sweep_pop_sizes[0] - 1) ) * 0.5)
                / (1.0 - allele_frequency[curr_step]) * sweep_dt / pop_size;
            p_rec_b = rec_rates[0] * sweep_dt;
            p_rec_B = rec_rates[1] * sweep_dt;
            sweep_pop_tot_rate = p_coal_b + p_coal_B + p_rec_b + p_rec_B;
            /* doing this to build in generality if we want >1 pop */
            total_rate = sweep_pop_tot_rate;
            event_prob *= 1.0 - total_rate;
            curr_step++;
            sweep_over = total_rate == 0;
        }
        if (sweep_over) {
            break;
        }
        /* passed check, choose event */
        tmp_rand = gsl_rng_uniform(self->rng);
        e_sum = p_coal_b;
        self->time = time[curr_step - 1];
        if (tmp_rand < e_sum / sweep_pop_tot_rate) {
            /* coalescent in b background */
            ret = self->common_ancestor_event(self, 0, 0);
        } else {
            e_sum += p_coal_B;
            if (tmp_rand < e_sum / sweep_pop_tot_rate) {
                /* coalescent in B background */
                ret = self->common_ancestor_event(self, 0, 1);
            } else {
                e_sum += rec_rates[0];
                if (tmp_rand < e_sum / sweep_pop_tot_rate) {
                    /* recomb in b background */
                    ret = msp_sweep_recombination_event(self, 0, sweep_locus,
                        (1.0 - allele_frequency[curr_step - 1]));
                } else {
                    /* recomb in B background */
                    ret = msp_sweep_recombination_event(self, 1, sweep_locus,
                        allele_frequency[curr_step - 1]);
                }
            }
        }
        if (ret != 0){
            goto out;
        }
        /*msp_print_state(self, stdout);*/
    }
    /* Check if any demographic events should have happened during the
     * event and raise an error if so. This is to keep computing population
     * sizes simple */
    if (self->next_demographic_event != NULL &&
            self->next_demographic_event->time <= self->time) {
        ret = MSP_ERR_EVENTS_DURING_SWEEP;
        goto out;
    }
    /* We could implement sampling events during a sweep if we wanted
     * to easily enough. Is this a sensible feature? */
    if (self->next_sampling_event < self->num_sampling_events
            && self->sampling_events[self->next_sampling_event].time <= self->time) {
        ret = MSP_ERR_EVENTS_DURING_SWEEP;
        goto out;
    }
    ret = msp_sweep_finalise(self);
    if (ret != 0) {
        goto out;
    }
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
    simulation_model_t *model = &self->model;
    double scaled_time = model->generations_to_model_time(model, max_time);

    if (self->state == MSP_STATE_INITIALISED) {
        self->state = MSP_STATE_SIMULATING;
    }
    if (self->state != MSP_STATE_SIMULATING) {
        ret = MSP_ERR_BAD_STATE;
        goto out;
    }
    if (self->store_full_arg && ! (
            self->model.type == MSP_MODEL_HUDSON
            || self->model.type == MSP_MODEL_SMC
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
        ret = msp_run_dtwf(self, scaled_time, max_events);
    } else if (self->model.type == MSP_MODEL_WF_PED) {
        if (self->pedigree == NULL ||
            self->pedigree->state != MSP_PED_STATE_UNCLIMBED) {
            ret = MSP_ERR_BAD_STATE;
            goto out;
        }
        ret = msp_pedigree_load_pop(self);
        if (ret != 0) {
            goto out;
        }
        ret = msp_pedigree_build_ind_queue(self);
        if (ret != 0) {
            goto out;
        }
        ret = msp_pedigree_climb(self);
        if (ret != 0) {
            goto out;
        }
    } else if (self->model.type == MSP_MODEL_SWEEP) {
        /* FIXME making sweep atomic for now as it's non-rentrant */
        ret = msp_run_sweep(self);
    } else {
        ret = msp_run_coalescent(self, scaled_time, max_events);
    }

    if (ret < 0) {
        goto out;
    }
    if (ret == MSP_EXIT_MAX_TIME) {
        /* Set the time to the max_time specified. If the tables are finalised
         * after this we will get unary edges on the end of each extant node
         * to this point so that the simulation can be resumed accurately.
         */
        self->time = scaled_time;
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
    node_id_t node;
    int64_t edge_start;
    tsk_node_table_t *nodes = &self->tables->nodes;
    const double current_time = self->model.model_time_to_generations(&self->model,
            self->time);
    tsk_bookmark_t bookmark;

    for (pop = 0; pop < (population_id_t) self->num_populations; pop++) {
        for (label = 0; label < (label_id_t) self->num_labels; label++) {
            for (a = self->populations[pop].ancestors[label].head; a != NULL; a = a->next) {
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
                    node = tsk_node_table_add_row(nodes, 0, current_time, pop,
                            TSK_NULL, NULL, 0);
                    if (node < 0) {
                        ret = msp_set_tsk_error(node);
                        goto out;
                    }
                }

                /* For every segment add an edge pointing to this new node */
                for (seg = (segment_t *) a->item; seg != NULL; seg = seg->next) {
                    if (seg->value != node) {
                        assert(nodes->time[node] > nodes->time[seg->value]);
                        ret = tsk_edge_table_add_row(&self->tables->edges,
                              seg->left, seg->right, node, seg->value);
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
    /* FIXME!!! This is a *horrible* ugly hack to work around sort not accepting
     * migrations. We know the migrations are sorted, so we shouldn't need to
     * sort them anyway, so we should be able to set bookmark.migrations =
     * tables.migrations.num_rows.
     * See https://github.com/tskit-dev/tskit/issues/117
     */
    tsk_size_t tmp_count = self->tables->migrations.num_rows;
    self->tables->migrations.num_rows = 0;
    ret = tsk_table_collection_sort(self->tables, &bookmark, 0);
    self->tables->migrations.num_rows = tmp_count;
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

    if (!msp_is_completed(self)) {
        ret = msp_insert_uncoalesced_edges(self);
        if (ret != 0) {
            goto out;
        }
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
    simulation_model_t *model = &self->model;
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
    if (! first_call && self->next_demographic_event != NULL) {
        de = self->next_demographic_event;

        /* Add in historical samples more recent than next demographic event */
        while (self->next_sampling_event < self->num_sampling_events
                && self->sampling_events[self->next_sampling_event].time <= de->time) {
            se = self->sampling_events + self->next_sampling_event;
            ret = msp_insert_sample(self, se->sample, se->population_id);
            if (ret != 0) {
                goto out;
            }
            self->next_sampling_event++;
        }

        ret = msp_apply_demographic_events(self);
        if (ret != 0) {
            goto out;
        }
    }
    if (self->next_demographic_event != NULL) {
        t = self->next_demographic_event->time;
    }
    *end_time = model->model_time_to_generations(model, t);
out:
    return ret;
}

/* Used for high-level debugging. */
int MSP_WARN_UNUSED
msp_compute_population_size(msp_t *self, size_t population_id, double time,
        double *pop_size)
{
    int ret = 0;
    population_t *pop;
    simulation_model_t *model = &self->model;
    double dt;

    if (population_id > self->num_populations) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    pop = &self->populations[population_id];
    if (pop->growth_rate == 0.0) {
        *pop_size = model->reference_size * pop->initial_size;
    } else {
        dt = model->generations_to_model_time(model, time) - pop->start_time;
        *pop_size = model->reference_size * pop->initial_size
            * exp(-pop->growth_rate * dt);
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
            ret = "wf_ped";
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

size_t
msp_get_num_samples(msp_t *self)
{
    return (size_t) self->num_samples;
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
msp_get_num_ancestors(msp_t *self)
{
    size_t n = 0;
    size_t j;
    label_id_t label;

    for (j = 0; j < self->num_populations; j++) {
        for (label = 0; label < (label_id_t) self->num_labels; label++) {
            n += avl_count(&self->populations[j].ancestors[label]);
        }
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
            for (node = population_ancestors->head;
                    node != NULL; node = node->next) {
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
        breakpoints[j] = (size_t) nm->left;
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
    simulation_model_t *model = &self->model;

    for (j = 0; j < N * N; j++) {
        migration_matrix[j] = model->model_rate_to_generation_rate(
                model, self->migration_matrix[j]);
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
msp_get_samples(msp_t *self, sample_t **samples)
{
    *samples = self->samples;
    return 0;
}

int MSP_WARN_UNUSED
msp_get_population_configuration(msp_t *self, size_t population_id, double *initial_size,
        double *growth_rate)
{
    int ret = 0;
    population_t *pop;
    simulation_model_t *model = &self->model;

    if (population_id > self->num_populations) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    pop = &self->populations[population_id];
    *initial_size = model->reference_size * pop->initial_size;
    *growth_rate = model->model_rate_to_generation_rate(model, pop->growth_rate);
out:
    return ret;
}

double
msp_get_time(msp_t *self)
{
    simulation_model_t *model = &self->model;
    return model->model_time_to_generations(model, self->time);
}

double
msp_get_gene_conversion_rate(msp_t *self)
{
    simulation_model_t *model = &self->model;
    return model->model_rate_to_generation_rate(model, self->gene_conversion_rate);
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
    ret_event->time = self->model.generations_to_model_time(&self->model, time);
    /* now insert this event at the end of the chain. */
    if (self->demographic_events_head == NULL) {
        self->demographic_events_head = ret_event;
        self->demographic_events_tail = ret_event;
    } else {
        assert(self->demographic_events_tail != NULL);
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
msp_change_single_population_parameters(msp_t *self, size_t population_id,
        double time, double initial_size, double growth_rate)
{
    int ret = 0;
    double dt;
    population_t *pop;
    simulation_model_t *model = &self->model;

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
        pop->initial_size = initial_size / model->reference_size;
    }
    /* Do not change the growth_rate unless it is specified */
    if (!gsl_isnan(growth_rate)) {
        pop->growth_rate = model->generation_rate_to_model_rate(model, growth_rate);
    }
    pop->start_time = time;
out:
    return ret;
}

static int
msp_change_population_parameters(msp_t *self, demographic_event_t *event)
{
    int ret = 0;
    population_id_t pid = event->params.population_parameters_change.population_id;
    double initial_size =
        event->params.population_parameters_change.initial_size;
    double growth_rate =
        event->params.population_parameters_change.growth_rate;

    if (pid == -1) {
        for (pid = 0; pid < (int) self->num_populations; pid++) {
            ret = msp_change_single_population_parameters(self, (size_t) pid,
                    event->time, initial_size, growth_rate);
            if (ret != 0) {
                goto out;
            }
        }
    } else {
        ret = msp_change_single_population_parameters(self, (size_t) pid,
                event->time, initial_size, growth_rate);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

static void
msp_print_population_parameters_change(msp_t * MSP_UNUSED(self),
        demographic_event_t *event, FILE *out)
{
    fprintf(out,
            "%f\tpopulation_parameters_change: %d -> initial_size=%f, growth_rate=%f\n",
            event->time,
            (int) event->params.population_parameters_change.population_id,
            event->params.population_parameters_change.initial_size,
            event->params.population_parameters_change.growth_rate);
}

/* Adds a population parameter change event. Time and growth_rate are measured in
 * units of generations, and the initial size is an absolute value. */
int
msp_add_population_parameters_change(msp_t *self, double time, int population_id,
        double initial_size, double growth_rate)
{
    int ret = -1;
    demographic_event_t *de;
    int N = (int) self->num_populations;

    if (population_id < -1 || population_id >= N) {
        ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
        goto out;
    }
    if (initial_size <= 0) {
        assert(! gsl_isnan(initial_size));
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
    de->params.population_parameters_change.population_id = population_id;
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
    simulation_model_t *model = &self->model;
    double rate = model->generation_rate_to_model_rate(model,
        event->params.migration_rate_change.migration_rate);

    if (index == -1) {
        for (index = 0; index < N * N; index++) {
            if (index % (N + 1) != 0) {
                ret = msp_change_migration_matrix_entry(self, (size_t) index,
                        rate);
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
msp_print_migration_rate_change(msp_t * MSP_UNUSED(self),
        demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tmigration_rate_change: %d -> %f\n",
            event->time,
            event->params.migration_rate_change.matrix_index,
            event->params.migration_rate_change.migration_rate);
}

/* Add a migration rate change event. Time and migration rate are measured in
 * units of generations. */
int MSP_WARN_UNUSED
msp_add_migration_rate_change(msp_t *self, double time, int matrix_index,
        double migration_rate)
{
    int ret = -1;
    demographic_event_t *de;
    int N = (int) self->num_populations;

    if (matrix_index < -1 || matrix_index >= N * N) {
        ret = MSP_ERR_BAD_MIGRATION_MATRIX_INDEX;
        goto out;
    }
    if (migration_rate < 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (matrix_index % (N + 1) == 0) {
        ret = MSP_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX;
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
msp_print_mass_migration(msp_t * MSP_UNUSED(self), demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tmass_migration: %d -> %d p = %f\n",
            event->time,
            (int) event->params.mass_migration.source,
            (int) event->params.mass_migration.destination,
            event->params.mass_migration.proportion);
}

/* Adds a mass migration event. Time is measured in units of generations */
int MSP_WARN_UNUSED
msp_add_mass_migration(msp_t *self, double time, int source, int destination,
        double proportion)
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

/* Simple bottlenecks. At some time we coalesce a fraction of the
 * extant lineages into a single ancestor. */

static int
msp_simple_bottleneck(msp_t *self, demographic_event_t *event)
{
    int ret = 0;
    population_id_t population_id = event->params.simple_bottleneck.population_id;
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
            assert(q_node != NULL);
        }
        node = next;
    }
    ret = msp_merge_ancestors(self, &Q, population_id, label, NULL, TSK_NULL);
out:
    return ret;
}

static void
msp_print_simple_bottleneck(msp_t * MSP_UNUSED(self), demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tsimple_bottleneck: %d I = %f\n",
            event->time,
            (int) event->params.simple_bottleneck.population_id,
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
    de->params.simple_bottleneck.population_id = population_id;
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
    population_id_t population_id = event->params.instantaneous_bottleneck.population_id;
    double T2 = event->params.instantaneous_bottleneck.strength;
    population_id_t N = (population_id_t) self->num_populations;
    node_id_t *lineages = NULL;
    node_id_t *pi = NULL;
    avl_node_t **avl_nodes = NULL;
    avl_tree_t *sets = NULL;
    node_id_t u, parent;
    uint32_t j, k, n, num_roots;
    double rate, t;
    avl_tree_t *pop;
    avl_node_t *node, *set_node;
    segment_t *individual;
    label_id_t label = 0; /* For now only support label 0 */

    /* This should have been caught on adding the event */
    if (population_id < 0 || population_id >= N) {
        ret = MSP_ERR_ASSERTION_FAILED;
        goto out;
    }
    if (self->model.type == MSP_MODEL_DTWF) {
        ret = MSP_ERR_DTWF_UNSUPPORTED_BOTTLENECK;
        goto out;
    }
    pop = &self->populations[population_id].ancestors[label];
    n = avl_count(pop);
    lineages = malloc(n * sizeof(node_id_t));
    avl_nodes = malloc(n * sizeof(avl_node_t *));
    pi = malloc(2 * n * sizeof(node_id_t));
    sets = malloc(2 * n * sizeof(avl_tree_t));
    if (lineages == NULL || avl_nodes == NULL || pi == NULL
            || sets == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (u = 0; u < (node_id_t) n; u++) {
        lineages[u] = u;
    }
    for (u = 0; u < (node_id_t) (2 * n); u++) {
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
    parent = (node_id_t) n;
    while (j > 0) {
        rate = j + 1;
        rate = rate * j;
        t += msp_get_common_ancestor_waiting_time_from_rate(self,
                &self->populations[population_id], rate);
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
        k = j > 0 ? (uint32_t) gsl_rng_uniform_int(self->rng, j): 0;
        pi[lineages[k]] = parent;
        lineages[k] = parent;
        parent++;
    }
    num_roots = j + 1;
    for (j = 0; j < num_roots; j++) {
        if (lineages[j] >= (node_id_t) n) {
            avl_init_tree(&sets[lineages[j]], cmp_segment_queue, NULL);
        }
    }

    /* Assign each lineage to the set corresponding to a given root.
     * For any root < n, this lineages has not been affected, so we
     * leave it alone.
     */
    for (j = 0; j < n; j++) {
        u = (node_id_t) j;
        while (pi[u] != TSK_NULL) {
            u = pi[u];
        }
        if (u >= (node_id_t) n) {
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
            assert(set_node != NULL);
        }
    }
    for (j = 0; j < num_roots; j++) {
        if (lineages[j] >= (node_id_t) n) {
            ret = msp_merge_ancestors(self, &sets[lineages[j]], population_id, label, NULL, TSK_NULL);
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
msp_print_instantaneous_bottleneck(msp_t *MSP_UNUSED(self),
        demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tinstantaneous_bottleneck: %d T2 = %f\n",
            event->time,
            (int) event->params.instantaneous_bottleneck.population_id,
            event->params.instantaneous_bottleneck.strength);
}

/* Add an instantaneous bottleneck event. Time and strength are measured in generations
 */
int MSP_WARN_UNUSED
msp_add_instantaneous_bottleneck(msp_t *self, double time, int population_id,
        double strength)
{
    int ret = 0;
    demographic_event_t *de;
    int N = (int) self->num_populations;
    simulation_model_t *model = &self->model;

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
    de->params.instantaneous_bottleneck.population_id = population_id;
    de->params.instantaneous_bottleneck.strength =
        model->generations_to_model_time(model, strength);
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
    double time = self->model.model_time_to_generations(&self->model, event->time);
    avl_tree_t *ancestors;
    avl_node_t *node;
    segment_t *seg;
    int i, j;
    node_id_t u;

    for (i = 0; i < (int) self->num_populations; i++) {
        for (j = 0; j < (int) self->num_labels; j++) {

            // Get segment from an ancestor in a population.
            ancestors = &self->populations[i].ancestors[j];
            node = ancestors->head;

            while (node != NULL){
                seg = (segment_t *) node->item;

                while (seg != NULL) {
                    // Add an edge to the edge table.
                    ret = msp_flush_edges(self);
                    if (ret != 0) {
                        goto out;
                    }
                    ret = tsk_node_table_add_row(&self->tables->nodes,
                            MSP_NODE_IS_CEN_EVENT, time, (population_id_t) i, TSK_NULL,
                            NULL, 0);
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
msp_print_census_event(msp_t * MSP_UNUSED(self), demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tcensus_event:\n",
            event->time);
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
std_model_time_to_generations(simulation_model_t *model, double t)
{
    return 4 * model->reference_size * t;
}

static double
std_generations_to_model_time(simulation_model_t *model, double g)
{
    return g / (4 * model->reference_size);
}

static double
std_generation_rate_to_model_rate(simulation_model_t *model, double rate)
{
    return rate * 4 * model->reference_size;
}

static double
std_model_rate_to_generation_rate(simulation_model_t *model, double rate)
{
    return rate / (4 * model->reference_size);
}

static double
msp_std_get_common_ancestor_waiting_time(msp_t *self, population_id_t pop_id,
        label_id_t label)
{
    population_t *pop = &self->populations[pop_id];
    double n = (double) avl_count(&pop->ancestors[label]);
    double lambda = n * (n - 1.0);

    return msp_get_common_ancestor_waiting_time_from_rate(self, pop, lambda);
}

static int MSP_WARN_UNUSED
msp_std_common_ancestor_event(msp_t *self, population_id_t population_id,
        label_id_t label)
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
    assert(x_node != NULL);
    x = (segment_t *) x_node->item;
    avl_unlink_node(ancestors, x_node);
    j = (uint32_t) gsl_rng_uniform_int(self->rng, n - 1);
    y_node = avl_at(ancestors, j);
    assert(y_node != NULL);
    y = (segment_t *) y_node->item;
    avl_unlink_node(ancestors, y_node);

    /* For SMC and SMC' models we reject some events to get the required
     * distribution. */
    if (msp_reject_ca_event(self, x, y)) {
        self->num_rejected_ca_events++;
        /* insert x and y back into the population */
        assert(x_node->item == x);
        node = avl_insert_node(ancestors, x_node);
        assert(node != NULL);
        assert(y_node->item == y);
        node = avl_insert_node(ancestors, y_node);
        assert(node != NULL);
    } else {
        self->num_ca_events++;
        msp_free_avl_node(self, x_node);
        msp_free_avl_node(self, y_node);
        ret = msp_merge_two_ancestors(self, population_id, label, x, y);
    }
    return ret;
}

/**************************************************************
 * Dirac coalescent
 *
 * TODO provide backround and documentation.
 **************************************************************/

static double
dirac_model_time_to_generations(simulation_model_t *model, double t)
{
    double x = 1 + (model->params.dirac_coalescent.c * gsl_pow_2(
               model->params.dirac_coalescent.psi));
    return 4 * gsl_pow_2(model->reference_size) * t / x;
}

static double
dirac_generations_to_model_time(simulation_model_t *model, double g)
{

    double x = 1 + (model->params.dirac_coalescent.c * gsl_pow_2(
               model->params.dirac_coalescent.psi));
    return g * x / (4 * gsl_pow_2(model->reference_size));
}

static double
dirac_generation_rate_to_model_rate(simulation_model_t *model, double rate)
{
    double x = 1 + (model->params.dirac_coalescent.c * gsl_pow_2(
               model->params.dirac_coalescent.psi));
    return rate * 4 * model->reference_size / x;
}

static double
dirac_model_rate_to_generation_rate(simulation_model_t *model, double rate)
{
    double x = 1 + (model->params.dirac_coalescent.c * gsl_pow_2(
               model->params.dirac_coalescent.psi));
    return rate * x / ( 4 * model->reference_size);
}

static double
msp_dirac_get_common_ancestor_waiting_time(msp_t *self, population_id_t pop_id,
        label_id_t label)
{
    population_t *pop = &self->populations[pop_id];
    unsigned int n = (unsigned int) avl_count(&pop->ancestors[label]);
    double c = self->model.params.dirac_coalescent.c;
    double lambda = 2 * (gsl_sf_choose(n, 2) + c);

    return msp_get_common_ancestor_waiting_time_from_rate(self, pop, lambda);
}

static int MSP_WARN_UNUSED
msp_dirac_common_ancestor_event(msp_t *self, population_id_t pop_id, label_id_t label)
{
    int ret = 0;
    uint32_t j, n, num_participants;
    avl_tree_t *ancestors, Q[4]; /* MSVC won't let us use num_pots here */
    avl_node_t *x_node, *y_node;
    segment_t *x, *y;
    double nC2, p;
    double psi = self->model.params.dirac_coalescent.psi;

    ancestors = &self->populations[pop_id].ancestors[label];
    n = avl_count(ancestors);
    nC2 = gsl_sf_choose(n, 2);
    p = (nC2 / (nC2 + self->model.params.dirac_coalescent.c));
    if (gsl_rng_uniform(self->rng) < p) {
        /* Choose x and y */
        n = avl_count(ancestors);
        j = (uint32_t) gsl_rng_uniform_int(self->rng, n);
        x_node = avl_at(ancestors, j);
        assert(x_node != NULL);
        x = (segment_t *) x_node->item;
        avl_unlink_node(ancestors, x_node);
        j = (uint32_t) gsl_rng_uniform_int(self->rng, n - 1);
        y_node = avl_at(ancestors, j);
        assert(y_node != NULL);
        y = (segment_t *) y_node->item;
        avl_unlink_node(ancestors, y_node);
        self->num_ca_events++;
        msp_free_avl_node(self, x_node);
        msp_free_avl_node(self, y_node);
        ret = msp_merge_two_ancestors(self, pop_id, label, x, y);
    } else {
        for (j = 0; j < 4; j++){
            avl_init_tree(&Q[j], cmp_segment_queue, NULL);
        }
        num_participants = gsl_ran_binomial(self->rng, psi, n);
        ret = msp_multi_merger_common_ancestor_event(self, ancestors, Q,
                                                     num_participants);
        if (ret < 0) {
            goto out;
        }
        /* All the lineages that have been assigned to the particular pots can now be
         * merged.
         */
        for (j = 0; j < 4; j++){
            ret = msp_merge_ancestors(self, &Q[j], pop_id, label, NULL, TSK_NULL);
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
beta_model_compute_timescale(simulation_model_t *model)
{
    double alpha = model->params.beta_coalescent.alpha;
    double truncation_point = model->params.beta_coalescent.truncation_point;
    double reference_size = model->reference_size;
    double m = 2.0 + exp(alpha * log(2) + (1 - alpha) * log(3) - log(alpha - 1));
    double timescale = exp(log(alpha) - alpha * log(m)
        - (alpha - 1) * log(reference_size))
        * gsl_sf_beta_inc(2 - alpha, alpha, truncation_point);
    return timescale;
}

static double
beta_model_time_to_generations(simulation_model_t *model, double t)
{
    return 4 * beta_model_compute_timescale(model) * t;
}

static double
beta_generations_to_model_time(simulation_model_t *model, double g)
{
    return g / (4 * beta_model_compute_timescale(model));
}

static double
beta_generation_rate_to_model_rate(simulation_model_t *model, double rate)
{
    return rate * 4 * beta_model_compute_timescale(model);
}

static double
beta_model_rate_to_generation_rate(simulation_model_t *model, double rate)
{
    return rate / (4 * beta_model_compute_timescale(model));
}

static double
msp_beta_get_common_ancestor_waiting_time(msp_t *self, population_id_t pop_id,
        label_id_t label)
{
    population_t *pop = &self->populations[pop_id];
    unsigned int n = (unsigned int) avl_count(&pop->ancestors[label]);
    /* Factor of 4 because only 1/4 of binary events result in a merger due to
     * diploidy, and 2 for consistency with the hudson model */
    double lambda = 8 * gsl_sf_choose(n, 2);
    double result = msp_get_common_ancestor_waiting_time_from_rate(self, pop, lambda);
    return result;
}

int MSP_WARN_UNUSED
msp_multi_merger_common_ancestor_event(msp_t *self, avl_tree_t *ancestors,
                                       avl_tree_t *Q, uint32_t k)
{
    int ret = 0;
    uint32_t j, i, l;
    avl_node_t  *node, *q_node;
    segment_t *u;
    uint32_t pot_size;
    uint32_t cumul_pot_size = 0;

    /* In the multiple merger regime we have four different 'pots' that
     * lineages get assigned to, where all lineages in a given pot are merged into
     * a common ancestor.
     */
    for (i = 0; i < 4; i++) {
        pot_size = gsl_ran_binomial(self->rng, 1.0 / (4.0 - i), k - cumul_pot_size);
        cumul_pot_size += pot_size;
        if (pot_size > 1) {
            for (l = 0; l < pot_size; l++) {
                j = (uint32_t) gsl_rng_uniform_int(self->rng, avl_count(ancestors));
                node = avl_at(ancestors, j);
                assert(node != NULL);

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
                assert(q_node != NULL);
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
    uint32_t j, n, num_participants;
    avl_tree_t *ancestors, Q[4]; /* MSVC won't let us use num_pots here */
    double beta_x, u, increment;

    for (j = 0; j < 4; j++){
        avl_init_tree(&Q[j], cmp_segment_queue, NULL);
    }
    ancestors = &self->populations[pop_id].ancestors[label];
    n = avl_count(ancestors);
    beta_x = ran_inc_beta(self->rng,
                 2.0 - self->model.params.beta_coalescent.alpha,
                 self->model.params.beta_coalescent.alpha,
                 self->model.params.beta_coalescent.truncation_point);

    /* We calculate the probability of accepting the event */
    if (beta_x > 1e-9) {
        u = (n - 1) * log(1 - beta_x)
            + log(1 + (n - 1) * beta_x);
        u = exp(log(1 - exp(u)) - 2 * log(beta_x)
            - gsl_sf_lnchoose(n, 2));
    } else {
        /* For very small values of beta_x we need a polynomial expansion
         * for numerical stability */
        u = 0;
        for (j = 2; j <= n; j += 2) {
            increment = (j - 1) * exp(gsl_sf_lnchoose(n, j)
                + (j - 2) * log(beta_x));
            if (increment / u < 1e-12) {
                /* We truncate the expansion adaptively once the increment
                 * becomes negligible. */
                break;
            }
            u += increment;
        }
        for (j = 3; j <= n; j += 2) {
            increment = (j - 1) * exp(gsl_sf_lnchoose(n, j)
                + (j - 2) * log(beta_x));
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
        } while (gsl_rng_uniform(self->rng) >
            1 / gsl_sf_choose(num_participants, 2));

        ret = msp_multi_merger_common_ancestor_event(self, ancestors, Q,
                                                     num_participants);
        if (ret < 0) {
            goto out;
        }

        /* All the lineages that have been assigned to the particular pots can now be
        * merged.
        */
        for (j = 0; j < 4; j++) {
            ret = msp_merge_ancestors(self, &Q[j], pop_id, label, NULL, TSK_NULL);
            if (ret < 0) {
                goto out;
            }
        }
    }

out:
    return ret;
}


/**************************************************************
 * Discrete Time Wright Fisher
 *
 * TODO provide background and documentation.
 **************************************************************/
static double
dtwf_model_time_to_generations(simulation_model_t *MSP_UNUSED(model), double t)
{
    return t;
}

static double
dtwf_generations_to_model_time(simulation_model_t *MSP_UNUSED(model), double g)
{
    return g;
}

static double
dtwf_generation_rate_to_model_rate(simulation_model_t *MSP_UNUSED(model), double rate)
{
    return rate;
}

static double
dtwf_model_rate_to_generation_rate(simulation_model_t *MSP_UNUSED(model), double rate)
{
    return rate;
}

/**************************************************************
 * Allele frequency trajectory simulation for genic selection
 *
 **************************************************************/

static double
genic_selection_stochastic_forwards(double dt, double freq, double alpha, double u)
{
    double ux = (alpha * freq * (1 - freq)) / tanh(alpha * freq);
    int sign = u < 0.5? 1: -1;
    return freq + (ux * dt) + sign * sqrt(freq * (1.0 - freq) * dt);
}

static int
genic_selection_generate_trajectory(sweep_t *self, msp_t *simulator,
        size_t *ret_num_steps, double **ret_time, double **ret_allele_frequency)
{
    int ret = 0;
    genic_selection_trajectory_t trajectory =
        self->trajectory_params.genic_selection_trajectory;
    gsl_rng *rng = simulator->rng;
    size_t max_steps = 64;
    double *time = malloc(max_steps * sizeof(*time));
    double *allele_frequency = malloc(max_steps * sizeof(*allele_frequency));
    double x, t, *tmp;
    size_t num_steps;
    double current_size = 1.0;

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
    t = simulator->time;
    num_steps = 0;
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
        time[num_steps] = t;
        allele_frequency[num_steps] = x;
        x = 1.0 - genic_selection_stochastic_forwards(
                trajectory.dt, 1.0 - x, trajectory.alpha * current_size,
                gsl_rng_uniform(rng));
        t += trajectory.dt;
        num_steps++;
    }
    assert(num_steps < max_steps); /* num_steps + 1 above guarantees this */
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
    genic_selection_trajectory_t *trajectory =
        &self->trajectory_params.genic_selection_trajectory;

    fprintf(out, "\tGenic selection trajectory\n");
    fprintf(out, "\t\tstart_frequency = %f\n", trajectory->start_frequency);
    fprintf(out, "\t\tend_frequency = %f\n", trajectory->end_frequency);
    fprintf(out, "\t\talpha = %f\n", trajectory->alpha);
    fprintf(out, "\t\tdt = %f\n", trajectory->dt);
}

/**************************************************************
 * Public API for setting simulation models.
 **************************************************************/

static void
msp_reset_segment_chain_masses(msp_t *self, segment_t *head)
{
    segment_t *seg, *prev;

    head->left_mass = recomb_map_position_to_mass(&self->recomb_map, head->left);
    head->right_mass = recomb_map_position_to_mass(&self->recomb_map, head->right);
    msp_set_single_segment_mass(self, head);

    prev = head;
    seg = head->next;
    while (seg != NULL) {
        if (seg->left == prev->right) {
            seg->left_mass = prev->right_mass;
        } else {
            seg->left_mass = recomb_map_position_to_mass(&self->recomb_map, seg->left);
        }
        seg->right_mass = recomb_map_position_to_mass(&self->recomb_map, seg->right);
        msp_set_segment_mass(self, seg, prev);

        prev = seg;
        seg = seg->next;
    }
}

static void
msp_reset_all_segment_masses(msp_t *self)
{
    uint32_t j;
    avl_tree_t *ancestors;
    avl_node_t *node;
    for (j = 0; j < self->num_populations; j++) {
        ancestors = self->populations[j].ancestors;
        for (node = ancestors->head; node != NULL; node = node->next) {
            msp_reset_segment_chain_masses(self, (segment_t *) node->item);
        }
    }
}

/* Unscale all times and rates from the current model time to generations. */
static int
msp_unscale_model_times(msp_t *self)
{
    uint32_t j;
    simulation_model_t *model = &self->model;
    demographic_event_t *de;

    self->start_time = self->model.model_time_to_generations(model, self->start_time);
    self->time = self->model.model_time_to_generations(model, self->time);
    self->gene_conversion_rate = self->model.model_rate_to_generation_rate(
            model, self->gene_conversion_rate);
    recomb_map_convert_rates(&self->recomb_map,
            (msp_convert_func) self->model.model_rate_to_generation_rate,
            &self->model);
    msp_reset_all_segment_masses(self);
    /* Samples */
    for (j = 0; j < self->num_samples; j++) {
        self->samples[j].time = model->model_time_to_generations(
                model, self->samples[j].time);
    }
    /* Sampling events */
    for (j = 0; j < self->num_sampling_events; j++) {
        self->sampling_events[j].time = model->model_time_to_generations(
                model, self->sampling_events[j].time);
    }
    /* Growth rates and start times for populations */
    for (j = 0; j < self->num_populations; j++) {
        self->populations[j].growth_rate = model->model_rate_to_generation_rate(
                model, self->populations[j].growth_rate);
        self->populations[j].start_time = model->model_time_to_generations(
                model, self->populations[j].start_time);
    }
    /* Migration rates */
    for (j = 0; j < gsl_pow_2(self->num_populations); j++) {
        self->migration_matrix[j] = model->model_rate_to_generation_rate(
                model, self->migration_matrix[j]);
    }
    /* Demographic events */
    for (de = self->demographic_events_head; de != NULL; de = de->next) {
        de->time = model->model_time_to_generations(model, de->time);
    }
    return 0;
}

/* Rescale all times and rates from generations back into the current model time */
static int
msp_rescale_model_times(msp_t *self)
{
    uint32_t j;
    simulation_model_t *model = &self->model;
    demographic_event_t *de;

    self->time = model->generations_to_model_time(model, self->time);
    self->start_time = model->generations_to_model_time(model, self->start_time);
    self->gene_conversion_rate = model->generation_rate_to_model_rate(
            model, self->gene_conversion_rate);
    recomb_map_convert_rates(&self->recomb_map,
            (msp_convert_func) self->model.generation_rate_to_model_rate,
            &self->model);
    msp_reset_all_segment_masses(self);
    /* Samples */
    for (j = 0; j < self->num_samples; j++) {
        self->samples[j].time = model->generations_to_model_time(
                model, self->samples[j].time);
    }
    /* Sampling events */
    for (j = 0; j < self->num_sampling_events; j++) {
        self->sampling_events[j].time = model->generations_to_model_time(
                model, self->sampling_events[j].time);
    }
    /* Growth rates and start times for populations */
    for (j = 0; j < self->num_populations; j++) {
        self->populations[j].growth_rate = model->generation_rate_to_model_rate(
                model, self->populations[j].growth_rate);
        self->populations[j].start_time = model->generations_to_model_time(
                model, self->populations[j].start_time);
    }
    /* Migration rates */
    for (j = 0; j < gsl_pow_2(self->num_populations); j++) {
        self->migration_matrix[j] = model->generation_rate_to_model_rate(
                model, self->migration_matrix[j]);
    }
    /* Demographic events */
    for (de = self->demographic_events_head; de != NULL; de = de->next) {
        de->time = model->generations_to_model_time(model, de->time);
    }
    return 0;
}

static int
msp_set_simulation_model(msp_t *self, int model, double reference_size)
{
    int ret = 0;

    if (model != MSP_MODEL_HUDSON && model != MSP_MODEL_SMC
            && model != MSP_MODEL_SMC_PRIME
            && model != MSP_MODEL_DIRAC
            && model != MSP_MODEL_BETA
            && model != MSP_MODEL_DTWF
            && model != MSP_MODEL_WF_PED
            && model != MSP_MODEL_SWEEP) {
        ret = MSP_ERR_BAD_MODEL;
        goto out;
    }
    if (reference_size <= 0) {
        ret = MSP_ERR_BAD_POPULATION_SIZE;
        goto out;
    }
    /* If this isn't the first time we've set the model, rescale times back
     * to generations so that we can scale them back into the appropriate values
     * after the model has been set */
    if (self->model.type != -1) {
        ret = msp_unscale_model_times(self);
        if (ret != 0) {
            goto out;
        }
        if (self->model.free != NULL) {
            self->model.free(&self->model);
        }
    }
    self->model.type = model;
    self->model.reference_size = reference_size;
    /* For convenience here we set these to what is needed for the standard
     * coalcescent. For other models, these functions must be overwritten
     * with the correct values *before* rescaling time. */
    self->model.model_time_to_generations = std_model_time_to_generations;
    self->model.generations_to_model_time = std_generations_to_model_time;
    self->model.generation_rate_to_model_rate = std_generation_rate_to_model_rate;
    self->model.model_rate_to_generation_rate = std_model_rate_to_generation_rate;
    self->get_common_ancestor_waiting_time = msp_std_get_common_ancestor_waiting_time;
    self->common_ancestor_event = msp_std_common_ancestor_event;
out:
    return ret;
}

int
msp_set_simulation_model_hudson(msp_t *self, double reference_size)
{
    int ret = msp_set_simulation_model(self, MSP_MODEL_HUDSON, reference_size);
    if (ret != 0) {
        goto out;
    }
    ret = msp_rescale_model_times(self);
out:
    return ret;
}

int
msp_set_simulation_model_smc(msp_t *self, double reference_size)
{
    int ret =  msp_set_simulation_model(self, MSP_MODEL_SMC, reference_size);
    if (ret != 0) {
        goto out;
    }
    ret = msp_rescale_model_times(self);
out:
    return ret;
}

int
msp_set_simulation_model_smc_prime(msp_t *self, double reference_size)
{
    int ret =  msp_set_simulation_model(self, MSP_MODEL_SMC_PRIME, reference_size);
    if (ret != 0) {
        goto out;
    }
    ret = msp_rescale_model_times(self);
out:
    return ret;
}

int
msp_set_simulation_model_dtwf(msp_t *self, double reference_size)
{
    int ret = 0;
    ret = msp_set_simulation_model(self, MSP_MODEL_DTWF, reference_size);
    if (ret != 0) {
        goto out;
    }
    self->model.model_time_to_generations = dtwf_model_time_to_generations;
    self->model.generations_to_model_time = dtwf_generations_to_model_time;
    self->model.model_rate_to_generation_rate = dtwf_model_rate_to_generation_rate;
    self->model.generation_rate_to_model_rate = dtwf_generation_rate_to_model_rate;

    ret = msp_rescale_model_times(self);
out:
    return ret;
}

int
msp_set_simulation_model_wf_ped(msp_t *self, double reference_size)
{
    int ret = 0;
    ret = msp_set_simulation_model(self, MSP_MODEL_WF_PED, reference_size);
    if (ret != 0) {
        goto out;
    }
    self->model.model_time_to_generations = dtwf_model_time_to_generations;
    self->model.generations_to_model_time = dtwf_generations_to_model_time;
    self->model.model_rate_to_generation_rate = dtwf_model_rate_to_generation_rate;
    self->model.generation_rate_to_model_rate = dtwf_generation_rate_to_model_rate;

    ret = msp_rescale_model_times(self);
out:
    return ret;
}

int
msp_set_simulation_model_dirac(msp_t *self, double reference_size, double psi, double c)
{
    int ret = 0;
/* We assume to be in the limit where N is infinite*/
    if (psi <= 0 || psi > 1.0) {
        ret = MSP_ERR_BAD_PSI;
        goto out;
    }

    if (c < 0.0 ) {
        ret = MSP_ERR_BAD_C;
        goto out;
    }
    ret = msp_set_simulation_model(self, MSP_MODEL_DIRAC, reference_size);
    if (ret != 0) {
        goto out;
    }
    self->model.params.dirac_coalescent.psi = psi;
    self->model.params.dirac_coalescent.c = c;
    self->model.model_time_to_generations = dirac_model_time_to_generations;
    self->model.generations_to_model_time = dirac_generations_to_model_time;
    self->model.generation_rate_to_model_rate = dirac_generation_rate_to_model_rate;
    self->model.model_rate_to_generation_rate = dirac_model_rate_to_generation_rate;
    self->get_common_ancestor_waiting_time = msp_dirac_get_common_ancestor_waiting_time;
    self->common_ancestor_event = msp_dirac_common_ancestor_event;
    ret = msp_rescale_model_times(self);
out:
    return ret;
}

int
msp_set_simulation_model_beta(msp_t *self, double reference_size, double alpha,
        double truncation_point)
{
    int ret = 0;

    if (alpha <= 1.0 || alpha >= 2.0) {
        ret = MSP_ERR_BAD_BETA_MODEL_ALPHA;
        goto out;
    }

    if (truncation_point <= 0.0 || truncation_point > 1.0) {
        ret = MSP_ERR_BAD_TRUNCATION_POINT;
        goto out;
    }

    ret = msp_set_simulation_model(self, MSP_MODEL_BETA, reference_size);
    if (ret != 0) {
        goto out;
    }

    self->model.params.beta_coalescent.alpha = alpha;
    self->model.params.beta_coalescent.truncation_point = truncation_point;

    self->model.model_time_to_generations = beta_model_time_to_generations;
    self->model.generations_to_model_time = beta_generations_to_model_time;
    self->model.generation_rate_to_model_rate = beta_generation_rate_to_model_rate;
    self->model.model_rate_to_generation_rate = beta_model_rate_to_generation_rate;
    self->get_common_ancestor_waiting_time = msp_beta_get_common_ancestor_waiting_time;
    self->common_ancestor_event = msp_beta_common_ancestor_event;
    ret = msp_rescale_model_times(self);
out:
    return ret;
}

int
msp_set_simulation_model_sweep_genic_selection(msp_t *self, double reference_size,
        double position, double start_frequency, double end_frequency,
        double alpha, double dt)
{
    int ret = 0;
    simulation_model_t *model = &self->model;
    genic_selection_trajectory_t *trajectory =
        &model->params.sweep.trajectory_params.genic_selection_trajectory;
    double L = recomb_map_get_sequence_length(&self->recomb_map);

    /* Check the inputs to make sure they make sense */
    if (position < 0 || position >= L) {
        ret = MSP_ERR_BAD_SWEEP_POSITION;
        goto out;
    }
    if (start_frequency <= 0.0 || start_frequency >= 1.0
            || end_frequency <= 0.0 || end_frequency >= 1.0) {
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
    if (alpha <= 0) {
        ret = MSP_ERR_BAD_SWEEP_GENIC_SELECTION_ALPHA;
        goto out;
    }

    ret = msp_set_simulation_model(self, MSP_MODEL_SWEEP, reference_size);
    if (ret != 0) {
        goto out;
    }
    model->params.sweep.locus = position;
    model->params.sweep.generate_trajectory = genic_selection_generate_trajectory;
    model->params.sweep.print_state = genic_selection_print_state;
    trajectory->start_frequency = start_frequency;
    trajectory->end_frequency = end_frequency;
    /* FIXME alpha must be rescaled here. See
     * https://github.com/tskit-dev/msprime/issues/941
     */
    trajectory->alpha = alpha;
    /* dt value is expressed in generations for consistency; translate to
     * model time. */
    trajectory->dt = self->model.generations_to_model_time(&self->model, dt);

    ret = msp_rescale_model_times(self);
out:
    return ret;
}
