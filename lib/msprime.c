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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>

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

static int
cmp_individual(const void *a, const void *b) {
    const segment_t *ia = (const segment_t *) a;
    const segment_t *ib = (const segment_t *) b;
    return (ia->id > ib->id) - (ia->id < ib->id);
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
    return self->segment_heap.num_blocks;
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
msp_set_event_time_file(msp_t *self, const char *event_time_file)
{
    int ret = 0;
    size_t len = strlen(event_time_file) + 1;

    self->event_time_file_name = malloc(len);
    if (self->event_time_file_name == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    strcpy(self->event_time_file_name, event_time_file);
out:
    return ret;
}

int
msp_set_num_populations(msp_t *self, size_t num_populations)
{
    int ret = 0;
    size_t j;

    if (num_populations < 1 || num_populations > UINT32_MAX) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->num_populations = (uint32_t) num_populations;
    /* Free any memory, if it has been allocated */
    if (self->initial_migration_matrix != NULL) {
        free(self->initial_migration_matrix);
    }
    if (self->migration_matrix != NULL) {
        free(self->migration_matrix);
    }
    if (self->num_migration_events != NULL) {
        free(self->num_migration_events);
    }
    if (self->initial_populations != NULL) {
        free(self->initial_populations);
    }
    if (self->populations != NULL) {
        free(self->populations);
    }
    /* Allocate storage for new num_populations */
    self->initial_migration_matrix = calloc(num_populations * num_populations,
            sizeof(double));
    self->migration_matrix = calloc(num_populations * num_populations,
            sizeof(double));
    self->num_migration_events = calloc(num_populations * num_populations,
            sizeof(size_t));
    self->initial_populations = calloc(num_populations, sizeof(population_t));
    self->populations = calloc(num_populations, sizeof(population_t));
    if (self->migration_matrix == NULL
            || self->initial_migration_matrix == NULL
            || self->num_migration_events == NULL
            || self->initial_populations == NULL
            || self->populations == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (j = 0; j < num_populations; j++) {
        avl_init_tree(&self->populations[j].ancestors, cmp_individual, NULL);
        /* Set the default sizes and growth rates. */
        self->initial_populations[j].growth_rate = 0.0;
        self->initial_populations[j].initial_size = 1.0;
        self->initial_populations[j].start_time = 0.0;
    }
out:
    return ret;
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
        initial_size / model->population_size;
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

static double
msp_genetic_to_phys(msp_t *self, double x)
{
    return recomb_map_genetic_to_phys(self->recomb_map, x);
}

static segment_t * WARN_UNUSED
msp_alloc_segment(msp_t *self, uint32_t left, uint32_t right, node_id_t value,
        population_id_t population_id, segment_t *prev, segment_t *next)
{
    segment_t *seg = NULL;

    if (object_heap_empty(&self->segment_heap)) {
        if (object_heap_expand(&self->segment_heap) != 0) {
            goto out;
        }
        if (fenwick_expand(&self->links, self->segment_block_size) != 0) {
            goto out;
        }
    }
    seg = (segment_t *) object_heap_alloc_object(&self->segment_heap);
    if (seg == NULL) {
        goto out;
    }
    seg->prev = prev;
    seg->next = next;
    seg->left = left;
    seg->right = right;
    seg->value = value;
    seg->population_id = population_id;
out:
    return seg;
}


/* Top level allocators and initialisation */

int
msp_alloc(msp_t *self,
        size_t num_samples, sample_t *samples,
        recomb_map_t *recomb_map, tree_sequence_t *from_ts, gsl_rng *rng) {
    int ret = -1;
    size_t j;

    memset(self, 0, sizeof(msp_t));
    if (rng == NULL || recomb_map == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (from_ts == NULL) {
        if (samples == NULL) {
            ret = MSP_ERR_BAD_PARAM_VALUE;
            goto out;
        }
        if (num_samples < 2) {
            ret = MSP_ERR_BAD_PARAM_VALUE;
            goto out;
        }
    } else {
        /* All samples must be specified in the tree sequence itself. */
        if (num_samples != 0) {
            ret = MSP_ERR_BAD_PARAM_VALUE;
            goto out;
        }
    }

    /* Use the standard coalescent with coalescent time units by default. */
    self->model.type = -1;
    ret = msp_set_simulation_model_hudson(self, 0.25);
    assert(ret == 0);
    self->rng = rng;
    self->recomb_map = recomb_map;
    self->num_loci = recomb_map_get_num_loci(self->recomb_map);

    self->recombination_rate = recomb_map_get_per_locus_recombination_rate(
            self->recomb_map);
    /* And rescale the recombination rate */
    self->recombination_rate = self->model.generation_rate_to_model_rate(
            &self->model, self->recombination_rate);

    ret = table_collection_alloc(&self->tables, MSP_ALLOC_TABLES);
    if (ret != 0) {
        goto out;
    }
    if (from_ts == NULL) {
        assert(samples != NULL);
        assert(num_samples > 1);
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
        assert(from_ts != NULL);

        /* Make a copy of the from_ts */
        ret = tree_sequence_dump_tables(from_ts, &self->tables, 0);
        if (ret != 0) {
            goto out;
        }
        self->from_ts = malloc(sizeof(*self->from_ts));
        if (self->from_ts == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        ret = tree_sequence_load_tables(self->from_ts, &self->tables, 0);
        if (ret != 0) {
            goto out;
        }
        ret = table_collection_record_position(&self->tables, &self->from_position);
        if (ret != 0) {
            goto out;
        }
        if (self->recomb_map->sequence_length != self->tables.sequence_length) {
            ret = MSP_ERR_INCOMPATIBLE_FROM_TS;
            goto out;
        }
    }

    self->start_time = -1;
    /* We have one population by default */
    ret = msp_set_num_populations(self, 1);
    if (ret != 0) {
        goto out;
    }
    /* Set sensible defaults for the sample_config and migration matrix */
    self->initial_migration_matrix[0] = 0.0;
    /* Set the memory defaults */
    self->store_migrations = false;
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

    /* Make sure GSL error handling is turned off */
    gsl_set_error_handler_off();
out:
    return ret;
}

static int
msp_alloc_memory_blocks(msp_t *self)
{
    int ret = 0;

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
    /* allocate the segments and Fenwick tree */
    ret = object_heap_init(&self->segment_heap, sizeof(segment_t),
           self->segment_block_size, segment_init);
    if (ret != 0) {
        goto out;
    }
    ret = fenwick_alloc(&self->links, self->segment_block_size);
    if (ret != 0) {
        goto out;
    }
    /* Allocate the edge records */
    self->num_buffered_edges = 0;
    self->max_buffered_edges = 128;
    self->buffered_edges = malloc(self->max_buffered_edges * sizeof(edge_t));
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
    demographic_event_t *de = self->demographic_events_head;
    demographic_event_t *tmp;

    while (de != NULL) {
        tmp = de->next;
        free(de);
        de = tmp;
    }
    msp_safe_free(self->initial_migration_matrix);
    msp_safe_free(self->migration_matrix);
    msp_safe_free(self->num_migration_events);
    msp_safe_free(self->initial_populations);
    msp_safe_free(self->populations);
    msp_safe_free(self->samples);
    msp_safe_free(self->sampling_events);
    msp_safe_free(self->buffered_edges);
    table_collection_free(&self->tables);
    /* free the object heaps */
    object_heap_free(&self->avl_node_heap);
    object_heap_free(&self->segment_heap);
    object_heap_free(&self->node_mapping_heap);
    fenwick_free(&self->links);
    if (self->from_ts != NULL) {
        tree_sequence_free(self->from_ts);
        free(self->from_ts);
    }
    if (self->model.free != NULL) {
        self->model.free(&self->model);
    }

    msp_safe_free(self->event_time_file_name);
    /* Ignore IO errors here --- should have been caught by earlier calls to flush. */
    if (self->event_time_file != NULL) {
        fclose(self->event_time_file);
        self->event_time_file = NULL;
    }
    ret = 0;
    return ret;
}

static inline avl_node_t * WARN_UNUSED
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
msp_get_segment(msp_t *self, size_t id)
{
    segment_t *u = object_heap_get_object(&self->segment_heap, id - 1);

    assert(u != NULL);
    assert(u->id == id);
    return u;
}

static void
msp_free_segment(msp_t *self, segment_t *seg)
{
    object_heap_free_object(&self->segment_heap, seg);
    fenwick_set_value(&self->links, seg->id, 0);
}

static inline int WARN_UNUSED
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
    node = avl_insert_node(
        &self->populations[u->population_id].ancestors, node);
    assert(node != NULL);
out:
    return ret;
}

static void
msp_print_segment_chain(msp_t * MSP_UNUSED(self), segment_t *head, FILE *out)
{
    segment_t *s = head;

    fprintf(out, "[%d]", (int) s->population_id);
    while (s != NULL) {
        fprintf(out, "[(%d-%d) %d] ", s->left, s->right, (int) s->value);
        s = s->next;
    }
    fprintf(out, "\n");
}

static void
msp_verify_segments(msp_t *self)
{
    int64_t s, ss, total_links, left, right, alt_total_links;
    size_t j;
    size_t total_segments = 0;
    size_t total_avl_nodes = 0;
    avl_node_t *node;
    segment_t *u;

    total_links = 0;
    alt_total_links = 0;
    for (j = 0; j < self->num_populations; j++) {
        node = (&self->populations[j].ancestors)->head;
        while (node != NULL) {
            u = (segment_t *) node->item;
            assert(u->prev == NULL);
            left = u->left;
            while (u != NULL) {
                total_segments++;
                assert(u->population_id == (population_id_t) j);
                assert(u->left < u->right);
                assert(u->right <= self->num_loci);
                if (u->prev != NULL) {
                    s = u->right - u->prev->right;
                } else {
                    s = u->right - u->left - 1;
                }
                ss = fenwick_get_value(&self->links, u->id);
                total_links += ss;
                assert(s == ss);
                if (s == ss) {
                    /* do nothing; just to keep compiler happy - see below also */
                }
                right = u->right;
                u = u->next;
            }
            alt_total_links += right - left - 1;
            node = node->next;
        }
    }
    assert(total_links == fenwick_get_total(&self->links));
    assert(total_links == alt_total_links);
    assert(total_segments == object_heap_get_num_allocated(&self->segment_heap));
    total_avl_nodes = msp_get_num_ancestors(self)
            + avl_count(&self->breakpoints)
            + avl_count(&self->overlap_counts);
    assert(total_avl_nodes == object_heap_get_num_allocated(
                &self->avl_node_heap));
    assert(total_avl_nodes - msp_get_num_ancestors(self)
            == object_heap_get_num_allocated(&self->node_mapping_heap));
    if (total_avl_nodes == total_segments) {
        /* do nothing - this is just to keep the compiler happy when
         * asserts are turned off.
         */
    }
}

static void
msp_verify_overlaps(msp_t *self)
{
    avl_node_t *node;
    node_mapping_t *nm;
    segment_t *u;
    double generations;
    uint32_t j, k, left, right, count;
    /* We check for every locus, so obviously this rules out large numbers
     * of loci. This code should never be called except during testing,
     * so we don't need to recover from malloc failure.
     */
    uint32_t *overlaps = calloc(self->num_loci, sizeof(uint32_t));

    assert(overlaps != NULL);
    /* Add in the counts for any historical samples that haven't been
     * included yet.
     */
    generations = self->model.model_time_to_generations(&self->model, self->time);
    for (j = 0; j < self->num_samples; j++) {
        if (self->samples[j].time > generations) {
            for (k = 0; k < self->num_loci; k++) {
                overlaps[k]++;
            }
        }
    }
    for (j = 0; j < self->num_populations; j++) {
        for (node = (&self->populations[j].ancestors)->head; node != NULL;
                node = node->next) {
            u = (segment_t *) node->item;
            while (u != NULL) {
                for (k = u->left; k < u->right; k++) {
                    overlaps[k]++;
                }
                u = u->next;
            }
        }
    }
    for (node = self->overlap_counts.head; node->next != NULL;
            node = node->next) {
        nm = (node_mapping_t *) node->item;
        left = nm->left;
        right = ((node_mapping_t *) node->next->item)->left;
        count = nm->value;
        for (k = left; k < right; k++) {
            assert(overlaps[k] == count);
        }
    }
    free(overlaps);
}

void
msp_verify(msp_t *self)
{
    msp_verify_segments(self);
    msp_verify_overlaps(self);
}

int
msp_print_state(msp_t *self, FILE *out)
{
    int ret = 0;
    avl_node_t *a;
    node_mapping_t *nm;
    segment_t *u;
    edge_t *edge;
    demographic_event_t *de;
    sampling_event_t *se;
    int64_t v;
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
    fprintf(out, "model population_size = %f\n", self->model.population_size);
    if (self->model.type == MSP_MODEL_BETA) {
        fprintf(out, "\tbeta coalescent parameters: alpha = %f, truncation_point = %f\n",
                self->model.params.beta_coalescent.alpha,
                self->model.params.beta_coalescent.truncation_point);
    } else if (self->model.type == MSP_MODEL_DIRAC) {
        fprintf(out, "\tdirac coalescent parameters: psi = %f, c = %f\n",
                self->model.params.dirac_coalescent.psi,
                self->model.params.dirac_coalescent.c);
    }
    fprintf(out, "n = %d\n", self->num_samples);
    fprintf(out, "m = %d\n", self->num_loci);
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

    fprintf(out, "num_links = %ld\n", (long) fenwick_get_total(&self->links));
    for (j = 0; j < self->num_populations; j++) {
        fprintf(out, "population[%d] = %d\n", j,
            avl_count(&self->populations[j].ancestors));
        fprintf(out, "\tstart_time = %f\n", self->populations[j].start_time);
        fprintf(out, "\tinitial_size = %f\n", self->populations[j].initial_size);
        fprintf(out, "\tgrowth_rate = %f\n", self->populations[j].growth_rate);
    }
    fprintf(out, "Time = %f\n", self->time);
    for (j = 0; j < msp_get_num_ancestors(self); j++) {
        fprintf(out, "\t");
        msp_print_segment_chain(self, ancestors[j], out);
    }
    fprintf(out, "Fenwick tree\n");
    for (j = 1; j <= (uint32_t) fenwick_get_size(&self->links); j++) {
        u = msp_get_segment(self, j);
        v = fenwick_get_value(&self->links, j);
        if (v != 0) {
            fprintf(out, "\t%ld\ti=%d l=%d r=%d v=%d prev=%p next=%p\n", (long) v,
                    (int) u->id, u->left,
                    u->right, (int) u->value, (void *) u->prev, (void *) u->next);
        }
    }
    fprintf(out, "Breakpoints = %d\n", avl_count(&self->breakpoints));
    for (a = self->breakpoints.head; a != NULL; a = a->next) {
        nm = (node_mapping_t *) a->item;
        fprintf(out, "\t%d -> %d\n", nm->left, nm->value);
    }
    fprintf(out, "Overlap count = %d\n", avl_count(&self->overlap_counts));
    for (a = self->overlap_counts.head; a != NULL; a = a->next) {
        nm = (node_mapping_t *) a->item;
        fprintf(out, "\t%d -> %d\n", nm->left, nm->value);
    }
    fprintf(out, "Tables = \n");
    table_collection_print_state(&self->tables, out);

    fprintf(out, "Buffered Edges = %ld\n", (long) self->num_buffered_edges);
    for (j = 0; j < self->num_buffered_edges; j++) {
        edge = &self->buffered_edges[j];
        fprintf(out, "\t%f\t%f\t%d\t%d\n", edge->left, edge->right, edge->parent,
                edge->child);
    }
    fprintf(out, "Memory heaps\n");
    fprintf(out, "avl_node_heap:");
    object_heap_print_state(&self->avl_node_heap, out);
    fprintf(out, "segment_heap:");
    object_heap_print_state(&self->segment_heap, out);
    fprintf(out, "node_mapping_heap:");
    object_heap_print_state(&self->node_mapping_heap, out);
    msp_verify(self);
out:
    if (ancestors != NULL) {
        free(ancestors);
    }
    return ret;
}

static int WARN_UNUSED
msp_record_migration(msp_t *self, uint32_t left, uint32_t right,
        node_id_t node, population_id_t source_pop, population_id_t dest_pop)
{
    int ret = 0;
    double scaled_time = self->model.model_time_to_generations(&self->model, self->time);

    ret = migration_table_add_row(self->tables.migrations,
            msp_genetic_to_phys(self, left),
            msp_genetic_to_phys(self, right),
            node, source_pop, dest_pop, scaled_time);
    if (ret < 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static int WARN_UNUSED
msp_move_individual(msp_t *self, avl_node_t *node, avl_tree_t *source,
        population_id_t dest_pop)
{
    int ret = 0;
    segment_t *ind, *x;

    ind = (segment_t *) node->item;
    avl_unlink_node(source, node);
    msp_free_avl_node(self, node);
    /* Need to set the population_id for each segment. */
    x = ind;
    while (x != NULL) {
        if (self->store_migrations) {
            ret = msp_record_migration(self, x->left, x->right, x->value,
                    x->population_id, dest_pop);
            if (ret != 0) {
                goto out;
            }
        }
        x->population_id = dest_pop;
        x = x->next;
    }
    ret = msp_insert_individual(self, ind);
out:
    return ret;
}

/*
 * Inserts a new breakpoint at the specified locus left.
 */
static int WARN_UNUSED
msp_insert_breakpoint(msp_t *self, uint32_t left)
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


/*
 * Inserts a new overlap_count at the specified locus left, mapping to the
 * specified number of overlapping segments b.
 */
static int WARN_UNUSED
msp_insert_overlap_count(msp_t *self, uint32_t left, uint32_t v)
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
static int WARN_UNUSED
msp_copy_overlap_count(msp_t *self, uint32_t k)
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

static int WARN_UNUSED
msp_store_edge(msp_t *self, double left, double right, node_id_t parent, node_id_t child)
{
    int ret = 0;
    edge_t *edge;
    const double *node_time = self->tables.nodes->time;

    assert(parent > child);
    assert(parent < (node_id_t) self->tables.nodes->num_rows);
    if (self->num_buffered_edges == self->max_buffered_edges - 1) {
        /* Grow the array */
        self->max_buffered_edges *= 2;
        edge = realloc(self->buffered_edges, self->max_buffered_edges * sizeof(edge_t));
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

static int WARN_UNUSED
msp_flush_edges(msp_t *self)
{
    int ret = 0;
    size_t j, num_edges;
    edge_t edge;

    if (self->num_buffered_edges > 0) {
        ret = squash_edges(self->buffered_edges, self->num_buffered_edges, &num_edges);
        if (ret != 0) {
            goto out;
        }
        for (j = 0; j < num_edges; j++) {
            edge = self->buffered_edges[j];
            ret = edge_table_add_row(self->tables.edges,
                    msp_genetic_to_phys(self, edge.left),
                    msp_genetic_to_phys(self, edge.right),
                    edge.parent, edge.child);
            if (ret < 0) {
                goto out;
            }

        }
        self->num_buffered_edges = 0;
    }
    ret = 0;
out:
    return ret;
}

static int WARN_UNUSED
msp_store_node(msp_t *self, uint32_t flags, double time, population_id_t population_id)
{
    int ret = 0;
    double scaled_time = self->model.model_time_to_generations(&self->model, time);

    ret = node_table_add_row(self->tables.nodes, flags, scaled_time, population_id,
            MSP_NULL_INDIVIDUAL, NULL, 0);
    if (ret < 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static int
msp_compress_overlap_counts(msp_t *self, uint32_t l, uint32_t r)
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

static int WARN_UNUSED
msp_conditional_compress_overlap_counts(msp_t *self, uint32_t l, uint32_t r)
{
    int ret = 0;
    double covered_fraction = (r - l) / (double) self->num_loci;

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
static int WARN_UNUSED
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
            fenwick_increment(&self->links, x->id, y->right - y->left);
            msp_free_segment(self, y);
        }
        y = x;
    }
    return 0;
}

static int WARN_UNUSED
msp_dtwf_recombine(msp_t *self, segment_t *x, segment_t **u, segment_t **v)
{
    int ret = 0;
    int ix;
    double mu;
    int64_t k;
    segment_t *y, *z;
    segment_t s1, s2;
    segment_t *seg_tails[] = {&s1, &s2};

    mu = 1.0 / self->recombination_rate;
    k = 1 + x->left + (int64_t) gsl_ran_exponential(self->rng, mu);
    assert(k > 0 && "Negative distance: Overflow in gsl_ran_exponential!");

    s1.next = NULL;
    s2.next = NULL;
    ix = (int) gsl_rng_uniform_int(self->rng, 2);
    seg_tails[ix]->next = x;
    x->prev = seg_tails[ix];

    while (x != NULL) {
        seg_tails[ix] = x;
        y = x->next;

        if (x->right > k) {
            // Make new segment
            assert(x->left <= k);
            self->num_re_events++;
            ix = (ix + 1) % 2;
            z = msp_alloc_segment(self, (uint32_t) k, x->right, x->value,
                    x->population_id, seg_tails[ix], x->next);
            assert(z->left < z->right);
            if (z == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            if (x->next != NULL) {
                x->next->prev = z;
            }
            seg_tails[ix]->next = z;
            seg_tails[ix] = z;
            x->next = NULL;
            x->right = (uint32_t) k;
            fenwick_increment(&self->links, x->id, k - z->right);
            assert(x->left < x->right);
            x = z;
            k = 1 + k + (int64_t) gsl_ran_exponential(self->rng, mu);
        }
        else if (x->right <= k && y != NULL && y->left >= k) {
            // Recombine in gap between segment and the next
            x->next = NULL;
            y->prev = NULL;
            while (y->left >= k) {
                self->num_re_events++;
                ix = (ix + 1) % 2;
                k = 1 + k + (int64_t) gsl_ran_exponential(self->rng, mu);
            }
            seg_tails[ix]->next = y;
            y->prev = seg_tails[ix];
            seg_tails[ix] = y;
            x = y;
            z = y;
            fenwick_set_value(&self->links, z->id, z->right - z->left - 1);

        } else {
            // Breakpoint in later segment
            x = y;
            z = y;
        }
    }
    // Remove sentinal segments
    *u = s1.next;
    *v = s2.next;

out:
    return ret;
}

static int WARN_UNUSED
msp_recombination_event(msp_t *self)
{
    int ret = 0;
    int64_t l, t, gap, k;
    size_t segment_id;
    node_mapping_t search;
    segment_t *x, *y, *z;
    int64_t num_links = fenwick_get_total(&self->links);

    self->num_re_events++;
    /* We can't use the GSL integer generator here as the range is too large */
    l = 1 + (int64_t) (gsl_rng_uniform(self->rng) * (double) num_links);
    //printf("l = %d, num_links = %d\n", (int)l, (int)num_links);
    assert(l > 0 && l <= num_links);
    segment_id = fenwick_find(&self->links, l);
    t = fenwick_get_cumulative_sum(&self->links, segment_id);
    gap = t - l;
    assert(gap >= 0 && gap < self->num_loci);
    y = msp_get_segment(self, segment_id);
    x = y->prev;
    k = y->right - gap - 1;
    assert(k >= 0 && k < self->num_loci);
    if (y->left < k) {
        z = msp_alloc_segment(self, (uint32_t) k, y->right, y->value,
                y->population_id, NULL, y->next);
        if (z == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        if (y->next != NULL) {
            y->next->prev = z;
        }
        y->next = NULL;
        y->right = (uint32_t) k;
        fenwick_increment(&self->links, y->id, k - z->right);
        search.left = (uint32_t) k;
        if (avl_search(&self->breakpoints, &search) == NULL) {
            ret = msp_insert_breakpoint(self, (uint32_t) k);
            if (ret != 0) {
                goto out;
            }
        } else {
            self->num_multiple_re_events++;
        }
    } else {
        assert(x != NULL);
        x->next = NULL;
        y->prev = NULL;
        z = y;
        self->num_trapped_re_events++;
    }
    fenwick_set_value(&self->links, z->id, z->right - z->left - 1);
    ret = msp_insert_individual(self, z);
out:
    return ret;
}

/* If we're implementing the SMC or SMC', discard this CA event if
 * there aren't any overlapping segments.
 */
static int WARN_UNUSED
msp_reject_ca_event(msp_t *self, segment_t *a, segment_t *b)
{
    int ret = 0;
    segment_t *x = a;
    segment_t *y = b;
    segment_t *beta;
    int64_t overlap, min_overlap;

    if (self->model.type == MSP_MODEL_SMC || self->model.type == MSP_MODEL_SMC_PRIME) {
        ret = 1;
        min_overlap = self->model.type == MSP_MODEL_SMC ? 1: 0;
        while (x != NULL && y != NULL) {
            if (y->left < x->left) {
                beta = x;
                x = y;
                y = beta;
            }
            overlap = ((int64_t) x->right) - ((int64_t) y->left);
            if (overlap >= min_overlap) {
                ret = 0;
                break;
            }
            x = x->next;
        }
    }
    return ret;
}

static int WARN_UNUSED
msp_merge_two_ancestors(msp_t *self, population_id_t population_id, segment_t *a,
        segment_t *b)
{
    int ret = 0;
    int coalescence = 0;
    int defrag_required = 0;
    node_id_t v;
    uint32_t l, r, l_min, r_max;
    avl_node_t *node;
    node_mapping_t *nm, search;
    segment_t *x, *y, *z, *alpha, *beta;

    x = a;
    y = b;
    /* Keep GCC happy */
    l_min = 0;
    r_max = 0;

    /* update num_links and get ready for loop */
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
                        x->population_id, NULL, NULL);
                if (alpha == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
                x->left = y->left;
            } else {
                l = x->left;
                r_max = GSL_MIN(x->right, y->right);
                if (!coalescence) {
                    coalescence = 1;
                    l_min = l;
                    ret = msp_flush_edges(self);
                    if (ret != 0) {
                        goto out;
                    }
                    ret = msp_store_node(self, 0, self->time, population_id);
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
                    alpha = msp_alloc_segment(self, l, r, v, population_id, NULL, NULL);
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
                fenwick_set_value(&self->links, alpha->id,
                        alpha->right - alpha->left - 1);
            } else {
                defrag_required |= z->right == alpha->left && z->value == alpha->value;
                z->next = alpha;
                fenwick_set_value(&self->links, alpha->id, alpha->right - z->right);
            }
            alpha->prev = z;
            z = alpha;
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

static int WARN_UNUSED
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

/* Merge the specified set of ancestors into a single ancestor. This is a
 * generalisation of the msp_common_ancestor_event method where we allow
 * any number of ancestors to merge. The AVL tree is a priority queue in
 * sorted by left coordinate.
 */
static int WARN_UNUSED
msp_merge_ancestors(msp_t *self, avl_tree_t *Q, population_id_t population_id)
{
    int ret = MSP_ERR_GENERIC;
    int coalescence = 0;
    int defrag_required = 0;
    node_id_t v;
    uint32_t j, l, r, h, r_max, next_l, l_min;
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
        r_max = self->num_loci;
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
                alpha = msp_alloc_segment(self, x->left, next_l, x->value,
                        x->population_id, NULL, NULL);
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
            if (!coalescence) {
                coalescence = 1;
                l_min = l;
                ret = msp_flush_edges(self);
                if (ret != 0) {
                    goto out;
                }
                ret = msp_store_node(self, 0, self->time, population_id);
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
                alpha = msp_alloc_segment(self, l, r, v, population_id,
                        NULL, NULL);
                if (alpha == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
            }
            /* Store the edges and update the priority queue */
            for (j = 0; j < h; j++) {
                x = H[j];
                ret = msp_store_edge(self, l, r, v, x->value);
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
                ret = msp_insert_individual(self, alpha);
                if (ret != 0) {
                    goto out;
                }
                fenwick_set_value(&self->links, alpha->id,
                        alpha->right - alpha->left - 1);
            } else {
                defrag_required |=
                    z->right == alpha->left && z->value == alpha->value;
                z->next = alpha;
                fenwick_set_value(&self->links, alpha->id,
                        alpha->right - z->right);
            }
            alpha->prev = z;
            z = alpha;
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

static int WARN_UNUSED
msp_migration_event(msp_t *self, population_id_t source_pop, population_id_t dest_pop)
{
    int ret = 0;
    uint32_t j;
    avl_node_t *node;
    avl_tree_t *source = &self->populations[source_pop].ancestors;
    size_t index = ((size_t) source_pop) * self->num_populations + (size_t) dest_pop;

    self->num_migration_events[index]++;
    j = (uint32_t) gsl_rng_uniform_int(self->rng, avl_count(source));
    node = avl_at(source, j);
    assert(node != NULL);
    ret = msp_move_individual(self, node, source, dest_pop);
    return ret;
}


static int WARN_UNUSED
msp_reset_memory_state(msp_t *self)
{
    int ret = 0;
    avl_node_t *node;
    node_mapping_t *nm;
    population_t *pop;
    segment_t *u, *v;
    size_t j;

    for (j = 0; j < self->num_populations; j++) {
        pop = &self->populations[j];
        for (node = pop->ancestors.head; node != NULL; node = node->next) {
            u = (segment_t *) node->item;
            while (u != NULL) {
                v = u->next;
                msp_free_segment(self, u);
                u = v;
            }
            avl_unlink_node(&pop->ancestors, node);
            msp_free_avl_node(self, node);
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

static int WARN_UNUSED
msp_insert_sample(msp_t *self, node_id_t sample, population_id_t population)
{
    int ret = MSP_ERR_GENERIC;
    segment_t *u;

    u = msp_alloc_segment(self, 0, self->num_loci, sample, population, NULL, NULL);
    if (u == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = msp_insert_individual(self, u);
    if (ret != 0) {
        goto out;
    }
    fenwick_set_value(&self->links, u->id, self->num_loci - 1);
out:
    return ret;
}

static inline int
msp_allocate_root_segments(msp_t *self, sparse_tree_t *tree,
        uint32_t left, uint32_t right,
        segment_t * restrict *root_segments_head,
        segment_t * restrict *root_segments_tail)
{
    int ret = 0;
    table_collection_t *tables = self->from_ts->tables;
    const population_id_t *node_population = tables->nodes->population;
    population_id_t population;
    node_id_t root;
    segment_t *seg, *tail;

    for (root = tree->left_root; root != MSP_NULL_NODE; root = tree->right_sib[root]) {
        population = node_population[root];
        /* Reference integrity has alreay been checked, but the null population
         * is still possibile */
        if (population == MSP_NULL_POPULATION) {
            ret = MSP_ERR_POPULATION_OUT_OF_BOUNDS;
            goto out;
        }
        if (root_segments_head[root] == NULL) {
            seg = msp_alloc_segment(self, left, right, root, population,
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
            } else {
                seg = msp_alloc_segment(self, left, right, root, population,
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
    sparse_tree_t t;
    node_id_t root;
    segment_t *seg;
    uint32_t left, right, num_roots, overlap, last_overlap;
    size_t num_nodes = self->tables.nodes->num_rows;
    segment_t **root_segments_head = calloc(num_nodes, sizeof(*root_segments_head));
    segment_t **root_segments_tail = calloc(num_nodes, sizeof(*root_segments_tail));

    ret = sparse_tree_alloc(&t, self->from_ts, 0);
    if (ret != 0) {
        goto out;
    }
    if (root_segments_head == NULL || root_segments_tail == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    /* Reset the tables to their correct position for replication */
    ret = table_collection_reset_position(&self->tables, &self->from_position);
    if (ret != 0) {
        goto out;
    }

    last_overlap = UINT32_MAX;
    for (t_iter = sparse_tree_first(&t); t_iter == 1; t_iter = sparse_tree_next(&t)) {
        ret = recomb_map_phys_to_discrete_genetic(self->recomb_map, t.left, &left);
        if (ret != 0) {
            goto out;
        }
        ret = recomb_map_phys_to_discrete_genetic(self->recomb_map, t.right, &right);
        if (ret != 0) {
            goto out;
        }
        if (left == right) {
            ret = MSP_ERR_RECOMB_MAP_TOO_COARSE;
            goto out;
        }
        assert(left < right);

        num_roots = (uint32_t) sparse_tree_get_num_roots(&t);
        overlap = 0;
        if (num_roots > 1) {
            overlap = num_roots;
            ret = msp_allocate_root_segments(self, &t,left, right,
                    root_segments_head, root_segments_tail);
            if (ret != 0) {
                goto out;
            }
        }
        if (overlap != last_overlap) {
            ret = msp_insert_overlap_count(self, left, overlap);
            if (ret != 0) {
                goto out;
            }
        }
    }
    if (t_iter != 0) {
        ret = t_iter;
        goto out;
    }

    ret = msp_insert_overlap_count(self, self->num_loci, UINT32_MAX);
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
            fenwick_set_value(&self->links, seg->id, seg->right - seg->left - 1);
            for (seg = seg->next; seg != NULL; seg = seg->next) {
                fenwick_set_value(&self->links, seg->id, seg->right - seg->prev->right);
            }
        }
    }

out:
    sparse_tree_free(&t);
    msp_safe_free(root_segments_head);
    msp_safe_free(root_segments_tail);
    return ret;
}

static int
msp_reset_from_samples(msp_t *self)
{
    int ret = 0;
    size_t j;
    node_id_t u;

    population_table_clear(self->tables.populations);
    edge_table_clear(self->tables.edges);
    node_table_clear(self->tables.nodes);
    migration_table_clear(self->tables.migrations);

    self->tables.sequence_length = self->recomb_map->sequence_length;
    for (j = 0; j < self->num_populations; j++) {
        ret = population_table_add_row(self->tables.populations, NULL, 0);
        if (ret < 0) {
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
        ret = msp_store_node(self, MSP_NODE_IS_SAMPLE,
                self->samples[u].time,
                self->samples[u].population_id);
        if (ret != 0) {
            goto out;
        }
    }
    ret = msp_insert_overlap_count(self, 0, self->num_samples);
    if (ret != 0) {
        goto out;
    }
    ret = msp_insert_overlap_count(self, self->num_loci, self->num_samples + 1);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

static int WARN_UNUSED
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

static int
msp_reopen_event_time_file(msp_t *self)
{
    int ret = MSP_ERR_IO;

    if (self->event_time_file != NULL) {
        if (fclose(self->event_time_file) != 0) {
            goto out;
        }
        self->event_time_file = NULL;
    }

    if (self->event_time_file_name != NULL) {
        self->event_time_file = fopen(self->event_time_file_name, "wb");
        if (self->event_time_file == NULL) {
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

int
msp_reset(msp_t *self)
{
    int ret = 0;
    size_t N = self->num_populations;
    population_id_t population_id;
    population_t *pop, *initial_pop;

    memcpy(&self->model, &self->initial_model, sizeof(self->model));
    ret = msp_reset_memory_state(self);
    if (ret != 0) {
        goto out;
    }
    /* Set up the initial segments and algorithm state */
    for (population_id = 0; population_id < (population_id_t) N; population_id++) {
        pop = self->populations + population_id;
        /* Set the initial population parameters */
        initial_pop = self->initial_populations + population_id;
        pop->growth_rate = initial_pop->growth_rate;
        pop->initial_size = initial_pop->initial_size;
        pop->start_time = self->time;
    }
    self->time = self->start_time;
    assert(self->time >= 0);
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
    self->num_ca_events = 0;
    self->num_rejected_ca_events = 0;
    self->num_trapped_re_events = 0;
    self->num_multiple_re_events = 0;

    ret = msp_reopen_event_time_file(self);
    if (ret != 0) {
        goto out;
    }
    memset(self->num_migration_events, 0, N * N * sizeof(size_t));
    self->state = MSP_STATE_INITIALISED;
out:
    return ret;
}

static int WARN_UNUSED
msp_initialise_from_samples(msp_t *self)
{
    int ret = 0;
    size_t j, k, initial_samples;
    population_id_t pop;

    if (self->start_time < 0) {
        self->start_time = 0;
    }
    if (self->num_samples < 2) {
        ret = MSP_ERR_INSUFFICIENT_SAMPLES;
        goto out;
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

static int WARN_UNUSED
msp_initialise_from_ts(msp_t *self)
{
    int ret = 0;
    double model_time, root_time;
    uint32_t num_samples;
    size_t num_nodes = self->tables.nodes->num_rows;
    population_id_t pop;
    size_t j;

    if (self->num_populations != self->tables.populations->num_rows) {
        ret = MSP_ERR_INCOMPATIBLE_FROM_TS;
        goto out;
    }
    /* Find the maximum time among the existing nodes */
    num_samples = 0;
    root_time = 0.0;
    for (j = 0; j < num_nodes; j++) {
        model_time = self->model.generations_to_model_time(
                &self->model, self->tables.nodes->time[j]);
        root_time = GSL_MAX(model_time, root_time);
        /* TODO we can catch ancient samples here and insert them as sampling
         * events, if we wish to support this. */
        if (self->tables.nodes->flags[j] & MSP_NODE_IS_SAMPLE) {
            num_samples++;
        }
        pop = self->tables.nodes->population[j];
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
int WARN_UNUSED
msp_initialise(msp_t *self)
{
    int ret = -1;

    /* These should really be proper checks with a return value */
    assert(self->num_loci >= 1);
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

static double
msp_get_population_size(msp_t *self, population_t *pop)
{
    double ret = 0;
    double alpha = pop->growth_rate;
    double t = self->time;
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
    }
    return ret;
}

static int WARN_UNUSED
msp_sanity_check(msp_t *self, int64_t num_links)
{
    int ret = 0;
    size_t n = msp_get_num_ancestors(self);

    /* TODO there are various worries when we're dealing with very
     * large population sizes and numbers of links. For example,
     * what if n is too big for the gsl integer random generator?
     * To at least determine if these are real things to worry
     * about, we put in some arbitrary (and hopefully safe) limits
     * and see if they are ever hit in practise.
     */
    if (n >= UINT32_MAX / 8) {
        ret = MSP_ERR_POPULATION_OVERFLOW;
        goto out;
    }
    if (num_links >= INT64_MAX / 8) {
        ret = MSP_ERR_LINKS_OVERFLOW;
        goto out;
    }
out:
    return ret;
}

static int
msp_record_event_time(msp_t *self, int32_t event)
{
    int ret = MSP_ERR_IO;
    double scaled_time;
    size_t written;

    if (self->event_time_file != NULL) {
        written = fwrite(&event, sizeof(event), 1, self->event_time_file);
        if (written != 1) {
            goto out;
        }
        scaled_time = self->model.model_time_to_generations(&self->model, self->time);
        written = fwrite(&scaled_time, sizeof(scaled_time), 1,
                self->event_time_file);
        if (written != 1) {
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

/* The main event loop for continuous time coalescent models. */
static int WARN_UNUSED
msp_run_coalescent(msp_t *self, double max_time, unsigned long max_events)
{
    int ret = 0;
    double lambda, t_temp, t_wait, ca_t_wait, re_t_wait, mig_t_wait,
           sampling_event_time, demographic_event_time;
    int64_t num_links;
    uint32_t j, k, n;
    population_id_t ca_pop_id, mig_source_pop, mig_dest_pop;
    unsigned long events = 0;
    sampling_event_t *se;

    while (msp_get_num_ancestors(self) > 0
            && self->time < max_time && events < max_events) {
        events++;
        num_links = fenwick_get_total(&self->links);
        ret = msp_sanity_check(self, num_links);
        if (ret != 0) {
            goto out;
        }
        /* Recombination */
        lambda = (double) num_links * self->recombination_rate;
        re_t_wait = DBL_MAX;
        if (lambda != 0.0) {
            re_t_wait = gsl_ran_exponential(self->rng, 1.0 / lambda);
        }
        /* Common ancestors */
        ca_t_wait = DBL_MAX;
        ca_pop_id = 0;
        for (j = 0; j < self->num_populations; j++) {
            t_temp = self->get_common_ancestor_waiting_time(self, (population_id_t) j);
            if (t_temp < 0) {
                /* Error condition signalled by negative waiting time */
                ret = (int) t_temp;
                goto out;
            }
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
            n = avl_count(&self->populations[j].ancestors);
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
        t_wait = GSL_MIN(GSL_MIN(re_t_wait, ca_t_wait), mig_t_wait);
        if (self->next_demographic_event == NULL
                && self->next_sampling_event == self->num_sampling_events
                && t_wait == DBL_MAX) {
            ret = MSP_ERR_INFINITE_WAITING_TIME;
            goto out;
        }
        t_temp = self->time + t_wait;
        /* printf("t = %f wait = %f\n", self->time, t_wait); */
        sampling_event_time = DBL_MAX;
        if (self->next_sampling_event < self->num_sampling_events) {
            sampling_event_time = self->sampling_events[
                self->next_sampling_event].time;
        }
        demographic_event_time = DBL_MAX;
        if (self->next_demographic_event != NULL) {
            demographic_event_time = self->next_demographic_event->time;
        }
        if (sampling_event_time < t_temp
                && sampling_event_time < demographic_event_time) {
            se = &self->sampling_events[self->next_sampling_event];
            self->time = se->time;
            ret = msp_insert_sample(self, se->sample, se->population_id);
            if (ret != 0) {
                goto out;
            }
            self->next_sampling_event++;
        } else if (demographic_event_time < t_temp) {
            ret = msp_apply_demographic_events(self);
            if (ret != 0) {
                goto out;
            }
        } else {
            self->time += t_wait;
            if (re_t_wait == t_wait) {
                ret = msp_record_event_time(self, 0);
                if (ret != 0) {
                    goto out;
                }
                ret = msp_recombination_event(self);
            } else if (ca_t_wait == t_wait) {
                ret = self->common_ancestor_event(self, ca_pop_id);
                if (ret == 1) {
                    /* The CA event has signalled that this event should be rejected */
                    self->time -= t_wait;
                    ret = 0;
                } else {
                    ret = msp_record_event_time(self, 1);
                    if (ret != 0) {
                        goto out;
                    }
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

/* List structure for collecting segments by parent */
typedef struct _segment_list_t {
    avl_node_t *node;
    struct _segment_list_t *next;
} segment_list_t;


/* Performs a single generation under the Wright Fisher model */
static int WARN_UNUSED
msp_dtwf_generation(msp_t *self)
{
    int ret = 0;
    int ix;
    uint32_t N, i, j, k, p;
    int64_t num_links;
    size_t segment_mem_offset;
    population_t *pop;
    segment_t *x;
    segment_t *u[2];
    segment_list_t **parents = NULL;
    segment_list_t *segment_mem = NULL;
    segment_list_t *s;
    avl_node_t *a, *node;
    avl_tree_t Q[2];

    for (i = 0; i < 2; i++) {
        avl_init_tree(&Q[i], cmp_segment_queue, NULL);
    }

    num_links = fenwick_get_total(&self->links);
    ret = msp_sanity_check(self, num_links);
    if (ret != 0) {
        goto out;
    }
    msp_verify_segments(self);

    for (j = 0; j < self->num_populations; j++) {

        pop = &self->populations[j];
        if (avl_count(&pop->ancestors) == 0) {
            continue;
        }
        /* For the DTWF, N for each population is the reference population size
         * from the model multiplied by the current population size, rounded
         * to the nearest integer. Thus, the population's size is always relative
         * to the reference model population size (which is also true for the
         * coalescent models. */
        N = (uint32_t) round(
            msp_get_population_size(self, pop) * self->model.population_size);
        /* printf("%u parents to choose from %f \n", N, msp_get_population_size(self, pop)); */
        if (N == 0) {
            ret = MSP_ERR_INFINITE_WAITING_TIME;
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
        for (a = pop->ancestors.head; a != NULL; a = a->next) {
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

        msp_verify_segments(self);
        // Iterate through offspring of parent k, adding to avl_tree
        for (k = 0; k < N; k++) {
            for (s = parents[k]; s != NULL; s = s->next) {
                node = s->node;
                x = (segment_t *) node->item;
                avl_unlink_node(&pop->ancestors, node);
                msp_free_avl_node(self, node);

                // Recombine ancestor

                if (self->recombination_rate > 0) {
                    ret = msp_dtwf_recombine(self, x, &u[0], &u[1]);
                    if (ret != 0) {
                        goto out;
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
                ret = msp_merge_ancestors(self, &Q[i], (population_id_t) j);
                if (ret != 0) {
                    goto out;
                }
            }
        }
        msp_verify_segments(self);
        free(parents);
        free(segment_mem);
        segment_mem = NULL;
        parents = NULL;
    }
out:
    if (parents != NULL) {
        free(parents);
    }
    if (segment_mem != NULL){
        free(segment_mem);
    }
    return ret;
}

/* The main event loop for the Wright Fisher model. */
static int WARN_UNUSED
msp_run_dtwf(msp_t *self, double max_time, unsigned long max_events)
{
    int ret = 0;
    unsigned long events = 0;
    int mig_source_pop, mig_dest_pop;
    sampling_event_t *se;
    double mu;
    uint32_t j, k, i, num_migrations, n;

    while (msp_get_num_ancestors(self) > 0
            && self->time < max_time && events < max_events) {
        events++;
        self->time++;
        ret = msp_dtwf_generation(self);
        if (ret != 0) {
            goto out;
        }

        /* TODO Need to reason more carefully about what order these events happen in!
         * Do migrations come before or after the actual breeding? Given a sampling
         * event at time t, does this happen before or after the breeding for time
         * t? Same for all demographic events */
        mig_source_pop = 0;
        mig_dest_pop = 0;
        for (j = 0; j < self->num_populations; j++) {
            for (k = 0; k < self->num_populations; k++) {
                n = avl_count(&self->populations[j].ancestors);
                mu = n * self->migration_matrix[j * self->num_populations + k];
                if (mu == 0) {
                    continue;
                }
                /* TODO this is brittle here, as it's easy to overflow gsl_ran_poisson.
                 * Also, is this the correct model? */
                num_migrations = GSL_MAX(n, gsl_ran_poisson(self->rng, mu));
                for (i = 0; i < num_migrations; i++) {
                    /* m[j, k] is the rate at which migrants move from
                     * population k to j forwards in time. Backwards
                     * in time, we move the individual from from
                     * population j into population k.
                     */
                    mig_source_pop = (population_id_t) j;
                    mig_dest_pop = (population_id_t) k;
                    ret = msp_migration_event(self, mig_source_pop, mig_dest_pop);
                    if (ret != 0) {
                        goto out;
                    }
                }
            }
        }
        if (self->next_sampling_event < self->num_sampling_events) {
            if (self->sampling_events[self->next_sampling_event].time >= self->time) {
                se = self->sampling_events + self->next_sampling_event;
                ret = msp_insert_sample(self, se->sample, se->population_id);
                if (ret != 0) {
                    goto out;
                }
                self->next_sampling_event++;
            }
        }
        if (self->next_demographic_event != NULL) {
            if (self->next_demographic_event->time >= self->time) {
                ret = msp_apply_demographic_events(self);
                if (ret != 0) {
                    goto out;
                }
            }
        }
    }
out:
    return ret;
}

/* Runs the simulation backwards in time until either the sample has coalesced,
 * or specified maximum simulation time has been reached or the specified maximum
 * number of events has been reached.
 */
int WARN_UNUSED
msp_run(msp_t *self, double max_time, unsigned long max_events)
{
    int ret = 0;
    simulation_model_t *model = &self->model;
    double scaled_time = model->generations_to_model_time(model, max_time);
    int err;

    /* printf("Running until %f generations = %f scaled time\n", max_time, */
    /*         scaled_time); */

    if (self->state == MSP_STATE_INITIALISED) {
        self->state = MSP_STATE_SIMULATING;
    }
    if (self->state != MSP_STATE_SIMULATING) {
        ret = MSP_ERR_BAD_STATE;
        goto out;
    }
    if (self->model.type == MSP_MODEL_DTWF) {
        ret = msp_run_dtwf(self, scaled_time, max_events);
    } else {
        ret = msp_run_coalescent(self, scaled_time, max_events);
    }
    if (ret != 0) {
        goto out;
    }
    ret = msp_flush_edges(self);
    if (ret != 0) {
        goto out;
    }
    if (msp_get_num_ancestors(self) != 0) {
        ret = 1;
        if (self->time >= scaled_time) {
            ret = 2;
        }
    }
    /* Flush the event time file */
    if (self->event_time_file != NULL) {
        err = fflush(self->event_time_file);
        if (err != 0) {
            ret = MSP_ERR_IO;
            goto out;
        }
    }
out:
    return ret;
}

/* Add in nodes and edges for the remaining segments to the output table. */
static int WARN_UNUSED
msp_insert_uncoalesced_edges(msp_t *self, table_collection_t *tables)
{
    int ret = 0;
    population_id_t pop;
    avl_node_t *a;
    segment_t *seg;
    node_id_t node;
    int64_t edge_start;
    node_table_t *nodes = tables->nodes;
    const double current_time = self->model.model_time_to_generations(&self->model,
            self->time);

    for (pop = 0; pop < (population_id_t) self->num_populations; pop++) {
        for (a = self->populations[pop].ancestors.head; a != NULL; a = a->next) {
            /* If there are any nodes in the segment chain with the current time,
             * then we don't make any unary edges for them. This is because (a)
             * we'd end up edges with the same parent and child time (if we didn't
             * hack an extra epsilon onto the parent time) and (b), this node
             * could only have arisen as the result of a coalescence and so this
             * node really does represent the current ancestor */
            node = MSP_NULL_NODE;
            for (seg = (segment_t *) a->item; seg != NULL; seg = seg->next) {
                if (nodes->time[seg->value] == current_time) {
                    node = seg->value;
                    break;
                }
            }
            if (node == MSP_NULL_NODE) {
                /* Add a node for this ancestor */
                node = node_table_add_row(nodes, 0, current_time, pop,
                        MSP_NULL_INDIVIDUAL, NULL, 0);
                if (node < 0) {
                    ret = node;
                    goto out;
                }
            }

            /* For every segment add an edge pointing to this new node */
            for (seg = (segment_t *) a->item; seg != NULL; seg = seg->next) {
                if (seg->value != node) {
                    assert(nodes->time[node] > nodes->time[seg->value]);
                    ret = edge_table_add_row(tables->edges,
                        msp_genetic_to_phys(self, seg->left),
                        msp_genetic_to_phys(self, seg->right),
                        node, seg->value);
                    if (ret < 0) {
                        goto out;
                    }
                }
            }
        }
    }

    /* Find the first edge with parent == current time */
    edge_start = ((int64_t) tables->edges->num_rows) - 1;
    while (edge_start >= 0
            && nodes->time[tables->edges->parent[edge_start]] == current_time) {
        edge_start--;
    }
    /* TODO This could be done more efficiently probably, but at least we only
     * end up sorting a handful of edges at the end, so it's probably not so
     * bad. */
    ret = table_collection_sort(tables, (size_t) (edge_start + 1), 0);
out:
    return ret;
}

int WARN_UNUSED
msp_populate_tables(msp_t *self, table_collection_t *tables)
{
    int ret = 0;

    ret = table_collection_copy(&self->tables, tables);
    if (ret != 0) {
        goto out;
    }
    if (!msp_is_completed(self)) {
        ret = msp_insert_uncoalesced_edges(self, tables);
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

    if (self->state == MSP_STATE_INITIALISED) {
        self->state = MSP_STATE_DEBUGGING;
        first_call = 1;
    }
    if (self->state != MSP_STATE_DEBUGGING) {
        ret = MSP_ERR_BAD_STATE;
        goto out;
    }
    if (! first_call && self->next_demographic_event != NULL) {
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
int WARN_UNUSED
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
        *pop_size = model->population_size * pop->initial_size;
    } else {
        dt = model->generations_to_model_time(model, time) - pop->start_time;
        *pop_size = model->population_size * pop->initial_size
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

size_t
msp_get_num_loci(msp_t *self)
{
    return (size_t) self->num_loci;
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
msp_get_num_ancestors(msp_t *self)
{
    size_t n = 0;
    size_t j;

    for (j = 0; j < self->num_populations; j++) {
        n += avl_count(&self->populations[j].ancestors);
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
    return (size_t) self->tables.nodes->num_rows;
}

size_t
msp_get_num_edges(msp_t *self)
{
    return (size_t) self->tables.edges->num_rows;
}

size_t
msp_get_num_migrations(msp_t *self)
{
    return (size_t) self->tables.migrations->num_rows;
}


int WARN_UNUSED
msp_get_ancestors(msp_t *self, segment_t **ancestors)
{
    int ret = -1;
    avl_node_t *node;
    avl_tree_t *population_ancestors;
    size_t j;
    size_t k = 0;

    for (j = 0; j < self->num_populations; j++) {
        population_ancestors = &self->populations[j].ancestors;
        for (node = population_ancestors->head; node != NULL; node = node->next) {
            ancestors[k] = (segment_t *) node->item;
            k++;
        }
    }
    ret = 0;
    return ret;
}

int WARN_UNUSED
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

int WARN_UNUSED
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

int WARN_UNUSED
msp_get_num_migration_events(msp_t *self, size_t *num_migration_events)
{
    size_t N = self->num_populations;

    memcpy(num_migration_events, self->num_migration_events, N * N * sizeof(size_t));
    return 0;
}

int WARN_UNUSED
msp_get_samples(msp_t *self, sample_t **samples)
{
    *samples = self->samples;
    return 0;
}

int WARN_UNUSED
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
    *initial_size = model->population_size * pop->initial_size;
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
msp_get_recombination_rate(msp_t *self)
{
    simulation_model_t *model = &self->model;
    return model->model_rate_to_generation_rate(model, self->recombination_rate);
}

/* Demographic events. All times and input parameters are specified in units
 * of generations. When we store these values, we must rescale them into
 * model time, as appropriate. */
static int WARN_UNUSED
msp_add_demographic_event(msp_t *self, double time, demographic_event_t **event)
{
    int ret = MSP_ERR_GENERIC;
    demographic_event_t *ret_event = NULL;

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
        pop->initial_size = initial_size / model->population_size;
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

static int WARN_UNUSED
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
int WARN_UNUSED
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

    /* This should have been caught on adding the event */
    if (source < 0 || source > N || dest < 0 || dest > N) {
        ret = MSP_ERR_ASSERTION_FAILED;
        goto out;
    }
    /*
     * Move lineages from source to dest with probability p.
     */
    pop = &self->populations[source].ancestors;
    node = pop->head;
    while (node != NULL) {
        next = node->next;
        if (gsl_rng_uniform(self->rng) < p) {
            ret = msp_move_individual(self, node, pop, dest);
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
int WARN_UNUSED
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
        ret = MSP_ERR_BAD_PARAM_VALUE;
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

    /* This should have been caught on adding the event */
    if (population_id < 0 || population_id > N) {
        ret = MSP_ERR_ASSERTION_FAILED;
        goto out;
    }
    avl_init_tree(&Q, cmp_segment_queue, NULL);
    /*
     * Find the individuals that descend from the common ancestor
     * during this simple_bottleneck.
     */
    pop = &self->populations[population_id].ancestors;
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
    ret = msp_merge_ancestors(self, &Q, population_id);
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
int WARN_UNUSED
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
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (self->model.type != MSP_MODEL_HUDSON) {
        ret = MSP_ERR_BAD_MODEL;
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

    /* This should have been caught on adding the event */
    if (population_id < 0 || population_id >= N) {
        ret = MSP_ERR_ASSERTION_FAILED;
        goto out;
    }
    pop = &self->populations[population_id].ancestors;
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
        pi[u] = MSP_NULL_NODE;
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
        while (pi[u] != MSP_NULL_NODE) {
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
            ret = msp_merge_ancestors(self, &sets[lineages[j]], population_id);
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
int WARN_UNUSED
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
    if (self->model.type != MSP_MODEL_HUDSON) {
        ret = MSP_ERR_BAD_MODEL;
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
    return 4 * model->population_size * t;
}

static double
std_generations_to_model_time(simulation_model_t *model, double g)
{
    return g / (4 * model->population_size);
}

static double
std_generation_rate_to_model_rate(simulation_model_t *model, double rate)
{
    return rate * 4 * model->population_size;
}

static double
std_model_rate_to_generation_rate(simulation_model_t *model, double rate)
{
    return rate / (4 * model->population_size);
}

static double
msp_std_get_common_ancestor_waiting_time(msp_t *self, population_id_t pop_id)
{
    population_t *pop = &self->populations[pop_id];
    double n = (double) avl_count(&pop->ancestors);
    double lambda = n * (n - 1.0);

    return msp_get_common_ancestor_waiting_time_from_rate(self, pop, lambda);
}

static int WARN_UNUSED
msp_std_common_ancestor_event(msp_t *self, population_id_t population_id)
{
    int ret = 0;
    uint32_t j, n;
    avl_tree_t *ancestors;
    avl_node_t *x_node, *y_node, *node;
    segment_t *x, *y;

    ancestors = &self->populations[population_id].ancestors;
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
        ret = msp_merge_two_ancestors(self, population_id, x, y);
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
    return 4 * model->population_size * t;
}

static double
dirac_generations_to_model_time(simulation_model_t *model, double g)
{
    return g /(4 * model->population_size);
}

static double
dirac_generation_rate_to_model_rate(simulation_model_t *model, double rate)
{
    return rate * 4 * model->population_size;
}

static double
dirac_model_rate_to_generation_rate(simulation_model_t *model, double rate)
{
    // This works for comparing dirac kingman case ... but why 4Ne^2
    return rate / ( 4 * gsl_pow_2(model->population_size));
}

double
compute_falling_factorial_log(unsigned int m)
{
    // (n)_m = n(n-1)(n-2)... (n-m + 1)  and  (n)_0 := 1
    // n = 4
    unsigned int l = 0;
    double ret = 1.0;
    while (l < m){
        l++;
        ret *= (4.0 - l + 1.0);
    }
    return gsl_sf_log(ret);
}

/* This calculates the rate given by Eq (2) in the notes
 */
double
compute_dirac_coalescence_rate(unsigned int num_ancestors, double psi, double c)
{
    unsigned int l, m;
    double r[5], r_max;
    double ret = 0;
    double b = num_ancestors;

    assert(b > 0);
    assert(psi > 0);
    assert(psi < 1);
    assert(c >= 0);

    m = GSL_MIN(num_ancestors, 4);
    /* An underflow error occurs because of the large exponent (b-l). We use the
     * LSE trick to approximate this calculation. For details, see at
     * https://en.wikipedia.org/wiki/LogSumExp
     */
    r[0] = b * gsl_sf_log(1 - psi);
    r_max = r[0];
    for (l = 1; l <= m; l++){
        r[l] = gsl_sf_lnchoose(num_ancestors, l)
                + compute_falling_factorial_log(l)
                + l * gsl_sf_log(psi / 4.0)
                + (b - l) * gsl_sf_log(1 - psi);
        r_max = GSL_MAX(r_max, r[l]);
    }

    for (l = 0; l <= m; l++)  {
        ret += exp(r[l] - r_max);
    }
    ret = exp(r_max + log(ret));

    ret = 1 - ret;
    ret *= c;
    ret += b * (b - 1) / 2;

    return ret;
}

static double
msp_dirac_compute_coalescence_rate(msp_t *self, unsigned int num_ancestors)
{
    double psi = self->model.params.dirac_coalescent.psi;
    double c = self->model.params.dirac_coalescent.c;

    return compute_dirac_coalescence_rate(num_ancestors, psi, c);
}

static double
msp_dirac_get_common_ancestor_waiting_time(msp_t *self, population_id_t pop_id)
{
    population_t *pop = &self->populations[pop_id];
    unsigned int n = (unsigned int) avl_count(&pop->ancestors);
    double lambda = msp_dirac_compute_coalescence_rate(self, n) * 2;

    return msp_get_common_ancestor_waiting_time_from_rate(self, pop, lambda);
}

static int WARN_UNUSED
msp_dirac_common_ancestor_event(msp_t *self, population_id_t pop_id)
{
    int ret = 0;
    uint32_t j, n, max_pot_size;
    const uint32_t num_pots = 4;
    avl_tree_t *ancestors, Q[4]; /* MSVC won't let us use num_pots here */
    avl_node_t *x_node, *y_node, *node, *next, *q_node;
    segment_t *x, *y, *u;

    ancestors = &self->populations[pop_id].ancestors;
    if (gsl_rng_uniform(self->rng) < (1 / (1.0 + self->model.params.dirac_coalescent.c))) {
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
        ret = msp_merge_two_ancestors(self, 0, x, y);
    } else {
        /* In the multiple merger regime we have four different 'pots' that
         * lineages get assigned to, where all lineages in a given pot are merged into
         * a common ancestor.
         */
        for (j = 0; j < num_pots; j++){
            avl_init_tree(&Q[j], cmp_segment_queue, NULL);
        }
        node = ancestors->head;
        while (node != NULL) {
            next = node->next;
            /* With probability psi / 4, a given lineage participates in this event. */
            if (gsl_rng_uniform(self->rng) < self->model.params.dirac_coalescent.psi / 4.0) {
                u = (segment_t *) node->item;
                avl_unlink_node(ancestors, node);
                msp_free_avl_node(self, node);
                q_node = msp_alloc_avl_node(self);
                if (q_node == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
                avl_init_node(q_node, u);
                /* Now assign this ancestor to a uniformly chosen pot */
                j = (uint32_t) gsl_rng_uniform_int(self->rng, num_pots);
                q_node = avl_insert_node(&Q[j], q_node);
                assert(q_node != NULL);
            }
            node = next;
        }
        /* All the lineages that have been assigned to the particular pots can now be
         * merged.
         */
        max_pot_size = 0;
        for (j = 0; j < num_pots; j++){
            max_pot_size = GSL_MAX(max_pot_size, avl_count(&Q[j]));
            ret = msp_merge_ancestors(self, &Q[j], 0);
            if (ret < 0) {
                goto out;
            }
        }
        /* If no coalescence has occured, we need to signal this back to the calling
         * function so that the event can be 'cancelled'. This is done by returning 1.
         */
        if (max_pot_size < 2){
            ret = 1;
        }
    }
out:
    return ret;
}

/**************************************************************
 * Beta coalescent
 **************************************************************/

static double
beta_model_time_to_generations(simulation_model_t *model, double t)
{
    return 4 * model->population_size * t;
}

static double
beta_generations_to_model_time(simulation_model_t *model, double g)
{
    return g / (4 * model->population_size);
}

static double
beta_generation_rate_to_model_rate(simulation_model_t *model, double rate)
{
    return rate * 4 * model->population_size;
}

static double
beta_model_rate_to_generation_rate(simulation_model_t *model, double rate)
{
    double alpha = model->params.beta_coalescent.alpha;
    double m = model->params.beta_coalescent.m;
    double phi = model->params.beta_coalescent.phi;
    double scalar = exp(log(alpha) - alpha * log(m)
            + gsl_sf_beta_inc (2 - alpha, alpha, phi)
            - (alpha - 1) * model->population_size);
    return rate / ( 4 * scalar );
}

static void
beta_model_free(simulation_model_t *model)
{
    if (model->params.beta_coalescent.integration_workspace != NULL) {
        gsl_integration_workspace_free(
            model->params.beta_coalescent.integration_workspace);
        model->params.beta_coalescent.integration_workspace = NULL;
    }
}

struct beta_integral_params {
    unsigned int num_ancestors;
    double alpha;
};

static double
beta_integrand(double x, void *p)
{
    struct beta_integral_params *params = (struct beta_integral_params *) p;
    unsigned int num_ancestors = params->num_ancestors;
    double alpha = params->alpha;
    unsigned int l, m;
    double r[5], r_max;
    double ret = 0;
    double b = num_ancestors;
    double log_x = gsl_sf_log(x);
    double log_1_minus_x = gsl_sf_log(1 - x);

    m = GSL_MIN(num_ancestors, 4);
    /* An underflow error occurs because of the large exponent (b-l). We use the
     * LSE trick to approximate this calculation. For details, see at
     * https://en.wikipedia.org/wiki/LogSumExp
     */
    r[0] = b * gsl_sf_log(1 - x);
    r_max = r[0];
    for (l = 1; l <= m; l++){
        r[l] = gsl_sf_lnchoose(num_ancestors, l)
                + compute_falling_factorial_log(l)
                - l * 1.386294      /* log(4) = 1.386294 */
                + l * log_x
                + (b - l) * log_1_minus_x;
        r_max = GSL_MAX(r_max, r[l]);
    }

    for (l = 0; l <= m; l++)  {
        ret += exp(r[l] - r_max);
    }
    ret = 1 - exp(r_max + log(ret));
    ret *= exp((-1 - alpha) * log_x + (alpha - 1) * log_1_minus_x);
    ret *= 4 / gsl_sf_beta(2 - alpha, alpha);
    return ret;
}

int
msp_compute_beta_integral(msp_t *self, unsigned int num_ancestors, double alpha, double *result)
{
    int ret = 0;
    double err;
    gsl_function F;
    gsl_integration_workspace *w = self->model.params.beta_coalescent.integration_workspace;
    size_t workspace_size = self->model.params.beta_coalescent.integration_workspace_size;
    double epsrel = self->model.params.beta_coalescent.integration_epsrel;
    double epsabs = self->model.params.beta_coalescent.integration_epsabs;
    struct beta_integral_params params = {num_ancestors, alpha};

    F.function = &beta_integrand;
    F.params = &params;
    ret = gsl_integration_qags(&F, 0, 1, epsabs, epsrel, workspace_size, w, result, &err);
    if (ret != 0) {
        /* It's ugly, but it should only happen while tuning models so we output
         * the GSL error to stderr */
        fprintf(stderr,  "GSL error: %s\n", gsl_strerror(ret));
        ret = MSP_ERR_INTEGRATION_FAILED;
    }
    return ret;
}

int
msp_beta_compute_coalescence_rate(msp_t *self, unsigned int num_ancestors, double *result)
{
    int ret = 0;
    double alpha = self->model.params.beta_coalescent.alpha;

    *result = 1;
    if (self->model.params.beta_coalescent.truncation_point != 0) {
        if (num_ancestors > 2) {
            ret = msp_compute_beta_integral(self, num_ancestors, alpha, result);
        }
    }
    return ret;
}

static double
msp_beta_get_common_ancestor_waiting_time(msp_t *self, population_id_t pop_id)
{
    int ret = 0;
    double result, lambda;
    population_t *pop = &self->populations[pop_id];
    unsigned int n = (unsigned int) avl_count(&pop->ancestors);

    ret = msp_beta_compute_coalescence_rate(self, n, &lambda);
    if (ret != 0) {
        assert(ret < 0);
        /* An error occured, and we signal this back using a negative waiting time */
        result = ret;
    } else {
        lambda *= 2; /* JK: Why do we multiply lambda by 2 here? */
        result = msp_get_common_ancestor_waiting_time_from_rate(self, pop, lambda);
    }
    return result;
}

static int WARN_UNUSED
msp_beta_common_ancestor_event(msp_t *self, population_id_t pop_id)
{
    int ret = 0;
    uint32_t j, n, max_pot_size;
    const uint32_t num_pots = 4;
    avl_tree_t *ancestors, Q[5]; /* MSVC won't let us use num_pots here */
    avl_node_t *x_node, *y_node, *node, *next, *q_node;
    segment_t *x, *y, *u;

    double beta_x;

    for (j = 0; j < 5; j++){
        avl_init_tree(&Q[j], cmp_segment_queue, NULL);
    }

    ancestors = &self->populations[pop_id].ancestors;
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

    q_node = msp_alloc_avl_node(self);
    if (q_node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    avl_init_node(q_node, x);
    q_node = avl_insert_node(&Q[4], q_node);

    q_node = msp_alloc_avl_node(self);
    if (q_node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    avl_init_node(q_node, y);
    q_node = avl_insert_node(&Q[4], q_node);
    beta_x = gsl_ran_beta(self->rng,
                 2.0 - self->model.params.beta_coalescent.alpha,
                 self->model.params.beta_coalescent.alpha);

    /* In the multiple merger regime we have four different 'pots' that
     * lineages get assigned to, where all lineages in a given pot are merged into
     * a common ancestor.
     */
    node = ancestors->head;
    while (node != NULL) {
        next = node->next;
        /* With probability beta_x / 4, a given lineage participates in this event. */
        if (gsl_rng_uniform(self->rng) < beta_x / 4.0) {
            u = (segment_t *) node->item;
            avl_unlink_node(ancestors, node);
            msp_free_avl_node(self, node);
            q_node = msp_alloc_avl_node(self);
            if (q_node == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            avl_init_node(q_node, u);
            /* Now assign this ancestor to a uniformly chosen pot */
            j = (uint32_t) gsl_rng_uniform_int(self->rng, num_pots);
            q_node = avl_insert_node(&Q[j], q_node);
            assert(q_node != NULL);
        }
        node = next;
    }
    /* All the lineages that have been assigned to the particular pots can now be
     * merged.
     */
    max_pot_size = 0;
    for (j = 0; j < 5; j++){
        max_pot_size = GSL_MAX(max_pot_size, avl_count(&Q[j]));
        ret = msp_merge_ancestors(self, &Q[j], 0);
        if (ret < 0) {
            goto out;
        }
    }
    /* If no coalescence has occured, we need to signal this back to the calling
     * function so that the event can be 'cancelled'. This is done by returning 1.
     */
    if (max_pot_size < 2){
        ret = 1;
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
 * Public API for setting simulation models.
 **************************************************************/

/* Unscale all times and rates from the current model time to generations. */
static int
msp_unscale_model_times(msp_t *self)
{
    uint32_t j;
    simulation_model_t *model = &self->model;
    demographic_event_t *de;

    self->start_time = self->model.model_time_to_generations(model, self->start_time);
    self->time = self->model.model_time_to_generations(model, self->time);
    self->recombination_rate = self->model.model_rate_to_generation_rate(
            model, self->recombination_rate);
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
    self->recombination_rate = model->generation_rate_to_model_rate(
            model, self->recombination_rate);
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
msp_set_simulation_model(msp_t *self, int model, double population_size)
{
    int ret = 0;

    if (model != MSP_MODEL_HUDSON && model != MSP_MODEL_SMC
            && model != MSP_MODEL_SMC_PRIME
            && model != MSP_MODEL_DIRAC
            && model != MSP_MODEL_BETA
            && model != MSP_MODEL_DTWF) {
        ret = MSP_ERR_BAD_MODEL;
        goto out;
    }
    if (population_size <= 0) {
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
    }
    self->model.type = model;
    self->model.population_size = population_size;
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
msp_set_simulation_model_hudson(msp_t *self, double population_size)
{
    int ret = msp_set_simulation_model(self, MSP_MODEL_HUDSON, population_size);
    if (ret != 0) {
        goto out;
    }
    ret = msp_rescale_model_times(self);
out:
    return ret;
}

int
msp_set_simulation_model_smc(msp_t *self, double population_size)
{
    int ret =  msp_set_simulation_model(self, MSP_MODEL_SMC, population_size);
    if (ret != 0) {
        goto out;
    }
    ret = msp_rescale_model_times(self);
out:
    return ret;
}

int
msp_set_simulation_model_smc_prime(msp_t *self, double population_size)
{
    int ret =  msp_set_simulation_model(self, MSP_MODEL_SMC_PRIME, population_size);
    if (ret != 0) {
        goto out;
    }
    ret = msp_rescale_model_times(self);
out:
    return ret;
}

int
msp_set_simulation_model_dtwf(msp_t *self, double population_size)
{
    int ret = 0;
    ret = msp_set_simulation_model(self, MSP_MODEL_DTWF, population_size);
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
msp_set_simulation_model_dirac(msp_t *self, double population_size, double psi, double c)
{
    int ret = 0;

    /* TODO bounds check psi: what are legal values? */
    ret = msp_set_simulation_model(self, MSP_MODEL_DIRAC, population_size);
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
msp_set_simulation_model_beta(msp_t *self, double population_size, double alpha,
        double truncation_point)
{

    int ret = 0;

    /* TODO bounds check alpha and truncation_point: what are legal values? */
    ret = msp_set_simulation_model(self, MSP_MODEL_BETA, population_size);
    if (ret != 0) {
        goto out;
    }
    self->model.params.beta_coalescent.alpha = alpha;
    self->model.params.beta_coalescent.truncation_point = truncation_point;
    self->model.params.beta_coalescent.m =
        2.0 + exp(alpha * 0.6931472 + (1 - alpha) * 1.098612 - log(alpha - 1));
    //self->model.params.beta_coalescent.phi = beta_compute_phi(population_size,
                        //truncation_point, self->model.params.beta_coalescent.m);
    //self->model.params.beta_coalescent.K = truncation_point;

    /* TODO we probably want to make these input parameters, as there will be situations
     * where integration fails and being able to tune them will be useful */
    self->model.params.beta_coalescent.integration_epsrel = 1e-3;
    self->model.params.beta_coalescent.integration_epsabs = 0;
    /* TODO Is 1000 a good size for the workspace here? Should it be a parameter? */
    self->model.params.beta_coalescent.integration_workspace_size = 1000;
    self->model.params.beta_coalescent.integration_workspace =
        gsl_integration_workspace_alloc(
            self->model.params.beta_coalescent.integration_workspace_size);
    if (self->model.params.beta_coalescent.integration_workspace == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    self->model.model_time_to_generations = beta_model_time_to_generations;
    self->model.generations_to_model_time = beta_generations_to_model_time;
    self->model.generation_rate_to_model_rate = beta_generation_rate_to_model_rate;
    self->model.model_rate_to_generation_rate = beta_model_rate_to_generation_rate;
    self->model.free = beta_model_free;
    self->get_common_ancestor_waiting_time = msp_beta_get_common_ancestor_waiting_time;
    self->common_ancestor_event = msp_beta_common_ancestor_event;
    ret = msp_rescale_model_times(self);
out:
    return ret;
}
