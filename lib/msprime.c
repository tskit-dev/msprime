/*
** Copyright (C) 2015-2016 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

#include <hdf5.h>

#include "err.h"
#include "avl.h"
#include "object_heap.h"
#include "fenwick.h"
#include "msprime.h"

#define MSP_HDF5_ERR_MSG_SIZE 1024

/* State machine for the simulator object. */
#define MSP_STATE_NEW 0
#define MSP_STATE_INITIALISED 1
#define MSP_STATE_SIMULATING 2
#define MSP_STATE_DEBUGGING 3

static char _hdf5_error[MSP_HDF5_ERR_MSG_SIZE];

static herr_t
hdf5_error_walker(unsigned n, const H5E_error2_t *err_desc, void *client_data)
{
    /* We only copy the message from the first element in the error stack */
    if (_hdf5_error[0] == '\0') {
        snprintf(_hdf5_error, MSP_HDF5_ERR_MSG_SIZE,
                "HDF5 Error: %d: %d:'%s'",
                (int) err_desc->maj_num, (int) err_desc->min_num,
                err_desc->desc);
    }
    return 0;
}

const char *
msp_strerror(int err)
{
    const char *ret = "Unknown error";

    if (err == MSP_ERR_NO_MEMORY) {
        ret = "Out of memory.";
    } else if (err == MSP_ERR_GENERIC) {
        ret = "Generic error; please file a bug report";
    } else if (err == MSP_ERR_FILE_FORMAT) {
        ret = "File format error";
    } else if (err == MSP_ERR_BAD_STATE) {
        ret = "Bad simulator state. initialise or reset must be called.";
    } else if (err == MSP_ERR_NEWICK_OVERFLOW) {
        ret = "Newick string generation overflow.";
    } else if (err == MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS) {
        ret = "Demographic events must be time sorted.";
    } else if (err == MSP_ERR_POPULATION_OVERFLOW) {
        ret = "Population Overflow occured.";
    } else if (err == MSP_ERR_LINKS_OVERFLOW) {
        ret = "Links Overflow occured.";
    } else if (err == MSP_ERR_OUT_OF_BOUNDS) {
        ret = "Array index out of bounds";
    } else if (err == MSP_ERR_BAD_ORDERING) {
        ret = "Bad record ordering requested";
    } else if (err == MSP_ERR_BAD_MUTATION) {
        ret = "Bad mutation provided";
    } else if (err == MSP_ERR_BAD_PARAM_VALUE) {
        ret = "Bad parameter value provided";
    } else if (err == MSP_ERR_UNSUPPORTED_OPERATION) {
        ret = "Operation cannot be performed in current configuration";
    } else if (err == MSP_ERR_BAD_POPULATION_CONFIGURATION) {
        ret = "Bad population configuration provided.";
    } else if (err == MSP_ERR_BAD_POPULATION_ID) {
        ret = "Bad population id provided.";
    } else if (err == MSP_ERR_BAD_MIGRATION_MATRIX) {
        ret = "Bad migration matrix provided.";
    } else if (err == MSP_ERR_BAD_MIGRATION_MATRIX_INDEX) {
        ret = "Bad migration matrix index provided.";
    } else if (err == MSP_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX) {
        ret = "Cannot set diagonal migration matrix elements.";
    } else if (err == MSP_ERR_INFINITE_WAITING_TIME) {
        ret = "Infinite waiting time until next simulation event.";
    } else if (err == MSP_ERR_ASSERTION_FAILED) {
        ret = "Internal error; please file a bug report.";
    } else if (err == MSP_ERR_SOURCE_DEST_EQUAL) {
        ret = "Source and destination populations equal.";
    } else if (err == MSP_ERR_BAD_RECOMBINATION_MAP) {
        ret = "Bad recombination map provided.";
    } else if (err == MSP_ERR_BAD_COALESCENCE_RECORDS) {
        ret = "Bad coalescence records in file.";
    } else if (err == MSP_ERR_BAD_SAMPLES) {
        ret = "Bad sample configuration provided.";
    } else if (err == MSP_ERR_NONBINARY_NEWICK) {
        ret = "Newick export not supported for non binary trees.";
    } else if (err == MSP_ERR_FILE_VERSION_TOO_OLD) {
        ret = "HDF5 file version too old. Please upgrade using the "
            "'msp upgrade' command";
    } else if (err == MSP_ERR_FILE_VERSION_TOO_NEW) {
        ret = "HDF5 file version is too new for this version of msprime. "
            "Please upgrade msprime to the latest version.";
    } else if (err == MSP_ERR_NONBINARY_NEWICK) {
        ret = "Newick export not supported for non binary trees.";
    } else if (err == MSP_ERR_DUPLICATE_TRACKED_LEAF) {
        ret = "Duplicate value provided in tracked leaf list.";
    } else if (err == MSP_ERR_REFCOUNT_NONZERO) {
        ret = "Cannot change the state of the tree sequence when "
            "other objects reference it. Make sure all trees are freed first.";
    } else if (err == MSP_ERR_PTHREAD) {
        ret = "Pthread error. Please file a bug report.";
    } else if (err == MSP_ERR_IO) {
        if (errno != 0) {
            ret = strerror(errno);
        } else {
            ret = "Unspecified IO error";
        }
    } else if (err == MSP_ERR_HDF5) {
        _hdf5_error[0] = '\0';
        if (H5Ewalk2(H5E_DEFAULT, H5E_WALK_UPWARD, hdf5_error_walker, NULL)
                != 0) {
            ret = "Eek! Error handling HDF5 error.";
            goto out;
        }
        ret = _hdf5_error;
    }
out:
    return ret;
}

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

static int
cmp_uint32_t(const void *a, const void *b) {
    const uint32_t *ia = (const uint32_t *) a;
    const uint32_t *ib = (const uint32_t *) b;
    return (*ia > *ib) - (*ia < *ib);
}

static void
segment_init(void **obj, size_t id)
{
    segment_t *seg = (segment_t *) obj;
    seg->id = id + 1;
}

static size_t
msp_get_avl_node_mem_increment(msp_t *self)
{
    return sizeof(void *) + self->avl_node_block_size
            * (sizeof(avl_node_t) + sizeof(void *));
}

static size_t
msp_get_binary_children_mem_increment(msp_t *self)
{
    return sizeof(void *) + self->coalescence_record_block_size
            * (2 * sizeof(uint32_t) + sizeof(void *));
}

static size_t
msp_get_segment_mem_increment(msp_t *self)
{
    /* we have a segment, a pointer to it and an entry in the Fenwick tree */
    size_t s = sizeof(segment_t) + sizeof(void *) + 2 * sizeof(int64_t);
    return sizeof(void *) + self->segment_block_size * s;
}

static size_t
msp_get_node_mapping_mem_increment(msp_t *self)
{
    return sizeof(void *) + self->node_mapping_block_size
            * sizeof(node_mapping_t);
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
msp_get_num_coalescence_record_blocks(msp_t *self)
{
    return self->num_coalescence_record_blocks;
}

size_t
msp_get_used_memory(msp_t *self)
{
    return self->used_memory;
}

size_t
msp_get_num_common_ancestor_events(msp_t *self)
{
    return self->num_ca_events;
}

size_t
msp_get_num_recombination_events(msp_t *self)
{
    return self->num_re_events;
}

int
msp_set_num_loci(msp_t *self, size_t num_loci)
{
    int ret = 0;

    if (num_loci < 1) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->num_loci = (uint32_t) num_loci;
out:
    return ret;
}

int
msp_set_num_populations(msp_t *self, size_t num_populations)
{
    int ret = 0;
    size_t j;

    if (num_populations < 1 || num_populations > UINT8_MAX) {
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
msp_set_scaled_recombination_rate(msp_t *self,
        double scaled_recombination_rate)
{
    int ret = 0;

    if (scaled_recombination_rate < 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->scaled_recombination_rate = scaled_recombination_rate;
out:
    return ret;
}

int
msp_set_population_configuration(msp_t *self, int population_id,
        double initial_size, double growth_rate)
{
    int ret = MSP_ERR_BAD_POPULATION_CONFIGURATION;

    if (population_id < 0 || population_id > (int) self->num_populations) {
        ret = MSP_ERR_BAD_POPULATION_ID;
        goto out;
    }
    self->initial_populations[population_id].initial_size = initial_size;
    self->initial_populations[population_id].growth_rate = growth_rate;
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
msp_set_max_memory(msp_t *self, size_t max_memory)
{
    int ret = 0;

    if (max_memory < 1) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->max_memory = max_memory;
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

int
msp_set_coalescence_record_block_size(msp_t *self, size_t block_size)
{
    int ret = 0;

    if (block_size < 1) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->coalescence_record_block_size = block_size;
out:
    return ret;
}

/* Top level allocators and initialisation */

int
msp_alloc(msp_t *self, size_t sample_size, sample_t *samples, gsl_rng *rng)
{
    int ret = -1;
    size_t j, k, initial_samples;

    memset(self, 0, sizeof(msp_t));
    if (sample_size < 2 || samples == NULL || rng == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->sample_size = (uint32_t) sample_size;
    self->rng = rng;
    self->num_loci = 1;
    self->scaled_recombination_rate = 0.0;
    self->samples = malloc(sample_size * sizeof(sample_t));
    if (self->samples == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    initial_samples = 0;
    for (j = 0; j < sample_size; j++) {
        self->samples[j].population_id = samples[j].population_id;
        self->samples[j].time = samples[j].time;
        if (self->samples[j].time < 0) {
            ret = MSP_ERR_BAD_PARAM_VALUE;
            goto out;
        }
        if (self->samples[j].time == 0) {
            initial_samples++;
        }
    }
    if (initial_samples == 0) {
        ret = MSP_ERR_BAD_SAMPLES;
        goto out;
    }
    /* Set up the historical sampling events */
    self->num_sampling_events = self->sample_size - initial_samples;
    self->sampling_events = NULL;
    if (self->num_sampling_events > 0) {
        self->sampling_events = malloc(self->num_sampling_events *
                sizeof(sampling_event_t));
        if (self->sampling_events == NULL) {
            goto out;
        }
        k = 0;
        for (j = 0; j < self->sample_size; j++) {
            if (self->samples[j].time > 0) {
                self->sampling_events[k].sample = (uint32_t) j;
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
    /* We have one population by default */
    ret = msp_set_num_populations(self, 1);
    if (ret != 0) {
        goto out;
    }
    /* Set sensible defaults for the sample_config and migration matrix */
    self->initial_migration_matrix[0] = 0.0;
    /* Set the memory defaults */
    self->avl_node_block_size = 1024;
    self->node_mapping_block_size = 1024;
    self->segment_block_size = 1024;
    self->max_memory = 1024 * 1024 * 1024; /* 1MiB */
    self->coalescence_record_block_size = 1024;
    /* set up the AVL trees */
    avl_init_tree(&self->breakpoints, cmp_node_mapping, NULL);
    avl_init_tree(&self->overlap_counts, cmp_node_mapping, NULL);
    /* Set up the demographic events */
    self->demographic_events_head = NULL;
    self->demographic_events_tail = NULL;
    self->next_demographic_event = NULL;
    self->state = MSP_STATE_NEW;
out:
    return ret;
}

static int
msp_alloc_memory_blocks(msp_t *self)
{
    int ret = 0;

    self->used_memory = msp_get_avl_node_mem_increment(self)
        + msp_get_segment_mem_increment(self)
        + msp_get_node_mapping_mem_increment(self);
    if (self->used_memory > self->max_memory) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
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
    ret = object_heap_init(&self->binary_children_heap, 2 * sizeof(uint32_t),
           self->coalescence_record_block_size, NULL);
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
    /* Allocate the coalescence records */
    self->coalescence_records = malloc(
            self->coalescence_record_block_size * sizeof(coalescence_record_t));
    self->max_coalescence_records = self->coalescence_record_block_size;
    self->num_coalescence_record_blocks = 1;
    if (self->coalescence_records == NULL) {
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

static uint32_t *
msp_alloc_children(msp_t *self, uint32_t num_children)
{
    uint32_t *ret = NULL;

    /* Binary children are by far the most common, so we optimise for this
     * case with a fixed memory heap. We malloc the rest directly. This
     * means that we must always iterate through the set of coalescence
     * records when freeing and be careful to always leave the records in
     * a consistent state. */
    if (num_children == 2) {
        if (object_heap_empty(&self->binary_children_heap)) {
            self->used_memory += msp_get_binary_children_mem_increment(self);
            if (self->used_memory > self->max_memory) {
                goto out;
            }
            if (object_heap_expand(&self->binary_children_heap) != 0) {
                goto out;
            }
        }
        ret = (uint32_t *) object_heap_alloc_object(
                &self->binary_children_heap);
    } else {
        ret = malloc(num_children * sizeof(uint32_t));
    }
out:
    return ret;
}

static void
msp_free_children(msp_t *self, uint32_t num_children, uint32_t *children)
{
    if (num_children == 2) {
        object_heap_free_object(&self->binary_children_heap, children);
    } else {
        free(children);
    }
}

int
msp_free(msp_t *self)
{
    int ret = -1;
    demographic_event_t *de = self->demographic_events_head;
    demographic_event_t *tmp;
    coalescence_record_t *cr;
    size_t j;

    /* Because we malloc the children for non-binary records
     * directly we must free these memory blocks */
    for (j = 0; j < self->num_coalescence_records; j++) {
        cr = &self->coalescence_records[j];
        if (cr->children != NULL) {
            msp_free_children(self, cr->num_children, cr->children);
        }
    }
    while (de != NULL) {
        tmp = de->next;
        free(de);
        de = tmp;
    }
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
    if (self->samples != NULL) {
        free(self->samples);
    }
    if (self->sampling_events != NULL) {
        free(self->sampling_events);
    }
    /* free the object heaps */
    object_heap_free(&self->avl_node_heap);
    object_heap_free(&self->segment_heap);
    object_heap_free(&self->node_mapping_heap);
    object_heap_free(&self->binary_children_heap);
    fenwick_free(&self->links);
    if (self->coalescence_records != NULL) {
        free(self->coalescence_records);
    }
    ret = 0;
    return ret;
}


static inline avl_node_t * WARN_UNUSED
msp_alloc_avl_node(msp_t *self)
{
    avl_node_t *ret = NULL;

    if (object_heap_empty(&self->avl_node_heap)) {
        self->used_memory += msp_get_avl_node_mem_increment(self);
        if (self->used_memory > self->max_memory) {
            goto out;
        }
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
        self->used_memory += msp_get_node_mapping_mem_increment(self);
        if (self->used_memory > self->max_memory) {
            goto out;
        }
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

static segment_t * WARN_UNUSED
msp_alloc_segment(msp_t *self, uint32_t left, uint32_t right, uint32_t value,
        uint32_t population_id, segment_t *prev, segment_t *next)
{
    segment_t *seg = NULL;

    if (object_heap_empty(&self->segment_heap)) {
        self->used_memory += msp_get_segment_mem_increment(self);
        if (self->used_memory > self->max_memory) {
            goto out;
        }
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
    assert(population_id < UINT8_MAX);
    seg->population_id = (uint8_t) population_id;
out:
    return seg;
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
msp_print_segment_chain(msp_t *self, segment_t *head, FILE *out)
{
    segment_t *s = head;

    fprintf(out, "[%d]", s->population_id);
    while (s != NULL) {
        fprintf(out, "[(%d-%d) %d] ", s->left, s->right, s->value);
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
                assert(u->population_id == j);
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
                ss = s; /* just to keep compiler happy - see below also */
                right = u->right;
                u = u->next;
            }
            alt_total_links += right - left - 1;
            node = node->next;
        }
    }
    assert(total_links == fenwick_get_total(&self->links));
    assert(total_links == alt_total_links);
    assert(total_segments == object_heap_get_num_allocated(
                &self->segment_heap));
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
    for (j = 0; j < self->sample_size; j++) {
        if (self->samples[j].time > self->time) {
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
    avl_node_t *node;
    node_mapping_t *nm;
    segment_t *u;
    coalescence_record_t *cr;
    demographic_event_t *de;
    sampling_event_t *se;
    int64_t v;
    uint32_t j, k;
    double gig = 1024.0 * 1024;
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
    fprintf(out, "used_memory = %f MiB\n", (double) self->used_memory / gig);
    fprintf(out, "max_memory  = %f MiB\n", (double) self->max_memory / gig);
    fprintf(out, "n = %d\n", self->sample_size);
    fprintf(out, "m = %d\n", self->num_loci);
    fprintf(out, "Samples = \n");
    for (j = 0; j < self->sample_size; j++) {
        fprintf(out, "\t%d\tpopulation=%d\ttime=%f\n", j, self->samples[j].population_id,
                self->samples[j].time);
    }
    fprintf(out, "Sampling events:\n");
    for (j = 0; j < self->num_sampling_events; j++) {
        if (j == self->next_sampling_event) {
            fprintf(out, "  ***");
        }
        se = &self->sampling_events[j];
        fprintf(out, "\t");
        fprintf(out, "%d @ %f in deme %d\n", se->sample, se->time, se->population_id);
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
    for (node = self->breakpoints.head; node != NULL; node = node->next) {
        nm = (node_mapping_t *) node->item;
        fprintf(out, "\t%d -> %d\n", nm->left, nm->value);
    }
    fprintf(out, "Overlap count = %d\n", avl_count(&self->overlap_counts));
    for (node = self->overlap_counts.head; node != NULL; node = node->next) {
        nm = (node_mapping_t *) node->item;
        fprintf(out, "\t%d -> %d\n", nm->left, nm->value);
    }
    fprintf(out, "Coalescence records = %ld\n",
            (long) self->num_coalescence_records);
    for (j = 0; j < self->num_coalescence_records; j++) {
        cr = &self->coalescence_records[j];
        fprintf(out, "\t%f\t%f\t%d\t(", cr->left, cr->right,
                cr->node);
        for (k = 0; k < cr->num_children; k++) {
            fprintf(out, "%d", cr->children[k]);
            if (k < cr->num_children - 1) {
                fprintf(out, ", ");
            }
        }
        fprintf(out, ")\t%f\t%d\n", cr->time, cr->population_id);
    }
    fprintf(out, "Memory heaps\n");
    fprintf(out, "avl_node_heap:");
    object_heap_print_state(&self->avl_node_heap, out);
    fprintf(out, "segment_heap:");
    object_heap_print_state(&self->segment_heap, out);
    fprintf(out, "node_mapping_heap:");
    object_heap_print_state(&self->node_mapping_heap, out);
    fprintf(out, "binary_children_heap:");
    object_heap_print_state(&self->binary_children_heap, out);
    msp_verify(self);
out:
    if (ancestors != NULL) {
        free(ancestors);
    }
    return ret;
}

static int WARN_UNUSED
msp_move_individual(msp_t *self, avl_node_t *node, avl_tree_t *source,
        uint32_t dest_pop)
{
    int ret = 0;
    segment_t *ind, *x;

    ind = (segment_t *) node->item;
    avl_unlink_node(source, node);
    msp_free_avl_node(self, node);
    /* Need to set the population_id for each segment. */
    x = ind;
    while (x != NULL) {
        x->population_id = (uint8_t) dest_pop;
        x = x->next;
    }
    ret = msp_insert_individual(self, ind);
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
msp_record_coalescence(msp_t *self, uint32_t left, uint32_t right,
        uint32_t num_children, uint32_t *children, uint32_t node,
        uint32_t population_id)
{
    int ret = 0;
    int equal;
    uint32_t j;
    coalescence_record_t *cr;
    coalescence_record_t *lcr;

    if (self->num_coalescence_records == self->max_coalescence_records - 1) {
        /* Grow the array */
        self->max_coalescence_records += self->coalescence_record_block_size;;
        cr = realloc(self->coalescence_records,
                self->max_coalescence_records * sizeof(coalescence_record_t));
        if (cr == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        self->coalescence_records = cr;
        self->num_coalescence_record_blocks++;
    }
    /* Sort the children */
    qsort(children, num_children, sizeof(uint32_t), cmp_uint32_t);

    cr = &self->coalescence_records[self->num_coalescence_records];
    if (self->num_coalescence_records != 0) {
        lcr = &self->coalescence_records[self->num_coalescence_records - 1];
        if (lcr->right == left && lcr->num_children == num_children
                && lcr->node == node) {
            /* Compare the children */
            equal = 1;
            for (j = 0; equal && j < num_children; j++) {
                equal = equal && (children[j] == lcr->children[j]);
            }
            if (equal) {
                /* squash this record into the last */
                lcr->right = right;
                cr = NULL;
                msp_free_children(self, num_children, children);
            }
        }
    }
    if (cr != NULL) {
        cr->num_children = num_children;
        cr->children = children;
        cr->left = (double) left;
        cr->right = (double) right;
        cr->node = node;
        cr->time = self->time;
        cr->population_id = (uint8_t) population_id;
        self->num_coalescence_records++;
    }
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

static int WARN_UNUSED
msp_common_ancestor_event(msp_t *self, uint32_t population_id)
{
    int ret = 0;
    int coalescence = 0;
    int defrag_required = 0;
    uint32_t j, n, l, r, l_min, r_max, v, *children;
    avl_tree_t *ancestors;
    avl_node_t *node;
    node_mapping_t *nm, search;
    segment_t *x, *y, *z, *alpha, *beta;

    ancestors = &self->populations[population_id].ancestors;
    self->num_ca_events++;
    /* Choose x and y */
    n = avl_count(ancestors);
    j = (uint32_t) gsl_rng_uniform_int(self->rng, n);
    node = avl_at(ancestors, j);
    assert(node != NULL);
    x = (segment_t *) node->item;
    avl_unlink_node(ancestors, node);
    msp_free_avl_node(self, node);
    j = (uint32_t) gsl_rng_uniform_int(self->rng, n - 1);
    node = avl_at(ancestors, j);
    assert(node != NULL);
    y = (segment_t *) node->item;
    avl_unlink_node(ancestors, node);
    msp_free_avl_node(self, node);
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
                    self->next_node++;
                    /* Check for overflow */
                    assert(self->next_node != 0);
                }
                v = self->next_node - 1;
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
                    alpha = msp_alloc_segment(self, l, r, v, population_id,
                            NULL, NULL);
                    if (alpha == NULL) {
                        ret = MSP_ERR_NO_MEMORY;
                        goto out;
                    }
                }
                children = msp_alloc_children(self, 2);
                if (children == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
                children[0] = x->value;
                children[1] = y->value;
                ret = msp_record_coalescence(self, l, r, 2, children, v,
                        population_id);
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
 * any number of ancestors to merge. The avl tree is a priority queue in
 * sorted by left coordinate.
 */
static int WARN_UNUSED
msp_merge_ancestors(msp_t *self, avl_tree_t *Q, uint32_t population_id)
{
    int ret = MSP_ERR_GENERIC;
    int coalescence = 0;
    int defrag_required = 0;
    uint32_t j, l, r, h, r_max, next_l, l_min, v, *children;
    avl_node_t *node;
    node_mapping_t *nm, search;
    segment_t *x, *z, *alpha;
    segment_t **H = NULL;

    H = malloc(avl_count(Q) * sizeof(segment_t *));
    if (H == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
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
                self->next_node++;
                /* Check for overflow */
                assert(self->next_node != 0);
            }
            v = self->next_node - 1;
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
            /* Create the record and update the priority queue */
            children = msp_alloc_children(self, h);
            if (children == NULL) {
                ret = MSP_ERR_NO_MEMORY;
                goto out;
            }
            for (j = 0; j < h; j++) {
                x = H[j];
                children[j] = x->value;
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
            ret = msp_record_coalescence(self, l, r, h, children, v,
                    population_id);
            if (ret != 0) {
                goto out;
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
msp_migration_event(msp_t *self, uint32_t source_pop, uint32_t dest_pop)
{
    int ret = 0;
    uint32_t j;
    avl_node_t *node;
    avl_tree_t *source = &self->populations[source_pop].ancestors;

    self->num_migration_events[
        source_pop * self->num_populations + dest_pop]++;
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
    coalescence_record_t *cr;
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
    for (j = 0; j < self->num_coalescence_records; j++) {
        cr = &self->coalescence_records[j];
        if (cr->children != NULL) {
            msp_free_children(self, cr->num_children, cr->children);
        }
        cr->children = NULL;
    }
    return ret;
}

static int WARN_UNUSED
msp_insert_sample(msp_t *self, uint32_t sample, uint8_t population)
{
    int ret = MSP_ERR_GENERIC;
    segment_t *u;

    u = msp_alloc_segment(self, 0, self->num_loci, sample, population,
            NULL, NULL);
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

int
msp_reset(msp_t *self)
{
    int ret = 0;
    uint32_t j, population_id;
    population_t *pop, *initial_pop;
    size_t N = self->num_populations;

    ret = msp_reset_memory_state(self);
    if (ret != 0) {
        goto out;
    }
    /* Set up the initial segments and algorithm state */
    for (population_id = 0; population_id < N; population_id++) {
        pop = &self->populations[population_id];
        assert(avl_count(&pop->ancestors) == 0);
        /* Set the initial population parameters */
        initial_pop = &self->initial_populations[population_id];
        pop->growth_rate = initial_pop->growth_rate;
        pop->initial_size = initial_pop->initial_size;
        pop->start_time = 0.0;
    }
    /* Set up the sample */
    for (j = 0; j < self->sample_size; j++) {
        if (self->samples[j].time == 0.0) {
            ret = msp_insert_sample(self, j, self->samples[j].population_id);
            if (ret != 0) {
                goto out;
            }
        }
    }
    self->next_node = self->sample_size;
    self->next_demographic_event = self->demographic_events_head;
    memcpy(self->migration_matrix, self->initial_migration_matrix,
            N * N * sizeof(double));
    ret = msp_insert_overlap_count(self, 0, self->sample_size);
    if (ret != 0) {
        goto out;
    }
    ret = msp_insert_overlap_count(self, self->num_loci,
            self->sample_size + 1);
    if (ret != 0) {
        goto out;
    }
    self->time = 0.0;
    self->next_sampling_event = 0;
    self->num_coalescence_records = 0;
    self->num_re_events = 0;
    self->num_ca_events = 0;
    self->num_trapped_re_events = 0;
    self->num_multiple_re_events = 0;
    memset(self->num_migration_events, 0, N * N * sizeof(size_t));
    self->state = MSP_STATE_INITIALISED;
out:
    return ret;
}

/*
 * Sets up the memory heaps, initial population and trees.
 */
int WARN_UNUSED
msp_initialise(msp_t *self)
{
    int ret = -1;
    uint32_t j;

    /* These should really be proper checks with a return value */
    assert(self->sample_size > 1);
    assert(self->num_loci >= 1);
    assert(self->num_populations >= 1);

    ret = msp_alloc_memory_blocks(self);
    if (ret != 0) {
        goto out;
    }
    /* First check that the sample configuration makes sense */
    for (j = 0; j < self->sample_size; j++) {
        if (self->samples[j].population_id >= self->num_populations) {
            ret = MSP_ERR_BAD_SAMPLES;
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
msp_get_common_ancestor_waiting_time(msp_t *self, uint32_t population_id)
{
    double ret = DBL_MAX;
    population_t *pop = &self->populations[population_id];
    /* Need to perform n * (n - 1) as a double due to overflow */
    double n = (double) avl_count(&pop->ancestors);
    double lambda = n * (n - 1.0);
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

int WARN_UNUSED
msp_run(msp_t *self, double max_time, unsigned long max_events)
{
    int ret = 0;
    double lambda, t_temp, t_wait, ca_t_wait, re_t_wait, mig_t_wait,
           sampling_event_time, demographic_event_time;
    int64_t num_links;
    uint32_t j, k, n;
    uint32_t ca_pop_id, mig_source_pop, mig_dest_pop;
    unsigned long events = 0;
    sampling_event_t *se;

    if (self->state == MSP_STATE_INITIALISED) {
        self->state = MSP_STATE_SIMULATING;
    }
    if (self->state != MSP_STATE_SIMULATING) {
        ret = MSP_ERR_BAD_STATE;
        goto out;
    }
    while (msp_get_num_ancestors(self) > 0
            && self->time < max_time && events < max_events) {
        events++;
        num_links = fenwick_get_total(&self->links);
        ret = msp_sanity_check(self, num_links);
        if (ret != 0) {
            goto out;
        }
        /* Recombination */
        lambda = (double) num_links * self->scaled_recombination_rate;
        re_t_wait = DBL_MAX;
        if (lambda != 0.0) {
            re_t_wait = gsl_ran_exponential(self->rng, 1.0 / lambda);
        }
        /* Common ancestors */
        ca_t_wait = DBL_MAX;
        ca_pop_id = 0;
        for (j = 0; j < self->num_populations; j++) {
            t_temp = msp_get_common_ancestor_waiting_time(self, j);
            if (t_temp < ca_t_wait) {
                ca_t_wait = t_temp;
                ca_pop_id = j;
            }
        }
        /* Migration */
        mig_t_wait = DBL_MAX;
        mig_source_pop = 0;
        mig_dest_pop = 0;
        for (j = 0; j < self->num_populations; j++) {
            n = avl_count(&self->populations[j].ancestors);
            for (k = 0; k < self->num_populations; k++) {
                lambda = n * self->migration_matrix[
                    j * self->num_populations + k];
                if (lambda != 0.0) {
                    t_temp = gsl_ran_exponential(self->rng, 1.0 / lambda);
                    if (t_temp < mig_t_wait) {
                        mig_t_wait = t_temp;
                        /* m[j, k] is the rate at which migrants move from
                         * population k to j forwards in time. Backwards
                         * in time, we move the individual from from
                         * population j into population k.
                         */
                        mig_source_pop = j;
                        mig_dest_pop = k;
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
                ret = msp_recombination_event(self);
            } else if (ca_t_wait == t_wait) {
                ret = msp_common_ancestor_event(self, ca_pop_id);
            } else {
                ret = msp_migration_event(self, mig_source_pop, mig_dest_pop);
            }
            if (ret != 0) {
                goto out;
            }
        }
    }
    if (msp_get_num_ancestors(self) != 0) {
        ret = 1;
        if (self->time >= max_time) {
            ret = 2;
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
    *end_time = t;
out:
    return ret;
}

size_t
msp_get_sample_size(msp_t *self)
{
    return (size_t) self->sample_size;
}

size_t
msp_get_num_loci(msp_t *self)
{
    return (size_t) self->num_loci;
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
msp_get_num_coalescence_records(msp_t *self)
{
    return self->num_coalescence_records;
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

    memcpy(migration_matrix, self->migration_matrix, N * N * sizeof(double));
    return 0;
}

int WARN_UNUSED
msp_get_num_migration_events(msp_t *self, size_t *num_migration_events)
{
    size_t N = self->num_populations;

    memcpy(num_migration_events, self->num_migration_events,
        N * N * sizeof(size_t));
    return 0;
}

int WARN_UNUSED
msp_get_coalescence_records(msp_t *self, coalescence_record_t **coalescence_records)
{
    *coalescence_records = self->coalescence_records;
    return 0;
}

int WARN_UNUSED
msp_get_samples(msp_t *self, sample_t **samples)
{
    *samples = self->samples;
    return 0;
}

int WARN_UNUSED
msp_get_population_configuration(msp_t *self, size_t population_id,
        double *initial_size, double *growth_rate)
{
    int ret = 0;
    population_t *pop;

    if (population_id > self->num_populations) {
        ret = MSP_ERR_BAD_POPULATION_ID;
        goto out;
    }
    pop = &self->populations[population_id];
    *initial_size = pop->initial_size;
    *growth_rate = pop->growth_rate;
out:
    return ret;
}

int WARN_UNUSED
msp_get_population(msp_t *self, size_t population_id,
        population_t **population)
{
    int ret = 0;

    if (population_id > self->num_populations) {
        ret = MSP_ERR_BAD_POPULATION_ID;
        goto out;
    }
    *population = &self->populations[population_id];
out:
    return ret;
}

/* Demographic events */
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
    ret_event->time = time;
    /* now insert this event at the end of the chain. */
    if (self->demographic_events_head == NULL) {
        self->demographic_events_head = ret_event;
        self->demographic_events_tail = ret_event;
    } else {
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

    if (population_id >= self->num_populations) {
        ret = MSP_ERR_BAD_POPULATION_ID;
        goto out;
    }
    pop = &self->populations[population_id];
    /* If initial_size is not specified, calculate the initial_size of the
     * population over the coming time period based on the growth rate over
     * the preceeding period.
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
    int pid = event->params.population_parameters_change.population_id;
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
msp_print_population_parameters_change(msp_t *self,
        demographic_event_t *event, FILE *out)
{
    fprintf(out,
            "%f\tpopulation_parameters_change: %d -> initial_size=%f, growth_rate=%f\n",
            event->time,
            event->params.population_parameters_change.population_id,
            event->params.population_parameters_change.initial_size,
            event->params.population_parameters_change.growth_rate);
}

int
msp_add_population_parameters_change(msp_t *self, double time, int population_id,
        double initial_size, double growth_rate)
{
    int ret = -1;
    demographic_event_t *de;
    int N = (int) self->num_populations;

    if (population_id < -1 || population_id >= N) {
        ret = MSP_ERR_BAD_POPULATION_ID;
        goto out;
    }
    if (initial_size < 0) {
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
    double rate  = event->params.migration_rate_change.migration_rate;
    int N = (int) self->num_populations;

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
msp_print_migration_rate_change(msp_t *self,
        demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tmigration_rate_change: %d -> %f\n",
            event->time,
            event->params.migration_rate_change.matrix_index,
            event->params.migration_rate_change.migration_rate);
}

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
    int source = event->params.mass_migration.source;
    int dest = event->params.mass_migration.destination;
    double p = event->params.mass_migration.proportion;
    int N = (int) self->num_populations;
    avl_node_t *node, *next;
    avl_tree_t *pop;

    /* This should have been caught on adding the event */
    if (source < 0 || source > N || dest < 0 || dest > N) {
        ret = MSP_ERR_ASSERTION_FAILED;
        goto out;
    }
    /*
     * Move lineages from source to dest with propabality p.
     */
    pop = &self->populations[source].ancestors;
    node = pop->head;
    while (node != NULL) {
        next = node->next;
        if (gsl_rng_uniform(self->rng) < p) {
            ret = msp_move_individual(self, node, pop, (uint32_t) dest);
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
msp_print_mass_migration(msp_t *self, demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tmass_migration: %d -> %d p = %f\n",
            event->time,
            event->params.mass_migration.source,
            event->params.mass_migration.destination,
            event->params.mass_migration.proportion);
}

int WARN_UNUSED
msp_add_mass_migration(msp_t *self, double time, int source, int destination,
        double proportion)
{
    int ret = 0;
    demographic_event_t *de;
    int N = (int) self->num_populations;

    if (source < 0 || source >= N || destination < 0 || destination >= N) {
        ret = MSP_ERR_BAD_POPULATION_ID;
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

/* Bottlenecks */

static int
msp_bottleneck(msp_t *self, demographic_event_t *event)
{
    int ret = 0;
    int population_id = event->params.bottleneck.population_id;
    double p = event->params.bottleneck.proportion;
    int N = (int) self->num_populations;
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
     * during this bottleneck.
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
    ret = msp_merge_ancestors(self, &Q, (uint32_t) population_id);
out:
    return ret;
}

static void
msp_print_bottleneck(msp_t *self, demographic_event_t *event, FILE *out)
{
    fprintf(out, "%f\tbottleneck: %d I = %f\n",
            event->time,
            event->params.bottleneck.population_id,
            event->params.bottleneck.proportion);
}

int WARN_UNUSED
msp_add_bottleneck(msp_t *self, double time, int population_id,
        double proportion)
{
    int ret = 0;
    demographic_event_t *de;
    int N = (int) self->num_populations;

    if (population_id < 0 || population_id >= N) {
        ret = MSP_ERR_BAD_POPULATION_ID;
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
    de->params.bottleneck.population_id = population_id;
    de->params.bottleneck.proportion = proportion;
    de->change_state = msp_bottleneck;
    de->print_state = msp_print_bottleneck;
    ret = 0;
out:
    return ret;
}


