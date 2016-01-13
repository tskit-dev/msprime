/*
** Copyright (C) 2015 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

static char _hdf5_error[MSP_HDF5_ERR_MSG_SIZE];

static herr_t
hdf5_error_walker(unsigned n, const H5E_error2_t *err_desc, void *client_data)
{
    /* We only copy the message from the first element in the error stack */
    if (_hdf5_error[0] == '\0') {
        snprintf(_hdf5_error, MSP_HDF5_ERR_MSG_SIZE,
                "HDF5 Error: %d: %d:'%s'",
                err_desc->maj_num, err_desc->min_num, err_desc->desc);
    }
    return 0;
}

const char *
msp_strerror(int err)
{
    const char *ret = "Unknown error";

    if (err == MSP_ERR_NO_MEMORY) {
        ret = "Out of memory. Try increasing the max_memory parameter";
    } else if (err == MSP_ERR_GENERIC) {
        ret = "Generic error; please file a bug report";
    } else if (err == MSP_ERR_FILE_FORMAT) {
        ret = "File format error";
    } else if (err == MSP_ERR_FILE_VERSION) {
        ret = "Unsupported file format version";
    } else if (err == MSP_ERR_BAD_MODE) {
        ret = "Bad tree file mode";
    } else if (err == MSP_ERR_BAD_POP_MODEL) {
        ret = "Bad population model type";
    } else if (err == MSP_ERR_NEWICK_OVERFLOW) {
        ret = "Newick string generation overflow.";
    } else if (err == MSP_ERR_UNSORTED_POP_MODELS) {
        ret = "Population models must sorted by start_time";
    } else if (err == MSP_ERR_POPULATION_OVERFLOW) {
        ret = "Population Overflow occured.";
    } else if (err == MSP_ERR_LINKS_OVERFLOW) {
        ret = "Links Overflow occured.";
    } else if (err == MSP_ERR_OUT_OF_BOUNDS) {
        ret = "Array index out of bounds";
    } else if (err == MSP_ERR_UNSUPPORTED_FILE_VERSION) {
        ret = "Unsupported file format version.";
    } else if (err == MSP_ERR_BAD_ORDERING) {
        ret = "Bad record ordering requested";
    } else if (err == MSP_ERR_BAD_MUTATION) {
        ret = "Bad mutation provided";
    } else if (err == MSP_ERR_BAD_PARAM_VALUE) {
        ret = "Bad parameter value provided";
    } else if (err == MSP_ERR_UNSUPPORTED_OPERATION) {
        ret = "Operation cannot be performed in current configuration";
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
    return (ia->index > ib->index) - (ia->index < ib->index);
}

static int
cmp_node_mapping(const void *a, const void *b) {
    const node_mapping_t *ia = (const node_mapping_t *) a;
    const node_mapping_t *ib = (const node_mapping_t *) b;
    return (ia->left > ib->left) - (ia->left < ib->left);
}


static void
segment_init(void **obj, size_t index)
{
    segment_t *seg = (segment_t *) obj;
    seg->index = index + 1;
}

/* population models */

static double
constant_population_model_get_size(population_model_t *self, double t)
{
    return self->param;
}

static double
constant_population_model_get_waiting_time(population_model_t *self,
        double lambda_coancestry, double t, gsl_rng* rng)
{
    return self->param * gsl_ran_exponential(rng, 1.0 / lambda_coancestry);
}

static double
exponential_population_model_get_size(population_model_t *self, double t)
{
    double alpha = self->param;
    return self->initial_size * exp(-alpha * (t - self->start_time));
}

static double
exponential_population_model_get_waiting_time(population_model_t *self,
        double lambda_coancestry, double t, gsl_rng* rng)
{
    double ret = DBL_MAX;
    double alpha = self->param;
    double u = gsl_ran_exponential(rng, 1.0 / lambda_coancestry);
    double dt, z;

    if (alpha == 0.0) {
        ret = self->initial_size * u;
    } else {
        dt = t - self->start_time;
        z = 1 + alpha * self->initial_size * exp(-alpha * dt) * u;
        /* if z is <= 0 no coancestry can occur */
        if (z > 0) {
            ret = log(z) / alpha;
        }
    }
    return ret;
}

/* Extends the list of population models by one and increments the
 * num_population_models counter. Also make sure the models are
 * inserted in time sorted order.
 */
static int WARN_UNUSED
msp_add_population_model(msp_t *self, double start_time)
{
    int ret = -1;
    population_model_t *m;

    self->num_population_models++;
    m = realloc(self->population_models,
            self->num_population_models * sizeof(population_model_t));
    if (m == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->population_models = m;
    if (self->num_population_models >= 2) {
        m = &self->population_models[self->num_population_models - 2];
        if (start_time <= m->start_time) {
            ret = MSP_ERR_BAD_PARAM_VALUE;
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

int
msp_add_constant_population_model(msp_t *self, double start_time, double size)
{
    int ret = -1;
    population_model_t *model;

    ret = msp_add_population_model(self, start_time);
    if (ret < 0) {
        goto out;
    }
    if (size == 0.0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    model = &self->population_models[self->num_population_models - 1];
    model->start_time = start_time;
    model->param = size;
    model->type = POP_MODEL_CONSTANT;
    model->get_size = constant_population_model_get_size;
    model->get_waiting_time = constant_population_model_get_waiting_time;
    ret = 0;
out:
    return ret;
}

int
msp_add_exponential_population_model(msp_t *self, double start_time, double alpha)
{
    int ret = -1;
    population_model_t *model;

    ret = msp_add_population_model(self, start_time);
    if (ret < 0) {
        goto out;
    }
    model = &self->population_models[self->num_population_models - 1];
    model->start_time = start_time;
    model->param = alpha;
    model->type = POP_MODEL_EXPONENTIAL;
    model->get_size = exponential_population_model_get_size;
    model->get_waiting_time = exponential_population_model_get_waiting_time;
    ret = 0;
out:
    return ret;
}

static size_t
msp_get_avl_node_mem_increment(msp_t *self)
{
    return sizeof(void *) + self->avl_node_block_size
            * (sizeof(avl_node_t) + sizeof(void *));
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

int
msp_set_random_seed(msp_t *self, unsigned long random_seed)
{
    self->random_seed = random_seed;
    gsl_rng_set(self->rng, self->random_seed);
    return 0;
}

int
msp_set_num_loci(msp_t *self, uint32_t num_loci)
{
    int ret = 0;

    if (num_loci < 1) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->num_loci = num_loci;
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

int msp_set_node_mapping_block_size(msp_t *self, size_t block_size)
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

int msp_set_segment_block_size(msp_t *self, size_t block_size)
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

int msp_set_avl_node_block_size(msp_t *self, size_t block_size)
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

int msp_set_coalescence_record_block_size(msp_t *self, size_t block_size)
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

int
msp_alloc(msp_t *self, uint32_t sample_size)
{
    int ret = -1;

    memset(self, 0, sizeof(msp_t));
    if (sample_size < 2) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->sample_size = sample_size;
    self->rng = gsl_rng_alloc(gsl_rng_default);
    if (self->rng == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    /* Set the parameter defaults */
    msp_set_random_seed(self, 1);
    self->num_loci = 1;
    self->scaled_recombination_rate = 0.0;
    /* Set the memory defaults */
    self->avl_node_block_size = 1024;
    self->node_mapping_block_size = 1024;
    self->segment_block_size = 1024;
    self->max_memory = 1024 * 1024 * 1024; /* 1MiB */
    self->coalescence_record_block_size = 1024;
    /* set up the AVL trees */
    avl_init_tree(&self->ancestral_population, cmp_individual, NULL);
    avl_init_tree(&self->breakpoints, cmp_node_mapping, NULL);
    avl_init_tree(&self->overlap_counts, cmp_node_mapping, NULL);
    /* Add the base population model */
    ret = msp_add_constant_population_model(self, -1.0, 1.0);
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

static int
msp_is_initialised(msp_t *self)
{
    return self->used_memory != 0;
}

/*
 * Returns true if the simulation has completed.
 */
int
msp_is_completed(msp_t *self)
{
    size_t n = avl_count(&self->ancestral_population);
    return msp_is_initialised(self) && n == 0;
}

int
msp_free(msp_t *self)
{
    int ret = -1;

    if (self->population_models != NULL) {
        free(self->population_models);
    }
    if (self->rng != NULL) {
        gsl_rng_free(self->rng);
    }
    /* free the object heaps */
    object_heap_free(&self->avl_node_heap);
    object_heap_free(&self->segment_heap);
    object_heap_free(&self->node_mapping_heap);
    fenwick_free(&self->links);
    if (self->coalescence_records != NULL) {
        free(self->coalescence_records);
    }
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
    if (ret == NULL) {
        goto out;
    }
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
        segment_t *prev, segment_t *next)
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
out:
    return seg;
}

/*
 * Returns the segment with the specified index.
 */
static segment_t *
msp_get_segment(msp_t *self, size_t index)
{
    segment_t *u = object_heap_get_object(&self->segment_heap, index - 1);

    assert(u != NULL);
    assert(u->index == index);
    return u;
}

static void
msp_free_segment(msp_t *self, segment_t *seg)
{
    object_heap_free_object(&self->segment_heap, seg);
    fenwick_set_value(&self->links, seg->index, 0);
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
    node = avl_insert_node(&self->ancestral_population, node);
    assert(node != NULL);
out:
    return ret;
}

static void
msp_print_segment_chain(msp_t *self, segment_t *head)
{
    unsigned int j;
    segment_t *s = head;
    j = 1;
    while (s != NULL) {
        while (j < s->left) {
            printf("%4c", '-');
            j++;
        }
        while (j <= s->right) {
            printf("%4d", s->value);
            j++;
        }
        s = s->next;
    }
    while (j <= self->num_loci) {
        printf("%4c", '-');
        j++;
    }
    printf("\n");
}


void
msp_verify(msp_t *self)
{
    int64_t s, ss, total_links, left, right, alt_total_links;
    size_t total_segments = 0;
    size_t total_avl_nodes = 0;
    avl_node_t *node;
    segment_t *u;

    total_links = 0;
    alt_total_links = 0;
    node = (&self->ancestral_population)->head;
    while (node != NULL) {
        u = (segment_t *) node->item;
        assert(u->prev == NULL);
        left = u->left;
        while (u != NULL) {
            total_segments++;
            assert(u->left < u->right);
            assert(u->right <= self->num_loci);
            /*
             TODO add back in checking for tree integrity
            for (j = u->left; j <= u->right; j++) {
                msp_get_tree(self, j, pi, tau);
                assert(pi[u->value] == 0);
            }
            */
            if (u->prev != NULL) {
                s = u->right - u->prev->right;
            } else {
                s = u->right - u->left - 1;
            }
            ss = fenwick_get_value(&self->links, u->index);
            total_links += ss;
            assert(s == ss);
            ss = s; /* just to keep compiler happy - see below also */
            right = u->right;
            u = u->next;
        }
        alt_total_links += right - left - 1;
        node = node->next;
    }
    assert(total_links == fenwick_get_total(&self->links));
    assert(total_links == alt_total_links);
    assert(total_segments == object_heap_get_num_allocated(
                &self->segment_heap));
    total_avl_nodes = avl_count(&self->ancestral_population)
            + avl_count(&self->breakpoints)
            + avl_count(&self->overlap_counts);
    assert(total_avl_nodes == object_heap_get_num_allocated(
                &self->avl_node_heap));
    assert(total_avl_nodes - avl_count(&self->ancestral_population)
            == object_heap_get_num_allocated(&self->node_mapping_heap));
    if (total_avl_nodes == total_segments) {
        /* do nothing - this is just to keep the compiler happy when
         * asserts are turned off.
         */
    }
}

int
msp_print_state(msp_t *self)
{
    int ret = 0;
    avl_node_t *node;
    node_mapping_t *nm;
    segment_t *u;
    coalescence_record_t *cr;
    population_model_t *m;
    int64_t v;
    uint32_t j;
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
    printf("used_memory = %f MiB\n", (double) self->used_memory / gig);
    printf("max_memory  = %f MiB\n", (double) self->max_memory / gig);
    printf("n = %d\n", self->sample_size);
    printf("m = %d\n", self->num_loci);
    printf("random seed = %ld\n", self->random_seed);
    printf("num_links = %ld\n", (long) fenwick_get_total(&self->links));
    printf("population = %d\n", avl_count(&self->ancestral_population));
    printf("population models = %d\n", (int) self->num_population_models - 1);
    for (j = 1; j < self->num_population_models; j++) {
        m = &self->population_models[j];
        printf("\t start_time=%f, type=%d, param=%f\n", m->start_time,
                m->type, m->param);
    }
    printf("time = %f\n", self->time);
    for (j = 0; j < msp_get_num_ancestors(self); j++) {
        printf("\t");
        msp_print_segment_chain(self, ancestors[j]);
    }
    printf("Fenwick tree\n");
    for (j = 1; j <= (uint32_t) fenwick_get_size(&self->links); j++) {
        u = msp_get_segment(self, j);
        v = fenwick_get_value(&self->links, j);
        if (v != 0) {
            printf("\t%ld\ti=%d l=%d r=%d v=%d prev=%p next=%p\n", (long) v,
                    (int) u->index, u->left,
                    u->right, (int) u->value, (void *) u->prev, (void *) u->next);
        }
    }
    printf("Breakpoints = %d\n", avl_count(&self->breakpoints));
    for (node = self->breakpoints.head; node != NULL; node = node->next) {
        nm = (node_mapping_t *) node->item;
        printf("\t%d -> %d\n", nm->left, nm->value);
    }
    printf("Overlap count = %d\n", avl_count(&self->overlap_counts));
    for (node = self->overlap_counts.head; node != NULL; node = node->next) {
        nm = (node_mapping_t *) node->item;
        printf("\t%d -> %d\n", nm->left, nm->value);
    }
    printf("Coalescence records = %ld\n", (long) self->num_coalescence_records);
    for (j = 0; j < self->num_coalescence_records; j++) {
        cr = &self->coalescence_records[j];
        printf("\t%d\t%d\t%d\t%d\t%d\t%f\n", cr->left, cr->right, cr->children[0],
                cr->children[1], cr->node, cr->time);
    }
    printf("Memory heaps\n");
    printf("avl_node_heap:");
    object_heap_print_state(&self->avl_node_heap);
    printf("segment_heap:");
    object_heap_print_state(&self->segment_heap);
    printf("node_mapping_heap:");
    object_heap_print_state(&self->node_mapping_heap);
    msp_verify(self);
out:
    if (ancestors != NULL) {
        free(ancestors);
    }
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
        uint32_t child1, uint32_t child2, uint32_t node)
{
    int ret = 0;
    uint32_t c1 = GSL_MIN(child1, child2);
    uint32_t c2 = GSL_MAX(child1, child2);
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
    cr = &self->coalescence_records[self->num_coalescence_records];
    if (self->num_coalescence_records != 0) {
        lcr = &self->coalescence_records[self->num_coalescence_records - 1];
        if (lcr->right == left
                && lcr->children[0] == c1
                && lcr->children[1] == c2
                && lcr->node == node) {
            /* squash this record into the last */
            lcr->right = right;
            cr = NULL;
        }
    }
    if (cr != NULL) {
        cr->left = left;
        cr->right = right;
        cr->children[0] = c1;
        cr->children[1] = c2;
        cr->node = node;
        cr->time = self->time;
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
msp_recombination_event(msp_t *self)
{
    int ret = 0;
    int64_t l, t, gap, k;
    size_t segment_index;
    node_mapping_t search;
    segment_t *x, *y, *z;
    int64_t num_links = fenwick_get_total(&self->links);

    self->num_re_events++;
    /* We can't use the GSL integer generator here as the range is too large */
    l = 1 + (int64_t) (gsl_rng_uniform(self->rng) * (double) num_links);
    assert(l > 0 && l <= num_links);
    segment_index = fenwick_find(&self->links, l);
    t = fenwick_get_cumulative_sum(&self->links, segment_index);
    gap = t - l;
    assert(gap >= 0 && gap < self->num_loci);
    y = msp_get_segment(self, segment_index);
    x = y->prev;
    k = y->right - gap - 1;
    assert(k >= 0 && k < self->num_loci);
    if (y->left < k) {
        z = msp_alloc_segment(self, (uint32_t) k, y->right, y->value, NULL, y->next);
        if (z == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        if (y->next != NULL) {
            y->next->prev = z;
        }
        y->next = NULL;
        y->right = (uint32_t) k;
        fenwick_increment(&self->links, y->index, k - z->right);
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
    fenwick_set_value(&self->links, z->index, z->right - z->left - 1);
    ret = msp_insert_individual(self, z);
out:
    return ret;
}

static int WARN_UNUSED
msp_coancestry_event(msp_t *self)
{
    int ret = 0;
    int coalescence = 0;
    int defrag_required = 0;
    double covered_fraction;
    uint32_t j, n, l, r, l_min, r_max, v;
    avl_node_t *node;
    node_mapping_t *nm, search;
    segment_t *x, *y, *z, *alpha, *beta;

    self->num_ca_events++;
    /* Choose x and y */
    n = avl_count(&self->ancestral_population);
    j = (uint32_t) gsl_rng_uniform_int(self->rng, n);
    node = avl_at(&self->ancestral_population, j);
    assert(node != NULL);
    x = (segment_t *) node->item;
    avl_unlink_node(&self->ancestral_population, node);
    msp_free_avl_node(self, node);
    j = (uint32_t) gsl_rng_uniform_int(self->rng, n - 1);
    node = avl_at(&self->ancestral_population, j);
    assert(node != NULL);
    y = (segment_t *) node->item;
    avl_unlink_node(&self->ancestral_population, node);
    msp_free_avl_node(self, node);

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
                        NULL, NULL);
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
                    alpha = msp_alloc_segment(self, l, r, v, NULL, NULL);
                    if (alpha == NULL) {
                        ret = MSP_ERR_NO_MEMORY;
                        goto out;
                    }
                }
                ret = msp_record_coalescence(self, l, r, x->value, y->value, v);
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
                fenwick_set_value(&self->links, alpha->index,
                        alpha->right - alpha->left - 1);
            } else {
                defrag_required |= z->right == alpha->left && z->value == alpha->value;
                z->next = alpha;
                fenwick_set_value(&self->links, alpha->index, alpha->right - z->right);
            }
            alpha->prev = z;
            z = alpha;
        }
    }
    if (defrag_required) {
        y = z;
        while (y->prev != NULL) {
            x = y->prev;
            if (x->right == y->left && x->value == y->value) {
                x->right = y->right;
                x->next = y->next;
                if (y->next != NULL) {
                    y->next->prev = x;
                }
                fenwick_increment(&self->links, x->index, y->right - y->left);
                msp_free_segment(self, y);
            }
            y = x;
        }
    }
    if (coalescence) {
        covered_fraction = (r_max - l_min) / (double) self->num_loci;
        /* This is a heuristic to prevent us spending a lot of time pointlessly
         * trying to defragment during the early stages of the simulation.
         * 5% of the overall length seems like a good value and leads to
         * a ~15% time reduction when doing large simulations.
         */
        if (covered_fraction < 0.05) {
            ret = msp_compress_overlap_counts(self, l_min, r_max);
            if (ret != 0) {
                goto out;
            }
        }
    }
out:
    return ret;
}

/* Given the value of the current_population_model, return the next one,
 * or NULL if none exists.
 */
static population_model_t *
msp_get_next_population_model(msp_t *self)
{
    population_model_t *ret = NULL;

    assert(self->current_population_model < self->num_population_models);
    if (self->current_population_model < self->num_population_models - 1) {
        ret = &self->population_models[self->current_population_model + 1];
    }
    return ret;
}

/*
 * Sets up the memory heaps, initial population and trees.
 */
static int WARN_UNUSED
msp_initialise(msp_t *self)
{
    int ret = -1;
    segment_t *u;
    unsigned int j;
    /* These should really be proper checks with a return value */
    assert(self->sample_size > 1);
    assert(self->num_loci >= 1);

    ret = msp_alloc_memory_blocks(self);
    if (ret != 0) {
        goto out;
    }
    /* Set up the initial segments and algorithm state */
    for (j = 1; j <= self->sample_size; j++) {
        u = msp_alloc_segment(self, 0, self->num_loci, j, NULL, NULL);
        if (u == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        ret = msp_insert_individual(self, u);
        if (ret != 0) {
            goto out;
        }
        fenwick_set_value(&self->links, u->index, self->num_loci - 1);
    }
    self->next_node = self->sample_size + 1;
    ret = msp_insert_overlap_count(self, 0, self->sample_size);
    if (ret != 0) {
        goto out;
    }
    ret = msp_insert_overlap_count(self, self->num_loci,
            self->sample_size + 1);
    if (ret != 0) {
        goto out;
    }
    self->current_population_model = 0;
    self->time = 0.0;
    self->num_coalescence_records = 0;
out:
    return ret;
}


int WARN_UNUSED
msp_run(msp_t *self, double max_time, unsigned long max_events)
{
    int ret = 0;
    double lambda_c, lambda_r, t_c, t_r, t_wait;
    int64_t num_links;
    uint32_t n;
    int (*event_method)(msp_t *);
    population_model_t *model, *next_model;
    unsigned long events = 0;

    if (!msp_is_initialised(self)) {
        ret = msp_initialise(self);
        if (ret != 0) {
            goto out;
        }
    }
    model = &self->population_models[self->current_population_model];
    assert(model != NULL);
    next_model = msp_get_next_population_model(self);
    n = avl_count(&self->ancestral_population);
    while (n > 1 && self->time <= max_time && events < max_events) {
        events++;
        num_links = fenwick_get_total(&self->links);
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
        lambda_r = (double) num_links * self->scaled_recombination_rate;
        t_r = DBL_MAX;
        if (lambda_r != 0.0) {
            t_r = gsl_ran_exponential(self->rng, 1.0 / lambda_r);
        }
        /* Need to perform n * (n - 1) as a double due to overflow */
        lambda_c = n * ((double) n - 1.0);
        t_c = model->get_waiting_time(model, lambda_c, self->time, self->rng);
        if (t_r < t_c) {
            t_wait = t_r;
            event_method = msp_recombination_event;
        } else {
            t_wait = t_c;
            event_method = msp_coancestry_event;
        }
        /* TODO check for infinite waiting time and return error. */
        assert(t_wait != DBL_MAX);
        if (next_model != NULL
                && self->time + t_wait >= next_model->start_time) {
            /* We skip ahead to the start time for the next demographic
             * model. The initial size of the population for the next
             * model is the size at this time for the current model.
             */
            self->time = next_model->start_time;
            next_model->initial_size = model->get_size(model, self->time);
            self->current_population_model++;
            model = next_model;
            assert(model != NULL);
            next_model = msp_get_next_population_model(self);
        } else {
            self->time += t_wait;
            ret = event_method(self);
            if (ret != 0) {
                goto out;
            }
            n = avl_count(&self->ancestral_population);
        }
    }
    if (n != 0) {
        ret = 1;
        if (self->time > max_time) {
            ret = 2;
        }
    }
out:
    return ret;
}

/* Returns the number of population models added by the user. */
size_t
msp_get_num_population_models(msp_t *self)
{
    return self->num_population_models - 1;
}

size_t
msp_get_num_ancestors(msp_t *self)
{
    return avl_count(&self->ancestral_population);
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

int
msp_get_ancestors(msp_t *self, segment_t **ancestors)
{
    int ret = -1;
    avl_node_t *node;
    size_t j = 0;

    for (node = (&self->ancestral_population)->head; node != NULL; node = node->next) {
        ancestors[j] = (segment_t *) node->item;
        j++;
    }
    ret = 0;
    return ret;
}

int
msp_get_breakpoints(msp_t *self, uint32_t *breakpoints)
{
    int ret = -1;
    avl_node_t *node;
    node_mapping_t *nm;
    size_t j = 0;

    for (node = (&self->breakpoints)->head; node != NULL; node = node->next) {
        nm = (node_mapping_t *) node->item;
        breakpoints[j] = nm->left;
        j++;
    }
    ret = 0;
    return ret;
}

int
msp_get_coalescence_records(msp_t *self, coalescence_record_t *coalescence_records)
{
    memcpy(coalescence_records, self->coalescence_records,
            self->num_coalescence_records * sizeof(coalescence_record_t));
    return 0;
}

int
msp_get_population_models(msp_t *self, population_model_t *population_models)
{
    /* We do not send back the base population model */
    memcpy(population_models, self->population_models + 1,
            (self->num_population_models - 1) * sizeof(population_model_t));
    return 0;
}
