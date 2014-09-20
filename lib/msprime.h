/*
** Copyright (C) 2014 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

#ifndef __MSPRIME_H__
#define __MSPRIME_H__

#include <stdio.h>

#include <avl.h>
#include <gsl/gsl_rng.h>

#include "fenwick.h"

#define POP_MODEL_CONSTANT 0
#define POP_MODEL_EXPONENTIAL 1

/* 2^32 * 32 bytes gives a maximum of 128 GiB of segments */
typedef struct segment_t_t {
    unsigned int left;
    unsigned int right;
    unsigned int value;
    unsigned int index;
    struct segment_t_t *prev;
    struct segment_t_t *next;
} segment_t;

/* int based oriented forests gives a maximum sample size of ~10^9 */
typedef struct {
    unsigned int left;
    unsigned int right;
    int children[2];
    int parent;
    float time;
} coalescence_record_t;

typedef struct {
    unsigned int left; /* TODO CHANGE THIS - not a good name! */
    int value;
} node_mapping_t;

typedef struct population_model_t_t {
    int type;
    double start_time;
    double initial_size;
    double param;
    double (*get_size)(struct population_model_t_t *, double);
    double (*get_waiting_time)(struct population_model_t_t *, double, double,
            gsl_rng*);
    struct population_model_t_t *next;
} population_model_t;

typedef struct memory_block_t_t {
    void *mem;
    struct memory_block_t_t *next;
} memory_block_t;

typedef struct {
    size_t object_size;
    size_t block_size; /* number of objects in a block */
    size_t top;
    size_t size;
    void **heap;
    size_t num_blocks;
    memory_block_t *mem_blocks;
    void (*init_object)(void **obj, size_t index);
} object_heap_t;

typedef struct {
    /* input parameters */
    int sample_size;
    int num_loci;
    double recombination_rate;
    long random_seed;
    char *coalescence_record_filename;
    size_t avl_node_block_size;
    size_t node_mapping_block_size;
    /* old */
    int max_segments;
    /* population models */
    population_model_t *population_models;
    population_model_t *current_population_model;
    /* Counters for statistics */
    unsigned int num_re_events;
    unsigned int num_ca_events;
    unsigned int num_trapped_re_events;
    unsigned int num_coalescence_records;
    /* state */
    float time;
    gsl_rng *rng;
    avl_tree_t *ancestral_population;
    avl_tree_t *breakpoints;
    fenwick_t *links;
    FILE *coalescence_record_file;
    /* memory management */
    object_heap_t avl_node_heap;
    object_heap_t segment_heap;
    /* node mappings are never freed so simpler requirements */
    memory_block_t *node_mapping_blocks;
    memory_block_t *current_node_mapping_block;
    size_t next_node_mapping;

    /* old memory heaps */
    /*
    avl_node_t **avl_node_heap;
    int avl_node_heap_top;
    avl_node_t *avl_node_mem;
    segment_t **segment_heap;
    int segment_heap_top;
    segment_t *segment_mem;
    */
} msp_t;

typedef struct {
    int sample_size;
    int num_loci;
    int num_trees;
    int *breakpoints;
    int *pi;
    float *tau;
    char *coalescence_record_filename;
} tree_viewer_t;

int msp_add_constant_population_model(msp_t *self, double time, double size);
int msp_add_exponential_population_model(msp_t *self, double time, double alpha);
int msp_alloc(msp_t *self);
int msp_free(msp_t *self);
int msp_initialise(msp_t *self);
int msp_run(msp_t *self, double max_time, unsigned long max_events);
int msp_print_state(msp_t *self);

tree_viewer_t * msp_get_tree_viewer(msp_t *self);
int tree_viewer_init(tree_viewer_t *self);
int tree_viewer_get_tree(tree_viewer_t *self, int j, int *breakpoint, int **pi,
        float **tau);
int tree_viewer_free(tree_viewer_t *self);

char * msp_strerror(int err);
#endif /*__BIT_H__*/
