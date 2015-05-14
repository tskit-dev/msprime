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

#include <gsl/gsl_rng.h>

#include "err.h"
#include "avl.h"
#include "fenwick.h"

#define POP_MODEL_CONSTANT 0
#define POP_MODEL_EXPONENTIAL 1

#define MSP_ALL_BREAKPOINTS 1
#define MSP_ZLIB_COMPRESSION 1
#define MSP_FILE_FORMAT_VERSION 1

#define MAX_BRANCH_LENGTH_STRING 24

/* Using a size_t for index allows us to have an effectively unlimited
 * number of segments. However, we end up wasting 4 bytes of space 
 * per segment because of alignments requirements. This means that 
 * we use 40 bytes instead of 32 (if we use a 32 bit index for a limit
 * of 4G segments), which is a 25% increase in space. It may be possible
 * to do something clever using the offsets of the pointers from the 
 * base address of the memory chunk, which might allow us to get this
 * memory back.
 */
typedef struct segment_t_t {
    uint32_t left;
    uint32_t right;
    uint32_t value;
    size_t index;
    struct segment_t_t *prev;
    struct segment_t_t *next;
} segment_t;

typedef struct {
    uint32_t left;
    uint32_t right;
    uint32_t node;
    double time;
    uint32_t children[2];
} coalescence_record_t;

typedef struct {
    uint32_t left; /* TODO CHANGE THIS - not a good name! */
    uint32_t value;
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

typedef struct {
    size_t object_size;
    size_t block_size; /* number of objects in a block */
    size_t top;
    size_t size;
    size_t num_blocks;
    void **heap;
    char **mem_blocks;
    void (*init_object)(void **obj, size_t index);
} object_heap_t;

typedef struct {
    /* input parameters */
    uint32_t sample_size;
    uint32_t num_loci;
    double scaled_recombination_rate;
    unsigned long random_seed;
    /* allocation block sizes */
    size_t avl_node_block_size;
    size_t node_mapping_block_size;
    size_t segment_block_size;
    size_t max_memory;
    /* population models */
    population_model_t *population_models;
    population_model_t *current_population_model;
    /* Counters for statistics */
    uint64_t num_re_events;
    uint64_t num_ca_events;
    uint64_t num_trapped_re_events;
    /* state */
    size_t used_memory;
    double time;
    uint32_t next_node;
    gsl_rng *rng;
    avl_tree_t ancestral_population;
    avl_tree_t breakpoints;
    fenwick_t links;
    /* memory management */
    object_heap_t avl_node_heap;
    object_heap_t segment_heap;
    /* node mappings are never freed, so simpler requirements */
    void **node_mapping_blocks;
    size_t num_node_mapping_blocks;
    size_t next_node_mapping;
    /* coalescence records are stored in an array */
    coalescence_record_t *coalescence_records;
    size_t num_coalescence_records;
    size_t max_coalescence_records;
    size_t coalescence_record_block_size;
    size_t num_coalescence_record_blocks;
} msp_t;

typedef struct {
    uint32_t sample_size;
    uint32_t num_loci;
    uint32_t *left;
    uint32_t *right;
    uint32_t *node;
    double *time;
    uint32_t *children;
    uint32_t *breakpoints;
    size_t num_records;
    size_t num_breakpoints;
} tree_sequence_t;

typedef struct tree_node {
    uint32_t id;
    uint32_t children[2];
    double time;
    struct tree_node *next;
} tree_node_t;

typedef struct {
    uint32_t key;
    tree_node_t *head;
    tree_node_t *tail;
} tree_node_list_t;

typedef struct {
    tree_sequence_t *tree_sequence;
    uint32_t current_left;
    uint32_t next_breakpoint;
    size_t current_breakpoint_index;
    size_t next_record_index;
    size_t num_records;
    int flags;
    tree_node_list_t nodes_in;
    avl_tree_t active_nodes;
    object_heap_t tree_node_heap;
    object_heap_t tree_node_list_heap;
    object_heap_t avl_node_heap;
} tree_diff_iterator_t;


typedef struct newick_tree_node {
    uint32_t id;
    double time;
    struct newick_tree_node *parent;
    struct newick_tree_node *children[2];
    char branch_length[MAX_BRANCH_LENGTH_STRING];
    char *subtree;
} newick_tree_node_t;

typedef struct {
    uint32_t sample_size;
    uint32_t num_loci;
    size_t precision;
    newick_tree_node_t *root;
    tree_diff_iterator_t diff_iterator;
    avl_tree_t tree;
    object_heap_t avl_node_heap;
} newick_converter_t;

typedef struct hapgen_tree_node {
    uint32_t id;
    double time;
    struct hapgen_tree_node *parent;
    struct hapgen_tree_node *children[2];
    double branch_length;
} hapgen_tree_node_t;

typedef struct {
    uint32_t sample_size;
    uint32_t num_loci;
    unsigned long random_seed;
    double mutation_rate;
    size_t max_haplotype_length;
    hapgen_tree_node_t *root;
    gsl_rng *rng;
    size_t num_segregating_sites;
    double total_branch_length;
    tree_diff_iterator_t diff_iterator;
    avl_tree_t tree;
    object_heap_t avl_node_heap;
    char **haplotypes;
    char *haplotype_mem;
    hapgen_tree_node_t **traversal_stack;
} hapgen_t;

int msp_alloc(msp_t *self);
int msp_add_constant_population_model(msp_t *self, double time, double size);
int msp_add_exponential_population_model(msp_t *self, double time, double alpha);
int msp_initialise(msp_t *self);
int msp_run(msp_t *self, double max_time, unsigned long max_events);
int msp_print_state(msp_t *self);
int msp_free(msp_t *self);

int msp_get_ancestors(msp_t *self, segment_t **ancestors);
int msp_get_breakpoints(msp_t *self, uint32_t *breakpoints);
int msp_get_coalescence_records(msp_t *self, coalescence_record_t *records);

size_t msp_get_num_ancestors(msp_t *self);
size_t msp_get_num_breakpoints(msp_t *self);
size_t msp_get_num_coalescence_records(msp_t *self);
size_t msp_get_num_avl_node_blocks(msp_t *self);
size_t msp_get_num_node_mapping_blocks(msp_t *self);
size_t msp_get_num_segment_blocks(msp_t *self);
size_t msp_get_num_coalescence_record_blocks(msp_t *self);
size_t msp_get_used_memory(msp_t *self);

int tree_sequence_create(tree_sequence_t *self, msp_t *sim);
int tree_sequence_load(tree_sequence_t *self, const char *filename);
int tree_sequence_free(tree_sequence_t *self);
int tree_sequence_dump(tree_sequence_t *self, const char *filename, 
        int flags);
size_t tree_sequence_get_num_breakpoints(tree_sequence_t *self);
size_t tree_sequence_get_num_coalescence_records(tree_sequence_t *self);
int tree_sequence_get_record(tree_sequence_t *self, size_t index, 
        coalescence_record_t *record);
int tree_sequence_get_breakpoints(tree_sequence_t *self, 
        uint32_t *breakpoints);

int tree_diff_iterator_alloc(tree_diff_iterator_t *self, 
        tree_sequence_t *tree_sequence, int flags);
int tree_diff_iterator_free(tree_diff_iterator_t *self);
int tree_diff_iterator_next(tree_diff_iterator_t *self, uint32_t *length,
        tree_node_t **nodes_out, tree_node_t **nodes_in);
void tree_diff_iterator_print_state(tree_diff_iterator_t *self);

int newick_converter_alloc(newick_converter_t *self, 
        tree_sequence_t *tree_sequence, size_t precision, 
        int all_breakpoints);
int newick_converter_next(newick_converter_t *self, uint32_t *length, 
        char **tree);
int newick_converter_free(newick_converter_t *self);
void newick_converter_print_state(newick_converter_t *self);

int hapgen_alloc(hapgen_t *self, tree_sequence_t *tree_sequence,
        double mutation_rate, unsigned long random_seed, 
        size_t max_haplotype_length);
int hapgen_generate(hapgen_t *self);
int hapgen_get_haplotype(hapgen_t *self, size_t j, char **haplotype);
int hapgen_free(hapgen_t *self);
void hapgen_print_state(hapgen_t *self);

const char * msp_strerror(int err);
const char * msp_gsl_version(void);
#endif /*__MSPRIME_H__*/
