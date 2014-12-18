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

#include "avl.h"
#include "fenwick.h"

#define POP_MODEL_CONSTANT 0
#define POP_MODEL_EXPONENTIAL 1

/* index: 2^32 * 32 bytes gives a maximum of 128 GiB of segments */
/* TODO should we change index to a size_t? This will make running out
 * of segment space irrelevant. However, we get an 8 byte penalty for 
 * doing so, because of the struct alignment requirements. */
typedef struct segment_t_t {
    uint32_t left;
    uint32_t right;
    uint32_t index;
    int32_t value;
    struct segment_t_t *prev;
    struct segment_t_t *next;
} segment_t;

/* int based oriented forests gives a maximum sample size of ~10^9 */
typedef struct {
    uint32_t left;
    uint32_t right;
    int32_t children[2];
    int32_t parent;
    float time;
} coalescence_record_t;

typedef struct {
    uint32_t left; /* TODO CHANGE THIS - not a good name! */
    int32_t value;
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
    void **mem_blocks;
    void (*init_object)(void **obj, size_t index);
} object_heap_t;

typedef struct {
    char *filename;
    char mode;
    FILE *file;
    uint32_t sample_size;
    uint32_t num_loci;
    uint32_t flags;
    char *metadata;
    size_t coalescence_record_offset;
    size_t metadata_offset;
} tree_file_t;

typedef struct {
    double mutation_rate;
    long random_seed;
    size_t max_haplotype_length;
    size_t haplotype_length;
    size_t num_trees;
    uint32_t sample_size;
    uint32_t num_loci;
    tree_file_t tree_file;
    gsl_rng *rng;
    char **haplotypes;
    char *haplotype_mem;
    int *pi;
    float *tau;
    /* tree traversal memory */
    int *child;
    int *sib;
    int *branch_mutations;
    int *mutation_sites;
} hapgen_t;

typedef struct {
    uint32_t sample_size;
    uint32_t num_loci;
    float *tau;
    float *branch_lengths;
    int **children;
    int *visited;
    int *stack;
    size_t stack_size;
    char *output_buffer;
    size_t output_buffer_size;
    int *children_mem;
    tree_file_t tree_file;
    uint32_t breakpoint;
    size_t num_trees;
    coalescence_record_t next_record;
    int completed;
} newick_t;


typedef struct {
    /* input parameters */
    /* TODO change these to uint32_t and run compiler with paranoid typing */
    int sample_size;
    int num_loci;
    double scaled_recombination_rate;
    long random_seed;
    char *tree_file_name;
    /* allocation block sizes */
    size_t avl_node_block_size;
    size_t node_mapping_block_size;
    size_t segment_block_size;
    size_t max_memory;
    /* population models */
    population_model_t *population_models;
    population_model_t *current_population_model;
    /* Counters for statistics */
    uint32_t num_re_events;
    uint32_t num_ca_events;
    uint32_t num_trapped_re_events;
    uint32_t num_coalescence_records;
    /* state */
    size_t used_memory;
    float time;
    gsl_rng *rng;
    avl_tree_t ancestral_population;
    avl_tree_t breakpoints;
    fenwick_t links;
    tree_file_t tree_file;
    /* memory management */
    object_heap_t avl_node_heap;
    object_heap_t segment_heap;
    /* node mappings are never freed, so simpler requirements */
    void **node_mapping_blocks;
    size_t num_node_mapping_blocks;
    size_t next_node_mapping;
    /* last coalescence record to enable squashing of adjacent records */
    coalescence_record_t last_coalesence_record;
} msp_t;

int msp_alloc(msp_t *self);
int msp_add_constant_population_model(msp_t *self, double time, double size);
int msp_add_exponential_population_model(msp_t *self, double time, double alpha);
int msp_initialise(msp_t *self);
int msp_run(msp_t *self, double max_time, unsigned long max_events);
int msp_get_ancestors(msp_t *self, segment_t **);
int msp_print_state(msp_t *self);
int msp_write_metadata(msp_t *self, FILE *f);
int msp_free(msp_t *self);
size_t msp_get_num_ancestors(msp_t *self);
size_t msp_get_num_breakpoints(msp_t *self);
size_t msp_get_num_avl_node_blocks(msp_t *self);
size_t msp_get_num_node_mapping_blocks(msp_t *self);
size_t msp_get_num_segment_blocks(msp_t *self);
size_t msp_get_used_memory(msp_t *self);

int tree_file_open(tree_file_t *self, const char *tree_file_name, char mode);
int tree_file_set_sample_size(tree_file_t *self, uint32_t sample_size);
int tree_file_set_num_loci(tree_file_t *self, uint32_t num_loci);
int tree_file_close(tree_file_t *self);
int tree_file_sort(tree_file_t *self);
int tree_file_finalise(tree_file_t *self, msp_t *msp);
int tree_file_append_record(tree_file_t *self, coalescence_record_t *r);
int tree_file_next_record(tree_file_t *self, coalescence_record_t *r);
int tree_file_print_state(tree_file_t *self);
int tree_file_print_records(tree_file_t *self);
int tree_file_iscomplete(tree_file_t *self);
int tree_file_issorted(tree_file_t *self);
int tree_file_isopen(tree_file_t *self);

int hapgen_alloc(hapgen_t *self, double mutation_rate, 
        const char *tree_file_name, long random_seed, 
        size_t max_haplotype_length);
int hapgen_generate(hapgen_t *self);
int hapgen_get_haplotypes(hapgen_t *self, char ***haplotypes, size_t *s);
int hapgen_free(hapgen_t *self);

int newick_alloc(newick_t *self, const char *tree_file_name);
int newick_next_tree(newick_t *self, uint32_t *l, char **tree);
int newick_free(newick_t *self);

char * msp_strerror(int err);
#endif /*__MSPRIME_H__*/
