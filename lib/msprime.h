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

#ifndef __MSPRIME_H__
#define __MSPRIME_H__

#ifndef MSP_LIBRARY_VERSION_STR
#define MSP_LIBRARY_VERSION_STR "undefined"
#endif

#include <stdio.h>

#include <gsl/gsl_rng.h>

#include "err.h"
#include "avl.h"
#include "fenwick.h"

/* Flags for tree sequence dump/load */
#define MSP_ZLIB_COMPRESSION 1
#define MSP_SKIP_H5CLOSE 2

#define MSP_FILE_FORMAT_VERSION_MAJOR 2
#define MSP_FILE_FORMAT_VERSION_MINOR 0

#define MSP_ORDER_TIME 0
#define MSP_ORDER_LEFT 1
#define MSP_ORDER_RIGHT 2

#define MSP_COUNT_LEAVES 1

#define MAX_BRANCH_LENGTH_STRING 24

/* printf pattern to print out double values losslessly;
 * TODO put this into a private header file somewhere so it's not exported */
#define MSP_LOSSLESS_DBL "%.17g"

/* The root node indicator */
#define MSP_NULL_NODE UINT32_MAX
/* Indicates the that the population ID has not been set. */
#define MSP_NULL_POPULATION_ID UINT8_MAX

typedef struct segment_t_t {
    uint8_t population_id;
    /* During simulation we use genetic coordinates */
    uint32_t left;
    uint32_t right;
    uint32_t value;
    size_t id;
    struct segment_t_t *prev;
    struct segment_t_t *next;
} segment_t;

typedef struct {
    uint8_t population_id;
    /* After simulation, all coordinates are converted to physical coordinates 
     * using a genetic map */
    double left;
    double right;
    uint32_t node;
    double time;
    uint32_t children[2];
} coalescence_record_t;

typedef struct {
    uint32_t left; /* TODO CHANGE THIS - not a good name! */
    uint32_t value;
} node_mapping_t;

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
    uint8_t population_id;
} sample_t;

typedef struct {
    uint32_t sample_size;
    double initial_size;
    double growth_rate;
    double start_time;
    avl_tree_t ancestors;
} population_t;

typedef struct {
    gsl_rng *rng;
    /* input parameters */
    uint32_t sample_size;
    uint32_t num_loci;
    double scaled_recombination_rate;
    uint32_t num_populations;
    sample_t *samples;
    double *initial_migration_matrix;
    population_t *initial_populations;
    /* allocation block sizes */
    size_t avl_node_block_size;
    size_t node_mapping_block_size;
    size_t segment_block_size;
    size_t max_memory;
    /* Counters for statistics */
    size_t num_re_events;
    size_t num_ca_events;
    size_t *num_migration_events;
    size_t num_trapped_re_events;
    size_t num_multiple_re_events;
    /* demographic events */
    struct demographic_event_t_t *demographic_events_head;
    struct demographic_event_t_t *demographic_events_tail;
    struct demographic_event_t_t *next_demographic_event;
    /* algorithm state */
    int state;
    size_t used_memory;
    double time;
    uint32_t next_node;
    double *migration_matrix;
    population_t *populations;
    avl_tree_t breakpoints;
    avl_tree_t overlap_counts;
    fenwick_t links;
    /* memory management */
    object_heap_t avl_node_heap;
    object_heap_t segment_heap;
    object_heap_t node_mapping_heap;
    /* coalescence records are stored in an array */
    coalescence_record_t *coalescence_records;
    size_t num_coalescence_records;
    size_t max_coalescence_records;
    size_t coalescence_record_block_size;
    size_t num_coalescence_record_blocks;
    /* JSON provenance string*/
    char *configuration_json;
} msp_t;

/* Demographic events */
typedef struct {
    int population_id;
    double initial_size;
    double growth_rate;
} population_parameters_change_t;

typedef struct {
    int matrix_index;
    double migration_rate;
} migration_rate_change_t;

typedef struct {
    int source;
    int destination;
    double proportion;
} mass_migration_t;

typedef struct demographic_event_t_t {
    double time;
    int (*change_state)(msp_t *, struct demographic_event_t_t *);
    void (*print_state)(msp_t *, struct demographic_event_t_t *);
    int (*json_snprintf)(struct demographic_event_t_t *, char *, size_t);
    union {
        mass_migration_t mass_migration;
        migration_rate_change_t migration_rate_change;
        population_parameters_change_t population_parameters_change;
    } params;
    struct demographic_event_t_t *next;
} demographic_event_t;

/* Recombination map */

typedef struct {
    uint32_t num_loci;      /* size of the genetic coordinate space  */
    double sequence_length; /* size of the physical coordinate space */
    double total_recombination_rate;
    size_t size;            /* the total number of values in the map */
    double *positions;
    double *rates;
} recomb_map_t;

/* Tree sequences */
typedef struct {
    uint32_t sample_size;
    double sequence_length;
    struct {
        uint8_t *population;
    } samples;
    struct {
        double *left;
        double *right;
        double *time;
        uint8_t *population;
        uint32_t *node;
        uint32_t *children;
        uint32_t *insertion_order;
        uint32_t *removal_order;
        char *parameters;
        char *environment;
    } trees;
    struct {
        uint32_t *node;
        double *position;
        char *parameters;
        char *environment;
    } mutations;
    uint32_t num_nodes;
    size_t num_records;
    size_t num_mutations;
} tree_sequence_t;

typedef struct node_record {
    uint32_t node;
    uint32_t children[2];
    double time;
    struct node_record *next;
} node_record_t;

typedef struct {
    double position;
    uint32_t node;
} mutation_t;

typedef struct leaf_list_node {
    uint32_t node;
    struct leaf_list_node *next;
} leaf_list_node_t;

typedef struct {
    uint32_t sample_size;
    double sequence_length;
    size_t num_nodes;
    size_t num_records;
    double tree_left;
    tree_sequence_t *tree_sequence;
    size_t insertion_index;
    size_t removal_index;
    node_record_t *node_records;
} tree_diff_iterator_t;

typedef struct {
    int flags;
    uint32_t sample_size;
    size_t num_nodes;
    uint32_t root;
    double left;
    double right;
    uint8_t *population;
    uint32_t *parent;
    uint32_t *children;
    double *time;
    uint32_t index;
    /* These are involved in the optional leaf tracking; num_leaves counts
     * all leaves below a give node, and num_tracked_leaves counts those 
     * from a specific subset. */
    uint32_t *num_leaves;
    uint32_t *num_tracked_leaves;
    leaf_list_node_t **leaf_list_head;
    leaf_list_node_t **leaf_list_tail;
    leaf_list_node_t *leaf_list_node_mem;
    /* traversal stacks */
    uint32_t *stack1;
    uint32_t *stack2;
    /* mutation storage */
    mutation_t *mutations;
    size_t num_mutations;
    size_t max_mutations;
} sparse_tree_t;

typedef struct {
    uint32_t sample_size;
    double sequence_length;
    size_t num_nodes;
    size_t num_records;
    tree_sequence_t *tree_sequence;
    sparse_tree_t *tree;
    size_t insertion_index;
    size_t removal_index;
    size_t mutation_index;
} sparse_tree_iterator_t;

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
    double sequence_length;
    size_t precision;
    double Ne;
    newick_tree_node_t *root;
    tree_diff_iterator_t diff_iterator;
    avl_tree_t tree;
    object_heap_t avl_node_heap;
} newick_converter_t;

typedef struct {
    uint32_t sample_size;
    double sequence_length;
    size_t num_mutations;
    tree_sequence_t *tree_sequence;
    /* the haplotype binary matrix */
    size_t words_per_row;
    uint64_t *haplotype_matrix;
    char *haplotype;
    sparse_tree_t tree;
    sparse_tree_iterator_t tree_iterator;
} hapgen_t;

typedef struct {
    uint32_t sample_size;
    double sequence_length;
    size_t num_mutations;
    tree_sequence_t *tree_sequence;
    char *variant;
    size_t tree_mutation_index;
    sparse_tree_t tree;
    sparse_tree_iterator_t tree_iterator;
} vargen_t;

typedef struct {
    gsl_rng *rng;
    double mutation_rate;
    double sequence_length;
    tree_sequence_t *tree_sequence;
    double *times;
    size_t num_mutations;
    size_t max_num_mutations;
    size_t mutation_block_size;
    mutation_t *mutations;
    char *parameters;
    char *environment;
} mutgen_t;

int msp_alloc(msp_t *self, size_t sample_size, gsl_rng *rng);
int msp_set_num_loci(msp_t *self, size_t num_loci);
int msp_set_num_populations(msp_t *self, size_t num_populations);
int msp_set_scaled_recombination_rate(msp_t *self, 
        double scaled_recombination_rate);
int msp_set_max_memory(msp_t *self, size_t max_memory);
int msp_set_node_mapping_block_size(msp_t *self, size_t block_size);
int msp_set_segment_block_size(msp_t *self, size_t block_size);
int msp_set_avl_node_block_size(msp_t *self, size_t block_size);
int msp_set_coalescence_record_block_size(msp_t *self, size_t block_size);
int msp_set_sample_configuration(msp_t *self, size_t num_populations,
        size_t *sample_configuration);
int msp_set_migration_matrix(msp_t *self, size_t size,
        double *migration_matrix);
int msp_set_population_configuration(msp_t *self, int population_id, 
        size_t sample_size, double initial_size, double growth_rate);

int msp_add_population_parameters_change(msp_t *self, double time, 
        int population_id, double size, double growth_rate);
int msp_add_migration_rate_change(msp_t *self, double time, int matrix_index,
        double migration_rate);
int msp_add_mass_migration(msp_t *self, double time, int source, int dest,
        double proportion);

int msp_initialise(msp_t *self);
int msp_run(msp_t *self, double max_time, unsigned long max_events);
int msp_debug_demography(msp_t *self, double *end_time);
int msp_reset(msp_t *self);
int msp_print_state(msp_t *self);
int msp_free(msp_t *self);
void msp_verify(msp_t *self);

int msp_get_ancestors(msp_t *self, segment_t **ancestors);
int msp_get_breakpoints(msp_t *self, size_t *breakpoints);
int msp_get_migration_matrix(msp_t *self, double *migration_matrix);
int msp_get_num_migration_events(msp_t *self, size_t *num_migration_events);
int msp_get_coalescence_records(msp_t *self, coalescence_record_t *records);
int msp_get_samples(msp_t *self, sample_t *samples);
int msp_get_population_configuration(msp_t *self, size_t population_id,
        size_t *sample_size, double *initial_size, double *growth_rate);
int msp_get_population(msp_t *self, size_t population_id, 
        population_t **population);
int msp_is_completed(msp_t *self);

size_t msp_get_sample_size(msp_t *self);
size_t msp_get_num_loci(msp_t *self);
size_t msp_get_num_populations(msp_t *self);
size_t msp_get_num_ancestors(msp_t *self);
size_t msp_get_num_breakpoints(msp_t *self);
size_t msp_get_num_coalescence_records(msp_t *self);
size_t msp_get_num_avl_node_blocks(msp_t *self);
size_t msp_get_num_node_mapping_blocks(msp_t *self);
size_t msp_get_num_segment_blocks(msp_t *self);
size_t msp_get_num_coalescence_record_blocks(msp_t *self);
size_t msp_get_used_memory(msp_t *self);
size_t msp_get_num_common_ancestor_events(msp_t *self);
size_t msp_get_num_recombination_events(msp_t *self);

char *msp_get_configuration_json(msp_t *self);

void tree_sequence_print_state(tree_sequence_t *self);
int tree_sequence_create(tree_sequence_t *self, msp_t *sim,
        recomb_map_t *recomb_map, double Ne);
int tree_sequence_load(tree_sequence_t *self, const char *filename, int flags);
int tree_sequence_free(tree_sequence_t *self);
int tree_sequence_dump(tree_sequence_t *self, const char *filename, int flags);
size_t tree_sequence_get_num_breakpoints(tree_sequence_t *self);
size_t tree_sequence_get_num_coalescence_records(tree_sequence_t *self);
size_t tree_sequence_get_num_mutations(tree_sequence_t *self);
uint32_t tree_sequence_get_num_nodes(tree_sequence_t *self);
uint32_t tree_sequence_get_sample_size(tree_sequence_t *self);
double tree_sequence_get_sequence_length(tree_sequence_t *self);

int tree_sequence_get_record(tree_sequence_t *self, size_t index, 
        coalescence_record_t *record, int order);
int tree_sequence_get_breakpoints(tree_sequence_t *self, 
        uint32_t *breakpoints);
int tree_sequence_get_mutations(tree_sequence_t *self, mutation_t *mutations);
int tree_sequence_get_population(tree_sequence_t *self, uint32_t u,
    uint32_t *population_id);
int tree_sequence_get_pairwise_diversity(tree_sequence_t *self,
    uint32_t *samples, uint32_t num_samples, double *pi);
char * tree_sequence_get_simulation_parameters(tree_sequence_t *self);
char * tree_sequence_get_mutation_parameters(tree_sequence_t *self);
int tree_sequence_alloc_sparse_tree(tree_sequence_t *self, 
        sparse_tree_t *tree, uint32_t *tracked_leaves, 
        uint32_t num_tracked_leaves, int flags);
int tree_sequence_set_mutations(tree_sequence_t *self, 
        size_t num_mutations, mutation_t *mutations, 
        const char *parameters, const char *environment);

int tree_diff_iterator_alloc(tree_diff_iterator_t *self, 
        tree_sequence_t *tree_sequence);
int tree_diff_iterator_free(tree_diff_iterator_t *self);
int tree_diff_iterator_next(tree_diff_iterator_t *self, double *length,
        node_record_t **nodes_out, node_record_t **nodes_in);
void tree_diff_iterator_print_state(tree_diff_iterator_t *self);

int sparse_tree_alloc(sparse_tree_t *self, uint32_t sample_size, 
        uint32_t num_nodes, size_t max_mutations, uint32_t *tracked_leaves, 
        uint32_t num_tracked_leaves, int flags);
int sparse_tree_free(sparse_tree_t *self);
int sparse_tree_clear(sparse_tree_t *self);
int sparse_tree_get_parent(sparse_tree_t *self, uint32_t u, uint32_t *parent);
int sparse_tree_get_mrca(sparse_tree_t *self, uint32_t u, uint32_t v,
        uint32_t *mrca);
int sparse_tree_get_num_leaves(sparse_tree_t *self, uint32_t u,
        uint32_t *num_leaves);
int sparse_tree_get_num_tracked_leaves(sparse_tree_t *self, uint32_t u,
        uint32_t *num_tracked_leaves);
int sparse_tree_get_leaf_list(sparse_tree_t *self, uint32_t u,
        leaf_list_node_t **head, leaf_list_node_t **tail);

int sparse_tree_iterator_alloc(sparse_tree_iterator_t *self, 
        tree_sequence_t *tree_sequence, sparse_tree_t *tree);
int sparse_tree_iterator_free(sparse_tree_iterator_t *self);
int sparse_tree_iterator_next(sparse_tree_iterator_t *self);
void sparse_tree_iterator_print_state(sparse_tree_iterator_t *self);

int newick_converter_alloc(newick_converter_t *self, 
        tree_sequence_t *tree_sequence, size_t precision, double Ne);
int newick_converter_next(newick_converter_t *self, double *length, 
        char **tree);
int newick_converter_free(newick_converter_t *self);
void newick_converter_print_state(newick_converter_t *self);

int hapgen_alloc(hapgen_t *self, tree_sequence_t *tree_sequence);
int hapgen_get_haplotype(hapgen_t *self, uint32_t j, char **haplotype);
size_t hapgen_get_num_segregating_sites(hapgen_t *self);
int hapgen_free(hapgen_t *self);
void hapgen_print_state(hapgen_t *self);

int vargen_alloc(vargen_t *self, tree_sequence_t *tree_sequence);
int vargen_next(vargen_t *self, double *position, char **variant);
int vargen_free(vargen_t *self);
void vargen_print_state(vargen_t *self);

int recomb_map_alloc(recomb_map_t *self, uint32_t num_loci, 
        double sequence_length, double *positions, double *rates, 
        size_t size);
int recomb_map_free(recomb_map_t *self);
uint32_t recomb_map_get_num_loci(recomb_map_t *self);
double recomb_map_get_sequence_length(recomb_map_t *self);
double recomb_map_get_per_locus_recombination_rate(recomb_map_t *self);
double recomb_map_get_total_recombination_rate(recomb_map_t *self);
double recomb_map_genetic_to_phys(recomb_map_t *self, double genetic_x);
double recomb_map_phys_to_genetic(recomb_map_t *self, double phys_x);
int recomb_map_genetic_to_phys_bulk(recomb_map_t *self, double *genetic_x, 
        size_t n);
size_t recomb_map_get_size(recomb_map_t *self);
int recomb_map_get_positions(recomb_map_t *self, double *positions);
int recomb_map_get_rates(recomb_map_t *self, double *rates);

void recomb_map_print_state(recomb_map_t *self);

int mutgen_alloc(mutgen_t *self, tree_sequence_t *tree_sequence, 
        double mutation_rate, gsl_rng *rng);
int mutgen_free(mutgen_t *self);
int mutgen_generate(mutgen_t *self);
void mutgen_print_state(mutgen_t *self);

int msp_encode_environment(char **destination);
const char * msp_strerror(int err);

#endif /*__MSPRIME_H__*/
