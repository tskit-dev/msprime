/*
** Copyright (C) 2015-2017 University of Oxford
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
#include <stdbool.h>

#include <gsl/gsl_rng.h>

#include "err.h"
#include "avl.h"
#include "fenwick.h"

/* Flags for tree sequence dump/load */
#define MSP_DUMP_ZLIB_COMPRESSION 1
#define MSP_LOAD_EXTENDED_CHECKS  1

#define MSP_FILE_FORMAT_VERSION_MAJOR 6
#define MSP_FILE_FORMAT_VERSION_MINOR 0

/* Flags for simplify() */
#define MSP_FILTER_INVARIANT_SITES 1

#define MSP_LEAF_COUNTS  1
#define MSP_LEAF_LISTS   2

#define MSP_DIR_FORWARD 1
#define MSP_DIR_REVERSE -1

#define MSP_GENOTYPES_AS_CHAR 1

#define MSP_ALPHABET_BINARY 0
#define MSP_ALPHABET_ASCII  1

#define MSP_MODEL_HUDSON 0
#define MSP_MODEL_SMC 1
#define MSP_MODEL_SMC_PRIME 2
#define MSP_MODEL_BETA 3
#define MSP_MODEL_DIRAC 4

#define MSP_NODE_IS_SAMPLE 1

#define MAX_BRANCH_LENGTH_STRING 24

/* The root node indicator */
#define MSP_NULL_NODE (-1)
/* Indicates the that the population ID has not been set. */
#define MSP_NULL_POPULATION_ID (-1)

#define MSP_INITIALISED_MAGIC 0x1234567

typedef int32_t node_id_t;
typedef int32_t population_id_t;
typedef int32_t site_id_t;
typedef int32_t mutation_id_t;
typedef uint32_t list_len_t;

typedef struct {
    size_t num_rows;
    size_t max_rows;
    size_t max_rows_increment;
    size_t total_ancestral_state_length;
    size_t max_total_ancestral_state_length;
    size_t max_total_ancestral_state_length_increment;
    char *ancestral_state;
    list_len_t *ancestral_state_length;
    double *position;
} site_table_t;

typedef struct {
    size_t num_rows;
    size_t max_rows;
    size_t max_rows_increment;
    size_t total_derived_state_length;
    size_t max_total_derived_state_length;
    size_t max_total_derived_state_length_increment;
    node_id_t *node;
    site_id_t *site;
    char *derived_state;
    list_len_t *derived_state_length;
} mutation_table_t;

typedef struct {
    size_t num_rows;
    size_t max_rows;
    size_t max_rows_increment;
    size_t total_name_length;
    size_t max_total_name_length;
    size_t max_total_name_length_increment;
    uint32_t *flags;
    double *time;
    population_id_t *population;
    char *name;
    list_len_t *name_length;
} node_table_t;

typedef struct {
    size_t num_rows;
    size_t max_rows;
    size_t max_rows_increment;
    size_t total_children_length;
    size_t max_total_children_length;
    size_t max_total_children_length_increment;
    double *left;
    double *right;
    node_id_t *parent;
    node_id_t *children;
    list_len_t *children_length;
} edgeset_table_t;

typedef struct {
    size_t num_rows;
    size_t max_rows;
    size_t max_rows_increment;
    population_id_t *source;
    population_id_t *dest;
    node_id_t *node;
    double *left;
    double *right;
    double *time;
} migration_table_t;

typedef struct segment_t_t {
    population_id_t population_id;
    /* During simulation we use genetic coordinates */
    uint32_t left;
    uint32_t right;
    node_id_t value;
    size_t id;
    struct segment_t_t *prev;
    struct segment_t_t *next;
} segment_t;

typedef struct {
    population_id_t population_id;
    uint32_t num_children;
    /* After simulation, all coordinates are converted to physical coordinates
     * using a genetic map */
    node_id_t node;
    double left;
    double right;
    double time;
    node_id_t *children;
} coalescence_record_t;

/* TODO remove migration record --- this is identical to the migration_t type.
 * Need to update the simulator code to use the new types and remove the
 * *_record types completely
 */
typedef struct {
    population_id_t source;
    population_id_t dest;
    node_id_t node;
    double left;
    double right;
    double time;
} migration_t;

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
    population_id_t population_id;
    double time;
} sample_t;

typedef struct {
    double initial_size;
    double growth_rate;
    double start_time;
    avl_tree_t ancestors;
} population_t;

typedef struct {
    double time;
    node_id_t sample;
    population_id_t population_id;
} sampling_event_t;

/* Simulation models */

typedef struct {
    double alpha;
    double truncation_point;
} beta_coalescent_t;

typedef struct {
    double psi;
    double c; // constant
} dirac_coalescent_t;

typedef struct {
    int type;
    union {
        beta_coalescent_t beta_coalescent;
        dirac_coalescent_t dirac_coalescent;
    } params;
} simulation_model_t;

typedef struct {
    gsl_rng *rng;
    /* input parameters */
    simulation_model_t model;
    bool store_migrations;
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
    size_t num_rejected_ca_events;
    size_t *num_migration_events;
    size_t num_trapped_re_events;
    size_t num_multiple_re_events;
    /* sampling events */
    sampling_event_t *sampling_events;
    size_t num_sampling_events;
    size_t next_sampling_event;
    /* Demographic events */
    struct demographic_event_t_t *demographic_events_head;
    struct demographic_event_t_t *demographic_events_tail;
    struct demographic_event_t_t *next_demographic_event;
    /* algorithm state */
    int state;
    size_t used_memory;
    double time;
    node_id_t next_node;
    double *migration_matrix;
    population_t *populations;
    avl_tree_t breakpoints;
    avl_tree_t overlap_counts;
    fenwick_t links;
    /* memory management */
    object_heap_t avl_node_heap;
    object_heap_t segment_heap;
    object_heap_t node_mapping_heap;
    object_heap_t binary_children_heap;
    /* coalescence records are stored in a flat array */
    coalescence_record_t *coalescence_records;
    size_t num_coalescence_records;
    size_t max_coalescence_records;
    size_t coalescence_record_block_size;
    size_t num_coalescence_record_blocks;
    /* migration records are stored in a flat array */
    migration_t *migrations;
    size_t num_migrations;
    size_t max_migrations;
    size_t migration_block_size;
    size_t num_migration_blocks;
} msp_t;

/* Demographic events */
typedef struct {
    population_id_t population_id;
    double initial_size;
    double growth_rate;
} population_parameters_change_t;

typedef struct {
    int matrix_index;
    double migration_rate;
} migration_rate_change_t;

typedef struct {
    population_id_t source;
    population_id_t destination;
    double proportion;
} mass_migration_t;

typedef struct {
    population_id_t population_id;
    double proportion;
} simple_bottleneck_t;

typedef struct {
    population_id_t population_id;
    double strength;
} instantaneous_bottleneck_t;

typedef struct demographic_event_t_t {
    double time;
    int (*change_state)(msp_t *, struct demographic_event_t_t *);
    void (*print_state)(msp_t *, struct demographic_event_t_t *, FILE *out);
    union {
        simple_bottleneck_t simple_bottleneck;
        instantaneous_bottleneck_t instantaneous_bottleneck;
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

/* Record definitions for tree sequence types. */

typedef struct {
    uint32_t flags;
    double time;
    population_id_t population;
    char *name;
} node_t;

typedef struct {
    mutation_id_t id;
    site_id_t site;
    node_id_t node;
    const char *derived_state;
    list_len_t derived_state_length;
    // TODO remove this and change to ID?
    size_t index;
} mutation_t;

typedef struct {
    site_id_t id;
    double position;
    const char *ancestral_state;
    list_len_t ancestral_state_length;
    mutation_t *mutations;
    list_len_t mutations_length;
} site_t;

typedef struct {
    node_id_t parent;
    node_id_t *children;
    list_len_t children_length;
    double left;
    double right;
    double time;
} edgeset_t;

/* Tree sequences */
typedef struct {
    uint32_t initialised_magic;
    size_t num_trees;
    double sequence_length;
    int alphabet;
    size_t sample_size;
    size_t max_sample_size;
    node_id_t *samples;
    struct {
        size_t num_records;
        size_t max_num_records;
        size_t total_name_length;
        size_t max_total_name_length;
        uint32_t *flags;
        population_id_t *population;
        double *time;
        list_len_t *name_length;
        char **name;
        char *name_mem;
        node_id_t *sample_index_map;
    } nodes;

    struct {
        size_t num_records;
        size_t max_num_records;
        size_t total_children_length;
        size_t max_total_children_length;
        double *left;
        double *right;
        node_id_t *parent;
        list_len_t *children_length;
        node_id_t **children;
        node_id_t *children_mem;
        struct {
            node_id_t *insertion_order;
            node_id_t *removal_order;
        } indexes;
    } edgesets;

    struct {
        size_t num_records;
        size_t max_num_records;
        size_t total_ancestral_state_length;
        size_t max_total_ancestral_state_length;
        char **ancestral_state;
        char *ancestral_state_mem;
        list_len_t *ancestral_state_length;
        double *position;
        site_t *tree_sites_mem;
        site_t **tree_sites;
        list_len_t *tree_sites_length;
        mutation_t *site_mutations_mem;
        mutation_t **site_mutations;
        list_len_t *site_mutations_length;
    } sites;

    struct {
        size_t num_records;
        size_t max_num_records;
        size_t total_derived_state_length;
        size_t max_total_derived_state_length;
        node_id_t *node;
        site_id_t *site;
        char **derived_state;
        char *derived_state_mem;
        list_len_t *derived_state_length;
    } mutations;

    struct {
        size_t num_records;
        size_t max_num_records;
        node_id_t *node;
        population_id_t *source;
        population_id_t *dest;
        double *left;
        double *right;
        double *time;
    } migrations;

    char **provenance_strings;
    size_t num_provenance_strings;
    size_t max_num_provenance_strings;
} tree_sequence_t;

/* TODO rename this struct. This is just used in the tree_diff iterator and
 * can easily be confused with the node_t type.
 */
typedef struct node_record {
    node_id_t node;
    size_t num_children;
    node_id_t *children;
    double time;
    struct node_record *next;
} node_record_t;

typedef struct leaf_list_node {
    node_id_t node;
    struct leaf_list_node *next;
} leaf_list_node_t;

typedef struct {
    size_t sample_size;
    double sequence_length;
    size_t num_nodes;
    size_t num_records;
    double tree_left;
    tree_sequence_t *tree_sequence;
    size_t insertion_index;
    size_t removal_index;
    size_t tree_index;
    node_record_t *node_records;
} tree_diff_iterator_t;

typedef struct {
    tree_sequence_t *tree_sequence;
    size_t sample_size;
    size_t num_nodes;
    int flags;
    node_id_t *samples;
    node_id_t root;
    /* Left and right physical coordinates of the tree */
    double left;
    double right;
    population_id_t *population;
    node_id_t *parent;
    /* TODO Change this to children_length for consistency with the tree sequence */
    list_len_t *num_children;
    node_id_t **children;
    double *time;
    size_t index;
    /* These are involved in the optional leaf tracking; num_leaves counts
     * all leaves below a give node, and num_tracked_leaves counts those
     * from a specific subset. */
    node_id_t *num_leaves;
    node_id_t *num_tracked_leaves;
    /* All nodes that are marked during a particular transition are marked
     * with a given value. */
    uint8_t *marked;
    uint8_t mark;
    /* These are for the optional leaf list tracking. */
    leaf_list_node_t **leaf_list_head;
    leaf_list_node_t **leaf_list_tail;
    leaf_list_node_t *leaf_list_node_mem;
    /* traversal stacks */
    node_id_t *stack1;
    node_id_t *stack2;
    /* The sites on this tree */
    site_t *sites;
    list_len_t sites_length;
    /* Counters needed for next() and prev() transformations. */
    int direction;
    node_id_t left_index;
    node_id_t right_index;
} sparse_tree_t;

typedef struct newick_tree_node {
    node_id_t id;
    double time;
    struct newick_tree_node *parent;
    struct newick_tree_node *children[2];
    char branch_length[MAX_BRANCH_LENGTH_STRING];
    char *subtree;
} newick_tree_node_t;

typedef struct {
    size_t sample_size;
    double sequence_length;
    size_t precision;
    double Ne;
    newick_tree_node_t *root;
    tree_diff_iterator_t diff_iterator;
    avl_tree_t tree;
    object_heap_t avl_node_heap;
} newick_converter_t;

typedef struct {
    size_t sample_size;
    double sequence_length;
    size_t num_sites;
    tree_sequence_t *tree_sequence;
    node_id_t *sample_index_map;
    /* The haplotype binary matrix. This is an optimised special case. */
    bool binary;
    size_t words_per_row;
    uint64_t *binary_haplotype_matrix;
    char *output_haplotype;
    /* The general haplotype matrix. */
    char *ascii_haplotype_matrix;
    sparse_tree_t tree;
} hapgen_t;

typedef struct {
    size_t sample_size;
    double sequence_length;
    size_t num_sites;
    tree_sequence_t *tree_sequence;
    node_id_t *sample_index_map;
    size_t tree_site_index;
    int finished;
    sparse_tree_t tree;
    int flags;
} vargen_t;

typedef struct {
    size_t sample_size;
    size_t num_vcf_samples;
    unsigned int ploidy;
    const char *chrom;
    char *genotypes;
    char *header;
    char *record;
    char *vcf_genotypes;
    size_t vcf_genotypes_size;
    size_t record_size;
    size_t num_sites;
    unsigned long contig_length;
    unsigned long *positions;
    vargen_t *vargen;
} vcf_converter_t;

typedef struct {
    sparse_tree_t *outer_tree;
    sparse_tree_t *inner_tree;
    size_t num_sites;
    int tree_changed;
    tree_sequence_t *tree_sequence;
} ld_calc_t;

typedef struct {
    double position;
    node_id_t node;
    const char *ancestral_state;
    const char *derived_state;
} infinite_sites_mutation_t;

typedef struct {
    int alphabet;
    gsl_rng *rng;
    double mutation_rate;
    size_t num_mutations;
    size_t max_num_mutations;
    size_t mutation_block_size;
    site_table_t *sites;
    infinite_sites_mutation_t *mutations;
    object_heap_t avl_node_heap;
} mutgen_t;

int msp_alloc(msp_t *self, size_t sample_size, sample_t *samples, gsl_rng *rng);
int msp_set_simulation_model_non_parametric(msp_t *self, int model);
int msp_set_simulation_model_dirac(msp_t *self, double psi, double c);
int msp_set_simulation_model_beta(msp_t *self, double alpha, double truncation_point);
int msp_set_num_loci(msp_t *self, size_t num_loci);
int msp_set_store_migrations(msp_t *self, bool store_migrations);
int msp_set_num_populations(msp_t *self, size_t num_populations);
int msp_set_scaled_recombination_rate(msp_t *self,
        double scaled_recombination_rate);
int msp_set_max_memory(msp_t *self, size_t max_memory);
int msp_set_node_mapping_block_size(msp_t *self, size_t block_size);
int msp_set_segment_block_size(msp_t *self, size_t block_size);
int msp_set_avl_node_block_size(msp_t *self, size_t block_size);
int msp_set_coalescence_record_block_size(msp_t *self, size_t block_size);
int msp_set_migration_block_size(msp_t *self, size_t block_size);
int msp_set_sample_configuration(msp_t *self, size_t num_populations,
        size_t *sample_configuration);
int msp_set_migration_matrix(msp_t *self, size_t size,
        double *migration_matrix);
int msp_set_population_configuration(msp_t *self, int population_id,
        double initial_size, double growth_rate);

int msp_add_population_parameters_change(msp_t *self, double time,
        int population_id, double size, double growth_rate);
int msp_add_migration_rate_change(msp_t *self, double time, int matrix_index,
        double migration_rate);
int msp_add_mass_migration(msp_t *self, double time, int source, int dest,
        double proportion);
int msp_add_simple_bottleneck(msp_t *self, double time, int population_id,
        double intensity);
int msp_add_instantaneous_bottleneck(msp_t *self, double time, int population_id,
        double strength);

int msp_initialise(msp_t *self);
int msp_run(msp_t *self, double max_time, unsigned long max_events);
int msp_debug_demography(msp_t *self, double *end_time);
int msp_populate_tables(msp_t *self, double Ne, recomb_map_t *recomb_map,
        node_table_t *node_table, edgeset_table_t *edgeset_table,
        migration_table_t *migration_table);
int msp_reset(msp_t *self);
int msp_print_state(msp_t *self, FILE *out);
int msp_free(msp_t *self);
void msp_verify(msp_t *self);

int msp_get_ancestors(msp_t *self, segment_t **ancestors);
int msp_get_breakpoints(msp_t *self, size_t *breakpoints);
int msp_get_migration_matrix(msp_t *self, double *migration_matrix);
int msp_get_num_migration_events(msp_t *self, size_t *num_migration_events);
int msp_get_coalescence_records(msp_t *self, coalescence_record_t **records);
int msp_get_migrations(msp_t *self, migration_t **records);
int msp_get_samples(msp_t *self, sample_t **samples);
int msp_get_population_configuration(msp_t *self, size_t population_id,
        double *initial_size, double *growth_rate);
int msp_get_population(msp_t *self, size_t population_id,
        population_t **population);
int msp_is_completed(msp_t *self);

simulation_model_t * msp_get_model(msp_t *self);
const char * msp_get_model_name(msp_t *self);
bool msp_get_store_migrations(msp_t *self);
size_t msp_get_sample_size(msp_t *self);
size_t msp_get_num_loci(msp_t *self);
size_t msp_get_num_populations(msp_t *self);
size_t msp_get_num_ancestors(msp_t *self);
size_t msp_get_num_breakpoints(msp_t *self);
size_t msp_get_num_coalescence_records(msp_t *self);
size_t msp_get_num_migrations(msp_t *self);
size_t msp_get_num_avl_node_blocks(msp_t *self);
size_t msp_get_num_node_mapping_blocks(msp_t *self);
size_t msp_get_num_segment_blocks(msp_t *self);
size_t msp_get_num_coalescence_record_blocks(msp_t *self);
size_t msp_get_num_migration_blocks(msp_t *self);
size_t msp_get_used_memory(msp_t *self);
size_t msp_get_num_common_ancestor_events(msp_t *self);
size_t msp_get_num_rejected_common_ancestor_events(msp_t *self);
size_t msp_get_num_recombination_events(msp_t *self);

void tree_sequence_print_state(tree_sequence_t *self, FILE *out);
int tree_sequence_initialise(tree_sequence_t *self);
/* Marking the x_tables API as tmp until we figure out what to do
 * with provenance. */
int tree_sequence_load_tables_tmp(tree_sequence_t *self,
        node_table_t *nodes, edgeset_table_t *edgesets, migration_table_t *migrations,
        site_table_t *sites, mutation_table_t *mutations,
        size_t num_provenance_strings, char **provenance_strings);
int tree_sequence_dump_tables_tmp(tree_sequence_t *self, node_table_t *node_table,
        edgeset_table_t *edgeset_table, migration_table_t *migration_table,
        site_table_t *sites, mutation_table_t *mutations,
        size_t *num_provenance_strings, char ***provenance_strings);
int tree_sequence_load(tree_sequence_t *self, const char *filename, int flags);
int tree_sequence_dump(tree_sequence_t *self, const char *filename, int flags);
int tree_sequence_free(tree_sequence_t *self);

size_t tree_sequence_get_num_nodes(tree_sequence_t *self);
size_t tree_sequence_get_num_migrations(tree_sequence_t *self);
size_t tree_sequence_get_num_edgesets(tree_sequence_t *self);
size_t tree_sequence_get_num_migrations(tree_sequence_t *self);
size_t tree_sequence_get_num_sites(tree_sequence_t *self);
size_t tree_sequence_get_num_mutations(tree_sequence_t *self);
size_t tree_sequence_get_num_trees(tree_sequence_t *self);
size_t tree_sequence_get_sample_size(tree_sequence_t *self);
double tree_sequence_get_sequence_length(tree_sequence_t *self);
int tree_sequence_get_alphabet(tree_sequence_t *self);
bool tree_sequence_is_sample(tree_sequence_t *self, node_id_t u);

int tree_sequence_get_node(tree_sequence_t *self, node_id_t index, node_t *node);
int tree_sequence_get_edgeset(tree_sequence_t *self, size_t index, edgeset_t *edgeset);
int tree_sequence_get_migration(tree_sequence_t *self, size_t index,
        migration_t *migration);
int tree_sequence_get_site(tree_sequence_t *self, site_id_t id, site_t *site);
int tree_sequence_get_mutation(tree_sequence_t *self, mutation_id_t id,
        mutation_t *mutation);
int tree_sequence_get_samples(tree_sequence_t *self, node_id_t **samples);
int tree_sequence_get_sample_index_map(tree_sequence_t *self,
        node_id_t **sample_index_map);

int tree_sequence_get_provenance_strings(tree_sequence_t *self,
        size_t *num_provenance_strings, char ***provenance_strings);
int tree_sequence_simplify(tree_sequence_t *self, node_id_t *samples,
        size_t sample_size, int flags, tree_sequence_t *output);
int tree_sequence_get_pairwise_diversity(tree_sequence_t *self,
    node_id_t *samples, size_t num_samples, double *pi);

int tree_diff_iterator_alloc(tree_diff_iterator_t *self,
        tree_sequence_t *tree_sequence);
int tree_diff_iterator_free(tree_diff_iterator_t *self);
int tree_diff_iterator_next(tree_diff_iterator_t *self, double *length,
        node_record_t **nodes_out, node_record_t **nodes_in);
void tree_diff_iterator_print_state(tree_diff_iterator_t *self, FILE *out);

int sparse_tree_alloc(sparse_tree_t *self, tree_sequence_t *tree_sequence,
        int flags);
int sparse_tree_free(sparse_tree_t *self);
int sparse_tree_copy(sparse_tree_t *self, sparse_tree_t *source);
int sparse_tree_equal(sparse_tree_t *self, sparse_tree_t *other);
int sparse_tree_set_tracked_leaves(sparse_tree_t *self,
        size_t num_tracked_leaves, node_id_t *tracked_leaves);
int sparse_tree_set_tracked_leaves_from_leaf_list(sparse_tree_t *self,
        leaf_list_node_t *head, leaf_list_node_t *tail);
int sparse_tree_get_root(sparse_tree_t *self, node_id_t *root);
int sparse_tree_get_parent(sparse_tree_t *self, node_id_t u, node_id_t *parent);
int sparse_tree_get_children(sparse_tree_t *self, node_id_t u,
        size_t *num_children, node_id_t **children);
int sparse_tree_get_time(sparse_tree_t *self, node_id_t u, double *t);
int sparse_tree_get_mrca(sparse_tree_t *self, node_id_t u, node_id_t v,
        node_id_t *mrca);
int sparse_tree_get_num_leaves(sparse_tree_t *self, node_id_t u,
        size_t *num_leaves);
int sparse_tree_get_num_tracked_leaves(sparse_tree_t *self, node_id_t u,
        size_t *num_tracked_leaves);
int sparse_tree_get_leaf_list(sparse_tree_t *self, node_id_t u,
        leaf_list_node_t **head, leaf_list_node_t **tail);
int sparse_tree_get_sites(sparse_tree_t *self, site_t **sites, list_len_t *sites_length);
void sparse_tree_print_state(sparse_tree_t *self, FILE *out);
/* Method for positioning the tree in the sequence. */
int sparse_tree_first(sparse_tree_t *self);
int sparse_tree_last(sparse_tree_t *self);
int sparse_tree_next(sparse_tree_t *self);
int sparse_tree_prev(sparse_tree_t *self);

int newick_converter_alloc(newick_converter_t *self,
        tree_sequence_t *tree_sequence, size_t precision, double Ne);
int newick_converter_next(newick_converter_t *self, double *length,
        char **tree);
int newick_converter_free(newick_converter_t *self);
void newick_converter_print_state(newick_converter_t *self, FILE *out);

int vcf_converter_alloc(vcf_converter_t *self,
        tree_sequence_t *tree_sequence, unsigned int ploidy, const char *chrom);
int vcf_converter_get_header(vcf_converter_t *self, char **header);
int vcf_converter_next(vcf_converter_t *self, char **record);
int vcf_converter_free(vcf_converter_t *self);
void vcf_converter_print_state(vcf_converter_t *self, FILE *out);

int ld_calc_alloc(ld_calc_t *self, tree_sequence_t *tree_sequence);
int ld_calc_free(ld_calc_t *self);
void ld_calc_print_state(ld_calc_t *self, FILE *out);
int ld_calc_get_r2(ld_calc_t *self, size_t a, size_t b, double *r2);
int ld_calc_get_r2_array(ld_calc_t *self, size_t a, int direction,
        size_t max_mutations, double max_distance,
        double *r2, size_t *num_r2_values);

int hapgen_alloc(hapgen_t *self, tree_sequence_t *tree_sequence);
int hapgen_get_haplotype(hapgen_t *self, node_id_t j, char **haplotype);
int hapgen_free(hapgen_t *self);
void hapgen_print_state(hapgen_t *self, FILE *out);

int vargen_alloc(vargen_t *self, tree_sequence_t *tree_sequence, int flags);
int vargen_next(vargen_t *self, site_t **site, char *genotypes);
int vargen_free(vargen_t *self);
void vargen_print_state(vargen_t *self, FILE *out);

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

void recomb_map_print_state(recomb_map_t *self, FILE *out);

int mutgen_alloc(mutgen_t *self, double mutation_rate, gsl_rng *rng,
        int alphabet, size_t mutation_block_size);
int mutgen_free(mutgen_t *self);
/* TODO finalise this interface */
int mutgen_generate_tables_tmp(mutgen_t *self, node_table_t *nodes,
        edgeset_table_t *edgesets);
int mutgen_populate_tables(mutgen_t *self, site_table_t *sites,
        mutation_table_t *mutations);
void mutgen_print_state(mutgen_t *self, FILE *out);

int node_table_alloc(node_table_t *self, size_t max_rows_increment,
        size_t max_total_name_length_increment);
int node_table_add_row(node_table_t *self, uint32_t flags, double time,
        population_id_t population, const char *name);
int node_table_set_columns(node_table_t *self, size_t num_rows, uint32_t *flags, double *time,
        population_id_t *population, char *name, list_len_t *name_length);
int node_table_reset(node_table_t *self);
int node_table_free(node_table_t *self);
void node_table_print_state(node_table_t *self, FILE *out);

int edgeset_table_alloc(edgeset_table_t *self, size_t max_rows_increment,
        size_t max_total_children_length_increment);
int edgeset_table_add_row(edgeset_table_t *self, double left, double right,
        node_id_t parent, node_id_t *children, list_len_t children_length);
int edgeset_table_set_columns(edgeset_table_t *self, size_t num_rows, double *left,
        double *right, node_id_t *parent, node_id_t *children,
        list_len_t *children_length);
int edgeset_table_reset(edgeset_table_t *self);
int edgeset_table_free(edgeset_table_t *self);
void edgeset_table_print_state(edgeset_table_t *self, FILE *out);

int site_table_alloc(site_table_t *self, size_t max_rows_increment,
        size_t max_total_ancestral_state_length_increment);
int site_table_add_row(site_table_t *self, double position, const char *ancestral_state,
        list_len_t ancestral_state_length);
int site_table_set_columns(site_table_t *self, size_t num_rows,
        double *position, const char *ancestral_state, list_len_t *ancestral_state_length);
bool site_table_equal(site_table_t *self, site_table_t *other);
int site_table_reset(site_table_t *self);
int site_table_free(site_table_t *self);
void site_table_print_state(site_table_t *self, FILE *out);

void mutation_table_print_state(mutation_table_t *self, FILE *out);
int mutation_table_alloc(mutation_table_t *self, size_t max_rows_increment,
        size_t max_total_derived_state_length_increment);
int mutation_table_add_row(mutation_table_t *self, site_id_t site, node_id_t node,
        const char *derived_state, list_len_t derived_state_length);
int mutation_table_set_columns(mutation_table_t *self, size_t num_rows,
        site_id_t *site, node_id_t *node, const char *derived_state,
        list_len_t *derived_state_length);
bool mutation_table_equal(mutation_table_t *self, mutation_table_t *other);
int mutation_table_reset(mutation_table_t *self);
int mutation_table_free(mutation_table_t *self);
void mutation_table_print_state(mutation_table_t *self, FILE *out);

int migration_table_alloc(migration_table_t *self, size_t max_rows_increment);
int migration_table_add_row(migration_table_t *self, double left, double right,
        node_id_t node, population_id_t source, population_id_t dest, double time);
int migration_table_set_columns(migration_table_t *self, size_t num_rows,
        double *left, double *right, node_id_t *node, population_id_t *source,
        population_id_t *dest, double *time);
int migration_table_reset(migration_table_t *self);
int migration_table_free(migration_table_t *self);
void migration_table_print_state(migration_table_t *self, FILE *out);

const char * msp_strerror(int err);
void __msp_safe_free(void **ptr);

double compute_falling_factorial_log(unsigned int  m);

#define msp_safe_free(pointer) __msp_safe_free((void **) &(pointer))


#endif /*__MSPRIME_H__*/
