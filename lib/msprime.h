/*
** Copyright (C) 2015-2020 University of Oxford
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
#include <stdbool.h>

#include <gsl/gsl_rng.h>
#include <tskit.h>

#include "util.h"
#include "avl.h"
#include "fenwick.h"
#include "object_heap.h"

#define MSP_MODEL_HUDSON 0
#define MSP_MODEL_SMC 1
#define MSP_MODEL_SMC_PRIME 2
#define MSP_MODEL_BETA 3
#define MSP_MODEL_DIRAC 4
#define MSP_MODEL_DTWF 5
#define MSP_MODEL_SWEEP 6
#define MSP_MODEL_WF_PED 7

/* Exit codes from msp_run to distinguish different reasons for exiting
 * before coalescence. */
#define MSP_EXIT_COALESCENCE 0
#define MSP_EXIT_MAX_EVENTS 1
#define MSP_EXIT_MAX_TIME 2

#define MSP_NODE_IS_RE_EVENT (1u << 17)
#define MSP_NODE_IS_CA_EVENT (1u << 18)
#define MSP_NODE_IS_MIG_EVENT (1u << 19)
#define MSP_NODE_IS_CEN_EVENT (1u << 20)

/* Flags for verify */
#define MSP_VERIFY_BREAKPOINTS (1 << 1)

/* Flags for mutgen */
#define MSP_KEEP_SITES 1
#define MSP_DISCRETE_SITES 2

/* Pedigree states */
#define MSP_PED_STATE_UNCLIMBED 0
#define MSP_PED_STATE_CLIMBING 1
#define MSP_PED_STATE_CLIMB_COMPLETE 2

/* FIXME: Using these typedefs to keep the diff size small on the initial
 * tskit transition. Can remove later. */
typedef tsk_id_t population_id_t;
typedef tsk_id_t node_id_t;
typedef tsk_id_t mutation_id_t;
typedef tsk_id_t site_id_t;

/* TODO int16 is surely fine, but it won't make the segment struct any smaller
 * because of alignment requirements. After tskit has been separarated,
 * we will probably want to make both population_id_t and label_id_t
 * int16_t, and then we'll reduce the size of the segment struct a
 * fair bit. It's not worth doing until then though.
 */
typedef tsk_id_t label_id_t;

typedef struct segment_t_t {
    /* TODO change to population */
    population_id_t population_id;
    label_id_t label;
    /* During simulation we use genetic coordinates */
    double left;
    double right;
    double left_mass;
    double right_mass;
    node_id_t value;
    size_t id;
    struct segment_t_t *prev;
    struct segment_t_t *next;
} segment_t;

typedef struct {
    double left; /* TODO CHANGE THIS - not a good name! */
    uint32_t value;
} node_mapping_t;

typedef struct {
    population_id_t population_id;
    double time;
} sample_t;

typedef struct {
    double initial_size;
    double growth_rate;
    double start_time;
    avl_tree_t *ancestors;
    tsk_size_t num_potential_destinations;
    tsk_id_t *potential_destinations;
} population_t;

typedef struct individual_t_t {
    tsk_id_t id;
    tsk_id_t tsk_id;
    struct individual_t_t **parents;
    avl_tree_t *segments;
    int sex;
    double time;
    bool queued;
    // For debugging, to ensure we only merge once.
    bool merged;
} individual_t;

typedef struct {
    individual_t *inds;
    size_t num_inds;
    size_t ploidy;
    individual_t **samples;
    size_t num_samples;
    avl_tree_t ind_heap;
    int state;
} pedigree_t;

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
    double c;
} dirac_coalescent_t;

/* Forward declaration */
struct _msp_t;

typedef struct {
    /* TODO document these parameters.*/
    double start_frequency;
    double end_frequency;
    double alpha;
    double dt;
} genic_selection_trajectory_t;

typedef struct _sweep_t {
    /* TODO change the name of this to position */
    double locus;
    union {
        /* Future trajectory simulation models would go here */
        genic_selection_trajectory_t genic_selection_trajectory;
    } trajectory_params;
    int (*generate_trajectory)(struct _sweep_t *self, struct _msp_t *simulator,
        size_t *num_steps, double **time, double **allele_frequency);
    void (*print_state)(struct _sweep_t *self, FILE *out);
} sweep_t;

typedef struct _simulation_model_t {
    int type;
    union {
        beta_coalescent_t beta_coalescent;
        dirac_coalescent_t dirac_coalescent;
        sweep_t sweep;
    } params;
    /* If the model allocates memory this function should be non-null. */
    void (*free)(struct _simulation_model_t *model);
} simulation_model_t;

typedef struct {
    double *position;
    double *value;
    size_t size;
} interval_map_t;

/* Recombination map */
typedef struct {
    interval_map_t map;
    double total_recombination_rate;
    double *cumulative;
    bool discrete;
} recomb_map_t;

typedef struct _msp_t {
    gsl_rng *rng;
    /* input parameters */
    simulation_model_t model;
    bool store_migrations;
    bool store_full_arg;
    uint32_t num_samples;
    double sequence_length;
    recomb_map_t recomb_map;
    double gene_conversion_rate;
    double gene_conversion_track_length;
    uint32_t num_populations;
    uint32_t num_labels;
    uint32_t ploidy;
    sample_t *samples;
    double start_time;
    tsk_treeseq_t *from_ts;
    simulation_model_t initial_model;
    double *initial_migration_matrix;
    population_t *initial_populations;
    pedigree_t *pedigree;
    /* allocation block sizes */
    size_t avl_node_block_size;
    size_t node_mapping_block_size;
    size_t segment_block_size;
    /* Counters for statistics */
    size_t num_re_events;
    size_t num_ca_events;
    size_t num_gc_events;
    size_t num_rejected_ca_events;
    size_t *num_migration_events;
    size_t num_trapped_re_events;
    size_t num_multiple_re_events;
    size_t num_noneffective_gc_events;
    size_t num_fenwick_rebuilds;
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
    double time;
    double *migration_matrix;
    population_t *populations;
    avl_tree_t non_empty_populations;
    avl_tree_t breakpoints;
    avl_tree_t overlap_counts;
    /* We keep an independent Fenwick tree for each label */
    fenwick_t *links;
    /* memory management */
    object_heap_t avl_node_heap;
    object_heap_t node_mapping_heap;
    /* We keep an independent segment heap for each label */
    object_heap_t *segment_heap;
    /* The tables used to store the simulation state */
    tsk_table_collection_t *tables;
    tsk_bookmark_t from_position;
    /* edges are buffered in a flat array until they are squashed and flushed */
    tsk_edge_t *buffered_edges;
    size_t num_buffered_edges;
    size_t max_buffered_edges;
    /* Methods for getting the waiting time until the next common ancestor
     * event and the event are defined by the simulation model */
    double (*get_common_ancestor_waiting_time)(
        struct _msp_t *self, population_id_t pop, label_id_t label);
    int (*common_ancestor_event)(
        struct _msp_t *selt, population_id_t pop, label_id_t label);
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

/* The site_t and mutation_t are similar the equivalent tsk_ types,
 * with some extra fields that we need for mutation generation. */
typedef struct mutation_t {
    tsk_id_t id;
    tsk_id_t site;
    tsk_id_t node;
    char *derived_state;
    tsk_size_t derived_state_length;
    char *metadata;
    tsk_size_t metadata_length;
    double time;
    struct mutation_t *parent;
    struct mutation_t *next;
    bool new;
    bool keep;
} mutation_t;

typedef struct {
    double position;
    char *ancestral_state;
    tsk_size_t ancestral_state_length;
    char *metadata;
    tsk_size_t metadata_length;
    mutation_t *mutations;
    size_t mutations_length;
    bool new;
} site_t;

typedef struct {
    size_t num_alleles;
    char **alleles;
    tsk_size_t *allele_length;
    double *root_distribution;
    double *transition_matrix;
} mutation_matrix_t;

typedef struct {
    int32_t mutation_type_id; // following SLiM's MutationMetadataRec
    int64_t next_mutation_id; // following SLiM's slim_mutationid_t
    tsk_blkalloc_t allocator;
} slim_mutator_t;

typedef struct {
    uint64_t start_allele;
    uint64_t next_allele;
    tsk_blkalloc_t allocator;
} infinite_alleles_t;

typedef struct _mutation_model_t {
    union {
        mutation_matrix_t mutation_matrix;
        slim_mutator_t slim_mutator;
        infinite_alleles_t infinite_alleles;
        /* Other known mutation models */
    } params;
    void (*print_state)(struct _mutation_model_t *model, FILE *out);
    int (*free)(struct _mutation_model_t *model);
    int (*choose_root_state)(
        struct _mutation_model_t *model, gsl_rng *rng, site_t *site);
    int (*transition)(struct _mutation_model_t *model, gsl_rng *rng,
        const char *parent_allele, tsk_size_t parent_allele_length,
        const char *parent_metadata, tsk_size_t parent_metadata_length,
        mutation_t *mutation);
} mutation_model_t;

typedef struct {
    gsl_rng *rng;
    interval_map_t *rate_map;
    double start_time;
    double end_time;
    size_t block_size;
    avl_tree_t sites;
    tsk_blkalloc_t allocator;
    mutation_model_t *model;
} mutgen_t;

int msp_alloc(msp_t *self, size_t num_samples, sample_t *samples,
    recomb_map_t *recomb_map, tsk_table_collection_t *from_ts_tables, gsl_rng *rng);
int msp_set_start_time(msp_t *self, double start_time);
int msp_set_simulation_model_hudson(msp_t *self);
int msp_set_simulation_model_smc(msp_t *self);
int msp_set_simulation_model_smc_prime(msp_t *self);
int msp_set_simulation_model_dtwf(msp_t *self);
int msp_set_simulation_model_wf_ped(msp_t *self);
int msp_set_simulation_model_dirac(msp_t *self, double psi, double c);
int msp_set_simulation_model_beta(msp_t *self, double alpha, double truncation_point);
int msp_set_simulation_model_sweep_genic_selection(msp_t *self, double position,
    double start_frequency, double end_frequency, double alpha, double dt);

int msp_set_store_migrations(msp_t *self, bool store_migrations);
int msp_set_store_full_arg(msp_t *self, bool store_full_arg);
int msp_set_ploidy(msp_t *self, int ploidy);
int msp_set_num_populations(msp_t *self, size_t num_populations);
int msp_set_dimensions(msp_t *self, size_t num_populations, size_t num_labels);
int msp_set_gene_conversion_rate(msp_t *self, double rate, double track_length);
int msp_set_node_mapping_block_size(msp_t *self, size_t block_size);
int msp_set_segment_block_size(msp_t *self, size_t block_size);
int msp_set_avl_node_block_size(msp_t *self, size_t block_size);
int msp_set_migration_matrix(msp_t *self, size_t size, double *migration_matrix);
int msp_set_population_configuration(
    msp_t *self, int population_id, double initial_size, double growth_rate);

int msp_add_population_parameters_change(
    msp_t *self, double time, int population_id, double size, double growth_rate);
int msp_add_migration_rate_change(
    msp_t *self, double time, int source_pop, int dest_pop, double migration_rate);
int msp_add_mass_migration(
    msp_t *self, double time, int source, int dest, double proportion);
int msp_add_simple_bottleneck(
    msp_t *self, double time, int population_id, double intensity);
int msp_add_instantaneous_bottleneck(
    msp_t *self, double time, int population_id, double strength);
int msp_add_census_event(msp_t *self, double time);

int alloc_individual(individual_t *ind, size_t ploidy);
int msp_alloc_pedigree(msp_t *self, size_t num_inds, size_t ploidy);
int msp_free_pedigree(msp_t *self);
int msp_set_pedigree(msp_t *self, size_t num_rows, int *inds, int *parents,
    double *times, uint32_t *is_sample);
int msp_pedigree_load_pop(msp_t *self);
void msp_check_samples(msp_t *self);
int msp_pedigree_build_ind_queue(msp_t *self);
int msp_pedigree_push_ind(msp_t *self, individual_t *ind);
int msp_pedigree_pop_ind(msp_t *self, individual_t **ind);
int msp_pedigree_add_individual_segment(
    msp_t *self, individual_t *ind, segment_t *segment, size_t parent_ix);
int msp_pedigree_climb(msp_t *self);
void msp_print_pedigree_inds(msp_t *self, FILE *out);

int msp_initialise(msp_t *self);
int msp_run(msp_t *self, double max_time, unsigned long max_events);
int msp_debug_demography(msp_t *self, double *end_time);
int msp_finalise_tables(msp_t *self);
int msp_reset(msp_t *self);
int msp_print_state(msp_t *self, FILE *out);
int msp_free(msp_t *self);
void msp_verify(msp_t *self, int options);

int msp_get_ancestors(msp_t *self, segment_t **ancestors);
int msp_get_breakpoints(msp_t *self, size_t *breakpoints);
int msp_get_migration_matrix(msp_t *self, double *migration_matrix);
int msp_get_num_migration_events(msp_t *self, size_t *num_migration_events);
int msp_get_samples(msp_t *self, sample_t **samples);
int msp_get_population_configuration(
    msp_t *self, size_t population_id, double *initial_size, double *growth_rate);
int msp_compute_population_size(
    msp_t *self, size_t population_id, double time, double *pop_size);
int msp_is_completed(msp_t *self);

simulation_model_t *msp_get_model(msp_t *self);
const char *msp_get_model_name(msp_t *self);
bool msp_get_store_migrations(msp_t *self);
double msp_get_recombination_rate(msp_t *self);
double msp_get_gene_conversion_rate(msp_t *self);
double msp_get_time(msp_t *self);
size_t msp_get_num_samples(msp_t *self);
size_t msp_get_num_loci(msp_t *self);
size_t msp_get_num_populations(msp_t *self);
size_t msp_get_num_labels(msp_t *self);
size_t msp_get_num_population_ancestors(msp_t *self, tsk_id_t population);
size_t msp_get_num_ancestors(msp_t *self);
size_t msp_get_num_breakpoints(msp_t *self);
size_t msp_get_num_nodes(msp_t *self);
size_t msp_get_num_edges(msp_t *self);
size_t msp_get_num_migrations(msp_t *self);
size_t msp_get_num_avl_node_blocks(msp_t *self);
size_t msp_get_num_node_mapping_blocks(msp_t *self);
size_t msp_get_num_segment_blocks(msp_t *self);
size_t msp_get_num_common_ancestor_events(msp_t *self);
size_t msp_get_num_rejected_common_ancestor_events(msp_t *self);
size_t msp_get_num_recombination_events(msp_t *self);
size_t msp_get_num_gene_conversion_events(msp_t *self);

int interval_map_alloc(
    interval_map_t *self, size_t size, double *position, double *value);
int interval_map_alloc_single(
    interval_map_t *self, double sequence_length, double value);
int interval_map_free(interval_map_t *self);
void interval_map_print_state(interval_map_t *self, FILE *out);
double interval_map_get_sequence_length(interval_map_t *self);
size_t interval_map_get_size(interval_map_t *self);
size_t interval_map_get_num_intervals(interval_map_t *self);
size_t interval_map_get_index(interval_map_t *self, double x);

typedef double (*msp_convert_func)(void *obj, double rate);

int recomb_map_alloc_uniform(
    recomb_map_t *self, double sequence_length, double rate, bool discrete);
int recomb_map_alloc(
    recomb_map_t *self, size_t size, double *position, double *rate, bool discrete);
int recomb_map_copy(recomb_map_t *to, recomb_map_t *from);
int recomb_map_free(recomb_map_t *self);
double recomb_map_get_sequence_length(recomb_map_t *self);
bool recomb_map_get_discrete(recomb_map_t *self);
double recomb_map_get_total_recombination_rate(recomb_map_t *self);
void recomb_map_convert_rates(recomb_map_t *self, msp_convert_func convert, void *obj);
size_t recomb_map_get_size(recomb_map_t *self);
int recomb_map_get_positions(recomb_map_t *self, double *positions);
int recomb_map_get_rates(recomb_map_t *self, double *rates);

void recomb_map_print_state(recomb_map_t *self, FILE *out);

double recomb_map_mass_between_left_exclusive(
    recomb_map_t *self, double left, double right);
double recomb_map_mass_between(recomb_map_t *self, double left, double right);
double recomb_map_mass_to_position(recomb_map_t *self, double mass);
double recomb_map_position_to_mass(recomb_map_t *self, double position);
double recomb_map_shift_by_mass(recomb_map_t *self, double pos, double mass);
double recomb_map_sample_poisson(recomb_map_t *self, gsl_rng *rng, double start);

int matrix_mutation_model_factory(mutation_model_t *self, int model);
int matrix_mutation_model_alloc(mutation_model_t *self, size_t num_alleles,
    char **alleles, size_t *allele_length, double *root_distribution,
    double *transition_matrix);
int slim_mutation_model_alloc(mutation_model_t *self, int32_t mutation_type_id,
    int64_t next_mutation_id, size_t block_size);
int infinite_alleles_mutation_model_alloc(
    mutation_model_t *self, uint64_t start_allele, tsk_flags_t options);
int mutation_model_free(mutation_model_t *self);
void mutation_model_print_state(mutation_model_t *self, FILE *out);

int mutgen_alloc(mutgen_t *self, gsl_rng *rng, interval_map_t *rate_map,
    mutation_model_t *model, size_t mutation_block_size);
int mutgen_set_time_interval(mutgen_t *self, double start_time, double end_time);
int mutgen_free(mutgen_t *self);
int mutgen_generate(mutgen_t *self, tsk_table_collection_t *tables, int flags);
void mutgen_print_state(mutgen_t *self, FILE *out);

/* Functions exposed here for unit testing. Not part of public API. */
int msp_multi_merger_common_ancestor_event(
    msp_t *self, avl_tree_t *ancestors, avl_tree_t *Q, uint32_t k);

#endif /*__MSPRIME_H__*/
