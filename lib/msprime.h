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
#include "rate_map.h"

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
#define MSP_KEEP_SITES (1 << 0)
#define MSP_DISCRETE_SITES (1 << 1)
#define MSP_KEPT_MUTATIONS_BEFORE_END_TIME (1 << 2)

/* Pedigree states */
#define MSP_PED_STATE_UNCLIMBED 0
#define MSP_PED_STATE_CLIMBING 1
#define MSP_PED_STATE_CLIMB_COMPLETE 2

typedef tsk_id_t population_id_t;
typedef tsk_id_t label_id_t;

typedef struct segment_t_t {
    population_id_t population;
    label_id_t label;
    double left;
    double right;
    tsk_id_t value;
    size_t id;
    struct segment_t_t *prev;
    struct segment_t_t *next;
} segment_t;

typedef struct {
    double position;
    tsk_size_t value;
} node_mapping_t;

typedef struct {
    double initial_size;
    double growth_rate;
    double start_time;
    avl_tree_t *ancestors;
    tsk_size_t num_potential_destinations;
    tsk_id_t *potential_destinations;
} population_t;

/* Note: we might want to make a distinction here between "individual"
 * and "pedigree individual". We might want to have a more generic
 * individual type at some point which just consists of an array of
 * ploidy segment pointers (one head for each strand). This struct
 * is quite specialised for the purpose of running the pedigree
 * simulation. */
typedef struct individual_t_t {
    tsk_id_t id;
    /* Note: probably simpler to make parents a list of tsk_id_ts,
     * save a bit of pointer fiddling */
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
    /* JK We don't need the samples any more, as we can derive this
     * from the tables. */
    individual_t **samples;
    size_t num_samples;
    avl_tree_t ind_heap;
    int state;
} pedigree_t;

typedef struct {
    double time;
    tsk_id_t sample;
    population_id_t population;
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
    double position;
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
    double left;
    tsk_size_t count;
} overlap_count_t;

typedef struct _msp_t {
    gsl_rng *rng;
    /* input parameters */
    simulation_model_t model;
    bool store_migrations;
    bool store_full_arg;
    double sequence_length;
    bool discrete_genome;
    rate_map_t recomb_map;
    rate_map_t gc_map;
    double gc_track_length;
    uint32_t num_populations;
    uint32_t num_labels;
    uint32_t ploidy;
    double start_time;
    pedigree_t *pedigree;
    /* Initial state for replication */
    segment_t **root_segments;
    overlap_count_t *initial_overlaps;
    simulation_model_t initial_model;
    double *initial_migration_matrix;
    population_t *initial_populations;
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
    fenwick_t *recomb_mass_index;
    fenwick_t *gc_mass_index;
    /* memory management */
    object_heap_t avl_node_heap;
    object_heap_t node_mapping_heap;
    /* We keep an independent segment heap for each label */
    object_heap_t *segment_heap;
    /* The tables used to store the simulation state */
    tsk_table_collection_t *tables;
    tsk_bookmark_t input_position;
    /* edges are buffered in a flat array until they are squashed and flushed */
    tsk_edge_t *buffered_edges;
    tsk_size_t num_buffered_edges;
    tsk_size_t max_buffered_edges;
    /* Methods for getting the waiting time until the next common ancestor
     * event and the event are defined by the simulation model */
    double (*get_common_ancestor_waiting_time)(
        struct _msp_t *self, population_id_t pop, label_id_t label);
    int (*common_ancestor_event)(
        struct _msp_t *selt, population_id_t pop, label_id_t label);
} msp_t;

/* Demographic events */
typedef struct {
    population_id_t population;
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
    population_id_t population;
    double proportion;
} simple_bottleneck_t;

typedef struct {
    population_id_t population;
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
    tsk_table_collection_t *tables;
    double start_time;
    double end_time;
    size_t block_size;
    rate_map_t rate_map;
    avl_tree_t sites;
    tsk_blkalloc_t allocator;
    mutation_model_t *model;
} mutgen_t;

int msp_alloc(msp_t *self, tsk_table_collection_t *tables, gsl_rng *rng);
int msp_set_simulation_model_hudson(msp_t *self);
int msp_set_simulation_model_smc(msp_t *self);
int msp_set_simulation_model_smc_prime(msp_t *self);
int msp_set_simulation_model_dtwf(msp_t *self);
int msp_set_simulation_model_wf_ped(msp_t *self);
int msp_set_simulation_model_dirac(msp_t *self, double psi, double c);
int msp_set_simulation_model_beta(msp_t *self, double alpha, double truncation_point);
int msp_set_simulation_model_sweep_genic_selection(msp_t *self, double position,
    double start_frequency, double end_frequency, double alpha, double dt);

int msp_set_start_time(msp_t *self, double start_time);
int msp_set_store_migrations(msp_t *self, bool store_migrations);
int msp_set_store_full_arg(msp_t *self, bool store_full_arg);
int msp_set_ploidy(msp_t *self, int ploidy);
int msp_set_recombination_map(msp_t *self, size_t size, double *position, double *rate);
int msp_set_recombination_rate(msp_t *self, double rate);
int msp_set_gene_conversion_map(
    msp_t *self, size_t size, double *position, double *rate);
int msp_set_gene_conversion_rate(msp_t *self, double rate);
int msp_set_gene_conversion_track_length(msp_t *self, double track_length);
int msp_set_discrete_genome(msp_t *self, bool is_discrete);
int msp_set_num_labels(msp_t *self, size_t num_labels);
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
int msp_get_population_configuration(
    msp_t *self, size_t population_id, double *initial_size, double *growth_rate);
int msp_compute_population_size(
    msp_t *self, size_t population_id, double time, double *pop_size);
int msp_is_completed(msp_t *self);

simulation_model_t *msp_get_model(msp_t *self);
const char *msp_get_model_name(msp_t *self);
bool msp_get_store_migrations(msp_t *self);
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

int mutgen_alloc(mutgen_t *self, gsl_rng *rng, tsk_table_collection_t *tables,
    mutation_model_t *model, size_t mutation_block_size);
int mutgen_set_time_interval(mutgen_t *self, double start_time, double end_time);
int mutgen_set_rate(mutgen_t *self, double rate);
int mutgen_set_rate_map(mutgen_t *self, size_t size, double *position, double *rate);
int mutgen_free(mutgen_t *self);
int mutgen_generate(mutgen_t *self, int flags);
void mutgen_print_state(mutgen_t *self, FILE *out);

/* Functions exposed here for unit testing. Not part of public API. */
int msp_multi_merger_common_ancestor_event(
    msp_t *self, avl_tree_t *ancestors, avl_tree_t *Q, uint32_t k, uint32_t num_pots);

#endif /*__MSPRIME_H__*/
