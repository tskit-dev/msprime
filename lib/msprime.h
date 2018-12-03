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


#include <stdio.h>
#include <stdbool.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>

#include "avl.h"
#include "fenwick.h"
#include "object_heap.h"
#include "tables.h"
#include "trees.h"

#define MSP_MODEL_HUDSON 0
#define MSP_MODEL_SMC 1
#define MSP_MODEL_SMC_PRIME 2
#define MSP_MODEL_BETA 3
#define MSP_MODEL_DIRAC 4
#define MSP_MODEL_DTWF 5

/* Alphabets for mutation generator */
#define MSP_ALPHABET_BINARY     0
#define MSP_ALPHABET_NUCLEOTIDE 1

/* Flags for mutgen */
#define MSP_KEEP_SITES  1


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
    uint32_t left; /* TODO CHANGE THIS - not a good name! */
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
    /* private state */
    double m;
    double phi;
    double integration_epsabs;
    double integration_epsrel;
    size_t integration_workspace_size;
    gsl_integration_workspace *integration_workspace;
} beta_coalescent_t;

typedef struct {
    double psi;
    double c; // constant
} dirac_coalescent_t;

typedef struct _simulation_model_t {
    int type;
    double population_size;
    union {
        beta_coalescent_t beta_coalescent;
        dirac_coalescent_t dirac_coalescent;
    } params;
    /* If the model allocates memory this function should be non-null. */
    void (*free)(struct _simulation_model_t *model);
    /* Time and rate conversions */
    double (*model_time_to_generations)(struct _simulation_model_t *model, double t);
    double (*generations_to_model_time)(struct _simulation_model_t *model, double g);
    double (*generation_rate_to_model_rate)(struct _simulation_model_t *model, double rg);
    double (*model_rate_to_generation_rate)(struct _simulation_model_t *model, double rm);
} simulation_model_t;

/* Recombination map */

typedef struct {
    uint32_t num_loci;      /* size of the genetic coordinate space  */
    double sequence_length; /* size of the physical coordinate space */
    double total_recombination_rate;
    size_t size;            /* the total number of values in the map */
    double *positions;
    double *rates;
} recomb_map_t;

typedef struct _msp_t {
    gsl_rng *rng;
    /* input parameters */
    simulation_model_t model;
    bool store_migrations;
    uint32_t num_samples;
    uint32_t num_loci;
    double recombination_rate;
    recomb_map_t *recomb_map;
    uint32_t num_populations;
    sample_t *samples;
    double start_time;
    tree_sequence_t *from_ts;
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
    double time;
    double *migration_matrix;
    population_t *populations;
    avl_tree_t breakpoints;
    avl_tree_t overlap_counts;
    fenwick_t links;
    /* memory management */
    object_heap_t avl_node_heap;
    object_heap_t segment_heap;
    object_heap_t node_mapping_heap;
    /* The tables used to store the simulation state */
    table_collection_t tables;
    table_collection_position_t from_position;
    /* edges are buffered in a flat array until they are squashed and flushed */
    edge_t *buffered_edges;
    size_t num_buffered_edges;
    size_t max_buffered_edges;
    /* Methods for getting the waiting time until the next common ancestor
     * event and the event are defined by the simulation model */
    double (*get_common_ancestor_waiting_time)(struct _msp_t *self, population_id_t pop);
    int (*common_ancestor_event)(struct _msp_t *selt, population_id_t pop);
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

typedef struct {
    double position;
    node_id_t node;
    const char *ancestral_state;
    const char *derived_state;
} infinite_sites_mutation_t;

typedef struct {
    int alphabet;
    gsl_rng *rng;
    double start_time;
    double end_time;
    double mutation_rate;
    avl_tree_t sites;
    block_allocator_t allocator;
} mutgen_t;

int msp_alloc(msp_t *self,
        size_t num_samples, sample_t *samples,
        recomb_map_t *recomb_map, tree_sequence_t *from, gsl_rng *rng);
int msp_set_start_time(msp_t *self, double start_time);
int msp_set_simulation_model_hudson(msp_t *self, double population_size);
int msp_set_simulation_model_smc(msp_t *self, double population_size);
int msp_set_simulation_model_smc_prime(msp_t *self, double population_size);
int msp_set_simulation_model_dtwf(msp_t *self, double population_size);
int msp_set_simulation_model_dirac(msp_t *self, double population_size, double psi,
    double c);
int msp_set_simulation_model_beta(msp_t *self, double population_size, double alpha,
        double truncation_point);
int msp_set_store_migrations(msp_t *self, bool store_migrations);
int msp_set_num_populations(msp_t *self, size_t num_populations);
int msp_set_node_mapping_block_size(msp_t *self, size_t block_size);
int msp_set_segment_block_size(msp_t *self, size_t block_size);
int msp_set_avl_node_block_size(msp_t *self, size_t block_size);
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
int msp_populate_tables(msp_t *self, table_collection_t *tables);
int msp_reset(msp_t *self);
int msp_print_state(msp_t *self, FILE *out);
int msp_free(msp_t *self);
void msp_verify(msp_t *self);

int msp_get_ancestors(msp_t *self, segment_t **ancestors);
int msp_get_breakpoints(msp_t *self, size_t *breakpoints);
int msp_get_migration_matrix(msp_t *self, double *migration_matrix);
int msp_get_num_migration_events(msp_t *self, size_t *num_migration_events);
int msp_get_samples(msp_t *self, sample_t **samples);
int msp_get_population_configuration(msp_t *self, size_t population_id,
        double *initial_size, double *growth_rate);
int msp_compute_population_size(msp_t *self, size_t population_id,
        double time, double *pop_size);
int msp_is_completed(msp_t *self);

simulation_model_t * msp_get_model(msp_t *self);
const char * msp_get_model_name(msp_t *self);
bool msp_get_store_migrations(msp_t *self);
double msp_get_recombination_rate(msp_t *self);
double msp_get_time(msp_t *self);
size_t msp_get_num_samples(msp_t *self);
size_t msp_get_num_loci(msp_t *self);
size_t msp_get_num_populations(msp_t *self);
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


int recomb_map_alloc_uniform(recomb_map_t *self, uint32_t num_loci,
        double sequence_length, double rate);
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
int recomb_map_phys_to_discrete_genetic(recomb_map_t *self, double phys_x,
        uint32_t *locus);
size_t recomb_map_get_size(recomb_map_t *self);
int recomb_map_get_positions(recomb_map_t *self, double *positions);
int recomb_map_get_rates(recomb_map_t *self, double *rates);

void recomb_map_print_state(recomb_map_t *self, FILE *out);

int mutgen_alloc(mutgen_t *self, double mutation_rate, gsl_rng *rng,
        int alphabet, size_t mutation_block_size);
int mutgen_set_time_interval(mutgen_t *self, double start_time, double end_time);
int mutgen_free(mutgen_t *self);
int mutgen_generate(mutgen_t *self, table_collection_t *tables, int flags);
void mutgen_print_state(mutgen_t *self, FILE *out);

/* Functions exposed here for unit testing. Not part of public API. */
double compute_falling_factorial_log(unsigned int  m);
double compute_dirac_coalescence_rate(unsigned int num_ancestors, double psi, double c);
int msp_compute_beta_integral(msp_t *self, unsigned int num_ancestors, double alpha, double *result);
int msp_beta_compute_coalescence_rate(msp_t *self, unsigned int num_ancestors, double *result);

int msp_multi_merger_common_ancestor_event(msp_t *self, double x, avl_tree_t *ancestors, avl_tree_t *Q);

#endif /*__MSPRIME_H__*/
