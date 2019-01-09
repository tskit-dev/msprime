/*
** Copyright (C) 2016-2018 University of Oxford
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

/*
 * Unit tests for the low-level msprime API.
 */

#include "msprime.h"

#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <unistd.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <CUnit/Basic.h>


/* Global variables used for test in state in the test suite */

char * _tmp_file_name;
FILE * _devnull;

static void
verify_simulator_tsk_treeseq_equality(msp_t *msp, tsk_treeseq_t *tree_seq,
        mutgen_t *mutgen, double scale)
{
    int ret;
    uint32_t num_samples = msp_get_num_samples(msp);
    uint32_t j;
    tsk_node_t node;
    sample_t *samples;
    node_id_t *sample_ids;

    CU_ASSERT_EQUAL_FATAL(
            tsk_treeseq_get_num_samples(tree_seq),
            msp_get_num_samples(msp));
    CU_ASSERT_EQUAL_FATAL(
            tsk_treeseq_get_num_edges(tree_seq),
            msp_get_num_edges(msp));
    CU_ASSERT_EQUAL_FATAL(
            tsk_treeseq_get_num_migrations(tree_seq),
            msp_get_num_migrations(msp));
    ret = msp_get_samples(msp, &samples);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_treeseq_get_num_nodes(tree_seq) >= num_samples);
    CU_ASSERT_TRUE(tsk_migration_tbl_equals(
                msp->tables->migrations, tree_seq->tables->migrations))

    for (j = 0; j < num_samples; j++) {
        ret = tsk_treeseq_get_node(tree_seq, j, &node);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(node.population, samples[j].population_id);
        CU_ASSERT_EQUAL(node.time, samples[j].time);
    }
    /* Samples should always be 0..n - 1 here for simulations */
    ret = tsk_treeseq_get_samples(tree_seq, &sample_ids);
    CU_ASSERT_FATAL(sample_ids != NULL);
    for (j = 0; j < num_samples; j++) {
        CU_ASSERT_EQUAL(j, sample_ids[j]);
    }
    mutgen_print_state(mutgen, _devnull);
    tsk_treeseq_print_state(tree_seq, _devnull);
}



/* Simple unit tests for the Fenwick tree API. */
static void
test_fenwick(void)
{
    fenwick_t t;
    int64_t s;
    size_t j, n;
    for (n = 1; n < 100; n++) {
        s = 0;
        CU_ASSERT(fenwick_alloc(&t, n) == 0);
        for (j = 1; j <= n; j++) {
            fenwick_increment(&t, j, (int64_t) j);
            s = s + (int64_t) j;
            CU_ASSERT(fenwick_get_value(&t, j) == (int64_t) j);
            CU_ASSERT(fenwick_get_cumulative_sum(&t, j) == s);
            CU_ASSERT(fenwick_get_total(&t) == s);
            CU_ASSERT(fenwick_find(&t, s) == j);
            fenwick_set_value(&t, j, 0);
            CU_ASSERT(fenwick_get_value(&t, j) == 0);
            CU_ASSERT(fenwick_get_cumulative_sum(&t, j) == s - (int64_t) j);
            fenwick_set_value(&t, j, (int64_t) j);
            CU_ASSERT(fenwick_get_value(&t, j) == (int64_t) j);
            /* Just make sure that we're seeing the same values even when
             * we expand.
             */
            CU_ASSERT(fenwick_expand(&t, 1) == 0);
        }
        CU_ASSERT(fenwick_free(&t) == 0);
    }
}

static void
test_fenwick_expand(void)
{
    fenwick_t t1, t2;
    int64_t s;
    size_t j, n;

    for (n = 1; n < 100; n++) {
        s = n;
        CU_ASSERT(fenwick_alloc(&t1, n) == 0);
        CU_ASSERT(fenwick_alloc(&t2, 3 * n) == 0);
        for (j = 1; j <= n; j++) {
            fenwick_increment(&t1, j, s);
            fenwick_increment(&t2, j, s);
            CU_ASSERT(fenwick_get_value(&t1, j) == s);
            CU_ASSERT(fenwick_get_value(&t2, j) == s);
        }
        /* After we expand, the internal tree values should be identical */
        CU_ASSERT(t1.size != t2.size);
        CU_ASSERT(fenwick_expand(&t1, 2 * n) == 0);
        CU_ASSERT_EQUAL(t1.size, t2.size);
        CU_ASSERT_EQUAL(memcmp(t1.tree, t2.tree, (t2.size + 1) * sizeof(int64_t)), 0);
        CU_ASSERT_EQUAL(memcmp(t1.values, t2.values, (t2.size + 1) * sizeof(int64_t)), 0);
        CU_ASSERT(fenwick_free(&t1) == 0);
        CU_ASSERT(fenwick_free(&t2) == 0);
    }
}

static void
test_single_locus_two_populations(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = {{0, 0.0}, {0, 0.0}, {1, 40.0}};
    tsk_tbl_collection_t tables;
    tsk_edge_tbl_t *edges;
    tsk_node_tbl_t *nodes;
    tsk_migration_tbl_t *migrations;
    size_t num_edges, num_migrations;
    uint32_t n = 3;
    double t0 = 30.0;
    double t1 = 30.5;
    double t2 = 40.5;
    recomb_map_t recomb_map;

    ret = recomb_map_alloc_uniform(&recomb_map, 1, 1.0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_num_populations(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_store_migrations(&msp, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_mass_migration(&msp, t0, 0, 1, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_mass_migration(&msp, t1, 1, 0, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_mass_migration(&msp, t2, 1, 0, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp);
    msp_print_state(&msp, _devnull);

    nodes = msp.tables->nodes;
    CU_ASSERT_EQUAL_FATAL(nodes->num_rows, 5);
    CU_ASSERT_EQUAL(nodes->time[0], 0);
    CU_ASSERT_EQUAL(nodes->population[0], 0);
    CU_ASSERT_EQUAL(nodes->time[1], 0);
    CU_ASSERT_EQUAL(nodes->population[1], 0);
    CU_ASSERT_EQUAL(nodes->time[2], 40.0);
    CU_ASSERT_EQUAL(nodes->population[2], 1);
    CU_ASSERT_TRUE(nodes->time[3] < 40);
    CU_ASSERT_TRUE(nodes->population[3] == 0);
    CU_ASSERT_TRUE(nodes->time[4] > 40.5);
    CU_ASSERT_TRUE(nodes->population[4] == 0);

    num_edges = msp_get_num_edges(&msp);
    CU_ASSERT_EQUAL_FATAL(num_edges, 4);
    edges = msp.tables->edges;
    CU_ASSERT_EQUAL(edges->parent[0], 3);
    CU_ASSERT_EQUAL(edges->parent[1], 3);
    CU_ASSERT_EQUAL(edges->parent[2], 4);
    CU_ASSERT_EQUAL(edges->parent[3], 4);

    num_migrations = msp_get_num_migrations(&msp);
    CU_ASSERT_EQUAL_FATAL(num_migrations, 3);

    migrations = msp.tables->migrations;
    CU_ASSERT_EQUAL(migrations->time[0], t0)
    CU_ASSERT_EQUAL(migrations->source[0], 0);
    CU_ASSERT_EQUAL(migrations->dest[0], 1);
    CU_ASSERT_EQUAL(migrations->time[1], t1);
    CU_ASSERT_EQUAL(migrations->source[1], 1);
    CU_ASSERT_EQUAL(migrations->dest[1], 0);
    CU_ASSERT_EQUAL(migrations->time[2], t2);
    CU_ASSERT_EQUAL(migrations->source[2], 1);
    CU_ASSERT_EQUAL(migrations->dest[2], 0);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

static void
test_single_locus_many_populations(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    uint32_t num_populations = 500;
    sample_t samples[] = {{0, 0.0}, {num_populations - 1, 0.0}};
    uint32_t n = 2;
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;

    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1, 1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_populations(&msp, num_populations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_add_mass_migration(&msp, 30.0, 0, num_populations - 1, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp);
    msp_print_state(&msp, _devnull);
    CU_ASSERT_EQUAL_FATAL(msp_get_num_edges(&msp), 2);
    CU_ASSERT_EQUAL_FATAL(msp_get_num_nodes(&msp), 3);
    CU_ASSERT_EQUAL(msp.tables->edges->parent[0], 2);
    CU_ASSERT_EQUAL(msp.tables->edges->parent[1], 2);

    CU_ASSERT_TRUE(msp.tables->nodes->time[2] > 30.0);
    CU_ASSERT_EQUAL(msp.tables->nodes->population[2], num_populations - 1);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

static void
test_single_locus_historical_sample(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = {{0, 0.0}, {0, 10.0}};
    tsk_edge_tbl_t *edges;
    tsk_node_tbl_t *nodes;
    uint32_t n = 2;
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;

    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1, 1.0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp);
    msp_print_state(&msp, _devnull);

    CU_ASSERT_EQUAL_FATAL(msp_get_num_nodes(&msp), 3);
    nodes = msp.tables->nodes;
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(nodes->time[0], 0);
    CU_ASSERT_EQUAL(nodes->time[1], 10);
    CU_ASSERT_TRUE(nodes->time[2] > 10);

    CU_ASSERT_EQUAL_FATAL(msp_get_num_edges(&msp), 2);
    edges = msp.tables->edges;
    CU_ASSERT_EQUAL(edges->left[0], 0);
    CU_ASSERT_EQUAL(edges->right[0], 1);
    CU_ASSERT_EQUAL(edges->parent[0], 2);
    CU_ASSERT_EQUAL(edges->left[1], 0);
    CU_ASSERT_EQUAL(edges->right[1], 1);
    CU_ASSERT_EQUAL(edges->parent[1], 2);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_tbl_collection_free(&tables);
    recomb_map_free(&recomb_map);
}

static void
test_single_locus_historical_sample_start_time(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = {{0, 0.0}, {0, 10.0}};
    tsk_edge_tbl_t *edges;
    tsk_node_tbl_t *nodes;
    uint32_t n = 2;
    size_t j;
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;
    double start_times[] = {0, 2, 10, 10.0001, 1000};
    /* double start_times[] = {0, 2, 9.99}; //10, 1000}; */
    /* double start_times[] = {10.00}; //10, 1000}; */

    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1, 1.0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < sizeof(start_times) / sizeof(double); j++) {

        ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_start_time(&msp, start_times[j]);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_initialise(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(msp.time, start_times[j]);

        msp_print_state(&msp, _devnull);
        ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        /* msp_print_state(&msp, stdout); */
        msp_verify(&msp);
        /* msp_print_state(&msp, _devnull); */

        CU_ASSERT_EQUAL_FATAL(msp_get_num_nodes(&msp), 3);
        nodes = msp.tables->nodes;
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(nodes->time[0], 0);
        CU_ASSERT_EQUAL(nodes->time[1], 10);
        CU_ASSERT_TRUE(nodes->time[2] > 10);
        CU_ASSERT_TRUE(nodes->time[2] > start_times[j]);

        CU_ASSERT_EQUAL_FATAL(msp_get_num_edges(&msp), 2);
        edges = msp.tables->edges;
        CU_ASSERT_EQUAL(edges->left[0], 0);
        CU_ASSERT_EQUAL(edges->right[0], 1);
        CU_ASSERT_EQUAL(edges->parent[0], 2);
        CU_ASSERT_EQUAL(edges->left[1], 0);
        CU_ASSERT_EQUAL(edges->right[1], 1);
        CU_ASSERT_EQUAL(edges->parent[1], 2);

        ret = msp_free(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tsk_tbl_collection_clear(&tables);
        CU_ASSERT_EQUAL(ret, 0);
    }
    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}


static void
test_simulator_getters_setters(void)
{
    int ret;
    uint32_t j;
    uint32_t n = 10;
    uint32_t m = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    double migration_matrix[] = {0, 0, 0, 0};
    double matrix[4], growth_rate, initial_size;
    double Ne = 4;
    size_t migration_events[4];
    size_t breakpoints[m];
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    ret = recomb_map_alloc_uniform(&recomb_map, m, 1.0, m - 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < n; j++) {
        samples[j].time = j;
        samples[j].population_id = j % 2;
    }
    CU_ASSERT_EQUAL(msp_alloc(&msp, 0, NULL, NULL, NULL, rng), MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);
    CU_ASSERT_EQUAL(msp_alloc(&msp, n, samples, NULL, NULL, rng), MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);
    ret = msp_alloc(&msp, 0, samples, &recomb_map, &tables, rng);
    msp_free(&msp);
    samples[0].time = 1.0;
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_SAMPLES);
    msp_free(&msp);
    samples[0].time = -1.0;
    CU_ASSERT_EQUAL(msp_alloc(&msp, n, samples, &recomb_map, &tables, rng),
            MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);
    samples[0].time = 0.0;

    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(msp_set_node_mapping_block_size(&msp, 0),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_segment_block_size(&msp, 0),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_avl_node_block_size(&msp, 0),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_num_populations(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
            msp_set_population_configuration(&msp, -1, 0, 0),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL(
            msp_set_population_configuration(&msp, 3, 0, 0),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);

    ret = msp_set_simulation_model_hudson(&msp, Ne);
    CU_ASSERT_EQUAL(msp_get_model(&msp)->type, MSP_MODEL_HUDSON);
    CU_ASSERT_EQUAL(msp_get_model(&msp)->population_size, Ne);

    ret = msp_set_num_populations(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(
            msp_get_population_configuration(&msp, 3, NULL, NULL),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    ret = msp_set_population_configuration(&msp, 0, 2 * Ne, 0.5);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(
            msp_set_migration_matrix(&msp, 0, NULL),
            MSP_ERR_BAD_MIGRATION_MATRIX);
    CU_ASSERT_EQUAL(
            msp_set_migration_matrix(&msp, 3, migration_matrix),
            MSP_ERR_BAD_MIGRATION_MATRIX);
    migration_matrix[0] = 1;
    CU_ASSERT_EQUAL(
            msp_set_migration_matrix(&msp, 4, migration_matrix),
            MSP_ERR_BAD_MIGRATION_MATRIX);
    migration_matrix[0] = 0;
    migration_matrix[1] = -1;
    CU_ASSERT_EQUAL(
            msp_set_migration_matrix(&msp, 4, migration_matrix),
            MSP_ERR_BAD_MIGRATION_MATRIX);

    migration_matrix[1] = 1;
    migration_matrix[2] = 1;
    ret = msp_set_migration_matrix(&msp, 4, migration_matrix);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_store_migrations(&msp, true);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_get_population_configuration(&msp, 0, &initial_size, &growth_rate);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(initial_size, 2 * Ne);
    CU_ASSERT_EQUAL(growth_rate, 0.5);
    CU_ASSERT_EQUAL(msp_get_recombination_rate(&msp), 1.0);

    CU_ASSERT_TRUE(msp_get_store_migrations(&msp));
    CU_ASSERT_EQUAL(msp_get_num_avl_node_blocks(&msp), 1);
    CU_ASSERT_EQUAL(msp_get_num_node_mapping_blocks(&msp), 1);
    CU_ASSERT_EQUAL(msp_get_num_segment_blocks(&msp), 1);
    CU_ASSERT_EQUAL(msp_get_num_populations(&msp), 2);

    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(msp_get_num_breakpoints(&msp), m - 1);
    ret = msp_get_breakpoints(&msp, breakpoints);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < m - 1; j++) {
        CU_ASSERT_EQUAL(breakpoints[j], j + 1);
    }
    ret = msp_get_num_migration_events(&msp, migration_events);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(migration_events[0], 0);
    CU_ASSERT(migration_events[1] > 0);
    CU_ASSERT(migration_events[2] > 0);
    CU_ASSERT_EQUAL(migration_events[3], 0);
    ret = msp_get_migration_matrix(&msp, matrix);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < 4; j++) {
        CU_ASSERT_EQUAL(matrix[j], migration_matrix[j]);
    }
    CU_ASSERT(msp_get_num_common_ancestor_events(&msp) > 0);
    CU_ASSERT(msp_get_num_recombination_events(&msp) > 0);

    free(samples);
    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

static void
test_simulator_model_errors(void)
{
    uint32_t n = 10;
    uint32_t j;
    int ret;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    msp_t msp;
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1, 1.0, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < n; j++) {
        samples[j].time = 0;
        samples[j].population_id = 0;
    }

    CU_ASSERT_EQUAL(msp_alloc(&msp, n, samples, &recomb_map, &tables, rng), 0);
    CU_ASSERT_EQUAL(msp_get_model(&msp)->type, MSP_MODEL_HUDSON);
    CU_ASSERT_EQUAL(msp_add_simple_bottleneck(&msp, 1, 0, 1), 0);
    CU_ASSERT_EQUAL(msp_add_instantaneous_bottleneck(&msp, 1, 0, 1), 0);
    CU_ASSERT_EQUAL(msp_initialise(&msp), 0);
    CU_ASSERT_EQUAL(msp_run(&msp, DBL_MAX, ULONG_MAX), 0);
    CU_ASSERT_EQUAL(msp_free(&msp), 0);

    for (j = 0; j < 2; j++) {
        tsk_tbl_collection_clear(&tables);
        CU_ASSERT_EQUAL(msp_alloc(&msp, n, samples, &recomb_map, &tables, rng), 0);
        if (j == 0) {
            CU_ASSERT_EQUAL(msp_set_simulation_model_smc(&msp, 0.25), 0);
        } else {
            CU_ASSERT_EQUAL(msp_set_simulation_model_smc_prime(&msp, 0.25), 0);
        }
        CU_ASSERT_EQUAL(msp_add_simple_bottleneck(&msp, 1, 0, 1), MSP_ERR_BAD_MODEL);
        CU_ASSERT_EQUAL(msp_add_instantaneous_bottleneck(&msp, 1, 0, 1),
                MSP_ERR_BAD_MODEL);
        CU_ASSERT_EQUAL(msp_free(&msp), 0);
    }

    free(samples);
    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

static void
test_demographic_events(void)
{
    int ret;
    uint32_t j, k;
    uint32_t n = 10;
    uint32_t m = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    double migration_matrix[] = {0, 1, 1, 0};
    double last_time, time, pop_size;
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    ret = recomb_map_alloc_uniform(&recomb_map, m, 1.0, m);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < n; j++) {
        samples[j].time = j;
        samples[j].population_id = j % 2;
    }
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    /* Zero or negative population sizes are not allowed */
    ret = msp_set_population_configuration(&msp, 0, -1, 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = msp_set_population_configuration(&msp, 0, 0, 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    ret = msp_set_num_populations(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 1, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 1, 2, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_migration_matrix(&msp, 4, migration_matrix);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(
        msp_add_mass_migration(&msp, 10, -1, 0, 1),
        MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL(
        msp_add_mass_migration(&msp, 10, 2, 0, 1),
        MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL(
        msp_add_mass_migration(&msp, 10, 0, 0, 1),
        MSP_ERR_SOURCE_DEST_EQUAL);
    CU_ASSERT_EQUAL(
        msp_add_mass_migration(&msp, 10, 0, 1, -5),
        MSP_ERR_BAD_PARAM_VALUE);

    CU_ASSERT_EQUAL(
        msp_add_migration_rate_change(&msp, 10, -2, 2.0),
        MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
    CU_ASSERT_EQUAL(
        msp_add_migration_rate_change(&msp, 10, -1, -2.0),
        MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
        msp_add_migration_rate_change(&msp, 10, 3, 2.0),
        MSP_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX);

    CU_ASSERT_EQUAL(
        msp_add_population_parameters_change(&msp, 10, -2, 0, 0),
        MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL(
        msp_add_population_parameters_change(&msp, 10, -1, -1, 0),
        MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
        msp_add_population_parameters_change(&msp, 10, -1, GSL_NAN, GSL_NAN),
        MSP_ERR_BAD_PARAM_VALUE);

    CU_ASSERT_EQUAL(
        msp_add_simple_bottleneck(&msp, 10, -1, 0),
        MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL(
        msp_add_simple_bottleneck(&msp, 10, 0, -1),
        MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
        msp_add_simple_bottleneck(&msp, 10, 0, 1.1),
        MSP_ERR_BAD_PARAM_VALUE);

    CU_ASSERT_EQUAL(
        msp_add_instantaneous_bottleneck(&msp, 10, 2, 0),
        MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL(
        msp_add_simple_bottleneck(&msp, 10, 0, -1),
        MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL_FATAL(
        msp_add_simple_bottleneck(&msp, 10, -1, 0),
        MSP_ERR_POPULATION_OUT_OF_BOUNDS);

    ret = msp_add_mass_migration(&msp, 0.1, 0, 1, 0.5);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_migration_rate_change(&msp, 0.2, 1, 2.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_migration_rate_change(&msp, 0.3, -1, 3.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.4, 0, 0.5, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.5, -1, 0.5, 2.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.6, 0, GSL_NAN, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.7, 1, 1, GSL_NAN);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_simple_bottleneck(&msp, 0.8, 0, 0.5);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_instantaneous_bottleneck(&msp, 0.9, 0, 2.0);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(
            msp_add_mass_migration(&msp, 0.1, 0, 1, 0.5),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
    CU_ASSERT_EQUAL(
            msp_add_migration_rate_change(&msp, 0.2, 1, 2.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
    CU_ASSERT_EQUAL(
            msp_add_population_parameters_change(&msp, 0.4, 0, 0.5, 1.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
    CU_ASSERT_EQUAL(
            msp_add_simple_bottleneck(&msp, 0.7, 0, 1.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
    CU_ASSERT_EQUAL(
            msp_add_instantaneous_bottleneck(&msp, 0.8, 0, 1.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);

    CU_ASSERT_EQUAL(
        msp_debug_demography(&msp, &time),
        MSP_ERR_BAD_STATE);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    j = 0;
    last_time = 0;
    do {
        ret = msp_debug_demography(&msp, &time);
        CU_ASSERT_EQUAL(ret, 0);
        msp_print_state(&msp, _devnull);
        ret = msp_compute_population_size(&msp, 10, last_time, &pop_size);
        CU_ASSERT_EQUAL(ret, MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        for (k = 0; k < 2; k++) {
            ret = msp_compute_population_size(&msp, k, last_time, &pop_size);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_TRUE(pop_size >= 0);
            ret = msp_compute_population_size(&msp, k, time, &pop_size);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_TRUE(pop_size >= 0);
            ret = msp_compute_population_size(&msp, k, last_time + (time - last_time) / 2,
                    &pop_size);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_TRUE(pop_size >= 0);
        }
        j++;
        last_time = time;
    } while (! gsl_isinf(time));
    CU_ASSERT_EQUAL(j, 10);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(
        msp_run(&msp, DBL_MAX, ULONG_MAX),
        MSP_ERR_BAD_STATE);
    ret = msp_reset(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);

    free(samples);
    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

static void
test_demographic_events_start_time(void)
{
    int ret;
    uint32_t n = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1, 1.0, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_start_time(msp, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(msp, 0.4, 0, 0.5, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_DEMOGRAPHIC_EVENT_TIME);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

static void
test_time_travel_error(void)
{
    int ret;
    uint32_t n = 100;
    sample_t *samples = calloc(n, sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1, 1.0, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_simple_bottleneck(&msp, 0.1, 0, 0.75);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_simple_bottleneck(&msp, 0.1, 0, 1.0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_TIME_TRAVEL);

    free(samples);
    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

static int
get_num_children(size_t node, tsk_edge_tbl_t *edges)
{
    int num_children = 0;
    size_t i;

    for (i = 0; i < edges->num_rows; i++) {
        if (edges->parent[i] == node) {
            num_children++;
        }
    }
    return num_children;
}

static void
test_dtwf_deterministic(void)
{
    int j, ret;
    uint32_t n = 10;
    uint32_t m = 2;
    unsigned long seed = 133;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_tbl_collection_t tables[2];
    recomb_map_t recomb_map;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, m, 1.0, m);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    for (j = 0; j < 2; j++) {
        ret = tsk_tbl_collection_alloc(&tables[j], 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        gsl_rng_set(rng, seed);
        ret = msp_alloc(msp, n, samples, &recomb_map, &tables[j], rng);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_simulation_model_dtwf(msp, n);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_population_configuration(msp, 0, n, 0);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_initialise(msp);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_run(msp, DBL_MAX, UINT32_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        msp_verify(msp);
        ret = msp_finalise_tables(msp);
        CU_ASSERT_EQUAL(ret, 0);
        msp_free(msp);
        CU_ASSERT_EQUAL(tables[j].migrations->num_rows, 0);
        CU_ASSERT(tables[j].nodes->num_rows > 0);
        CU_ASSERT(tables[j].edges->num_rows > 0);

    }
    CU_ASSERT_TRUE(tsk_node_tbl_equals(tables[0].nodes, tables[1].nodes));
    CU_ASSERT_TRUE(tsk_edge_tbl_equals(tables[0].edges, tables[1].edges));

    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    for (j = 0; j < 2; j++) {
        tsk_tbl_collection_free(&tables[j]);
    }
    recomb_map_free(&recomb_map);
}

static void
test_mixed_model_simulation(void)
{
    int ret;
    uint32_t j, k;
    uint32_t n = 10;
    double N = 100;
    int model;
    double initial_size, growth_rate;
    const char *model_name;
    tsk_tbl_collection_t tables;
    tsk_treeseq_t ts;
    recomb_map_t recomb_map;
    double g = -1.0 / 8192;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    memset(samples, 0, n * sizeof(sample_t));
    ret = recomb_map_alloc_uniform(&recomb_map, 10, 1.0, 10);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_populations(msp, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(msp, N);
    CU_ASSERT_EQUAL(ret, 0);
    /* Set the populations to 1, 2, and 3N. We don't simulate them,
     * but they should be equal at the end */
    ret = msp_set_population_configuration(msp, 0, N, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(msp, 1, 2 * N, g);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(msp, 2, 3 * N, 2 * g);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    /* Run for 10 generations each for the different models, alternating */
    j = 0;
    while ((ret = msp_run(msp, j * 10, UINT32_MAX)) == 2) {
        msp_verify(msp);
        /* Check that our populations and growth rates are still correct */
        for (k = 0; k < 3; k++) {
            ret = msp_get_population_configuration(msp, k, &initial_size, &growth_rate);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(initial_size, (k + 1) * N);
            CU_ASSERT_EQUAL_FATAL(growth_rate, k * g);
        }
        CU_ASSERT_FALSE(msp_is_completed(msp));
        if (j % 2 == 1) {
            model = msp_get_model(msp)->type;
            CU_ASSERT_EQUAL(model, MSP_MODEL_HUDSON);
            model_name = msp_get_model_name(msp);
            CU_ASSERT_STRING_EQUAL(model_name, "hudson");
            ret = msp_set_simulation_model_dtwf(msp, N);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            model = msp_get_model(msp)->type;
            CU_ASSERT_EQUAL(model, MSP_MODEL_DTWF);
            model_name = msp_get_model_name(msp);
            CU_ASSERT_STRING_EQUAL(model_name, "dtwf");
            ret = msp_set_simulation_model_hudson(msp, N);
            CU_ASSERT_EQUAL(ret, 0);
        }
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(msp_is_completed(msp));
    CU_ASSERT_TRUE(j > 5);

    ret = msp_finalise_tables(msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* TODO remove this when populate tables takes table_collection as arg */
    tables.sequence_length = msp->num_loci;
    CU_ASSERT_EQUAL_FATAL(tables.sequence_length, msp->num_loci);
    ret = tsk_treeseq_alloc(&ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_print_state(&ts, _devnull);
    tsk_treeseq_free(&ts);
    tsk_tbl_collection_free(&tables);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

static void
test_dtwf_single_locus_simulation(void)
{
    int ret;
    const char *model_name;
    int i;
    uint32_t n = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    int num_coalescent_events = 0;
    int num_children;
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1, 1.0, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(msp, n);
    CU_ASSERT_EQUAL(ret, 0);
    model_name = msp_get_model_name(msp);
    CU_ASSERT_STRING_EQUAL(model_name, "dtwf");

    ret = msp_run(msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(msp);

    /* For the single locus sim we should have n-1 coalescent events,
     * counting multiple mergers as multiple coalescent events */
    for (i = 0; i < msp->tables->nodes->num_rows; i++) {
        num_children = get_num_children(i, msp->tables->edges);
        if (num_children > 0) {
            num_coalescent_events += num_children - 1;
        }
    }
    CU_ASSERT_EQUAL(num_coalescent_events, n-1);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

static void
test_single_locus_simulation(void)
{
    int ret;
    int model;
    const char *model_name;
    uint32_t j;
    uint32_t n = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1, 1.0, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    /* For the single locus sim we should have exactly n - 1 events */
    for (j = 0; j < n - 2; j++) {
        ret = msp_run(msp, DBL_MAX, 1);
        CU_ASSERT_EQUAL(ret, 1);
        msp_verify(msp);
    }
    ret = msp_run(msp, DBL_MAX, 1);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(msp);

    model = msp_get_model(msp)->type;
    CU_ASSERT_EQUAL(model, MSP_MODEL_HUDSON);
    model_name = msp_get_model_name(msp);
    CU_ASSERT_STRING_EQUAL(model_name, "hudson");

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

static void
test_dtwf_multi_locus_simulation(void)
{
    int ret;
    uint32_t n = 100;
    uint32_t m = 100;
    long seed = 10;
    double migration_matrix[] = {0, 0.1, 0.1, 0};
    const char *model_name;
    size_t num_ca_events, num_re_events;
    double t;
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;

    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, m, 1.0, m);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_set(rng, seed);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(msp, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_populations(msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_migration_matrix(msp, 4, migration_matrix);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_store_migrations(msp, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);
    model_name = msp_get_model_name(msp);
    CU_ASSERT_STRING_EQUAL(model_name, "dtwf");

    ret = msp_run(msp, DBL_MAX, ULONG_MAX);
    msp_verify(msp);
    num_ca_events = msp_get_num_common_ancestor_events(msp);
    num_re_events = msp_get_num_recombination_events(msp);
    CU_ASSERT_TRUE(num_ca_events > 0);
    CU_ASSERT_TRUE(num_re_events > 0);
    CU_ASSERT_EQUAL(ret, 0);
    msp_free(msp);

    /* Realloc the simulator under different memory params to see if
     * we get the same result. */
    gsl_rng_set(rng, seed);
    tsk_tbl_collection_clear(&tables);
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(msp, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_populations(msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_migration_matrix(msp, 4, migration_matrix);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);
    t = 1;
    /* We should be able to step backward here generation-by-generation until
     * coalescence.
     */
    while ((ret = msp_run(msp, t, ULONG_MAX)) > 0) {
        msp_verify(msp);
        CU_ASSERT_EQUAL_FATAL(msp->time, t);
        t++;
    }
    msp_verify(msp);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(num_ca_events == msp_get_num_common_ancestor_events(msp));
    CU_ASSERT_TRUE(num_re_events == msp_get_num_recombination_events(msp));

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

static void
test_multi_locus_simulation(void)
{
    int ret;
    uint32_t num_events;
    uint32_t n = 100;
    uint32_t m = 100;
    long seed = 10;
    int model;
    int models[] = {MSP_MODEL_HUDSON, MSP_MODEL_SMC, MSP_MODEL_SMC_PRIME};
    double migration_matrix[] = {0, 1, 1, 0};
    size_t migration_events[4];
    const char *model_names[] = {"hudson", "smc", "smc_prime"};
    const char *model_name;
    size_t j;
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;

    ret = recomb_map_alloc_uniform(&recomb_map, m, 1.0, m);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < sizeof(models) / sizeof(int); j++) {
        sample_t *samples = malloc(n * sizeof(sample_t));
        msp_t *msp = malloc(sizeof(msp_t));
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

        CU_ASSERT_FATAL(msp != NULL);
        CU_ASSERT_FATAL(samples != NULL);
        CU_ASSERT_FATAL(rng != NULL);
        gsl_rng_set(rng, seed);
        tsk_tbl_collection_clear(&tables);
        memset(samples, 0, n * sizeof(sample_t));
        ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_num_populations(msp, 2);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_migration_matrix(msp, 4, migration_matrix);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_store_migrations(msp, true);
        CU_ASSERT_EQUAL(ret, 0);
        /* set all the block sizes to something small to provoke the memory
         * expansions. */
        ret = msp_set_avl_node_block_size(msp, 1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_node_mapping_block_size(msp, 1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_segment_block_size(msp, 1);
        CU_ASSERT_EQUAL(ret, 0);
        switch (j) {
            case 0:
                ret = msp_set_simulation_model_hudson(msp, 0.25);
                break;
            case 1:
                ret = msp_set_simulation_model_smc(msp, 0.25);
                break;
            case 2:
                ret = msp_set_simulation_model_smc_prime(msp, 0.25);
                break;
        }
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_initialise(msp);
        CU_ASSERT_EQUAL(ret, 0);

        num_events = 0;
        while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
            msp_verify(msp);
            num_events++;
        }
        CU_ASSERT_EQUAL(ret, 0);
        msp_verify(msp);
        ret = msp_get_num_migration_events(msp, migration_events);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT(num_events > n - 1);
        CU_ASSERT_EQUAL(1 + num_events,
                /* diagonals must be zero here */
                migration_events[1] + migration_events[2] +
                msp_get_num_recombination_events(msp) +
                msp_get_num_common_ancestor_events(msp) +
                msp_get_num_rejected_common_ancestor_events(msp));
        if (models[j] == MSP_MODEL_HUDSON) {
            CU_ASSERT_EQUAL(msp_get_num_rejected_common_ancestor_events(msp), 0);
        }

        gsl_rng_set(rng, seed);
        ret = msp_reset(msp);
        num_events = 0;
        while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
            msp_verify(msp);
            num_events++;
        }
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_get_num_migration_events(msp, migration_events);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(1 + num_events,
                migration_events[1] + migration_events[2] +
                msp_get_num_recombination_events(msp) +
                msp_get_num_common_ancestor_events(msp) +
                msp_get_num_rejected_common_ancestor_events(msp));
        if (models[j] == MSP_MODEL_HUDSON) {
            CU_ASSERT_EQUAL(msp_get_num_rejected_common_ancestor_events(msp), 0);
        }

        model = msp_get_model(msp)->type;
        CU_ASSERT_EQUAL(model, models[j]);
        model_name = msp_get_model_name(msp);
        CU_ASSERT_STRING_EQUAL(model_name, model_names[j]);

        ret = msp_free(msp);
        CU_ASSERT_EQUAL(ret, 0);
        gsl_rng_free(rng);
        free(msp);
        free(samples);
    }
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

static void
test_simulation_replicates(void)
{
    int ret;
    uint32_t n = 100;
    uint32_t m = 100;
    double mutation_rate = 2;
    size_t num_replicates = 10;
    long seed = 10;
    double migration_matrix[] = {0, 1, 1, 0};
    size_t j;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    msp_t msp;
    tsk_treeseq_t ts;
    mutgen_t mutgen;
    tsk_tbl_collection_t tables;
    recomb_map_t recomb_map;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    ret = recomb_map_alloc_uniform(&recomb_map, m, 0.5, 1.0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    gsl_rng_set(rng, seed);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_populations(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_migration_matrix(&msp, 4, migration_matrix);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_store_migrations(&msp, true);
    CU_ASSERT_EQUAL(ret, 0);
    /* set all the block sizes to something small to provoke the memory
     * expansions. */
    ret = msp_set_avl_node_block_size(&msp, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_node_mapping_block_size(&msp, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_segment_block_size(&msp, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = mutgen_alloc(&mutgen, mutation_rate, rng, 0, 3);
    CU_ASSERT_EQUAL(ret, 0);

    for (j = 0; j < num_replicates; j++) {
        ret = msp_run(&msp, DBL_MAX, SIZE_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        msp_verify(&msp);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = msp_finalise_tables(&msp);
        ret = mutgen_generate(&mutgen, &tables, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        tables.sequence_length = m;
        ret = tsk_treeseq_alloc(&ts, &tables, TSK_BUILD_INDEXES);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        verify_simulator_tsk_treeseq_equality(&msp, &ts, &mutgen, 1.0);
        tsk_treeseq_print_state(&ts, _devnull);
        ret = msp_reset(&msp);
        CU_ASSERT_EQUAL_FATAL(msp_get_num_edges(&msp), 0);
        CU_ASSERT_EQUAL_FATAL(msp_get_num_migrations(&msp), 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        tsk_treeseq_free(&ts);
    }
    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(samples);
    tsk_tbl_collection_free(&tables);
    recomb_map_free(&recomb_map);
}

static void
test_bottleneck_simulation(void)
{
    int ret;
    uint32_t n = 100;
    uint32_t m = 100;
    long seed = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    double t1 = 0.1;
    double t2 = 0.5;
    int t1_found = 0;
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    ret = recomb_map_alloc_uniform(&recomb_map, m, 1.0, m);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    gsl_rng_set(rng, seed);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    /* set all the block sizes to something small to provoke the memory
     * expansions. */
    ret = msp_set_avl_node_block_size(msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_node_mapping_block_size(msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_segment_block_size(msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    /* Add a bottleneck that does nothing at time t1 */
    ret = msp_add_simple_bottleneck(msp, t1, 0, 0);
    CU_ASSERT_EQUAL(ret, 0);
    /* Add a bottleneck that coalesces everything at t2 */
    ret = msp_add_simple_bottleneck(msp, t2, 0, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(msp, t1, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 2);
    CU_ASSERT_FALSE(msp_is_completed(msp));
    msp_print_state(msp, _devnull);

    ret = msp_run(msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(msp_is_completed(msp));
    CU_ASSERT_EQUAL(msp->time, t2);
    msp_verify(msp);

    msp_reset(msp);
    while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
        msp_verify(msp);
        if (msp->time == 0.1) {
            t1_found = 1;
        }
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(t1_found);
    CU_ASSERT_EQUAL(msp->time, t2);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

#define SIMPLE_BOTTLENECK 0
#define INSTANTANEOUS_BOTTLENECK 1

typedef struct {
    int type;
    double time;
    uint32_t population_id;
    double parameter;
} bottleneck_desc_t;

static void
test_large_bottleneck_simulation(void)
{
    int ret;
    uint32_t j;
    uint32_t n = 10000;
    uint32_t m = 100;
    long seed = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    uint32_t num_bottlenecks = 10;
    bottleneck_desc_t bottlenecks[num_bottlenecks];
    double t;
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    ret = recomb_map_alloc_uniform(&recomb_map, m, 1.0, m);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    t = 0.1;
    for (j = 0; j < num_bottlenecks; j++) {
        bottlenecks[j].type = SIMPLE_BOTTLENECK;
        bottlenecks[j].time = t;
        bottlenecks[j].parameter = 0.1;
        t += 0.01;
    }
    /* Set the last bottleneck to be full intensity */
    bottlenecks[num_bottlenecks - 1].parameter = 1.0;

    gsl_rng_set(rng, seed);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_bottlenecks; j++) {
        ret = msp_add_simple_bottleneck(msp, bottlenecks[j].time, 0,
                bottlenecks[j].parameter);
        CU_ASSERT_EQUAL(ret, 0);
    }
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    for (j = 0; j < num_bottlenecks - 1; j++) {
        ret = msp_run(msp, bottlenecks[j].time, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 2);
        CU_ASSERT_FALSE(msp_is_completed(msp));
        CU_ASSERT_EQUAL(msp->time, bottlenecks[j].time);
        msp_verify(msp);
    }
    ret = msp_run(msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(msp_is_completed(msp));
    CU_ASSERT_EQUAL(msp->time, bottlenecks[num_bottlenecks - 1].time);
    msp_verify(msp);

    /* Test out resets on partially completed simulations. */
    ret = msp_reset(msp);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_bottlenecks - 1; j++) {
        ret = msp_run(msp, bottlenecks[j].time, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 2);
    }
    ret = msp_reset(msp);
    msp_verify(msp);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}

static void
test_compute_falling_factorial(void)
{
    CU_ASSERT_DOUBLE_EQUAL(compute_falling_factorial_log(0), 0, 0.000000);
    CU_ASSERT_DOUBLE_EQUAL(compute_falling_factorial_log(1), 1.386294, 0.000001);
    CU_ASSERT_DOUBLE_EQUAL(compute_falling_factorial_log(2), 2.484907, 0.000001);
    CU_ASSERT_DOUBLE_EQUAL(compute_falling_factorial_log(3), 3.178054, 0.000001);
    CU_ASSERT_DOUBLE_EQUAL(compute_falling_factorial_log(4), 3.178054, 0.000001);
}


static void
test_compute_dirac_coalescence_rate(void)
{
    // compute_dirac_coalescence_rate(unsigned int num_ancestors, double psi, double c)

    // Falls to Kingman coalescent
    CU_ASSERT_DOUBLE_EQUAL(compute_dirac_coalescence_rate(2, 0.1, 0), 1.0, 0.000000);
    CU_ASSERT_DOUBLE_EQUAL(compute_dirac_coalescence_rate(3, 0.1, 0), 3.0, 0.000000);
    CU_ASSERT_DOUBLE_EQUAL(compute_dirac_coalescence_rate(4, 0.1, 0), 6.0, 0.000000);

    // Pairwise coalescent, λ_2 = 1 + cψ^2 /4
    CU_ASSERT_DOUBLE_EQUAL(compute_dirac_coalescence_rate(2, 0.1, 0.1), 1.00025, 0.000001);
    CU_ASSERT_DOUBLE_EQUAL(compute_dirac_coalescence_rate(2, 0.1, 1.0), 1.0025, 0.000001);
    CU_ASSERT_DOUBLE_EQUAL(compute_dirac_coalescence_rate(2, 0.1, 10.0), 1.025, 0.000001);

    // Other cases, check against lambdab r code
    // PASS when e = 1e-5, FAIL when e = 1e-6, in particular with low psi
    CU_ASSERT_DOUBLE_EQUAL(compute_dirac_coalescence_rate(10, 0.001, 0.5), 45.00001, 0.00001);
    CU_ASSERT_DOUBLE_EQUAL(compute_dirac_coalescence_rate(10, 0.001, 50), 45.00056, 0.00001);
    CU_ASSERT_DOUBLE_EQUAL(compute_dirac_coalescence_rate(10, 0.01, 0.5), 45.00055, 0.00001);
    CU_ASSERT_DOUBLE_EQUAL(compute_dirac_coalescence_rate(100, 0.99, 50.0), 5000.0, 0.00001);
    CU_ASSERT_DOUBLE_EQUAL(compute_dirac_coalescence_rate(101, 0.51, 5.0), 5055.0, 0.00001);
}

/* Because the beta coalescent uses allocated memory we must allocate a simulator.
 * This function hides this detail away for convenience. */
static double
compute_beta_coalescence_rate(unsigned int num_ancestors, double alpha)
{
    int ret;
    msp_t msp;
    unsigned int n = 10;
    double value;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;

    ret = recomb_map_alloc_uniform(&recomb_map, 1, 1.0, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* TODO is truncation_point arbitrary here? */
    ret = msp_set_simulation_model_beta(&msp, 1, alpha, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_beta_compute_coalescence_rate(&msp, num_ancestors, &value);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    msp_free(&msp);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
    free(samples);
    gsl_rng_free(rng);
    return value;
}

static void
test_compute_beta_coalescence_rate(void)
{
    CU_ASSERT_DOUBLE_EQUAL(compute_beta_coalescence_rate(100, 1.01), 225.6396, 0.001);
    CU_ASSERT_DOUBLE_EQUAL(compute_beta_coalescence_rate(100, 1.5), 1140.782, 1140.782*1e-3);
    CU_ASSERT_DOUBLE_EQUAL(compute_beta_coalescence_rate(100, 1.8), 2815.267, 0.1);

    // Pairwise coalescent, alpha is irrelevant
    CU_ASSERT_DOUBLE_EQUAL(compute_beta_coalescence_rate(2, 1.1), 1.0, 0.000000);
    CU_ASSERT_DOUBLE_EQUAL(compute_beta_coalescence_rate(2, 1.5), 1.0, 0.000000);
    CU_ASSERT_DOUBLE_EQUAL(compute_beta_coalescence_rate(2, 1.9), 1.0, 0.000000);
}

static int
compute_beta_coalescence_rate_fails(unsigned int num_ancestors, double alpha)
{
    int ret;
    msp_t msp;
    unsigned int n = 10;
    double value;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;

    ret = recomb_map_alloc_uniform(&recomb_map, 1, 1.0, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_simulation_model_beta(&msp, 1, alpha, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_beta_compute_coalescence_rate(&msp, num_ancestors, &value);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INTEGRATION_FAILED);

    msp_free(&msp);
    recomb_map_free(&recomb_map);
    free(samples);
    gsl_rng_free(rng);
    tsk_tbl_collection_free(&tables);
    return 0;
}

static void
test_gsl_error_handling_beta_coalescent(void)
{
    FILE *old_stderr = stderr;
    gsl_error_handler_t *old_handler;
    old_handler = gsl_set_error_handler_off();
    /* Redirect stderr to avoid spamming output */
    stderr = _devnull;
    compute_beta_coalescence_rate_fails(1000, 1e20);
    gsl_set_error_handler(old_handler);
    stderr = old_stderr;
}

static void
test_multiple_mergers_simulation(void)
{
    int ret;
    size_t j;
    uint32_t n = 100;
    uint32_t m = 100;
    long seed = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, m, 1.0, m);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < 2; j++) {
        gsl_rng_set(rng, seed);
        /* TODO check non-zero sample times here to make sure they fail. */
        memset(samples, 0, n * sizeof(sample_t));
        tsk_tbl_collection_clear(&tables);
        ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
        CU_ASSERT_EQUAL(ret, 0);
        /* TODO what are good parameters here?? */
        if (j == 0) {
            // Use psi = 0.5 for now, but should definitely test for 0 and 1 cases
            ret = msp_set_simulation_model_dirac(msp, 1, 0.5, 1);
        } else {
            ret = msp_set_simulation_model_beta(msp, 1, 1.5, 10.0);
        }
        CU_ASSERT_EQUAL(ret, 0);
        /* TODO check for adding various complications like multiple populations etc
         * to ensure they fail.
         */
        ret = msp_initialise(msp);
        CU_ASSERT_EQUAL(ret, 0);
        msp_print_state(msp, _devnull);

        ret = msp_run(msp, DBL_MAX, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_TRUE(msp_is_completed(msp));
        CU_ASSERT_TRUE(msp->time > 0);
        msp_verify(msp);

        msp_reset(msp);
        while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
            msp_verify(msp);
        }
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_TRUE(msp_is_completed(msp));

        ret = msp_free(msp);
        CU_ASSERT_EQUAL(ret, 0);
    }
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&tables);
}


static void
test_simple_recomb_map(void)
{
    int ret;
    recomb_map_t recomb_map;
    double positions[] = {0.0, 1.0};
    double rates[] = {0.0, 1.0};
    uint32_t num_loci[] = {1, 100, UINT32_MAX};
    size_t j;

    for (j = 0; j < sizeof(num_loci) / sizeof(uint32_t); j++) {
        ret = recomb_map_alloc(&recomb_map, num_loci[j], 1.0,
                positions, rates, 2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        recomb_map_print_state(&recomb_map, _devnull);
        CU_ASSERT_EQUAL(recomb_map_get_num_loci(&recomb_map), num_loci[j]);
        CU_ASSERT_EQUAL(recomb_map_get_size(&recomb_map), 2);
        CU_ASSERT_EQUAL(recomb_map_get_sequence_length(&recomb_map), 1.0);
        CU_ASSERT_EQUAL(
                recomb_map_genetic_to_phys(&recomb_map, num_loci[j]), 1.0);
        CU_ASSERT_EQUAL(
                recomb_map_phys_to_genetic(&recomb_map, 1.0), num_loci[j]);
        ret = recomb_map_free(&recomb_map);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
}

static void
test_recomb_map_errors(void)
{
    int ret;
    recomb_map_t recomb_map;
    uint32_t locus;
    double positions[] = {0.0, 1.0, 2.0};
    double rates[] = {1.0, 2.0, 0.0};

    ret = recomb_map_alloc(&recomb_map, 10, 1.0, positions, rates, 0);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);

    ret = recomb_map_alloc(&recomb_map, 10, 1.0, positions, rates, 1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);

    ret = recomb_map_alloc(&recomb_map, 0, 1.0, positions, rates, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);

    ret = recomb_map_alloc(&recomb_map, 10, 2.0, positions, rates, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);

    positions[0] = 1;
    ret = recomb_map_alloc(&recomb_map, 10, 1.0, positions, rates, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);
    positions[0] = 0;

    positions[1] = 3.0;
    ret = recomb_map_alloc(&recomb_map, 10, 2.0, positions, rates, 3);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);
    positions[1] = 1.0;

    ret = recomb_map_alloc(&recomb_map, 10, 2.0, positions, rates, 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Physical coordinates must be 0 <= x <= L */
    ret = recomb_map_phys_to_discrete_genetic(&recomb_map, -1, &locus);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    ret = recomb_map_phys_to_discrete_genetic(&recomb_map, 2.1, &locus);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    recomb_map_free(&recomb_map);
}

static void
verify_recomb_map(uint32_t num_loci, double length, double *positions,
        double *rates, size_t size)
{

    int ret;
    recomb_map_t recomb_map;
    double total_rate, x, y, z;
    uint32_t locus;
    size_t j;
    size_t num_checks = 1000;
    double eps = 1e-6;
    double *ret_rates, *ret_positions;

    ret_rates = malloc(size * sizeof(double));
    ret_positions = malloc(size * sizeof(double));

    CU_ASSERT_FATAL(ret_rates != NULL);
    CU_ASSERT_FATAL(ret_positions != NULL);

    ret = recomb_map_alloc(&recomb_map, num_loci, length,
           positions, rates, size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    recomb_map_print_state(&recomb_map, _devnull);

    CU_ASSERT_EQUAL(recomb_map_get_num_loci(&recomb_map), num_loci);
    CU_ASSERT_EQUAL(recomb_map_get_size(&recomb_map), size);
    CU_ASSERT_EQUAL(recomb_map_genetic_to_phys(&recomb_map, 0), 0);
    CU_ASSERT_EQUAL(
            recomb_map_genetic_to_phys(&recomb_map, num_loci), length);
    total_rate = recomb_map_get_total_recombination_rate(&recomb_map);
    CU_ASSERT_DOUBLE_EQUAL(
            total_rate / (num_loci - 1),
            recomb_map_get_per_locus_recombination_rate(&recomb_map), eps)

    for (j = 0; j < num_checks; j++) {
        x = j * (num_loci / num_checks);
        y = recomb_map_genetic_to_phys(&recomb_map, x);
        CU_ASSERT_TRUE(0 <= y && y <= length);
        z = recomb_map_phys_to_genetic(&recomb_map, y);
        CU_ASSERT_DOUBLE_EQUAL(x, z, eps);
        ret = recomb_map_phys_to_discrete_genetic(&recomb_map, y, &locus);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(locus, (uint32_t) round(x));
    }
    ret = recomb_map_get_positions(&recomb_map, ret_positions);
    CU_ASSERT_EQUAL(ret, 0);
    ret = recomb_map_get_rates(&recomb_map, ret_rates);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < size; j++) {
        CU_ASSERT_EQUAL(ret_rates[j], rates[j]);
        CU_ASSERT_EQUAL(ret_positions[j], positions[j]);
    }
    ret = recomb_map_free(&recomb_map);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    free(ret_rates);
    free(ret_positions);
}

static void
test_recomb_map_examples(void)
{
    double p1[] =  {0, 0.05019838393314813, 0.36933662489552865, 1};
    double r1[] = {3.5510784169955434, 4.184964179610539, 3.800808140657212, 0};
    double p2[] =  {0, 0.125, 0.875, 1, 4, 8, 16};
    double r2[] = {0.1, 6.0, 3.333, 2.1, 0.0, 2.2, 0};

    verify_recomb_map(2, 1.0, p1, r1, 4);
    verify_recomb_map(1000, 1.0, p1, r1, 4);
    verify_recomb_map(UINT32_MAX, 1.0, p1, r1, 4);

    verify_recomb_map(2, 16.0, p2, r2, 7);
    verify_recomb_map(100, 16.0, p2, r2, 7);
    verify_recomb_map(UINT32_MAX, 16.0, p2, r2, 7);
}

static void
verify_simulate_from(int model, recomb_map_t *recomb_map,
        tsk_tbl_collection_t *from_tables, size_t num_replicates)
{
    int ret;
    size_t j;
    tsk_tbl_collection_position_t pos;
    tsk_tbl_collection_t tables;
    tsk_treeseq_t final;
    tsk_tree_t tree;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_copy(from_tables, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_alloc(&msp, 0, NULL, recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    if (model == MSP_MODEL_DTWF) {
        ret = msp_set_simulation_model_dtwf(&msp, 10);
        CU_ASSERT_EQUAL(ret, 0);
    }
    /* TODO add dirac and other models */
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < num_replicates; j++) {
        msp_verify(&msp);
        ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_TRUE(msp_is_completed(&msp));
        ret = msp_finalise_tables(&msp);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_treeseq_alloc(&final, &tables, TSK_BUILD_INDEXES);
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        ret = tsk_tree_alloc(&tree, &final, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        for (ret = tsk_tree_first(&tree); ret == 1; ret = tsk_tree_next(&tree)) {
            CU_ASSERT_EQUAL_FATAL(tsk_tree_get_num_roots(&tree), 1);
        }
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        ret = tsk_tbl_collection_record_position(from_tables, &pos);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tbl_collection_reset_position(&tables, &pos);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_TRUE(tsk_tbl_collection_equals(from_tables, &tables));

        tsk_treeseq_free(&final);

        tsk_tree_free(&tree);
        ret = msp_reset(&msp);
        /* printf("ret = %s\n", msp_strerror(ret)); */
        CU_ASSERT_EQUAL(ret, 0);
    }
    msp_free(&msp);
    gsl_rng_free(rng);
    tsk_tbl_collection_free(&tables);
}

/* Verify that the initial state we get in a new simulator from calling
 * with from_ts is equivalent to the state in the original simulator */
static void
verify_initial_simulate_from_state(msp_t *msp_source, recomb_map_t *recomb_map,
        tsk_tbl_collection_t *from_tables)
{
    int ret;
    msp_t msp_dest;
    size_t num_ancestors;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    ret = msp_alloc(&msp_dest, 0, NULL, recomb_map, from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp_dest);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    msp_print_state(&msp_dest, _devnull);

    /* In principle we can examine the generated ancestors to ensure
     * that we have correctly translated the segments. In practise this
     * is quite tricky, so we content ourselves with checking the number
     * is the same. */

    num_ancestors = msp_get_num_ancestors(msp_source);
    CU_ASSERT_EQUAL_FATAL(num_ancestors, msp_get_num_ancestors(&msp_dest));

    msp_free(&msp_dest);
    gsl_rng_free(rng);
}

static void
verify_simple_simulate_from(int model, uint32_t n, size_t num_loci, double sequence_length,
        double recombination_rate, size_t num_events, size_t num_replicates)
{
    int ret;
    tsk_tbl_collection_t tables;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t msp;
    recomb_map_t recomb_map;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, num_loci, sequence_length,
            recombination_rate);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    /* Partially run the simulation */
    ret = msp_run(&msp, DBL_MAX, num_events);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_FALSE(msp_is_completed(&msp));
    ret = msp_finalise_tables(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Add some provenances */
    ret = tsk_provenance_tbl_add_row(tables.provenances, "time", 4, "record", 6);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_initial_simulate_from_state(&msp, &recomb_map, &tables);
    verify_simulate_from(model, &recomb_map, &tables, num_replicates);

    msp_free(&msp);
    gsl_rng_free(rng);
    tsk_tbl_collection_free(&tables);
    recomb_map_free(&recomb_map);
    free(samples);
}

static void
test_simulate_from_single_locus(void)
{
    verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 1, 1.0, 0, 5, 1);
    verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 1, 1.0, 0, 5, 1);

    verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 1, 1.0, 1, 5, 1);
    verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 1, 1.0, 1, 5, 1);
}

static void
test_simulate_from_single_locus_sequence_length(void)
{
    double L[] = {0.1, 0.99, 2, 3.33333333, 10, 1e6};
    double rate[] = {1e-7, 0.1, 1, 5};
    size_t j, k;

    for (j = 0; j < sizeof(L) / sizeof(*L); j++) {
        for (k = 0; k < sizeof(rate) / sizeof(*rate); k++) {
            verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 1, L[j], rate[k], 5, 1);
            verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 1, L[j], rate[k], 5, 1);
        }
    }
}

static void
test_simulate_from_multi_locus_sequence_length(void)
{
    double L[] = {0.1, 0.99, 2, 3.33333333, 10, 1e6};
    double rate[] = {1e-8, 0.1, 1, 5};
    size_t j, k;

    for (j = 0; j < sizeof(L) / sizeof(*L); j++) {
        for (k = 0; k < sizeof(rate) / sizeof(*rate); k++) {
            verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 100, L[j], rate[k], 5, 1);
            verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 100, L[j], rate[k], 5, 1);
        }
    }
}

static void
test_simulate_from_single_locus_replicates(void)
{
    verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 1, 1.0, 0, 5, 10);
    verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 1, 1.0, 0, 5, 10);
}

static void
test_simulate_from_multi_locus(void)
{
    verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 100, 1.0, 10.0, 20, 1);
    verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 100, 1.0, 10.0, 20, 1);
}

static void
test_simulate_from_multi_locus_replicates(void)
{
    verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 100, 1.0, 10.0, 20, 10);
    verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 100, 1.0, 10.0, 20, 10);
}

static void
test_simulate_from_empty(void)
{
    verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 1, 1.0, 0, 0, 1);
    verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 1, 1.0, 0, 0, 1);
}

static void
test_simulate_from_completed(void)
{
    int ret;
    uint32_t n = 25;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t msp;
    recomb_map_t recomb_map;
    tsk_tbl_collection_t tables;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    size_t num_loci = 10;
    double recombination_rate = 2;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, num_loci, 1.0, recombination_rate);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables,0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(msp_is_completed(&msp));
    ret = msp_finalise_tables(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    verify_simulate_from(MSP_MODEL_HUDSON, &recomb_map, &tables, 1);

    msp_free(&msp);
    gsl_rng_free(rng);
    tsk_tbl_collection_free(&tables);
    recomb_map_free(&recomb_map);
    free(samples);
}

static void
test_simulate_from_incompatible(void)
{
    int ret;
    msp_t msp;
    recomb_map_t recomb_map;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_tbl_collection_t from_tables;

    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1, 10.0, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&from_tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    from_tables.sequence_length = 1.0;
    /* Add a node so that we get past the first check */
    ret = tsk_population_tbl_add_row(from_tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_tbl_add_row(from_tables.nodes, 0, 0.0,
            0, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Malformed tree sequence */
    from_tables.nodes->individual[0] = 100;
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_FATAL(msp_is_tsk_error(ret));
    CU_ASSERT_EQUAL_FATAL(ret ^ (1 << MSP_TSK_ERR_BIT), TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    CU_ASSERT_STRING_EQUAL(msp_strerror(ret), tsk_strerror(TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS));
    from_tables.nodes->individual[0] = -1;
    msp_free(&msp);

    /* Sequence length mismatch */
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INCOMPATIBLE_FROM_TS);
    msp_free(&msp);

    /* zero samples */
    from_tables.sequence_length = 10.0;
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, MSP_ERR_INSUFFICIENT_SAMPLES);
    msp_free(&msp);

    /* older nodes */
    from_tables.sequence_length = 10.0;
    ret = tsk_node_tbl_add_row(from_tables.nodes, TSK_NODE_IS_SAMPLE, 1.0, 0,
            TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_tbl_add_row(from_tables.nodes, TSK_NODE_IS_SAMPLE, 2.0, 0,
            TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_start_time(&msp, 1.999);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_START_TIME_FROM_TS);
    msp_free(&msp);

    /* Num populations should be 1 */
    tsk_population_tbl_add_row(from_tables.populations, NULL, 0);
    ret = msp_set_start_time(&msp, 0);
    CU_ASSERT_EQUAL(ret, 0);
    from_tables.sequence_length = 10.0;
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, MSP_ERR_INCOMPATIBLE_FROM_TS);
    msp_free(&msp);

    tsk_population_tbl_clear(from_tables.populations);
    tsk_population_tbl_add_row(from_tables.populations, NULL, 0);
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_start_time(&msp, 1.999);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_START_TIME_FROM_TS);
    msp_free(&msp);

    /* Must have legitimate population references */
    from_tables.nodes->population[0] = -1;
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_start_time(&msp, 2.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    msp_free(&msp);

    /* Check to make sure we can run this correctly */
    from_tables.nodes->population[0] = 0;
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_start_time(&msp, 2.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    msp_free(&msp);

    /* Make a tree sequence that we cannot recover trees from. This only happens
     * at initialisation time. */
    ret = tsk_edge_tbl_add_row(from_tables.edges, 0, 1, 1, 0);
    CU_ASSERT_FATAL(ret >= 0);
    tsk_edge_tbl_add_row(from_tables.edges, 0, 1, 2, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_FATAL(msp_is_tsk_error(ret));
    CU_ASSERT_EQUAL_FATAL(ret ^ (1 << MSP_TSK_ERR_BIT),
            TSK_ERR_BAD_EDGES_CONTRADICTORY_CHILDREN);
    CU_ASSERT_STRING_EQUAL(msp_strerror(ret),
            tsk_strerror(TSK_ERR_BAD_EDGES_CONTRADICTORY_CHILDREN));
    msp_free(&msp);

    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    tsk_tbl_collection_free(&from_tables);
}

static void
test_simulate_init_errors(void)
{
    int ret;
    uint32_t n = 25;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t msp;
    recomb_map_t recomb_map;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_tbl_collection_t tables;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1, 1.0, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1.0;

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(&msp, 0, samples, &recomb_map, NULL, rng);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = msp_alloc(&msp, n, samples, &recomb_map, NULL, rng);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = msp_alloc(&msp, n, NULL, &recomb_map, NULL, rng);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = msp_alloc(&msp, n, NULL, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    tsk_tbl_collection_clear(&tables);
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_set_start_time(&msp, -1);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_START_TIME);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    msp_free(&msp);
    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    free(samples);
    tsk_tbl_collection_free(&tables);
}

static void
insert_single_tree(tsk_tbl_collection_t *tables)
{
    /*
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  2   0\n"
        "0  3   0\n";
    const char *edges =
        "0  1   4   0,1\n"
        "0  1   5   2,3\n"
        "0  1   6   4,5\n";
    */
    int ret;
    ret = tsk_node_tbl_add_row(tables->nodes, TSK_NODE_IS_SAMPLE, 0.0, 0,
            TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_tbl_add_row(tables->nodes, TSK_NODE_IS_SAMPLE, 0.0, 0,
            TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_tbl_add_row(tables->nodes, TSK_NODE_IS_SAMPLE, 0.0, 0,
            TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_tbl_add_row(tables->nodes, TSK_NODE_IS_SAMPLE, 0.0, 0,
            TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_tbl_add_row(tables->nodes, 0, 1.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_tbl_add_row(tables->nodes, 0, 2.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_tbl_add_row(tables->nodes, 0, 3.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_tbl_add_row(tables->edges, 0, 1, 4, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_tbl_add_row(tables->edges, 0, 1, 4, 1);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_tbl_add_row(tables->edges, 0, 1, 5, 2);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_tbl_add_row(tables->edges, 0, 1, 5, 3);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_tbl_add_row(tables->edges, 0, 1, 6, 4);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_tbl_add_row(tables->edges, 0, 1, 6, 5);
    CU_ASSERT_FATAL(ret >= 0);

    /* Add a site and a mutation */
    ret = tsk_site_tbl_add_row(tables->sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_tbl_add_row(tables->mutations, 0, 0, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
}

static void
test_single_tree_mutgen(void)
{
    int ret = 0;
    size_t j;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_tbl_collection_t tables1, tables2;

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_tbl_collection_alloc(&tables1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&tables2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables1);
    insert_single_tree(&tables2);

    ret = mutgen_alloc(&mutgen, 0.0, rng, 0, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables1.mutations->num_rows == 0);
    mutgen_print_state(&mutgen, _devnull);
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    gsl_rng_set(rng, 1);
    ret = mutgen_alloc(&mutgen, 10.0, rng, 0, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    mutgen_print_state(&mutgen, _devnull);
    CU_ASSERT_TRUE(tables1.mutations->num_rows > 0);
    CU_ASSERT_TRUE(tables1.mutations->num_rows == tables1.sites->num_rows);
    for (j = 0; j < tables1.mutations->num_rows; j++) {
        CU_ASSERT_TRUE(tables1.mutations->site[j] == j);
        CU_ASSERT_TRUE(tables1.sites->position[j] <= 1.0);
        CU_ASSERT_TRUE(tables1.mutations->node[j] < 6);
        CU_ASSERT_EQUAL(tables1.sites->ancestral_state[j], '0');
        CU_ASSERT_EQUAL(tables1.mutations->derived_state[j], '1');
    }
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Test the reallocing behavior by setting a very small
     * block size.
     */
    gsl_rng_set(rng, 1);
    ret = mutgen_alloc(&mutgen, 10.0, rng, 0, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_tbl_collection_equals(&tables1, &tables2));
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_tbl_collection_free(&tables1);
    tsk_tbl_collection_free(&tables2);
    gsl_rng_free(rng);
}

static void
test_single_tree_mutgen_keep_sites(void)
{
    int ret = 0;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_tbl_collection_t tables;
    tsk_tbl_collection_t copy;
    mutgen_t mutgen;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables);
    ret = tsk_tbl_collection_alloc(&copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&copy);
    CU_ASSERT_TRUE(tsk_tbl_collection_equals(&tables, &copy));

    /* With a mutation rate of 0, we should keep exactly the same set
     * of mutations */
    ret = mutgen_alloc(&mutgen, 0.0, rng, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_tbl_collection_equals(&tables, &copy));
    mutgen_free(&mutgen);

    gsl_rng_set(rng, 2);
    ret = mutgen_alloc(&mutgen, 10.0, rng, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites->num_rows > copy.sites->num_rows);
    CU_ASSERT_TRUE(tables.mutations->num_rows > copy.mutations->num_rows);
    mutgen_free(&mutgen);

    /* If we run precisely the same mutations again we should rejection
     * sample away all of the original positions */
    gsl_rng_set(rng, 2);
    ret = mutgen_alloc(&mutgen, 10.0, rng, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites->num_rows > copy.sites->num_rows);
    CU_ASSERT_TRUE(tables.mutations->num_rows > copy.mutations->num_rows);

    /* add a duplicate site to the original */
    ret = tsk_site_tbl_add_row(tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_TRUE(ret > 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_DUPLICATE_SITE_POSITION);

    mutgen_free(&mutgen);
    tsk_tbl_collection_free(&tables);
    tsk_tbl_collection_free(&copy);
    gsl_rng_free(rng);
}

static void
test_single_tree_mutgen_interval(void)
{
    int ret = 0;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_tbl_collection_t tables1;
    size_t j;
    tsk_id_t node;

    CU_ASSERT_FATAL(rng != NULL);

    ret = tsk_tbl_collection_alloc(&tables1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables1);

    ret = mutgen_alloc(&mutgen, 10.0, rng, 0, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables1.mutations->num_rows > 0);
    mutgen_print_state(&mutgen, _devnull);

    /* End before start is an error */
    ret = mutgen_set_time_interval(&mutgen, 0, -1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    /* Setting start and end == 0 should give 0 mutations */
    ret = mutgen_set_time_interval(&mutgen, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables1.sites->num_rows == 0);
    CU_ASSERT_TRUE(tables1.mutations->num_rows == 0);

    /* Setting start = 3 should give 0 mutations */
    ret = mutgen_set_time_interval(&mutgen, 3, DBL_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables1.sites->num_rows == 0);
    CU_ASSERT_TRUE(tables1.mutations->num_rows == 0);

    /* Setting start = 2 should give mutations only above 4 and 5 */
    ret = mutgen_set_time_interval(&mutgen, 2, DBL_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables1.sites->num_rows > 0);
    CU_ASSERT_TRUE(tables1.mutations->num_rows > 0);
    for (j = 0; j < tables1.sites->num_rows; j++) {
        node = tables1.mutations->node[j];
        CU_ASSERT_TRUE(node == 4 || node == 5);
    }

    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_tbl_collection_free(&tables1);
    gsl_rng_free(rng);
}

static void
test_strerror(void)
{
    int j;
    const char *msg;
    int max_error_code = 1024; /* totally arbitrary */

    for (j = 1; j > -max_error_code; j--) {
        msg = msp_strerror(-j);
        CU_ASSERT_FATAL(msg != NULL);
        CU_ASSERT(strlen(msg) > 0);
    }
}

static void
test_strerror_tskit(void)
{
    int tskit_errors[] = {TSK_ERR_NO_MEMORY,
        TSK_ERR_NODE_OUT_OF_BOUNDS, TSK_ERR_EDGE_OUT_OF_BOUNDS};
    size_t j;
    int err;

    for (j = 0; j < sizeof(tskit_errors) / sizeof(*tskit_errors); j++) {
        err = msp_set_tsk_error(tskit_errors[j]);
        CU_ASSERT_TRUE(msp_is_tsk_error(err));
        CU_ASSERT_STRING_EQUAL(msp_strerror(err), tsk_strerror(tskit_errors[j]));
    }
}

static int
msprime_suite_init(void)
{
    int fd;
    static char template[] = "/tmp/msp_sim_c_test_XXXXXX";

    _tmp_file_name = NULL;
    _devnull = NULL;

    _tmp_file_name = malloc(sizeof(template));
    if (_tmp_file_name == NULL) {
        return CUE_NOMEMORY;
    }
    strcpy(_tmp_file_name, template);
    fd = mkstemp(_tmp_file_name);
    if (fd == -1) {
        return CUE_SINIT_FAILED;
    }
    close(fd);
    _devnull = fopen("/dev/null", "w");
    if (_devnull == NULL) {
        return CUE_SINIT_FAILED;
    }
    return CUE_SUCCESS;
}

static int
msprime_suite_cleanup(void)
{
    if (_tmp_file_name != NULL) {
        unlink(_tmp_file_name);
        free(_tmp_file_name);
    }
    if (_devnull != NULL) {
        fclose(_devnull);
    }
    return CUE_SUCCESS;
}

static void
handle_cunit_error()
{
    fprintf(stderr, "CUnit error occured: %d: %s\n",
            CU_get_error(), CU_get_error_msg());
    exit(EXIT_FAILURE);
}

int
main(int argc, char **argv)
{
    int ret;
    CU_pTest test;
    CU_pSuite suite;
    CU_TestInfo tests[] = {
        {"test_fenwick", test_fenwick},
        {"test_fenwick_expand", test_fenwick_expand},
        {"test_single_locus_two_populations", test_single_locus_two_populations},
        {"test_single_locus_many_populations", test_single_locus_many_populations},
        {"test_single_locus_historical_sample", test_single_locus_historical_sample},
        {"test_single_locus_historical_sample_start_time",
            test_single_locus_historical_sample_start_time},
        {"test_simulator_getters_setters", test_simulator_getters_setters},
        {"test_simulator_model_errors", test_simulator_model_errors},
        {"test_demographic_events", test_demographic_events},
        {"test_demographic_events_start_time", test_demographic_events_start_time},
        {"test_time_travel_error", test_time_travel_error},
        {"test_single_locus_simulation", test_single_locus_simulation},
        {"test_mixed_model_simulation", test_mixed_model_simulation},
        {"test_dtwf_deterministic", test_dtwf_deterministic},
        {"test_dtwf_single_locus_simulation", test_dtwf_single_locus_simulation},
        {"test_multi_locus_simulation", test_multi_locus_simulation},
        {"test_dtwf_multi_locus_simulation", test_dtwf_multi_locus_simulation},
        {"test_simulation_replicates", test_simulation_replicates},
        {"test_bottleneck_simulation", test_bottleneck_simulation},
        {"test_compute_falling_factorial", test_compute_falling_factorial},
        {"test_compute_dirac_coalescence_rate", test_compute_dirac_coalescence_rate},
        {"test_compute_beta_coalescence_rate", test_compute_beta_coalescence_rate},
        {"test_gsl_error_handling_beta_coalescent", test_gsl_error_handling_beta_coalescent},
        {"test_multiple_mergers_simulation", test_multiple_mergers_simulation},
        {"test_large_bottleneck_simulation", test_large_bottleneck_simulation},
        {"test_simple_recombination_map", test_simple_recomb_map},
        {"test_recombination_map_errors", test_recomb_map_errors},
        {"test_recombination_map_examples", test_recomb_map_examples},
        {"test_simulate_from_single_locus", test_simulate_from_single_locus},
        {"test_simulate_from_single_locus_sequence_length",
            test_simulate_from_single_locus_sequence_length},
        {"test_simulate_from_multi_locus_sequence_length",
            test_simulate_from_multi_locus_sequence_length},
        {"test_simulate_from_single_locus_replicates",
            test_simulate_from_single_locus_replicates},
        {"test_simulate_from_multi_locus", test_simulate_from_multi_locus},
        {"test_simulate_from_multi_locus_replicates",
            test_simulate_from_multi_locus_replicates},
        {"test_simulate_from_empty", test_simulate_from_empty},
        {"test_simulate_from_completed", test_simulate_from_completed},
        {"test_simulate_from_incompatible", test_simulate_from_incompatible},
        {"test_simulate_init_errors", test_simulate_init_errors},

        {"test_single_tree_mutgen", test_single_tree_mutgen},
        {"test_single_tree_mutgen_keep_sites", test_single_tree_mutgen_keep_sites},
        {"test_single_tree_mutgen_interval", test_single_tree_mutgen_interval},

        {"test_strerror", test_strerror},
        {"test_strerror_tskit", test_strerror_tskit},
        CU_TEST_INFO_NULL,
    };

    /* We use initialisers here as the struct definitions change between
     * versions of CUnit */
    CU_SuiteInfo suites[] = {
        {
            .pName = "msprime",
            .pInitFunc = msprime_suite_init,
            .pCleanupFunc = msprime_suite_cleanup,
            .pTests = tests
        },
        CU_SUITE_INFO_NULL,
    };
    if (CUE_SUCCESS != CU_initialize_registry()) {
        handle_cunit_error();
    }
    if (CUE_SUCCESS != CU_register_suites(suites)) {
        handle_cunit_error();
    }
    CU_basic_set_mode(CU_BRM_VERBOSE);

    if (argc == 1) {
        CU_basic_run_tests();
    } else if (argc == 2) {
        suite = CU_get_suite_by_name("msprime", CU_get_registry());
        if (suite == NULL) {
            printf("Suite not found\n");
            return EXIT_FAILURE;
        }
        test = CU_get_test_by_name(argv[1], suite);
        if (test == NULL) {
            printf("Test '%s' not found\n", argv[1]);
            return EXIT_FAILURE;
        }
        CU_basic_run_test(suite, test);
    } else {
        printf("usage: ./simulation_tests <test_name>\n");
        return EXIT_FAILURE;
    }

    ret = EXIT_SUCCESS;
    if (CU_get_number_of_tests_failed() != 0) {
        printf("Test failed!\n");
        ret = EXIT_FAILURE;
    }
    CU_cleanup_registry();
    return ret;
}
