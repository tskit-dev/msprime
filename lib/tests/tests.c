/*
** Copyright (C) 2016-2020 University of Oxford
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
#include "likelihood.h"

#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <CUnit/Basic.h>

#define ALPHABET_BINARY 0
#define ALPHABET_NUCLEOTIDE 1

/* Global variables used for test in state in the test suite */

char *_tmp_file_name;
FILE *_devnull;

static void
verify_simulator_tsk_treeseq_equality(msp_t *msp, tsk_treeseq_t *tree_seq, double scale)
{
    int ret;
    uint32_t num_samples = msp_get_num_samples(msp);
    uint32_t j;
    tsk_node_t node;
    sample_t *samples;
    node_id_t *sample_ids;

    CU_ASSERT_EQUAL_FATAL(
        tsk_treeseq_get_num_samples(tree_seq), msp_get_num_samples(msp));
    CU_ASSERT_EQUAL_FATAL(tsk_treeseq_get_num_edges(tree_seq), msp_get_num_edges(msp));
    CU_ASSERT_EQUAL_FATAL(
        tsk_treeseq_get_num_migrations(tree_seq), msp_get_num_migrations(msp));
    ret = msp_get_samples(msp, &samples);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_treeseq_get_num_nodes(tree_seq) >= num_samples);
    CU_ASSERT_TRUE(tsk_migration_table_equals(
        &msp->tables->migrations, &tree_seq->tables->migrations))

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
    tsk_treeseq_print_state(tree_seq, _devnull);
}

/* Simple unit tests for the Fenwick tree API. */
static void
test_fenwick(void)
{
    fenwick_t t;
    double s;
    size_t j, n;

    for (n = 1; n < 100; n++) {
        s = 0;
        CU_ASSERT(fenwick_alloc(&t, n) == 0);
        for (j = 1; j <= n; j++) {
            fenwick_increment(&t, j, j);
            s = s + j;
            CU_ASSERT(fenwick_get_value(&t, j) == j);
            CU_ASSERT(fenwick_get_cumulative_sum(&t, j) == s);
            CU_ASSERT(fenwick_get_total(&t) == s);
            CU_ASSERT(fenwick_get_numerical_drift(&t) == 0.0);
            CU_ASSERT(fenwick_find(&t, s) == j);
            fenwick_set_value(&t, j, 0);
            CU_ASSERT(fenwick_get_value(&t, j) == 0);
            CU_ASSERT(fenwick_get_cumulative_sum(&t, j) == s - j);
            fenwick_set_value(&t, j, j);
            CU_ASSERT(fenwick_get_value(&t, j) == j);
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
        CU_ASSERT_EQUAL(memcmp(t1.tree, t2.tree, (t2.size + 1) * sizeof(*t1.tree)), 0);
        CU_ASSERT_EQUAL(
            memcmp(t1.values, t2.values, (t2.size + 1) * sizeof(*t2.values)), 0);
        CU_ASSERT(fenwick_free(&t1) == 0);
        CU_ASSERT(fenwick_free(&t2) == 0);
    }
}

static void
test_fenwick_zero_values(void)
{
    fenwick_t t;
    size_t n = 10;
    size_t j;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(rng != 0);
    gsl_rng_set(rng, 42);
    CU_ASSERT(fenwick_alloc(&t, n) == 0);

    /* Adding in lots of small values is fine, and we can recover these
     * to a high degree of precision */
    for (j = 0; j < 1000; j++) {
        fenwick_set_value(&t, j % n + 1, gsl_ran_flat(rng, 0, 1e-12));
        fenwick_verify(&t, 1e-9);
    }
    fenwick_print_state(&t, _devnull);

    /* Set everything before 4 to zero */
    fenwick_set_value(&t, 1, 0);
    fenwick_set_value(&t, 2, 0);
    fenwick_set_value(&t, 3, 0);
    fenwick_set_value(&t, 4, 0);
    /* Because of numerical precision issues, the internal node 4 will not
     * compute to *exactly* zero in the tree. */
    CU_ASSERT_FATAL(t.tree[4] > 0);

    /* 5 is the first non-zero node in the tree, so any values smaller
     * than this should search to it. */
    CU_ASSERT_EQUAL(fenwick_find(&t, t.values[5]), 5);
    CU_ASSERT_EQUAL(fenwick_find(&t, DBL_EPSILON), 5);
    CU_ASSERT_EQUAL(fenwick_find(&t, DBL_MIN), 5);
    CU_ASSERT_EQUAL(fenwick_find(&t, 0), 5);

    /* Set the remaining values to zero and search */
    for (j = 5; j <= n; j++) {
        fenwick_set_value(&t, j, 0);
    }
    CU_ASSERT_EQUAL(fenwick_find(&t, 1), n + 1);
    CU_ASSERT_EQUAL(fenwick_find(&t, 0), n + 1);

    fenwick_free(&t);
    gsl_rng_free(rng);
}

static void
test_fenwick_drift(void)
{
    fenwick_t t;
    double s;
    size_t j, n;

    for (n = 1; n < 100; n++) {
        s = 0;
        CU_ASSERT(fenwick_alloc(&t, n) == 0);
        for (j = 1; j <= n; j++) {
            fenwick_increment(&t, j, j);
            s = s + j;
            CU_ASSERT(fenwick_get_value(&t, j) == j);
            CU_ASSERT(fenwick_get_cumulative_sum(&t, j) == s);
            CU_ASSERT(fenwick_get_total(&t) == s);
            CU_ASSERT(fenwick_get_numerical_drift(&t) == 0.0);
        }
        /* put some drift into the tree */
        for (j = 1; j <= n; j++) {
            t.tree[j] += 1e-9;
        }
        fenwick_print_state(&t, _devnull);
        CU_ASSERT(fenwick_get_numerical_drift(&t) > 0);
        fenwick_rebuild(&t);
        CU_ASSERT(fenwick_get_numerical_drift(&t) == 0);

        s = 0;
        for (j = 1; j <= n; j++) {
            s = s + j;
            CU_ASSERT(fenwick_get_value(&t, j) == j);
            CU_ASSERT(fenwick_get_cumulative_sum(&t, j) == s);
            CU_ASSERT(fenwick_get_numerical_drift(&t) == 0.0);
        }

        CU_ASSERT(fenwick_free(&t) == 0);
    }
}

static void
test_single_locus_two_populations(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = { { 0, 0.0 }, { 0, 0.0 }, { 1, 40.0 } };
    tsk_table_collection_t tables;
    tsk_edge_table_t *edges;
    tsk_node_table_t *nodes;
    tsk_migration_table_t *migrations;
    size_t num_edges, num_migrations;
    uint32_t n = 3;
    double t0 = 30.0;
    double t1 = 30.5;
    double t2 = 40.5;
    recomb_map_t recomb_map;

    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 0, false);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
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
    msp_verify(&msp, 0);
    msp_print_state(&msp, _devnull);

    nodes = &msp.tables->nodes;
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
    edges = &msp.tables->edges;
    CU_ASSERT_EQUAL(edges->parent[0], 3);
    CU_ASSERT_EQUAL(edges->parent[1], 3);
    CU_ASSERT_EQUAL(edges->parent[2], 4);
    CU_ASSERT_EQUAL(edges->parent[3], 4);

    num_migrations = msp_get_num_migrations(&msp);
    CU_ASSERT_EQUAL_FATAL(num_migrations, 3);

    migrations = &msp.tables->migrations;
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
    tsk_table_collection_free(&tables);
}

static void
test_single_locus_many_populations(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    uint32_t num_populations = 500;
    sample_t samples[] = { { 0, 0.0 }, { num_populations - 1, 0.0 } };
    uint32_t n = 2;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1, 0, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
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
    ret = msp_run(&msp, DBL_MAX, 1);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_EVENTS);
    msp_verify(&msp, 0);
    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);
    msp_print_state(&msp, _devnull);
    CU_ASSERT_EQUAL_FATAL(msp_get_num_edges(&msp), 2);
    CU_ASSERT_EQUAL_FATAL(msp_get_num_nodes(&msp), 3);
    CU_ASSERT_EQUAL(msp.tables->edges.parent[0], 2);
    CU_ASSERT_EQUAL(msp.tables->edges.parent[1], 2);

    CU_ASSERT_TRUE(msp.tables->nodes.time[2] > 30.0);
    CU_ASSERT_EQUAL(msp.tables->nodes.population[2], num_populations - 1);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_dtwf_simultaneous_historical_samples(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = { { 0, 0 }, { 0, 1.1 }, { 0, 1.2 } };
    tsk_node_table_t *nodes;
    uint32_t n = 3;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    gsl_rng_set(rng, 5);

    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 0, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_table_collection_clear(&tables);
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 100, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);
    msp_print_state(&msp, _devnull);

    nodes = &msp.tables->nodes;
    CU_ASSERT_EQUAL(nodes->time[0], 0);
    CU_ASSERT_EQUAL(nodes->time[1], 1.1);
    CU_ASSERT_EQUAL(nodes->time[2], 1.2);
    CU_ASSERT_EQUAL_FATAL(msp_get_num_nodes(&msp), 5);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
    recomb_map_free(&recomb_map);
}

static void
test_single_locus_historical_sample(void)
{
    int ret, j;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = { { 0, 0.0 }, { 0, 101.0 } };
    tsk_edge_table_t *edges;
    tsk_node_table_t *nodes;
    uint32_t n = 2;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    gsl_rng_set(rng, 5);

    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 0, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < 2; j++) {
        tsk_table_collection_clear(&tables);
        ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
        CU_ASSERT_EQUAL(ret, 0);

        if (j == 0) {
            ret = msp_set_simulation_model_hudson(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            ret = msp_set_simulation_model_dtwf(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        }
        ret = msp_set_population_configuration(&msp, 0, 100, 0);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_initialise(&msp);
        CU_ASSERT_EQUAL(ret, 0);

        msp_print_state(&msp, _devnull);
        ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        msp_verify(&msp, 0);
        msp_print_state(&msp, _devnull);

        CU_ASSERT_EQUAL_FATAL(msp_get_num_nodes(&msp), 3);
        nodes = &msp.tables->nodes;
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(nodes->time[0], 0);
        CU_ASSERT_EQUAL(nodes->time[1], 101.0);
        CU_ASSERT_TRUE(nodes->time[2] > 101.0);

        CU_ASSERT_EQUAL_FATAL(msp_get_num_edges(&msp), 2);
        edges = &msp.tables->edges;
        CU_ASSERT_EQUAL(edges->left[0], 0);
        CU_ASSERT_EQUAL(edges->right[0], 1);
        CU_ASSERT_EQUAL(edges->parent[0], 2);
        CU_ASSERT_EQUAL(edges->left[1], 0);
        CU_ASSERT_EQUAL(edges->right[1], 1);
        CU_ASSERT_EQUAL(edges->parent[1], 2);

        ret = msp_free(&msp);
        CU_ASSERT_EQUAL(ret, 0);
    }
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
    recomb_map_free(&recomb_map);
}

static void
test_single_locus_multiple_historical_samples(void)
{
    int ret, j;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = { { 0, 0.0 }, { 0, 10.0 }, { 0, 10.0 }, { 0, 10.0 } };
    /* tsk_edge_table_t *edges; */
    /* tsk_node_table_t *nodes; */
    uint32_t n = 4;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    gsl_rng_set(rng, 5);

    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 0, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < 2; j++) {
        tsk_table_collection_clear(&tables);
        ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
        CU_ASSERT_EQUAL(ret, 0);

        if (j == 0) {
            ret = msp_set_simulation_model_hudson(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            ret = msp_set_simulation_model_dtwf(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        }
        ret = msp_initialise(&msp);
        CU_ASSERT_EQUAL(ret, 0);

        msp_print_state(&msp, _devnull);
        ret = msp_run(&msp, 10, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_TIME);
        CU_ASSERT_EQUAL(msp_get_num_ancestors(&msp), 1);
        msp_verify(&msp, 0);
        ret = msp_run(&msp, DBL_MAX, 1);
        CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_EVENTS);
        /* All the samples should be added in now. */
        CU_ASSERT_EQUAL(msp_get_num_ancestors(&msp), 4);
        msp_verify(&msp, 0);

        ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        msp_verify(&msp, 0);

        ret = msp_free(&msp);
        CU_ASSERT_EQUAL(ret, 0);
    }
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
    recomb_map_free(&recomb_map);
}

static void
test_single_locus_historical_sample_start_time(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = { { 0, 0.0 }, { 0, 10.0 } };
    tsk_edge_table_t *edges;
    tsk_node_table_t *nodes;
    uint32_t n = 2;
    size_t j, model;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;
    double start_times[] = { 0, 2, 10, 10.0001, 1000 };
    /* double start_times[] = {0, 2, 9.99}; //10, 1000}; */
    /* double start_times[] = {10.00}; //10, 1000}; */

    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 0, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (model = 0; model < 2; model++) {
        for (j = 0; j < sizeof(start_times) / sizeof(double); j++) {

            ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
            CU_ASSERT_EQUAL(ret, 0);
            if (model == 0) {
                ret = msp_set_simulation_model_hudson(&msp);
                CU_ASSERT_EQUAL(ret, 0);
            } else {
                ret = msp_set_simulation_model_dtwf(&msp);
                CU_ASSERT_EQUAL(ret, 0);
            }
            ret = msp_set_population_configuration(&msp, 0, 100, 0);
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
            msp_verify(&msp, 0);
            /* msp_print_state(&msp, _devnull); */

            CU_ASSERT_EQUAL_FATAL(msp_get_num_nodes(&msp), 3);
            nodes = &msp.tables->nodes;
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(nodes->time[0], 0);
            CU_ASSERT_EQUAL(nodes->time[1], 10);
            CU_ASSERT_TRUE(nodes->time[2] > 10);
            CU_ASSERT_TRUE(nodes->time[2] > start_times[j]);

            CU_ASSERT_EQUAL_FATAL(msp_get_num_edges(&msp), 2);
            edges = &msp.tables->edges;
            CU_ASSERT_EQUAL(edges->left[0], 0);
            CU_ASSERT_EQUAL(edges->right[0], 1);
            CU_ASSERT_EQUAL(edges->parent[0], 2);
            CU_ASSERT_EQUAL(edges->left[1], 0);
            CU_ASSERT_EQUAL(edges->right[1], 1);
            CU_ASSERT_EQUAL(edges->parent[1], 2);

            ret = msp_free(&msp);
            CU_ASSERT_EQUAL(ret, 0);
            ret = tsk_table_collection_clear(&tables);
            CU_ASSERT_EQUAL(ret, 0);
        }
    }
    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_single_locus_historical_sample_end_time(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = { { 0, 0.0 }, { 0, 10.0 } };
    tsk_node_table_t *nodes;
    uint32_t n = 2;
    size_t model;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 0, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (model = 0; model < 2; model++) {
        ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
        CU_ASSERT_EQUAL(ret, 0);
        if (model == 0) {
            ret = msp_set_simulation_model_hudson(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            ret = msp_set_simulation_model_dtwf(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        }
        ret = msp_set_population_configuration(&msp, 0, 100, 0);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_initialise(&msp);
        CU_ASSERT_EQUAL(ret, 0);

        msp_print_state(&msp, _devnull);
        ret = msp_run(&msp, 10.0, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_TIME);
        /* msp_verify(&msp); */

        /* The sampling event should *not* have been applied */
        CU_ASSERT_EQUAL_FATAL(msp.next_sampling_event, 0);
        CU_ASSERT_EQUAL_FATAL(msp_get_num_ancestors(&msp), 1);
        CU_ASSERT_EQUAL_FATAL(msp_get_num_nodes(&msp), 2);

        nodes = &msp.tables->nodes;
        CU_ASSERT_EQUAL(nodes->time[0], 0);
        CU_ASSERT_EQUAL(nodes->time[1], 10);
        CU_ASSERT_EQUAL_FATAL(msp_get_num_edges(&msp), 0);

        /* msp_print_state(&msp, stdout); */

        ret = msp_finalise_tables(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(tables.nodes.num_rows, 3);
        CU_ASSERT_EQUAL(tables.edges.num_rows, 1);

        ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 0);

        ret = msp_finalise_tables(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(tables.nodes.num_rows, 4);
        CU_ASSERT_EQUAL(tables.edges.num_rows, 3);

        ret = msp_free(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tsk_table_collection_clear(&tables);
        CU_ASSERT_EQUAL(ret, 0);
    }
    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
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
    double migration_matrix[] = { 0, 0, 0, 0 };
    double matrix[4], growth_rate, initial_size;
    double Ne = 4;
    size_t migration_events[4];
    size_t breakpoints[m];
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    ret = recomb_map_alloc_uniform(&recomb_map, m, 1, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < n; j++) {
        samples[j].time = j;
        samples[j].population_id = j % 2;
    }
    CU_ASSERT_EQUAL(msp_alloc(&msp, 0, NULL, NULL, NULL, rng), MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);
    CU_ASSERT_EQUAL(
        msp_alloc(&msp, n, samples, NULL, NULL, rng), MSP_ERR_BAD_PARAM_VALUE);
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
    CU_ASSERT_EQUAL(
        msp_alloc(&msp, n, samples, &recomb_map, &tables, rng), MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);
    samples[0].time = 0.0;

    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(msp_set_dimensions(&msp, 0, 1), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_dimensions(&msp, 1, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_node_mapping_block_size(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_segment_block_size(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_avl_node_block_size(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_num_populations(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_population_configuration(&msp, -1, 0, 0),
        MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL(msp_set_population_configuration(&msp, 3, 0, 0),
        MSP_ERR_POPULATION_OUT_OF_BOUNDS);

    ret = msp_set_simulation_model_hudson(&msp);
    CU_ASSERT_EQUAL(msp_get_model(&msp)->type, MSP_MODEL_HUDSON);
    ret = msp_set_population_configuration(&msp, 0, Ne, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_set_num_populations(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(msp_get_population_configuration(&msp, 3, NULL, NULL),
        MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    ret = msp_set_population_configuration(&msp, 0, 2 * Ne, 0.5);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(
        msp_set_migration_matrix(&msp, 0, NULL), MSP_ERR_BAD_MIGRATION_MATRIX);
    CU_ASSERT_EQUAL(msp_set_migration_matrix(&msp, 3, migration_matrix),
        MSP_ERR_BAD_MIGRATION_MATRIX);
    migration_matrix[0] = 1;
    CU_ASSERT_EQUAL(msp_set_migration_matrix(&msp, 4, migration_matrix),
        MSP_ERR_BAD_MIGRATION_MATRIX);
    migration_matrix[0] = 0;
    migration_matrix[1] = -1;
    CU_ASSERT_EQUAL(msp_set_migration_matrix(&msp, 4, migration_matrix),
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
    // TODO CU_ASSERT_EQUAL(msp_get_recombination_rate(&msp), 1.0);

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
    tsk_table_collection_free(&tables);
}

static void
test_demographic_events(void)
{
    int ret;
    uint32_t j, k, model;
    uint32_t n = 10;
    uint32_t m = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    double migration_matrix[] = { 0, 1, 1, 0 };
    double last_time, time, pop_size;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, m, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (model = 0; model < 2; model++) {
        for (j = 0; j < n; j++) {
            samples[j].time = j;
            samples[j].population_id = j % 2;
        }

        tsk_table_collection_clear(&tables);
        ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
        CU_ASSERT_EQUAL(ret, 0);

        /* Zero or negative population sizes are not allowed */
        ret = msp_set_population_configuration(&msp, 0, -1, 0);
        CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
        ret = msp_set_population_configuration(&msp, 0, 0, 0);
        CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

        ret = msp_set_num_populations(&msp, 2);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_population_configuration(&msp, 0, 1, 0.001);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_population_configuration(&msp, 1, 2, 0.002);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_migration_matrix(&msp, 4, migration_matrix);
        CU_ASSERT_EQUAL(ret, 0);

        if (model == 0) {
            ret = msp_set_simulation_model_hudson(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            ret = msp_set_simulation_model_dtwf(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        }

        ret = msp_set_population_configuration(&msp, 0, 100, 0);
        CU_ASSERT_EQUAL(ret, 0);

        CU_ASSERT_EQUAL(msp_add_mass_migration(&msp, 10, -1, 0, 1),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        CU_ASSERT_EQUAL(
            msp_add_mass_migration(&msp, 10, 2, 0, 1), MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        CU_ASSERT_EQUAL(
            msp_add_mass_migration(&msp, 10, 0, 0, 1), MSP_ERR_SOURCE_DEST_EQUAL);
        CU_ASSERT_EQUAL(
            msp_add_mass_migration(&msp, 10, 0, 1, -5), MSP_ERR_BAD_PROPORTION);

        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, -1, 0, 2.0),
            MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, 0, -1, 2.0),
            MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, 0, 2, 2.0),
            MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, 2, 0, 2.0),
            MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, -1, 0, 2.0),
            MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, 0, -1, 2.0),
            MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, -1, -1, -2.0),
            MSP_ERR_BAD_PARAM_VALUE);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, 0, 0, 2.0),
            MSP_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX);

        CU_ASSERT_EQUAL(msp_add_population_parameters_change(&msp, 10, -2, 0, 0),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        CU_ASSERT_EQUAL(msp_add_population_parameters_change(&msp, 10, -1, -1, 0),
            MSP_ERR_BAD_PARAM_VALUE);
        CU_ASSERT_EQUAL(
            msp_add_population_parameters_change(&msp, 10, -1, GSL_NAN, GSL_NAN),
            MSP_ERR_BAD_PARAM_VALUE);

        CU_ASSERT_EQUAL(msp_add_simple_bottleneck(&msp, 10, -1, 0),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        CU_ASSERT_EQUAL(
            msp_add_simple_bottleneck(&msp, 10, 0, -1), MSP_ERR_BAD_PROPORTION);
        CU_ASSERT_EQUAL(
            msp_add_simple_bottleneck(&msp, 10, 0, 1.1), MSP_ERR_BAD_PROPORTION);

        CU_ASSERT_EQUAL(msp_add_instantaneous_bottleneck(&msp, 10, 2, 0),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        CU_ASSERT_EQUAL(
            msp_add_simple_bottleneck(&msp, 10, 0, -1), MSP_ERR_BAD_PROPORTION);
        CU_ASSERT_EQUAL_FATAL(msp_add_simple_bottleneck(&msp, 10, -1, 0),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);

        CU_ASSERT_EQUAL(msp_add_census_event(&msp, -0.5), MSP_ERR_BAD_PARAM_VALUE);

        ret = msp_add_census_event(&msp, 0.05);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_add_mass_migration(&msp, 0.1, 0, 1, 0.5);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_add_migration_rate_change(&msp, 0.2, 0, 1, 2.0);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_add_migration_rate_change(&msp, 0.3, -1, -1, 3.0);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_add_population_parameters_change(&msp, 0.4, 1, 1, GSL_NAN);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_add_population_parameters_change(&msp, 0.5, 0, 0.5, 0.002);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_add_population_parameters_change(&msp, 0.6, -1, 0.5, 0.001);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_add_population_parameters_change(&msp, 0.7, -1, 1, GSL_NAN);
        CU_ASSERT_EQUAL(ret, 0);

        if (model == 0) {
            ret = msp_add_population_parameters_change(&msp, 0.7, 0, GSL_NAN, 0);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_add_simple_bottleneck(&msp, 1.5, 0, 0.5);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_add_instantaneous_bottleneck(&msp, 2.5, 0, 2.0);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            ret = msp_add_population_parameters_change(&msp, 0.8, -1, 1, 0);
            CU_ASSERT_EQUAL(ret, 0);
            /* Need to lower final migration rate for DTWF or else lineages will
             * alternate pops every generation and miss each other - need to let
             * one lineage migrate while the others stay put */
            ret = msp_add_migration_rate_change(&msp, 1.5, -1, -1, 0.3);
            CU_ASSERT_EQUAL(ret, 0);
            /* Bottleneck events only supported in Hudson model so we add
             * another mass migration to have the same number of events */
            ret = msp_add_mass_migration(&msp, 2.5, 1, 0, 0.6);
            CU_ASSERT_EQUAL(ret, 0);
        }

        CU_ASSERT_EQUAL(msp_add_mass_migration(&msp, 0.1, 0, 1, 0.5),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 0.2, 0, 1, 2.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
        CU_ASSERT_EQUAL(msp_add_population_parameters_change(&msp, 0.4, 0, 0.5, 1.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);

        CU_ASSERT_EQUAL(msp_add_simple_bottleneck(&msp, 0.7, 0, 1.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
        CU_ASSERT_EQUAL(msp_add_instantaneous_bottleneck(&msp, 0.8, 0, 1.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);

        CU_ASSERT_EQUAL(msp_debug_demography(&msp, &time), MSP_ERR_BAD_STATE);

        ret = msp_initialise(&msp);
        CU_ASSERT_EQUAL(ret, 0);

        j = 0;
        last_time = 0;
        do {
            ret = msp_debug_demography(&msp, &time);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
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
                ret = msp_compute_population_size(
                    &msp, k, last_time + (time - last_time) / 2, &pop_size);
                CU_ASSERT_EQUAL(ret, 0);
                CU_ASSERT_TRUE(pop_size >= 0);
            }
            j++;
            last_time = time;
        } while (!gsl_isinf(time));
        CU_ASSERT_TRUE(j >= 9);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(msp_run(&msp, DBL_MAX, ULONG_MAX), MSP_ERR_BAD_STATE);
        ret = msp_reset(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        msp_print_state(&msp, _devnull);
        ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_free(&msp);
        CU_ASSERT_EQUAL(ret, 0);
    }

    free(samples);
    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_census_event(void)
{
    int ret;
    uint32_t n = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;
    int num_census_nodes = 0;
    int i;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);

    /* Add a census event in at 0.5 generations. */
    ret = msp_add_census_event(msp, 0.5);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(msp, 0);
    msp_print_state(msp, _devnull);

    /* Check there is more than 1 node at the census time. */
    for (i = 0; i < tables.nodes.num_rows; i++) {
        if (tables.nodes.time[i] == 0.5) {
            num_census_nodes++;
        }
    }
    CU_ASSERT_TRUE(num_census_nodes > 1);

    /* Free things. */
    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_dtwf_events_between_generations(void)
{
    int ret, i, j;
    uint32_t n = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t msp;
    population_t pop;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    double migration_matrix[] = { 0, 0, 0, 0 };
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    for (j = 0; j < n; j++) {
        samples[j].time = 0;
        samples[j].population_id = j % 2;
    }

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_set_num_populations(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 10, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 1, 10, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_migration_matrix(&msp, 4, migration_matrix);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_add_population_parameters_change(&msp, 0.1, 0, 5, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.2, 1, 5, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 1.1, 0, 10, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 1.2, 1, 10, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_migration_rate_change(&msp, 3, -1, -1, 0.3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Check that we stay on integer-valued times, and that both demographic
     * events between generations occurred  */
    ret = msp_run(&msp, DBL_MAX, 1);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(msp.time, 1);
    for (i = 0; i < 2; i++) {
        pop = msp.populations[i];
        CU_ASSERT_EQUAL(pop.initial_size, 5);
    }

    ret = msp_run(&msp, DBL_MAX, 1);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(msp.time, 2);
    for (i = 0; i < 2; i++) {
        pop = msp.populations[i];
        CU_ASSERT_EQUAL(pop.initial_size, 10);
    }

    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    gsl_rng_free(rng);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_dtwf_zero_pop_size(void)
{
    int ret;
    uint32_t n = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    memset(samples, 0, n * sizeof(sample_t));

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* DTWF population size must round to >= 1 */
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 0.4, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_DTWF_ZERO_POPULATION_SIZE);
    msp_free(&msp);

    /* With no high growth rate, population sizes crash to zero. */
    tsk_table_collection_clear(&tables);
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 10, 100);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_DTWF_ZERO_POPULATION_SIZE);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    gsl_rng_free(rng);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_demographic_events_start_time(void)
{
    int ret;
    uint32_t n = 10;
    uint32_t model;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (model = 0; model < 1; model++) {
        memset(samples, 0, n * sizeof(sample_t));
        ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
        CU_ASSERT_EQUAL(ret, 0);

        if (model == 0) {
            ret = msp_set_start_time(msp, 1.0);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_add_population_parameters_change(msp, 0.4, 0, 0.5, 1.0);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_set_simulation_model_hudson(msp);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            /* Times need to be bumped to integer values for DTWF */
            ret = msp_set_start_time(msp, 2.0);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_add_population_parameters_change(msp, 1, 0, 0.5, 1.0);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_set_simulation_model_dtwf(msp);
            CU_ASSERT_EQUAL(ret, 0);
        }

        ret = msp_initialise(msp);
        CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_DEMOGRAPHIC_EVENT_TIME);
        ret = msp_free(msp);
        CU_ASSERT_EQUAL(ret, 0);
    }

    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_dtwf_unsupported_bottleneck(void)
{
    int ret;
    uint32_t n = 10;
    uint32_t j;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < n; j++) {
        samples[j].time = 0;
        samples[j].population_id = 0;
    }

    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_add_simple_bottleneck(&msp, 0.8, 0, 0.5);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_DTWF_UNSUPPORTED_BOTTLENECK);

    ret = msp_reset(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_hudson(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_instantaneous_bottleneck(&msp, 0.8, 0, 0.5);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_DTWF_UNSUPPORTED_BOTTLENECK);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    gsl_rng_free(rng);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_time_travel_error(void)
{
    int ret;
    uint32_t n = 100;
    sample_t *samples = calloc(n, sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
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
    tsk_table_collection_free(&tables);
}

static int
get_num_children(size_t node, tsk_edge_table_t *edges)
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
test_pedigree_multi_locus_simulation(void)
{
    int ret;
    const char *model_name;
    /* int i; */
    /* int num_coalescent_events = 0; */
    /* int num_children; */
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;
    int num_inds = 4;
    int ploidy = 2;
    tsk_id_t inds[4] = { 1, 2, 3, 4 };
    tsk_id_t parents[8] = { 2, 3, 2, 3, -1, -1, -1, -1 }; // size num_inds * ploidy
    double times[4] = { 0, 0, 1, 1 };
    tsk_flags_t is_sample[4] = { 1, 1, 0, 0 };
    uint32_t N = 4;
    uint32_t n = 2 * ploidy;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 100.0, 10.0, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_alloc_pedigree(msp, num_inds, ploidy);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_pedigree(msp, num_inds, inds, parents, times, is_sample);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(msp, 0, N, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_wf_ped(msp);
    CU_ASSERT_EQUAL(ret, 0);
    model_name = msp_get_model_name(msp);
    CU_ASSERT_STRING_EQUAL(model_name, "wf_ped");

    ret = msp_run(msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(msp, 0);
    ret = msp_set_simulation_model_dtwf(msp);
    CU_ASSERT_EQUAL(ret, 0);
    model_name = msp_get_model_name(msp);
    CU_ASSERT_STRING_EQUAL(model_name, "dtwf");
    msp_print_pedigree_inds(msp, _devnull);
    ret = msp_run(msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_finalise_tables(msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tables.migrations.num_rows, 0);
    CU_ASSERT(tables.nodes.num_rows > 0);
    CU_ASSERT(tables.edges.num_rows > 0);

    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_pedigree_specification(void)
{
    int ret;
    const char *model_name;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;
    int num_inds = 4;
    int ploidy = 2;
    tsk_id_t inds[4] = { 1, 2, 3, 4 };
    tsk_id_t parents[8] = { 2, 3, 2, 3, -1, -1, -1, -1 }; // size num_inds * ploidy
    double times_good[4] = { 0, 0, 1, 1 };
    tsk_flags_t is_sample[4] = { 1, 1, 0, 0 };
    uint32_t n = 2 * ploidy;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_alloc_pedigree(msp, num_inds, ploidy);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_pedigree(msp, num_inds + 1, inds, parents, times_good, is_sample);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = msp_set_pedigree(msp, num_inds, inds, parents, times_good, is_sample);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_wf_ped(msp);
    CU_ASSERT_EQUAL(ret, 0);
    model_name = msp_get_model_name(msp);
    CU_ASSERT_STRING_EQUAL(model_name, "wf_ped");

    ret = msp_run(msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(msp, 0);

    ret = msp_finalise_tables(msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);

    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_pedigree_single_locus_simulation(void)
{
    int ret;
    const char *model_name;
    int i;
    int num_coalescent_events = 0;
    int num_children;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;
    int num_inds = 4;
    int ploidy = 2;
    tsk_id_t inds[4] = { 1, 2, 3, 4 };
    tsk_id_t bad_inds[4] = { 0, 1, 2, 3 };
    tsk_id_t parents[8] = { 2, 3, 2, 3, -1, -1, -1, -1 }; // size num_inds * ploidy
    double times[4] = { 0, 0, 1, 1 };
    tsk_flags_t is_sample[4] = { 1, 1, 0, 0 };
    tsk_flags_t bad_is_sample[4] = { 1, 0, 0, 0 };
    uint32_t n = 2 * ploidy;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_alloc_pedigree(msp, num_inds, ploidy);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_pedigree(msp, num_inds, bad_inds, parents, times, is_sample);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PEDIGREE_ID);
    ret = msp_set_pedigree(msp, num_inds, inds, parents, times, bad_is_sample);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PEDIGREE_ID);
    ret = msp_set_pedigree(msp, num_inds, inds, parents, times, is_sample);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_wf_ped(msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(msp, 0, n, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);
    model_name = msp_get_model_name(msp);
    CU_ASSERT_STRING_EQUAL(model_name, "wf_ped");

    ret = msp_run(msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(msp, 0);
    ret = msp_run(msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_STATE);
    ret = msp_set_simulation_model_dtwf(msp);
    CU_ASSERT_EQUAL(ret, 0);
    model_name = msp_get_model_name(msp);
    CU_ASSERT_STRING_EQUAL(model_name, "dtwf");
    ret = msp_run(msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(msp, 0);

    /* For the single locus sim we should have n-1 coalescent events,
     * counting multiple mergers as multiple coalescent events */
    for (i = 0; i < msp->tables->nodes.num_rows; i++) {
        num_children = get_num_children(i, &msp->tables->edges);
        if (num_children > 0) {
            num_coalescent_events += num_children - 1;
        }
    }
    CU_ASSERT_EQUAL(num_coalescent_events, n - 1);

    ret = msp_finalise_tables(msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tables.migrations.num_rows, 0);
    CU_ASSERT(tables.nodes.num_rows > 0);
    CU_ASSERT(tables.edges.num_rows > 0);

    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
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
    tsk_table_collection_t tables[2];
    recomb_map_t recomb_map;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, m, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    for (j = 0; j < 2; j++) {
        ret = tsk_table_collection_init(&tables[j], 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        gsl_rng_set(rng, seed);
        ret = msp_alloc(msp, n, samples, &recomb_map, &tables[j], rng);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_simulation_model_dtwf(msp);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_population_configuration(msp, 0, n, 0);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_initialise(msp);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_run(msp, DBL_MAX, UINT32_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        msp_verify(msp, 0);
        ret = msp_finalise_tables(msp);
        CU_ASSERT_EQUAL(ret, 0);
        msp_free(msp);
        CU_ASSERT_EQUAL(tables[j].migrations.num_rows, 0);
        CU_ASSERT(tables[j].nodes.num_rows > 0);
        CU_ASSERT(tables[j].edges.num_rows > 0);
    }
    CU_ASSERT_TRUE(tsk_node_table_equals(&tables[0].nodes, &tables[1].nodes));
    CU_ASSERT_TRUE(tsk_edge_table_equals(&tables[0].edges, &tables[1].edges));

    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    for (j = 0; j < 2; j++) {
        tsk_table_collection_free(&tables[j]);
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
    tsk_table_collection_t tables;
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
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 10, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_populations(msp, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(msp);
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
        msp_verify(msp, 0);
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
            ret = msp_set_simulation_model_dtwf(msp);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            model = msp_get_model(msp)->type;
            CU_ASSERT_EQUAL(model, MSP_MODEL_DTWF);
            model_name = msp_get_model_name(msp);
            CU_ASSERT_STRING_EQUAL(model_name, "dtwf");
            ret = msp_set_simulation_model_hudson(msp);
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
    tables.sequence_length = msp->sequence_length;
    CU_ASSERT_EQUAL_FATAL(tables.sequence_length, msp->sequence_length);
    ret = tsk_treeseq_init(&ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_print_state(&ts, _devnull);
    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
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
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(msp);
    CU_ASSERT_EQUAL(ret, 0);
    model_name = msp_get_model_name(msp);
    CU_ASSERT_STRING_EQUAL(model_name, "dtwf");

    ret = msp_run(msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(msp, 0);

    /* For the single locus sim we should have n-1 coalescent events,
     * counting multiple mergers as multiple coalescent events */
    for (i = 0; i < msp->tables->nodes.num_rows; i++) {
        num_children = get_num_children(i, &msp->tables->edges);
        if (num_children > 0) {
            num_coalescent_events += num_children - 1;
        }
    }
    CU_ASSERT_EQUAL(num_coalescent_events, n - 1);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_dtwf_low_recombination(void)
{
    int ret;
    uint32_t n = 2;

    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1e-9, false);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(msp, 0, n, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    /* msp_verify(msp); */

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
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
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    /* For the single locus sim we should have exactly n - 1 events */
    for (j = 0; j < n - 2; j++) {
        ret = msp_run(msp, DBL_MAX, 1);
        CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_EVENTS);
        msp_verify(msp, 0);
    }
    ret = msp_run(msp, DBL_MAX, 1);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(msp, 0);

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
    tsk_table_collection_free(&tables);
}

static void
test_single_locus_gene_conversion(void)
{
    int ret;
    uint32_t j;
    uint32_t n = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_gene_conversion_rate(msp, 10, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(msp_get_gene_conversion_rate(msp), 10);

    /* For the single locus sim we should have exactly n - 1 events */
    for (j = 0; j < n - 2; j++) {
        ret = msp_run(msp, DBL_MAX, 1);
        CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_EVENTS);
        msp_verify(msp, 0);
    }
    ret = msp_run(msp, DBL_MAX, 1);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(msp, 0);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_multi_locus_bottleneck_arg(void)
{
    int ret;
    uint32_t n = 2;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 100.0, 10.0, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_instantaneous_bottleneck(msp, 1.0, 0, 10);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_store_full_arg(msp, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(msp, 0);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_floating_point_extremes(void)
{
    int ret;
    uint32_t n = 2;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;
    double denormal_min = pow(2, -52) * pow(2, -1022);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    memset(samples, 0, n * sizeof(sample_t));
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = recomb_map_alloc_uniform(&recomb_map, denormal_min, DBL_MAX, false);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(msp, 0, DBL_MAX, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BREAKPOINT_RESAMPLE_OVERFLOW);
    msp_print_state(msp, _devnull);
    recomb_map_free(&recomb_map);
    msp_free(msp);

    tsk_table_collection_clear(&tables);
    /* A long sequence length and high recombination should overflow */
    ret = recomb_map_alloc_uniform(&recomb_map, DBL_MAX, DBL_MAX, false);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);
    msp_print_state(msp, _devnull);
    ret = msp_run(msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BREAKPOINT_MASS_NON_FINITE);

    msp_free(msp);
    recomb_map_free(&recomb_map);

    gsl_rng_free(rng);
    free(msp);
    free(samples);
    tsk_table_collection_free(&tables);
}

static void
test_dtwf_multi_locus_simulation(void)
{
    int ret;
    uint32_t n = 100;
    uint32_t m = 10;
    long seed = 10;
    double migration_matrix[] = { 0, 0.1, 0.1, 0 };
    const char *model_name;
    size_t num_ca_events, num_re_events;
    double t;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, m, .01, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_set(rng, seed);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(msp, 0, n, 0);
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
    msp_verify(msp, 0);
    num_ca_events = msp_get_num_common_ancestor_events(msp);
    num_re_events = msp_get_num_recombination_events(msp);
    CU_ASSERT_TRUE(num_ca_events > 0);
    CU_ASSERT_TRUE(num_re_events > 0);
    CU_ASSERT_EQUAL(ret, 0);
    msp_free(msp);

    /* Realloc the simulator under different memory params to see if
     * we get the same result. */
    gsl_rng_set(rng, seed);
    tsk_table_collection_clear(&tables);
    ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(msp, 0, n, 0);
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
    while ((ret = msp_run(msp, DBL_MAX, 1)) > 0) {
        msp_verify(msp, 0);
        CU_ASSERT_EQUAL_FATAL(msp->time, t);
        t++;
    }
    msp_verify(msp, 0);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(num_ca_events == msp_get_num_common_ancestor_events(msp));
    CU_ASSERT_TRUE(num_re_events == msp_get_num_recombination_events(msp));

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_gene_conversion_simulation(void)
{
    int ret;
    uint32_t n = 10;
    uint32_t m = 10;
    double track_lengths[] = { 1.0, 1.3333, 5, 10 };
    long seed = 10;
    size_t j, num_events, num_ca_events, num_re_events, num_gc_events;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, m, 0.1, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_set(rng, seed);

    for (j = 0; j < sizeof(track_lengths) / sizeof(double); j++) {
        memset(samples, 0, n * sizeof(sample_t));
        ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_gene_conversion_rate(msp, 1.0, track_lengths[j]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = msp_initialise(msp);
        CU_ASSERT_EQUAL(ret, 0);

        num_events = 0;
        while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
            msp_verify(msp, MSP_VERIFY_BREAKPOINTS);
            num_events++;
        }
        CU_ASSERT_EQUAL(ret, 0);
        msp_verify(msp, MSP_VERIFY_BREAKPOINTS);
        num_ca_events = msp_get_num_common_ancestor_events(msp);
        num_re_events = msp_get_num_recombination_events(msp);
        num_gc_events = msp_get_num_gene_conversion_events(msp);
        CU_ASSERT_TRUE(num_ca_events > 0);
        CU_ASSERT_TRUE(num_re_events > 0);
        CU_ASSERT_TRUE(num_gc_events > 0);
        msp_free(msp);

        /* Make sure we can build a tree sequence out of the result */
        ret = tsk_treeseq_init(&ts, &tables, TSK_BUILD_INDEXES);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        tsk_treeseq_free(&ts);

        tsk_table_collection_clear(&tables);
    }
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_likelihood_errors(void)
{
    int ret;
    double lik = 0;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 1.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 1.0, 2, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 1.0, 2, 1);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    /* Site has two mutations */
    ret = tsk_mutation_table_add_row(&tables.mutations, 0, 0, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 0, 1, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Ne <= 0 is an error */
    ret = msp_log_likelihood_arg(&ts, 1, 0, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_POPULATION_SIZE);
    ret = msp_log_likelihood_arg(&ts, -1, 0, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_POPULATION_SIZE);

    ret = msp_unnormalised_log_likelihood_mut(&ts, 0, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(lik, -DBL_MAX);

    ret = msp_unnormalised_log_likelihood_mut(&ts, 1, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_likelihood_zero_edges(void)
{
    int ret;
    double lik = 0;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_treeseq_init(&ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Zero edges should give a likelihood of zero */
    ret = msp_log_likelihood_arg(&ts, 1, 1, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(lik, 0);

    /* Zero mutations, so we get zero. */
    ret = msp_unnormalised_log_likelihood_mut(&ts, 0, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(lik, 0);
    ret = msp_unnormalised_log_likelihood_mut(&ts, 1, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(lik, 0);

    tsk_treeseq_free(&ts);

    /* Add in a site and a mutation */
    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 0, 0, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Zero edges should give a likelihood of zero */
    ret = msp_log_likelihood_arg(&ts, 1, 1, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(lik, 0);
    ret = msp_log_likelihood_arg(&ts, 0, 1, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(lik, 0);

    /* Zero mutation rate gives a likelihood of zero when there are mutations in
     * the ARG */
    ret = msp_unnormalised_log_likelihood_mut(&ts, 0, &lik);
    CU_ASSERT_EQUAL_FATAL(lik, -DBL_MAX);
    ret = msp_unnormalised_log_likelihood_mut(&ts, 1, &lik);
    CU_ASSERT_EQUAL_FATAL(lik, -DBL_MAX);

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_likelihood_three_leaves(void)
{
    int i;
    int ret;
    double rho[] = { 0.1, 1, 10 };
    double theta[] = { 0.1, 1, 10 };
    double ll_exact;
    double lik = 0;
    double tree_length = 19.0 / 8;
    double tol = 1e-9;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 3, 2);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 4, 2);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 1, 5, 1);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 5, 3);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.25, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 1, 6, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 6, 4);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.5, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 1, 7, 5);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 1, 7, 6);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_mutation_table_add_row(&tables.mutations, 0, 0, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 1, 1, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 2, 3, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 3, 0, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 4, 5, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.2, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.3, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.4, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.45, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    tsk_population_table_add_row(&tables.populations, NULL, 0);

    ret = tsk_treeseq_init(&ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (i = 0; i < 3; i++) {
        ll_exact = log(rho[i]) - (3 + 3 * rho[i]) * 0.1;
        ll_exact -= (6 + 3 * rho[i]) * 0.15;
        ll_exact -= (3 + 2.5 * rho[i]) * 0.25;
        ll_exact -= (1 + 2 * rho[i]) * 0.5;
        ret = msp_log_likelihood_arg(&ts, rho[i], 0.5, &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }
    for (i = 0; i < 3; i++) {
        ll_exact = (5 * log(tree_length * theta[i]) - tree_length * theta[i]);
        ll_exact -= 2 * log(4 * tree_length);
        ll_exact -= 2 * log(tree_length);
        ll_exact += log(3 / (4 * tree_length));
        ret = msp_unnormalised_log_likelihood_mut(&ts, theta[i], &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_likelihood_two_mrcas(void)
{
    int i;
    int ret;
    double rho[] = { 0.1, 1, 10 };
    double theta[] = { 0.1, 1, 10 };
    double ll_exact;
    double lik = 0;
    double tree_length = 1.5;
    double tol = 1e-9;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 2, 1);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 3, 1);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 4, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 5, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.15, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.15, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 6, 2);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 6, 4);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.5, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 7, 3);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 7, 5);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_mutation_table_add_row(&tables.mutations, 0, 1, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 1, 4, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 2, 3, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.2, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.7, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    tsk_population_table_add_row(&tables.populations, NULL, 0);

    ret = tsk_treeseq_init(&ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (i = 0; i < 3; i++) {
        ll_exact = log(rho[i]) - (1 + 2 * rho[i]) * 0.1;
        ll_exact += log(rho[i]) - (3 + 2 * rho[i]) * 0.05;
        ll_exact -= (6 + 2 * rho[i]) * 0.35;
        ll_exact -= (1 + rho[i]) * 0.5;
        ret = msp_log_likelihood_arg(&ts, rho[i], 0.5, &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }
    for (i = 0; i < 3; i++) {
        ll_exact = 3 * log(tree_length * theta[i]) - tree_length * theta[i];
        ll_exact -= log(tree_length);
        ll_exact -= 2 * log(2 * tree_length);
        ret = msp_unnormalised_log_likelihood_mut(&ts, theta[i], &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_likelihood_material_overhang(void)
{
    int i;
    int ret;
    double rho[] = { 0.1, 1, 10 };
    double theta[] = { 0.1, 1, 10 };
    double ll_exact;
    double lik = 0;
    double tree_length = 81.0 / 50;
    double tol = 1e-9;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 2, 1);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 3, 1);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.7, 4, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.7, 1, 5, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.15, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.15, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 6, 2);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.7, 6, 4);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.5, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 7, 3);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.7, 1, 7, 5);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 0.7, 8, 6);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 0.7, 8, 7);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 1.3, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_mutation_table_add_row(&tables.mutations, 0, 1, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 1, 4, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 2, 3, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.6, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.75, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    tsk_population_table_add_row(&tables.populations, NULL, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (i = 0; i < 3; i++) {
        ll_exact = log(rho[i]) - (1 + 2 * rho[i]) * 0.1;
        ll_exact += log(rho[i]) - (3 + 2 * rho[i]) * 0.05;
        ll_exact -= (6 + 2 * rho[i]) * 0.35;
        ll_exact -= (3 + rho[i]) * 0.5;
        ll_exact -= (1 + 0.4 * rho[i]) * 0.3;
        ret = msp_log_likelihood_arg(&ts, rho[i], 0.5, &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }
    for (i = 0; i < 3; i++) {
        ll_exact = 3 * log(tree_length * theta[i]) - tree_length * theta[i];
        ll_exact -= log(2 * tree_length);
        ll_exact += log(1.3 / tree_length);
        ll_exact -= log(tree_length);
        ret = msp_unnormalised_log_likelihood_mut(&ts, theta[i], &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_likelihood_material_gap(void)
{
    int i;
    int ret;
    double rho[] = { 0.1, 1, 10 };
    double theta[] = { 0.1, 1, 10 };
    double ll_exact;
    double lik = 0;
    double tree_length = 0.84;
    double tol = 1e-9;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 2, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 1, 3, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 0.5, 4, 3);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 5, 3);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.2, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.2, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 6, 2);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 6, 5);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.3, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 1, 7, 1);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 7, 6);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 7, 6);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.4, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 0.5, 8, 4);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 0.5, 8, 7);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.5, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_mutation_table_add_row(&tables.mutations, 0, 2, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 1, 3, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 2, 7, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 3, 6, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 4, 5, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.35, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.4, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.6, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.8, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    tsk_population_table_add_row(&tables.populations, NULL, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (i = 0; i < 3; i++) {
        ll_exact = log(rho[i]) - (1 + 2 * rho[i]) * 0.1;
        ll_exact += log(rho[i]) - (3 + 2 * rho[i]) * 0.1;
        ll_exact -= (6 + 2 * rho[i]) * 0.1;
        ll_exact -= (3 + 2.2 * rho[i]) * 0.1;
        ll_exact -= (1 + 0.4 * rho[i]) * 0.1;
        ret = msp_log_likelihood_arg(&ts, rho[i], 0.5, &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }
    for (i = 0; i < 3; i++) {
        ll_exact = 5 * log(tree_length * theta[i]) - tree_length * theta[i];
        ll_exact += 3 * log(0.4 / tree_length);
        ll_exact += 2 * log(0.5 / tree_length);
        ret = msp_unnormalised_log_likelihood_mut(&ts, theta[i], &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_likelihood_recombination_in_material_gap(void)
{
    int i;
    int ret;
    double rho[] = { 0.1, 1, 10 };
    double theta[] = { 0.1, 1, 10 };
    double ll_exact;
    double lik = 0;
    double tree_length = 1.34;
    double tol = 1e-9;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 2, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 1, 3, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 0.5, 4, 3);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 5, 3);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.2, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.2, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 6, 2);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 6, 5);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.3, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 7, 6);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 8, 6);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.4, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.4, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 0.5, 9, 4);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 9, 8);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.5, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 1, 10, 1);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 10, 7);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.6, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 1, 11, 9);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 1, 11, 10);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.7, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_mutation_table_add_row(&tables.mutations, 0, 2, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 1, 3, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 2, 6, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 3, 5, -1, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.35, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.6, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.8, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    tsk_population_table_add_row(&tables.populations, NULL, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (i = 0; i < 3; i++) {
        ll_exact = log(rho[i]) - (1 + 2 * rho[i]) * 0.1;
        ll_exact += log(rho[i]) - (3 + 2 * rho[i]) * 0.1;
        ll_exact -= (6 + 2 * rho[i]) * 0.1;
        ll_exact += log(0.2 * rho[i]) - (3 + 2.2 * rho[i]) * 0.1;
        ll_exact -= (6 + 2 * rho[i]) * 0.1;
        ll_exact -= (3 + 2 * rho[i]) * 0.1;
        ll_exact -= (1 + 1.4 * rho[i]) * 0.1;
        ret = msp_log_likelihood_arg(&ts, rho[i], 0.5, &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }
    for (i = 0; i < 3; i++) {
        ll_exact = 4 * log(tree_length * theta[i]) - tree_length * theta[i];
        ll_exact += log(0.6 / tree_length);
        ll_exact += 3 * log(0.7 / tree_length);
        ret = msp_unnormalised_log_likelihood_mut(&ts, theta[i], &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_multi_locus_simulation(void)
{
    int ret;
    uint32_t num_events = 0;
    uint32_t n = 100;
    uint32_t m = 100;
    long seed = 10;
    int model;
    double migration_matrix[] = { 0, 1, 1, 0 };
    size_t migration_events[4];
    int models[] = { MSP_MODEL_HUDSON, MSP_MODEL_SMC, MSP_MODEL_SMC_PRIME };
    const char *model_names[] = { "hudson", "smc", "smc_prime" };
    const char *model_name;
    bool store_full_arg[] = { true, false };
    size_t j, k;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    ret = recomb_map_alloc_uniform(&recomb_map, m, 1 / m, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < sizeof(models) / sizeof(int); j++) {
        for (k = 0; k < sizeof(store_full_arg) / sizeof(bool); k++) {
            sample_t *samples = malloc(n * sizeof(sample_t));
            msp_t *msp = malloc(sizeof(msp_t));
            gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

            CU_ASSERT_FATAL(msp != NULL);
            CU_ASSERT_FATAL(samples != NULL);
            CU_ASSERT_FATAL(rng != NULL);
            gsl_rng_set(rng, seed);
            tsk_table_collection_clear(&tables);
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
                    ret = msp_set_simulation_model_hudson(msp);
                    break;
                case 1:
                    ret = msp_set_simulation_model_smc(msp);
                    break;
                case 2:
                    ret = msp_set_simulation_model_smc_prime(msp);
                    break;
            }
            ret = msp_set_store_full_arg(msp, store_full_arg[k]);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            ret = msp_initialise(msp);
            CU_ASSERT_EQUAL(ret, 0);

            num_events = 0;
            while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
                msp_verify(msp, MSP_VERIFY_BREAKPOINTS);
                num_events++;
            }
            CU_ASSERT_EQUAL(ret, 0);
            msp_verify(msp, MSP_VERIFY_BREAKPOINTS);
            ret = msp_get_num_migration_events(msp, migration_events);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT(num_events > n - 1);
            CU_ASSERT_EQUAL(1 + num_events,
                /* diagonals must be zero here */
                migration_events[1] + migration_events[2]
                    + msp_get_num_recombination_events(msp)
                    + msp_get_num_common_ancestor_events(msp)
                    + msp_get_num_rejected_common_ancestor_events(msp));
            if (models[j] == MSP_MODEL_HUDSON) {
                CU_ASSERT_EQUAL(msp_get_num_rejected_common_ancestor_events(msp), 0);
            }

            gsl_rng_set(rng, seed);
            ret = msp_reset(msp);
            num_events = 0;
            while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
                msp_verify(msp, MSP_VERIFY_BREAKPOINTS);
                num_events++;
            }
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_get_num_migration_events(msp, migration_events);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(
                1 + num_events, migration_events[1] + migration_events[2]
                                    + msp_get_num_recombination_events(msp)
                                    + msp_get_num_common_ancestor_events(msp)
                                    + msp_get_num_rejected_common_ancestor_events(msp));
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
    }
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
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
    double migration_matrix[] = { 0, 1, 1, 0 };
    size_t j;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    msp_t msp;
    tsk_treeseq_t ts;
    mutgen_t mutgen;
    mutation_model_t mut_model;
    tsk_table_collection_t tables;
    recomb_map_t recomb_map;
    interval_map_t mut_map;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    ret = recomb_map_alloc_uniform(&recomb_map, m / 2, 1.0, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = interval_map_alloc_single(&mut_map, m / 2, mutation_rate);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_BINARY);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_init(&tables, 0);
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
    ret = mutgen_alloc(&mutgen, rng, &mut_map, &mut_model, 3);
    CU_ASSERT_EQUAL(ret, 0);

    mutgen_print_state(&mutgen, _devnull);

    for (j = 0; j < num_replicates; j++) {
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_run(&msp, DBL_MAX, SIZE_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        msp_verify(&msp, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = msp_finalise_tables(&msp);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = mutgen_generate(&mutgen, &tables, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        tables.sequence_length = m;
        ret = tsk_treeseq_init(&ts, &tables, TSK_BUILD_INDEXES);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        verify_simulator_tsk_treeseq_equality(&msp, &ts, 1.0);
        tsk_treeseq_print_state(&ts, _devnull);
        mutgen_print_state(&mutgen, _devnull);
        ret = msp_reset(&msp);
        CU_ASSERT_EQUAL_FATAL(msp_get_num_edges(&msp), 0);
        CU_ASSERT_EQUAL_FATAL(msp_get_num_migrations(&msp), 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        tsk_treeseq_free(&ts);
    }
    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    mutgen_free(&mutgen);
    mutation_model_free(&mut_model);
    gsl_rng_free(rng);
    free(samples);
    tsk_table_collection_free(&tables);
    recomb_map_free(&recomb_map);
    interval_map_free(&mut_map);
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
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    ret = recomb_map_alloc_uniform(&recomb_map, m, 1 / m, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
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
    msp_verify(msp, 0);

    msp_reset(msp);
    while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
        msp_verify(msp, 0);
        if (msp->time == 0.1) {
            t1_found = 1;
        }
        CU_ASSERT_EQUAL(msp_get_time(msp), msp->time);
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
    tsk_table_collection_free(&tables);
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
    uint32_t n = 1000;
    uint32_t m = 10;
    long seed = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    uint32_t num_bottlenecks = 10;
    bottleneck_desc_t bottlenecks[num_bottlenecks];
    double t;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    ret = recomb_map_alloc_uniform(&recomb_map, m, m, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
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
        ret = msp_add_simple_bottleneck(
            msp, bottlenecks[j].time, 0, bottlenecks[j].parameter);
        CU_ASSERT_EQUAL(ret, 0);
    }
    ret = msp_set_population_configuration(msp, 0, 0.25, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    for (j = 0; j < num_bottlenecks - 1; j++) {
        ret = msp_run(msp, bottlenecks[j].time + 1e-6, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_TIME);
        CU_ASSERT_FALSE(msp_is_completed(msp));
        CU_ASSERT_EQUAL(msp->time, bottlenecks[j].time + 1e-6);
        msp_verify(msp, 0);
    }
    ret = msp_run(msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(msp_is_completed(msp));
    CU_ASSERT_EQUAL(msp->time, bottlenecks[num_bottlenecks - 1].time);
    msp_verify(msp, 0);

    /* Test out resets on partially completed simulations. */
    ret = msp_reset(msp);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_bottlenecks - 1; j++) {
        ret = msp_run(msp, bottlenecks[j].time, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 2);
    }
    ret = msp_reset(msp);
    msp_verify(msp, 0);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_dirac_coalescent_bad_parameters(void)
{
    int j;
    int ret;
    msp_t msp;
    unsigned int n = 10;
    double cs[] = { -1 };
    double psis[] = { -1e6, 1e6 };
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < sizeof(cs) / sizeof(*cs); j++) {
        ret = msp_set_simulation_model_dirac(&msp, 0.1, cs[j]);
        CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_C);
    }
    for (j = 0; j < sizeof(psis) / sizeof(*psis); j++) {
        ret = msp_set_simulation_model_dirac(&msp, psis[j], 10);
        CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PSI);
    }

    msp_free(&msp);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
    free(samples);
    gsl_rng_free(rng);
}

static void
test_beta_coalescent_bad_parameters(void)
{
    int j;
    int ret;
    msp_t msp;
    unsigned int n = 10;
    double alphas[] = { -1e6, 0, 0.001, 1.0, 2.0, 1e6 };
    double truncation_points[] = { -1e6, 0, 1.001, 1e6 };
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < sizeof(alphas) / sizeof(*alphas); j++) {
        ret = msp_set_simulation_model_beta(&msp, alphas[j], 1);
        CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_BETA_MODEL_ALPHA);
    }
    for (j = 0; j < sizeof(truncation_points) / sizeof(*truncation_points); j++) {
        ret = msp_set_simulation_model_beta(&msp, 1.5, truncation_points[j]);
        CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_TRUNCATION_POINT);
    }

    msp_free(&msp);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
    free(samples);
    gsl_rng_free(rng);
}

static void
test_multiple_mergers_simulation(void)
{
    int ret;
    size_t j, k, o, p;
    uint32_t n = 10;
    uint32_t m = 10;
    long seed = 10;
    bool store_full_arg[] = { true, false };
    /* These simulations can be slow, so just choose a few param combinations */
    double beta_params[][2] = { { 1.1, 0.5 }, { 1.99, 1 } };
    /* TODO what are good psi parameters here? */
    double psi_params[][2] = { { 0.9, 10 }, { 0.1, 1 } };
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, m, m, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < 2; j++) {
        if (j == 0) {
            o = sizeof(psi_params) / sizeof(*psi_params);
        } else if (j == 1) {
            o = sizeof(beta_params) / sizeof(*beta_params);
        }
        for (k = 0; k < sizeof(store_full_arg) / sizeof(bool); k++) {
            for (p = 0; p < o; p++) {
                gsl_rng_set(rng, seed);
                /* TODO check non-zero sample times here to make sure they fail. */
                memset(samples, 0, n * sizeof(sample_t));
                tsk_table_collection_clear(&tables);
                ret = msp_alloc(msp, n, samples, &recomb_map, &tables, rng);
                CU_ASSERT_EQUAL(ret, 0);
                if (j == 0) {
                    ret = msp_set_simulation_model_dirac(
                        msp, psi_params[p][0], psi_params[p][1]);
                } else {
                    ret = msp_set_simulation_model_beta(
                        msp, beta_params[p][0], beta_params[p][1]);
                }
                CU_ASSERT_EQUAL(ret, 0);
                /* TODO check for adding various complications like multiple
                 * populations etc to ensure they fail.
                 */
                ret = msp_set_store_full_arg(msp, store_full_arg[k]);
                CU_ASSERT_EQUAL(ret, 0);
                ret = msp_initialise(msp);
                CU_ASSERT_EQUAL(ret, 0);
                msp_print_state(msp, _devnull);

                ret = msp_run(msp, DBL_MAX, ULONG_MAX);
                CU_ASSERT_EQUAL_FATAL(ret, 0);
                CU_ASSERT_TRUE(msp_is_completed(msp));
                CU_ASSERT_TRUE(msp->time > 0);
                msp_verify(msp, 0);

                msp_reset(msp);
                while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
                    msp_verify(msp, 0);
                }
                CU_ASSERT_EQUAL_FATAL(ret, 0);
                CU_ASSERT_TRUE(msp_is_completed(msp));

                ret = msp_free(msp);
                CU_ASSERT_EQUAL(ret, 0);
            }
        }
    }
    gsl_rng_free(rng);
    free(msp);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_multiple_mergers_growth_rate(void)
{
    int ret;
    size_t n = 10;
    size_t m = 1;
    int j;
    sample_t *samples = calloc(n, sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    msp_t msp;
    recomb_map_t recomb_map;
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, m, m, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < 2; j++) {

        ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
        CU_ASSERT_EQUAL(ret, 0);
        if (j == 0) {
            ret = msp_set_simulation_model_dirac(&msp, 0.9, 100);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            ret = msp_set_simulation_model_beta(&msp, 1.9, 1);
            CU_ASSERT_EQUAL(ret, 0);
        }

        /* Set to a nonzero growth_rate */
        ret = msp_set_population_configuration(&msp, 0, 1, -0.01);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = msp_initialise(&msp);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 0);

        msp_print_state(&msp, _devnull);

        msp_free(&msp);

        tsk_table_collection_clear(&tables);
    }
    gsl_rng_free(rng);
    free(samples);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
}

static void
test_simple_recomb_map(void)
{
    int ret;
    recomb_map_t recomb_map;
    double seq_length;
    double positions[][2] = { { 0.0, 1.0 }, { 0.0, 100.0 }, { 0.0, 10000.0 } };
    double rates[] = { 0.0, 1.0 };
    size_t j;

    for (j = 0; j < 3; j++) {
        seq_length = positions[j][1];
        ret = recomb_map_alloc(&recomb_map, 2, positions[j], rates, true);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        recomb_map_print_state(&recomb_map, _devnull);
        CU_ASSERT_EQUAL(recomb_map_get_size(&recomb_map), 2);
        CU_ASSERT_EQUAL(recomb_map_get_sequence_length(&recomb_map), seq_length);
        ret = recomb_map_free(&recomb_map);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
}

static void
verify_recomb_maps_equal(recomb_map_t *actual, recomb_map_t *expected)
{
    size_t i;

    CU_ASSERT_EQUAL(actual->map.size, expected->map.size);
    CU_ASSERT_EQUAL(actual->discrete, expected->discrete);
    for (i = 0; i < actual->map.size; i++) {
        CU_ASSERT_DOUBLE_EQUAL(actual->map.position[i], expected->map.position[i], 0.0);
        CU_ASSERT_DOUBLE_EQUAL(actual->map.value[i], expected->map.value[i], 0.0);
        CU_ASSERT_DOUBLE_EQUAL(actual->cumulative[i], expected->cumulative[i], 1e-6);
    }
}

static double
double_rate(void *_, double rate)
{
    return 2 * rate;
}

static void
test_recomb_map_copy(void)
{
    int ret = 0;
    recomb_map_t map, other, copy;
    double positions[] = { 0.0, 1.0, 2.0 };
    double rates[] = { 1.0, 2.0, 0.0 };
    double other_rates[] = { 2.0, 4.0, 0.0 };

    ret = recomb_map_alloc(&map, 3, positions, rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = recomb_map_copy(&copy, &map);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_recomb_maps_equal(&copy, &map);
    recomb_map_free(&map);
    recomb_map_free(&copy);

    ret = recomb_map_alloc(&map, 3, positions, rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = recomb_map_alloc(&other, 3, positions, other_rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    recomb_map_convert_rates(&map, double_rate, NULL);
    verify_recomb_maps_equal(&map, &other);
    recomb_map_free(&map);
    recomb_map_free(&other);
}

static void
test_recomb_map_errors(void)
{
    int ret;
    recomb_map_t recomb_map;
    double positions[] = { 0.0, 1.0, 2.0 };
    double rates[] = { 1.0, 2.0, 0.0 };
    double short_positions[] = { 0.0, 0.25, 0.5 };

    ret = recomb_map_alloc(&recomb_map, 0, positions, rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INSUFFICIENT_INTERVALS);
    recomb_map_free(&recomb_map);

    ret = recomb_map_alloc(&recomb_map, 3, short_positions, rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);

    ret = recomb_map_alloc(&recomb_map, 1, positions, rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INSUFFICIENT_INTERVALS);
    recomb_map_free(&recomb_map);

    positions[0] = 1;
    ret = recomb_map_alloc(&recomb_map, 2, positions, rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INTERVAL_MAP_START_NON_ZERO);
    recomb_map_free(&recomb_map);
    positions[0] = 0;

    positions[1] = 3.0;
    ret = recomb_map_alloc(&recomb_map, 3, positions, rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INTERVAL_POSITIONS_UNSORTED);
    recomb_map_free(&recomb_map);
    positions[1] = 1.0;

    positions[0] = -1;
    ret = recomb_map_alloc(&recomb_map, 3, positions, rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INTERVAL_MAP_START_NON_ZERO);
    recomb_map_free(&recomb_map);
    positions[0] = 0.0;

    rates[0] = -1;
    ret = recomb_map_alloc(&recomb_map, 3, positions, rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);
    rates[0] = 1.0;

    rates[0] = INFINITY;
    ret = recomb_map_alloc(&recomb_map, 3, positions, rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);
    rates[0] = 1.0;

    rates[0] = NAN;
    ret = recomb_map_alloc(&recomb_map, 3, positions, rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);
    rates[0] = 1.0;

    ret = recomb_map_alloc(&recomb_map, 3, positions, rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    recomb_map_free(&recomb_map);

    ret = recomb_map_alloc(&recomb_map, 3, short_positions, rates, false);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    recomb_map_free(&recomb_map);
}

static void
verify_recomb_map(double length, double *positions, double *rates, size_t size)
{

    int ret;
    recomb_map_t recomb_map;
    size_t j;
    double *ret_rates, *ret_positions;

    ret_rates = malloc(size * sizeof(double));
    ret_positions = malloc(size * sizeof(double));

    CU_ASSERT_FATAL(ret_rates != NULL);
    CU_ASSERT_FATAL(ret_positions != NULL);

    ret = recomb_map_alloc(&recomb_map, size, positions, rates, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    recomb_map_print_state(&recomb_map, _devnull);
    CU_ASSERT_EQUAL(recomb_map_get_size(&recomb_map), size);

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
    double p1[] = { 0, 0.05019838393314813, 0.36933662489552865, 1 };
    double r1[] = { 3.5510784169955434, 4.184964179610539, 3.800808140657212, 0 };
    double p2[] = { 0, 0.125, 0.875, 1, 4, 8, 16 };
    double r2[] = { 0.1, 6.0, 3.333, 2.1, 0.0, 2.2, 0 };

    verify_recomb_map(1.0, p1, r1, 4);
    verify_recomb_map(1.0, p1, r1, 4);
    verify_recomb_map(1.0, p1, r1, 4);

    verify_recomb_map(16.0, p2, r2, 7);
    verify_recomb_map(16.0, p2, r2, 7);
    verify_recomb_map(16.0, p2, r2, 7);
}

static void
test_translate_position_and_recomb_mass(void)
{
    recomb_map_t map;
    double p1[] = { 0, 6, 13, 20 };
    double r1[] = { 3, 0, 1, 0 };
    recomb_map_alloc(&map, 4, p1, r1, true);

    /* interval edges */
    CU_ASSERT_EQUAL(recomb_map_position_to_mass(&map, 0), 0);
    CU_ASSERT_EQUAL(recomb_map_position_to_mass(&map, 6), 18);
    CU_ASSERT_EQUAL(recomb_map_position_to_mass(&map, 13), 18);
    CU_ASSERT_EQUAL(recomb_map_position_to_mass(&map, 19), 24);
    CU_ASSERT_EQUAL(recomb_map_mass_to_position(&map, 0), 0);
    CU_ASSERT_EQUAL(recomb_map_mass_to_position(&map, 18), 6);
    CU_ASSERT_EQUAL(recomb_map_mass_to_position(&map, 24), 19);

    /* intervals with recombination */
    CU_ASSERT_EQUAL(recomb_map_position_to_mass(&map, 4), 12);
    CU_ASSERT_EQUAL(recomb_map_position_to_mass(&map, 14), 19);
    CU_ASSERT_EQUAL(recomb_map_position_to_mass(&map, 16), 21);
    CU_ASSERT_EQUAL(recomb_map_mass_to_position(&map, 12), 4);
    CU_ASSERT_EQUAL(recomb_map_mass_to_position(&map, 19), 14);
    CU_ASSERT_EQUAL(recomb_map_mass_to_position(&map, 21), 16);

    /* inside recombination interval */
    CU_ASSERT_EQUAL(recomb_map_position_to_mass(&map, 8), 18);

    recomb_map_free(&map);
}

static void
test_recomb_map_mass_between(void)
{
    recomb_map_t discrete_map;
    recomb_map_t cont_map;
    double p1[] = { 0, 6, 13, 20 };
    double r1[] = { 3, 0, 1, 0 };
    double tol = 1e-9;

    recomb_map_alloc(&discrete_map, 4, p1, r1, true);
    recomb_map_alloc(&cont_map, 4, p1, r1, false);

    CU_ASSERT_DOUBLE_EQUAL_FATAL(recomb_map_mass_between(&discrete_map, 0, 2), 6, tol);
    CU_ASSERT_DOUBLE_EQUAL_FATAL(recomb_map_mass_between(&cont_map, 0, 2), 6, tol);
    CU_ASSERT_DOUBLE_EQUAL_FATAL(
        recomb_map_mass_between_left_exclusive(&discrete_map, 0, 2), 3, tol);
    CU_ASSERT_DOUBLE_EQUAL_FATAL(
        recomb_map_mass_between_left_exclusive(&cont_map, 0, 2), 6, tol);

    recomb_map_free(&discrete_map);
    recomb_map_free(&cont_map);
}

static void
test_msp_binary_interval_search(void)
{
    double values[] = { -10, 10, 20, 30 };
    size_t size = 4;
    size_t idx;

    // Search from bottom
    idx = msp_binary_interval_search(-11, values, size);
    CU_ASSERT_EQUAL(idx, 0);
    // Exact match returns index of value
    idx = msp_binary_interval_search(-10, values, size);
    CU_ASSERT_EQUAL(idx, 0);
    // values[index-1] < query <= values[index]
    idx = msp_binary_interval_search(9, values, size);
    CU_ASSERT_EQUAL(idx, 1);
    // exact match
    idx = msp_binary_interval_search(10, values, size);
    CU_ASSERT_EQUAL(idx, 1);
    // Within mid interval
    idx = msp_binary_interval_search(11, values, size);
    CU_ASSERT_EQUAL(idx, 2);
    idx = msp_binary_interval_search(19, values, size);
    CU_ASSERT_EQUAL(idx, 2);
    // Exact
    idx = msp_binary_interval_search(20, values, size);
    CU_ASSERT_EQUAL(idx, 2);
    // Within
    idx = msp_binary_interval_search(21, values, size);
    CU_ASSERT_EQUAL(idx, 3);
    // Exact
    idx = msp_binary_interval_search(30, values, size);
    CU_ASSERT_EQUAL(idx, 3);
    // from the top - return last element
    idx = msp_binary_interval_search(31, values, size);
    CU_ASSERT_EQUAL(idx, 3);
    // way above - same
    idx = msp_binary_interval_search(300, values, size);
    CU_ASSERT_EQUAL(idx, 3);
}

static void
test_msp_binary_interval_search_repeating(void)
{
    double values_repeating[] = { 0, 10, 10, 30 };
    size_t size = 4;
    size_t idx;

    // Same as above
    idx = msp_binary_interval_search(-1, values_repeating, size);
    CU_ASSERT_EQUAL(idx, 0);
    // Want leftmost interval
    idx = msp_binary_interval_search(10, values_repeating, size);
    CU_ASSERT_EQUAL(idx, 1);
    // Same as above
    idx = msp_binary_interval_search(11, values_repeating, size);
    CU_ASSERT_EQUAL(idx, 3);
}

static void
test_msp_binary_interval_search_edge_cases(void)
{
    double values_empty[] = {};
    size_t idx;

    // Empty list
    idx = msp_binary_interval_search(0, values_empty, 0);
    CU_ASSERT_EQUAL(idx, 0);

    // Size 1 list
    double values_one[] = { 10 };

    // below
    idx = msp_binary_interval_search(9, values_one, 1);
    CU_ASSERT_EQUAL(idx, 0);
    // exact
    idx = msp_binary_interval_search(10, values_one, 1);
    CU_ASSERT_EQUAL(idx, 0);
    // above
    idx = msp_binary_interval_search(11, values_one, 1);
    CU_ASSERT_EQUAL(idx, 0);

    // Size 2 list
    double values_two[] = { 10, 20 };
    idx = msp_binary_interval_search(9, values_two, 2);
    CU_ASSERT_EQUAL(idx, 0);
    idx = msp_binary_interval_search(10, values_two, 2);
    CU_ASSERT_EQUAL(idx, 0);
    idx = msp_binary_interval_search(15, values_two, 2);
    CU_ASSERT_EQUAL(idx, 1);
    idx = msp_binary_interval_search(20, values_two, 2);
    CU_ASSERT_EQUAL(idx, 1);
    idx = msp_binary_interval_search(21, values_two, 2);
    CU_ASSERT_EQUAL(idx, 1);

    // All zeros
    double values_zeros[] = { 0, 0, 0 };
    idx = msp_binary_interval_search(-1, values_zeros, 3);
    CU_ASSERT_EQUAL(idx, 0);
    idx = msp_binary_interval_search(0, values_zeros, 3);
    CU_ASSERT_EQUAL(idx, 0);
    idx = msp_binary_interval_search(1, values_zeros, 3);
    CU_ASSERT_EQUAL(idx, 2);
}

static void
verify_simulate_from(int model, recomb_map_t *recomb_map,
    tsk_table_collection_t *from_tables, size_t num_replicates)
{
    int ret;
    size_t j;
    tsk_bookmark_t pos;
    tsk_table_collection_t tables;
    tsk_treeseq_t final;
    tsk_tree_t tree;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    ret = tsk_table_collection_copy(from_tables, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_alloc(&msp, 0, NULL, recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    if (model == MSP_MODEL_DTWF) {
        ret = msp_set_simulation_model_dtwf(&msp);
        CU_ASSERT_EQUAL(ret, 0);
    }
    ret = msp_set_population_configuration(&msp, 0, 10, 0);
    CU_ASSERT_EQUAL(ret, 0);
    /* TODO add dirac and other models */
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < num_replicates; j++) {
        msp_verify(&msp, 0);
        ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_TRUE(msp_is_completed(&msp));
        ret = msp_finalise_tables(&msp);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_treeseq_init(&final, &tables, TSK_BUILD_INDEXES);
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        ret = tsk_tree_init(&tree, &final, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        for (ret = tsk_tree_first(&tree); ret == 1; ret = tsk_tree_next(&tree)) {
            CU_ASSERT_EQUAL_FATAL(tsk_tree_get_num_roots(&tree), 1);
        }
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        ret = tsk_table_collection_record_num_rows(from_tables, &pos);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_table_collection_truncate(&tables, &pos);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_TRUE(tsk_table_collection_equals(from_tables, &tables));

        tsk_treeseq_free(&final);

        tsk_tree_free(&tree);
        ret = msp_reset(&msp);
        /* printf("ret = %s\n", msp_strerror(ret)); */
        CU_ASSERT_EQUAL(ret, 0);
    }
    msp_free(&msp);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

/* Verify that the initial state we get in a new simulator from calling
 * with from_ts is equivalent to the state in the original simulator */
static void
verify_initial_simulate_from_state(
    msp_t *msp_source, recomb_map_t *recomb_map, tsk_table_collection_t *from_tables)
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
verify_simple_simulate_from(int model, uint32_t n, size_t num_loci,
    double sequence_length, double recombination_rate, size_t num_events,
    size_t num_replicates)
{
    int ret;
    tsk_table_collection_t tables;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t msp;
    recomb_map_t recomb_map;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(
        &recomb_map, sequence_length, recombination_rate, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
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
    ret = tsk_provenance_table_add_row(&tables.provenances, "time", 4, "record", 6);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_initial_simulate_from_state(&msp, &recomb_map, &tables);
    verify_simulate_from(model, &recomb_map, &tables, num_replicates);

    msp_free(&msp);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
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
test_simulate_from_single_locus_replicates(void)
{
    verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 1, 1.0, 0, 5, 10);
    verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 1, 1.0, 0, 5, 10);
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
    tsk_table_collection_t tables;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    double recombination_rate = 2;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, recombination_rate, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
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
    tsk_table_collection_free(&tables);
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
    tsk_table_collection_t from_tables;

    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 10.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&from_tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    from_tables.sequence_length = 1.0;
    /* Add a node so that we get past the first check */
    ret = tsk_population_table_add_row(&from_tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_add_row(&from_tables.nodes, 0, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Malformed tree sequence */
    from_tables.nodes.individual[0] = 100;
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_FATAL(msp_is_tsk_error(ret));
    CU_ASSERT_EQUAL_FATAL(
        ret ^ (1 << MSP_TSK_ERR_BIT), TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    CU_ASSERT_STRING_EQUAL(
        msp_strerror(ret), tsk_strerror(TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS));
    from_tables.nodes.individual[0] = -1;
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
    ret = tsk_node_table_add_row(
        &from_tables.nodes, TSK_NODE_IS_SAMPLE, 1.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &from_tables.nodes, TSK_NODE_IS_SAMPLE, 2.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_start_time(&msp, 1.999);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_START_TIME_FROM_TS);
    msp_free(&msp);

    /* Num populations should be 1 */
    tsk_population_table_add_row(&from_tables.populations, NULL, 0);
    ret = msp_set_start_time(&msp, 0);
    CU_ASSERT_EQUAL(ret, 0);
    from_tables.sequence_length = 10.0;
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, MSP_ERR_INCOMPATIBLE_FROM_TS);
    msp_free(&msp);

    tsk_population_table_clear(&from_tables.populations);
    tsk_population_table_add_row(&from_tables.populations, NULL, 0);
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_start_time(&msp, 1.999);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_START_TIME_FROM_TS);
    msp_free(&msp);

    /* Must have legitimate population references */
    from_tables.nodes.population[0] = -1;
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_start_time(&msp, 2.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    msp_free(&msp);

    /* Check to make sure we can run this correctly */
    from_tables.nodes.population[0] = 0;
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_start_time(&msp, 2.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    msp_free(&msp);

    /* Make a tree sequence that we cannot recover trees from. This only happens
     * at initialisation time. */
    ret = tsk_edge_table_add_row(&from_tables.edges, 0, 1, 1, 0);
    CU_ASSERT_FATAL(ret >= 0);
    tsk_edge_table_add_row(&from_tables.edges, 0, 1, 2, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = msp_alloc(&msp, 0, NULL, &recomb_map, &from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_FATAL(msp_is_tsk_error(ret));
    CU_ASSERT_EQUAL_FATAL(
        ret ^ (1 << MSP_TSK_ERR_BIT), TSK_ERR_BAD_EDGES_CONTRADICTORY_CHILDREN);
    CU_ASSERT_STRING_EQUAL(
        msp_strerror(ret), tsk_strerror(TSK_ERR_BAD_EDGES_CONTRADICTORY_CHILDREN));
    msp_free(&msp);

    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&from_tables);
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
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1.0;

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(&msp, 0, samples, &recomb_map, NULL, rng);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);
    ret = msp_alloc(&msp, n, samples, &recomb_map, NULL, rng);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);
    ret = msp_alloc(&msp, n, NULL, &recomb_map, NULL, rng);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);
    ret = msp_alloc(&msp, n, NULL, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);

    tsk_table_collection_clear(&tables);
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_set_start_time(&msp, -1);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_START_TIME);

    ret = msp_set_gene_conversion_rate(&msp, -1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = msp_set_gene_conversion_rate(&msp, 1, -1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = msp_set_gene_conversion_rate(&msp, 1, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    msp_free(&msp);
    gsl_rng_free(rng);
    recomb_map_free(&recomb_map);
    free(samples);
    tsk_table_collection_free(&tables);
}

static void
insert_single_tree(tsk_table_collection_t *tables, int alphabet)
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
    tables->sequence_length = 1.0;
    ret = tsk_node_table_add_row(
        &tables->nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_add_row(
        &tables->nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables->nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables->nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables->nodes, 0, 1.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables->nodes, 0, 2.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables->nodes, 0, 3.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables->edges, 0, 1, 4, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables->edges, 0, 1, 4, 1);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables->edges, 0, 1, 5, 2);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables->edges, 0, 1, 5, 3);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables->edges, 0, 1, 6, 4);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables->edges, 0, 1, 6, 5);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_population_table_add_row(&tables->populations, NULL, 0);
    CU_ASSERT_FATAL(ret == 0);

    /* Add a site and a mutation */
    if (alphabet == ALPHABET_BINARY) {
        ret = tsk_site_table_add_row(&tables->sites, 0.1, "0", 1, NULL, 0);
        CU_ASSERT_FATAL(ret >= 0);
        ret = tsk_mutation_table_add_row(&tables->mutations, 0, 0, -1, "1", 1, NULL, 0);
        CU_ASSERT_FATAL(ret >= 0);
    } else if (alphabet == ALPHABET_NUCLEOTIDE) {
        ret = tsk_site_table_add_row(&tables->sites, 0.1, "A", 1, NULL, 0);
        CU_ASSERT_FATAL(ret >= 0);
        ret = tsk_mutation_table_add_row(&tables->mutations, 0, 0, -1, "C", 1, NULL, 0);
        CU_ASSERT_FATAL(ret >= 0);
    }

    ret = tsk_table_collection_check_integrity(tables, TSK_CHECK_OFFSETS);

    CU_ASSERT_FATAL(ret == 0);
}

static void
test_mutgen_simple_map(void)
{
    int ret = 0;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    mutation_model_t mut_model;
    tsk_table_collection_t tables;
    interval_map_t rate_map;
    double pos[] = { 0, 0.1, 0.2, 0.3, 0.4, 1.0 };
    double rate[] = { 0, 0.01, 0.02, 0.03, 0.04, 0.0 };

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);

    ret = interval_map_alloc(&rate_map, 6, pos, rate);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map->size, 6);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map->position[0], 0);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map->position[1], 0.1);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map->position[2], 0.2);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map->position[3], 0.3);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map->position[4], 0.4);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map->position[5], 1.0);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map->value[0], 0);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map->value[1], 0.01);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map->value[2], 0.02);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map->value[3], 0.03);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map->value[4], 0.04);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map->value[5], 0);

    ret = mutgen_generate(&mutgen, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 2.0;
    ret = mutgen_generate(&mutgen, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INCOMPATIBLE_MUTATION_MAP);

    mutgen_free(&mutgen);
    interval_map_free(&rate_map);
    mutation_model_free(&mut_model);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_mutgen_errors(void)
{
    int ret = 0;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    interval_map_t rate_map;
    mutation_model_t mut_model, mut_model_binary;

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);

    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = matrix_mutation_model_factory(&mut_model_binary, ALPHABET_BINARY);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = interval_map_alloc_single(&rate_map, 1, -1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_MUTATION_MAP_RATE);
    interval_map_free(&rate_map);
    mutgen_free(&mutgen);

    ret = interval_map_alloc_single(&rate_map, 5, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INCOMPATIBLE_MUTATION_MAP);

    tables.sequence_length = 0.1;
    ret = mutgen_generate(&mutgen, &tables, 0);
    CU_ASSERT_FATAL(msp_is_tsk_error(ret));
    tables.sequence_length = 1.0;
    interval_map_free(&rate_map);
    mutgen_free(&mutgen);

    /* mix of binary and nucleotide alleles */
    ret = interval_map_alloc_single(&rate_map, 1, 20);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model_binary, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_time_interval(&mutgen, 0.5, 10.0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* we shouldn't error the first time since existing site is nucleotide
     * but not at an integer loction */
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    mutgen_free(&mutgen);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_time_interval(&mutgen, 0.0, 0.5);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_UNKNOWN_ALLELE);

    mutgen_free(&mutgen);
    interval_map_free(&rate_map);
    mutation_model_free(&mut_model);
    mutation_model_free(&mut_model_binary);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_mutgen_bad_mutation_order(void)
{
    int ret = 0;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    interval_map_t rate_map;
    mutation_model_t mut_model;

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);
    ret = interval_map_alloc_single(&rate_map, 1, 20);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* bad mutation generation order */
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_time_interval(&mutgen, 0.0, 0.5);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_time_interval(&mutgen, 0.5, 1.0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_MUTATION_GENERATION_OUT_OF_ORDER);

    mutgen_free(&mutgen);
    interval_map_free(&rate_map);
    mutation_model_free(&mut_model);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_single_tree_mutgen(void)
{
    int ret = 0;
    size_t j;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables1, tables2;
    mutation_model_t mut_model;
    interval_map_t rate_map;

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables1, ALPHABET_BINARY);
    insert_single_tree(&tables2, ALPHABET_BINARY);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_BINARY);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = interval_map_alloc_single(&rate_map, 1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables1.mutations.num_rows == 0);
    mutgen_print_state(&mutgen, _devnull);
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    gsl_rng_set(rng, 1);
    /* Cheat a bit to avoid reallocating a new rate_map */
    rate_map.value[0] = 10;
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    mutgen_print_state(&mutgen, _devnull);
    CU_ASSERT_TRUE(tables1.mutations.num_rows > 0);
    CU_ASSERT_TRUE(tables1.mutations.num_rows == tables1.sites.num_rows);
    for (j = 0; j < tables1.mutations.num_rows; j++) {
        CU_ASSERT_TRUE(tables1.mutations.site[j] == j);
        CU_ASSERT_TRUE(tables1.sites.position[j] <= 1.0);
        CU_ASSERT_TRUE(tables1.mutations.node[j] < 6);
        CU_ASSERT_EQUAL(tables1.sites.ancestral_state[j], '0');
        CU_ASSERT_EQUAL(tables1.mutations.derived_state[j], '1');
    }
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Test the reallocing behavior by setting a very small
     * block size.
     */
    gsl_rng_set(rng, 1);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables1, &tables2));
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_table_collection_free(&tables1);
    tsk_table_collection_free(&tables2);
    mutation_model_free(&mut_model);
    interval_map_free(&rate_map);
    gsl_rng_free(rng);
}

static void
test_single_tree_mutgen_keep_sites(void)
{
    int ret = 0;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    tsk_table_collection_t copy, copy2;
    interval_map_t rate_map;
    mutation_model_t mut_model;
    mutgen_t mutgen;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);
    ret = tsk_table_collection_init(&copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&copy, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy));
    ret = tsk_table_collection_init(&copy2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&copy2, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy2));
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* With a mutation rate of 0, we should keep exactly the same set
     * of mutations */
    ret = interval_map_alloc_single(&rate_map, 1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy));
    mutgen_free(&mutgen);

    /* and, with discrete sites */
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_KEEP_SITES | MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy));
    mutgen_free(&mutgen);

    /* Turn up the mutation rate */
    rate_map.value[0] = 10;
    gsl_rng_set(rng, 2);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites.num_rows > copy.sites.num_rows);
    CU_ASSERT_TRUE(tables.mutations.num_rows > copy.mutations.num_rows);
    mutgen_free(&mutgen);
    /* and, discrete sites */
    gsl_rng_set(rng, 2);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &copy2, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy2));
    ret = mutgen_generate(&mutgen, &tables, MSP_KEEP_SITES | MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites.num_rows > copy2.sites.num_rows);
    CU_ASSERT_TRUE(tables.mutations.num_rows > copy2.mutations.num_rows);
    mutgen_free(&mutgen);

    /* If we run precisely the same mutations again we should rejection
     * sample away all of the original positions */
    gsl_rng_set(rng, 2);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites.num_rows > copy.sites.num_rows);
    CU_ASSERT_TRUE(tables.mutations.num_rows > copy.mutations.num_rows);

    /* add a duplicate site to the original */
    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_TRUE(ret > 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_DUPLICATE_SITE_POSITION);
    /* and, discrete sites */
    ret = mutgen_generate(&mutgen, &tables, MSP_KEEP_SITES | MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_DUPLICATE_SITE_POSITION);

    mutgen_free(&mutgen);
    mutation_model_free(&mut_model);
    tsk_table_collection_free(&tables);
    tsk_table_collection_free(&copy);
    tsk_table_collection_free(&copy2);
    gsl_rng_free(rng);
    interval_map_free(&rate_map);
}

static void
test_single_tree_mutgen_discrete_sites(void)
{
    int ret = 0;
    int j;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    mutgen_t mutgen;
    interval_map_t rate_map;
    mutation_model_t mut_model;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* in the single tree, node 4, at time 1, is ancestral to node 0
     * so if we later add more mutations on the interval [0,1) then
     * all should be good */
    ret = tsk_site_table_add_row(&tables.sites, 0.0, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 0, 4, -1, "C", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    ret = tsk_mutation_table_add_row(&tables.mutations, 0, 4, 1, "G", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 2);

    ret = interval_map_alloc_single(&rate_map, 1, 10);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    gsl_rng_set(rng, 2);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_time_interval(&mutgen, 0.0, 1.0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.num_rows, 1);
    CU_ASSERT_FATAL(tables.mutations.num_rows > 1);
    for (j = 0; j < tables.sites.num_rows; j++) {
        CU_ASSERT_EQUAL_FATAL(tables.sites.position[j], ceil(tables.sites.position[j]));
    }
    ret = tsk_table_collection_check_integrity(&tables, TSK_CHECK_ALL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_table_collection_clear(&tables);

    /* now, keep: the single tree also has a mutation at position 0.1 */
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);

    ret = tsk_site_table_add_row(&tables.sites, 0.0, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(&tables.mutations, 0, 4, -1, "C", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    ret = tsk_mutation_table_add_row(&tables.mutations, 0, 4, 1, "G", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 2);

    ret = mutgen_generate(&mutgen, &tables, MSP_KEEP_SITES | MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.num_rows, 2);
    CU_ASSERT_FATAL(tables.mutations.num_rows > 3);
    ret = tsk_table_collection_check_integrity(&tables, TSK_CHECK_ALL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    mutgen_free(&mutgen);
    mutation_model_free(&mut_model);
    interval_map_free(&rate_map);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_single_tree_mutgen_keep_sites_many_mutations(void)
{
    int ret = 0;
    int j;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    mutgen_t mutgen;
    interval_map_t rate_map;
    mutation_model_t mut_model;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < 8192; j++) {
        ret = tsk_mutation_table_add_row(&tables.mutations, 0, 0, -1, "C", 1, NULL, 0);
        CU_ASSERT_EQUAL_FATAL(ret, j + 1);
    }

    ret = interval_map_alloc_single(&rate_map, 1, 10);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    gsl_rng_set(rng, 2);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < 10; j++) {
        ret = mutgen_generate(&mutgen, &tables, MSP_KEEP_SITES);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    CU_ASSERT_TRUE(tables.sites.num_rows > 2);
    mutgen_free(&mutgen);

    interval_map_free(&rate_map);
    mutation_model_free(&mut_model);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_single_tree_mutgen_interval(void)
{
    int ret = 0;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    interval_map_t rate_map;
    mutation_model_t mut_model;
    size_t j;
    tsk_id_t node;

    CU_ASSERT_FATAL(rng != NULL);

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = interval_map_alloc_single(&rate_map, 1, 10);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.mutations.num_rows > 0);
    mutgen_print_state(&mutgen, _devnull);

    /* End before start is an error */
    ret = mutgen_set_time_interval(&mutgen, 0, -1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    /* Setting start and end == 0 should give 0 mutations */
    ret = mutgen_set_time_interval(&mutgen, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites.num_rows == 0);
    CU_ASSERT_TRUE(tables.mutations.num_rows == 0);

    /* Setting start = 3 should give 0 mutations */
    ret = mutgen_set_time_interval(&mutgen, 3, DBL_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites.num_rows == 0);
    CU_ASSERT_TRUE(tables.mutations.num_rows == 0);

    /* Setting start = 2 should give mutations only above 4 and 5 */
    ret = mutgen_set_time_interval(&mutgen, 2, DBL_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites.num_rows > 0);
    CU_ASSERT_TRUE(tables.mutations.num_rows > 0);
    for (j = 0; j < tables.sites.num_rows; j++) {
        node = tables.mutations.node[j];
        CU_ASSERT_TRUE(node == 4 || node == 5);
    }

    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    mutation_model_free(&mut_model);
    interval_map_free(&rate_map);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_single_tree_mutgen_empty_site(void)
{
    int ret = 0;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    interval_map_t rate_map;
    mutation_model_t mut_model;

    CU_ASSERT_FATAL(rng != NULL);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = interval_map_alloc_single(&rate_map, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);
    ret = tsk_site_table_add_row(&tables.sites, 0.5, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 1);

    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    mutgen_free(&mutgen);
    interval_map_free(&rate_map);
    mutation_model_free(&mut_model);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_single_tree_mutgen_do_nothing_mutations(void)
{
    int ret = 0;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables, copy;
    interval_map_t rate_map;
    mutation_model_t mut_model;
    size_t lengths[] = { 1, 1 };
    const char *binary_alleles[] = { "0", "1" };
    double root_distribution[] = { 0.5, 0.5 };
    double transition_matrix[] = { 1.0, 0.0, 0.0, 1.0 };

    CU_ASSERT_FATAL(rng != NULL);
    ret = matrix_mutation_model_alloc(&mut_model, 2,
        (char **) (uintptr_t *) binary_alleles, lengths, root_distribution,
        transition_matrix);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = interval_map_alloc_single(&rate_map, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);
    insert_single_tree(&copy, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy));

    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy));
    mutgen_free(&mutgen);

    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.num_rows, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.num_rows, 0);
    mutgen_free(&mutgen);

    interval_map_free(&rate_map);
    mutation_model_free(&mut_model);
    tsk_table_collection_free(&tables);
    tsk_table_collection_free(&copy);
    gsl_rng_free(rng);
}

static void
test_single_tree_mutgen_many_mutations(void)
{
    int ret = 0;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    interval_map_t rate_map;
    mutation_model_t mut_model;

    CU_ASSERT_FATAL(rng != NULL);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_BINARY);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = interval_map_alloc_single(&rate_map, 1, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_BINARY);

    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_MUTATION_GENERATION_OUT_OF_ORDER);
    mutgen_free(&mutgen);

    interval_map_free(&rate_map);
    mutation_model_free(&mut_model);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static int
cmp_int64(const void *a, const void *b)
{
    int64_t *x, *y;
    x = (int64_t *) a;
    y = (int64_t *) b;
    return (int) (*x - *y);
}

static int
parse_text_int64(char *ds, tsk_size_t n)
{
    int64_t mut_id = 0;

    // not null-terminated, so we convert by hand
    while ((n > 0) && (isdigit(ds[0]))) {
        mut_id = mut_id * 10 + ds[0] - '0';
        n--;
        ds++;
    }
    CU_ASSERT_FATAL(n == 0);
    return mut_id;
}

static void
verify_slim_mutation_ids(int64_t *mut_ids, size_t mut_ids_length, int64_t min_mut_id)
{
    int j;
    CU_ASSERT_FATAL(mut_ids_length > 0);
    qsort(mut_ids, mut_ids_length, sizeof(*mut_ids), &cmp_int64);
    CU_ASSERT_EQUAL_FATAL(mut_ids[0], min_mut_id);
    if (mut_ids_length > 1) {
        for (j = 1; j < mut_ids_length; j++) {
            CU_ASSERT_EQUAL_FATAL(mut_ids[j] - mut_ids[j - 1], 1);
        }
    }
}

static void
verify_slim_metadata(char *metadata, size_t metadata_length)
{
    size_t n;
    int32_t *mutation_type_id_;
    float *selection_coeff_;
    int32_t *subpop_index_;
    int32_t *origin_generation_;
    int8_t *nucleotide_;

    CU_ASSERT_EQUAL_FATAL(metadata_length, 17);

    n = 0;
    mutation_type_id_ = (int32_t *) (metadata + n);
    CU_ASSERT_FATAL(*mutation_type_id_ >= 0);
    n += sizeof(int32_t);
    selection_coeff_ = (float *) (metadata + n);
    CU_ASSERT_EQUAL_FATAL(*selection_coeff_, 0.0);
    n += sizeof(float);
    subpop_index_ = (int32_t *) (metadata + n);
    CU_ASSERT_EQUAL_FATAL(*subpop_index_, TSK_NULL);
    n += sizeof(int32_t);
    origin_generation_ = (int32_t *) (metadata + n);
    CU_ASSERT_EQUAL_FATAL(*origin_generation_, 0.0);
    n += sizeof(int32_t);
    nucleotide_ = (int8_t *) (metadata + n);
    CU_ASSERT_EQUAL_FATAL(*nucleotide_, -1);
}

static void
test_mutgen_slim_mutations(void)
{
    int ret = 0;
    int j, k;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    interval_map_t rate_map;
    mutation_model_t mut_model;
    size_t len, parent_len;
    char *ds;
    int64_t mut_id, *all_mut_ids;
    int32_t mutation_type_id = 10;
    int64_t next_mutation_id = 23;

    CU_ASSERT_FATAL(rng != NULL);
    ret = slim_mutation_model_alloc(&mut_model, mutation_type_id, next_mutation_id, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = interval_map_alloc_single(&rate_map, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, -1);

    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.mutations.num_rows > 0);

    // should have empty ancestral states
    for (j = 0; j < tables.sites.num_rows; j++) {
        CU_ASSERT_EQUAL_FATAL(tables.sites.ancestral_state_offset[j], 0);
    }
    // check that derived states append unique integers,
    // counting up from next_mutation_id
    all_mut_ids = malloc(tables.mutations.num_rows * sizeof(int64_t));
    CU_ASSERT_FATAL(all_mut_ids != NULL);
    for (j = 0; j < tables.mutations.num_rows; j++) {
        ds = tables.mutations.derived_state + tables.mutations.derived_state_offset[j];
        len = (tables.mutations.derived_state_offset[j + 1]
               - tables.mutations.derived_state_offset[j]);
        k = tables.mutations.parent[j];
        if (k == TSK_NULL) {
            parent_len = 0;
        } else {
            parent_len = (tables.mutations.derived_state_offset[k + 1]
                          - tables.mutations.derived_state_offset[k]);
            CU_ASSERT_EQUAL_FATAL(memcmp(ds,
                                      tables.mutations.derived_state
                                          + tables.mutations.derived_state_offset[k],
                                      parent_len),
                0);
            CU_ASSERT_EQUAL_FATAL((ds + parent_len)[0], 44); // 44 is ',' in ascii
            parent_len++;
        }
        CU_ASSERT_FATAL(len > parent_len);
        mut_id = parse_text_int64(ds + parent_len, len - parent_len);
        all_mut_ids[j] = mut_id;
    }
    verify_slim_mutation_ids(all_mut_ids, tables.mutations.num_rows, next_mutation_id);

    // check that metadata is formed by appending slim_metadata to the previous one
    for (j = 0; j < tables.mutations.num_rows; j++) {
        len = (tables.mutations.metadata_offset[j + 1]
               - tables.mutations.metadata_offset[j]);
        k = tables.mutations.parent[j];
        if (k == TSK_NULL) {
            parent_len = 0;
        } else {
            parent_len = (tables.mutations.metadata_offset[k + 1]
                          - tables.mutations.metadata_offset[k]);
            CU_ASSERT_EQUAL_FATAL(
                memcmp(tables.mutations.metadata + tables.mutations.metadata_offset[j],
                    tables.mutations.metadata + tables.mutations.metadata_offset[k],
                    parent_len),
                0);
        }
        verify_slim_metadata(
            tables.mutations.metadata + tables.mutations.metadata_offset[j] + parent_len,
            tables.mutations.metadata_offset[j + 1] - tables.mutations.metadata_offset[j]
                - parent_len);
    }

    mutgen_print_state(&mutgen, _devnull);

    mutgen_free(&mutgen);
    free(all_mut_ids);
    interval_map_free(&rate_map);
    mutation_model_free(&mut_model);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_mutgen_slim_mutation_large_values(void)
{
    int ret = 0;
    mutgen_t mutgen;
    tsk_mutation_t mut;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    interval_map_t rate_map;
    mutation_model_t mut_model;
    char value[21]; /* longest 64 bit int has 20 digits */
    int64_t next_mutation_id;

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, -1);
    ret = interval_map_alloc_single(&rate_map, 1, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Trying to generate mutations that overflow raises an error */
    ret = slim_mutation_model_alloc(&mut_model, 0, INT64_MAX, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_MUTATION_ID_OVERFLOW);
    mutgen_free(&mutgen);
    mutation_model_free(&mut_model);

    tsk_mutation_table_clear(&tables.mutations);
    tsk_site_table_clear(&tables.sites);
    /* Try out with a large value that doesn't hit the ceiling */
    next_mutation_id = INT64_MAX - 100;
    ret = slim_mutation_model_alloc(&mut_model, 0, next_mutation_id, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.mutations.num_rows > 10);
    mutgen_free(&mutgen);
    mutation_model_free(&mut_model);

    ret = tsk_mutation_table_get_row(&tables.mutations, 0, &mut);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(mut.derived_state_length, 19);
    sprintf(value, "%" PRId64, next_mutation_id);
    CU_ASSERT_NSTRING_EQUAL(value, mut.derived_state, mut.derived_state_length);

    interval_map_free(&rate_map);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_mutgen_infinite_alleles(void)
{
    int ret = 0;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    interval_map_t rate_map;
    mutation_model_t mut_model;
    tsk_site_t site;
    tsk_mutation_t mutation;
    tsk_size_t j;
    char buff[100];

    CU_ASSERT_FATAL(rng != NULL);
    ret = interval_map_alloc_single(&rate_map, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, -1);

    ret = infinite_alleles_mutation_model_alloc(&mut_model, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.mutations.num_rows > 5);
    CU_ASSERT_TRUE(tables.sites.num_rows == 1);

    ret = tsk_site_table_get_row(&tables.sites, 0, &site);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_NSTRING_EQUAL("0", site.ancestral_state, site.ancestral_state_length);

    for (j = 0; j < tables.mutations.num_rows; j++) {
        ret = tsk_mutation_table_get_row(&tables.mutations, j, &mutation);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        sprintf(buff, "%d", j + 1);
        CU_ASSERT_NSTRING_EQUAL(
            buff, mutation.derived_state, mutation.derived_state_length);
    }
    mutgen_print_state(&mutgen, _devnull);

    mutgen_free(&mutgen);
    interval_map_free(&rate_map);
    mutation_model_free(&mut_model);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_mutgen_infinite_alleles_large_values(void)
{
    int ret = 0;
    mutgen_t mutgen;
    tsk_mutation_t mut;
    tsk_site_t site;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    interval_map_t rate_map;
    mutation_model_t mut_model;
    char value[21]; /* longest 64 bit uint has 20 digits */

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, -1);
    ret = interval_map_alloc_single(&rate_map, 1, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = infinite_alleles_mutation_model_alloc(&mut_model, UINT64_MAX, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &rate_map, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, &tables, MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites.num_rows == 1);
    CU_ASSERT_TRUE(tables.mutations.num_rows > 10);
    mutgen_print_state(&mutgen, _devnull);

    ret = tsk_site_table_get_row(&tables.sites, 0, &site);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(site.ancestral_state_length, 20);
    sprintf(value, "%" PRIu64, UINT64_MAX);
    CU_ASSERT_NSTRING_EQUAL(value, site.ancestral_state, site.ancestral_state_length);

    ret = tsk_mutation_table_get_row(&tables.mutations, 0, &mut);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_NSTRING_EQUAL("0", mut.derived_state, mut.derived_state_length);

    mutgen_free(&mutgen);
    mutation_model_free(&mut_model);
    interval_map_free(&rate_map);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
verify_simple_genic_selection_trajectory(
    double start_frequency, double end_frequency, double alpha, double dt)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = { { 0, 0.0 }, { 0, 0.0 } };
    tsk_table_collection_t tables;
    recomb_map_t recomb_map;
    size_t j, num_steps;
    double *allele_frequency, *time;

    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 0, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_alloc(&msp, 2, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_simulation_model_sweep_genic_selection(
        &msp, 0.5, start_frequency, end_frequency, alpha, dt);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* compute the trajectory */
    ret = msp.model.params.sweep.generate_trajectory(
        &msp.model.params.sweep, &msp, &num_steps, &time, &allele_frequency);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(num_steps > 1);
    CU_ASSERT_EQUAL(time[0], 0);
    CU_ASSERT_EQUAL(allele_frequency[0], end_frequency);
    CU_ASSERT_TRUE(allele_frequency[num_steps - 1] == start_frequency);

    for (j = 0; j < num_steps; j++) {
        CU_ASSERT_TRUE(allele_frequency[j] >= 0);
        CU_ASSERT_TRUE(allele_frequency[j] <= 1);
        if (j > 0) {
            CU_ASSERT_DOUBLE_EQUAL_FATAL(time[j], time[j - 1] + dt, 1e-9);
        }
    }

    free(time);
    free(allele_frequency);
    msp_free(&msp);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_genic_selection_trajectory(void)
{
    verify_simple_genic_selection_trajectory(0.1, 0.9, 0.1, 0.0125);
    verify_simple_genic_selection_trajectory(0.1, 0.9, 0.01, 0.00125);
    verify_simple_genic_selection_trajectory(0.8999, 0.9, 0.1, 0.2);
    verify_simple_genic_selection_trajectory(0.1, 0.9, 100, 0.1);
    verify_simple_genic_selection_trajectory(0.1, 0.9, 1, 10);
}

static void
test_sweep_genic_selection_bad_parameters(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = { { 0, 0.0 }, { 0, 0.0 } };
    tsk_table_collection_t tables;
    recomb_map_t recomb_map;

    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 0, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_alloc(&msp, 2, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_set_simulation_model_sweep_genic_selection(
        &msp, 0.5, -0.01, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_ALLELE_FREQUENCY);
    ret = msp_set_simulation_model_sweep_genic_selection(
        &msp, 0.5, 0.01, -0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_ALLELE_FREQUENCY);
    ret = msp_set_simulation_model_sweep_genic_selection(
        &msp, 0.5, 10.01, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_ALLELE_FREQUENCY);
    ret = msp_set_simulation_model_sweep_genic_selection(
        &msp, 0.5, 0.01, 10.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_ALLELE_FREQUENCY);

    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.01, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_TRAJECTORY_START_END);
    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.1, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_TRAJECTORY_START_END);

    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.9, 0.1, 0.0);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_TIME_DELTA);
    ret = msp_set_simulation_model_sweep_genic_selection(
        &msp, 0.5, 0.1, 0.9, 0.1, -0.01);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_TIME_DELTA);

    ret = msp_set_simulation_model_sweep_genic_selection(&msp, -0.5, 0.1, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_SWEEP_POSITION);
    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 5.0, 0.1, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_SWEEP_POSITION);

    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.9, -666, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_SWEEP_GENIC_SELECTION_ALPHA);
    /* The incorrect number of populations was specified */
    ret = msp_set_dimensions(&msp, 2, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_UNSUPPORTED_OPERATION);

    msp_free(&msp);
    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_sweep_genic_selection_events(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = { { 0, 0.0 }, { 0, 0.0 }, { 0, 0.0 } };
    tsk_table_collection_t tables;
    recomb_map_t recomb_map;

    ret = recomb_map_alloc_uniform(&recomb_map, 1.0, 0, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_alloc(&msp, 2, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_dimensions(&msp, 1, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.1, 0, 1, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_EVENTS_DURING_SWEEP);
    msp_free(&msp);

    tsk_table_collection_clear(&tables);
    samples[1].time = 0.1;
    ret = msp_alloc(&msp, 3, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_dimensions(&msp, 1, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_EVENTS_DURING_SWEEP);
    msp_free(&msp);

    recomb_map_free(&recomb_map);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
verify_sweep_genic_selection(uint32_t num_loci, double growth_rate)
{
    int j, ret;
    uint32_t n = 10;
    unsigned long seed = 133;
    msp_t msp;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables[2];
    recomb_map_t recomb_map;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, num_loci, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    memset(samples, 0, n * sizeof(sample_t));

    for (j = 0; j < 2; j++) {
        ret = tsk_table_collection_init(&tables[j], 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        gsl_rng_set(rng, seed);
        ret = msp_alloc(&msp, n, samples, &recomb_map, &tables[j], rng);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_dimensions(&msp, 1, 2);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_simulation_model_sweep_genic_selection(
            &msp, num_loci / 2, 0.1, 0.9, 0.1, 0.01);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_population_configuration(&msp, 0, 1.0, growth_rate);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_initialise(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
        CU_ASSERT_TRUE(ret >= 0);
        msp_print_state(&msp, _devnull);
        ret = msp_finalise_tables(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        msp_free(&msp);
        CU_ASSERT_EQUAL(tables[j].migrations.num_rows, 0);
        CU_ASSERT(tables[j].nodes.num_rows > 0);
        CU_ASSERT(tables[j].edges.num_rows > 0);
    }
    CU_ASSERT_TRUE(tsk_node_table_equals(&tables[0].nodes, &tables[1].nodes));
    CU_ASSERT_TRUE(tsk_edge_table_equals(&tables[0].edges, &tables[1].edges));
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(samples);
    for (j = 0; j < 2; j++) {
        tsk_table_collection_free(&tables[j]);
    }
    recomb_map_free(&recomb_map);
}

static void
test_sweep_genic_selection_single_locus(void)
{
    verify_sweep_genic_selection(1, 0.0);
    verify_sweep_genic_selection(1, 1.0);
    verify_sweep_genic_selection(1, -1.0);
}

static void
test_sweep_genic_selection_recomb(void)
{
    verify_sweep_genic_selection(100, 0.0);
    verify_sweep_genic_selection(100, 1.0);
    verify_sweep_genic_selection(100, -1.0);
}

static void
test_sweep_genic_selection_time_change(void)
{
    int j, ret;
    uint32_t n = 10;
    uint32_t num_loci = 10;
    unsigned long seed = 133234;
    double t;
    msp_t msp;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    recomb_map_t recomb_map;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    ret = recomb_map_alloc_uniform(&recomb_map, num_loci, 1, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    memset(samples, 0, n * sizeof(sample_t));

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    gsl_rng_set(rng, seed);
    ret = msp_alloc(&msp, n, samples, &recomb_map, &tables, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_dimensions(&msp, 1, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    /* Run some time and then change the model */
    ret = msp_run(&msp, 0.125, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_TIME);
    t = msp_get_time(&msp);

    for (j = 0; j < 10; j++) {

        ret = msp_set_simulation_model_sweep_genic_selection(
            &msp, num_loci / 2, 0.1, 0.9, 0.1, 0.01);
        CU_ASSERT_EQUAL(ret, 0);

        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(msp_get_time(&msp), t);

        ret = msp_set_simulation_model_hudson(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(msp_get_time(&msp), t);
    }
    msp_free(&msp);

    gsl_rng_free(rng);
    free(samples);
    tsk_table_collection_free(&tables);
    recomb_map_free(&recomb_map);
}

static void
test_interval_map(void)
{
    int ret;
    size_t j;
    interval_map_t imap;
    double position[] = { 0, 1, 2, 3, 4, 5 };
    double value[] = { 0.1, 1.1, 2.1, 3.1, 4.1, 5.1 };

    ret = interval_map_alloc(&imap, 0, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INSUFFICIENT_INTERVALS);
    interval_map_free(&imap);

    ret = interval_map_alloc(&imap, 1, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INSUFFICIENT_INTERVALS);
    interval_map_free(&imap);

    for (j = 2; j < 6; j++) {
        ret = interval_map_alloc(&imap, j, position, value);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(interval_map_get_size(&imap), j);
        CU_ASSERT_EQUAL(interval_map_get_num_intervals(&imap), j - 1);
        CU_ASSERT_EQUAL(interval_map_get_sequence_length(&imap), j - 1);
        interval_map_print_state(&imap, _devnull);
        interval_map_free(&imap);
    }

    position[1] = -1;
    ret = interval_map_alloc(&imap, 6, position, value);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_NEGATIVE_INTERVAL_POSITION);
    interval_map_free(&imap);
    position[1] = 0;

    position[0] = 1;
    ret = interval_map_alloc(&imap, 6, position, value);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INTERVAL_MAP_START_NON_ZERO);
    interval_map_free(&imap);
}

static void
test_matrix_mutation_model_errors(void)
{
    int ret;
    mutation_model_t model;
    const char *alleles[] = { "A", "B" };
    size_t lengths[] = { 1, 1 };
    double dist[] = { 0.5, 0.5 };
    double matrix[] = { 0, 1.0, 1.0, 0 };

    ret = matrix_mutation_model_alloc(&model, 0, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INSUFFICIENT_ALLELES);
    mutation_model_free(&model);

    /* Bad probability values */
    dist[0] = -1;
    ret = matrix_mutation_model_alloc(
        &model, 2, (char **) (uintptr_t) alleles, lengths, dist, matrix);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_ROOT_PROBABILITIES);
    mutation_model_free(&model);

    /* Sums to 1, but bad values */
    dist[0] = -1;
    dist[1] = 2;
    ret = matrix_mutation_model_alloc(
        &model, 2, (char **) (uintptr_t) alleles, lengths, dist, matrix);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_ROOT_PROBABILITIES);
    mutation_model_free(&model);

    dist[0] = 1.1;
    dist[1] = 0.5;
    ret = matrix_mutation_model_alloc(
        &model, 2, (char **) (uintptr_t) alleles, lengths, dist, matrix);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_ROOT_PROBABILITIES);
    mutation_model_free(&model);

    /* Must sum to 1 */
    dist[0] = 0.9;
    ret = matrix_mutation_model_alloc(
        &model, 2, (char **) (uintptr_t) alleles, lengths, dist, matrix);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_ROOT_PROBABILITIES);
    mutation_model_free(&model);

    dist[0] = 0.1;
    ret = matrix_mutation_model_alloc(
        &model, 2, (char **) (uintptr_t) alleles, lengths, dist, matrix);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_ROOT_PROBABILITIES);
    mutation_model_free(&model);

    /* Make sure dist is still OK . */
    dist[0] = 0.5;
    ret = matrix_mutation_model_alloc(
        &model, 2, (char **) (uintptr_t) alleles, lengths, dist, matrix);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    mutation_model_free(&model);

    /* Matrix values must be probablilities */
    matrix[0] = -1;
    ret = matrix_mutation_model_alloc(
        &model, 2, (char **) (uintptr_t) alleles, lengths, dist, matrix);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_TRANSITION_MATRIX);
    mutation_model_free(&model);

    /* Sums to 1, but bad values */
    matrix[0] = -1;
    matrix[1] = 2;
    ret = matrix_mutation_model_alloc(
        &model, 2, (char **) (uintptr_t) alleles, lengths, dist, matrix);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_TRANSITION_MATRIX);
    mutation_model_free(&model);

    matrix[0] = 1.1;
    matrix[1] = 1.0;
    ret = matrix_mutation_model_alloc(
        &model, 2, (char **) (uintptr_t) alleles, lengths, dist, matrix);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_TRANSITION_MATRIX);
    mutation_model_free(&model);

    matrix[0] = 0.6;
    matrix[1] = 0.5;
    ret = matrix_mutation_model_alloc(
        &model, 2, (char **) (uintptr_t) alleles, lengths, dist, matrix);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_TRANSITION_MATRIX);
    mutation_model_free(&model);

    matrix[0] = 0.5;
    matrix[1] = 0.5;
    ret = matrix_mutation_model_alloc(
        &model, 2, (char **) (uintptr_t) alleles, lengths, dist, matrix);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    mutation_model_free(&model);
}

static void
test_matrix_mutation_model_properties(void)
{
    int ret;
    mutation_model_t model;
    const char *alleles[] = { "", "BBBBB" };
    size_t lengths[] = { 0, 5 };
    double dist[] = { 0.5, 0.5 };
    double matrix[] = { 0, 1.0, 1.0, 0 };

    ret = matrix_mutation_model_alloc(
        &model, 2, (char **) (uintptr_t) alleles, lengths, dist, matrix);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(model.params.mutation_matrix.num_alleles, 2);
    CU_ASSERT_EQUAL_FATAL(model.params.mutation_matrix.allele_length[0], 0);
    CU_ASSERT_EQUAL_FATAL(model.params.mutation_matrix.allele_length[1], 5);
    CU_ASSERT_NSTRING_EQUAL(model.params.mutation_matrix.alleles[0], "", 0);
    CU_ASSERT_NSTRING_EQUAL(model.params.mutation_matrix.alleles[1], "BBBBB", 5);

    mutation_model_free(&model);
}

static void
test_slim_mutation_model_errors(void)
{
    int ret;
    mutation_model_t model;
    int32_t mutation_type_id = 0;
    int64_t next_mutation_id = 0;

    next_mutation_id--;
    ret = slim_mutation_model_alloc(&model, mutation_type_id, next_mutation_id, 0);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_SLIM_PARAMETERS);
    mutation_model_free(&model);

    mutation_type_id--;
    next_mutation_id++;
    ret = slim_mutation_model_alloc(&model, mutation_type_id, next_mutation_id, 0);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_SLIM_PARAMETERS);
    mutation_model_free(&model);
}

static void
test_slim_mutation_model_properties(void)
{
    int ret;
    mutation_model_t model;
    int32_t mutation_type_id = 1;
    int64_t next_mutation_id = 2;

    ret = slim_mutation_model_alloc(&model, mutation_type_id, next_mutation_id, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(model.params.slim_mutator.next_mutation_id, 2);
    CU_ASSERT_EQUAL_FATAL(model.params.slim_mutator.mutation_type_id, 1);

    mutation_model_free(&model);
}

static void
test_strerror(void)
{
    int j;
    const char *msg;
    int max_error_code = 1024; /* totally arbitrary */

    for (j = 1; j > -max_error_code; j--) {
        msg = msp_strerror(j);
        CU_ASSERT_FATAL(msg != NULL);
        CU_ASSERT(strlen(msg) > 0);
    }
}

static void
test_strerror_tskit(void)
{
    int tskit_errors[]
        = { TSK_ERR_NO_MEMORY, TSK_ERR_NODE_OUT_OF_BOUNDS, TSK_ERR_EDGE_OUT_OF_BOUNDS };
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
    fprintf(stderr, "CUnit error occured: %d: %s\n", CU_get_error(), CU_get_error_msg());
    exit(EXIT_FAILURE);
}

int
main(int argc, char **argv)
{
    int ret;
    CU_pTest test;
    CU_pSuite suite;
    CU_TestInfo tests[] = {
        { "test_fenwick", test_fenwick },
        { "test_fenwick_expand", test_fenwick_expand },
        { "test_fenwick_zero_values", test_fenwick_zero_values },
        { "test_fenwick_drift", test_fenwick_drift },
        { "test_single_locus_two_populations", test_single_locus_two_populations },
        { "test_single_locus_many_populations", test_single_locus_many_populations },
        { "test_single_locus_historical_sample", test_single_locus_historical_sample },
        { "test_single_locus_multiple_historical_samples",
            test_single_locus_multiple_historical_samples },
        { "test_dtwf_simultaneous_historical_samples",
            test_dtwf_simultaneous_historical_samples },
        { "test_single_locus_historical_sample_start_time",
            test_single_locus_historical_sample_start_time },
        { "test_single_locus_historical_sample_end_time",
            test_single_locus_historical_sample_end_time },
        { "test_simulator_getters_setters", test_simulator_getters_setters },
        { "test_demographic_events", test_demographic_events },
        { "test_demographic_events_start_time", test_demographic_events_start_time },
        { "test_census_event", test_census_event },
        { "test_dtwf_unsupported_bottleneck", test_dtwf_unsupported_bottleneck },
        { "test_time_travel_error", test_time_travel_error },
        { "test_single_locus_simulation", test_single_locus_simulation },
        { "test_single_locus_gene_conversion", test_single_locus_gene_conversion },
        { "test_multi_locus_bottleneck_arg", test_multi_locus_bottleneck_arg },
        { "test_floating_point_extremes", test_floating_point_extremes },
        { "test_mixed_model_simulation", test_mixed_model_simulation },
        { "test_dtwf_deterministic", test_dtwf_deterministic },
        { "test_dtwf_zero_pop_size", test_dtwf_zero_pop_size },
        { "test_dtwf_events_between_generations", test_dtwf_events_between_generations },
        { "test_dtwf_single_locus_simulation", test_dtwf_single_locus_simulation },
        { "test_dtwf_low_recombination", test_dtwf_low_recombination },
        { "test_pedigree_single_locus_simulation",
            test_pedigree_single_locus_simulation },
        { "test_pedigree_multi_locus_simulation", test_pedigree_multi_locus_simulation },
        { "test_pedigree_specification", test_pedigree_specification },

        /* TODO we should move the likelihood tests to their own file */
        { "test_likelihood_errors", test_likelihood_errors },
        { "test_likelihood_zero_edges", test_likelihood_zero_edges },
        { "test_likelihood_three_leaves", test_likelihood_three_leaves },
        { "test_likelihood_two_mrcas", test_likelihood_two_mrcas },
        { "test_likelihood_material_overhang", test_likelihood_material_overhang },
        { "test_likelihood_material_gap", test_likelihood_material_gap },
        { "test_likelihood_recombination_in_material_gap",
            test_likelihood_recombination_in_material_gap },

        { "test_multi_locus_simulation", test_multi_locus_simulation },
        { "test_dtwf_multi_locus_simulation", test_dtwf_multi_locus_simulation },
        { "test_gene_conversion_simulation", test_gene_conversion_simulation },
        { "test_simulation_replicates", test_simulation_replicates },
        { "test_bottleneck_simulation", test_bottleneck_simulation },
        { "test_dirac_coalescent_bad_parameters", test_dirac_coalescent_bad_parameters },
        { "test_beta_coalescent_bad_parameters", test_beta_coalescent_bad_parameters },
        { "test_multiple_mergers_simulation", test_multiple_mergers_simulation },
        { "test_multiple_mergers_growth_rate", test_multiple_mergers_growth_rate },
        { "test_large_bottleneck_simulation", test_large_bottleneck_simulation },
        { "test_simple_recombination_map", test_simple_recomb_map },
        { "test_recombination_map_copy", test_recomb_map_copy },
        { "test_recombination_map_errors", test_recomb_map_errors },
        { "test_recombination_map_examples", test_recomb_map_examples },

        { "test_translate_position_and_recomb_mass",
            test_translate_position_and_recomb_mass },
        { "test_recomb_map_mass_between", test_recomb_map_mass_between },

        { "test_binary_search", test_msp_binary_interval_search },
        { "test_binary_search_repeating", test_msp_binary_interval_search_repeating },
        { "test_binary_search_edge_cases", test_msp_binary_interval_search_edge_cases },

        { "test_simulate_from_single_locus", test_simulate_from_single_locus },
        { "test_simulate_from_single_locus_replicates",
            test_simulate_from_single_locus_replicates },
        { "test_simulate_from_empty", test_simulate_from_empty },
        { "test_simulate_from_completed", test_simulate_from_completed },
        { "test_simulate_from_incompatible", test_simulate_from_incompatible },
        { "test_simulate_init_errors", test_simulate_init_errors },

        { "test_mutgen_simple_map", test_mutgen_simple_map },
        { "test_mutgen_errors", test_mutgen_errors },
        { "test_mutgen_bad_mutation_order", test_mutgen_bad_mutation_order },
        { "test_single_tree_mutgen", test_single_tree_mutgen },
        { "test_single_tree_mutgen_keep_sites", test_single_tree_mutgen_keep_sites },
        { "test_single_tree_mutgen_discrete_sites",
            test_single_tree_mutgen_discrete_sites },
        { "test_single_tree_mutgen_keep_sites_many_mutations",
            test_single_tree_mutgen_keep_sites_many_mutations },
        { "test_single_tree_mutgen_interval", test_single_tree_mutgen_interval },
        { "test_single_tree_mutgen_empty_site", test_single_tree_mutgen_empty_site },
        { "test_single_tree_mutgen_do_nothing_mutations",
            test_single_tree_mutgen_do_nothing_mutations },
        { "test_single_tree_mutgen_many_mutations",
            test_single_tree_mutgen_many_mutations },
        { "test_mutgen_slim_mutations", test_mutgen_slim_mutations },
        { "test_mutgen_slim_mutation_large_values",
            test_mutgen_slim_mutation_large_values },
        { "test_mutgen_infinite_alleles", test_mutgen_infinite_alleles },
        { "test_mutgen_infinite_alleles_large_values",
            test_mutgen_infinite_alleles_large_values },

        { "test_genic_selection_trajectory", test_genic_selection_trajectory },
        { "test_sweep_genic_selection_bad_parameters",
            test_sweep_genic_selection_bad_parameters },
        { "test_sweep_genic_selection_events", test_sweep_genic_selection_events },
        { "test_sweep_genic_selection_single_locus",
            test_sweep_genic_selection_single_locus },
        { "test_sweep_genic_selection_recomb", test_sweep_genic_selection_recomb },
        { "test_sweep_genic_selection_time_change",
            test_sweep_genic_selection_time_change },

        { "test_interval_map", test_interval_map },

        { "test_matrix_mutation_model_errors", test_matrix_mutation_model_errors },
        { "test_matrix_mutation_model_properties",
            test_matrix_mutation_model_properties },
        { "test_slim_mutation_model_errors", test_slim_mutation_model_errors },
        { "test_slim_mutation_model_properties", test_slim_mutation_model_properties },

        { "test_strerror", test_strerror },
        { "test_strerror_tskit", test_strerror_tskit },
        CU_TEST_INFO_NULL,
    };

    /* We use initialisers here as the struct definitions change between
     * versions of CUnit */
    CU_SuiteInfo suites[] = {
        { .pName = "msprime",
            .pInitFunc = msprime_suite_init,
            .pCleanupFunc = msprime_suite_cleanup,
            .pTests = tests },
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
