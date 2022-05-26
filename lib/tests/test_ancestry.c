/*
** Copyright (C) 2016-2021 University of Oxford
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

#include "testlib.h"

static void
test_single_locus_simulation(void)
{
    int ret;
    int model;
    const char *model_name;
    uint32_t j;
    uint32_t n = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(&msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);

    memset(samples, 0, n * sizeof(sample_t));
    ret = build_sim(&msp, &tables, rng, 1, 1, samples, n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* For the single locus sim we should have exactly n - 1 events */
    for (j = 0; j < n - 2; j++) {
        ret = msp_run(&msp, DBL_MAX, 1);
        CU_ASSERT_EQUAL_FATAL(ret, MSP_EXIT_MAX_EVENTS);
        msp_verify(&msp, 0);
    }
    ret = msp_run(&msp, DBL_MAX, 1);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);

    model = msp_get_model(&msp)->type;
    CU_ASSERT_EQUAL(model, MSP_MODEL_HUDSON);
    model_name = msp_get_model_name(&msp);
    CU_ASSERT_STRING_EQUAL(model_name, "hudson");

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(samples);
    tsk_table_collection_free(&tables);
}

static void
test_single_locus_two_populations(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
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

    ret = build_sim(&msp, &tables, rng, 1, 2, samples, n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_discrete_genome(&msp, false), 0);
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
    tsk_table_collection_free(&tables);
}

static void
test_single_locus_many_populations(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    uint32_t num_populations = 500;
    sample_t samples[] = { { 0, 0.0 }, { num_populations - 1, 0.0 } };
    uint32_t n = 2;
    tsk_table_collection_t tables;

    ret = build_sim(&msp, &tables, rng, 1, num_populations, samples, n);
    CU_ASSERT_EQUAL(ret, 0);
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
    tsk_table_collection_free(&tables);
}

static void
test_single_locus_labels(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = { { 0, 0.0 }, { 0, 0.0 }, { 0, 0.0 } };
    tsk_size_t num_labels = 10;
    uint32_t n = 3;
    tsk_table_collection_t tables;

    /* Make sure the logic around reallocating the mass indexes
     * is robust. We need to have a non-zero recombination rate to
     * do this */
    ret = build_sim(&msp, &tables, rng, 1, 1, samples, n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_recombination_rate(&msp, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_simulation_model_hudson(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_num_labels(&msp, num_labels);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_simulation_model_hudson(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, 1);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_EVENTS);
    msp_verify(&msp, 0);
    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL_FATAL(msp_get_num_edges(&msp), 4);
    CU_ASSERT_EQUAL_FATAL(msp_get_num_nodes(&msp), 5);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_single_locus_historical_sample(void)
{
    int ret, j;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = { { 0, 0.0 }, { 0, 101.0 } };
    tsk_edge_table_t *edges;
    tsk_node_table_t *nodes;
    uint32_t n = 2;
    tsk_table_collection_t tables;

    gsl_rng_set(rng, 5);

    for (j = 0; j < 2; j++) {
        ret = build_sim(&msp, &tables, rng, 1, 1, samples, n);
        CU_ASSERT_EQUAL(ret, 0);

        if (j == 0) {
            ret = msp_set_simulation_model_hudson(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            ret = msp_set_simulation_model_dtwf(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        }
        ret = msp_set_population_configuration(&msp, 0, 100, 0, true);
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
        tsk_table_collection_free(&tables);
    }
    gsl_rng_free(rng);
}

static void
test_single_locus_all_historical(void)
{
    int ret, j;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = { { 0, 10.0 }, { 0, 101.0 } };
    tsk_edge_table_t *edges;
    tsk_node_table_t *nodes;
    uint32_t n = 2;
    tsk_table_collection_t tables;

    gsl_rng_set(rng, 5);

    for (j = 0; j < 2; j++) {
        ret = build_sim(&msp, &tables, rng, 1, 1, samples, n);
        CU_ASSERT_EQUAL(ret, 0);

        if (j == 0) {
            ret = msp_set_simulation_model_hudson(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            ret = msp_set_simulation_model_dtwf(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        }
        ret = msp_set_population_configuration(&msp, 0, 100, 0, true);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_initialise(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(msp.time, 10.0);

        msp_print_state(&msp, _devnull);
        ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        msp_verify(&msp, 0);
        msp_print_state(&msp, _devnull);

        CU_ASSERT_EQUAL_FATAL(msp_get_num_nodes(&msp), 3);
        nodes = &msp.tables->nodes;
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(nodes->time[0], 10.0);
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
        tsk_table_collection_free(&tables);
    }
    gsl_rng_free(rng);
}

static void
test_single_locus_multiple_historical_samples(void)
{
    int ret, j;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = { { 0, 0.0 }, { 0, 10.0 }, { 0, 10.0 }, { 0, 10.0 } };
    uint32_t n = 4;
    tsk_table_collection_t tables;

    gsl_rng_set(rng, 5);
    for (j = 0; j < 2; j++) {
        ret = build_sim(&msp, &tables, rng, 1, 1, samples, n);
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
        ret = msp_run(&msp, 10 - 1E-6, ULONG_MAX);
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
        tsk_table_collection_free(&tables);
    }
    gsl_rng_free(rng);
}

static void
test_single_locus_historical_sample_start_time(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = { { 0, 0.0 }, { 0, 10.0 } };
    tsk_edge_table_t *edges;
    tsk_node_table_t *nodes;
    uint32_t n = 2;
    size_t j, model;
    tsk_table_collection_t tables;
    double start_times[] = { -1, 0, 2, 10, 10.0001, 1000 };

    for (model = 0; model < 2; model++) {
        for (j = 0; j < sizeof(start_times) / sizeof(double); j++) {

            ret = build_sim(&msp, &tables, rng, 1, 1, samples, n);
            CU_ASSERT_EQUAL(ret, 0);
            if (model == 0) {
                ret = msp_set_simulation_model_hudson(&msp);
                CU_ASSERT_EQUAL(ret, 0);
            } else {
                ret = msp_set_simulation_model_dtwf(&msp);
                CU_ASSERT_EQUAL(ret, 0);
            }
            ret = msp_set_population_configuration(&msp, 0, 100, 0, true);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_set_start_time(&msp, start_times[j]);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_initialise(&msp);
            CU_ASSERT_EQUAL(ret, 0);
            /* Start time cannot be less than the time of the first sample */
            CU_ASSERT_EQUAL(msp.time, GSL_MAX(start_times[j], 0));

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
            tsk_table_collection_free(&tables);
        }
    }
    gsl_rng_free(rng);
}

static void
test_single_locus_historical_sample_end_time(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = { { 0, 0.0 }, { 0, 10.0 } };
    tsk_node_table_t *nodes;
    uint32_t n = 2;
    size_t model;
    tsk_table_collection_t tables;

    for (model = 0; model < 2; model++) {
        ret = build_sim(&msp, &tables, rng, 1, 1, samples, n);
        CU_ASSERT_EQUAL(ret, 0);
        if (model == 0) {
            ret = msp_set_simulation_model_hudson(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            ret = msp_set_simulation_model_dtwf(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        }
        ret = msp_set_population_configuration(&msp, 0, 100, 0, true);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_initialise(&msp);
        CU_ASSERT_EQUAL(ret, 0);

        msp_print_state(&msp, _devnull);
        ret = msp_run(&msp, 10.0 - 1E-6, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_TIME);
        /* msp_verify(&msp); */

        /* The second sampling event should *not* have been applied */
        CU_ASSERT_EQUAL_FATAL(msp.next_sampling_event, 1);
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
        tsk_table_collection_free(&tables);
    }
    gsl_rng_free(rng);
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
    tsk_table_collection_t tables;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();

    gsl_rng_set(rng, seed);

    for (j = 0; j < sizeof(models) / sizeof(int); j++) {
        for (k = 0; k < sizeof(store_full_arg) / sizeof(bool); k++) {

            ret = build_sim(&msp, &tables, rng, m, 2, NULL, n);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 1.0 / m), 0);
            ret = msp_set_migration_matrix(&msp, 4, migration_matrix);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_set_store_migrations(&msp, true);
            CU_ASSERT_EQUAL(ret, 0);
            /* set all the block sizes to something small to provoke the memory
             * expansions. */
            ret = msp_set_avl_node_block_size(&msp, 1);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_set_node_mapping_block_size(&msp, 1);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_set_segment_block_size(&msp, 1);
            CU_ASSERT_EQUAL(ret, 0);
            switch (j) {
                case 0:
                    ret = msp_set_simulation_model_hudson(&msp);
                    break;
                case 1:
                    ret = msp_set_simulation_model_smc(&msp);
                    break;
                case 2:
                    ret = msp_set_simulation_model_smc_prime(&msp);
                    break;
            }
            ret = msp_set_store_full_arg(&msp, store_full_arg[k]);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            ret = msp_initialise(&msp);
            CU_ASSERT_EQUAL(ret, 0);

            num_events = 0;
            while ((ret = msp_run(&msp, DBL_MAX, 1)) == 1) {
                msp_verify(&msp, MSP_VERIFY_BREAKPOINTS);
                num_events++;
            }
            CU_ASSERT_EQUAL(ret, 0);
            msp_verify(&msp, MSP_VERIFY_BREAKPOINTS);
            ret = msp_get_num_migration_events(&msp, migration_events);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT(num_events > n - 1);
            CU_ASSERT_EQUAL(1 + num_events,
                /* diagonals must be zero here */
                migration_events[1] + migration_events[2]
                    + msp_get_num_recombination_events(&msp)
                    + msp_get_num_common_ancestor_events(&msp)
                    + msp_get_num_rejected_common_ancestor_events(&msp));
            if (models[j] == MSP_MODEL_HUDSON) {
                CU_ASSERT_EQUAL(msp_get_num_rejected_common_ancestor_events(&msp), 0);
            }

            gsl_rng_set(rng, seed);
            ret = msp_reset(&msp);
            num_events = 0;
            while ((ret = msp_run(&msp, DBL_MAX, 1)) == 1) {
                msp_verify(&msp, MSP_VERIFY_BREAKPOINTS);
                num_events++;
            }
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_get_num_migration_events(&msp, migration_events);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(
                1 + num_events, migration_events[1] + migration_events[2]
                                    + msp_get_num_recombination_events(&msp)
                                    + msp_get_num_common_ancestor_events(&msp)
                                    + msp_get_num_rejected_common_ancestor_events(&msp));
            if (models[j] == MSP_MODEL_HUDSON) {
                CU_ASSERT_EQUAL(msp_get_num_rejected_common_ancestor_events(&msp), 0);
            }
            model = msp_get_model(&msp)->type;
            CU_ASSERT_EQUAL(model, models[j]);
            model_name = msp_get_model_name(&msp);
            CU_ASSERT_STRING_EQUAL(model_name, model_names[j]);

            ret = msp_free(&msp);
            CU_ASSERT_EQUAL(ret, 0);
            tsk_table_collection_free(&tables);
        }
    }
    gsl_rng_free(rng);
}

static void
test_multi_locus_bottleneck_arg(void)
{
    int ret;
    uint32_t n = 2;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    ret = build_sim(&msp, &tables, rng, 100, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 10), 0);
    ret = msp_add_instantaneous_bottleneck(&msp, 1.0, 0, 10);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_store_full_arg(&msp, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
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
test_dtwf_single_locus_simulation(void)
{
    int ret;
    const char *model_name;
    int i;
    uint32_t n = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    int num_coalescent_events = 0;
    int num_children;
    tsk_table_collection_t tables;

    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    model_name = msp_get_model_name(&msp);
    CU_ASSERT_STRING_EQUAL(model_name, "dtwf");

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);

    /* For the single locus sim we should have n-1 coalescent events,
     * counting multiple mergers as multiple coalescent events */
    for (i = 0; i < msp.tables->nodes.num_rows; i++) {
        num_children = get_num_children(i, &msp.tables->edges);
        if (num_children > 0) {
            num_coalescent_events += num_children - 1;
        }
    }
    CU_ASSERT_EQUAL(num_coalescent_events, n - 1);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_dtwf_multi_locus_simulation(void)
{
    int ret;
    uint32_t n = 100;
    long seed = 10;
    double migration_matrix[] = { 0, 0.1, 0.1, 0 };
    const char *model_name;
    size_t num_ca_events, num_re_events;
    double t;
    tsk_table_collection_t tables;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();

    gsl_rng_set(rng, seed);
    ret = build_sim(&msp, &tables, rng, 10, 2, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 0.1), 0);
    ret = msp_set_population_configuration(&msp, 0, n, 0, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_migration_matrix(&msp, 4, migration_matrix);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_store_migrations(&msp, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    model_name = msp_get_model_name(&msp);
    CU_ASSERT_STRING_EQUAL(model_name, "dtwf");

    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    msp_verify(&msp, 0);
    num_ca_events = msp_get_num_common_ancestor_events(&msp);
    num_re_events = msp_get_num_recombination_events(&msp);
    CU_ASSERT_TRUE(num_ca_events > 0);
    CU_ASSERT_TRUE(num_re_events > 0);
    CU_ASSERT_EQUAL(ret, 0);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* Realloc the simulator under different memory params to see if
     * we get the same result. */
    gsl_rng_set(rng, seed);
    ret = build_sim(&msp, &tables, rng, 10, 2, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 0.1), 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, n, 0, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_migration_matrix(&msp, 4, migration_matrix);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    t = 1;
    /* We should be able to step backward here generation-by-generation until
     * coalescence.
     */
    while ((ret = msp_run(&msp, DBL_MAX, 1)) > 0) {
        msp_verify(&msp, 0);
        CU_ASSERT_EQUAL_FATAL(msp.time, t);
        t++;
    }
    msp_verify(&msp, 0);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(num_ca_events == msp_get_num_common_ancestor_events(&msp));
    CU_ASSERT_TRUE(num_re_events == msp_get_num_recombination_events(&msp));

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_dtwf_deterministic(void)
{
    int j, ret;
    uint32_t n = 10;
    unsigned long seed = 133;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables[2];

    for (j = 0; j < 2; j++) {
        gsl_rng_set(rng, seed);
        ret = build_sim(&msp, &tables[j], rng, 2, 1, NULL, n);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_simulation_model_dtwf(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_recombination_rate(&msp, 1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_population_configuration(&msp, 0, n, 0, true);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_initialise(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        msp_verify(&msp, 0);
        ret = msp_finalise_tables(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        msp_free(&msp);
        CU_ASSERT_EQUAL(tables[j].migrations.num_rows, 0);
        CU_ASSERT(tables[j].nodes.num_rows > 0);
        CU_ASSERT(tables[j].edges.num_rows > 0);
    }
    CU_ASSERT_TRUE(tsk_node_table_equals(&tables[0].nodes, &tables[1].nodes, 0));
    CU_ASSERT_TRUE(tsk_edge_table_equals(&tables[0].edges, &tables[1].edges, 0));

    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    for (j = 0; j < 2; j++) {
        tsk_table_collection_free(&tables[j]);
    }
}

static void
test_dtwf_simultaneous_historical_samples(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = { { 0, 0 }, { 0, 1.1 }, { 0, 1.2 } };
    tsk_node_table_t *nodes;
    uint32_t n = 3;
    tsk_table_collection_t tables;

    gsl_rng_set(rng, 5);

    ret = build_sim(&msp, &tables, rng, 1, 1, samples, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 100, 0, true);
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
}

static void
test_dtwf_low_recombination(void)
{
    int ret;
    uint32_t n = 2;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_discrete_genome(&msp, false), 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 1e-9), 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, n, 0, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
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
    gsl_rng *rng = safe_rng_alloc();
    double migration_matrix[] = { 0, 0, 0, 0 };
    tsk_table_collection_t tables;

    for (j = 0; j < n; j++) {
        samples[j].time = 0;
        samples[j].population = j % 2;
    }

    CU_ASSERT_FATAL(samples != NULL);

    ret = build_sim(&msp, &tables, rng, 5, 2, samples, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_recombination_rate(&msp, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 10, 0, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 1, 10, 0, true);
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
    tsk_table_collection_free(&tables);
}

static void
test_dtwf_unsupported_bottleneck(void)
{
    int ret;
    uint32_t n = 10;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    msp_t msp;

    ret = build_sim(&msp, &tables, rng, 1.0, 1, NULL, n);
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
    tsk_table_collection_free(&tables);
}

static void
test_dtwf_zero_pop_size(void)
{
    int ret;
    uint32_t n = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    /* DTWF population size must round to >= 1 */
    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_recombination_rate(&msp, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 0.4, 0, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_DTWF_ZERO_POPULATION_SIZE);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* With no high growth rate, population sizes crash to zero. */
    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_recombination_rate(&msp, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 10, 100, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_DTWF_ZERO_POPULATION_SIZE);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    gsl_rng_free(rng);
}

static void
test_dtwf_migration_matrix_not_stochastic(void)
{
    int ret;
    uint32_t n = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    double migration_matrix[] = { 0, .1, .1, .1, 0, .1, .1, .1, 0 };

    /* Rows of migration matrix must sum to <=1 in DTWF */
    ret = build_sim(&msp, &tables, rng, 1, 3, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_recombination_rate(&msp, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 1, 10, 0, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 2, 10, 0, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_migration_matrix(&msp, 9, migration_matrix);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_migration_rate_change(&msp, 0, 0, 1, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_DTWF_MIGRATION_MATRIX_NOT_STOCHASTIC);

    /* With a second migration rate change in same generation to satisfy constraint */
    ret = msp_add_migration_rate_change(&msp, 0, 0, 2, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_reset(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_dtwf_errors(void)
{
    int ret;
    uint32_t n = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    /* GC is supported */
    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_gene_conversion_rate(&msp, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_gene_conversion_tract_length(&msp, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_DTWF_GC_NOT_SUPPORTED);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* Any ploidy other than 2 is an error */
    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_ploidy(&msp, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_DTWF_DIPLOID_ONLY);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    gsl_rng_free(rng);
}

static void
test_mixed_hudson_dtwf(void)
{
    int ret;
    uint32_t j, k;
    uint32_t n = 10;
    double N = 100;
    int model;
    double initial_size, growth_rate;
    int state;
    const char *model_name;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    double g = -1.0 / 8192;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();

    ret = build_sim(&msp, &tables, rng, 10, 3, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_recombination_rate(&msp, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    /* Set the populations to 1, 2, and 3N. We don't simulate them,
     * but they should be equal at the end */
    ret = msp_set_population_configuration(&msp, 0, N, 0, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 1, 2 * N, g, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 2, 3 * N, 2 * g, true);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    /* Run for 10 generations each for the different models, alternating */
    j = 0;
    while ((ret = msp_run(&msp, j * 10, UINT32_MAX)) == 2) {
        msp_verify(&msp, 0);
        /* Check that our populations and growth rates are still correct */
        for (k = 0; k < 3; k++) {
            ret = msp_get_population_configuration(
                &msp, k, &initial_size, &growth_rate, &state);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(initial_size, (k + 1) * N);
            CU_ASSERT_EQUAL_FATAL(growth_rate, k * g);
            CU_ASSERT_EQUAL(state, MSP_POP_STATE_ACTIVE);
        }
        CU_ASSERT_FALSE(msp_is_completed(&msp));
        if (j % 2 == 1) {
            model = msp_get_model(&msp)->type;
            CU_ASSERT_EQUAL(model, MSP_MODEL_HUDSON);
            model_name = msp_get_model_name(&msp);
            CU_ASSERT_STRING_EQUAL(model_name, "hudson");
            ret = msp_set_simulation_model_dtwf(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            model = msp_get_model(&msp)->type;
            CU_ASSERT_EQUAL(model, MSP_MODEL_DTWF);
            model_name = msp_get_model_name(&msp);
            CU_ASSERT_STRING_EQUAL(model_name, "dtwf");
            ret = msp_set_simulation_model_hudson(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        }
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(msp_is_completed(&msp));
    CU_ASSERT_TRUE(j > 10);

    ret = msp_finalise_tables(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_print_state(&ts, _devnull);
    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_mixed_hudson_smc(void)
{
    int ret;
    uint32_t j;
    uint32_t n = 10;
    int model;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();

    ret = build_sim(&msp, &tables, rng, 10, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_segment_block_size(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_recombination_rate(&msp, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_gene_conversion_rate(&msp, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_gene_conversion_tract_length(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    /* Run for 1 event each for the models, interleaving */
    j = 0;
    while ((ret = msp_run(&msp, DBL_MAX, 1)) == MSP_EXIT_MAX_EVENTS) {
        msp_verify(&msp, 0);
        CU_ASSERT_FALSE(msp_is_completed(&msp));
        model = msp_get_model(&msp)->type;
        if (j % 2 == 1) {
            CU_ASSERT_EQUAL(model, MSP_MODEL_SMC);
            ret = msp_set_simulation_model_hudson(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            CU_ASSERT_EQUAL(model, MSP_MODEL_HUDSON);
            ret = msp_set_simulation_model_smc(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        }
        if (j == 10) {
            msp_print_state(&msp, _devnull);
        }
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(msp_is_completed(&msp));
    CU_ASSERT_TRUE(j > 10);
    CU_ASSERT_TRUE(msp_get_num_recombination_events(&msp) > 1);
    CU_ASSERT_TRUE(msp_get_num_gene_conversion_events(&msp) > 1);

    ret = msp_finalise_tables(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_print_state(&ts, _devnull);
    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
run_gc_simulation(double sequence_length, double gc_rate, double tract_length,
    double recombination_rate, bool discrete_genome, bool full_arg)
{
    int ret;
    uint32_t n = 10;
    long seed = 10;
    size_t num_events, num_ca_events, num_re_events, num_gc_events,
        num_internal_gc_events;
    double sum_internal_gc_tract_lengths;
    bool single_locus = sequence_length == 1 && discrete_genome;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();

    gsl_rng_set(rng, seed);

    ret = build_sim(&msp, &tables, rng, sequence_length, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_discrete_genome(&msp, discrete_genome), 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, recombination_rate), 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_gene_conversion_rate(&msp, gc_rate), 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_gene_conversion_tract_length(&msp, tract_length), 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_store_full_arg(&msp, full_arg), 0);
    /* Set a very small block size to force lots of fenwick tree rebuilds */
    CU_ASSERT_EQUAL_FATAL(msp_set_segment_block_size(&msp, 16), 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    num_events = 0;
    while ((ret = msp_run(&msp, DBL_MAX, 1)) == 1) {
        msp_verify(&msp, MSP_VERIFY_BREAKPOINTS);
        num_events++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    msp_verify(&msp, MSP_VERIFY_BREAKPOINTS);
    num_ca_events = msp_get_num_common_ancestor_events(&msp);
    CU_ASSERT_TRUE(num_ca_events > 0);
    num_re_events = msp_get_num_recombination_events(&msp);
    if (recombination_rate == 0 || single_locus) {
        CU_ASSERT_EQUAL(num_re_events, 0);
    } else {
        CU_ASSERT_TRUE(num_re_events > 0);
    }
    num_gc_events = msp_get_num_gene_conversion_events(&msp);
    if (gc_rate == 0 || single_locus) {
        CU_ASSERT_EQUAL(num_gc_events, 0);
    } else {
        CU_ASSERT_TRUE(num_gc_events > 0);
    }
    num_internal_gc_events = msp_get_num_internal_gene_conversion_events(&msp);
    CU_ASSERT_TRUE(num_internal_gc_events >= num_gc_events);
    sum_internal_gc_tract_lengths = msp_get_sum_internal_gc_tract_lengths(&msp);
    if (discrete_genome) {
        CU_ASSERT_TRUE(sum_internal_gc_tract_lengths >= num_internal_gc_events);
    }
    msp_free(&msp);

    /* Make sure we can build a tree sequence out of the result */
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);

    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_gc_single_locus(void)
{
    run_gc_simulation(10, 0.0, 1.0, 0.0, false, false);
    run_gc_simulation(1, 1.0, 1.0, 0.0, true, false);
    run_gc_simulation(10, 0.0, 1.0, 0.0, false, true);
    run_gc_simulation(1, 1.0, 1.0, 0.0, true, true);
}

static void
test_gc_tract_lengths(void)
{
    double tract_lengths[] = { 1.0, 1.3333, 5, 10 };
    size_t j;

    for (j = 0; j < sizeof(tract_lengths) / sizeof(double); j++) {
        run_gc_simulation(10, 1.0, tract_lengths[j], 0.1, true, false);
        run_gc_simulation(10, 1.0, tract_lengths[j], 0.1, false, false);
        run_gc_simulation(10, 1.0, tract_lengths[j], 0.1, true, true);
        run_gc_simulation(10, 1.0, tract_lengths[j], 0.1, false, true);
    }
}

static void
test_gc_zero_recombination(void)
{
    run_gc_simulation(10, 1.0, 5, 0.0, true, false);
    run_gc_simulation(10, 1.0, 1, 0.0, false, false);
    run_gc_simulation(10, 1.0, 5, 0.0, true, true);
    run_gc_simulation(10, 1.0, 1, 0.0, false, true);
}

static void
test_gc_rates(void)
{
    size_t full_arg;
    for (full_arg = 0; full_arg < 2; full_arg++) {
        run_gc_simulation(10, 0.1, 1.5, 1.0, false, full_arg);
        run_gc_simulation(5, 5.0, 2.5, 0.01, false, full_arg);
        run_gc_simulation(10, 1, 1, 1.0, false, full_arg);
        run_gc_simulation(30, 1.0, 6.5, 1.0, true, full_arg);
    }
}

static void
test_multiple_mergers_simulation(void)
{
    int ret;
    size_t j, k, o, p, q;
    uint32_t n = 10;
    long seed = 10;
    bool store_full_arg[] = { true, false };
    /* These simulations can be slow, so just choose a few param combinations */
    double beta_params[][2] = { { 1.1, 0.5 }, { 1.99, 1 } };
    /* TODO what are good psi parameters here? */
    double psi_params[][2] = { { 0.9, 10 }, { 0.1, 1 } };
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    for (j = 0; j < 2; j++) {
        if (j == 0) {
            o = sizeof(psi_params) / sizeof(*psi_params);
        } else if (j == 1) {
            o = sizeof(beta_params) / sizeof(*beta_params);
        }
        for (k = 0; k < sizeof(store_full_arg) / sizeof(bool); k++) {
            for (p = 0; p < o; p++) {
                for (q = 1; q < 3; q++) {
                    gsl_rng_set(rng, seed);
                    ret = build_sim(&msp, &tables, rng, 10, 1, NULL, n);
                    CU_ASSERT_EQUAL(ret, 0);
                    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 1), 0);
                    if (j == 0) {
                        ret = msp_set_simulation_model_dirac(
                            &msp, psi_params[p][0], psi_params[p][1]);
                    } else {
                        ret = msp_set_simulation_model_beta(
                            &msp, beta_params[p][0], beta_params[p][1]);
                    }
                    msp_set_ploidy(&msp, q);
                    CU_ASSERT_EQUAL(ret, 0);
                    /* TODO check for adding various complications like multiple
                     * populations etc to ensure they fail.
                     */
                    ret = msp_set_store_full_arg(&msp, store_full_arg[k]);
                    CU_ASSERT_EQUAL(ret, 0);
                    ret = msp_initialise(&msp);
                    CU_ASSERT_EQUAL(ret, 0);
                    msp_print_state(&msp, _devnull);

                    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                    CU_ASSERT_TRUE(msp_is_completed(&msp));
                    CU_ASSERT_TRUE(msp.time > 0);
                    msp_verify(&msp, 0);

                    msp_reset(&msp);
                    while ((ret = msp_run(&msp, DBL_MAX, 1)) == 1) {
                        msp_verify(&msp, 0);
                    }
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                    CU_ASSERT_TRUE(msp_is_completed(&msp));

                    ret = msp_free(&msp);
                    CU_ASSERT_EQUAL(ret, 0);
                    tsk_table_collection_free(&tables);
                }
            }
        }
    }
    gsl_rng_free(rng);
}

static void
test_multiple_mergers_growth_rate(void)
{
    int ret;
    size_t n = 10;
    int j, k;
    gsl_rng *rng = safe_rng_alloc();
    msp_t msp;
    tsk_table_collection_t tables;

    for (j = 0; j < 2; j++) {
        for (k = 1; k < 3; k++) {
            ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
            CU_ASSERT_EQUAL(ret, 0);
            if (j == 0) {
                ret = msp_set_simulation_model_dirac(&msp, 0.9, 100);
                CU_ASSERT_EQUAL(ret, 0);
            } else {
                ret = msp_set_simulation_model_beta(&msp, 1.9, 1);
                CU_ASSERT_EQUAL(ret, 0);
            }
            msp_set_ploidy(&msp, k);

            /* Set to a nonzero growth_rate */
            ret = msp_set_population_configuration(&msp, 0, 1, -0.01, true);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            ret = msp_initialise(&msp);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
            CU_ASSERT_EQUAL(ret, 0);

            msp_print_state(&msp, _devnull);
            msp_free(&msp);
            tsk_table_collection_free(&tables);
        }
    }
    gsl_rng_free(rng);
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
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
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
    tsk_table_collection_free(&tables);
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
    double truncation_points[] = { -1e6, 0 };
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
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
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_simulator_getters_setters(void)
{
    int ret;
    uint32_t j;
    uint32_t n = 10;
    uint32_t m = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = safe_rng_alloc();
    double migration_matrix[] = { 0, 0, 0, 0 };
    double matrix[4], growth_rate, initial_size;
    int state;
    double Ne = 4;
    size_t migration_events[4];
    size_t breakpoints[m];
    double position[] = { 0, 1.0 };
    double rate = 0;
    tsk_table_collection_t tables;
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);

    for (j = 0; j < n; j++) {
        samples[j].time = j;
        samples[j].population = j % 2;
    }
    CU_ASSERT_EQUAL(msp_alloc(&msp, NULL, rng), MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);
    CU_ASSERT_EQUAL(msp_alloc(&msp, &tables, NULL), MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);

    samples[0].time = 1.0;
    ret = build_sim(&msp, &tables, rng, m, 2, samples, n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    samples[0].population = TSK_NULL;
    ret = build_sim(&msp, &tables, rng, m, 2, samples, n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    msp_free(&msp);
    tsk_table_collection_free(&tables);
    samples[0].population = 0;

    ret = build_sim(&msp, &tables, rng, m, 2, samples, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(msp_set_ploidy(&msp, -1), MSP_ERR_BAD_PLOIDY);
    CU_ASSERT_EQUAL(msp_set_ploidy(&msp, 0), MSP_ERR_BAD_PLOIDY);
    CU_ASSERT_EQUAL(msp_set_num_labels(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_node_mapping_block_size(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_segment_block_size(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_avl_node_block_size(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_population_configuration(&msp, -1, 0, 0, true),
        MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL(msp_set_population_configuration(&msp, 3, 0, 0, true),
        MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL(msp_set_recombination_rate(&msp, -1), MSP_ERR_BAD_RATE_VALUE);
    CU_ASSERT_EQUAL(
        msp_set_recombination_map(&msp, 1, position, &rate), MSP_ERR_BAD_RATE_MAP);

    ret = msp_set_gene_conversion_rate(&msp, -1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RATE_VALUE);
    ret = msp_set_gene_conversion_tract_length(&msp, -1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = msp_set_gene_conversion_tract_length(&msp, tables.sequence_length + 1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
        msp_set_gene_conversion_map(&msp, 1, position, &rate), MSP_ERR_BAD_RATE_MAP);

    ret = msp_set_start_time(&msp, NAN);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_START_TIME);

    ret = msp_set_simulation_model_hudson(&msp);
    CU_ASSERT_EQUAL(msp_get_model(&msp)->type, MSP_MODEL_HUDSON);
    ret = msp_set_population_configuration(&msp, 0, Ne, 0, true);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(msp_get_population_configuration(&msp, 3, NULL, NULL, NULL),
        MSP_ERR_POPULATION_OUT_OF_BOUNDS);
    ret = msp_set_population_configuration(&msp, 0, 2 * Ne, 0.5, true);
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

    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 1), 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_get_population_configuration(&msp, 0, &initial_size, &growth_rate, &state);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(initial_size, 2 * Ne);
    CU_ASSERT_EQUAL(growth_rate, 0.5);
    CU_ASSERT_EQUAL(state, MSP_POP_STATE_ACTIVE);

    CU_ASSERT_TRUE(msp_get_store_migrations(&msp));
    CU_ASSERT_EQUAL(msp_get_num_avl_node_blocks(&msp), 1);
    CU_ASSERT_EQUAL(msp_get_num_node_mapping_blocks(&msp), 1);
    CU_ASSERT_EQUAL(msp_get_num_segment_blocks(&msp), 1);
    CU_ASSERT_EQUAL(msp_get_num_populations(&msp), 2);

    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(msp_get_num_breakpoints(&msp), m - 1);
    ret = msp_get_breakpoints(&msp, breakpoints);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
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
    gsl_rng *rng = safe_rng_alloc();
    double migration_matrix[] = { 0, 0.1, 0.1, 0 };
    double last_time, time, pop_size;
    int pops[2];
    double probs[2] = { 0, 1 };
    tsk_table_collection_t tables;
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);

    for (j = 0; j < n; j++) {
        samples[j].time = j;
        samples[j].population = j % 2;
    }

    for (model = 0; model < 2; model++) {
        ret = build_sim(&msp, &tables, rng, m, 2, samples, n);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 1), 0);

        /* Negative population sizes are not allowed */
        ret = msp_set_population_configuration(&msp, 0, -1, 0, true);
        CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

        ret = msp_set_population_configuration(&msp, 0, 1, 0.001, true);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_population_configuration(&msp, 1, 2, 0.002, true);
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

        ret = msp_set_population_configuration(&msp, 0, 0, 0, true);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_population_configuration(&msp, 0, 100, 0, true);
        CU_ASSERT_EQUAL(ret, 0);

        CU_ASSERT_EQUAL(msp_add_mass_migration(&msp, 10, -1, 0, 1),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        CU_ASSERT_EQUAL(
            msp_add_mass_migration(&msp, 10, 2, 0, 1), MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        CU_ASSERT_EQUAL(
            msp_add_mass_migration(&msp, 10, 0, 0, 1), MSP_ERR_SOURCE_DEST_EQUAL);
        CU_ASSERT_EQUAL(
            msp_add_mass_migration(&msp, 10, 0, 1, -5), MSP_ERR_BAD_PROPORTION);

        pops[0] = 0;
        CU_ASSERT_EQUAL(msp_add_population_split(&msp, 10, 100, pops, 1),
            MSP_ERR_TOO_MANY_EVENT_POPULATIONS);
        CU_ASSERT_EQUAL(msp_add_population_split(&msp, 10, 1, pops, 2),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        CU_ASSERT_EQUAL_FATAL(
            msp_add_population_split(&msp, 10, 1, pops, 0), MSP_ERR_SOURCE_DEST_EQUAL);
        pops[0] = -1;
        CU_ASSERT_EQUAL(msp_add_population_split(&msp, 10, 1, pops, 0),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        pops[0] = 0;
        pops[1] = 0;
        CU_ASSERT_EQUAL(msp_add_population_split(&msp, 10, 2, pops, 1),
            MSP_ERR_DUPLICATE_POPULATION);
        CU_ASSERT_EQUAL(
            msp_add_population_split(&msp, -1, 1, pops, 1), MSP_ERR_BAD_PARAM_VALUE);

        pops[0] = 0;
        CU_ASSERT_EQUAL(msp_add_admixture(&msp, 10.0, 1, 100, pops, probs),
            MSP_ERR_TOO_MANY_EVENT_POPULATIONS);
        CU_ASSERT_EQUAL(msp_add_admixture(&msp, 10.0, 2, 1, pops, probs),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        CU_ASSERT_EQUAL(msp_add_admixture(&msp, 10.0, -1, 1, pops, probs),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        CU_ASSERT_EQUAL_FATAL(
            msp_add_admixture(&msp, 10, 0, 1, pops, probs), MSP_ERR_SOURCE_DEST_EQUAL);
        pops[0] = -1;
        CU_ASSERT_EQUAL_FATAL(msp_add_admixture(&msp, 10, 0, 1, pops, probs),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        pops[0] = 2;
        CU_ASSERT_EQUAL_FATAL(msp_add_admixture(&msp, 10, 0, 1, pops, probs),
            MSP_ERR_POPULATION_OUT_OF_BOUNDS);
        pops[0] = 1;
        pops[1] = 1;
        CU_ASSERT_EQUAL_FATAL(msp_add_admixture(&msp, 10, 0, 2, pops, probs),
            MSP_ERR_DUPLICATE_POPULATION);
        probs[0] = -1;
        CU_ASSERT_EQUAL_FATAL(
            msp_add_admixture(&msp, 10, 0, 1, pops, probs), MSP_ERR_BAD_PROPORTION);
        probs[0] = 1.1;
        CU_ASSERT_EQUAL_FATAL(
            msp_add_admixture(&msp, 10, 0, 1, pops, probs), MSP_ERR_BAD_PROPORTION);
        probs[0] = 1;
        CU_ASSERT_EQUAL_FATAL(
            msp_add_admixture(&msp, -1, 0, 1, pops, probs), MSP_ERR_BAD_PARAM_VALUE);

        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, -1, 0, 0.2),
            MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, 0, -1, 0.2),
            MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, 0, 2, 0.2),
            MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, 2, 0, 0.2),
            MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, -1, 0, 0.2),
            MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, 0, -1, 0.2),
            MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, -1, -1, -0.2),
            MSP_ERR_BAD_PARAM_VALUE);
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 10, 0, 0, 0.2),
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
        ret = msp_add_migration_rate_change(&msp, 0.2, 0, 1, 0.2);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_add_migration_rate_change(&msp, 0.3, -1, -1, 0.3);
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
        CU_ASSERT_EQUAL(msp_add_migration_rate_change(&msp, 0.2, 0, 1, 0.2),
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
        tsk_table_collection_free(&tables);
    }

    free(samples);
    gsl_rng_free(rng);
}

static void
test_demographic_events_start_time(void)
{
    int ret;
    uint32_t n = 10;
    uint32_t model;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    for (model = 0; model < 2; model++) {
        ret = build_sim(&msp, &tables, rng, 2, 1, NULL, n);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_recombination_rate(&msp, 1);
        CU_ASSERT_EQUAL(ret, 0);

        if (model == 0) {
            ret = msp_set_start_time(&msp, 1.0);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_add_population_parameters_change(&msp, 0.4, 0, 0.5, 1.0);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_set_simulation_model_hudson(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            /* Times need to be bumped to integer values for DTWF */
            ret = msp_set_start_time(&msp, 2.0);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_add_population_parameters_change(&msp, 1, 0, 0.5, 1.0);
            CU_ASSERT_EQUAL(ret, 0);
            ret = msp_set_simulation_model_dtwf(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        }

        ret = msp_initialise(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        msp_run(&msp, DBL_MAX, UINT32_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_free(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        tsk_table_collection_free(&tables);
    }

    gsl_rng_free(rng);
}

static void
test_population_split(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = {
        { .time = 0, .population = 0 },
        { .time = 0, .population = 0 },
        { .time = 0, .population = 1 },
        { .time = 0, .population = 1 },
        { .time = 0, .population = 2 },
        { .time = 0, .population = 2 },
    };
    int source[] = { 0, 1, 2 };
    tsk_table_collection_t tables;

    ret = build_sim(&msp, &tables, rng, 1, 4, samples, 6);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_split(&msp, 0.5, 3, source, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[3].state, MSP_POP_STATE_ACTIVE);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);
    msp_print_state(&msp, _devnull);

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[3].state, MSP_POP_STATE_ACTIVE);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_population_split_debug(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    int32_t derived[2];
    double time;
    /*
            6
          ┏━┻━┓
          ┃   5
          ┃ ┏━┻┓
          ┃ ┃  4
          ┃ ┃ ┏┻┓
          0 1 2 3
    */
    ret = build_sim(&msp, &tables, rng, 1, 7, NULL, 0);
    CU_ASSERT_EQUAL(ret, 0);

    derived[0] = 2;
    derived[1] = 3;
    ret = msp_add_population_split(&msp, 4, 2, derived, 4);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 4, 1, 0, false);
    CU_ASSERT_EQUAL(ret, 0);
    derived[0] = 1;
    derived[1] = 4;
    ret = msp_add_population_split(&msp, 5, 2, derived, 5);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 5, 1, 0, false);
    CU_ASSERT_EQUAL(ret, 0);
    derived[0] = 0;
    derived[1] = 5;
    ret = msp_add_population_split(&msp, 6, 2, derived, 6);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 6, 1, 0, false);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[3].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[4].state, MSP_POP_STATE_INACTIVE);
    CU_ASSERT_EQUAL(msp.populations[5].state, MSP_POP_STATE_INACTIVE);
    CU_ASSERT_EQUAL(msp.populations[6].state, MSP_POP_STATE_INACTIVE);

    ret = msp_debug_demography(&msp, &time);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(time, 4);

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[3].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[4].state, MSP_POP_STATE_INACTIVE);
    CU_ASSERT_EQUAL(msp.populations[5].state, MSP_POP_STATE_INACTIVE);
    CU_ASSERT_EQUAL(msp.populations[6].state, MSP_POP_STATE_INACTIVE);

    ret = msp_debug_demography(&msp, &time);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(time, 5);

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[3].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[4].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[5].state, MSP_POP_STATE_INACTIVE);
    CU_ASSERT_EQUAL(msp.populations[6].state, MSP_POP_STATE_INACTIVE);

    ret = msp_debug_demography(&msp, &time);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(time, 6);

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[3].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[4].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[5].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[6].state, MSP_POP_STATE_INACTIVE);

    ret = msp_debug_demography(&msp, &time);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(isinf(time));

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[3].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[4].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[5].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[6].state, MSP_POP_STATE_ACTIVE);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_population_split_replicates(void)
{
    int ret;
    int j;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    int32_t derived[] = { 0, 1 };
    sample_t samples[]
        = { { .time = 0, .population = 0 }, { .time = 0, .population = 1 } };
    /*
        2
       ┏┻┓
       0 1
    */
    /* Nominal case. */
    ret = build_sim(&msp, &tables, rng, 1, 3, samples, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_split(&msp, 2, 2, derived, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 2, 1, 0, false);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    for (j = 0; j < 10; j++) {
        ret = msp_reset(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_ACTIVE);
        CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_ACTIVE);
        CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_INACTIVE);

        ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
        CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
        CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_ACTIVE);
    }

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_inactive_populations(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    int32_t derived[] = { 0, 1 };
    sample_t samples[]
        = { { .time = 0, .population = 0 }, { .time = 0, .population = 1 } };
    /*
        2
       ┏┻┓
       0 1
    */
    /* Nominal case. */
    ret = build_sim(&msp, &tables, rng, 1, 3, samples, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_split(&msp, 2, 2, derived, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 2, 1, 0, false);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* Sample from inactive at time zero. */
    samples[0].population = 2;
    ret = build_sim(&msp, &tables, rng, 1, 3, samples, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_split(&msp, 2, 2, derived, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 2, 1, 0, false);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, MSP_ERR_POPULATION_INACTIVE_SAMPLE);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* Sample from inactive as event. */
    samples[0].population = 2;
    samples[0].time = 1;
    ret = build_sim(&msp, &tables, rng, 1, 3, samples, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_split(&msp, 2, 2, derived, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 2, 1, 0, false);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_POPULATION_INACTIVE_SAMPLE);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* Migrate into inactive via mass migration */
    samples[0].population = 0;
    samples[0].time = 0;
    ret = build_sim(&msp, &tables, rng, 1, 3, samples, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_split(&msp, 2, 2, derived, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_mass_migration(&msp, 2.0001, 2, 1, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_POPULATION_INACTIVE_MOVE);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* Migrate into inactive via continuous migration */
    ret = build_sim(&msp, &tables, rng, 1, 3, samples, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_split(&msp, 2, 2, derived, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_migration_rate_change(&msp, 2.0001, 2, 1, 10.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_POPULATION_INACTIVE_MOVE);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    gsl_rng_free(rng);
}

static void
test_admixture(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = {
        { .time = 0, .population = 3 },
        { .time = 0, .population = 3 },
        { .time = 0, .population = 3 },
        { .time = 0, .population = 3 },
        { .time = 0, .population = 3 },
        { .time = 0, .population = 3 },
    };
    int admix_source[] = { 0, 1, 2 };
    double probs[] = { 0.25, 0.25, 0.5 };
    tsk_table_collection_t tables;

    /* Build a 5 pop model, where 3 is admixed from 0, 1, 2, and then
     * 0, 1, 2 all split from 4. */
    ret = build_sim(&msp, &tables, rng, 1, 5, samples, 6);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_admixture(&msp, 0.5, 3, 3, admix_source, probs);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_split(&msp, 1.5, 3, admix_source, 4);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 4, 1, 0, false);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[3].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[4].state, MSP_POP_STATE_INACTIVE);

    ret = msp_run(&msp, 0.51, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_TIME);
    msp_verify(&msp, 0);
    msp_print_state(&msp, _devnull);

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[3].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[4].state, MSP_POP_STATE_INACTIVE);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);
    msp_print_state(&msp, _devnull);

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[3].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[4].state, MSP_POP_STATE_ACTIVE);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_admixture_excess_probability(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = {
        { .time = 0, .population = 3 },
        { .time = 0, .population = 3 },
        { .time = 0, .population = 3 },
        { .time = 0, .population = 3 },
        { .time = 0, .population = 3 },
        { .time = 0, .population = 1 },
    };
    int admix_source[] = { 0, 1, 2 };
    /* Put all the weight on population 1, but don't make the total sum = 1 */
    double probs[] = { 0.00001, 1, 0.5 };
    tsk_table_collection_t tables;

    /* Build a 4 pop model, where 3 is admixed from 0, 1, 2. Since all
     * the probability is from pop 1, the lineages should all migrate there
     * and coalesce */
    ret = build_sim(&msp, &tables, rng, 1, 4, samples, 6);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_admixture(&msp, 0.5, 3, 3, admix_source, probs);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[3].state, MSP_POP_STATE_ACTIVE);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);
    msp_print_state(&msp, _devnull);

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[3].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_population_state_machine_errors(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    double time;
    double proba[] = { 1 };
    int32_t derived[] = { 0 };
    int32_t ancestral[] = { 0 };

    /* PREVIOUSLY_ACTIVE -> ACTIVE */
    ret = build_sim(&msp, &tables, rng, 1, 3, NULL, 0);
    CU_ASSERT_EQUAL(ret, 0);
    /* 0 -> 1 */
    derived[0] = 0;
    ret = msp_add_population_split(&msp, 2, 1, derived, 1);
    CU_ASSERT_EQUAL(ret, 0);
    /* 1 -> 2 */
    derived[0] = 1;
    ret = msp_add_population_split(&msp, 2, 1, derived, 2);
    CU_ASSERT_EQUAL(ret, 0);
    /* 2 -> 1 */
    derived[0] = 2;
    ret = msp_add_population_split(&msp, 2, 1, derived, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_debug_demography(&msp, &time);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(time, 2);
    ret = msp_debug_demography(&msp, &time);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_POPULATION_PREVIOUSLY_ACTIVE);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* Inactive population as admixture derived.*/
    ret = build_sim(&msp, &tables, rng, 1, 2, NULL, 0);
    CU_ASSERT_EQUAL(ret, 0);
    /* 0 -> 1 */
    ancestral[0] = 1;
    ret = msp_add_admixture(&msp, 2, 0, 1, ancestral, proba);
    CU_ASSERT_EQUAL(ret, 0);
    /* 1 -> 0 */
    derived[0] = 1;
    ret = msp_add_population_split(&msp, 2, 1, derived, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 1, 0, false);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_debug_demography(&msp, &time);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(time, 2);
    ret = msp_debug_demography(&msp, &time);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_ADMIX_DERIVED_NOT_ACTIVE);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* derived not active in split */
    ret = build_sim(&msp, &tables, rng, 1, 2, NULL, 0);
    CU_ASSERT_EQUAL(ret, 0);
    /* 0 -> 1 */
    derived[0] = 0;
    ret = msp_add_population_split(&msp, 2, 1, derived, 1);
    CU_ASSERT_EQUAL(ret, 0);
    /* 1 -> 0 */
    derived[0] = 1;
    ret = msp_add_population_split(&msp, 2, 1, derived, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 1, 0, false);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_debug_demography(&msp, &time);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(time, 2);
    ret = msp_debug_demography(&msp, &time);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_SPLIT_DERIVED_NOT_ACTIVE);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    gsl_rng_free(rng);
}

static void
test_census_event(void)
{
    int ret;
    uint32_t n = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    int num_census_nodes = 0;
    int i;

    ret = build_sim(&msp, &tables, rng, 2, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_recombination_rate(&msp, 1);

    /* Add a census event in at 0.5 generations. */
    ret = msp_add_census_event(&msp, 0.5);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);
    msp_print_state(&msp, _devnull);

    /* Check there is more than 1 node at the census time. */
    for (i = 0; i < tables.nodes.num_rows; i++) {
        if (tables.nodes.time[i] == 0.5) {
            num_census_nodes++;
        }
    }
    CU_ASSERT_TRUE(num_census_nodes > 1);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_activate_population_event(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = {
        { .time = 2, .population = 0 },
        { .time = 2, .population = 0 },
        { .time = 0, .population = 1 },
        { .time = 0, .population = 1 },
    };
    tsk_table_collection_t tables;

    ret = build_sim(&msp, &tables, rng, 1, 3, samples, 4);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_set_population_configuration(&msp, 0, 1, 0, false);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_add_activate_population_event(&msp, 1, 0);
    CU_ASSERT_EQUAL(ret, 0);

    int derived[] = { 0, 1 };
    ret = msp_add_population_split(&msp, 4, 2, derived, 2);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_INACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_ACTIVE);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);
    msp_print_state(&msp, _devnull);

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[2].state, MSP_POP_STATE_ACTIVE);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_activate_population_event_errors(void)
{
    int ret;
    uint32_t n = 2;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    /* Activate already active population. */
    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_activate_population_event(&msp, 1, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_POPULATION_CURRENTLY_ACTIVE);

    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* Activate nonexistent populations. */
    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_activate_population_event(&msp, 1, -1);
    CU_ASSERT_EQUAL(ret, MSP_ERR_POPULATION_OUT_OF_BOUNDS);

    msp_free(&msp);
    tsk_table_collection_free(&tables);

    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_activate_population_event(&msp, 1, 1);
    CU_ASSERT_EQUAL(ret, MSP_ERR_POPULATION_OUT_OF_BOUNDS);

    msp_free(&msp);
    tsk_table_collection_free(&tables);

    gsl_rng_free(rng);
}

static void
test_deactivate_population_event(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = {
        { .time = 0, .population = 0 },
        { .time = 0, .population = 0 },
    };
    tsk_table_collection_t tables;

    ret = build_sim(&msp, &tables, rng, 1, 2, samples, 2);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_add_deactivate_population_event(&msp, 1, 1);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_ACTIVE);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);
    msp_print_state(&msp, _devnull);

    CU_ASSERT_EQUAL(msp.populations[0].state, MSP_POP_STATE_ACTIVE);
    CU_ASSERT_EQUAL(msp.populations[1].state, MSP_POP_STATE_PREVIOUSLY_ACTIVE);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_deactivate_population_event_errors(void)
{
    int ret;
    uint32_t n = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    /* Deactivate inactive population. */
    ret = build_sim(&msp, &tables, rng, 1, 2, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 1, 1, 0, false);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_deactivate_population_event(&msp, 1, 1);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_DEACTIVATE_INACTIVE_POPULATION);

    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* Deactivate previously active population. */
    ret = build_sim(&msp, &tables, rng, 1, 2, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_deactivate_population_event(&msp, 0.1, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_deactivate_population_event(&msp, 0.2, 1); /* deactivate again */
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_DEACTIVATE_PREVIOUSLY_ACTIVE_POPULATION);

    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* Activate nonexistent populations. */
    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_deactivate_population_event(&msp, 1, -1);
    CU_ASSERT_EQUAL(ret, MSP_ERR_POPULATION_OUT_OF_BOUNDS);

    msp_free(&msp);
    tsk_table_collection_free(&tables);

    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_deactivate_population_event(&msp, 1, 1);
    CU_ASSERT_EQUAL(ret, MSP_ERR_POPULATION_OUT_OF_BOUNDS);

    msp_free(&msp);
    tsk_table_collection_free(&tables);

    gsl_rng_free(rng);
}

static void
test_time_travel_error(void)
{
    int ret;
    uint32_t n = 100;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    msp_t msp;

    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_add_simple_bottleneck(&msp, 0.1, 0, 0.75);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_simple_bottleneck(&msp, 0.1, 0, 1.0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_TIME_TRAVEL);

    msp_free(&msp);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_floating_point_extremes(void)
{
    int ret;
    uint32_t n = 2;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    double denormal_min = pow(2, -52) * pow(2, -1022);

    ret = build_sim(&msp, &tables, rng, denormal_min, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, DBL_MAX), 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_discrete_genome(&msp, false), 0);
    ret = msp_set_population_configuration(&msp, 0, DBL_MAX, 0, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BREAKPOINT_RESAMPLE_OVERFLOW);
    msp_print_state(&msp, _devnull);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* A long sequence length and high recombination should overflow */
    ret = build_sim(&msp, &tables, rng, DBL_MAX, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_discrete_genome(&msp, false), 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, DBL_MAX), 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BREAKPOINT_MASS_NON_FINITE);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    gsl_rng_free(rng);
}

static void
test_simulation_replicates(void)
{
    int ret;
    uint32_t n = 100;
    uint32_t m = 10;
    double mutation_rate = 2;
    size_t num_replicates = 10;
    long seed = 10;
    double migration_matrix[] = { 0, 1, 1, 0 };
    size_t j;
    gsl_rng *rng = safe_rng_alloc();
    msp_t msp;
    tsk_treeseq_t ts;
    mutgen_t mutgen;
    mutation_model_t mut_model;
    tsk_table_collection_t tables;

    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_BINARY);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    gsl_rng_set(rng, seed);
    ret = build_sim(&msp, &tables, rng, m, 2, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 1), 0);
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
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, mutation_rate);
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
        ret = mutgen_generate(&mutgen, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
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
    tsk_table_collection_free(&tables);
}

static void
test_bottleneck_simulation(void)
{
    int ret;
    uint32_t n = 100;
    uint32_t m = 10;
    long seed = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    double t1 = 0.1;
    double t2 = 0.5;
    int t1_found = 0;
    tsk_table_collection_t tables;

    gsl_rng_set(rng, seed);
    ret = build_sim(&msp, &tables, rng, m, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 1.0 / m), 0);
    /* set all the block sizes to something small to provoke the memory
     * expansions. */
    ret = msp_set_avl_node_block_size(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_node_mapping_block_size(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_segment_block_size(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    /* Add a bottleneck that does nothing at time t1 */
    ret = msp_add_simple_bottleneck(&msp, t1, 0, 0);
    CU_ASSERT_EQUAL(ret, 0);
    /* Add a bottleneck that coalesces everything at t2 */
    ret = msp_add_simple_bottleneck(&msp, t2, 0, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(&msp, t1, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 2);
    CU_ASSERT_FALSE(msp_is_completed(&msp));
    msp_print_state(&msp, _devnull);

    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(msp_is_completed(&msp));
    CU_ASSERT_EQUAL(msp.time, t2);
    msp_verify(&msp, 0);

    msp_reset(&msp);
    while ((ret = msp_run(&msp, DBL_MAX, 1)) == 1) {
        msp_verify(&msp, 0);
        if (msp.time == 0.1) {
            t1_found = 1;
        }
        CU_ASSERT_EQUAL(msp_get_time(&msp), msp.time);
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(t1_found);
    CU_ASSERT_EQUAL(msp.time, t2);

    msp_free(&msp);
    gsl_rng_free(rng);
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
    long seed = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    uint32_t num_bottlenecks = 10;
    bottleneck_desc_t bottlenecks[num_bottlenecks];
    double t;
    tsk_table_collection_t tables;

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
    ret = build_sim(&msp, &tables, rng, 10, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 1), 0);
    for (j = 0; j < num_bottlenecks; j++) {
        ret = msp_add_simple_bottleneck(
            &msp, bottlenecks[j].time, 0, bottlenecks[j].parameter);
        CU_ASSERT_EQUAL(ret, 0);
    }
    ret = msp_set_population_configuration(&msp, 0, 0.25, 0, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    for (j = 0; j < num_bottlenecks - 1; j++) {
        ret = msp_run(&msp, bottlenecks[j].time + 1e-6, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_TIME);
        CU_ASSERT_FALSE(msp_is_completed(&msp));
        CU_ASSERT_DOUBLE_EQUAL(msp.time, bottlenecks[j].time + 1e-6, 1e-9);
        msp_verify(&msp, 0);
    }
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(msp_is_completed(&msp));
    CU_ASSERT_EQUAL(msp.time, bottlenecks[num_bottlenecks - 1].time);
    msp_verify(&msp, 0);

    /* Test out resets on partially completed simulations. */
    ret = msp_reset(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_bottlenecks - 1; j++) {
        ret = msp_run(&msp, bottlenecks[j].time, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 2);
    }
    ret = msp_reset(&msp);
    msp_verify(&msp, 0);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
verify_simulate_from(int model, rate_map_t *recomb_map,
    tsk_table_collection_t *from_tables, size_t num_replicates)
{
    int ret;
    size_t j;
    tsk_bookmark_t pos;
    tsk_table_collection_t tables;
    tsk_treeseq_t final;
    tsk_tree_t tree;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();

    ret = tsk_table_collection_copy(from_tables, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_alloc(&msp, &tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_recombination_map(
        &msp, recomb_map->size, recomb_map->position, recomb_map->rate);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    if (model == MSP_MODEL_DTWF) {
        ret = msp_set_simulation_model_dtwf(&msp);
        CU_ASSERT_EQUAL(ret, 0);
    }
    ret = msp_set_population_configuration(&msp, 0, 10, 0, true);
    CU_ASSERT_EQUAL(ret, 0);
    /* TODO add dirac and other models */
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < num_replicates; j++) {
        msp_verify(&msp, MSP_VERIFY_BREAKPOINTS);
        /* Run for a little bit and check the state again */
        ret = msp_run(&msp, DBL_MAX, 10);
        CU_ASSERT_FATAL(ret >= 0);
        msp_verify(&msp, MSP_VERIFY_BREAKPOINTS);
        ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_TRUE(msp_is_completed(&msp));
        ret = msp_finalise_tables(&msp);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_treeseq_init(&final, &tables, TSK_TS_INIT_BUILD_INDEXES);
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
        CU_ASSERT_TRUE(tsk_table_collection_equals(from_tables, &tables, 0));

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
    msp_t *msp_source, rate_map_t *recomb_map, tsk_table_collection_t *from_tables)
{
    int ret;
    msp_t msp_dest;
    size_t num_ancestors;
    gsl_rng *rng = safe_rng_alloc();

    ret = msp_alloc(&msp_dest, from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_recombination_map(
        &msp_dest, recomb_map->size, recomb_map->position, recomb_map->rate);
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
verify_simple_simulate_from(int model, uint32_t n, double sequence_length,
    double recombination_rate, size_t num_events, size_t num_replicates)
{
    int ret;
    tsk_table_collection_t tables;
    msp_t msp;
    rate_map_t recomb_map;
    gsl_rng *rng = safe_rng_alloc();

    ret = rate_map_alloc_single(&recomb_map, sequence_length, recombination_rate);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = build_sim(&msp, &tables, rng, sequence_length, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_recombination_map(
        &msp, recomb_map.size, recomb_map.position, recomb_map.rate);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
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
    rate_map_free(&recomb_map);
}

static void
test_simulate_from_single_locus(void)
{
    verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 1.0, 0, 5, 1);
    verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 1.0, 0, 5, 1);

    verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 1.0, 1, 5, 1);
    verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 1.0, 1, 5, 1);
}

static void
test_simulate_from_single_locus_replicates(void)
{
    verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 1.0, 0, 5, 10);
    verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 1.0, 0, 5, 10);
}

static void
test_simulate_from_multilocus(void)
{
    verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 10.0, 1, 5, 1);
    verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 1.0, 0, 5, 1);
    verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 1.0, 1, 5, 1);
    verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 1.0, 1, 5, 1);
}

static void
test_simulate_from_large_example(void)
{
    verify_simple_simulate_from(MSP_MODEL_HUDSON, 100, 100.0, 1, 500, 2);
}

static void
test_simulate_from_empty(void)
{
    verify_simple_simulate_from(MSP_MODEL_HUDSON, 10, 1.0, 0, 0, 1);
    verify_simple_simulate_from(MSP_MODEL_DTWF, 10, 1.0, 0, 0, 1);
}

static void
test_simulate_from_ancient_samples(void)
{
    int ret;
    uint32_t n = 25;
    msp_t msp;
    tsk_table_collection_t tables;
    rate_map_t recomb_map;
    gsl_rng *rng = safe_rng_alloc();

    /* We don't simulate any recombination here because if we
     * do there's a chance we can get malformed tables in the
     * output. This happens because edges can be in not-quite
     * the right order in the output when we have these unsorted
     * edges. It's not a supported input type anyway, so maybe
     * not worth worrying about. */
    ret = rate_map_alloc_single(&recomb_map, 100.0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = build_sim(&msp, &tables, rng, 100, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_recombination_rate(&msp, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, 10);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_EXIT_MAX_EVENTS);
    CU_ASSERT_FALSE(msp_is_completed(&msp));
    /* We don't finalise the tables so we don't have all the
     * unary edges. We should therefore get lots of roots
     * at different times. */
    msp_free(&msp);

    verify_simulate_from(MSP_MODEL_HUDSON, &recomb_map, &tables, 1);

    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
    rate_map_free(&recomb_map);
}

static void
test_simulate_from_completed(void)
{
    int ret;
    uint32_t n = 25;
    msp_t msp;
    rate_map_t recomb_map;
    tsk_table_collection_t tables;
    gsl_rng *rng = safe_rng_alloc();
    double recombination_rate = 2;

    ret = rate_map_alloc_single(&recomb_map, 1.0, recombination_rate);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_recombination_map(
        &msp, recomb_map.size, recomb_map.position, recomb_map.rate);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
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
    rate_map_free(&recomb_map);
}

static void
test_simulate_from_incompatible(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t from_tables;

    ret = tsk_table_collection_init(&from_tables, 0);

    CU_ASSERT_EQUAL_FATAL(ret, 0);
    from_tables.sequence_length = 1.0;
    /* Add two nodes so that we meet the minimum requirements */
    ret = tsk_population_table_add_row(&from_tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_add_row(
        &from_tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_add_row(
        &from_tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 1);

    /* Malformed tree sequence */
    from_tables.nodes.individual[0] = 100;
    ret = msp_alloc(&msp, &from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_FATAL(msp_is_tsk_error(ret));
    CU_ASSERT_EQUAL_FATAL(
        ret ^ (1 << MSP_TSK_ERR_BIT), TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    CU_ASSERT_STRING_EQUAL(
        msp_strerror(ret), tsk_strerror(TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS));
    from_tables.nodes.individual[0] = -1;
    msp_free(&msp);

    /* older nodes */
    from_tables.sequence_length = 10.0;
    ret = tsk_node_table_add_row(&from_tables.nodes, 0, 1.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&from_tables.nodes, 0, 2.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    /* Make a tree sequence that we cannot recover trees from */
    ret = tsk_edge_table_add_row(&from_tables.edges, 0, 1, 2, 0, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    tsk_edge_table_add_row(&from_tables.edges, 0, 1, 3, 0, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = msp_alloc(&msp, &from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_FATAL(msp_is_tsk_error(ret));
    CU_ASSERT_EQUAL_FATAL(
        ret ^ (1 << MSP_TSK_ERR_BIT), TSK_ERR_BAD_EDGES_CONTRADICTORY_CHILDREN);
    CU_ASSERT_STRING_EQUAL(
        msp_strerror(ret), tsk_strerror(TSK_ERR_BAD_EDGES_CONTRADICTORY_CHILDREN));
    msp_free(&msp);
    tsk_edge_table_clear(&from_tables.edges);

    tsk_table_collection_clear(&from_tables, 0);
    from_tables.sequence_length = 0.0;
    ret = msp_alloc(&msp, &from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_SEQUENCE_LENGTH);
    msp_free(&msp);

    from_tables.sequence_length = 1.0;
    ret = msp_alloc(&msp, &from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_ZERO_POPULATIONS);
    msp_free(&msp);

    tsk_population_table_add_row(&from_tables.populations, NULL, 0);
    ret = msp_alloc(&msp, &from_tables, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    msp_free(&msp);

    gsl_rng_free(rng);
    tsk_table_collection_free(&from_tables);
}

static void
test_simulate_empty_tables(void)
{
    int ret;
    tsk_table_collection_t tables;
    gsl_rng *rng = safe_rng_alloc();
    msp_t msp;

    ret = build_sim(&msp, &tables, rng, 1, 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_initialise(&msp), 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 0);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 0);
    CU_ASSERT_EQUAL_FATAL(tables.populations.num_rows, 1);
    msp_print_state(&msp, _devnull);

    msp_free(&msp);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_simulate_init_errors(void)
{
    int ret;
    uint32_t n = 25;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    ret = msp_alloc(&msp, NULL, rng);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);

    ret = msp_alloc(&msp, &tables, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);

    /* No populations */
    ret = build_sim(&msp, &tables, rng, 1, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_ZERO_POPULATIONS);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    ret = build_sim(&msp, &tables, rng, -1, 1, NULL, n);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_SEQUENCE_LENGTH);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    gsl_rng_free(rng);
}

static void
check_zero_population_size(
    sample_t *samples, size_t n, double N, double T, int init_ret, int run_ret)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    int model;

    for (model = 0; model < 2; model++) {
        ret = build_sim(&msp, &tables, rng, 1, 3, samples, n);
        CU_ASSERT_EQUAL(ret, 0);

        if (model == 0) {
            ret = msp_set_simulation_model_hudson(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            ret = msp_set_simulation_model_dtwf(&msp);
            CU_ASSERT_EQUAL(ret, 0);
        }

        ret = msp_set_population_configuration(&msp, 0, 0, 0, true);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_population_configuration(&msp, 1, N, 0, true);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_population_configuration(&msp, 2, N, 0, true);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_add_mass_migration(&msp, T, 1, 0, 1.0);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_add_mass_migration(&msp, T, 2, 0, 1.0);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_add_population_parameters_change(&msp, T, 0, N, 0);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_add_population_parameters_change(&msp, T, 1, 0, 0);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_add_population_parameters_change(&msp, T, 2, 0, 0);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_initialise(&msp);
        CU_ASSERT_EQUAL_FATAL(ret, init_ret);

        if (init_ret == 0) {
            msp_print_state(&msp, _devnull);
            ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
            CU_ASSERT_EQUAL_FATAL(ret, run_ret);
            if (ret == 0) {
                msp_verify(&msp, 0);
                msp_print_state(&msp, _devnull);
            }
        }
        CU_ASSERT_EQUAL(msp_free(&msp), 0);
        tsk_table_collection_free(&tables);
    }

    gsl_rng_free(rng);
}

static void
test_zero_population_size(void)
{
    const size_t n = 2;
    const double N = 1000;
    const double T = 500;
    sample_t samples_good_1[] = { { 1, 0.0 }, { 1, 0.0 } };
    sample_t samples_good_2[] = { { 1, 0.0 }, { 2, 0.0 } };
    sample_t samples_good_3[] = { { 1, 0.0 }, { 0, T } };
    sample_t samples_bad_1[] = { { 0, 0.0 }, { 1, 0.0 } };
    sample_t samples_bad_2[] = { { 1, T }, { 1, 0.0 } };

    check_zero_population_size(samples_good_1, n, N, T, 0, 0);
    check_zero_population_size(samples_good_2, n, N, T, 0, 0);
    check_zero_population_size(samples_good_3, n, N, T, 0, 0);
    check_zero_population_size(samples_bad_1, n, N, T, MSP_ERR_POP_SIZE_ZERO_SAMPLE, 0);
    check_zero_population_size(samples_bad_2, n, N, T, 0, MSP_ERR_POP_SIZE_ZERO_SAMPLE);
}

static void
test_population_size_start_time(void)
{

    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    double t = 10;
    double size1, size2;

    ret = build_sim(&msp, &tables, rng, 1, 2, NULL, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 1, 0.1, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_compute_population_size(&msp, 0, t, &size1);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* The size of the population should not be affected by start_time */
    ret = build_sim(&msp, &tables, rng, 1, 2, NULL, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 1, 0.1, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_start_time(&msp, 5);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_compute_population_size(&msp, 0, t, &size2);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    CU_ASSERT_EQUAL_FATAL(size1, size2);
    gsl_rng_free(rng);
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_single_locus_simulation", test_single_locus_simulation },
        { "test_single_locus_two_populations", test_single_locus_two_populations },
        { "test_single_locus_many_populations", test_single_locus_many_populations },
        { "test_single_locus_labels", test_single_locus_labels },
        { "test_single_locus_historical_sample", test_single_locus_historical_sample },
        { "test_single_locus_all_historical", test_single_locus_all_historical },
        { "test_single_locus_multiple_historical_samples",
            test_single_locus_multiple_historical_samples },
        { "test_single_locus_historical_sample_start_time",
            test_single_locus_historical_sample_start_time },
        { "test_single_locus_historical_sample_end_time",
            test_single_locus_historical_sample_end_time },

        { "test_multi_locus_simulation", test_multi_locus_simulation },
        { "test_multi_locus_bottleneck_arg", test_multi_locus_bottleneck_arg },

        { "test_dtwf_single_locus_simulation", test_dtwf_single_locus_simulation },
        { "test_dtwf_multi_locus_simulation", test_dtwf_multi_locus_simulation },
        { "test_dtwf_deterministic", test_dtwf_deterministic },
        { "test_dtwf_simultaneous_historical_samples",
            test_dtwf_simultaneous_historical_samples },
        { "test_dtwf_low_recombination", test_dtwf_low_recombination },
        { "test_dtwf_events_between_generations", test_dtwf_events_between_generations },
        { "test_dtwf_unsupported_bottleneck", test_dtwf_unsupported_bottleneck },
        { "test_dtwf_zero_pop_size", test_dtwf_zero_pop_size },
        { "test_dtwf_migration_matrix_not_stochastic",
            test_dtwf_migration_matrix_not_stochastic },
        { "test_dtwf_errors", test_dtwf_errors },

        { "test_mixed_hudson_dtwf", test_mixed_hudson_dtwf },
        { "test_mixed_hudson_smc", test_mixed_hudson_smc },

        { "test_gc_single_locus", test_gc_single_locus },
        { "test_gc_tract_lengths", test_gc_tract_lengths },
        { "test_gc_zero_recombination", test_gc_zero_recombination },
        { "test_gc_rates", test_gc_rates },

        { "test_multiple_mergers_simulation", test_multiple_mergers_simulation },
        { "test_multiple_mergers_growth_rate", test_multiple_mergers_growth_rate },
        { "test_dirac_coalescent_bad_parameters", test_dirac_coalescent_bad_parameters },
        { "test_beta_coalescent_bad_parameters", test_beta_coalescent_bad_parameters },

        { "test_simulator_getters_setters", test_simulator_getters_setters },
        { "test_demographic_events", test_demographic_events },
        { "test_demographic_events_start_time", test_demographic_events_start_time },
        { "test_population_split", test_population_split },
        { "test_population_split_debug", test_population_split_debug },
        { "test_population_split_replicates", test_population_split_replicates },
        { "test_inactive_populations", test_inactive_populations },
        { "test_admixture", test_admixture },
        { "test_admixture_excess_probability", test_admixture_excess_probability },
        { "test_population_state_machine_errors", test_population_state_machine_errors },
        { "test_census_event", test_census_event },
        { "test_activate_population_event", test_activate_population_event },
        { "test_activate_population_event_errors",
            test_activate_population_event_errors },
        { "test_deactivate_population_event", test_deactivate_population_event },
        { "test_deactivate_population_event_errors",
            test_deactivate_population_event_errors },
        { "test_time_travel_error", test_time_travel_error },
        { "test_floating_point_extremes", test_floating_point_extremes },
        { "test_simulation_replicates", test_simulation_replicates },
        { "test_bottleneck_simulation", test_bottleneck_simulation },
        { "test_large_bottleneck_simulation", test_large_bottleneck_simulation },

        { "test_simulate_from_single_locus", test_simulate_from_single_locus },
        { "test_simulate_from_single_locus_replicates",
            test_simulate_from_single_locus_replicates },
        { "test_simulate_from_multilocus", test_simulate_from_multilocus },
        { "test_simulate_from_large_example", test_simulate_from_large_example },
        { "test_simulate_from_empty", test_simulate_from_empty },
        { "test_simulate_from_ancient_samples", test_simulate_from_ancient_samples },
        { "test_simulate_from_completed", test_simulate_from_completed },
        { "test_simulate_from_incompatible", test_simulate_from_incompatible },
        { "test_simulate_empty_tables", test_simulate_empty_tables },

        { "test_simulate_init_errors", test_simulate_init_errors },
        { "test_zero_population_size", test_zero_population_size },
        { "test_population_size_start_time", test_population_size_start_time },

        CU_TEST_INFO_NULL,
    };

    return test_main(tests, argc, argv);
}
