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
#include "msprime.c"

static void
test_pedigree_trio(void)
{
    int ret;
    tsk_table_collection_t tables;
    int num_inds = 3;
    int ploidy = 2;
    tsk_id_t parents[] = { -1, -1, -1, -1, 0, 1 }; // size num_inds * ploidy
    double time[] = { 1, 1, 0 };
    tsk_flags_t is_sample[] = { 0, 0, 1 };
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();

    ret = build_pedigree_sim(
        &msp, &tables, rng, 1, ploidy, num_inds, parents, time, is_sample);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);

    // msp_print_pedigree_inds(&msp, stdout);
    tsk_node_table_print_state(&msp.tables->nodes, stdout);

    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_pedigree_three_generations(void)
{
    int ret;
    tsk_table_collection_t tables;
    int num_inds = 7;
    int ploidy = 2;
    tsk_id_t parents[]
        = { -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5 }; // size num_inds * ploidy
    double time[] = { 2, 2, 2, 2, 1, 1, 0 };
    tsk_flags_t is_sample[] = { 0, 0, 0, 0, 0, 0, 1 };
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();

    ret = build_pedigree_sim(
        &msp, &tables, rng, 1, ploidy, num_inds, parents, time, is_sample);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);

    msp_print_pedigree_inds(&msp, stdout);
    tsk_node_table_print_state(&msp.tables->nodes, stdout);

    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_pedigree_single_locus_simulation(void)
{
    int ret;
    tsk_table_collection_t tables;
    int num_inds = 4;
    int ploidy = 2;
    tsk_id_t parents[] = { -1, -1, -1, -1, 0, 0, 1, 1 }; // size num_inds * ploidy
    double time[] = { 1, 1, 0, 0 };
    tsk_flags_t is_sample[] = { 0, 0, 1, 1 };
    tsk_size_t num_nodes;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();

    ret = build_pedigree_sim(
        &msp, &tables, rng, 1, ploidy, num_inds, parents, time, is_sample);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    num_nodes = msp.tables->nodes.num_rows;
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);
    CU_ASSERT_EQUAL(msp.tables->nodes.num_rows, num_nodes);
    msp_verify(&msp, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_STATE);
    /* TODO put in some meaningful tests of the WF pedigree */

    /* msp_print_state(&msp, stdout); */
    /* tsk_table_collection_print_state(&tables, stdout); */

    /* Complete the simulation */
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);

    ret = msp_finalise_tables(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tables.migrations.num_rows, 0);
    CU_ASSERT(tables.nodes.num_rows > 0);
    CU_ASSERT(tables.edges.num_rows > 0);

    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_pedigree_multi_locus_simulation(void)
{
    int ret;
    const char *model_name;
    tsk_table_collection_t tables;
    int num_inds = 4;
    int ploidy = 2;
    tsk_id_t parents[] = { -1, -1, -1, -1, 0, 0, 1, 1 }; // size num_inds * ploidy
    double time[] = { 1, 1, 0, 0 };
    tsk_flags_t is_sample[] = { 0, 0, 1, 1 };
    tsk_size_t num_nodes;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();

    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    num_nodes = msp.tables->nodes.num_rows;
    ret = msp_set_recombination_rate(&msp, 10);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);
    CU_ASSERT_EQUAL(msp.tables->nodes.num_rows, num_nodes);
    msp_verify(&msp, 0);
    /* TODO put in some meaningful tests of the pedigree */

    /* Complete the simulation */
    ret = msp_set_simulation_model_dtwf(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    model_name = msp_get_model_name(&msp);
    CU_ASSERT_STRING_EQUAL(model_name, "dtwf");
    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp, 0);

    ret = msp_finalise_tables(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tables.migrations.num_rows, 0);
    CU_ASSERT(tables.nodes.num_rows > 0);
    CU_ASSERT(tables.edges.num_rows > 0);

    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_pedigree_errors(void)
{
    int ret;
    size_t num_inds = 4;
    size_t ploidy = 2;
    tsk_id_t parents[] = { -1, -1, -1, -1, 0, 0, 1, 1 }; // size num_inds * ploidy
    double time[] = { 1, 1, 0, 0 };
    tsk_flags_t is_sample[] = { 0, 0, 1, 1 };
    msp_t msp;
    tsk_table_collection_t tables;
    gsl_rng *rng = safe_rng_alloc();

    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, 1, num_inds, parents, time, is_sample);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PLOIDY);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, 0, parents, time, is_sample);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    /* Any demographic events during the pedigree sim are errors */
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.5, 0, 1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_UNSUPPORTED_OPERATION);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    time[2] = 1;
    time[3] = 1;
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_TIME_TRAVEL);
    tsk_table_collection_free(&tables);
    msp_free(&msp);
    time[2] = 0;
    time[3] = 0;

    parents[0] = -2;
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample);
    CU_ASSERT_EQUAL_FATAL(ret, msp_set_tsk_error(TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS));
    ret = msp_initialise(&msp);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    parents[0] = 100;
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample);
    CU_ASSERT_EQUAL_FATAL(ret, msp_set_tsk_error(TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS));
    ret = msp_initialise(&msp);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    /* TODO lots more tests when we update the pedigree simulation. */

    gsl_rng_free(rng);
}
int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {

        { "test_pedigree_trio", test_pedigree_trio },
        { "test_pedigree_three_generations", test_pedigree_three_generations },
        { "test_pedigree_single_locus_simulation",
            test_pedigree_single_locus_simulation },
        { "test_pedigree_multi_locus_simulation", test_pedigree_multi_locus_simulation },
        { "test_pedigree_errors", test_pedigree_errors },

        CU_TEST_INFO_NULL,
    };

    return test_main(tests, argc, argv);
}
