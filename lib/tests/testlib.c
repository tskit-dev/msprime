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

#include "testlib.h"

/* Global variables used for test in state in the test suite */
char *_tmp_file_name;
FILE *_devnull;

/* Utility function to create a simulation from tables and a set
 * of samples */
int
build_sim(msp_t *msp, tsk_table_collection_t *tables, gsl_rng *rng,
    double sequence_length, size_t num_populations, sample_t *samples,
    size_t num_samples)
{
    size_t j;
    int ret;
    double time = 0.0;
    tsk_id_t population = 0;

    ret = tsk_table_collection_init(tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables->sequence_length = sequence_length;

    for (j = 0; j < num_samples; j++) {
        /* If the input samples are NULL, default to n samples from first
         * population */
        if (samples != NULL) {
            time = samples[j].time;
            population = samples[j].population;
        }
        ret = tsk_node_table_add_row(
            &tables->nodes, TSK_NODE_IS_SAMPLE, time, population, TSK_NULL, NULL, 0);
        CU_ASSERT_EQUAL_FATAL(ret, j);
    }
    for (j = 0; j < num_populations; j++) {
        ret = tsk_population_table_add_row(&tables->populations, NULL, 0);
        CU_ASSERT_EQUAL_FATAL(ret, j);
    }
    ret = msp_alloc(msp, tables, rng);
    return ret;
}

/* Utility function to create a pedigree simulation */
int
build_pedigree_sim(msp_t *msp, tsk_table_collection_t *tables, gsl_rng *rng,
    double sequence_length, size_t ploidy, size_t num_individuals, tsk_id_t *parents,
    double *time, tsk_flags_t *is_sample, tsk_id_t *population)
{
    int ret;
    size_t j, k;
    tsk_id_t pop_id, ind_id, max_pop;
    tsk_flags_t flags;

    ret = tsk_table_collection_init(tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables->sequence_length = sequence_length;

    max_pop = 0;
    /* Insert the pedigree individuals */
    for (j = 0; j < num_individuals; j++) {
        ret = tsk_individual_table_add_row(
            &tables->individuals, 0, NULL, 0, parents + j * ploidy, ploidy, NULL, 0);
        CU_ASSERT_EQUAL_FATAL(ret, j);
        ind_id = ret;
        pop_id = population == NULL ? 0 : population[j];
        if (is_sample == NULL) {
            flags = time[j] == 0;
        } else {
            flags = is_sample[j] ? TSK_NODE_IS_SAMPLE : 0;
        }
        for (k = 0; k < ploidy; k++) {
            ret = tsk_node_table_add_row(
                &tables->nodes, flags, time[j], pop_id, ind_id, NULL, 0);
            CU_ASSERT_FATAL(ret >= 0);
        }
        max_pop = TSK_MAX(max_pop, pop_id);
    }

    for (j = 0; j <= max_pop; j++) {
        ret = tsk_population_table_add_row(&tables->populations, NULL, 0);
        CU_ASSERT_EQUAL_FATAL(ret, j);
    }
    ret = msp_alloc(msp, tables, rng);
    if (ret != 0) {
        goto out;
    }
    ret = msp_set_ploidy(msp, ploidy);
    if (ret != 0) {
        goto out;
    }
    ret = msp_set_simulation_model_fixed_pedigree(msp);
out:
    return ret;
}

gsl_rng *
safe_rng_alloc(void)
{
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    CU_ASSERT_FATAL(rng != NULL);
    return rng;
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
test_main(CU_TestInfo *tests, int argc, char **argv)
{
    int ret;
    CU_pTest test;
    CU_pSuite suite;

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
