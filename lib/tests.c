/*
** Copyright (C) 2016 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

#include <CUnit/Basic.h>


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

        }
        CU_ASSERT(fenwick_free(&t) == 0);
    }
}

/* Utility function to return a tree sequence for testing. It is the
 * callers responsilibility to free all memory.
 */
static tree_sequence_t *
get_example_tree_sequence(uint32_t sample_size, uint32_t num_loci,
        double scaled_recombination_rate, double mutation_rate)
{
    int ret;
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tree_sequence_t *tree_seq = malloc(sizeof(tree_sequence_t));
    recomb_map_t *recomb_map = malloc(sizeof(recomb_map_t));
    mutgen_t *mutgen = malloc(sizeof(mutgen_t));
    double positions[] = {0.0, 0.0};
    double rates[] = {0.0, 0.0};

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    CU_ASSERT_FATAL(recomb_map != NULL);

    ret = msp_alloc(msp, sample_size, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(msp, num_loci);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(msp, scaled_recombination_rate);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);

    rates[0] = scaled_recombination_rate;
    positions[1] = num_loci;
    recomb_map_alloc(recomb_map, num_loci, num_loci, positions, rates, 2);

    /* Create the tree_sequence from the state of the simulator.
     * We want to use coalescent time here, so use an Ne of 1/4
     * to cancel scaling factor. */
    ret = tree_sequence_create(tree_seq, msp, recomb_map, 0.25);
    CU_ASSERT_EQUAL(ret, 0);

    ret = mutgen_alloc(mutgen, tree_seq, mutation_rate, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = mutgen_generate(mutgen);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_mutations(tree_seq, mutgen->num_mutations,
            mutgen->mutations, mutgen->parameters, mutgen->environment);
    CU_ASSERT_EQUAL(ret, 0);

    gsl_rng_free(rng);
    msp_free(msp);
    free(msp);
    recomb_map_free(recomb_map);
    free(recomb_map);
    mutgen_free(mutgen);
    free(mutgen);
    return tree_seq;
}

static void
test_vcf(void)
{
    int ret;
    char *str = NULL;
    unsigned int ploidy, num_variants;
    vcf_converter_t *vc = malloc(sizeof(vcf_converter_t));
    tree_sequence_t *ts = get_example_tree_sequence(10, 100, 1.0, 1.0);
    FILE *devnull = fopen("/dev/null", "w");

    CU_ASSERT_FATAL(ts != NULL);
    CU_ASSERT_FATAL(vc != NULL);
    CU_ASSERT_FATAL(devnull != NULL);

    ret = vcf_converter_alloc(vc, ts, 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = vcf_converter_alloc(vc, ts, 3);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = vcf_converter_alloc(vc, ts, 11);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    for (ploidy = 1; ploidy < 3; ploidy++) {
        ret = vcf_converter_alloc(vc, ts, ploidy);
        CU_ASSERT_EQUAL(ret, 0);
        vcf_converter_print_state(vc, devnull);
        ret = vcf_converter_get_header(vc, &str);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_NSTRING_EQUAL("##", str, 2);
        num_variants = 0;
        while ((ret = vcf_converter_next(vc, &str)) == 1) {
            CU_ASSERT_NSTRING_EQUAL("1\t", str, 2);
            num_variants++;
        }
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_TRUE(num_variants > 0);
        vcf_converter_free(vc);
    }

    fclose(devnull);
    free(vc);
    tree_sequence_free(ts);
    free(ts);
}

static void
test_single_locus_simulation(void)
{
    int ret;
    uint32_t j;
    uint32_t n = 10;
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    ret = msp_alloc(msp, n, rng);
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

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
}

static void
test_simplest_records(void)
{
    int ret = 0;
    coalescence_record_t records[] = {
        {0, 0, 1, 2, 1.0, {0, 1}},
    };
    size_t num_records = sizeof(records) / sizeof(coalescence_record_t);
    tree_sequence_t ts;

    ret = tree_sequence_load_records(&ts, 2, 1, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    tree_sequence_free(&ts);
}

static void
test_single_tree_good_records(void)
{
    int ret = 0;
    coalescence_record_t records[] = {
        {0, 0, 1, 4, 1.0, {0, 1}},
        {0, 0, 1, 5, 2.0, {2, 3}},
        {0, 0, 1, 6, 3.0, {4, 5}}
    };
    size_t num_records = sizeof(records) / sizeof(coalescence_record_t);
    tree_sequence_t ts;

    ret = tree_sequence_load_records(&ts, 4, 1, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 7);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    tree_sequence_free(&ts);
}

int
main(void)
{
    CU_pSuite pSuite = NULL;

    /* initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* add a suite to the registry */
    pSuite = CU_add_suite("msprime", NULL, NULL);
    if (NULL == pSuite) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* add the tests to the suite */
    if (
        (NULL == CU_add_test(pSuite, "Fenwick tree", test_fenwick)) ||
        (NULL == CU_add_test(pSuite, "VCF", test_vcf)) ||
        (NULL == CU_add_test(
             pSuite, "Simplest records", test_simplest_records)) ||
        (NULL == CU_add_test(
             pSuite, "Single tree good records", test_single_tree_good_records)) ||
        (NULL == CU_add_test(
             pSuite, "Single locus simulation", test_single_locus_simulation))) {
        CU_cleanup_registry();
        return CU_get_error();
    }
    /* Run all tests using the CUnit Basic interface */
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return CU_get_error();
}
