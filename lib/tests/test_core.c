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

static void
test_probability_list_select(void)
{
    double short_one = 1.0 - 0.5 * DBL_EPSILON;
    CU_TEST(short_one == nexttoward(1.0, 0));
    double short_half = 0.5 - 0.25 * DBL_EPSILON;
    CU_TEST(short_half == nexttoward(0.5, 0));
    double long_half = 0.5 + 0.5 * DBL_EPSILON;
    CU_TEST(long_half == nexttoward(0.5, 1));
    {
        double probs[2] = { 0.5, 0.5 };
        CU_ASSERT_EQUAL(0, probability_list_select(0.0, 2, probs))
        CU_ASSERT_EQUAL(0, probability_list_select(short_half, 2, probs))
        CU_ASSERT_EQUAL(1, probability_list_select(0.5, 2, probs))
        CU_ASSERT_EQUAL(1, probability_list_select(1.0, 2, probs))
        CU_ASSERT_EQUAL(1, probability_list_select(1.1, 2, probs))
    }
    {
        double probs[2] = { 1.0, 0.0 };
        CU_ASSERT_EQUAL(0, probability_list_select(short_one, 2, probs))
        /* Note: gsl_ran_flat does not return 1.0 per documentation */
        CU_ASSERT_EQUAL(1, probability_list_select(1.0, 2, probs))
    }
    {
        double probs[2] = { 0.0, 1.0 };
        CU_ASSERT_EQUAL(0, probability_list_select(0.0, 0, probs))
        CU_ASSERT_EQUAL(0, probability_list_select(1.1, 0, probs))
        CU_ASSERT_EQUAL(0, probability_list_select(0.0, 1, probs))
        CU_ASSERT_EQUAL(0, probability_list_select(1.1, 1, probs))
        CU_ASSERT_EQUAL(1, probability_list_select(0.0, 2, probs))
        CU_ASSERT_EQUAL(1, probability_list_select(1.1, 2, probs))
    }
    {
        double probs[10] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
        CU_TEST(1.0 != 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1);
        /* Note that 0.1 == 0x1.999999999999ap-4 in hexadecimal */
        CU_ASSERT_EQUAL(0, probability_list_select(0.0, 10, probs))
        CU_ASSERT_EQUAL(4, probability_list_select(short_half, 10, probs))
        CU_ASSERT_EQUAL(5, probability_list_select(0.5, 10, probs))
        CU_ASSERT_EQUAL(5, probability_list_select(long_half, 10, probs))
        CU_ASSERT_EQUAL(9, probability_list_select(0.9, 10, probs))
        CU_ASSERT_EQUAL(9, probability_list_select(short_one, 10, probs))
        CU_ASSERT_EQUAL(9, probability_list_select(1.0, 10, probs))
    }
    {
        double probs[5] = { 0.2, 0.2, 0.2, 0.2, 0.2 };
        CU_ASSERT_EQUAL(0, probability_list_select(0.0, 5, probs))
        CU_ASSERT_EQUAL(4, probability_list_select(0.8, 5, probs))
        CU_ASSERT_EQUAL(4, probability_list_select(short_one, 5, probs))
        CU_ASSERT_EQUAL(4, probability_list_select(1.0, 5, probs))
    }
    {
        double probs[3] = { short_half, DBL_EPSILON / 8.0, 0.5 };
#ifdef MSP_TEST_EXACT_FLOAT_COMPARISONS
        CU_TEST(1.0 == probs[0] + probs[1] + probs[2]);
#endif
        CU_ASSERT_EQUAL(0, probability_list_select(0.0, 3, probs))
        CU_ASSERT_EQUAL(1, probability_list_select(short_half, 3, probs))
        CU_ASSERT_EQUAL(2, probability_list_select(0.5, 3, probs))
        CU_ASSERT_EQUAL(2, probability_list_select(1.0, 3, probs))
    }
    {
        double probs[3] = { short_half, DBL_EPSILON / 16.0, 0.5 };
#ifdef MSP_TEST_EXACT_FLOAT_COMPARISONS
        CU_TEST(1.0 == probs[0] + probs[1] + probs[2]);
#endif
        CU_ASSERT_EQUAL(0, probability_list_select(0.0, 3, probs))
        CU_ASSERT_EQUAL(2, probability_list_select(short_half, 3, probs))
        CU_ASSERT_EQUAL(2, probability_list_select(0.5, 3, probs))
        CU_ASSERT_EQUAL(2, probability_list_select(1.0, 3, probs))
    }
}

static void
test_tskit_version(void)
{
    /* Make sure we don't have any accidental changes to the tskit submodule. */
    CU_ASSERT_EQUAL(TSK_VERSION_MAJOR, 1);
    CU_ASSERT_EQUAL(TSK_VERSION_MINOR, 0);
    CU_ASSERT_EQUAL(TSK_VERSION_PATCH, 0);
}

static void
test_gsl_ran_flat_patch(void)
{
    double left, right, x;
    int i;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    right = 85190021.000000;
    left = nexttoward(right, 0.0);
    gsl_rng_set(rng, 12345);
    x = gsl_ran_flat(rng, left, right);
    CU_ASSERT(!(right > x));
    gsl_rng_set(rng, 12345);
    x = msp_gsl_ran_flat(rng, left, right);
    CU_ASSERT(right > x);
    CU_ASSERT(left <= x);

    left = 0.0;
    right = 1.0;
    for (i = 0; i < 1000; i++) {
        x = msp_gsl_ran_flat(rng, left, right);
        CU_ASSERT(left <= x);
        CU_ASSERT(right > x);
    }

    left = 2.0;
    right = 2.0;
    x = msp_gsl_ran_flat(rng, left, right);
    CU_ASSERT(left == x && right == x);

    gsl_rng_free(rng);
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_strerror", test_strerror },
        { "test_strerror_tskit", test_strerror_tskit },
        { "test_probability_list_select", test_probability_list_select },
        { "test_tskit_version", test_tskit_version },
        { "test_gsl_ran_flat_patch", test_gsl_ran_flat_patch },
        CU_TEST_INFO_NULL,
    };

    return test_main(tests, argc, argv);
}
