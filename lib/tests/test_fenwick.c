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
    CU_ASSERT_EQUAL(fenwick_find(&t, 1), 0);
    CU_ASSERT_EQUAL(fenwick_find(&t, 0), 0);

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
test_fenwick_rebuild(void)
{
    fenwick_t t;
    size_t n = 14;
    size_t j;
    double drift_before, drift_after;

    CU_ASSERT(fenwick_alloc(&t, n) == 0);
    for (j = 1; j <= n; j++) {
        fenwick_set_value(&t, j, 0.1);
    }
    drift_before = fenwick_get_numerical_drift(&t);
    CU_ASSERT_TRUE(fenwick_rebuild_required(&t));
    fenwick_print_state(&t, _devnull);
    fenwick_rebuild(&t);
    drift_after = fenwick_get_numerical_drift(&t);
    /* After we've rebuilt we signal that a rebuild is not required */
    CU_ASSERT_FALSE(fenwick_rebuild_required(&t));
    /* even though the drift values are identical (this is as good as
     * we can do with these numbers) */
    CU_ASSERT_EQUAL(drift_before, drift_after);

    CU_ASSERT(fenwick_free(&t) == 0);
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_fenwick", test_fenwick },
        { "test_fenwick_expand", test_fenwick_expand },
        { "test_fenwick_zero_values", test_fenwick_zero_values },
        { "test_fenwick_drift", test_fenwick_drift },
        { "test_fenwick_rebuild", test_fenwick_rebuild },
        CU_TEST_INFO_NULL,
    };

    return test_main(tests, argc, argv);
}
