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

tsk_id_t deep_n10_parents[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, 1, 8, 1, 8, 0, 5, 0, 5, 2, 7, 2, 7, 4, 9, 4, 9, 3, 6, 3, 6,
    12, 18, 12, 18, 10, 16, 10, 16, 17, 19, 17, 19, 11, 13, 11, 13, 14, 15, 14, 15, 21,
    25, 21, 25, 24, 28, 24, 28, 20, 27, 20, 27, 23, 26, 23, 26, 22, 29, 22, 29, 30, 35,
    30, 35, 32, 36, 32, 36, 33, 37, 33, 37, 31, 34, 31, 34, 38, 39, 38, 39, 40, 42, 40,
    42, 44, 48, 44, 48, 46, 49, 46, 49, 41, 45, 41, 45, 43, 47, 43, 47 };
double deep_n10_time[] = { 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 4.0, 4.0,
    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
    3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

tsk_id_t large_family_parents[]
    = { -1, -1, -1, -1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
          0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
          0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
          0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1 };
double large_family_time[] = { 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

tsk_id_t very_deep_n2_parents[] = { -1, -1, -1, -1, 0, 1, 0, 1, 2, 3, 2, 3, 4, 5, 4, 5,
    6, 7, 6, 7, 8, 9, 8, 9, 10, 11, 10, 11, 12, 13, 12, 13, 14, 15, 14, 15, 16, 17, 16,
    17, 18, 19, 18, 19, 20, 21, 20, 21, 22, 23, 22, 23, 24, 25, 24, 25, 26, 27, 26, 27,
    28, 29, 28, 29, 30, 31, 30, 31, 32, 33, 32, 33, 34, 35, 34, 35, 36, 37, 36, 37, 38,
    39, 38, 39, 40, 41, 40, 41, 42, 43, 42, 43, 44, 45, 44, 45, 46, 47, 46, 47, 48, 49,
    48, 49, 50, 51, 50, 51, 52, 53, 52, 53, 54, 55, 54, 55, 56, 57, 56, 57, 58, 59, 58,
    59, 60, 61, 60, 61 };
double very_deep_n2_time[] = { 31.0, 31.0, 30.0, 30.0, 29.0, 29.0, 28.0, 28.0, 27.0,
    27.0, 26.0, 26.0, 25.0, 25.0, 24.0, 24.0, 23.0, 23.0, 22.0, 22.0, 21.0, 21.0, 20.0,
    20.0, 19.0, 19.0, 18.0, 18.0, 17.0, 17.0, 16.0, 16.0, 15.0, 15.0, 14.0, 14.0, 13.0,
    13.0, 12.0, 12.0, 11.0, 11.0, 10.0, 10.0, 9.0, 9.0, 8.0, 8.0, 7.0, 7.0, 6.0, 6.0,
    5.0, 5.0, 4.0, 4.0, 3.0, 3.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0 };

/* Two independent pedigrees */
tsk_id_t two_pedigrees_parents[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 4, 1, 4,
    0, 2, 0, 2, 5, 8, 5, 8, 6, 7, 6, 7, 9, 10, 9, 10, 11, 12, 11, 12, 13, 15, 13, 15, 14,
    16, 14, 16, -1, -1, -1, -1, -1, -1, 21, 22, 21, 22, 24, 25, 24, 25, 26, 27, 26, 27,
    28, 29, 28, 29, 30, 31, 30, 31, 32, 33, 32, 33, 34, 35, 34, 35 };
double two_pedigrees_time[] = { 4.0, 4.0, 4.0, 4.0, 4.0, 3.0, 3.0, 3.0, 3.0, 2.0, 2.0,
    2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 7.0, 7.0, 7.0, 6.0, 6.0, 5.0, 5.0,
    4.0, 4.0, 3.0, 3.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0 };

static void
verify_complete_pedigree_simulation(
    tsk_table_collection_t *tables, double recombination_rate)
{
    int ret = 0;
    tsk_table_collection_t tables_copy;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_treeseq_t ts;
    tsk_tree_t tree;

    ret = tsk_table_collection_copy(tables, &tables_copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_alloc(&msp, &tables_copy, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    msp_verify(&msp, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_EXIT_COALESCENCE);
    msp_verify(&msp, 0);
    ret = msp_finalise_tables(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_init(&ts, &tables_copy, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (ret = tsk_tree_first(&tree); ret == 1; ret = tsk_tree_next(&tree)) {
        CU_ASSERT(tsk_tree_get_num_roots(&tree) == 1);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);

    msp_free(&msp);
    tsk_table_collection_free(&tables_copy);
    gsl_rng_free(rng);
}

static void
verify_pedigree(double recombination_rate, unsigned long seed,
    tsk_size_t num_individuals, tsk_id_t *parents, double *time, tsk_flags_t *is_sample,
    tsk_id_t *population)
{
    int ret;
    int ploidy = 2;
    tsk_size_t num_nodes;
    msp_t msp;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    gsl_rng *rng = safe_rng_alloc();
    bool coalescence = false;

    ret = build_pedigree_sim(&msp, &tables, rng, 100, ploidy, num_individuals, parents,
        time, is_sample, population);
    /* tsk_table_collection_print_state(&tables, stdout); */
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    num_nodes = msp.tables->nodes.num_rows;
    ret = msp_set_recombination_rate(&msp, recombination_rate);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* msp_print_state(&msp, stdout); */

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    /* printf("ret = %d\n", ret); */
    if (ret == MSP_EXIT_COALESCENCE) {
        coalescence = true;
    } else {
        CU_ASSERT_EQUAL_FATAL(ret, MSP_EXIT_MODEL_COMPLETE);
    }
    CU_ASSERT_EQUAL(msp.tables->nodes.num_rows, num_nodes);
    msp_verify(&msp, 0);
    /* msp_print_state(&msp, stdout); */

    ret = msp_finalise_tables(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tables.migrations.num_rows, 0);
    CU_ASSERT(tables.nodes.num_rows == num_nodes);
    /* msp_print_state(&msp, stdout); */

    verify_complete_pedigree_simulation(&tables, recombination_rate);

    /* Is this a valid tree sequence? */
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* tsk_table_collection_print_state(&tables, stdout); */

    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (ret = tsk_tree_first(&tree); ret == 1; ret = tsk_tree_next(&tree)) {
        if (coalescence) {
            CU_ASSERT(tsk_tree_get_num_roots(&tree) == 1);
        }
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
    gsl_rng_free(rng);
    msp_free(&msp);
    tsk_table_collection_free(&tables);
}

static void
test_trio(void)
{

    tsk_id_t parents[] = { -1, -1, -1, -1, 0, 1 };
    double time[] = { 1, 1, 0 };

    verify_pedigree(0, 1, 3, parents, time, NULL, NULL);
}

static void
test_three_generations(void)
{
    tsk_id_t parents[] = { -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5 };
    double time[] = { 2, 2, 2, 2, 1, 1, 0 };

    verify_pedigree(0, 1234, 7, parents, time, NULL, NULL);
    verify_pedigree(0.1, 1234, 7, parents, time, NULL, NULL);
}

static void
test_sibs(void)
{
    tsk_id_t parents[] = { -1, -1, -1, -1, 0, 1, 0, 1 };
    double time[] = { 1, 1, 0, 0 };
    int seed;

    for (seed = 1; seed < 10; seed++) {
        verify_pedigree(1, seed, 4, parents, time, NULL, NULL);
        verify_pedigree(0, seed, 4, parents, time, NULL, NULL);
    }
}

static void
test_ancient_sample(void)
{
    tsk_id_t parents[] = { -1, -1, -1, -1, 0, 1, 0, 1 };
    double time[] = { 2.0, 2.0, 0.0, 1.0 };
    tsk_flags_t is_sample[] = { false, false, true, true };
    int seed;

    for (seed = 1; seed < 10; seed++) {
        verify_pedigree(0.1, seed, 4, parents, time, is_sample, NULL);
        verify_pedigree(0, seed, 4, parents, time, is_sample, NULL);
    }
}

static void
test_large_family(void)
{
    int seed;
    size_t n = sizeof(large_family_time) / sizeof(*large_family_time);

    for (seed = 1; seed < 10; seed++) {
        verify_pedigree(1, seed, n, large_family_parents, large_family_time, NULL, NULL);
        verify_pedigree(0, seed, n, large_family_parents, large_family_time, NULL, NULL);
    }
}

static void
test_unrelated_n3(void)
{
    tsk_id_t parents[] = { -1, -1, -1, -1, -1, -1 };
    double time[] = { 0.0, 0.0, 0.0 };

    verify_pedigree(0, 1, 3, parents, time, NULL, NULL);
    verify_pedigree(0.1, 1, 3, parents, time, NULL, NULL);
}

static void
test_very_deep_n2(void)
{
    size_t n = sizeof(very_deep_n2_time) / sizeof(*very_deep_n2_time);

    verify_pedigree(0, 1, n, very_deep_n2_parents, very_deep_n2_time, NULL, NULL);
    verify_pedigree(0.1, 1, n, very_deep_n2_parents, very_deep_n2_time, NULL, NULL);
}

static void
test_two_pedigrees(void)
{
    size_t n = sizeof(two_pedigrees_time) / sizeof(*two_pedigrees_time);

    verify_pedigree(0, 1, n, two_pedigrees_parents, two_pedigrees_time, NULL, NULL);
    verify_pedigree(0.1, 1, n, two_pedigrees_parents, two_pedigrees_time, NULL, NULL);
}

static void
test_deep_n2(void)
{
    tsk_id_t parents[]
        = { -1, -1, -1, -1, 0, 1, 0, 1, 2, 3, 2, 3, 4, 5, 4, 5, 6, 7, 6, 7, 8, 9, 8, 9 };
    double time[] = { 5.0, 5.0, 4.0, 4.0, 3.0, 3.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0 };

    verify_pedigree(0, 1, 12, parents, time, NULL, NULL);
    verify_pedigree(0.1, 1, 12, parents, time, NULL, NULL);
}

static void
test_deep_n10(void)
{
    size_t n = sizeof(deep_n10_time) / sizeof(*deep_n10_time);

    verify_pedigree(0, 1, n, deep_n10_parents, deep_n10_time, NULL, NULL);
    verify_pedigree(0.1, 1, n, deep_n10_parents, deep_n10_time, NULL, NULL);
}

static void
test_shallow_n100(void)
{
    tsk_id_t parents[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, 53, 83, 53, 83, 45, 70, 45, 70, 39, 44, 39, 44, 22, 80, 22,
        80, 0, 10, 0, 10, 18, 30, 18, 30, 33, 73, 33, 73, 4, 90, 4, 90, 76, 77, 76, 77,
        12, 31, 12, 31, 55, 88, 55, 88, 26, 42, 26, 42, 15, 69, 15, 69, 40, 96, 40, 96,
        9, 72, 9, 72, 11, 47, 11, 47, 28, 85, 28, 85, 5, 93, 5, 93, 65, 66, 65, 66, 16,
        35, 16, 35, 34, 49, 34, 49, 7, 95, 7, 95, 19, 27, 19, 27, 25, 81, 25, 81, 13, 62,
        13, 62, 3, 24, 3, 24, 17, 38, 17, 38, 8, 78, 8, 78, 6, 64, 6, 64, 36, 89, 36, 89,
        56, 99, 56, 99, 43, 54, 43, 54, 50, 67, 50, 67, 46, 68, 46, 68, 61, 97, 61, 97,
        41, 79, 41, 79, 48, 58, 48, 58, 57, 98, 57, 98, 32, 75, 32, 75, 59, 94, 59, 94,
        63, 84, 63, 84, 29, 37, 29, 37, 1, 52, 1, 52, 2, 21, 2, 21, 23, 87, 23, 87, 74,
        91, 74, 91, 82, 86, 82, 86, 20, 60, 20, 60, 14, 71, 14, 71, 51, 92, 51, 92 };

    double time[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    size_t n = sizeof(time) / sizeof(*time);

    for (int seed = 1; seed < 5; seed++) {
        verify_pedigree(0, seed, n, parents, time, NULL, NULL);
        verify_pedigree(0.1, seed, n, parents, time, NULL, NULL);
    }
}

static void
test_event_by_event(void)
{
    size_t num_inds = sizeof(deep_n10_time) / sizeof(*deep_n10_time);
    size_t ploidy = 2;
    msp_t msp1, msp2;
    tsk_table_collection_t tables1, tables2;
    gsl_rng *rng1 = safe_rng_alloc();
    gsl_rng *rng2 = safe_rng_alloc();
    int ret;

    ret = build_pedigree_sim(&msp1, &tables1, rng1, 100, ploidy, num_inds,
        deep_n10_parents, deep_n10_time, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = build_pedigree_sim(&msp2, &tables2, rng2, 100, ploidy, num_inds,
        deep_n10_parents, deep_n10_time, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_initialise(&msp1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_run(&msp1, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);

    ret = MSP_EXIT_MAX_EVENTS;
    while (ret == MSP_EXIT_MAX_EVENTS) {
        ret = msp_run(&msp2, DBL_MAX, 1);
        CU_ASSERT_FATAL(ret >= 0);
        msp_verify(&msp2, 0);
    }
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);

    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables1, &tables2, 0));

    tsk_table_collection_free(&tables1);
    msp_free(&msp1);
    tsk_table_collection_free(&tables2);
    msp_free(&msp2);
    gsl_rng_free(rng1);
    gsl_rng_free(rng2);
}

static void
test_generation_by_generation(void)
{
    size_t num_inds = sizeof(deep_n10_time) / sizeof(*deep_n10_time);
    size_t ploidy = 2;
    double time;
    msp_t msp1, msp2;
    tsk_table_collection_t tables1, tables2;
    gsl_rng *rng1 = safe_rng_alloc();
    gsl_rng *rng2 = safe_rng_alloc();
    int ret;

    ret = build_pedigree_sim(&msp1, &tables1, rng1, 100, ploidy, num_inds,
        deep_n10_parents, deep_n10_time, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = build_pedigree_sim(&msp2, &tables2, rng2, 100, ploidy, num_inds,
        deep_n10_parents, deep_n10_time, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_initialise(&msp1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_run(&msp1, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);

    for (time = 1; time < 5; time++) {
        ret = msp_run(&msp2, time, UINT32_MAX);
        CU_ASSERT_FATAL(ret >= 0);
        CU_ASSERT_EQUAL(msp2.time, time);
        CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_TIME);
        msp_verify(&msp2, 0);
    }
    ret = msp_run(&msp2, 5, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);
    msp_verify(&msp2, 0);

    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables1, &tables2, 0));

    tsk_table_collection_free(&tables1);
    msp_free(&msp1);
    tsk_table_collection_free(&tables2);
    msp_free(&msp2);
    gsl_rng_free(rng1);
    gsl_rng_free(rng2);
}

static void
test_replicates(void)
{
    size_t num_inds = sizeof(large_family_time) / sizeof(*large_family_time);
    size_t ploidy = 2;
    msp_t msp;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    gsl_rng *rng = safe_rng_alloc();
    int j, ret;

    ret = build_pedigree_sim(&msp, &tables, rng, 100, ploidy, num_inds,
        large_family_parents, large_family_time, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < 20; j++) {
        msp_verify(&msp, 0);
        ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
        CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);
        msp_verify(&msp, 0);

        ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tree_init(&tree, &ts, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        for (ret = tsk_tree_first(&tree); ret == 1; ret = tsk_tree_next(&tree)) {
            /* All the progeny should have found a common ancestor in one of the
             * ploids of each of the two parents. */
            CU_ASSERT(tsk_tree_get_num_roots(&tree) == 4);
        }
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        tsk_treeseq_free(&ts);
        tsk_tree_free(&tree);

        /* msp_print_state(&msp, stdout); */
        ret = msp_reset(&msp);
        CU_ASSERT_EQUAL(ret, 0);
    }

    tsk_table_collection_free(&tables);
    msp_free(&msp);
    gsl_rng_free(rng);
}

static void
test_replicates_ancient_samples(void)
{

    /* Large single family with parents at different times and each child at a
     * different time */
    tsk_id_t parents[] = { -1, -1, -1, -1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
        1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1 };
    double time[] = { 21.0, 22.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
        11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0 };
    tsk_flags_t is_sample[]
        = { 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    size_t num_inds = sizeof(time) / sizeof(*time);
    assert(sizeof(is_sample) / sizeof(*is_sample) == num_inds);
    size_t ploidy = 2;
    msp_t msp;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    gsl_rng *rng = safe_rng_alloc();
    int j, ret;

    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < 20; j++) {
        msp_verify(&msp, 0);
        ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
        CU_ASSERT_EQUAL_FATAL(ret, MSP_EXIT_MODEL_COMPLETE);
        msp_verify(&msp, 0);
        ret = msp_finalise_tables(&msp);
        CU_ASSERT_EQUAL(ret, 0);

        ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tree_init(&tree, &ts, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        for (ret = tsk_tree_first(&tree); ret == 1; ret = tsk_tree_next(&tree)) {
            /* All the progeny should have found a common ancestor in one of the
             * ploids of each of the two parents. */
            CU_ASSERT(tsk_tree_get_num_roots(&tree) == 4);
        }
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        tsk_treeseq_free(&ts);
        tsk_tree_free(&tree);

        /* msp_print_state(&msp, stdout); */
        ret = msp_reset(&msp);
        CU_ASSERT_EQUAL(ret, 0);
    }

    tsk_table_collection_free(&tables);
    msp_free(&msp);
    gsl_rng_free(rng);
}

static void
test_replicates_early_exit(void)
{
    size_t num_inds = sizeof(large_family_time) / sizeof(*large_family_time);
    size_t ploidy = 2;
    msp_t msp;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    gsl_rng *rng = safe_rng_alloc();
    size_t total_offspring;
    int j, k, ploid, ret;

    ret = build_pedigree_sim(&msp, &tables, rng, 100, ploidy, num_inds,
        large_family_parents, large_family_time, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < 20; j++) {
        msp_verify(&msp, 0);
        /* Simulate until the end of generation 0 */
        ret = msp_run(&msp, 0, UINT32_MAX);
        CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_TIME);
        msp_print_state(&msp, _devnull);

        for (k = 0; k < 2; k++) {
            total_offspring = 0;
            for (ploid = 0; ploid < 2; ploid++) {
                total_offspring
                    += avl_count(&msp.pedigree.individuals[k].common_ancestors[ploid]);
            }
            /* The two parents had 50 offspring, so each one should have this many
             * descendents across their two ploids */
            CU_ASSERT_EQUAL(total_offspring, 50);
        }

        ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tree_init(&tree, &ts, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        for (ret = tsk_tree_first(&tree); ret == 1; ret = tsk_tree_next(&tree)) {
            CU_ASSERT(tsk_tree_get_num_roots(&tree) == tsk_treeseq_get_num_samples(&ts));
        }
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        tsk_treeseq_free(&ts);
        tsk_tree_free(&tree);

        ret = msp_reset(&msp);
        CU_ASSERT_EQUAL(ret, 0);
    }

    tsk_table_collection_free(&tables);
    msp_free(&msp);
    gsl_rng_free(rng);
}

static void
test_replicates_exit_coalescence(void)
{
    size_t num_inds = sizeof(very_deep_n2_time) / sizeof(*very_deep_n2_time);
    size_t ploidy = 2;
    msp_t msp;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    gsl_rng *rng = safe_rng_alloc();
    int j, ret;

    ret = build_pedigree_sim(&msp, &tables, rng, 100, ploidy, num_inds,
        very_deep_n2_parents, very_deep_n2_time, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < 20; j++) {
        msp_verify(&msp, 0);
        ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
        CU_ASSERT_EQUAL(ret, MSP_EXIT_COALESCENCE);
        msp_print_state(&msp, _devnull);

        ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tree_init(&tree, &ts, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        for (ret = tsk_tree_first(&tree); ret == 1; ret = tsk_tree_next(&tree)) {
            CU_ASSERT(tsk_tree_get_num_roots(&tree) == 1);
        }
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        tsk_treeseq_free(&ts);
        tsk_tree_free(&tree);

        ret = msp_reset(&msp);
        CU_ASSERT_EQUAL(ret, 0);
    }

    tsk_table_collection_free(&tables);
    msp_free(&msp);
    gsl_rng_free(rng);
}

static void
test_errors(void)
{
    int ret;
    size_t num_inds = 4;
    size_t ploidy = 2;
    tsk_id_t parents[] = { -1, -1, -1, -1, 0, 0, 1, 1 };
    double time[] = { 1, 1, 0, 0 };
    tsk_flags_t is_sample[] = { 0, 0, 1, 1 };
    tsk_id_t population[] = { 0, 0, 0, 0 };
    msp_t msp;
    tsk_table_collection_t tables;
    gsl_rng *rng = safe_rng_alloc();

    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample, population);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, 1, num_inds, parents, time, is_sample, population);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_DTWF_DIPLOID_ONLY);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, 0, parents, time, is_sample, population);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_EMPTY_PEDIGREE);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    /* Any demographic events during the pedigree sim are errors */
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample, population);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.5, 0, 1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_UNSUPPORTED_OPERATION);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    /* Record full ARG is an error */
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample, population);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_store_full_arg(&msp, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_UNSUPPORTED_OPERATION);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    /* non-zero GC rate is an error */
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample, population);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_gene_conversion_rate(&msp, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_DTWF_GC_NOT_SUPPORTED);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    parents[0] = -2;
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample, population);
    CU_ASSERT_EQUAL_FATAL(ret, msp_set_tsk_error(TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS));
    ret = msp_initialise(&msp);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    parents[0] = 100;
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample, population);
    CU_ASSERT_EQUAL_FATAL(ret, msp_set_tsk_error(TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS));
    ret = msp_initialise(&msp);
    tsk_table_collection_free(&tables);
    msp_free(&msp);
    parents[0] = -1;

    is_sample[2] = 0;
    is_sample[3] = 0;
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample, population);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INSUFFICIENT_SAMPLES);
    tsk_table_collection_free(&tables);
    msp_free(&msp);
    is_sample[2] = 1;
    is_sample[3] = 1;

    /* Different times for two nodes in an individual */
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample, population);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    msp_free(&msp);
    tables.nodes.time[0] = 0.001;
    CU_ASSERT_EQUAL_FATAL(msp_alloc(&msp, &tables, rng), 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_ploidy(&msp, 2), 0);
    ret = msp_set_simulation_model_fixed_pedigree(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_PEDIGREE_IND_NODE_TIME_DISAGREE);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    /* Different populations for two nodes in an individual */
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample, population);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    msp_free(&msp);
    tsk_population_table_add_row(&tables.populations, NULL, 0);
    tables.nodes.population[0] = 1;
    CU_ASSERT_EQUAL_FATAL(msp_alloc(&msp, &tables, rng), 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_ploidy(&msp, 2), 0);
    ret = msp_set_simulation_model_fixed_pedigree(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_PEDIGREE_IND_NODE_POPULATION_DISAGREE);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    /* non diploid individuals is an error */
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample, population);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    msp_free(&msp);
    tables.nodes.individual[0] = TSK_NULL;
    CU_ASSERT_EQUAL_FATAL(msp_alloc(&msp, &tables, rng), 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_ploidy(&msp, 2), 0);
    ret = msp_set_simulation_model_fixed_pedigree(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_PEDIGREE_IND_NOT_DIPLOID);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    /* not having two parents is an error */
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample, population);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    msp_free(&msp);
    tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(msp_alloc(&msp, &tables, rng), 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_ploidy(&msp, 2), 0);
    ret = msp_set_simulation_model_fixed_pedigree(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_PEDIGREE_IND_NOT_TWO_PARENTS);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    time[2] = 1;
    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, ploidy, num_inds, parents, time, is_sample, population);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_PEDIGREE_TIME_TRAVEL);
    time[2] = 0;
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    gsl_rng_free(rng);
}

static void
test_internal_samples(void)
{
    int ret;
    tsk_id_t parents[] = { -1, -1, -1, -1, 0, 1 };
    double time[] = { 1.0, 1.0, 0.0 };
    tsk_flags_t is_sample[] = { true, true, true };
    msp_t msp;
    tsk_table_collection_t tables;
    gsl_rng *rng = safe_rng_alloc();

    ret = build_pedigree_sim(
        &msp, &tables, rng, 100, 2, 3, parents, time, is_sample, NULL);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_PEDIGREE_INTERNAL_SAMPLE);
    tsk_table_collection_free(&tables);
    msp_free(&msp);

    gsl_rng_free(rng);
}

static void
test_combined_with_other_models(void)
{
    size_t num_inds = sizeof(large_family_time) / sizeof(*large_family_time);
    size_t ploidy = 2;
    msp_t msp;
    tsk_table_collection_t tables;
    gsl_rng *rng = safe_rng_alloc();
    int ret;

    ret = build_pedigree_sim(&msp, &tables, rng, 100, ploidy, num_inds,
        large_family_parents, large_family_time, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp.model.type, MSP_MODEL_WF_PED);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp.model.type, MSP_MODEL_WF_PED);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);
    ret = msp_set_simulation_model_hudson(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_OTHER_MODELS_WITH_PED);
    CU_ASSERT_EQUAL_FATAL(msp.model.type, MSP_MODEL_WF_PED);

    /* This is cheating slightly because we should be testing if we throw the
     * error after simulating other models, but the effect is the same */
    ret = msp_set_simulation_model_fixed_pedigree(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_OTHER_MODELS_WITH_PED);

    tsk_table_collection_free(&tables);
    msp_free(&msp);
    gsl_rng_free(rng);
}

static void
test_trio_same_pop(void)
{

    tsk_id_t parents[] = { -1, -1, -1, -1, 0, 1 };
    double time[] = { 1, 1, 0 };
    tsk_id_t population[] = { 2, 2, 2 };

    verify_pedigree(0, 1, 3, parents, time, NULL, population);
}

static void
test_trio_child_different_pop(void)
{

    tsk_id_t parents[] = { -1, -1, -1, -1, 0, 1 };
    double time[] = { 1, 1, 0 };
    tsk_id_t population[] = { 2, 2, 1 };

    verify_pedigree(0, 1, 3, parents, time, NULL, population);
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {

        { "test_trio", test_trio },
        { "test_three_generations", test_three_generations },
        { "test_sibs", test_sibs },
        { "test_ancient_sample", test_ancient_sample },
        { "test_internal_samples", test_internal_samples },
        { "test_large_family", test_large_family },
        { "test_unrelated_n3", test_unrelated_n3 },
        { "test_very_deep_n2", test_very_deep_n2 },
        { "test_two_pedigrees", test_two_pedigrees },
        { "test_deep_n2", test_deep_n2 },
        { "test_deep_n10", test_deep_n10 },
        { "test_shallow_n100", test_shallow_n100 },
        { "test_event_by_event", test_event_by_event },
        { "test_generation_by_generation", test_generation_by_generation },
        { "test_replicates", test_replicates },
        { "test_replicates_ancient_samples", test_replicates_ancient_samples },
        { "test_replicates_early_exit", test_replicates_early_exit },
        { "test_replicates_exit_coalescence", test_replicates_exit_coalescence },
        { "test_errors", test_errors },
        { "test_combined_with_other_models", test_combined_with_other_models },

        { "test_trio_same_pop", test_trio_same_pop },
        { "test_trio_child_different_pop", test_trio_child_different_pop },
        CU_TEST_INFO_NULL,
    };

    return test_main(tests, argc, argv);
}
