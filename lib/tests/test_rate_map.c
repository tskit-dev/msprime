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
test_simple_rate_map(void)
{
    int ret;
    rate_map_t recomb_map;
    double seq_length;
    double positions[][2] = { { 0.0, 1.0 }, { 0.0, 100.0 }, { 0.0, 10000.0 } };
    double rates[] = { 0.0, 1.0 };
    size_t j;

    for (j = 0; j < 3; j++) {
        seq_length = positions[j][1];
        ret = rate_map_alloc(&recomb_map, 1, positions[j], rates);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        rate_map_print_state(&recomb_map, _devnull);
        CU_ASSERT_EQUAL(rate_map_get_size(&recomb_map), 1);
        CU_ASSERT_EQUAL(rate_map_get_sequence_length(&recomb_map), seq_length);
        ret = rate_map_free(&recomb_map);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
}

static void
test_rate_map_errors(void)
{
    int ret;
    rate_map_t rate_map;
    double positions[] = { 0.0, 1.0, 2.0 };
    double rates[] = { 1.0, 2.0, 0.0 };
    double short_positions[] = { 0.0, 0.25, 0.5 };

    ret = rate_map_alloc(&rate_map, 0, positions, rates);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INSUFFICIENT_INTERVALS);
    rate_map_free(&rate_map);

    positions[0] = 1;
    ret = rate_map_alloc(&rate_map, 2, positions, rates);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INTERVAL_MAP_START_NON_ZERO);
    rate_map_free(&rate_map);
    positions[0] = 0;

    positions[1] = 3.0;
    ret = rate_map_alloc(&rate_map, 2, positions, rates);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INTERVAL_POSITIONS_UNSORTED);
    rate_map_free(&rate_map);
    positions[1] = 1.0;

    positions[0] = -1;
    ret = rate_map_alloc(&rate_map, 2, positions, rates);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INTERVAL_MAP_START_NON_ZERO);
    rate_map_free(&rate_map);
    positions[0] = 0.0;

    rates[0] = -1;
    ret = rate_map_alloc(&rate_map, 2, positions, rates);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RATE_VALUE);
    rate_map_free(&rate_map);
    rates[0] = 1.0;

    rates[0] = INFINITY;
    ret = rate_map_alloc(&rate_map, 2, positions, rates);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RATE_VALUE);
    rate_map_free(&rate_map);
    rates[0] = 1.0;

    rates[0] = NAN;
    ret = rate_map_alloc(&rate_map, 2, positions, rates);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RATE_VALUE);
    rate_map_free(&rate_map);
    rates[0] = 1.0;

    ret = rate_map_alloc(&rate_map, 2, positions, rates);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    rate_map_free(&rate_map);

    ret = rate_map_alloc(&rate_map, 2, short_positions, rates);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    rate_map_free(&rate_map);
}

static void
verify_rate_map(double length, double *positions, double *rates, size_t size)
{

    int ret;
    rate_map_t rate_map;
    size_t j;
    double *ret_rates, *ret_positions;

    ret_rates = malloc(size * sizeof(double));
    ret_positions = malloc(size * sizeof(double));

    CU_ASSERT_FATAL(ret_rates != NULL);
    CU_ASSERT_FATAL(ret_positions != NULL);

    ret = rate_map_alloc(&rate_map, size, positions, rates);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    rate_map_print_state(&rate_map, _devnull);
    CU_ASSERT_EQUAL(rate_map_get_size(&rate_map), size);

    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j <= size; j++) {
        CU_ASSERT_EQUAL(rate_map.position[j], positions[j]);
        if (j < size) {
            CU_ASSERT_EQUAL(rate_map.rate[j], rates[j]);
        }
    }
    ret = rate_map_free(&rate_map);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    free(ret_rates);
    free(ret_positions);
}

static void
test_rate_map_examples(void)
{
    double p1[] = { 0, 0.05019838393314813, 0.36933662489552865, 1 };
    double r1[] = { 3.5510784169955434, 4.184964179610539, 3.800808140657212 };
    double p2[] = { 0, 0.125, 0.875, 1, 4, 8, 16 };
    double r2[] = { 0.1, 6.0, 3.333, 2.1, 0.0, 2.2 };

    verify_rate_map(1.0, p1, r1, 3);
    verify_rate_map(1.0, p1, r1, 3);
    verify_rate_map(1.0, p1, r1, 3);

    verify_rate_map(16.0, p2, r2, 6);
    verify_rate_map(16.0, p2, r2, 6);
    verify_rate_map(16.0, p2, r2, 6);
}

static void
test_translate_position_and_recomb_mass(void)
{
    rate_map_t map;
    double p1[] = { 0, 6, 13, 20 };
    double r1[] = { 3, 0, 1 };
    rate_map_alloc(&map, 3, p1, r1);

    /* interval edges */
    CU_ASSERT_EQUAL(rate_map_position_to_mass(&map, 0), 0);
    CU_ASSERT_EQUAL(rate_map_position_to_mass(&map, 6), 18);
    CU_ASSERT_EQUAL(rate_map_position_to_mass(&map, 13), 18);
    CU_ASSERT_EQUAL(rate_map_position_to_mass(&map, 19), 24);
    CU_ASSERT_EQUAL(rate_map_mass_to_position(&map, 0), 0);
    CU_ASSERT_EQUAL(rate_map_mass_to_position(&map, 18), 6);
    CU_ASSERT_EQUAL(rate_map_mass_to_position(&map, 24), 19);

    /* intervals with recombination */
    CU_ASSERT_EQUAL(rate_map_position_to_mass(&map, 4), 12);
    CU_ASSERT_EQUAL(rate_map_position_to_mass(&map, 14), 19);
    CU_ASSERT_EQUAL(rate_map_position_to_mass(&map, 16), 21);
    CU_ASSERT_EQUAL(rate_map_mass_to_position(&map, 12), 4);
    CU_ASSERT_EQUAL(rate_map_mass_to_position(&map, 19), 14);
    CU_ASSERT_EQUAL(rate_map_mass_to_position(&map, 21), 16);

    /* inside recombination interval */
    CU_ASSERT_EQUAL(rate_map_position_to_mass(&map, 8), 18);

    rate_map_free(&map);
}

static void
test_rate_map_mass_between(void)
{
    rate_map_t discrete_map;
    double p1[] = { 0, 6, 13, 20 };
    double r1[] = { 3, 0, 1 };
    double tol = 1e-9;

    rate_map_alloc(&discrete_map, 3, p1, r1);

    CU_ASSERT_DOUBLE_EQUAL_FATAL(rate_map_mass_between(&discrete_map, 0, 2), 6, tol);

    rate_map_free(&discrete_map);
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
    CU_ASSERT_EQUAL(idx, 4);
    // way above - one past, like numpy.searchsorted and C++ std::lower_bound
    idx = msp_binary_interval_search(300, values, size);
    CU_ASSERT_EQUAL(idx, 4);
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
    CU_ASSERT_EQUAL(idx, 1);

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
    CU_ASSERT_EQUAL(idx, 2);

    // All zeros
    double values_zeros[] = { 0, 0, 0 };
    idx = msp_binary_interval_search(-1, values_zeros, 3);
    CU_ASSERT_EQUAL(idx, 0);
    idx = msp_binary_interval_search(0, values_zeros, 3);
    CU_ASSERT_EQUAL(idx, 0);
    idx = msp_binary_interval_search(1, values_zeros, 3);
    CU_ASSERT_EQUAL(idx, 3);
}

static void
test_interval_map(void)
{
    int ret;
    size_t j;
    rate_map_t imap;
    double position[] = { 0, 1, 2, 3, 4, 5, 6 };
    double value[] = { 0.1, 1.1, 2.1, 3.1, 4.1, 5.1 };

    ret = rate_map_alloc(&imap, 0, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INSUFFICIENT_INTERVALS);
    rate_map_free(&imap);

    for (j = 1; j < 7; j++) {
        ret = rate_map_alloc(&imap, j, position, value);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(rate_map_get_size(&imap), j);
        CU_ASSERT_EQUAL(rate_map_get_sequence_length(&imap), j);
        rate_map_print_state(&imap, _devnull);
        rate_map_free(&imap);
    }

    position[0] = 1;
    ret = rate_map_alloc(&imap, 6, position, value);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INTERVAL_MAP_START_NON_ZERO);
    rate_map_free(&imap);

    position[0] = 0;
    position[1] = -1;
    ret = rate_map_alloc(&imap, 6, position, value);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INTERVAL_POSITIONS_UNSORTED);
    rate_map_free(&imap);
    position[1] = 0;

    position[1] = NAN;
    ret = rate_map_alloc(&imap, 6, position, value);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_NONFINITE_INTERVAL_POSITION);
    rate_map_free(&imap);
    position[1] = 0;
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_simple_rate_map", test_simple_rate_map },
        { "test_rate_map_errors", test_rate_map_errors },
        { "test_rate_map_examples", test_rate_map_examples },
        { "test_translate_position_and_recomb_mass",
            test_translate_position_and_recomb_mass },
        { "test_rate_map_mass_between", test_rate_map_mass_between },
        { "test_binary_search", test_msp_binary_interval_search },
        { "test_binary_search_repeating", test_msp_binary_interval_search_repeating },
        { "test_binary_search_edge_cases", test_msp_binary_interval_search_edge_cases },
        { "test_interval_map", test_interval_map },
        CU_TEST_INFO_NULL,
    };

    return test_main(tests, argc, argv);
}
