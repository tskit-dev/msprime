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
        CU_ASSERT_FATAL(ret ==  0);
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

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    tree_sequence_free(&ts);
}

static void
test_simplest_bad_records(void)
{
    int ret = 0;
    coalescence_record_t records[] = {
        {0, 0, 1, 2, 1.0, {0, 1}},
    };
    size_t num_records = sizeof(records) / sizeof(coalescence_record_t);
    tree_sequence_t ts;

    /* An empty sequence should be an error */
    ret = tree_sequence_load_records(&ts, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);

    /* Bad sequence length */
    records[0].right = 0.0;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].right = 1.0;

    /* Equal nodes in the children */
    records[0].children[0] = 1;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].children[0] = 0;

    /* children node >= parent */
    records[0].children[0] = 2;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].children[0] = 0;

    /* Unsorted nodes in the children */
    records[0].children[0] = 1;
    records[0].children[1] = 0;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].children[0] = 0;
    records[0].children[1] = 1;

    /* A gap in the input nodes */
    records[0].node = 3;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].node = 2;

    /* Null parent */
    records[0].node = MSP_NULL_NODE;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].node = 2;

    /* Null child */
    records[0].children[0] = MSP_NULL_NODE;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].children[0] = 0;

    /* Make sure we've preserved a good tree sequence */
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    tree_sequence_free(&ts);
    records[0].node = 2;
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

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 7);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    tree_sequence_free(&ts);
}

static void
test_single_tree_bad_records(void)
{
    int ret = 0;
    coalescence_record_t records[] = {
        {0, 0, 1, 4, 1.0, {0, 1}},
        {0, 0, 1, 5, 2.0, {2, 3}},
        {0, 0, 1, 6, 3.0, {4, 5}}
    };
    size_t num_records = sizeof(records) / sizeof(coalescence_record_t);
    tree_sequence_t ts;

    /* Not sorted in time order */
    records[2].time = 0.5;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[2].time = 3;

    /* Missing node 6 */
    records[2].node = 7;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[2].node = 6;

    /* Left value greater than sequence length */
    records[2].left = 2.0;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[2].left = 0.0;

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    tree_sequence_free(&ts);
}

static void
test_single_tree_good_mutations(void)
{
    int ret = 0;
    coalescence_record_t records[] = {
        {0, 0, 6, 4, 1.0, {0, 1}},
        {0, 0, 6, 5, 2.0, {2, 3}},
        {0, 0, 6, 6, 3.0, {4, 5}}
    };
    mutation_t *mutations = NULL;
    mutation_t *other_mutations = NULL;
    size_t num_records = sizeof(records) / sizeof(coalescence_record_t);
    size_t num_mutations = 6;
    size_t j;
    tree_sequence_t ts;

    mutations = malloc(num_mutations * sizeof(mutation_t));
    CU_ASSERT_FATAL(mutations != NULL);
    for (j = 0; j < num_mutations; j++) {
        mutations[j].position = (double) j;
        mutations[j].node = (uint32_t) j;
    }
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 6.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 7);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);

    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations, "", "");
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), num_mutations);
    other_mutations = malloc(num_mutations * sizeof(mutation_t));
    CU_ASSERT_FATAL(other_mutations != NULL);
    ret = tree_sequence_get_mutations(&ts, other_mutations);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_mutations; j++) {
        CU_ASSERT_EQUAL(mutations[j].position, other_mutations[j].position);
        CU_ASSERT_EQUAL(mutations[j].node, other_mutations[j].node);
    }
    free(mutations);
    free(other_mutations);
    tree_sequence_free(&ts);
}

static void
test_single_tree_bad_mutations(void)
{
    int ret = 0;
    coalescence_record_t records[] = {
        {0, 0, 1, 4, 1.0, {0, 1}},
        {0, 0, 1, 5, 2.0, {2, 3}},
        {0, 0, 1, 6, 3.0, {4, 5}}
    };
    mutation_t mutations[] = {{0, 0}, {0, 1}};
    size_t num_records = sizeof(records) / sizeof(coalescence_record_t);
    size_t num_mutations = 2;
    tree_sequence_t ts;

    /* negative coordinate */
    mutations[0].position = -1.0;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations, "", "");
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_MUTATION);
    tree_sequence_free(&ts);
    mutations[0].position = 0.0;

    /* coordinate > sequence length */
    mutations[0].position = 1.1;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations, "", "");
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_MUTATION);
    tree_sequence_free(&ts);
    mutations[0].position = 0.0;

    /* node = NULL */
    mutations[0].node = MSP_NULL_NODE;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations, "", "");
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_MUTATION);
    tree_sequence_free(&ts);
    mutations[0].node = 0;

    /* node >= num_nodes */
    mutations[0].node = 7;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations, "", "");
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_MUTATION);
    tree_sequence_free(&ts);
    mutations[0].node = 0;

    /* Check to make sure we've maintained legal mutations */
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations, "", "");
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), num_mutations);
    tree_sequence_free(&ts);
}

static void
test_single_tree_iter(void)
{
    int ret = 0;
    coalescence_record_t records[] = {
        {0, 0, 1, 4, 1.0, {0, 1}},
        {0, 0, 1, 5, 2.0, {2, 3}},
        {0, 0, 1, 6, 3.0, {4, 5}}
    };
    uint32_t parents[] = {4, 4, 5, 5, 6, 6, MSP_NULL_NODE};
    size_t num_records = 3;
    uint32_t num_nodes = 7;
    uint32_t u, v, num_leaves, w;
    tree_sequence_t ts;
    sparse_tree_t tree;
    sparse_tree_iterator_t iter;
    FILE *devnull = fopen("/dev/null", "w");

    CU_ASSERT_FATAL(devnull != NULL);

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_alloc_sparse_tree(&ts, &tree, NULL, 0, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_iterator_alloc(&iter, &ts, &tree);
    CU_ASSERT_EQUAL(ret, 0);

    ret = sparse_tree_iterator_next(&iter);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), num_nodes);
    sparse_tree_iterator_print_state(&iter, devnull);

    for (u = 0; u < num_nodes; u++) {
        ret = sparse_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
    }
    ret = sparse_tree_get_num_leaves(&tree, 0, &num_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_leaves, 1);
    ret = sparse_tree_get_num_leaves(&tree, 4, &num_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_leaves, 2);
    ret = sparse_tree_get_num_leaves(&tree, 6, &num_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_leaves, 4);
    ret = sparse_tree_get_mrca(&tree, 0, 1, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 4);
    ret = sparse_tree_get_mrca(&tree, 0, 2, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 6);

    ret = sparse_tree_iterator_next(&iter);
    CU_ASSERT_EQUAL(ret, 0);

    sparse_tree_iterator_free(&iter);
    sparse_tree_free(&tree);
    tree_sequence_free(&ts);
    fclose(devnull);
}

static void
verify_trees(size_t num_records, coalescence_record_t *records,
        size_t num_trees, size_t num_nodes, uint32_t* parents,
        size_t num_mutations, mutation_t *mutations)
{
    int ret;
    uint32_t u, v, j, k, mutation_index;
    tree_sequence_t ts;
    sparse_tree_t tree;
    sparse_tree_iterator_t iter;
    mutation_t *tree_mutations;
    size_t num_tree_mutations;
    FILE *devnull = fopen("/dev/null", "w");

    CU_ASSERT_FATAL(devnull != NULL);

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), num_nodes);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations, "", "");
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), num_mutations);
    ret = tree_sequence_alloc_sparse_tree(&ts, &tree, NULL, 0, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_iterator_alloc(&iter, &ts, &tree);
    CU_ASSERT_EQUAL(ret, 0);

    mutation_index = 0;
    for (j = 0; j < num_trees; j++) {
        ret = sparse_tree_iterator_next(&iter);
        CU_ASSERT_EQUAL(ret, 1);
        sparse_tree_iterator_print_state(&iter, devnull);
        for (u = 0; u < num_nodes; u++) {
            ret = sparse_tree_get_parent(&tree, u, &v);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(v, parents[j * num_nodes + u]);
        }
        ret = sparse_tree_get_mutations(&tree, &num_tree_mutations,
                &tree_mutations);
        CU_ASSERT_EQUAL(ret, 0);
        for (k = 0; k < num_tree_mutations; k++) {
            CU_ASSERT_EQUAL(
                tree_mutations[k].position, mutations[mutation_index].position);
            CU_ASSERT_EQUAL(
                tree_mutations[k].node, mutations[mutation_index].node);
            mutation_index++;
        }
    }
    CU_ASSERT_EQUAL(mutation_index, num_mutations);

    ret = sparse_tree_iterator_next(&iter);
    CU_ASSERT_EQUAL(ret, 0);

    sparse_tree_iterator_free(&iter);
    sparse_tree_free(&tree);
    tree_sequence_free(&ts);
    fclose(devnull);
}

static void
verify_tree_iter_fails(size_t num_records, coalescence_record_t *records,
        size_t num_mutations, mutation_t *mutations,
        uint32_t tree_index, int error_code)
{
    int ret;
    uint32_t index = 0;
    tree_sequence_t ts;
    sparse_tree_t tree;
    sparse_tree_iterator_t iter;

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations, "", "");
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_alloc_sparse_tree(&ts, &tree, NULL, 0, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_iterator_alloc(&iter, &ts, &tree);
    CU_ASSERT_EQUAL(ret, 0);

    while ((ret = sparse_tree_iterator_next(&iter)) == 1) {
        index++;
    }
    CU_ASSERT_EQUAL(index, tree_index);
    CU_ASSERT_EQUAL(ret, error_code);

    sparse_tree_iterator_free(&iter);
    sparse_tree_free(&tree);
    tree_sequence_free(&ts);
}

static void
test_tree_sequence_iter(void)
{
    coalescence_record_t records[] = {
        {0, 2, 10, 4, 0.071, {2, 3}},
        {0, 0, 2,  5, 0.090, {1, 3}},
        {0, 2, 10, 5, 0.090, {1, 4}},
        {0, 0, 7,  6, 0.170, {0, 5}},
        {0, 7, 10, 7, 0.202, {0, 5}},
        {0, 0, 2,  8, 0.253, {2, 6}},
    };
    /* We make one mutation for each tree */
    mutation_t mutations[] = {{1, 2}, {4.5, 0}, {8.5, 5}};
    uint32_t parents[] = {
        6, 5, 8, 5, MSP_NULL_NODE, 6, 8, MSP_NULL_NODE, MSP_NULL_NODE,
        6, 5, 4, 4, 5, 6, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
        7, 5, 4, 4, 5, 7, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
    };
    size_t num_records = 6;
    uint32_t num_nodes = 9;
    uint32_t num_trees = 3;
    uint32_t num_mutations = 3;

    verify_trees(num_records, records, num_trees, num_nodes, parents,
            num_mutations, mutations);
}

static void
test_tree_sequence_bad_records(void)
{
    int ret = 0;
    coalescence_record_t records[] = {
        {0, 2, 10, 4, 0.071, {2, 3}},
        {0, 0, 2,  5, 0.090, {1, 3}},
        {0, 2, 10, 5, 0.090, {1, 4}},
        {0, 0, 7,  6, 0.170, {0, 5}},
        {0, 7, 10, 7, 0.202, {0, 5}},
        {0, 0, 2,  8, 0.253, {2, 6}},
    };
    size_t num_records = sizeof(records) / sizeof(coalescence_record_t);
    tree_sequence_t ts;

    /* Not sorted in time order */
    records[2].time = 0.5;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[2].time = 0.090;

    /* Missing node 6 */
    records[3].node = 9;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[3].node = 6;

    /* Left value greater than right */
    records[0].left = 10.0;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].left = 2.0;

    /* Child value greater than parent*/
    records[3].children[1] = 7;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[3].children[1] = 5;

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    tree_sequence_free(&ts);
}

static void
test_single_tree_iter_failure(void)
{
    int ret = 0;
    coalescence_record_t records[] = {
        {0, 0, 1, 4, 1.0, {0, 1}},
        {0, 0, 1, 5, 2.0, {2, 3}},
        {0, 0, 1, 6, 3.0, {4, 5}}
    };
    size_t num_records = sizeof(records) / sizeof(coalescence_record_t);
    tree_sequence_t ts;
    sparse_tree_t tree;
    sparse_tree_iterator_t iter;

    /* change the left coordinate of one record so we can't build the tree */
    records[0].left = 0.5;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_alloc_sparse_tree(&ts, &tree, NULL, 0, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_iterator_alloc(&iter, &ts, &tree);
    CU_ASSERT_EQUAL(ret, 0);

    ret = sparse_tree_iterator_next(&iter);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);

    sparse_tree_iterator_free(&iter);
    sparse_tree_free(&tree);
    tree_sequence_free(&ts);
    records[0].left = 0.0;
}

static void
test_tree_sequence_iter_failure(void)
{
    coalescence_record_t records[] = {
        {0, 2, 10, 4, 0.071, {2, 3}},
        {0, 0, 2,  5, 0.090, {1, 3}},
        {0, 2, 10, 5, 0.090, {1, 4}},
        {0, 0, 7,  6, 0.170, {0, 5}},
        {0, 7, 10, 7, 0.202, {0, 5}},
        {0, 0, 2,  8, 0.253, {2, 6}},
    };
    uint32_t parents[] = {
        6, 5, 8, 5, MSP_NULL_NODE, 6, 8, MSP_NULL_NODE, MSP_NULL_NODE,
        6, 5, 4, 4, 5, 6, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
        7, 5, 4, 4, 5, 7, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
    };
    size_t num_records = 6;
    uint32_t num_nodes = 9;
    uint32_t num_trees = 3;

    /* The first tree is missing a record */
    records[5].left = 1;
    verify_tree_iter_fails(num_records, records, 0, NULL, 0,
            MSP_ERR_BAD_COALESCENCE_RECORDS);
    records[5].left = 0;
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

    /* Make a gap between adjacent records */
    records[1].right = 1;
    verify_tree_iter_fails(num_records, records, 0, NULL, 1,
            MSP_ERR_BAD_COALESCENCE_RECORDS);
    records[1].right = 2;
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

    /* Make a gap in the middle of the sequence */
    records[0].left = 7;
    records[2].left = 7;
    records[3].right = 2;
    verify_tree_iter_fails(num_records, records, 0, NULL, 1,
            MSP_ERR_BAD_COALESCENCE_RECORDS);
    records[0].left = 2;
    records[2].left = 2;
    records[3].right = 7;
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

    /* Make a gap before the last tree */
    records[4].left = 8;
    verify_tree_iter_fails(num_records, records, 0, NULL, 2,
            MSP_ERR_BAD_COALESCENCE_RECORDS);
    records[4].left = 7;
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

    /* Add an extra record to the first tree */
    records[4].left = 2;
    verify_tree_iter_fails(num_records, records, 0, NULL, 1,
            MSP_ERR_BAD_COALESCENCE_RECORDS);
    records[4].left = 7;
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

    /* Add an extra record to the second tree */
    records[0].left = 0;
    verify_tree_iter_fails(num_records, records, 0, NULL, 0,
            MSP_ERR_BAD_COALESCENCE_RECORDS);
    records[0].left = 2;
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

    /* Remove the last record */
    verify_tree_iter_fails(num_records - 1, records, 0, NULL, 0,
            MSP_ERR_BAD_COALESCENCE_RECORDS);
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

}

static void
test_tree_sequence_mutations_iter_failure(void)
{
    coalescence_record_t records[] = {
        {0, 2, 10, 4, 0.071, {2, 3}},
        {0, 0, 2,  5, 0.090, {1, 3}},
        {0, 2, 10, 5, 0.090, {1, 4}},
        {0, 0, 7,  6, 0.170, {0, 5}},
        {0, 7, 10, 7, 0.202, {0, 5}},
        {0, 0, 2,  8, 0.253, {2, 6}},
    };
    uint32_t parents[] = {
        6, 5, 8, 5, MSP_NULL_NODE, 6, 8, MSP_NULL_NODE, MSP_NULL_NODE,
        6, 5, 4, 4, 5, 6, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
        7, 5, 4, 4, 5, 7, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
    };
    mutation_t mutations[] = {{0, 0}};
    size_t num_records = 6;
    size_t num_mutations = 1;
    uint32_t num_nodes = 9;
    uint32_t num_trees = 3;

    /* Mutation over the root in the first tree */
    mutations[0].node = 8;
    verify_tree_iter_fails(num_records, records, num_mutations, mutations, 0,
            MSP_ERR_BAD_MUTATION);
    mutations[0].node = 0;
    verify_trees(num_records, records, num_trees, num_nodes, parents,
            num_mutations, mutations);

    /* Mutation at a node that does not exist in the first tree */
    mutations[0].node = 7;
    verify_tree_iter_fails(num_records, records, num_mutations, mutations, 0,
            MSP_ERR_BAD_MUTATION);
    mutations[0].node = 0;
    verify_trees(num_records, records, num_trees, num_nodes, parents,
            num_mutations, mutations);

    /* Mutation over the root in the first tree */
    mutations[0].node = 6;
    mutations[0].position = 2;
    verify_tree_iter_fails(num_records, records, num_mutations, mutations, 1,
            MSP_ERR_BAD_MUTATION);
    mutations[0].node = 0;
    verify_trees(num_records, records, num_trees, num_nodes, parents,
            num_mutations, mutations);

    /* Mutation at a node that does not exist in the second tree */
    mutations[0].node = 8;
    mutations[0].position = 2;
    verify_tree_iter_fails(num_records, records, num_mutations, mutations, 1,
            MSP_ERR_BAD_MUTATION);
    mutations[0].node = 0;
    verify_trees(num_records, records, num_trees, num_nodes, parents,
            num_mutations, mutations);

    /* Mutation over the root in the third tree */
    mutations[0].node = 7;
    mutations[0].position = 7;
    verify_tree_iter_fails(num_records, records, num_mutations, mutations, 2,
            MSP_ERR_BAD_MUTATION);
    mutations[0].node = 0;
    verify_trees(num_records, records, num_trees, num_nodes, parents,
            num_mutations, mutations);

    /* Mutation at a node that does not exist in the third tree */
    mutations[0].node = 6;
    mutations[0].position = 7;
    verify_tree_iter_fails(num_records, records, num_mutations, mutations, 2,
            MSP_ERR_BAD_MUTATION);
    mutations[0].node = 0;
    verify_trees(num_records, records, num_trees, num_nodes, parents,
            num_mutations, mutations);
}

static void
test_records_equivalent(void)
{
    int ret;
    tree_sequence_t *ts1 = get_example_tree_sequence(10, 100, 1.0, 1.0);
    tree_sequence_t ts2;
    coalescence_record_t *records, r1, r2;
    size_t j, num_records;

    CU_ASSERT_FATAL(ts1 != NULL);
    num_records = tree_sequence_get_num_coalescence_records(ts1);
    records = malloc(num_records * sizeof(coalescence_record_t));
    CU_ASSERT_FATAL(records != NULL);
    for (j = 0; j < num_records; j++) {
        ret = tree_sequence_get_record(ts1, j, &records[j], MSP_ORDER_TIME);
        CU_ASSERT_EQUAL(ret, 0);
    }
    ret = tree_sequence_load_records(&ts2, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(
        tree_sequence_get_sample_size(ts1),
        tree_sequence_get_sample_size(&ts2));
    CU_ASSERT_EQUAL(
        tree_sequence_get_sequence_length(ts1),
        tree_sequence_get_sequence_length(&ts2));
    for (j = 0; j < num_records; j++) {
        ret = tree_sequence_get_record(ts1, j, &records[j], MSP_ORDER_TIME);
        CU_ASSERT_EQUAL(ret, 0);
    }
    for (j = 0; j < num_records; j++) {
        ret = tree_sequence_get_record(ts1, j, &r1, MSP_ORDER_TIME);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_get_record(&ts2, j, &r2, MSP_ORDER_TIME);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(r1.left, r1.left);
        CU_ASSERT_EQUAL(r1.right, r1.right);
        CU_ASSERT_EQUAL(r1.node, r1.node);
        CU_ASSERT_EQUAL(r1.children[0], r1.children[0]);
        CU_ASSERT_EQUAL(r1.children[1], r1.children[1]);
        CU_ASSERT_EQUAL(r1.time, r1.time);
        CU_ASSERT_EQUAL(r1.population_id, r1.population_id);
    }
    tree_sequence_free(&ts2);
    tree_sequence_free(ts1);
    free(ts1);
    free(records);
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
             pSuite, "Simplest bad records", test_simplest_bad_records)) ||
        (NULL == CU_add_test(
             pSuite, "Single tree good records",
             test_single_tree_good_records)) ||
        (NULL == CU_add_test(
             pSuite, "Single tree bad records",
             test_single_tree_bad_records)) ||
        (NULL == CU_add_test(
             pSuite, "Single tree good mutations",
             test_single_tree_good_mutations)) ||
        (NULL == CU_add_test(
             pSuite, "Single tree bad mutations",
             test_single_tree_bad_mutations)) ||
        (NULL == CU_add_test(
             pSuite, "Single tree iterator",
             test_single_tree_iter)) ||
        (NULL == CU_add_test(
             pSuite, "Tree sequence iterator",
             test_tree_sequence_iter)) ||
        (NULL == CU_add_test(
             pSuite, "Tree sequence bad records",
             test_tree_sequence_bad_records)) ||
        (NULL == CU_add_test(
             pSuite, "Single tree iterator failure",
             test_single_tree_iter_failure)) ||
        (NULL == CU_add_test(
             pSuite, "Tree sequence iterator failure",
             test_tree_sequence_iter_failure)) ||
        (NULL == CU_add_test(
             pSuite, "Tree sequence mutation iterator failure",
             test_tree_sequence_mutations_iter_failure)) ||
        (NULL == CU_add_test(
             pSuite, "Test records equivalent after import",
             test_records_equivalent)) ||
        (NULL == CU_add_test(
             pSuite, "Single locus simulation",
             test_single_locus_simulation))) {
        CU_cleanup_registry();
        return CU_get_error();
    }
    /* Run all tests using the CUnit Basic interface */
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return CU_get_error();
}
