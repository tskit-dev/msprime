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

    ret = tsk_edge_table_add_row(&tables->edges, 0, 1, 4, 0, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables->edges, 0, 1, 4, 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables->edges, 0, 1, 5, 2, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables->edges, 0, 1, 5, 3, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables->edges, 0, 1, 6, 4, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables->edges, 0, 1, 6, 5, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_population_table_add_row(&tables->populations, NULL, 0);
    CU_ASSERT_FATAL(ret == 0);

    /* Add a site and a mutation */
    if (alphabet == ALPHABET_BINARY) {
        ret = tsk_site_table_add_row(&tables->sites, 0.1, "0", 1, NULL, 0);
        CU_ASSERT_FATAL(ret >= 0);
        ret = tsk_mutation_table_add_row(
            &tables->mutations, 0, 0, -1, TSK_UNKNOWN_TIME, "1", 1, NULL, 0);
        CU_ASSERT_FATAL(ret >= 0);
    } else if (alphabet == ALPHABET_NUCLEOTIDE) {
        ret = tsk_site_table_add_row(&tables->sites, 0.1, "A", 1, NULL, 0);
        CU_ASSERT_FATAL(ret >= 0);
        ret = tsk_mutation_table_add_row(
            &tables->mutations, 0, 0, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
        CU_ASSERT_FATAL(ret >= 0);
    }

    ret = tsk_table_collection_check_integrity(tables, 0);

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
    double pos[] = { 0, 0.1, 0.2, 0.3, 0.4, 1.0 };
    double rate[] = { 0, 0.01, 0.02, 0.03, 0.04 };

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);

    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate_map(&mutgen, 5, pos, rate);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map.size, 5);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map.position[0], 0);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map.position[1], 0.1);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map.position[2], 0.2);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map.position[3], 0.3);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map.position[4], 0.4);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map.position[5], 1.0);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map.rate[0], 0);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map.rate[1], 0.01);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map.rate[2], 0.02);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map.rate[3], 0.03);
    CU_ASSERT_EQUAL_FATAL(mutgen.rate_map.rate[4], 0.04);

    ret = mutgen_generate(&mutgen, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    mutgen_free(&mutgen);
    mutation_model_free(&mut_model);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_mutgen_errors(void)
{
    int ret = 0;
    mutgen_t mutgen;
    double pos[] = { 0, 10 };
    double rate = 0;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    mutation_model_t mut_model, mut_model_binary;

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);

    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = matrix_mutation_model_factory(&mut_model_binary, ALPHABET_BINARY);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate_map(&mutgen, 1, pos, &rate);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INCOMPATIBLE_MUTATION_MAP);
    mutgen_free(&mutgen);

    tables.sequence_length = 10;
    rate = -1;
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate_map(&mutgen, 1, pos, &rate);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RATE_VALUE);
    mutgen_free(&mutgen);

    tables.sequence_length = 0.1;
    ret = mutgen_generate(&mutgen, 0);
    CU_ASSERT_FATAL(msp_is_tsk_error(ret));
    tables.sequence_length = 1.0;
    mutgen_free(&mutgen);

    /* mix of binary and nucleotide alleles */
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model_binary, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 20);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_time_interval(&mutgen, 0.5, 10.0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* we shouldn't error the first time since existing site is nucleotide
     * but not at an integer loction */
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    mutgen_free(&mutgen);
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_time_interval(&mutgen, 0.0, 0.5);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 20);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_UNKNOWN_ALLELE);

    mutgen_free(&mutgen);

    ret = mutgen_alloc(&mutgen, NULL, &tables, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    mutgen_free(&mutgen);

    ret = mutgen_alloc(&mutgen, rng, NULL, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    mutgen_free(&mutgen);

    tables.sequence_length = -1;
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    mutgen_free(&mutgen);

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
    mutation_model_t mut_model;

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* bad mutation generation order */
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_time_interval(&mutgen, 0.0, 0.5);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 20);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_time_interval(&mutgen, 0.5, 1.0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_MUTATION_GENERATION_OUT_OF_ORDER);

    mutgen_free(&mutgen);
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

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables1, ALPHABET_BINARY);
    insert_single_tree(&tables2, ALPHABET_BINARY);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_BINARY);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = mutgen_alloc(&mutgen, rng, &tables1, &mut_model, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables1.mutations.num_rows == 0);
    mutgen_print_state(&mutgen, _devnull);
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    gsl_rng_set(rng, 1);
    ret = mutgen_alloc(&mutgen, rng, &tables1, &mut_model, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 10);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, 0);
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
    ret = mutgen_alloc(&mutgen, rng, &tables2, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 10);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables1, &tables2));
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_table_collection_free(&tables1);
    tsk_table_collection_free(&tables2);
    mutation_model_free(&mut_model);
    gsl_rng_free(rng);
}

static void
test_single_tree_mutgen_keep_sites(void)
{
    int ret = 0;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    tsk_table_collection_t copy, copy2;
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
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy));
    mutgen_free(&mutgen);

    /* and, with discrete sites */
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_KEEP_SITES | MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy));
    mutgen_free(&mutgen);

    /* Turn up the mutation rate */
    gsl_rng_set(rng, 2);
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 10);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites.num_rows > copy.sites.num_rows);
    CU_ASSERT_TRUE(tables.mutations.num_rows > copy.mutations.num_rows);
    mutgen_free(&mutgen);

    /* and, discrete sites */
    gsl_rng_set(rng, 2);
    ret = mutgen_alloc(&mutgen, rng, &copy2, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 10);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy2));
    mutgen_free(&mutgen);

    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 10);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_KEEP_SITES | MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites.num_rows > copy2.sites.num_rows);
    CU_ASSERT_TRUE(tables.mutations.num_rows > copy2.mutations.num_rows);
    mutgen_free(&mutgen);

    /* If we run precisely the same mutations again we should rejection
     * sample away all of the original positions */
    gsl_rng_set(rng, 2);
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 10);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites.num_rows > copy.sites.num_rows);
    CU_ASSERT_TRUE(tables.mutations.num_rows > copy.mutations.num_rows);

    /* add a duplicate site to the original */
    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_TRUE(ret > 0);
    ret = mutgen_generate(&mutgen, MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_DUPLICATE_SITE_POSITION);
    /* and, discrete sites */
    ret = mutgen_generate(&mutgen, MSP_KEEP_SITES | MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_DUPLICATE_SITE_POSITION);

    mutgen_free(&mutgen);
    mutation_model_free(&mut_model);
    tsk_table_collection_free(&tables);
    tsk_table_collection_free(&copy);
    tsk_table_collection_free(&copy2);
    gsl_rng_free(rng);
}

static void
test_single_tree_mutgen_discrete_sites(void)
{
    int ret = 0;
    int j;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tsk_table_collection_t tables;
    mutgen_t mutgen;
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
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 0, 4, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 0, 4, 1, TSK_UNKNOWN_TIME, "G", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 2);

    gsl_rng_set(rng, 2);
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_time_interval(&mutgen, 0.0, 1.0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 10);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.num_rows, 1);
    CU_ASSERT_FATAL(tables.mutations.num_rows > 1);
    for (j = 0; j < tables.sites.num_rows; j++) {
        CU_ASSERT_EQUAL_FATAL(tables.sites.position[j], ceil(tables.sites.position[j]));
    }
    ret = tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_table_collection_clear(&tables);

    /* now, keep: the single tree also has a mutation at position 0.1 */
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);

    ret = tsk_site_table_add_row(&tables.sites, 0.0, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 0, 4, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 0, 4, 1, TSK_UNKNOWN_TIME, "G", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 2);

    ret = mutgen_generate(&mutgen, MSP_KEEP_SITES | MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.num_rows, 2);
    CU_ASSERT_FATAL(tables.mutations.num_rows > 3);
    ret = tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    mutgen_free(&mutgen);
    mutation_model_free(&mut_model);
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
    mutation_model_t mut_model;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < 8192; j++) {
        ret = tsk_mutation_table_add_row(
            &tables.mutations, 0, 0, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
        CU_ASSERT_EQUAL_FATAL(ret, j + 1);
    }

    gsl_rng_set(rng, 2);
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 10);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < 10; j++) {
        ret = mutgen_generate(&mutgen, MSP_KEEP_SITES);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
    CU_ASSERT_TRUE(tables.sites.num_rows > 2);

    mutgen_free(&mutgen);
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
    mutation_model_t mut_model;
    size_t j;
    tsk_id_t node;

    CU_ASSERT_FATAL(rng != NULL);

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 10);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.mutations.num_rows > 0);
    mutgen_print_state(&mutgen, _devnull);

    /* End before start is an error */
    ret = mutgen_set_time_interval(&mutgen, 0, -1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    /* Setting start and end == 0 should give 0 mutations */
    ret = mutgen_set_time_interval(&mutgen, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites.num_rows == 0);
    CU_ASSERT_TRUE(tables.mutations.num_rows == 0);

    /* Setting start = 3 should give 0 mutations */
    ret = mutgen_set_time_interval(&mutgen, 3, DBL_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.sites.num_rows == 0);
    CU_ASSERT_TRUE(tables.mutations.num_rows == 0);

    /* Setting start = 2 should give mutations only above 4 and 5 */
    ret = mutgen_set_time_interval(&mutgen, 2, DBL_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, 0);
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
    mutation_model_t mut_model;

    CU_ASSERT_FATAL(rng != NULL);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);
    ret = tsk_site_table_add_row(&tables.sites, 0.5, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 1);

    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    mutgen_free(&mutgen);
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

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_NUCLEOTIDE);
    insert_single_tree(&copy, ALPHABET_NUCLEOTIDE);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy));

    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy));
    mutgen_free(&mutgen);

    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.num_rows, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.num_rows, 0);
    mutgen_free(&mutgen);

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
    mutation_model_t mut_model;

    CU_ASSERT_FATAL(rng != NULL);
    ret = matrix_mutation_model_factory(&mut_model, ALPHABET_BINARY);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, ALPHABET_BINARY);

    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES | MSP_KEEP_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_MUTATION_GENERATION_OUT_OF_ORDER);
    mutgen_free(&mutgen);

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
    mutation_model_t mut_model;
    size_t len, parent_len;
    char *ds;
    int64_t mut_id, *all_mut_ids;
    int32_t mutation_type_id = 10;
    int64_t next_mutation_id = 23;

    CU_ASSERT_FATAL(rng != NULL);
    ret = slim_mutation_model_alloc(&mut_model, mutation_type_id, next_mutation_id, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, -1);

    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES);
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
    mutation_model_t mut_model;
    char value[21]; /* longest 64 bit int has 20 digits */
    int64_t next_mutation_id;

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, -1);

    /* Trying to generate mutations that overflow raises an error */
    ret = slim_mutation_model_alloc(&mut_model, 0, INT64_MAX, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_MUTATION_ID_OVERFLOW);
    mutgen_free(&mutgen);
    mutation_model_free(&mut_model);

    tsk_mutation_table_clear(&tables.mutations);
    tsk_site_table_clear(&tables.sites);
    /* Try out with a large value that doesn't hit the ceiling */
    next_mutation_id = INT64_MAX - 100;
    ret = slim_mutation_model_alloc(&mut_model, 0, next_mutation_id, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tables.mutations.num_rows > 10);
    mutgen_free(&mutgen);
    mutation_model_free(&mut_model);

    ret = tsk_mutation_table_get_row(&tables.mutations, 0, &mut);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(mut.derived_state_length, 19);
    sprintf(value, "%" PRId64, next_mutation_id);
    CU_ASSERT_NSTRING_EQUAL(value, mut.derived_state, mut.derived_state_length);

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
    mutation_model_t mut_model;
    tsk_site_t site;
    tsk_mutation_t mutation;
    tsk_size_t j;
    char buff[100];

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, -1);

    ret = infinite_alleles_mutation_model_alloc(&mut_model, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES);
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
    mutation_model_t mut_model;
    char value[21]; /* longest 64 bit uint has 20 digits */

    CU_ASSERT_FATAL(rng != NULL);
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_single_tree(&tables, -1);

    ret = infinite_alleles_mutation_model_alloc(&mut_model, UINT64_MAX, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_set_rate(&mutgen, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate(&mutgen, MSP_DISCRETE_SITES);
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
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
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

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
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
        { "test_matrix_mutation_model_errors", test_matrix_mutation_model_errors },
        { "test_matrix_mutation_model_properties",
            test_matrix_mutation_model_properties },
        { "test_slim_mutation_model_errors", test_slim_mutation_model_errors },
        { "test_slim_mutation_model_properties", test_slim_mutation_model_properties },
        CU_TEST_INFO_NULL,
    };

    return test_main(tests, argc, argv);
}
