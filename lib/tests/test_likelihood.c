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
test_likelihood_errors(void)
{
    int ret;
    double lik = 0;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 1.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 1.0, 2, 0, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 1.0, 2, 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    /* Site has two mutations */
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 0, 1, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Ne <= 0 is an error */
    ret = msp_log_likelihood_arg(&ts, 1, 0, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_POPULATION_SIZE);
    ret = msp_log_likelihood_arg(&ts, -1, 0, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_POPULATION_SIZE);

    ret = msp_unnormalised_log_likelihood_mut(&ts, 0, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(lik, -DBL_MAX);

    ret = msp_unnormalised_log_likelihood_mut(&ts, 1, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_likelihood_zero_edges(void)
{
    int ret;
    double lik = 0;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Zero edges should give a likelihood of zero */
    ret = msp_log_likelihood_arg(&ts, 1, 1, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(lik, 0);

    /* Zero mutations, so we get zero. */
    ret = msp_unnormalised_log_likelihood_mut(&ts, 0, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(lik, 0);
    ret = msp_unnormalised_log_likelihood_mut(&ts, 1, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(lik, 0);

    tsk_treeseq_free(&ts);

    /* Add in a site and a mutation */
    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Zero edges should give a likelihood of zero */
    ret = msp_log_likelihood_arg(&ts, 1, 1, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(lik, 0);
    ret = msp_log_likelihood_arg(&ts, 0, 1, &lik);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(lik, 0);

    /* Zero mutation rate gives a likelihood of zero when there are mutations in
     * the ARG */
    ret = msp_unnormalised_log_likelihood_mut(&ts, 0, &lik);
    CU_ASSERT_EQUAL_FATAL(lik, -DBL_MAX);
    ret = msp_unnormalised_log_likelihood_mut(&ts, 1, &lik);
    CU_ASSERT_EQUAL_FATAL(lik, -DBL_MAX);

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_likelihood_three_leaves(void)
{
    int i;
    int ret;
    double rho[] = { 0.1, 1, 10 };
    double theta[] = { 0.1, 1, 10 };
    double ll_exact;
    double lik = 0;
    double tree_length = 19.0 / 8;
    double tol = 1e-9;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 3, 2, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 4, 2, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 1, 5, 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 5, 3, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.25, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 1, 6, 0, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 6, 4, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.5, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 1, 7, 5, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 1, 7, 6, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 1, 1, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 2, 3, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 3, 0, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 4, 5, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.2, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.3, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.4, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.45, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    tsk_population_table_add_row(&tables.populations, NULL, 0);

    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (i = 0; i < 3; i++) {
        ll_exact = log(rho[i]) - (3 + 3 * rho[i]) * 0.1;
        ll_exact -= (6 + 3 * rho[i]) * 0.15;
        ll_exact -= (3 + 2.5 * rho[i]) * 0.25;
        ll_exact -= (1 + 2 * rho[i]) * 0.5;
        ret = msp_log_likelihood_arg(&ts, rho[i], 0.5, &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }
    for (i = 0; i < 3; i++) {
        ll_exact = (5 * log(tree_length * theta[i]) - tree_length * theta[i]);
        ll_exact -= 2 * log(4 * tree_length);
        ll_exact -= 2 * log(tree_length);
        ll_exact += log(3 / (4 * tree_length));
        ret = msp_unnormalised_log_likelihood_mut(&ts, theta[i], &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_likelihood_two_mrcas(void)
{
    int i;
    int ret;
    double rho[] = { 0.1, 1, 10 };
    double theta[] = { 0.1, 1, 10 };
    double ll_exact;
    double lik = 0;
    double tree_length = 1.5;
    double tol = 1e-9;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 2, 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 3, 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 4, 0, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 5, 0, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.15, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.15, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 6, 2, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 6, 4, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.5, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 7, 3, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 7, 5, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_mutation_table_add_row(
        &tables.mutations, 0, 1, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 1, 4, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 2, 3, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.2, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.7, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    tsk_population_table_add_row(&tables.populations, NULL, 0);

    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (i = 0; i < 3; i++) {
        ll_exact = log(rho[i]) - (1 + 2 * rho[i]) * 0.1;
        ll_exact += log(rho[i]) - (3 + 2 * rho[i]) * 0.05;
        ll_exact -= (6 + 2 * rho[i]) * 0.35;
        ll_exact -= (1 + rho[i]) * 0.5;
        ret = msp_log_likelihood_arg(&ts, rho[i], 0.5, &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }
    for (i = 0; i < 3; i++) {
        ll_exact = 3 * log(tree_length * theta[i]) - tree_length * theta[i];
        ll_exact -= log(tree_length);
        ll_exact -= 2 * log(2 * tree_length);
        ret = msp_unnormalised_log_likelihood_mut(&ts, theta[i], &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_likelihood_material_overhang(void)
{
    int i;
    int ret;
    double rho[] = { 0.1, 1, 10 };
    double theta[] = { 0.1, 1, 10 };
    double ll_exact;
    double lik = 0;
    double tree_length = 81.0 / 50;
    double tol = 1e-9;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 2, 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 3, 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.7, 4, 0, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.7, 1, 5, 0, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.15, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.15, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.5, 6, 2, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.7, 6, 4, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.5, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 7, 3, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.7, 1, 7, 5, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 0.7, 8, 6, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 0.7, 8, 7, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 1.3, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_mutation_table_add_row(
        &tables.mutations, 0, 1, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 1, 4, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 2, 3, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.6, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.75, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    tsk_population_table_add_row(&tables.populations, NULL, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (i = 0; i < 3; i++) {
        ll_exact = log(rho[i]) - (1 + 2 * rho[i]) * 0.1;
        ll_exact += log(rho[i]) - (3 + 2 * rho[i]) * 0.05;
        ll_exact -= (6 + 2 * rho[i]) * 0.35;
        ll_exact -= (3 + rho[i]) * 0.5;
        ll_exact -= (1 + 0.4 * rho[i]) * 0.3;
        ret = msp_log_likelihood_arg(&ts, rho[i], 0.5, &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }
    for (i = 0; i < 3; i++) {
        ll_exact = 3 * log(tree_length * theta[i]) - tree_length * theta[i];
        ll_exact -= log(2 * tree_length);
        ll_exact += log(1.3 / tree_length);
        ll_exact -= log(tree_length);
        ret = msp_unnormalised_log_likelihood_mut(&ts, theta[i], &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_likelihood_material_gap(void)
{
    int i;
    int ret;
    double rho[] = { 0.1, 1, 10 };
    double theta[] = { 0.1, 1, 10 };
    double ll_exact;
    double lik = 0;
    double tree_length = 0.84;
    double tol = 1e-9;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 2, 0, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 1, 3, 0, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 0.5, 4, 3, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 5, 3, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.2, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.2, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 6, 2, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 6, 5, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.3, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 1, 7, 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 7, 6, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 7, 6, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.4, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 0.5, 8, 4, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 0.5, 8, 7, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.5, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_mutation_table_add_row(
        &tables.mutations, 0, 2, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 1, 3, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 2, 7, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 3, 6, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 4, 5, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.35, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.4, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.6, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.8, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    tsk_population_table_add_row(&tables.populations, NULL, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (i = 0; i < 3; i++) {
        ll_exact = log(rho[i]) - (1 + 2 * rho[i]) * 0.1;
        ll_exact += log(rho[i]) - (3 + 2 * rho[i]) * 0.1;
        ll_exact -= (6 + 2 * rho[i]) * 0.1;
        ll_exact -= (3 + 2.2 * rho[i]) * 0.1;
        ll_exact -= (1 + 0.4 * rho[i]) * 0.1;
        ret = msp_log_likelihood_arg(&ts, rho[i], 0.5, &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }
    for (i = 0; i < 3; i++) {
        ll_exact = 5 * log(tree_length * theta[i]) - tree_length * theta[i];
        ll_exact += 3 * log(0.4 / tree_length);
        ll_exact += 2 * log(0.5 / tree_length);
        ret = msp_unnormalised_log_likelihood_mut(&ts, theta[i], &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_likelihood_recombination_in_material_gap(void)
{
    int i;
    int ret;
    double rho[] = { 0.1, 1, 10 };
    double theta[] = { 0.1, 1, 10 };
    double ll_exact;
    double lik = 0;
    double tree_length = 1.34;
    double tol = 1e-9;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 2, 0, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 1, 3, 0, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.1, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 0.5, 4, 3, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 5, 3, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.2, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.2, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 6, 2, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 6, 5, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.3, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 7, 6, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 8, 6, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.4, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(
        &tables.nodes, MSP_NODE_IS_RE_EVENT, 0.4, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 0.5, 9, 4, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.5, 1, 9, 8, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.5, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0, 1, 10, 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0, 0.3, 10, 7, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.6, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 1, 11, 9, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_edge_table_add_row(&tables.edges, 0.3, 1, 11, 10, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 0.7, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_mutation_table_add_row(
        &tables.mutations, 0, 2, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 1, 3, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 2, 6, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_mutation_table_add_row(
        &tables.mutations, 3, 5, -1, TSK_UNKNOWN_TIME, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    ret = tsk_site_table_add_row(&tables.sites, 0.1, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.35, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.6, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_table_add_row(&tables.sites, 0.8, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    tsk_population_table_add_row(&tables.populations, NULL, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (i = 0; i < 3; i++) {
        ll_exact = log(rho[i]) - (1 + 2 * rho[i]) * 0.1;
        ll_exact += log(rho[i]) - (3 + 2 * rho[i]) * 0.1;
        ll_exact -= (6 + 2 * rho[i]) * 0.1;
        ll_exact += log(0.2 * rho[i]) - (3 + 2.2 * rho[i]) * 0.1;
        ll_exact -= (6 + 2 * rho[i]) * 0.1;
        ll_exact -= (3 + 2 * rho[i]) * 0.1;
        ll_exact -= (1 + 1.4 * rho[i]) * 0.1;
        ret = msp_log_likelihood_arg(&ts, rho[i], 0.5, &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }
    for (i = 0; i < 3; i++) {
        ll_exact = 4 * log(tree_length * theta[i]) - tree_length * theta[i];
        ll_exact += log(0.6 / tree_length);
        ll_exact += 3 * log(0.7 / tree_length);
        ret = msp_unnormalised_log_likelihood_mut(&ts, theta[i], &lik);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL(ll_exact, lik, tol);
    }

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_likelihood_errors", test_likelihood_errors },
        { "test_likelihood_zero_edges", test_likelihood_zero_edges },
        { "test_likelihood_three_leaves", test_likelihood_three_leaves },
        { "test_likelihood_two_mrcas", test_likelihood_two_mrcas },
        { "test_likelihood_material_overhang", test_likelihood_material_overhang },
        { "test_likelihood_material_gap", test_likelihood_material_gap },
        { "test_likelihood_recombination_in_material_gap",
            test_likelihood_recombination_in_material_gap },
        CU_TEST_INFO_NULL,
    };

    return test_main(tests, argc, argv);
}
