/*
** Copyright (C) 2016-2017 University of Oxford
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
#include <stdio.h>
#include <unistd.h>

#include <hdf5.h>
#include <gsl/gsl_math.h>
#include <CUnit/Basic.h>

/* Global variables used for test in state in the test suite */

char * _tmp_file_name;
FILE * _devnull;

#define SIMPLE_BOTTLENECK 0
#define INSTANTANEOUS_BOTTLENECK 1

typedef struct {
    int type;
    double time;
    uint32_t population_id;
    double parameter;
} bottleneck_desc_t;

/* Example tree sequences used in some of the tests. */


/* Simple single tree example. */
const char *single_tree_ex_nodes =
    "1  0   0\n"
    "1  0   0\n"
    "1  0   0\n"
    "1  0   0\n"
    "0  1   0\n"
    "0  2   0\n"
    "0  3   0\n";
const char *single_tree_ex_edgesets =
    "0  1   4   0,1\n"
    "0  1   5   2,3\n"
    "0  1   6   4,5\n";
const char *single_tree_ex_sites =
    "0.1  0\n"
    "0.2  0\n"
    "0.3  0\n";
const char *single_tree_ex_mutations =
    "0    2     1\n"
    "1    4     1\n"
    "1    0     0\n"  /* Back mutation over 0 */
    "2    0     1\n"  /* recurrent mutations over leaves */
    "2    1     1\n"
    "2    2     1\n"
    "2    3     1\n";

/* Example from the PLOS paper */
const char *paper_ex_nodes =
    "1  0   0\n"
    "1  0   0\n"
    "1  0   0\n"
    "1  0   0\n"
    "0  0.071   0\n"
    "0  0.090   0\n"
    "0  0.170   0\n"
    "0  0.202   0\n"
    "0  0.253   0\n";
const char *paper_ex_edgesets =
    "2 10 4 2,3\n"
    "0 2  5 1,3\n"
    "2 10 5 1,4\n"
    "0 7  6 0,5\n"
    "7 10 7 0,5\n"
    "0 2  8 2,6\n";
/* We make one mutation for each tree */
const char *paper_ex_sites =
    "1      0\n"
    "4.5    0\n"
    "8.5    0\n";
const char *paper_ex_mutations =
    "0      2   1\n"
    "1      0   1\n"
    "2      5   1\n";

/* An example of a nonbinary tree sequence */
const char *nonbinary_ex_nodes =
    "1  0   0\n"
    "1  0   0\n"
    "1  0   0\n"
    "1  0   0\n"
    "1  0   0\n"
    "1  0   0\n"
    "1  0   0\n"
    "1  0   0\n"
    "0  0.01    0\n"
    "0  0.068   0\n"
    "0  0.130   0\n"
    "0  0.279   0\n"
    "0  0.405   0\n";
const char *nonbinary_ex_edgesets =
    "0	100	8	0,1,2,3\n"
    "0	100	9	6,8\n"
    "0	17	10	4,5,7\n"
    "17	100	10	4,7\n"
    "17	100	11	5,9\n"
    "0	17	12	9,10\n"
    "17	100	12	10,11";
const char *nonbinary_ex_sites =
        "1  0\n"
        "18 0\n";
const char *nonbinary_ex_mutations =
    "0    2   1\n"
    "1    11  1";

/* An example of a tree sequence with unary nodes. */

const char *unary_ex_nodes =
    "1  0   0\n"
    "1  0   0\n"
    "1  0   0\n"
    "1  0   0\n"
    "0  0.071   0\n"
    "0  0.090   0\n"
    "0  0.170   0\n"
    "0  0.202   0\n"
    "0  0.253   0\n";
const char *unary_ex_edgesets =
    "2 10 4 2,3\n"
    "0 2  5 1,3\n"
    "2 10 5 1,4\n"
    "0 7  6 0,5\n"
    "7 10 7 0,5\n"
    "0 2  7 2\n"
    "0 2  8 6,7\n"
    "2 7  8 6\n";
/* We make one mutation for each tree, over unary nodes if this exist */
const char *unary_ex_sites =
    "1.0    0\n"
    "4.5    0\n"
    "8.5    0\n";
const char *unary_ex_mutations =
    "0    2   1\n"
    "1    6   1\n"
    "2    5   1\n";

/* Simple utilities to parse text so we can write declaritive
 * tests. This is not intended as a robust general input mechanism.
 */

static void
parse_nodes(const char *text, node_table_t *node_table)
{
    int ret;
    size_t c, k;
    size_t MAX_LINE = 1024;
    char line[MAX_LINE];
    const char *whitespace = " \t";
    char *p;
    double time;
    int flags, population;
    char *name;

    c = 0;
    while (text[c] != '\0') {
        /* Fill in the line */
        k = 0;
        while (text[c] != '\n' && text[c] != '\0') {
            CU_ASSERT_FATAL(k < MAX_LINE - 1);
            line[k] = text[c];
            c++;
            k++;
        }
        if (text[c] == '\n') {
            c++;
        }
        line[k] = '\0';
        p = strtok(line, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        flags = atoi(p);
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        time = atof(p);
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        population = atoi(p);
        p = strtok(NULL, whitespace);
        if (p == NULL) {
            name = "";
        } else {
            name = p;
        }
        ret = node_table_add_row(node_table, flags, time, population, name);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
}

static void
parse_edgesets(const char *text, edgeset_table_t *edgeset_table)
{
    int ret;
    size_t c, k;
    size_t MAX_LINE = 1024;
    size_t MAX_CHILDREN = 1024;
    char line[MAX_LINE], sub_line[MAX_LINE];
    const char *whitespace = " \t";
    char *p, *q;
    double left, right;
    node_id_t parent, children[MAX_CHILDREN];
    uint32_t num_children;

    c = 0;
    while (text[c] != '\0') {
        /* Fill in the line */
        k = 0;
        while (text[c] != '\n' && text[c] != '\0') {
            CU_ASSERT_FATAL(k < MAX_LINE - 1);
            line[k] = text[c];
            c++;
            k++;
        }
        if (text[c] == '\n') {
            c++;
        }
        line[k] = '\0';
        p = strtok(line, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        left = atof(p);
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        right = atof(p);
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        parent = atoi(p);
        num_children = 0;
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p != NULL);

        num_children = 1;
        q = p;
        while (*q != '\0') {
            if (*q == ',') {
                num_children++;
            }
            q++;
        }
        CU_ASSERT_FATAL(num_children >= 1);
        CU_ASSERT_FATAL(num_children < MAX_CHILDREN);
        strncpy(sub_line, p, MAX_LINE);
        q = strtok(sub_line, ",");
        for (k = 0; k < num_children; k++) {
            CU_ASSERT_FATAL(q != NULL);
            children[k] = atoi(q);
            q = strtok(NULL, ",");
        }
        CU_ASSERT_FATAL(q == NULL);
        ret = edgeset_table_add_row(edgeset_table, left, right, parent,
                children, num_children);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
}

static void
parse_sites(const char *text, site_table_t *site_table)
{
    int ret;
    size_t c, k;
    size_t MAX_LINE = 1024;
    char line[MAX_LINE];
    double position;
    char ancestral_state[MAX_LINE];
    const char *whitespace = " \t";
    char *p;

    c = 0;
    while (text[c] != '\0') {
        /* Fill in the line */
        k = 0;
        while (text[c] != '\n' && text[c] != '\0') {
            CU_ASSERT_FATAL(k < MAX_LINE - 1);
            line[k] = text[c];
            c++;
            k++;
        }
        if (text[c] == '\n') {
            c++;
        }
        line[k] = '\0';
        p = strtok(line, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        position = atof(p);
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        strncpy(ancestral_state, p, MAX_LINE);
        ret = site_table_add_row(site_table, position, ancestral_state,
                strlen(ancestral_state));
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
}

static void
parse_mutations(const char *text, mutation_table_t *mutation_table)
{
    int ret;
    size_t c, k;
    size_t MAX_LINE = 1024;
    char line[MAX_LINE];
    const char *whitespace = " \t";
    char *p;
    node_id_t node;
    site_id_t site;
    char derived_state[MAX_LINE];

    c = 0;
    while (text[c] != '\0') {
        /* Fill in the line */
        k = 0;
        while (text[c] != '\n' && text[c] != '\0') {
            CU_ASSERT_FATAL(k < MAX_LINE - 1);
            line[k] = text[c];
            c++;
            k++;
        }
        if (text[c] == '\n') {
            c++;
        }
        line[k] = '\0';
        p = strtok(line, whitespace);
        site = atoi(p);
        CU_ASSERT_FATAL(p != NULL);
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        node = atoi(p);
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        strncpy(derived_state, p, MAX_LINE);
        ret = mutation_table_add_row(mutation_table, site, node, derived_state,
                strlen(derived_state));
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
}

static void
tree_sequence_from_text(tree_sequence_t *ts, const char *nodes, const char *edgesets,
        const char *migrations, const char *sites, const char *mutations,
        const char *provenance)
{
    int ret;
    node_table_t node_table;
    edgeset_table_t edgeset_table;
    mutation_table_t mutation_table;
    site_table_t site_table;
    migration_table_t migration_table;
    size_t default_size_increment = 1024;

    CU_ASSERT_FATAL(ts != NULL);
    CU_ASSERT_FATAL(nodes != NULL);
    CU_ASSERT_FATAL(edgesets != NULL);

    ret = node_table_alloc(&node_table, default_size_increment,
            default_size_increment);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = edgeset_table_alloc(&edgeset_table, default_size_increment,
            default_size_increment);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = site_table_alloc(&site_table, default_size_increment, default_size_increment);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutation_table_alloc(&mutation_table, default_size_increment,
            default_size_increment);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = migration_table_alloc(&migration_table, default_size_increment);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    parse_nodes(nodes, &node_table);
    parse_edgesets(edgesets, &edgeset_table);
    if (sites != NULL) {
        parse_sites(sites, &site_table);
    }
    if (mutations != NULL) {
        parse_mutations(mutations, &mutation_table);
    }
    ret = tree_sequence_initialise(ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(ts, &node_table, &edgeset_table,
            &migration_table, &site_table, &mutation_table, 0, NULL);
    /* tree_sequence_print_state(ts, stdout); */
    /* printf("ret = %s\n", msp_strerror(ret)); */
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    node_table_free(&node_table);
    edgeset_table_free(&edgeset_table);
    site_table_free(&site_table);
    mutation_table_free(&mutation_table);
    migration_table_free(&migration_table);
}

static void
verify_nodes_equal(node_t *n1, node_t *n2)
{
    double eps = 1e-6;

    CU_ASSERT_DOUBLE_EQUAL_FATAL(n1->time, n1->time, eps);
    CU_ASSERT_EQUAL_FATAL(n1->population, n2->population);
    CU_ASSERT_EQUAL_FATAL(n1->flags, n2->flags);
    CU_ASSERT_FATAL(n1->name != NULL);
    CU_ASSERT_FATAL(n2->name != NULL);
    CU_ASSERT_STRING_EQUAL_FATAL(n1->name, n2->name);
}

static void
verify_edgesets_equal(edgeset_t *r1, edgeset_t *r2, double scale)
{
    uint32_t j;
    double eps = 1e-6;

    CU_ASSERT_DOUBLE_EQUAL_FATAL(r1->left * scale, r2->left, eps);
    CU_ASSERT_DOUBLE_EQUAL_FATAL(r1->right * scale, r2->right, eps);
    CU_ASSERT_EQUAL_FATAL(r1->children_length, r2->children_length);
    for (j = 0; j < r1->children_length; j++) {
        CU_ASSERT_EQUAL(r1->children[j], r2->children[j]);
    }
}

static void
verify_migrations_equal(migration_t *r1, migration_t *r2, double scale)
{
    double eps = 1e-6;

    CU_ASSERT_DOUBLE_EQUAL_FATAL(r1->left * scale, r2->left, eps);
    CU_ASSERT_DOUBLE_EQUAL_FATAL(r1->right * scale, r2->right, eps);
    CU_ASSERT_DOUBLE_EQUAL_FATAL(r1->time, r2->time, eps);
    CU_ASSERT_EQUAL_FATAL(r1->node, r2->node);
    CU_ASSERT_EQUAL_FATAL(r1->source, r2->source);
    CU_ASSERT_EQUAL_FATAL(r1->dest, r2->dest);
}

static void
verify_hapgen(tree_sequence_t *ts)
{
    int ret;
    hapgen_t hapgen;
    char *haplotype;
    size_t sample_size = tree_sequence_get_sample_size(ts);
    size_t num_sites = tree_sequence_get_num_sites(ts);
    size_t j;

    ret = hapgen_alloc(&hapgen, ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    hapgen_print_state(&hapgen, _devnull);

    for (j = 0; j < sample_size; j++) {
        ret = hapgen_get_haplotype(&hapgen, j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(strlen(haplotype), num_sites);
    }
    for (j = sample_size; j < sample_size + 10; j++) {
        ret = hapgen_get_haplotype(&hapgen, j, &haplotype);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
    }
    ret = hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
}

static void
verify_vargen(tree_sequence_t *ts)
{
    int ret;
    vargen_t vargen;
    site_t *site;
    size_t sample_size = tree_sequence_get_sample_size(ts);
    size_t num_sites = tree_sequence_get_num_sites(ts);
    char *genotypes = malloc((sample_size + 1) * sizeof(char));
    size_t j, k;

    CU_ASSERT_FATAL(genotypes != NULL);
    ret = vargen_alloc(&vargen, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    vargen_print_state(&vargen, _devnull);
    j = 0;
    while ((ret = vargen_next(&vargen, &site, genotypes)) == 1) {
        CU_ASSERT_EQUAL(site->id, j);
        for (k = 0; k < sample_size; k++) {
            CU_ASSERT(genotypes[k] == 0 || genotypes[k] == 1);
        }
        j++;
    }
    if (ret != 0) {
        tree_sequence_print_state(ts, stdout);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(j, num_sites);
    CU_ASSERT_EQUAL_FATAL(vargen_next(&vargen, &site, genotypes), 0);
    ret = vargen_free(&vargen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = vargen_alloc(&vargen, ts, MSP_GENOTYPES_AS_CHAR);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    vargen_print_state(&vargen, _devnull);
    j = 0;
    genotypes[sample_size] = '\0';
    while ((ret = vargen_next(&vargen, &site, genotypes)) == 1) {
        CU_ASSERT_EQUAL(site->id, j);
        for (k = 0; k < sample_size; k++) {
            CU_ASSERT(genotypes[k] == '0' || genotypes[k] == '1');
        }
        CU_ASSERT_EQUAL_FATAL(genotypes[sample_size], '\0');
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(j, num_sites);
    CU_ASSERT_EQUAL_FATAL(vargen_next(&vargen, &site, genotypes), 0);
    ret = vargen_free(&vargen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    free(genotypes);
}

static void
verify_stats(tree_sequence_t *ts)
{
    int ret;
    uint32_t sample_size = tree_sequence_get_sample_size(ts);
    node_id_t *samples = malloc(sample_size * sizeof(node_id_t));
    uint32_t j;
    double pi;

    CU_ASSERT_FATAL(samples != NULL);

    ret = tree_sequence_get_pairwise_diversity(ts, NULL, 0, &pi);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = tree_sequence_get_pairwise_diversity(ts, NULL, 1, &pi);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = tree_sequence_get_pairwise_diversity(ts, NULL, sample_size + 1, &pi);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    for (j = 0; j < sample_size; j++) {
        samples[j] = j;
    }
    for (j = 2; j < sample_size; j++) {
        ret = tree_sequence_get_pairwise_diversity(ts, samples, j, &pi);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_TRUE(pi >= 0);
    }
    free(samples);
}

static void
verify_simplify_properties(tree_sequence_t *ts, tree_sequence_t *subset,
        node_id_t *samples, uint32_t num_samples)
{
    int ret;
    node_t n1, n2;
    sparse_tree_t full_tree, subset_tree;
    site_t *tree_sites;
    list_len_t tree_sites_length;
    uint32_t j, k;
    node_id_t u, mrca1, mrca2;
    double tmrca1, tmrca2;
    size_t total_sites;

    CU_ASSERT_EQUAL(
        tree_sequence_get_sequence_length(ts),
        tree_sequence_get_sequence_length(subset));
    CU_ASSERT(
        tree_sequence_get_num_nodes(ts) >= tree_sequence_get_num_nodes(subset));
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(subset), num_samples);

    /* Check the sample properties */
    for (j = 0; j < num_samples; j++) {
        ret = tree_sequence_get_node(ts, samples[j], &n1);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tree_sequence_get_node(subset, j, &n2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(n1.population, n2.population);
        CU_ASSERT_EQUAL_FATAL(n1.time, n2.time);
        CU_ASSERT_EQUAL_FATAL(n1.flags, n2.flags);
        /* TODO compare name */
    }
    /* Check the pairwise MRCAs */
    ret = sparse_tree_alloc(&full_tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_alloc(&subset_tree, subset, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_first(&full_tree);
    CU_ASSERT_EQUAL(ret, 1);
    ret = sparse_tree_first(&subset_tree);
    CU_ASSERT_EQUAL(ret, 1);

    total_sites = 0;
    while (1) {
        while (full_tree.right <= subset_tree.right) {
            for (j = 0; j < num_samples; j++) {
                for (k = j + 1; k < num_samples; k++) {
                    ret = sparse_tree_get_mrca(&full_tree, samples[j], samples[k], &mrca1);
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                    ret = sparse_tree_get_time(&full_tree, mrca1, &tmrca1);
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                    ret = sparse_tree_get_mrca(&subset_tree, j, k, &mrca2);
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                    ret = sparse_tree_get_time(&subset_tree, mrca2, &tmrca2);
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                    CU_ASSERT_EQUAL(tmrca1, tmrca2);
                }
            }
            ret = sparse_tree_next(&full_tree);
            CU_ASSERT_FATAL(ret >= 0);
            if (ret != 1) {
                break;
            }
        }
        /* Check the sites in this tree */
        ret = sparse_tree_get_sites(&subset_tree, &tree_sites, &tree_sites_length);
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        for (j = 0; j < tree_sites_length; j++) {
            CU_ASSERT(subset_tree.left <= tree_sites[j].position);
            CU_ASSERT(tree_sites[j].position < subset_tree.right);
            for (k = 0; k < tree_sites[j].mutations_length; k++) {
                ret = sparse_tree_get_parent(&subset_tree,
                        tree_sites[j].mutations[k].node, &u);
                CU_ASSERT_EQUAL(ret, 0);
                CU_ASSERT_FATAL(u != MSP_NULL_NODE);
            }
            total_sites++;
        }
        ret = sparse_tree_next(&subset_tree);
        if (ret != 1) {
            break;
        }
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_sites(subset), total_sites);

    sparse_tree_free(&subset_tree);
    sparse_tree_free(&full_tree);

    /* Check some other operations on the subset */
    if (tree_sequence_get_alphabet(subset) == MSP_ALPHABET_BINARY) {
        /* TODO the binary alphabet is currently equivalent to strict
         * infinite sites, so if have a site with a 1 ancestral state
         * it will have the CHAR alphabet, which vargen doesn't support.
         */
        verify_vargen(subset);
    }
    verify_hapgen(subset);
}

static void
verify_simplify(tree_sequence_t *ts)
{
    int ret;
    uint32_t n = tree_sequence_get_sample_size(ts);
    uint32_t sample_sizes[] = {2, 3, n / 2, n - 1, n};
    size_t j;
    node_id_t *sample = malloc(n * sizeof(node_id_t));
    tree_sequence_t subset;
    int flags = MSP_FILTER_INVARIANT_SITES;

    CU_ASSERT_FATAL(sample != NULL);
    for (j = 0; j < n; j++) {
        sample[j] = j;
    }
    ret = tree_sequence_initialise(&subset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < sizeof(sample_sizes) / sizeof(uint32_t); j++) {
        if (sample_sizes[j] > 1 && sample_sizes[j] <= n) {
            ret = tree_sequence_simplify(ts, sample, sample_sizes[j], flags, &subset);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            verify_simplify_properties(ts, &subset, sample, sample_sizes[j]);
        }
    }
    tree_sequence_free(&subset);
    free(sample);
}

static void
verify_simulator_tree_sequence_equality(msp_t *msp, tree_sequence_t *tree_seq,
        mutgen_t *mutgen, double scale)
{
    int ret;
    uint32_t sample_size = msp_get_sample_size(msp);
    edgeset_t ts_edgeset, sim_edgeset;
    coalescence_record_t *sim_records;
    migration_t *sim_mig_records, ts_mig_record;
    uint32_t j;
    size_t num_edgesets, num_migrations;
    node_t node;
    sample_t *samples;

    CU_ASSERT_EQUAL_FATAL(
            tree_sequence_get_num_edgesets(tree_seq),
            msp_get_num_coalescence_records(msp));
    CU_ASSERT_EQUAL_FATAL(
            tree_sequence_get_sample_size(tree_seq),
            msp_get_sample_size(msp));
    CU_ASSERT_EQUAL_FATAL(
            tree_sequence_get_num_migrations(tree_seq),
            msp_get_num_migrations(msp));
    CU_ASSERT_FATAL(tree_sequence_get_num_nodes(tree_seq) >= sample_size);
    /* CU_ASSERT_EQUAL_FATAL( */
    /*         tree_sequence_get_num_mutations(tree_seq), */
    /*         mutgen_get_num_mutations(mutgen)); */
    ret = msp_get_coalescence_records(msp, &sim_records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    num_edgesets = msp_get_num_coalescence_records(msp);
    ret = msp_get_samples(msp, &samples);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < num_edgesets; j++) {
        ret = tree_sequence_get_edgeset(tree_seq, j, &ts_edgeset);
        CU_ASSERT_EQUAL(ret, 0);
        sim_edgeset.left = sim_records[j].left;
        sim_edgeset.right = sim_records[j].right;
        sim_edgeset.parent = sim_records[j].node;
        sim_edgeset.children_length = sim_records[j].num_children;
        sim_edgeset.children = sim_records[j].children;
        verify_edgesets_equal(&sim_edgeset, &ts_edgeset, scale);
    }
    for (j = num_edgesets; j < num_edgesets + 10; j++) {
        ret = tree_sequence_get_edgeset(tree_seq, j, &ts_edgeset);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
    }

    ret = msp_get_migrations(msp, &sim_mig_records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    num_migrations = msp_get_num_migrations(msp);
    for (j = 0; j < num_migrations; j++) {
        ret = tree_sequence_get_migration(tree_seq, j, &ts_mig_record);
        CU_ASSERT_EQUAL(ret, 0);
        verify_migrations_equal(&sim_mig_records[j], &ts_mig_record, scale);
    }
    for (j = num_migrations; j < num_migrations + 10; j++) {
        ret = tree_sequence_get_migration(tree_seq, j, &ts_mig_record);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
    }
    for (j = 0; j < sample_size; j++) {
        ret = tree_sequence_get_node(tree_seq, j, &node);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(node.population, samples[j].population_id);
        CU_ASSERT_EQUAL(node.time, samples[j].time);
    }
    mutgen_print_state(mutgen, _devnull);
    tree_sequence_print_state(tree_seq, _devnull);
}

/* Utility function to return a tree sequence for testing. It is the
 * callers responsilibility to free all memory.
 */
static tree_sequence_t *
get_example_tree_sequence(uint32_t sample_size,
        uint32_t num_historical_samples, uint32_t num_loci,
        double sequence_length, double scaled_recombination_rate,
        double mutation_rate, uint32_t num_bottlenecks,
        bottleneck_desc_t *bottlenecks, int alphabet)
{
    int ret;
    msp_t *msp = malloc(sizeof(msp_t));
    sample_t *samples = malloc(sample_size * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tree_sequence_t *tree_seq = malloc(sizeof(tree_sequence_t));
    recomb_map_t *recomb_map = malloc(sizeof(recomb_map_t));
    mutgen_t *mutgen = malloc(sizeof(mutgen_t));
    node_table_t *nodes = malloc(sizeof(node_table_t));
    edgeset_table_t *edgesets = malloc(sizeof(edgeset_table_t));
    migration_table_t *migrations= malloc(sizeof(migration_table_t));
    site_table_t *sites = malloc(sizeof(site_table_t));
    mutation_table_t *mutations = malloc(sizeof(mutation_table_t));
    char *provenance[] = {"get_example_tree_sequence"};
    uint32_t j;
    size_t num_populations = 3;
    double migration_matrix[] = {
        0.0, 1.0, 1.0,
        1.0, 0.0, 1.0,
        1.0, 1.0, 0.0
    };
    double positions[] = {0.0, 0.0};
    double rates[] = {0.0, 0.0};

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    CU_ASSERT_FATAL(recomb_map != NULL);
    CU_ASSERT_FATAL(nodes != NULL);
    CU_ASSERT_FATAL(edgesets != NULL);
    CU_ASSERT_FATAL(migrations != NULL);
    CU_ASSERT_FATAL(sites != NULL);
    CU_ASSERT_FATAL(mutations != NULL);
    gsl_rng_set(rng, 1);

    ret = node_table_alloc(nodes, 10, 10);
    CU_ASSERT_EQUAL(ret, 0);
    ret = edgeset_table_alloc(edgesets, 10, 10);
    CU_ASSERT_EQUAL(ret, 0);
    ret = migration_table_alloc(migrations, 10);
    CU_ASSERT_EQUAL(ret, 0);
    ret = site_table_alloc(sites, 10, 10);
    CU_ASSERT_EQUAL(ret, 0);
    ret = mutation_table_alloc(mutations, 10, 10);
    CU_ASSERT_EQUAL(ret, 0);

    ret = mutgen_alloc(mutgen, mutation_rate, rng, alphabet, 10);
    CU_ASSERT_EQUAL(ret, 0);
    /* initialise the samples to zero for the default configuration */
    memset(samples, 0, sample_size * sizeof(sample_t));
    for (j = 0; j < num_historical_samples; j++) {
        samples[j].time = 0.1 * (j + 1);
    }
    ret = msp_alloc(msp, sample_size, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(msp, num_loci);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(msp, scaled_recombination_rate);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_bottlenecks; j++) {
        if (bottlenecks[j].type == SIMPLE_BOTTLENECK) {
            ret = msp_add_simple_bottleneck(msp, bottlenecks[j].time,
                    bottlenecks[j].population_id, bottlenecks[j].parameter);
            CU_ASSERT_EQUAL(ret, 0);
        } else if(bottlenecks[j].type == INSTANTANEOUS_BOTTLENECK) {
            ret = msp_add_instantaneous_bottleneck(msp, bottlenecks[j].time,
                    bottlenecks[j].population_id, bottlenecks[j].parameter);
            CU_ASSERT_EQUAL(ret, 0);
        } else {
            CU_ASSERT_FATAL(0 == 1);
        }
    }
    ret = msp_set_num_populations(msp, num_populations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_migration_matrix(msp, num_populations * num_populations,
            migration_matrix);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_store_migrations(msp, true);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_run(msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    rates[0] = scaled_recombination_rate;
    positions[1] = sequence_length;
    ret = recomb_map_alloc(recomb_map, num_loci, sequence_length,
            positions, rates, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Create the tree_sequence from the state of the simulator.
     * We want to use coalescent time here, so use an Ne of 1/4
     * to cancel scaling factor. */
    ret = tree_sequence_initialise(tree_seq);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_populate_tables(msp, 0.25, recomb_map, nodes, edgesets, migrations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate_tables_tmp(mutgen, nodes, edgesets);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_populate_tables(mutgen, sites, mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(tree_seq, nodes, edgesets, migrations,
            sites, mutations, 1, provenance);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_simulator_tree_sequence_equality(msp, tree_seq, mutgen,
            sequence_length / num_loci);

    gsl_rng_free(rng);
    free(samples);
    msp_free(msp);
    free(msp);
    recomb_map_free(recomb_map);
    free(recomb_map);
    mutgen_free(mutgen);
    free(mutgen);
    node_table_free(nodes);
    edgeset_table_free(edgesets);
    mutation_table_free(mutations);
    site_table_free(sites);
    migration_table_free(migrations);
    free(nodes);
    free(edgesets);
    free(migrations);
    free(mutations);
    free(sites);
    return tree_seq;
}

tree_sequence_t **
get_example_nonbinary_tree_sequences(void)
{
    size_t max_examples = 1024;
    tree_sequence_t **ret = malloc(max_examples * sizeof(tree_sequence_t *));
    bottleneck_desc_t bottlenecks[] = {
        {SIMPLE_BOTTLENECK, 0.1, 0, 0.5},
        {INSTANTANEOUS_BOTTLENECK, 0.4, 0, 10.0},
    };
    bottleneck_desc_t other_bottlenecks[] = {
        {SIMPLE_BOTTLENECK, 0.1, 0, 0.1},
        {SIMPLE_BOTTLENECK, 0.15, 0, 0.75},
        {INSTANTANEOUS_BOTTLENECK, 0.2, 0, 0.1},
    };

    CU_ASSERT_FATAL(ret != NULL);
    ret[0] = get_example_tree_sequence(100, 0, 100, 100.0, 10.0, 1.0,
            1, bottlenecks, MSP_ALPHABET_BINARY);
    ret[1] = get_example_tree_sequence(10, 2, 100, 10.0, 1.0, 2.0,
            1, bottlenecks, MSP_ALPHABET_BINARY);
    ret[2] = get_example_tree_sequence(500, 10, 10, 1000.0, 0.5, 3.0,
            2, bottlenecks, MSP_ALPHABET_BINARY);
    ret[3] = get_example_tree_sequence(100, 0, 100, 1.0, 1.0, 0.0,
            3, other_bottlenecks, MSP_ALPHABET_BINARY);
    ret[4] = NULL;
    return ret;
}

tree_sequence_t *
make_recurrent_and_back_mutations_copy(tree_sequence_t *ts)
{
    int ret;
    size_t num_provenance_strings;
    size_t alloc_size = 8192;
    char **provenance_strings;
    tree_sequence_t *new_ts = malloc(sizeof(tree_sequence_t));
    sparse_tree_t tree;
    node_table_t nodes;
    edgeset_table_t edgesets;
    migration_table_t migrations;
    mutation_table_t mutations;
    site_table_t sites;
    node_id_t *stack;
    node_id_t u;
    site_id_t site_id;
    char *state = NULL;
    int stack_top = 0;

    CU_ASSERT_FATAL(new_ts != NULL);
    ret = node_table_alloc(&nodes, alloc_size, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = edgeset_table_alloc(&edgesets, alloc_size, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = migration_table_alloc(&migrations, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = site_table_alloc(&sites, alloc_size, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutation_table_alloc(&mutations, alloc_size, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = sparse_tree_alloc(&tree, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    state = malloc(tree_sequence_get_num_nodes(ts) * sizeof(char));
    CU_ASSERT_FATAL(state != NULL);

    stack = tree.stack1;
    site_id = 0;
    for (ret = sparse_tree_first(&tree); ret == 1; ret = sparse_tree_next(&tree)) {
        ret = site_table_add_row(&sites, tree.left, "0", 1);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        /* Traverse down the tree */
        stack_top = 0;
        stack[0] = tree.root;
        state[tree.root] = 0;
        while (stack_top >= 0) {
            u = stack[stack_top];
            stack_top--;
            if (u != tree.root) {
                state[u] = (state[u] + 1) % 2;
                ret = mutation_table_add_row(&mutations, site_id, u,
                        state[u] == 0? "0": "1", 1);
                CU_ASSERT_EQUAL_FATAL(ret, 0);
            }
            /* To ensure that the mutations are sorted in time order, we only
             * traverse down the left-most path to the leaf. This means we
             * don't really need a stack at all, but we may extend this to a full
             * traversal in the future. */
            if (tree.num_children[u] > 0) {
                stack_top++;
                stack[stack_top] = tree.children[u][0];
                state[tree.children[u][0]] = state[u];
            }
        }
        site_id++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tree_sequence_dump_tables_tmp(ts, &nodes, &edgesets,
            &migrations, NULL, NULL, &num_provenance_strings,
            &provenance_strings);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_initialise(new_ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(new_ts, &nodes, &edgesets, &migrations,
            &sites, &mutations, num_provenance_strings, provenance_strings);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    node_table_free(&nodes);
    edgeset_table_free(&edgesets);
    migration_table_free(&migrations);
    site_table_free(&sites);
    mutation_table_free(&mutations);
    sparse_tree_free(&tree);
    free(state);
    return new_ts;
}

tree_sequence_t **
get_example_tree_sequences(int include_nonbinary)
{
    size_t max_examples = 1024;
    int j, k;
    tree_sequence_t **ret = malloc(max_examples * sizeof(tree_sequence_t *));
    tree_sequence_t **nonbinary = NULL;
    CU_ASSERT_FATAL(ret != NULL);

    ret[0] = get_example_tree_sequence(10, 0, 100, 100.0, 1.0, 1.0, 0, NULL,
            MSP_ALPHABET_BINARY);
    ret[1] = get_example_tree_sequence(2, 0, 1, 0.1, 1.0, 1.0, 0, NULL,
            MSP_ALPHABET_BINARY);
    ret[2] = get_example_tree_sequence(3, 0, 3, 10.0, 10.0, 0.0, 0, NULL,
            MSP_ALPHABET_BINARY);
    ret[3] = get_example_tree_sequence(10, 0, UINT32_MAX, 10.0,
            9.31322575049e-08, 10.0, 0, NULL, MSP_ALPHABET_BINARY);
    ret[4] = make_recurrent_and_back_mutations_copy(ret[0]);
    k = 5;
    if (include_nonbinary) {
        nonbinary = get_example_nonbinary_tree_sequences();
        for (j = 0; nonbinary[j] != NULL; j++) {
            ret[k] = nonbinary[j];
            k++;
        }
        free(nonbinary);
    }
    ret[k] = NULL;
    return ret;
}


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
            /* Just make sure that we're seeing the same values even when
             * we expand.
             */
            CU_ASSERT(fenwick_expand(&t, 1) == 0);
        }
        CU_ASSERT(fenwick_free(&t) == 0);
    }
}

static void
verify_vcf_converter(tree_sequence_t *ts, unsigned int ploidy)
{
    int ret;
    char *str = NULL;
    vcf_converter_t vc;
    unsigned int num_variants;

    ret = vcf_converter_alloc(&vc, ts, ploidy);
    CU_ASSERT_FATAL(ret ==  0);
    vcf_converter_print_state(&vc, _devnull);
    ret = vcf_converter_get_header(&vc, &str);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_NSTRING_EQUAL("##", str, 2);
    num_variants = 0;
    while ((ret = vcf_converter_next(&vc, &str)) == 1) {
        CU_ASSERT_NSTRING_EQUAL("1\t", str, 2);
        num_variants++;
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(num_variants == tree_sequence_get_num_mutations(ts));
    vcf_converter_free(&vc);
}

static void
test_vcf(void)
{
    int ret;
    unsigned int ploidy;
    vcf_converter_t *vc = malloc(sizeof(vcf_converter_t));
    tree_sequence_t *ts = get_example_tree_sequence(10, 0, 100, 100.0, 1.0, 1.0,
            0, NULL, MSP_ALPHABET_BINARY);

    CU_ASSERT_FATAL(ts != NULL);
    CU_ASSERT_FATAL(vc != NULL);

    ret = vcf_converter_alloc(vc, ts, 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = vcf_converter_alloc(vc, ts, 3);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = vcf_converter_alloc(vc, ts, 11);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    for (ploidy = 1; ploidy < 3; ploidy++) {
        verify_vcf_converter(ts, ploidy);
    }

    free(vc);
    tree_sequence_free(ts);
    free(ts);
}

static void
test_vcf_no_mutations(void)
{
    int ret;
    char *str = NULL;
    vcf_converter_t *vc = malloc(sizeof(vcf_converter_t));
    tree_sequence_t *ts = get_example_tree_sequence(100, 0, 1, 1.0, 0.0, 0.0, 0, NULL,
            MSP_ALPHABET_BINARY);

    CU_ASSERT_FATAL(ts != NULL);
    CU_ASSERT_FATAL(vc != NULL);
    CU_ASSERT_EQUAL_FATAL(tree_sequence_get_num_mutations(ts), 0);

    ret = vcf_converter_alloc(vc, ts, 1);
    CU_ASSERT_FATAL(ret ==  0);
    vcf_converter_print_state(vc, _devnull);
    ret = vcf_converter_get_header(vc, &str);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_NSTRING_EQUAL("##", str, 2);
    ret = vcf_converter_next(vc, &str);
    CU_ASSERT_EQUAL(ret, 0);
    vcf_converter_free(vc);

    free(vc);
    tree_sequence_free(ts);
    free(ts);
}

static void
test_simple_recomb_map(void)
{
    int ret;
    recomb_map_t recomb_map;
    double positions[] = {0.0, 1.0};
    double rates[] = {0.0, 1.0};
    uint32_t num_loci[] = {1, 100, UINT32_MAX};
    size_t j;

    for (j = 0; j < sizeof(num_loci) / sizeof(uint32_t); j++) {
        ret = recomb_map_alloc(&recomb_map, num_loci[j], 1.0,
                positions, rates, 2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        recomb_map_print_state(&recomb_map, _devnull);
        CU_ASSERT_EQUAL(recomb_map_get_num_loci(&recomb_map), num_loci[j]);
        CU_ASSERT_EQUAL(recomb_map_get_size(&recomb_map), 2);
        CU_ASSERT_EQUAL(recomb_map_get_sequence_length(&recomb_map), 1.0);
        CU_ASSERT_EQUAL(
                recomb_map_genetic_to_phys(&recomb_map, num_loci[j]), 1.0);
        CU_ASSERT_EQUAL(
                recomb_map_phys_to_genetic(&recomb_map, 1.0), num_loci[j]);
        ret = recomb_map_free(&recomb_map);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
}

static void
test_recomb_map_errors(void)
{
    int ret;
    recomb_map_t recomb_map;
    double positions[] = {0.0, 1.0, 2.0};
    double rates[] = {1.0, 2.0, 0.0};
    double values[] = {0.0, 1.0};

    ret = recomb_map_alloc(&recomb_map, 10, 1.0, positions, rates, 0);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);

    ret = recomb_map_alloc(&recomb_map, 10, 1.0, positions, rates, 1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);

    ret = recomb_map_alloc(&recomb_map, 0, 1.0, positions, rates, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);

    ret = recomb_map_alloc(&recomb_map, 10, 2.0, positions, rates, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);

    positions[0] = 1;
    ret = recomb_map_alloc(&recomb_map, 10, 1.0, positions, rates, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);
    positions[0] = 0;

    positions[1] = 3.0;
    ret = recomb_map_alloc(&recomb_map, 10, 2.0, positions, rates, 3);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);
    positions[1] = 1.0;

    ret = recomb_map_alloc(&recomb_map, 10, 2.0, positions, rates, 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = recomb_map_genetic_to_phys_bulk(&recomb_map, values, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* values must be increasing */
    values[0] = 2.0;
    ret = recomb_map_genetic_to_phys_bulk(&recomb_map, values, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_GENERIC);

    /* values must be <= num_loci */
    values[0] = 1000;
    ret = recomb_map_genetic_to_phys_bulk(&recomb_map, values, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_GENERIC);

    recomb_map_free(&recomb_map);
}

static void
verify_recomb_map(uint32_t num_loci, double length, double *positions,
        double *rates, size_t size)
{

    int ret;
    recomb_map_t recomb_map;
    double total_rate, x, y, z;
    size_t j;
    size_t num_checks = 1000;
    double eps = 1e-6;
    double *ret_rates, *ret_positions, *bulk_x;

    ret_rates = malloc(size * sizeof(double));
    ret_positions = malloc(size * sizeof(double));
    bulk_x = malloc(num_checks * sizeof(double));

    CU_ASSERT_FATAL(ret_rates != NULL);
    CU_ASSERT_FATAL(ret_positions != NULL);

    ret = recomb_map_alloc(&recomb_map, num_loci, length,
           positions, rates, size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    recomb_map_print_state(&recomb_map, _devnull);

    CU_ASSERT_EQUAL(recomb_map_get_num_loci(&recomb_map), num_loci);
    CU_ASSERT_EQUAL(recomb_map_get_size(&recomb_map), size);
    CU_ASSERT_EQUAL(recomb_map_genetic_to_phys(&recomb_map, 0), 0);
    CU_ASSERT_EQUAL(
            recomb_map_genetic_to_phys(&recomb_map, num_loci), length);
    total_rate = recomb_map_get_total_recombination_rate(&recomb_map);
    CU_ASSERT_DOUBLE_EQUAL(
            total_rate / (num_loci - 1),
            recomb_map_get_per_locus_recombination_rate(&recomb_map), eps)

    for (j = 0; j < num_checks; j++) {
        x = j * (num_loci / num_checks);
        bulk_x[j] = x;
        y = recomb_map_genetic_to_phys(&recomb_map, x);
        CU_ASSERT_TRUE(0 <= y && y <= length);
        z = recomb_map_phys_to_genetic(&recomb_map, y);
        CU_ASSERT_DOUBLE_EQUAL(x, z, eps);
    }
    ret = recomb_map_genetic_to_phys_bulk(&recomb_map, bulk_x, num_checks);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_checks; j++) {
        x = j * (num_loci / num_checks);
        y = recomb_map_genetic_to_phys(&recomb_map, x);
        CU_ASSERT_DOUBLE_EQUAL(bulk_x[j], y, 1e-12);
    }
    ret = recomb_map_get_positions(&recomb_map, ret_positions);
    CU_ASSERT_EQUAL(ret, 0);
    ret = recomb_map_get_rates(&recomb_map, ret_rates);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < size; j++) {
        CU_ASSERT_EQUAL(ret_rates[j], rates[j]);
        CU_ASSERT_EQUAL(ret_positions[j], positions[j]);
    }
    ret = recomb_map_free(&recomb_map);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    free(ret_rates);
    free(ret_positions);
    free(bulk_x);
}

static void
test_recomb_map_examples(void)
{
    double p1[] =  {0, 0.05019838393314813, 0.36933662489552865, 1};
    double r1[] = {3.5510784169955434, 4.184964179610539, 3.800808140657212, 0};
    double p2[] =  {0, 0.125, 0.875, 1, 4, 8, 16};
    double r2[] = {0.1, 6.0, 3.333, 2.1, 0.0, 2.2, 0};

    verify_recomb_map(2, 1.0, p1, r1, 4);
    verify_recomb_map(1000, 1.0, p1, r1, 4);
    verify_recomb_map(UINT32_MAX, 1.0, p1, r1, 4);

    verify_recomb_map(2, 16.0, p2, r2, 7);
    verify_recomb_map(100, 16.0, p2, r2, 7);
    verify_recomb_map(UINT32_MAX, 16.0, p2, r2, 7);
}

static void
test_single_locus_two_populations(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = {{0, 0.0}, {0, 0.0}, {1, 40.0}};
    coalescence_record_t *coalescence_records;
    migration_t *migrations;
    size_t num_coalescence_records, num_migrations;
    uint32_t n = 3;
    double t0 = 30.0;
    double t1 = 30.5;
    double t2 = 40.5;

    CU_ASSERT_FATAL(rng != NULL);

    ret = msp_alloc(&msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_populations(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_store_migrations(&msp, true);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_mass_migration(&msp, t0, 0, 1, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_mass_migration(&msp, t1, 1, 0, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_mass_migration(&msp, t2, 1, 0, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp);
    msp_print_state(&msp, _devnull);
    num_coalescence_records = msp_get_num_coalescence_records(&msp);
    CU_ASSERT_EQUAL_FATAL(num_coalescence_records, 2);
    ret = msp_get_coalescence_records(&msp, &coalescence_records);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(coalescence_records[0].node, 3);
    CU_ASSERT_TRUE(coalescence_records[0].time < 40.0);
    CU_ASSERT_EQUAL(coalescence_records[0].population_id, 0);
    CU_ASSERT_EQUAL(coalescence_records[1].node, 4);
    CU_ASSERT_TRUE(coalescence_records[1].time > 40.5);
    CU_ASSERT_EQUAL(coalescence_records[1].population_id, 0);
    num_migrations = msp_get_num_migrations(&msp);
    CU_ASSERT_EQUAL_FATAL(num_migrations, 3);
    ret = msp_get_migrations(&msp, &migrations);
    CU_ASSERT_EQUAL(migrations[0].time, t0)
    CU_ASSERT_EQUAL(migrations[0].source, 0);
    CU_ASSERT_EQUAL(migrations[0].dest, 1);
    CU_ASSERT_EQUAL(migrations[1].time, t1);
    CU_ASSERT_EQUAL(migrations[1].source, 1);
    CU_ASSERT_EQUAL(migrations[1].dest, 0);
    CU_ASSERT_EQUAL(migrations[2].time, t2);
    CU_ASSERT_EQUAL(migrations[2].source, 1);
    CU_ASSERT_EQUAL(migrations[2].dest, 0);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
}

static void
test_single_locus_many_populations(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    uint32_t num_populations = 500;
    sample_t samples[] = {{0, 0.0}, {num_populations - 1, 0.0}};
    coalescence_record_t *records;
    size_t num_records;
    uint32_t n = 2;

    CU_ASSERT_FATAL(rng != NULL);

    ret = msp_alloc(&msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_populations(&msp, num_populations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_add_mass_migration(&msp, 30.0, 0, num_populations - 1, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp);
    msp_print_state(&msp, _devnull);
    num_records = msp_get_num_coalescence_records(&msp);
    CU_ASSERT_EQUAL_FATAL(num_records, 1);
    ret = msp_get_coalescence_records(&msp, &records);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(records[0].node, 2);
    CU_ASSERT_TRUE(records[0].time > 30.0);
    CU_ASSERT_EQUAL(records[0].population_id, num_populations - 1);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
}

static void
test_single_locus_historical_sample(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = {{0, 0.0}, {0, 10.0}};
    coalescence_record_t *record;
    size_t num_records;
    uint32_t n = 2;

    CU_ASSERT_FATAL(rng != NULL);

    ret = msp_alloc(&msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp);
    msp_print_state(&msp, _devnull);
    num_records = msp_get_num_coalescence_records(&msp);
    CU_ASSERT_EQUAL_FATAL(num_records, 1);
    ret = msp_get_coalescence_records(&msp, &record);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(record->left, 0);
    CU_ASSERT_EQUAL(record->right, 1);
    CU_ASSERT_EQUAL(record->node, 2);
    CU_ASSERT_TRUE(record->time > 10.0);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
}

static void
test_simulator_getters_setters(void)
{
    int ret;
    uint32_t j;
    uint32_t n = 10;
    uint32_t m = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    double migration_matrix[] = {0, 0, 0, 0};
    double matrix[4], growth_rate, initial_size;
    size_t migration_events[4];
    size_t breakpoints[m];
    population_t *population;
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    for (j = 0; j < n; j++) {
        samples[j].time = j;
        samples[j].population_id = j % 2;
    }
    CU_ASSERT_EQUAL(msp_alloc(&msp, 0, NULL, rng), MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);
    samples[0].time = 1.0;
    CU_ASSERT_EQUAL(msp_alloc(&msp, n, samples, rng), MSP_ERR_BAD_SAMPLES);
    msp_free(&msp);
    samples[0].time = -1.0;
    CU_ASSERT_EQUAL(msp_alloc(&msp, n, samples, rng), MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);
    samples[0].time = 0.0;

    ret = msp_alloc(&msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(msp_set_max_memory(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_node_mapping_block_size(&msp, 0),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_segment_block_size(&msp, 0),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_avl_node_block_size(&msp, 0),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_coalescence_record_block_size(&msp, 0),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_migration_block_size(&msp, 0),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_num_loci(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_num_populations(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
            msp_set_scaled_recombination_rate(&msp, -1),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
            msp_set_population_configuration(&msp, -1, 0, 0),
            MSP_ERR_BAD_POPULATION_ID);
    CU_ASSERT_EQUAL(
            msp_set_population_configuration(&msp, 3, 0, 0),
            MSP_ERR_BAD_POPULATION_ID);

    CU_ASSERT_EQUAL(msp_get_model(&msp)->type, MSP_MODEL_HUDSON);

    ret = msp_set_num_populations(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(
            msp_get_population_configuration(&msp, 3, NULL, NULL),
            MSP_ERR_BAD_POPULATION_ID);
    ret = msp_set_population_configuration(&msp, 0, 2, 0.5);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(
            msp_set_migration_matrix(&msp, 0, NULL),
            MSP_ERR_BAD_MIGRATION_MATRIX);
    CU_ASSERT_EQUAL(
            msp_set_migration_matrix(&msp, 3, migration_matrix),
            MSP_ERR_BAD_MIGRATION_MATRIX);
    migration_matrix[0] = 1;
    CU_ASSERT_EQUAL(
            msp_set_migration_matrix(&msp, 4, migration_matrix),
            MSP_ERR_BAD_MIGRATION_MATRIX);
    migration_matrix[0] = 0;
    migration_matrix[1] = -1;
    CU_ASSERT_EQUAL(
            msp_set_migration_matrix(&msp, 4, migration_matrix),
            MSP_ERR_BAD_MIGRATION_MATRIX);

    migration_matrix[1] = 1;
    migration_matrix[2] = 1;
    ret = msp_set_migration_matrix(&msp, 4, migration_matrix);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(&msp, m);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(&msp, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_store_migrations(&msp, true);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_get_population_configuration(&msp, 0,
            &initial_size, &growth_rate);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(initial_size, 2);
    CU_ASSERT_EQUAL(growth_rate, 0.5);
    CU_ASSERT_EQUAL(
            msp_get_population(&msp, 3, NULL),
            MSP_ERR_BAD_POPULATION_ID);
    ret = msp_get_population(&msp, 0, &population);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(population->initial_size, 2);
    CU_ASSERT_EQUAL(population->growth_rate, 0.5);
    CU_ASSERT_EQUAL(population->start_time, 0.0);

    CU_ASSERT_TRUE(msp_get_store_migrations(&msp));
    CU_ASSERT_EQUAL(msp_get_num_avl_node_blocks(&msp), 1);
    CU_ASSERT_EQUAL(msp_get_num_node_mapping_blocks(&msp), 1);
    CU_ASSERT_EQUAL(msp_get_num_segment_blocks(&msp), 1);
    CU_ASSERT_EQUAL(msp_get_num_coalescence_record_blocks(&msp), 1);
    CU_ASSERT_EQUAL(msp_get_num_migration_blocks(&msp), 1);
    CU_ASSERT(msp_get_used_memory(&msp) > 0);
    CU_ASSERT_EQUAL(msp_get_num_populations(&msp), 2);

    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(msp_get_num_breakpoints(&msp), m - 1);
    ret = msp_get_breakpoints(&msp, breakpoints);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < m - 1; j++) {
        CU_ASSERT_EQUAL(breakpoints[j], j + 1);
    }
    ret = msp_get_num_migration_events(&msp, migration_events);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(migration_events[0], 0);
    CU_ASSERT(migration_events[1] > 0);
    CU_ASSERT(migration_events[2] > 0);
    CU_ASSERT_EQUAL(migration_events[3], 0);
    ret = msp_get_migration_matrix(&msp, matrix);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < 4; j++) {
        CU_ASSERT_EQUAL(matrix[j], migration_matrix[j]);
    }
    CU_ASSERT(msp_get_num_common_ancestor_events(&msp) > 0);
    CU_ASSERT(msp_get_num_recombination_events(&msp) > 0);

    free(samples);
    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
}

static void
test_simulator_model_errors(void)
{
    uint32_t n = 10;
    uint32_t j;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    int models[] = {MSP_MODEL_SMC, MSP_MODEL_SMC_PRIME};
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    for (j = 0; j < n; j++) {
        samples[j].time = 0;
        samples[j].population_id = 0;
    }

    CU_ASSERT_EQUAL(msp_alloc(&msp, n, samples, rng), 0);
    CU_ASSERT_EQUAL(msp_get_model(&msp)->type, MSP_MODEL_HUDSON);
    CU_ASSERT_EQUAL(msp_add_simple_bottleneck(&msp, 1, 0, 1), 0);
    CU_ASSERT_EQUAL(msp_add_instantaneous_bottleneck(&msp, 1, 0, 1), 0);
    CU_ASSERT_EQUAL(msp_initialise(&msp), 0);
    CU_ASSERT_EQUAL(msp_run(&msp, DBL_MAX, ULONG_MAX), 0);
    CU_ASSERT_EQUAL(msp_free(&msp), 0);

    CU_ASSERT_EQUAL(msp_alloc(&msp, n, samples, rng), 0);
    CU_ASSERT_EQUAL(msp_add_simple_bottleneck(&msp, 1, 0, 1), 0);
    CU_ASSERT_EQUAL(msp_add_instantaneous_bottleneck(&msp, 1, 0, 1), 0);
    CU_ASSERT_EQUAL(msp_set_simulation_model_non_parametric(&msp, MSP_MODEL_HUDSON),
            MSP_ERR_UNSUPPORTED_OPERATION);
    CU_ASSERT_EQUAL(msp_free(&msp), 0);

    for (j = 0; j < sizeof(models) / sizeof(int); j++) {
        CU_ASSERT_EQUAL(msp_alloc(&msp, n, samples, rng), 0);
        CU_ASSERT_EQUAL(msp_set_simulation_model_non_parametric(&msp, models[j]), 0);
        CU_ASSERT_EQUAL(msp_add_simple_bottleneck(&msp, 1, 0, 1), MSP_ERR_BAD_MODEL);
        CU_ASSERT_EQUAL(msp_add_instantaneous_bottleneck(&msp, 1, 0, 1),
                MSP_ERR_BAD_MODEL);
        CU_ASSERT_EQUAL(msp_free(&msp), 0);
    }

    free(samples);
    gsl_rng_free(rng);
}

static void
test_simulator_demographic_events(void)
{
    int ret;
    uint32_t j;
    uint32_t n = 10;
    uint32_t m = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    double migration_matrix[] = {0, 1, 1, 0};
    double time;
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    for (j = 0; j < n; j++) {
        samples[j].time = j;
        samples[j].population_id = j % 2;
    }
    ret = msp_alloc(&msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_populations(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 1, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 1, 2, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_migration_matrix(&msp, 4, migration_matrix);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(&msp, m);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(&msp, 1.0);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(
        msp_add_mass_migration(&msp, 10, -1, 0, 1),
        MSP_ERR_BAD_POPULATION_ID);
    CU_ASSERT_EQUAL(
        msp_add_mass_migration(&msp, 10, 2, 0, 1),
        MSP_ERR_BAD_POPULATION_ID);
    CU_ASSERT_EQUAL(
        msp_add_mass_migration(&msp, 10, 0, 0, 1),
        MSP_ERR_SOURCE_DEST_EQUAL);
    CU_ASSERT_EQUAL(
        msp_add_mass_migration(&msp, 10, 0, 1, -5),
        MSP_ERR_BAD_PARAM_VALUE);

    CU_ASSERT_EQUAL(
        msp_add_migration_rate_change(&msp, 10, -2, 2.0),
        MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
    CU_ASSERT_EQUAL(
        msp_add_migration_rate_change(&msp, 10, -1, -2.0),
        MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
        msp_add_migration_rate_change(&msp, 10, 3, 2.0),
        MSP_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX);

    CU_ASSERT_EQUAL(
        msp_add_population_parameters_change(&msp, 10, -2, 0, 0),
        MSP_ERR_BAD_POPULATION_ID);
    CU_ASSERT_EQUAL(
        msp_add_population_parameters_change(&msp, 10, -1, -1, 0),
        MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
        msp_add_population_parameters_change(&msp, 10, -1, GSL_NAN, GSL_NAN),
        MSP_ERR_BAD_PARAM_VALUE);

    CU_ASSERT_EQUAL(
        msp_add_simple_bottleneck(&msp, 10, -1, 0),
        MSP_ERR_BAD_POPULATION_ID);
    CU_ASSERT_EQUAL(
        msp_add_simple_bottleneck(&msp, 10, 0, -1),
        MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
        msp_add_simple_bottleneck(&msp, 10, 0, 1.1),
        MSP_ERR_BAD_PARAM_VALUE);

    CU_ASSERT_EQUAL(
        msp_add_instantaneous_bottleneck(&msp, 10, 2, 0),
        MSP_ERR_BAD_POPULATION_ID);
    CU_ASSERT_EQUAL(
        msp_add_simple_bottleneck(&msp, 10, 0, -1),
        MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL_FATAL(
        msp_add_simple_bottleneck(&msp, 10, -1, 0),
        MSP_ERR_BAD_POPULATION_ID);

    ret = msp_add_mass_migration(&msp, 0.1, 0, 1, 0.5);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_migration_rate_change(&msp, 0.2, 1, 2.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_migration_rate_change(&msp, 0.3, -1, 3.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.4, 0, 0.5, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.5, -1, 0.5, 2.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.6, 0, GSL_NAN, 2.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.7, 1, 1, GSL_NAN);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_simple_bottleneck(&msp, 0.8, 0, 0.5);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_instantaneous_bottleneck(&msp, 0.9, 0, 2.0);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(
            msp_add_mass_migration(&msp, 0.1, 0, 1, 0.5),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
    CU_ASSERT_EQUAL(
            msp_add_migration_rate_change(&msp, 0.2, 1, 2.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
    CU_ASSERT_EQUAL(
            msp_add_population_parameters_change(&msp, 0.4, 0, 0.5, 1.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
    CU_ASSERT_EQUAL(
            msp_add_simple_bottleneck(&msp, 0.7, 0, 1.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
    CU_ASSERT_EQUAL(
            msp_add_instantaneous_bottleneck(&msp, 0.8, 0, 1.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);

    CU_ASSERT_EQUAL(
        msp_debug_demography(&msp, &time),
        MSP_ERR_BAD_STATE);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    j = 0;
    do {
        ret = msp_debug_demography(&msp, &time);
        CU_ASSERT_EQUAL(ret, 0);
        msp_print_state(&msp, _devnull);
        j++;
    } while (! gsl_isinf(time));
    CU_ASSERT_EQUAL(j, 10);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(
        msp_run(&msp, DBL_MAX, ULONG_MAX),
        MSP_ERR_BAD_STATE);
    ret = msp_reset(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);

    free(samples);
    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
}

static void
test_single_locus_simulation(void)
{
    int ret;
    int model;
    const char *model_name;
    uint32_t j;
    uint32_t n = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, rng);
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

    model = msp_get_model(msp)->type;
    CU_ASSERT_EQUAL(model, MSP_MODEL_HUDSON);
    model_name = msp_get_model_name(msp);
    CU_ASSERT_STRING_EQUAL(model_name, "hudson");

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
}

static void
test_simulation_memory_limit(void)
{
    int ret;
    uint32_t n = 1000;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_max_memory(msp, 1024 * 1024);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(msp, 10000);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(msp, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_NO_MEMORY);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
}

static void
test_multi_locus_simulation(void)
{
    int ret;
    uint32_t num_events;
    uint32_t n = 100;
    uint32_t m = 100;
    long seed = 10;
    int model;
    int models[] = {MSP_MODEL_HUDSON, MSP_MODEL_SMC, MSP_MODEL_SMC_PRIME};
    double migration_matrix[] = {0, 1, 1, 0};
    size_t migration_events[4];
    const char *model_names[] = {"hudson", "smc", "smc_prime"};
    const char *model_name;
    size_t j;

    for (j = 0; j < sizeof(models) / sizeof(int); j++) {
        sample_t *samples = malloc(n * sizeof(sample_t));
        msp_t *msp = malloc(sizeof(msp_t));
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

        CU_ASSERT_FATAL(msp != NULL);
        CU_ASSERT_FATAL(samples != NULL);
        CU_ASSERT_FATAL(rng != NULL);
        gsl_rng_set(rng, seed);
        memset(samples, 0, n * sizeof(sample_t));
        ret = msp_alloc(msp, n, samples, rng);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_num_populations(msp, 2);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_migration_matrix(msp, 4, migration_matrix);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_store_migrations(msp, true);
        CU_ASSERT_EQUAL(ret, 0);
        /* set all the block sizes to something small to provoke the memory
         * expansions. */
        ret = msp_set_avl_node_block_size(msp, 1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_node_mapping_block_size(msp, 1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_segment_block_size(msp, 1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_coalescence_record_block_size(msp, 1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_migration_block_size(msp, 1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_num_loci(msp, m);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_scaled_recombination_rate(msp, 1.0);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_simulation_model_non_parametric(msp, models[j]);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_initialise(msp);
        CU_ASSERT_EQUAL(ret, 0);

        num_events = 0;
        while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
            msp_verify(msp);
            num_events++;
        }
        CU_ASSERT_EQUAL(ret, 0);
        msp_verify(msp);
        ret = msp_get_num_migration_events(msp, migration_events);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT(num_events > n - 1);
        CU_ASSERT_EQUAL(1 + num_events,
                /* diagonals must be zero here */
                migration_events[1] + migration_events[2] +
                msp_get_num_recombination_events(msp) +
                msp_get_num_common_ancestor_events(msp) +
                msp_get_num_rejected_common_ancestor_events(msp));
        if (models[j] == MSP_MODEL_HUDSON) {
            CU_ASSERT_EQUAL(msp_get_num_rejected_common_ancestor_events(msp), 0);
        }

        gsl_rng_set(rng, seed);
        ret = msp_reset(msp);
        num_events = 0;
        while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
            msp_verify(msp);
            num_events++;
        }
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_get_num_migration_events(msp, migration_events);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(1 + num_events,
                migration_events[1] + migration_events[2] +
                msp_get_num_recombination_events(msp) +
                msp_get_num_common_ancestor_events(msp) +
                msp_get_num_rejected_common_ancestor_events(msp));
        if (models[j] == MSP_MODEL_HUDSON) {
            CU_ASSERT_EQUAL(msp_get_num_rejected_common_ancestor_events(msp), 0);
        }

        model = msp_get_model(msp)->type;
        CU_ASSERT_EQUAL(model, models[j]);
        model_name = msp_get_model_name(msp);
        CU_ASSERT_STRING_EQUAL(model_name, model_names[j]);

        ret = msp_free(msp);
        CU_ASSERT_EQUAL(ret, 0);
        gsl_rng_free(rng);
        free(msp);
        free(samples);
    }
}

static void
test_simulation_replicates(void)
{
    int ret;
    uint32_t n = 100;
    uint32_t m = 100;
    double mutation_rate = 2;
    size_t num_replicates = 10;
    long seed = 10;
    double migration_matrix[] = {0, 1, 1, 0};
    size_t j;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    msp_t msp;
    tree_sequence_t ts;
    mutgen_t mutgen;
    node_table_t nodes;
    edgeset_table_t edgesets;
    site_table_t sites;
    mutation_table_t mutations;
    migration_table_t migrations;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    /* Set all the table block sizes to 1 to force reallocs */
    ret = node_table_alloc(&nodes, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = edgeset_table_alloc(&edgesets, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = site_table_alloc(&sites, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutation_table_alloc(&mutations, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = migration_table_alloc(&migrations, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    gsl_rng_set(rng, seed);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(&msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_populations(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_migration_matrix(&msp, 4, migration_matrix);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_store_migrations(&msp, true);
    CU_ASSERT_EQUAL(ret, 0);
    /* set all the block sizes to something small to provoke the memory
     * expansions. */
    ret = msp_set_avl_node_block_size(&msp, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_node_mapping_block_size(&msp, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_segment_block_size(&msp, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_coalescence_record_block_size(&msp, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_migration_block_size(&msp, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(&msp, m);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(&msp, 0.5);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = mutgen_alloc(&mutgen, mutation_rate, rng, MSP_ALPHABET_BINARY, 3);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL(ret, 0);

    for (j = 0; j < num_replicates; j++) {
        ret = msp_run(&msp, DBL_MAX, SIZE_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        msp_verify(&msp);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = msp_populate_tables(&msp, 0.25, NULL, &nodes, &edgesets, &migrations);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = mutgen_generate_tables_tmp(&mutgen, &nodes, &edgesets);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = mutgen_populate_tables(&mutgen, &sites, &mutations);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tree_sequence_load_tables_tmp(&ts, &nodes, &edgesets, &migrations,
                &sites, &mutations, 0, NULL);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        verify_simulator_tree_sequence_equality(&msp, &ts, &mutgen, 1.0);
        tree_sequence_print_state(&ts, _devnull);
        ret = msp_reset(&msp);
        CU_ASSERT_EQUAL_FATAL(msp_get_num_coalescence_records(&msp), 0);
        CU_ASSERT_EQUAL_FATAL(msp_get_num_migrations(&msp), 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);

    }
    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(samples);
    node_table_free(&nodes);
    edgeset_table_free(&edgesets);
    site_table_free(&sites);
    mutation_table_free(&mutations);
    migration_table_free(&migrations);
}

static void
test_bottleneck_simulation(void)
{
    int ret;
    uint32_t n = 100;
    uint32_t m = 100;
    long seed = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    double t1 = 0.1;
    double t2 = 0.5;
    int t1_found = 0;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    gsl_rng_set(rng, seed);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    /* set all the block sizes to something small to provoke the memory
     * expansions. */
    ret = msp_set_avl_node_block_size(msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_node_mapping_block_size(msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_segment_block_size(msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_coalescence_record_block_size(msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(msp, m);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(msp, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    /* Add a bottleneck that does nothing at time t1 */
    ret = msp_add_simple_bottleneck(msp, t1, 0, 0);
    CU_ASSERT_EQUAL(ret, 0);
    /* Add a bottleneck that coalesces everything at t2 */
    ret = msp_add_simple_bottleneck(msp, t2, 0, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(msp, t1, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 2);
    CU_ASSERT_FALSE(msp_is_completed(msp));
    msp_print_state(msp, _devnull);

    ret = msp_run(msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(msp_is_completed(msp));
    CU_ASSERT_EQUAL(msp->time, t2);
    msp_verify(msp);

    msp_reset(msp);
    while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
        msp_verify(msp);
        if (msp->time == 0.1) {
            t1_found = 1;
        }
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(t1_found);
    CU_ASSERT_EQUAL(msp->time, t2);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
}

static void
test_large_bottleneck_simulation(void)
{
    int ret;
    uint32_t j;
    uint32_t n = 10000;
    uint32_t m = 100;
    long seed = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    uint32_t num_bottlenecks = 10;
    bottleneck_desc_t bottlenecks[num_bottlenecks];
    double t;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    t = 0.1;
    for (j = 0; j < num_bottlenecks; j++) {
        bottlenecks[j].type = SIMPLE_BOTTLENECK;
        bottlenecks[j].time = t;
        bottlenecks[j].parameter = 0.1;
        t += 0.01;
    }
    /* Set the last bottleneck to be full intensity */
    bottlenecks[num_bottlenecks - 1].parameter = 1.0;

    gsl_rng_set(rng, seed);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(msp, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(msp, m);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_bottlenecks; j++) {
        ret = msp_add_simple_bottleneck(msp, bottlenecks[j].time, 0,
                bottlenecks[j].parameter);
        CU_ASSERT_EQUAL(ret, 0);
    }
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    for (j = 0; j < num_bottlenecks - 1; j++) {
        ret = msp_run(msp, bottlenecks[j].time, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 2);
        CU_ASSERT_FALSE(msp_is_completed(msp));
        CU_ASSERT_EQUAL(msp->time, bottlenecks[j].time);
        msp_verify(msp);
    }
    ret = msp_run(msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(msp_is_completed(msp));
    CU_ASSERT_EQUAL(msp->time, bottlenecks[num_bottlenecks - 1].time);
    msp_verify(msp);

    /* Test out resets on partially completed simulations. */
    ret = msp_reset(msp);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_bottlenecks - 1; j++) {
        ret = msp_run(msp, bottlenecks[j].time, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 2);
    }
    ret = msp_reset(msp);
    msp_verify(msp);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
}

static void
test_multiple_mergers_simulation(void)
{
    int ret;
    size_t j;
    uint32_t n = 100;
    uint32_t m = 100;
    long seed = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    for (j = 0; j < 2; j++) {
        gsl_rng_set(rng, seed);
        /* TODO check non-zero sample times here to make sure they fail. */
        memset(samples, 0, n * sizeof(sample_t));
        ret = msp_alloc(msp, n, samples, rng);
        CU_ASSERT_EQUAL(ret, 0);
        /* TODO what are good parameters here?? */
        if (j == 0) {
            // Use psi = 0.5 for now, but should definitely test for 0 and 1 cases
            ret = msp_set_simulation_model_dirac(msp, 0.5, 1);
        } else {
            ret = msp_set_simulation_model_beta(msp, 1.0, 10.0);
        }
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_num_loci(msp, m);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_scaled_recombination_rate(msp, 10.0);
        CU_ASSERT_EQUAL(ret, 0);
        /* TODO check for adding various complications like multiple populations etc
         * to ensure they fail.
         */
        ret = msp_initialise(msp);
        CU_ASSERT_EQUAL(ret, 0);
        msp_print_state(msp, _devnull);

        ret = msp_run(msp, DBL_MAX, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_TRUE(msp_is_completed(msp));
        CU_ASSERT_TRUE(msp->time > 0);
        msp_verify(msp);

        msp_reset(msp);
        while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
            msp_verify(msp);
        }
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_TRUE(msp_is_completed(msp));

        ret = msp_free(msp);
        CU_ASSERT_EQUAL(ret, 0);
    }
    gsl_rng_free(rng);
    free(msp);
    free(samples);
}

static void
test_node_names(void)
{
    const char *nodes =
        "1  0   0   n1\n"
        "1  0   0   n2\n"
        "0  1   0   A_much_longer_name\n"
        "0  1   0\n"
        "0  1   0   n4";
    const char *edgesets =
        "0  1   2   0,1\n";
    tree_sequence_t ts;
    int ret;
    node_t node;

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 5);

    ret = tree_sequence_get_node(&ts, 0, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(node.name, "n1");

    ret = tree_sequence_get_node(&ts, 1, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(node.name, "n2");

    ret = tree_sequence_get_node(&ts, 2, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(node.name, "A_much_longer_name");

    ret = tree_sequence_get_node(&ts, 3, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(node.name, "");

    ret = tree_sequence_get_node(&ts, 4, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(node.name, "n4");

    tree_sequence_free(&ts);
}

static void
test_simplest_records(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0";
    const char *edgesets =
        "0  1   2   0,1\n";
    tree_sequence_t ts;

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    tree_sequence_free(&ts);
}

static void
test_simplest_nonbinary_records(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0";
    const char *edgesets =
        "0  1   4   0,1,2,3\n";
    tree_sequence_t ts;

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 5);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    tree_sequence_free(&ts);
}

static void
test_simplest_unary_records(void)
{
    int ret;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  1   0\n"
        "0  2   0";
    const char *edgesets =
        "0  1   2   0\n"
        "0  1   3   1\n"
        "0  1   4   2,3\n";
    tree_sequence_t ts, simplified;
    node_id_t sample_ids[] = {0, 1};

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 5);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);

    ret = tree_sequence_initialise(&simplified);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_simplify(&ts, sample_ids, 2, 0, &simplified);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&simplified), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&simplified), 3);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&simplified), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&simplified), 1);

    tree_sequence_free(&ts);
    tree_sequence_free(&simplified);
}

static void
test_simplest_non_sample_leaf_records(void)
{
    int ret;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  0   0\n"
        "0  0   0";
    const char *edgesets =
        "0  1   2   0,1,3,4\n";
    const char *sites =
        "0.1  0\n"
        "0.2  0\n"
        "0.3  0\n"
        "0.4  0\n";
    const char *mutations =
        "0    0     1\n"
        "1    1     1\n"
        "2    3     1\n"
        "3    4     1";
    tree_sequence_t ts, simplified;
    node_id_t sample_ids[] = {0, 1};
    hapgen_t hapgen;
    vargen_t vargen;
    char *haplotype, genotypes[64];
    site_t *site;

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, sites, mutations, NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 5);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);

    ret = hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = hapgen_get_haplotype(&hapgen, 0, &haplotype);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "1000");
    ret = hapgen_get_haplotype(&hapgen, 1, &haplotype);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "0100");
    hapgen_free(&hapgen);

    genotypes[2] = '\0';
    ret = vargen_alloc(&vargen, &ts, MSP_GENOTYPES_AS_CHAR);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    vargen_print_state(&vargen, _devnull);
    ret = vargen_next(&vargen, &site, genotypes);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_STRING_EQUAL(genotypes, "10");
    ret = vargen_next(&vargen, &site, genotypes);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_STRING_EQUAL(genotypes, "01");
    ret = vargen_next(&vargen, &site, genotypes);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_STRING_EQUAL(genotypes, "00");
    ret = vargen_next(&vargen, &site, genotypes);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_STRING_EQUAL(genotypes, "00");
    ret = vargen_next(&vargen, &site, genotypes);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    vargen_free(&vargen);

    ret = tree_sequence_initialise(&simplified);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_simplify(&ts, sample_ids, 2, 0, &simplified);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&simplified), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&simplified), 3);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&simplified), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&simplified), 1);

    tree_sequence_free(&ts);
    tree_sequence_free(&simplified);
}

static void
test_simplest_degenerate_multiple_root_records(void)
{

    int ret;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  1   0\n";
    const char *edgesets =
        "0  1   2   0\n"
        "0  1   3   1\n";
    tree_sequence_t ts, simplified;
    node_id_t sample_ids[] = {0, 1};

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);

    ret = tree_sequence_initialise(&simplified);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_simplify(&ts, sample_ids, 2, 0, &simplified);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_CANNOT_SIMPLIFY);
    tree_sequence_free(&ts);
}

static void
test_simplest_multiple_root_records(void)
{
    int ret;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  1   0\n";
    const char *edgesets =
        "0  1   4   0,1\n"
        "0  1   5   2,3\n";
    tree_sequence_t ts, simplified;
    node_id_t sample_ids[] = {0, 1, 2, 3};

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 6);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);

    ret = tree_sequence_initialise(&simplified);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_simplify(&ts, sample_ids, 4, 0, &simplified);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 6);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);

    tree_sequence_free(&ts);
    tree_sequence_free(&simplified);
}

static void
test_simplest_root_mutations(void)
{
    int ret;
    uint32_t j;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n";
    const char *edgesets =
        "0  1   2   0,1\n";
    const char *sites =
        "0.1 0";
    const char *mutations =
        "0    2     1";
    hapgen_t hapgen;
    char *haplotype;
    int flags = 0;
    node_id_t sample_ids[] = {0, 1};
    tree_sequence_t ts, simplified;

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, sites, mutations, NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tree_sequence_get_num_sites(&ts), 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);

    ret = hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < 2; j++) {
        ret = hapgen_get_haplotype(&hapgen, j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, "1");
    }
    hapgen_free(&hapgen);

    ret = tree_sequence_initialise(&simplified);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_simplify(&ts, sample_ids, 2, flags, &simplified);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&simplified), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&simplified), 3);
    CU_ASSERT_EQUAL(tree_sequence_get_num_sites(&simplified), 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&simplified), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&simplified), 1);
    tree_sequence_free(&simplified);

    flags = MSP_FILTER_INVARIANT_SITES;
    ret = tree_sequence_initialise(&simplified);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_simplify(&ts, sample_ids, 2, flags, &simplified);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&simplified), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&simplified), 3);
    CU_ASSERT_EQUAL(tree_sequence_get_num_sites(&simplified), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&simplified), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&simplified), 1);
    tree_sequence_free(&simplified);

    tree_sequence_free(&ts);
}

static void
test_simplest_back_mutations(void)
{
    int ret;
    uint32_t j;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  2   0\n";
    const char *edgesets =
        "0  1   3   0,1\n"
        "0  1   4   2,3\n";
    const char *sites =
        "0.5 0";
    const char *mutations =
        "0    3     1\n"
        "0    0     0";
    hapgen_t hapgen;
    const char *haplotypes[] = {"0", "1", "0"};
    char *haplotype;
    char genotypes[4];
    tree_sequence_t ts;
    vargen_t vargen;
    site_t *site;

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, sites, mutations, NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 3);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 5);
    CU_ASSERT_EQUAL(tree_sequence_get_num_sites(&ts), 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);

    ret = hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < 2; j++) {
        ret = hapgen_get_haplotype(&hapgen, j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, haplotypes[j]);
    }
    hapgen_free(&hapgen);

    ret = vargen_alloc(&vargen, &ts, MSP_GENOTYPES_AS_CHAR);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    genotypes[3] = '\0';
    ret = vargen_next(&vargen, &site, genotypes);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_STRING_EQUAL(genotypes, "010");
    CU_ASSERT_EQUAL(site->id, 0);
    CU_ASSERT_EQUAL(site->mutations_length, 2);
    vargen_free(&vargen);

    tree_sequence_free(&ts);
}


static void
test_simplest_bad_records(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0";
    const char *edgesets =
        "0  1   2   0,1\n";
    tree_sequence_t ts;
    node_table_t node_table;
    edgeset_table_t edgeset_table;
    int ret;

    ret = node_table_alloc(&node_table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = edgeset_table_alloc(&edgeset_table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    parse_nodes(nodes, &node_table);
    CU_ASSERT_EQUAL_FATAL(node_table.num_rows, 3);
    parse_edgesets(edgesets, &edgeset_table);
    CU_ASSERT_EQUAL_FATAL(edgeset_table.num_rows, 1);

    /* Make sure we have a good set of records */
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table,
            NULL, NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    tree_sequence_free(&ts);

    /* NULL for nodes or edges should be an error */
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, NULL, NULL, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    tree_sequence_free(&ts);

    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, NULL, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    tree_sequence_free(&ts);

    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, NULL, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    tree_sequence_free(&ts);

    /* Bad interval */
    edgeset_table.right[0] = 0.0;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_RECORD_INTERVAL);
    tree_sequence_free(&ts);
    edgeset_table.right[0]= 1.0;

    /* Equal nodes in the children */
    edgeset_table.children[0] = 1;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_UNSORTED_CHILDREN);
    tree_sequence_free(&ts);
    edgeset_table.children[0] = 0;

    /* children node == parent */
    edgeset_table.children[1] = 2;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_NODE_TIME_ORDERING);
    tree_sequence_free(&ts);
    edgeset_table.children[1] = 1;

    /* Unsorted nodes in the children */
    edgeset_table.children[0] = 1;
    edgeset_table.children[1] = 0;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_UNSORTED_CHILDREN);
    tree_sequence_free(&ts);
    edgeset_table.children[0] = 0;
    edgeset_table.children[1] = 1;

    /* Null parent */
    edgeset_table.parent[0] = MSP_NULL_NODE;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_NULL_NODE_IN_RECORD);
    tree_sequence_free(&ts);
    edgeset_table.parent[0] = 2;

    /* parent not in nodes list */
    node_table.num_rows = 2;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_NODE_OUT_OF_BOUNDS);
    tree_sequence_free(&ts);
    node_table.num_rows = 3;

    /* parent negative */
    edgeset_table.parent[0] = -2;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_NODE_OUT_OF_BOUNDS);
    tree_sequence_free(&ts);
    edgeset_table.parent[0] = 2;

    /* Null child */
    edgeset_table.children[0] = MSP_NULL_NODE;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_NULL_NODE_IN_RECORD);
    tree_sequence_free(&ts);
    edgeset_table.children[0] = 0;

    /* child node reference out of bounds */
    edgeset_table.children[0] = 3;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_NODE_OUT_OF_BOUNDS);
    tree_sequence_free(&ts);
    edgeset_table.children[0] = 0;

    /* child node reference negative */
    edgeset_table.children[0] = -2;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_NODE_OUT_OF_BOUNDS);
    tree_sequence_free(&ts);
    edgeset_table.children[0] = 0;

    /* 0 children */
    edgeset_table.total_children_length = 0;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_CHILDREN_ARRAY);
    tree_sequence_free(&ts);
    edgeset_table.total_children_length = 2;

    /* Make sure we've preserved a good tree sequence */
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    tree_sequence_free(&ts);

    node_table_free(&node_table);
    edgeset_table_free(&edgeset_table);
}

static void
test_alphabet_detection(void)
{
    tree_sequence_t ts;

    /* Default to binary */
    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets, NULL,
            NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_num_sites(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_alphabet(&ts), MSP_ALPHABET_BINARY);
    tree_sequence_free(&ts);

    /* All 0->1 mutations are binary */
    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets, NULL,
            "0  0", "0 0 1", NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_num_sites(&ts), 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 1);
    CU_ASSERT_EQUAL(tree_sequence_get_alphabet(&ts), MSP_ALPHABET_BINARY);
    tree_sequence_free(&ts);

    /* A non-zero ancestral state means ASCII */
    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets, NULL,
            "0  1", "0 0 0", NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_num_sites(&ts), 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 1);
    CU_ASSERT_EQUAL(tree_sequence_get_alphabet(&ts), MSP_ALPHABET_ASCII);
    tree_sequence_free(&ts);

    /* Back mutations are still binary */
    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets, NULL,
            "0  0", "0 0 1\n0 1 0", NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_num_sites(&ts), 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_alphabet(&ts), MSP_ALPHABET_BINARY);
    tree_sequence_free(&ts);

    /* Any non-0 or 1 chars make it ASCII */
    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets, NULL,
            "0  0", "0 0 1\n0 1 A", NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_num_sites(&ts), 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_alphabet(&ts), MSP_ALPHABET_ASCII);
    tree_sequence_free(&ts);

}

static void
test_single_tree_good_records(void)
{
    tree_sequence_t ts;

    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets,
            NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 7);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    tree_sequence_free(&ts);
}


static void
test_single_nonbinary_tree_good_records(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  2   0\n"
        "0  3   0\n";
    const char *edgesets =
        "0 1 7 0,1,2,3\n"
        "0 1 8 4,5\n"
        "0 1 9 6,7,8";
    tree_sequence_t ts;

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 7);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 10);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    tree_sequence_free(&ts);
}

static void
test_single_tree_bad_records(void)
{
    int ret = 0;
    tree_sequence_t ts;
    node_table_t node_table;
    edgeset_table_t edgeset_table;

    ret = node_table_alloc(&node_table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = edgeset_table_alloc(&edgeset_table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    parse_nodes(single_tree_ex_nodes, &node_table);
    CU_ASSERT_EQUAL_FATAL(node_table.num_rows, 7);
    parse_edgesets(single_tree_ex_edgesets, &edgeset_table);
    CU_ASSERT_EQUAL_FATAL(edgeset_table.num_rows, 3);

    /* Not sorted in time order */
    node_table.time[5] = 0.5;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_RECORDS_NOT_TIME_SORTED);
    tree_sequence_free(&ts);
    node_table.time[5] = 2.0;

    /* Left value greater than sequence right */
    edgeset_table.left[2] = 2.0;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_RECORD_INTERVAL);
    tree_sequence_free(&ts);
    edgeset_table.left[2] = 0.0;

    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    tree_sequence_free(&ts);

    edgeset_table_free(&edgeset_table);
    node_table_free(&node_table);
}


static void
test_single_tree_good_mutations(void)
{
    tree_sequence_t ts;
    size_t j;
    size_t num_sites = 3;
    size_t num_mutations = 7;
    site_t other_sites[num_sites];
    mutation_t other_mutations[num_mutations];
    int ret;

    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets,
            NULL, single_tree_ex_sites, single_tree_ex_mutations, NULL);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 7);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_sites(&ts), num_sites);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), num_mutations);

    for (j = 0; j < num_sites; j++) {
        ret = tree_sequence_get_site(&ts, j, other_sites + j);
        CU_ASSERT_EQUAL(ret, 0);
    }
    for (j = 0; j < num_mutations; j++) {
        ret = tree_sequence_get_mutation(&ts, j, other_mutations + j);
        CU_ASSERT_EQUAL(ret, 0);
    }
    CU_ASSERT_EQUAL(other_sites[0].position, 0.1);
    CU_ASSERT_STRING_EQUAL(other_sites[0].ancestral_state, "0");
    CU_ASSERT_EQUAL(other_sites[1].position, 0.2);
    CU_ASSERT_STRING_EQUAL(other_sites[1].ancestral_state, "0");
    CU_ASSERT_EQUAL(other_sites[2].position, 0.3);
    CU_ASSERT_STRING_EQUAL(other_sites[2].ancestral_state, "0");

    CU_ASSERT_EQUAL(other_mutations[0].index, 0);
    CU_ASSERT_EQUAL(other_mutations[0].node, 2);
    CU_ASSERT_STRING_EQUAL(other_mutations[0].derived_state, "1");
    CU_ASSERT_EQUAL(other_mutations[1].index, 1);
    CU_ASSERT_EQUAL(other_mutations[1].node, 4);
    CU_ASSERT_STRING_EQUAL(other_mutations[1].derived_state, "1");
    CU_ASSERT_EQUAL(other_mutations[2].index, 2);
    CU_ASSERT_EQUAL(other_mutations[2].node, 0);
    CU_ASSERT_STRING_EQUAL(other_mutations[2].derived_state, "0");
    CU_ASSERT_EQUAL(other_mutations[3].index, 3);
    CU_ASSERT_EQUAL(other_mutations[3].node, 0);
    CU_ASSERT_STRING_EQUAL(other_mutations[3].derived_state, "1");
    CU_ASSERT_EQUAL(other_mutations[4].index, 4);
    CU_ASSERT_EQUAL(other_mutations[4].node, 1);
    CU_ASSERT_STRING_EQUAL(other_mutations[4].derived_state, "1");
    CU_ASSERT_EQUAL(other_mutations[5].index, 5);
    CU_ASSERT_EQUAL(other_mutations[5].node, 2);
    CU_ASSERT_STRING_EQUAL(other_mutations[5].derived_state, "1");
    CU_ASSERT_EQUAL(other_mutations[6].index, 6);
    CU_ASSERT_EQUAL(other_mutations[6].node, 3);
    CU_ASSERT_STRING_EQUAL(other_mutations[6].derived_state, "1");

    tree_sequence_free(&ts);
}

static void
test_single_tree_bad_mutations(void)
{
    int ret = 0;
    const char *sites =
        "0       0\n"
        "0.1     0\n"
        "0.2     0\n";
    const char *mutations =
        "0   0  1\n"
        "1   1  1\n"
        "2   0  1\n"
        "2   1  1\n";
    tree_sequence_t ts;
    node_table_t node_table;
    edgeset_table_t edgeset_table;
    site_table_t site_table;
    mutation_table_t mutation_table;

    ret = node_table_alloc(&node_table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = edgeset_table_alloc(&edgeset_table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = site_table_alloc(&site_table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutation_table_alloc(&mutation_table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    parse_nodes(single_tree_ex_nodes, &node_table);
    CU_ASSERT_EQUAL_FATAL(node_table.num_rows, 7);
    parse_edgesets(single_tree_ex_edgesets, &edgeset_table);
    CU_ASSERT_EQUAL_FATAL(edgeset_table.num_rows, 3);
    parse_sites(sites, &site_table);
    parse_mutations(mutations, &mutation_table);
    CU_ASSERT_EQUAL_FATAL(site_table.num_rows, 3);
    CU_ASSERT_EQUAL_FATAL(mutation_table.num_rows, 4);

    /* negative coordinate */
    site_table.position[0] = -1.0;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            &site_table, &mutation_table, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_SITE_POSITION);
    tree_sequence_free(&ts);
    site_table.position[0] = 0.0;

    /* coordinate == sequence length */
    site_table.position[2] = 1.0;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            &site_table, &mutation_table, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_SITE_POSITION);
    tree_sequence_free(&ts);
    site_table.position[2] = 0.2;

    /* coordinate > sequence length */
    site_table.position[2] = 1.1;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            &site_table, &mutation_table, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_SITE_POSITION);
    tree_sequence_free(&ts);
    site_table.position[2] = 0.2;

    /* Unsorted positions */
    site_table.position[0] = 0.3;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            &site_table, &mutation_table, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_UNSORTED_SITES);
    tree_sequence_free(&ts);
    site_table.position[0] = 0.0;

    /* site < 0 */
    mutation_table.site[0] = -2;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            &site_table, &mutation_table, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_SITE_OUT_OF_BOUNDS);
    tree_sequence_free(&ts);
    mutation_table.site[0] = 0;

    /* site == num_sites */
    mutation_table.site[0] = 3;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            &site_table, &mutation_table, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_SITE_OUT_OF_BOUNDS);
    tree_sequence_free(&ts);
    mutation_table.site[0] = 0;

    /* node = NULL */
    mutation_table.node[0] = MSP_NULL_NODE;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            &site_table, &mutation_table, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_NODE_OUT_OF_BOUNDS);
    tree_sequence_free(&ts);
    mutation_table.node[0] = 0;

    /* node >= num_nodes */
    mutation_table.node[0] = 7;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            &site_table, &mutation_table, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_NODE_OUT_OF_BOUNDS);
    tree_sequence_free(&ts);
    mutation_table.node[0] = 0;

    /* Two mutations at the same node for a given site */
    /* FIXME: this condition was relaxed because of the difficulty in maintaining
     * it for the outputs of simplify(). This should be fixed and reinstated though.
     */
    /* mutation_table.node[2] = 1; */
    /* ret = tree_sequence_initialise(&ts); */
    /* CU_ASSERT_EQUAL_FATAL(ret, 0); */
    /* ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL, */
    /*         &site_table, &mutation_table, 0, NULL); */
    /* CU_ASSERT_EQUAL(ret, MSP_ERR_DUPLICATE_MUTATION_NODES); */
    /* tree_sequence_free(&ts); */
    /* mutation_table.node[2] = 0; */

    /* Unsorted nodes */
    /* FIXME: this condition was relaxed because it's too difficult to maintain
     * on the output of simplify.
     */
    /* mutation_table.node[3] = 5; */
    /* ret = tree_sequence_initialise(&ts); */
    /* CU_ASSERT_EQUAL_FATAL(ret, 0); */
    /* ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL, */
    /*         &site_table, &mutation_table, 0, NULL); */
    /* CU_ASSERT_EQUAL(ret, MSP_ERR_UNSORTED_MUTATION_NODES); */
    /* tree_sequence_free(&ts); */
    /* mutation_table.node[3] = 1; */

    /* Check to make sure we've maintained legal mutations */
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            &site_table, &mutation_table, 0, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 4);
    tree_sequence_free(&ts);

    node_table_free(&node_table);
    edgeset_table_free(&edgeset_table);
    site_table_free(&site_table);
    mutation_table_free(&mutation_table);
}

static void
test_single_tree_iter(void)
{
    int ret;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  2   0\n"
        "0  3   0\n";
    const char *edgesets =
        "0  6   4   0,1\n"
        "0  6   5   2,3\n"
        "0  6   6   4,5\n";
    node_id_t parents[] = {4, 4, 5, 5, 6, 6, MSP_NULL_NODE};
    tree_sequence_t ts;
    sparse_tree_t tree;
    node_id_t u, v, w;
    size_t num_leaves;
    uint32_t num_nodes = 7;

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, NULL, NULL, NULL);
    ret = sparse_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    sparse_tree_print_state(&tree, _devnull);

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

    ret = sparse_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    sparse_tree_free(&tree);
    tree_sequence_free(&ts);
}

static void
test_single_nonbinary_tree_iter(void)
{

    int ret;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  2   0\n"
        "0  3   0\n";
    const char *edgesets =
        "0  1   7   0,1,2,3\n"
        "0  1   8   4,5\n"
        "0  1   9   6,7,8\n";
    node_id_t parents[] = {7, 7, 7, 7, 8, 8, 9, 9, 9, MSP_NULL_NODE};
    tree_sequence_t ts;
    sparse_tree_t tree;
    node_id_t u, v, w, *children;
    size_t num_leaves, num_children;
    size_t num_nodes = 10;
    size_t num_samples = 7;

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, NULL, NULL, NULL);
    ret = sparse_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    sparse_tree_print_state(&tree, _devnull);

    for (u = 0; u < num_nodes; u++) {
        ret = sparse_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
    }
    for (u = 0; u < num_samples; u++) {
        ret = sparse_tree_get_num_leaves(&tree, u, &num_leaves);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(num_leaves, 1);
        ret = sparse_tree_get_children(&tree, u, &num_children, &children);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(children, NULL);
        CU_ASSERT_EQUAL(num_children, 0);
    }

    u = 7;
    ret = sparse_tree_get_num_leaves(&tree, u, &num_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_leaves, 4);
    ret = sparse_tree_get_children(&tree, u, &num_children, &children);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_children, 4);
    CU_ASSERT_EQUAL(children[0], 0);
    CU_ASSERT_EQUAL(children[1], 1);
    CU_ASSERT_EQUAL(children[2], 2);
    CU_ASSERT_EQUAL(children[3], 3);

    u = 8;
    ret = sparse_tree_get_num_leaves(&tree, u, &num_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_leaves, 2);
    ret = sparse_tree_get_children(&tree, u, &num_children, &children);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_children, 2);
    CU_ASSERT_EQUAL(children[0], 4);
    CU_ASSERT_EQUAL(children[1], 5);

    u = 9;
    ret = sparse_tree_get_num_leaves(&tree, u, &num_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_leaves, 7);
    ret = sparse_tree_get_children(&tree, u, &num_children, &children);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_children, 3);
    CU_ASSERT_EQUAL(children[0], 6);
    CU_ASSERT_EQUAL(children[1], 7);
    CU_ASSERT_EQUAL(children[2], 8);

    ret = sparse_tree_get_root(&tree, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 9);

    ret = sparse_tree_get_mrca(&tree, 0, 1, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 7);
    ret = sparse_tree_get_mrca(&tree, 0, 4, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 9);

    ret = sparse_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    sparse_tree_free(&tree);
    tree_sequence_free(&ts);
}

static void
test_single_tree_iter_times(void)
{
    int ret = 0;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  2   0\n"
        "1  3   0\n"
        "0  1   0\n"
        "0  4   0\n"
        "0  5   0\n";
    const char *edgesets =
        "0  6   4   0,1\n"
        "0  6   5   2,3\n"
        "0  6   6   4,5\n";
    node_id_t parents[] = {4, 4, 5, 5, 6, 6, MSP_NULL_NODE};
    double times[] = {0.0, 0.0, 2.0, 3.0, 1.0, 4.0, 5.0};
    double t;
    tree_sequence_t ts;
    sparse_tree_t tree;
    node_id_t u, v;
    uint32_t num_nodes = 7;

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, NULL, NULL, NULL);
    ret = sparse_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    sparse_tree_print_state(&tree, _devnull);

    for (u = 0; u < num_nodes; u++) {
        ret = sparse_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
        ret = sparse_tree_get_time(&tree, u, &t);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(t, times[u]);
    }
    ret = sparse_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    sparse_tree_free(&tree);
    tree_sequence_free(&ts);
}

static void
test_single_tree_hapgen_char_alphabet(void)
{
    int ret = 0;
    const char *sites =
        "0.0    A\n"
        "0.1    A\n"
        "0.2    C\n"
        "0.4    A\n";
    const char *mutations =
        "0    0     T\n"
        "1    1     T\n"
        "2    0     G\n"
        "2    1     A\n"
        "2    2     T\n"  // A bunch of different leaf mutations
        "3    4     T\n"
        "3    0     A\n"; // A back mutation from T -> A
    uint32_t num_samples = 4;
    tree_sequence_t ts;
    char *haplotype;
    size_t j;
    hapgen_t hapgen;

    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets, NULL,
            NULL, NULL, NULL);
    ret = hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < num_samples; j++) {
        ret = hapgen_get_haplotype(&hapgen, j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, "");
    }
    ret = hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tree_sequence_free(&ts);

    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets, NULL,
            sites, mutations, NULL);
    ret = hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    hapgen_print_state(&hapgen, _devnull);

    ret = hapgen_get_haplotype(&hapgen, 0, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "TAGA");
    ret = hapgen_get_haplotype(&hapgen, 1, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "ATAT");
    ret = hapgen_get_haplotype(&hapgen, 2, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "AATA");
    ret = hapgen_get_haplotype(&hapgen, 3, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "AACA");

    ret = hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tree_sequence_free(&ts);
}

static void
test_single_tree_vargen_char_alphabet(void)
{
    int ret = 0;
    const char *sites =
        "0.0    A\n"
        "0.1    A\n"
        "0.2    C\n"
        "0.4    A\n";
    const char *mutations =
        "0    0     T\n"
        "1    1     T\n"
        "2    0     G\n"
        "2    1     A\n"
        "2    2     T\n"  // A bunch of different leaf mutations
        "3    4     T\n"
        "3    0     A\n"; // A back mutation from T -> A
    tree_sequence_t ts;
    vargen_t vargen;

    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets, NULL,
            sites, mutations, NULL);
    ret = vargen_alloc(&vargen, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_NONBINARY_MUTATIONS_UNSUPPORTED);
    ret = vargen_alloc(&vargen, &ts, MSP_GENOTYPES_AS_CHAR);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_NONBINARY_MUTATIONS_UNSUPPORTED);

    tree_sequence_free(&ts);
}

static void
test_single_tree_hapgen_binary_alphabet(void)
{
    int ret = 0;
    uint32_t num_samples = 4;
    tree_sequence_t ts;
    char *haplotype;
    size_t j;
    hapgen_t hapgen;

    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets, NULL,
            NULL, NULL, NULL);
    ret = hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < num_samples; j++) {
        ret = hapgen_get_haplotype(&hapgen, j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, "");
    }
    ret = hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tree_sequence_free(&ts);

    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets, NULL,
            single_tree_ex_sites, single_tree_ex_mutations, NULL);
    ret = hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    hapgen_print_state(&hapgen, _devnull);

    ret = hapgen_get_haplotype(&hapgen, 0, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "001");
    ret = hapgen_get_haplotype(&hapgen, 1, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "011");
    ret = hapgen_get_haplotype(&hapgen, 2, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "101");
    ret = hapgen_get_haplotype(&hapgen, 3, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "001");

    ret = hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tree_sequence_free(&ts);
}


static void
test_single_tree_vargen_binary_alphabet(void)
{
    int ret = 0;
    tree_sequence_t ts;
    char genotypes[5];
    site_t *site;
    vargen_t vargen;

    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets, NULL,
            single_tree_ex_sites, single_tree_ex_mutations, NULL);
    ret = vargen_alloc(&vargen, &ts, MSP_GENOTYPES_AS_CHAR);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    vargen_print_state(&vargen, _devnull);

    genotypes[4] = '\0';
    ret = vargen_next(&vargen, &site, genotypes);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_STRING_EQUAL(genotypes, "0010");
    CU_ASSERT_EQUAL(site->id, 0);
    CU_ASSERT_EQUAL(site->mutations_length, 1);

    ret = vargen_next(&vargen, &site, genotypes);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_STRING_EQUAL(genotypes, "0100");
    CU_ASSERT_EQUAL(site->id, 1);
    CU_ASSERT_EQUAL(site->mutations_length, 2);

    ret = vargen_next(&vargen, &site, genotypes);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_STRING_EQUAL(genotypes, "1111");
    CU_ASSERT_EQUAL(site->id, 2);
    CU_ASSERT_EQUAL(site->mutations_length, 4);

    ret = vargen_next(&vargen, &site, genotypes);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = vargen_free(&vargen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tree_sequence_free(&ts);
}

static void
test_single_tree_simplify(void)
{
    tree_sequence_t ts;

    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets, NULL,
            single_tree_ex_sites, single_tree_ex_mutations, NULL);
    verify_simplify(&ts);

    tree_sequence_free(&ts);
}

static void
test_single_tree_inconsistent_mutations(void)
{
    const char *sites =
        "0.0     0\n"
        "0.1     0\n"
        "0.2     0\n";
    const char *mutations =
        "0    0     1\n"
        "1    1     1\n"
        "2    4     1\n"
        "2    0     1\n";
    tree_sequence_t ts;
    char genotypes[5];
    site_t *site;
    vargen_t vargen;
    hapgen_t hapgen;
    int ret;

    tree_sequence_from_text(&ts, single_tree_ex_nodes, single_tree_ex_edgesets, NULL,
            sites, mutations, NULL);

    ret = hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INCONSISTENT_MUTATIONS);
    ret = hapgen_free(&hapgen);

    ret = vargen_alloc(&vargen, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = vargen_next(&vargen, &site, genotypes);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    ret = vargen_next(&vargen, &site, genotypes);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    ret = vargen_next(&vargen, &site, genotypes);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_INCONSISTENT_MUTATIONS);
    ret = vargen_free(&vargen);

    tree_sequence_free(&ts);
}

static void
test_single_unary_tree_hapgen(void)
{
    int ret = 0;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  1   0\n"
        "0  2   0\n"
        "0  3   0\n"
        "0  4   0\n";
    const char *edgesets =
        "0 1 2 0\n"
        "0 1 3 1\n"
        "0 1 4 2,3\n"
        "0 1 5 4\n"
        "0 1 6 5\n";
    const char *sites =
        "0     0\n"
        "0.1   0\n"
        "0.2   0\n"
        "0.3   0\n";
    const char *mutations =
        "0    0     1\n"
        "1    1     1\n"
        "2    4     1\n"
        "3    5     1\n";
    tree_sequence_t ts;
    size_t num_samples = 2;
    size_t j;
    hapgen_t hapgen;
    char *haplotype;

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, NULL, NULL, NULL);

    ret = hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < num_samples; j++) {
        ret = hapgen_get_haplotype(&hapgen, j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, "");
    }
    ret = hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tree_sequence_free(&ts);

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, sites, mutations, NULL);
    ret = hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    hapgen_print_state(&hapgen, _devnull);

    ret = hapgen_get_haplotype(&hapgen, 0, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "1011");
    ret = hapgen_get_haplotype(&hapgen, 1, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "0111");

    ret = hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tree_sequence_free(&ts);
}

static void
test_single_tree_mutgen(void)
{
    int ret = 0;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  2   0\n"
        "0  3   0\n";
    const char *edgesets =
        "0  1   4   0,1\n"
        "0  1   5   2,3\n"
        "0  1   6   4,5\n";
    size_t j;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    node_table_t node_table;
    edgeset_table_t edgeset_table;
    site_table_t sites, sites_after;
    mutation_table_t mutations, mutations_after;

    CU_ASSERT_FATAL(rng != NULL);
    ret = node_table_alloc(&node_table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = edgeset_table_alloc(&edgeset_table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    parse_nodes(nodes, &node_table);
    parse_edgesets(edgesets, &edgeset_table);
    ret = site_table_alloc(&sites, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutation_table_alloc(&mutations, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = site_table_alloc(&sites_after, 100, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutation_table_alloc(&mutations_after, 100, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = mutgen_alloc(&mutgen, 0.0, rng, MSP_ALPHABET_BINARY, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate_tables_tmp(&mutgen, &node_table, &edgeset_table);
    CU_ASSERT_EQUAL(ret, 0);
    ret = mutgen_populate_tables(&mutgen, &sites, &mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(mutations.num_rows == 0);
    mutgen_print_state(&mutgen, _devnull);
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    gsl_rng_set(rng, 1);
    ret = mutgen_alloc(&mutgen, 10.0, rng, MSP_ALPHABET_BINARY, 100);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate_tables_tmp(&mutgen, &node_table, &edgeset_table);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    mutgen_print_state(&mutgen, _devnull);
    ret = mutgen_populate_tables(&mutgen, &sites, &mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(mutations.num_rows > 0);
    CU_ASSERT_TRUE(mutations.num_rows == sites.num_rows);
    for (j = 0; j < mutations.num_rows; j++) {
        CU_ASSERT_TRUE(mutations.site[j] == j);
        CU_ASSERT_TRUE(sites.position[j] <= 1.0);
        CU_ASSERT_TRUE(mutations.node[j] < 6);
        CU_ASSERT_EQUAL(sites.ancestral_state[j], '0');
        CU_ASSERT_EQUAL(mutations.derived_state[j], '1');
    }
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Test the reallocing behavior by setting a very small
     * block size.
     */
    gsl_rng_set(rng, 1);
    ret = mutgen_alloc(&mutgen, 10.0, rng, MSP_ALPHABET_BINARY, 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = mutgen_alloc(&mutgen, 10.0, rng, MSP_ALPHABET_BINARY, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_generate_tables_tmp(&mutgen, &node_table, &edgeset_table);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutgen_populate_tables(&mutgen, &sites_after, &mutations_after);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(mutation_table_equal(&mutations, &mutations_after));
    CU_ASSERT_TRUE(site_table_equal(&sites, &sites_after));
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    edgeset_table_free(&edgeset_table);
    node_table_free(&node_table);
    mutation_table_free(&mutations);
    site_table_free(&sites);
    mutation_table_free(&mutations_after);
    site_table_free(&sites_after);
    gsl_rng_free(rng);
}

static void
verify_trees(tree_sequence_t *ts, uint32_t num_trees, node_id_t* parents)
{
    int ret;
    node_id_t u, v;
    uint32_t j, mutation_index, site_index;
    list_len_t k, l, tree_sites_length;
    site_t *sites = NULL;
    sparse_tree_t tree;
    size_t num_nodes = tree_sequence_get_num_nodes(ts);
    size_t num_sites = tree_sequence_get_num_sites(ts);
    size_t num_mutations = tree_sequence_get_num_mutations(ts);

    ret = sparse_tree_alloc(&tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(ts), num_trees);

    site_index = 0;
    mutation_index = 0;
    j = 0;
    for (ret = sparse_tree_first(&tree); ret == 1; ret = sparse_tree_next(&tree)) {
        CU_ASSERT_EQUAL(j, tree.index);
        sparse_tree_print_state(&tree, _devnull);
        for (u = 0; u < num_nodes; u++) {
            ret = sparse_tree_get_parent(&tree, u, &v);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(v, parents[j * num_nodes + u]);
        }
        ret = sparse_tree_get_sites(&tree, &sites, &tree_sites_length);
        CU_ASSERT_EQUAL(ret, 0);
        for (k = 0; k < tree_sites_length; k++) {
            CU_ASSERT_EQUAL(sites[k].id, site_index);
            for (l = 0; l < sites[k].mutations_length; l++) {
                CU_ASSERT_EQUAL(sites[k].mutations[l].index, mutation_index);
                CU_ASSERT_EQUAL(sites[k].mutations[l].site, site_index);
                mutation_index++;
            }
            site_index++;
        }
        j++;
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(site_index, num_sites);
    CU_ASSERT_EQUAL(mutation_index, num_mutations);

    ret = sparse_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    sparse_tree_free(&tree);
}

static void
verify_trees_consistent(tree_sequence_t *ts)
{
    int ret, found;
    size_t sample_size, num_trees;
    node_id_t u, v, root, *children;
    size_t j, k, num_children;
    sparse_tree_t tree;

    sample_size = tree_sequence_get_sample_size(ts);
    ret = sparse_tree_alloc(&tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    num_trees = 0;
    for (ret = sparse_tree_first(&tree); ret == 1;
            ret = sparse_tree_next(&tree)) {
        ret = sparse_tree_get_root(&tree, &root);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(tree.index, num_trees);
        num_trees++;
        for (j = 0; j < sample_size; j++) {
            v = j;
            while (v != MSP_NULL_NODE) {
                u = v;
                ret = sparse_tree_get_parent(&tree, u, &v);
                CU_ASSERT_EQUAL(ret, 0);
                if (v != MSP_NULL_NODE) {
                    ret = sparse_tree_get_children(&tree, v, &num_children, &children);
                    CU_ASSERT_EQUAL(ret, 0);
                    CU_ASSERT(num_children >= 1);
                    found = 0;
                    for (k = 0; k < num_children; k++) {
                        if (children[k] == u) {
                            found = 1;
                        }
                    }
                    CU_ASSERT(found);
                }
            }
            CU_ASSERT_EQUAL(u, root);
        }
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(ts), num_trees);

    sparse_tree_free(&tree);
}

static void
test_sparse_tree_errors(void)
{
    int ret;
    size_t j;
    uint32_t num_nodes = 9;
    uint32_t u;
    tree_sequence_t ts, other_ts;
    sparse_tree_t t, other_t;
    node_id_t bad_nodes[] = {num_nodes, num_nodes + 1, -1};
    node_id_t tracked_leaves[] = {0, 0, 0};

    tree_sequence_from_text(&ts, paper_ex_nodes, paper_ex_edgesets, NULL, NULL, NULL, NULL);

    ret = sparse_tree_alloc(&t, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = sparse_tree_alloc(&t, &ts, MSP_LEAF_COUNTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = sparse_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 1);

    /* Out-of-bounds queries */
    for (j = 0; j < sizeof(bad_nodes) / sizeof(node_id_t); j++) {
        u = bad_nodes[j];
        ret = sparse_tree_get_children(&t, u, NULL, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_parent(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_time(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_mrca(&t, u, 0, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_mrca(&t, 0, u, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_num_leaves(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_num_tracked_leaves(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_leaf_list(&t, u, NULL, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
    }

    tracked_leaves[0] = 0;
    tracked_leaves[1] = tree_sequence_get_sample_size(&ts);
    ret = sparse_tree_set_tracked_leaves(&t, 2, tracked_leaves);
    CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
    tracked_leaves[1] = 0;
    ret = sparse_tree_set_tracked_leaves(&t, 2, tracked_leaves);
    CU_ASSERT_EQUAL(ret, MSP_ERR_DUPLICATE_SAMPLE);
    ret = sparse_tree_set_tracked_leaves_from_leaf_list(&t, NULL, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    tree_sequence_from_text(&other_ts, paper_ex_nodes, paper_ex_edgesets, NULL, NULL, NULL, NULL);
    ret = sparse_tree_alloc(&other_t, &other_ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_copy(&t, &t);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = sparse_tree_copy(&t, &other_t);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    /* TODO run checks for the various unsupported operations with flags */

    sparse_tree_free(&t);
    sparse_tree_free(&other_t);
    tree_sequence_free(&other_ts);
    tree_sequence_free(&ts);
}

static void
test_tree_sequence_iter(void)
{
    node_id_t parents[] = {
        6, 5, 8, 5, MSP_NULL_NODE, 6, 8, MSP_NULL_NODE, MSP_NULL_NODE,
        6, 5, 4, 4, 5, 6, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
        7, 5, 4, 4, 5, 7, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
    };
    uint32_t num_trees = 3;
    tree_sequence_t ts;

    tree_sequence_from_text(&ts, paper_ex_nodes, paper_ex_edgesets, NULL,
            paper_ex_sites, paper_ex_mutations, NULL);
    verify_trees(&ts, num_trees, parents);
    tree_sequence_free(&ts);
}

static void
test_unary_tree_sequence_iter(void)
{
    node_id_t parents[] = {
        6, 5, 7, 5, MSP_NULL_NODE, 6, 8, 8, MSP_NULL_NODE,
        6, 5, 4, 4, 5, 6, 8, MSP_NULL_NODE, MSP_NULL_NODE,
        7, 5, 4, 4, 5, 7, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
    };
    tree_sequence_t ts;
    uint32_t num_trees = 3;

    tree_sequence_from_text(&ts, unary_ex_nodes, unary_ex_edgesets, NULL,
            unary_ex_sites, unary_ex_mutations, NULL);
    verify_trees(&ts, num_trees, parents);
    tree_sequence_free(&ts);
}

static void
test_nonbinary_tree_sequence_iter(void)
{
    /* We make one mutation for each tree */
    node_id_t parents[] = {
        8, 8, 8, 8, 10, 10, 9, 10, 9, 12, 12, MSP_NULL_NODE, MSP_NULL_NODE,
        8, 8, 8, 8, 10, 11, 9, 10, 9, 11, 12, 12, MSP_NULL_NODE,
    };

    tree_sequence_t ts;
    uint32_t num_trees = 2;

    tree_sequence_from_text(&ts, nonbinary_ex_nodes, nonbinary_ex_edgesets, NULL,
            nonbinary_ex_sites, nonbinary_ex_mutations, NULL);
    verify_trees(&ts, num_trees, parents);
    tree_sequence_free(&ts);
}

static void
test_left_to_right_tree_sequence_iter(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  0.090   0\n"
        "0  0.170   0\n"
        "0  0.253   0\n"
        "0  0.071   0\n"
        "0  0.202   0\n";
    const char *edgesets =
        "2 10 7 2,3\n"
        "0 2  4 1,3\n"
        "2 10 4 1,7\n"
        "0 7  5 0,4\n"
        "7 10 8 0,4\n"
        "0 2  6 2,5\n";
    const char *sites =
        "1      0\n"
        "4.5    0\n"
        "8.5    0\n";
    const char *mutations =
        "0    2    1\n"
        "1    0    1\n"
        "2    4    1\n";

    node_id_t parents[] = {
        5, 4, 6, 4, 5, 6, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
        5, 4, 7, 7, 5, MSP_NULL_NODE, MSP_NULL_NODE, 4, MSP_NULL_NODE,
        8, 4, 7, 7, 8, MSP_NULL_NODE, MSP_NULL_NODE, 4, MSP_NULL_NODE,
    };
    tree_sequence_t ts;
    uint32_t num_trees = 3;

    tree_sequence_from_text(&ts, nodes, edgesets, NULL, sites, mutations, NULL);
    verify_trees(&ts, num_trees, parents);
    tree_sequence_free(&ts);
}

typedef struct {
    uint32_t tree_index;
    uint32_t node;
    uint32_t count;
} leaf_count_test_t;

static void
verify_leaf_counts(tree_sequence_t *ts, size_t num_tests,
        leaf_count_test_t *tests)
{
    int ret;
    size_t j, num_leaves, n, k;
    node_id_t *tracked_leaves = NULL;
    sparse_tree_t tree;
    leaf_list_node_t *u, *head, *tail;

    n = tree_sequence_get_sample_size(ts);

    /* First run without the MSP_LEAF_COUNTS feature */
    ret = sparse_tree_alloc(&tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = sparse_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        ret = sparse_tree_get_num_leaves(&tree, tests[j].node, &num_leaves);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_leaves);
        /* all operations depending on tracked leaves should fail. */
        ret = sparse_tree_get_num_tracked_leaves(&tree, 0, &num_leaves);
        CU_ASSERT_EQUAL(ret, MSP_ERR_UNSUPPORTED_OPERATION);
        ret = sparse_tree_get_leaf_list(&tree, 0, NULL, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_UNSUPPORTED_OPERATION);
    }
    sparse_tree_free(&tree);

    /* Now run with MSP_LEAF_COUNTS but with no leaves tracked. */
    ret = sparse_tree_alloc(&tree, ts, MSP_LEAF_COUNTS);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = sparse_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        ret = sparse_tree_get_num_leaves(&tree, tests[j].node, &num_leaves);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_leaves);
        /* all operations depending on tracked leaves should fail. */
        ret = sparse_tree_get_num_tracked_leaves(&tree, 0, &num_leaves);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(num_leaves, 0);
        /* Getting leaf lists should still fail, as it's not enabled. */
        ret = sparse_tree_get_leaf_list(&tree, 0, NULL, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_UNSUPPORTED_OPERATION);
    }
    sparse_tree_free(&tree);

    /* Run with MSP_LEAF_LISTS, but without MSP_LEAF_COUNTS */
    ret = sparse_tree_alloc(&tree, ts, MSP_LEAF_LISTS);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = sparse_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        ret = sparse_tree_get_num_leaves(&tree, tests[j].node, &num_leaves);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_leaves);
        /* all operations depending on tracked leaves should fail. */
        ret = sparse_tree_get_num_tracked_leaves(&tree, 0, &num_leaves);
        CU_ASSERT_EQUAL(ret, MSP_ERR_UNSUPPORTED_OPERATION);
        /* Getting leaf lists should still fail, as it's not enabled. */
        ret = sparse_tree_get_leaf_list(&tree, tests[j].node, &head, &tail);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        u = head;
        k = 0;
        while (1) {
            k++;
            if (u == tail) {
                break;
            }
            u = u->next;
        }
        CU_ASSERT_EQUAL(tests[j].count, k);
    }
    sparse_tree_free(&tree);

    /* Now use MSP_LEAF_COUNTS|MSP_LEAF_LISTS */
    tracked_leaves = malloc(n * sizeof(node_id_t));
    for (j = 0; j < n; j++) {
        tracked_leaves[j] = j;
    }
    ret = sparse_tree_alloc(&tree, ts, MSP_LEAF_COUNTS|MSP_LEAF_LISTS);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_set_tracked_leaves(&tree, n, tracked_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = sparse_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        ret = sparse_tree_get_num_leaves(&tree, tests[j].node, &num_leaves);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_leaves);
        /* We're tracking all leaves, so the count should be the same */
        ret = sparse_tree_get_num_tracked_leaves(&tree, tests[j].node,
                &num_leaves);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_leaves);
        ret = sparse_tree_get_leaf_list(&tree, tests[j].node, &head, &tail);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        u = head;
        k = 0;
        while (1) {
            k++;
            if (u == tail) {
                break;
            }
            u = u->next;
        }
        CU_ASSERT_EQUAL(tests[j].count, k);
    }
    sparse_tree_free(&tree);
    free(tracked_leaves);
}


static void
verify_leaf_sets_for_tree(sparse_tree_t *tree)
{
    int ret, stack_top, j;
    node_id_t u, v, n, num_nodes, num_leaves;
    size_t tmp;
    node_id_t *stack, *leaves;
    leaf_list_node_t *z, *head, *tail;
    tree_sequence_t *ts = tree->tree_sequence;

    n = tree_sequence_get_sample_size(ts);
    num_nodes = tree_sequence_get_num_nodes(ts);
    stack = malloc(n * sizeof(node_id_t));
    leaves = malloc(n * sizeof(node_id_t));
    CU_ASSERT_FATAL(stack != NULL);
    CU_ASSERT_FATAL(leaves != NULL);
    for (u = 0; u < num_nodes; u++) {
        if (tree->num_children[u] == 0 && u >= n) {
            ret = sparse_tree_get_leaf_list(tree, u, &head, &tail);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(head, NULL);
            CU_ASSERT_EQUAL(tail, NULL);
        } else {
            stack_top = 0;
            num_leaves = 0;
            stack[stack_top] = u;
            while (stack_top >= 0) {
                v = stack[stack_top];
                stack_top--;
                if (v < n) {
                    leaves[num_leaves] = v;
                    num_leaves++;
                }
                for (j = tree->num_children[v] - 1; j >= 0; j--) {
                    stack_top++;
                    stack[stack_top] = tree->children[v][j];
                }
            }
            ret = sparse_tree_get_num_leaves(tree, u, &tmp);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(num_leaves, tmp);
            ret = sparse_tree_get_leaf_list(tree, u, &head, &tail);
            CU_ASSERT_EQUAL(ret, 0);
            z = head;
            j = 0;
            while (1) {
                CU_ASSERT_TRUE_FATAL(j < num_leaves);
                CU_ASSERT_EQUAL(z->node, leaves[j]);
                j++;
                if (z == tail) {
                    break;
                }
                z = z->next;
            }
            CU_ASSERT_EQUAL(j, num_leaves);
        }
    }
    free(stack);
    free(leaves);
}

static void
verify_leaf_sets(tree_sequence_t *ts)
{
    int ret;
    sparse_tree_t t;

    ret = sparse_tree_alloc(&t, ts, MSP_LEAF_COUNTS|MSP_LEAF_LISTS);
    CU_ASSERT_EQUAL(ret, 0);

    for (ret = sparse_tree_first(&t); ret == 1; ret = sparse_tree_next(&t)) {
        verify_leaf_sets_for_tree(&t);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (ret = sparse_tree_last(&t); ret == 1; ret = sparse_tree_prev(&t)) {
        verify_leaf_sets_for_tree(&t);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    sparse_tree_free(&t);
}

static sparse_tree_t *
get_tree_list(tree_sequence_t *ts)
{
    int ret;
    sparse_tree_t t, *trees;
    size_t num_trees;

    num_trees = tree_sequence_get_num_trees(ts);
    ret = sparse_tree_alloc(&t, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    trees = malloc(num_trees * sizeof(sparse_tree_t));
    CU_ASSERT_FATAL(trees != NULL);
    for (ret = sparse_tree_first(&t); ret == 1; ret = sparse_tree_next(&t)) {
        CU_ASSERT_FATAL(t.index < num_trees);
        ret = sparse_tree_alloc(&trees[t.index], ts, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = sparse_tree_copy(&trees[t.index], &t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = sparse_tree_equal(&trees[t.index], &t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        /* Make sure the left and right coordinates are also OK */
        CU_ASSERT_DOUBLE_EQUAL(trees[t.index].left, t.left, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(trees[t.index].right, t.right, 1e-6);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = sparse_tree_free(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    return trees;
}

static void
verify_tree_equals(tree_sequence_t *ts)
{
    int ret;
    sparse_tree_t *trees, t;
    size_t j, k;
    tree_sequence_t *other_ts = get_example_tree_sequence(
            10, 0, 100, 100.0, 1.0, 1.0, 0, NULL, MSP_ALPHABET_BINARY);
    int flags[] = {0, MSP_LEAF_LISTS, MSP_LEAF_COUNTS,
        MSP_LEAF_LISTS | MSP_LEAF_COUNTS};

    trees = get_tree_list(ts);
    ret = sparse_tree_alloc(&t, other_ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < tree_sequence_get_num_trees(ts); j++) {
        ret = sparse_tree_equal(&t, &trees[j]);
        CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
        for (k = 0; k < tree_sequence_get_num_trees(ts); k++) {
            ret = sparse_tree_equal(&trees[j], &trees[k]);
            if (j == k) {
                CU_ASSERT_EQUAL_FATAL(ret, 0);
            } else {
                CU_ASSERT_EQUAL_FATAL(ret, 1);
            }
        }
    }
    ret = sparse_tree_free(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < sizeof(flags) / sizeof(int); j++) {
        ret = sparse_tree_alloc(&t, ts, flags[j]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        for (ret = sparse_tree_first(&t); ret == 1;
                ret = sparse_tree_next(&t)) {
            for (k = 0; k < tree_sequence_get_num_trees(ts); k++) {
                ret = sparse_tree_equal(&t, &trees[k]);
                if (t.index == k) {
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                } else {
                    CU_ASSERT_EQUAL_FATAL(ret, 1);
                }
            }
        }
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = sparse_tree_free(&t);
        CU_ASSERT_EQUAL(ret, 0);
    }
    for (j = 0; j < tree_sequence_get_num_trees(ts); j++) {
        ret = sparse_tree_free(&trees[j]);
    }
    free(trees);
    tree_sequence_free(other_ts);
    free(other_ts);
}

static void
verify_tree_next_prev(tree_sequence_t *ts)
{
    int ret;
    sparse_tree_t *trees, t;
    size_t j;
    size_t num_trees = tree_sequence_get_num_trees(ts);

    trees = get_tree_list(ts);
    ret = sparse_tree_alloc(&t, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Single forward pass */
    j = 0;
    for (ret = sparse_tree_first(&t); ret == 1; ret = sparse_tree_next(&t)) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);

    /* Single reverse pass */
    j = num_trees;
    for (ret = sparse_tree_last(&t); ret == 1; ret = sparse_tree_prev(&t)) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);

    /* Full forward, then reverse */
    j = 0;
    for (ret = sparse_tree_first(&t); ret == 1; ret = sparse_tree_next(&t)) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);
    j--;
    while ((ret = sparse_tree_prev(&t)) == 1) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);
    CU_ASSERT_EQUAL_FATAL(t.index, 0);
    /* Calling prev should return 0 and have no effect. */
    for (j = 0; j < 10; j++) {
        ret = sparse_tree_prev(&t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(t.index, 0);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    /* Full reverse then forward */
    j = num_trees;
    for (ret = sparse_tree_last(&t); ret == 1; ret = sparse_tree_prev(&t)) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);
    j++;
    while ((ret = sparse_tree_next(&t)) == 1) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);
    CU_ASSERT_EQUAL_FATAL(t.index, num_trees - 1);
    /* Calling next should return 0 and have no effect. */
    for (j = 0; j < 10; j++) {
        ret = sparse_tree_next(&t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(t.index, num_trees - 1);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    /* Do a zigzagging traversal */
    ret = sparse_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 1; j < GSL_MIN(10, num_trees / 2); j++) {
        while (t.index < num_trees - j) {
            ret = sparse_tree_next(&t);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        CU_ASSERT_EQUAL_FATAL(t.index, num_trees - j);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        while (t.index > j) {
            ret = sparse_tree_prev(&t);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        CU_ASSERT_EQUAL_FATAL(t.index, j);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    /* Free the trees. */
    ret = sparse_tree_free(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < tree_sequence_get_num_trees(ts); j++) {
        ret = sparse_tree_free(&trees[j]);
    }
    free(trees);
}

static void
test_leaf_sets(void)
{
    leaf_count_test_t tests[] = {
        {0, 0, 1}, {0, 5, 2}, {0, 6, 3},
        {1, 4, 2}, {1, 5, 3}, {1, 6, 4}};
    uint32_t num_tests = 6;
    tree_sequence_t ts;

    tree_sequence_from_text(&ts, paper_ex_nodes, paper_ex_edgesets, NULL, NULL, NULL, NULL);
    verify_leaf_counts(&ts, num_tests, tests);
    verify_leaf_sets(&ts);

    tree_sequence_free(&ts);
}

static void
test_nonbinary_leaf_sets(void)
{
    leaf_count_test_t tests[] = {
        {0, 0, 1}, {0, 8, 4}, {0, 9, 5}, {0, 10, 3}, {0, 12, 8},
        {1, 5, 1}, {1, 8, 4}, {1, 9, 5}, {0, 10, 2}, {0, 11, 1}};
    uint32_t num_tests = 8;
    tree_sequence_t ts;

    tree_sequence_from_text(&ts, nonbinary_ex_nodes, nonbinary_ex_edgesets, NULL,
            NULL, NULL, NULL);
    verify_leaf_counts(&ts, num_tests, tests);
    verify_leaf_sets(&ts);

    tree_sequence_free(&ts);
}

static void
test_tree_sequence_bad_records(void)
{
    int ret = 0;
    tree_sequence_t ts;
    node_table_t node_table;
    edgeset_table_t edgeset_table;
    uint32_t num_trees = 3;
    node_id_t parents[] = {
        6, 5, 8, 5, MSP_NULL_NODE, 6, 8, MSP_NULL_NODE, MSP_NULL_NODE,
        6, 5, 4, 4, 5, 6, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
        7, 5, 4, 4, 5, 7, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
    };

    ret = node_table_alloc(&node_table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = edgeset_table_alloc(&edgeset_table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    parse_nodes(paper_ex_nodes, &node_table);
    parse_edgesets(paper_ex_edgesets, &edgeset_table);

    /* Make sure we have a good set of records */
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    verify_trees(&ts, num_trees, parents);
    tree_sequence_free(&ts);

    /* Left value greater than right */
    edgeset_table.left[0] = 10.0;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_RECORD_INTERVAL);
    tree_sequence_free(&ts);
    edgeset_table.left[0] = 2.0;

    /* Children equal */
    edgeset_table.children[7] = 0;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_UNSORTED_CHILDREN);
    tree_sequence_free(&ts);
    edgeset_table.children[7] = 5;

    /* Children not sorted */
    edgeset_table.children[6] = 5;
    edgeset_table.children[7] = 0;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_UNSORTED_CHILDREN);
    tree_sequence_free(&ts);
    edgeset_table.children[6] = 0;
    edgeset_table.children[7] = 5;

    /* Make a gap between adjacent records */
    edgeset_table.right[1] = 1.0;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORD_NONMATCHING_RIGHT);
    tree_sequence_free(&ts);
    edgeset_table.right[1] = 2.0;

    /* Make a gap in the middle of the sequence */
    edgeset_table.left[0] = 7;
    edgeset_table.left[2] = 7;
    edgeset_table.right[3] = 2;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORD_NONMATCHING_RIGHT);
    tree_sequence_free(&ts);
    edgeset_table.left[0] = 2;
    edgeset_table.left[2] = 2;
    edgeset_table.right[3] = 7;

    /* Make a gap before the last tree */
    edgeset_table.left[4] = 8;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORD_NONMATCHING_RIGHT);
    tree_sequence_free(&ts);
    edgeset_table.left[4] = 7;

    /* Add an extra record to the first tree */
    edgeset_table.left[4] = 2;
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORD_NONMATCHING_RIGHT);
    tree_sequence_free(&ts);
    edgeset_table.left[4] = 7;

    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_tables_tmp(&ts, &node_table, &edgeset_table, NULL,
            NULL, NULL, 0, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    verify_trees(&ts, num_trees, parents);
    tree_sequence_free(&ts);

    edgeset_table_free(&edgeset_table);
    node_table_free(&node_table);
}

static void
verify_tree_diffs(tree_sequence_t *ts)
{
    int ret;
    tree_diff_iterator_t iter;
    sparse_tree_t tree;
    node_record_t *record, *records_out, *records_in;
    size_t num_nodes = tree_sequence_get_num_nodes(ts);
    size_t j, k, num_in, num_out, num_trees;
    double length, t, x;
    node_t node;
    node_id_t u;
    node_id_t *pi = malloc(num_nodes * sizeof(node_id_t));
    double *tau = malloc(num_nodes * sizeof(double));
    int first_tree;

    CU_ASSERT_FATAL(pi != NULL);
    CU_ASSERT_FATAL(tau != NULL);
    for (j = 0; j < num_nodes; j++) {
        pi[j] = MSP_NULL_NODE;
        tau[j] = 0.0;
    }
    ret = tree_diff_iterator_alloc(&iter, ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = sparse_tree_alloc(&tree, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    tree_diff_iterator_print_state(&iter, _devnull);
    /* FIXME general samples will break this */
    for (j = 0; j < tree_sequence_get_sample_size(ts); j++) {
        ret = tree_sequence_get_node(ts, j, &node);
        CU_ASSERT_EQUAL(ret, 0);
        tau[j] = node.time;
    }

    first_tree = 1;
    x = 0.0;
    num_trees = 0;
    while ((ret = tree_diff_iterator_next(
                &iter, &length, &records_out, &records_in)) == 1) {
        tree_diff_iterator_print_state(&iter, _devnull);
        num_trees++;
        record = records_out;
        num_out = 0;
        while (record != NULL) {
            for (k = 0; k < record->num_children; k++) {
                pi[record->children[k]] = MSP_NULL_NODE;
            }
            tau[record->node] = 0;
            num_out += record->num_children - 1;
            record = record->next;
        }
        record = records_in;
        num_in = 0;
        while (record != NULL) {
            for (k = 0; k < record->num_children; k++) {
                pi[record->children[k]] = record->node;
            }
            tau[record->node] = record->time;
            num_in += record->num_children - 1;
            record = record->next;
        }
        if (first_tree) {
            CU_ASSERT_EQUAL(num_in, tree_sequence_get_sample_size(ts) - 1);
        } else {
            CU_ASSERT_EQUAL(num_in, num_out);
        }
        /* Now check against the sparse tree iterator. */
        for (j = 0; j < num_nodes; j++) {
            ret = sparse_tree_get_parent(&tree, j, &u);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(pi[j], u);
            ret = sparse_tree_get_time(&tree, j, &t);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(tau[j], t);
        }
        CU_ASSERT_EQUAL(tree.left, x);
        x += length;
        CU_ASSERT_EQUAL(tree.right, x);
        first_tree = 0;
        ret = sparse_tree_next(&tree);
        if (num_trees < tree_sequence_get_num_trees(ts)) {
            CU_ASSERT_EQUAL(ret, 1);
        } else {
            CU_ASSERT_EQUAL(ret, 0);
        }
    }
    if (num_trees != tree_sequence_get_num_trees(ts)) {
        printf("ERROR: %d %d\n", (int) num_trees, (int) tree_sequence_get_num_trees(ts));
    }
    CU_ASSERT_EQUAL(num_trees, tree_sequence_get_num_trees(ts));
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = sparse_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_diff_iterator_free(&iter);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = sparse_tree_free(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    free(pi);
    free(tau);
}

static void
test_tree_sequence_diff_iter(void)
{
    int ret;
    tree_sequence_t ts;

    tree_sequence_from_text(&ts, paper_ex_nodes, paper_ex_edgesets, NULL, NULL, NULL, NULL);
    verify_tree_diffs(&ts);

    ret = tree_sequence_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_nonbinary_tree_sequence_diff_iter(void)
{
    int ret;
    tree_sequence_t ts;

    tree_sequence_from_text(&ts, nonbinary_ex_nodes, nonbinary_ex_edgesets, NULL,
            NULL, NULL, NULL);
    verify_tree_diffs(&ts);

    ret = tree_sequence_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_unary_tree_sequence_diff_iter(void)
{
    int ret;
    tree_sequence_t ts;

    tree_sequence_from_text(&ts, unary_ex_nodes, unary_ex_edgesets, NULL,
            NULL, NULL, NULL);
    verify_tree_diffs(&ts);

    ret = tree_sequence_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_diff_iter_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_tree_diffs(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_tree_iter_from_examples(void)
{

    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_trees_consistent(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_leaf_sets_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_leaf_sets(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_tree_equals_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_tree_equals(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_next_prev_from_examples(void)
{

    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_tree_next_prev(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}


static void
test_hapgen_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_hapgen(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
verify_ld(tree_sequence_t *ts)
{
    int ret;
    size_t num_sites = tree_sequence_get_num_sites(ts);
    site_t *sites = malloc(num_sites * sizeof(site_t));
    ld_calc_t ld_calc;
    double *r2, *r2_prime, x;
    size_t j, num_r2_values;
    double eps = 1e-6;

    r2 = malloc(num_sites * sizeof(double));
    r2_prime = malloc(num_sites * sizeof(double));
    CU_ASSERT_FATAL(r2 != NULL);

    ret = ld_calc_alloc(&ld_calc, ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ld_calc_print_state(&ld_calc, _devnull);


    for (j = 0; j < num_sites; j++) {
        ret = tree_sequence_get_site(ts, j, sites + j);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = ld_calc_get_r2(&ld_calc, j, j, &x);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL_FATAL(x, 1.0, eps);
    }

    if (num_sites > 0) {
        /* Some checks in the forward direction */
        ret = ld_calc_get_r2_array(&ld_calc, 0, MSP_DIR_FORWARD,
                num_sites, DBL_MAX, r2, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, num_sites - 1);
        ld_calc_print_state(&ld_calc, _devnull);

        ret = ld_calc_get_r2_array(&ld_calc, num_sites - 2, MSP_DIR_FORWARD,
                num_sites, DBL_MAX, r2_prime, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, 1);
        ld_calc_print_state(&ld_calc, _devnull);

        ret = ld_calc_get_r2_array(&ld_calc, 0, MSP_DIR_FORWARD,
                num_sites, DBL_MAX, r2_prime, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, num_sites - 1);
        ld_calc_print_state(&ld_calc, _devnull);

        for (j = 0; j < num_r2_values; j++) {
            CU_ASSERT_EQUAL_FATAL(r2[j], r2_prime[j]);
            ret = ld_calc_get_r2(&ld_calc, 0, j + 1, &x);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_DOUBLE_EQUAL_FATAL(r2[j], x, eps);
        }

        /* Some checks in the reverse direction */
        ret = ld_calc_get_r2_array(&ld_calc, num_sites - 1,
                MSP_DIR_REVERSE, num_sites, DBL_MAX,
                r2, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, num_sites - 1);
        ld_calc_print_state(&ld_calc, _devnull);

        ret = ld_calc_get_r2_array(&ld_calc, 1, MSP_DIR_REVERSE,
                num_sites, DBL_MAX, r2_prime, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, 1);
        ld_calc_print_state(&ld_calc, _devnull);

        ret = ld_calc_get_r2_array(&ld_calc, num_sites - 1,
                MSP_DIR_REVERSE, num_sites, DBL_MAX,
                r2_prime, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, num_sites - 1);
        ld_calc_print_state(&ld_calc, _devnull);

        for (j = 0; j < num_r2_values; j++) {
            CU_ASSERT_EQUAL_FATAL(r2[j], r2_prime[j]);
            ret = ld_calc_get_r2(&ld_calc, num_sites - 1,
                    num_sites - j - 2, &x);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_DOUBLE_EQUAL_FATAL(r2[j], x, eps);
        }

        /* Check some error conditions */
        ret = ld_calc_get_r2_array(&ld_calc, 0, 0, num_sites, DBL_MAX,
            r2, &num_r2_values);
        CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    }

    if (num_sites > 3) {
        /* Check for some basic distance calculations */
        j = num_sites / 2;
        x = sites[j + 1].position - sites[j].position;
        ret = ld_calc_get_r2_array(&ld_calc, j, MSP_DIR_FORWARD, num_sites,
                x, r2, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, 1);

        x = sites[j].position - sites[j - 1].position;
        ret = ld_calc_get_r2_array(&ld_calc, j, MSP_DIR_REVERSE, num_sites,
                x, r2, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, 1);
    }

    /* Check some error conditions */
    for (j = num_sites; j < num_sites + 2; j++) {
        ret = ld_calc_get_r2_array(&ld_calc, j, MSP_DIR_FORWARD,
                num_sites, DBL_MAX, r2, &num_r2_values);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = ld_calc_get_r2(&ld_calc, j, 0, r2);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = ld_calc_get_r2(&ld_calc, 0, j, r2);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
    }

    ld_calc_free(&ld_calc);
    free(r2);
    free(r2_prime);
    free(sites);
}

static void
test_ld_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        /* j = 4 corresponds to the recurrent mutation case where we
         * trigger an assert */
        if (j != 4) {
            verify_ld(examples[j]);
        }
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_vargen_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        if (j == 4) {
            printf("\n\nFIXME multiple mutation vargen\n");
            tree_sequence_free(examples[j]);
            free(examples[j]);
            continue;
        }
        verify_vargen(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_stats_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        if (j == 4) {
            printf("\n\nFIXME multiple mutation PI\n");
            tree_sequence_free(examples[j]);
            free(examples[j]);
            continue;
        }
        verify_stats(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
verify_simplify_errors(tree_sequence_t *ts)
{
    int ret;
    uint32_t n = tree_sequence_get_sample_size(ts);
    tree_sequence_t subset;
    node_id_t sample[] = {0, 1, 2, 3};

    ret = tree_sequence_simplify(ts, sample, 0, 0, &subset);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = tree_sequence_simplify(ts, sample, 1, 0, &subset);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    sample[1] = n;
    ret = tree_sequence_simplify(ts, sample, 2, 0, &subset);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_SAMPLES);
    sample[0] = 0;
    sample[1] = 0;
    ret = tree_sequence_simplify(ts, sample, 2, 0, &subset);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_DUPLICATE_SAMPLE);

}

static void
test_simplify_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_simplify(examples[j]);
        verify_simplify_errors(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
verify_newick(tree_sequence_t *ts, bool should_fail)
{
    newick_converter_t nc;
    double length;
    char *tree;
    int ret;

    ret = newick_converter_alloc(&nc, ts, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    while ((ret = newick_converter_next(&nc, &length, &tree)) == 1) {
        CU_ASSERT(length > 0);
        newick_converter_print_state(&nc, _devnull);
        CU_ASSERT_FATAL(tree != NULL);
        CU_ASSERT(strlen(tree) > 0);
    }
    if (should_fail) {
        CU_ASSERT_EQUAL(ret, MSP_ERR_NONBINARY_NEWICK);
    } else {
        CU_ASSERT_EQUAL(ret, 0);
    }
    newick_converter_free(&nc);
}

static void
test_newick_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(0);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_newick(examples[j], false);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);

    examples = get_example_nonbinary_tree_sequences();
    for (j = 0; examples[j] != NULL; j++) {
        verify_newick(examples[j], true);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
verify_tree_sequences_equal(tree_sequence_t *ts1, tree_sequence_t *ts2,
        bool check_migrations, bool check_mutations,
        bool check_provenance_strings)
{
    int ret, err1, err2;
    size_t j, nps1, nps2;
    edgeset_t r1, r2;
    node_t n1, n2;
    migration_t m1, m2;
    char **ps1, **ps2;
    size_t num_mutations = tree_sequence_get_num_mutations(ts1);
    site_t site_1, site_2;
    mutation_t mutation_1, mutation_2;
    sparse_tree_t t1, t2;

    CU_ASSERT_EQUAL(
        tree_sequence_get_sample_size(ts1),
        tree_sequence_get_sample_size(ts2))
    CU_ASSERT_EQUAL(
        tree_sequence_get_sequence_length(ts1),
        tree_sequence_get_sequence_length(ts2))
    CU_ASSERT_EQUAL(
        tree_sequence_get_num_edgesets(ts1),
        tree_sequence_get_num_edgesets(ts2));
    CU_ASSERT_EQUAL(
        tree_sequence_get_num_nodes(ts1),
        tree_sequence_get_num_nodes(ts2));
    CU_ASSERT_EQUAL(
        tree_sequence_get_num_trees(ts1),
        tree_sequence_get_num_trees(ts2));

    for (j = 0; j < tree_sequence_get_num_nodes(ts1); j++) {
        ret = tree_sequence_get_node(ts1, j, &n1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_get_node(ts2, j, &n2);
        CU_ASSERT_EQUAL(ret, 0);
        verify_nodes_equal(&n1, &n2);
    }
    for (j = 0; j < tree_sequence_get_num_edgesets(ts1); j++) {
        ret = tree_sequence_get_edgeset(ts1, j, &r1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_get_edgeset(ts2, j, &r2);
        CU_ASSERT_EQUAL(ret, 0);
        verify_edgesets_equal(&r1, &r2, 1.0);
    }
    if (check_mutations) {
        CU_ASSERT_EQUAL_FATAL(
            tree_sequence_get_num_sites(ts1),
            tree_sequence_get_num_sites(ts2));
        for (j = 0; j < tree_sequence_get_num_sites(ts1); j++) {
            ret = tree_sequence_get_site(ts1, j, &site_1);
            CU_ASSERT_EQUAL(ret, 0);
            ret = tree_sequence_get_site(ts2, j, &site_2);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(site_1.position, site_2.position);
            CU_ASSERT_STRING_EQUAL(site_1.ancestral_state, site_2.ancestral_state);
        }
        CU_ASSERT_EQUAL_FATAL(
            tree_sequence_get_num_mutations(ts1),
            tree_sequence_get_num_mutations(ts2));
        for (j = 0; j < num_mutations; j++) {
            ret = tree_sequence_get_mutation(ts1, j, &mutation_1);
            CU_ASSERT_EQUAL(ret, 0);
            ret = tree_sequence_get_mutation(ts2, j, &mutation_2);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(mutation_1.id, j);
            CU_ASSERT_EQUAL(mutation_1.id, mutation_2.id);
            CU_ASSERT_EQUAL(mutation_1.index, j);
            CU_ASSERT_EQUAL(mutation_1.index, mutation_2.index);
            CU_ASSERT_EQUAL(mutation_1.site, mutation_2.site);
            CU_ASSERT_EQUAL(mutation_1.node, mutation_2.node);
            CU_ASSERT_STRING_EQUAL(mutation_1.derived_state, mutation_2.derived_state);
        }
    }
    if (check_migrations) {
        CU_ASSERT_EQUAL_FATAL(
            tree_sequence_get_num_migrations(ts1),
            tree_sequence_get_num_migrations(ts2));
        for (j = 0; j < tree_sequence_get_num_migrations(ts1); j++) {
            ret = tree_sequence_get_migration(ts1, j, &m1);
            CU_ASSERT_EQUAL(ret, 0);
            ret = tree_sequence_get_migration(ts2, j, &m2);
            CU_ASSERT_EQUAL(ret, 0);
            verify_migrations_equal(&m1, &m2, 1.0);
        }
    }
    if (check_provenance_strings) {
        ret = tree_sequence_get_provenance_strings(ts1, &nps1, &ps1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_get_provenance_strings(ts2, &nps2, &ps2);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(nps1, nps2);
        for (j = 0; j < nps1; j++) {
            CU_ASSERT_FATAL(ps1[j] != NULL);
            CU_ASSERT_FATAL(ps2[j] != NULL);
            CU_ASSERT_STRING_EQUAL(ps1[j], ps2[j]);
        }
    }
    ret = sparse_tree_alloc(&t1, ts1, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_alloc(&t2, ts2, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_first(&t1);
    CU_ASSERT_EQUAL(ret, 1);
    ret = sparse_tree_first(&t2);
    CU_ASSERT_EQUAL(ret, 1);
    while (1) {
        err1 = sparse_tree_next(&t1);
        err2 = sparse_tree_next(&t2);
        CU_ASSERT_EQUAL_FATAL(err1, err2);
        if (err1 != 1) {
            break;
        }
    }
    sparse_tree_free(&t1);
    sparse_tree_free(&t2);
}

static void
verify_empty_tree_sequence(tree_sequence_t *ts)
{
    CU_ASSERT_EQUAL(tree_sequence_get_num_edgesets(ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_migrations(ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(ts), 0);
    verify_trees_consistent(ts);
    verify_ld(ts);
    verify_stats(ts);
    verify_hapgen(ts);
    verify_vargen(ts);
    verify_newick(ts, false);
    verify_vcf_converter(ts, 1);
}

static void
test_save_empty_hdf5(void)
{
    int ret;
    tree_sequence_t ts1, ts2;

    ret = tree_sequence_initialise(&ts1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_initialise(&ts2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_dump(&ts1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_empty_tree_sequence(&ts1);
    ret = tree_sequence_load(&ts2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_empty_tree_sequence(&ts2);
}

static void
test_save_hdf5(void)
{
    int ret;
    size_t j, k;
    tree_sequence_t **examples = get_example_tree_sequences(1);
    tree_sequence_t ts2;
    tree_sequence_t *ts1;
    int dump_flags[] = {0, MSP_DUMP_ZLIB_COMPRESSION};

    CU_ASSERT_FATAL(examples != NULL);

    for (j = 0; examples[j] != NULL; j++) {
        ts1 = examples[j];
        for (k = 0; k < sizeof(dump_flags) / sizeof(int); k++) {
            ret = tree_sequence_dump(ts1, _tmp_file_name, dump_flags[k]);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            ret = tree_sequence_initialise(&ts2);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            ret = tree_sequence_load(&ts2, _tmp_file_name, MSP_LOAD_EXTENDED_CHECKS);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            /* FIXME storing migrations */
            verify_tree_sequences_equal(ts1, &ts2, false, true, true);
            tree_sequence_print_state(&ts2, _devnull);
            verify_hapgen(&ts2);
            verify_vargen(&ts2);
            tree_sequence_free(&ts2);
        }
        tree_sequence_free(ts1);
        free(ts1);
    }
    free(examples);
}

static void
test_dump_tables(void)
{
    int ret;
    tree_sequence_t **examples = get_example_tree_sequences(1);
    tree_sequence_t ts2;
    tree_sequence_t *ts1;
    size_t j, num_provenance_strings;
    size_t alloc_size = 8192;
    char **provenance_strings;
    node_table_t nodes;
    edgeset_table_t edgesets;
    migration_table_t migrations;
    site_table_t sites;
    mutation_table_t mutations;

    ret = node_table_alloc(&nodes, alloc_size, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = edgeset_table_alloc(&edgesets, alloc_size, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = migration_table_alloc(&migrations, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = site_table_alloc(&sites, alloc_size, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutation_table_alloc(&mutations, alloc_size, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_FATAL(examples != NULL);

    for (j = 0; examples[j] != NULL; j++) {
        ts1 = examples[j];

        ret = tree_sequence_dump_tables_tmp(ts1, NULL, &edgesets,
                &migrations, &sites, &mutations, &num_provenance_strings,
                &provenance_strings);
        CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
        ret = tree_sequence_dump_tables_tmp(ts1, &nodes, NULL,
                &migrations, &sites, &mutations, &num_provenance_strings,
                &provenance_strings);
        CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
        ret = tree_sequence_load_tables_tmp(&ts2, NULL, &edgesets,
                &migrations, &sites, &mutations, num_provenance_strings,
                provenance_strings);
        CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
        ret = tree_sequence_load_tables_tmp(&ts2, &nodes, NULL,
                &migrations, &sites, &mutations, num_provenance_strings,
                provenance_strings);
        CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);

        ret = tree_sequence_dump_tables_tmp(ts1, &nodes, &edgesets,
                &migrations, &sites, &mutations, &num_provenance_strings,
                &provenance_strings);
        ret = tree_sequence_initialise(&ts2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tree_sequence_load_tables_tmp(&ts2, &nodes, &edgesets,
                &migrations, &sites, &mutations, num_provenance_strings,
                provenance_strings);
        verify_tree_sequences_equal(ts1, &ts2, true, true, true);
        tree_sequence_print_state(&ts2, _devnull);
        tree_sequence_free(&ts2);

        ret = tree_sequence_dump_tables_tmp(ts1, &nodes, &edgesets,
                NULL, &sites, &mutations, &num_provenance_strings,
                &provenance_strings);
        ret = tree_sequence_initialise(&ts2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tree_sequence_load_tables_tmp(&ts2, &nodes, &edgesets,
                NULL, &sites, &mutations, num_provenance_strings,
                provenance_strings);
        verify_tree_sequences_equal(ts1, &ts2, false, true, true);
        CU_ASSERT_EQUAL_FATAL(
            tree_sequence_get_num_migrations(&ts2), 0);
        tree_sequence_free(&ts2);

        ret = tree_sequence_dump_tables_tmp(ts1, &nodes, &edgesets,
                &migrations, NULL, NULL, &num_provenance_strings,
                &provenance_strings);
        ret = tree_sequence_initialise(&ts2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tree_sequence_load_tables_tmp(&ts2, &nodes, &edgesets,
                &migrations, NULL, NULL, num_provenance_strings,
                provenance_strings);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        verify_tree_sequences_equal(ts1, &ts2, true, false, true);
        CU_ASSERT_EQUAL_FATAL(
            tree_sequence_get_num_mutations(&ts2), 0);
        tree_sequence_free(&ts2);

        tree_sequence_free(ts1);
        free(ts1);
    }
    free(examples);
    node_table_free(&nodes);
    edgeset_table_free(&edgesets);
    migration_table_free(&migrations);
    site_table_free(&sites);
    mutation_table_free(&mutations);
}

static void
test_dump_tables_hdf5(void)
{
    int ret;
    size_t k, num_provenance_strings;
    tree_sequence_t *ts1, ts2, ts3, **examples;
    size_t alloc_size = 8192;
    char **provenance_strings;
    node_table_t nodes;
    edgeset_table_t edgesets;
    migration_table_t migrations;
    site_table_t sites;
    mutation_table_t mutations;

    ret = node_table_alloc(&nodes, alloc_size, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = edgeset_table_alloc(&edgesets, alloc_size, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = migration_table_alloc(&migrations, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = site_table_alloc(&sites, alloc_size, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = mutation_table_alloc(&mutations, alloc_size, alloc_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    examples = get_example_tree_sequences(1);
    for (k = 0; examples[k] != NULL; k++) {
        ts1 = examples[k];
        CU_ASSERT_FATAL(ts1 != NULL);
        ret = tree_sequence_dump_tables_tmp(ts1, &nodes, &edgesets,
                &migrations, &sites, &mutations, &num_provenance_strings,
                &provenance_strings);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tree_sequence_initialise(&ts2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tree_sequence_load_tables_tmp(&ts2, &nodes, &edgesets,
                &migrations, &sites, &mutations, num_provenance_strings,
                provenance_strings);
        ret = tree_sequence_dump(&ts2, _tmp_file_name, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tree_sequence_initialise(&ts3);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tree_sequence_load(&ts3, _tmp_file_name, MSP_LOAD_EXTENDED_CHECKS);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        /* FIXME storing migrations */
        verify_tree_sequences_equal(ts1, &ts3, false, true, true);
        tree_sequence_print_state(&ts2, _devnull);

        tree_sequence_free(&ts2);
        tree_sequence_free(&ts3);
        tree_sequence_free(ts1);
        free(ts1);
    }
    free(examples);
    node_table_free(&nodes);
    edgeset_table_free(&edgesets);
    migration_table_free(&migrations);
    site_table_free(&sites);
    mutation_table_free(&mutations);
}


static void
test_strerror(void)
{
    int ret;
    int j;
    const char *msg;
    int max_error_code = 1024; /* totally arbitrary */
    FILE *f;
    tree_sequence_t ts;

    for (j = 0; j < max_error_code; j++) {
        msg = msp_strerror(-j);
        CU_ASSERT_FATAL(msg != NULL);
        if (-j == MSP_ERR_HDF5) {
            /* There is no HDF5 error, so... */
            CU_ASSERT_EQUAL(strlen(msg), 0);
        } else {
            CU_ASSERT(strlen(msg) > 0);
        }
    }
    /* Provoke an HDF5 error */
    ret = tree_sequence_initialise(&ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load(&ts, "/file/does/not/exist", 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_HDF5);
    msg = msp_strerror(ret);
    CU_ASSERT_FATAL(msg != NULL);
    CU_ASSERT(strlen(msg) > 0);
    /* Provoke an IO error */
    f = fopen("/file/does/not/exist", "r");
    CU_ASSERT_EQUAL_FATAL(f, NULL);
    msg = msp_strerror(MSP_ERR_IO);
    CU_ASSERT_FATAL(msg != NULL);
    CU_ASSERT_STRING_EQUAL(msg, strerror(errno));
}

static void
test_node_table(void)
{
    int ret;
    node_table_t table;
    size_t num_rows = 100;
    size_t j;
    uint32_t *flags;
    population_id_t *population;
    double *time;
    char *name;
    uint32_t *name_length;
    const char *test_name = "test";
    size_t test_name_length = 4;
    char name_copy[test_name_length + 1];

    name_copy[test_name_length] = '\0';
    ret = node_table_alloc(&table, 0, 1);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = node_table_alloc(&table, 1, 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    ret = node_table_alloc(&table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    node_table_print_state(&table, _devnull);
    ret = node_table_add_row(&table, 0, 0, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    for (j = 0; j < num_rows; j++) {
        ret = node_table_add_row(&table, j, j, j, test_name);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(table.flags[j], j);
        CU_ASSERT_EQUAL(table.time[j], j);
        CU_ASSERT_EQUAL(table.population[j], j);
        CU_ASSERT_EQUAL(table.num_rows, j + 1);
        CU_ASSERT_EQUAL(table.total_name_length, (j + 1) * test_name_length);
        /* check the name */
        memcpy(name_copy, table.name + j * test_name_length, test_name_length);
        CU_ASSERT_STRING_EQUAL(name_copy, test_name);
    }
    node_table_print_state(&table, _devnull);
    node_table_reset(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.total_name_length, 0);

    num_rows *= 2;
    flags = malloc(num_rows * sizeof(uint32_t));
    CU_ASSERT_FATAL(flags != NULL);
    memset(flags, 1, num_rows * sizeof(uint32_t));
    population = malloc(num_rows * sizeof(uint32_t));
    CU_ASSERT_FATAL(population != NULL);
    memset(population, 2, num_rows * sizeof(uint32_t));
    time = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(time != NULL);
    memset(time, 0, num_rows * sizeof(double));
    name = malloc(num_rows * sizeof(char));
    memset(name, 'a', num_rows * sizeof(char));
    CU_ASSERT_FATAL(name != NULL);
    name_length = malloc(num_rows * sizeof(uint32_t));
    CU_ASSERT_FATAL(name_length != NULL);
    for (j = 0; j < num_rows; j++) {
        name_length[j] = 1;
    }
    ret = node_table_set_columns(&table, num_rows, flags, time, population,
            name, name_length);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.flags, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.population, population, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.name, name, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.name_length, name_length, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.total_name_length, num_rows);

    /* If population is NULL it should be set the -1. If name is NULL all names
     * should be set to the empty string. */
    num_rows = 10;
    memset(population, 0xff, num_rows * sizeof(uint32_t));
    ret = node_table_set_columns(&table, num_rows, flags, time, NULL, name, name_length);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.flags, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.population, population, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.name, name, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.name_length, name_length, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.total_name_length, num_rows);

    /* flags and time cannot be NULL */
    ret = node_table_set_columns(&table, num_rows, NULL, time, population, name, name_length);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = node_table_set_columns(&table, num_rows, flags, NULL, population, name, name_length);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = node_table_set_columns(&table, num_rows, flags, time, population, NULL, name_length);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = node_table_set_columns(&table, num_rows, flags, time, population, name, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    /* if name and name_length are both null, all names are zero length */
    num_rows = 10;
    memset(name_length, 0, num_rows * sizeof(uint32_t));
    ret = node_table_set_columns(&table, num_rows, flags, time, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.flags, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.name_length, name_length, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.total_name_length, 0);

    node_table_free(&table);
    free(flags);
    free(population);
    free(time);
    free(name);
    free(name_length);
}

static void
test_edgeset_table(void)
{
    int ret;
    edgeset_table_t table;
    size_t num_rows = 100;
    size_t max_children = 10;
    size_t j, k, total_children_length;
    node_id_t *parent, *children;
    list_len_t *children_length;
    double *left, *right;
    node_id_t c[max_children];

    ret = edgeset_table_alloc(&table, 0, 1);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = edgeset_table_alloc(&table, 1, 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    ret = edgeset_table_alloc(&table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    edgeset_table_print_state(&table, _devnull);

    /* Adding 0 children is an error */
    ret = edgeset_table_add_row(&table, 0, 0, 0, c, 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    memset(c, 0, max_children * sizeof(node_id_t));

    total_children_length = 0;
    for (j = 0; j < num_rows; j++) {
        k = GSL_MIN(j + 1, max_children);
        total_children_length += k;
        ret = edgeset_table_add_row(&table, j, j, j, c, k);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(table.left[j], j);
        CU_ASSERT_EQUAL(table.right[j], j);
        CU_ASSERT_EQUAL(table.parent[j], j);
        CU_ASSERT_EQUAL(table.children_length[j], k);
        CU_ASSERT_EQUAL(table.num_rows, j + 1);
        CU_ASSERT_EQUAL(table.total_children_length, total_children_length);
    }
    edgeset_table_print_state(&table, _devnull);
    edgeset_table_reset(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.total_children_length, 0);

    num_rows *= 2;
    left = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(left != NULL);
    memset(left, 0, num_rows * sizeof(double));
    right = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(right != NULL);
    memset(right, 0, num_rows * sizeof(double));
    parent = malloc(num_rows * sizeof(node_id_t));
    CU_ASSERT_FATAL(parent != NULL);
    memset(parent, 1, num_rows * sizeof(node_id_t));
    children = malloc(2 * num_rows * sizeof(node_id_t));
    children_length = malloc(num_rows * sizeof(list_len_t));
    CU_ASSERT_FATAL(children_length != NULL);
    for (j = 0; j < num_rows; j++) {
        children[2 * j] = j;
        children[2 * j + 1] = j + 1;
        children_length[j] = 2;
    }

    ret = edgeset_table_set_columns(&table, num_rows, left, right, parent,
            children, children_length);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.parent, parent, num_rows * sizeof(node_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.children, children, 2 * num_rows * sizeof(node_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.children_length, children_length,
                num_rows * sizeof(list_len_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.total_children_length, 2 * num_rows);

    /* Inputs cannot be NULL */
    ret = edgeset_table_set_columns(&table, num_rows, NULL, right, parent,
            children, children_length);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = edgeset_table_set_columns(&table, num_rows, left, NULL, parent,
            children, children_length);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = edgeset_table_set_columns(&table, num_rows, left, right, NULL,
            children, children_length);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = edgeset_table_set_columns(&table, num_rows, left, right, parent,
            NULL, children_length);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = edgeset_table_set_columns(&table, num_rows, left, right, parent,
            children, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    edgeset_table_free(&table);
    free(left);
    free(right);
    free(parent);
    free(children);
    free(children_length);
}

static void
test_site_table(void)
{
    int ret;
    site_table_t table;
    size_t num_rows, j;
    char *ancestral_state;
    double *position;
    uint32_t *ancestral_state_length;

    ret = site_table_alloc(&table, 0, 1);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = site_table_alloc(&table, 1, 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    ret = site_table_alloc(&table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    site_table_print_state(&table, _devnull);

    ret = site_table_add_row(&table, 0, "A", 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.position[0], 0);
    CU_ASSERT_EQUAL(table.ancestral_state_length[0], 1);
    CU_ASSERT_EQUAL(table.num_rows, 1);

    ret = site_table_add_row(&table, 1, "AA", 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.position[1], 1);
    CU_ASSERT_EQUAL(table.ancestral_state_length[1], 2);
    CU_ASSERT_EQUAL(table.num_rows, 2);

    ret = site_table_add_row(&table, 2, "A", 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.position[1], 1);
    CU_ASSERT_EQUAL(table.ancestral_state_length[2], 1);
    CU_ASSERT_EQUAL(table.num_rows, 3);
    CU_ASSERT_EQUAL(table.total_ancestral_state_length, 4);

    site_table_print_state(&table, _devnull);
    site_table_reset(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.total_ancestral_state_length, 0);

    num_rows = 100;
    position = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(position != NULL);
    ancestral_state = malloc(num_rows * sizeof(char));
    CU_ASSERT_FATAL(ancestral_state != NULL);
    ancestral_state_length = malloc(num_rows * sizeof(uint32_t));
    CU_ASSERT_FATAL(ancestral_state_length != NULL);

    for (j = 0; j < num_rows; j++) {
        position[j] = j;
        ancestral_state[j] = j;
        ancestral_state_length[j] = 1;
    }
    ret = site_table_set_columns(&table, num_rows, position, ancestral_state,
            ancestral_state_length);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.position, position,
                num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.ancestral_state, ancestral_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.total_ancestral_state_length, num_rows);

    /* Inputs cannot be NULL */
    ret = site_table_set_columns(&table, num_rows, NULL, ancestral_state,
            ancestral_state_length);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = site_table_set_columns(&table, num_rows, position, NULL, ancestral_state_length);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = site_table_set_columns(&table, num_rows, position, ancestral_state, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    site_table_free(&table);
    free(position);
    free(ancestral_state);
    free(ancestral_state_length);
}

static void
test_mutation_table(void)
{
    int ret;
    mutation_table_t table;
    size_t num_rows = 100;
    size_t max_len = 20;
    size_t j, k, len;
    node_id_t *node;
    site_id_t *site;
    char *derived_state;
    char c[max_len + 1];
    uint32_t *derived_state_length;

    for (j = 0; j < max_len; j++) {
        c[j] = j + 1;
    }
    ret = mutation_table_alloc(&table, 0, 1);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = mutation_table_alloc(&table, 1, 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    ret = mutation_table_alloc(&table, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    mutation_table_print_state(&table, _devnull);

    len = 0;
    for (j = 0; j < num_rows; j++) {
        k = GSL_MIN(j + 1, max_len);
        len += k;
        ret = mutation_table_add_row(&table, j, j, c, k);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(table.site[j], j);
        CU_ASSERT_EQUAL(table.node[j], j);
        CU_ASSERT_EQUAL(table.derived_state_length[j], k);
        CU_ASSERT_EQUAL(table.num_rows, j + 1);
        CU_ASSERT_EQUAL(table.total_derived_state_length, len);
    }
    mutation_table_print_state(&table, _devnull);
    mutation_table_reset(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.total_derived_state_length, 0);

    num_rows *= 2;
    site = malloc(num_rows * sizeof(site_id_t));
    CU_ASSERT_FATAL(site != NULL);
    node = malloc(num_rows * sizeof(node_id_t));
    CU_ASSERT_FATAL(node != NULL);
    derived_state = malloc(num_rows * sizeof(char));
    CU_ASSERT_FATAL(derived_state != NULL);
    derived_state_length = malloc(num_rows * sizeof(uint32_t));
    CU_ASSERT_FATAL(derived_state_length != NULL);

    for (j = 0; j < num_rows; j++) {
        node[j] = j;
        site[j] = j + 1;
        derived_state[j] = j;
        derived_state_length[j] = 1;
    }
    ret = mutation_table_set_columns(&table, num_rows, site, node, derived_state,
            derived_state_length);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.site, site, num_rows * sizeof(site_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.node, node, num_rows * sizeof(node_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.derived_state, derived_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.derived_state_length, derived_state_length,
                num_rows * sizeof(list_len_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.total_derived_state_length, num_rows);

    /* Inputs cannot be NULL */
    ret = mutation_table_set_columns(&table, num_rows, NULL, node, derived_state,
            derived_state_length);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = mutation_table_set_columns(&table, num_rows, site, NULL, derived_state,
            derived_state_length);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = mutation_table_set_columns(&table, num_rows, site, node, NULL,
            derived_state_length);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = mutation_table_set_columns(&table, num_rows, site, node, derived_state,
            NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    mutation_table_free(&table);
    free(site);
    free(node);
    free(derived_state);
    free(derived_state_length);
}

static void
test_migration_table(void)
{
    int ret;
    migration_table_t table;
    size_t num_rows = 100;
    size_t j;
    node_id_t *node;
    population_id_t *source, *dest;
    double *left, *right, *time;

    ret = migration_table_alloc(&table, 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    ret = migration_table_alloc(&table, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    migration_table_print_state(&table, _devnull);

    for (j = 0; j < num_rows; j++) {
        ret = migration_table_add_row(&table, j, j, j, j, j, j);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(table.left[j], j);
        CU_ASSERT_EQUAL(table.right[j], j);
        CU_ASSERT_EQUAL(table.node[j], j);
        CU_ASSERT_EQUAL(table.source[j], j);
        CU_ASSERT_EQUAL(table.dest[j], j);
        CU_ASSERT_EQUAL(table.time[j], j);
        CU_ASSERT_EQUAL(table.num_rows, j + 1);
    }
    migration_table_print_state(&table, _devnull);
    migration_table_reset(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);

    num_rows *= 2;
    left = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(left != NULL);
    memset(left, 1, num_rows * sizeof(double));
    right = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(right != NULL);
    memset(right, 2, num_rows * sizeof(double));
    time = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(time != NULL);
    memset(time, 3, num_rows * sizeof(double));
    node = malloc(num_rows * sizeof(node_id_t));
    CU_ASSERT_FATAL(node != NULL);
    memset(node, 4, num_rows * sizeof(node_id_t));
    source = malloc(num_rows * sizeof(population_id_t));
    CU_ASSERT_FATAL(source != NULL);
    memset(source, 5, num_rows * sizeof(population_id_t));
    dest = malloc(num_rows * sizeof(population_id_t));
    CU_ASSERT_FATAL(dest != NULL);
    memset(dest, 6, num_rows * sizeof(population_id_t));

    ret = migration_table_set_columns(&table, num_rows, left, right, node, source,
            dest, time);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.node, node, num_rows * sizeof(node_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.source, source, num_rows * sizeof(population_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.dest, dest, num_rows * sizeof(population_id_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);

    /* inputs cannot be NULL */
    ret = migration_table_set_columns(&table, num_rows, NULL, right, node, source,
            dest, time);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = migration_table_set_columns(&table, num_rows, left, NULL, node, source,
            dest, time);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = migration_table_set_columns(&table, num_rows, left, right, NULL, source,
            dest, time);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = migration_table_set_columns(&table, num_rows, left, right, node, NULL,
            dest, time);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = migration_table_set_columns(&table, num_rows, left, right, node, source,
            NULL, time);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = migration_table_set_columns(&table, num_rows, left, right, node, source,
            dest, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    migration_table_free(&table);
    free(left);
    free(right);
    free(time);
    free(node);
    free(source);
    free(dest);
}


static int
msprime_suite_init(void)
{
    int fd;
    static char template[] = "/tmp/msp_c_test_XXXXXX";

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
    /* Silence HDF5 errors */
    if (H5Eset_auto(H5E_DEFAULT, NULL, NULL) < 0) {
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
    fprintf(stderr, "CUnit error occured: %d: %s\n",
            CU_get_error(), CU_get_error_msg());
    exit(EXIT_FAILURE);
}

int
main(int argc, char **argv)
{
    int ret;
    CU_pTest test;
    CU_pSuite suite;
    CU_TestInfo tests[] = {
        {"test_fenwick_tree", test_fenwick},
        {"test_vcf", test_vcf},
        {"test_vcf_no_mutations", test_vcf_no_mutations},
        {"test_simple_recombination_map", test_simple_recomb_map},
        {"test_recombination_map_errors", test_recomb_map_errors},
        {"test_recombination_map_examples", test_recomb_map_examples},
        {"test_node_names", test_node_names},
        {"test_simplest_records", test_simplest_records},
        {"test_simplest_nonbinary_records", test_simplest_nonbinary_records},
        {"test_simplest_unary_records", test_simplest_unary_records},
        {"test_simplest_non_sample_leaf_records", test_simplest_non_sample_leaf_records},
        {"test_simplest_degenerate_multiple_root_records",
            test_simplest_degenerate_multiple_root_records},
        {"test_simplest_multiple_root_records", test_simplest_multiple_root_records},
        {"test_simplest_root_mutations", test_simplest_root_mutations},
        {"test_simplest_back_mutations", test_simplest_back_mutations},
        {"test_simplest_bad_records", test_simplest_bad_records},
        {"test_alphabet_detection", test_alphabet_detection},
        {"test_single_tree_good_records", test_single_tree_good_records},
        {"test_single_nonbinary_tree_good_records",
            test_single_nonbinary_tree_good_records},
        {"test_single_tree_bad_records", test_single_tree_bad_records},
        {"test_single_tree_good_mutations", test_single_tree_good_mutations},
        {"test_single_tree_bad_mutations", test_single_tree_bad_mutations},
        {"test_single_tree_iter", test_single_tree_iter},
        {"test_single_nonbinary_tree_iter", test_single_nonbinary_tree_iter},
        {"test_single_tree_iter_times", test_single_tree_iter_times},
        {"test_single_tree_hapgen_char_alphabet", test_single_tree_hapgen_char_alphabet},
        {"test_single_tree_hapgen_binary_alphabet", test_single_tree_hapgen_binary_alphabet},
        {"test_single_tree_vargen_char_alphabet", test_single_tree_vargen_char_alphabet},
        {"test_single_tree_vargen_binary_alphabet", test_single_tree_vargen_binary_alphabet},
        {"test_single_tree_simplify", test_single_tree_simplify},
        {"test_single_tree_inconsistent_mutations", test_single_tree_inconsistent_mutations},
        {"test_single_unary_tree_hapgen", test_single_unary_tree_hapgen},
        {"test_single_tree_mutgen", test_single_tree_mutgen},
        {"test_sparse_tree_errors", test_sparse_tree_errors},
        {"test_tree_sequence_iter", test_tree_sequence_iter},
        {"test_leaf_sets", test_leaf_sets},
        {"test_nonbinary_leaf_sets", test_nonbinary_leaf_sets},
        {"test_nonbinary_tree_sequence_iter", test_nonbinary_tree_sequence_iter},
        {"test_unary_tree_sequence_iter", test_unary_tree_sequence_iter},
        {"test_left_to_right_tree_sequence_iter", test_left_to_right_tree_sequence_iter},
        {"test_tree_sequence_bad_records", test_tree_sequence_bad_records},
        {"test_tree_sequence_diff_iter", test_tree_sequence_diff_iter},
        {"test_nonbinary_tree_sequence_diff_iter",
            test_nonbinary_tree_sequence_diff_iter},
        {"test_unary_tree_sequence_diff_iter",
            test_unary_tree_sequence_diff_iter},
        {"test_diff_iter_from_examples", test_diff_iter_from_examples},
        {"test_tree_iter_from_examples", test_tree_iter_from_examples},
        {"test_tree_equals_from_examples", test_tree_equals_from_examples},
        {"test_tree_next_and_prev_from_examples", test_next_prev_from_examples},
        {"test_leaf_sets_from_examples", test_leaf_sets_from_examples},
        {"test_hapgen_from_examples", test_hapgen_from_examples},
        {"test_vargen_from_examples", test_vargen_from_examples},
        {"test_newick_from_examples", test_newick_from_examples},
        {"test_stats_from_examples", test_stats_from_examples},
        {"test_ld_from_examples", test_ld_from_examples},
        {"test_simplify_from_examples", test_simplify_from_examples},
        {"test_save_empty_hdf5", test_save_empty_hdf5},
        {"test_save_hdf5", test_save_hdf5},
        {"test_dump_tables", test_dump_tables},
        {"test_dump_tables_hdf5", test_dump_tables_hdf5},
        {"test_single_locus_two_populations", test_single_locus_two_populations},
        {"test_many_populations", test_single_locus_many_populations},
        {"test_historical_samples", test_single_locus_historical_sample},
        {"test_simulator_getters/setters", test_simulator_getters_setters},
        {"test_model_errors", test_simulator_model_errors},
        {"test_demographic_events", test_simulator_demographic_events},
        {"test_single_locus_simulation", test_single_locus_simulation},
        {"test_simulation_memory_limit", test_simulation_memory_limit},
        {"test_multi_locus_simulation", test_multi_locus_simulation},
        {"test_simulation_replicates", test_simulation_replicates},
        {"test_bottleneck_simulation", test_bottleneck_simulation},
        {"test_multiple_mergers_simulation", test_multiple_mergers_simulation},
        {"test_large_bottleneck_simulation", test_large_bottleneck_simulation},
        {"test_error_messages", test_strerror},
        {"test_node_table", test_node_table},
        {"test_edgeset_table", test_edgeset_table},
        {"test_site_table", test_site_table},
        {"test_mutation_table", test_mutation_table},
        {"test_migration_table", test_migration_table},
        CU_TEST_INFO_NULL,
    };

    /* We use initialisers here as the struct definitions change between
     * versions of CUnit */
    CU_SuiteInfo suites[] = {
        {
            .pName = "msprime",
            .pInitFunc = msprime_suite_init,
            .pCleanupFunc = msprime_suite_cleanup,
            .pTests = tests
        },
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
        printf("usage: ./tests <test_name>\n");
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
