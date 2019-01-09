/*
** Copyright (C) 2016-2018 University of Oxford
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

#define _GNU_SOURCE
/*
 * Unit tests for the low-level msprime API.
 */
#include "tsk_genotypes.h"
#include "tsk_convert.h"
#include "tsk_stats.h"

#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <unistd.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
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
const char *single_tree_ex_nodes =/*          6          */
    "1  0   -1   -1\n"             /*         / \         */
    "1  0   -1   -1\n"             /*        /   \        */
    "1  0   -1   -1\n"             /*       /     \       */
    "1  0   -1   -1\n"             /*      /       5      */
    "0  1   -1   -1\n"             /*     4       / \     */
    "0  2   -1   -1\n"             /*    / \     /   \    */
    "0  3   -1   -1\n";            /*   0   1   2     3   */
const char *single_tree_ex_edges =
    "0  1   4   0,1\n"
    "0  1   5   2,3\n"
    "0  1   6   4,5\n";
const char *single_tree_ex_sites =
    "0.1  0\n"
    "0.2  0\n"
    "0.3  0\n";
const char *single_tree_ex_mutations =
    "0    2     1   -1\n"
    "1    4     1   -1\n"
    "1    0     0   1\n"  /* Back mutation over 0 */
    "2    0     1   -1\n"  /* recurrent mutations over samples */
    "2    1     1   -1\n"
    "2    2     1   -1\n"
    "2    3     1   -1\n";

/* Example from the PLOS paper */
const char *paper_ex_nodes =
    "1  0       -1   0\n"
    "1  0       -1   0\n"
    "1  0       -1   1\n"
    "1  0       -1   1\n"
    "0  0.071   -1   -1\n"
    "0  0.090   -1   -1\n"
    "0  0.170   -1   -1\n"
    "0  0.202   -1   -1\n"
    "0  0.253   -1   -1\n";
const char *paper_ex_edges =
    "2 10 4 2\n"
    "2 10 4 3\n"
    "0 10 5 1\n"
    "0 2  5 3\n"
    "2 10 5 4\n"
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
/* Two (diploid) indivduals */
const char *paper_ex_individuals =
    "0      0.2,1.5\n"
    "0      0.0,0.0\n";

/* An example of a nonbinary tree sequence */
const char *nonbinary_ex_nodes =
    "1  0       0   -1\n"
    "1  0       0   -1\n"
    "1  0       0   -1\n"
    "1  0       0   -1\n"
    "1  0       0   -1\n"
    "1  0       0   -1\n"
    "1  0       0   -1\n"
    "1  0       0   -1\n"
    "0  0.01    0   -1\n"
    "0  0.068   0   -1\n"
    "0  0.130   0   -1\n"
    "0  0.279   0   -1\n"
    "0  0.405   0   -1\n";
const char *nonbinary_ex_edges =
    "0	100	8	0,1,2,3\n"
    "0	100	9	6,8\n"
    "0  100 10  4\n"
    "0  17  10  5\n"
    "0  100 10  7\n"
    "17	100	11	5,9\n"
    "0	17	12	9\n"
    "0  100 12  10\n"
    "17	100	12	11";
const char *nonbinary_ex_sites =
        "1  0\n"
        "18 0\n";
const char *nonbinary_ex_mutations =
    "0    2   1\n"
    "1    11  1";

/* An example of a tree sequence with unary nodes. */

const char *unary_ex_nodes =
    "1  0       0  -1\n"
    "1  0       0  -1\n"
    "1  0       0  -1\n"
    "1  0       0  -1\n"
    "0  0.071   0  -1\n"
    "0  0.090   0  -1\n"
    "0  0.170   0  -1\n"
    "0  0.202   0  -1\n"
    "0  0.253   0  -1\n";
const char *unary_ex_edges =
    "2 10 4 2,3\n"
    "0 10 5 1\n"
    "0 2  5 3\n"
    "2 10 5 4\n"
    "0 7  6 0,5\n"
    "7 10 7 0\n"
    "0 2  7 2\n"
    "7 10 7 5\n"
    "0 7  8 6\n"
    "0 2  8 7\n";

/* We make one mutation for each tree, over unary nodes if this exist */
const char *unary_ex_sites =
    "1.0    0\n"
    "4.5    0\n"
    "8.5    0\n";
const char *unary_ex_mutations =
    "0    2   1\n"
    "1    6   1\n"
    "2    5   1\n";

/* An example of a tree sequence with internally sampled nodes. */

/* TODO: find a way to draw these side-by-side */
/*
  7
+-+-+
|   5
| +-++
| |  4
| | +++
| | | 3
| | |
| 1 2
|
0

  8
+-+-+
|   5
| +-++
| |  4
| | +++
3 | | |
  | | |
  1 2 |
      |
      0

  6
+-+-+
|   5
| +-++
| |  4
| | +++
| | | 3
| | |
| 1 2
|
0
*/

const char *internal_sample_ex_nodes =
    "1  0.0   0   -1\n"
    "1  0.1   0   -1\n"
    "1  0.1   0   -1\n"
    "1  0.2   0   -1\n"
    "0  0.4   0   -1\n"
    "1  0.5   0   -1\n"
    "0  0.7   0   -1\n"
    "0  1.0   0   -1\n"
    "0  1.2   0   -1\n";
const char *internal_sample_ex_edges =
    "2 8  4 0\n"
    "0 10 4 2\n"
    "0 2  4 3\n"
    "8 10 4 3\n"
    "0 10 5 1,4\n"
    "8 10 6 0,5\n"
    "0 2  7 0,5\n"
    "2 8  8 3,5\n";
/* We make one mutation for each tree, some above the internal node */
const char *internal_sample_ex_sites =
    "1.0    0\n"
    "4.5    0\n"
    "8.5    0\n";
const char *internal_sample_ex_mutations =
    "0    2   1\n"
    "1    5   1\n"
    "2    5   1\n";


/* Simple utilities to parse text so we can write declaritive
 * tests. This is not intended as a robust general input mechanism.
 */

static void
parse_nodes(const char *text, tsk_node_tbl_t *node_table)
{
    int ret;
    size_t c, k;
    size_t MAX_LINE = 1024;
    char line[MAX_LINE];
    const char *whitespace = " \t";
    char *p;
    double time;
    int flags, population, individual;
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
            individual = -1;
        } else {
            individual = atoi(p);
            p = strtok(NULL, whitespace);
        }
        if (p == NULL) {
            name = "";
        } else {
            name = p;
        }
        ret = tsk_node_tbl_add_row(node_table, flags, time, population,
                individual, name, strlen(name));
        CU_ASSERT_FATAL(ret >= 0);
    }
}

static void
parse_edges(const char *text, tsk_edge_tbl_t *edge_table)
{
    int ret;
    size_t c, k;
    size_t MAX_LINE = 1024;
    char line[MAX_LINE], sub_line[MAX_LINE];
    const char *whitespace = " \t";
    char *p, *q;
    double left, right;
    tsk_id_t parent, child;
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
        strncpy(sub_line, p, MAX_LINE);
        q = strtok(sub_line, ",");
        for (k = 0; k < num_children; k++) {
            CU_ASSERT_FATAL(q != NULL);
            child = atoi(q);
            ret = tsk_edge_tbl_add_row(edge_table, left, right, parent, child);
            CU_ASSERT_FATAL(ret >= 0);
            q = strtok(NULL, ",");
        }
        CU_ASSERT_FATAL(q == NULL);
    }
}

static void
parse_sites(const char *text, tsk_site_tbl_t *site_table)
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
        ret = tsk_site_tbl_add_row(site_table, position, ancestral_state,
                strlen(ancestral_state), NULL, 0);
        CU_ASSERT_FATAL(ret >= 0);
    }
}

static void
parse_mutations(const char *text, tsk_mutation_tbl_t *mutation_table)
{
    int ret;
    size_t c, k;
    size_t MAX_LINE = 1024;
    char line[MAX_LINE];
    const char *whitespace = " \t";
    char *p;
    tsk_id_t node;
    tsk_id_t site;
    tsk_id_t parent;
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
        parent = TSK_NULL;
        p = strtok(NULL, whitespace);
        if (p != NULL) {
            parent = atoi(p);
        }
        ret = tsk_mutation_tbl_add_row(mutation_table, site, node, parent,
                derived_state, strlen(derived_state), NULL, 0);
        CU_ASSERT_FATAL(ret >= 0);
    }
}

static void
parse_individuals(const char *text, tsk_individual_tbl_t *individual_table)
{
    int ret;
    size_t c, k;
    size_t MAX_LINE = 1024;
    char line[MAX_LINE];
    char sub_line[MAX_LINE];
    const char *whitespace = " \t";
    char *p, *q;
    double location[MAX_LINE];
    int location_len;
    int flags;
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
        // the locations are comma-separated
        location_len = 1;
        q = p;
        while (*q != '\0') {
            if (*q == ',') {
                location_len++;
            }
            q++;
        }
        CU_ASSERT_FATAL(location_len >= 1);
        strncpy(sub_line, p, MAX_LINE);
        q = strtok(sub_line, ",");
        for (k = 0; k < location_len; k++) {
            CU_ASSERT_FATAL(q != NULL);
            location[k] = atof(q);
            q = strtok(NULL, ",");
        }
        CU_ASSERT_FATAL(q == NULL);
        p = strtok(NULL, whitespace);
        if (p == NULL) {
            name = "";
        } else {
            name = p;
        }
        ret = tsk_individual_tbl_add_row(individual_table, flags, location, location_len,
                name, strlen(name));
        CU_ASSERT_FATAL(ret >= 0);
    }
}

static void
tsk_treeseq_from_text(tsk_treeseq_t *ts, double sequence_length,
        const char *nodes, const char *edges,
        const char *migrations, const char *sites, const char *mutations,
        const char *individuals, const char *provenance)
{
    int ret;
    tsk_tbl_collection_t tables;
    tsk_id_t max_population_id;
    tsk_tbl_size_t j;

    CU_ASSERT_FATAL(ts != NULL);
    CU_ASSERT_FATAL(nodes != NULL);
    CU_ASSERT_FATAL(edges != NULL);
    /* Not supporting provenance here for now */
    CU_ASSERT_FATAL(provenance == NULL);

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = sequence_length;
    parse_nodes(nodes, tables.nodes);
    parse_edges(edges, tables.edges);
    if (sites != NULL) {
        parse_sites(sites, tables.sites);
    }
    if (mutations != NULL) {
        parse_mutations(mutations, tables.mutations);
    }
    if (individuals != NULL) {
        parse_individuals(individuals, tables.individuals);
    }
    /* We need to add in populations if they are referenced */
    max_population_id = -1;
    for (j = 0; j < tables.nodes->num_rows; j++) {
        max_population_id = TSK_MAX(max_population_id, tables.nodes->population[j]);
    }
    if (max_population_id >= 0) {
        for (j = 0; j <= (tsk_tbl_size_t) max_population_id; j++) {
            ret = tsk_population_tbl_add_row(tables.populations, NULL, 0);
            CU_ASSERT_EQUAL_FATAL(ret, j);
        }
    }

    ret = tsk_treeseq_alloc(ts, &tables, TSK_BUILD_INDEXES);
    /* tsk_treeseq_print_state(ts, stdout); */
    /* printf("ret = %s\n", tsk_strerror(ret)); */
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_tbl_collection_free(&tables);
}

static int
get_max_site_mutations(tsk_treeseq_t *ts)
{
    int ret;
    int max_mutations = 0;
    size_t j;
    tsk_site_t site;

    for (j = 0; j < tsk_treeseq_get_num_sites(ts); j++) {
        ret = tsk_treeseq_get_site(ts, j, &site);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        max_mutations = TSK_MAX(max_mutations, site.mutations_length);
    }
    return max_mutations;
}

static bool
multi_mutations_exist(tsk_treeseq_t *ts, size_t start, size_t end)
{
    int ret;
    size_t j;
    tsk_site_t site;

    for (j = 0; j < TSK_MIN(tsk_treeseq_get_num_sites(ts), end); j++) {
        ret = tsk_treeseq_get_site(ts, j, &site);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        if (site.mutations_length > 1) {
            return true;
        }
    }
    return false;
}

static void
unsort_edges(tsk_edge_tbl_t *edges, size_t start)
{
    size_t j, k;
    size_t n = edges->num_rows - start;
    tsk_edge_t *buff = malloc(n * sizeof(tsk_edge_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(edges != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    gsl_rng_set(rng, 1);

    for (j = 0; j < n; j++) {
        k = start + j;
        buff[j].left = edges->left[k];
        buff[j].right = edges->right[k];
        buff[j].parent = edges->parent[k];
        buff[j].child = edges->child[k];
    }
    gsl_ran_shuffle(rng, buff, n, sizeof(tsk_edge_t));
    for (j = 0; j < n; j++) {
        k = start + j;
        edges->left[k] = buff[j].left;
        edges->right[k] = buff[j].right;
        edges->parent[k] = buff[j].parent;
        edges->child[k] = buff[j].child;
    }
    free(buff);
    gsl_rng_free(rng);
}

static void
unsort_sites(tsk_site_tbl_t *sites, tsk_mutation_tbl_t *mutations)
{
    double position;
    char *ancestral_state = NULL;
    size_t j, k, length;

    if (sites->num_rows > 1) {
        /* Swap the first two sites */
        CU_ASSERT_EQUAL_FATAL(sites->ancestral_state_offset[0], 0);

        position = sites->position[0];
        length = sites->ancestral_state_offset[1];
        /* Save a copy of the first ancestral state */
        ancestral_state = malloc(length);
        CU_ASSERT_FATAL(ancestral_state != NULL);
        memcpy(ancestral_state, sites->ancestral_state, length);
        /* Now write the ancestral state for the site 1 here */
        k = 0;
        for (j = sites->ancestral_state_offset[1]; j < sites->ancestral_state_offset[2];
                j++) {
            sites->ancestral_state[k] = sites->ancestral_state[j];
            k++;
        }
        sites->ancestral_state_offset[1] = k;
        memcpy(sites->ancestral_state + k, ancestral_state, length);
        sites->position[0] = sites->position[1];
        sites->position[1] = position;

        /* Update the mutations for these sites */
        j = 0;
        while (j < mutations->num_rows && mutations->site[j] == 0) {
            mutations->site[j] = 1;
            j++;
        }
        while (j < mutations->num_rows && mutations->site[j] == 1) {
            mutations->site[j] = 0;

            j++;
        }
    }
    tsk_safe_free(ancestral_state);
}

static void
add_individuals(tsk_treeseq_t *ts)
{
    int ret;
    int max_inds = 20;
    tsk_id_t j;
    int k = 0;
    int ploidy = 2;
    tsk_tbl_collection_t tables;
    char *metadata = "abc";
    size_t metadata_length = 3;
    tsk_id_t *samples;
    tsk_tbl_size_t num_samples = tsk_treeseq_get_num_samples(ts);

    ret = tsk_treeseq_get_samples(ts, &samples);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_copy_tables(ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_individual_tbl_clear(tables.individuals);
    memset(tables.nodes->individual, 0xff, tables.nodes->num_rows * sizeof(tsk_id_t));

    k = 0;
    for (j = 0; j < num_samples; j++) {
        if ((k % ploidy) == 0) {
            tsk_individual_tbl_add_row(tables.individuals, (uint32_t) k,
                    NULL, 0, metadata, metadata_length);
            CU_ASSERT_TRUE(ret >= 0)
        }
        tables.nodes->individual[samples[j]] = k / ploidy;
        k += 1;
        if (k >= ploidy * max_inds) {
            break;
        }
    }
    ret = tsk_treeseq_free(ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_alloc(ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_tbl_collection_free(&tables);
}

static void
verify_nodes_equal(tsk_node_t *n1, tsk_node_t *n2)
{
    double eps = 1e-6;

    CU_ASSERT_DOUBLE_EQUAL_FATAL(n1->time, n1->time, eps);
    CU_ASSERT_EQUAL_FATAL(n1->population, n2->population);
    CU_ASSERT_EQUAL_FATAL(n1->flags, n2->flags);
    CU_ASSERT_FATAL(n1->metadata_length == n2->metadata_length);
    CU_ASSERT_NSTRING_EQUAL_FATAL(n1->metadata, n2->metadata, n1->metadata_length);
}

static void
verify_edges_equal(tsk_edge_t *r1, tsk_edge_t *r2, double scale)
{
    double eps = 1e-6;

    CU_ASSERT_DOUBLE_EQUAL_FATAL(r1->left * scale, r2->left, eps);
    CU_ASSERT_DOUBLE_EQUAL_FATAL(r1->right * scale, r2->right, eps);
    CU_ASSERT_EQUAL_FATAL(r1->parent, r2->parent);
    CU_ASSERT_EQUAL_FATAL(r1->child, r2->child);
}

static void
verify_migrations_equal(tsk_migration_t *r1, tsk_migration_t *r2, double scale)
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
verify_provenances_equal(tsk_provenance_t *p1, tsk_provenance_t *p2)
{
    CU_ASSERT_FATAL(p1->timestamp_length == p2->timestamp_length);
    CU_ASSERT_NSTRING_EQUAL_FATAL(p1->timestamp, p2->timestamp, p1->timestamp_length);
    CU_ASSERT_FATAL(p1->record_length == p2->record_length);
    CU_ASSERT_NSTRING_EQUAL_FATAL(p1->record, p2->record, p1->record_length);
}

static void
verify_individuals_equal(tsk_individual_t *i1, tsk_individual_t *i2)
{
    tsk_tbl_size_t j;

    CU_ASSERT_FATAL(i1->id == i2->id);
    CU_ASSERT_FATAL(i1->flags == i2->flags);
    CU_ASSERT_FATAL(i1->metadata_length == i2->metadata_length);
    CU_ASSERT_NSTRING_EQUAL_FATAL(i1->metadata, i2->metadata, i1->metadata_length);
    CU_ASSERT_FATAL(i1->location_length == i2->location_length);
    for (j = 0; j < i1->location_length; j++) {
        CU_ASSERT_EQUAL_FATAL(i1->location[j], i2->location[j]);
    }
}

static void
verify_populations_equal(tsk_population_t *p1, tsk_population_t *p2)
{
    CU_ASSERT_FATAL(p1->id == p2->id);
    CU_ASSERT_FATAL(p1->metadata_length == p2->metadata_length);
    CU_ASSERT_NSTRING_EQUAL_FATAL(p1->metadata, p2->metadata, p1->metadata_length);
}

static tsk_tree_t *
get_tree_list(tsk_treeseq_t *ts)
{
    int ret;
    tsk_tree_t t, *trees;
    size_t num_trees;

    num_trees = tsk_treeseq_get_num_trees(ts);
    ret = tsk_tree_alloc(&t, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    trees = malloc(num_trees * sizeof(tsk_tree_t));
    CU_ASSERT_FATAL(trees != NULL);
    for (ret = tsk_tree_first(&t); ret == 1; ret = tsk_tree_next(&t)) {
        CU_ASSERT_FATAL(t.index < num_trees);
        ret = tsk_tree_alloc(&trees[t.index], ts, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tree_copy(&trees[t.index], &t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tree_equal(&trees[t.index], &t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        /* Make sure the left and right coordinates are also OK */
        CU_ASSERT_DOUBLE_EQUAL(trees[t.index].left, t.left, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(trees[t.index].right, t.right, 1e-6);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_free(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    return trees;
}

static void
verify_tree_next_prev(tsk_treeseq_t *ts)
{
    int ret;
    tsk_tree_t *trees, t;
    size_t j;
    size_t num_trees = tsk_treeseq_get_num_trees(ts);

    trees = get_tree_list(ts);
    ret = tsk_tree_alloc(&t, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Single forward pass */
    j = 0;
    for (ret = tsk_tree_first(&t); ret == 1; ret = tsk_tree_next(&t)) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);

    /* Single reverse pass */
    j = num_trees;
    for (ret = tsk_tree_last(&t); ret == 1; ret = tsk_tree_prev(&t)) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        if (ret != 0) {
            printf("trees differ\n");
            printf("REVERSE tree::\n");
            tsk_tree_print_state(&t, stdout);
            printf("FORWARD tree::\n");
            tsk_tree_print_state(&trees[t.index], stdout);
        }
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);

    /* Full forward, then reverse */
    j = 0;
    for (ret = tsk_tree_first(&t); ret == 1; ret = tsk_tree_next(&t)) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);
    j--;
    while ((ret = tsk_tree_prev(&t)) == 1) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);
    CU_ASSERT_EQUAL_FATAL(t.index, 0);
    /* Calling prev should return 0 and have no effect. */
    for (j = 0; j < 10; j++) {
        ret = tsk_tree_prev(&t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(t.index, 0);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    /* Full reverse then forward */
    j = num_trees;
    for (ret = tsk_tree_last(&t); ret == 1; ret = tsk_tree_prev(&t)) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);
    j++;
    while ((ret = tsk_tree_next(&t)) == 1) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);
    CU_ASSERT_EQUAL_FATAL(t.index, num_trees - 1);
    /* Calling next should return 0 and have no effect. */
    for (j = 0; j < 10; j++) {
        ret = tsk_tree_next(&t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(t.index, num_trees - 1);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    /* Do a zigzagging traversal */
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 1; j < TSK_MIN(10, num_trees / 2); j++) {
        while (t.index < num_trees - j) {
            ret = tsk_tree_next(&t);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        CU_ASSERT_EQUAL_FATAL(t.index, num_trees - j);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        while (t.index > j) {
            ret = tsk_tree_prev(&t);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        CU_ASSERT_EQUAL_FATAL(t.index, j);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    /* Free the trees. */
    ret = tsk_tree_free(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < tsk_treeseq_get_num_trees(ts); j++) {
        ret = tsk_tree_free(&trees[j]);
    }
    free(trees);
}

static void
verify_hapgen(tsk_treeseq_t *ts)
{
    int ret;
    tsk_hapgen_t hapgen;
    char *haplotype;
    size_t num_samples = tsk_treeseq_get_num_samples(ts);
    size_t num_sites = tsk_treeseq_get_num_sites(ts);
    tsk_site_t site;
    size_t j;
    int k;
    bool single_char = true;

    for (j = 0; j < num_sites; j++) {
        ret = tsk_treeseq_get_site(ts, j, &site);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        if (site.ancestral_state_length != 1) {
            single_char = false;
        }
        for (k = 0; k < site.mutations_length; k++) {
            if (site.mutations[k].derived_state_length != 1) {
                single_char = false;
            }
        }
    }

    ret = tsk_hapgen_alloc(&hapgen, ts);
    if (single_char) {
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        tsk_hapgen_print_state(&hapgen, _devnull);

        for (j = 0; j < num_samples; j++) {
            ret = tsk_hapgen_get_haplotype(&hapgen, j, &haplotype);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(strlen(haplotype), num_sites);
        }
        for (j = num_samples; j < num_samples + 10; j++) {
            ret = tsk_hapgen_get_haplotype(&hapgen, j, &haplotype);
            CU_ASSERT_EQUAL(ret, TSK_ERR_OUT_OF_BOUNDS);
        }
    } else {
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NON_SINGLE_CHAR_MUTATION);
    }
    ret = tsk_hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
}

static void
verify_vargen(tsk_treeseq_t *ts)
{
    int ret;
    tsk_vargen_t vargen;
    size_t num_samples = tsk_treeseq_get_num_samples(ts);
    size_t num_sites = tsk_treeseq_get_num_sites(ts);
    tsk_variant_t *var;
    size_t j, k, f, s;
    int flags[] = {0, TSK_16_BIT_GENOTYPES};
    tsk_id_t *samples[] = {NULL, NULL};

    ret = tsk_treeseq_get_samples(ts, samples);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (s = 0; s < 2; s++) {
        for (f = 0; f < sizeof(flags) / sizeof(*flags); f++) {
            ret = tsk_vargen_alloc(&vargen, ts, samples[s], num_samples, flags[f]);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            tsk_vargen_print_state(&vargen, _devnull);
            j = 0;
            while ((ret = tsk_vargen_next(&vargen, &var)) == 1) {
                CU_ASSERT_EQUAL(var->site->id, j);
                if (var->site->mutations_length == 0) {
                    CU_ASSERT_EQUAL(var->num_alleles, 1);
                } else {
                    CU_ASSERT_TRUE(var->num_alleles > 1);
                }
                CU_ASSERT_EQUAL(var->allele_lengths[0], var->site->ancestral_state_length);
                CU_ASSERT_NSTRING_EQUAL_FATAL(var->alleles[0], var->site->ancestral_state,
                        var->allele_lengths[0]);
                for (k = 0; k < var->num_alleles; k++) {
                    CU_ASSERT_TRUE(var->allele_lengths[k] >= 0);
                }
                for (k = 0; k < num_samples; k++) {
                    if (flags[f] == TSK_16_BIT_GENOTYPES) {
                        CU_ASSERT(var->genotypes.u16[k] <= var->num_alleles);
                    } else  {
                        CU_ASSERT(var->genotypes.u8[k] <= var->num_alleles);
                    }
                }
                j++;
            }
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_EQUAL(j, num_sites);
            CU_ASSERT_EQUAL_FATAL(tsk_vargen_next(&vargen, &var), 0);
            ret = tsk_vargen_free(&vargen);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
        }
    }
}

static void
verify_stats(tsk_treeseq_t *ts)
{
    int ret;
    uint32_t num_samples = tsk_treeseq_get_num_samples(ts);
    tsk_id_t *samples;
    uint32_t j;
    double pi;
    int max_site_mutations = get_max_site_mutations(ts);

    ret = tsk_treeseq_get_pairwise_diversity(ts, NULL, 0, &pi);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_treeseq_get_pairwise_diversity(ts, NULL, 1, &pi);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_treeseq_get_pairwise_diversity(ts, NULL, num_samples + 1, &pi);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    ret = tsk_treeseq_get_samples(ts, &samples);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 2; j < num_samples; j++) {
        ret = tsk_treeseq_get_pairwise_diversity(ts, samples, j, &pi);
        if (max_site_mutations <= 1) {
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_TRUE_FATAL(pi >= 0);
        } else {
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_ONLY_INFINITE_SITES);
        }
    }
}

/* FIXME: this test is weak and should check the return value somehow.
 * We should also have simplest and single tree tests along with separate
 * tests for the error conditions. This should be done as part of the general
 * stats framework.
 */
static void
verify_genealogical_nearest_neighbours(tsk_treeseq_t *ts)
{
    int ret;
    tsk_id_t *samples;
    tsk_id_t *sample_sets[2];
    size_t sample_set_size[2];
    size_t num_samples = tsk_treeseq_get_num_samples(ts);
    double *A = malloc(2 * num_samples * sizeof(double));
    CU_ASSERT_FATAL(A != NULL);

    ret = tsk_treeseq_get_samples(ts, &samples);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    sample_sets[0] = samples;
    sample_set_size[0] = num_samples / 2;
    sample_sets[1] = samples + sample_set_size[0];
    sample_set_size[1] = num_samples - sample_set_size[0];

    ret = tsk_treeseq_genealogical_nearest_neighbours(ts,
        samples, num_samples, sample_sets, sample_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    free(A);
}

/* FIXME: this test is weak and should check the return value somehow.
 * We should also have simplest and single tree tests along with separate
 * tests for the error conditions. This should be done as part of the general
 * stats framework.
 */
static void
verify_mean_descendants(tsk_treeseq_t *ts)
{
    int ret;
    tsk_id_t *samples;
    tsk_id_t *sample_sets[2];
    size_t sample_set_size[2];
    size_t num_samples = tsk_treeseq_get_num_samples(ts);
    double *C = malloc(2 * tsk_treeseq_get_num_nodes(ts) * sizeof(double));
    CU_ASSERT_FATAL(C != NULL);

    ret = tsk_treeseq_get_samples(ts, &samples);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    sample_sets[0] = samples;
    sample_set_size[0] = num_samples / 2;
    sample_sets[1] = samples + sample_set_size[0];
    sample_set_size[1] = num_samples - sample_set_size[0];

    ret = tsk_treeseq_mean_descendants(ts, sample_sets, sample_set_size, 2, 0, C);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Check some error conditions */
    ret = tsk_treeseq_mean_descendants(ts, sample_sets, sample_set_size, 0, 0, C);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    samples[0] = -1;
    ret = tsk_treeseq_mean_descendants(ts, sample_sets, sample_set_size, 2, 0, C);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_OUT_OF_BOUNDS);
    samples[0] = tsk_treeseq_get_num_nodes(ts) + 1;
    ret = tsk_treeseq_mean_descendants(ts, sample_sets, sample_set_size, 2, 0, C);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_OUT_OF_BOUNDS);

    free(C);
}


static void
verify_compute_mutation_parents(tsk_treeseq_t *ts)
{
    int ret;
    size_t size = tsk_treeseq_get_num_mutations(ts) * sizeof(tsk_id_t);
    tsk_id_t *parent = malloc(size);
    tsk_tbl_collection_t tables;

    CU_ASSERT_FATAL(parent != NULL);
    ret = tsk_treeseq_copy_tables(ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    memcpy(parent, tables.mutations->parent, size);
    /* tsk_tbl_collection_print_state(&tables, stdout); */
    /* Make sure the tables are actually updated */
    memset(tables.mutations->parent, 0xff, size);

    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(memcmp(parent, tables.mutations->parent, size), 0);
    /* printf("after\n"); */
    /* tsk_tbl_collection_print_state(&tables, stdout); */

    free(parent);
    tsk_tbl_collection_free(&tables);
}

static void
verify_individual_nodes(tsk_treeseq_t *ts)
{
    int ret;
    tsk_individual_t individual;
    tsk_id_t k;
    size_t num_nodes = tsk_treeseq_get_num_nodes(ts);
    size_t num_individuals = tsk_treeseq_get_num_individuals(ts);
    size_t j;

    for (k = 0; k < (tsk_id_t) num_individuals; k++) {
        ret = tsk_treeseq_get_individual(ts, k, &individual);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_FATAL(individual.nodes_length >= 0);
        for (j = 0; j < individual.nodes_length; j++) {
            CU_ASSERT_FATAL(individual.nodes[j] < num_nodes);
            CU_ASSERT_EQUAL_FATAL(k,
                    ts->tables->nodes->individual[individual.nodes[j]]);
        }
    }
}

/* When we keep all sites in simplify, the genotypes for the subset of the
 * samples should be the same as the original */
static void
verify_simplify_genotypes(tsk_treeseq_t *ts, tsk_treeseq_t *subset,
        tsk_id_t *samples, uint32_t num_samples, tsk_id_t *node_map)
{
    int ret;
    size_t m = tsk_treeseq_get_num_sites(ts);
    tsk_vargen_t vargen, subset_vargen;
    tsk_variant_t *variant, *subset_variant;
    size_t j, k;
    tsk_id_t *all_samples;
    uint8_t a1, a2;
    tsk_id_t *sample_index_map;

    tsk_treeseq_get_sample_index_map(ts, &sample_index_map);

    /* tsk_treeseq_print_state(ts, stdout); */
    /* tsk_treeseq_print_state(subset, stdout); */

    ret = tsk_vargen_alloc(&vargen, ts, NULL, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_vargen_alloc(&subset_vargen, subset, NULL, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(m, tsk_treeseq_get_num_sites(subset));
    tsk_treeseq_get_samples(ts, &all_samples);

    for (j = 0; j < m; j++) {
        ret = tsk_vargen_next(&vargen, &variant);
        CU_ASSERT_EQUAL_FATAL(ret, 1);
        ret = tsk_vargen_next(&subset_vargen, &subset_variant);
        CU_ASSERT_EQUAL_FATAL(ret, 1);
        CU_ASSERT_EQUAL(variant->site->id, j)
        CU_ASSERT_EQUAL(subset_variant->site->id, j)
        CU_ASSERT_EQUAL(variant->site->position, subset_variant->site->position);
        for (k = 0; k < num_samples; k++) {
            CU_ASSERT_FATAL(sample_index_map[samples[k]] < ts->num_samples);
            a1 = variant->genotypes.u8[sample_index_map[samples[k]]];
            a2 = subset_variant->genotypes.u8[k];
            /* printf("a1 = %d, a2 = %d\n", a1, a2); */
            /* printf("k = %d original node = %d " */
            /*         "original_index = %d a1=%.*s a2=%.*s\n", */
            /*         (int) k, samples[k], sample_index_map[samples[k]], */
            /*         variant->allele_lengths[a1], variant->alleles[a1], */
            /*         subset_variant->allele_lengths[a2], subset_variant->alleles[a2]); */
            CU_ASSERT_FATAL(a1 < variant->num_alleles);
            CU_ASSERT_FATAL(a2 < subset_variant->num_alleles);
            CU_ASSERT_EQUAL_FATAL(variant->allele_lengths[a1],
                    subset_variant->allele_lengths[a2]);
            CU_ASSERT_NSTRING_EQUAL_FATAL(
                variant->alleles[a1], subset_variant->alleles[a2],
                variant->allele_lengths[a1]);
        }
    }
    tsk_vargen_free(&vargen);
    tsk_vargen_free(&subset_vargen);
}


static void
verify_simplify_properties(tsk_treeseq_t *ts, tsk_treeseq_t *subset,
        tsk_id_t *samples, uint32_t num_samples, tsk_id_t *node_map)
{
    int ret;
    tsk_node_t n1, n2;
    tsk_tree_t full_tree, subset_tree;
    tsk_site_t *tree_sites;
    tsk_tbl_size_t tree_sites_length;
    uint32_t j, k;
    tsk_id_t u, mrca1, mrca2;
    size_t total_sites;

    CU_ASSERT_EQUAL(
        tsk_treeseq_get_sequence_length(ts),
        tsk_treeseq_get_sequence_length(subset));
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(subset), num_samples);
    CU_ASSERT(
        tsk_treeseq_get_num_nodes(ts) >= tsk_treeseq_get_num_nodes(subset));
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(subset), num_samples);

    /* Check the sample properties */
    for (j = 0; j < num_samples; j++) {
        ret = tsk_treeseq_get_node(ts, samples[j], &n1);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(node_map[samples[j]], j);
        ret = tsk_treeseq_get_node(subset, node_map[samples[j]], &n2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(n1.population, n2.population);
        CU_ASSERT_EQUAL_FATAL(n1.time, n2.time);
        CU_ASSERT_EQUAL_FATAL(n1.flags, n2.flags);
        CU_ASSERT_EQUAL_FATAL(n1.metadata_length, n2.metadata_length);
        CU_ASSERT_NSTRING_EQUAL(n1.metadata, n2.metadata, n2.metadata_length);
    }
    /* Check that node mappings are correct */
    for (j = 0; j < tsk_treeseq_get_num_nodes(ts); j++) {
        ret = tsk_treeseq_get_node(ts, j, &n1);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        if (node_map[j] != TSK_NULL) {
            ret = tsk_treeseq_get_node(subset, node_map[j], &n2);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(n1.population, n2.population);
            CU_ASSERT_EQUAL_FATAL(n1.time, n2.time);
            CU_ASSERT_EQUAL_FATAL(n1.flags, n2.flags);
            CU_ASSERT_EQUAL_FATAL(n1.metadata_length, n2.metadata_length);
            CU_ASSERT_NSTRING_EQUAL(n1.metadata, n2.metadata, n2.metadata_length);
        }
    }
    if (num_samples == 0) {
        CU_ASSERT_EQUAL(tsk_treeseq_get_num_edges(subset), 0);
        CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(subset), 0);
    } else if (num_samples == 1) {
        CU_ASSERT_EQUAL(tsk_treeseq_get_num_edges(subset), 0);
        CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(subset), 1);
    }
    /* Check the pairwise MRCAs */
    ret = tsk_tree_alloc(&full_tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_alloc(&subset_tree, subset, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&full_tree);
    CU_ASSERT_EQUAL(ret, 1);
    ret = tsk_tree_first(&subset_tree);
    CU_ASSERT_EQUAL(ret, 1);

    total_sites = 0;
    while (1) {
        while (full_tree.right <= subset_tree.right) {
            for (j = 0; j < num_samples; j++) {
                for (k = j + 1; k < num_samples; k++) {
                    ret = tsk_tree_get_mrca(&full_tree, samples[j], samples[k], &mrca1);
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                    ret = tsk_tree_get_mrca(&subset_tree,
                            node_map[samples[j]], node_map[samples[k]], &mrca2);
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                    if (mrca1 == TSK_NULL) {
                        CU_ASSERT_EQUAL_FATAL(mrca2, TSK_NULL);
                    } else {
                        CU_ASSERT_EQUAL(node_map[mrca1], mrca2);
                    }
                }
            }
            ret = tsk_tree_next(&full_tree);
            CU_ASSERT_FATAL(ret >= 0);
            if (ret != 1) {
                break;
            }
        }
        /* Check the sites in this tree */
        ret = tsk_tree_get_sites(&subset_tree, &tree_sites, &tree_sites_length);
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        for (j = 0; j < tree_sites_length; j++) {
            CU_ASSERT(subset_tree.left <= tree_sites[j].position);
            CU_ASSERT(tree_sites[j].position < subset_tree.right);
            for (k = 0; k < tree_sites[j].mutations_length; k++) {
                ret = tsk_tree_get_parent(&subset_tree,
                        tree_sites[j].mutations[k].node, &u);
                CU_ASSERT_EQUAL(ret, 0);
            }
            total_sites++;
        }
        ret = tsk_tree_next(&subset_tree);
        if (ret != 1) {
            break;
        }
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(subset), total_sites);

    tsk_tree_free(&subset_tree);
    tsk_tree_free(&full_tree);
    verify_vargen(subset);
    verify_hapgen(subset);
}

static void
verify_simplify(tsk_treeseq_t *ts)
{
    int ret;
    uint32_t n = tsk_treeseq_get_num_samples(ts);
    uint32_t num_samples[] = {0, 1, 2, 3, n / 2, n - 1, n};
    size_t j;
    tsk_id_t *sample;
    tsk_id_t *node_map = malloc(tsk_treeseq_get_num_nodes(ts) * sizeof(tsk_id_t));
    tsk_treeseq_t subset;
    int flags = TSK_FILTER_SITES;

    CU_ASSERT_FATAL(node_map != NULL);
    ret = tsk_treeseq_get_samples(ts, &sample);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    if (tsk_treeseq_get_num_migrations(ts) > 0) {
        ret = tsk_treeseq_simplify(ts, sample, 2, 0, &subset, NULL);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SIMPLIFY_MIGRATIONS_NOT_SUPPORTED);
        /* Exiting early here because simplify isn't supported with migrations. */
        goto out;
    }

    for (j = 0; j < sizeof(num_samples) / sizeof(uint32_t); j++) {
        if (num_samples[j] <= n) {
            ret = tsk_treeseq_simplify(ts, sample, num_samples[j], flags, &subset,
                    node_map);
            /* printf("ret = %s\n", tsk_strerror(ret)); */
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            verify_simplify_properties(ts, &subset, sample, num_samples[j], node_map);
            tsk_treeseq_free(&subset);

            /* Keep all sites */
            ret = tsk_treeseq_simplify(ts, sample, num_samples[j], 0, &subset,
                    node_map);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            verify_simplify_properties(ts, &subset, sample, num_samples[j], node_map);
            verify_simplify_genotypes(ts, &subset, sample, num_samples[j], node_map);
            tsk_treeseq_free(&subset);
        }
    }
out:
    free(node_map);
}

static void
verify_reduce_topology(tsk_treeseq_t *ts)
{
    int ret;
    size_t j;
    tsk_id_t *sample;
    tsk_treeseq_t reduced;
    tsk_edge_t edge;
    double *X;
    size_t num_sites;
    size_t n = tsk_treeseq_get_num_samples(ts);
    int flags = TSK_REDUCE_TO_SITE_TOPOLOGY;

    ret = tsk_treeseq_get_samples(ts, &sample);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    if (tsk_treeseq_get_num_migrations(ts) > 0) {
        ret = tsk_treeseq_simplify(ts, sample, 2, flags, &reduced, NULL);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SIMPLIFY_MIGRATIONS_NOT_SUPPORTED);
        return;
    }

    ret = tsk_treeseq_simplify(ts, sample, n, flags, &reduced, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    X = reduced.tables->sites->position;
    num_sites = reduced.tables->sites->num_rows;
    if (num_sites == 0) {
        CU_ASSERT_EQUAL_FATAL(tsk_treeseq_get_num_edges(&reduced), 0);
    }
    for (j = 0; j < tsk_treeseq_get_num_edges(&reduced); j++) {
        ret = tsk_treeseq_get_edge(&reduced, j, &edge);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        if (edge.left != 0) {
            CU_ASSERT_EQUAL_FATAL(edge.left,
                X[tsk_search_sorted(X, num_sites, edge.left)]);
        }
        if (edge.right != tsk_treeseq_get_sequence_length(&reduced)) {
            CU_ASSERT_EQUAL_FATAL(edge.right,
                X[tsk_search_sorted(X, num_sites, edge.right)]);
        }
    }
    tsk_treeseq_free(&reduced);
}

/* Utility function to return a tree sequence for testing. It is the
 * callers responsilibility to free all memory.
 */
static tsk_treeseq_t *
get_example_tree_sequence(uint32_t num_samples,
        uint32_t num_historical_samples, uint32_t num_loci,
        double sequence_length, double recombination_rate,
        double mutation_rate, uint32_t num_bottlenecks,
        bottleneck_desc_t *bottlenecks, int alphabet)
{
    return NULL;
}

tsk_treeseq_t **
get_example_nonbinary_tree_sequences(void)
{
    return NULL;
}

tsk_treeseq_t *
make_recurrent_and_back_mutations_copy(tsk_treeseq_t *ts)
{
    return NULL;
}

tsk_treeseq_t *
make_permuted_nodes_copy(tsk_treeseq_t *ts)
{
    return NULL;
}

/* Insert some gaps into the specified tree sequence, i.e., positions
 * that no edge covers. */
tsk_treeseq_t *
make_gappy_copy(tsk_treeseq_t *ts)
{
    return NULL;
}

/* Return a copy of the tree sequence after deleting half of its edges.
 */
tsk_treeseq_t *
make_decapitated_copy(tsk_treeseq_t *ts)
{
    return NULL;
}

tsk_treeseq_t *
make_multichar_mutations_copy(tsk_treeseq_t *ts)
{
    return NULL;
}

tsk_treeseq_t **
get_example_tree_sequences(int include_nonbinary)
{
    size_t max_examples = 1024;
    tsk_treeseq_t **ret = malloc(max_examples * sizeof(tsk_treeseq_t *));
    ret[0] = NULL;
    return ret;
}

static void
verify_vcf_converter(tsk_treeseq_t *ts, unsigned int ploidy)
{
    int ret;
    char *str = NULL;
    tsk_vcf_converter_t vc;
    unsigned int num_variants;

    ret = tsk_vcf_converter_alloc(&vc, ts, ploidy, "chr1234");
    CU_ASSERT_FATAL(ret ==  0);
    tsk_vcf_converter_print_state(&vc, _devnull);
    ret = tsk_vcf_converter_get_header(&vc, &str);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_NSTRING_EQUAL("##", str, 2);
    num_variants = 0;
    while ((ret = tsk_vcf_converter_next(&vc, &str)) == 1) {
        CU_ASSERT_NSTRING_EQUAL("chr1234\t", str, 2);
        num_variants++;
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(num_variants == tsk_treeseq_get_num_mutations(ts));
    tsk_vcf_converter_free(&vc);
}

static void
test_vcf(void)
{
    int ret;
    unsigned int ploidy;
    tsk_vcf_converter_t *vc = malloc(sizeof(tsk_vcf_converter_t));
    tsk_treeseq_t *ts = get_example_tree_sequence(10, 0, 100, 100.0, 1.0, 1.0,
            0, NULL, 0);

    CU_ASSERT_FATAL(ts != NULL);
    CU_ASSERT_FATAL(vc != NULL);

    ret = tsk_vcf_converter_alloc(vc, ts, 0, "1");
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_vcf_converter_alloc(vc, ts, 3, "1");
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_vcf_converter_alloc(vc, ts, 11, "1");
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    for (ploidy = 1; ploidy < 3; ploidy++) {
        verify_vcf_converter(ts, ploidy);
    }

    free(vc);
    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_vcf_no_mutations(void)
{
    int ret;
    char *str = NULL;
    tsk_vcf_converter_t *vc = malloc(sizeof(tsk_vcf_converter_t));
    tsk_treeseq_t *ts = get_example_tree_sequence(100, 0, 1, 1.0, 0.0, 0.0, 0, NULL,
            0);

    CU_ASSERT_FATAL(ts != NULL);
    CU_ASSERT_FATAL(vc != NULL);
    CU_ASSERT_EQUAL_FATAL(tsk_treeseq_get_num_mutations(ts), 0);

    ret = tsk_vcf_converter_alloc(vc, ts, 1, "1");
    CU_ASSERT_FATAL(ret ==  0);
    tsk_vcf_converter_print_state(vc, _devnull);
    ret = tsk_vcf_converter_get_header(vc, &str);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_NSTRING_EQUAL("##", str, 2);
    ret = tsk_vcf_converter_next(vc, &str);
    CU_ASSERT_EQUAL(ret, 0);
    tsk_vcf_converter_free(vc);

    free(vc);
    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_node_metadata(void)
{
    const char *nodes =
        "1  0   0   -1   n1\n"
        "1  0   0   -1   n2\n"
        "0  1   0   -1   A_much_longer_name\n"
        "0  1   0   -1\n"
        "0  1   0   -1   n4";
    const char *edges =
        "0  1   2   0,1\n";
    tsk_treeseq_t ts;
    int ret;
    tsk_node_t node;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 5);

    ret = tsk_treeseq_get_node(&ts, 0, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_NSTRING_EQUAL(node.metadata, "n1", 2);

    ret = tsk_treeseq_get_node(&ts, 1, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_NSTRING_EQUAL(node.metadata, "n2", 2);

    ret = tsk_treeseq_get_node(&ts, 2, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_NSTRING_EQUAL(node.metadata, "A_much_longer_name", 18);

    ret = tsk_treeseq_get_node(&ts, 3, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_NSTRING_EQUAL(node.metadata, "", 0);

    ret = tsk_treeseq_get_node(&ts, 4, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_NSTRING_EQUAL(node.metadata, "n4", 2);

    tsk_treeseq_free(&ts);
}

static void
verify_trees_consistent(tsk_treeseq_t *ts)
{
    int ret;
    size_t num_trees;
    tsk_tree_t tree;

    ret = tsk_tree_alloc(&tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    num_trees = 0;
    for (ret = tsk_tree_first(&tree); ret == 1; ret = tsk_tree_next(&tree)) {
        tsk_tree_print_state(&tree, _devnull);
        CU_ASSERT_EQUAL(tree.index, num_trees);
        num_trees++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(ts), num_trees);

    tsk_tree_free(&tree);
}

static void
verify_ld(tsk_treeseq_t *ts)
{
    int ret;
    size_t num_sites = tsk_treeseq_get_num_sites(ts);
    tsk_site_t *sites = malloc(num_sites * sizeof(tsk_site_t));
    int *num_site_mutations = malloc(num_sites * sizeof(int));
    tsk_ld_calc_t ld_calc;
    double *r2, *r2_prime, x;
    size_t j, num_r2_values;
    double eps = 1e-6;

    r2 = calloc(num_sites, sizeof(double));
    r2_prime = calloc(num_sites, sizeof(double));
    CU_ASSERT_FATAL(r2 != NULL);
    CU_ASSERT_FATAL(r2_prime != NULL);
    CU_ASSERT_FATAL(sites != NULL);
    CU_ASSERT_FATAL(num_site_mutations != NULL);

    ret = tsk_ld_calc_alloc(&ld_calc, ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_ld_calc_print_state(&ld_calc, _devnull);

    for (j = 0; j < num_sites; j++) {
        ret = tsk_treeseq_get_site(ts, j, sites + j);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        num_site_mutations[j] = sites[j].mutations_length;
        ret = tsk_ld_calc_get_r2(&ld_calc, j, j, &x);
        if (num_site_mutations[j] <= 1) {
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_DOUBLE_EQUAL_FATAL(x, 1.0, eps);
        } else {
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_ONLY_INFINITE_SITES);
        }
    }

    if (num_sites > 0) {
        /* Some checks in the forward direction */
        ret = tsk_ld_calc_get_r2_array(&ld_calc, 0, TSK_DIR_FORWARD,
                num_sites, DBL_MAX, r2, &num_r2_values);
        if (multi_mutations_exist(ts, 0, num_sites)) {
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_ONLY_INFINITE_SITES);
        } else {
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(num_r2_values, num_sites - 1);
        }
        tsk_ld_calc_print_state(&ld_calc, _devnull);

        ret = tsk_ld_calc_get_r2_array(&ld_calc, num_sites - 2, TSK_DIR_FORWARD,
                num_sites, DBL_MAX, r2_prime, &num_r2_values);
        if (multi_mutations_exist(ts, num_sites - 2, num_sites)) {
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_ONLY_INFINITE_SITES);
        } else {
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(num_r2_values, 1);
        }
        tsk_ld_calc_print_state(&ld_calc, _devnull);

        ret = tsk_ld_calc_get_r2_array(&ld_calc, 0, TSK_DIR_FORWARD,
                num_sites, DBL_MAX, r2_prime, &num_r2_values);
        if (multi_mutations_exist(ts, 0, num_sites)) {
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_ONLY_INFINITE_SITES);
        } else {
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(num_r2_values, num_sites - 1);
            tsk_ld_calc_print_state(&ld_calc, _devnull);
            for (j = 0; j < num_r2_values; j++) {
                CU_ASSERT_EQUAL_FATAL(r2[j], r2_prime[j]);
                ret = tsk_ld_calc_get_r2(&ld_calc, 0, j + 1, &x);
                CU_ASSERT_EQUAL_FATAL(ret, 0);
                CU_ASSERT_DOUBLE_EQUAL_FATAL(r2[j], x, eps);
            }

        }

        /* Some checks in the reverse direction */
        ret = tsk_ld_calc_get_r2_array(&ld_calc, num_sites - 1,
                TSK_DIR_REVERSE, num_sites, DBL_MAX,
                r2, &num_r2_values);
        if (multi_mutations_exist(ts, 0, num_sites)) {
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_ONLY_INFINITE_SITES);
        } else {
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(num_r2_values, num_sites - 1);
        }
        tsk_ld_calc_print_state(&ld_calc, _devnull);

        ret = tsk_ld_calc_get_r2_array(&ld_calc, 1, TSK_DIR_REVERSE,
                num_sites, DBL_MAX, r2_prime, &num_r2_values);
        if (multi_mutations_exist(ts, 0, 1)) {
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_ONLY_INFINITE_SITES);
        } else {
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(num_r2_values, 1);
        }
        tsk_ld_calc_print_state(&ld_calc, _devnull);

        ret = tsk_ld_calc_get_r2_array(&ld_calc, num_sites - 1,
                TSK_DIR_REVERSE, num_sites, DBL_MAX,
                r2_prime, &num_r2_values);
        if (multi_mutations_exist(ts, 0, num_sites)) {
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_ONLY_INFINITE_SITES);
        } else {
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(num_r2_values, num_sites - 1);
            tsk_ld_calc_print_state(&ld_calc, _devnull);

            for (j = 0; j < num_r2_values; j++) {
                CU_ASSERT_EQUAL_FATAL(r2[j], r2_prime[j]);
                ret = tsk_ld_calc_get_r2(&ld_calc, num_sites - 1,
                        num_sites - j - 2, &x);
                CU_ASSERT_EQUAL_FATAL(ret, 0);
                CU_ASSERT_DOUBLE_EQUAL_FATAL(r2[j], x, eps);
            }
        }

        /* Check some error conditions */
        ret = tsk_ld_calc_get_r2_array(&ld_calc, 0, 0, num_sites, DBL_MAX,
            r2, &num_r2_values);
        CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    }

    if (num_sites > 3) {
        /* Check for some basic distance calculations */
        j = num_sites / 2;
        x = sites[j + 1].position - sites[j].position;
        ret = tsk_ld_calc_get_r2_array(&ld_calc, j, TSK_DIR_FORWARD, num_sites,
                x, r2, &num_r2_values);
        if (multi_mutations_exist(ts, j, num_sites)) {
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_ONLY_INFINITE_SITES);
        } else {
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(num_r2_values, 1);
        }

        x = sites[j].position - sites[j - 1].position;
        ret = tsk_ld_calc_get_r2_array(&ld_calc, j, TSK_DIR_REVERSE, num_sites,
                x, r2, &num_r2_values);
        if (multi_mutations_exist(ts, 0, j + 1)) {
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_ONLY_INFINITE_SITES);
        } else {
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(num_r2_values, 1);
        }
    }

    /* Check some error conditions */
    for (j = num_sites; j < num_sites + 2; j++) {
        ret = tsk_ld_calc_get_r2_array(&ld_calc, j, TSK_DIR_FORWARD,
                num_sites, DBL_MAX, r2, &num_r2_values);
        CU_ASSERT_EQUAL(ret, TSK_ERR_OUT_OF_BOUNDS);
        ret = tsk_ld_calc_get_r2(&ld_calc, j, 0, r2);
        CU_ASSERT_EQUAL(ret, TSK_ERR_OUT_OF_BOUNDS);
        ret = tsk_ld_calc_get_r2(&ld_calc, 0, j, r2);
        CU_ASSERT_EQUAL(ret, TSK_ERR_OUT_OF_BOUNDS);
    }

    tsk_ld_calc_free(&ld_calc);
    free(r2);
    free(r2_prime);
    free(sites);
    free(num_site_mutations);
}

static void
verify_empty_tree_sequence(tsk_treeseq_t *ts, double sequence_length)
{
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_edges(ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_migrations(ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(ts), sequence_length);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(ts), 1);
    verify_trees_consistent(ts);
    verify_ld(ts);
    verify_stats(ts);
    verify_hapgen(ts);
    verify_vargen(ts);
    verify_vcf_converter(ts, 1);
}
static void
test_single_tree_newick(void)
{
    /* int ret; */
    /* tsk_treeseq_t ts; */
    /* tsk_tree_t t; */
    /* size_t buffer_size = 1024; */
    /* char newick[buffer_size]; */

    /* tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, */
    /*         NULL, NULL, NULL, NULL, NULL); */
    /* CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 4); */
    /* CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1); */

    /* ret = tsk_tree_alloc(&t, &ts, 0); */
    /* CU_ASSERT_EQUAL_FATAL(ret, 0) */
    /* ret = tsk_tree_first(&t); */
    /* CU_ASSERT_EQUAL_FATAL(ret, 1) */


    /* ret = tsk_tree_get_newick(&t, -1, 1, 0, buffer_size, newick); */
    /* CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS); */
    /* ret = tsk_tree_get_newick(&t, 7, 1, 0, buffer_size, newick); */
    /* CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS); */

    /* ret = tsk_tree_get_newick(&t, 0, 0, 0, buffer_size, newick); */
    /* CU_ASSERT_EQUAL_FATAL(ret, 0); */
    /* /1* Seems odd, but this is what a single node newick tree looks like. */
    /*  * Newick parsers seems to accept it in any case *1/ */
    /* CU_ASSERT_STRING_EQUAL(newick, "1;"); */

    /* ret = tsk_tree_get_newick(&t, 4, 0, 0, buffer_size, newick); */
    /* CU_ASSERT_EQUAL_FATAL(ret, 0); */
    /* CU_ASSERT_STRING_EQUAL(newick, "(1:1,2:1);"); */

    /* ret = tsk_tree_get_newick(&t, 6, 0, 0, buffer_size, newick); */
    /* CU_ASSERT_EQUAL_FATAL(ret, 0); */
    /* CU_ASSERT_STRING_EQUAL(newick, "((1:1,2:1):2,(3:2,4:2):1);"); */

    /* tsk_tree_free(&t); */
    /* tsk_treeseq_free(&ts); */
}


static void
verify_sample_sets_for_tree(tsk_tree_t *tree)
{
    int ret, stack_top, j;
    tsk_id_t u, v, n, num_nodes, num_samples;
    size_t tmp;
    tsk_id_t *stack, *samples;
    tsk_treeseq_t *ts = tree->tree_sequence;
    tsk_id_t *sample_index_map = ts->sample_index_map;
    const tsk_id_t *list_left = tree->left_sample;
    const tsk_id_t *list_right = tree->right_sample;
    const tsk_id_t *list_next = tree->next_sample;
    tsk_id_t stop, sample_index;

    n = tsk_treeseq_get_num_samples(ts);
    num_nodes = tsk_treeseq_get_num_nodes(ts);
    stack = malloc(n * sizeof(tsk_id_t));
    samples = malloc(n * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(stack != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    for (u = 0; u < num_nodes; u++) {
        if (tree->left_child[u] == TSK_NULL && !tsk_treeseq_is_sample(ts, u)) {
            CU_ASSERT_EQUAL(list_left[u], TSK_NULL);
            CU_ASSERT_EQUAL(list_right[u], TSK_NULL);
        } else {
            stack_top = 0;
            num_samples = 0;
            stack[stack_top] = u;
            while (stack_top >= 0) {
                v = stack[stack_top];
                stack_top--;
                if (tsk_treeseq_is_sample(ts, v)) {
                    samples[num_samples] = v;
                    num_samples++;
                }
                for (v = tree->right_child[v]; v != TSK_NULL; v = tree->left_sib[v]) {
                    stack_top++;
                    stack[stack_top] = v;
                }
            }
            ret = tsk_tree_get_num_samples(tree, u, &tmp);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(num_samples, tmp);

            j = 0;
            sample_index = list_left[u];
            if (sample_index != TSK_NULL) {
                stop = list_right[u];
                while (true) {
                    CU_ASSERT_TRUE_FATAL(j < n);
                    CU_ASSERT_EQUAL_FATAL(sample_index, sample_index_map[samples[j]]);
                    j++;
                    if (sample_index == stop) {
                        break;
                    }
                    sample_index = list_next[sample_index];
                }
            }
            CU_ASSERT_EQUAL_FATAL(j, num_samples);
        }
    }
    free(stack);
    free(samples);
}

static void
verify_sample_sets(tsk_treeseq_t *ts)
{
    int ret;
    tsk_tree_t t;

    ret = tsk_tree_alloc(&t, ts, TSK_SAMPLE_COUNTS|TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL(ret, 0);

    for (ret = tsk_tree_first(&t); ret == 1; ret = tsk_tree_next(&t)) {
        verify_sample_sets_for_tree(&t);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (ret = tsk_tree_last(&t); ret == 1; ret = tsk_tree_prev(&t)) {
        verify_sample_sets_for_tree(&t);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_tree_free(&t);
}

static void
verify_tree_equals(tsk_treeseq_t *ts)
{
    int ret;
    tsk_tree_t *trees, t;
    size_t j, k;
    tsk_treeseq_t *other_ts = get_example_tree_sequence(
            10, 0, 100, 100.0, 1.0, 1.0, 0, NULL, 0);
    int flags[] = {0, TSK_SAMPLE_LISTS, TSK_SAMPLE_COUNTS,
        TSK_SAMPLE_LISTS | TSK_SAMPLE_COUNTS};

    trees = get_tree_list(ts);
    ret = tsk_tree_alloc(&t, other_ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < tsk_treeseq_get_num_trees(ts); j++) {
        ret = tsk_tree_equal(&t, &trees[j]);
        CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
        for (k = 0; k < tsk_treeseq_get_num_trees(ts); k++) {
            ret = tsk_tree_equal(&trees[j], &trees[k]);
            if (j == k) {
                CU_ASSERT_EQUAL_FATAL(ret, 0);
            } else {
                CU_ASSERT_EQUAL_FATAL(ret, 1);
            }
        }
    }
    ret = tsk_tree_free(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < sizeof(flags) / sizeof(int); j++) {
        ret = tsk_tree_alloc(&t, ts, flags[j]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        for (ret = tsk_tree_first(&t); ret == 1;
                ret = tsk_tree_next(&t)) {
            for (k = 0; k < tsk_treeseq_get_num_trees(ts); k++) {
                ret = tsk_tree_equal(&t, &trees[k]);
                if (t.index == k) {
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                } else {
                    CU_ASSERT_EQUAL_FATAL(ret, 1);
                }
            }
        }
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tree_free(&t);
        CU_ASSERT_EQUAL(ret, 0);
    }
    for (j = 0; j < tsk_treeseq_get_num_trees(ts); j++) {
        ret = tsk_tree_free(&trees[j]);
    }
    free(trees);
    tsk_treeseq_free(other_ts);
    free(other_ts);
}

static void
test_individual_nodes_from_examples(void)
{
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_individual_nodes(examples[j]);
        add_individuals(examples[j]);
        verify_individual_nodes(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_diff_iter_from_examples(void)
{
    /* tsk_treeseq_t **examples = get_example_tree_sequences(1); */
    /* uint32_t j; */

    /* CU_ASSERT_FATAL(examples != NULL); */
    /* for (j = 0; examples[j] != NULL; j++) { */
    /*     verify_tree_diffs(examples[j]); */
    /*     tsk_treeseq_free(examples[j]); */
    /*     free(examples[j]); */
    /* } */
    /* free(examples); */
}

static void
test_tree_iter_from_examples(void)
{

    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_trees_consistent(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_sample_sets_from_examples(void)
{
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_sample_sets(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_tree_equals_from_examples(void)
{
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_tree_equals(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_next_prev_from_examples(void)
{

    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_tree_next_prev(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_tsk_hapgen_from_examples(void)
{
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_hapgen(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_ld_from_examples(void)
{
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_ld(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_tsk_vargen_from_examples(void)
{
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_vargen(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_stats_from_examples(void)
{
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_stats(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_genealogical_nearest_neighbours_from_examples(void)
{
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_genealogical_nearest_neighbours(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_mean_descendants_from_examples(void)
{
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_mean_descendants(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}



static void
test_compute_mutation_parents_from_examples(void)
{
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_compute_mutation_parents(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
verify_simplify_errors(tsk_treeseq_t *ts)
{
    int ret;
    tsk_id_t *s;
    tsk_id_t u;
    tsk_treeseq_t subset;
    tsk_id_t sample[2];

    ret = tsk_treeseq_get_samples(ts, &s);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    memcpy(sample, s, 2 * sizeof(tsk_id_t));

    for (u = 0; u < (tsk_id_t) tsk_treeseq_get_num_nodes(ts); u++) {
        if (! tsk_treeseq_is_sample(ts, u)) {
            sample[1] = u;
            ret = tsk_treeseq_simplify(ts, sample, 2, 0, &subset, NULL);
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_SAMPLES);
        }
    }
    sample[0] = -1;
    ret = tsk_treeseq_simplify(ts, sample, 2, 0, &subset, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    sample[0] = s[0];
    sample[1] = s[0];
    ret = tsk_treeseq_simplify(ts, sample, 2, 0, &subset, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SAMPLE);
}

static void
test_simplify_from_examples(void)
{
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_simplify(examples[j]);
        if (tsk_treeseq_get_num_migrations(examples[j]) == 0) {
            /* Migrations are not supported at the moment, so skip these tests
             * rather than complicate them */
            verify_simplify_errors(examples[j]);
        }
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_reduce_topology_from_examples(void)
{
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_reduce_topology(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
verify_newick(tsk_treeseq_t *ts)
{
    /* int ret, err; */
    /* tsk_tree_t t; */
    /* tsk_id_t root; */
    /* size_t precision = 4; */
    /* size_t buffer_size = 1024 * 1024; */
    /* char *newick = malloc(buffer_size); */
    /* size_t j, size; */

    /* CU_ASSERT_FATAL(newick != NULL); */

    /* ret = tsk_tree_alloc(&t, ts, 0); */
    /* CU_ASSERT_EQUAL_FATAL(ret, 0); */
    /* ret = tsk_tree_first(&t); */
    /* CU_ASSERT_FATAL(ret == 1); */
    /* for (root = t.left_root; root != TSK_NULL; root = t.right_sib[root]) { */
    /*     err = tsk_tree_get_newick(&t, root, precision, 0, buffer_size, newick); */
    /*     CU_ASSERT_EQUAL_FATAL(err, 0); */
    /*     size = strlen(newick); */
    /*     CU_ASSERT_TRUE(size > 0); */
    /*     CU_ASSERT_TRUE(size < buffer_size); */
    /*     for (j = 0; j <= size; j++) { */
    /*         err = tsk_tree_get_newick(&t, root, precision, 0, j, newick); */
    /*         CU_ASSERT_EQUAL_FATAL(err, TSK_ERR_BUFFER_OVERFLOW); */
    /*     } */
    /*     err = tsk_tree_get_newick(&t, root, precision, 0, size + 1, newick); */
    /*     CU_ASSERT_EQUAL_FATAL(err, 0); */
    /* } */

    /* for (ret = tsk_tree_first(&t); ret == 1; ret = tsk_tree_next(&t)) { */
    /*     for (root = t.left_root; root != TSK_NULL; root = t.right_sib[root]) { */
    /*         err = tsk_tree_get_newick(&t, root, precision, 0, 0, NULL); */
    /*         CU_ASSERT_EQUAL_FATAL(err, TSK_ERR_BAD_PARAM_VALUE); */
    /*         err = tsk_tree_get_newick(&t, root, precision, 0, buffer_size, newick); */
    /*         CU_ASSERT_EQUAL_FATAL(err, 0); */
    /*         size = strlen(newick); */
    /*         CU_ASSERT_EQUAL(newick[size - 1], ';'); */
    /*     } */
    /* } */
    /* CU_ASSERT_EQUAL_FATAL(ret, 0); */

    /* tsk_tree_free(&t); */
    /* free(newick); */
}

static void
test_newick_from_examples(void)
{
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_newick(examples[j]);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
verify_tree_sequences_equal(tsk_treeseq_t *ts1, tsk_treeseq_t *ts2,
        bool check_migrations, bool check_mutations,
        bool check_provenance)
{
    int ret, err1, err2;
    size_t j;
    tsk_edge_t r1, r2;
    tsk_node_t n1, n2;
    tsk_migration_t m1, m2;
    tsk_provenance_t p1, p2;
    tsk_individual_t i1, i2;
    tsk_population_t pop1, pop2;
    size_t num_mutations = tsk_treeseq_get_num_mutations(ts1);
    tsk_site_t site_1, site_2;
    tsk_mutation_t mutation_1, mutation_2;
    tsk_tree_t t1, t2;

    /* tsk_treeseq_print_state(ts1, stdout); */
    /* tsk_treeseq_print_state(ts2, stdout); */

    CU_ASSERT_EQUAL(
        tsk_treeseq_get_num_samples(ts1),
        tsk_treeseq_get_num_samples(ts2));
    CU_ASSERT_EQUAL(
        tsk_treeseq_get_sequence_length(ts1),
        tsk_treeseq_get_sequence_length(ts2));
    CU_ASSERT_EQUAL(
        tsk_treeseq_get_num_edges(ts1),
        tsk_treeseq_get_num_edges(ts2));
    CU_ASSERT_EQUAL(
        tsk_treeseq_get_num_nodes(ts1),
        tsk_treeseq_get_num_nodes(ts2));
    CU_ASSERT_EQUAL(
        tsk_treeseq_get_num_trees(ts1),
        tsk_treeseq_get_num_trees(ts2));

    for (j = 0; j < tsk_treeseq_get_num_nodes(ts1); j++) {
        ret = tsk_treeseq_get_node(ts1, j, &n1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tsk_treeseq_get_node(ts2, j, &n2);
        CU_ASSERT_EQUAL(ret, 0);
        verify_nodes_equal(&n1, &n2);
    }
    for (j = 0; j < tsk_treeseq_get_num_edges(ts1); j++) {
        ret = tsk_treeseq_get_edge(ts1, j, &r1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tsk_treeseq_get_edge(ts2, j, &r2);
        CU_ASSERT_EQUAL(ret, 0);
        verify_edges_equal(&r1, &r2, 1.0);
    }
    if (check_mutations) {
        CU_ASSERT_EQUAL_FATAL(
            tsk_treeseq_get_num_sites(ts1),
            tsk_treeseq_get_num_sites(ts2));
        for (j = 0; j < tsk_treeseq_get_num_sites(ts1); j++) {
            ret = tsk_treeseq_get_site(ts1, j, &site_1);
            CU_ASSERT_EQUAL(ret, 0);
            ret = tsk_treeseq_get_site(ts2, j, &site_2);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(site_1.position, site_2.position);
            CU_ASSERT_EQUAL(site_1.ancestral_state_length, site_2.ancestral_state_length);
            CU_ASSERT_NSTRING_EQUAL(site_1.ancestral_state, site_2.ancestral_state,
                    site_1.ancestral_state_length);
            CU_ASSERT_EQUAL(site_1.metadata_length, site_2.metadata_length);
            CU_ASSERT_NSTRING_EQUAL(site_1.metadata, site_2.metadata,
                    site_1.metadata_length);
        }
        CU_ASSERT_EQUAL_FATAL(
            tsk_treeseq_get_num_mutations(ts1),
            tsk_treeseq_get_num_mutations(ts2));
        for (j = 0; j < num_mutations; j++) {
            ret = tsk_treeseq_get_mutation(ts1, j, &mutation_1);
            CU_ASSERT_EQUAL(ret, 0);
            ret = tsk_treeseq_get_mutation(ts2, j, &mutation_2);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(mutation_1.id, j);
            CU_ASSERT_EQUAL(mutation_1.id, mutation_2.id);
            CU_ASSERT_EQUAL(mutation_1.site, mutation_2.site);
            CU_ASSERT_EQUAL(mutation_1.node, mutation_2.node);
            CU_ASSERT_EQUAL_FATAL(mutation_1.parent, mutation_2.parent);
            CU_ASSERT_EQUAL_FATAL(mutation_1.derived_state_length,
                    mutation_2.derived_state_length);
            CU_ASSERT_NSTRING_EQUAL(mutation_1.derived_state,
                    mutation_2.derived_state, mutation_1.derived_state_length);
            CU_ASSERT_EQUAL_FATAL(mutation_1.metadata_length,
                    mutation_2.metadata_length);
            CU_ASSERT_NSTRING_EQUAL(mutation_1.metadata,
                    mutation_2.metadata, mutation_1.metadata_length);
        }
    }
    if (check_migrations) {
        CU_ASSERT_EQUAL_FATAL(
            tsk_treeseq_get_num_migrations(ts1),
            tsk_treeseq_get_num_migrations(ts2));
        for (j = 0; j < tsk_treeseq_get_num_migrations(ts1); j++) {
            ret = tsk_treeseq_get_migration(ts1, j, &m1);
            CU_ASSERT_EQUAL(ret, 0);
            ret = tsk_treeseq_get_migration(ts2, j, &m2);
            CU_ASSERT_EQUAL(ret, 0);
            verify_migrations_equal(&m1, &m2, 1.0);
        }
    }
    if (check_provenance) {
        CU_ASSERT_EQUAL_FATAL(
            tsk_treeseq_get_num_provenances(ts1),
            tsk_treeseq_get_num_provenances(ts2));
        for (j = 0; j < tsk_treeseq_get_num_provenances(ts1); j++) {
            ret = tsk_treeseq_get_provenance(ts1, j, &p1);
            CU_ASSERT_EQUAL(ret, 0);
            ret = tsk_treeseq_get_provenance(ts2, j, &p2);
            CU_ASSERT_EQUAL(ret, 0);
            verify_provenances_equal(&p1, &p2);
        }
    }

    CU_ASSERT_EQUAL_FATAL(
        tsk_treeseq_get_num_individuals(ts1),
        tsk_treeseq_get_num_individuals(ts2));
    for (j = 0; j < tsk_treeseq_get_num_individuals(ts1); j++) {
        ret = tsk_treeseq_get_individual(ts1, j, &i1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tsk_treeseq_get_individual(ts2, j, &i2);
        CU_ASSERT_EQUAL(ret, 0);
        verify_individuals_equal(&i1, &i2);
    }

    CU_ASSERT_EQUAL_FATAL(
        tsk_treeseq_get_num_populations(ts1),
        tsk_treeseq_get_num_populations(ts2));
    for (j = 0; j < tsk_treeseq_get_num_populations(ts1); j++) {
        ret = tsk_treeseq_get_population(ts1, j, &pop1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tsk_treeseq_get_population(ts2, j, &pop2);
        CU_ASSERT_EQUAL(ret, 0);
        verify_populations_equal(&pop1, &pop2);
    }

    ret = tsk_tree_alloc(&t1, ts1, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_alloc(&t2, ts2, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&t1);
    CU_ASSERT_EQUAL(ret, 1);
    ret = tsk_tree_first(&t2);
    CU_ASSERT_EQUAL(ret, 1);
    while (1) {
        err1 = tsk_tree_next(&t1);
        err2 = tsk_tree_next(&t2);
        CU_ASSERT_EQUAL_FATAL(err1, err2);
        if (err1 != 1) {
            break;
        }
    }
    tsk_tree_free(&t1);
    tsk_tree_free(&t2);
}

static void
test_save_empty_kas(void)
{
    int ret;
    tsk_treeseq_t ts1, ts2;
    double sequence_length = 1234.00;
    tsk_tbl_collection_t tables;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = sequence_length;

    ret = tsk_treeseq_alloc(&ts1, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_dump(&ts1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_empty_tree_sequence(&ts1, sequence_length);
    ret = tsk_treeseq_load(&ts2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_empty_tree_sequence(&ts2, sequence_length);

    tsk_treeseq_free(&ts1);
    tsk_treeseq_free(&ts2);
    tsk_tbl_collection_free(&tables);
}

static void
test_save_kas(void)
{
    int ret;
    size_t j, k;
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    tsk_treeseq_t ts2;
    tsk_treeseq_t *ts1;
    char *file_uuid;
    int dump_flags[] = {0};

    CU_ASSERT_FATAL(examples != NULL);

    for (j = 0; examples[j] != NULL; j++) {
        ts1 = examples[j];
        file_uuid = tsk_treeseq_get_file_uuid(ts1);
        CU_ASSERT_EQUAL_FATAL(file_uuid, NULL);
        for (k = 0; k < sizeof(dump_flags) / sizeof(int); k++) {
            ret = tsk_treeseq_dump(ts1, _tmp_file_name, dump_flags[k]);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            ret = tsk_treeseq_load(&ts2, _tmp_file_name, TSK_LOAD_EXTENDED_CHECKS);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            verify_tree_sequences_equal(ts1, &ts2, true, true, true);
            tsk_treeseq_print_state(&ts2, _devnull);
            verify_hapgen(&ts2);
            verify_vargen(&ts2);
            file_uuid = tsk_treeseq_get_file_uuid(&ts2);
            CU_ASSERT_NOT_EQUAL_FATAL(file_uuid, NULL);
            CU_ASSERT_EQUAL(strlen(file_uuid), TSK_UUID_SIZE);
            tsk_treeseq_free(&ts2);
        }
        tsk_treeseq_free(ts1);
        free(ts1);
    }
    free(examples);
}

static void
test_save_kas_tables(void)
{
    int ret;
    size_t j, k;
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    tsk_treeseq_t *ts1;
    tsk_tbl_collection_t t1, t2;
    int dump_flags[] = {0};

    CU_ASSERT_FATAL(examples != NULL);

    for (j = 0; examples[j] != NULL; j++) {
        ts1 = examples[j];
        ret = tsk_tbl_collection_alloc(&t1, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_treeseq_copy_tables(ts1, &t1);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(t1.file_uuid, NULL);
        for (k = 0; k < sizeof(dump_flags) / sizeof(int); k++) {
            ret = tsk_tbl_collection_dump(&t1, _tmp_file_name, dump_flags[k]);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            ret = tsk_tbl_collection_alloc(&t2, 0);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            ret = tsk_tbl_collection_load(&t2, _tmp_file_name, 0);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_TRUE(tsk_tbl_collection_equals(&t1, &t2));
            CU_ASSERT_EQUAL_FATAL(t1.file_uuid, NULL);
            CU_ASSERT_NOT_EQUAL_FATAL(t2.file_uuid, NULL);
            CU_ASSERT_EQUAL(strlen(t2.file_uuid), TSK_UUID_SIZE);
            tsk_tbl_collection_free(&t2);
        }
        tsk_tbl_collection_free(&t1);
        tsk_treeseq_free(ts1);
        free(ts1);
    }
    free(examples);
}

static void
test_sort_tables(void)
{
    int ret;
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    tsk_treeseq_t ts2;
    tsk_treeseq_t *ts1;
    size_t j, k, start, starts[3];
    tsk_tbl_collection_t tables;
    int load_flags = TSK_BUILD_INDEXES;
    tsk_id_t tmp_node;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(examples != NULL);

    for (j = 0; examples[j] != NULL; j++) {
        ts1 = examples[j];

        ret = tsk_treeseq_copy_tables(ts1, &tables);
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        /* Check the input validation */
        ret = tsk_tbl_collection_sort(NULL, 0, 0);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
        /* Check edge sorting */
        if (tables.edges->num_rows == 2) {
            starts[0] = 0;
            starts[1] = 0;
            starts[2] = 0;
        } else {
            starts[0] = 0;
            starts[1] = tables.edges->num_rows / 2;
            starts[2] = tables.edges->num_rows - 2;
        }
        for (k = 0; k < 3; k++) {
            start = starts[k];
            unsort_edges(tables.edges, start);
            ret = tsk_treeseq_alloc(&ts2, &tables, load_flags);
            CU_ASSERT_NOT_EQUAL_FATAL(ret, 0);
            tsk_treeseq_free(&ts2);

            ret = tsk_tbl_collection_sort(&tables, 0, 0);
            CU_ASSERT_EQUAL_FATAL(ret, 0);

            ret = tsk_treeseq_alloc(&ts2, &tables, load_flags);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            verify_tree_sequences_equal(ts1, &ts2, true, true, false);
            tsk_treeseq_free(&ts2);
        }

        /* A start value of num_tables.edges should have no effect */
        ret = tsk_tbl_collection_sort(&tables, tables.edges->num_rows, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_treeseq_alloc(&ts2, &tables, load_flags);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        verify_tree_sequences_equal(ts1, &ts2, true, true, false);
        tsk_treeseq_free(&ts2);

        if (tables.sites->num_rows > 1) {
            /* Check site sorting */
            unsort_sites(tables.sites, tables.mutations);
            ret = tsk_treeseq_alloc(&ts2, &tables, load_flags);
            CU_ASSERT_NOT_EQUAL(ret, 0);
            tsk_treeseq_free(&ts2);

            ret = tsk_tbl_collection_sort(&tables, 0, 0);
            CU_ASSERT_EQUAL_FATAL(ret, 0);

            ret = tsk_treeseq_alloc(&ts2, &tables, load_flags);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            verify_tree_sequences_equal(ts1, &ts2, true, true, false);
            tsk_treeseq_free(&ts2);

            /* Check for site bounds error */
            tables.mutations->site[0] = tables.sites->num_rows;
            ret = tsk_tbl_collection_sort(&tables, 0, 0);
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
            tables.mutations->site[0] = 0;
            ret = tsk_tbl_collection_sort(&tables, 0, 0);
            CU_ASSERT_EQUAL_FATAL(ret, 0);

            /* Check for edge node bounds error */
            tmp_node = tables.edges->parent[0];
            tables.edges->parent[0] = tables.nodes->num_rows;
            ret = tsk_tbl_collection_sort(&tables, 0, 0);
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
            tables.edges->parent[0] = tmp_node;
            ret = tsk_tbl_collection_sort(&tables, 0, 0);
            CU_ASSERT_EQUAL_FATAL(ret, 0);

            /* Check for mutation node bounds error */
            tmp_node = tables.mutations->node[0];
            tables.mutations->node[0] = tables.nodes->num_rows;
            ret = tsk_tbl_collection_sort(&tables, 0, 0);
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
            tables.mutations->node[0] = tmp_node;

            /* Check for mutation parent bounds error */
            tables.mutations->parent[0] = tables.mutations->num_rows;
            ret = tsk_tbl_collection_sort(&tables, 0, 0);
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);
            tables.mutations->parent[0] = TSK_NULL;
            ret = tsk_tbl_collection_sort(&tables, 0, 0);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
        }
        tsk_treeseq_free(ts1);
        free(ts1);
    }
    free(examples);
    tsk_tbl_collection_free(&tables);
}
static void
test_dump_tables(void)
{
    int ret;
    tsk_treeseq_t **examples = get_example_tree_sequences(1);
    tsk_treeseq_t ts2;
    tsk_treeseq_t *ts1;
    tsk_tbl_collection_t tables;
    size_t j;
    int load_flags = TSK_BUILD_INDEXES;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(examples != NULL);

    for (j = 0; examples[j] != NULL; j++) {
        ts1 = examples[j];

        ret = tsk_treeseq_copy_tables(ts1, NULL);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);

        ret = tsk_treeseq_copy_tables(ts1, &tables);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_treeseq_alloc(&ts2, &tables, load_flags);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        verify_tree_sequences_equal(ts1, &ts2, true, true, true);
        tsk_treeseq_print_state(&ts2, _devnull);
        tsk_treeseq_free(&ts2);
        tsk_treeseq_free(ts1);
        free(ts1);
    }

    free(examples);
    tsk_tbl_collection_free(&tables);
}

static void
test_dump_tables_kas(void)
{
    int ret;
    size_t k;
    tsk_treeseq_t *ts1, ts2, ts3, **examples;
    tsk_tbl_collection_t tables;
    int load_flags = TSK_BUILD_INDEXES;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    examples = get_example_tree_sequences(1);
    for (k = 0; examples[k] != NULL; k++) {
        ts1 = examples[k];
        CU_ASSERT_FATAL(ts1 != NULL);
        ret = tsk_treeseq_copy_tables(ts1, &tables);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_treeseq_alloc(&ts2, &tables, load_flags);
        ret = tsk_treeseq_dump(&ts2, _tmp_file_name, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_treeseq_load(&ts3, _tmp_file_name, TSK_LOAD_EXTENDED_CHECKS);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        verify_tree_sequences_equal(ts1, &ts3, true, true, true);
        tsk_treeseq_print_state(&ts2, _devnull);

        tsk_treeseq_free(&ts2);
        tsk_treeseq_free(&ts3);
        tsk_treeseq_free(ts1);
        free(ts1);
    }
    free(examples);
    tsk_tbl_collection_free(&tables);
}

void
test_tsk_tbl_collection_simplify_errors(void)
{
    int ret;
    tsk_tbl_collection_t tables;
    tsk_id_t samples[] = {0, 1};


    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_site_tbl_add_row(tables.sites, 0, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_tbl_add_row(tables.sites, 0, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_tbl_collection_simplify(&tables, samples, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SITE_POSITION);

    /* Out of order positions */
    tables.sites->position[0] = 0.5;
    ret = tsk_tbl_collection_simplify(&tables, samples, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNSORTED_SITES);

    /* Position out of bounds */
    tables.sites->position[0] = 1.5;
    ret = tsk_tbl_collection_simplify(&tables, samples, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_SITE_POSITION);

    /* TODO More tests for this: see
     * https://github.com/tskit-dev/msprime/issues/517 */

    tsk_tbl_collection_free(&tables);

}
void
test_tsk_tbl_collection_position_errors(void)
{
    int ret;
    int j;
    tsk_tbl_collection_t t1, t2;
    tsk_tbl_collection_position_t pos1, pos2;
    tsk_treeseq_t **examples = get_example_tree_sequences(1);

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        // set-up
        ret = tsk_tbl_collection_alloc(&t1, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tbl_collection_alloc(&t2, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_treeseq_copy_tables(examples[j], &t1);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tbl_collection_copy(&t1, &t2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        tsk_tbl_collection_record_position(&t1, &pos1);

        // for each table, add a new row to t2, bookmark that location,
        // then try to reset t1 to this illegal location

        // individuals
        tsk_individual_tbl_add_row(t2.individuals, 0, NULL, 0, NULL, 0);
        tsk_tbl_collection_record_position(&t2, &pos2);
        ret = tsk_tbl_collection_reset_position(&t1, &pos2);
        CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_TABLE_POSITION);
        ret = tsk_tbl_collection_reset_position(&t2, &pos1);
        CU_ASSERT_EQUAL(ret, 0);

        // nodes
        tsk_node_tbl_add_row(t2.nodes, 0, 1.2, 0, -1, NULL, 0);
        tsk_tbl_collection_record_position(&t2, &pos2);
        ret = tsk_tbl_collection_reset_position(&t1, &pos2);
        CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_TABLE_POSITION);
        ret = tsk_tbl_collection_reset_position(&t2, &pos1);
        CU_ASSERT_EQUAL(ret, 0);

        // edges
        tsk_edge_tbl_add_row(t2.edges, 0.1, 0.4, 0, 3);
        tsk_tbl_collection_record_position(&t2, &pos2);
        ret = tsk_tbl_collection_reset_position(&t1, &pos2);
        CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_TABLE_POSITION);
        ret = tsk_tbl_collection_reset_position(&t2, &pos1);
        CU_ASSERT_EQUAL(ret, 0);

        // migrations
        tsk_migration_tbl_add_row(t2.migrations, 0.1, 0.2, 2, 1, 2, 1.2);
        tsk_tbl_collection_record_position(&t2, &pos2);
        ret = tsk_tbl_collection_reset_position(&t1, &pos2);
        CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_TABLE_POSITION);
        ret = tsk_tbl_collection_reset_position(&t2, &pos1);
        CU_ASSERT_EQUAL(ret, 0);

        // sites
        tsk_site_tbl_add_row(t2.sites, 0.3, "A", 1, NULL, 0);
        tsk_tbl_collection_record_position(&t2, &pos2);
        ret = tsk_tbl_collection_reset_position(&t1, &pos2);
        CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_TABLE_POSITION);
        ret = tsk_tbl_collection_reset_position(&t2, &pos1);
        CU_ASSERT_EQUAL(ret, 0);

        // mutations
        tsk_mutation_tbl_add_row(t2.mutations, 0, 1, -1, "X", 1, NULL, 0);
        tsk_tbl_collection_record_position(&t2, &pos2);
        ret = tsk_tbl_collection_reset_position(&t1, &pos2);
        CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_TABLE_POSITION);
        ret = tsk_tbl_collection_reset_position(&t2, &pos1);
        CU_ASSERT_EQUAL(ret, 0);

        // populations
        tsk_population_tbl_add_row(t2.populations, NULL, 0);
        tsk_tbl_collection_record_position(&t2, &pos2);
        ret = tsk_tbl_collection_reset_position(&t1, &pos2);
        CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_TABLE_POSITION);
        ret = tsk_tbl_collection_reset_position(&t2, &pos1);
        CU_ASSERT_EQUAL(ret, 0);

        // provenance
        tsk_provenance_tbl_add_row(t2.provenances, "abc", 3, NULL, 0);
        tsk_tbl_collection_record_position(&t2, &pos2);
        ret = tsk_tbl_collection_reset_position(&t1, &pos2);
        CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_TABLE_POSITION);
        ret = tsk_tbl_collection_reset_position(&t2, &pos1);
        CU_ASSERT_EQUAL(ret, 0);

        tsk_tbl_collection_free(&t1);
        tsk_tbl_collection_free(&t2);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

void
test_tsk_tbl_collection_position(void)
{
    int ret;
    int j, k;
    tsk_treeseq_t **examples;
    tsk_tbl_collection_t t1, t2, t3;
    tsk_tbl_collection_position_t pos1, pos2;

    examples = get_example_tree_sequences(1);
    CU_ASSERT_FATAL(examples != NULL);

    for (j = 0; examples[j] != NULL; j++) {
        ret = tsk_tbl_collection_alloc(&t1, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tbl_collection_alloc(&t2, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tbl_collection_alloc(&t3, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        ret = tsk_treeseq_copy_tables(examples[j], &t1);

        // bookmark at pos1
        tsk_tbl_collection_record_position(&t1, &pos1);
        // copy to t2
        ret = tsk_tbl_collection_copy(&t1, &t2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        // resetting position should do nothing
        ret = tsk_tbl_collection_reset_position(&t2, &pos1);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_TRUE(tsk_tbl_collection_equals(&t1, &t2));
        // add more rows to t2
        // (they don't have to make sense for this test)
        for (k = 0; k < 3; k++) {
            tsk_node_tbl_add_row(t2.nodes, 0, 1.2, 0, -1, NULL, 0);
            tsk_node_tbl_add_row(t2.nodes, 0, 1.2, k, -1, NULL, 0);
            tsk_edge_tbl_add_row(t2.edges, 0.1, 0.5, k, k+1);
            tsk_edge_tbl_add_row(t2.edges, 0.3, 0.8, k, k+2);
        }
        // bookmark at pos2
        tsk_tbl_collection_record_position(&t2, &pos2);
        // copy to t3
        ret = tsk_tbl_collection_copy(&t2, &t3);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        // add more rows to t3
        for (k = 0; k < 3; k++) {
            tsk_node_tbl_add_row(t3.nodes, 0, 1.2, k+5, -1, NULL, 0);
            tsk_site_tbl_add_row(t3.sites, 0.2, "A", 1, NULL, 0);
            tsk_site_tbl_add_row(t3.sites, 0.2, "C", 1, NULL, 0);
            tsk_mutation_tbl_add_row(t3.mutations, 0, k, -1, "T", 1, NULL, 0);
            tsk_migration_tbl_add_row(t3.migrations, 0.0, 0.5, 1, 0, 1, 1.2);
            tsk_individual_tbl_add_row(t3.individuals, k, NULL, 0, NULL, 0);
            tsk_population_tbl_add_row(t3.populations, "X", 1);
            tsk_provenance_tbl_add_row(t3.provenances, "abc", 3, NULL, 0);
        }
        // now resetting t3 to pos2 should equal t2
        ret = tsk_tbl_collection_reset_position(&t3, &pos2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_TRUE(tsk_tbl_collection_equals(&t2, &t3));
        // and resetting to pos1 should equal t1
        ret = tsk_tbl_collection_reset_position(&t3, &pos1);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_TRUE(tsk_tbl_collection_equals(&t1, &t3));

        ret = tsk_tbl_collection_clear(&t1);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(t1.individuals->num_rows, 0);
        CU_ASSERT_EQUAL(t1.populations->num_rows, 0);
        CU_ASSERT_EQUAL(t1.nodes->num_rows, 0);
        CU_ASSERT_EQUAL(t1.edges->num_rows, 0);
        CU_ASSERT_EQUAL(t1.migrations->num_rows, 0);
        CU_ASSERT_EQUAL(t1.sites->num_rows, 0);
        CU_ASSERT_EQUAL(t1.mutations->num_rows, 0);
        CU_ASSERT_EQUAL(t1.provenances->num_rows, 0);

        tsk_tbl_collection_free(&t1);
        tsk_tbl_collection_free(&t2);
        tsk_tbl_collection_free(&t3);
        tsk_treeseq_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static int
msprime_suite_init(void)
{
    int fd;
    static char template[] = "/tmp/tsk_c_test_XXXXXX";

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
        {"test_vcf", test_vcf},
        {"test_vcf_no_mutations", test_vcf_no_mutations},
        {"test_node_metadata", test_node_metadata},

        {"test_single_tree_newick", test_single_tree_newick},


        {"test_diff_iter_from_examples", test_diff_iter_from_examples},
        {"test_tree_iter_from_examples", test_tree_iter_from_examples},
        {"test_tree_equals_from_examples", test_tree_equals_from_examples},
        {"test_next_prev_from_examples", test_next_prev_from_examples},
        {"test_sample_sets_from_examples", test_sample_sets_from_examples},
        {"test_tsk_hapgen_from_examples", test_tsk_hapgen_from_examples},
        {"test_tsk_vargen_from_examples", test_tsk_vargen_from_examples},
        {"test_newick_from_examples", test_newick_from_examples},
        {"test_stats_from_examples", test_stats_from_examples},
        {"test_compute_mutation_parents_from_examples",
            test_compute_mutation_parents_from_examples},
        {"test_individual_nodes_from_examples",
            test_individual_nodes_from_examples},
        {"test_ld_from_examples", test_ld_from_examples},
        {"test_simplify_from_examples", test_simplify_from_examples},
        {"test_reduce_topology_from_examples", test_reduce_topology_from_examples},
        {"test_save_empty_kas", test_save_empty_kas},
        {"test_save_kas", test_save_kas},
        {"test_save_kas_tables", test_save_kas_tables},
        {"test_dump_tables", test_dump_tables},
        {"test_sort_tables", test_sort_tables},
        {"test_dump_tables_kas", test_dump_tables_kas},

        {"test_tsk_tbl_collection_position", test_tsk_tbl_collection_position},
        {"test_tsk_tbl_collection_position_errors", test_tsk_tbl_collection_position_errors},

        {"test_genealogical_nearest_neighbours_from_examples",
            test_genealogical_nearest_neighbours_from_examples},
        {"test_mean_descendants_from_examples", test_mean_descendants_from_examples},
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
