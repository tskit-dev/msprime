#define _GNU_SOURCE
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "testlib.h"

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

void
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

void
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

void
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

void
parse_mutations(const char *text, tsk_mutation_tbl_t *mutation_table)
{
    int ret;
    size_t c, k;
    size_t MAX_LINE = 1024;
    char line[MAX_LINE];
    const char *whitespace = " \t";
    char *p;
    tsk_id_t node, site, parent;
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

void
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

void
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

void
unsort_edges(tsk_edge_tbl_t *edges, size_t start)
{
    size_t j, k;
    size_t n = edges->num_rows - start;
    tsk_edge_t *buff = malloc(n * sizeof(tsk_edge_t));
    CU_ASSERT_FATAL(buff != NULL);

    for (j = 0; j < n; j++) {
        k = start + j;
        buff[j].left = edges->left[k];
        buff[j].right = edges->right[k];
        buff[j].parent = edges->parent[k];
        buff[j].child = edges->child[k];
    }
    for (j = 0; j < n; j++) {
        k = start + j;
        edges->left[k] = buff[n - j - 1].left;
        edges->right[k] = buff[n - j - 1].right;
        edges->parent[k] = buff[n - j - 1].parent;
        edges->child[k] = buff[n - j - 1].child;
    }
    free(buff);
}

static int
tskit_suite_init(void)
{
    int fd = -1;
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
tskit_suite_cleanup(void)
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
handle_cunit_error(void)
{
    fprintf(stderr, "CUnit error occured: %d: %s\n", CU_get_error(), CU_get_error_msg());
    exit(EXIT_FAILURE);
}

int
test_main(CU_TestInfo *tests, int argc, char **argv)
{
    int ret;
    CU_pTest test;
    CU_pSuite suite;
    CU_SuiteInfo suites[] = {
        {
            .pName = "tskit",
            .pInitFunc = tskit_suite_init,
            .pCleanupFunc = tskit_suite_cleanup,
            .pTests = tests,
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
        suite = CU_get_suite_by_name("tskit", CU_get_registry());
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
        printf("usage: %s <test_name>\n", argv[0]);
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
