#ifndef __TESTLIB_H__
#define __TESTLIB_H__

#include <stdio.h>

#include <CUnit/Basic.h>
#include "tsk_trees.h"

/* Global variables used in the test suite */

char * _tmp_file_name;
FILE * _devnull;

int test_main(CU_TestInfo *tests, int argc, char **argv);

void tsk_treeseq_from_text(tsk_treeseq_t *ts,
        double sequence_length,
        const char *nodes, const char *edges, const char *migrations,
        const char *sites, const char *mutations,
        const char *individuals, const char *provenance);

void parse_nodes(const char *text, tsk_node_tbl_t *node_table);
void parse_edges(const char *text, tsk_edge_tbl_t *edge_table);
void parse_sites(const char *text, tsk_site_tbl_t *site_table);
void parse_mutations(const char *text, tsk_mutation_tbl_t *mutation_table);
void parse_individuals(const char *text, tsk_individual_tbl_t *individual_table);

void unsort_edges(tsk_edge_tbl_t *edges, size_t start);

extern const char *single_tree_ex_nodes;
extern const char *single_tree_ex_edges;
extern const char *single_tree_ex_sites;
extern const char *single_tree_ex_mutations;

extern const char *nonbinary_ex_nodes;
extern const char *nonbinary_ex_edges;
extern const char *nonbinary_ex_sites;
extern const char *nonbinary_ex_mutations;

extern const char *unary_ex_nodes;
extern const char *unary_ex_edges;
extern const char *unary_ex_sites;
extern const char *unary_ex_mutations;

extern const char *internal_sample_ex_nodes;
extern const char *internal_sample_ex_edges;
extern const char *internal_sample_ex_sites;
extern const char *internal_sample_ex_mutations;

extern const char *paper_ex_nodes;
extern const char *paper_ex_edges;
extern const char *paper_ex_sites;
extern const char *paper_ex_mutations;
extern const char *paper_ex_individuals;

#endif
