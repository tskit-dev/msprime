/*
 * MIT License
 *
 * Copyright (c) 2019-2022 Tskit Developers
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "testlib.h"
#include <tskit/trees.h>
#include <tskit/genotypes.h>

#include <unistd.h>
#include <stdlib.h>

/*=======================================================
 * Verification utilities.
 *======================================================*/

/* Checks if the specified trees are topologically equivalent, i.e, represent
 * the same tree without checking state specific to seeking.*/
static void
check_trees_equal(tsk_tree_t *self, tsk_tree_t *other)
{
    tsk_size_t N = self->num_nodes;

    CU_ASSERT_FATAL(self->tree_sequence == other->tree_sequence);
    CU_ASSERT_FATAL(self->index == other->index);
    CU_ASSERT_FATAL(self->interval.left == other->interval.left);
    CU_ASSERT_FATAL(self->interval.right == other->interval.right);
    CU_ASSERT_FATAL(self->sites_length == other->sites_length);
    CU_ASSERT_FATAL(self->sites == other->sites);
    CU_ASSERT_FATAL(self->samples == other->samples);
    CU_ASSERT_FATAL(self->num_edges == other->num_edges);
    CU_ASSERT_FATAL(tsk_memcmp(self->parent, other->parent, N * sizeof(tsk_id_t)) == 0);
    CU_ASSERT_FATAL(tsk_tree_equals(self, other));
}

static void
check_trees_identical(tsk_tree_t *self, tsk_tree_t *other)
{
    tsk_size_t N = self->num_nodes;

    check_trees_equal(self, other);
    CU_ASSERT_FATAL(self->left_index == other->left_index);
    CU_ASSERT_FATAL(self->right_index == other->right_index);
    CU_ASSERT_FATAL(self->direction == other->direction);

    CU_ASSERT_FATAL(
        tsk_memcmp(self->left_child, other->left_child, N * sizeof(tsk_id_t)) == 0);
    CU_ASSERT_FATAL(
        tsk_memcmp(self->right_child, other->right_child, N * sizeof(tsk_id_t)) == 0);
    CU_ASSERT_FATAL(
        tsk_memcmp(self->left_sib, other->left_sib, N * sizeof(tsk_id_t)) == 0);
    CU_ASSERT_FATAL(
        tsk_memcmp(self->right_sib, other->right_sib, N * sizeof(tsk_id_t)) == 0);

    CU_ASSERT_EQUAL_FATAL(self->num_samples == NULL, other->num_samples == NULL)
    CU_ASSERT_EQUAL_FATAL(
        self->num_tracked_samples == NULL, other->num_tracked_samples == NULL)
    if (self->num_samples != NULL) {
        CU_ASSERT_FATAL(tsk_memcmp(self->num_samples, other->num_samples,
                            N * sizeof(*self->num_samples))
                        == 0);
        CU_ASSERT_FATAL(tsk_memcmp(self->num_tracked_samples, other->num_tracked_samples,
                            N * sizeof(*self->num_tracked_samples))
                        == 0);
    }

    CU_ASSERT_EQUAL_FATAL(self->left_sample == NULL, other->left_sample == NULL)
    CU_ASSERT_EQUAL_FATAL(self->right_sample == NULL, other->left_sample == NULL)
    CU_ASSERT_EQUAL_FATAL(self->next_sample == NULL, other->next_sample == NULL)
    if (self->left_sample != NULL) {
        CU_ASSERT_FATAL(tsk_memcmp(self->left_sample, other->left_sample,
                            N * sizeof(*self->left_sample))
                        == 0);
        CU_ASSERT_FATAL(tsk_memcmp(self->right_sample, other->right_sample,
                            N * sizeof(*self->right_sample))
                        == 0);
        CU_ASSERT_FATAL(
            tsk_memcmp(self->next_sample, other->next_sample,
                self->tree_sequence->num_samples * sizeof(*self->next_sample))
            == 0);
    }
}

static void
verify_compute_mutation_parents(tsk_treeseq_t *ts)
{
    int ret;
    tsk_size_t size = tsk_treeseq_get_num_mutations(ts) * sizeof(tsk_id_t);
    tsk_id_t *parent = tsk_malloc(size);
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(parent != NULL);
    ret = tsk_treeseq_copy_tables(ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_memcpy(parent, tables.mutations.parent, size);
    /* tsk_table_collection_print_state(&tables, stdout); */
    /* Make sure the tables are actually updated */
    tsk_memset(tables.mutations.parent, 0xff, size);

    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tsk_memcmp(parent, tables.mutations.parent, size), 0);
    /* printf("after\n"); */
    /* tsk_table_collection_print_state(&tables, stdout); */

    free(parent);
    tsk_table_collection_free(&tables);
}

static void
verify_compute_mutation_times(tsk_treeseq_t *ts)
{
    int ret;
    tsk_size_t j;
    tsk_size_t size = tsk_treeseq_get_num_mutations(ts) * sizeof(tsk_id_t);
    tsk_id_t *time = tsk_malloc(size);
    tsk_table_collection_t tables;

    CU_ASSERT_FATAL(time != NULL);
    ret = tsk_treeseq_copy_tables(ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_memcpy(time, tables.mutations.time, size);
    /* Time should be set to TSK_UNKNOWN_TIME before computing */
    for (j = 0; j < size; j++) {
        tables.mutations.time[j] = TSK_UNKNOWN_TIME;
    }

    ret = tsk_table_collection_compute_mutation_times(&tables, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tsk_memcmp(time, tables.mutations.time, size), 0);

    free(time);
    tsk_table_collection_free(&tables);
}

static void
verify_individual_nodes(tsk_treeseq_t *ts)
{
    int ret;
    tsk_individual_t individual;
    tsk_id_t k;
    tsk_size_t num_nodes = tsk_treeseq_get_num_nodes(ts);
    tsk_size_t num_individuals = tsk_treeseq_get_num_individuals(ts);
    tsk_size_t j;

    for (k = 0; k < (tsk_id_t) num_individuals; k++) {
        ret = tsk_treeseq_get_individual(ts, k, &individual);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        for (j = 0; j < individual.nodes_length; j++) {
            CU_ASSERT_FATAL(individual.nodes[j] < (tsk_id_t) num_nodes);
            CU_ASSERT_EQUAL_FATAL(k, ts->tables->nodes.individual[individual.nodes[j]]);
        }
    }
}

static void
verify_trees(tsk_treeseq_t *ts, tsk_size_t num_trees, tsk_id_t *parents)
{
    int ret;
    tsk_id_t u, j, v;
    uint32_t mutation_index, site_index;
    tsk_size_t k, l, tree_sites_length;
    const tsk_site_t *sites = NULL;
    tsk_tree_t tree;
    tsk_size_t num_edges;
    tsk_size_t num_nodes = tsk_treeseq_get_num_nodes(ts);
    tsk_size_t num_sites = tsk_treeseq_get_num_sites(ts);
    tsk_size_t num_mutations = tsk_treeseq_get_num_mutations(ts);
    const double *breakpoints = tsk_treeseq_get_breakpoints(ts);

    ret = tsk_tree_init(&tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(ts), num_trees);

    CU_ASSERT_EQUAL(tree.index, -1);
    site_index = 0;
    mutation_index = 0;
    j = 0;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)) {
        CU_ASSERT_EQUAL(j, (tsk_id_t) tree.index);
        tsk_tree_print_state(&tree, _devnull);
        /* tsk_tree_print_state(&tree, stdout); */
        CU_ASSERT_EQUAL(tree.interval.left, breakpoints[j]);
        num_edges = 0;
        for (u = 0; u < (tsk_id_t) num_nodes; u++) {
            ret = tsk_tree_get_parent(&tree, u, &v);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(v, parents[j * (tsk_id_t) num_nodes + u]);
            if (v != TSK_NULL) {
                num_edges++;
            }
        }
        CU_ASSERT_EQUAL(num_edges, tree.num_edges);
        ret = tsk_tree_get_sites(&tree, &sites, &tree_sites_length);
        CU_ASSERT_EQUAL(ret, 0);
        for (k = 0; k < tree_sites_length; k++) {
            CU_ASSERT_EQUAL(sites[k].id, (tsk_id_t) site_index);
            for (l = 0; l < sites[k].mutations_length; l++) {
                CU_ASSERT_EQUAL(sites[k].mutations[l].id, (tsk_id_t) mutation_index);
                CU_ASSERT_EQUAL(sites[k].mutations[l].site, (tsk_id_t) site_index);
                mutation_index++;
            }
            site_index++;
        }
        j++;
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(site_index, num_sites);
    CU_ASSERT_EQUAL(mutation_index, num_mutations);
    CU_ASSERT_EQUAL(tree.index, -1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(ts), breakpoints[j]);

    tsk_tree_free(&tree);
}

static tsk_tree_t *
get_tree_list(tsk_treeseq_t *ts)
{
    int ret;
    tsk_tree_t t, *trees;
    tsk_size_t num_trees;

    num_trees = tsk_treeseq_get_num_trees(ts);
    ret = tsk_tree_init(&t, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    trees = tsk_malloc(num_trees * sizeof(tsk_tree_t));
    CU_ASSERT_FATAL(trees != NULL);
    for (ret = tsk_tree_first(&t); ret == TSK_TREE_OK; ret = tsk_tree_next(&t)) {
        CU_ASSERT_FATAL(t.index < (tsk_id_t) num_trees);
        ret = tsk_tree_copy(&t, &trees[t.index], 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        check_trees_equal(&trees[t.index], &t);
        /* Make sure the left and right coordinates are also OK */
        CU_ASSERT_EQUAL(trees[t.index].interval.left, t.interval.left);
        CU_ASSERT_EQUAL(trees[t.index].interval.right, t.interval.right);
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
    tsk_id_t j;
    tsk_id_t num_trees = (tsk_id_t) tsk_treeseq_get_num_trees(ts);

    trees = get_tree_list(ts);
    ret = tsk_tree_init(&t, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Single forward pass */
    j = 0;
    for (ret = tsk_tree_first(&t); ret == TSK_TREE_OK; ret = tsk_tree_next(&t)) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        check_trees_equal(&t, &trees[t.index]);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);

    /* Single reverse pass */
    j = num_trees;
    for (ret = tsk_tree_last(&t); ret == TSK_TREE_OK; ret = tsk_tree_prev(&t)) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        check_trees_equal(&t, &trees[t.index]);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);

    /* Full forward, then reverse */
    j = 0;
    for (ret = tsk_tree_first(&t); ret == TSK_TREE_OK; ret = tsk_tree_next(&t)) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        check_trees_equal(&t, &trees[t.index]);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);
    while ((ret = tsk_tree_prev(&t)) == TSK_TREE_OK) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        check_trees_equal(&t, &trees[t.index]);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);
    CU_ASSERT_EQUAL_FATAL(t.index, -1);

    /* Full reverse then forward */
    j = num_trees;
    for (ret = tsk_tree_last(&t); ret == TSK_TREE_OK; ret = tsk_tree_prev(&t)) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        check_trees_equal(&t, &trees[t.index]);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);
    while ((ret = tsk_tree_next(&t)) == TSK_TREE_OK) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        check_trees_equal(&t, &trees[t.index]);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);
    CU_ASSERT_EQUAL_FATAL(t.index, -1);

    /* Do a zigzagging traversal */
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    for (j = 1; j < TSK_MIN(10, num_trees / 2); j++) {
        while (t.index < num_trees - j) {
            ret = tsk_tree_next(&t);
            CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
        }
        CU_ASSERT_EQUAL_FATAL(t.index, num_trees - j);
        check_trees_equal(&t, &trees[t.index]);
        while (t.index > j) {
            ret = tsk_tree_prev(&t);
            CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
        }
        CU_ASSERT_EQUAL_FATAL(t.index, j);
        check_trees_equal(&t, &trees[t.index]);
    }

    ret = tsk_tree_clear(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Calling next() on a cleared tree should be the same as first() */
    j = 0;
    while ((ret = tsk_tree_next(&t)) == TSK_TREE_OK) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        check_trees_equal(&t, &trees[t.index]);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);

    ret = tsk_tree_free(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_init(&t, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Calling prev() on an uninitialised tree should be the same as last() */

    j = num_trees;
    while ((ret = tsk_tree_prev(&t)) == TSK_TREE_OK) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        check_trees_equal(&t, &trees[t.index]);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);

    /* Free the trees. */
    ret = tsk_tree_free(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < (tsk_id_t) tsk_treeseq_get_num_trees(ts); j++) {
        tsk_tree_free(&trees[j]);
    }
    free(trees);
}

static void
verify_tree_diffs(tsk_treeseq_t *ts, tsk_flags_t options)
{
    int ret, valid_tree;
    tsk_diff_iter_t iter;
    tsk_tree_t tree;
    tsk_edge_list_node_t *record;
    tsk_edge_list_t records_out, records_in;
    tsk_size_t num_nodes = tsk_treeseq_get_num_nodes(ts);
    tsk_size_t j, num_trees;
    double lft, rgt;
    tsk_id_t *parent = tsk_malloc(num_nodes * sizeof(tsk_id_t));
    tsk_id_t *child = tsk_malloc(num_nodes * sizeof(tsk_id_t));
    tsk_id_t *sib = tsk_malloc(num_nodes * sizeof(tsk_id_t));

    CU_ASSERT_FATAL(parent != NULL);
    CU_ASSERT_FATAL(child != NULL);
    CU_ASSERT_FATAL(sib != NULL);
    for (j = 0; j < num_nodes; j++) {
        parent[j] = TSK_NULL;
        child[j] = TSK_NULL;
        sib[j] = TSK_NULL;
    }
    ret = tsk_diff_iter_init(&iter, ts, options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_init(&tree, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    valid_tree = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(valid_tree, TSK_TREE_OK);
    tsk_diff_iter_print_state(&iter, _devnull);

    num_trees = 0;
    while ((ret = tsk_diff_iter_next(&iter, &lft, &rgt, &records_out, &records_in))
           == TSK_TREE_OK) {
        tsk_diff_iter_print_state(&iter, _devnull);
        num_trees++;
        /* Update forwards */
        for (record = records_out.head; record != NULL; record = record->next) {
            parent[record->edge.child] = TSK_NULL;
        }
        for (record = records_in.head; record != NULL; record = record->next) {
            parent[record->edge.child] = record->edge.parent;
        }
        if (valid_tree) {
            /* Now check against the sparse tree iterator. */
            for (j = 0; j < num_nodes; j++) {
                CU_ASSERT_EQUAL(parent[j], tree.parent[j]);
            }
        }
        /* Update backwards */
        for (record = records_out.tail; record != NULL; record = record->prev) {
            parent[record->edge.child] = TSK_NULL;
        }
        for (record = records_in.tail; record != NULL; record = record->prev) {
            parent[record->edge.child] = record->edge.parent;
        }
        if (valid_tree) {
            /* Now check against the sparse tree iterator. */
            for (j = 0; j < num_nodes; j++) {
                CU_ASSERT_EQUAL(parent[j], tree.parent[j]);
            }
            CU_ASSERT_EQUAL(tree.interval.left, lft);
            CU_ASSERT_EQUAL(tree.interval.right, rgt);
            valid_tree = tsk_tree_next(&tree);
            if (num_trees < tsk_treeseq_get_num_trees(ts)) {
                CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
            } else {
                CU_ASSERT_EQUAL(valid_tree, 0);
            }
        } else {
            CU_ASSERT_TRUE_FATAL(options & TSK_INCLUDE_TERMINAL);
            for (j = 0; j < num_nodes; j++) {
                CU_ASSERT_EQUAL(parent[j], -1);
            }
            CU_ASSERT_EQUAL(lft, tsk_treeseq_get_sequence_length(ts));
            CU_ASSERT_EQUAL(rgt, tsk_treeseq_get_sequence_length(ts));
        }
    }
    if (options & TSK_INCLUDE_TERMINAL) {
        CU_ASSERT_EQUAL(num_trees, tsk_treeseq_get_num_trees(ts) + 1);
    } else {
        CU_ASSERT_EQUAL(num_trees, tsk_treeseq_get_num_trees(ts));
    }
    CU_ASSERT_EQUAL_FATAL(valid_tree, 0);
    ret = tsk_diff_iter_free(&iter);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_free(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    free(parent);
    free(child);
    free(sib);
}

/* When we keep all sites in simplify, the genotypes for the subset of the
 * samples should be the same as the original */
static void
verify_simplify_genotypes(tsk_treeseq_t *ts, tsk_treeseq_t *subset,
    const tsk_id_t *samples, tsk_size_t num_samples)
{
    int ret;
    tsk_size_t m = tsk_treeseq_get_num_sites(ts);
    tsk_vargen_t vargen, subset_vargen;
    tsk_variant_t *variant, *subset_variant;
    tsk_size_t j, k;
    int32_t a1, a2;
    const tsk_id_t *sample_index_map;

    sample_index_map = tsk_treeseq_get_sample_index_map(ts);

    /* tsk_treeseq_print_state(ts, stdout); */
    /* tsk_treeseq_print_state(subset, stdout); */

    ret = tsk_vargen_init(&vargen, ts, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_vargen_init(
        &subset_vargen, subset, NULL, 0, NULL, TSK_ISOLATED_NOT_MISSING);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(m, tsk_treeseq_get_num_sites(subset));

    for (j = 0; j < m; j++) {
        ret = tsk_vargen_next(&vargen, &variant);
        CU_ASSERT_EQUAL_FATAL(ret, 1);
        ret = tsk_vargen_next(&subset_vargen, &subset_variant);
        CU_ASSERT_EQUAL_FATAL(ret, 1);
        CU_ASSERT_EQUAL(variant->site.id, (tsk_id_t) j)
        CU_ASSERT_EQUAL(subset_variant->site.id, (tsk_id_t) j)
        CU_ASSERT_EQUAL(variant->site.position, subset_variant->site.position);
        for (k = 0; k < num_samples; k++) {
            CU_ASSERT_FATAL(sample_index_map[samples[k]] < (tsk_id_t) ts->num_samples);
            a1 = variant->genotypes[sample_index_map[samples[k]]];
            a2 = subset_variant->genotypes[k];
            /* printf("a1 = %d, a2 = %d\n", a1, a2); */
            /* printf("k = %d original node = %d " */
            /*         "original_index = %d a1=%.*s a2=%.*s\n", */
            /*         (int) k, samples[k], sample_index_map[samples[k]], */
            /*         variant->allele_lengths[a1], variant->alleles[a1], */
            /*         subset_variant->allele_lengths[a2], subset_variant->alleles[a2]);
             */
            CU_ASSERT_FATAL(a1 < (int) variant->num_alleles);
            CU_ASSERT_FATAL(a2 < (int) subset_variant->num_alleles);
            CU_ASSERT_EQUAL_FATAL(
                variant->allele_lengths[a1], subset_variant->allele_lengths[a2]);
            CU_ASSERT_NSTRING_EQUAL_FATAL(variant->alleles[a1],
                subset_variant->alleles[a2], variant->allele_lengths[a1]);
        }
    }
    tsk_vargen_free(&vargen);
    tsk_vargen_free(&subset_vargen);
}

static void
verify_simplify_properties(tsk_treeseq_t *ts, tsk_treeseq_t *subset,
    const tsk_id_t *samples, tsk_size_t num_samples, tsk_id_t *node_map)
{
    int ret;
    tsk_node_t n1, n2;
    tsk_tree_t full_tree, subset_tree;
    const tsk_site_t *tree_sites;
    tsk_size_t tree_sites_length;
    uint32_t j, k;
    tsk_id_t u, mrca1, mrca2;
    tsk_size_t total_sites;

    CU_ASSERT_EQUAL(
        tsk_treeseq_get_sequence_length(ts), tsk_treeseq_get_sequence_length(subset));
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(subset), num_samples);
    CU_ASSERT(tsk_treeseq_get_num_nodes(ts) >= tsk_treeseq_get_num_nodes(subset));
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(subset), num_samples);

    /* Check the sample properties */
    for (j = 0; j < num_samples; j++) {
        ret = tsk_treeseq_get_node(ts, samples[j], &n1);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(node_map[samples[j]], (tsk_id_t) j);
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
        ret = tsk_treeseq_get_node(ts, (tsk_id_t) j, &n1);
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
    ret = tsk_tree_init(&full_tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_init(&subset_tree, subset, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&full_tree);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    ret = tsk_tree_first(&subset_tree);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);

    total_sites = 0;
    while (1) {
        while (full_tree.interval.right <= subset_tree.interval.right) {
            for (j = 0; j < num_samples; j++) {
                for (k = j + 1; k < num_samples; k++) {
                    ret = tsk_tree_get_mrca(&full_tree, samples[j], samples[k], &mrca1);
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                    ret = tsk_tree_get_mrca(&subset_tree, node_map[samples[j]],
                        node_map[samples[k]], &mrca2);
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
            CU_ASSERT(subset_tree.interval.left <= tree_sites[j].position);
            CU_ASSERT(tree_sites[j].position < subset_tree.interval.right);
            for (k = 0; k < tree_sites[j].mutations_length; k++) {
                ret = tsk_tree_get_parent(
                    &subset_tree, tree_sites[j].mutations[k].node, &u);
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
}

static void
verify_simplify(tsk_treeseq_t *ts)
{
    int ret;
    tsk_size_t n = tsk_treeseq_get_num_samples(ts);
    tsk_size_t num_samples[] = { 0, 1, 2, 3, n / 2, n - 1, n };
    tsk_size_t j;
    const tsk_id_t *sample;
    tsk_id_t *node_map = tsk_malloc(tsk_treeseq_get_num_nodes(ts) * sizeof(tsk_id_t));
    tsk_treeseq_t subset;
    tsk_flags_t options = TSK_SIMPLIFY_FILTER_SITES;

    CU_ASSERT_FATAL(node_map != NULL);
    sample = tsk_treeseq_get_samples(ts);
    if (tsk_treeseq_get_num_migrations(ts) > 0) {
        ret = tsk_treeseq_simplify(ts, sample, 2, 0, &subset, NULL);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SIMPLIFY_MIGRATIONS_NOT_SUPPORTED);
        /* Exiting early here because simplify isn't supported with migrations. */
        goto out;
    }

    for (j = 0; j < sizeof(num_samples) / sizeof(*num_samples); j++) {
        if (num_samples[j] <= n) {
            ret = tsk_treeseq_simplify(
                ts, sample, num_samples[j], options, &subset, node_map);
            /* printf("ret = %s\n", tsk_strerror(ret)); */
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            verify_simplify_properties(ts, &subset, sample, num_samples[j], node_map);
            tsk_treeseq_free(&subset);

            /* Keep all sites */
            ret = tsk_treeseq_simplify(ts, sample, num_samples[j], 0, &subset, node_map);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            verify_simplify_properties(ts, &subset, sample, num_samples[j], node_map);
            verify_simplify_genotypes(ts, &subset, sample, num_samples[j]);
            tsk_treeseq_free(&subset);
        }
    }
out:
    free(node_map);
}

typedef struct {
    tsk_id_t tree_index;
    tsk_id_t node;
    tsk_size_t count;
} sample_count_test_t;

static void
verify_sample_counts(tsk_treeseq_t *ts, tsk_size_t num_tests, sample_count_test_t *tests)
{
    int ret;
    tsk_size_t j, num_samples, n, k;
    tsk_id_t stop, sample_index;
    tsk_tree_t tree;
    const tsk_id_t *samples;

    n = tsk_treeseq_get_num_samples(ts);
    samples = tsk_treeseq_get_samples(ts);

    /* First run with the TSK_NO_SAMPLE_COUNTS feature */

    ret = tsk_tree_init(&tree, ts, TSK_NO_SAMPLE_COUNTS);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = tsk_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
        }
        ret = tsk_tree_get_num_samples(&tree, tests[j].node, &num_samples);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_samples);
        /* all operations depending on tracked samples should fail. */
        ret = tsk_tree_get_num_tracked_samples(&tree, 0, &num_samples);
        CU_ASSERT_EQUAL(ret, TSK_ERR_UNSUPPORTED_OPERATION);
        /* The root should be NULL */
        CU_ASSERT_EQUAL(tsk_tree_get_left_root(&tree), TSK_NULL);
    }
    tsk_tree_free(&tree);

    /* Now run with TSK_SAMPLE_COUNTS but with no samples tracked. */
    ret = tsk_tree_init(&tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = tsk_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
        }
        ret = tsk_tree_get_num_samples(&tree, tests[j].node, &num_samples);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_samples);
        /* all operations depending on tracked samples should fail. */
        ret = tsk_tree_get_num_tracked_samples(&tree, 0, &num_samples);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(num_samples, 0);
        /* The root should not be NULL */
        CU_ASSERT_NOT_EQUAL(tree.virtual_root, TSK_NULL);
    }
    tsk_tree_free(&tree);

    /* Run with TSK_SAMPLE_LISTS and TSK_NO_SAMPLE_COUNTS */
    ret = tsk_tree_init(&tree, ts, TSK_SAMPLE_LISTS | TSK_NO_SAMPLE_COUNTS);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = tsk_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
        }
        ret = tsk_tree_get_num_samples(&tree, tests[j].node, &num_samples);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_samples);
        /* all operations depending on tracked samples should fail. */
        ret = tsk_tree_get_num_tracked_samples(&tree, 0, &num_samples);
        CU_ASSERT_EQUAL(ret, TSK_ERR_UNSUPPORTED_OPERATION);

        sample_index = tree.left_sample[tests[j].node];
        k = 0;
        if (sample_index != TSK_NULL) {
            stop = tree.right_sample[tests[j].node];
            while (true) {
                k++;
                CU_ASSERT_FATAL(k <= tests[j].count);
                if (sample_index == stop) {
                    break;
                }
                sample_index = tree.next_sample[sample_index];
            }
        }
        CU_ASSERT_EQUAL(tests[j].count, k);
    }
    tsk_tree_free(&tree);

    /* Now use TSK_SAMPLE_LISTS */
    ret = tsk_tree_init(&tree, ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_set_tracked_samples(&tree, n, samples);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = tsk_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
        }
        ret = tsk_tree_get_num_samples(&tree, tests[j].node, &num_samples);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_samples);

        /* We're tracking all samples, so the count should be the same */
        ret = tsk_tree_get_num_tracked_samples(&tree, tests[j].node, &num_samples);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_samples);

        sample_index = tree.left_sample[tests[j].node];
        k = 0;
        if (sample_index != TSK_NULL) {
            stop = tree.right_sample[tests[j].node];
            while (true) {
                k++;
                if (sample_index == stop) {
                    break;
                }
                sample_index = tree.next_sample[sample_index];
            }
        }
        CU_ASSERT_EQUAL(tests[j].count, k);
    }
    tsk_tree_free(&tree);
}

static void
verify_sample_sets_for_tree(tsk_tree_t *tree)
{
    int ret, stack_top, j;
    tsk_id_t u, v;
    tsk_size_t tmp, n, num_nodes, num_samples;
    tsk_id_t *stack, *samples;
    const tsk_treeseq_t *ts = tree->tree_sequence;
    tsk_id_t *sample_index_map = ts->sample_index_map;
    const tsk_id_t *list_left = tree->left_sample;
    const tsk_id_t *list_right = tree->right_sample;
    const tsk_id_t *list_next = tree->next_sample;
    tsk_id_t stop, sample_index;

    n = tsk_treeseq_get_num_samples(ts);
    num_nodes = tsk_treeseq_get_num_nodes(ts);
    stack = tsk_malloc(n * sizeof(tsk_id_t));
    samples = tsk_malloc(n * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(stack != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    for (u = 0; u < (tsk_id_t) num_nodes; u++) {
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
                    CU_ASSERT_TRUE_FATAL(j < (tsk_id_t) n);
                    CU_ASSERT_EQUAL_FATAL(sample_index, sample_index_map[samples[j]]);
                    j++;
                    if (sample_index == stop) {
                        break;
                    }
                    sample_index = list_next[sample_index];
                }
            }
            CU_ASSERT_EQUAL_FATAL(j, (int) num_samples);
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

    ret = tsk_tree_init(&t, ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL(ret, 0);

    for (ret = tsk_tree_first(&t); ret == TSK_TREE_OK; ret = tsk_tree_next(&t)) {
        verify_sample_sets_for_tree(&t);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (ret = tsk_tree_last(&t); ret == TSK_TREE_OK; ret = tsk_tree_prev(&t)) {
        verify_sample_sets_for_tree(&t);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_tree_free(&t);
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
}

/*=======================================================
 * Simplest test cases.
 *======================================================*/

static void
test_simplest_discrete_genome(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0";
    const char *edges = "0  1   2   0,1\n";
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_id_t ret_id;
    int ret;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    CU_ASSERT_TRUE(tsk_treeseq_get_discrete_genome(&ts));

    ret = tsk_table_collection_copy(ts.tables, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);

    tables.sequence_length = 1.001;
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_treeseq_get_discrete_genome(&ts));
    tsk_treeseq_free(&ts);
    tables.sequence_length = 1;

    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_treeseq_get_discrete_genome(&ts));
    tsk_treeseq_free(&ts);

    tables.edges.right[0] = 0.999;
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_treeseq_get_discrete_genome(&ts));
    tsk_treeseq_free(&ts);
    tables.edges.right[0] = 1.0;

    tables.edges.left[0] = 0.999;
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_treeseq_get_discrete_genome(&ts));
    tsk_treeseq_free(&ts);
    tables.edges.left[0] = 0;

    ret_id = tsk_site_table_add_row(&tables.sites, 0, "A", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_treeseq_get_discrete_genome(&ts));
    tsk_treeseq_free(&ts);

    tables.sites.position[0] = 0.001;
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_treeseq_get_discrete_genome(&ts));
    tsk_treeseq_free(&ts);
    tables.sites.position[0] = 0;

    /* Need another population for a migration */
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);

    ret_id
        = tsk_migration_table_add_row(&tables.migrations, 0, 1, 0, 0, 1, 1.0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_treeseq_get_discrete_genome(&ts));
    tsk_treeseq_free(&ts);

    tables.migrations.left[0] = 0.001;
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_treeseq_get_discrete_genome(&ts));
    tsk_treeseq_free(&ts);
    tables.migrations.left[0] = 0;

    tables.migrations.right[0] = 0.999;
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_treeseq_get_discrete_genome(&ts));
    tsk_treeseq_free(&ts);
    tables.migrations.right[0] = 1;

    /* An empty tree sequence is has a discrete genome. */
    tsk_table_collection_clear(&tables, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_treeseq_get_discrete_genome(&ts));
    tsk_treeseq_free(&ts);

    tsk_table_collection_free(&tables);
}

static void
test_simplest_discrete_time(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  0   0\n"
                        "0  0   0";
    const char *edges = "0  1   2   0,1,3,4\n";
    const char *sites = "0.1  0\n"
                        "0.2  0\n"
                        "0.3  0\n"
                        "0.4  0\n";
    const char *mutations = "0    0     1\n"
                            "1    1     1\n"
                            "2    3     1\n"
                            "3    4     1";
    const char *migrations = "0  1  0  0  1  1";

    tsk_treeseq_from_text(
        &ts, 1, nodes, edges, migrations, sites, mutations, NULL, NULL, 0);
    CU_ASSERT_TRUE(tsk_treeseq_get_discrete_time(&ts));

    ret = tsk_table_collection_copy(ts.tables, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);

    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_treeseq_get_discrete_time(&ts));
    tsk_treeseq_free(&ts);

    tables.nodes.time[0] = 0.0001;
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_treeseq_get_discrete_time(&ts));
    tsk_treeseq_free(&ts);
    tables.nodes.time[0] = 0;

    tables.mutations.time[0] = 0.001;
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_treeseq_get_discrete_time(&ts));
    tsk_treeseq_free(&ts);
    tables.mutations.time[0] = 0;

    tables.migrations.time[0] = 0.001;
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_treeseq_get_discrete_time(&ts));
    tsk_treeseq_free(&ts);
    tables.migrations.time[0] = 0;

    tables.mutations.time[0] = TSK_UNKNOWN_TIME;
    tables.mutations.time[1] = TSK_UNKNOWN_TIME;
    tables.mutations.time[2] = TSK_UNKNOWN_TIME;
    tables.mutations.time[3] = TSK_UNKNOWN_TIME;
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_treeseq_get_discrete_time(&ts));
    tsk_treeseq_free(&ts);

    /* An empty tree sequence is has a discrete time. */
    tsk_table_collection_clear(&tables, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_treeseq_get_discrete_time(&ts));
    tsk_treeseq_free(&ts);

    tsk_table_collection_free(&tables);
}

static void
test_simplest_records(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0";
    const char *edges = "0  1   2   0,1\n";
    tsk_treeseq_t ts, simplified;
    tsk_id_t sample_ids[] = { 0, 1 };
    int ret;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 2, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, simplified.tables, 0));
    tsk_treeseq_free(&simplified);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 2,
        TSK_SIMPLIFY_KEEP_UNARY | TSK_SIMPLIFY_KEEP_UNARY_IN_INDIVIDUALS, &simplified,
        NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_KEEP_UNARY_MUTUALLY_EXCLUSIVE);
    tsk_treeseq_free(&simplified);

    ret = tsk_treeseq_simplify(
        &ts, sample_ids, 2, TSK_SIMPLIFY_KEEP_UNARY, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, simplified.tables, 0));
    tsk_treeseq_free(&simplified);

    ret = tsk_treeseq_simplify(
        &ts, sample_ids, 2, TSK_SIMPLIFY_KEEP_UNARY_IN_INDIVIDUALS, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, simplified.tables, 0));
    tsk_treeseq_free(&simplified);

    tsk_treeseq_free(&ts);
}

static void
test_simplest_nonbinary_records(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0";
    const char *edges = "0  1   4   0,1,2,3\n";
    tsk_treeseq_t ts, simplified;
    tsk_id_t sample_ids[] = { 0, 1, 2, 3 };
    int ret;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 5);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 4, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, simplified.tables, 0));
    tsk_treeseq_free(&simplified);

    ret = tsk_treeseq_simplify(
        &ts, sample_ids, 4, TSK_SIMPLIFY_KEEP_UNARY, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, simplified.tables, 0));
    tsk_treeseq_free(&simplified);

    ret = tsk_treeseq_simplify(
        &ts, sample_ids, 4, TSK_SIMPLIFY_KEEP_UNARY_IN_INDIVIDUALS, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, simplified.tables, 0));
    tsk_treeseq_free(&simplified);

    tsk_treeseq_free(&ts);
}

static void
test_simplest_unary_records(void)
{
    int ret;
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  1   0\n"
                        "0  2   0";
    const char *edges = "0  1   2   0\n"
                        "0  1   3   1\n"
                        "0  1   4   2,3\n";
    tsk_treeseq_t ts, simplified, simplified_other;
    tsk_id_t sample_ids[] = { 0, 1 };

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 5);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_populations(&ts), 1);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 2, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&simplified), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&simplified), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_edges(&simplified), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&simplified), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&simplified), 1);
    tsk_treeseq_free(&simplified);

    ret = tsk_treeseq_simplify(
        &ts, sample_ids, 2, TSK_SIMPLIFY_KEEP_UNARY, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, simplified.tables, 0));
    tsk_treeseq_free(&simplified);

    ret = tsk_treeseq_simplify(
        &ts, sample_ids, 2, TSK_SIMPLIFY_KEEP_UNARY_IN_INDIVIDUALS, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_equals(ts.tables, simplified.tables, 0));
    ret = tsk_treeseq_simplify(&ts, sample_ids, 2, 0, &simplified_other, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(
        tsk_table_collection_equals(simplified.tables, simplified_other.tables, 0));
    tsk_treeseq_free(&simplified);
    tsk_treeseq_free(&simplified_other);

    tsk_treeseq_free(&ts);
}

static void
test_simplest_unary_with_individuals(void)
{
    int ret;
    const char *nodes = "1  0   0   -1\n"
                        "1  0   0   0\n"
                        "0  1   0   -1\n"
                        "0  1   0   1\n"
                        "0  2   0   -1\n"
                        "0  3   0   -1\n"
                        "0  3   0   2\n"
                        "0  1   0   -1\n"
                        "0  1   0   3\n"
                        "0  0   0   -1\n"
                        "0  0   0   4\n"
                        "0  1   0   3\n";
    const char *edges = "0  2   2   0\n"
                        "0  2   3   1\n"
                        "2  3   7   0\n"
                        "2  3   8   1,9\n"
                        "2  3   11   10\n"
                        "0  2   4   2,3\n"
                        "0  1   5   4\n"
                        "1  2   6   4\n";
    const char *individuals = "0    0.5     -1,-1\n"
                              "0    1.5,3.1 -1,-1\n"
                              "0    2.1     0,1\n"
                              "0    3.2     1,2\n"
                              "0    4.2     2,3\n";
    const char *nodes_expect = "1  0   0   -1\n"
                               "1  0   0   0\n"
                               "0  1   0   1\n"
                               "0  1   0   3\n"
                               "0  2   0   -1\n"
                               "0  3   0   2\n";
    const char *edges_expect = "0  2   2   1\n"
                               "2  3   3   1\n"
                               "0  2   4   0,2\n"
                               "1  2   5   4\n";
    const char *individuals_expect = "0    0.5     -1,-1\n"
                                     "0    1.5,3.1 -1,-1\n"
                                     "0    2.1     0,1\n"
                                     "0    3.2     1,2\n";
    tsk_treeseq_t ts, simplified, expected;
    tsk_id_t sample_ids[] = { 0, 1 };

    tsk_treeseq_from_text(&ts, 3, nodes, edges, NULL, NULL, NULL, individuals, NULL, 0);
    tsk_treeseq_from_text(&expected, 3, nodes_expect, edges_expect, NULL, NULL, NULL,
        individuals_expect, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 3.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 12);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_individuals(&ts), 5);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_populations(&ts), 1);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 2,
        TSK_SIMPLIFY_KEEP_UNARY_IN_INDIVIDUALS | TSK_SIMPLIFY_FILTER_INDIVIDUALS,
        &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(simplified.tables, expected.tables, 0));
    tsk_treeseq_free(&simplified);

    tsk_treeseq_free(&expected);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_non_sample_leaf_records(void)
{
    int ret;
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  0   0\n"
                        "0  0   0";
    const char *edges = "0  1   2   0,1,3,4\n";
    const char *sites = "0.1  0\n"
                        "0.2  0\n"
                        "0.3  0\n"
                        "0.4  0\n";
    const char *mutations = "0    0     1\n"
                            "1    1     1\n"
                            "2    3     1\n"
                            "3    4     1";
    tsk_treeseq_t ts, simplified;
    tsk_id_t sample_ids[] = { 0, 1 };
    tsk_vargen_t vargen;
    tsk_variant_t *var;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, sites, mutations, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 5);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_vargen_init(&vargen, &ts, NULL, 0, NULL, 0);
    tsk_vargen_print_state(&vargen, _devnull);
    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->genotypes[0], 1);
    CU_ASSERT_EQUAL(var->genotypes[1], 0);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->genotypes[0], 0);
    CU_ASSERT_EQUAL(var->genotypes[1], 1);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->genotypes[0], 0);
    CU_ASSERT_EQUAL(var->genotypes[1], 0);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->genotypes[0], 0);
    CU_ASSERT_EQUAL(var->genotypes[1], 0);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_vargen_free(&vargen);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 2, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&simplified), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&simplified), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&simplified), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&simplified), 1);

    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&simplified);
}

static void
test_simplest_degenerate_multiple_root_records(void)
{

    int ret;
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  1   0\n";
    const char *edges = "0  1   2   0\n"
                        "0  1   3   1\n";
    tsk_treeseq_t ts, simplified;
    tsk_tree_t t;
    tsk_id_t sample_ids[] = { 0, 1 };

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&t), 2);
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&t), 2);
    CU_ASSERT_EQUAL(tsk_tree_get_right_root(&t), 3);
    CU_ASSERT_EQUAL(t.num_edges, 2);
    CU_ASSERT_EQUAL(t.right_sib[2], 3);
    CU_ASSERT_EQUAL(t.right_sib[3], TSK_NULL);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 2, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&simplified), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&simplified), 2);
    tsk_treeseq_free(&simplified);

    ret = tsk_treeseq_simplify(
        &ts, sample_ids, 2, TSK_SIMPLIFY_KEEP_UNARY, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, simplified.tables, 0));
    tsk_treeseq_free(&simplified);

    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_simplest_multiple_root_records(void)
{
    int ret;
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  1   0\n";
    const char *edges = "0  1   4   0,1\n"
                        "0  1   5   2,3\n";
    tsk_treeseq_t ts, simplified;
    tsk_id_t sample_ids[] = { 0, 1, 2, 3 };

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 6);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 4, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&simplified), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&simplified), 6);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&simplified), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&simplified), 1);
    tsk_treeseq_free(&simplified);

    /* Make one tree degenerate */
    ret = tsk_treeseq_simplify(&ts, sample_ids, 3, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&simplified), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&simplified), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&simplified), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&simplified), 1);
    tsk_treeseq_free(&simplified);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_zero_root_tree(void)
{
    int ret;
    const char *nodes = "0  0   0\n"
                        "0  0   0\n"
                        "0  0   0\n"
                        "0  0   0\n"
                        "0  1   0\n"
                        "0  1   0\n";
    const char *edges = "0  1   4   0,1\n"
                        "0  1   5   2,3\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 6);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&t), 0);
    CU_ASSERT_EQUAL(t.num_edges, 4);
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&t), TSK_NULL);
    CU_ASSERT_EQUAL(tsk_tree_get_right_root(&t), TSK_NULL);
    CU_ASSERT_EQUAL(t.right_sib[2], 3);
    CU_ASSERT_EQUAL(t.right_sib[3], TSK_NULL);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_multi_root_tree(void)
{
    int ret;
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n";
    const char *edges = "0  1   3   1,2\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_tree_init(&t, &ts, 0);

    tsk_tree_print_state(&t, _devnull);

    /* Make sure the initial roots are set correctly */
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&t), 0);
    CU_ASSERT_EQUAL(t.left_sib[0], TSK_NULL);
    CU_ASSERT_EQUAL(t.right_sib[0], 1);
    CU_ASSERT_EQUAL(t.left_sib[1], 0);
    CU_ASSERT_EQUAL(t.right_sib[1], 2);
    CU_ASSERT_EQUAL(t.left_sib[2], 1);
    CU_ASSERT_EQUAL(t.right_sib[2], TSK_NULL);

    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&t), 2);
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&t), 0);
    CU_ASSERT_EQUAL(t.right_sib[0], 3);
    CU_ASSERT_EQUAL(t.num_edges, 2);

    tsk_tree_print_state(&t, _devnull);

    CU_ASSERT_EQUAL(tsk_tree_set_root_threshold(&t, 1), TSK_ERR_UNSUPPORTED_OPERATION);
    ret = tsk_tree_next(&t);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_tree_set_root_threshold(&t, 0), TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_tree_set_root_threshold(&t, 2);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_tree_get_root_threshold(&t), 2);

    ret = tsk_tree_next(&t);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&t), 1);
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&t), 3);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_tree_mrca(void)
{
    int ret;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    tsk_tree_t t;
    tsk_id_t mrca, ret_id;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_tree_get_mrca(&t, 0, 0, &mrca);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(mrca, 0);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_simplest_root_mutations(void)
{
    int ret;
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n";
    const char *edges = "0  1   2   0,1\n";
    const char *sites = "0.1 0";
    const char *mutations = "0    2     1";
    tsk_flags_t options = 0;
    tsk_id_t sample_ids[] = { 0, 1 };
    tsk_treeseq_t ts, simplified;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, sites, mutations, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 2, options, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&simplified), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&simplified), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&simplified), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&simplified), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&simplified), 1);
    tsk_treeseq_free(&simplified);

    tsk_treeseq_free(&ts);
}

static void
test_simplest_back_mutations(void)
{
    int ret;
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  2   0\n";
    const char *edges = "0  1   3   0,1\n"
                        "0  1   4   2,3\n";
    const char *sites = "0.5 0";
    const char *mutations = "0    3     1   -1\n"
                            "0    0     0   0";
    tsk_treeseq_t ts;
    tsk_vargen_t vargen;
    tsk_variant_t *var;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, sites, mutations, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 5);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_vargen_init(&vargen, &ts, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->genotypes[0], 0);
    CU_ASSERT_EQUAL(var->genotypes[1], 1);
    CU_ASSERT_EQUAL(var->genotypes[2], 0);
    CU_ASSERT_EQUAL(var->site.id, 0);
    CU_ASSERT_EQUAL(var->site.mutations_length, 2);
    tsk_vargen_free(&vargen);

    tsk_treeseq_free(&ts);
}

static void
test_simplest_general_samples(void)
{
    const char *nodes = "1  0   0\n"
                        "0  1   0\n"
                        "1  0   0";
    const char *edges = "0  1   1   0,2\n";
    const char *sites = "0.5  0\n"
                        "0.75 0\n";
    const char *mutations = "0    2     1\n"
                            "1    0     1";
    const tsk_id_t samples[2] = { 0, 2 };
    const tsk_id_t *s;
    int ret;

    tsk_treeseq_t ts, simplified;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, sites, mutations, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    s = tsk_treeseq_get_samples(&ts);
    CU_ASSERT_FATAL(s != NULL);
    CU_ASSERT_EQUAL(s[0], 0);
    CU_ASSERT_EQUAL(s[1], 2);

    ret = tsk_treeseq_simplify(&ts, samples, 2, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    s = tsk_treeseq_get_samples(&simplified);
    CU_ASSERT_FATAL(s != NULL);
    CU_ASSERT_EQUAL(s[0], 0);
    CU_ASSERT_EQUAL(s[1], 1);

    tsk_treeseq_free(&simplified);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_holey_tree_sequence(void)
{
    const char *nodes_txt = "1  0   0\n"
                            "1  0   0\n"
                            "0  1   0";
    const char *edges_txt = "0  1   2   0\n"
                            "2  3   2   0\n"
                            "0  1   2   1\n"
                            "2  3   2   1\n";
    const char *sites_txt = "0.5  0\n"
                            "1.5  0\n"
                            "2.5  0\n";
    const char *mutations_txt = "0    0     1\n"
                                "1    1     1\n"
                                "2    2     1\n";
    int ret;
    tsk_treeseq_t ts, simplified;
    tsk_id_t sample_ids[] = { 0, 1 };

    tsk_treeseq_from_text(
        &ts, 3, nodes_txt, edges_txt, NULL, sites_txt, mutations_txt, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 3.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 3);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 2, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, simplified.tables, 0));
    tsk_treeseq_free(&simplified);

    ret = tsk_treeseq_simplify(
        &ts, sample_ids, 2, TSK_SIMPLIFY_KEEP_UNARY, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, simplified.tables, 0));
    tsk_treeseq_free(&simplified);

    tsk_treeseq_free(&ts);
}

static void
test_simplest_holey_tsk_treeseq_mutation_parents(void)
{
    const char *nodes_txt = "1  0   0\n"
                            "1  0   0\n"
                            "0  1   0";
    const char *edges_txt = "0  1   2   0\n"
                            "2  3   2   0\n"
                            "0  1   2   1\n"
                            "2  3   2   1\n";
    const char *sites_txt = "0.5  0\n"
                            "1.5  0\n"
                            "2.5  0\n";
    const char *mutations_txt = "0    0     1\n"
                                "0    0     1\n"
                                "1    1     1\n"
                                "1    1     1\n"
                                "2    2     1\n"
                                "2    2     1\n";
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    int ret;

    tsk_treeseq_from_text(
        &ts, 3, nodes_txt, edges_txt, NULL, sites_txt, mutations_txt, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 3);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.mutations.parent[0], -1);
    CU_ASSERT_EQUAL(tables.mutations.parent[1], 0);
    CU_ASSERT_EQUAL(tables.mutations.parent[2], -1);
    CU_ASSERT_EQUAL(tables.mutations.parent[3], 2);
    CU_ASSERT_EQUAL(tables.mutations.parent[4], -1);
    CU_ASSERT_EQUAL(tables.mutations.parent[5], 4);
    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_initial_gap_tree_sequence(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0";
    const char *edges = "2  3   2   0,1\n";
    const char *sites = "0.5  0\n"
                        "1.5  0\n"
                        "2.5  0\n";
    const char *mutations = "0    0     1\n"
                            "1    1     1\n"
                            "2    2     1";
    int ret;
    tsk_treeseq_t ts, simplified;
    const tsk_id_t z = TSK_NULL;
    tsk_id_t parents[] = {
        z,
        z,
        z,
        2,
        2,
        z,
    };
    tsk_size_t num_trees = 2;
    tsk_id_t sample_ids[] = { 0, 1 };

    tsk_treeseq_from_text(&ts, 3, nodes, edges, NULL, sites, mutations, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 3.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 2);

    verify_trees(&ts, num_trees, parents);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 2, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, simplified.tables, 0));
    tsk_treeseq_free(&simplified);

    ret = tsk_treeseq_simplify(
        &ts, sample_ids, 2, TSK_SIMPLIFY_KEEP_UNARY, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, simplified.tables, 0));
    tsk_treeseq_free(&simplified);

    tsk_treeseq_free(&ts);
}

static void
test_simplest_initial_gap_zero_roots(void)
{
    const char *nodes = "0  0   0\n"
                        "0  0   0\n"
                        "0  1   0";
    const char *edges = "2  3   2   0,1\n";
    int ret;
    tsk_treeseq_t ts;
    const tsk_id_t z = TSK_NULL;
    tsk_id_t parents[] = {
        z,
        z,
        z,
        2,
        2,
        z,
    };
    uint32_t num_trees = 2;
    tsk_tree_t tree;

    tsk_treeseq_from_text(&ts, 3, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 3.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 2);

    verify_trees(&ts, num_trees, parents);

    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&tree), TSK_NULL);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&tree), 0);
    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&tree), TSK_NULL);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&tree), 0);
    CU_ASSERT_EQUAL(tree.parent[0], 2);
    CU_ASSERT_EQUAL(tree.parent[1], 2);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_holey_tsk_treeseq_zero_roots(void)
{
    const char *nodes_txt = "0  0   0\n"
                            "0  0   0\n"
                            "0  1   0";
    const char *edges_txt = "0  1   2   0\n"
                            "2  3   2   0\n"
                            "0  1   2   1\n"
                            "2  3   2   1\n";
    int ret;
    tsk_treeseq_t ts;
    const tsk_id_t z = TSK_NULL;
    tsk_id_t parents[] = {
        2,
        2,
        z,
        z,
        z,
        z,
        2,
        2,
        z,
    };
    uint32_t num_trees = 3;
    tsk_tree_t tree;

    tsk_treeseq_from_text(&ts, 3, nodes_txt, edges_txt, NULL, NULL, NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 3.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 3);

    verify_trees(&ts, num_trees, parents);

    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&tree), TSK_NULL);
    CU_ASSERT_EQUAL(tree.parent[0], 2);
    CU_ASSERT_EQUAL(tree.parent[1], 2);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&tree), 0);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&tree), TSK_NULL);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&tree), 0);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&tree), TSK_NULL);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&tree), 0);
    CU_ASSERT_EQUAL(tree.parent[0], 2);
    CU_ASSERT_EQUAL(tree.parent[1], 2);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_initial_gap_tsk_treeseq_mutation_parents(void)
{
    const char *nodes_txt = "1  0   0\n"
                            "1  0   0\n"
                            "0  1   0";
    const char *edges_txt = "2  3   2   0,1\n";
    const char *sites_txt = "0.5  0\n"
                            "1.5  0\n"
                            "2.5  0\n";
    const char *mutations_txt = "0    0     1\n"
                                "0    0     1\n"
                                "1    1     1\n"
                                "1    1     1\n"
                                "2    2     1\n"
                                "2    2     1\n";
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    int ret;

    tsk_treeseq_from_text(
        &ts, 3, nodes_txt, edges_txt, NULL, sites_txt, mutations_txt, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 2);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.mutations.parent[0], -1);
    CU_ASSERT_EQUAL(tables.mutations.parent[1], 0);
    CU_ASSERT_EQUAL(tables.mutations.parent[2], -1);
    CU_ASSERT_EQUAL(tables.mutations.parent[3], 2);
    CU_ASSERT_EQUAL(tables.mutations.parent[4], -1);
    CU_ASSERT_EQUAL(tables.mutations.parent[5], 4);
    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_final_gap_tree_sequence(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0";
    const char *edges = "0  2   2   0,1\n";
    const char *sites = "0.5  0\n"
                        "1.5  0\n"
                        "2.5  0\n";
    const char *mutations = "0    0     1\n"
                            "1    1     1\n"
                            "2    0     1";
    tsk_treeseq_t ts;
    const tsk_id_t z = TSK_NULL;
    tsk_id_t parents[] = {
        2,
        2,
        z,
        z,
        z,
        z,
    };
    uint32_t num_trees = 2;

    tsk_treeseq_from_text(&ts, 3, nodes, edges, NULL, sites, mutations, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 3.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 2);

    verify_trees(&ts, num_trees, parents);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_final_gap_tsk_treeseq_mutation_parents(void)
{
    const char *nodes_txt = "1  0   0\n"
                            "1  0   0\n"
                            "0  1   0";
    const char *edges_txt = "0  2   2   0,1\n";
    const char *sites_txt = "0.5  0\n"
                            "1.5  0\n"
                            "2.5  0\n";
    const char *mutations_txt = "0    0     1\n"
                                "0    0     1\n"
                                "1    1     1\n"
                                "1    1     1\n"
                                "2    0     1\n"
                                "2    0     1\n";
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    int ret;

    tsk_treeseq_from_text(
        &ts, 3, nodes_txt, edges_txt, NULL, sites_txt, mutations_txt, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 2);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.mutations.parent[0], -1);
    CU_ASSERT_EQUAL(tables.mutations.parent[1], 0);
    CU_ASSERT_EQUAL(tables.mutations.parent[2], -1);
    CU_ASSERT_EQUAL(tables.mutations.parent[3], 2);
    CU_ASSERT_EQUAL(tables.mutations.parent[4], -1);
    CU_ASSERT_EQUAL(tables.mutations.parent[5], 4);
    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_individuals(void)
{
    const char *individuals = "1      0.25     -1,-1\n"
                              "2      0.5,0.25 -1,-1\n"
                              "3      0.75     0,1\n";
    const char *nodes = "1  0   -1  -1\n"
                        "1  0   -1  1\n"
                        "0  0   -1  -1\n"
                        "1  0   -1  0\n"
                        "0  0   -1  1\n"
                        "0  0   -1  2\n";
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    tsk_node_t node;
    tsk_individual_t individual;
    tsk_flags_t load_flags = TSK_TS_INIT_BUILD_INDEXES;
    int ret;
    tsk_id_t pat_id, mat_id;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1.0;
    parse_individuals(individuals, &tables.individuals);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.num_rows, 3);

    parse_nodes(nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 6);

    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_get_node(&ts, 0, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(node.individual, TSK_NULL);

    ret = tsk_treeseq_get_node(&ts, 1, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(node.individual, 1);

    ret = tsk_treeseq_get_individual(&ts, 0, &individual);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(individual.id, 0);
    CU_ASSERT_EQUAL_FATAL(individual.flags, 1);
    CU_ASSERT_EQUAL_FATAL(individual.location_length, 1);
    CU_ASSERT_EQUAL_FATAL(individual.location[0], 0.25);
    CU_ASSERT_EQUAL_FATAL(individual.parents_length, 2);
    CU_ASSERT_EQUAL_FATAL(individual.parents[0], -1);
    CU_ASSERT_EQUAL_FATAL(individual.parents[1], -1);
    pat_id = individual.id;
    CU_ASSERT_EQUAL_FATAL(individual.nodes_length, 1);
    CU_ASSERT_EQUAL_FATAL(individual.nodes[0], 3);

    ret = tsk_treeseq_get_individual(&ts, 1, &individual);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(individual.id, 1);
    CU_ASSERT_EQUAL_FATAL(individual.flags, 2);
    CU_ASSERT_EQUAL_FATAL(individual.location_length, 2);
    CU_ASSERT_EQUAL_FATAL(individual.location[0], 0.5);
    CU_ASSERT_EQUAL_FATAL(individual.location[1], 0.25);
    CU_ASSERT_EQUAL_FATAL(individual.parents_length, 2);
    CU_ASSERT_EQUAL_FATAL(individual.parents[0], -1);
    CU_ASSERT_EQUAL_FATAL(individual.parents[1], -1);
    mat_id = individual.id;
    CU_ASSERT_EQUAL_FATAL(individual.nodes_length, 2);
    CU_ASSERT_EQUAL_FATAL(individual.nodes[0], 1);
    CU_ASSERT_EQUAL_FATAL(individual.nodes[1], 4);

    ret = tsk_treeseq_get_individual(&ts, 2, &individual);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(individual.id, 2);
    CU_ASSERT_EQUAL_FATAL(individual.flags, 3);
    CU_ASSERT_EQUAL_FATAL(individual.location_length, 1);
    CU_ASSERT_EQUAL_FATAL(individual.location[0], 0.75);
    CU_ASSERT_EQUAL_FATAL(individual.parents_length, 2);
    CU_ASSERT_EQUAL_FATAL(individual.parents[0], pat_id);
    CU_ASSERT_EQUAL_FATAL(individual.parents[1], mat_id);
    CU_ASSERT_EQUAL_FATAL(individual.nodes_length, 1);
    CU_ASSERT_EQUAL_FATAL(individual.nodes[0], 5);

    ret = tsk_treeseq_get_individual(&ts, 3, &individual);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);

    /* NaN/ifinity values are allowed in locations they do not
     * affect the integrity of the model. */
    tables.individuals.location[0] = NAN;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_treeseq_get_individual(&ts, 0, &individual);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT(!tsk_isfinite(individual.location[0]));
    tsk_treeseq_free(&ts);

    tsk_table_collection_free(&tables);
}

static void
test_simplest_bad_individuals(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "1  0   0\n"
                        "0  1   0\n";
    const char *edges = "0  1   2   0\n"
                        "0  1   2   1\n"
                        "0  1   4   3\n";
    const char *individuals = "1      0.25     -1\n"
                              "2      0.5,0.25 0\n";
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_flags_t load_flags = TSK_TS_INIT_BUILD_INDEXES;
    tsk_id_t ret_id;
    int ret;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1.0;
    parse_nodes(nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 5);
    parse_edges(edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 3);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    /* Make sure we have a good set of records */
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);

    /* Bad individual ID */
    tables.nodes.individual[0] = -2;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.nodes.individual[0] = TSK_NULL;

    /* Bad individual ID */
    tables.nodes.individual[0] = 0;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.nodes.individual[0] = TSK_NULL;

    /* Add two individuals */
    parse_individuals(individuals, &tables.individuals);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.num_rows, 2);

    /* Make sure we have a good set of records */
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);

    /* Bad individual ID */
    tables.nodes.individual[0] = 2;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.nodes.individual[0] = TSK_NULL;

    /* Bad parent ID */
    tables.individuals.parents[0] = -2;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.individuals.parents[0] = 42;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.individuals.parents[0] = TSK_NULL;

    /* Parent is self */
    tables.individuals.parents[0] = 0;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_INDIVIDUAL_SELF_PARENT);
    tsk_treeseq_free(&ts);
    tables.individuals.parents[0] = TSK_NULL;

    /* Unsorted individuals are OK*/
    tables.individuals.parents[0] = 1;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, 0);
    tsk_treeseq_free(&ts);
    tables.individuals.parents[0] = TSK_NULL;

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_simplest_bad_edges(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "1  0   0\n"
                        "0  1   0\n";
    const char *edges = "0  1   2   0\n"
                        "0  1   2   1\n"
                        "0  1   4   3\n";
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    int ret;
    tsk_id_t ret_id;
    tsk_flags_t load_flags = TSK_TS_INIT_BUILD_INDEXES;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1.0;
    parse_nodes(nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 5);
    parse_edges(edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 3);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    /* Make sure we have a good set of records */
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);

    /* Bad population ID */
    tables.nodes.population[0] = -2;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.nodes.population[0] = 0;

    /* Bad population ID */
    tables.nodes.population[0] = 1;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.nodes.population[0] = 0;

    /* Bad interval */
    tables.edges.right[0] = 0.0;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_EDGE_INTERVAL);
    tsk_treeseq_free(&ts);
    tables.edges.right[0] = 1.0;

    /* Nonfinite coords */
    tables.edges.left[0] = NAN;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_GENOME_COORDS_NONFINITE);
    tsk_treeseq_free(&ts);
    tables.edges.left[0] = 1.0;

    tables.edges.left[0] = INFINITY;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_GENOME_COORDS_NONFINITE);
    tsk_treeseq_free(&ts);
    tables.edges.left[0] = 1.0;

    tables.edges.right[0] = NAN;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_GENOME_COORDS_NONFINITE);
    tsk_treeseq_free(&ts);
    tables.edges.right[0] = 1.0;

    tables.edges.right[0] = -INFINITY;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_GENOME_COORDS_NONFINITE);
    tsk_treeseq_free(&ts);
    tables.edges.right[0] = 1.0;

    /* Left coordinate < 0. */
    tables.edges.left[0] = -1;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_LEFT_LESS_ZERO);
    tsk_treeseq_free(&ts);
    tables.edges.left[0] = 0.0;

    /* Right coordinate > sequence length. */
    tables.edges.right[0] = 2.0;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_RIGHT_GREATER_SEQ_LENGTH);
    tsk_treeseq_free(&ts);
    tables.edges.right[0] = 1.0;

    /* Duplicate records */
    tables.edges.child[0] = 1;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_DUPLICATE_EDGES);
    tsk_treeseq_free(&ts);
    tables.edges.child[0] = 0;

    /* Duplicate records */
    tables.edges.child[0] = 1;
    tables.edges.left[0] = 0.5;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_EDGES_NOT_SORTED_LEFT);
    tsk_treeseq_free(&ts);
    tables.edges.child[0] = 0;
    tables.edges.left[0] = 0.0;

    /* child node == parent */
    tables.edges.child[1] = 2;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_NODE_TIME_ORDERING);
    tsk_treeseq_free(&ts);
    tables.edges.child[1] = 1;

    /* Unsorted child nodes */
    tables.edges.child[0] = 1;
    tables.edges.child[1] = 0;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_EDGES_NOT_SORTED_CHILD);
    tsk_treeseq_free(&ts);
    tables.edges.child[0] = 0;
    tables.edges.child[1] = 1;

    /* discontinuous parent nodes */
    /* Swap rows 1 and 2 */
    tables.edges.parent[1] = 4;
    tables.edges.child[1] = 3;
    tables.edges.parent[2] = 2;
    tables.edges.child[2] = 1;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_EDGES_NONCONTIGUOUS_PARENTS);
    tsk_treeseq_free(&ts);
    tables.edges.parent[2] = 4;
    tables.edges.child[2] = 3;
    tables.edges.parent[1] = 2;
    tables.edges.child[1] = 1;

    /* Null parent */
    tables.edges.parent[0] = TSK_NULL;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NULL_PARENT);
    tsk_treeseq_free(&ts);
    tables.edges.parent[0] = 2;

    /* parent not in nodes list */
    tables.nodes.num_rows = 2;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.nodes.num_rows = 5;

    /* parent negative */
    tables.edges.parent[0] = -2;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.edges.parent[0] = 2;

    /* Null child */
    tables.edges.child[0] = TSK_NULL;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NULL_CHILD);
    tsk_treeseq_free(&ts);
    tables.edges.child[0] = 0;

    /* child node reference out of bounds */
    tables.edges.child[0] = 100;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.edges.child[0] = 0;

    /* child node reference negative */
    tables.edges.child[0] = -2;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.edges.child[0] = 0;

    /* Make sure we've preserved a good tree sequence */
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, 0);
    tsk_treeseq_free(&ts);

    tsk_table_collection_free(&tables);
}

static void
test_simplest_bad_indexes(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "1  0   0\n"
                        "0  1   0\n";
    const char *edges = "0  1   2   0\n"
                        "0  1   2   1\n"
                        "0  1   4   3\n";
    tsk_table_collection_t tables;
    tsk_id_t bad_indexes[] = { -1, 3, 4, 1000 };
    tsk_size_t j;
    tsk_id_t ret_id;
    tsk_id_t ret_num_trees;
    int ret;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1.0;
    parse_nodes(nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 5);
    parse_edges(edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 3);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    /* Make sure we have a good set of records */
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = (int) tsk_table_collection_check_integrity(&tables, TSK_CHECK_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLES_NOT_INDEXED);
    ret = tsk_table_collection_build_index(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_num_trees = tsk_table_collection_check_integrity(&tables, TSK_CHECK_TREES);
    /* TSK_CHECK_TREES returns the number of trees */
    CU_ASSERT_EQUAL_FATAL(ret_num_trees, 1);

    for (j = 0; j < sizeof(bad_indexes) / sizeof(*bad_indexes); j++) {
        tables.indexes.edge_insertion_order[0] = bad_indexes[j];
        ret_num_trees = tsk_table_collection_check_integrity(&tables, TSK_CHECK_TREES);
        CU_ASSERT_EQUAL_FATAL(ret_num_trees, TSK_ERR_EDGE_OUT_OF_BOUNDS);
        tables.indexes.edge_insertion_order[0] = 0;

        tables.indexes.edge_removal_order[0] = bad_indexes[j];
        ret_num_trees = tsk_table_collection_check_integrity(&tables, TSK_CHECK_TREES);
        CU_ASSERT_EQUAL_FATAL(ret_num_trees, TSK_ERR_EDGE_OUT_OF_BOUNDS);
        tables.indexes.edge_removal_order[0] = 0;
    }

    ret = tsk_table_collection_drop_index(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, TSK_CHECK_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLES_NOT_INDEXED);

    tsk_table_collection_free(&tables);
}

static void
test_simplest_bad_migrations(void)
{
    tsk_table_collection_t tables;
    int ret;
    tsk_id_t ret_id;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    /* insert two populations and one node to refer to. */
    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    /* One migration, node 0 goes from population 0 to 1. */
    ret_id
        = tsk_migration_table_add_row(&tables.migrations, 0, 1, 0, 0, 1, 1.0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    /* We only need basic intregity checks for migrations */
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Bad node reference */
    tables.migrations.node[0] = -1;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tables.migrations.node[0] = 0;

    /* Bad node reference */
    tables.migrations.node[0] = 1;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tables.migrations.node[0] = 0;

    /* Bad population reference */
    tables.migrations.source[0] = -1;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tables.migrations.source[0] = 0;

    /* Bad population reference */
    tables.migrations.source[0] = 2;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tables.migrations.source[0] = 0;

    /* Bad population reference */
    tables.migrations.dest[0] = -1;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tables.migrations.dest[0] = 1;

    /* Bad population reference */
    tables.migrations.dest[0] = 2;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tables.migrations.dest[0] = 1;

    /* Bad time values */
    tables.migrations.time[0] = NAN;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_TIME_NONFINITE);
    tables.migrations.time[0] = 1.0;

    tables.migrations.time[0] = INFINITY;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_TIME_NONFINITE);
    tables.migrations.time[0] = 1.0;

    /* Bad left coordinate */
    tables.migrations.left[0] = -1;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_LEFT_LESS_ZERO);
    tables.migrations.left[0] = 0;

    tables.migrations.left[0] = NAN;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_GENOME_COORDS_NONFINITE);
    tables.migrations.left[0] = 0;

    tables.migrations.left[0] = -INFINITY;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_GENOME_COORDS_NONFINITE);
    tables.migrations.left[0] = 0;

    /* Bad right coordinate */
    tables.migrations.right[0] = 2;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_RIGHT_GREATER_SEQ_LENGTH);
    tables.migrations.right[0] = 1;

    tables.migrations.right[0] = NAN;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_GENOME_COORDS_NONFINITE);
    tables.migrations.right[0] = 1;

    tables.migrations.right[0] = INFINITY;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_GENOME_COORDS_NONFINITE);
    tables.migrations.right[0] = 1;

    /* Bad interval coordinate */
    tables.migrations.right[0] = 0;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_EDGE_INTERVAL);
    tables.migrations.right[0] = 1;

    tsk_table_collection_free(&tables);
}

static void
test_simplest_migration_simplify(void)
{
    tsk_table_collection_t tables;
    int ret;
    tsk_id_t ret_id;
    tsk_id_t samples[] = { 0, 1 };

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    /* insert two populations and one node to refer to. */
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    /* One migration, node 0 goes from population 0 to 1. */
    ret_id
        = tsk_migration_table_add_row(&tables.migrations, 0, 1, 0, 0, 1, 1.0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    ret = tsk_table_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SIMPLIFY_MIGRATIONS_NOT_SUPPORTED);

    tsk_table_collection_free(&tables);
}

static void
test_simplest_overlapping_parents(void)
{
    const char *nodes = "1  0   -1\n"
                        "1  0   -1\n"
                        "0  1   -1\n";
    const char *edges = "0  1   2   0\n"
                        "0  1   2   1\n";
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_tree_t tree;
    int ret;
    tsk_flags_t load_flags = TSK_TS_INIT_BUILD_INDEXES;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    parse_nodes(nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 3);
    parse_edges(edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 2);

    tables.edges.left[0] = 0;
    tables.edges.parent[0] = 2;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tree.parent[0], 2);
    CU_ASSERT_EQUAL(tree.parent[1], 2);
    CU_ASSERT_EQUAL(tree.left_sib[2], TSK_NULL);
    CU_ASSERT_EQUAL(tree.right_sib[2], TSK_NULL);
    CU_ASSERT_EQUAL(tree.left_child[2], 0);
    CU_ASSERT_EQUAL(tree.right_child[2], 1);
    CU_ASSERT_EQUAL(tree.left_sib[0], TSK_NULL);
    CU_ASSERT_EQUAL(tree.right_sib[0], 1);
    CU_ASSERT_EQUAL(tree.left_sib[1], 0);
    CU_ASSERT_EQUAL(tree.right_sib[1], TSK_NULL);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_simplest_contradictory_children(void)
{
    const char *nodes = "1  0   -1\n"
                        "1  1   -1\n"
                        "0  1   -1\n";
    const char *edges = "0  1   1   0\n"
                        "0  1   2   0\n";
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    int ret;
    tsk_flags_t load_flags = TSK_TS_INIT_BUILD_INDEXES;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    parse_nodes(nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 3);
    parse_edges(edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 2);
    tables.sequence_length = 1.0;

    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_EDGES_CONTRADICTORY_CHILDREN);

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_simplest_overlapping_edges_simplify(void)
{
    const char *nodes = "1  0   -1\n"
                        "1  0   -1\n"
                        "1  0   -1\n"
                        "0  1   -1";
    const char *edges = "0  2   3   0\n"
                        "1  3   3   1\n"
                        "0  3   3   2\n";
    tsk_id_t samples[] = { 0, 1, 2 };
    tsk_table_collection_t tables;
    int ret;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 3;
    parse_nodes(nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 4);
    parse_edges(edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 3);

    ret = tsk_table_collection_simplify(&tables, samples, 3, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL(tables.nodes.num_rows, 4);
    CU_ASSERT_EQUAL(tables.edges.num_rows, 3);

    /* Identical to the input.
    0  2   3   0
    1  3   3   1
    0  3   3   2
    */
    CU_ASSERT_EQUAL(tables.edges.left[0], 0);
    CU_ASSERT_EQUAL(tables.edges.left[1], 1);
    CU_ASSERT_EQUAL(tables.edges.left[2], 0);
    CU_ASSERT_EQUAL(tables.edges.right[0], 2);
    CU_ASSERT_EQUAL(tables.edges.right[1], 3);
    CU_ASSERT_EQUAL(tables.edges.right[2], 3);
    CU_ASSERT_EQUAL(tables.edges.parent[0], 3);
    CU_ASSERT_EQUAL(tables.edges.parent[1], 3);
    CU_ASSERT_EQUAL(tables.edges.parent[2], 3);
    CU_ASSERT_EQUAL(tables.edges.child[0], 0);
    CU_ASSERT_EQUAL(tables.edges.child[1], 1);
    CU_ASSERT_EQUAL(tables.edges.child[2], 2);

    tsk_table_collection_free(&tables);
}

static void
test_simplest_overlapping_unary_edges_simplify(void)
{
    const char *nodes = "1  0   -1\n"
                        "1  0   -1\n"
                        "0  1   -1";
    const char *edges = "0  2   2   0\n"
                        "1  3   2   1\n";
    tsk_id_t samples[] = { 0, 1 };
    tsk_table_collection_t tables;
    int ret;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 3;
    parse_nodes(nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 3);
    parse_edges(edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 2);

    ret = tsk_table_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL(tables.nodes.num_rows, 3);
    CU_ASSERT_EQUAL(tables.edges.num_rows, 2);

    /* Because we only sample 0 and 1, the flanking unary edges are removed
     1       2       2       0
     1       2       2       1
     */
    CU_ASSERT_EQUAL(tables.edges.left[0], 1);
    CU_ASSERT_EQUAL(tables.edges.right[0], 2);
    CU_ASSERT_EQUAL(tables.edges.parent[0], 2);
    CU_ASSERT_EQUAL(tables.edges.child[0], 0);
    CU_ASSERT_EQUAL(tables.edges.left[1], 1);
    CU_ASSERT_EQUAL(tables.edges.right[1], 2);
    CU_ASSERT_EQUAL(tables.edges.parent[1], 2);
    CU_ASSERT_EQUAL(tables.edges.child[1], 1);

    tsk_table_collection_free(&tables);
}

static void
test_simplest_overlapping_unary_edges_internal_samples_simplify(void)
{
    const char *nodes = "1  0   -1\n"
                        "1  0   -1\n"
                        "1  1   -1";
    const char *edges = "0  2   2   0\n"
                        "1  3   2   1\n";
    tsk_id_t samples[] = { 0, 1, 2 };
    tsk_table_collection_t tables;
    int ret;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 3;
    parse_nodes(nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 3);
    parse_edges(edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 2);

    ret = tsk_table_collection_simplify(&tables, samples, 3, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL(tables.nodes.num_rows, 3);
    CU_ASSERT_EQUAL(tables.edges.num_rows, 2);
    /* Identical to the input.
        0  2   2   0
        1  3   2   1
     */
    CU_ASSERT_EQUAL(tables.edges.left[0], 0);
    CU_ASSERT_EQUAL(tables.edges.left[1], 1);
    CU_ASSERT_EQUAL(tables.edges.right[0], 2);
    CU_ASSERT_EQUAL(tables.edges.right[1], 3);
    CU_ASSERT_EQUAL(tables.edges.parent[0], 2);
    CU_ASSERT_EQUAL(tables.edges.parent[1], 2);
    CU_ASSERT_EQUAL(tables.edges.child[0], 0);
    CU_ASSERT_EQUAL(tables.edges.child[1], 1);

    tsk_table_collection_free(&tables);
}

static void
test_simplest_reduce_site_topology(void)
{
    /* Two trees side by side, with a site on the second one. The first
     * tree should disappear. */
    const char *nodes = "1  0   -1\n"
                        "1  0   -1\n"
                        "0  1   -1\n"
                        "0  2   -1\n";
    const char *edges = "0  1   2   0\n"
                        "0  1   2   1\n"
                        "1  2   3   0\n"
                        "1  2   3   1\n";
    const char *sites = "1.0  0\n";
    tsk_id_t samples[] = { 0, 1 };
    tsk_table_collection_t tables;
    int ret;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 2;
    parse_nodes(nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 4);
    parse_edges(edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 4);
    parse_sites(sites, &tables.sites);
    CU_ASSERT_EQUAL_FATAL(tables.sites.num_rows, 1);

    ret = tsk_table_collection_simplify(
        &tables, samples, 2, TSK_SIMPLIFY_REDUCE_TO_SITE_TOPOLOGY, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL(tables.nodes.num_rows, 3);
    CU_ASSERT_EQUAL(tables.edges.num_rows, 2);
    CU_ASSERT_EQUAL(tables.edges.left[0], 0);
    CU_ASSERT_EQUAL(tables.edges.left[1], 0);
    CU_ASSERT_EQUAL(tables.edges.right[0], 2);
    CU_ASSERT_EQUAL(tables.edges.right[1], 2);
    CU_ASSERT_EQUAL(tables.edges.parent[0], 2);
    CU_ASSERT_EQUAL(tables.edges.parent[1], 2);
    CU_ASSERT_EQUAL(tables.edges.child[0], 0);
    CU_ASSERT_EQUAL(tables.edges.child[1], 1);

    tsk_table_collection_free(&tables);
}

static void
test_simplest_simplify_defragment(void)
{
    const char *nodes = "0        2     -1\n"
                        "0        2     -1\n"
                        "0        2     -1\n"
                        "0        2     -1\n"
                        "0        2     -1\n"
                        "0        2     -1\n"
                        "0        1     -1\n"
                        "0        1     -1\n"
                        "0        1     -1\n"
                        "0        1     -1\n"
                        "0        1     -1\n"
                        "0        1     -1\n"
                        "1        0     -1\n"
                        "1        0     -1\n"
                        "1        0     -1\n"
                        "1        0     -1\n"
                        "1        0     -1\n"
                        "1        0     -1\n";
    const char *edges = "0.00000000      0.20784841      8       12\n"
                        "0.00000000      0.42202433      8       15\n"
                        "0.00000000      0.63541014      8       16\n"
                        "0.42202433      1.00000000      9       15\n"
                        "0.00000000      1.00000000      9       17\n"
                        "0.00000000      1.00000000      10      14\n"
                        "0.20784841      1.00000000      11      12\n"
                        "0.00000000      1.00000000      11      13\n"
                        "0.63541014      1.00000000      11      16\n"
                        "0.00000000      1.00000000      0       10\n"
                        "0.62102072      1.00000000      1       9\n"
                        "0.00000000      1.00000000      1       11\n"
                        "0.00000000      0.26002984      2       6\n"
                        "0.26002984      1.00000000      2       6\n"
                        "0.00000000      0.62102072      2       9\n"
                        "0.55150554      1.00000000      3       8\n"
                        "0.00000000      1.00000000      4       7\n"
                        "0.00000000      0.55150554      5       8\n";

    tsk_id_t samples[] = { 12, 13, 14, 15, 16, 17 };
    tsk_table_collection_t tables;
    int ret;

    /* This was the simplest example I could find that exercised the
     * inner loops of the simplifier_extract_ancestry function */
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;
    parse_nodes(nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 18);
    parse_edges(edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 18);

    ret = tsk_table_collection_simplify(&tables, samples, 6, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL(tables.nodes.num_rows, 10);
    CU_ASSERT_EQUAL(tables.edges.num_rows, 10);

    tsk_table_collection_free(&tables);
}

static void
test_simplest_population_filter(void)
{
    tsk_table_collection_t tables;
    tsk_id_t samples[] = { 0, 1 };
    int ret;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    tsk_population_table_add_row(&tables.populations, "0", 1);
    tsk_population_table_add_row(&tables.populations, "1", 1);
    tsk_population_table_add_row(&tables.populations, "2", 1);
    /* Two nodes referring to population 1 */
    tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 1, TSK_NULL, NULL, 0);
    tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 1, TSK_NULL, NULL, 0);

    ret = tsk_table_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.nodes.num_rows, 2);
    CU_ASSERT_EQUAL(tables.populations.num_rows, 3);
    CU_ASSERT_EQUAL(tables.populations.metadata[0], '0');
    CU_ASSERT_EQUAL(tables.populations.metadata[1], '1');
    CU_ASSERT_EQUAL(tables.populations.metadata[2], '2');

    ret = tsk_table_collection_simplify(
        &tables, samples, 2, TSK_SIMPLIFY_FILTER_POPULATIONS, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.nodes.num_rows, 2);
    CU_ASSERT_EQUAL(tables.nodes.population[0], 0);
    CU_ASSERT_EQUAL(tables.nodes.population[1], 0);
    CU_ASSERT_EQUAL(tables.populations.num_rows, 1);
    CU_ASSERT_EQUAL(tables.populations.metadata[0], '1');

    tsk_table_collection_free(&tables);
}

static void
test_simplest_individual_filter(void)
{
    tsk_table_collection_t tables;
    tsk_id_t samples[] = { 0, 1 };
    int ret;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0, NULL, 0, "0", 1);
    tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0, NULL, 0, "1", 1);
    tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0, NULL, 0, "2", 1);
    /* Two nodes referring to individual 1 */
    tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, 1, NULL, 0);
    tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, 1, NULL, 0);

    ret = tsk_table_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.nodes.num_rows, 2);
    CU_ASSERT_EQUAL(tables.individuals.num_rows, 3);
    CU_ASSERT_EQUAL(tables.individuals.metadata[0], '0');
    CU_ASSERT_EQUAL(tables.individuals.metadata[1], '1');
    CU_ASSERT_EQUAL(tables.individuals.metadata[2], '2');

    ret = tsk_table_collection_simplify(
        &tables, samples, 2, TSK_SIMPLIFY_FILTER_INDIVIDUALS, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.nodes.num_rows, 2);
    CU_ASSERT_EQUAL(tables.nodes.individual[0], 0);
    CU_ASSERT_EQUAL(tables.nodes.individual[1], 0);
    CU_ASSERT_EQUAL(tables.individuals.num_rows, 1);
    CU_ASSERT_EQUAL(tables.individuals.metadata[0], '1');

    tsk_table_collection_free(&tables);
}

static void
test_simplest_map_mutations(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0";
    const char *edges = "0  1   2   0,1\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int32_t genotypes[] = { 0, 0 };
    tsk_size_t num_transitions;
    tsk_state_transition_t *transitions;
    int32_t ancestral_state;
    int ret;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_tree_next(&t));

    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 0);
    free(transitions);

    genotypes[0] = 1;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[0].node, 0);
    CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[0].state, 1);
    free(transitions);

    genotypes[0] = -1;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 0);
    free(transitions);

    /* Check the null tree */
    genotypes[0] = 1;
    CU_ASSERT_FALSE(tsk_tree_next(&t));
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[0].node, 0);
    CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[0].state, 1);
    free(transitions);

    /* Assign the ancestral_state */
    genotypes[0] = 1;
    genotypes[1] = 1;
    ancestral_state = 0;
    ret = tsk_tree_map_mutations(&t, genotypes, NULL, TSK_MM_FIXED_ANCESTRAL_STATE,
        &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 2);
    CU_ASSERT_EQUAL_FATAL(transitions[0].node, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[0].state, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[1].node, 0);
    CU_ASSERT_EQUAL_FATAL(transitions[1].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[1].state, 1);
    free(transitions);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_nonbinary_map_mutations(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0";
    const char *edges = "0  1   4   0,1,2,3\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int32_t genotypes[] = { 0, 0, 0, 0 };
    tsk_size_t num_transitions;
    tsk_state_transition_t *transitions;
    int32_t ancestral_state;
    int ret;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_tree_next(&t));

    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 0);
    free(transitions);

    genotypes[0] = 1;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[0].node, 0);
    CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[0].state, 1);
    free(transitions);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_unary_map_mutations(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  1   0\n"
                        "0  2   0";
    const char *edges = "0  1   2   0\n"
                        "0  1   3   1\n"
                        "0  1   4   2,3\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int32_t genotypes[] = { 0, 0 };
    tsk_size_t num_transitions;
    tsk_state_transition_t *transitions;
    int32_t ancestral_state;
    int ret;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_tree_next(&t));

    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 0);
    free(transitions);

    genotypes[0] = 1;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[0].node, 2);
    CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[0].state, 1);
    free(transitions);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_non_sample_leaf_map_mutations(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  0   0\n"
                        "0  0   0";
    const char *edges = "0  1   2   0,1,3,4\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int32_t genotypes[] = { 0, 0 };
    tsk_size_t num_transitions;
    tsk_state_transition_t *transitions;
    int32_t ancestral_state;
    int ret;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_tree_next(&t));

    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 0);
    free(transitions);

    genotypes[0] = 1;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[0].node, 0);
    CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[0].state, 1);
    free(transitions);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_internal_sample_map_mutations(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  1   0";
    const char *edges = "0  1   2   0,1\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int32_t genotypes[] = { 0, 0, 0 };
    tsk_size_t num_transitions;
    tsk_state_transition_t *transitions;
    int32_t ancestral_state;
    int ret;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_tree_next(&t));

    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 0);
    free(transitions);

    genotypes[0] = 1;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[0].node, 0);
    CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[0].state, 1);
    free(transitions);

    genotypes[2] = 1;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 1);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[0].node, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[0].state, 0);
    free(transitions);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_multiple_root_map_mutations(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  1   0\n";
    const char *edges = "0  1   4   0,1\n"
                        "0  1   5   2,3\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int32_t genotypes[] = { 0, 0, 0, 0 };
    tsk_size_t num_transitions;
    tsk_state_transition_t *transitions;
    int32_t ancestral_state;
    int ret;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_tree_next(&t));

    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 0);
    free(transitions);

    genotypes[0] = 1;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[0].node, 0);
    CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[0].state, 1);
    free(transitions);

    genotypes[1] = 1;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[0].node, 4);
    CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[0].state, 1);
    free(transitions);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_chained_map_mutations(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  1   0\n"
                        "1  1   0\n"
                        "0  2   0";
    const char *edges = "0  1   2   0\n"
                        "0  1   3   1\n"
                        "0  1   4   2,3\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int32_t genotypes[] = { 0, 0, 0, 0 };
    tsk_size_t num_transitions;
    tsk_state_transition_t *transitions;
    int32_t ancestral_state;
    int ret;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_tree_next(&t));

    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 0);
    free(transitions);

    genotypes[2] = 1;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 2);
    CU_ASSERT_EQUAL_FATAL(transitions[0].node, 2);
    CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[0].state, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[1].node, 0);
    CU_ASSERT_EQUAL_FATAL(transitions[1].parent, 0);
    CU_ASSERT_EQUAL_FATAL(transitions[1].state, 0);
    free(transitions);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_mutation_edges(void)
{
    const char *nodes = "1  0   0\n"
                        "0  1   0\n"
                        "0  1   0";
    const char *edges = "0  1   1   0\n"
                        "1  2   2   0\n";
    const char *sites = "0.5  0\n"
                        "1.5  0\n";
    const char *mutations = "0    1     1\n"
                            "0    0     1\n"
                            "0    2     1\n"
                            "1    0     1\n"
                            "1    1     1\n"
                            "1    2     1\n";
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    /* We have mutations over roots, samples and just isolated nodes */
    tsk_id_t mutation_edges[] = { -1, 0, -1, 1, -1, -1 };
    tsk_size_t i, j, k, t;
    tsk_mutation_t mut;
    tsk_site_t site;
    int ret;

    tsk_treeseq_from_text(&ts, 2, nodes, edges, NULL, sites, mutations, NULL, NULL, 0);

    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);

    for (j = 0; j < tsk_treeseq_get_num_mutations(&ts); j++) {
        ret = tsk_treeseq_get_mutation(&ts, (tsk_id_t) j, &mut);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(mut.edge, mutation_edges[j]);
    }

    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    i = 0;
    for (t = 0; t < 2; t++) {
        ret = tsk_tree_next(&tree);
        CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
        for (j = 0; j < tree.sites_length; j++) {
            site = tree.sites[j];
            for (k = 0; k < site.mutations_length; k++) {
                CU_ASSERT_EQUAL(site.mutations[k].edge, mutation_edges[i]);
                i++;
            }
        }
    }
    CU_ASSERT_EQUAL(i, 6);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

/*=======================================================
 * Single tree tests.
 *======================================================*/

static void
test_single_tree_good_records(void)
{
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 7);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_treeseq_free(&ts);
}

static void
test_single_nonbinary_tree_good_records(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  2   0\n"
                        "0  3   0\n";
    const char *edges = "0 1 7 0,1,2,3\n"
                        "0 1 8 4,5\n"
                        "0 1 9 6,7,8";
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 7);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 10);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_bad_records(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_flags_t load_flags = TSK_TS_INIT_BUILD_INDEXES;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    parse_nodes(single_tree_ex_nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 7);
    parse_edges(single_tree_ex_edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 6);

    /* Not sorted in time order */
    tables.nodes.time[5] = 0.5;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_EDGES_NOT_SORTED_PARENT_TIME);
    tsk_treeseq_free(&ts);
    tables.nodes.time[5] = 2.0;

    /* Left value greater than sequence right */
    tables.edges.left[2] = 2.0;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_EDGE_INTERVAL);
    tsk_treeseq_free(&ts);
    tables.edges.left[2] = 0.0;

    /* Non finite */
    tables.nodes.time[5] = INFINITY;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_TIME_NONFINITE);
    tsk_treeseq_free(&ts);
    tables.nodes.time[5] = 2.0;

    tables.nodes.time[5] = NAN;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_TIME_NONFINITE);
    tsk_treeseq_free(&ts);
    tables.nodes.time[5] = 2.0;

    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, 0);
    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_single_tree_good_mutations(void)
{
    tsk_treeseq_t ts;
    tsk_size_t j;
    tsk_size_t num_sites = 3;
    tsk_size_t num_mutations = 7;
    tsk_site_t other_sites[num_sites];
    tsk_mutation_t other_mutations[num_mutations];
    int ret;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 7);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), num_sites);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), num_mutations);

    for (j = 0; j < num_sites; j++) {
        ret = tsk_treeseq_get_site(&ts, (tsk_id_t) j, other_sites + j);
        CU_ASSERT_EQUAL(ret, 0);
    }
    for (j = 0; j < num_mutations; j++) {
        ret = tsk_treeseq_get_mutation(&ts, (tsk_id_t) j, other_mutations + j);
        CU_ASSERT_EQUAL(ret, 0);
    }
    CU_ASSERT_EQUAL(other_sites[0].position, 0.125);
    CU_ASSERT_NSTRING_EQUAL(other_sites[0].ancestral_state, "0", 1);
    CU_ASSERT_EQUAL(other_sites[1].position, 0.25);
    CU_ASSERT_NSTRING_EQUAL(other_sites[1].ancestral_state, "0", 1);
    CU_ASSERT_EQUAL(other_sites[2].position, 0.5);
    CU_ASSERT_NSTRING_EQUAL(other_sites[2].ancestral_state, "0", 1);

    CU_ASSERT_EQUAL(other_mutations[0].id, 0);
    CU_ASSERT_EQUAL(other_mutations[0].node, 2);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[0].derived_state, "1", 1);
    CU_ASSERT_EQUAL(other_mutations[1].id, 1);
    CU_ASSERT_EQUAL(other_mutations[1].node, 4);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[1].derived_state, "1", 1);
    CU_ASSERT_EQUAL(other_mutations[2].id, 2);
    CU_ASSERT_EQUAL(other_mutations[2].node, 0);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[2].derived_state, "0", 1);
    CU_ASSERT_EQUAL(other_mutations[3].id, 3);
    CU_ASSERT_EQUAL(other_mutations[3].node, 0);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[3].derived_state, "1", 1);
    CU_ASSERT_EQUAL(other_mutations[4].id, 4);
    CU_ASSERT_EQUAL(other_mutations[4].node, 1);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[4].derived_state, "1", 1);
    CU_ASSERT_EQUAL(other_mutations[5].id, 5);
    CU_ASSERT_EQUAL(other_mutations[5].node, 2);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[5].derived_state, "1", 1);
    CU_ASSERT_EQUAL(other_mutations[6].id, 6);
    CU_ASSERT_EQUAL(other_mutations[6].node, 3);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[6].derived_state, "1", 1);

    tsk_treeseq_free(&ts);
}

static void
test_single_tree_bad_mutations(void)
{
    int ret = 0;
    const char *sites = "0       0\n"
                        "0.1     0\n"
                        "0.2     0\n";
    const char *mutations = "0   0  1  -1  0\n"
                            "1   1  1  -1  0\n"
                            "2   4  1  -1  1\n"
                            "2   1  0  2   0\n"
                            "2   1  1  3   0\n"
                            "2   2  1  -1  0\n";
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_flags_t load_flags = TSK_TS_INIT_BUILD_INDEXES;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    parse_nodes(single_tree_ex_nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 7);
    parse_edges(single_tree_ex_edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 6);
    parse_sites(sites, &tables.sites);
    parse_mutations(mutations, &tables.mutations);
    CU_ASSERT_EQUAL_FATAL(tables.sites.num_rows, 3);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.num_rows, 6);
    tables.sequence_length = 1.0;

    /* Check to make sure we have legal mutations */
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);
    tsk_treeseq_free(&ts);

    /* negative coordinate */
    tables.sites.position[0] = -1.0;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_SITE_POSITION);
    tsk_treeseq_free(&ts);
    tables.sites.position[0] = 0.0;

    /* non finite coordinates */
    tables.sites.position[0] = NAN;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_SITE_POSITION);
    tsk_treeseq_free(&ts);
    tables.sites.position[0] = 0.0;

    tables.sites.position[0] = INFINITY;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_SITE_POSITION);
    tsk_treeseq_free(&ts);
    tables.sites.position[0] = 0.0;

    /* coordinate == sequence length */
    tables.sites.position[2] = 1.0;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_SITE_POSITION);
    tsk_treeseq_free(&ts);
    tables.sites.position[2] = 0.2;

    /* coordinate > sequence length */
    tables.sites.position[2] = 1.1;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_SITE_POSITION);
    tsk_treeseq_free(&ts);
    tables.sites.position[2] = 0.2;

    /* Duplicate positions */
    tables.sites.position[0] = 0.1;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_DUPLICATE_SITE_POSITION);
    tsk_treeseq_free(&ts);
    tables.sites.position[0] = 0.0;

    /* Unsorted positions */
    tables.sites.position[0] = 0.3;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_UNSORTED_SITES);
    tsk_treeseq_free(&ts);
    tables.sites.position[0] = 0.0;

    /* site < 0 */
    tables.mutations.site[0] = -2;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.mutations.site[0] = 0;

    /* site == num_sites */
    tables.mutations.site[0] = 3;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.mutations.site[0] = 0;

    /* node = NULL */
    tables.mutations.node[0] = TSK_NULL;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.mutations.node[0] = 0;

    /* node >= num_nodes */
    tables.mutations.node[0] = 7;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.mutations.node[0] = 0;

    /* parent < -1 */
    tables.mutations.parent[0] = -2;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.mutations.parent[0] = TSK_NULL;

    /* parent >= num_mutations */
    tables.mutations.parent[0] = 7;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.mutations.parent[0] = TSK_NULL;

    /* parent on a different site */
    tables.mutations.parent[1] = 0;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_PARENT_DIFFERENT_SITE);
    tsk_treeseq_free(&ts);
    tables.mutations.parent[1] = TSK_NULL;

    /* parent is the same mutation */
    tables.mutations.parent[0] = 0;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_PARENT_EQUAL);
    tsk_treeseq_free(&ts);
    tables.mutations.parent[0] = TSK_NULL;

    /* parent_id > mutation id */
    tables.mutations.parent[3] = 4;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_PARENT_AFTER_CHILD);
    tsk_treeseq_free(&ts);
    tables.mutations.parent[3] = TSK_NULL;

    /* time < node time */
    tables.mutations.time[2] = 0;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_TIME_YOUNGER_THAN_NODE);
    tsk_treeseq_free(&ts);
    tables.mutations.time[2] = 1;

    /* time > parent mutation */
    tables.mutations.time[4] = 0.5;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_TIME_OLDER_THAN_PARENT_MUTATION);
    tsk_treeseq_free(&ts);
    tables.mutations.time[4] = 0;

    /* time > parent node */
    tables.mutations.time[0] = 1.5;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_TIME_OLDER_THAN_PARENT_NODE);
    tsk_treeseq_free(&ts);
    tables.mutations.time[0] = 0;

    /* Check to make sure we've maintained legal mutations */
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);
    tsk_treeseq_free(&ts);

    tsk_table_collection_free(&tables);
}

static void
test_single_tree_iter(void)
{
    int ret;
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  2   0\n"
                        "0  3   0\n";
    const char *edges = "0  6   4   0,1\n"
                        "0  6   5   2,3\n"
                        "0  6   6   4,5\n";
    tsk_id_t parents[] = { 4, 4, 5, 5, 6, 6, TSK_NULL };
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    tsk_id_t u, v, w;
    tsk_size_t num_samples;
    tsk_size_t num_nodes = 7;

    tsk_treeseq_from_text(&ts, 6, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_tree_print_state(&tree, _devnull);

    for (u = 0; u < (tsk_id_t) num_nodes; u++) {
        ret = tsk_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
    }
    ret = tsk_tree_get_num_samples(&tree, 0, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 1);
    ret = tsk_tree_get_num_samples(&tree, 4, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 2);
    ret = tsk_tree_get_num_samples(&tree, 6, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 4);
    ret = tsk_tree_get_mrca(&tree, 0, 1, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 4);
    ret = tsk_tree_get_mrca(&tree, 0, 2, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 6);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_single_nonbinary_tree_iter(void)
{
    int ret;
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  2   0\n"
                        "0  3   0\n";
    const char *edges = "0  1   7   0,1,2,3\n"
                        "0  1   8   4,5\n"
                        "0  1   9   6,7,8\n";
    tsk_id_t parents[] = { 7, 7, 7, 7, 8, 8, 9, 9, 9, TSK_NULL };
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    tsk_id_t u, v, w;
    tsk_size_t num_samples;
    tsk_size_t num_nodes = 10;
    tsk_size_t total_samples = 7;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_tree_print_state(&tree, _devnull);

    for (u = 0; u < (tsk_id_t) num_nodes; u++) {
        ret = tsk_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
    }
    for (u = 0; u < (tsk_id_t) total_samples; u++) {
        ret = tsk_tree_get_num_samples(&tree, u, &num_samples);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(num_samples, 1);
        CU_ASSERT_EQUAL(tree.left_child[u], TSK_NULL);
    }

    u = 7;
    ret = tsk_tree_get_num_samples(&tree, u, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 4);
    CU_ASSERT_EQUAL(tree.right_child[u], 3);
    CU_ASSERT_EQUAL(tree.left_sib[3], 2);
    CU_ASSERT_EQUAL(tree.left_sib[2], 1);
    CU_ASSERT_EQUAL(tree.left_sib[1], 0);
    CU_ASSERT_EQUAL(tree.left_sib[0], TSK_NULL);

    u = 8;
    ret = tsk_tree_get_num_samples(&tree, u, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 2);
    CU_ASSERT_EQUAL(tree.right_child[u], 5);
    CU_ASSERT_EQUAL(tree.left_sib[5], 4);
    CU_ASSERT_EQUAL(tree.left_sib[4], TSK_NULL);

    u = 9;
    ret = tsk_tree_get_num_samples(&tree, u, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 7);
    CU_ASSERT_EQUAL(tree.right_child[u], 8);
    CU_ASSERT_EQUAL(tree.left_sib[8], 7);
    CU_ASSERT_EQUAL(tree.left_sib[7], 6);
    CU_ASSERT_EQUAL(tree.left_sib[6], TSK_NULL);

    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&tree), 1);
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&tree), 9);

    ret = tsk_tree_get_mrca(&tree, 0, 1, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 7);
    ret = tsk_tree_get_mrca(&tree, 0, 4, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 9);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_general_samples_iter(void)
{
    int ret;
    const char *nodes = "0  3   0\n"
                        "0  2   0\n"
                        "0  1   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n";
    const char *edges = "0  6   2   3,4\n"
                        "0  6   1   5,6\n"
                        "0  6   0   1,2\n";
    tsk_id_t parents[] = { TSK_NULL, 0, 0, 2, 2, 1, 1 };
    const tsk_id_t *samples;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    tsk_id_t u, v, w;
    tsk_size_t num_samples;
    tsk_size_t num_nodes = 7;

    tsk_treeseq_from_text(&ts, 6, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    samples = tsk_treeseq_get_samples(&ts);
    CU_ASSERT_EQUAL(samples[0], 3);
    CU_ASSERT_EQUAL(samples[1], 4);
    CU_ASSERT_EQUAL(samples[2], 5);
    CU_ASSERT_EQUAL(samples[3], 6);

    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_tree_print_state(&tree, _devnull);

    for (u = 0; u < (tsk_id_t) num_nodes; u++) {
        ret = tsk_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
    }
    ret = tsk_tree_get_num_samples(&tree, 3, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 1);
    ret = tsk_tree_get_num_samples(&tree, 2, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 2);
    ret = tsk_tree_get_num_samples(&tree, 0, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 4);
    ret = tsk_tree_get_mrca(&tree, 3, 4, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 2);
    ret = tsk_tree_get_mrca(&tree, 3, 6, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 0);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_iter_times(void)
{
    int ret = 0;
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  2   0\n"
                        "1  3   0\n"
                        "0  1   0\n"
                        "0  4   0\n"
                        "0  5   0\n";
    const char *edges = "0  6   4   0,1\n"
                        "0  6   5   2,3\n"
                        "0  6   6   4,5\n";
    tsk_id_t parents[] = { 4, 4, 5, 5, 6, 6, TSK_NULL };
    double times[] = { 0.0, 0.0, 2.0, 3.0, 1.0, 4.0, 5.0 };
    double t;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    tsk_id_t u, v;
    uint32_t num_nodes = 7;

    tsk_treeseq_from_text(&ts, 6, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_tree_print_state(&tree, _devnull);

    for (u = 0; u < (tsk_id_t) num_nodes; u++) {
        ret = tsk_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
        ret = tsk_tree_get_time(&tree, u, &t);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(t, times[u]);
    }
    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_iter_depths(void)
{
    int ret = 0;
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  2   0\n"
                        "0  3   0\n";
    const char *edges = "0  6   4   0,1\n"
                        "0  6   5   2,3\n"
                        "0  6   6   4,5\n";
    int depths[] = { 2, 2, 2, 2, 1, 1, 0 };
    int depth;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    tsk_id_t u;
    uint32_t num_nodes = 7;

    tsk_treeseq_from_text(&ts, 6, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), num_nodes);

    for (u = 0; u < (tsk_id_t) num_nodes; u++) {
        ret = tsk_tree_get_depth(&tree, u, &depth);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(depth, depths[u]);
    }

    ret = tsk_tree_get_depth(&tree, (tsk_id_t) num_nodes + 1, &depth);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    ret = tsk_tree_get_depth(&tree, TSK_NULL, &depth);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_simplify(void)
{
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_id_t samples[] = { 0, 1 };
    int ret;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);
    verify_simplify(&ts);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.nodes.num_rows, 3);
    CU_ASSERT_EQUAL(tables.edges.num_rows, 2);

    /* Zero samples gives us the empty table collection */
    ret = tsk_table_collection_simplify(&tables, samples, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.nodes.num_rows, 0);
    CU_ASSERT_EQUAL(tables.edges.num_rows, 0);

    /* Make sure we detect unsorted edges */
    ret = tsk_treeseq_copy_tables(&ts, &tables, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    unsort_edges(&tables.edges, 0);
    ret = tsk_table_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGES_NOT_SORTED_CHILD);

    /* detect bad parents */
    ret = tsk_treeseq_copy_tables(&ts, &tables, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.edges.parent[0] = -1;
    ret = tsk_table_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NULL_PARENT);

    /* detect bad children */
    ret = tsk_treeseq_copy_tables(&ts, &tables, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.edges.child[0] = -1;
    ret = tsk_table_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NULL_CHILD);

    /* detect loops */
    ret = tsk_treeseq_copy_tables(&ts, &tables, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.edges.child[0] = tables.edges.parent[0];
    ret = tsk_table_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_NODE_TIME_ORDERING);

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_single_tree_simplify_debug(void)
{
    tsk_treeseq_t ts, simplified;
    tsk_id_t samples[] = { 0, 1 };
    int ret;
    FILE *tmp = fopen(_tmp_file_name, "w");

    CU_ASSERT_FATAL(tmp != NULL);
    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);

    tsk_set_debug_stream(tmp);
    ret = tsk_treeseq_simplify(&ts, samples, 2, TSK_DEBUG, &simplified, NULL);
    tsk_set_debug_stream(stdout);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(ftell(tmp) > 0);

    fclose(tmp);
    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&simplified);
}

static void
test_single_tree_simplify_keep_input_roots(void)
{
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_id_t samples[] = { 0, 1 };
    int ret;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);
    verify_simplify(&ts);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_simplify(
        &tables, samples, 2, TSK_SIMPLIFY_KEEP_INPUT_ROOTS, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.nodes.num_rows, 4);
    CU_ASSERT_EQUAL(tables.edges.num_rows, 3);
    CU_ASSERT_EQUAL(tables.sites.num_rows, 3);
    CU_ASSERT_EQUAL(tables.mutations.num_rows, 4);

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_single_tree_simplify_no_sample_nodes(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t t1, t2;
    tsk_id_t samples[] = { 0, 1, 2, 3 };

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_copy_tables(&ts, &t2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* We zero out the sample column in t1, and run simplify. We should
     * get back the same table */

    tsk_memset(t1.nodes.flags, 0, sizeof(*t1.nodes.flags) * t1.nodes.num_rows);

    ret = tsk_table_collection_simplify(&t1, samples, 4, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));

    tsk_table_collection_free(&t1);
    tsk_table_collection_free(&t2);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_simplify_null_samples(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t t1, t2;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_copy_tables(&ts, &t2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_simplify(&t1, NULL, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));

    tsk_table_collection_free(&t1);
    tsk_table_collection_free(&t2);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_compute_mutation_parents(void)
{
    int ret = 0;
    const char *sites = "0       0\n"
                        "0.1     0\n"
                        "0.2     0\n";
    const char *mutations = "0   0  1  -1\n"
                            "1   1  1  -1\n"
                            "2   4  1  -1\n"
                            "2   1  0  2 \n"
                            "2   1  1  3 \n"
                            "2   2  1  -1\n";
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    parse_nodes(single_tree_ex_nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 7);
    parse_edges(single_tree_ex_edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 6);
    parse_sites(sites, &tables.sites);
    parse_mutations(mutations, &tables.mutations);
    CU_ASSERT_EQUAL_FATAL(tables.sites.num_rows, 3);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.num_rows, 6);
    tables.sequence_length = 1.0;

    ret = tsk_table_collection_build_index(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Check to make sure we have legal mutations */
    ret = tsk_treeseq_init(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);

    /* Compute the mutation parents */
    verify_compute_mutation_parents(&ts);
    tsk_treeseq_free(&ts);

    /* Bad site reference */
    tables.mutations.site[0] = -1;
    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tables.mutations.site[0] = 0;

    /* Bad site reference */
    tables.mutations.site[0] = -1;
    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tables.mutations.site[0] = 0;

    /* mutation sites out of order */
    tables.mutations.site[0] = 2;
    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_UNSORTED_MUTATIONS);
    tables.mutations.site[0] = 0;

    /* sites out of order */
    tables.sites.position[0] = 0.11;
    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_UNSORTED_SITES);
    tables.sites.position[0] = 0;

    /* Bad node reference */
    tables.mutations.node[0] = -1;
    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tables.mutations.node[0] = 0;

    /* Bad node reference */
    tables.mutations.node[0] = (tsk_id_t) tables.nodes.num_rows;
    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tables.mutations.node[0] = 0;

    /* Mutations not ordered by tree */
    tables.mutations.node[2] = 1;
    tables.mutations.node[3] = 4;
    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_PARENT_AFTER_CHILD);
    tables.mutations.node[2] = 4;
    tables.mutations.node[3] = 1;

    /* Need to reset the parent field here */
    tsk_memset(
        tables.mutations.parent, 0xff, tables.mutations.num_rows * sizeof(tsk_id_t));
    /* Mutations not ordered by site */
    tables.mutations.site[3] = 1;
    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNSORTED_MUTATIONS);
    tables.mutations.site[3] = 2;

    /* Check to make sure we still have legal mutations */
    ret = tsk_table_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_init(&ts, &tables, 0);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);
    tsk_treeseq_free(&ts);

    tsk_table_collection_free(&tables);
}

static void
test_single_tree_compute_mutation_times(void)
{
    int ret = 0;
    const char *sites = "0       0\n"
                        "0.1     0\n"
                        "0.2     0\n"
                        "0.3     0\n";
    const char *mutations = "0   0  1  -1  3\n"
                            "1   1  1  -1  3\n"
                            "2   4  1  -1  8\n"
                            "2   1  0  2   4\n"
                            "2   2  1  -1  4\n"
                            "2   1  1  3   2\n"
                            "3   6  1  -1  10\n";
    /*          6          */
    /*          6          */
    /*         / \         */
    /*        /   \        */
    /*       2     \       */
    /*      /       5      */
    /*     4       / \     */
    /*    0 1,3,4 5   \    */
    /*   0   1   2     3   */

    tsk_treeseq_t ts;
    tsk_table_collection_t tables;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    parse_nodes(single_tree_ex_nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 7);
    tables.nodes.time[4] = 6;
    tables.nodes.time[5] = 8;
    tables.nodes.time[6] = 10;
    parse_edges(single_tree_ex_edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 6);
    parse_sites(sites, &tables.sites);
    parse_mutations(mutations, &tables.mutations);
    CU_ASSERT_EQUAL_FATAL(tables.sites.num_rows, 4);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.num_rows, 7);
    tables.sequence_length = 1.0;

    ret = tsk_table_collection_build_index(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Check to make sure we have legal mutations */
    ret = tsk_treeseq_init(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 7);

    /* Compute the mutation times */
    verify_compute_mutation_times(&ts);

    /* Verify consistency of individuals */
    verify_individual_nodes(&ts);
    tsk_treeseq_free(&ts);

    /* Bad random param */
    ret = tsk_table_collection_compute_mutation_times(&tables, (double *) 1, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Bad site reference */
    tables.mutations.site[0] = -1;
    ret = tsk_table_collection_compute_mutation_times(&tables, NULL, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tables.mutations.site[0] = 0;

    /* Bad site reference */
    tables.mutations.site[0] = -1;
    ret = tsk_table_collection_compute_mutation_times(&tables, NULL, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tables.mutations.site[0] = 0;

    /* mutation sites out of order */
    tables.mutations.site[0] = 2;
    ret = tsk_table_collection_compute_mutation_times(&tables, NULL, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_UNSORTED_MUTATIONS);
    tables.mutations.site[0] = 0;

    /* sites out of order */
    tables.sites.position[0] = 0.11;
    ret = tsk_table_collection_compute_mutation_times(&tables, NULL, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_UNSORTED_SITES);
    tables.sites.position[0] = 0;

    /* Bad node reference */
    tables.mutations.node[0] = -1;
    ret = tsk_table_collection_compute_mutation_times(&tables, NULL, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tables.mutations.node[0] = 0;

    /* Bad node reference */
    tables.mutations.node[0] = (tsk_id_t) tables.nodes.num_rows;
    ret = tsk_table_collection_compute_mutation_times(&tables, NULL, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tables.mutations.node[0] = 0;

    /* Mutations not ordered by site */
    tables.mutations.site[2] = 0;
    ret = tsk_table_collection_compute_mutation_times(&tables, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNSORTED_MUTATIONS);
    tables.mutations.site[2] = 2;

    ret = tsk_treeseq_init(&ts, &tables, 0);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 7);
    tsk_treeseq_free(&ts);

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_single_tree_mutation_edges(void)
{
    int ret = 0;
    tsk_size_t i, j, k;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    tsk_mutation_t mut;
    tsk_site_t site;
    tsk_id_t mutation_edges[] = { 2, 4, 0, 0, 1, 2, 3 };

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);

    for (j = 0; j < 7; j++) {
        ret = tsk_treeseq_get_mutation(&ts, (tsk_id_t) j, &mut);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(mut.edge, mutation_edges[j]);
    }

    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, TSK_TREE_OK);
    i = 0;
    for (j = 0; j < tree.sites_length; j++) {
        site = tree.sites[j];
        for (k = 0; k < site.mutations_length; k++) {
            CU_ASSERT_EQUAL(site.mutations[k].edge, mutation_edges[i]);
            i++;
        }
    }
    CU_ASSERT_EQUAL(i, 7);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_is_descendant(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t tree;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    CU_ASSERT_TRUE(tsk_tree_is_descendant(&tree, 0, 4));
    CU_ASSERT_TRUE(tsk_tree_is_descendant(&tree, 1, 4));
    CU_ASSERT_TRUE(tsk_tree_is_descendant(&tree, 0, 6));
    CU_ASSERT_TRUE(tsk_tree_is_descendant(&tree, 1, 6));
    CU_ASSERT_TRUE(tsk_tree_is_descendant(&tree, 4, 6));
    CU_ASSERT_TRUE(tsk_tree_is_descendant(&tree, 2, 5));
    CU_ASSERT_TRUE(tsk_tree_is_descendant(&tree, 3, 5));
    CU_ASSERT_TRUE(tsk_tree_is_descendant(&tree, 2, 6));
    CU_ASSERT_TRUE(tsk_tree_is_descendant(&tree, 3, 6));
    CU_ASSERT_TRUE(tsk_tree_is_descendant(&tree, 5, 6));
    /* Nodes are descendents of themselves. */
    CU_ASSERT_TRUE(tsk_tree_is_descendant(&tree, 0, 0));
    CU_ASSERT_TRUE(tsk_tree_is_descendant(&tree, 1, 1));

    CU_ASSERT_FALSE(tsk_tree_is_descendant(&tree, 0, 1));
    CU_ASSERT_FALSE(tsk_tree_is_descendant(&tree, 0, 2));
    CU_ASSERT_FALSE(tsk_tree_is_descendant(&tree, 0, 5));

    /* Out of bounds nodes always return false.*/
    CU_ASSERT_FALSE(tsk_tree_is_descendant(&tree, -1, 5));
    CU_ASSERT_FALSE(tsk_tree_is_descendant(&tree, 100, 5));
    CU_ASSERT_FALSE(tsk_tree_is_descendant(&tree, -1, -1));

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_total_branch_length(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    double length;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_total_branch_length(&tree, TSK_NULL, &length), 0);
    CU_ASSERT_EQUAL_FATAL(length, 9);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_total_branch_length(&tree, 7, &length), 0);
    CU_ASSERT_EQUAL_FATAL(length, 9);
    CU_ASSERT_EQUAL_FATAL(
        tsk_tree_get_total_branch_length(&tree, tree.virtual_root, &length), 0);
    CU_ASSERT_EQUAL_FATAL(length, 9);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_total_branch_length(&tree, 4, &length), 0);
    CU_ASSERT_EQUAL_FATAL(length, 2);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_total_branch_length(&tree, 0, &length), 0);
    CU_ASSERT_EQUAL_FATAL(length, 0);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_total_branch_length(&tree, 5, &length), 0);
    CU_ASSERT_EQUAL_FATAL(length, 4);

    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_total_branch_length(&tree, -2, &length),
        TSK_ERR_NODE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(
        tsk_tree_get_total_branch_length(&tree, 8, &length), TSK_ERR_NODE_OUT_OF_BOUNDS);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_map_mutations(void)
{
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int32_t genotypes[] = { 0, 1, 1, 1 };
    int ret = 0;
    tsk_size_t num_transitions;
    tsk_state_transition_t *transitions;
    int32_t ancestral_state, j;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 4);
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_tree_next(&t));

    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 1);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[0].node, 0);
    CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[0].state, 0);
    free(transitions);

    genotypes[0] = 1;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 1);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 0);
    free(transitions);

    genotypes[0] = 0;
    genotypes[1] = 0;
    genotypes[2] = 0;
    genotypes[3] = 0;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 0);
    free(transitions);

    for (j = 1; j < 64; j++) {
        genotypes[0] = j;
        genotypes[1] = 0;
        genotypes[2] = 0;
        genotypes[3] = 0;
        ret = tsk_tree_map_mutations(
            &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
        CU_ASSERT_EQUAL_FATAL(num_transitions, 1);
        CU_ASSERT_EQUAL_FATAL(transitions[0].node, 0);
        CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
        CU_ASSERT_EQUAL_FATAL(transitions[0].state, j);
        free(transitions);
    }

    genotypes[0] = 2;
    genotypes[1] = 1;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 2);
    CU_ASSERT_EQUAL_FATAL(transitions[0].node, 4);
    CU_ASSERT_EQUAL_FATAL(transitions[0].parent, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(transitions[0].state, 1);
    CU_ASSERT_EQUAL_FATAL(transitions[1].node, 0);
    CU_ASSERT_EQUAL_FATAL(transitions[1].parent, 0);
    CU_ASSERT_EQUAL_FATAL(transitions[1].state, 2);
    free(transitions);

    genotypes[0] = 1;
    genotypes[1] = 2;
    genotypes[2] = 3;
    genotypes[3] = 4;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 3);
    free(transitions);

    ancestral_state = 5;
    ret = tsk_tree_map_mutations(&t, genotypes, NULL, TSK_MM_FIXED_ANCESTRAL_STATE,
        &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 4);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 5);
    free(transitions);

    ancestral_state = -1;
    ret = tsk_tree_map_mutations(&t, genotypes, NULL, TSK_MM_FIXED_ANCESTRAL_STATE,
        &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_ANCESTRAL_STATE);

    ancestral_state = 64;
    ret = tsk_tree_map_mutations(&t, genotypes, NULL, TSK_MM_FIXED_ANCESTRAL_STATE,
        &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_ANCESTRAL_STATE);

    genotypes[0] = 64;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_GENOTYPE);

    genotypes[0] = -2;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_GENOTYPE);

    genotypes[0] = -1;
    genotypes[1] = -1;
    genotypes[2] = -1;
    genotypes[3] = -1;
    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_GENOTYPES_ALL_MISSING);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_map_mutations_internal_samples(void)
{
    /* Example derived from test case provoking a segfault */
    const char *nodes = "0       0.00000000000000   0\n"
                        "0       0.00000000000000   0\n"
                        "1       0.00000000000000   0\n"
                        "1       0.00000000000000   0\n"
                        "1       0.00000000000000   0\n"
                        "0       0.10792116530237   0\n"
                        "1       1.00674711128465   0\n"
                        "1       1.24675560985525   0\n"
                        "0       1.78536352520779   0\n";
    const char *edges = "0.00000000      1.00000000      5       0\n"
                        "0.00000000      1.00000000      5       2\n"
                        "0.00000000      1.00000000      6       4\n"
                        "0.00000000      1.00000000      6       5\n"
                        "0.00000000      1.00000000      7       1\n"
                        "0.00000000      1.00000000      7       3\n"
                        "0.00000000      1.00000000      8       6\n"
                        "0.00000000      1.00000000      8       7\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int32_t genotypes[] = { 0, 2, 2, 1, 0 };
    int ret = 0;
    tsk_size_t num_transitions;
    tsk_state_transition_t *transitions;
    int32_t ancestral_state;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 5);
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_tree_next(&t));

    ret = tsk_tree_map_mutations(
        &t, genotypes, NULL, 0, &ancestral_state, &num_transitions, &transitions);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ancestral_state, 0);
    CU_ASSERT_EQUAL_FATAL(num_transitions, 4);
    free(transitions);

    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_single_tree_tracked_samples(void)
{
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    tsk_id_t samples[] = { 0, 1 };
    tsk_size_t n;
    int ret;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);

    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_tree_set_tracked_samples(&tree, 2, samples);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_get_num_tracked_samples(&tree, 0, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 1);
    ret = tsk_tree_get_num_tracked_samples(&tree, 4, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 0);
    ret = tsk_tree_get_num_tracked_samples(&tree, tree.virtual_root, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 2);

    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    ret = tsk_tree_get_num_tracked_samples(&tree, 0, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 1);
    ret = tsk_tree_get_num_tracked_samples(&tree, 4, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 2);
    ret = tsk_tree_get_num_tracked_samples(&tree, 5, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 0);
    ret = tsk_tree_get_num_tracked_samples(&tree, 6, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 2);
    ret = tsk_tree_get_num_tracked_samples(&tree, tree.virtual_root, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 2);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_get_num_tracked_samples(&tree, 0, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 1);
    ret = tsk_tree_get_num_tracked_samples(&tree, 4, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 0);
    ret = tsk_tree_get_num_tracked_samples(&tree, tree.virtual_root, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 2);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    ret = tsk_tree_get_num_tracked_samples(&tree, 0, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 1);
    ret = tsk_tree_get_num_tracked_samples(&tree, 4, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 2);
    ret = tsk_tree_get_num_tracked_samples(&tree, tree.virtual_root, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 2);

    ret = tsk_tree_set_tracked_samples(&tree, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_get_num_tracked_samples(&tree, 0, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 0);
    ret = tsk_tree_get_num_tracked_samples(&tree, tree.virtual_root, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 0);

    tsk_treeseq_free(&ts);
    tsk_tree_free(&tree);
}

/*=======================================================
 * Multi tree tests.
 *======================================================*/

static void
test_simple_multi_tree(void)
{
    // clang-format off
    tsk_id_t parents[] = {
        6, 5, 8, 5, TSK_NULL, 6, 8, TSK_NULL, TSK_NULL,
        6, 5, 4, 4, 5, 6, TSK_NULL, TSK_NULL, TSK_NULL,
        7, 5, 4, 4, 5, 7, TSK_NULL, TSK_NULL, TSK_NULL,
    };
    // clang-format on
    uint32_t num_trees = 3;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, paper_ex_sites,
        paper_ex_mutations, paper_ex_individuals, NULL, 0);
    verify_trees(&ts, num_trees, parents);
    tsk_treeseq_free(&ts);
}

static void
test_unary_multi_tree(void)
{
    // clang-format off
    tsk_id_t parents[] = {
        6, 5, 7, 5, TSK_NULL, 6, 8, 8, TSK_NULL,
        6, 5, 4, 4, 5, 6, 8, TSK_NULL, TSK_NULL,
        7, 5, 4, 4, 5, 7, TSK_NULL, TSK_NULL, TSK_NULL,
    };
    // clang-format on
    tsk_treeseq_t ts;
    uint32_t num_trees = 3;

    tsk_treeseq_from_text(&ts, 10, unary_ex_nodes, unary_ex_edges, NULL, unary_ex_sites,
        unary_ex_mutations, NULL, NULL, 0);
    verify_trees(&ts, num_trees, parents);
    tsk_treeseq_free(&ts);
}

static void
test_internal_sample_multi_tree(void)
{
    // clang-format off
    tsk_id_t parents[] = {
        7, 5, 4, 4, 5, 7, TSK_NULL, TSK_NULL, TSK_NULL,
        4, 5, 4, 8, 5, 8, TSK_NULL, TSK_NULL, TSK_NULL,
        6, 5, 4, 4, 5, 6, TSK_NULL, TSK_NULL, TSK_NULL,
    };
    // clang-format on
    tsk_treeseq_t ts;
    uint32_t num_trees = 3;

    tsk_treeseq_from_text(&ts, 10, internal_sample_ex_nodes, internal_sample_ex_edges,
        NULL, internal_sample_ex_sites, internal_sample_ex_mutations, NULL, NULL, 0);
    verify_trees(&ts, num_trees, parents);
    tsk_treeseq_free(&ts);
}

static void
test_internal_sample_simplified_multi_tree(void)
{
    int ret;
    tsk_treeseq_t ts, simplified;
    tsk_id_t samples[] = { 2, 3, 5 };
    tsk_id_t node_map[9];
    tsk_id_t z = TSK_NULL;
    // clang-format off
    tsk_id_t parents[] = {
    /*  0  1  2  3  4 */
        3, 3, z, 2, z,
        2, 4, 4, z, z,
        3, 3, z, 2, z,
    };
    // clang-format on
    uint32_t num_trees = 3;

    tsk_treeseq_from_text(&ts, 10, internal_sample_ex_nodes, internal_sample_ex_edges,
        NULL, internal_sample_ex_sites, internal_sample_ex_mutations, NULL, NULL, 0);
    ret = tsk_treeseq_simplify(&ts, samples, 3, 0, &simplified, node_map);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(node_map[2], 0);
    CU_ASSERT_EQUAL(node_map[3], 1);
    CU_ASSERT_EQUAL(node_map[5], 2);

    verify_trees(&simplified, num_trees, parents);
    tsk_treeseq_free(&simplified);
    tsk_treeseq_free(&ts);
}

static void
test_nonbinary_multi_tree(void)
{
    /* We make one mutation for each tree */
    // clang-format off
    tsk_id_t parents[] = {
        8, 8, 8, 8, 10, 10, 9, 10, 9, 12, 12, TSK_NULL, TSK_NULL,
        8, 8, 8, 8, 10, 11, 9, 10, 9, 11, 12, 12, TSK_NULL,
    };
    // clang-format on

    tsk_treeseq_t ts;
    uint32_t num_trees = 2;

    tsk_treeseq_from_text(&ts, 100, nonbinary_ex_nodes, nonbinary_ex_edges, NULL,
        nonbinary_ex_sites, nonbinary_ex_mutations, NULL, NULL, 0);
    verify_trees(&ts, num_trees, parents);
    tsk_treeseq_free(&ts);
}

static void
test_simplify_keep_input_roots_multi_tree(void)
{

    /*
    0.25┊     8   ┊         ┊         ┊
        ┊   ┏━┻━┓ ┊         ┊         ┊
    0.20┊   ┃   ┃ ┊         ┊   7     ┊
        ┊   ┃   ┃ ┊         ┊ ┏━┻━┓   ┊
    0.17┊   6   ┃ ┊   6     ┊ ┃   ┃   ┊
        ┊ ┏━┻┓  ┃ ┊ ┏━┻━┓   ┊ ┃   ┃   ┊
    0.09┊ ┃  5  ┃ ┊ ┃   5   ┊ ┃   5   ┊
        ┊ ┃ ┏┻┓ ┃ ┊ ┃ ┏━┻┓  ┊ ┃ ┏━┻┓  ┊
    0.07┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃  4  ┊ ┃ ┃  4  ┊
        ┊ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┏┻┓ ┊ ┃ ┃ ┏┻┓ ┊
    0.00┊ 0 1 3 2 ┊ 0 1 2 3 ┊ 0 1 2 3 ┊
      0.00      2.00      7.00      10.00

    Simplifies to

    0.25┊  4  ┊     ┊     ┊
        ┊  ┃  ┊     ┊     ┊
    0.20┊  ┃  ┊     ┊  3  ┊
        ┊  ┃  ┊     ┊ ┏┻┓ ┊
    0.17┊  2  ┊  2  ┊ ┃ ┃ ┊
        ┊ ┏┻┓ ┊ ┏┻┓ ┊ ┃ ┃ ┊
    0.00┊ 0 1 ┊ 0 1 ┊ 0 1 ┊
      0.00  2.00  7.00  10.00

    */
    int ret = 0;
    // clang-format off
    tsk_id_t parents[] = {
        2, 2, 4, -1, -1,
        2, 2, -1, -1, -1,
        3, 3, -1, -1, -1,
    };
    // clang-format on
    uint32_t num_trees = 3;

    tsk_id_t samples[] = { 0, 3 };
    tsk_treeseq_t ts, simplified;

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, paper_ex_sites,
        paper_ex_mutations, paper_ex_individuals, NULL, 0);
    tsk_treeseq_dump(&ts, "tmp.trees", 0);
    ret = tsk_treeseq_simplify(
        &ts, samples, 2, TSK_SIMPLIFY_KEEP_INPUT_ROOTS, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_trees(&simplified, num_trees, parents);

    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&simplified);
}

static void
test_left_to_right_multi_tree(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  0.090   0\n"
                        "0  0.170   0\n"
                        "0  0.253   0\n"
                        "0  0.071   0\n"
                        "0  0.202   0\n";
    const char *edges = "2 10 7 2,3\n"
                        "0 2  4 1\n"
                        "2 10 4 1\n"
                        "0 2  4 3\n"
                        "2 10 4 7\n"
                        "0 7  5 0,4\n"
                        "7 10 8 0,4\n"
                        "0 2  6 2,5\n";
    const char *sites = "1      0\n"
                        "4.5    0\n"
                        "8.5    0\n";
    const char *mutations = "0    2    1\n"
                            "1    0    1\n"
                            "2    4    1\n";

    // clang-format off
    tsk_id_t parents[] = {
        5, 4, 6, 4, 5, 6, TSK_NULL, TSK_NULL, TSK_NULL,
        5, 4, 7, 7, 5, TSK_NULL, TSK_NULL, 4, TSK_NULL,
        8, 4, 7, 7, 8, TSK_NULL, TSK_NULL, 4, TSK_NULL,
    };
    // clang-format on
    tsk_treeseq_t ts;
    uint32_t num_trees = 3;

    tsk_treeseq_from_text(&ts, 10, nodes, edges, NULL, sites, mutations, NULL, NULL, 0);
    verify_trees(&ts, num_trees, parents);
    verify_tree_next_prev(&ts);
    tsk_treeseq_free(&ts);
}

static void
test_gappy_multi_tree(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  0.090   0\n"
                        "0  0.170   0\n"
                        "0  0.253   0\n"
                        "0  0.071   0\n"
                        "0  0.202   0\n";
    const char *edges = "2 7  7 2\n"
                        "8 10 7 2\n"
                        "2 7  7 3\n"
                        "8 10 7 3\n"
                        "1 2  4 1\n"
                        "2 7  4 1\n"
                        "8 10 4 1\n"
                        "1 2  4 3\n"
                        "2 7  4 7\n"
                        "8 10 4 7\n"
                        "1 7  5 0,4\n"
                        "8 10 8 0,4\n"
                        "1 2  6 2,5\n";
    tsk_id_t z = TSK_NULL;
    // clang-format off
    tsk_id_t parents[] = {
        z, z, z, z, z, z, z, z, z,
        5, 4, 6, 4, 5, 6, z, z, z,
        5, 4, 7, 7, 5, z, z, 4, z,
        z, z, z, z, z, z, z, z, z,
        8, 4, 7, 7, 8, z, z, 4, z,
        z, z, z, z, z, z, z, z, z,
    };
    // clang-format on
    tsk_treeseq_t ts;
    uint32_t num_trees = 6;

    tsk_treeseq_from_text(&ts, 12, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    verify_trees(&ts, num_trees, parents);
    verify_tree_next_prev(&ts);
    tsk_treeseq_free(&ts);
}

static void
test_tsk_treeseq_bad_records(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    uint32_t num_trees = 3;
    // clang-format off
    tsk_id_t parents[] = {
        6, 5, 8, 5, TSK_NULL, 6, 8, TSK_NULL, TSK_NULL,
        6, 5, 4, 4, 5, 6, TSK_NULL, TSK_NULL, TSK_NULL,
        7, 5, 4, 4, 5, 7, TSK_NULL, TSK_NULL, TSK_NULL,
    };
    // clang-format on
    tsk_flags_t load_flags = TSK_TS_INIT_BUILD_INDEXES;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 10;
    parse_nodes(paper_ex_nodes, &tables.nodes);
    parse_edges(paper_ex_edges, &tables.edges);
    parse_individuals(paper_ex_individuals, &tables.individuals);

    /* Make sure we have a good set of records */
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ts.num_trees, 3);
    verify_trees(&ts, num_trees, parents);
    tsk_treeseq_free(&ts);

    /* Left value greater than right */
    tables.edges.left[0] = 10.0;
    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_EDGE_INTERVAL);
    tsk_treeseq_free(&ts);
    tables.edges.left[0] = 2.0;

    ret = tsk_treeseq_init(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, 0);
    verify_trees(&ts, num_trees, parents);
    tsk_treeseq_free(&ts);

    tsk_table_collection_free(&tables);
}

/*=======================================================
 * Diff iter tests.
 *======================================================*/

static void
test_simple_diff_iter(void)
{
    int ret;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, NULL, NULL,
        paper_ex_individuals, NULL, 0);

    verify_tree_diffs(&ts, 0);
    verify_tree_diffs(&ts, TSK_INCLUDE_TERMINAL);

    ret = tsk_treeseq_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_nonbinary_diff_iter(void)
{
    int ret;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 100, nonbinary_ex_nodes, nonbinary_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    verify_tree_diffs(&ts, 0);
    verify_tree_diffs(&ts, TSK_INCLUDE_TERMINAL);

    ret = tsk_treeseq_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_unary_diff_iter(void)
{
    int ret;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(
        &ts, 10, unary_ex_nodes, unary_ex_edges, NULL, NULL, NULL, NULL, NULL, 0);
    verify_tree_diffs(&ts, 0);
    verify_tree_diffs(&ts, TSK_INCLUDE_TERMINAL);

    ret = tsk_treeseq_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_internal_sample_diff_iter(void)
{
    int ret;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 10, internal_sample_ex_nodes, internal_sample_ex_edges,
        NULL, NULL, NULL, NULL, NULL, 0);
    verify_tree_diffs(&ts, 0);
    verify_tree_diffs(&ts, TSK_INCLUDE_TERMINAL);

    ret = tsk_treeseq_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_multiroot_mrca(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    tsk_id_t mrca;

    tsk_treeseq_from_text(&ts, 10, multiroot_ex_nodes, multiroot_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&tree, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_mrca(&tree, 0, 0, &mrca), 0);
    CU_ASSERT_EQUAL(mrca, 0);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_mrca(&tree, 0, 1, &mrca), 0);
    CU_ASSERT_EQUAL(mrca, 10);
    /* MRCA of two nodes in different subtrees is TSK_NULL */
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_mrca(&tree, 0, 2, &mrca), 0);
    CU_ASSERT_EQUAL(mrca, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_mrca(&tree, 2, 0, &mrca), 0);
    CU_ASSERT_EQUAL(mrca, TSK_NULL);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_multiroot_diff_iter(void)
{
    int ret;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 10, multiroot_ex_nodes, multiroot_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    verify_tree_diffs(&ts, 0);
    verify_tree_diffs(&ts, TSK_INCLUDE_TERMINAL);

    ret = tsk_treeseq_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_empty_diff_iter(void)
{
    int ret;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(
        &ts, 10, empty_ex_nodes, empty_ex_edges, NULL, NULL, NULL, NULL, NULL, 0);
    verify_tree_diffs(&ts, 0);
    verify_tree_diffs(&ts, TSK_INCLUDE_TERMINAL);

    ret = tsk_treeseq_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
}

/*=======================================================
 * Sample sets
 *======================================================*/

static void
test_simple_sample_sets(void)
{
    // clang-format off
    sample_count_test_t tests[] = {
        {0, 0, 1}, {0, 5, 2}, {0, 6, 3},
        {1, 4, 2}, {1, 5, 3}, {1, 6, 4}};
    // clang-format on
    uint32_t num_tests = 6;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, NULL, NULL,
        paper_ex_individuals, NULL, 0);
    verify_sample_counts(&ts, num_tests, tests);
    verify_sample_sets(&ts);

    tsk_treeseq_free(&ts);
}

static void
test_nonbinary_sample_sets(void)
{
    // clang-format off
    sample_count_test_t tests[] = {
        {0, 0, 1}, {0, 8, 4}, {0, 9, 5}, {0, 10, 3}, {0, 12, 8},
        {1, 5, 1}, {1, 8, 4}, {1, 9, 5}, {0, 10, 2}, {0, 11, 1}};
    // clang-format on
    uint32_t num_tests = 8;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 100, nonbinary_ex_nodes, nonbinary_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    verify_sample_counts(&ts, num_tests, tests);
    verify_sample_sets(&ts);

    tsk_treeseq_free(&ts);
}

static void
test_internal_sample_sample_sets(void)
{
    // clang-format off
    sample_count_test_t tests[] = {
        {0, 0, 1}, {0, 5, 4}, {0, 4, 2}, {0, 7, 5},
        {1, 4, 2}, {1, 5, 4}, {1, 8, 5},
        {2, 5, 4}, {2, 6, 5}};
    // clang-format on
    uint32_t num_tests = 9;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 10, internal_sample_ex_nodes, internal_sample_ex_edges,
        NULL, NULL, NULL, NULL, NULL, 0);
    verify_sample_counts(&ts, num_tests, tests);
    verify_sample_sets(&ts);

    tsk_treeseq_free(&ts);
}

static void
test_non_sample_leaf_sample_lists(void)
{
    const char *nodes = "1  0   0\n"
                        "0  0   0\n"
                        "1  2   0\n";
    const char *edges = "0 1  2 0,1\n";
    const tsk_id_t left_sample[3] = { 0, -1, 1 };
    const tsk_id_t right_sample[3] = { 0, -1, 0 };
    const tsk_id_t next_sample[2] = { -1, 0 };
    const tsk_id_t samples[2] = { 0, 2 };
    const tsk_id_t sample_index_map[3] = { 0, -1, 1 };
    tsk_treeseq_t ts;
    tsk_tree_t t;
    tsk_id_t i;
    int ret;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);

    ret = tsk_tree_init(&t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    for (i = 0; i < 3; i++) {
        CU_ASSERT_EQUAL_FATAL(left_sample[i], t.left_sample[i]);
        CU_ASSERT_EQUAL_FATAL(right_sample[i], t.right_sample[i]);
        CU_ASSERT_EQUAL_FATAL(sample_index_map[i], ts.sample_index_map[i]);
    }
    for (i = 0; i < 2; i++) {
        CU_ASSERT_EQUAL_FATAL(next_sample[i], t.next_sample[i]);
        CU_ASSERT_EQUAL_FATAL(samples[i], t.samples[i]);
    }

    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_virtual_root_properties(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int depth;
    double time, length;
    tsk_id_t node;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_depth(&t, t.virtual_root, &depth), 0)
    CU_ASSERT_EQUAL_FATAL(depth, -1);

    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_time(&t, t.virtual_root, &time), 0)
    /* Workaround problems in IEEE floating point macros. We may want to
     * add tsk_isinf (like tsk_isnan) at some point, but not worth it just
     * for this test case */
    CU_ASSERT_TRUE(isinf((float) time));

    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_mrca(&t, t.virtual_root, 0, &node), 0)
    CU_ASSERT_EQUAL(node, t.virtual_root);

    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_mrca(&t, 0, t.virtual_root, &node), 0)
    CU_ASSERT_EQUAL(node, t.virtual_root);

    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_parent(&t, t.virtual_root, &node), 0)
    CU_ASSERT_EQUAL(node, TSK_NULL);

    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_branch_length(&t, t.virtual_root, &length), 0)
    CU_ASSERT_EQUAL(length, 0);

    /* The definition of "descendant" is that node v is on the path from
     * u to a root. Since there is no parent link from roots to the
     * virtual_root, it's consistent with this definition to return false
     * for every node. */
    CU_ASSERT_FALSE(tsk_tree_is_descendant(&t, 0, t.virtual_root));
    CU_ASSERT_FALSE(
        tsk_tree_is_descendant(&t, t.left_child[t.virtual_root], t.virtual_root));
    CU_ASSERT_FALSE(tsk_tree_is_descendant(&t, t.virtual_root, 0));
    /* The virtual_root *is* a descendent of itself, though. This is
     * consistent with other nodes that are not "in" the tree being
     * descendents of themselves, despite not being roots in the tree. */
    CU_ASSERT_TRUE(tsk_tree_is_descendant(&t, t.virtual_root, t.virtual_root));

    CU_ASSERT_FALSE(tsk_tree_is_sample(&t, t.virtual_root));

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_no_sample_count_semantics(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t t;
    tsk_id_t nodes;
    tsk_size_t n;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);

    ret = tsk_tree_init(&t, &ts, TSK_NO_SAMPLE_COUNTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&t), 0);
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&t), TSK_NULL);
    CU_ASSERT_EQUAL(tsk_tree_preorder(&t, &nodes, &n), TSK_ERR_UNSUPPORTED_OPERATION);
    CU_ASSERT_EQUAL(tsk_tree_postorder(&t, &nodes, &n), TSK_ERR_UNSUPPORTED_OPERATION);
    CU_ASSERT_EQUAL(tsk_tree_preorder_samples_from(&t, -1, &nodes, &n),
        TSK_ERR_UNSUPPORTED_OPERATION);

    CU_ASSERT_EQUAL(tsk_tree_preorder_from(&t, t.virtual_root, &nodes, &n),
        TSK_ERR_UNSUPPORTED_OPERATION);
    CU_ASSERT_EQUAL(tsk_tree_postorder_from(&t, t.virtual_root, &nodes, &n),
        TSK_ERR_UNSUPPORTED_OPERATION);
    CU_ASSERT_EQUAL(tsk_tree_preorder_samples_from(&t, t.virtual_root, &nodes, &n),
        TSK_ERR_UNSUPPORTED_OPERATION);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

/*=======================================================
 * Tree traversals
 *=======================================================*/

static void
verify_node_lists(tsk_size_t n, tsk_id_t *l1, tsk_id_t *l2)
{
    tsk_size_t j;

    for (j = 0; j < n; j++) {
        /* printf("%d %d\n", l1[j], l2[j]); */
        CU_ASSERT_EQUAL(l1[j], l2[j]);
    }
}

static void
test_single_tree_traversal(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t t;
    tsk_size_t num_nodes = 7;
    tsk_id_t preorder[] = { 6, 4, 0, 1, 5, 2, 3 };
    tsk_id_t preorder_vr[] = { 7, 6, 4, 0, 1, 5, 2, 3 };
    tsk_id_t preorder_samples[] = { 0, 1, 2, 3 };
    tsk_id_t postorder[] = { 0, 1, 4, 2, 3, 5, 6 };
    tsk_id_t postorder_vr[] = { 0, 1, 4, 2, 3, 5, 6, 7 };
    tsk_id_t nodes[num_nodes + 1];
    tsk_size_t n;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    ret = tsk_tree_preorder(&t, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, num_nodes);
    verify_node_lists(n, nodes, preorder);

    ret = tsk_tree_preorder_from(&t, -1, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, num_nodes);
    verify_node_lists(n, nodes, preorder);

    ret = tsk_tree_preorder_from(&t, t.virtual_root, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, num_nodes + 1);
    verify_node_lists(n, nodes, preorder_vr);

    ret = tsk_tree_preorder_samples_from(&t, -1, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 4);
    verify_node_lists(n, nodes, preorder_samples);

    ret = tsk_tree_preorder_samples_from(&t, t.virtual_root, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 4);
    verify_node_lists(n, nodes, preorder_samples);

    ret = tsk_tree_preorder_from(&t, 5, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 3);
    verify_node_lists(n, nodes, preorder + 4);

    ret = tsk_tree_preorder_samples_from(&t, 5, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 2);
    verify_node_lists(n, nodes, preorder_samples + 2);

    ret = tsk_tree_postorder(&t, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, num_nodes);
    verify_node_lists(n, nodes, postorder);

    ret = tsk_tree_postorder_from(&t, -1, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, num_nodes);
    verify_node_lists(n, nodes, postorder);

    ret = tsk_tree_postorder_from(&t, t.virtual_root, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, num_nodes + 1);
    verify_node_lists(n, nodes, postorder_vr);

    ret = tsk_tree_postorder_from(&t, 4, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 3);
    verify_node_lists(n, nodes, postorder);

    /* Check errors */
    ret = tsk_tree_preorder_from(&t, -2, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    ret = tsk_tree_preorder_from(&t, 8, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    ret = tsk_tree_preorder_samples_from(&t, -2, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    ret = tsk_tree_preorder_samples_from(&t, 8, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    ret = tsk_tree_postorder_from(&t, -2, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    ret = tsk_tree_postorder_from(&t, 8, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

/* printed out in tree order.
0.90┊             ┊         11  ┊             ┊
    ┊             ┊         ┏┻┓ ┊             ┊
0.80┊         10  ┊         ┃ ┃ ┊             ┊
    ┊         ┏┻┓ ┊         ┃ ┃ ┊             ┊
0.40┊     9   ┃ ┃ ┊    9    ┃ ┃ ┊     9       ┊
    ┊   ┏━┻┓  ┃ ┃ ┊  ┏━┻━┓  ┃ ┃ ┊   ┏━┻━━┓    ┊
0.30┊   ┃  ┃  ┃ ┃ ┊  ┃   8  ┃ ┃ ┊   ┃    8    ┊
    ┊   ┃  ┃  ┃ ┃ ┊  ┃  ┏┻┓ ┃ ┃ ┊   ┃   ┏┻┓   ┊
0.20┊   ┃  7  ┃ ┃ ┊  7  ┃ ┃ ┃ ┃ ┊   7   ┃ ┃   ┊
    ┊   ┃ ┏┻┓ ┃ ┃ ┊ ┏┻┓ ┃ ┃ ┃ ┃ ┊ ┏━┻┓  ┃ ┃   ┊
0.10┊   ┃ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┃ ┃ ┊ ┃  6  ┃ ┃   ┊
    ┊   ┃ ┃ ┃ ┃ ┃ ┊ ┃ ┃ ┃ ┃ ┃ ┃ ┊ ┃ ┏┻┓ ┃ ┃   ┊
0.00┊ 5 2 3 4 0 1 ┊ 3 4 1 2 0 5 ┊ 4 0 3 1 2 5 ┊
    0             4             8            10
*/

static void
test_multiroot_tree_traversal(void)
{
    int ret;
    tsk_treeseq_t ts;

    tsk_tree_t t;
    tsk_id_t preorder[] = { 5, 9, 2, 7, 3, 4, 10, 0, 1 };
    tsk_id_t preorder_vr[] = { 12, 5, 9, 2, 7, 3, 4, 10, 0, 1 };
    tsk_id_t preorder_samples[] = { 5, 2, 3, 4, 0, 1 };
    tsk_id_t postorder[] = { 5, 2, 3, 4, 7, 9, 0, 1, 10 };
    tsk_id_t postorder_vr[] = { 5, 2, 3, 4, 7, 9, 0, 1, 10, 12 };
    tsk_id_t nodes[13];
    tsk_size_t n;

    tsk_treeseq_from_text(&ts, 10, multiroot_ex_nodes, multiroot_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    ret = tsk_tree_preorder(&t, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 9);
    verify_node_lists(n, nodes, preorder);

    ret = tsk_tree_preorder_from(&t, -1, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 9);
    verify_node_lists(n, nodes, preorder);

    ret = tsk_tree_preorder_from(&t, t.virtual_root, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 10);
    verify_node_lists(n, nodes, preorder_vr);

    ret = tsk_tree_preorder_from(&t, 10, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 3);
    verify_node_lists(n, nodes, preorder + 6);

    ret = tsk_tree_preorder_samples_from(&t, -1, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 6);
    verify_node_lists(n, nodes, preorder_samples);

    ret = tsk_tree_preorder_samples_from(&t, t.virtual_root, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 6);
    verify_node_lists(n, nodes, preorder_samples);

    ret = tsk_tree_preorder_samples_from(&t, 5, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 1);
    verify_node_lists(n, nodes, preorder_samples);

    ret = tsk_tree_preorder_samples_from(&t, 10, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 2);
    verify_node_lists(n, nodes, preorder_samples + 4);

    ret = tsk_tree_postorder(&t, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 9);
    verify_node_lists(n, nodes, postorder);

    ret = tsk_tree_postorder_from(&t, -1, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 9);
    verify_node_lists(n, nodes, postorder);

    ret = tsk_tree_postorder_from(&t, t.virtual_root, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 10);
    verify_node_lists(n, nodes, postorder_vr);

    ret = tsk_tree_postorder_from(&t, 10, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 3);
    verify_node_lists(n, nodes, postorder + 6);

    /* Nodes that aren't "in" the tree have singleton traversal lists and
     * connect to no samples */

    ret = tsk_tree_preorder_from(&t, 11, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 1);
    CU_ASSERT_EQUAL_FATAL(nodes[0], 11);

    ret = tsk_tree_postorder_from(&t, 11, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 1);
    CU_ASSERT_EQUAL_FATAL(nodes[0], 11);

    ret = tsk_tree_preorder_samples_from(&t, 11, nodes, &n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(n, 0);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_seek_multi_tree(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t t;
    double breakpoints[] = { 0, 2, 7, 10 };
    tsk_id_t num_trees = 3;
    tsk_id_t j, k;

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, NULL, NULL,
        paper_ex_individuals, NULL, 0);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < num_trees; j++) {
        ret = tsk_tree_seek(&t, breakpoints[j], 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(t.index, j);
        for (k = 0; k < num_trees; k++) {
            ret = tsk_tree_seek(&t, breakpoints[k], 0);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(t.index, k);
        }
    }

    ret = tsk_tree_seek(&t, 1.99999, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(t.index, 0);
    ret = tsk_tree_seek(&t, 6.99999, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(t.index, 1);
    ret = tsk_tree_seek(&t, 9.99999, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(t.index, 2);

    tsk_tree_free(&t);

    /* Seek to all positions from a new tree. */
    for (j = 0; j < num_trees; j++) {
        ret = tsk_tree_init(&t, &ts, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tree_seek(&t, breakpoints[j], 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(t.index, j);
        tsk_tree_free(&t);
    }

    /* Seek to all positions from a non-new tree in the null state*/
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < num_trees; j++) {
        ret = tsk_tree_seek(&t, 0, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tree_prev(&t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(t.index, -1);
        ret = tsk_tree_seek(&t, breakpoints[j], 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(t.index, j);
    }
    tsk_tree_free(&t);

    tsk_treeseq_free(&ts);
}

static void
test_seek_errors(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t t;

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, NULL, NULL,
        paper_ex_individuals, NULL, 0);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_seek(&t, -1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SEEK_OUT_OF_BOUNDS);
    ret = tsk_tree_seek(&t, 10, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SEEK_OUT_OF_BOUNDS);
    ret = tsk_tree_seek(&t, 11, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SEEK_OUT_OF_BOUNDS);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

/*=======================================================
 * KC Distance tests.
 *=======================================================*/

static void
test_isolated_node_kc(void)
{
    const char *single_leaf = "1 0 0";
    const char *single_internal = "0 0 0";
    const char *edges = "";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int ret;
    double result = 0;

    tsk_treeseq_from_text(&ts, 1, single_leaf, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_treeseq_kc_distance(&ts, &ts, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(result, 0);
    ret = tsk_tree_init(&t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    ret = tsk_tree_kc_distance(&t, &t, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(result, 0);
    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);

    tsk_treeseq_from_text(
        &ts, 1, single_internal, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_treeseq_kc_distance(&ts, &ts, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MULTIPLE_ROOTS);
    ret = tsk_tree_init(&t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_left_root(&t), TSK_NULL);
    ret = tsk_tree_kc_distance(&t, &t, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MULTIPLE_ROOTS);
    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_single_tree_kc(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t t, other_t;
    double result = 0;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_kc_distance(&ts, &ts, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(result, 0);
    ret = tsk_treeseq_kc_distance(&ts, &ts, 1, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(result, 0);

    ret = tsk_tree_init(&t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    ret = tsk_tree_init(&other_t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&other_t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    ret = tsk_tree_copy(&t, &other_t, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    check_trees_identical(&t, &other_t);

    ret = tsk_tree_kc_distance(&t, &other_t, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(result, 0);
    ret = tsk_tree_kc_distance(&t, &other_t, 1, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(result, 0);

    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
    tsk_tree_free(&other_t);
}

static void
test_two_trees_kc(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  2   0\n"
                        "0  3   0\n";
    const char *nodes_other = "1  0   0\n"
                              "1  0   0\n"
                              "1  0   0\n"
                              "0  4   0\n"
                              "0  6   0\n";
    const char *edges = "0 1  3 0,1\n"
                        "0 1  4 2,3\n";
    int ret;
    tsk_treeseq_t ts, other_ts;
    tsk_tree_t t, other_t;
    double result = 0;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    tsk_treeseq_from_text(
        &other_ts, 1, nodes_other, edges, NULL, NULL, NULL, NULL, NULL, 0);

    ret = tsk_treeseq_kc_distance(&ts, &other_ts, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(result, 0);
    ret = tsk_treeseq_kc_distance(&ts, &other_ts, 1, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_DOUBLE_EQUAL_FATAL(result, 4.243, 1e-2);

    ret = tsk_tree_init(&other_t, &other_ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&other_t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    ret = tsk_tree_kc_distance(&t, &other_t, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(result, 0);
    ret = tsk_tree_kc_distance(&t, &other_t, 1, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_DOUBLE_EQUAL_FATAL(result, 4.243, 1e-2);

    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&other_ts);
    tsk_tree_free(&t);
    tsk_tree_free(&other_t);
}

static void
test_empty_tree_kc(void)
{
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_tree_t t;
    tsk_id_t v;
    int ret;
    double result = 0;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_SEQUENCE_LENGTH);
    tsk_treeseq_free(&ts);
    tables.sequence_length = 1.0;
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    verify_empty_tree_sequence(&ts, 1.0);

    ret = tsk_treeseq_kc_distance(&ts, &ts, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MULTIPLE_ROOTS);

    ret = tsk_tree_init(&t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_left_root(&t), TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(t.interval.left, 0);
    CU_ASSERT_EQUAL_FATAL(t.interval.right, 1);
    CU_ASSERT_EQUAL_FATAL(t.parent[0], TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(t.left_child[0], TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(t.right_child[0], TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(t.left_sib[0], TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(t.right_sib[0], TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_parent(&t, 1, &v), TSK_ERR_NODE_OUT_OF_BOUNDS);

    ret = tsk_tree_kc_distance(&t, &t, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MULTIPLE_ROOTS);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_nonbinary_tree_kc(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0";
    const char *edges = "0  1   4   0,1,2,3\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int ret;
    double result = 0;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);

    tsk_treeseq_kc_distance(&ts, &ts, 0, &result);
    CU_ASSERT_EQUAL_FATAL(result, 0);

    ret = tsk_tree_init(&t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    tsk_tree_kc_distance(&t, &t, 0, &result);
    CU_ASSERT_EQUAL_FATAL(result, 0);
    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_nonzero_samples_kc(void)
{
    const char *nodes = "0  0   0\n" /* unused node at the start */
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0";
    const char *edges = "0  1   3   1,2\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int ret;
    double result = 0;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);

    ret = tsk_treeseq_kc_distance(&ts, &ts, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(result, 0);

    ret = tsk_tree_init(&t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    ret = tsk_tree_kc_distance(&t, &t, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(result, 0);
    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_internal_samples_kc(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  1   0";
    const char *edges = "0  1   2   0,1\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int ret;
    double result = 0;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);

    /* Permitted in tree sequences */
    ret = tsk_treeseq_kc_distance(&ts, &ts, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(result, 0.0);

    ret = tsk_tree_init(&t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    ret = tsk_tree_kc_distance(&t, &t, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_non_sample_leaf_kc(void)
{
    const char *nodes = "1  0   0\n"
                        "0  0   0\n"
                        "0  1   0\n";
    const char *edges = "0 1  2 0,1\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int ret;
    double result = 0;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);

    ret = tsk_treeseq_kc_distance(&ts, &ts, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(result, 0.0);

    ret = tsk_tree_init(&t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    ret = tsk_tree_kc_distance(&t, &t, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(result, 0.0);

    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_unequal_sample_size_kc(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  2   0\n"
                        "0  3   0\n";
    const char *nodes_other = "1  0   0\n"
                              "1  0   0\n"
                              "0  1   0\n";
    const char *edges = "0 1  3 0,1\n"
                        "0 1  4 2,3\n";
    const char *edges_other = "0 1  2 0,1\n";
    int ret;
    tsk_treeseq_t ts, other_ts;
    tsk_tree_t t, other_t;
    double result = 0;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    tsk_treeseq_from_text(
        &other_ts, 1, nodes_other, edges_other, NULL, NULL, NULL, NULL, NULL, 0);

    ret = tsk_treeseq_kc_distance(&ts, &other_ts, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SAMPLE_SIZE_MISMATCH);

    ret = tsk_tree_init(&t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    ret = tsk_tree_init(&other_t, &other_ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&other_t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    ret = tsk_tree_kc_distance(&t, &other_t, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SAMPLE_SIZE_MISMATCH);
    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&other_ts);
    tsk_tree_free(&t);
    tsk_tree_free(&other_t);
}

static void
test_unequal_samples_kc(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  2   0\n"
                        "0  3   0\n";
    const char *nodes_other = "0  0   0\n" /* Unused node at the start */
                              "1  0   0\n"
                              "1  0   0\n"
                              "1  0   0\n"
                              "0  2   0\n"
                              "0  3   0\n";
    const char *edges = "0 1  3 0,1\n"
                        "0 1  4 2,3\n";
    const char *edges_other = "0 1  4 1,2\n"
                              "0 1  5 3,4\n";
    int ret;
    tsk_treeseq_t ts, other_ts;
    tsk_tree_t t, other_t;
    double result = 0;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    tsk_treeseq_from_text(
        &other_ts, 1, nodes_other, edges_other, NULL, NULL, NULL, NULL, NULL, 0);

    ret = tsk_treeseq_kc_distance(&ts, &other_ts, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SAMPLES_NOT_EQUAL);

    ret = tsk_tree_init(&t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    ret = tsk_tree_init(&other_t, &other_ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&other_t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    ret = tsk_tree_kc_distance(&t, &other_t, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SAMPLES_NOT_EQUAL);

    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&other_ts);
    tsk_tree_free(&t);
    tsk_tree_free(&other_t);
}

static void
test_unary_nodes_kc(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  2   0";
    const char *edges = "0  1   2   0,1\n"
                        "0  1   3   2";
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int ret;
    double result = 0;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&t, &ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    ret = tsk_tree_kc_distance(&t, &t, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNARY_NODES);

    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_no_sample_lists_kc(void)
{
    tsk_treeseq_t ts;
    tsk_tree_t t;
    int ret = 0;
    double result = 0;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    ret = tsk_tree_kc_distance(&t, &t, 9, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NO_SAMPLE_LISTS);

    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_unequal_sequence_lengths_kc(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  2   0\n"
                        "0  3   0\n";
    const char *edges_1 = "0 1  3 0,1\n"
                          "0 1  4 2,3\n";
    const char *edges_2 = "0 2  3 0,1\n"
                          "0 2  4 2,3\n";

    tsk_treeseq_t ts, other;
    int ret;
    double result = 0;

    tsk_treeseq_from_text(&ts, 1, nodes, edges_1, NULL, NULL, NULL, NULL, NULL, 0);
    tsk_treeseq_from_text(&other, 2, nodes, edges_2, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_treeseq_kc_distance(&ts, &other, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SEQUENCE_LENGTH_MISMATCH);

    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&other);
}

static void
test_different_number_trees_kc(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  1   0\n"
                        "0  2   0\n"
                        "0  3   0\n"
                        "0  4   0\n"
                        "0  5   0\n";
    const char *edges = "0 10  5 0,1\n"
                        "0 10  6 3,4\n"
                        "5 10  7 2,5\n"
                        "0 5   8 2\n"
                        "0 10  8 6\n"
                        "5 10  8 7\n"
                        "0 5   9 5,8\n";

    const char *other_nodes = "1  0   0\n"
                              "1  0   0\n"
                              "1  0   0\n"
                              "1  0   0\n"
                              "1  0   0\n"
                              "0  1   0\n"
                              "0  2   0\n"
                              "0  3   0\n"
                              "0  4   0\n";
    const char *other_edges = "0 10  5 0,1\n"
                              "0 10  6 2,3\n"
                              "0 10  7 4,5\n"
                              "0 10  8 6,7\n";
    tsk_treeseq_t ts, other;
    double result, expected;
    int ret = 0;

    tsk_treeseq_from_text(&ts, 10, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    tsk_treeseq_from_text(
        &other, 10, other_nodes, other_edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_treeseq_kc_distance(&ts, &other, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    expected = (sqrt(8.0) * 5.0 + sqrt(6.0) * 5.0) / 10.0;
    CU_ASSERT_DOUBLE_EQUAL_FATAL(result, expected, 1e-2);

    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&other);
}

static void
test_offset_trees_with_errors_kc(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "1  0   0\n"
                        "0  2   0\n"
                        "0  3   0\n"
                        "0  4   0\n";
    const char *edges = "0 10  4 0,1\n"
                        "0 10  5 2,3\n"
                        "0 10  6 4,5\n";
    tsk_treeseq_t ts, other;
    double result;
    int ret = 0;

    tsk_treeseq_from_text(
        &ts, 10, unary_ex_nodes, unary_ex_edges, NULL, NULL, NULL, NULL, NULL, 0);
    tsk_treeseq_from_text(&other, 10, nodes, edges, NULL, NULL, NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 10);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&other), 10);

    ret = tsk_treeseq_kc_distance(&ts, &other, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNARY_NODES);

    ret = tsk_treeseq_kc_distance(&other, &ts, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNARY_NODES);

    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&other);
}

/*=======================================================
 * Miscellaneous tests.
 *======================================================*/

static void
test_genealogical_nearest_neighbours_errors(void)
{
    int ret;
    tsk_treeseq_t ts;
    const tsk_id_t *reference_sets[2];
    tsk_id_t reference_set_0[4], reference_set_1[4];
    tsk_id_t focal[] = { 0, 1, 2, 3 };
    tsk_size_t reference_set_size[2];
    tsk_size_t num_focal = 4;
    double *A = tsk_malloc(2 * num_focal * sizeof(double));
    CU_ASSERT_FATAL(A != NULL);

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_treeseq_genealogical_nearest_neighbours(
        &ts, focal, num_focal, reference_sets, reference_set_size, 0, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_treeseq_genealogical_nearest_neighbours(
        &ts, focal, num_focal, reference_sets, reference_set_size, INT16_MAX, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Overlapping sample sets */
    reference_sets[0] = focal;
    reference_set_size[0] = 1;
    reference_sets[1] = focal;
    reference_set_size[1] = num_focal;
    ret = tsk_treeseq_genealogical_nearest_neighbours(
        &ts, focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SAMPLE);

    /* bad values in the sample sets */
    reference_set_0[0] = 0;
    reference_set_0[1] = 1;
    reference_set_1[0] = 2;
    reference_set_1[1] = 3;
    reference_set_size[0] = 2;
    reference_set_size[1] = 2;
    reference_sets[0] = reference_set_0;
    reference_sets[1] = reference_set_1;
    ret = tsk_treeseq_genealogical_nearest_neighbours(
        &ts, focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    reference_set_0[0] = -1;
    ret = tsk_treeseq_genealogical_nearest_neighbours(
        &ts, focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    reference_set_0[0] = (tsk_id_t) tsk_treeseq_get_num_nodes(&ts);
    ret = tsk_treeseq_genealogical_nearest_neighbours(
        &ts, focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    reference_set_0[0] = (tsk_id_t) tsk_treeseq_get_num_nodes(&ts) + 1;
    ret = tsk_treeseq_genealogical_nearest_neighbours(
        &ts, focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    /* Duplicate values in the focal sets */
    reference_set_0[0] = 1;
    ret = tsk_treeseq_genealogical_nearest_neighbours(
        &ts, focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SAMPLE);
    reference_set_0[0] = 3;
    ret = tsk_treeseq_genealogical_nearest_neighbours(
        &ts, focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SAMPLE);

    /* Bad sample ID */
    reference_sets[0] = focal;
    reference_set_size[0] = 1;
    reference_sets[1] = focal + 1;
    reference_set_size[1] = num_focal - 1;
    focal[0] = -1;
    ret = tsk_treeseq_genealogical_nearest_neighbours(
        &ts, focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    focal[0] = (tsk_id_t) tsk_treeseq_get_num_nodes(&ts);
    ret = tsk_treeseq_genealogical_nearest_neighbours(
        &ts, focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    focal[0] = (tsk_id_t) tsk_treeseq_get_num_nodes(&ts) + 100;
    ret = tsk_treeseq_genealogical_nearest_neighbours(
        &ts, focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    tsk_treeseq_free(&ts);
    free(A);
}

static void
test_single_tree_balance(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t t;
    tsk_size_t sackin;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    /* Balanced binary tree with 4 leaves */
    CU_ASSERT_EQUAL_FATAL(tsk_tree_sackin_index(&t, &sackin), 0);
    CU_ASSERT_EQUAL(sackin, 8);

    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_multiroot_balance(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t t;
    tsk_size_t sackin;

    tsk_treeseq_from_text(&ts, 10, multiroot_ex_nodes, multiroot_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    /* 0.80┊         10 */
    /*     ┊         ┏┻┓ */
    /* 0.40┊     9   ┃ ┃ */
    /*     ┊   ┏━┻┓  ┃ ┃ */
    /* 0.30┊   ┃  ┃  ┃ ┃ */
    /*     ┊   ┃  ┃  ┃ ┃ */
    /* 0.20┊   ┃  7  ┃ ┃ */
    /*     ┊   ┃ ┏┻┓ ┃ ┃ */
    /* 0.10┊   ┃ ┃ ┃ ┃ ┃ */
    /*     ┊   ┃ ┃ ┃ ┃ ┃ */
    /* 0.00┊ 5 2 3 4 0 1 */

    CU_ASSERT_EQUAL_FATAL(tsk_tree_sackin_index(&t, &sackin), 0);
    CU_ASSERT_EQUAL(sackin, 7);

    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_empty_tree_balance(void)
{
    int ret;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    tsk_tree_t t;
    tsk_size_t sackin;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1.0;
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    CU_ASSERT_EQUAL_FATAL(tsk_tree_sackin_index(&t, &sackin), 0);
    CU_ASSERT_EQUAL(sackin, 0);

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_tree_errors(void)
{
    int ret;
    tsk_size_t j;
    tsk_id_t num_nodes = 9;
    tsk_id_t u;
    tsk_node_t node;
    tsk_treeseq_t ts, other_ts;
    tsk_tree_t t, other_t;
    tsk_id_t bad_nodes[] = { num_nodes + 1, num_nodes + 2, -1 };
    tsk_id_t tracked_samples[] = { 0, 0, 0 };

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, NULL, NULL,
        paper_ex_individuals, NULL, 0);

    ret = tsk_tree_init(&t, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);

    /* Out-of-bounds queries */
    for (j = 0; j < sizeof(bad_nodes) / sizeof(tsk_id_t); j++) {
        u = bad_nodes[j];
        ret = tsk_tree_get_parent(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        ret = tsk_tree_get_time(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        ret = tsk_tree_get_branch_length(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        ret = tsk_tree_get_mrca(&t, u, 0, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        ret = tsk_tree_get_mrca(&t, 0, u, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        ret = tsk_tree_get_num_samples(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        ret = tsk_tree_get_num_tracked_samples(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        /* Also check tree sequence methods */
        ret = tsk_treeseq_get_node(&ts, (tsk_id_t) u, &node);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        CU_ASSERT(!tsk_treeseq_is_sample(&ts, u));
        CU_ASSERT(!tsk_tree_is_sample(&t, u));
    }

    tracked_samples[0] = 0;
    tracked_samples[1] = (tsk_id_t) tsk_treeseq_get_num_samples(&ts);
    ret = tsk_tree_set_tracked_samples(&t, 2, tracked_samples);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_SAMPLES);
    tracked_samples[1] = (tsk_id_t) tsk_treeseq_get_num_nodes(&ts);
    ret = tsk_tree_set_tracked_samples(&t, 2, tracked_samples);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tracked_samples[1] = 0;
    ret = tsk_tree_set_tracked_samples(&t, 2, tracked_samples);
    CU_ASSERT_EQUAL(ret, TSK_ERR_DUPLICATE_SAMPLE);

    tsk_treeseq_from_text(&other_ts, 10, paper_ex_nodes, paper_ex_edges, NULL, NULL,
        NULL, paper_ex_individuals, NULL, 0);

    ret = tsk_tree_init(&other_t, &other_ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_copy(&t, &other_t, TSK_NO_INIT);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    tsk_tree_free(&t);
    tsk_tree_free(&other_t);

    ret = tsk_tree_init(&t, &other_ts, TSK_NO_SAMPLE_COUNTS);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_copy(&t, &other_t, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_UNSUPPORTED_OPERATION);
    tsk_tree_free(&other_t);
    ret = tsk_tree_copy(&t, &other_t, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL(ret, TSK_ERR_UNSUPPORTED_OPERATION);
    tsk_tree_free(&other_t);

    tsk_tree_free(&t);
    tsk_treeseq_free(&other_ts);
    tsk_treeseq_free(&ts);
}

static void
test_treeseq_row_access_errors(void)
{
    int ret;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_get_individual(&ts, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    ret = tsk_treeseq_get_node(&ts, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    ret = tsk_treeseq_get_edge(&ts, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);
    ret = tsk_treeseq_get_migration(&ts, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATION_OUT_OF_BOUNDS);
    ret = tsk_treeseq_get_site(&ts, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    ret = tsk_treeseq_get_mutation(&ts, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);
    ret = tsk_treeseq_get_population(&ts, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    ret = tsk_treeseq_get_provenance(&ts, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_PROVENANCE_OUT_OF_BOUNDS);

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_tree_copy_flags(void)
{
    int iret, ret;
    tsk_size_t j;
    tsk_treeseq_t ts;
    tsk_tree_t t, other_t;
    tsk_flags_t options[] = { 0, TSK_NO_SAMPLE_COUNTS, TSK_SAMPLE_LISTS,
        TSK_NO_SAMPLE_COUNTS | TSK_SAMPLE_LISTS };

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, NULL, NULL,
        paper_ex_individuals, NULL, 0);

    for (j = 0; j < sizeof(options) / sizeof(*options); j++) {
        ret = tsk_tree_init(&t, &ts, options[j]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tree_init(&other_t, &ts, options[j]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tree_copy(&t, &other_t, TSK_NO_INIT);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        check_trees_identical(&t, &other_t);
        tsk_tree_free(&other_t);

        while ((iret = tsk_tree_next(&t)) == TSK_TREE_OK) {
            ret = tsk_tree_copy(&t, &other_t, options[j]);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            check_trees_identical(&t, &other_t);
            tsk_tree_free(&other_t);
        }
        CU_ASSERT_EQUAL_FATAL(iret, 0);

        ret = tsk_tree_first(&t);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
        ret = tsk_tree_copy(&t, &other_t, options[j]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        while (true) {
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            check_trees_identical(&t, &other_t);
            CU_ASSERT_EQUAL_FATAL(tsk_tree_next(&t), tsk_tree_next(&other_t));
            if (t.index == -1) {
                break;
            }
        }

        ret = tsk_tree_last(&t);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
        ret = tsk_tree_copy(&t, &other_t, TSK_NO_INIT | options[j]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        while (true) {
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            check_trees_identical(&t, &other_t);
            CU_ASSERT_EQUAL_FATAL(tsk_tree_prev(&t), tsk_tree_prev(&other_t));
            if (t.index == -1) {
                break;
            }
        }

        tsk_tree_free(&other_t);
        tsk_tree_free(&t);
    }
    tsk_treeseq_free(&ts);
}

static void
test_deduplicate_sites(void)
{
    int ret;
    // Modified from paper_ex
    const char *tidy_sites = "1      0\n"
                             "4.5    0\n"
                             "8.5    0\n";
    const char *tidy_mutations = "0      2   1\n"
                                 "0      1   2\n"
                                 "0      6   3\n"
                                 "0      3   4\n"
                                 "1      0   1\n"
                                 "1      2   2\n"
                                 "1      4   3\n"
                                 "1      5   4\n"
                                 "2      5   1\n"
                                 "2      7   2\n"
                                 "2      1   3\n"
                                 "2      0   4\n";
    const char *messy_sites = "1      0\n"
                              "1      0\n"
                              "1      0\n"
                              "1      0\n"
                              "4.5    0\n"
                              "4.5    0\n"
                              "4.5    0\n"
                              "4.5    0\n"
                              "8.5    0\n"
                              "8.5    0\n"
                              "8.5    0\n"
                              "8.5    0\n";
    const char *messy_mutations = "0      2   1\n"
                                  "1      1   2\n"
                                  "2      6   3\n"
                                  "3      3   4\n"
                                  "4      0   1\n"
                                  "5      2   2\n"
                                  "6      4   3\n"
                                  "7      5   4\n"
                                  "8      5   1\n"
                                  "9      7   2\n"
                                  "10     1   3\n"
                                  "11     0   4\n";
    tsk_table_collection_t tidy, messy;

    ret = tsk_table_collection_init(&tidy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&messy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    messy.sequence_length = 10;
    tidy.sequence_length = 10;
    parse_individuals(paper_ex_individuals, &tidy.individuals);
    parse_nodes(paper_ex_nodes, &tidy.nodes);
    parse_sites(tidy_sites, &tidy.sites);
    parse_mutations(tidy_mutations, &tidy.mutations);
    // test cleaning doesn't mess up the tidy one
    parse_individuals(paper_ex_individuals, &messy.individuals);
    parse_nodes(paper_ex_nodes, &messy.nodes);
    parse_sites(tidy_sites, &messy.sites);
    parse_mutations(tidy_mutations, &messy.mutations);

    ret = tsk_table_collection_deduplicate_sites(&messy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_site_table_equals(&tidy.sites, &messy.sites, 0));
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&tidy.mutations, &messy.mutations, 0));

    tsk_site_table_clear(&messy.sites);
    tsk_mutation_table_clear(&messy.mutations);

    // test with the actual messy one
    parse_sites(messy_sites, &messy.sites);
    parse_mutations(messy_mutations, &messy.mutations);

    ret = tsk_table_collection_deduplicate_sites(&messy, 0);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(tsk_site_table_equals(&tidy.sites, &messy.sites, 0));
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&tidy.mutations, &messy.mutations, 0));

    tsk_table_collection_free(&tidy);
    tsk_table_collection_free(&messy);
}

static void
test_deduplicate_sites_errors(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 10;
    ret_id = tsk_site_table_add_row(&tables.sites, 2, "A", 1, "m", 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 2, "TT", 2, "MM", 2);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    ret_id = tsk_mutation_table_add_row(&tables.mutations, 0, 0, -1, 0, "T", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    /* Negative position */
    tables.sites.position[0] = -1;
    ret = tsk_table_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_SITE_POSITION);
    tables.sites.position[0] = 2;

    /* unsorted position */
    tables.sites.position[1] = 0.5;
    ret = tsk_table_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_UNSORTED_SITES);
    tables.sites.position[1] = 2;

    /* negative site ID */
    tables.mutations.site[0] = -1;
    ret = tsk_table_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tables.mutations.site[0] = 0;

    /* site ID out of bounds */
    tables.mutations.site[0] = 2;
    ret = tsk_table_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tables.mutations.site[0] = 0;

    /* Bad offset in metadata */
    tables.sites.metadata_offset[0] = 2;
    ret = tsk_table_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    tables.sites.metadata_offset[0] = 0;

    /* Bad length in metadata */
    tables.sites.metadata_offset[2] = 100;
    ret = tsk_table_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    tables.sites.metadata_offset[2] = 3;

    /* Bad offset in ancestral_state */
    tables.sites.ancestral_state_offset[0] = 2;
    ret = tsk_table_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    tables.sites.ancestral_state_offset[0] = 0;

    /* Bad length in ancestral_state */
    tables.sites.ancestral_state_offset[2] = 100;
    ret = tsk_table_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    tables.sites.ancestral_state_offset[2] = 3;

    ret = tsk_table_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, 0);

    tsk_table_collection_free(&tables);
}

static void
test_deduplicate_sites_zero_rows(void)
{

    int ret;
    tsk_table_collection_t tables;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    ret = tsk_table_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tables.sites.num_rows, 0)

    tsk_table_collection_free(&tables);
}

static void
test_deduplicate_sites_multichar(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 10;
    ret_id = tsk_site_table_add_row(&tables.sites, 0, "AA", 1, "M", 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0, "0", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    ret_id = tsk_site_table_add_row(&tables.sites, 1, "BBBBB", 5, "NNNNN", 5);
    CU_ASSERT_EQUAL_FATAL(ret_id, 2);
    ret_id = tsk_site_table_add_row(&tables.sites, 1, "0", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 3);

    ret = tsk_table_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.num_rows, 2);
    CU_ASSERT_EQUAL_FATAL(tables.sites.position[0], 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.position[1], 1);
    CU_ASSERT_EQUAL_FATAL(tables.sites.ancestral_state[0], 'A');
    CU_ASSERT_EQUAL_FATAL(tables.sites.ancestral_state_offset[1], 1);
    CU_ASSERT_EQUAL_FATAL(tables.sites.metadata[0], 'M');
    CU_ASSERT_EQUAL_FATAL(tables.sites.metadata_offset[1], 1);

    CU_ASSERT_NSTRING_EQUAL(tables.sites.ancestral_state + 1, "BBBBB", 5);
    CU_ASSERT_EQUAL_FATAL(tables.sites.ancestral_state_offset[2], 6);
    CU_ASSERT_NSTRING_EQUAL(tables.sites.metadata + 1, "NNNNN", 5);
    CU_ASSERT_EQUAL_FATAL(tables.sites.metadata_offset[2], 6);

    tsk_table_collection_free(&tables);
}

static void
test_empty_tree_sequence(void)
{
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_tree_t t;
    tsk_id_t v;
    int ret;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_SEQUENCE_LENGTH);
    tsk_treeseq_free(&ts);
    tables.sequence_length = 1.0;
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    verify_empty_tree_sequence(&ts, 1.0);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_left_root(&t), TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(t.interval.left, 0);
    CU_ASSERT_EQUAL_FATAL(t.interval.right, 1);
    CU_ASSERT_EQUAL_FATAL(t.num_edges, 0);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_parent(&t, 0, &v), 0);
    CU_ASSERT_EQUAL_FATAL(v, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_parent(&t, 1, &v), TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_tree_free(&t);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_last(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_left_root(&t), TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(t.interval.left, 0);
    CU_ASSERT_EQUAL_FATAL(t.interval.right, 1);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_parent(&t, 1, &v), TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_tree_free(&t);

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tables);
}

static void
test_zero_edges(void)
{
    const char *nodes = "1  0   0\n"
                        "1  0   0\n";
    const char *edges = "";
    const char *sites = "0.1  0\n"
                        "0.2  0\n";
    const char *mutations = "0    0     1\n"
                            "1    1     1\n";
    tsk_treeseq_t ts, tss;
    tsk_tree_t t;
    tsk_id_t samples, node_map;
    const tsk_id_t z = TSK_NULL;
    tsk_id_t parents[] = {
        z,
        z,
    };
    int ret;

    tsk_treeseq_from_text(&ts, 2, nodes, edges, NULL, sites, mutations, NULL, NULL, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 2.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_treeseq_print_state(&ts, _devnull);

    verify_trees(&ts, 1, parents);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(t.interval.left, 0);
    CU_ASSERT_EQUAL(t.interval.right, 2);
    CU_ASSERT_EQUAL(t.num_edges, 0);
    CU_ASSERT_EQUAL(t.parent[0], TSK_NULL);
    CU_ASSERT_EQUAL(t.parent[1], TSK_NULL);
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&t), 0);
    CU_ASSERT_EQUAL(t.left_sib[0], TSK_NULL);
    CU_ASSERT_EQUAL(t.right_sib[0], 1);
    tsk_tree_print_state(&t, _devnull);
    tsk_tree_free(&t);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_last(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK);
    CU_ASSERT_EQUAL(t.interval.left, 0);
    CU_ASSERT_EQUAL(t.interval.right, 2);
    CU_ASSERT_EQUAL(t.parent[0], TSK_NULL);
    CU_ASSERT_EQUAL(t.parent[1], TSK_NULL);
    CU_ASSERT_EQUAL(tsk_tree_get_left_root(&t), 0);
    CU_ASSERT_EQUAL(t.left_sib[0], TSK_NULL);
    CU_ASSERT_EQUAL(t.right_sib[0], 1);
    tsk_tree_print_state(&t, _devnull);
    tsk_tree_free(&t);

    /* We give pointers ot samples and node_map here as they must be non null */
    ret = tsk_treeseq_simplify(&ts, &samples, 0, 0, &tss, &node_map);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&tss), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&tss), 2.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&tss), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&tss), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&tss), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&tss), 1);
    tsk_treeseq_print_state(&ts, _devnull);

    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&tss);
}

static void
test_tree_sequence_metadata(void)
{
    int ret;
    tsk_table_collection_t tc;
    tsk_treeseq_t ts;

    char example_metadata[100] = "An example of metadata with unicode 🎄🌳🌴🌲🎋";
    char example_metadata_schema[100]
        = "An example of metadata schema with unicode 🎄🌳🌴🌲🎋";
    char example_time_units[100] = "An example of time units ⏰";
    tsk_size_t example_metadata_length = (tsk_size_t) strlen(example_metadata);
    tsk_size_t example_time_units_length = (tsk_size_t) strlen(example_metadata_schema);
    tsk_size_t example_metadata_schema_length = (tsk_size_t) strlen(example_time_units);

    ret = tsk_table_collection_init(&tc, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tc.sequence_length = 1.0;
    ret = tsk_table_collection_build_index(&tc, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_set_metadata(
        &tc, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_set_metadata_schema(
        &tc, example_metadata_schema, example_metadata_schema_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_set_time_units(
        &tc, example_time_units, example_time_units_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_init(&ts, &tc, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL(tsk_treeseq_get_metadata_length(&ts), example_metadata_length);
    CU_ASSERT_EQUAL(
        tsk_treeseq_get_metadata_schema_length(&ts), example_metadata_schema_length);
    CU_ASSERT_EQUAL(tsk_memcmp(tsk_treeseq_get_metadata(&ts), example_metadata,
                        example_metadata_length),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(tsk_treeseq_get_metadata_schema(&ts),
                        example_metadata_schema, example_metadata_schema_length),
        0);

    CU_ASSERT_EQUAL(tsk_treeseq_get_time_units_length(&ts), example_time_units_length);
    CU_ASSERT_EQUAL(tsk_memcmp(tsk_treeseq_get_time_units(&ts), example_time_units,
                        example_time_units_length),
        0);

    tsk_treeseq_free(&ts);
    tsk_table_collection_free(&tc);
}

static int
dummy_stat(tsk_size_t K, const double *X, tsk_size_t M, double *Y, void *params)
{
    tsk_size_t k;
    CU_ASSERT_FATAL(M == K);
    CU_ASSERT_FATAL(params == NULL);

    for (k = 0; k < K; k++) {
        Y[k] = X[k];
    }
    return 0;
}

static void
test_time_uncalibrated(void)
{
    int ret;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    tsk_treeseq_t ts2;
    tsk_size_t sample_set_sizes[] = { 2, 2 };
    tsk_id_t samples[] = { 0, 1, 2, 3 };
    tsk_size_t num_samples;
    double result[10];
    double *W;
    double *sigma;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ts.time_uncalibrated, false);
    tsk_treeseq_free(&ts);

    ret = tsk_table_collection_set_time_units(
        &tables, TSK_TIME_UNITS_UNCALIBRATED, strlen(TSK_TIME_UNITS_UNCALIBRATED));
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ts.time_uncalibrated, true);
    tsk_treeseq_free(&ts);

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, paper_ex_sites,
        paper_ex_mutations, paper_ex_individuals, NULL, 0);
    ret = tsk_table_collection_set_time_units(
        ts.tables, TSK_TIME_UNITS_UNCALIBRATED, strlen(TSK_TIME_UNITS_UNCALIBRATED));
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_init(&ts2, ts.tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_allele_frequency_spectrum(
        &ts2, 2, sample_set_sizes, samples, 0, NULL, TSK_STAT_SITE, result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_allele_frequency_spectrum(
        &ts2, 2, sample_set_sizes, samples, 0, NULL, TSK_STAT_BRANCH, result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TIME_UNCALIBRATED);
    ret = tsk_treeseq_allele_frequency_spectrum(&ts2, 2, sample_set_sizes, samples, 0,
        NULL, TSK_STAT_BRANCH | TSK_STAT_ALLOW_TIME_UNCALIBRATED, result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    sigma = tsk_calloc(tsk_treeseq_get_num_nodes(&ts2), sizeof(double));
    num_samples = tsk_treeseq_get_num_samples(&ts2);
    W = tsk_calloc(num_samples, sizeof(double));

    ret = tsk_treeseq_general_stat(&ts2, 1, W, 1, dummy_stat, NULL,
        tsk_treeseq_get_num_trees(&ts2), tsk_treeseq_get_breakpoints(&ts2),
        TSK_STAT_SITE, sigma);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_general_stat(&ts2, 1, W, 1, dummy_stat, NULL,
        tsk_treeseq_get_num_trees(&ts2), tsk_treeseq_get_breakpoints(&ts2),
        TSK_STAT_BRANCH, sigma);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TIME_UNCALIBRATED);
    ret = tsk_treeseq_general_stat(&ts2, 1, W, 1, dummy_stat, NULL,
        tsk_treeseq_get_num_trees(&ts2), tsk_treeseq_get_breakpoints(&ts2),
        TSK_STAT_BRANCH | TSK_STAT_ALLOW_TIME_UNCALIBRATED, sigma);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_safe_free(W);
    tsk_safe_free(sigma);
    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&ts2);
    tsk_table_collection_free(&tables);
}

static void
test_reference_sequence(void)
{
    int ret;
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_treeseq_has_reference_sequence(&ts));
    tsk_treeseq_free(&ts);

    ret = tsk_reference_sequence_set_data(&tables.reference_sequence, "abc", 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_treeseq_has_reference_sequence(&ts));
    tsk_treeseq_free(&ts);

    tsk_table_collection_free(&tables);
}

static void
test_init_take_ownership_no_edge_metadata(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t *tables = tsk_malloc(sizeof(tsk_table_collection_t));

    CU_ASSERT_NOT_EQUAL_FATAL(tables, NULL);

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, paper_ex_sites,
        paper_ex_mutations, paper_ex_individuals, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, tables, TSK_TC_NO_EDGE_METADATA);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);

    ret = tsk_treeseq_init(&ts, tables, TSK_TAKE_OWNERSHIP);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_CANT_TAKE_OWNERSHIP_NO_EDGE_METADATA);

    tsk_treeseq_free(&ts);
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        /* simplest example tests */
        { "test_simplest_discrete_genome", test_simplest_discrete_genome },
        { "test_simplest_discrete_time", test_simplest_discrete_time },
        { "test_simplest_records", test_simplest_records },
        { "test_simplest_nonbinary_records", test_simplest_nonbinary_records },
        { "test_simplest_unary_records", test_simplest_unary_records },
        { "test_simplest_unary_with_individuals", test_simplest_unary_with_individuals },
        { "test_simplest_non_sample_leaf_records",
            test_simplest_non_sample_leaf_records },
        { "test_simplest_degenerate_multiple_root_records",
            test_simplest_degenerate_multiple_root_records },
        { "test_simplest_multiple_root_records", test_simplest_multiple_root_records },
        { "test_simplest_zero_root_tree", test_simplest_zero_root_tree },
        { "test_simplest_multi_root_tree", test_simplest_multi_root_tree },
        { "test_simplest_tree_mrca", test_simplest_tree_mrca },
        { "test_simplest_root_mutations", test_simplest_root_mutations },
        { "test_simplest_back_mutations", test_simplest_back_mutations },
        { "test_simplest_general_samples", test_simplest_general_samples },
        { "test_simplest_holey_tree_sequence", test_simplest_holey_tree_sequence },
        { "test_simplest_holey_tsk_treeseq_zero_roots",
            test_simplest_holey_tsk_treeseq_zero_roots },
        { "test_simplest_holey_tsk_treeseq_mutation_parents",
            test_simplest_holey_tsk_treeseq_mutation_parents },
        { "test_simplest_initial_gap_tree_sequence",
            test_simplest_initial_gap_tree_sequence },
        { "test_simplest_initial_gap_zero_roots", test_simplest_initial_gap_zero_roots },
        { "test_simplest_initial_gap_tsk_treeseq_mutation_parents",
            test_simplest_initial_gap_tsk_treeseq_mutation_parents },
        { "test_simplest_final_gap_tree_sequence",
            test_simplest_final_gap_tree_sequence },
        { "test_simplest_final_gap_tsk_treeseq_mutation_parents",
            test_simplest_final_gap_tsk_treeseq_mutation_parents },
        { "test_simplest_individuals", test_simplest_individuals },
        { "test_simplest_bad_individuals", test_simplest_bad_individuals },
        { "test_simplest_bad_edges", test_simplest_bad_edges },
        { "test_simplest_bad_indexes", test_simplest_bad_indexes },
        { "test_simplest_bad_migrations", test_simplest_bad_migrations },
        { "test_simplest_migration_simplify", test_simplest_migration_simplify },
        { "test_simplest_overlapping_parents", test_simplest_overlapping_parents },
        { "test_simplest_contradictory_children", test_simplest_contradictory_children },
        { "test_simplest_overlapping_edges_simplify",
            test_simplest_overlapping_edges_simplify },
        { "test_simplest_overlapping_unary_edges_simplify",
            test_simplest_overlapping_unary_edges_simplify },
        { "test_simplest_overlapping_unary_edges_internal_samples_simplify",
            test_simplest_overlapping_unary_edges_internal_samples_simplify },
        { "test_simplest_reduce_site_topology", test_simplest_reduce_site_topology },
        { "test_simplest_simplify_defragment", test_simplest_simplify_defragment },
        { "test_simplest_population_filter", test_simplest_population_filter },
        { "test_simplest_individual_filter", test_simplest_individual_filter },
        { "test_simplest_map_mutations", test_simplest_map_mutations },
        { "test_simplest_nonbinary_map_mutations",
            test_simplest_nonbinary_map_mutations },
        { "test_simplest_unary_map_mutations", test_simplest_unary_map_mutations },
        { "test_simplest_non_sample_leaf_map_mutations",
            test_simplest_non_sample_leaf_map_mutations },
        { "test_simplest_internal_sample_map_mutations",
            test_simplest_internal_sample_map_mutations },
        { "test_simplest_multiple_root_map_mutations",
            test_simplest_multiple_root_map_mutations },
        { "test_simplest_chained_map_mutations", test_simplest_chained_map_mutations },
        { "test_simplest_mutation_edges", test_simplest_mutation_edges },

        /* Single tree tests */
        { "test_single_tree_good_records", test_single_tree_good_records },
        { "test_single_nonbinary_tree_good_records",
            test_single_nonbinary_tree_good_records },
        { "test_single_tree_bad_records", test_single_tree_bad_records },
        { "test_single_tree_good_mutations", test_single_tree_good_mutations },
        { "test_single_tree_bad_mutations", test_single_tree_bad_mutations },
        { "test_single_tree_iter", test_single_tree_iter },
        { "test_single_tree_general_samples_iter",
            test_single_tree_general_samples_iter },
        { "test_single_nonbinary_tree_iter", test_single_nonbinary_tree_iter },
        { "test_single_tree_iter_times", test_single_tree_iter_times },
        { "test_single_tree_iter_depths", test_single_tree_iter_depths },
        { "test_single_tree_simplify", test_single_tree_simplify },
        { "test_single_tree_simplify_debug", test_single_tree_simplify_debug },
        { "test_single_tree_simplify_keep_input_roots",
            test_single_tree_simplify_keep_input_roots },
        { "test_single_tree_simplify_no_sample_nodes",
            test_single_tree_simplify_no_sample_nodes },
        { "test_single_tree_simplify_null_samples",
            test_single_tree_simplify_null_samples },
        { "test_single_tree_compute_mutation_parents",
            test_single_tree_compute_mutation_parents },
        { "test_single_tree_compute_mutation_times",
            test_single_tree_compute_mutation_times },
        { "test_single_tree_mutation_edges", test_single_tree_mutation_edges },
        { "test_single_tree_is_descendant", test_single_tree_is_descendant },
        { "test_single_tree_total_branch_length", test_single_tree_total_branch_length },
        { "test_single_tree_map_mutations", test_single_tree_map_mutations },
        { "test_single_tree_map_mutations_internal_samples",
            test_single_tree_map_mutations_internal_samples },
        { "test_single_tree_tracked_samples", test_single_tree_tracked_samples },

        /* Multi tree tests */
        { "test_simple_multi_tree", test_simple_multi_tree },
        { "test_nonbinary_multi_tree", test_nonbinary_multi_tree },
        { "test_unary_multi_tree", test_unary_multi_tree },
        { "test_internal_sample_multi_tree", test_internal_sample_multi_tree },
        { "test_internal_sample_simplified_multi_tree",
            test_internal_sample_simplified_multi_tree },
        { "test_simplify_keep_input_roots_multi_tree",
            test_simplify_keep_input_roots_multi_tree },
        { "test_left_to_right_multi_tree", test_left_to_right_multi_tree },
        { "test_gappy_multi_tree", test_gappy_multi_tree },

        { "test_tsk_treeseq_bad_records", test_tsk_treeseq_bad_records },

        /* multiroot tests */
        { "test_multiroot_mrca", test_multiroot_mrca },
        { "test_multiroot_diff_iter", test_multiroot_diff_iter },

        /* Diff iter tests */
        { "test_simple_diff_iter", test_simple_diff_iter },
        { "test_nonbinary_diff_iter", test_nonbinary_diff_iter },
        { "test_unary_diff_iter", test_unary_diff_iter },
        { "test_internal_sample_diff_iter", test_internal_sample_diff_iter },
        { "test_empty_diff_iter", test_empty_diff_iter },

        /* Sample sets */
        { "test_simple_sample_sets", test_simple_sample_sets },
        { "test_nonbinary_sample_sets", test_nonbinary_sample_sets },
        { "test_internal_sample_sample_sets", test_internal_sample_sample_sets },
        { "test_non_sample_leaf_sample_lists", test_non_sample_leaf_sample_lists },

        { "test_no_sample_count_semantics", test_no_sample_count_semantics },
        { "test_virtual_root_properties", test_virtual_root_properties },

        /* tree traversal orders */
        { "test_single_tree_traversal", test_single_tree_traversal },
        { "test_multiroot_tree_traversal", test_multiroot_tree_traversal },

        /* Seek */
        { "test_seek_multi_tree", test_seek_multi_tree },
        { "test_seek_errors", test_seek_errors },

        /* KC distance tests */
        { "test_single_tree_kc", test_single_tree_kc },
        { "test_isolated_node_kc", test_isolated_node_kc },
        { "test_two_trees_kc", test_two_trees_kc },
        { "test_empty_tree_kc", test_empty_tree_kc },
        { "test_nonbinary_tree_kc", test_nonbinary_tree_kc },
        { "test_nonzero_samples_kc", test_nonzero_samples_kc },
        { "test_internal_samples_kc", test_internal_samples_kc },
        { "test_non_sample_leaf_kc", test_non_sample_leaf_kc },
        { "test_unequal_sample_size_kc", test_unequal_sample_size_kc },
        { "test_unequal_samples_kc", test_unequal_samples_kc },
        { "test_unary_nodes_kc", test_unary_nodes_kc },
        { "test_no_sample_lists_kc", test_no_sample_lists_kc },
        { "test_unequal_sequence_lengths_kc", test_unequal_sequence_lengths_kc },
        { "test_different_number_trees_kc", test_different_number_trees_kc },
        { "test_offset_trees_with_errors_kc", test_offset_trees_with_errors_kc },

        /* Tree balance/imbalance index tests */
        { "test_single_tree_balance", test_single_tree_balance },
        { "test_multiroot_balance", test_multiroot_balance },
        { "test_empty_tree_balance", test_empty_tree_balance },

        /* Misc */
        { "test_tree_errors", test_tree_errors },
        { "test_treeseq_row_access_errors", test_treeseq_row_access_errors },
        { "test_tree_copy_flags", test_tree_copy_flags },
        { "test_genealogical_nearest_neighbours_errors",
            test_genealogical_nearest_neighbours_errors },
        { "test_deduplicate_sites", test_deduplicate_sites },
        { "test_deduplicate_sites_errors", test_deduplicate_sites_errors },
        { "test_deduplicate_sites_zero_rows", test_deduplicate_sites_zero_rows },
        { "test_deduplicate_sites_multichar", test_deduplicate_sites_multichar },
        { "test_empty_tree_sequence", test_empty_tree_sequence },
        { "test_zero_edges", test_zero_edges },
        { "test_tree_sequence_metadata", test_tree_sequence_metadata },
        { "test_time_uncalibrated", test_time_uncalibrated },
        { "test_reference_sequence", test_reference_sequence },
        { "test_init_take_ownership_no_edge_metadata",
            test_init_take_ownership_no_edge_metadata },
        { NULL, NULL },
    };

    return test_main(tests, argc, argv);
}
