/*
** Copyright (C) 2014 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#include "err.h"
#include "object_heap.h"
#include "msprime.h"

static int
cmp_hapgen_tree_node(const void *a, const void *b) {
    const hapgen_tree_node_t *ia = (const hapgen_tree_node_t *) a;
    const hapgen_tree_node_t *ib = (const hapgen_tree_node_t *) b;
    return (ia->id > ib->id) - (ia->id < ib->id);
}

/* Ensure the tree is in a consistent state */
static void
hapgen_check_state(hapgen_t *self)
{
    uint32_t j;
    avl_node_t *avl_node;
    hapgen_tree_node_t *node = NULL;
    hapgen_tree_node_t search;
    double tbl = 0.0;

    assert(avl_count(&self->tree) == 2 * self->sample_size - 1);
    for (j = 1; j <= self->sample_size; j++) {
        search.id = j;
        avl_node = avl_search(&self->tree, &search);
        assert(avl_node != NULL);
        node = (hapgen_tree_node_t *) avl_node->item;
        while (node->parent != NULL) {
            node = node->parent;
        }
        assert(node != NULL);
        assert(node == self->root);
    }
    for (avl_node = self->tree.head; avl_node != NULL; avl_node = avl_node->next) {
        node = (hapgen_tree_node_t *) avl_node->item;
        tbl += node->branch_length;
    }
    if (gsl_fcmp(tbl, self->total_branch_length, 1e-8) != 0) {
        printf("tbl %f\t%f\n", tbl, self->total_branch_length);
        assert(tbl == self->total_branch_length);
    }
}

void
hapgen_print_state(hapgen_t *self)
{
    avl_node_t *avl_node;
    hapgen_tree_node_t *node;
    size_t j;

    printf("Hapgen state\n");
    printf("total branch length = %f\n", self->total_branch_length);
    printf("segregating sites = %d\n", self->num_segregating_sites);
    printf("num_nodes = %d\n", avl_count(&self->tree));
    printf("root = %d\n", self->root == NULL? 0: self->root->id);
    for (avl_node = self->tree.head; avl_node != NULL; avl_node = avl_node->next) {
        node = (hapgen_tree_node_t *) avl_node->item;
        printf("%d\t", node->id);
        for (j = 0; j < 2; j++) {
            printf("%d\t", node->children[j] == NULL? 0 :
                    node->children[j]->id);
        }
        printf("%d\t%f\t%f\n",
                node->parent == NULL ? 0: node->parent->id, node->time,
                node->branch_length);
    }
    printf("avl_node_heap\n");
    object_heap_print_state(&self->avl_node_heap);
    hapgen_check_state(self);
}

static inline avl_node_t * WARN_UNUSED
hapgen_alloc_avl_node(hapgen_t *self, uint32_t node_id, double time)
{
    avl_node_t *ret = NULL;
    hapgen_tree_node_t *node;
    char *p;

    if (object_heap_empty(&self->avl_node_heap)) {
        goto out;
    }
    p = object_heap_alloc_object(&self->avl_node_heap);
    if (p == NULL) {
        goto out;
    }
    /* We have alloced two objects at once: the avl_node and the
     * hapgen_tree_node that goes with it.
     */
    ret = (avl_node_t *) p;
    node = (hapgen_tree_node_t *) (p + sizeof(avl_node_t));
    node->id = node_id;
    node->children[0] = NULL;
    node->children[1] = NULL;
    node->time = time;
    node->parent = NULL;
    node->branch_length = 0.0;
    avl_init_node(ret, node);
out:
    return ret;
}

static inline void
hapgen_free_avl_node(hapgen_t *self, avl_node_t *avl_node)
{
    object_heap_free_object(&self->avl_node_heap, avl_node);
}

int
hapgen_alloc(hapgen_t *self, tree_sequence_t *tree_sequence,
        double mutation_rate, unsigned long random_seed,
        size_t max_haplotype_length)
{
    int ret = -1;
    uint32_t j;
    uint32_t n = tree_sequence->sample_size;
    avl_node_t *avl_node;

    memset(self, 0, sizeof(hapgen_t));
    self->mutation_rate = mutation_rate;
    self->random_seed = random_seed;
    self->max_haplotype_length = max_haplotype_length;
    self->total_branch_length = 0.0;
    self->num_segregating_sites = 0;
    self->rng = gsl_rng_alloc(gsl_rng_default);
    if (self->rng == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    gsl_rng_set(self->rng, self->random_seed);
    self->sample_size = tree_sequence->sample_size;
    self->num_loci = tree_sequence->num_loci;
    memset(&self->diff_iterator, 0, sizeof(tree_diff_iterator_t));
    ret = tree_diff_iterator_alloc(&self->diff_iterator, tree_sequence, 0);
    if (ret != 0) {
        goto out;
    }
    ret = object_heap_init(&self->avl_node_heap,
            sizeof(avl_node_t) + sizeof(hapgen_tree_node_t), 3 * n, NULL);
    if (ret != 0) {
        goto out;
    }
    avl_init_tree(&self->tree, cmp_hapgen_tree_node, NULL);
    /* Add in the leaf nodes */
    for (j = 1; j <= n; j++) {
        avl_node = hapgen_alloc_avl_node(self, j, 0.0);
        if (avl_node == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        avl_node = avl_insert_node(&self->tree, avl_node);
        assert(avl_node != NULL);
    }
    /* set up the haplotypes */
    self->haplotypes = malloc((n + 1) * sizeof(char *));
    self->haplotype_mem = malloc(n * self->max_haplotype_length);
    if (self->haplotypes == NULL || self->haplotype_mem == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->haplotypes[0] = NULL;
    for (j = 1; j <= n; j++) {
        self->haplotypes[j] = self->haplotype_mem
            + (j - 1) * self->max_haplotype_length;
    }
    /* initialise all haplotypes to '0' */
    memset(self->haplotype_mem, '0', n * self->max_haplotype_length);
    self->traversal_stack = malloc(2 * n * sizeof(hapgen_tree_node_t));
    if (self->traversal_stack == NULL) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

int
hapgen_free(hapgen_t *self)
{
    if (self->rng != NULL) {
        gsl_rng_free(self->rng);
    }
    tree_diff_iterator_free(&self->diff_iterator);
    object_heap_free(&self->avl_node_heap);
    if (self->haplotypes != NULL) {
        free(self->haplotypes);
    }
    if (self->haplotype_mem != NULL) {
        free(self->haplotype_mem);
    }
    if (self->traversal_stack != NULL) {
        free(self->traversal_stack);
    }
    return 0;
}

static int
hapgen_delete_node(hapgen_t *self, uint32_t node_id)
{
    int ret = 0;
    hapgen_tree_node_t search, *node;
    avl_node_t *avl_node;

    search.id = node_id;
    avl_node = avl_search(&self->tree, &search);
    assert(avl_node != NULL);
    node = (hapgen_tree_node_t *) avl_node->item;
    /* If the first child is NULL, we're no longer using this node
     * and can free it.
     */
    if (node->children[0] == NULL) {
        avl_unlink_node(&self->tree, avl_node);
        hapgen_free_avl_node(self, avl_node);
    }
    return ret;
}

/* Update the node to indicate that it has been removed. This is done
 * by setting the first child to NULL. If this has not been reset to a
 * different node after the in records have been applied, we know that
 * this node can be removed from the tree.
 */
static int
hapgen_update_out_node(hapgen_t *self, uint32_t node_id)
{
    int ret = 0;
    hapgen_tree_node_t search, *node;
    avl_node_t *avl_node;
    unsigned int j;

    search.id = node_id;
    avl_node = avl_search(&self->tree, &search);
    assert(avl_node != NULL);
    node = (hapgen_tree_node_t *) avl_node->item;
    for (j = 0; j < 2; j++) {
        self->total_branch_length -= node->children[j]->branch_length;
        node->children[j] = NULL;
    }
    return ret;
}

static int
hapgen_insert_node(hapgen_t *self, uint32_t node_id, uint32_t *children,
        double time)
{
    int ret = 0;
    unsigned int j;
    avl_node_t *avl_node;
    hapgen_tree_node_t search, *node, *child_node;

    search.id = node_id;
    avl_node = avl_search(&self->tree, &search);
    if (avl_node == NULL) {
        avl_node = hapgen_alloc_avl_node(self, node_id, time);
        if (avl_node == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        avl_node = avl_insert_node(&self->tree, avl_node);
        assert(avl_node != NULL);
    }
    node = (hapgen_tree_node_t *) avl_node->item;
    /* Update the parent pointers for the children */
    for (j = 0; j < 2; j++) {
        search.id = children[j];
        avl_node = avl_search(&self->tree, &search);
        assert(avl_node != 0);
        child_node = (hapgen_tree_node_t *) avl_node->item;
        node->children[j] = child_node;
        child_node->parent = node;
    }
out:
    return ret;
}

static int
hapgen_update_branch_lengths(hapgen_t *self, uint32_t node_id)
{
    int ret = 0;
    unsigned int j;
    double bl;
    avl_node_t *avl_node;
    hapgen_tree_node_t search, *node;

    search.id = node_id;
    avl_node = avl_search(&self->tree, &search);
    assert(avl_node != NULL);
    node = (hapgen_tree_node_t *) avl_node->item;
    /* Update the parent pointers for the children */
    for (j = 0; j < 2; j++) {
        bl = node->time - node->children[j]->time;
        node->children[j]->branch_length = bl;
        self->total_branch_length += bl;
    }
    return ret;
}


static int
hapgen_update_root(hapgen_t *self)
{
    avl_node_t *avl_node;
    hapgen_tree_node_t *node = NULL;
    hapgen_tree_node_t search;

    search.id = 1;
    while ((avl_node = avl_search(&self->tree, &search)) != NULL) {
        node = (hapgen_tree_node_t *) avl_node->item;
        search.id = node->parent == NULL? 0 : node->parent->id;
    }
    if (node->parent != NULL) {
        node->parent = NULL;
    }
    node->branch_length = 0.0;
    self->root = node;
    return 0;
}

/* Returns a branch chosen with probability equal to its relative length. */
static hapgen_tree_node_t *
hapgen_select_branch(hapgen_t *self)
{
    avl_node_t *avl_node = self->tree.head;
    hapgen_tree_node_t *node = NULL;
    hapgen_tree_node_t *ret = NULL;
    double lhs = 0.0;
    double rhs = 0.0;
    double u = gsl_rng_uniform(self->rng);

    /* Note that this might do better if we performed a top-down traversal
     * of the tree rather than a random-ish walk around it, since the
     * longest branches (and hence the most likely to be chosen) are near
     * the root.
     */
    assert(avl_node != NULL);
    do {
        node = (hapgen_tree_node_t *) avl_node->item;
        lhs = rhs;
        rhs += node->branch_length / self->total_branch_length;
        if (lhs < u && u < rhs) {
            ret = node;
        }
        avl_node = avl_node->next;
    } while (avl_node != NULL && ret == NULL);
    return ret;
}

/* Adds a mutation on the branch above this node. This means adding
 * a 1 to the end of all haplotypes below this node.
 */
static int
hapgen_add_mutation(hapgen_t *self, hapgen_tree_node_t *node)
{
    int ret = -1;
    hapgen_tree_node_t **stack = self->traversal_stack;
    int stack_top = 0;
    unsigned int j;
    hapgen_tree_node_t *u;

    if (self->num_segregating_sites >= self->max_haplotype_length - 1) {
        ret = MSP_ERR_TOO_MANY_SEG_SITES;
        goto out;
    }
    stack[0] = node;
    while (stack_top >= 0) {
        u = stack[stack_top];
        stack_top--;
        if (u->children[0] == NULL) {
            /* leaf node */
            assert(u->id >= 1 && u->id <= self->sample_size);
            self->haplotypes[u->id][self->num_segregating_sites] = '1';
        } else {
            for (j = 0; j < 2; j++) {
                stack_top++;
                stack[stack_top] = u->children[j];
            }
        }
    }
    self->num_segregating_sites++;
    ret = 0;
out:
    return ret;
}

static int
hapgen_generate_mutations(hapgen_t *self, uint32_t length)
{
    int ret = -1;
    unsigned int j;
    double mu = (self->total_branch_length * self->mutation_rate * length)
        / self->num_loci;
    unsigned int num_mutations = gsl_ran_poisson(self->rng, 1.0 / mu);
    hapgen_tree_node_t *node;

    for (j = 0; j < num_mutations; j++) {
        node = hapgen_select_branch(self);
        if (node == NULL) {
            goto out;
        }
        ret = hapgen_add_mutation(self, node);
        if (ret < 0) {
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

static int
hapgen_process_tree(hapgen_t *self, uint32_t length,
        tree_node_t *nodes_out, tree_node_t *nodes_in)
{
    int ret = 0;
    tree_node_t *tree_node;

    /* mark these nodes as removed. */
    tree_node = nodes_out;
    while (tree_node != NULL) {
        ret = hapgen_update_out_node(self, tree_node->id);
        if (ret != 0) {
            goto out;
        }
        tree_node = tree_node->next;
    }
    /* insert the new records */
    tree_node = nodes_in;
    while (tree_node != NULL) {
        ret = hapgen_insert_node(self, tree_node->id,
                tree_node->children, tree_node->time);
        if (ret != 0) {
            goto out;
        }
        tree_node = tree_node->next;
    }
    /* now, delete any nodes that we need to clear out of the tree */
    tree_node = nodes_out;
    while (tree_node != NULL) {
        ret = hapgen_delete_node(self, tree_node->id);
        if (ret != 0) {
            goto out;
        }
        tree_node = tree_node->next;
    }
    /* update the root */
    ret = hapgen_update_root(self);
    if (ret != 0) {
        goto out;
    }
    /* Update the new branch lengths */
    tree_node = nodes_in;
    while (tree_node != NULL) {
        ret = hapgen_update_branch_lengths(self, tree_node->id);
        if (ret != 0) {
            goto out;
        }
        tree_node = tree_node->next;
    }
    ret = hapgen_generate_mutations(self, length);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int
hapgen_generate(hapgen_t *self)
{
    int ret = -1;
    uint32_t j, length;
    tree_node_t *nodes_out, *nodes_in;

    while ((ret = tree_diff_iterator_next(&self->diff_iterator,
            &length, &nodes_out, &nodes_in)) == 1) {
        ret = hapgen_process_tree(self, length, nodes_out, nodes_in);
        if (ret != 0) {
            goto out;
        }
        assert(avl_count(&self->tree) == 2 * self->sample_size - 1);
    }
    /* NULL terminate the haplotype strings. */
    assert(self->num_segregating_sites < self->max_haplotype_length);
    for (j = 1; j <= self->sample_size; j++) {
        self->haplotypes[j][self->num_segregating_sites] = '\0';
    }
    ret = 0;
out:
    return ret;
}

int
hapgen_get_haplotype(hapgen_t *self, size_t j, char **haplotype)
{
    int ret = 0;
    if (j < 1 || j > self->sample_size) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    *haplotype = self->haplotypes[j];
out:
    return ret;
}
