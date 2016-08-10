/*
** Copyright (C) 2015 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

#include "err.h"
#include "object_heap.h"
#include "msprime.h"

static int
cmp_newick_tree_node(const void *a, const void *b) {
    const newick_tree_node_t *ia = (const newick_tree_node_t *) a;
    const newick_tree_node_t *ib = (const newick_tree_node_t *) b;
    return (ia->id > ib->id) - (ia->id < ib->id);
}

/* Ensure the tree is in a consistent state */
static void
newick_converter_check_state(newick_converter_t *self)
{
    uint32_t j;
    avl_node_t *avl_node;
    newick_tree_node_t *node = NULL;
    newick_tree_node_t search;

    assert(avl_count(&self->tree) == 2 * self->sample_size - 1);
    for (j = 0; j < self->sample_size; j++) {
        search.id = j;
        avl_node = avl_search(&self->tree, &search);
        assert(avl_node != NULL);
        node = (newick_tree_node_t *) avl_node->item;
        while (node->parent != NULL) {
            node = node->parent;
        }
        assert(node != NULL);
        assert(node == self->root);
    }
}

void
newick_converter_print_state(newick_converter_t *self, FILE *out)
{
    avl_node_t *avl_node;
    newick_tree_node_t *node;
    size_t j;

    fprintf(out, "Newick converter state\n");
    fprintf(out, "num_nodes = %d\n", avl_count(&self->tree));
    fprintf(out, "root = %d\n", self->root == NULL? 0: self->root->id);
    for (avl_node = self->tree.head; avl_node != NULL; avl_node = avl_node->next) {
        node = (newick_tree_node_t *) avl_node->item;
        fprintf(out, "%d\t", node->id);
        for (j = 0; j < 2; j++) {
            fprintf(out, "%d\t", node->children[j] == NULL? 0 :
                    node->children[j]->id);
        }
        fprintf(out, "%d\t%f\t%s\t",
                node->parent == NULL ? 0: node->parent->id, node->time,
                node->branch_length);
        if (node->subtree == NULL) {
            fprintf(out, "NULL\n");
        } else {
            fprintf(out, "%s\n", node->subtree);
        }
    }
    fprintf(out, "avl_node_heap\n");
    object_heap_print_state(&self->avl_node_heap, out);
    newick_converter_check_state(self);
}

static inline avl_node_t * WARN_UNUSED
newick_converter_alloc_avl_node(newick_converter_t *self, uint32_t node_id,
        double time)
{
    avl_node_t *ret = NULL;
    newick_tree_node_t *node;
    char *p;

    if (object_heap_empty(&self->avl_node_heap)) {
        goto out;
    }
    p = object_heap_alloc_object(&self->avl_node_heap);
    if (p == NULL) {
        goto out;
    }
    /* We have alloced two objects at once: the avl_node and the
     * newick_tree_node that goes with it.
     */
    /* cast to void * here to avoid alignment warnings from clang. */
    ret = (void *) p;
    node = (void *) (p + sizeof(avl_node_t));

    node->id = node_id;
    node->children[0] = NULL;
    node->children[1] = NULL;
    node->time = time;
    node->parent = NULL;
    node->branch_length[0] = '\0';
    node->subtree = NULL;
    avl_init_node(ret, node);
out:
    return ret;
}

static inline void
newick_converter_free_avl_node(newick_converter_t *self, avl_node_t *avl_node)
{
    object_heap_free_object(&self->avl_node_heap, avl_node);
}

static int
newick_converter_delete_node(newick_converter_t *self, uint32_t node_id)
{
    int ret = 0;
    newick_tree_node_t search, *node;
    avl_node_t *avl_node;

    search.id = node_id;
    avl_node = avl_search(&self->tree, &search);
    assert(avl_node != NULL);
    node = (newick_tree_node_t *) avl_node->item;
    /* If the first child is NULL, we're no longer using this node
     * and can free it.
     */
    if (node->children[0] == NULL) {
        avl_unlink_node(&self->tree, avl_node);
        newick_converter_free_avl_node(self, avl_node);
    }
    return ret;
}

/* Update the node to indicate that it has been removed. This is done
 * by setting the first child to NULL. If this has not been reset to a
 * different node after the in records have been applied, we know that
 * this node can be removed from the tree.
 */
static int
newick_converter_update_out_node(newick_converter_t *self, uint32_t node_id)
{
    int ret = 0;
    newick_tree_node_t search, *node, *u;;
    avl_node_t *avl_node;
    unsigned int j;

    search.id = node_id;
    avl_node = avl_search(&self->tree, &search);
    assert(avl_node != NULL);
    node = (newick_tree_node_t *) avl_node->item;
    for (j = 0; j < 2; j++) {
        u = node->children[j];
        node->children[j] = NULL;
        /* Free the subtree and propagate this up the tree */
        while (u != NULL && u->subtree != NULL) {
            free(u->subtree);
            u->subtree = NULL;
            u = u->parent;
        }
    }
    return ret;
}

static int
newick_tree_set_branch_length(newick_converter_t *self,
        newick_tree_node_t *node, double parent_time)
{
    int ret = 0;
    double length = parent_time - node->time;
    int r;

    /* We rescale branch lengths to be in coalescent time units. */
    length /= 4 * self->Ne;
    r = snprintf(node->branch_length, MAX_BRANCH_LENGTH_STRING,
            "%.*f", (int) self->precision, length);
    if (r >= MAX_BRANCH_LENGTH_STRING) {
        ret = MSP_ERR_NEWICK_OVERFLOW;
        goto out;
    }
out:
    return ret;
}

static int
newick_converter_insert_node(newick_converter_t *self, uint32_t node_id,
        uint32_t *children, double time)
{
    int ret = 0;
    unsigned int j;
    avl_node_t *avl_node;
    newick_tree_node_t search, *node, *child_node;

    search.id = node_id;
    avl_node = avl_search(&self->tree, &search);
    if (avl_node == NULL) {
        avl_node = newick_converter_alloc_avl_node(self, node_id, time);
        if (avl_node == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        avl_node = avl_insert_node(&self->tree, avl_node);
        assert(avl_node != NULL);
    }
    node = (newick_tree_node_t *) avl_node->item;
    /* Update the parent pointers for the children */
    for (j = 0; j < 2; j++) {
        search.id = children[j];
        avl_node = avl_search(&self->tree, &search);
        assert(avl_node != 0);
        child_node = (newick_tree_node_t *) avl_node->item;
        node->children[j] = child_node;
        child_node->parent = node;
        ret = newick_tree_set_branch_length(self, child_node, time);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

static int
newick_converter_update_root(newick_converter_t *self)
{
    avl_node_t *avl_node;
    newick_tree_node_t *node = NULL;
    newick_tree_node_t search;

    search.id = 0;
    while ((avl_node = avl_search(&self->tree, &search)) != NULL) {
        node = (newick_tree_node_t *) avl_node->item;
        search.id = node->parent == NULL? MSP_NULL_NODE : node->parent->id;
    }
    if (node->parent != NULL) {
        node->parent = NULL;
        node->branch_length[0] = '\0';
    }
    self->root = node;
    return 0;
}

static int
newick_converter_generate_subtree(newick_converter_t *self,
        newick_tree_node_t *node)
{
    int ret = 0;
    size_t size, s1_len, s2_len;
    const char *leaf_format = "%d:%s";
    char sep, *s, *s1, *s2;
    int label;

    if (node->children[0] == NULL) {
        /* leaf node */
        /* TODO For ms compatablility we set the ID to 1 here. We should make
         * this a configurable behaviour.
         */
        label = (int) node->id + 1;
        size = (size_t) snprintf(NULL, 0, leaf_format, label,
                node->branch_length);
        node->subtree = malloc(size + 1);
        if (node->subtree == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        sprintf(node->subtree, leaf_format, label, node->branch_length);
    } else {
        s1 = node->children[0]->subtree;
        assert(s1 != NULL);
        s1_len = strlen(s1);
        s2 = node->children[1]->subtree;
        assert(s2 != NULL);
        s2_len = strlen(s2);
        size = s1_len + s2_len + 5;
        sep = ';';
        if (node != self->root) {
            size += strlen(node->branch_length);
            sep = ':';
        }
        node->subtree = malloc(size);
        if (node->subtree == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        s = node->subtree;
        s[size - 1] = '\0';
        *s = '(';
        s++;
        strcpy(s, s1);
        s += s1_len;
        *s = ',';
        s++;
        strcpy(s, s2);
        s += s2_len;
        *s = ')';
        s++;
        *s = sep;
        if (node != self->root) {
            s++;
            strcpy(s, node->branch_length);
        }
    }
out:
    return ret;
}

static int
newick_converter_update_subtrees(newick_converter_t *self)
{
    int ret = 0;
    newick_tree_node_t **traversal_stack = NULL;
    newick_tree_node_t **visit_list = NULL;
    int traversal_stack_top = 0;
    int visit_list_top = -1;
    size_t j;
    newick_tree_node_t *node;

    /* We use a two-stack iterative method to perform a post-order
     * traversal of the nodes which have NULL subtrees.
     */
    traversal_stack = malloc(self->sample_size * sizeof(newick_tree_node_t *));
    visit_list = malloc(2 * self->sample_size * sizeof(newick_tree_node_t *));
    if (traversal_stack == NULL || visit_list == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    traversal_stack[0] = self->root;
    while (traversal_stack_top >= 0) {
        node = traversal_stack[traversal_stack_top];
        traversal_stack_top--;
        if (node->subtree == NULL) {
            visit_list_top++;
            visit_list[visit_list_top] = node;
            for (j = 0; j < 2; j++) {
                if (node->children[j] != NULL) {
                    traversal_stack_top++;
                    traversal_stack[traversal_stack_top] = node->children[j];
                }
            }
        }
    }
    while (visit_list_top >= 0) {
        node = visit_list[visit_list_top];
        ret = newick_converter_generate_subtree(self, node);
        if (ret != 0) {
            goto out;
        }
        visit_list_top--;
    }
out:
    if (traversal_stack != NULL) {
        free(traversal_stack);
    }
    if (visit_list != NULL) {
        free(visit_list);
    }
    return ret;
}

static int
newick_converter_process_tree(newick_converter_t *self, node_record_t *nodes_out,
        node_record_t *nodes_in)
{
    int ret = 0;
    node_record_t *tree_node;

    /* mark these nodes as removed. */
    tree_node = nodes_out;
    while (tree_node != NULL) {
        if (tree_node->num_children > 2) {
            ret = MSP_ERR_NONBINARY_NEWICK;
            goto out;
        }
        ret = newick_converter_update_out_node(self, tree_node->node);
        if (ret != 0) {
            goto out;
        }
        tree_node = tree_node->next;
    }
    /* insert the new records */
    tree_node = nodes_in;
    while (tree_node != NULL) {
        if (tree_node->num_children > 2) {
            ret = MSP_ERR_NONBINARY_NEWICK;
            goto out;
        }
        ret = newick_converter_insert_node(self, tree_node->node,
                tree_node->children, tree_node->time);
        if (ret != 0) {
            goto out;
        }
        tree_node = tree_node->next;
    }
    /* now, delete any nodes that we need to clear out of the tree */
    tree_node = nodes_out;
    while (tree_node != NULL) {
        ret = newick_converter_delete_node(self, tree_node->node);
        if (ret != 0) {
            goto out;
        }
        tree_node = tree_node->next;
    }
    /* update the root */
    ret = newick_converter_update_root(self);
    if (ret != 0) {
        goto out;
    }
    ret = newick_converter_update_subtrees(self);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int
newick_converter_next(newick_converter_t *self, double *length, char **tree)
{
    int ret = -1;
    int err;
    node_record_t *nodes_out, *nodes_in;

    ret = tree_diff_iterator_next(&self->diff_iterator, length, &nodes_out,
            &nodes_in);
    if (ret < 0) {
        goto out;
    }
    if (ret == 1) {
        err = newick_converter_process_tree(self, nodes_out, nodes_in);
        if (err != 0) {
            ret = err;
            goto out;
        }
        assert(self->root->subtree != NULL);
        assert(avl_count(&self->tree) == 2 * self->sample_size - 1);
        *tree = self->root->subtree;
    }
out:
    return ret;
}

int
newick_converter_alloc(newick_converter_t *self,
        tree_sequence_t *tree_sequence, size_t precision, double Ne)
{
    int ret = -1;
    uint32_t j;
    avl_node_t *avl_node;

    memset(self, 0, sizeof(newick_converter_t));
    self->sample_size = tree_sequence_get_sample_size(tree_sequence);
    self->sequence_length = tree_sequence_get_sequence_length(tree_sequence);
    self->precision = precision;
    self->Ne = Ne;
    memset(&self->diff_iterator, 0, sizeof(tree_diff_iterator_t));
    ret = tree_diff_iterator_alloc(&self->diff_iterator, tree_sequence);
    if (ret != 0) {
        goto out;
    }
    ret = object_heap_init(&self->avl_node_heap,
            sizeof(avl_node_t) + sizeof(newick_tree_node_t),
            3 * self->sample_size, NULL);
    if (ret != 0) {
        goto out;
    }
    avl_init_tree(&self->tree, cmp_newick_tree_node, NULL);
    /* Add in the leaf nodes */
    for (j = 0; j < self->sample_size; j++) {
        avl_node = newick_converter_alloc_avl_node(self, j, 0.0);
        if (avl_node == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        avl_node = avl_insert_node(&self->tree, avl_node);
        assert(avl_node != NULL);
    }
out:
    return ret;
}

int
newick_converter_free(newick_converter_t *self)
{
    avl_node_t *avl_node;
    newick_tree_node_t *node;

    /* Clear out any dangling subtree strings. */
    for (avl_node = self->tree.head; avl_node != NULL;
            avl_node = avl_node->next) {
        node = (newick_tree_node_t *) avl_node->item;
        if (node->subtree != NULL) {
            free(node->subtree);
        }
    }
    tree_diff_iterator_free(&self->diff_iterator);
    object_heap_free(&self->avl_node_heap);
    return 0;
}
