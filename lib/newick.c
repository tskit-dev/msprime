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

#include "err.h"
#include "object_heap.h"
#include "msprime.h"


typedef struct {
    uint32_t id;
    uint32_t parent;
    uint32_t children[2];
    double time;
    double branch_length;
} newick_tree_node_t;

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
    for (j = 1; j <= self->sample_size; j++) {
        search.id = j;
        while ((avl_node = avl_search(&self->tree, &search)) != NULL) {
            node = (newick_tree_node_t *) avl_node->item;
            search.id = node->parent;
        }
        assert(node != NULL);
        assert(node->id = self->root);
    }
}

void
newick_converter_print_state(newick_converter_t *self)
{
    avl_node_t *avl_node;
    newick_tree_node_t *node;

    printf("Newick converter state\n");
    for (avl_node = self->tree.head; avl_node != NULL; avl_node = avl_node->next) {
        node = (newick_tree_node_t *) avl_node->item;
        printf("%d\t%d\t%d\t%d\t%f\n", node->id, node->children[0], node->children[1],
                node->parent, node->time);
    }

    printf("avl_node_heap\n");
    object_heap_print_state(&self->avl_node_heap);
    newick_converter_check_state(self);
}
static inline avl_node_t * WARN_UNUSED
newick_converter_alloc_avl_node(newick_converter_t *self, uint32_t node_id,
        uint32_t *children, double time)
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
    ret = (avl_node_t *) p;
    node = (newick_tree_node_t *) (p + sizeof(avl_node_t));

    node->id = node_id;
    node->children[0] = children[0];
    node->children[1] = children[1];
    node->time = time;
    node->parent = 0;
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
    if (node->children[0] == 0) {
        avl_unlink_node(&self->tree, avl_node);
        newick_converter_free_avl_node(self, avl_node);
    }
    return ret;
}

/* Update the node to indicate that it has been removed. This is done
 * by setting the first child to 0. If this has not been reset to a
 * different node after the in records have been applied, we know that
 * this node can be removed from the tree.
 */
static int
newick_converter_update_out_node(newick_converter_t *self, uint32_t node_id)
{
    int ret = 0;
    newick_tree_node_t search, *node;
    avl_node_t *avl_node;

    search.id = node_id;
    avl_node = avl_search(&self->tree, &search);
    assert(avl_node != NULL);
    node = (newick_tree_node_t *) avl_node->item;
    node->children[0] = 0;
    return ret;
}


static int
newick_converter_insert_node(newick_converter_t *self, uint32_t node_id,
        uint32_t *children, double time)
{
    int ret = 0;
    unsigned int j;
    avl_node_t *avl_node;
    newick_tree_node_t search, *node;

    search.id = node_id;
    avl_node = avl_search(&self->tree, &search);
    if (avl_node == NULL) {
        avl_node = newick_converter_alloc_avl_node(self, node_id, children, time);
        if (avl_node == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        avl_node = avl_insert_node(&self->tree, avl_node);
        assert(avl_node != NULL);
    } else {
        /* node already exists so we update it */
        node = (newick_tree_node_t *) avl_node->item;
        node->children[0] = children[0];
        node->children[1] = children[1];
    }
    /* Update the parent pointers for the children */
    for (j = 0; j < 2; j++) {
        search.id = children[j];
        avl_node = avl_search(&self->tree, &search);
        assert(avl_node != NULL);
        node = (newick_tree_node_t *) avl_node->item;
        node->parent = node_id;
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

    search.id = 1;
    while ((avl_node = avl_search(&self->tree, &search)) != NULL) {
        node = (newick_tree_node_t *) avl_node->item;
        search.id = node->parent;
    }
    assert(node != NULL);
    node->parent = 0;
    self->root = node->id;
    return 0;
}

static int
newick_converter_process_tree(newick_converter_t *self, tree_node_t *nodes_out,
        tree_node_t *nodes_in)
{
    int ret = 0;
    tree_node_t *tree_node;

    /* mark these nodes as removed. */
    tree_node = nodes_out;
    while (tree_node != NULL) {
        ret = newick_converter_update_out_node(self, tree_node->id);
        if (ret != 0) {
            goto out;
        }
        tree_node = tree_node->next;
    }
    /* insert the new records */
    tree_node = nodes_in;
    while (tree_node != NULL) {
        ret = newick_converter_insert_node(self, tree_node->id,
                tree_node->children, tree_node->time);
        if (ret != 0) {
            goto out;
        }
        tree_node = tree_node->next;
    }
    /* now, delete any nodes that we need to clear out of the tree */
    tree_node = nodes_out;
    while (tree_node != NULL) {
        ret = newick_converter_delete_node(self, tree_node->id);
        if (ret != 0) {
            goto out;
        }
        tree_node = tree_node->next;
    }
    /* update the root */
    ret = newick_converter_update_root(self);
out:
    return ret;
}

int
newick_converter_next(newick_converter_t *self, uint32_t *length, char **tree)
{
    int ret = -1;
    int err;
    tree_node_t *nodes_out, *nodes_in;

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
        newick_converter_print_state(self);
    }
out:
    return ret;
}

int
newick_converter_alloc(newick_converter_t *self,
        tree_sequence_t *tree_sequence, size_t precision, int all_breakpoints)
{
    int ret = -1;
    int flags = 0;
    uint32_t j;
    uint32_t empty_children[] = {0, 0};
    uint32_t n = tree_sequence->sample_size;
    avl_node_t *avl_node;

    self->sample_size = tree_sequence->sample_size;
    self->num_loci = tree_sequence->num_loci;
    memset(&self->diff_iterator, 0, sizeof(tree_diff_iterator_t));
    if (all_breakpoints) {
        flags = MSP_ALL_BREAKPOINTS;
    }
    ret = tree_diff_iterator_alloc(&self->diff_iterator, tree_sequence, flags);
    if (ret != 0) {
        goto out;
    }
    ret = object_heap_init(&self->avl_node_heap,
            sizeof(avl_node_t) + sizeof(newick_tree_node_t), 3 * n, NULL);
    if (ret != 0) {
        goto out;
    }
    avl_init_tree(&self->tree, cmp_newick_tree_node, NULL);
    /* Add in the leaf nodes */
    for (j = 1; j <= n; j++) {
        avl_node = newick_converter_alloc_avl_node(self, j, empty_children,
                0.0);
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
    tree_diff_iterator_free(&self->diff_iterator);
    object_heap_free(&self->avl_node_heap);
    return 0;
}
