/*
** Copyright (C) 2016 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

#include "err.h"
#include "object_heap.h"
#include "msprime.h"

static void
ld_calc_check_state(ld_calc_t *self)
{
    uint32_t u;
    uint32_t num_nodes = tree_sequence_get_num_nodes(self->tree_sequence);
    sparse_tree_t *tA = self->outer_tree;
    sparse_tree_t *tB = self->inner_tree;

    assert(tA->index == tB->index);

    /* The inner tree's mark values should all be zero. */
    for (u = 0; u < num_nodes; u++) {
        assert(tA->marked[u] == 0);
        assert(tB->marked[u] == 0);
    }
}

void
ld_calc_print_state(ld_calc_t *self, FILE *out)
{
    fprintf(out, "tree_sequence = %p\n", (void *) self->tree_sequence);
    fprintf(out, "outer tree index = %d\n", (int) self->outer_tree->index);
    fprintf(out, "outer tree interval = (%f, %f)\n",
            self->outer_tree->left, self->outer_tree->right);
    fprintf(out, "inner tree index = %d\n", (int) self->inner_tree->index);
    fprintf(out, "inner tree interval = (%f, %f)\n",
            self->inner_tree->left, self->inner_tree->right);
    ld_calc_check_state(self);
}

int WARN_UNUSED
ld_calc_alloc(ld_calc_t *self, tree_sequence_t *tree_sequence)
{
    int ret = MSP_ERR_GENERIC;

    memset(self, 0, sizeof(ld_calc_t));
    if (pthread_mutex_init(&self->work_mutex, NULL) != 0) {
        ret = MSP_ERR_PTHREAD;
        goto out;
    }
    self->tree_sequence = tree_sequence;
    self->num_mutations = tree_sequence_get_num_mutations(tree_sequence);
    self->outer_tree = malloc(sizeof(sparse_tree_t));
    self->inner_tree = malloc(sizeof(sparse_tree_t));
    if (self->outer_tree == NULL || self->inner_tree == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = sparse_tree_alloc(self->outer_tree, self->tree_sequence,
            MSP_LEAF_COUNTS|MSP_LEAF_LISTS);
    if (ret != 0) {
        goto out;
    }
    ret = sparse_tree_alloc(self->inner_tree, self->tree_sequence,
            MSP_LEAF_COUNTS);
    if (ret != 0) {
        goto out;
    }
    ret = sparse_tree_first(self->outer_tree);
    if (ret < 0) {
        goto out;
    }
    assert(ret != 0);
    ret = sparse_tree_first(self->inner_tree);
    if (ret < 0) {
        goto out;
    }
    assert(ret != 0);
    ret = tree_sequence_get_mutations(self->tree_sequence, &self->mutations);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int
ld_calc_free(ld_calc_t *self)
{
    if (self->inner_tree != NULL) {
        sparse_tree_free(self->inner_tree);
        free(self->inner_tree);
    }
    if (self->outer_tree != NULL) {
        sparse_tree_free(self->outer_tree);
        free(self->outer_tree);
    }
    pthread_mutex_destroy(&self->work_mutex);
    return 0;
}

/* Position the two trees so that the specified mutation is within their
 * interval.
 */
static int WARN_UNUSED
ld_calc_position_trees(ld_calc_t *self, size_t mutation_index)
{
    int ret = MSP_ERR_GENERIC;
    double x = self->mutations[mutation_index].position;
    sparse_tree_t *tA = self->outer_tree;
    sparse_tree_t *tB = self->inner_tree;

    assert(tA->index == tB->index);
    while (x > tA->right) {
        ret = sparse_tree_next(tA);
        if (ret < 0) {
            goto out;
        }
        assert(ret == 1);
        ret = sparse_tree_next(tB);
        if (ret < 0) {
            goto out;
        }
        assert(ret == 1);
    }
    while (x < tA->left) {
        ret = sparse_tree_prev(tA);
        if (ret < 0) {
            goto out;
        }
        assert(ret == 1);
        ret = sparse_tree_prev(tB);
        if (ret < 0) {
            goto out;
        }
        assert(ret == 1);
    }
    ret = 0;
    assert(x >= tA->left && x < tB->right);
    assert(tA->index == tB->index);
out:
    return ret;
}

static uint32_t
ld_calc_overlap_within_tree(ld_calc_t *self, mutation_t mA, mutation_t mB)
{
    sparse_tree_t *t = self->inner_tree;
    uint32_t u, v, nAB;

    u = mA.node;
    v = mB.node;
    if (t->time[u] > t->time[v]) {
        v = mA.node;
        u = mB.node;
    }
    while (u != v && u != MSP_NULL_NODE) {
        u = t->parent[u];
    }
    nAB = 0;
    if (u == v) {
        nAB = GSL_MIN(t->num_leaves[mA.node], t->num_leaves[mB.node]);
    }
    return nAB;
}

static inline int WARN_UNUSED
ld_calc_set_tracked_leaves(ld_calc_t *self, mutation_t mA)
{
    int ret = 0;
    leaf_list_node_t *head, *tail;

    ret = sparse_tree_get_leaf_list(self->outer_tree, mA.node, &head, &tail);
    if (ret != 0) {
        goto out;
    }
    ret = sparse_tree_set_tracked_leaves_from_leaf_list(self->inner_tree,
            head, tail);
out:
    return ret;
}

static int WARN_UNUSED
ld_calc_get_r2_array_forward(ld_calc_t *self, size_t source_index,
        size_t max_mutations, double max_distance, double *r2,
        size_t *num_r2_values)
{
    int ret = MSP_ERR_GENERIC;
    mutation_t mA, mB;
    double fA, fB, fAB, D;
    int tracked_leaves_set = 0;
    sparse_tree_t *tA, *tB;
    double n = tree_sequence_get_sample_size(self->tree_sequence);
    uint32_t nAB;
    size_t j;

    tA = self->outer_tree;
    tB = self->inner_tree;
    mA = self->mutations[source_index];
    assert(tA->parent[mA.node] != MSP_NULL_NODE);
    fA = tA->num_leaves[mA.node] / n;
    assert(fA > 0);
    tB->mark = 1;
    for (j = 0; j < max_mutations; j++) {
        if (source_index + j + 1 >= self->num_mutations) {
            break;
        }
        mB = self->mutations[source_index + j + 1];
        if (mB.position - mA.position > max_distance) {
            break;
        }
        while (mB.position >= tB->right) {
            ret = sparse_tree_next(tB);
            if (ret < 0) {
                goto out;
            }
            assert(ret == 1);
        }
        assert(tB->parent[mB.node] != MSP_NULL_NODE);
        fB = tB->num_leaves[mB.node] / n;
        assert(fB > 0);
        if (mB.position < tA->right) {
            nAB = ld_calc_overlap_within_tree(self, mA, mB);
        } else {
            if (!tracked_leaves_set && tB->marked[mA.node] == 1) {
                tracked_leaves_set = 1;
                ret = ld_calc_set_tracked_leaves(self, mA);
                if (ret != 0) {
                    goto out;
                }
            }
            if (tracked_leaves_set) {
                nAB = tB->num_tracked_leaves[mB.node];
            } else {
                nAB = ld_calc_overlap_within_tree(self, mA, mB);
            }
        }
        fAB = nAB / n;
        D = fAB - fA * fB;
        r2[j] = D * D / (fA * fB * (1 - fA) * (1 - fB));
    }

    /* Now rewind back the inner iterator and unmark all nodes that
     * were set to 1 as we moved forward. */
    tB->mark = 0;
    while (tB->index > tA->index) {
        ret = sparse_tree_prev(tB);
        if (ret < 0) {
            goto out;
        }
        assert(ret == 1);
    }
    *num_r2_values = j;
    ret = 0;
out:
    return ret;
}

static int WARN_UNUSED
ld_calc_get_r2_array_reverse(ld_calc_t *self, size_t source_index,
        size_t max_mutations, double max_distance, double *r2,
        size_t *num_r2_values)
{
    int ret = MSP_ERR_GENERIC;
    mutation_t mA, mB;
    double fA, fB, fAB, D;
    int tracked_leaves_set = 0;
    sparse_tree_t *tA, *tB;
    double n = tree_sequence_get_sample_size(self->tree_sequence);
    uint32_t nAB;
    size_t j;
    ssize_t mutation_index;

    tA = self->outer_tree;
    tB = self->inner_tree;
    mA = self->mutations[source_index];
    assert(tA->parent[mA.node] != MSP_NULL_NODE);
    fA = tA->num_leaves[mA.node] / n;
    assert(fA > 0);
    tB->mark = 1;
    for (j = 0; j < max_mutations; j++) {
        mutation_index = ((ssize_t) source_index) - ((ssize_t) j) - 1;
        if (mutation_index < 0) {
            break;
        }
        mB = self->mutations[mutation_index];
        if (mA.position - mB.position > max_distance) {
            break;
        }
        while (mB.position < tB->left) {
            ret = sparse_tree_prev(tB);
            if (ret < 0) {
                goto out;
            }
            assert(ret == 1);
        }
        assert(tB->parent[mB.node] != MSP_NULL_NODE);
        fB = tB->num_leaves[mB.node] / n;
        assert(fB > 0);
        if (mB.position >= tA->left) {
            nAB = ld_calc_overlap_within_tree(self, mA, mB);
        } else {
            if (!tracked_leaves_set && tB->marked[mA.node] == 1) {
                tracked_leaves_set = 1;
                ret = ld_calc_set_tracked_leaves(self, mA);
                if (ret != 0) {
                    goto out;
                }
            }
            if (tracked_leaves_set) {
                nAB = tB->num_tracked_leaves[mB.node];
            } else {
                nAB = ld_calc_overlap_within_tree(self, mA, mB);
            }
        }
        fAB = nAB / n;
        D = fAB - fA * fB;
        r2[j] = D * D / (fA * fB * (1 - fA) * (1 - fB));
    }

    /* Now fast forward the inner iterator and unmark all nodes that
     * were set to 1 as we moved back. */
    tB->mark = 0;
    while (tB->index < tA->index) {
        ret = sparse_tree_next(tB);
        if (ret < 0) {
            goto out;
        }
        assert(ret == 1);
    }
    *num_r2_values = j;
    ret = 0;
out:
    return ret;
}

int WARN_UNUSED
ld_calc_get_r2_array(ld_calc_t *self, size_t a, int direction,
        size_t max_mutations, double max_distance, double *r2,
        size_t *num_r2_values)
{
    int ret = MSP_ERR_GENERIC;
    int lock_acquired = 0;

    if (pthread_mutex_lock(&self->work_mutex) != 0) {
        ret = MSP_ERR_PTHREAD;
        goto out;
    }
    lock_acquired = 1;

    if (a >= self->num_mutations) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    ret = ld_calc_position_trees(self, a);
    if (ret != 0) {
        goto out;
    }
    if (direction == MSP_DIR_FORWARD) {
        ret = ld_calc_get_r2_array_forward(self, a, max_mutations, max_distance,
                r2, num_r2_values);
    } else if (direction == MSP_DIR_REVERSE) {
        ret = ld_calc_get_r2_array_reverse(self, a, max_mutations, max_distance,
                r2, num_r2_values);
    } else {
        ret = MSP_ERR_BAD_PARAM_VALUE;
    }
out:
    if (lock_acquired) {
        if (pthread_mutex_unlock(&self->work_mutex) != 0) {
            if (ret != 0) {
                ret = MSP_ERR_PTHREAD;
            }
        }
    }
    return ret;
}

int WARN_UNUSED
ld_calc_get_r2(ld_calc_t *self, size_t a, size_t b, double *r2)
{
    int ret = MSP_ERR_GENERIC;
    mutation_t mA, mB;
    double fA, fB, fAB, D;
    sparse_tree_t *tA, *tB;
    double n = tree_sequence_get_sample_size(self->tree_sequence);
    uint32_t nAB;
    size_t tmp;
    int lock_acquired = 0;

    if (pthread_mutex_lock(&self->work_mutex) != 0) {
        ret = MSP_ERR_PTHREAD;
        goto out;
    }
    lock_acquired = 1;

    if (a >= self->num_mutations || b >= self->num_mutations) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    if (a > b) {
        tmp = a;
        a = b;
        b = tmp;
    }
    ret = ld_calc_position_trees(self, a);
    if (ret != 0) {
        goto out;
    }
    /* We can probably do a lot better than this implementation... */
    tA = self->outer_tree;
    tB = self->inner_tree;
    mA = self->mutations[a];
    mB = self->mutations[b];
    assert(tA->parent[mA.node] != MSP_NULL_NODE);
    fA = tA->num_leaves[mA.node] / n;
    assert(fA > 0);
    ret = ld_calc_set_tracked_leaves(self, mA);
    if (ret != 0) {
        goto out;
    }

    while (mB.position >= tB->right) {
        ret = sparse_tree_next(tB);
        if (ret < 0) {
            goto out;
        }
        assert(ret == 1);
    }
    assert(tB->parent[mB.node] != MSP_NULL_NODE);
    fB = tB->num_leaves[mB.node] / n;
    assert(fB > 0);
    nAB = tB->num_tracked_leaves[mB.node];
    fAB = nAB / n;
    D = fAB - fA * fB;
    *r2 = D * D / (fA * fB * (1 - fA) * (1 - fB));

    /* Now rewind the inner iterator back. */
    while (tB->index > tA->index) {
        ret = sparse_tree_prev(tB);
        if (ret < 0) {
            goto out;
        }
        assert(ret == 1);
    }
    ret = 0;
out:
    if (lock_acquired) {
        if (pthread_mutex_unlock(&self->work_mutex) != 0) {
            if (ret != 0) {
                ret = MSP_ERR_PTHREAD;
            }
        }
    }
    return ret;
}
