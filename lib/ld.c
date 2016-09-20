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

#define MAX_LABEL_LEN 16

static int WARN_UNUSED
ld_calc_make_position_labels(ld_calc_t *self)
{
    int ret = 0;
    size_t j, size, offset;
    size_t buffer_size = MAX_LABEL_LEN * self->num_mutations;
    const char *pattern = "%d\t";
    char *label;
    int pos, written;

    offset = 0;
    /* TODO put in the proper rounding algorithm so we don't have collisions */
    for (j = 0; j < self->num_mutations; j++) {
        pos = (int) round(self->mutations[j].position);
        size = (size_t) snprintf(NULL, 0, pattern, pos);
        if (offset + size >= buffer_size) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        label = self->label_mem + offset;
        written = snprintf(label, size + 1, pattern, pos);
        if (written < 0) {
            ret = MSP_ERR_IO;
            goto out;
        }
        assert(written == (int) size);
        self->position_labels[j] = label;
        offset += size + 1;
    }
out:
    return ret;
}

int WARN_UNUSED
ld_calc_alloc(ld_calc_t *self, tree_sequence_t *tree_sequence, size_t max_sites,
        double max_distance, double r2_threshold)
{

    int ret = MSP_ERR_GENERIC;

    memset(self, 0, sizeof(ld_calc_t));
    self->tree_sequence = tree_sequence;
    self->max_sites = max_sites;
    self->max_distance = max_distance;
    self->r2_threshold = r2_threshold;
    self->num_mutations = tree_sequence_get_num_mutations(tree_sequence);
    self->outer_tree = malloc(sizeof(sparse_tree_t));
    self->inner_tree = malloc(sizeof(sparse_tree_t));
    self->position_labels = malloc(self->num_mutations * sizeof(char *));
    self->label_mem = malloc(self->num_mutations * MAX_LABEL_LEN);
    if (self->outer_tree == NULL || self->inner_tree == NULL
            || self->position_labels == NULL || self->label_mem == NULL) {
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
    ret = tree_sequence_get_mutations(self->tree_sequence, &self->mutations);
    if (ret != 0) {
        goto out;
    }
    ret = ld_calc_make_position_labels(self);
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
    if (self->position_labels != NULL) {
        free(self->position_labels);
    }
    if (self->label_mem != NULL) {
        free(self->label_mem);
    }
    return 0;
}

void
ld_calc_print_state(ld_calc_t *self, FILE *out)
{
    fprintf(out, "tree_sequence = %p\n", (void *) self->tree_sequence);
    fprintf(out, "max_sites = %d\n", (int) self->max_sites);
    fprintf(out, "max_distance = %f\n", self->max_distance);
    fprintf(out, "r2_threshold = %f\n", self->r2_threshold);
    fprintf(out, "outer tree index = %d\n", (int) self->outer_tree->index);
    fprintf(out, "inner tree index = %d\n", (int) self->inner_tree->index);
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

int WARN_UNUSED
ld_calc_write_table(ld_calc_t *self, FILE *out)
{
    int ret = MSP_ERR_GENERIC;
    int err, tracked_leaves_set;
    size_t j, k;
    uint32_t nAB;
    mutation_t mA, mB;
    double r2, fA, fB, fAB, D;
    sparse_tree_t *tA, *tB;
    double n = tree_sequence_get_sample_size(self->tree_sequence);

    /* set the trees to the first tree */
    tA = self->outer_tree;
    tB = self->inner_tree;
    ret = sparse_tree_first(tA);
    if (ret < 0) {
        goto out;
    }
    ret = sparse_tree_next(tB);
    if (ret < 0) {
        goto out;
    }

    for (j = 0; j < self->num_mutations; j++) {
        assert(tA->index == tB->index);
        mA = self->mutations[j];
        while (mA.position >= tA->right) {
            ret = sparse_tree_next(tA);
            if (ret < 0) {
                goto out;
            }
            ret = sparse_tree_next(tB);
            if (ret < 0) {
                goto out;
            }
        }
        assert(tA->parent[mA.node] != MSP_NULL_NODE);
        fA = tA->num_leaves[mA.node] / n;
        assert(fA > 0);
        tB->mark = 1;
        tracked_leaves_set = 0;
        for (k = j + 1; k < self->num_mutations; k++) {
            mB = self->mutations[k];
            if (k - j + 1 > self->max_sites ||
                    mB.position - mA.position > self->max_distance) {
                break;
            }
            while (mB.position >= tB->right) {
                ret = sparse_tree_next(tB);
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
            r2 = D * D / (fA * fB * (1 - fA) * (1 - fB));
            if (r2 >= self->r2_threshold) {
                err = fputs(self->position_labels[j], out);
                if (err < 0) {
                    ret = MSP_ERR_IO;
                    goto out;
                }
                err = fputs(self->position_labels[k], out);
                if (err < 0) {
                    ret = MSP_ERR_IO;
                    goto out;
                }
                err = fprintf(out, "%G\n", r2);
                if (err < 0) {
                    ret = MSP_ERR_IO;
                    goto out;
                }
            }
        }
        tB->mark = 0;
        /* Now rewind back the inner iterator */
        while (tB->index > tA->index) {
            ret = sparse_tree_prev(tB);
            if (ret < 0) {
                goto out;
            }
        }
    }
    ret = 0;
out:
    return ret;
}
