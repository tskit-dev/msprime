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

#include "err.h"
#include "object_heap.h"
#include "msprime.h"

void
vargen_print_state(vargen_t *self)
{
    printf("vargen state\n");
    printf("tree_mutation_index = %d\n", (int) self->tree_mutation_index);
    printf("variant = '%s'\n", self->variant);
    sparse_tree_iterator_print_state(&self->tree_iterator);
}

static int
vargen_next_tree(vargen_t *self)
{
    int ret = 0;

    ret = sparse_tree_iterator_next(&self->tree_iterator);
    if (ret < 0) {
        goto out;
    }
    self->tree_mutation_index = 0;
out:
    return ret;
}

int
vargen_alloc(vargen_t *self, tree_sequence_t *tree_sequence)
{
    int ret = MSP_ERR_NO_MEMORY;

    assert(tree_sequence != NULL);
    memset(self, 0, sizeof(vargen_t));
    self->sample_size = tree_sequence_get_sample_size(tree_sequence);
    self->sequence_length = tree_sequence_get_sequence_length(tree_sequence);
    self->num_mutations = tree_sequence_get_num_mutations(tree_sequence);
    self->tree_sequence = tree_sequence;
    self->variant = calloc(self->sample_size + 1, sizeof(char));

    if (self->variant == NULL) {
        goto out;
    }
    ret = tree_sequence_alloc_sparse_tree(tree_sequence, &self->tree,
            NULL, 0, MSP_COUNT_LEAVES);
    if (ret != 0) {
        goto out;
    }
    ret = sparse_tree_iterator_alloc(&self->tree_iterator,
            self->tree_sequence, &self->tree);
    if (ret != 0) {
        goto out;
    }
    ret = vargen_next_tree(self);
    /* We must have at least one tree in the iterator */
    assert(ret == 1);
    ret = 0;
out:
    return ret;
}

int
vargen_free(vargen_t *self)
{
    if (self->variant != NULL) {
        free(self->variant);
    }
    sparse_tree_free(&self->tree);
    sparse_tree_iterator_free(&self->tree_iterator);
    return 0;
}

static int
vargen_apply_tree_mutation(vargen_t *self, mutation_t *mut)
{
    int ret = 0;
    leaf_list_node_t *w, *tail;
    uint32_t parent;
    int not_done = 1;

    ret = sparse_tree_get_parent(&self->tree, mut->node, &parent);
    if (ret != 0) {
        goto out;
    }
    if (parent == MSP_NULL_NODE) {
        ret = MSP_ERR_BAD_MUTATION;
        goto out;
    }
    ret = sparse_tree_get_leaf_list(&self->tree, mut->node, &w, &tail);
    if (ret != 0) {
        goto out;
    }
    while (not_done) {
        assert(w != NULL);
        assert(w->node < self->sample_size);
        self->variant[w->node] = '1';
        not_done = w != tail;
        w = w->next;
    }
out:
    return ret;
}

int
vargen_next(vargen_t *self, double *position, char **variant)
{
    int ret = 0;
    int not_done = 1;
    mutation_t *mutation;

    while (not_done && self->tree_mutation_index == self->tree.num_mutations) {
        ret = vargen_next_tree(self);
        if (ret < 0) {
            goto out;
        }
        not_done = ret == 1;
    }
    if (not_done) {
        memset(self->variant, '0', self->sample_size);
        mutation = &self->tree.mutations[self->tree_mutation_index];
        ret = vargen_apply_tree_mutation(self, mutation);
        if (ret != 0) {
            goto out;
        }
        self->tree_mutation_index++;
        *variant = self->variant;
        *position = mutation->position;
        ret = 1;
    }
out:
    return ret;
}
