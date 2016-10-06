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
vargen_print_state(vargen_t *self, FILE *out)
{
    fprintf(out, "vargen state\n");
    fprintf(out, "tree_mutation_index = %d\n", (int) self->tree_mutation_index);
}

static int
vargen_next_tree(vargen_t *self)
{
    int ret = 0;

    ret = sparse_tree_next(&self->tree);
    if (ret == 0) {
        self->finished = 1;
    } else if (ret < 0) {
        goto out;
    }
    self->tree_mutation_index = 0;
out:
    return ret;
}

int
vargen_alloc(vargen_t *self, tree_sequence_t *tree_sequence, int flags)
{
    int ret = MSP_ERR_NO_MEMORY;

    assert(tree_sequence != NULL);
    memset(self, 0, sizeof(vargen_t));
    self->sample_size = tree_sequence_get_sample_size(tree_sequence);
    self->sequence_length = tree_sequence_get_sequence_length(tree_sequence);
    self->num_mutations = tree_sequence_get_num_mutations(tree_sequence);
    self->tree_sequence = tree_sequence;
    self->flags = flags;

    ret = sparse_tree_alloc(&self->tree, tree_sequence, MSP_LEAF_LISTS);
    if (ret != 0) {
        goto out;
    }
    self->finished = 0;
    self->tree_mutation_index = 0;
    ret = sparse_tree_first(&self->tree);
    if (ret < 0) {
        goto out;
    }
    assert(ret == 1);
    ret = 0;
out:
    return ret;
}

int
vargen_free(vargen_t *self)
{
    sparse_tree_free(&self->tree);
    return 0;
}

static int
vargen_apply_tree_mutation(vargen_t *self, mutation_t *mut, char *genotypes)
{
    int ret = 0;
    leaf_list_node_t *w, *tail;
    uint32_t parent;
    int not_done = 1;
    char one = self->flags & MSP_GENOTYPES_AS_CHAR? '1': 1;

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
        genotypes[w->node] = one;
        not_done = w != tail;
        w = w->next;
    }
out:
    return ret;
}

int
vargen_next(vargen_t *self, mutation_t **mutation, char *genotypes)
{
    int ret = 0;
    int not_done = 1;
    mutation_t *m;
    char zero = self->flags & MSP_GENOTYPES_AS_CHAR? '0': 0;

    if (!self->finished) {
        while (not_done && self->tree_mutation_index == self->tree.num_mutations) {
            ret = vargen_next_tree(self);
            if (ret < 0) {
                goto out;
            }
            not_done = ret == 1;
        }
        if (not_done) {
            memset(genotypes, zero, self->sample_size);
            m = &self->tree.mutations[self->tree_mutation_index];
            ret = vargen_apply_tree_mutation(self, m, genotypes);
            if (ret != 0) {
                goto out;
            }
            self->tree_mutation_index++;
            *mutation = m;
            ret = 1;
        }
    }
out:
    return ret;
}
