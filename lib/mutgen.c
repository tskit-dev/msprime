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

#include <gsl/gsl_randist.h>

#include "err.h"
#include "msprime.h"
#include "object_heap.h"

static int
cmp_mutation(const void *a, const void *b) {
    const mutation_t *ia = (const mutation_t *) a;
    const mutation_t *ib = (const mutation_t *) b;
    return (ia->position > ib->position) - (ia->position < ib->position);
}

static void
mutgen_check_state(mutgen_t *self)
{
    /* TODO some checks! */
}

void
mutgen_print_state(mutgen_t *self, FILE *out)
{
    size_t j, k;

    fprintf(out, "Mutgen state\n");
    fprintf(out, "\tmutation_rate = %f\n", (double) self->mutation_rate);
    fprintf(out, "\tsequence_length = %f\n", (double) self->sequence_length);
    fprintf(out, "\tmutation_block_size = %d\n", (int) self->mutation_block_size);
    fprintf(out, "\tmax_num_mutations  = %d\n", (int) self->max_num_mutations);
    fprintf(out, "\tMUTATIONS\t%d\n", (int) self->num_mutations);
    for (j = 0; j < self->num_mutations; j++) {
        fprintf(out, "\t%f\t", self->mutations[j].position);
        for (k = 0; k < self->mutations[j].num_nodes; k++) {
            fprintf(out, "%d,", self->mutations[j].nodes[k]);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\tnode_heap  = \n");
    object_heap_print_state(&self->node_heap, out);

    mutgen_check_state(self);
}


int WARN_UNUSED
mutgen_alloc(mutgen_t *self, tree_sequence_t *tree_sequence,
        double mutation_rate, gsl_rng *rng)
{
    int ret = MSP_ERR_NO_MEMORY;
    uint32_t j;
    sample_t sample;

    assert(tree_sequence != NULL);
    assert(rng != NULL);
    memset(self, 0, sizeof(mutgen_t));
    self->mutation_rate = mutation_rate;
    self->tree_sequence = tree_sequence;
    self->sequence_length = tree_sequence_get_sequence_length(tree_sequence);
    self->rng = rng;
    self->num_mutations = 0;
    self->mutation_block_size = 1024 * 1024;
    /* Avoid potential portability issues with realloc(NULL, newsize)
     * by mallocing enough space for 1 mutation initiall. This gives the user
     * control over the overall malloc behavior.
     */
    self->max_num_mutations = 1;
    self->mutations = malloc(sizeof(mutation_t));
    if (self->mutations == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = object_heap_init(&self->node_heap, sizeof(uint32_t),
            self->mutation_block_size, NULL);
    if (ret != 0) {
        goto out;
    }
    self->times = calloc(
        tree_sequence_get_num_nodes(tree_sequence) + 1, sizeof(double));
    if (self->times == NULL) {
        goto out;
    }
    /* Set the initial times for samples */
    for (j = 0; j < tree_sequence_get_sample_size(tree_sequence); j++) {
        ret = tree_sequence_get_sample(tree_sequence, j, &sample);
        if (ret != 0) {
            goto out;
        }
        self->times[j] = sample.time;
    }
out:
    return ret;
}

int
mutgen_free(mutgen_t *self)
{
    if (self->mutations != NULL) {
        free(self->mutations);
    }
    if (self->times != NULL) {
        free(self->times);
    }
    object_heap_free(&self->node_heap);
    return 0;
}

int WARN_UNUSED
mutgen_set_mutation_block_size(mutgen_t *self, size_t mutation_block_size)
{
    int ret = 0;
    if (mutation_block_size == 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->mutation_block_size = mutation_block_size;
out:
    return ret;
}

static int WARN_UNUSED
mutgen_add_mutation(mutgen_t *self, uint32_t node, double position)
{
    int ret = 0;
    mutation_t *tmp_buffer;
    uint32_t *p;

    assert(self->num_mutations <= self->max_num_mutations);

    if (self->num_mutations == self->max_num_mutations) {
        self->max_num_mutations += self->mutation_block_size;
        tmp_buffer = realloc(self->mutations,
            self->max_num_mutations * sizeof(mutation_t));
        if (tmp_buffer == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        self->mutations = tmp_buffer;
    }
    if (object_heap_empty(&self->node_heap)) {
        ret = object_heap_expand(&self->node_heap);
        if ( ret != 0) {
            goto out;
        }
    }
    p = (uint32_t *) object_heap_alloc_object(&self->node_heap);
    self->mutations[self->num_mutations].nodes = p;
    self->mutations[self->num_mutations].nodes[0] = node;
    self->mutations[self->num_mutations].num_nodes = 1;
    self->mutations[self->num_mutations].position = position;
    self->num_mutations++;
out:
    return ret;
}

static int WARN_UNUSED
mutgen_generate_record_mutations(mutgen_t *self, coalescence_record_t *cr)
{
    int ret = -1;
    size_t k, l, branch_mutations;
    double branch_length, position, mu;
    double distance = cr->right - cr->left;
    uint32_t child;

    self->times[cr->node] = cr->time;
    for (k = 0; k < cr->num_children; k++) {
        child = cr->children[k];
        branch_length = cr->time - self->times[child];
        mu = branch_length * distance * self->mutation_rate;
        branch_mutations = gsl_ran_poisson(self->rng, mu);
        for (l = 0; l < branch_mutations; l++) {
            position = gsl_ran_flat(self->rng, cr->left, cr->right);
            assert(cr->left <= position && position < cr->right);
            ret = mutgen_add_mutation(self, child, position);
            if (ret != 0) {
                goto out;
            }
        }
    }
    ret = 0;
out:
    return ret;
}

int WARN_UNUSED
mutgen_generate(mutgen_t *self)
{
    int ret = -1;
    tree_sequence_t *ts = self->tree_sequence;
    coalescence_record_t cr;
    size_t j;

    for (j = 0; j < tree_sequence_get_num_coalescence_records(ts); j++) {
        ret = tree_sequence_get_coalescence_record(ts, j, &cr, MSP_ORDER_TIME);
        if (ret != 0) {
            goto out;
        }
        ret = mutgen_generate_record_mutations(self, &cr);
        if (ret != 0) {
            goto out;
        }
    }
    qsort(self->mutations, self->num_mutations, sizeof(mutation_t),
                cmp_mutation);
    ret = 0;
out:
    return ret;
}

size_t
mutgen_get_num_mutations(mutgen_t *self)
{
    return self->num_mutations;
}

int  WARN_UNUSED
mutgen_get_mutations(mutgen_t *self, mutation_t **mutations)
{
    *mutations = self->mutations;
    return 0;
}

