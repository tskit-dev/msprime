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

#define HG_WORD_SIZE 64

static int
cmp_mutation(const void *a, const void *b) {
    const mutation_t *ia = (const mutation_t *) a;
    const mutation_t *ib = (const mutation_t *) b;
    return (ia->position > ib->position) - (ia->position < ib->position);
}

static int
cmp_mutation_pointer(const void *a, const void *b) {
    mutation_t *const*ia = (mutation_t *const*) a;
    mutation_t *const*ib = (mutation_t *const*) b;
    return cmp_mutation(*ia, *ib);
}

/* Ensure the tree is in a consistent state */
static void
hapgen_check_state(hapgen_t *self)
{
    /* TODO some checks! */
}

void
hapgen_print_state(hapgen_t *self)
{
    size_t j, k;
    mutation_t *u;

    printf("Hapgen state\n");
    printf("num_mutations = %d\n", (int) self->num_mutations);
    printf("words_per_row = %d\n", (int) self->words_per_row);
    for (j = 0; j < self->num_mutations; j++) {
        u = self->mutations[j];
        printf("\t\t%f @ %d\n", u->position, (int) u->site);
    }
    printf("haplotype matrix\n");
    for (j = 0; j < self->sample_size; j++) {
        for (k = 0; k < self->words_per_row; k++) {
            printf("%llu ", (unsigned long long)
                    self->haplotype_matrix[j * self->words_per_row + k]);
        }
        printf("\n");
    }
    hapgen_check_state(self);
}


static inline int
hapgen_set_bit(hapgen_t *self, size_t row, size_t column)
{
    /* get the word that column falls in */
    size_t word = column / HG_WORD_SIZE;
    size_t bit = column % HG_WORD_SIZE;

    assert(word < self->words_per_row);
    self->haplotype_matrix[row * self->words_per_row + word] |= 1ULL << bit;
    return 0;
}

static int
hapgen_initialise_mutations(hapgen_t *self)
{
    int ret = -1;
    size_t j;
    mutation_t *u;
    uint32_t *nodes = NULL;
    double *positions = NULL;

    nodes = malloc(self->num_mutations * sizeof(uint32_t));
    positions = malloc(self->num_mutations * sizeof(double));
    if (nodes == NULL || positions == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = tree_sequence_get_mutations(self->tree_sequence, nodes,
            positions);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < self->num_mutations; j++) {
        u = &self->mutation_mem[j];
        self->mutations[j] = u;
        self->mutations[j]->node = nodes[j];
        self->mutations[j]->position = positions[j];
    }
    /* Now we sort by position so we can assign site numbers. */
    qsort(self->mutations, self->num_mutations, sizeof(mutation_t *),
            cmp_mutation_pointer);
    for (j = 0; j < self->num_mutations; j++) {
        self->mutations[j]->site = j;
    }
    ret = 0;
out:
    if (nodes != NULL) {
        free(nodes);
    }
    if (positions != NULL) {
        free(positions);
    }
    return ret;
}


static int
hapgen_apply_tree_mutation(hapgen_t *self, sparse_tree_t *tree,
        mutation_t *mut)
{
    int ret = 0;
    uint32_t *stack = self->traversal_stack;
    uint32_t u, c;
    int stack_top = 0;

    stack[0] = mut->node;
    while (stack_top >= 0) {
        u = stack[stack_top];
        stack_top--;
        if (tree->children[2 * u] == 0) {
            hapgen_set_bit(self, u - 1, mut->site);
        } else {
            for (c = 0; c < 2; c++) {
                stack_top++;
                stack[stack_top] = tree->children[2 * u + c];
            }
        }
    }
    return ret;
}

static int
hapgen_generate_all_haplotypes(hapgen_t *self)
{
    int ret = 0;
    sparse_tree_t *tree;
    size_t j;

    j = 0;
    while ((ret = sparse_tree_iterator_next(
            &self->tree_iterator, &tree)) == 1) {
        while (j < self->num_mutations
                && self->mutations[j]->position < tree->right) {
            ret = hapgen_apply_tree_mutation(self, tree,
                    self->mutations[j]);
            if (ret != 0) {
                goto out;
            }
            j++;
        }
    }
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int
hapgen_alloc(hapgen_t *self, tree_sequence_t *tree_sequence)
{
    int ret = MSP_ERR_NO_MEMORY;

    assert(tree_sequence != NULL);
    memset(self, 0, sizeof(hapgen_t));
    self->sample_size = tree_sequence_get_sample_size(tree_sequence);
    self->num_loci = tree_sequence_get_num_loci(tree_sequence);
    self->num_mutations = tree_sequence_get_num_mutations(tree_sequence);
    self->tree_sequence = tree_sequence;

    ret = sparse_tree_iterator_alloc(&self->tree_iterator,
            self->tree_sequence);
    if (ret != 0) {
        goto out;
    }
    self->traversal_stack = malloc(self->sample_size * sizeof(uint32_t));
    if (self->traversal_stack == NULL) {
        goto out;
    }
    self->mutations = malloc(self->num_mutations * sizeof(mutation_t *));
    self->mutation_mem = malloc(self->num_mutations * sizeof(mutation_t));
    if (self->mutations == NULL || self->mutation_mem == NULL) {
        goto out;
    }
    /* set up the haplotype binary matrix */
    /* The number of words per row is the number of mutations divided by 64 */
    self->words_per_row = (self->num_mutations / HG_WORD_SIZE) + 1;
    self->haplotype_matrix = calloc(self->words_per_row * self->sample_size,
            sizeof(uint64_t));
    self->haplotype = malloc(self->words_per_row * HG_WORD_SIZE + 1);
    if (self->haplotype_matrix == NULL || self->haplotype == NULL) {
        goto out;
    }
    ret = hapgen_initialise_mutations(self);
    if (ret != 0) {
        goto out;
    }
    ret = hapgen_generate_all_haplotypes(self);
    if (ret != 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

int
hapgen_free(hapgen_t *self)
{
    if (self->mutations != NULL) {
        free(self->mutations);
    }
    if (self->mutation_mem != NULL) {
        free(self->mutation_mem);
    }
    if (self->haplotype_matrix != NULL) {
        free(self->haplotype_matrix);
    }
    if (self->haplotype != NULL) {
        free(self->haplotype);
    }
    if (self->traversal_stack != NULL) {
        free(self->traversal_stack);
    }
    sparse_tree_iterator_free(&self->tree_iterator);
    return 0;
}


int
hapgen_get_haplotype(hapgen_t *self, uint32_t sample_id, char **haplotype)
{
    int ret = 0;
    size_t j, k, l;
    uint64_t word;

    if (sample_id < 1 || sample_id > self->sample_size) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    l = 0;
    for (j = 0; j < self->words_per_row; j++) {
        word = self->haplotype_matrix[(sample_id - 1) * self->words_per_row + j];
        for (k = 0; k < HG_WORD_SIZE; k++) {
            self->haplotype[l] = (word >> k) & 1ULL ? '1': '0';
            l++;
        }
    }
    self->haplotype[self->num_mutations] = '\0';
    *haplotype = self->haplotype;
out:
    return ret;
}

size_t
hapgen_get_num_segregating_sites(hapgen_t *self)
{
    return self->num_mutations;
}
