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
    /* uint32_t j; */
    /* avl_node_t *avl_node; */
    /* hapgen_tree_node_t *node = NULL; */
    /* hapgen_tree_node_t search; */
    /* double tbl = 0.0; */

    /* assert(avl_count(&self->tree) == 2 * self->sample_size - 1); */
    /* for (j = 1; j <= self->sample_size; j++) { */
    /*     search.id = j; */
    /*     avl_node = avl_search(&self->tree, &search); */
    /*     assert(avl_node != NULL); */
    /*     node = (hapgen_tree_node_t *) avl_node->item; */
    /*     while (node->parent != NULL) { */
    /*         node = node->parent; */
    /*     } */
    /*     assert(node != NULL); */
    /*     assert(node == self->root); */
    /* } */
    /* for (avl_node = self->tree.head; avl_node != NULL; avl_node = avl_node->next) { */
    /*     node = (hapgen_tree_node_t *) avl_node->item; */
    /*     tbl += node->branch_length; */
    /* } */
    /* if (gsl_fcmp(tbl, self->total_branch_length, 1e-8) != 0) { */
    /*     printf("tbl %f\t%f\n", tbl, self->total_branch_length); */
    /*     assert(tbl == self->total_branch_length); */
    /* } */
}

void
hapgen_print_state(hapgen_t *self)
{
    size_t j;
    avl_node_t *avl_node;
    avl_tree_t *mutations;
    mutation_t *u;

    printf("Hapgen state\n");
    printf("num_nodes = %d\n", (int) self->num_nodes);
    printf("num_mutations = %d\n", (int) self->num_mutations);
    for (j = 1; j <= self->num_nodes; j++) {
        printf("\t%d\n", (int) j);
        mutations = &self->mutations[j];
        for (avl_node = mutations->head; avl_node != NULL; avl_node = avl_node->next) {
            u = (mutation_t *) avl_node->item;
            printf("\t\t%f @ %d\n", u->position, (int) u->site);
        }
    }
    printf("avl_node_heap\n");
    object_heap_print_state(&self->avl_node_heap);
    hapgen_check_state(self);
}

static inline avl_node_t * WARN_UNUSED
hapgen_alloc_avl_node(hapgen_t *self, double position)
{
    avl_node_t *ret = NULL;
    mutation_t *u;
    char *p;

    if (object_heap_empty(&self->avl_node_heap)) {
        goto out;
    }
    p = object_heap_alloc_object(&self->avl_node_heap);
    if (p == NULL) {
        goto out;
    }
    /* We have alloced two objects at once: the avl_node and the
     * mutation that goes with it.
     */
    ret = (avl_node_t *) p;
    u = (mutation_t *) (p + sizeof(avl_node_t));
    u->position = position;
    u->site = 0;
    avl_init_node(ret, u);
out:
    return ret;
}

static int
hapgen_add_mutation(hapgen_t *self, uint32_t node, double position)
{
    int ret = -1;
    avl_node_t *avl_node = hapgen_alloc_avl_node(self, position);

    if (avl_node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    avl_node = avl_insert_node(&self->mutations[node], avl_node);
    /* We assume that we never generate conflicting positions. */
    assert(avl_node != NULL);
    ret = 0;
out:
    return ret;
}

static int
hapgen_initialise_mutations(hapgen_t *self)
{
    int ret = -1;
    size_t j;
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
        hapgen_add_mutation(self, nodes[j], positions[j]);
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
hapgen_assign_mutation_sites(hapgen_t *self)
{
    int ret = -1;
    mutation_t **all_mutations = NULL;
    size_t j, k;
    avl_node_t *avl_node;
    avl_tree_t *mutations;

    /* sort the mutations so we can assign site numbers. */
    all_mutations = malloc(self->num_mutations * sizeof(mutation_t *));
    if (all_mutations == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    k = 0;
    for (j = 1; j <= self->num_nodes; j++) {
        mutations = &self->mutations[j];
        for (avl_node = mutations->head; avl_node != NULL;
                avl_node = avl_node->next) {
            all_mutations[k] = (mutation_t *) avl_node->item;
            k++;
        }
    }
    assert(k == self->num_mutations);
    qsort(all_mutations, k, sizeof(mutation_t *), cmp_mutation_pointer);
    /* now assign the sites and fill the positions array */
    for (j = 0; j < self->num_mutations; j++) {
        all_mutations[j]->site = j;
        self->positions[j] = all_mutations[j]->position;
    }

    ret = 0;
out:
    if (all_mutations != NULL) {
        free(all_mutations);
    }
    return ret;
}

static inline double
_get_position(avl_node_t *node)
{
    mutation_t *mut;
    assert(node != NULL);
    mut = (mutation_t *) node->item;
    assert(mut != NULL);
    return mut->position;
}


static int
hapgen_apply_node_mutations(hapgen_t *self, uint32_t sample_id, uint32_t node,
        uint32_t left, uint32_t right)
{
    int ret = 0;
    int where;
    avl_node_t *avl_node, *found;
    mutation_t *mut, search;
    uint32_t haplotype_index;

    haplotype_index = 0;
    if (self->mode == MSP_HAPGEN_MODE_ALL) {
        haplotype_index = sample_id - 1;
    }
    search.position = (double) left;
    where = avl_search_closest(&self->mutations[node], &search, &found);
    assert(found != NULL);
    avl_node = found;
    mut = (mutation_t *) avl_node->item;
    /* printf("found %d: %f\n", where, mut->position); */
    if (where < 0) {
        /* The closest position is > left, so we start from there */
        assert(_get_position(avl_node) >= left);
    } else {
        /* the closest position is <= left, so we start from either the
         * current node or the next one. */
        assert(_get_position(avl_node) <= left);
        if (_get_position(avl_node) < left) {
            avl_node = avl_node->next;
        }
    }
    while (avl_node != NULL
            && _get_position(avl_node) >= left
            && _get_position(avl_node) < right) {
        mut = (mutation_t *) avl_node->item;
        /* printf("\tapplying %f \n", mut->position); */
        self->haplotypes[haplotype_index][mut->site] = '1';
        assert(left <= mut->position && mut->position < right);
        avl_node = avl_node->next;
    }

/*     /1* TMP: check that this has been done properly *1/ */
/*     for (avl_node = self->mutations[node].head; avl_node != NULL; */
/*             avl_node = avl_node->next) { */
/*         mut = (mutation_t *) avl_node->item; */
/*         assert(mut != NULL); */
/*         if (left <= mut->position && mut->position < right) { */
/*             if (self->haplotype[mut->site] != '1') { */
/*                 printf("Error at site %d: %f\n", (int) mut->site, mut->position); */
/*                 assert(1 == 0); */
/*             } */
/*         } */
/*     } */
    return ret;
}

static int
hapgen_apply_tree_mutations(hapgen_t *self, uint32_t sample_id,
        sparse_tree_t *tree, uint32_t left, uint32_t right)
{
    int ret = 0;
    uint32_t u;

    u = sample_id;
    while (u != 0) {
        if (avl_count(&self->mutations[u]) > 0) {
            ret = hapgen_apply_node_mutations(self, sample_id, u, left, right);
            if (ret != 0) {
                goto out;
            }
        }
        u = tree->parent[u];
    }
out:
    return ret;
}

static int
hapgen_apply_all_tree_mutations(hapgen_t *self, sparse_tree_t *tree,
        uint32_t left, uint32_t right)
{
    int ret = 0;
    uint32_t j;

    for (j = 1; j <= self->sample_size; j++) {
        ret = hapgen_apply_tree_mutations(self, j, tree, left, right);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}


static int
hapgen_generate_all_haplotypes(hapgen_t *self)
{
    int ret = 0;
    uint32_t length, left;
    sparse_tree_t *tree;

    sparse_tree_iterator_reset(&self->tree_iterator);
    left = 0;
    while ((ret = sparse_tree_iterator_next(
                    &self->tree_iterator, &length, &tree)) == 1) {
        ret = hapgen_apply_all_tree_mutations(self, tree, left,
                left + length);
        if (ret != 0) {
            goto out;
        }
        left += length;
    }
    if (ret != 0) {
        goto out;
    }

out:
    return ret;
}

int
hapgen_alloc(hapgen_t *self, tree_sequence_t *tree_sequence, int mode)
{
    int ret = MSP_ERR_NO_MEMORY;
    uint32_t j, num_to_alloc;

    assert(tree_sequence != NULL);
    memset(self, 0, sizeof(hapgen_t));
    if (mode != MSP_HAPGEN_MODE_SINGLE || mode != MSP_HAPGEN_MODE_ALL) {
        ret = MSP_ERR_BAD_MODE;
    }
    self->mode = mode;
    self->sample_size = tree_sequence_get_sample_size(tree_sequence);
    self->num_nodes = tree_sequence_get_num_nodes(tree_sequence);
    self->num_loci = tree_sequence_get_num_loci(tree_sequence);
    self->num_mutations = tree_sequence_get_num_mutations(tree_sequence);
    self->tree_sequence = tree_sequence;

    ret = object_heap_init(&self->avl_node_heap,
            sizeof(avl_node_t) + sizeof(mutation_t), self->num_mutations,
            NULL);
    if (ret != 0) {
        goto out;
    }
    ret = sparse_tree_iterator_alloc(&self->tree_iterator,
            self->tree_sequence);
    if (ret != 0) {
        goto out;
    }
    self->mutations = malloc((self->num_nodes + 1) * sizeof(avl_tree_t));
    if (self->mutations == NULL) {
        goto out;
    }
    for (j = 1; j <= self->num_nodes; j++) {
        avl_init_tree(&self->mutations[j], cmp_mutation, NULL);
    }
    ret = hapgen_initialise_mutations(self);
    if (ret != 0) {
        goto out;
    }
    /* allocate the list of positions and the haplotype(s) */
    self->positions = malloc(self->num_mutations * sizeof(double));
    num_to_alloc = 1;
    if (self->mode == MSP_HAPGEN_MODE_ALL) {
        num_to_alloc = self->sample_size;
    }
    self->haplotype_mem = malloc(
            num_to_alloc * (self->num_mutations + 1) * sizeof(char));
    self->haplotypes = malloc(num_to_alloc * sizeof(char *));
    if (self->haplotype_mem == NULL || self->haplotypes == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memset(self->haplotype_mem, '0',
            num_to_alloc * (self->num_mutations + 1) * sizeof(char));
    for (j = 0; j < num_to_alloc; j++) {
        self->haplotypes[j] = self->haplotype_mem
            + j * (self->num_mutations + 1);
        self->haplotypes[j][self->num_mutations] = '\0';
    }
    ret = hapgen_assign_mutation_sites(self);
    if (ret != 0) {
        goto out;
    }
    if (self->mode == MSP_HAPGEN_MODE_ALL) {
        hapgen_generate_all_haplotypes(self);
    }
    /* hapgen_print_state(self); */
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
    if (self->positions != NULL) {
        free(self->positions);
    }
    if (self->haplotypes != NULL) {
        free(self->haplotypes);
    }
    if (self->haplotype_mem != NULL) {
        free(self->haplotype_mem);
    }
    object_heap_free(&self->avl_node_heap);
    sparse_tree_iterator_free(&self->tree_iterator);
    return 0;
}

static int
hapgen_get_haplotype_mode_single(hapgen_t *self, uint32_t sample_id,
        char **haplotype)
{
    int ret = 0;
    uint32_t length, left;
    sparse_tree_t *tree;

    memset(self->haplotypes[0], '0', self->num_mutations);
    sparse_tree_iterator_reset(&self->tree_iterator);
    left = 0;
    while ((ret = sparse_tree_iterator_next(
                    &self->tree_iterator, &length, &tree)) == 1) {
        ret = hapgen_apply_tree_mutations(self, sample_id, tree,
                left, left + length);
        if (ret != 0) {
            goto out;
        }
        left += length;
    }
    if (ret != 0) {
        goto out;
    }

    *haplotype = self->haplotypes[0];
out:
    return ret;
}

int
hapgen_get_haplotype(hapgen_t *self, uint32_t sample_id, char **haplotype)
{
    int ret = 0;

    if (sample_id < 1 || sample_id> self->sample_size) {
        ret = MSP_ERR_OUT_OF_BOUNDS;
        goto out;
    }
    if (self->mode == MSP_HAPGEN_MODE_SINGLE) {
        ret = hapgen_get_haplotype_mode_single(self, sample_id, haplotype);
    } else {
        *haplotype = self->haplotypes[sample_id - 1];
    }
out:
    return ret;
}

size_t
hapgen_get_num_segregating_sites(hapgen_t *self)
{
    return self->num_mutations;
}
