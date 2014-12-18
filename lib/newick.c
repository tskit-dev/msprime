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
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "err.h"
#include "msprime.h"

int
newick_alloc(newick_t *self, const char *tree_file_name)
{
    int ret = -1;
    size_t j, N;
    int *p;

    memset(self, 0, sizeof(newick_t));
    ret = tree_file_open(&self->tree_file, tree_file_name, 'r');
    if (ret != 0) {
        goto out;
    }
    if (!tree_file_issorted(&self->tree_file)) {
        ret = MSP_ERR_TREE_FILE_NOT_SORTED;
        goto out;
    }
    self->sample_size = self->tree_file.sample_size;
    N  = 2 * self->sample_size;
    self->tau = malloc(N * sizeof(float));
    self->branch_lengths = malloc(N * sizeof(float));
    self->children = malloc(N * sizeof(int *));
    self->visited = malloc(N * sizeof(int));
    self->stack = malloc(self->sample_size * sizeof(int));
    /* TODO fix this! should be function of n */
    self->output_buffer_size = 16 * 1024 * 1024;
    self->output_buffer = malloc(self->output_buffer_size);
    self->children_mem = malloc(N * sizeof(int));
    if (self->tau == NULL || self->branch_lengths == NULL
            || self->children == NULL || self->visited == NULL
            || self->stack == NULL || self->output_buffer == NULL
            || self->children_mem == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    /* Set the pointers for the child storage. We could save some memory
     * here by not keeping the null pointers for 1 to n.
     */
    p = self->children_mem;
    for (j = self->sample_size; j < N; j++) {
        self->children[j] = p;
        p += 2;
    }
out:
    return ret;
}

int
newick_free(newick_t *self)
{
    /* TODO set all these to NULL to protect against double free. */
    if (self->tau != NULL) {
        free(self->tau);
    }
    if (self->branch_lengths != NULL) {
        free(self->branch_lengths);
    }
    if (self->children != NULL) {
        free(self->children);
    }
    if (self->visited != NULL) {
        free(self->visited);
    }
    if (self->stack != NULL) {
        free(self->stack);
    }
    if (self->output_buffer != NULL) {
        free(self->output_buffer);
    }
    if (self->children_mem != NULL) {
        free(self->children_mem);
    }
    tree_file_close(&self->tree_file);
    return 0;
}

int
newick_next_tree(newick_t *self, uint32_t *l, char **tree)
{
    int ret = -1;
    return ret;
}
