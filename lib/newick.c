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
    self->num_loci = self->tree_file.num_loci;
    N  = 2 * self->sample_size;
    self->tau = malloc(N * sizeof(float));
    self->branch_lengths = malloc(N * sizeof(float));
    self->children = malloc(N * sizeof(int *));
    self->visited = malloc(N * sizeof(int));
    self->stack_size = self->sample_size;
    self->stack = malloc(self->stack_size * sizeof(int));
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
    memset(self->children, 0, N * sizeof(int *));
    memset(self->tau, 0, N * sizeof(float));
    memset(self->branch_lengths, 0, N * sizeof(float));
    p = self->children_mem;
    for (j = self->sample_size + 1; j < N; j++) {
        self->children[j] = p;
        p += 2;
    }
    self->breakpoint = 1;
    self->num_trees = 0;
    self->completed = 0;
    /* start off the process by getting the first record */
    ret = tree_file_next_record(&self->tree_file, &self->next_record);
    if (ret < 0) {
        goto out;
    }
    /* We must always have 1 record, so this is an error condition */
    assert(ret == 1);
    ret = 0;
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

/* Get the next tree from the treefile */
static int newick_update_tree(newick_t *self, uint32_t *length)
{
    int ret = 0;
    int tree_ret = 1;
    int32_t c1, c2, p;

    while (self->next_record.left == self->breakpoint && tree_ret == 1) {
        /*
        printf("examining: %d\t(%d, %d)->%d @ %f\n",
                self->next_record.left, self->next_record.children[0],
                self->next_record.children[1], self->next_record.parent, self->next_record.time);
        */
        self->tau[self->next_record.parent] = self->next_record.time;
        c1 = self->next_record.children[0];
        c2 = self->next_record.children[1];
        p = self->next_record.parent;
        if (c1 < c2) {
            self->children[p][0] = c1;
            self->children[p][1] = c2;
        } else {
            self->children[p][1] = c1;
            self->children[p][0] = c2;
        }
        self->branch_lengths[c1] = self->next_record.time - self->tau[c1];
        self->branch_lengths[c2] = self->next_record.time - self->tau[c2];
        tree_ret = tree_file_next_record(&self->tree_file, &self->next_record);
    }
    if (tree_ret < 0) {
        goto out;
    }
    self->num_trees++;
    if (tree_ret == 1) {
        *length = self->next_record.left - self->breakpoint;
        self->breakpoint = self->next_record.left;
    } else {
        self->completed = 1;
        *length = self->num_loci - self->breakpoint + 1;
    }
    ret = 1;
out:
    return ret;
}

/*
 * Generates the newick string from the tree in memory
 */
static int
newick_generate_string(newick_t *self, size_t *output_length)
{
    int ret = 0;
    uint32_t u;
    int32_t n = self->sample_size;
    int32_t root = 2 * n - 1;
    int **c = self->children;
    float *branch_lengths = self->branch_lengths;
    int *stack = self->stack;
    int *visited = self->visited;
    int stack_top = 0;
    char *s = self->output_buffer;
    size_t max_length = self->output_buffer_size;
    size_t length = 0;
    size_t r;

    /* TODO this needs test cases to make sure all the corner cases
     * work properly!
     */
    stack[0] = root;
    stack_top = 0;
    memset(visited, 0, 2 * n * sizeof(int32_t));
    while (stack_top >= 0) {
        u = stack[stack_top];
        stack_top--;
        if (c[u] == NULL) {
            /* leaf node */
            r = snprintf(s + length, max_length - length, "%d:%.3f",
                    u, branch_lengths[u]);
            length += r;
            if (length >= max_length) {
                ret = MSP_ERR_NEWICK_OVERFLOW;
                goto out;
            }
        } else {
            if (visited[u] == 0) {
                s[length] = '(';
                length++;
                if (length >= max_length) {
                    ret = MSP_ERR_NEWICK_OVERFLOW;
                    goto out;
                }
                stack_top++;
                assert(stack_top < self->stack_size);
                stack[stack_top] = u;
                stack_top++;
                assert(stack_top < self->stack_size);
                stack[stack_top] = c[u][0];
            } else if (visited[u] == 1) {
                s[length] = ',';
                length++;
                if (length >= max_length) {
                    ret = MSP_ERR_NEWICK_OVERFLOW;
                    goto out;
                }
                stack_top++;
                assert(stack_top < self->stack_size);
                stack[stack_top] = u;
                stack_top++;
                assert(stack_top < self->stack_size);
                stack[stack_top] = c[u][1];
            } else {
                s[length] = ')';
                length++;
                if (length >= max_length) {
                    ret = MSP_ERR_NEWICK_OVERFLOW;
                    goto out;
                }
                if (u == root) {
                    s[length] = ';';
                    length++;
                } else {
                    r = snprintf(s + length, max_length - length, ":%.3f",
                            branch_lengths[u]);
                    length += r;
                }
                if (length >= max_length) {
                    ret = MSP_ERR_NEWICK_OVERFLOW;
                    goto out;
                }
            }
            visited[u]++;
        }
    }
    s[length] = '\0';
    *output_length = length;
out:
    return ret;
}

int
newick_next_tree(newick_t *self, uint32_t *tree_length, char **tree,
        size_t *str_length)
{
    int ret = -1;
    int newick_ret;

    if (self->completed) {
        ret = 0;
    } else {
        ret = newick_update_tree(self, tree_length);
        if (ret < 0) {
            goto out;
        }
        newick_ret = newick_generate_string(self, str_length);
        if (newick_ret != 0) {
            ret = newick_ret;
            goto out;
        }
        *tree = self->output_buffer;
    }
out:
    return ret;
}

int
newick_output_ms_format(newick_t *self, FILE *out)
{
    int ret = -1;
    int io_ret;
    char *tree = NULL;
    uint32_t l;
    size_t str_len;

    while ((ret = newick_next_tree(self, &l, &tree, &str_len)) == 1) {
        io_ret = fprintf(out, "[%d]", l);
        if (io_ret < 0) {
            ret = MSP_ERR_IO;
            goto out;
        }
        io_ret = fprintf(out, tree);
        if (io_ret < 0) {
            ret = MSP_ERR_IO;
            goto out;
        }
        io_ret = fprintf(out, "\n");
        if (io_ret < 0) {
            ret = MSP_ERR_IO;
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}
