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

/*
 * Binary index tree (also known as a Fenwick tree) implementation.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "err.h"
#include "fenwick.h"

static void
fenwick_set_log_size(fenwick_t *self)
{
    size_t u = self->size;

    while (u != 0) {
        self->log_size = u;
        u -= (u & -u);
    }
}

static int WARN_UNUSED
fenwick_alloc_buffers(fenwick_t *self)
{
    int ret = -1;

    self->tree = NULL;
    self->values = NULL;
    self->tree = calloc((1 + self->size), sizeof(int64_t));
    if (self->tree == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->values = calloc((1 + self->size), sizeof(int64_t));
    if (self->values == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    fenwick_set_log_size(self);
    ret = 0;
out:
    return ret;
}

int WARN_UNUSED
fenwick_alloc(fenwick_t *self, size_t initial_size)
{
    self->size = initial_size;
    return fenwick_alloc_buffers(self);
}

int WARN_UNUSED
fenwick_expand(fenwick_t *self, size_t increment)
{
    int ret = -1;
    size_t j;
    int64_t *values = self->values;
    size_t old_size = self->size;

    if (self->tree == NULL || self->values == NULL) {
        goto out;
    }
    self->size += increment;
    free(self->tree);
    /* we still need values until we have rebuilt the tree; we need to
     * make sure we don't double free in an error condition though.
     */
    self->values = NULL;
    ret = fenwick_alloc_buffers(self);
    if (ret != 0) {
        goto out;
    }
    /* now insert all of the values into the new tree. There is probably a
     * better way to do this...
     */
    for (j = 1; j <= old_size; j++) {
        fenwick_set_value(self, j, values[j]);
    }
out:
    if (values != NULL) {
        free(values);
    }
    return ret;
}

int
fenwick_free(fenwick_t *self)
{
    int ret = -1;

    if (self->tree != NULL) {
        free(self->tree);
        self->tree = NULL;
    }
    if (self->values != NULL) {
        free(self->values);
        self->values = NULL;
    }
    ret = 0;
    return ret;
}

size_t
fenwick_get_size(fenwick_t *self)
{
    return self->size;
}

int64_t
fenwick_get_total(fenwick_t *self)
{
    return fenwick_get_cumulative_sum(self, self->size);
}

void
fenwick_increment(fenwick_t *self, size_t index, int64_t value)
{
    size_t j = index;

    assert(0 < index && index <= self->size);
    self->values[index] += value;
    while (j <= self->size) {
        self->tree[j] += value;
        j += (j & -j);
    }
}

void
fenwick_set_value(fenwick_t *self, size_t index, int64_t value)
{
    int64_t v = value - self->values[index];

    fenwick_increment(self, index, v);
}

int64_t
fenwick_get_cumulative_sum(fenwick_t *self, size_t index)
{
    int64_t ret = 0;
    size_t j = index;

    assert(0 < index && index <= self->size);
    while (j > 0) {
        ret += self->tree[j];
        j -= (j & -j);
    }
    return ret;
}

int64_t
fenwick_get_value(fenwick_t *self, size_t index)
{
    assert(0 < index && index <= self->size);
    return self->values[index];
}


size_t
fenwick_find(fenwick_t *self, int64_t sum)
{
    size_t j = 0;
    size_t k;
    int64_t s = sum;
    size_t half = self->log_size;

    while (half > 0) {
        /* Skip non-existent entries */
        while (j + half > self->size) {
            half >>= 1;
        }
        k = j + half;
        if (s > self->tree[k]) {
            j = k;
            s -= self->tree[j];
        }
        half >>= 1;
    }
    return j + 1;
}
