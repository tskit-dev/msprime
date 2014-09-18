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

int WARN_UNUSED
fenwick_alloc(fenwick_t *self)
{
    int ret = -1;
    unsigned int u = self->max_index;
    while (u != 0) {
        self->log_max_index = u;
        u -= (u & -u);
    }
    self->tree = NULL;
    self->values = NULL;
    self->tree = calloc((1 + self->max_index), sizeof(long long));
    if (self->tree == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->values = calloc((1 + self->max_index), sizeof(long long));
    if (self->values == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = 0;
out:
    return ret;
}

int
fenwick_free(fenwick_t *self)
{
    int ret = -1;
    if (self->tree != NULL) {
        free(self->tree);
    }
    if (self->values != NULL) {
        free(self->values);
    }
    ret = 0;
    return ret;
}

long long
fenwick_get_total(fenwick_t *self)
{
    long long ret = fenwick_get_cumulative_sum(self, self->max_index);
    return ret;
}

void
fenwick_increment(fenwick_t *self, unsigned int index, long long value)
{
    unsigned int j = index;
    assert(0 < index && index <= self->max_index);
    self->values[index] += value;
    while (j <= self->max_index) {
        self->tree[j] += value;
        j += (j & -j);
    }
}

void
fenwick_set_value(fenwick_t *self, unsigned int index, long long value)
{
    long long v = value - self->values[index];
    fenwick_increment(self, index, v);
}

long long fenwick_get_cumulative_sum(fenwick_t *self, unsigned int index)
{
    long long ret = 0;
    unsigned int j = index;
    assert(0 < index && index <= self->max_index);
    while (j > 0) {
        ret += self->tree[j];
        j -= (j & -j);
    }
    return ret;

}

long long fenwick_get_value(fenwick_t *self, unsigned int index)
{
    assert(0 < index && index <= self->max_index);
    return self->values[index];
    /*
    long long ret = 0;
    unsigned int j = index;
    unsigned int p = j & (j - 1);
    assert(0 < index && index <= self->max_index);
    ret = self->tree[j];
    j--;
    while (p != j) {
        ret -= self->tree[j];
        j = j & (j - 1);
    }
    return ret;
    */
}


unsigned int
fenwick_find(fenwick_t *self, long long sum)
{
    unsigned int j = 0;
    unsigned int k;
    long long s = sum;
    unsigned int half = self->log_max_index;
    while (half > 0) {
        /* Skip non-existent entries */
        while (j + half > self->max_index) {
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



