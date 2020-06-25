/*
** Copyright (C) 2015-2020 University of Oxford
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
#include <string.h>

#include <gsl/gsl_math.h>

#include "util.h"
#include "fenwick.h"

/* Return the value for the specified index as represented within
 * the tree. */
static double
fenwick_compute_tree_value(fenwick_t *self, size_t index)
{
    double ret = fenwick_get_cumulative_sum(self, index);
    if (index > 1) {
        ret -= fenwick_get_cumulative_sum(self, index - 1);
    }
    return ret;
}

void
fenwick_verify(fenwick_t *self, double eps)
{
    size_t j;
    double computed_value;

    for (j = 1; j <= self->size; j++) {
        computed_value = fenwick_compute_tree_value(self, j);
        assert(gsl_fcmp(computed_value, self->values[j], eps) == 0);
    }
}

void
fenwick_print_state(fenwick_t *self, FILE *out)
{
    size_t j;

    fprintf(out, "Fenwick tree @%p\n", (void *) self);

    for (j = 1; j <= self->size; j++) {
        fprintf(out, "%d\t%.16g\t%.16g\t%.16g\n", (int) j, self->values[j],
            self->tree[j], fabs(self->values[j] - fenwick_compute_tree_value(self, j)));
    }
}

static void
fenwick_set_log_size(fenwick_t *self)
{
    size_t u = self->size;

    while (u != 0) {
        self->log_size = u;
        u -= (u & -u);
    }
}

int MSP_WARN_UNUSED
fenwick_alloc(fenwick_t *self, size_t initial_size)
{
    int ret = 0;

    memset(self, 0, sizeof(*self));
    self->size = initial_size;
    self->tree = calloc((1 + self->size), sizeof(*self->tree));
    self->values = calloc((1 + self->size), sizeof(*self->tree));
    if (self->tree == NULL || self->values == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    fenwick_set_log_size(self);
out:
    return ret;
}

int MSP_WARN_UNUSED
fenwick_expand(fenwick_t *self, size_t increment)
{
    int ret = MSP_ERR_NO_MEMORY;
    size_t j, n, k;
    void *p;

    p = realloc(self->tree, (1 + self->size + increment) * sizeof(*self->tree));
    if (p == NULL) {
        goto out;
    }
    self->tree = p;
    p = realloc(self->values, (1 + self->size + increment) * sizeof(*self->tree));
    if (p == NULL) {
        goto out;
    }
    self->values = p;

    self->size += increment;
    fenwick_set_log_size(self);
    for (j = self->size - increment + 1; j <= self->size; j++) {
        self->values[j] = 0;
        self->tree[j] = 0;
        n = j;
        k = 1;
        while (n % 2 == 0) {
            self->tree[j] += self->tree[j - k];
            k *= 2;
            n >>= 1;
        }
    }
    ret = 0;
out:
    return ret;
}

int
fenwick_free(fenwick_t *self)
{
    msp_safe_free(self->tree);
    msp_safe_free(self->values);
    return 0;
}

size_t
fenwick_get_size(fenwick_t *self)
{
    return self->size;
}

double
fenwick_get_total(fenwick_t *self)
{
    return fenwick_get_cumulative_sum(self, self->size);
}

void
fenwick_increment(fenwick_t *self, size_t index, double value)
{
    size_t j;
    const size_t size = self->size;
    double *restrict tree = self->tree;

    assert(0 < index && index <= size);
    self->values[index] += value;
    for (j = index; j <= size; j += (j & -j)) {
        tree[j] += value;
    }
}

void
fenwick_set_value(fenwick_t *self, size_t index, double value)
{
    /* propagate loss of the old value up the tree */
    fenwick_increment(self, index, -self->values[index]);
    /* Ensure the new value is exactly 0 */
    self->values[index] = 0;
    /* propagate gain of the new value up the tree */
    fenwick_increment(self, index, value);
}

double
fenwick_get_cumulative_sum(fenwick_t *self, size_t index)
{
    double ret = 0;
    const double *restrict tree = self->tree;
    size_t j;

    assert(0 < index && index <= self->size);
    for (j = index; j > 0; j -= (j & -j)) {
        ret += tree[j];
    }
    return ret;
}

double
fenwick_get_value(fenwick_t *self, size_t index)
{
    assert(0 < index && index <= self->size);
    return self->values[index];
}

size_t
fenwick_find(fenwick_t *self, double sum)
{
    size_t j = 0;
    size_t k, index;
    double s = sum;
    const double *restrict tree = self->tree;
    const double *restrict values = self->values;
    const size_t size = self->size;
    size_t half = self->log_size;

    while (half > 0) {
        /* Skip non-existent entries */
        while (j + half > size) {
            half >>= 1;
        }
        k = j + half;
        if (s > tree[k]) {
            j = k;
            s -= tree[j];
        }
        half >>= 1;
    }
    index = j + 1;
    /* We can have situations due to numerical imprecision where
     * the sum points to an index that's actually mapped to zero,
     * so skip ahead until we find a non-zero value. */
    while (index <= size && values[index] == 0) {
        index++;
    }
    return index;
}
