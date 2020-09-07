/*
** Copyright (C) 2015 University of Oxford
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

#ifndef __FENWICK_H__
#define __FENWICK_H__

#include <stdlib.h>
#include <inttypes.h>

typedef struct {
    size_t size;
    size_t log_size;
    double rebuild_threshold;
    /* Variables used for Kahan summation of the running total */
    double total_sum;
    double total_c;
    double *tree;
    double *values;
} fenwick_t;

void fenwick_print_state(fenwick_t *self, FILE *out);
void fenwick_verify(fenwick_t *self, double eps);
int fenwick_alloc(fenwick_t *, size_t);
int fenwick_expand(fenwick_t *, size_t);
int fenwick_free(fenwick_t *);
double fenwick_get_total(fenwick_t *);
void fenwick_rebuild(fenwick_t *);
bool fenwick_rebuild_required(fenwick_t *);
double fenwick_get_numerical_drift(fenwick_t *self);
void fenwick_increment(fenwick_t *, size_t, double);
void fenwick_set_value(fenwick_t *, size_t, double);
double fenwick_get_cumulative_sum(fenwick_t *, size_t);
double fenwick_get_value(fenwick_t *, size_t);
size_t fenwick_find(fenwick_t *, double);
size_t fenwick_get_size(fenwick_t *);

#endif /*__FENWICK_H__*/
