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

#ifndef __FENWICK_H__
#define __FENWICK_H__

typedef struct {
    unsigned int max_index;
    unsigned int log_max_index;
    long long *tree;
    long long *values;
} fenwick_t;


void fenwick_alloc(fenwick_t *);
void fenwick_free(fenwick_t *);
long long fenwick_get_total(fenwick_t *);
void fenwick_increment(fenwick_t *, unsigned int, long long);
void fenwick_set_value(fenwick_t *, unsigned int, long long);
long long fenwick_get_cumulative_sum(fenwick_t *, unsigned int);
long long fenwick_get_value(fenwick_t *, unsigned int);
unsigned int fenwick_find(fenwick_t *, long long);

#endif /*__FENWICK_H__*/
