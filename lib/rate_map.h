/*
** Copyright (C) 2020 University of Oxford
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

#ifndef __RATE_MAP_H__
#define __RATE_MAP_H__

#include <stddef.h>
#include <stdio.h>

#include "util.h"

typedef struct {
    size_t size;
    double *position;
    double *rate;
    double *cumulative_mass;
    fast_search_t position_lookup;
} rate_map_t;

int rate_map_alloc(rate_map_t *self, size_t size, double *position, double *value);
int rate_map_alloc_single(rate_map_t *self, double sequence_length, double value);
int rate_map_copy(rate_map_t *to, rate_map_t *from);
int rate_map_free(rate_map_t *self);
void rate_map_print_state(rate_map_t *self, FILE *out);
double rate_map_get_sequence_length(rate_map_t *self);
size_t rate_map_get_size(rate_map_t *self);
size_t rate_map_get_num_intervals(rate_map_t *self);
size_t rate_map_get_index(rate_map_t *self, double x);
double rate_map_get_total_mass(rate_map_t *self);
double rate_map_mass_between(rate_map_t *self, double left, double right);
double rate_map_mass_to_position(rate_map_t *self, double mass);
double rate_map_position_to_mass(rate_map_t *self, double position);
double rate_map_shift_by_mass(rate_map_t *self, double pos, double mass);

#endif /*__RATE_MAP_H__*/
