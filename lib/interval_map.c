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
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "util.h"
#include "msprime.h"

void
interval_map_print_state(interval_map_t *self, FILE *out)
{
    size_t j;

    fprintf(out, "interval_map (%p):: size = %d\n", (void *) self, (int) self->size);
    fprintf(out, "\tsequence_length = %f\n", interval_map_get_sequence_length(self));
    fprintf(out, "\tindex\tposition\tvalue\n");
    for (j = 0; j < self->size; j++) {
        fprintf(out, "\t%d\t%f\t%f\n", (int) j, self->position[j], self->value[j]);
    }
}

int MSP_WARN_UNUSED
interval_map_alloc(interval_map_t *self, size_t size, double *position, double *value)
{
    int ret = 0;
    size_t j;

    memset(self, 0, sizeof(interval_map_t));
    if (size < 2) {
        ret = MSP_ERR_INSUFFICIENT_INTERVALS;
        goto out;
    }
    /* Check the framing position */
    if (position[0] != 0.0) {
        ret = MSP_ERR_INTERVAL_MAP_START_NON_ZERO;
        goto out;
    }
    self->position = malloc(size * sizeof(*self->position));
    self->value = malloc(size * sizeof(*self->value));
    if (self->position == NULL || self->value == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->size = size;
    for (j = 0; j < size; j++) {
        if (position[j] < 0) {
            ret = MSP_ERR_NEGATIVE_INTERVAL_POSITION;
            goto out;
        }
        if (j > 0) {
            if (position[j] <= position[j - 1]) {
                ret = MSP_ERR_INTERVAL_POSITIONS_UNSORTED;
                goto out;
            }
        }
        self->value[j] = value[j];
        self->position[j] = position[j];
    }
out:
    return ret;
}

int MSP_WARN_UNUSED
interval_map_alloc_single(interval_map_t *self, double sequence_length, double value)
{
    double position[2] = {0, sequence_length};
    return interval_map_alloc(self, 2, position, &value);
}

int
interval_map_free(interval_map_t *self)
{
    msp_safe_free(self->position);
    msp_safe_free(self->value);
    return 0;
}

double
interval_map_get_sequence_length(interval_map_t *self)
{
    return self->position[self->size - 1];
}

size_t
interval_map_get_size(interval_map_t *self)
{
    return self->size;
}

size_t
interval_map_get_num_intervals(interval_map_t *self)
{
    return self->size - 1;
}

size_t
interval_map_get_index(interval_map_t *self, double x)
{
    size_t index = tsk_search_sorted(self->position, self->size, x);

    if (self->position[index] > x) {
        index--;
    }
    return index;
}
