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
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "err.h"
#include "msprime.h"

void
recomb_map_print_state(recomb_map_t *self)
{
    size_t j;

    printf("recombination_map::\n");
    printf("size = %d\n", (int) self->size);
    for (j = 0; j < self->size; j++) {
        printf("\t%d\t%f\n", self->coordinates[j], self->rates[j]);
    }
}

int WARN_UNUSED
recomb_map_alloc(recomb_map_t *self, uint32_t *coordinates, double *rates,
        size_t size)
{
    int ret = MSP_ERR_BAD_RECOMBINATION_MAP;
    size_t j;

    memset(self, 0, sizeof(recomb_map_t));
    if (size < 2) {
        goto out;
    }
    self->coordinates = malloc(size * sizeof(uint32_t));
    self->rates = malloc(size * sizeof(double));
    if (self->coordinates == NULL || self->rates == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->size = size;
    if (coordinates[0] != 0) {
        goto out;
    }
    for (j = 0; j < size; j++) {
        if (rates[j] < 0) {
            goto out;
        }
        if (j > 0) {
            /* Coordinates must be sorted */
            if (coordinates[j] <= coordinates[j - 1]) {
                goto out;
            }
        }
        self->rates[j] = rates[j];
        self->coordinates[j] = coordinates[j];
    }
    ret = 0;
out:
    return ret;
}

int
recomb_map_free(recomb_map_t *self)
{
    if (self->coordinates != NULL) {
        free(self->coordinates);
    }
    if (self->rates != NULL) {
        free(self->rates);
    }
    return 0;
}

double
recomb_map_get_effective_rate(recomb_map_t *self)
{
    double effective_rate = 0.0;
    double length;
    size_t j;

    for (j = 0; j < self->size - 1; j++) {
        length = self->coordinates[j + 1] - self->coordinates[j];
        effective_rate += self->rates[j] * length;
    }
    effective_rate /= self->coordinates[self->size - 1];
    return effective_rate;
}
