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

#include <gsl/gsl_math.h>

#include "err.h"
#include "msprime.h"

void
recomb_map_print_state(recomb_map_t *self)
{
    size_t j;

    printf("recombination_map:: size = %d\n", (int) self->size);
    printf("\tsequence_length = %f\n", recomb_map_get_sequence_length(self));
    printf("\ttotal_rate = %f\n", recomb_map_get_total_recombination_rate(self));
    printf("\tindex\tlocation\trate\n");
    for (j = 0; j < self->size; j++) {
        printf("\t%d\t%f\t%f\n", (int) j, self->positions[j], self->rates[j]);
    }
}

int WARN_UNUSED
recomb_map_alloc(recomb_map_t *self, double *positions, double *rates,
        size_t size)
{
    int ret = MSP_ERR_BAD_RECOMBINATION_MAP;
    double length;
    size_t j;

    memset(self, 0, sizeof(recomb_map_t));
    if (size < 2) {
        goto out;
    }
    self->positions = malloc(size * sizeof(double));
    self->rates = malloc(size * sizeof(double));
    if (self->positions == NULL || self->rates == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->total_mass = 0.0;
    self->size = size;
    if (positions[0] != 0.0) {
        goto out;
    }
    for (j = 0; j < size; j++) {
        if (rates[j] < 0 || positions[j] < 0) {
            goto out;
        }
        if (j > 0) {
            /* Coordinates must be sorted */
            if (positions[j] <= positions[j - 1]) {
                goto out;
            }
            length = positions[j] - positions[j - 1];
            self->total_mass += length * rates[j - 1];
        }
        self->rates[j] = rates[j];
        self->positions[j] = positions[j];
    }
    if (self->total_mass == 0.0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

int
recomb_map_free(recomb_map_t *self)
{
    if (self->positions != NULL) {
        free(self->positions);
    }
    if (self->rates != NULL) {
        free(self->rates);
    }
    return 0;
}

/* Returns the total recombination rate along the entire sequence.
 */
double
recomb_map_get_total_recombination_rate(recomb_map_t *self)
{
    return self->total_mass;
}

double
recomb_map_get_sequence_length(recomb_map_t *self)
{
    return self->positions[self->size - 1];
}

/* Remaps the specified genetic coordinate in the range {0, 1} to
 * the physical coordinate space in the range (0, sequence_length)
 */
double
recomb_map_phys_to_genetic(recomb_map_t *self, double x)
{
    size_t j;
    double s = 0.0;
    double ret = 0.0;
    double rate = 1.0;
    double last_phys_x, phys_x;

    if (self->total_mass > 0) {
        last_phys_x = 0;
        for (j = 1; j < self->size && x > self->positions[j]; j++) {
            phys_x = self->positions[j];
            rate = self->rates[j - 1];
            s += (phys_x - last_phys_x) * rate;
            last_phys_x = phys_x;
        }
        rate = self->rates[j - 1];
        s += (x - last_phys_x) * rate;
        assert(s >= 0 && s <= self->total_mass);
        ret = s / self->total_mass;
    }
    return ret;
}

double
recomb_map_genetic_to_phys(recomb_map_t *self, double x)
{
    size_t j;
    double s = 0.0;
    double ret = 0.0;
    double rate = 1.0;
    double last_phys_x, phys_x;
    /* the coordinate x is provided as a fraction from 0 to 1 so we
     * rescale into the rate [0, total_mass). */
    double genetic_x = x * self->total_mass;

    assert(x >= 0 && x <= 1.0);
    if (self->total_mass == 0.0) {
        if (x != 0.0) {
            ret = GSL_NAN;
        }
    } else {
        last_phys_x = 0;
        for (j = 1; j < self->size && s < genetic_x; j++) {
            phys_x = self->positions[j];
            rate = self->rates[j - 1];
            s += (phys_x - last_phys_x) * rate;
            last_phys_x = phys_x;
        }
        ret = last_phys_x - (s - genetic_x) / rate;
        assert(ret >= 0 && ret <= self->positions[self->size - 1]);
    }
    return ret;
}
