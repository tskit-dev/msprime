/*
** Copyright (C) 2014 Jerome Kelleher <jerome.kelleher@ed.ac.uk>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
** 
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software 
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/
/*
 * Binary index tree (also known as a Fenwick tree) implementation.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <assert.h>

#include "util.h"
#include "bit.h"



void 
bit_alloc(bit_t *self)
{
    unsigned int u = self->max_index;
    while (u != 0) {
        self->log_max_index = u;
        u -= (u & -u);
    }
    self->tree = xcalloc((1 + self->max_index), sizeof(long long));
    self->values = xcalloc((1 + self->max_index), sizeof(long long));
}

void bit_free(bit_t *self)
{
    free(self->tree);
    free(self->values);
}

long long 
bit_get_total(bit_t *self)
{
    long long ret = bit_get_cumulative_sum(self, self->max_index);
    return ret;
}

void 
bit_increment(bit_t *self, unsigned int index, long long value) 
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
bit_set_value(bit_t *self, unsigned int index, long long value) 
{
    long long v = value - self->values[index];
    bit_increment(self, index, v);
}

long long bit_get_cumulative_sum(bit_t *self, unsigned int index)
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

long long bit_get_value(bit_t *self, unsigned int index)
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
bit_find(bit_t *self, long long sum)
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



