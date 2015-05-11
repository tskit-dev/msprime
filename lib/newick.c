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
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "err.h"
#include "msprime.h"

int
newick_converter_alloc(newick_converter_t *self, tree_sequence_t *tree_sequence, size_t precision,
        int all_breakpoints)
{
    int ret = -1;
    int flags = 0;

    self->sample_size = tree_sequence->sample_size;
    self->num_loci = tree_sequence->num_loci;
    memset(&self->diff_iterator, 0, sizeof(tree_diff_iterator_t));
    if (all_breakpoints) {
        flags = MSP_ALL_BREAKPOINTS;
    }
    ret = tree_diff_iterator_alloc(&self->diff_iterator, tree_sequence, flags);
    if (ret != 0) {
        goto out;
    }

out:
    return ret;
}

int
newick_converter_free(newick_converter_t *self)
{
    tree_diff_iterator_free(&self->diff_iterator);
    return 0;
}
