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
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "err.h"
#include "msprime.h"

int
newick_alloc(newick_t *self, const char *tree_file_name)
{
    int ret = -1;
    size_t N;

    memset(self, 0, sizeof(newick_t));
    ret = tree_file_open(&self->tree_file, tree_file_name, 'r');
    if (ret != 0) {
        goto out;
    }
    if (!tree_file_issorted(&self->tree_file)) {
        ret = MSP_ERR_TREE_FILE_NOT_SORTED;
        goto out;
    }
    self->sample_size = self->tree_file.sample_size;
    N  = 2 * self->sample_size;
    printf("N = %d\n", (int) N);
out:
    return ret;
}

int
newick_free(newick_t *self)
{
    tree_file_close(&self->tree_file);
    return 0;
}

int
newick_next(newick_t *self)
{
    int ret = -1;
    return ret;
}
