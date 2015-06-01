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

/*
 * Binary index tree (also known as a Fenwick tree) implementation.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "fenwick.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

static void
test_simple_cases(void)
{
    fenwick_t t;
    int64_t s;
    size_t j, n;
    for (n = 1; n < 100; n++) {
        s = 0;
        assert(fenwick_alloc(&t, n) == 0);
        for (j = 1; j <= n; j++) {
            fenwick_increment(&t, j, (int64_t) j);
            s = s + (int64_t) j;
            assert(fenwick_get_value(&t, j) == (int64_t) j);
            assert(fenwick_get_cumulative_sum(&t, j) == s);
            assert(fenwick_get_total(&t) == s);
            assert(fenwick_find(&t, s) == j);
            fenwick_set_value(&t, j, 0);
            assert(fenwick_get_value(&t, j) == 0);
            assert(fenwick_get_cumulative_sum(&t, j) == s - (int64_t) j);
            fenwick_set_value(&t, j, (int64_t) j);
            assert(fenwick_get_value(&t, j) == (int64_t) j);

        }
        assert(fenwick_free(&t) == 0);
    }
}


int
main(int argc, char** argv)
{
    test_simple_cases();
    /* TODO add some proper test cases */
    return EXIT_SUCCESS;

}
