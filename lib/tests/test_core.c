/*
** Copyright (C) 2016-2020 University of Oxford
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

#include "testlib.h"

static void
test_strerror(void)
{
    int j;
    const char *msg;
    int max_error_code = 1024; /* totally arbitrary */

    for (j = 1; j > -max_error_code; j--) {
        msg = msp_strerror(j);
        CU_ASSERT_FATAL(msg != NULL);
        CU_ASSERT(strlen(msg) > 0);
    }
}

static void
test_strerror_tskit(void)
{
    int tskit_errors[]
        = { TSK_ERR_NO_MEMORY, TSK_ERR_NODE_OUT_OF_BOUNDS, TSK_ERR_EDGE_OUT_OF_BOUNDS };
    size_t j;
    int err;

    for (j = 0; j < sizeof(tskit_errors) / sizeof(*tskit_errors); j++) {
        err = msp_set_tsk_error(tskit_errors[j]);
        CU_ASSERT_TRUE(msp_is_tsk_error(err));
        CU_ASSERT_STRING_EQUAL(msp_strerror(err), tsk_strerror(tskit_errors[j]));
    }
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_strerror", test_strerror },
        { "test_strerror_tskit", test_strerror_tskit },
        CU_TEST_INFO_NULL,
    };

    return test_main(tests, argc, argv);
}
