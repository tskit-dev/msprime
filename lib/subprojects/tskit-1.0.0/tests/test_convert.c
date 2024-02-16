/*
 * MIT License
 *
 * Copyright (c) 2019-2022 Tskit Developers
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "testlib.h"
#include <tskit/convert.h>

#include <unistd.h>
#include <stdlib.h>

static void
test_single_tree_newick(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t t;
    size_t buffer_size = 1024;
    char newick[buffer_size];

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0)
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK)

    ret = tsk_convert_newick(&t, 0, 0, TSK_NEWICK_LEGACY_MS_LABELS, buffer_size, newick);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Seems odd, but this is what a single node newick tree looks like.
     * Newick parsers seems to accept it in any case */
    CU_ASSERT_STRING_EQUAL(newick, "1;");

    ret = tsk_convert_newick(&t, 0, 0, 0, buffer_size, newick);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(newick, "n0;");

    ret = tsk_convert_newick(&t, 4, 0, TSK_NEWICK_LEGACY_MS_LABELS, buffer_size, newick);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(newick, "(1:1,2:1);");
    ret = tsk_convert_newick(&t, 4, 0, 0, buffer_size, newick);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(newick, "(n0:1,n1:1);");

    ret = tsk_convert_newick(&t, 6, 0, TSK_NEWICK_LEGACY_MS_LABELS, buffer_size, newick);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(newick, "((1:1,2:1):2,(3:2,4:2):1);");

    ret = tsk_convert_newick(&t, 6, 0, 0, buffer_size, newick);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(newick, "((n0:1,n1:1):2,(n2:2,n3:2):1);");

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_newick_errors(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t t;
    size_t j, len;
    size_t buffer_size = 1024;
    char newick[buffer_size];

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);

    ret = tsk_tree_init(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0)
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_TREE_OK)

    ret = tsk_convert_newick(&t, -1, 1, 0, buffer_size, newick);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    ret = tsk_convert_newick(&t, 7, 1, 0, buffer_size, newick);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    ret = tsk_convert_newick(&t, 6, 0, 0, buffer_size, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_convert_newick(&t, 6, 0, 0, buffer_size, newick);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    len = 1 + strlen(newick);
    for (j = 0; j < len; j++) {
        ret = tsk_convert_newick(&t, 6, 0, 0, j, newick);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BUFFER_OVERFLOW);
    }
    ret = tsk_convert_newick(&t, 6, 0, TSK_NEWICK_LEGACY_MS_LABELS, len, newick);

    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(newick, "((1:1,2:1):2,(3:2,4:2):1);");

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_single_tree_newick", test_single_tree_newick },
        { "test_single_tree_newick_errors", test_single_tree_newick_errors },
        { NULL, NULL },
    };
    return test_main(tests, argc, argv);
}
