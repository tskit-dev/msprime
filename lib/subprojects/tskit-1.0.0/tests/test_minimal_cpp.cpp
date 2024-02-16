/* * MIT License
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

/* Minimal tests to make sure that tskit at least compiles and links
 * in a simple C++ program */

#include <iostream>
#include <cassert>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstring>

#include <tskit.h>

using namespace std;

void
test_kas_strerror()
{
    std::cout << "test_kas_strerror" << endl;
    std::ostringstream o;
    o << kas_strerror(KAS_ERR_NO_MEMORY);
    assert(std::string("Out of memory").compare(o.str()) == 0);
}

void
test_strerror()
{
    std::cout << "test_strerror" << endl;
    std::ostringstream o;
    o << tsk_strerror(TSK_ERR_NO_MEMORY);
    assert(std::string("Out of memory. (TSK_ERR_NO_MEMORY)").compare(o.str()) == 0);
}

void
test_load_error()
{
    std::cout << "test_open_error" << endl;
    tsk_treeseq_t ts;
    int ret = tsk_treeseq_load(&ts, "no such file", 0);
    assert(ret == TSK_ERR_IO);
    tsk_treeseq_free(&ts);
}

void
test_table_basics()
{
    std::cout << "test_table_basics" << endl;
    tsk_table_collection_t tables;
    int ret = tsk_table_collection_init(&tables, 0);
    assert(ret == 0);

    ret = tsk_node_table_add_row(&tables.nodes, 0, 1.0, TSK_NULL, TSK_NULL, NULL, 0);
    assert(ret == 0);
    ret = tsk_node_table_add_row(&tables.nodes, 0, 2.0, TSK_NULL, TSK_NULL, NULL, 0);
    assert(ret == 1);
    assert(tables.nodes.num_rows == 2);

    tsk_table_collection_free(&tables);
}

/* A definition of sort_edges that uses C++ std::sort and inlining of the
 * comparison function to achieve significantly better performance than
 * the builtin method in tskit.
 */
int
cpp_sort_edges(tsk_table_sorter_t *sorter, tsk_size_t start)
{
    struct _edge {
        double left, right;
        tsk_id_t parent, child;

        _edge(double l, double r, tsk_id_t p, tsk_id_t c)
            : left{ l }, right{ r }, parent{ p }, child{ c }
        {
        }
    };
    tsk_edge_table_t *edges = &sorter->tables->edges;
    const double *node_time = sorter->tables->nodes.time;
    std::vector<_edge> sorted_edges;
    size_t num_edges = edges->num_rows;
    size_t j;

    /* This is the comparison function.  We cannot define an
     * operator < for _edge because we need to bind the node times
     * so we have to use a functional method. This is a copy of the cmp
     * from fwdpp.  Only difference is the final time comparison
     * (fwdpp table times go forwards). */
    const auto cmp = [&node_time](const _edge &lhs, const _edge &rhs) {
        auto tl = node_time[lhs.parent];
        auto tr = node_time[rhs.parent];
        if (tl == tr) {
            if (lhs.parent == rhs.parent) {
                if (lhs.child == rhs.child) {
                    return lhs.left < rhs.left;
                }
                return lhs.child < rhs.child;
            }
            return lhs.parent < rhs.parent;
        }
        return tl < tr;
    };

    assert(start == 0);
    /* Let's not bother with metadata */
    assert(edges->metadata_length == 0);

    sorted_edges.reserve(num_edges);
    for (j = 0; j < num_edges; j++) {
        sorted_edges.emplace_back(
            edges->left[j], edges->right[j], edges->parent[j], edges->child[j]);
    }

    std::sort(begin(sorted_edges), end(sorted_edges), cmp);

    for (j = 0; j < num_edges; j++) {
        edges->left[j] = sorted_edges[j].left;
        edges->right[j] = sorted_edges[j].right;
        edges->parent[j] = sorted_edges[j].parent;
        edges->child[j] = sorted_edges[j].child;
    }
    return 0;
}

void
test_edge_sorting()
{
    std::cout << "test_edge_sorting" << endl;
    tsk_table_collection_t tables;
    tsk_id_t n = 10;
    tsk_id_t j;
    int ret = tsk_table_collection_init(&tables, 0);
    assert(ret == 0);

    tables.sequence_length = 1.0;
    /* Make a stick tree */
    /* Add nodes and edges */
    for (j = 0; j < n; j++) {
        ret = tsk_node_table_add_row(
            &tables.nodes, TSK_NODE_IS_SAMPLE, j + 1, TSK_NULL, TSK_NULL, NULL, 0);
        assert(ret == j);
    }
    for (j = n - 1; j > 0; j--) {
        tsk_edge_table_add_row(&tables.edges, 0, 1, j, j - 1, NULL, 0);
    }
    assert(tables.nodes.num_rows == (tsk_size_t) n);
    assert(tables.edges.num_rows == (tsk_size_t) n - 1);

    /* Make sure the edges are unsorted */
    /* Not calling TSK_CHECK_TREES so casting is safe */
    ret = (int) tsk_table_collection_check_integrity(&tables, TSK_CHECK_EDGE_ORDERING);
    assert(ret == TSK_ERR_EDGES_NOT_SORTED_PARENT_TIME);

    /* Sort the tables */
    tsk_table_sorter_t sorter;
    ret = tsk_table_sorter_init(&sorter, &tables, 0);
    assert(ret == 0);
    /* Set the sort_edges to our local C++ version. We could also set some
     * persistent state in sorter.params if we wanted to. */
    sorter.sort_edges = cpp_sort_edges;
    ret = tsk_table_sorter_run(&sorter, NULL);
    assert(ret == 0);
    tsk_table_sorter_free(&sorter);

    /* Make sure the edges are now sorted */
    ret = (int) tsk_table_collection_check_integrity(&tables, TSK_CHECK_EDGE_ORDERING);
    assert(ret == 0);

    tsk_table_collection_free(&tables);
}

int
sort_edges_raises_exception(tsk_table_sorter_t *sorter, tsk_size_t start)
{
    throw std::exception();
    return 0;
}

int
sort_edges_raises_non_exception(tsk_table_sorter_t *sorter, tsk_size_t start)
{
    throw 42;
    return 0;
}

int
safe_sort_edges(tsk_table_sorter_t *sorter, tsk_size_t start)
{
    int ret = 0;
    if (sorter->user_data == NULL) {
        try {
            ret = sort_edges_raises_exception(sorter, start);
        } catch (...) {
            ret = -12345;
        }
    } else {
        try {
            ret = sort_edges_raises_non_exception(sorter, start);
        } catch (...) {
            ret = -12346;
        }
    }
    return ret;
}

void
test_edge_sorting_errors()
{
    /* Some inexplicable error happened here on 32 bit Windows where the
     * exceptions were not being caught as expected. This seems much
     * more likely to be a platform quirk that a real bug in our code,
     * so just disabling the test there.
     *
     * https://github.com/tskit-dev/tskit/issues/1790
     * https://github.com/tskit-dev/tskit/pull/1791
     */
#if !defined(_WIN32)
    std::cout << "test_edge_sorting_errors" << endl;
    tsk_table_collection_t tables;
    tsk_table_sorter_t sorter;
    tsk_id_t ret = tsk_table_collection_init(&tables, 0);

    assert(ret == 0);
    tables.sequence_length = 1.0;

    ret = tsk_table_sorter_init(&sorter, &tables, 0);
    assert(ret == 0);
    sorter.sort_edges = safe_sort_edges;
    ret = tsk_table_sorter_run(&sorter, NULL);
    assert(ret == -12345);

    /* Use the user_data as a way to communicate with the sorter
     * function. Here, we want to try out two different types
     * of exception that get thrown. */
    sorter.user_data = &tables;
    ret = tsk_table_sorter_run(&sorter, NULL);
    assert(ret == -12346);

    tsk_table_sorter_free(&sorter);
    tsk_table_collection_free(&tables);
#endif
}

int
main()
{
    test_kas_strerror();
    test_strerror();
    test_load_error();
    test_table_basics();
    test_edge_sorting();
    test_edge_sorting_errors();
    return 0;
}
