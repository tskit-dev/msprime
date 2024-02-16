#include <cstddef>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <type_traits>
#include <tskit.h>

static void
handle_tskit_return_code(int code)
{
    if (code != 0) {
        std::ostringstream o;
        o << tsk_strerror(code);
        throw std::runtime_error(o.str());
    }
}

struct edge_plus_time {
    double time;
    tsk_id_t parent, child;
    double left, right;
};

int
sort_edges(tsk_table_sorter_t *sorter, tsk_size_t start)
{
    if (sorter->tables->edges.metadata_length != 0) {
        throw std::invalid_argument(
            "the sorter does not currently handle edge metadata");
    }
    if (start != 0) {
        throw std::invalid_argument("the sorter requires start==0");
    }

    std::vector<edge_plus_time> temp;
    temp.reserve(static_cast<std::size_t>(sorter->tables->edges.num_rows));

    auto edges = &sorter->tables->edges;
    auto nodes = &sorter->tables->nodes;

    for (tsk_size_t i = 0; i < sorter->tables->edges.num_rows; ++i) {
        temp.push_back(edge_plus_time{ nodes->time[edges->parent[i]], edges->parent[i],
            edges->child[i], edges->left[i], edges->right[i] });
    }

    std::sort(begin(temp), end(temp),
        [](const edge_plus_time &lhs, const edge_plus_time &rhs) {
            if (lhs.time == rhs.time) {
                if (lhs.parent == rhs.parent) {
                    if (lhs.child == rhs.child) {
                        return lhs.left < rhs.left;
                    }
                    return lhs.child < rhs.child;
                }
                return lhs.parent < rhs.parent;
            }
            return lhs.time < rhs.time;
        });

    for (std::size_t i = 0; i < temp.size(); ++i) {
        edges->left[i] = temp[i].left;
        edges->right[i] = temp[i].right;
        edges->parent[i] = temp[i].parent;
        edges->child[i] = temp[i].child;
    }

    return 0;
}

int
main(int argc, char **argv)
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input.trees output.trees\n";
        std::exit(0);
    }
    const char *infile = argv[1];
    const char *outfile = argv[2];

    tsk_table_collection_t tables;
    auto ret = tsk_table_collection_load(&tables, infile, 0);
    handle_tskit_return_code(ret);

    tsk_table_sorter_t sorter;
    ret = tsk_table_sorter_init(&sorter, &tables, 0);
    handle_tskit_return_code(ret);
    sorter.sort_edges = sort_edges;
    try {
        ret = tsk_table_sorter_run(&sorter, NULL);
    } catch (std::exception &e) {
        std::cerr << e.what() << '\n';
        std::exit(1);
    }
    handle_tskit_return_code(ret);
    ret = tsk_table_collection_dump(&tables, outfile, 0);
    handle_tskit_return_code(ret);
    ret = tsk_table_collection_free(&tables);
    handle_tskit_return_code(ret);
}

