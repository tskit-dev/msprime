#include "tables.h"

#include <iostream>
#include <stdlib.h>

using namespace std;


static void
handle_error(string msg, int err)
{
    cout << "Error:" << msg << ":" << msp_strerror(err) << endl;
    exit(1);
}

int
main(int argc, char **argv)
{
    int j, ret;
    table_collection_t tables;
    node_id_t samples[] = {0};

    /* Allocate a table collection, with all the internal tables initialised
     * using default alloc sizes */
    ret = table_collection_alloc(&tables, MSP_ALLOC_TABLES);
    if (ret != 0) {
        handle_error("table_collection_alloc", ret);
    }
    /* NB: must set the sequence_length !! */
    tables.sequence_length = 1.0;

    /* Create a simple chain of nodes, with 0 as the only sample. */
    for (j = 0; j < 10; j++) {
        ret = node_table_add_row(&tables.nodes, j == 0, j, 0, NULL, 0);
        if (ret < 0) {
            handle_error("add_node", ret);
        }
        if (j > 0) {
            ret = edge_table_add_row(&tables.edges, 0, 1, j, j - 1);
            if (ret < 0) {
                handle_error("add_edge", ret);
            }
        }
    }

    /* Write the state out to file */
    ret = table_collection_dump(&tables, "tmp.hdf5", 0);
    if (ret != 0) {
        handle_error("dump", ret);
    }

    /* Useful debugging feature */
    table_collection_print_state(&tables, stdout);

    ret = table_collection_simplify(&tables, samples, 1, 0, NULL);
    if (ret != 0) {
        handle_error("simplify", ret);
    }
    /* After simplify, we only have 1 node left and no edges */
    node_table_print_state(&tables.nodes, stdout);
    edge_table_print_state(&tables.edges, stdout);

    /* Clean up. This should usually also be done in the error handling case,
     * but since this is a simple standalone program and we can exit on
     * error. */
    table_collection_free(&tables);
    return 0;
}
