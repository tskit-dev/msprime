#include <stdio.h>
#include <stdlib.h>
#include <tskit/tables.h>

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        fprintf(stderr, "line %d: %s", __LINE__, tsk_strerror(val));                    \
        exit(EXIT_FAILURE);                                                             \
    }

int
main(int argc, char **argv)
{
    int j, ret;
    tsk_edge_table_t edges;

    ret = tsk_edge_table_init(&edges, 0);
    check_tsk_error(ret);
    for (j = 0; j < 5; j++) {
        ret = tsk_edge_table_add_row(&edges, 0, 1, j + 1, j, NULL, 0);
        check_tsk_error(ret);
    }
    tsk_edge_table_print_state(&edges, stdout);
    tsk_edge_table_free(&edges);

    return EXIT_SUCCESS;
}
