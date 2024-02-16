#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <err.h>

#include <tskit/tables.h>

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        errx(EXIT_FAILURE, "line %d: %s", __LINE__, tsk_strerror(val));                 \
    }

void
simulate(
    tsk_table_collection_t *tables, int N, int T, int simplify_interval)
{
    tsk_id_t *buffer, *parents, *children, child, left_parent, right_parent;
    double breakpoint;
    int ret, j, t, b;

    assert(simplify_interval != 0); // leads to division by zero
    buffer = malloc(2 * N * sizeof(tsk_id_t));
    if (buffer == NULL) {
        errx(EXIT_FAILURE, "Out of memory");
    }
    tables->sequence_length = 1.0;
    parents = buffer;
    for (j = 0; j < N; j++) {
        parents[j]
            = tsk_node_table_add_row(&tables->nodes, 0, T, TSK_NULL, TSK_NULL, NULL, 0);
        check_tsk_error(parents[j]);
    }
    b = 0;
    for (t = T - 1; t >= 0; t--) {
        /* Alternate between using the first and last N values in the buffer */
        parents = buffer + (b * N);
        b = (b + 1) % 2;
        children = buffer + (b * N);
        for (j = 0; j < N; j++) {
            child = tsk_node_table_add_row(
                &tables->nodes, 0, t, TSK_NULL, TSK_NULL, NULL, 0);
            check_tsk_error(child);
            /* NOTE: the use of rand() is discouraged for
             * research code and proper random number generator
             * libraries should be preferred.
             */
            left_parent = parents[(size_t)((rand()/(1.+RAND_MAX))*N)];
            right_parent = parents[(size_t)((rand()/(1.+RAND_MAX))*N)];
            do {
                breakpoint = rand()/(1.+RAND_MAX);
            } while (breakpoint == 0); /* tiny proba of breakpoint being 0 */
            ret = tsk_edge_table_add_row(
                &tables->edges, 0, breakpoint, left_parent, child, NULL, 0);
            check_tsk_error(ret);
            ret = tsk_edge_table_add_row(
                &tables->edges, breakpoint, 1, right_parent, child, NULL, 0);
            check_tsk_error(ret);
            children[j] = child;
        }
        if (t % simplify_interval == 0) {
            printf("Simplify at generation %lld: (%lld nodes %lld edges)",
                (long long) t,
                (long long) tables->nodes.num_rows,
                (long long) tables->edges.num_rows);
            /* Note: Edges must be sorted for simplify to work, and we use a brute force
             * approach of sorting each time here for simplicity. This is inefficient. */
            ret = tsk_table_collection_sort(tables, NULL, 0);
            check_tsk_error(ret);
            ret = tsk_table_collection_simplify(tables, children, N, 0, NULL);
            check_tsk_error(ret);
            printf(" -> (%lld nodes %lld edges)\n",
                (long long) tables->nodes.num_rows,
                (long long) tables->edges.num_rows);
            for (j = 0; j < N; j++) {
                children[j] = j;
            }
        }
    }
    free(buffer);
}

int
main(int argc, char **argv)
{
    int ret;
    tsk_table_collection_t tables;

    if (argc != 6) {
        errx(EXIT_FAILURE, "usage: N T simplify-interval output-file seed");
    }
    ret = tsk_table_collection_init(&tables, 0);
    check_tsk_error(ret);
    srand((unsigned)atoi(argv[5]));
    simulate(&tables, atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
    ret = tsk_table_collection_dump(&tables, argv[4], 0);
    check_tsk_error(ret);

    tsk_table_collection_free(&tables);
    return 0;
}
