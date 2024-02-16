#include <stdio.h>
#include <stdlib.h>
#include <tskit.h>

int
main(int argc, char **argv)
{
    int ret;
    tsk_treeseq_t ts;

    if (argc != 2) {
        fprintf(stderr, "usage: <tree sequence file>");
        exit(EXIT_FAILURE);
    }
    ret = tsk_treeseq_load(&ts, argv[1], 0);
    if (ret < 0) {
        /* Error condition. Free and exit */
        tsk_treeseq_free(&ts);
        fprintf(stderr, "%s", tsk_strerror(ret));
        exit(EXIT_FAILURE);
    }
    printf("Loaded tree sequence with %lld nodes and %lld edges from %s\n",
        (long long) tsk_treeseq_get_num_nodes(&ts),
        (long long) tsk_treeseq_get_num_edges(&ts),
        argv[1]);
    tsk_treeseq_free(&ts);

    return EXIT_SUCCESS;
}
