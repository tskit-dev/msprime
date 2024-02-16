#include <stdio.h>
#include <stdlib.h>
#include <tskit/tables.h>

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        fprintf(stderr, "Error: line %d: %s\n", __LINE__, tsk_strerror(val));           \
        exit(EXIT_FAILURE);                                                             \
    }

int
main(int argc, char **argv)
{
    int ret;
    int j = 0;
    tsk_table_collection_t tables;

    ret = tsk_table_collection_init(&tables, 0);
    check_tsk_error(ret);

    while (true) {
        ret = tsk_table_collection_loadf(&tables, stdin, TSK_NO_INIT);
        if (ret == TSK_ERR_EOF) {
            break;
        }
        check_tsk_error(ret);
        fprintf(stderr, "Tree sequence %d had %lld mutations\n", j,
            (long long) tables.mutations.num_rows);
        ret = tsk_mutation_table_truncate(&tables.mutations, 0);
        check_tsk_error(ret);
        ret = tsk_table_collection_dumpf(&tables, stdout, 0);
        check_tsk_error(ret);
        j++;
    }
    tsk_table_collection_free(&tables);
    return EXIT_SUCCESS;
}
