#include <stdio.h>
#include <stdlib.h>
#include <err.h>

#include <tskit.h>

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        errx(EXIT_FAILURE, "line %d: %s", __LINE__, tsk_strerror(val));                 \
    }

static void
traverse_standard(const tsk_tree_t *tree)
{
    int ret;
    tsk_size_t num_nodes, j;
    tsk_id_t *nodes = malloc(tsk_tree_get_size_bound(tree) * sizeof(*nodes));

    if (nodes == NULL) {
        errx(EXIT_FAILURE, "Out of memory");
    }
    ret = tsk_tree_preorder(tree, nodes, &num_nodes);
    check_tsk_error(ret);
    for (j = 0; j < num_nodes; j++) {
        printf("Visit preorder %lld\n", (long long) nodes[j]);
    }

    ret = tsk_tree_postorder(tree, nodes, &num_nodes);
    check_tsk_error(ret);
    for (j = 0; j < num_nodes; j++) {
        printf("Visit postorder %lld\n", (long long) nodes[j]);
    }

    free(nodes);
}

static void
_traverse(const tsk_tree_t *tree, tsk_id_t u, int depth)
{
    tsk_id_t v;
    int j;

    for (j = 0; j < depth; j++) {
        printf("    ");
    }
    printf("Visit recursive %lld\n", (long long) u);
    for (v = tree->left_child[u]; v != TSK_NULL; v = tree->right_sib[v]) {
        _traverse(tree, v, depth + 1);
    }
}

static void
traverse_recursive(const tsk_tree_t *tree)
{
    _traverse(tree, tree->virtual_root, -1);
}

static void
traverse_stack(const tsk_tree_t *tree)
{
    int stack_top;
    tsk_id_t u, v;
    tsk_id_t *stack = malloc(tsk_tree_get_size_bound(tree) * sizeof(*stack));

    if (stack == NULL) {
        errx(EXIT_FAILURE, "Out of memory");
    }
    stack_top = 0;
    stack[stack_top] = tree->virtual_root;
    while (stack_top >= 0) {
        u = stack[stack_top];
        stack_top--;
        printf("Visit stack %lld\n", (long long) u);
        /* Put nodes on the stack right-to-left, so we visit in left-to-right */
        for (v = tree->right_child[u]; v != TSK_NULL; v = tree->left_sib[v]) {
            stack_top++;
            stack[stack_top] = v;
        }
    }
    free(stack);
}

static void
traverse_upwards(const tsk_tree_t *tree)
{
    const tsk_id_t *samples = tsk_treeseq_get_samples(tree->tree_sequence);
    tsk_size_t num_samples = tsk_treeseq_get_num_samples(tree->tree_sequence);
    tsk_size_t j;
    tsk_id_t u;

    for (j = 0; j < num_samples; j++) {
        u = samples[j];
        while (u != TSK_NULL) {
            printf("Visit upwards: %lld\n", (long long) u);
            u = tree->parent[u];
        }
    }
}

int
main(int argc, char **argv)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tree_t tree;

    if (argc != 2) {
        errx(EXIT_FAILURE, "usage: <tree sequence file>");
    }
    ret = tsk_treeseq_load(&ts, argv[1], 0);
    check_tsk_error(ret);
    ret = tsk_tree_init(&tree, &ts, 0);
    check_tsk_error(ret);
    ret = tsk_tree_first(&tree);
    check_tsk_error(ret);

    traverse_standard(&tree);

    traverse_recursive(&tree);

    traverse_stack(&tree);

    traverse_upwards(&tree);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
    return 0;
}
