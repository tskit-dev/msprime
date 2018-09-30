#ifndef TSK_STATS_H
#define TSK_STATS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "tsk_trees.h"

typedef struct {
    tsk_tree_t *outer_tree;
    tsk_tree_t *inner_tree;
    size_t num_sites;
    int tree_changed;
    tsk_treeseq_t *tree_sequence;
} tsk_ld_calc_t;

int tsk_ld_calc_alloc(tsk_ld_calc_t *self, tsk_treeseq_t *tree_sequence);
int tsk_ld_calc_free(tsk_ld_calc_t *self);
void tsk_ld_calc_print_state(tsk_ld_calc_t *self, FILE *out);
int tsk_ld_calc_get_r2(tsk_ld_calc_t *self, size_t a, size_t b, double *r2);
int tsk_ld_calc_get_r2_array(tsk_ld_calc_t *self, size_t a, int direction,
        size_t max_mutations, double max_distance,
        double *r2, size_t *num_r2_values);


#ifdef __cplusplus
}
#endif
#endif
