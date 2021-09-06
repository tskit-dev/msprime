/*
** Copyright (C) 2019-2020 University of Oxford
**
** This file is part of msprime.
**
** msprime is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** msprime is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with msprime.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Implementation of log-likelihood functions for full-ARG.
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "msprime.h"
#include "likelihood.h"

static double
get_total_material(tsk_treeseq_t *ts)
{
    double ret = 0;
    tsk_size_t j;
    const tsk_edge_table_t *edges = &ts->tables->edges;
    const tsk_node_table_t *nodes = &ts->tables->nodes;
    tsk_id_t parent;
    tsk_id_t child;

    for (j = 0; j < edges->num_rows; j++) {
        parent = edges->parent[j];
        child = edges->child[j];
        ret += (edges->right[j] - edges->left[j])
               * (nodes->time[parent] - nodes->time[child]);
    }
    return ret;
}

int
msp_unnormalised_log_likelihood_mut(tsk_treeseq_t *ts, double mu, double *r_lik)
{
    int ret = 0;
    tsk_size_t j;
    tsk_mutation_t mut;
    double num_mutations = (double) tsk_treeseq_get_num_mutations(ts);
    const double total_material = get_total_material(ts);
    const double *node_time = ts->tables->nodes.time;
    double branch_length, lik;
    tsk_tree_t tree;
    tsk_id_t parent, child, it;

    ret = tsk_tree_init(&tree, ts, 0);
    if (ret != 0) {
        goto out;
    }
    if (total_material > 0 && mu > 0) {
        lik = num_mutations * log(total_material * mu) - total_material * mu;
        for (it = tsk_tree_first(&tree); it == 1; it = tsk_tree_next(&tree)) {
            for (j = 0; j < tree.sites_length; j++) {
                if (tree.sites[j].mutations_length != 1) {
                    ret = MSP_ERR_BAD_PARAM_VALUE;
                    goto out;
                }
                mut = tree.sites[j].mutations[0];
                child = mut.node;
                parent = tree.parent[child];
                branch_length = node_time[parent] - node_time[child];
                child = parent;
                parent = tree.parent[parent];
                while (parent != TSK_NULL
                       && tree.left_child[child] == tree.right_child[child]) {
                    branch_length += node_time[parent] - node_time[child];
                    child = parent;
                    parent = tree.parent[parent];
                }
                child = mut.node;
                while (tree.left_child[child] != TSK_NULL
                       && tree.left_child[child] == tree.right_child[child]) {
                    // unary nodes have left_child[node] == right_child[node]
                    // so this measures the valid leafwards branch length on
                    // which mutation mut could have taken place
                    parent = child;
                    child = tree.left_child[child];
                    branch_length += node_time[parent] - node_time[child];
                }
                lik += log(branch_length / total_material);
            }
        }
        if (it < 0) {
            ret = it;
            goto out;
        }
    } else if (num_mutations > 0) {
        lik = -DBL_MAX;
    } else {
        lik = 0;
    }
    *r_lik = lik;
out:
    tsk_tree_free(&tree);
    return ret;
}

int
msp_log_likelihood_arg(tsk_treeseq_t *ts, double r, double Ne, double *r_lik)
{
    int ret = 0;
    tsk_id_t i;
    double lineages = (double) tsk_treeseq_get_num_samples(ts);
    double sim_time = 0;
    double material = lineages * tsk_treeseq_get_sequence_length(ts);
    double material_in_children, material_in_parent, rate, gap;
    double lik = 0;
    const tsk_edge_table_t *edges = &ts->tables->edges;
    const tsk_node_table_t *nodes = &ts->tables->nodes;
    tsk_id_t *first_parent_edge = NULL;
    tsk_id_t *last_parent_edge = NULL;
    tsk_id_t edge = 0;
    tsk_id_t parent;

    if (Ne <= 0) {
        ret = MSP_ERR_BAD_POPULATION_SIZE;
        goto out;
    }

    first_parent_edge = malloc(nodes->num_rows * sizeof(tsk_id_t));
    last_parent_edge = malloc(nodes->num_rows * sizeof(tsk_id_t));
    if (first_parent_edge == NULL || last_parent_edge == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memset(first_parent_edge, TSK_NULL, sizeof(tsk_id_t) * nodes->num_rows);
    memset(last_parent_edge, TSK_NULL, sizeof(tsk_id_t) * nodes->num_rows);
    for (i = 0; i < (tsk_id_t) edges->num_rows; i++) {
        if (first_parent_edge[edges->child[i]] == TSK_NULL) {
            first_parent_edge[edges->child[i]] = i;
        }
        last_parent_edge[edges->child[i]] = i;
    }
    while (edge < (tsk_id_t) edges->num_rows && lineages > 0) {
        rate = lineages * (lineages - 1) / (4 * Ne) + material * r;
        parent = edges->parent[edge];
        lik -= rate * (nodes->time[parent] - sim_time);
        sim_time = nodes->time[parent];
        if (nodes->flags[parent] & MSP_NODE_IS_RE_EVENT) {
            if (r <= 0) {
                *r_lik = -DBL_MAX;
                ret = 0;
                goto out;
            }
            while (edge < (tsk_id_t) edges->num_rows && edges->parent[edge] == parent) {
                edge++;
            }
            gap = edges->left[edge] - edges->right[edge - 1];
            edge = last_parent_edge[edges->child[edge]];
            material -= gap;
            lineages++;
            if (gap <= 0) {
                // we evaluate the density rather than probability
                gap = 1;
            }
            lik += log(r * gap);
        } else {
            material_in_children = -edges->left[edge];
            edge = last_parent_edge[edges->child[edge]];
            material_in_children += edges->right[edge];
            edge++;
            material_in_children -= edges->left[edge];
            edge = last_parent_edge[edges->child[edge]];
            material_in_children += edges->right[edge];
            if (first_parent_edge[parent] == -1) {
                lineages -= 2;
                material -= material_in_children;
            } else {
                material_in_parent = edges->right[last_parent_edge[parent]]
                                     - edges->left[first_parent_edge[parent]];
                lineages--;
                material -= material_in_children - material_in_parent;
            }
            lik -= log(2 * Ne);
        }
        if (lineages > 0) {
            edge++;
        }
    }
    *r_lik = lik;
    ret = 0;
out:
    msp_safe_free(first_parent_edge);
    msp_safe_free(last_parent_edge);
    return ret;
}
