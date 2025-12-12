/*
** Copyright (C) 2016-2020 University of Oxford
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

#include "testlib.h"

static void
verify_simple_genic_selection_trajectory(
    double start_frequency, double end_frequency, double alpha, double dt)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = { { 0, 0.0 }, { 0, 0.0 } };
    tsk_table_collection_t tables;
    size_t j, num_steps;
    double *allele_frequency, *time;

    ret = build_sim(&msp, &tables, rng, 1.0, 1, samples, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_simulation_model_sweep_genic_selection(
        &msp, 0.5, start_frequency, end_frequency, alpha, dt);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* compute the trajectory */
    ret = msp.model.params.sweep.generate_trajectory(
        &msp.model.params.sweep, &msp, &num_steps, &time, &allele_frequency);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(num_steps > 1);
    CU_ASSERT_EQUAL(time[0], 0);
    CU_ASSERT_EQUAL(allele_frequency[0], end_frequency);
    CU_ASSERT_TRUE(allele_frequency[num_steps - 1] == start_frequency);

    for (j = 0; j < num_steps; j++) {
        CU_ASSERT_TRUE(allele_frequency[j] >= 0);
        CU_ASSERT_TRUE(allele_frequency[j] <= 1);
        if (j > 0) {
            CU_ASSERT_DOUBLE_EQUAL_FATAL(time[j], time[j - 1] + dt, 1e-9);
        }
    }

    free(time);
    free(allele_frequency);
    msp_free(&msp);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_genic_selection_trajectory(void)
{
    verify_simple_genic_selection_trajectory(0.1, 0.9, 200, 0.0001);
    verify_simple_genic_selection_trajectory(0.1, 0.9, 101, 0.000001);
    verify_simple_genic_selection_trajectory(0.8, 0.9, 500, 0.000002);
    verify_simple_genic_selection_trajectory(0.1, 0.7, 100, 0.000001);
    verify_simple_genic_selection_trajectory(0.1, 0.4, 50, 0.0001);
}

static void
test_sweep_genic_selection_bad_parameters(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = { { 0, 0.0 }, { 0, 0.0 } };
    tsk_table_collection_t tables;

    ret = build_sim(&msp, &tables, rng, 1.0, 1, samples, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_set_simulation_model_sweep_genic_selection(
        &msp, 0.5, -0.01, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_ALLELE_FREQUENCY);
    ret = msp_set_simulation_model_sweep_genic_selection(
        &msp, 0.5, 0.01, -0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_ALLELE_FREQUENCY);
    ret = msp_set_simulation_model_sweep_genic_selection(
        &msp, 0.5, 10.01, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_ALLELE_FREQUENCY);
    ret = msp_set_simulation_model_sweep_genic_selection(
        &msp, 0.5, 0.01, 10.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_ALLELE_FREQUENCY);

    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.01, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_TRAJECTORY_START_END);
    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.1, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_TRAJECTORY_START_END);

    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.9, 0.1, 0.0);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_TIME_DELTA);
    ret = msp_set_simulation_model_sweep_genic_selection(
        &msp, 0.5, 0.1, 0.9, 0.1, -0.01);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_TIME_DELTA);

    ret = msp_set_simulation_model_sweep_genic_selection(&msp, -0.5, 0.1, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_SWEEP_POSITION);
    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 5.0, 0.1, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_SWEEP_POSITION);
    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.9, -666, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_SWEEP_GENIC_SELECTION_S);

    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    /* The incorrect number of populations was specified */
    ret = build_sim(&msp, &tables, rng, 1.0, 2, samples, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_UNSUPPORTED_OPERATION);

    msp_free(&msp);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_sweep_genic_selection_events(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    sample_t samples[] = { { 0, 0.0 }, { 0, 0.0 }, { 0, 0.0 } };
    tsk_table_collection_t tables;

    ret = build_sim(&msp, &tables, rng, 1.0, 1, samples, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_num_labels(&msp, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.1, 0, 1, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_EVENTS_DURING_SWEEP);
    msp_free(&msp);
    tsk_table_collection_free(&tables);

    samples[1].time = 0.1;
    ret = build_sim(&msp, &tables, rng, 1, 1, samples, 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 0.5, 0.1, 0.9, 0.1, 0.1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_num_labels(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_EVENTS_DURING_SWEEP);
    msp_free(&msp);

    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
verify_sweep_genic_selection(double sequence_length, double s)
{
    int j, ret;
    uint32_t n = 10;
    unsigned long seed = 133;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables[2];

    for (j = 0; j < 2; j++) {
        gsl_rng_set(rng, seed);
        ret = build_sim(&msp, &tables[j], rng, sequence_length, 1, NULL, n);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 1), 0);
        ret = msp_set_num_labels(&msp, 2);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_set_simulation_model_sweep_genic_selection(
            &msp, sequence_length / 2, 0.1, 0.9, s, 1e-6);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = msp_initialise(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
        CU_ASSERT_TRUE(ret >= 0);
        msp_print_state(&msp, _devnull);
        ret = msp_finalise_tables(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        msp_free(&msp);
        CU_ASSERT_EQUAL(tables[j].migrations.num_rows, 0);
        CU_ASSERT(tables[j].nodes.num_rows > 0);
        CU_ASSERT(tables[j].edges.num_rows > 0);
    }
    CU_ASSERT_TRUE(tsk_node_table_equals(&tables[0].nodes, &tables[1].nodes, 0));
    CU_ASSERT_TRUE(tsk_edge_table_equals(&tables[0].edges, &tables[1].edges, 0));
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    for (j = 0; j < 2; j++) {
        tsk_table_collection_free(&tables[j]);
    }
}

static void
test_sweep_genic_selection_single_locus(void)
{
    verify_sweep_genic_selection(1, 0.02);
    verify_sweep_genic_selection(1, 0.01);
    verify_sweep_genic_selection(1, 0.06);
}

static void
test_sweep_genic_selection_recomb(void)
{
    verify_sweep_genic_selection(10, 0.2);
    verify_sweep_genic_selection(10, 1.0);
    verify_sweep_genic_selection(10, 1e-2);
}

static void
test_sweep_genic_selection_gc(void)
{
    int ret;
    uint32_t n = 100;
    unsigned long seed = 133;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    gsl_rng_set(rng, seed);
    ret = build_sim(&msp, &tables, rng, 10.0, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 1), 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_gene_conversion_rate(&msp, 1), 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_gene_conversion_tract_length(&msp, 1), 0);
    ret = msp_set_num_labels(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, 100);
    CU_ASSERT_TRUE(msp_get_num_gene_conversion_events(&msp) > 0);
    CU_ASSERT_TRUE(ret >= 0);
    ret = msp_set_simulation_model_sweep_genic_selection(&msp, 5, 0.1, 0.9, 0.1, 0.01);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    /* GC rate in sweep model is not implemented */
    CU_ASSERT_TRUE(ret == MSP_ERR_SWEEPS_GC_NOT_SUPPORTED);
    msp_print_state(&msp, _devnull);
    ret = msp_finalise_tables(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    msp_free(&msp);

    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_sweep_genic_selection_time_change(void)
{
    int j, ret;
    uint32_t n = 10;
    double num_loci = 10;
    unsigned long seed = 133234;
    double t;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    gsl_rng_set(rng, seed);
    ret = build_sim(&msp, &tables, rng, num_loci, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, 1), 0);
    ret = msp_set_num_labels(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    /* Run some time and then change the model */
    ret = msp_run(&msp, 0.125, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MAX_TIME);
    t = msp_get_time(&msp);

    for (j = 0; j < 10; j++) {

        ret = msp_set_simulation_model_sweep_genic_selection(
            &msp, num_loci / 2, 0.1, 0.9, 0.1, 0.01);
        CU_ASSERT_EQUAL(ret, 0);

        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(msp_get_time(&msp), t);

        ret = msp_set_simulation_model_hudson(&msp);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(msp_get_time(&msp), t);
    }
    msp_free(&msp);

    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
sweep_genic_selection_mimic_msms_single_run(unsigned long int seed)
{

    /* Try to mimic the msms parameters used in verification.py
           "100 300 -t 200 -r 200 500000"
           " -SF 0 0.9 -Sp 0.5 -SaA 5000 -SAA 10000 -N 10000"
     */
    int ret;
    uint32_t n = 10;
    double num_loci = 500001;
    double position = num_loci / 2;
    double s = 10000;
    double recom_rate = 0.0004;
    double start_frequency = 0.5 / 1000;
    double end_frequency = 1 - 0.5 / 1000;
    double dt = 1.0 / 400000;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;

    // Test over differnt seeds
    gsl_rng_set(rng, seed);

    ret = build_sim(&msp, &tables, rng, num_loci, 1, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, recom_rate), 0);
    ret = msp_set_num_labels(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);

    // To mimic the verfication.py call
    msp_set_discrete_genome(&msp, 0);
    msp_set_gene_conversion_rate(&msp, 0);
    msp_set_gene_conversion_tract_length(&msp, 1);
    msp_set_avl_node_block_size(&msp, 65536);
    msp_set_node_mapping_block_size(&msp, 65536);
    msp_set_segment_block_size(&msp, 65536);

    ret = msp_set_simulation_model_sweep_genic_selection(
        &msp, position, start_frequency, end_frequency, s, dt);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);

    msp_verify(&msp, 0);

    msp_free(&msp);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

static void
test_sweep_genic_selection_mimic_msms(void)
{
    /* To mimic the nrepeats = 300  parameter in msms cmdline arguments*/
    for (int i = 0; i < 1; i++) {
        sweep_genic_selection_mimic_msms_single_run(i + 1);
    }
}

static void
verify_sweep_reverse_wm_no_recomb(unsigned long int seed)
{
    int ret;
    uint32_t n = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    int num_demes = 1;
    tsk_id_t mut_pop[] = { 10001 };
    int num_events = mut_pop[0] - 1;
    double *time_of_ev;
    int *ev_type;
    int *start_deme;
    int *end_deme;
    int r, i;
    double p;

    tsk_id_t final_mut_pop[] = { 1 };

    gsl_rng_set(rng, seed);

    time_of_ev = (double *) malloc(sizeof(double) * num_events);
    ev_type = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);
    end_deme = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);
    start_deme = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);

    r = rand() % RAND_MAX;
    p = (double) r / RAND_MAX;
    time_of_ev[0] = 1000 - p;
    ev_type[0] = 0;
    start_deme[0] = 0;
    end_deme[0] = 0;

    for (i = 1; i < num_events; i++) {
        r = rand() % RAND_MAX;
        p = (double) r / RAND_MAX;
        time_of_ev[i] = time_of_ev[i - 1] - p;
        ev_type[i] = 0;
        start_deme[i] = 0;
        end_deme[i] = 0;
    }

    ret = build_sim(&msp, &tables, rng, 1.0, num_demes, NULL, n);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = msp_set_num_labels(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);

    msp_set_discrete_genome(&msp, 0);
    msp_set_gene_conversion_rate(&msp, 0);
    msp_set_gene_conversion_tract_length(&msp, 1);
    msp_set_avl_node_block_size(&msp, 65536);
    msp_set_node_mapping_block_size(&msp, 65536);
    msp_set_segment_block_size(&msp, 65536);

    ret = msp_set_simulation_model_sweep_genic_selection_reverse(&msp, 0.5, num_events,
        num_demes, num_events + 1, 0.0, mut_pop, final_mut_pop, time_of_ev, ev_type,
        start_deme, end_deme);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);

    msp_verify(&msp, 0);

    free(time_of_ev);
    free(ev_type);
    free(start_deme);
    free(end_deme);
    msp_free(&msp);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_sweep_reverse_wm_no_recomb(void)
{
    /* To mimic the nrepeats = 300  parameter in msms cmdline arguments*/
    for (int i = 0; i < 300; i++) {
        verify_sweep_reverse_wm_no_recomb(i + 256329651209);
    }
}

static void
verify_sweep_reverse_wm(unsigned long int seed)
{
    int ret;
    uint32_t n = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    tsk_id_t mut_pop[] = { 10001 };
    int num_events = mut_pop[0] - 1;
    int num_demes = 1;
    double *time_of_ev;
    int *ev_type;
    int *start_deme;
    int *end_deme;
    int r, i;
    double p;
    double num_loci = 5001;
    double position = num_loci / 2;
    double recom_rate = 0.00004;
    tsk_id_t final_mut_pop[] = { 1 };

    time_of_ev = (double *) malloc(sizeof(double) * num_events);
    ev_type = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);
    end_deme = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);
    start_deme = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);

    // Test over differnt seeds
    gsl_rng_set(rng, seed);

    r = rand() % RAND_MAX;
    p = (double) r / RAND_MAX;
    time_of_ev[0] = num_events - p;
    ev_type[0] = 0;
    start_deme[0] = 0;
    end_deme[0] = 0;

    for (i = 1; i < num_events; i++) {
        r = rand() % RAND_MAX;
        p = (double) r / RAND_MAX;
        time_of_ev[i] = time_of_ev[i - 1] - p;
        ev_type[i] = 0;
        start_deme[i] = 0;
        end_deme[i] = 0;
    }

    ret = build_sim(&msp, &tables, rng, num_loci, num_demes, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, recom_rate), 0);
    ret = msp_set_num_labels(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);

    // To mimic the verfication.py call
    msp_set_discrete_genome(&msp, 0);
    msp_set_gene_conversion_rate(&msp, 0);
    msp_set_gene_conversion_tract_length(&msp, 1);
    msp_set_avl_node_block_size(&msp, 65536);
    msp_set_node_mapping_block_size(&msp, 65536);
    msp_set_segment_block_size(&msp, 65536);

    ret = msp_set_simulation_model_sweep_genic_selection_reverse(&msp, position,
        num_events, num_demes, num_events + 1, 0.0, mut_pop, final_mut_pop, time_of_ev,
        ev_type, start_deme, end_deme);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);

    msp_verify(&msp, 0);

    free(time_of_ev);
    free(ev_type);
    free(start_deme);
    free(end_deme);
    msp_free(&msp);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_sweep_reverse_wm(void)
{
    /* To mimic the nrepeats = 300  parameter in msms cmdline arguments*/
    for (int i = 0; i < 300; i++) {
        verify_sweep_reverse_wm(i + 6584367292);
    }
}

static void
verify_sweep_reverse_wm_back_and_forth(unsigned long int seed)
{
    int ret;
    uint32_t n = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    tsk_id_t mut_pop[] = { 1001 };
    int num_events = 3 * mut_pop[0] - 5;
    int num_demes = 1;
    double *time_of_ev;
    int *ev_type;
    int *start_deme;
    int *end_deme;
    int r, i;
    double p;
    double num_loci = 5001;
    double position = num_loci / 2;
    double recom_rate = 0.00004;
    tsk_id_t final_mut_pop[] = { 1 };

    gsl_rng_set(rng, seed);

    time_of_ev = (double *) malloc(sizeof(double) * num_events);
    ev_type = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);
    end_deme = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);
    start_deme = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);

    r = rand() % RAND_MAX;
    p = (double) r / RAND_MAX;
    time_of_ev[0] = num_events - p;
    ev_type[0] = 0;
    start_deme[0] = 0;
    end_deme[0] = 0;

    for (i = 1; i < num_events; i++) {
        r = rand() % RAND_MAX;
        p = (double) r / RAND_MAX;
        time_of_ev[i] = time_of_ev[i - 1] - p;
        if (i % 3 == 2) {
            ev_type[i] = 1;
        } else {
            ev_type[i] = 0;
        }
        start_deme[i] = 0;
        end_deme[i] = 0;
    }

    ret = build_sim(&msp, &tables, rng, num_loci, num_demes, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, recom_rate), 0);
    ret = msp_set_num_labels(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);

    msp_set_discrete_genome(&msp, 0);
    msp_set_gene_conversion_rate(&msp, 0);
    msp_set_gene_conversion_tract_length(&msp, 1);
    msp_set_avl_node_block_size(&msp, 65536);
    msp_set_node_mapping_block_size(&msp, 65536);
    msp_set_segment_block_size(&msp, 65536);

    ret = msp_set_simulation_model_sweep_genic_selection_reverse(&msp, position,
        num_events, num_demes, mut_pop[0], 0.0, mut_pop, final_mut_pop, time_of_ev,
        ev_type, start_deme, end_deme);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);

    msp_verify(&msp, 0);

    free(time_of_ev);
    free(ev_type);
    free(start_deme);
    free(end_deme);
    msp_free(&msp);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_sweep_reverse_wm_back_and_forth(void)
{
    /* To mimic the nrepeats = 300  parameter in msms cmdline arguments*/
    for (int i = 0; i < 300; i++) {
        verify_sweep_reverse_wm_back_and_forth(i + 462523296209);
    }
}

static void
verify_sweep_reverse_wm_oscillating(unsigned long int seed)
{
    int ret;
    uint32_t n = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    tsk_id_t mut_pop[] = { 10001 };
    int num_events = 3 * mut_pop[0] - 9;
    int num_demes = 1;
    double *time_of_ev;
    int *ev_type;
    int *start_deme;
    int *end_deme;
    int r, i;
    double p;
    double num_loci = 5001;
    double position = num_loci / 2;
    double recom_rate = 0.00004;
    tsk_id_t final_mut_pop[] = { 1 };

    gsl_rng_set(rng, seed);

    time_of_ev = (double *) malloc(sizeof(double) * num_events);
    ev_type = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);
    end_deme = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);
    start_deme = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);

    r = rand() % RAND_MAX;
    p = (double) r / RAND_MAX;
    time_of_ev[0] = num_events - p;
    ev_type[0] = 0;
    start_deme[0] = 0;
    end_deme[0] = 0;

    for (i = 1; i < num_events; i++) {
        r = rand() % RAND_MAX;
        p = (double) r / RAND_MAX;
        time_of_ev[i] = time_of_ev[i - 1] - p;
        if ((i < 2 * mut_pop[0] - 6) && (i > mut_pop[0] - 3)) {
            ev_type[i] = 1;
        } else {
            ev_type[i] = 0;
        }
        start_deme[i] = 0;
        end_deme[i] = 0;
    }

    ret = build_sim(&msp, &tables, rng, num_loci, num_demes, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, recom_rate), 0);
    ret = msp_set_num_labels(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);

    msp_set_discrete_genome(&msp, 0);
    msp_set_gene_conversion_rate(&msp, 0);
    msp_set_gene_conversion_tract_length(&msp, 1);
    msp_set_avl_node_block_size(&msp, 65536);
    msp_set_node_mapping_block_size(&msp, 65536);
    msp_set_segment_block_size(&msp, 65536);

    ret = msp_set_simulation_model_sweep_genic_selection_reverse(&msp, position,
        num_events, num_demes, mut_pop[0], 0.0, mut_pop, final_mut_pop, time_of_ev,
        ev_type, start_deme, end_deme);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);

    msp_verify(&msp, 0);

    free(time_of_ev);
    free(ev_type);
    free(start_deme);
    free(end_deme);
    msp_free(&msp);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_sweep_reverse_wm_oscillating(void)
{
    /* To mimic the nrepeats = 300  parameter in msms cmdline arguments*/
    for (int i = 0; i < 300; i++) {
        verify_sweep_reverse_wm_oscillating(i + 462941996216);
    }
}

static void
verify_sweep_reverse_wm_sudden_drop(unsigned long int seed)
{
    int ret;
    uint32_t n = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    tsk_id_t mut_pop[] = { 10001 };
    int num_events = mut_pop[0] - 1;
    int num_demes = 1;
    double *time_of_ev;
    int *ev_type;
    int *start_deme;
    int *end_deme;
    int r, i;
    double p;
    double num_loci = 5001;
    double position = num_loci / 2;
    double recom_rate = 0.0004;
    tsk_id_t final_mut_pop[] = { 1 };

    time_of_ev = (double *) malloc(sizeof(double) * num_events);
    ev_type = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);
    end_deme = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);
    start_deme = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);

    // Test over differnt seeds
    gsl_rng_set(rng, seed);

    r = rand() % RAND_MAX;
    p = (double) r / RAND_MAX;
    time_of_ev[0] = num_events - p / 10000;
    ev_type[0] = 0;
    start_deme[0] = 0;
    end_deme[0] = 0;

    for (i = 1; i < num_events; i++) {
        r = rand() % RAND_MAX;
        p = (double) r / RAND_MAX;
        time_of_ev[i] = time_of_ev[i - 1] - p / 10000;
        ev_type[i] = 0;
        start_deme[i] = 0;
        end_deme[i] = 0;
    }

    ret = build_sim(&msp, &tables, rng, num_loci, num_demes, NULL, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, recom_rate), 0);
    ret = msp_set_num_labels(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);

    // To mimic the verfication.py call
    msp_set_discrete_genome(&msp, 0);
    msp_set_gene_conversion_rate(&msp, 0);
    msp_set_gene_conversion_tract_length(&msp, 1);
    msp_set_avl_node_block_size(&msp, 65536);
    msp_set_node_mapping_block_size(&msp, 65536);
    msp_set_segment_block_size(&msp, 65536);

    ret = msp_set_simulation_model_sweep_genic_selection_reverse(&msp, position,
        num_events, num_demes, num_events + 1, 0.0, mut_pop, final_mut_pop, time_of_ev,
        ev_type, start_deme, end_deme);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);

    msp_verify(&msp, 0);

    free(time_of_ev);
    free(ev_type);
    free(start_deme);
    free(end_deme);
    msp_free(&msp);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_sweep_reverse_wm_sudden_drop(void)
{
    /* To mimic the nrepeats = 300  parameter in msms cmdline arguments*/
    for (int i = 0; i < 300; i++) {
        verify_sweep_reverse_wm_sudden_drop(i + 6584367292);
    }
}

static void
verify_sweep_reverse_1D_deme_wise(unsigned long int seed)
{
    int ret;
    uint32_t n = 10;
    msp_t msp;
    gsl_rng *rng = safe_rng_alloc();
    tsk_table_collection_t tables;
    tsk_id_t *mut_pop;
    tsk_id_t tot_pop = 100;
    int num_demes = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    int num_events = num_demes * tot_pop - 1;
    double migration_rate = 0.25;
    size_t deme_index;
    double *time_of_ev;
    int *ev_type;
    int *start_deme;
    int *end_deme;
    int r, i;
    double p;
    double num_loci = 5001;
    double position = num_loci / 2;
    double recom_rate = 0.00004;
    tsk_id_t *final_mut_pop;

    time_of_ev = (double *) malloc(sizeof(double) * num_events);
    ev_type = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);
    end_deme = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);
    start_deme = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_events);
    mut_pop = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_demes);
    final_mut_pop = (tsk_id_t *) malloc(sizeof(tsk_id_t) * num_demes);

    for (i = 0; i < n; i++) {
        samples[i].time = 0;
        samples[i].population = i % num_demes;
    }

    for (i = 0; i < num_demes; i++) {
        mut_pop[i] = tot_pop;
        final_mut_pop[i] = 0;
    }
    final_mut_pop[0] = 1;

    // Test over differnt seeds
    gsl_rng_set(rng, seed);

    r = rand() % RAND_MAX;
    p = (double) r / RAND_MAX;
    time_of_ev[0] = num_events - p;
    ev_type[0] = 0;
    deme_index = num_demes - 1;
    start_deme[0] = deme_index;
    end_deme[0] = deme_index;

    for (i = 1; i < num_events; i++) {
        r = rand() % RAND_MAX;
        p = (double) r / RAND_MAX;
        time_of_ev[i] = time_of_ev[i - 1] - p;
        if (i % tot_pop == (tot_pop - 1)) {
            ev_type[i] = 2;
            end_deme[i] = deme_index - 1;
            start_deme[i] = deme_index;
            deme_index -= 1;
        } else {
            ev_type[i] = 0;
            start_deme[i] = deme_index;
            end_deme[i] = deme_index;
        }
    }

    ret = build_sim(&msp, &tables, rng, num_loci, num_demes, samples, n);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(msp_set_recombination_rate(&msp, recom_rate), 0);
    ret = msp_set_num_labels(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);

    // To mimic the verfication.py call
    msp_set_discrete_genome(&msp, 0);
    msp_set_gene_conversion_rate(&msp, 0);
    msp_set_gene_conversion_tract_length(&msp, 1);
    msp_set_avl_node_block_size(&msp, 65536);
    msp_set_node_mapping_block_size(&msp, 65536);
    msp_set_segment_block_size(&msp, 65536);

    ret = msp_set_simulation_model_sweep_genic_selection_reverse(&msp, position,
        num_events, num_demes, num_events + 1, migration_rate, mut_pop, final_mut_pop,
        time_of_ev, ev_type, start_deme, end_deme);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
    CU_ASSERT_EQUAL(ret, MSP_EXIT_MODEL_COMPLETE);

    msp_verify(&msp, 0);

    free(time_of_ev);
    free(ev_type);
    free(start_deme);
    free(end_deme);
    free(mut_pop);
    free(final_mut_pop);
    free(samples);
    msp_free(&msp);
    tsk_table_collection_free(&tables);
    gsl_rng_free(rng);
}

static void
test_sweep_reverse_1D_deme_wise(void)
{
    /* To mimic the nrepeats = 300  parameter in msms cmdline arguments*/
    for (int i = 0; i < 300; i++) {
        verify_sweep_reverse_1D_deme_wise(i + 6584367292);
    }
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_genic_selection_trajectory", test_genic_selection_trajectory },
        { "test_sweep_genic_selection_bad_parameters",
         test_sweep_genic_selection_bad_parameters },
        { "test_sweep_genic_selection_events", test_sweep_genic_selection_events },
        { "test_sweep_genic_selection_single_locus",
         test_sweep_genic_selection_single_locus },
        { "test_sweep_genic_selection_recomb", test_sweep_genic_selection_recomb },
        { "test_sweep_genic_selection_gc", test_sweep_genic_selection_gc },
        { "test_sweep_genic_selection_time_change",
         test_sweep_genic_selection_time_change },
        { "test_sweep_genic_selection_mimic_msms",
            test_sweep_genic_selection_mimic_msms },
        { "test_sweep_reverse_wm_no_recomb",
            test_sweep_reverse_wm_no_recomb },
        { "test_sweep_reverse_wm",
            test_sweep_reverse_wm },
        { "test_sweep_reverse_wm_back_and_forth",
            test_sweep_reverse_wm_back_and_forth },
        { "test_sweep_reverse_wm_oscillating",
            test_sweep_reverse_wm_oscillating },
        { "test_sweep_reverse_wm_sudden_drop",
            test_sweep_reverse_wm_sudden_drop },
        //{ "test_sweep_reverse_1D_deme_wise", test_sweep_reverse_1D_deme_wise },
        CU_TEST_INFO_NULL,
    };

    return test_main(tests, argc, argv);
}
