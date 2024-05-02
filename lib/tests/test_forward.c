#include "forward.h"
#include <stdio.h>
#include <stdlib.h>

int test_genic_selection_from_forwards_1d(size_t *ret_num_steps, double **ret_time, double **ret_allele_frequency,
    unsigned int pastward_time, unsigned int N, unsigned int L, double s, double m, long SEED)
{
    printf("running forward\n");    
    sweep_freqs *freqs;
    freqs = forward_1d(L, N, s, m, pastward_time, SEED);
    printf("running forward done\n");    

    /* Initialization */
    for (unsigned int i = 0; i < L; i++){
        ret_num_steps[i] = 0;
        for (int j = 0; j < freqs->length; j++) {
            ret_allele_frequency[i][j] = 0;
            ret_time[i][j] = 0;
        }
    }
    printf("test_genic_selection_from_forwards_1d finish\n");

    double t = 0;
    int current_gen = 0;
    while (current_gen < freqs->length) {
        size_t pop_id = freqs->pop_id_records[freqs->length - current_gen - 1];
        double x = (double) freqs->sweep_frequencies[freqs->length - current_gen - 1];
        
        ret_allele_frequency[pop_id][current_gen] = x / N;

        t += 0.000001;
        ret_time[pop_id][current_gen] = t;   
        
        ret_num_steps[pop_id]++;
        current_gen++;
  }

      for (int j = 0; j < L; j++) {
        printf("j:%zu, %zu\n", j, ret_num_steps[j]);
        for (int k = 0; k < freqs->length; k++) {
            printf("[%zu][%zu] time: %f, allele_frequency: %f\n", j, k, ret_time[j][k], ret_allele_frequency[j][k]);
        }
    }
    return 0;
}

int
main(int argc, char **argv)
{
    
    int num_populations = 3; // L
    
    int N = 1000;
    
    double s = 1;
    double mig = 0.15;
    int tfinal = 1000;
    long SEED = 1;


    /* Compute trajactory */
    double **allele_frequency = calloc(num_populations, sizeof(double*));
    double **time = calloc(num_populations, sizeof(double*));
    size_t *num_steps = calloc(num_populations, sizeof(size_t));  

    for (int i = 0; i < num_populations; i++) {
        allele_frequency[i] = calloc(tfinal, sizeof(double));
        time[i] = calloc(tfinal, sizeof(double));
    }

    test_genic_selection_from_forwards_1d(num_steps, time, allele_frequency, tfinal, N, num_populations, s, mig, SEED);

    for (int j = 0; j < num_populations; j++) {
        printf("j:%zu, %f\n", j, num_steps[j]);
        for (int k = 0; k < tfinal; k++) {
            printf("[%zu][%zu] time: %f, allele_frequency: %f\n", j, k, time[j][k], allele_frequency[j][k]);
        }
    }
}
