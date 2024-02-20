#include <stdio.h>
#include "forward.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

const gsl_rng_type *T;
gsl_rng *R;

/*function roundint to round a double to an integer*/
int roundint(double x) {
	int result;
	
	result = floor(x);
	if (x - result >= 0.5) result++;
	
	return result;
}

/**
 * 1D Forward Simulation
*/

void next_gen_1d(unsigned int n[], double s, double mig, int L, unsigned int N) {

    double x[L]; // frequencies
    double xmig[L]; // after migration
    int i;

	// Initialize frequencies:
	for (i = 0; i < L; i++){
		x[i] = (double)n[i] / N;
	}
	// Migration:
	for (i = 1; i < L - 1; i++){
		xmig[i] = x[i] + mig * (0.5 * (x[i-1] + x[i+1]) - x[i]);
	}
	xmig[0] = x[0] + mig * 0.5 * (x[1] - x[0]);
	xmig[L-1] = x[L-1] + mig * 0.5 * (x[L-2] - x[L-1]);

    // Sampling and selection within demes
    for (i = 0; i < L; i++) {
		n[i] = gsl_ran_binomial(R, xmig[i] + s * xmig[i] * (1 - xmig[i]), N);
	}

}

void print_sweep_freqs(sweep_freqs* freqs) {
    fprintf(stderr, "freqs->length: %d\n", freqs->length);
    for (int i = 0; i < freqs->length; i++) {
        if (freqs->sweep_frequencies != NULL) {
            fprintf(stderr, "freqs->sweep_frequencies[%u]: %u\n", i, freqs->sweep_frequencies[i]);
        } else {
            fprintf(stderr, "freqs->sweep_frequencies NULL \n");
        }
        
        if (freqs->pop_id_records != NULL) {
            fprintf(stderr, "freqs->pop_id_records[%u]: %u\n", i, freqs->pop_id_records[i]);
        } else {
            fprintf(stderr, "freqs->pop_id_records NULL \n");
        }
    }
    return;
}

sweep_freqs* forward_1d(unsigned int L, unsigned int N, double s, double mig, unsigned int tfinal, long SEED) {
    gsl_rng_env_setup();
    T = gsl_rng_default;
    R = gsl_rng_alloc (T);
    gsl_rng_set(R,SEED);

    // Initialize population:
    unsigned int* n = malloc(sizeof(unsigned int));

    if (L * N * s < 100){
    	fprintf(stderr, "Error: meta-population is too small for selection to be strong.\n");
    	return NULL;
    }

    for (unsigned int i = 0; i < L; i++){
        n[i] = 0;
    }
    
    unsigned int t; //Current time. Cannot exceed tfinal.
    unsigned int ntot; //Frequency count

    unsigned int length = 0;
    unsigned int* sweep_frequencies = NULL;
    unsigned int* pop_id_records = NULL;
    //sweep_freqs* result = (sweep_freqs*) malloc(sizeof(sweep_freqs));
    while (n[0] == 0){
        
        length = 0;
        
        if (sweep_frequencies != NULL) {
            free(sweep_frequencies);
        }
        sweep_frequencies = NULL;
        
        if (pop_id_records != NULL) {
            free(pop_id_records);
        }
        pop_id_records = NULL;

        for (unsigned int i = 0; i < L; i++){
            n[i] = 0;
        }
        n[0] = 1;
        
        sweep_frequencies = (unsigned int*) malloc(tfinal * L * sizeof(unsigned int));    
        pop_id_records = (unsigned int*) malloc(tfinal * L * sizeof(unsigned int));    
        for (t = 0; t < tfinal; t++) {
            // Record the status of the population and check for fixation:
            ntot = 0;
            for (unsigned int i = 0; i < L; i++){
                ntot += n[i]; 
                sweep_frequencies[t * L + i] = n[i];
                length = (t * L) + i;
                pop_id_records[t * L + i] = i;
            }
        
            // Stop if one of the alleles is fixed or extinct:
            if ((ntot == 0) || (ntot == N * L)) {                
                break;
            }
            
            // Evolve for one generation
            next_gen_1d(n, s, mig, L, N);

        }
    }
    
    sweep_freqs* result = (sweep_freqs*) malloc(sizeof(sweep_freqs));
    result->length = length;
    result->sweep_frequencies = sweep_frequencies;
    result->pop_id_records = pop_id_records;
    print_sweep_freqs(result);
    
    if (t == tfinal){
		fprintf(stderr, "Error: Simulation finished without fixation.\n");
		return NULL;
	}

    return result;
}



/**
 * 2D Forward Simulation
*/

void next_gen_2d(int L, unsigned int n[L][L], double s, double mig, unsigned int N) {

    double x[L][L]; // frequencies
    double xmig[L][L]; // after migration
    int i;
    int j;

    // Initialize frequencies:
    
    for (i = 0; i < L; i++)
    {
        for(j = 0; j < L; j++)
        {
            x[i][j] = (double)n[i][j] / N;
        }
    }

    // Migration:
    for (i = 1; i < L - 1; i++)
    {
        for(j = 1; j < L - 1; j++)
        {
            xmig[i][j] = x[i][j] + mig * (0.25 * (x[i-1][j] + x[i+1][j] + x[i][j-1] + x[i][j+1] ) - x[i][j]);
        }
    }

    //Boundary Conditions (PBC)
    for (i = 1; i < L - 1; i++)
    {
        xmig[0][i] = x[0][i] + mig * 0.25 * (x[1][i] + x[L-1][i] + x[0][i-1] + x[0][i+1]) - mig*x[0][i];
        xmig[L-1][i] = x[L-1][i] + mig * 0.25 * (x[L-2][i] + x[0][i] + x[L-1][i-1] + x[L-1][i+1]) - mig*x[L-1][i];
        xmig[i][0] = x[i][0] + mig * 0.25 * (x[i][1] + x[i][L-1] + x[i-1][0] + x[i+1][0]) - mig*x[i][0];
        xmig[i][L-1] = x[i][L-1] + mig * 0.25 * (x[i][0] + x[i][L-2] + x[i-1][L-1] + x[i+1][L-1]) - mig*x[i][L-1];

    }
    
    xmig[0][0] = x[0][0] + mig * 0.25 * (x[1][0] + x[L-1][0]+ x[0][1] + x[0][L-1]) - mig*x[0][0];  // Origin at bottom right
    xmig[L-1][0] = x[L-1][0] + mig * 0.25 * (x[L-2][0] + x[0][0] + x[L-1][1] + x[L-1][L-1]) - mig*x[L-1][0];   //Top right
    xmig[0][L-1] = x[0][L-1] + mig * 0.25 * (x[0][L-2] + x[0][0] + x[1][L-1] + x[L-1][L-1])- mig*x[0][L-1];   //Bottom Left
    xmig[L-1][L-1] = x[L-1][L-1] + mig * 0.25 * (x[L-2][L-1] + x[0][L-1] +  x[L-1][L-2] + x[L-1][0]) - mig*x[L-1][L-1];  //Top Left

    // Sampling and selection within demes
    for (i = 0; i < L; i++) 
    {

    for(j = 0; j < L; j++)
        {
            n[i][j] = gsl_ran_binomial(R, xmig[i][j] + s * xmig[i][j] * (1 - xmig[i][j]), N);  // Logistic Growth and Sampling within a deme 
        }
        
    }
    
}

sweep_freqs *forward_2d (unsigned int L, unsigned int N, double s, double mig, unsigned int tfinal, unsigned int l0, long SEED) {
    gsl_rng_env_setup();
    T = gsl_rng_default;
    R = gsl_rng_alloc(T);
    gsl_rng_set(R, 0);

    unsigned int n[L][L];

    if (L * N * s < 100) {
        fprintf(stderr, "Error: meta-population is too small for selection to be strong.\n");
        return NULL;
    }

    for (unsigned int i = 0; i < L; i++) {
        for (unsigned int k = 0; k < L; k++) {
            n[i][k] = 0;
        }
    }
    unsigned int ntot;
    unsigned int t;
    sweep_freqs* result = (sweep_freqs*) malloc(sizeof(sweep_freqs));
    while (n[l0][l0] == 0) {

        result->length = 0;
        if (result->sweep_frequencies != NULL) {
            free(result->sweep_frequencies);
        }
        result->sweep_frequencies = NULL;
        for (unsigned int i = 0; i < L; i++) {
            for (unsigned int k = 0; k < L; k++) {
                n[i][k] = 0;
            }
        }
        n[l0][l0] = 1;
        result->sweep_frequencies = (unsigned int*) malloc(tfinal * L * L * sizeof(unsigned int));

        for (t = 0; t < tfinal; t++) {
            // Record the status of the population and check for fixation:
            ntot = 0;
            for (unsigned int i = 0; i < L; i++) {
                for (unsigned int k = 0; k < L; k++) {
                    ntot += n[i][k];
                    result->sweep_frequencies[t*L*L + i*L + k] = n[i][k];
                }
            }

            // Stop if one of the alleles is fixed or extinct:
            if ((ntot == 0) || (ntot == N * L * L)) {
                result->length = (t * L * L) + 1;
                break;
            }

            // Evolve for one generation
            next_gen_2d(L, n, s, mig, N);
        }
    }

    if (t == tfinal) {
		fprintf(stderr, "Error: Simulation finished without fixation.\n");
		return NULL;
	}

    return result;
}

void free_freqs_struct(sweep_freqs* freqs) {
    if (freqs != NULL) {
        if (freqs->sweep_frequencies != NULL) {
            free(freqs->sweep_frequencies);
        }
        free(freqs);
    }
}
