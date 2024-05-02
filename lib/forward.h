#ifndef __FORWARD_H__
#define __FORWARD_H__

//Contains the sweep data.
//The list has length L * t, where t is the sweep time and L is the number of demes.
//To access the number of mutants at time i in deme k, query a[i * k].
//Both sweep timer and the deme count start at 0.
typedef struct sweep_frequencies {
    unsigned int *sweep_frequencies;
    unsigned int *pop_id_records;
    int length;
} sweep_freqs;

//Create the forward simulation data structures.
sweep_freqs *forward_1d(unsigned int L, unsigned int N, double s, double mig, unsigned int tfinal, long SEED);
sweep_freqs *forward_2d (unsigned int L, unsigned int N, double s, double mig, unsigned int tfinal, unsigned int l0, long SEED);

//Memory management.
void free_freqs_struct(sweep_freqs* freqs);

#endif
