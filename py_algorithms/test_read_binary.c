#include <stdio.h>
#include <stdlib.h>

int main() {
    FILE *file;
    const char *filename = "sweep_trajectories.bin";
    double demes_dbl, num_steps_dbl;
    double **allele_frequency;
    double *time;
    size_t demes, num_steps;

    // Open the binary file for reading
    file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Read the number of demes and number of steps as doubles
    fread(&demes_dbl, sizeof(double), 1, file);
    fread(&num_steps_dbl, sizeof(double), 1, file);

    demes = (size_t)demes_dbl;
    num_steps = (size_t)num_steps_dbl;

    // Allocate memory for time
    time = (double *)malloc(sizeof(double) * num_steps);

    // Allocate memory for allele frequencies
    allele_frequency = (double **)malloc(sizeof(double *) * num_steps);
    for (size_t step = 0; step < num_steps; step++) {
        allele_frequency[step] = (double *)malloc(sizeof(double) * demes);
    }

    // Read time and allele frequencies
    for (size_t step = 0; step < num_steps; step++) {
        fread(&time[step], sizeof(double), 1, file); // Read time
        for (size_t i = 0; i < demes; i++) {
            fread(&allele_frequency[step][i], sizeof(double), 1, file); // Read frequency for each deme
        }
    }

    // New Part: Read the last time element
    double last_time_element;
    fread(&last_time_element, sizeof(double), 1, file);


    // Close the file
    fclose(file);

    // Print the read data for verification
    printf("Demes: %zu, Timesteps: %zu\n", demes, num_steps);
    for (size_t step = 0; step < num_steps; step++) {
        printf("Time[%zu]: %f\n", step, time[step]);
        for (size_t i = 0; i < demes; i++) {
            printf("Allele Frequency[%zu][%zu]: %f\n", step, i, allele_frequency[step][i]);
        }
    }
    printf("Last Time Element: %f\n", last_time_element);
    

    // Free allocated memory
    for (size_t step = 0; step < num_steps; step++) {
        free(allele_frequency[step]);
    }
    free(allele_frequency);
    free(time);

    return 0;
}
