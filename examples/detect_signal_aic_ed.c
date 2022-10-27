/*
 * Copyright (C) 2022 Aniol Martí
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted under the
 * terms of the BSD 3-clause license.
 *
 * The software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * BSD 3-clause license for more details.
 *
 * You should have received a copy of the BSD 3-clause license
 * along with the software. If not, see <https://opensource.org/licenses/BSD-3-Clause>.
 */

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <kissfft/kiss_fft.h>
#include "spectral_opts.h"
#include "windows.h"
#include "periodogram_methods.h"
#include "noise_power_est.h"
#include "signal_detection.h"

/**
 * This example opens a file passed as an argument, computes its PSD
 * using Welch's method, estimates the noise power using the AIC and
 * finally it detects signal presence with an energy detector.
 */
int main(int argc, char **argv) {
    /* Open file */
    char *filename;
    if (argc >= 2) {
        filename = argv[1];
    } else {
        fprintf(stderr, "error opening file\n");
        exit(EXIT_FAILURE);
    }
    FILE *fp = fopen(filename, "rb");
    /* Get filesize and allocate memory */
    fseek(fp, 0, SEEK_END);
    size_t size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    printf("FILE SIZE: %lu\n", (unsigned long)size);
    int16_t *data = malloc(size);
    int samples = size / 2;
    /* Read bytes */
    if (fread(data, size, 1, fp)*size != size) {
        fprintf(stderr, "error reading file\n");
        exit(EXIT_FAILURE);
    }
    double *datad = malloc(samples*sizeof(double));
    for (int i = 0; i < samples; i++) {
        datad[i] = (double)data[i];
    }
    int nperseg = 2048;
    double *freqs = calloc(nperseg, sizeof(*freqs));
    double *power = calloc(nperseg, sizeof(*power));
    bool *signal_presence = malloc(samples/2 * sizeof(bool));
    spectralOpts wopts = new_spectral_opts_basic(samples, HANN, nperseg, DENSITY);
    welch(datad, freqs, power, samples, &wopts);
    double noise_power = noise_power_aic(power, samples, &wopts);
    printf("Estimated noise power = %lf\n", noise_power);
    energy_detector(datad, signal_presence, samples, noise_power, 0.05, 4096);
    /* Free memory */
    fclose(fp);
    free(data);
    free(freqs);
    free(power);
    free(signal_presence);
}
