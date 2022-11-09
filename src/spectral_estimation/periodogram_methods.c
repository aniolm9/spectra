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
#include "spectral_common.h"
#include "periodogram_methods.h"

/**
 * Estimate power spectral density using Welch's method. The method consists in splitting
 * the data in overlapping segments, computing the modified periodogram for each segment
 * and finally averaging the periodograms.
 * 
 * @param data Input array with the data samples (as doubles).
 * @param freqs Output array containing the sample frequencies.
 * @param power Output array with the power spectral density of data.
 * @param N Number of samples in data.
 * @param welchOpts Struct containing the internal settings of the estimators.
 */
void welch(double *data, double *freqs, double *power, int N, spectralOpts *welchOpts) {
    /* Check the options passed by the user */
    check_spectral_opts(*welchOpts);
    /* Declare some useful variables from the welchOpts struct */
    int nperseg = welchOpts->nperseg;
    int num_frames = welchOpts->nframes;
    int average = welchOpts->average;
    /* Allocate a 2-D array for the PSD */
    double **psd = (double **) malloc(num_frames * sizeof(double*));
    for(int i = 0; i < num_frames; i++) {
        psd[i] = (double *)malloc(nperseg * sizeof(double));
    }
    /* Use the spectral helper function to compute the PSD */
    spectral_helper(data, data, freqs, psd, N, N, welchOpts);
    /* Decide which averaging method we should use */
    switch (average) {
        case MEDIAN:
            fprintf(stderr, "median averaging method not supported\n");
            exit(EXIT_FAILURE);
            break;

        case MEAN:
            for (int k = 0; k < num_frames; k++) {
                for (int j = 0; j < nperseg; j++) {
                    power[j] += psd[k][j] / num_frames;
                }
            }
            break;

        default:
            fprintf(stderr, "averaging method not supported\n");
            break;
    }
    /* Free the PSD 2-D array */
    for(int i = 0; i < num_frames; i++) {
        free(psd[i]);
    }
    free(psd);
}

/**
 * Estimate power spectral density using a periodogram. This implementation
 * performs a temporal average. This could be avoided by setting nperseg=N.
 *
 * @param data Input array with the data samples (as doubles).
 * @param freqs Output array containing the sample frequencies.
 * @param power Output array with the power spectral density of data.
 * @param N Number of samples in data.
 * @param welchOpts Struct containing the internal settings of the estimators.
 */
void periodogram(double *data, double *freqs, double *power, int N, spectralOpts *welchOpts) {
    /* The periodogram is like the Welch method but for noverlap = 0 */
    spectralOpts periodogramOpts = *welchOpts;
    periodogramOpts.noverlap = 0;
    welch(data, freqs, power, N, &periodogramOpts);
}
