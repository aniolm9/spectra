#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <kissfft/kiss_fft.h>
#include <gsl/gsl_cdf.h>
#include "opts.h"
#include "windows.h"
#include "spectral_common.h"
#include "spectrum.h"

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
void welch(double *data, double *freqs, double *power, int N, opts *welchOpts) {
    /* TODO: check welch options */
    //check_welch_opts(welchOpts);
    /* Declare some useful variables from the welchOpts struct */
    int nperseg = welchOpts->nperseg;
    int num_frames = compute_num_frames(N, welchOpts);
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
            break;
        
        default:
            fprintf(stderr, "averaging method not supported, falling back to mean method\n");
        
        case MEAN:
            for (int k = 0; k < num_frames; k++) {
                for (int j = 0; j < nperseg; j++) {
                    power[j] += psd[k][j] / num_frames;
                }
            }
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
void periodogram(double *data, double *freqs, double *power, int N, opts *welchOpts) {
    /* TODO: check welch options */
    //check_welch_opts(welchOpts);
    /* The periodogram is like the Welch method but for noverlap = 0 */
    opts periodogramOpts = *welchOpts;
    periodogramOpts.noverlap = 0;
    welch(data, freqs, power, N, &periodogramOpts);
}
