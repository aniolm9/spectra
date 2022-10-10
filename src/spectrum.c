#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <kissfft/kiss_fft.h>
#include <gsl/gsl_cdf.h>
#include "opts.h"
#include "spectrum.h"
#include "windows.h"

void new_opts() {
    /* Some comments about the input parameters:
     * nfft < nperseg => error
     * noverlap >= nperseg => error
     * nperseg = N => We must take nperseg=2^k < N
    **/
}

static void check_welch_opts(opts *welchOpts) {
    if (!welchOpts) {

    }
}

/**
 * Transform an array of doubles (in IQ format) into a
 * kiss_fft_cpx array.
 *
 * @param iq Input array of IQ samples (as doubles).
 * @param cpx Output array as kiss_fft_cpx scalars.
 * @param N Number of samples in the input array.
 */
void IQ2fftcpx(double *iq, kiss_fft_cpx *cpx, int N) {
    int k = 0;
    for (int j = 0; j < N; j+=2) {
        cpx[k].r = iq[j];
        cpx[k].i = iq[j+1];
        k++;
    }
}

/**
 * Compute the number of frames in which the data array is divided. It uses
 * the total number of samples, the frame length and the number of overlapping samples.
 *
 * @param N Number of samples in data.
 * @param welchOpts Struct containing the internal settings of the estimators.
 * @return The number of frames in which the data array is divided.
 */
int compute_num_frames(int N, opts *spectralOpts) {
    /* TODO: support padding. The -1 will go away if we pad */
    int nperseg = spectralOpts->nperseg;
    int noverlap = spectralOpts->noverlap;
    int nstep = nperseg - noverlap;
    N = N/2;
    int num_frames = N / nstep + (N % nstep != 0) - 1;
    return num_frames;
}

/**
 * Compute the Discrete Fourier Transform (DFT) sample frequencies.
 *
 * @param freqs Output array containing the sample frequencies.
 * @param M Window length (in general, nfft or nperseg).
 * @param ts Inverse of the sampling frequency.
 */
void fftfreq(double *freqs, int M, float ts) {
    int i;
    int cut_index = (int)((M-1)/2) + 1;
    for (i = 0; i < cut_index; i++) {
        freqs[i] = i/(ts*M);
    }
    freqs[cut_index] = (-M/2)/(ts*M);
    for (i = cut_index + 1; i < M; i++) {
        freqs[i] = freqs[i-1] + 1/(ts*M);
    }
}

/**
 * Helper function that can be used in different spectral estimation tools, such as
 * the periodogram or the averaged periodogram.
 *
 * @param data_x First input data array (as doubles).
 * @param data_y Second input data array (as doubles).
 * @param freqs Output array containing the sample frequencies.
 * @param psd Output matrix where each row is a frame and each column a sample.
 * @param N_x Number of samples in data_x.
 * @param N_y Number of samples in data_y.
 * @param spectralOpts Struct containing the internal settings of the estimators.
 */
static void spectral_helper(double *data_x, double *data_y, double *freqs, double **psd, int N_x, int N_y, opts *spectralOpts) {
    /* TODO: We do not support cross power spectral density estimations yet */
    if (data_x != data_y) {
        fprintf(stderr, "different input arrays not supported\n");
    } else {
        int nperseg = spectralOpts->nperseg;
        int noverlap = spectralOpts->noverlap;
        int nstep = nperseg - noverlap;
        kiss_fft_cfg cfg = kiss_fft_alloc(spectralOpts->nfft, 0, 0, 0);
        kiss_fft_cpx *in = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * spectralOpts->nfft);
        kiss_fft_cpx *out = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * spectralOpts->nfft);
        kiss_fft_cpx *data_cpx = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * N_x/2);
        /* TODO: check if data is complex or not */
        IQ2fftcpx(data_x, data_cpx, N_x);
        int N = N_x/2;
        int num_frames = N / nstep + (N % nstep != 0) - 1 - 1;
        kiss_fft_cpx *frame = malloc(nperseg * sizeof(kiss_fft_cpx));
        kiss_fft_cpx *windowed_frame = malloc(nperseg * sizeof(kiss_fft_cpx));
        double powVal;
        int k;
        for (k = 0; k < num_frames; k++) {
            /* Create a frame of size nperseg with nstep new samples */
            memcpy(frame, data_cpx+k*nstep, nperseg*sizeof(kiss_fft_cpx));
            /* Apply a window */
            kiss_fft_scalar scale = windowing(frame, windowed_frame, nperseg, spectralOpts->window, spectralOpts->fs, spectralOpts->scaling);
            /* Copy the windowed frame into a new array and compute the FFT */
            memcpy(in, windowed_frame, nperseg*sizeof(kiss_fft_cpx));
            kiss_fft(cfg, in, out);
            /* Compute PSD */
            for (int j = 0; j < nperseg; j++) {
                powVal = (out[j].r*out[j].r + out[j].i*out[j].i) * scale;
                psd[k][j] = powVal;
            }
        }
        /* Compute the frequency values */
        fftfreq(freqs, spectralOpts->nfft, 1/spectralOpts->fs);
        /* Last frame may be smaller if it is not padded, thus we ignore it
        * TODO: implement a padding to avoid ignoring it.
        **/
        /* Free memory */
        free(in);
        free(out);
        free(frame);
        free(windowed_frame);
        free(data_cpx);
        kiss_fft_free(cfg);
    }
}

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
