#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <kissfft/kiss_fft.h>
#include <gsl/gsl_cdf.h>
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
static void IQ2fftcpx(double *iq, kiss_fft_cpx *cpx, int N) {
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
static int compute_num_frames(int N, opts *spectralOpts) {
    /* TODO: support padding. The -1 will go away if we pad */
    int nperseg = spectralOpts->nperseg;
    int noverlap = spectralOpts->noverlap;
    int nstep = nperseg - noverlap;
    N = N/2;
    int num_frames = N / nstep + (N % nstep != 0) - 1;
    printf("%d\n", num_frames);
    return num_frames;
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
 * Compare two akaike structs.
 *
 * @param a First akaike struct.
 * @param b Second akaike struct.
 * @return -1 if the value of a is greater than that of b, 1 if smaller, 0 otherwise.
 */
static int cmp_akaike(const void *a, const void *b) {
    akaike *a1 = (akaike *)a;
    akaike *a2 = (akaike *)b;
    if ((*a1).value > (*a2).value)
        return -1;
    else if ((*a1).value < (*a2).value)
        return 1;
    else
        return 0;
}

/**
 * Estimate the noise power of a signal using the Akaike Information Criterion (AIC).
 *
 * @param psd Output matrix where each row is a frame and each column a sample.
 * @param N Number of samples in data.
 * @param welchOpts Struct containing the internal settings of the estimators.
 * @return A double with the noise power estimate.
 */
double noise_power_aic(const double *psd, int N, opts *welchOpts) {
    int nperseg = welchOpts->nperseg;
    int num_frames = compute_num_frames(N, welchOpts);
    akaike *aic = calloc(nperseg, sizeof(*aic));
    /* Compute Akaike Information Criterion (AIC) */
    double powSum;
    double powProduct;
    double aic_n;
    for (int k = 0; k < nperseg; k++) {
        powSum = 0.0;
        powProduct = 1.0;
        for (int i=k; i < nperseg; i++) {
            powSum += psd[i] / (nperseg-k);
            powProduct *= pow(psd[i], 1.0/(nperseg-k));
        }
        aic_n = (nperseg-k)*num_frames*log(powSum/powProduct)+k*(2*nperseg-k);
        aic[k].index = k;
        aic[k].value = (int32_t)aic_n;
    }
    /* Sort the AIC array and find the index of the minimum value */
    qsort(aic, nperseg, sizeof(akaike), cmp_akaike);
    int kmin = aic[nperseg-1].index;
    double noise_power = 0.0;
    for (int i = kmin; i < nperseg; i++) {
        noise_power += psd[i]/(nperseg-kmin);
    }
    return noise_power;
}

/**
 * Detects signal presence based on a classical energy detector.
 *
 * @param data Input array with the data samples (as doubles).
 * @param signal_presense Output array with "true" in the positions where a signal is detected.
 * @param N Number of samples in data.
 * @param noise_power Power of the noise.
 * @param Pfa Probability of false alarm.
 * @param df Degrees of freedom.
 */
void energy_detector(double *data, bool *signal_presence, int N, double noise_power, float Pfa, int df) {
    kiss_fft_cpx *data_cpx = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * N/2);
    /* TODO: check if data is complex or not */
    IQ2fftcpx(data, data_cpx, N);
    N = N/2;
    double gamma = gsl_cdf_chisq_Qinv(Pfa, df) * noise_power;
    double energy;
    for (int j = 0; j < N; j++) {
        energy = 0.0;
        for (int k = j; k < j+df; k++) {
            energy += data_cpx[k].r*data_cpx[k].r + data_cpx[k].i*data_cpx[k].i;
        }
        signal_presence[j] = energy > gamma;
    }
    /* Free memory */
    free(data_cpx);
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
