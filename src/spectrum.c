#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <kissfft/kiss_fft.h>
#include "spectrum.h"
#include "opts.h"
#include "windows.h"


static void spectral_helper(double *data_x, double *data_y, int N_x, int N_y, opts *spectralOpts) {
    /* We do not support cross power spectral density estimations yet */
    if (data_x != data_y) {
        fprintf(stderr, "different input arrays not supported\n");
    }
}

void new_opts() {

}

void IQ2fftcpx(double *iq, kiss_fft_cpx *cpx, int N) {
    int k = 0;
    for (int j = 0; j < N; j+=2) {
        cpx[k].r = iq[j];
        cpx[k].i = iq[j+1];
        k++;
    }
}

void check_welch_opts(opts *welchOpts) {
    if (!welchOpts) {

    }
}


typedef struct {
    int32_t value;
    uint16_t index;
} akaike;

int cmp(const void *a, const void *b) {
    akaike *a1 = (akaike *)a;
    akaike *a2 = (akaike *)b;
    if ((*a1).value > (*a2).value)
        return -1;
    else if ((*a1).value < (*a2).value)
        return 1;
    else
        return 0;
}

double noise_power_aic(const double *psd, int N, opts *welchOpts) {
    int nperseg = welchOpts->nperseg;
    int noverlap = welchOpts->noverlap;
    int nstep = nperseg - noverlap;
    N = N/2;
    int num_frames = N / nstep + (N % nstep != 0) - 1;
    akaike *aic = calloc(nperseg, sizeof(*aic));
    /* Compute Akaike Information Criterion (AIC) */
    float powSum;
    float powProduct;
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
    qsort(aic, nperseg, sizeof(akaike), cmp);
    int kmin = aic[nperseg-1].index;
    double noise_power = 0.0;
    for (int i = kmin; i < nperseg; i++) {
        noise_power += psd[i]/(nperseg-kmin);
    }
    return noise_power;
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
 */
void welch(double *data, double *freqs, double *power, int N, opts *welchOpts) {
    //check_welch_opts(welchOpts);
    /* Part of the code below should be moved to spectral_helper */
    int nperseg = welchOpts->nperseg;
    int noverlap = welchOpts->noverlap;
    int nstep = nperseg - noverlap;
    kiss_fft_cfg cfg = kiss_fft_alloc(welchOpts->nfft, 0, 0, 0);
    kiss_fft_cpx *in = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * welchOpts->nfft);
    kiss_fft_cpx *out = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * welchOpts->nfft);
    kiss_fft_cpx *data_cpx = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * N/2);
    /* TODO: check if data is complex or not */
    IQ2fftcpx(data, data_cpx, N);
    N = N/2;
    int num_frames = N / nstep + (N % nstep != 0) - 1;
    printf("nperseg: %d, noverlap: %d, nstep: %d, nframes: %d\n", nperseg, noverlap, nstep, num_frames);
    kiss_fft_cpx *frame = malloc(nperseg * sizeof(kiss_fft_cpx));
    kiss_fft_cpx *windowed_frame = malloc(nperseg * sizeof(kiss_fft_cpx));
    double powVal;
    int k;
    for (k = 0; k < num_frames; k++) {
        memcpy(frame, data_cpx+k*nstep, nperseg*sizeof(kiss_fft_cpx));
        /* Apply a window
         * TODO: check non-rectangular windows.
        **/
        windowing(frame, windowed_frame, nperseg, welchOpts->window);
        for (int j = 0; j < nperseg; j++) {
            printf("%lf %lf\n", windowed_frame[j].r, windowed_frame[j].i);
        }
        memcpy(in, windowed_frame, nperseg*sizeof(kiss_fft_cpx));
        kiss_fft(cfg, in, out);

        /* Compute PSD */
        for (int j = 0; j < nperseg; j++) {
            powVal = (out[j].r*out[j].r + out[j].i*out[j].i) / nperseg;
            power[j] += powVal / num_frames;
        }
    }
    /* Last frame may be smaller if it is not padded, thus we ignore it
     * TODO: implement a padding to avoid ignoring it.
    **/
    printf("memcpy ok\n");
    for (int j = 0; j < nperseg; j++) {
        //printf("%d\t %lf\n", j, power[j]);
    }

    free(in);
    free(out);
    free(frame);
    free(windowed_frame);
    free(data_cpx);
    kiss_fft_free(cfg);
}

int main() {
    /* Open file */
    const char *filename = "sdr_iq_data_0G433000_60dB-20210505_184759.iqdat";
    FILE *fp = fopen(filename, "rb");
    /* Get filesize and allocate memory */
    fseek(fp, 0, SEEK_END);
    size_t size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    printf("FILE SIZE: %lu\n", (unsigned long)size);
    int16_t *data = malloc(size);
    int samples = size / 2;
    /* Read bytes */
    fread(data, size, 1, fp);
    double *datad = malloc(samples*sizeof(double));
    for (int i = 0; i < samples; i++) {
        datad[i] = (double)data[i];
    }
    int nperseg = 2048;
    double *freqs;
    double *power = calloc(nperseg, sizeof(*power));
    opts wopts = {1.0, 0, nperseg, 1024, nperseg, true, 1, 1, 1};
    welch(datad, freqs, power, samples, &wopts);
    double noise_power = noise_power_aic(power, samples, &wopts);
    printf("Estimated noise power = %lf\n", noise_power);
    /* Free memory */
    fclose(fp);
    free(data);
    free(power);
}
