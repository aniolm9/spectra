#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>
#include "spectrum.h"
#include "opts.h"

static void fft_initializer(fftw_complex *in, fftw_complex *out, fftw_plan plan, int nfft, int sides) {
    /* TODO
     * in complex => sides = 2
     * in real => sides = 1 | sides = 2
     * If data are real but sides == 2, we proceed as in the complex case.
     */
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nfft);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nfft);
    plan = fftw_plan_dft_1d(nfft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
}

static void fft_cleaner(fftw_complex *in, fftw_complex *out, fftw_plan plan) {
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
}

static void spectral_helper(double *data_x, double *data_y, int N_x, int N_y, opts *spectralOpts) {
    /* We do not support cross power spectral density estimations yet */
    if (data_x != data_y) {
        fprintf(stderr, "different input arrays not supported\n");
    }
}

void new_opts() {

}

void check_welch_opts(opts *welchOpts) {
    if (!welchOpts) {

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
 */
void welch(double *data, double *freqs, double *power, int N, opts *welchOpts) {
    fftw_plan plan;
    fftw_complex *in, *out;
    //check_welch_opts(welchOpts);
    /* Part of the code below should be moved to spectral_helper */
    int nperseg = welchOpts->nperseg;
    int noverlap = welchOpts->noverlap;
    int nstep = nperseg - noverlap;
    printf("nperseg: %d, noverlap: %d, nstep: %d\n", nperseg, noverlap, nstep);
    //fft_initializer(in, out, plan, welchOpts->nfft, welchOpts->sides);
    printf("fft ok\n");
    int num_frames = N / nstep + (N % nstep != 0);
    double *frame = malloc(nperseg * sizeof(double));
    int i;
    for (i = 0; i < num_frames-1; i++) {
        memcpy(frame, data+i*nstep, nperseg*sizeof(double));
        for (int j = 0; j < nperseg; j++) {
            printf("%d\t%lf\t%lf\n", j, data[j], frame[j]);
        }
    }
    /* Last frame may be smaller if it is not padded, thus we ignore it
     * TODO: implement a padding to avoid ignoring it.
    **/
    printf("memcpy ok\n");
    
    /* */
    //fft_cleaner(in, out, plan); // There is a segfault here.
    //printf("fft clean ok\n");
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
    double *freqs;
    double *power;
    opts wopts = {1.0, 0, 2048, 1024, 2048, 1, 1, 1};
    welch(datad, freqs, power, samples, &wopts);
    /* Free memory */
    fclose(fp);
    free(data);
}
