#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <complex.h>
#include <fftw3.h>
#include "spectrum.h"

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

void welch(double *data, double *output, float fs, int window, int nperseg, int noverlap, int nfft, int sides) {
    fftw_plan plan;
    fftw_complex *in, *out;
    fft_initializer(in, out, plan, nfft, sides);
    //...
    fft_cleaner(in, out, plan);
}

int main() {
    double complex z = 5.0;
    printf(
        "z = %.1f% + .1fi\n",
        creal(z), cimag(z));
    if (cimag(z) == 0.0) {
        printf("Real number\n");
    }
}
