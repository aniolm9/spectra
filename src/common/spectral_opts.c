#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <kissfft/kiss_fft.h>
#include "spectral_opts.h"
#include "tools.h"
#include "windows.h"

spectralOpts new_spectral_opts_basic(int window, int nperseg, int scaling) {
    /* Some comments about the input parameters:
     * nperseg = N => We must take nperseg=2^k < N
    **/
    spectralOpts opts = {(float)1.0,
                         window,
                         nperseg,
                         nperseg/2,
                         nperseg,
                         true,
                         1,
                         scaling,
                         MEAN};
    check_spectral_opts(opts);
    return opts;
}

void check_spectral_opts(spectralOpts spectralOpts) {
    if (spectralOpts.window < 0 || spectralOpts.window > MAX_WINDOW) {
        fprintf(stderr, "invalid window\n");
        exit(EXIT_FAILURE);
    }
    if (spectralOpts.nfft < spectralOpts.nperseg) {
        fprintf(stderr, "nfft can't be smaller than nperseg\n");
        exit(EXIT_FAILURE);
    }
    if (spectralOpts.noverlap >= spectralOpts.nperseg) {
        fprintf(stderr, "noverlap must be smaller than nperseg\n");
        exit(EXIT_FAILURE);
    }
    if (spectralOpts.scaling != SPECTRUM && spectralOpts.scaling != DENSITY) {
        fprintf(stderr, "scaling method not available\n");
        exit(EXIT_FAILURE);
    }
    if (spectralOpts.average != MEAN) {
        fprintf(stderr, "averaging method not supported\n");
        exit(EXIT_FAILURE);
    }
}
