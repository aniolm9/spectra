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
#include "tools.h"
#include "windows.h"

/**
 * Create a new spectralOpts struct with the options passed
 * by the user. This function only accepts some basic tuning.
 *
 * @param window Type of the window.
 * @param nperseg Number of samples per segment.
 * @param scaling Type of scaling: PSD or spectrum.
 * @return A spectralOpts struct.
 */
spectralOpts new_spectral_opts_basic(int window, int nperseg, int scaling) {
    /* Some comments about the input parameters:
     * nperseg = N => We must take nperseg=2^k < N
     */
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

/**
 * Checks that the parameters in a spectralOpts struct are valid.
 *
 * @param spectralOpts The struct to be checked.
 */
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
