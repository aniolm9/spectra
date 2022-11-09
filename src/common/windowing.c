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
#include <math.h>
#include <kissfft/kiss_fft.h>

#include "windowing.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

/**
 * Hann window.
 *
 * @param n Index of current sample.
 * @param M Window size.
 * @return The value of the window of size M at position n.
 */
kiss_fft_scalar hann(int n, int M) {
    return 0.5 - 0.5*cos(2*M_PI*n/(M-1));
}

/**
 * Hamming window.
 *
 * @param n Index of current sample.
 * @param M Window size.
 * @return The value of the window of size M at position n.
 */
kiss_fft_scalar hamming(int n, int M) {
    return 0.54 - 0.46*cos(2*M_PI*n/(M-1));
}

/**
 * Blackman window.
 *
 * @param n Index of current sample.
 * @param M Window size.
 * @return The value of the window of size M at position n.
 */
kiss_fft_scalar blackman(int n, int M) {
    return 0.42 - 0.5*cos(2*M_PI*n/M) + 0.08*cos(4*M_PI*n/M);
}

/**
 * Apply a supported window to an array of samples.
 *
 * @param data Array of input samples as kiss_fft_cpx scalars.
 * @param output Array of output samples after applying the window.
 * @param M Size of the window.
 * @param window Type of the window.
 * @param fs Sampling frequency.
 * @param scaling Type of scaling: PSD or spectrum.
 * @return The scaling factor.
 */
kiss_fft_scalar windowing(const kiss_fft_cpx *data, kiss_fft_cpx *output, int M, int window, float fs, int scaling) {
    kiss_fft_scalar scale_density = 0.0;
    kiss_fft_scalar scale_spectrum = 0.0;
    for (int n = 0; n < M; n++) {
        switch (window) {
            case RECTANGULAR:
                output[n] = data[n];
                scale_density += 1;
                scale_spectrum += 1;
                break;
            case HANN:
                output[n].r = data[n].r * hann(n, M);
                output[n].i = data[n].i * hann(n, M);
                scale_density += hann(n, M) * hann(n, M);
                scale_spectrum += hann(n, M);
                break;
            case HAMMING:
                output[n].r = data[n].r * hamming(n, M);
                output[n].i = data[n].i * hamming(n, M);
                scale_density += hamming(n, M) * hamming(n, M);
                scale_spectrum += hamming(n, M);
                break;
            case BLACKMAN:
                output[n].r = data[n].r * blackman(n, M);
                output[n].i = data[n].i * blackman(n, M);
                scale_density += blackman(n, M) * blackman(n, M);
                scale_spectrum += blackman(n, M);
                break;
            default:
                fprintf(stderr, "invalid window\n");
                exit(EXIT_FAILURE);
        }
    }
    if (scaling == DENSITY) {
        return 1.0 / (fs * scale_density);
    } else if (scaling == SPECTRUM) {
        return 1.0 / (scale_spectrum * scale_spectrum);
    } else {
        fprintf(stderr, "invalid scaling method\n");
        exit(EXIT_FAILURE);
    }
}
