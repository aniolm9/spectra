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
#include "spectral_common.h"
#include "windowing.h"

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
void spectral_helper(double *data_x, double *data_y, double *freqs, double **psd, int N_x, int N_y, spectralOpts *spectralOpts) {
    /* TODO: We do not support cross power spectral density estimations yet */
    if (data_x != data_y) {
        fprintf(stderr, "different input arrays not supported\n");
        exit(EXIT_FAILURE);
    } else {
        int nperseg = spectralOpts->nperseg;
        int noverlap = spectralOpts->noverlap;
        int nstep = nperseg - noverlap;
        kiss_fft_cfg cfg = kiss_fft_alloc(spectralOpts->nfft, 0, NULL, NULL);
        kiss_fft_cpx *in = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * spectralOpts->nfft);
        kiss_fft_cpx *out = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * spectralOpts->nfft);
        kiss_fft_cpx *data_cpx = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * N_x/2);
        /* TODO: check if data is complex or not */
        IQ2fftcpx(data_x, data_cpx, N_x);
        int num_frames = spectralOpts->nframes;
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
