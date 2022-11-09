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

#include <stdint.h>
#include <stdbool.h>
#include <kissfft/kiss_fft.h>
#include <gsl/gsl_cdf.h>
#include "spectral_opts.h"
#include "tools.h"
#include "signal_detection.h"

/**
 * Detects signal presence based on a classical energy detector.
 *
 * @param data Input array with the data samples (as doubles).
 * @param signal_presence Output array with "true" in the positions where a signal is detected.
 * @param N Number of samples in data.
 * @param noise_power Power of the noise.
 * @param Pfa Probability of false alarm.
 * @param df Degrees of freedom.
 * @param nstep Step size between energy windows, between 1 (highest overlap) and df (no overlap).
 */
void energy_detector(double *data, bool *signal_presence, int N, double noise_power, float Pfa, int df, int nstep) {
    /* Check that 1 <= nstep <= df */
    if (nstep > df) {
        nstep = df;
    } else if (nstep < 1) {
        nstep = 1;
    }
    kiss_fft_cpx *data_cpx = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * N/2);
    /* TODO: check if data is complex or not */
    IQ2fftcpx(data, data_cpx, N);
    N = N/2;
    double gamma = gsl_cdf_chisq_Qinv(Pfa, df) * noise_power;
    double energy;
    int j, k;
    /* We decide signal presence in blocks of samples */
    for (j = 0; j < N-nstep; j += nstep) {
        int looplim = (j+df < N) ? (j+df) : (N);
        energy = 0.0;
        for (k = j; k < looplim; k++) {
            energy += data_cpx[k].r*data_cpx[k].r + data_cpx[k].i*data_cpx[k].i;
        }
        looplim = (j+nstep < N) ? (j+nstep) : (N);
        for (k = j; k < j+nstep; k++) {
            signal_presence[k] = energy > gamma;
        }
    }
    /* Detect the last samples */
    for (; j < N; j++) {
        int looplim = (j+df < N) ? (j+df) : (N);
        energy = 0.0;
        for (k = j; k < looplim; k++) {
            energy += data_cpx[k].r*data_cpx[k].r + data_cpx[k].i*data_cpx[k].i;
        }
        signal_presence[j] = energy > gamma;
    }
    /* Free memory */
    free(data_cpx);
}
