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

#ifndef SPECTRAL_OPTS_H
#define SPECTRAL_OPTS_H

    #define MEAN 1
    #define MEDIAN 2

    /** @struct spectralOpts
     *  Struct with the options for the spectral estimators.
     *
     *  @var spectralOpts::fs
     *    Samplig frequency.
     *  @var spectralOpts::window
     *    Window type. See windows.h for the available types.
     *  @var spectralOpts::nperseg
     *    Length of each segment.
     *  @var spectralOpts::noverlap
     *    Number of points to overlap between segments.
     *  @var spectralOpts::nfft
     *    Length of the FFT used.
     *  @var spectralOpts::nframes
     *    Number of frames. It is automatically computed if new_spectral_opts_basic() is used.
     *  @var spectralOpts::isComplex
     *    Real or complex input sequence?
     *  @var spectralOpts::sides
     *    Sides of the FFT for real data (1 or 2). For complex data this parameter is not taken into account.
     *  @var spectralOpts::scaling
     *    Compute the PSD (V^2/Hz) or the spectrum (V^2) if data is measured in V and fs in Hz.
     *  @var spectralOpts::average
     *    Method to use when averaging periodograms (mean or median).
     */
    typedef struct spectralOpts {
        float fs;
        int window;
        int nperseg;
        int noverlap;
        int nfft;
        int nframes;
        bool isComplex;
        int sides;
        int scaling;
        int average;
    } spectralOpts;

    int compute_num_frames(int N, int nperseg, int noverlap);
    spectralOpts new_spectral_opts_basic(int N, int window, int nperseg, int scaling);
    void check_spectral_opts(spectralOpts spectralOpts);
#endif
