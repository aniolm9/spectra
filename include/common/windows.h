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

#ifndef WINDOWS_H
#define WINDOWS_H
    /* Types of window */
    #define RECTANGULAR 0
    #define HANN 1
    #define HAMMING 2
    #define BLACKMAN 3

    #define MAX_WINDOW 3

    /* Scaling methods */
    #define DENSITY 1
    #define SPECTRUM 2

    kiss_fft_scalar hann(int n, int M);
    kiss_fft_scalar hamming(int n, int M);
    kiss_fft_scalar blackman(int n, int M);
    kiss_fft_scalar windowing(const kiss_fft_cpx *data, kiss_fft_cpx *output, int M, int window, float fs, int scaling);
#endif
