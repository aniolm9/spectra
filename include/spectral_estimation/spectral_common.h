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

#ifndef SPECTRAL_COMMON_H
#define SPECTRAL_COMMON_H
    void fftfreq(double *freqs, int M, float ts);
    void spectral_helper(double *data_x, double *data_y, double *freqs, double **psd, int N_x, int N_y, spectralOpts *spectralOpts);
#endif
