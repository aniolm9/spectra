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

#ifndef TOOLS_H
#define TOOLS_H
    /** @struct akaike
     *  Struct for the pairs of AIC values and index.
     *
     *  @var akaike::value
     *    Value of the AIC for a given index.
     *  @var akaike::index
     *    The index of the AIC value.
     */
    typedef struct {
        int32_t value;
        uint16_t index;
    } akaike;

    void IQ2fftcpx(double *iq, kiss_fft_cpx *cpx, int N);
    int cmp_akaike(const void *a, const void *b);
#endif
