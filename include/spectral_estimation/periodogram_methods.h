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

#ifndef PERIODOGRAM_METHODS_H
#define PERIODOGRAM_METHODS_H
    void periodogram(double *data, double *freqs, double *power, int N, spectralOpts *welchOpts);
    void welch(double *data, double *freqs, double *power, int N, spectralOpts *welchOpts);
#endif
