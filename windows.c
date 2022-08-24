#include <math.h>
#include <stdio.h>
#include "windows.h"

double hann(int n, int M) {
    return 0.5 - 0.5*cos(2*M_PI*n/(M-1));
}

double hamming(int n, int M) {
    return 0.54 - 0.46*cos(2*M_PI*n/(M-1));
}

double blackman(int n, int M) {
    return 0.42 - 0.5*cos(2*M_PI*n/M) + 0.08*cos(4*M_PI*n/M);
}

void windowing(int M, const double *data, int window, double *output) {
    for (int n = 0; n < M; n++) {
        switch (window) {
            case RECTANGULAR:
                output[n] = data[n];
                break;
            case HANN:
                output[n] = data[n] * hann(n, M);
                break;
            case HAMMING:
                output[n] = data[n] * hamming(n, M);
                break;
            case BLACKMAN:
                output[n] = data[n] * blackman(n, M);
                break;
            default:
                fprintf(stderr, "invalid window\n");
                break;
        }
    }
}
