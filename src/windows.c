#include <math.h>
#include <stdio.h>
#include <kissfft/kiss_fft.h>
#include "windows.h"

kiss_fft_scalar hann(int n, int M) {
    return 0.5 - 0.5*cos(2*M_PI*n/(M-1));
}

kiss_fft_scalar hamming(int n, int M) {
    return 0.54 - 0.46*cos(2*M_PI*n/(M-1));
}

kiss_fft_scalar blackman(int n, int M) {
    return 0.42 - 0.5*cos(2*M_PI*n/M) + 0.08*cos(4*M_PI*n/M);
}

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
                break;
        }
    }
    if (scaling == DENSITY) {
        return 1.0 / (fs * scale_density);
    } else if (scaling == SPECTRUM) {
        return 1.0 / (scale_spectrum * scale_spectrum);
    } else {
        fprintf(stderr, "invalid scaling method\n");
        return -1;
    }
}
