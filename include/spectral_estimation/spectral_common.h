#ifndef SPECTRAL_COMMON_H
#define SPECTRAL_COMMON_H
    void fftfreq(double *freqs, int M, float ts);
    void spectral_helper(double *data_x, double *data_y, double *freqs, double **psd, int N_x, int N_y, spectralOpts *spectralOpts);
#endif
