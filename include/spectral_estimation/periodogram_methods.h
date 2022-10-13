#ifndef PERIODOGRAM_METHODS_H
#define PERIODOGRAM_METHODS_H
    void periodogram(double *data, double *freqs, double *power, int N, spectralOpts *welchOpts);
    void welch(double *data, double *freqs, double *power, int N, spectralOpts *welchOpts);
#endif
