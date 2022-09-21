#define MEAN 1
#define MEDIAN 2

#include "opts.h"

typedef struct {
    int32_t value;
    uint16_t index;
} akaike;

void welch(double *data, double *freqs, double *power, int N, opts *welchOpts);
double noise_power_aic(const double *psd, int N, opts *welchOpts);
void energy_detector(double *data, bool *signal_presence, int N, double noise_power, float Pfa, int df);
