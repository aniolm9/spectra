#include <stdint.h>
#include <stdbool.h>
#include <kissfft/kiss_fft.h>
#include "spectral_opts.h"
#include "tools.h"
#include "noise_power_est.h"

/**
 * Estimate the noise power of a signal using the Akaike Information Criterion (AIC).
 *
 * @param psd Output matrix where each row is a frame and each column a sample.
 * @param N Number of samples in data.
 * @param welchOpts Struct containing the internal settings of the estimators.
 * @return A double with the noise power estimate.
 */
double noise_power_aic(const double *psd, int N, spectralOpts *welchOpts) {
    int nperseg = welchOpts->nperseg;
    int num_frames = compute_num_frames(N, welchOpts);
    akaike *aic = calloc(nperseg, sizeof(*aic));
    /* Compute Akaike Information Criterion (AIC) */
    double powSum;
    double powProduct;
    double aic_n;
    for (int k = 0; k < nperseg; k++) {
        powSum = 0.0;
        powProduct = 1.0;
        for (int i=k; i < nperseg; i++) {
            powSum += psd[i] / (nperseg-k);
            powProduct *= pow(psd[i], 1.0/(nperseg-k));
        }
        aic_n = (nperseg-k)*num_frames*log(powSum/powProduct)+k*(2*nperseg-k);
        aic[k].index = k;
        aic[k].value = (int32_t)aic_n;
    }
    /* Sort the AIC array and find the index of the minimum value */
    qsort(aic, nperseg, sizeof(akaike), cmp_akaike);
    int kmin = aic[nperseg-1].index;
    double noise_power = 0.0;
    for (int i = kmin; i < nperseg; i++) {
        noise_power += psd[i]/(nperseg-kmin);
    }
    return noise_power;
}
