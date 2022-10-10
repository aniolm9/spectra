#include <stdbool.h>
#include <math.h>
#include <kissfft/kiss_fft.h>
#include <gsl/gsl_cdf.h>
#include "opts.h"
#include "spectrum.h"
#include "signal_detection.h"

/**
 * Detects signal presence based on a classical energy detector.
 *
 * @param data Input array with the data samples (as doubles).
 * @param signal_presense Output array with "true" in the positions where a signal is detected.
 * @param N Number of samples in data.
 * @param noise_power Power of the noise.
 * @param Pfa Probability of false alarm.
 * @param df Degrees of freedom.
 */
void energy_detector(double *data, bool *signal_presence, int N, double noise_power, float Pfa, int df) {
    kiss_fft_cpx *data_cpx = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * N/2);
    /* TODO: check if data is complex or not */
    IQ2fftcpx(data, data_cpx, N);
    N = N/2;
    double gamma = gsl_cdf_chisq_Qinv(Pfa, df) * noise_power;
    double energy;
    for (int j = 0; j < N; j++) {
        energy = 0.0;
        for (int k = j; k < j+df; k++) {
            energy += data_cpx[k].r*data_cpx[k].r + data_cpx[k].i*data_cpx[k].i;
        }
        signal_presence[j] = energy > gamma;
    }
    /* Free memory */
    free(data_cpx);
}