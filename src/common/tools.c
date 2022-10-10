#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <kissfft/kiss_fft.h>
#include "opts.h"
#include "tools.h"

/**
 * Transform an array of doubles (in IQ format) into a
 * kiss_fft_cpx array.
 *
 * @param iq Input array of IQ samples (as doubles).
 * @param cpx Output array as kiss_fft_cpx scalars.
 * @param N Number of samples in the input array.
 */
void IQ2fftcpx(double *iq, kiss_fft_cpx *cpx, int N) {
    int k = 0;
    for (int j = 0; j < N; j+=2) {
        cpx[k].r = iq[j];
        cpx[k].i = iq[j+1];
        k++;
    }
}

/**
 * Compute the number of frames in which the data array is divided. It uses
 * the total number of samples, the frame length and the number of overlapping samples.
 *
 * @param N Number of samples in data.
 * @param welchOpts Struct containing the internal settings of the estimators.
 * @return The number of frames in which the data array is divided.
 */
int compute_num_frames(int N, opts *spectralOpts) {
    /* TODO: support padding. The -1 will go away if we pad */
    int nperseg = spectralOpts->nperseg;
    int noverlap = spectralOpts->noverlap;
    int nstep = nperseg - noverlap;
    N = N/2;
    int num_frames = N / nstep + (N % nstep != 0) - 1;
    return num_frames;
}

/**
 * Compare two akaike structs.
 *
 * @param a First akaike struct.
 * @param b Second akaike struct.
 * @return -1 if the value of a is greater than that of b, 1 if smaller, 0 otherwise.
 */
int cmp_akaike(const void *a, const void *b) {
    akaike *a1 = (akaike *)a;
    akaike *a2 = (akaike *)b;
    if ((*a1).value > (*a2).value)
        return -1;
    else if ((*a1).value < (*a2).value)
        return 1;
    else
        return 0;
}
