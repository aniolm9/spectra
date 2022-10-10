#define MEAN 1
#define MEDIAN 2

void periodogram(double *data, double *freqs, double *power, int N, opts *welchOpts);
void welch(double *data, double *freqs, double *power, int N, opts *welchOpts);
void IQ2fftcpx(double *iq, kiss_fft_cpx *cpx, int N);
int compute_num_frames(int N, opts *spectralOpts);
