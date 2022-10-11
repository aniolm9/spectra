typedef struct {
    int32_t value;
    uint16_t index;
} akaike;

void IQ2fftcpx(double *iq, kiss_fft_cpx *cpx, int N);
int compute_num_frames(int N, spectralOpts *spectralOpts);
int cmp_akaike(const void *a, const void *b);
