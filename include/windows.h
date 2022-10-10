/* Types of window */
#define RECTANGULAR 0
#define HANN 1
#define HAMMING 2
#define BLACKMAN 3

/* Scaling methods */
#define DENSITY 1
#define SPECTRUM 2

kiss_fft_scalar hann(int n, int M);
kiss_fft_scalar hamming(int n, int M);
kiss_fft_scalar blackman(int n, int M);
kiss_fft_scalar windowing(const kiss_fft_cpx *data, kiss_fft_cpx *output, int M, int window, float fs, int scaling);
