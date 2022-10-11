#define MEAN 1
#define MEDIAN 2

typedef struct spectralOpts {
    float fs;       /* Sampling frequency of the data. */
    int window;     /* Type of window. See windows.h for the available options. */
    int nperseg;    /* Length of each segment. */
    int noverlap;   /* Number of points to overlap between segments. */
    int nfft;       /* Length of the FFT used. */
    bool isComplex; /* Real or complex input sequence? */
    int sides;      /* Sides of the FFT for real data (1 or 2). For complex data this parameter is not taken into account. */
    int scaling;    /* Compute the PSD (V^2/Hz) or the spectrum (V^2) if data is measured in V and fs in Hz. */
    int average;    /* Method to use when averaging periodograms (mean or median). */
} spectralOpts;

spectralOpts new_spectral_opts_basic(int window, int nperseg, int scaling);
void check_spectral_opts(spectralOpts spectralOpts);
