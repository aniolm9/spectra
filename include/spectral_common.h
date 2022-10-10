void new_spectral_opts();
void check_spectral_opts(opts *spectralOpts);
void fftfreq(double *freqs, int M, float ts);
void spectral_helper(double *data_x, double *data_y, double *freqs, double **psd, int N_x, int N_y, opts *spectralOpts);
