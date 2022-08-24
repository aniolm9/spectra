#define RECTANGULAR 0
#define HANN 1
#define HAMMING 2
#define BLACKMAN 3

double hann(int n, int M);
double hamming(int n, int M);
double blackman(int n, int M);
void windowing(int M, const double *data, int window, double *output);
