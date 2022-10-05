#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <kissfft/kiss_fft.h>
#include "spectrum.h"
#include "windows.h"

int main() {
    /* Open file */
    const char *filename = "sdr_iq_data_0G433000_60dB-20210505_184759.iqdat";
    FILE *fp = fopen(filename, "rb");
    /* Get filesize and allocate memory */
    fseek(fp, 0, SEEK_END);
    size_t size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    printf("FILE SIZE: %lu\n", (unsigned long)size);
    int16_t *data = malloc(size);
    int samples = size / 2;
    /* Read bytes */
    fread(data, size, 1, fp);
    double *datad = malloc(samples*sizeof(double));
    for (int i = 0; i < samples; i++) {
        datad[i] = (double)data[i];
    }
    int nperseg = 2048;
    double *freqs;
    double *power = calloc(nperseg, sizeof(*power));
    bool *signal_presence = malloc(samples/2 * sizeof(bool));
    opts wopts = {1.0, 1, nperseg, 1024, nperseg, true, 1, 1, MEAN};
    periodogram(datad, freqs, power, samples, &wopts);
    double noise_power = noise_power_aic(power, samples, &wopts);
    printf("Estimated noise power = %lf\n", noise_power);
    energy_detector(datad, signal_presence, samples, noise_power, 0.05, 4096);
    /* Free memory */
    fclose(fp);
    free(data);
    free(power);
    free(signal_presence);
}