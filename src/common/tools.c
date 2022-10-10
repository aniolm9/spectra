#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <kissfft/kiss_fft.h>
#include <gsl/gsl_cdf.h>
#include "opts.h"
#include "tools.h"

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
