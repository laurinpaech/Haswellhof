#pragma once 

#include <math.h>
#include <stdlib.h>
#include "integral_image.h"
#include "interest_point.h"

#define PATCH_SIZE 20

static inline float* get_gaussian(float sigma) {
    /* computes matrix of shape (size x size) containing prob values of
    2d gaussian with std=sigma and mean=(size/2, size/2) */
    int size = PATCH_SIZE;
    float* GW = (float*) malloc(size*size * sizeof(float));

    float variance = sigma*sigma;
    float normalization = 1/(2*M_PI*variance);
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            GW[j*size + i] = normalization * exp(-((i-size/2)*(i-size/2)+(j-size/2)*(j-size/2))/(2*variance));
        }
    }

    return GW;
}

void get_descriptor(struct integral_image* iimage, struct interest_point* ipoint, float* GW);