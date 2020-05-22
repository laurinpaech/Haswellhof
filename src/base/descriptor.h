#pragma once 

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include "integral_image.h"
#include "interest_point.h"
#include <vector>


#define PATCH_SIZE 20

inline float* get_gaussian(float sigma) {
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

inline float gaussian(float x, float y, float sig) {
    return 1.0f / (2.0f * M_PI * sig*sig) * exp(-(x*x+y*y)/(2.0f*sig*sig)); // 2 div + 8 x mul + 1 x add + 1 x exp
}

inline float haarX(struct integral_image *iimage, int row, int col, int s) {
    return box_integral(iimage, row-s/2, col, s, s/2) - box_integral(iimage, row-s/2, col-s/2, s, s/2); // 1 x add + 2 x box_integral(3 x add + 1 x fmax) => Total: 7 x add + 2 x fmax
}

inline float haarY(struct integral_image *iimage, int row, int col, int s) {
    return box_integral(iimage, row, col-s/2, s/2, s) - box_integral(iimage, row-s/2, col-s/2, s/2, s); // 1 x add + 2 x box_integral(3 x add + 1 x fmax) => Total: 7 x add + 2 x fmax
}

void get_msurf_descriptor(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);
