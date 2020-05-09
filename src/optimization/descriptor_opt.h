#pragma once 

#include <math.h>
#include <stdlib.h>
#include "integral_image.h"
#include "interest_point.h"
#include "helper.h"

// #define PATCH_SIZE 20

// inline float* get_gaussian(float sigma) {
//     /* computes matrix of shape (size x size) containing prob values of
//     2d gaussian with std=sigma and mean=(size/2, size/2) */
//     int size = PATCH_SIZE;
//     float* GW = (float*) malloc(size*size * sizeof(float));

//     float variance = sigma*sigma;
//     float normalization = 1/(2*M_PI*variance);
//     for (int i=0; i<size; i++) {
//         for (int j=0; j<size; j++) {
//             GW[j*size + i] = normalization * exp(-((i-size/2)*(i-size/2)+(j-size/2)*(j-size/2))/(2*variance));
//         }
//     }

//     return GW;
// }

inline float gaussianf(float x, float y, float sig) {
    return 1.0f / (2.0f * M_PI * sig*sig) * expf(-(x*x+y*y)/(2.0f*sig*sig));
}

void get_descriptor_inlinedHaarWavelets(struct integral_image* iimage, struct interest_point* ipoint, float* GW);

void get_msurf_descriptor_improved(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptor_inlined(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptor_inlinedHaarWavelets(struct integral_image* iimage, struct interest_point* ipoint);

inline void haarXY(float* ii_data, int height, int width, int row, int col, int scale, float* haarX, float* haarY) {
    // subtracting by one for row/col because row/col is inclusive.
    int r0 = MIN(row, height) - 1;         
    int c0 = MIN(col, width) - 1;         
    int r1 = MIN(row + scale, height) - 1;  
    int c1 = MIN(col + scale, width) - 1;   
    int r2 = MIN(row + 2*scale, height) - 1;
    int c2 = MIN(col + 2*scale, width) - 1;

    float r0c0 = 0.0f; // A
    float r0c1 = 0.0f; // B
    float r0c2 = 0.0f; // C
    float r1c0 = 0.0f; // D
    float r1c2 = 0.0f; // E
    float r2c0 = 0.0f; // F
    float r2c1 = 0.0f; // G
    float r2c2 = 0.0f; // H

    if (r0 >= 0 && c0 >= 0)  {
        r0c0 = ii_data[r0 * width + c0];
    }
    if (r0 >= 0 && c1 >= 0) {
        r0c1 = ii_data[r0 * width + c1];
    }
    if (r0 >= 0 && c2 >= 0) {
        r0c2 = ii_data[r0 * width + c2];
    }
    if (r1 >= 0 && c0 >= 0) {
        r1c0 = ii_data[r1 * width + c0];
    }
    if (r1 >= 0 && c2 >= 0) {
        r1c2 = ii_data[r1 * width + c2];
    }
    if (r2 >= 0 && c0 >= 0) {
        r2c0 = ii_data[r2 * width + c0];
    }
    if (r2 >= 0 && c1 >= 0) {
        r2c1 = ii_data[r2 * width + c1];
    }
    if (r2 >= 0 && c2 >= 0) {
        r2c2 = ii_data[r2 * width + c2];
    }

    float r0c0_sub_r2c2 = r0c0 - r2c2;
    float r0c2_sub_r2c0 = r0c2 - r2c0;

    *haarX = -1*(r0c0_sub_r2c2 - 2*(r0c1 - r2c1) + r0c2_sub_r2c0);
    *haarY = -1*(r0c0_sub_r2c2 - 2*(r1c0 - r1c2) - r0c2_sub_r2c0);
}

inline void haarXY_precheck_boundaries(float* ii_data, int height, int width, int col, int row, int scale, float* haarX, float* haarY) {
    // (row,col) is upper left corner of haar wavelet filter
    if (row <= 0 
        || col <= 0 
        || (row + 2*scale) > height 
        || (col + 2*scale) > width) {
        // wavelet filters that can not be applied completely will be skipped
        // I think it even makes sense to skip the whole keypoint if its 
        // descriptor goes over boundaries as it then can not be invariant 
        return;
    }
    
    // subtracting by one for row/col because row/col is inclusive.
    int r0 = row - 1;         
    int c0 = col - 1;         
    int r1 = row + scale - 1;  
    int c1 = col + scale - 1;   
    int r2 = row + 2*scale - 1;
    int c2 = col + 2*scale - 1;

    float r0c0 = ii_data[r0 * width + c0];
    float r0c1 = ii_data[r0 * width + c1];
    float r0c2 = ii_data[r0 * width + c2];
    float r1c0 = ii_data[r1 * width + c0];
    float r1c2 = ii_data[r1 * width + c2];
    float r2c0 = ii_data[r2 * width + c0];
    float r2c1 = ii_data[r2 * width + c1];
    float r2c2 = ii_data[r2 * width + c2];

    float r0c0_sub_r2c2 = r0c0 - r2c2;
    float r0c2_sub_r2c0 = r0c2 - r2c0;

    *haarX = r0c0_sub_r2c2 - 2*(r0c1 - r2c1) + r0c2_sub_r2c0;
    *haarY = r0c0_sub_r2c2 - 2*(r1c0 - r1c2) - r0c2_sub_r2c0;
}

void get_msurf_descriptor_precompute_gauss_s2(struct integral_image* iimage, struct interest_point* ipoint);
