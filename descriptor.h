#pragma once 

#include <math.h>
#include <stdlib.h>
#include "integral_image.h"
#include "interest_point.h"

// initializing patch size
int PATCH_SIZE = 20;

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

static inline void get_descriptor(struct integral_image* iimage, struct interest_point* ipoint, float* GW) {

    int width = iimage->width;
    int height = iimage->height;

    float scale = ipoint->scale;
    float ipoint_x = ipoint->x;
    float ipoint_y = ipoint->y;

    float col_offset = ipoint_x-PATCH_SIZE*scale/2;
    float row_offset = ipoint_y-PATCH_SIZE*scale/2;

    // store wavelet responses
    float dx[PATCH_SIZE][PATCH_SIZE];
    float dy[PATCH_SIZE][PATCH_SIZE];

    // compute haar wavelet responses
    for (int i=0; i<PATCH_SIZE; i++) { // x coordinate
        for (int j=0; j<PATCH_SIZE; j++) { // y coordinate
            float gw = GW[i*PATCH_SIZE + j];

            // compute needed corners of the [2*scale x 2*scale] patch 
            // centered at (ipoint_x, ipoint_y)

            // c1 XX c2 XX 
            // XX XX XX XX
            // c4 XX XX XX
            // XX XX XX XX

            int c1_row = round(col_offset + (j)*scale);
            int c1_col = round(row_offset + (i)*scale);

            int c2_row = round(col_offset + (j)*scale);
            int c2_col = round(row_offset + (i+1)*scale);

            int c4_row = round(col_offset + (j+1)*scale);
            int c4_col = round(row_offset + (i)*scale);

            dx[i][j] = gw * (box_integral(iimage, c1_row, c1_col, round(scale), round(2*scale)) 
                            - box_integral(iimage, c2_row, c2_col, round(scale), round(2*scale)));

            dx[i][j] = gw * (box_integral(iimage, c1_row, c1_col, round(2*scale), round(scale)) 
                            - box_integral(iimage, c4_row, c4_col, round(2*scale), round(scale)));
        }
    }

    // build descriptor
    float* descriptor = ipoint->descriptor; //(float *) malloc(64 * sizeof(float));
    int desc_idx = 0;
    float sum_of_squares = 0;

    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) { // iterate over 4x4 sub_patches 
            descriptor[desc_idx] = 0;
            descriptor[desc_idx+1] = 0;
            descriptor[desc_idx+2] = 0;
            descriptor[desc_idx+3] = 0;
            for (int k=i*5; k<i*5+5; k++) {
                for (int l=j*5; l<j*5+5; l++) { // iterate over 5x5 sample points sub_patch[i][j]
                    float x = dx[k][l]; 
                    float y = dy[k][l]; 

                    descriptor[desc_idx] += x; // sum(x)
                    descriptor[desc_idx+1] += y; // sum(y)
                    descriptor[desc_idx+2] += (float)fabs(x); // sum(abs(x))
                    descriptor[desc_idx+3] += (float)fabs(y); // sum(abs(y))
                }   
            }
            
            // precompute for normaliztion
            for (int m=0; m<4; m++) 
                sum_of_squares += descriptor[m]*descriptor[m];

            desc_idx += 4;
        }
    }

    // rescale to unit vector
    float norm_factor = 1./sqrt(sum_of_squares);
    for (int i=0; i<64; i++) 
        descriptor[i] *= norm_factor;

}