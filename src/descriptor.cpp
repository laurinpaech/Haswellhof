#include "descriptor.h"
#include "integral_image.h"
#include "interest_point.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define MAX(a,b) (((a)>(b))?(a):(b))


void get_descriptor(struct integral_image* iimage, struct interest_point* ipoint, float* GW) {

    float scale = ipoint->scale;
    int ipoint_x = (int) (ipoint->x - 0.5);
    int ipoint_y = (int) (ipoint->y - 0.5);

    int step = MAX((int)(scale/2 + 0.5),1); // rounding is done this way in the original implementaion

    int col_offset = ipoint_x-step*10; //10 = PATCH_SIZE/2;
    int row_offset = ipoint_y-step*10;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0;

    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) { // iterate over 4x4 sub_patches 
            descriptor[desc_idx] = 0;
            descriptor[desc_idx+1] = 0;
            descriptor[desc_idx+2] = 0;
            descriptor[desc_idx+3] = 0;
            for (int k=i*5; k<i*5+5; ++k) {
                for (int l=j*5; l<j*5+5; ++l) { 
                    // iterate over 5x5 sample points of sub_patch[i][j]
                    // and compute haar wavelet responses

                    float gw = GW[k*PATCH_SIZE + l];

                    // compute needed corners of the [2*scale x 2*scale] patch:
                    // c1 XX c2 XX 
                    // XX XX XX XX
                    // c4 XX XX XX
                    // XX XX XX XX

                    int c1_row = row_offset +  l   *step;
                    int c1_col = col_offset +  k   *step;

                    int c2_col = col_offset + (k+1)*step;
                    int c4_row = row_offset + (l+1)*step;

                    // haarX
                    float x = gw * (box_integral(iimage, c1_row, c1_col, 2*step, step)
                                    - box_integral(iimage, c1_row, c2_col, 2*step, step));

                    // haarY
                    float y = gw * (box_integral(iimage, c1_row, c1_col, step, 2*step) 
                                    - box_integral(iimage, c4_row, c1_col, step, 2*step));


                    descriptor[desc_idx] += x; // sum(x)
                    descriptor[desc_idx+1] += y; // sum(y)
                    descriptor[desc_idx+2] += (float)fabs(x); // sum(abs(x))
                    descriptor[desc_idx+3] += (float)fabs(y); // sum(abs(y))
                }   
            }
            
            // precompute for normaliztion
            for (int m=0; m<4; ++m) 
                sum_of_squares += descriptor[desc_idx+m]*descriptor[desc_idx+m];

            desc_idx += 4;
        }
    }

    // rescale to unit vector
    float norm_factor = 1./sqrt(sum_of_squares);

    for (int i=0; i<64; ++i) 
        descriptor[i] *= norm_factor;

}