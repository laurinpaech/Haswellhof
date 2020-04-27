#include "descriptor.h"
#include "integral_image.h"
#include "interest_point.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>



void get_descriptor(struct integral_image* iimage, struct interest_point* ipoint, float* GW) {

    // int width = iimage->width;
    // int height = iimage->height;

    float scale = ipoint->scale;
    float ipoint_x = ipoint->x;
    float ipoint_y = ipoint->y;

    float col_offset = ipoint_x-PATCH_SIZE*scale/2;
    float row_offset = ipoint_y-PATCH_SIZE*scale/2;

    // printf("%f %f %f %f\n", col_offset, row_offset, ipoint_x, ipoint_y);

    // build descriptor
    float* descriptor = ipoint->descriptor; //(float *) malloc(64 * sizeof(float));
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

                    int c1_row = round(row_offset +  l   *scale);
                    int c1_col = round(col_offset +  k   *scale);

                    // int c2_row = round(row_offset +  l   *scale);
                    int c2_col = round(col_offset + (k+1)*scale);

                    int c4_row = round(row_offset + (l+1)*scale);
                    // int c4_col = round(col_offset +  k   *scale);

                    // haarX
                    // printf("haarX\n");
                    float x = gw * (box_integral(iimage, c1_row, c1_col, round(2*scale), round(scale))
                                    - box_integral(iimage, c1_row, c2_col, round(2*scale), round(scale)));

                    // haarY
                    // printf("haarY\n");
                    float y = gw * (box_integral(iimage, c1_row, c1_col, round(scale), round(2*scale)) 
                                    - box_integral(iimage, c4_row, c1_col, round(scale), round(2*scale)));

                    // printf("c1:(%d %d) c2_col:%d c4_row:%d x:%f y:%f gw:%f\n", c1_row, c1_col, c2_col, c4_row, x, y, gw);

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
    // printf("%f %f \n", sum_of_squares, norm_factor);

    for (int i=0; i<64; ++i) 
        descriptor[i] *= norm_factor;

}