#include "descriptor.h"
#include "integral_image.h"
#include "interest_point.h"

#define _USE_MATH_DEFINES

#include <math.h>
#include <stdlib.h>

#define MAX(a,b) (((a)>(b))?(a):(b))


void get_descriptor(struct integral_image* iimage, struct interest_point* ipoint, float* GW) {

    float scale = ipoint->scale;
    // TODO: (Sebastian) Is this correct with  "- 0.5"?
    int ipoint_x = (int) (ipoint->x + 0.5);
    int ipoint_y = (int) (ipoint->y + 0.5);

    int step = MAX((int)(scale/2 + 0.5),1); // rounding is done this way in the original implementaion

    int col_offset = ipoint_x-step*11; //10 - 1 = PATCH_SIZE/2 + (shift to obtain upper left);
    int row_offset = ipoint_y-step*11;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) { // iterate over 4x4 sub_patches 
            float sum_x = 0.0f;
            float sum_y = 0.0f;
            float abs_x = 0.0f;
            float abs_y = 0.0f;
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

                    sum_x += x; // sum(x)
                    sum_y += y; // sum(y)
                    // TODO: (Sebastian) Why cast here?
                    abs_x += (float)fabs(x); // sum(abs(x))
                    abs_y += (float)fabs(y); // sum(abs(y))
                }
            }
            descriptor[desc_idx] = sum_x;
            descriptor[desc_idx+1] = sum_y;   
            descriptor[desc_idx+2] = abs_x;
            descriptor[desc_idx+3] = abs_y;
            
            // precompute for normaliztion
            //for (int m=0; m<4; ++m) 
            //    sum_of_squares += descriptor[desc_idx+m]*descriptor[desc_idx+m];
            sum_of_squares += sum_x * sum_x + sum_y * sum_y + abs_x * abs_x + abs_y * abs_y;


            desc_idx += 4;
        }
    }

    // rescale to unit vector
    float norm_factor = 1./sqrt(sum_of_squares);

    for (int i=0; i<64; ++i) {
        descriptor[i] *= norm_factor;
    }

}


void get_msurf_descriptor(struct integral_image* iimage, struct interest_point* ipoint) {

    float scale = ipoint->scale;
    int ipoint_x = (int) round(ipoint->x);
    int ipoint_y = (int) round(ipoint->y);

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;
    
    // subregion centers for the 4x4 gaussian weighting
    float cx = -0.5f;
    float cy = 0.0f;

    //float co = 1;
    //float si = 0;

    int i = -8;

    // calculate descriptor for this interest point
    while(i < 12) {

        int j = -8;
        i = i - 4;

        cx += 1.0f;
        cy = -0.5f;

        while(j < 12) {
            
            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            cy += 1.0f;

            j = j - 4;

            //int xs = (int) round(ipoint_x + ( -jx*scale*si + ix*scale*co));
            //int ys = (int) round(ipoint_y + ( jx*scale*co + ix*scale*si));
            int xs = (int) round(ipoint_x + (i + 5) * scale);
            int ys = (int) round(ipoint_y + (j + 5) * scale);

            for (int k = i; k < i + 9; ++k) {
                for (int l = j; l < j + 9; ++l) {

                    //Get coords of sample point on the rotated axis
                    //int sample_x = (int) round(ipoint_x + (-l*scale*si + k*scale*co));
                    //int sample_y = (int) round(ipoint_y + ( l*scale*co + k*scale*si));
                    int sample_x = (int) round(ipoint_x + k * scale);
                    int sample_y = (int) round(ipoint_y + l * scale);

                    //Get the gaussian weighted x and y responses
                    // TODO: (Sebastian) Precompute this...
                    float gauss_s1 = gaussian((float) xs-sample_x, (float) ys-sample_y, 2.5f * scale);
                                        
                    float rx = haarX(iimage, sample_y, sample_x, (int) 2.0 * round(scale));
                    float ry = haarY(iimage, sample_y, sample_x, (int) 2.0 * round(scale));
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    //float rrx = gauss_s1*(-rx*si + ry*co);
                    //float rry = gauss_s1*(rx*co + ry*si);
                    // TODO: (Sebastian) This seems weird...
                    float rrx = gauss_s1 * (ry);
                    float rry = gauss_s1 * (rx);

                    dx += rrx;
                    dy += rry;
                    mdx += fabs(rrx);
                    mdy += fabs(rry);

                }
            }

            // TODO: (Sebastian) Precompute this...
            float gauss_s2 = gaussian(cx-2.0f, cy-2.0f, 1.5f);

            // add the values to the descriptor vector
            descriptor[desc_idx] = dx * gauss_s2;
            descriptor[desc_idx+1] = dy * gauss_s2;
            descriptor[desc_idx+2] = mdx * gauss_s2;
            descriptor[desc_idx+3] = mdy * gauss_s2;

            // precompute for normaliztion
            sum_of_squares += (dx*dx + dy*dy + mdx*mdx + mdy*mdy) * gauss_s2*gauss_s2;

            desc_idx += 4;

            j += 9;

        }

        i += 9;

    }

    // rescale to unit vector
    float norm_factor = 1.0f / sqrt(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }

}
