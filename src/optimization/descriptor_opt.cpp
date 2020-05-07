#include "descriptor.h"
#include "descriptor_opt.h"
#include "integral_image.h"
#include "interest_point.h"

#define _USE_MATH_DEFINES

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PRECHECK_BOUNDARIES 1

void get_descriptor_inlinedHaarWavelets(struct integral_image* iimage, struct interest_point* ipoint, float* GW) {

    float scale = ipoint->scale;
    // TODO: (Sebastian) Is this correct with  "- 0.5"?
    int ipoint_x = (int) (ipoint->x + 0.5);
    int ipoint_y = (int) (ipoint->y + 0.5);

    int step = MAX((int)(scale + 0.5),1); // rounding is done this way in the original implementaion

    int col_offset = ipoint_x-step*10 - step/2;
    int row_offset = ipoint_y-step*10 - step/2;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    float *data = (float *)iimage->data;
    int width = iimage->width;
    int height = iimage->height;

    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) { // iterate over 4x4 sub_patches 
            float sum_x = 0.0f;
            float sum_y = 0.0f;
            float abs_x = 0.0f;
            float abs_y = 0.0f;
            // float missed_scaling = 1;
            for (int k=i*5; k<i*5+5; ++k) {
                for (int l=j*5; l<j*5+5; ++l) { 
                    // iterate over 5x5 sample points of sub_patch[i][j]
                    // and compute haar wavelet responses

                    float gw = GW[k*PATCH_SIZE + l];

                    // compute left corner of the [2*scale x 2*scale] patch:
                    int row = row_offset + l*step;
                    int col = col_offset + k*step;

                    float haarX = 0.0f;
                    float haarY = 0.0f;

                    // TODO: (valentin) can parameter passing be done better?
#if (!PRECHECK_BOUNDARIES)
                    haarXY(data, height, width, row, col, step, &haarX, &haarY);
#else
                    haarXY_precheck_boundaries(data, height, width, row, col, step, &haarX, &haarY);
#endif
                    // haarX
                    float x = gw * haarX;
                    float y = gw * haarY;

                    sum_x += x; // sum(x)
                    sum_y += y; // sum(y)
                    abs_x += fabs(x); // sum(abs(x))
                    abs_y += fabs(y); // sum(abs(y))

                    if (row <= 0 
                        || col <= 0 
                        || (row + 2*step) > height 
                        || (col + 2*step) > width) {
                        // wavelets that can not be applied completely will be skipped
                        // this is not just more efficeint but I believe may also be more accurate
                        // as cut off wavelet responses can distort features
                        // printf("%f %f\n", x, y);
                        continue;
                    }

                }
            }

            descriptor[desc_idx] = sum_x;
            descriptor[desc_idx+1] = sum_y;   
            descriptor[desc_idx+2] = abs_x;
            descriptor[desc_idx+3] = abs_y;
            
            // precompute for normaliztion
            sum_of_squares += sum_x * sum_x + sum_y * sum_y + abs_x * abs_x + abs_y * abs_y;

            desc_idx += 4;
        }
    }

    // rescale to unit vector
    float norm_factor = 1./sqrt(sum_of_squares);

    for (int i=0; i<64; i+=4) {
        descriptor[i] *= norm_factor;
        descriptor[i+1] *= norm_factor;
        descriptor[i+2] *= norm_factor;
        descriptor[i+3] *= norm_factor;
    }

}


void get_msurf_descriptor_inlinedHaarWavelets(struct integral_image* iimage, struct interest_point* ipoint) {

    float scale = ipoint->scale;
    int ipoint_x = (int) round(ipoint->x);
    int ipoint_y = (int) round(ipoint->y);

    float *data = (float *)iimage->data;
    int width = iimage->width;
    int height = iimage->height;

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
                    //float gauss_s1 = gaussian(xs-sample_x, ys-sample_y, 2.5f * scale);
                    //(1.0f/(2.0f*pi*sig*sig)) * exp( -(x*x+y*y)/(2.0f*sig*sig))
                    float g_factor = 0.08f / (scale*scale); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
                    float g_x = xs - sample_x;
                    float g_y = ys - sample_y;
                    float gauss_s1 = M_1_PI * g_factor * exp(-g_factor * (g_x*g_x + g_y*g_y));
                    
                    //float rx = haarX(sample_y, sample_x, (int) 2.0 * round(scale));
                    //float ry = haarY(sample_y, sample_x, (int) 2.0 * round(scale));
                    int s = (int) round(scale);
                    float rx = 0.0f; //box_integral(iimage, sample_y-s, sample_x, 2*s, s) - box_integral(iimage, sample_y-s, sample_x-s, 2*s, s);
                    float ry = 0.0f; //box_integral(iimage, sample_y, sample_x-s, s, 2*s) - box_integral(iimage, sample_y-s, sample_x-s, s, 2*s);

#if (!PRECHECK_BOUNDARIES)
                    haarXY(data, height, width, sample_y-s, sample_x-s, s, &rx, &ry);
#else
                    haarXY_precheck_boundaries(data, height, width, sample_y-s, sample_x-s, s, &rx, &ry);
#endif
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    //float rrx = gauss_s1*(-rx*si + ry*co);
                    //float rry = gauss_s1*(rx*co + ry*si);
                    float rrx = gauss_s1 * (ry);
                    float rry = gauss_s1 * (rx);

                    dx += rrx;
                    dy += rry;
                    mdx += fabs(rrx);
                    mdy += fabs(rry);

                }
            }

            // TODO: (Sebastian) Precompute this...
            //float gauss_s2 = gaussian(cx-2.0f,cy-2.0f,1.5f);
            //(1.0f/(2.0f*pi*sig*sig)) * exp( -(x*x+y*y)/(2.0f*sig*sig));
            float g_factor = 1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)
            float g_x = cx - 2.0f;
            float g_y = cy - 2.0f;
            float gauss_s2 = M_1_PI * g_factor * exp(-g_factor * (g_x*g_x + g_y*g_y));

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
