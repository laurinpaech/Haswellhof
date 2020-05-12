#include "descriptor.h"
#include "descriptor_opt.h"
#include "integral_image.h"
#include "interest_point.h"

#define _USE_MATH_DEFINES

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>

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
                    abs_x += fabsf(x); // sum(abs(x))
                    abs_y += fabsf(y); // sum(abs(y))

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
    float norm_factor = 1./sqrtf(sum_of_squares);

    for (int i=0; i<64; i+=4) {
        descriptor[i] *= norm_factor;
        descriptor[i+1] *= norm_factor;
        descriptor[i+2] *= norm_factor;
        descriptor[i+3] *= norm_factor;
    }

}


void get_msurf_descriptor_improved(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - replaced outer wihle loops with for loops, simplifying index calculation
        - scalar replacement & FLOP reduction by computing in the outer most (possible) loop
        - changed math functions to their float counterparts (??? or is that part of base)
    
    ideas:
        - flip outer for loops (this will change how desc_idx hat to be incremeted to yield the same feature vector)
        otherwise the result is a mere permutation, doubt this helps locality much)
    */

    float scale = ipoint->scale;
    float scale_mul_25f = 2.5f*scale;
    int int_scale_mul_2 = (int) 2 * roundf(scale);

    int ipoint_x = (int) roundf(ipoint->x);
    int ipoint_y = (int) roundf(ipoint->y);

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;
    
    // subregion centers for the 4x4 gaussian weighting
    float cx = -2.5f;
    float cy;

    // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        cx += 1.0f;
        cy = -2.5f;

        for (int j=-8; j<8; j+=5) {

            cy += 1.0f;
            
            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int xs = (int) roundf(ipoint_x + i * scale);
            int ys = (int) roundf(ipoint_y + j * scale);

            for (int k = i-4; k < i + 5; ++k) {

                //Get x coords of sample point
                int sample_x = (int) roundf(ipoint_x + k * scale);
                float xs_sub_sample_x = (float) xs-sample_x;

                for (int l = j-4; l < j + 5; ++l) {

                    //Get y coords of sample point
                    int sample_y = (int) roundf(ipoint_y + l * scale);
                    float ys_sub_sample_y = (float) ys-sample_y;

                

                    //Get the gaussian weighted x and y responses
                    // TODO: (Sebastian) Precompute this...
                    float gauss_s1 = gaussianf(xs_sub_sample_x, ys_sub_sample_y, scale_mul_25f);
                    
                    float rx = haarX(iimage, sample_y, sample_x, int_scale_mul_2);
                    float ry = haarY(iimage, sample_y, sample_x, int_scale_mul_2);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // TODO: (Sebastian) Precompute this...
            float gauss_s2 = gaussianf(cx, cy, 1.5f);

            // add the values to the descriptor vector
            float d1 = dx * gauss_s2;
            float d2 = dy * gauss_s2;
            float d3 = mdx * gauss_s2;
            float d4 = mdy * gauss_s2;

            descriptor[desc_idx] = d1;
            descriptor[desc_idx+1] = d2;
            descriptor[desc_idx+2] = d3;
            descriptor[desc_idx+3] = d4;

            // precompute for normaliztion
            sum_of_squares += (d1*d1 + d2*d2 + d3*d3 + d4*d4);

            desc_idx += 4;

        }
    }

    // rescale to unit vector
    // NOTE: using sqrtf() for floats
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}


void get_msurf_descriptors_improved(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_improved(iimage, &interest_points->at(i));
	}
}


void get_msurf_descriptor_improved_flip(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - same as get_msurf_descriptor_improved
        - switched order of inner loops to go along x direction for better locality
    */

    float scale = ipoint->scale;
    float scale_mul_25f = 2.5f*scale;
    int int_scale_mul_2 = (int) 2 * roundf(scale);

    int ipoint_x = (int) roundf(ipoint->x);
    int ipoint_y = (int) roundf(ipoint->y);

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;
    
    // subregion centers for the 4x4 gaussian weighting
    float cx = -2.5f;
    float cy;

    // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        cx += 1.0f;
        cy = -2.5f;

        for (int j=-8; j<8; j+=5) {

            cy += 1.0f;
            
            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int xs = (int) roundf(ipoint_x + i * scale);
            int ys = (int) roundf(ipoint_y + j * scale);

            for (int l = j-4; l < j + 5; ++l) {

                //Get y coords of sample point
                int sample_y = (int) roundf(ipoint_y + l * scale);
                float ys_sub_sample_y = (float) ys-sample_y;

                for (int k = i-4; k < i + 5; ++k) {

                    //Get x coords of sample point
                    int sample_x = (int) roundf(ipoint_x + k * scale);
                    float xs_sub_sample_x = (float) xs-sample_x;

                    //Get the gaussian weighted x and y responses
                    // TODO: (Sebastian) Precompute this...
                    float gauss_s1 = gaussianf(xs_sub_sample_x, ys_sub_sample_y, scale_mul_25f);
                    
                    float rx = haarX(iimage, sample_y, sample_x, int_scale_mul_2);
                    float ry = haarY(iimage, sample_y, sample_x, int_scale_mul_2);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // TODO: (Sebastian) Precompute this...
            float gauss_s2 = gaussianf(cx, cy, 1.5f);

            // add the values to the descriptor vector
            float d1 = dx * gauss_s2;
            float d2 = dy * gauss_s2;
            float d3 = mdx * gauss_s2;
            float d4 = mdy * gauss_s2;

            descriptor[desc_idx] = d1;
            descriptor[desc_idx+1] = d2;
            descriptor[desc_idx+2] = d3;
            descriptor[desc_idx+3] = d4;

            // precompute for normaliztion
            sum_of_squares += (d1*d1 + d2*d2 + d3*d3 + d4*d4);

            desc_idx += 4;

        }
    }

    // rescale to unit vector
    // NOTE: using sqrtf() for floats
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}


void get_msurf_descriptors_improved_flip(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_improved_flip(iimage, &interest_points->at(i));
	}
}


void get_msurf_descriptor_improved_flip_flip(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - same as get_msurf_descriptor_improved
        - switched order of inner and outer loops to go along x direction for better locality
    */

    float scale = ipoint->scale;
    float scale_mul_25f = 2.5f*scale;
    int int_scale_mul_2 = (int) 2 * roundf(scale);

    int ipoint_x = (int) roundf(ipoint->x);
    int ipoint_y = (int) roundf(ipoint->y);

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;
    
    // subregion centers for the 4x4 gaussian weighting
    float cy = -2.5f;
    float cx;

    // calculate descriptor for this interest point
    for (int j=-8; j<8; j+=5) {

        cy += 1.0f;
        cx = -2.5f;

        for (int i=-8; i<8; i+=5) {

            cx += 1.0f;
            
            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int xs = (int) roundf(ipoint_x + i * scale);
            int ys = (int) roundf(ipoint_y + j * scale);

            for (int l = j-4; l < j + 5; ++l) {

                //Get y coords of sample point
                int sample_y = (int) roundf(ipoint_y + l * scale);
                float ys_sub_sample_y = (float) ys-sample_y;

                for (int k = i-4; k < i + 5; ++k) {

                    //Get x coords of sample point
                    int sample_x = (int) roundf(ipoint_x + k * scale);
                    float xs_sub_sample_x = (float) xs-sample_x;

                    //Get the gaussian weighted x and y responses
                    // TODO: (Sebastian) Precompute this...
                    float gauss_s1 = gaussianf(xs_sub_sample_x, ys_sub_sample_y, scale_mul_25f);
                    
                    float rx = haarX(iimage, sample_y, sample_x, int_scale_mul_2);
                    float ry = haarY(iimage, sample_y, sample_x, int_scale_mul_2);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // TODO: (Sebastian) Precompute this...
            float gauss_s2 = gaussianf(cx, cy, 1.5f);

            // add the values to the descriptor vector
            float d1 = dx * gauss_s2;
            float d2 = dy * gauss_s2;
            float d3 = mdx * gauss_s2;
            float d4 = mdy * gauss_s2;

            descriptor[desc_idx] = d1;
            descriptor[desc_idx+1] = d2;
            descriptor[desc_idx+2] = d3;
            descriptor[desc_idx+3] = d4;

            // precompute for normaliztion
            sum_of_squares += (d1*d1 + d2*d2 + d3*d3 + d4*d4);

            desc_idx += 4;

        }
    }

    // rescale to unit vector
    // NOTE: using sqrtf() for floats
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}


void get_msurf_descriptors_improved_flip_flip(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_improved_flip_flip(iimage, &interest_points->at(i));
	}
}


void get_msurf_descriptor_inlined(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - all of get_msurf_descriptor_improved_flip
        - inlined gaussian funktion and removed its noralization (as we normalize in the end)

    ideas:
        - keep xs_sub_sample_x and xs_sub_sample_x_squared as int and cast in expf (same for y)
    */

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale_mul_2 = (int) 2 * roundf(scale);
    float scale_squared = scale*scale;
    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    int ipoint_x = (int) roundf(ipoint->x);
    int ipoint_y = (int) roundf(ipoint->y);

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;
    
    // subregion centers for the 4x4 gaussian weighting
    float cx = -2.5f;
    float cy;

    // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        cx += 1.0f;
        float cx_squared = cx*cx;
        
        cy = -2.5f;

        for (int j=-8; j<8; j+=5) {

            cy += 1.0f;
            float cy_squared = cy*cy;
            
            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int xs = (int) roundf(ipoint_x + i * scale);
            int ys = (int) roundf(ipoint_y + j * scale);

            for (int l = j-4; l < j + 5; ++l) {

                //Get y coords of sample point
                int sample_y = (int) roundf(ipoint_y + l * scale);
                float ys_sub_sample_y = (float) ys-sample_y;
                float ys_sub_sample_y_squared = ys_sub_sample_y*ys_sub_sample_y;

                for (int k = i-4; k < i + 5; ++k) {

                    //Get x coords of sample point
                    int sample_x = (int) roundf(ipoint_x + k * scale);
                    float xs_sub_sample_x = (float) xs-sample_x;
                    float xs_sub_sample_x_squared = xs_sub_sample_x*xs_sub_sample_x;

                    //Get the gaussian weighted x and y responses
                    // NOTE: We use expf here
                    float gauss_s1 = expf(g1_factor * (xs_sub_sample_x_squared + ys_sub_sample_y_squared));
                    
                    float rx = haarX(iimage, sample_y, sample_x, int_scale_mul_2);
                    float ry = haarY(iimage, sample_y, sample_x, int_scale_mul_2);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // TODO: (Sebastian) Precompute this...
            // NOTE: We use expf here
            float gauss_s2 = expf(g2_factor * (cx_squared + cy_squared)); 

            // add the values to the descriptor vector
            float d1 = dx * gauss_s2;
            float d2 = dy * gauss_s2;
            float d3 = mdx * gauss_s2;
            float d4 = mdy * gauss_s2;

            descriptor[desc_idx] = d1;
            descriptor[desc_idx+1] = d2;
            descriptor[desc_idx+2] = d3;
            descriptor[desc_idx+3] = d4;

            // precompute for normaliztion
            sum_of_squares += (d1*d1 + d2*d2 + d3*d3 + d4*d4);

            desc_idx += 4;

        }
    }

    // rescale to unit vector
    // NOTE: using sqrtf() for floats
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}


void get_msurf_descriptors_inlined(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptors_inlined(iimage, &interest_points->at(i));
	}
}


void get_msurf_descriptor_inlinedHaarWavelets(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - all of get_msurf_descriptor_inlined
        - using haarXY function to compute Haar Wavelet responses at once
    */

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale = (int) roundf(scale);
    // int int_scale_mul_2 =  2 * int_scale;
    float scale_squared = scale*scale;

    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    int ipoint_x = (int) roundf(ipoint->x);
    int ipoint_y = (int) roundf(ipoint->y);

    float *data = (float *)iimage->data;
    int width = iimage->width;
    int height = iimage->height;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;
    
    // subregion centers for the 4x4 gaussian weighting
    float cx = -2.5f;
    float cy;

    // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        cx += 1.0f;
        float cx_squared = cx*cx;
        
        cy = -2.5f;

        for (int j=-8; j<8; j+=5) {

            cy += 1.0f;
            float cy_squared = cy*cy;
            
            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int xs = (int) roundf(ipoint_x + i * scale);
            int ys = (int) roundf(ipoint_y + j * scale);

            for (int l = j-4; l < j + 5; ++l) {

                //Get y coords of sample point
                int sample_y = (int) roundf(ipoint_y + l * scale);
                float ys_sub_sample_y = (float) ys-sample_y;
                float ys_sub_sample_y_squared = ys_sub_sample_y*ys_sub_sample_y;

                int sample_y_sub_int_scale = sample_y-int_scale;

                for (int k = i-4; k < i + 5; ++k) {

                    //Get x coords of sample point
                    int sample_x = (int) roundf(ipoint_x + k * scale);
                    float xs_sub_sample_x = (float) xs-sample_x;
                    float xs_sub_sample_x_squared = xs_sub_sample_x*xs_sub_sample_x;
                    int sample_x_sub_int_scale = sample_x-int_scale;

                    //Get the gaussian weighted x and y responses
                    // NOTE: We use expf here
                    float gauss_s1 = expf(g1_factor * (xs_sub_sample_x_squared + ys_sub_sample_y_squared));
                    
                    float rx = 0.0f;
                    float ry = 0.0f;
                    haarXY(data, height, width, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &rx, &ry);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // TODO: (Sebastian) Precompute this...
            // NOTE: We use expf here
            float gauss_s2 = expf(g2_factor * (cx_squared + cy_squared)); 

            // add the values to the descriptor vector
            float d1 = dx * gauss_s2;
            float d2 = dy * gauss_s2;
            float d3 = mdx * gauss_s2;
            float d4 = mdy * gauss_s2;

            descriptor[desc_idx] = d1;
            descriptor[desc_idx+1] = d2;
            descriptor[desc_idx+2] = d3;
            descriptor[desc_idx+3] = d4;

            // precompute for normaliztion
            sum_of_squares += (d1*d1 + d2*d2 + d3*d3 + d4*d4);

            desc_idx += 4;

        }
    }

    // rescale to unit vector
    // NOTE: using sqrtf() for floats
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}


void get_msurf_descriptors_inlinedHaarWavelets(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_inlinedHaarWavelets(iimage, &interest_points->at(i));
	}
}


void get_msurf_descriptor_inlinedHaarWavelets_precheck_boundaries(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - all of get_msurf_descriptor_inlined
        - using haarXY_precheck_boundaries function to compute Haar Wavelet responses at once with very naive reduction of if/MIN usage

    ideas:
        - more sophisticated boundary checking (before the haarXY_precheck_boundaries function call)
    */

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale = (int) roundf(scale);
    // int int_scale_mul_2 =  2 * int_scale;
    float scale_squared = scale*scale;

    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    int ipoint_x = (int) roundf(ipoint->x);
    int ipoint_y = (int) roundf(ipoint->y);

    float *data = (float *)iimage->data;
    int width = iimage->width;
    int height = iimage->height;


    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;
    
    // subregion centers for the 4x4 gaussian weighting
    float cx = -2.5f;
    float cy;


    // check if we ever hit a boundary
    if (((int) roundf(ipoint_x - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_y - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_x + 11*scale)) + int_scale > width 
        || ((int) roundf(ipoint_y + 11*scale)) + int_scale > height) 
    { // some outside
    for (int i=-8; i<8; i+=5) {

        cx += 1.0f;
        float cx_squared = cx*cx;
        
        cy = -2.5f;

        for (int j=-8; j<8; j+=5) {

            cy += 1.0f;
            float cy_squared = cy*cy;
            
            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int xs = (int) roundf(ipoint_x + i * scale);
            int ys = (int) roundf(ipoint_y + j * scale);

            for (int l = j-4; l < j + 5; ++l) {

                //Get y coords of sample point
                int sample_y = (int) roundf(ipoint_y + l * scale);
                float ys_sub_sample_y = (float) ys-sample_y;
                float ys_sub_sample_y_squared = ys_sub_sample_y*ys_sub_sample_y;

                int sample_y_sub_int_scale = sample_y-int_scale;

                for (int k = i-4; k < i + 5; ++k) {

                    //Get x coords of sample point
                    int sample_x = (int) roundf(ipoint_x + k * scale);
                    float xs_sub_sample_x = (float) xs-sample_x;
                    float xs_sub_sample_x_squared = xs_sub_sample_x*xs_sub_sample_x;
                    int sample_x_sub_int_scale = sample_x-int_scale;

                    //Get the gaussian weighted x and y responses
                    // NOTE: We use expf here
                    float gauss_s1 = expf(g1_factor * (xs_sub_sample_x_squared + ys_sub_sample_y_squared));
                    
                    float rx = 0.0f;
                    float ry = 0.0f;
                    haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &rx, &ry);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // TODO: (Sebastian) Precompute this...
            // NOTE: We use expf here
            float gauss_s2 = expf(g2_factor * (cx_squared + cy_squared)); 

            // add the values to the descriptor vector
            float d1 = dx * gauss_s2;
            float d2 = dy * gauss_s2;
            float d3 = mdx * gauss_s2;
            float d4 = mdy * gauss_s2;

            descriptor[desc_idx] = d1;
            descriptor[desc_idx+1] = d2;
            descriptor[desc_idx+2] = d3;
            descriptor[desc_idx+3] = d4;

            // precompute for normaliztion
            sum_of_squares += (d1*d1 + d2*d2 + d3*d3 + d4*d4);

            desc_idx += 4;

        }
    }
    

    } else { // everything in borders
    // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        cx += 1.0f;
        float cx_squared = cx*cx;
        
        cy = -2.5f;

        for (int j=-8; j<8; j+=5) {

            cy += 1.0f;
            float cy_squared = cy*cy;
            
            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int xs = (int) roundf(ipoint_x + i * scale);
            int ys = (int) roundf(ipoint_y + j * scale);

            for (int l = j-4; l < j + 5; ++l) {

                //Get y coords of sample point
                int sample_y = (int) roundf(ipoint_y + l * scale);
                float ys_sub_sample_y = (float) ys-sample_y;
                float ys_sub_sample_y_squared = ys_sub_sample_y*ys_sub_sample_y;

                int sample_y_sub_int_scale = sample_y-int_scale;

                for (int k = i-4; k < i + 5; ++k) {

                    //Get x coords of sample point
                    int sample_x = (int) roundf(ipoint_x + k * scale);
                    float xs_sub_sample_x = (float) xs-sample_x;
                    float xs_sub_sample_x_squared = xs_sub_sample_x*xs_sub_sample_x;
                    int sample_x_sub_int_scale = sample_x-int_scale;

                    //Get the gaussian weighted x and y responses
                    // NOTE: We use expf here
                    float gauss_s1 = expf(g1_factor * (xs_sub_sample_x_squared + ys_sub_sample_y_squared));
                    
                    float rx = 0.0f;
                    float ry = 0.0f;
                    haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &rx, &ry);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // TODO: (Sebastian) Precompute this...
            // NOTE: We use expf here
            float gauss_s2 = expf(g2_factor * (cx_squared + cy_squared)); 

            // add the values to the descriptor vector
            float d1 = dx * gauss_s2;
            float d2 = dy * gauss_s2;
            float d3 = mdx * gauss_s2;
            float d4 = mdy * gauss_s2;

            descriptor[desc_idx] = d1;
            descriptor[desc_idx+1] = d2;
            descriptor[desc_idx+2] = d3;
            descriptor[desc_idx+3] = d4;

            // precompute for normaliztion
            sum_of_squares += (d1*d1 + d2*d2 + d3*d3 + d4*d4);

            desc_idx += 4;

        }
    }
    }
    // rescale to unit vector
    // NOTE: using sqrtf() for floats
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}


void get_msurf_descriptors_inlinedHaarWavelets_precheck_boundaries(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_inlinedHaarWavelets_precheck_boundaries(iimage, &interest_points->at(i));
	}
}


static const float gauss_s2_precomputed[] = {0.026022f, 0.040585f, 0.040585f, 0.026022f, 0.040585f, 0.063297f, 0.063297f, 0.040585f, 0.040585f, 0.063297f, 0.063297f, 0.040585f, 0.026022f, 0.040585f, 0.040585f, 0.026022f};



/*
static const float gauss_s2_precomputed[] = { 
    {0.026022f, 0.040585f, 0.040585f, 0.026022f}, 
    {0.040585f, 0.063297f, 0.063297f, 0.040585f}, 
    {0.040585f, 0.063297f, 0.063297f, 0.040585f}, 
    {0.026022f, 0.040585f, 0.040585f, 0.026022f}
};
*/

void get_msurf_descriptor_precompute_gauss_s2(struct integral_image* iimage, struct interest_point* ipoint) {

    float scale = ipoint->scale;
    int ipoint_x = (int) roundf(ipoint->x);
    int ipoint_y = (int) roundf(ipoint->y);

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    int i = -8;

    // calculate descriptor for this interest point
    while(i < 12) {

        int j = -8;
        i = i - 4;

        while(j < 12) {
            
            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            j = j - 4;

            int xs = (int) roundf(ipoint_x + (i + 4) * scale);
            int ys = (int) roundf(ipoint_y + (j + 4) * scale);

            for (int k = i; k < i + 9; ++k) {
                for (int l = j; l < j + 9; ++l) {

                    //Get coords of sample point on the rotated axis
                    int sample_x = (int) roundf(ipoint_x + k * scale);
                    int sample_y = (int) roundf(ipoint_y + l * scale);

                    //Get the gaussian weighted x and y responses
                    // TODO: (Sebastian) Precompute this...
                    float gauss_s1 = gaussianf((float) xs-sample_x, (float) ys-sample_y, 2.5f * scale);

                    float rx = haarX(iimage, sample_y, sample_x, (int) 2 * roundf(scale));
                    float ry = haarY(iimage, sample_y, sample_x, (int) 2 * roundf(scale));
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * (ry);
                    float rry = gauss_s1 * (rx);

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);

                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_precomputed[gauss_s2_index++];

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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }

}


void get_msurf_descriptors_precompute_gauss_s2(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_precompute_gauss_s2(iimage, &interest_points->at(i));
	}
}


void get_msurf_descriptor_separable_gauss_s1(struct integral_image* iimage, struct interest_point* ipoint) {

    float scale = ipoint->scale;
    int ipoint_x = (int) roundf(ipoint->x);
    int ipoint_y = (int) roundf(ipoint->y);

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    int i = -8;

    // calculate descriptor for this interest point
    while(i < 12) {

        int j = -8;
        i = i - 4;

        while(j < 12) {
            
            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            j = j - 4;

            int xs = (int) roundf(ipoint_x + (i + 4) * scale);
            int ys = (int) roundf(ipoint_y + (j + 4) * scale);

    

            for (int k = i; k < i + 9; ++k) {
                for (int l = j; l < j + 9; ++l) {

                    //Get coords of sample point on the rotated axis
                    int sample_x = (int) roundf(ipoint_x + k * scale);
                    int sample_y = (int) roundf(ipoint_y + l * scale);

                    //Get the gaussian weighted x and y responses
                    // TODO: (Sebastian) Precompute this...
                    float gauss_s1 = gaussianf((float) xs-sample_x, (float) ys-sample_y, 2.5f * scale);

                    float rx = haarX(iimage, sample_y, sample_x, (int) 2 * roundf(scale));
                    float ry = haarY(iimage, sample_y, sample_x, (int) 2 * roundf(scale));
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * (ry);
                    float rry = gauss_s1 * (rx);

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);

                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_precomputed[gauss_s2_index++];

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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }

}


void get_msurf_descriptors_separable_gauss_s1(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_separable_gauss_s1(iimage, &interest_points->at(i));
	}
}