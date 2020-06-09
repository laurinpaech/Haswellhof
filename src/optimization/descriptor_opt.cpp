#include "descriptor.h"
#include "descriptor_opt.h"
#include "integral_image.h"
#include "interest_point.h"

#define _USE_MATH_DEFINES

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <immintrin.h>

#include <vector>


void get_msurf_descriptor_improved(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - replaced outer wihle loops with for loops, simplifying index calculation
        - scalar replacement & FLOP reduction by computing everything in the most outer loop
        - changed math functions to their correct float counterparts
    
    ideas:
        - flip outer for loops (this will change how desc_idx hat to be incremeted to yield the same feature vector)
        otherwise the result is a mere permutation, doubt this helps locality much)
    */

    float scale = ipoint->scale;
    float scale_mul_25f = 2.5f*scale;
    int int_scale_mul_2 = (int) 2 * roundf(scale);

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

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
                    
                    float rx = haarX_improved(iimage, sample_y, sample_x, int_scale_mul_2);
                    float ry = haarY_improved(iimage, sample_y, sample_x, int_scale_mul_2);
                    
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

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

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
                    
                    float rx = haarX_improved(iimage, sample_y, sample_x, int_scale_mul_2);
                    float ry = haarY_improved(iimage, sample_y, sample_x, int_scale_mul_2);
                    
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

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

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
                    
                    float rx = haarX_improved(iimage, sample_y, sample_x, int_scale_mul_2);
                    float ry = haarY_improved(iimage, sample_y, sample_x, int_scale_mul_2);
                    
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
        - inlined gaussian function and removed its noralization constatnt (as we normalize in the end)

    ideas:
        - keep xs_sub_sample_x and xs_sub_sample_x_squared as int and cast in expf (same for y)
    */

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale_mul_2 = (int) 2 * roundf(scale);
    float scale_squared = scale*scale;
    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

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
                    
                    float rx = haarX_improved(iimage, sample_y, sample_x, int_scale_mul_2);
                    float ry = haarY_improved(iimage, sample_y, sample_x, int_scale_mul_2);
                    
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
        get_msurf_descriptor_inlined(iimage, &interest_points->at(i));
	}
}


void get_msurf_descriptor_gauss_s1_separable_test(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - all of get_msurf_descriptor_inlined
    ideas:
        - test if gauss_s1 can be separated in x and y direction gauss_s1 = gauss_s1_x * gauss_s1_y
    */

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale_mul_2 = (int) 2 * roundf(scale);
    float scale_squared = scale*scale;
    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

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

                    // Get the separable gaussian kernels in x and y direction
                    float gauss_s1_x = expf(g1_factor * (xs_sub_sample_x_squared));
                    float gauss_s1_y = expf(g1_factor * (ys_sub_sample_y_squared));

                    // Computing gauss_s1 2d kernel weight with separable gaussian in x and y direction
                    float gauss_s1 = gauss_s1_x * gauss_s1_y;
                    
                    float rx = haarX_improved(iimage, sample_y, sample_x, int_scale_mul_2);
                    float ry = haarY_improved(iimage, sample_y, sample_x, int_scale_mul_2);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Computing gauss_s2 2d kernel weight
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}


void get_msurf_descriptors_gauss_s1_separable_test(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_gauss_s1_separable_test(iimage, &interest_points->at(i));
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


void get_msurf_descriptor_gauss_s2_precomputed(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - all of get_msurf_descriptor_inlined
        - use array gauss_s2_precomputed and index gauss_s2_index for precomputed gauss_s2 values

    ideas:
        - 
    */

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale_mul_2 = (int) 2 * roundf(scale);
    float scale_squared = scale*scale;
    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;
    
     // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        for (int j=-8; j<8; j+=5) {
            
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
                    float gauss_s1 = expf(g1_factor * (xs_sub_sample_x_squared + ys_sub_sample_y_squared));
                    
                    float rx = haarX_improved(iimage, sample_y, sample_x, int_scale_mul_2);
                    float ry = haarY_improved(iimage, sample_y, sample_x, int_scale_mul_2);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_precomputed[gauss_s2_index++];

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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}


void get_msurf_descriptors_gauss_s2_precomputed(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_gauss_s2_precomputed(iimage, &interest_points->at(i));
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

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

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
                    haarXY(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &rx, &ry);
                    // printf("sample_x,y:%i %i int_scale:%i rx:%f ry:%f\n",sample_x_sub_int_scale,sample_y_sub_int_scale, int_scale, rx,ry);

                    // float rx = haarX_improved(iimage, sample_y, sample_x, int_scale*2);
                    // float ry = haarY_improved(iimage, sample_y, sample_x, int_scale*2);
                    
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

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

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
                    haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &rx, &ry);
                    
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
                    haarXY_unconditional(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &rx, &ry);
                    
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


void get_msurf_descriptor_precompute_gauss_case(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - all of get_msurf_inlinedHaarWavelets
        - precompute gauss_s2 like get_msurf_descriptor_gauss_s2_precomputed 
          (changed to use case statements here)
        - precompute gauss_s1 by computing them once 
          with separable kernels and using case statements
    */

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale = (int) roundf(scale);
    // int int_scale_mul_2 =  2 * int_scale;
    float scale_squared = scale*scale;

    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    //float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    //r(2.5*s) + r(1.5*s) == -(r(2.5*s) - r(6.5*s))

    float s0  = roundf( 0.5 * scale);
    float s1  = roundf( 1.5 * scale);
    float s2  = roundf( 2.5 * scale);
    float s3  = roundf( 3.5 * scale);
    float s4  = roundf( 4.5 * scale);
    float s5  = roundf( 5.5 * scale);
    float s6  = roundf( 6.5 * scale);
    float s7  = roundf( 7.5 * scale);
    float s8  = roundf( 8.5 * scale);
    float s9  = roundf( 9.5 * scale);
    float s10 = roundf(10.5 * scale);
    float s11 = roundf(11.5 * scale);

    float e_c0_m4 = s2 + s1; // CAREFUL HERE!
    float e_c0_m3 = s2 + s0; // CAREFUL HERE!
    float e_c0_m2 = s2 - s0;
    float e_c0_m1 = s2 - s1;
    //float e_c0_z0 = s2 - s2;
    float e_c0_p1 = s2 - s3;
    float e_c0_p2 = s2 - s4;
    float e_c0_p3 = s2 - s5;
    float e_c0_p4 = s2 - s6;

    float e_c1_m4 = s7 - s3;
    float e_c1_m3 = s7 - s4;
    float e_c1_m2 = s7 - s5;
    float e_c1_m1 = s7 - s6;
    //float e_c1_z0 = s7 - s7;
    float e_c1_p1 = s7 - s8;
    float e_c1_p2 = s7 - s9;
    float e_c1_p3 = s7 - s10;
    float e_c1_p4 = s7 - s11;

    float gauss_s1_c0_m4 = expf(g1_factor * (e_c0_m4 * e_c0_m4));
    float gauss_s1_c0_m3 = expf(g1_factor * (e_c0_m3 * e_c0_m3));
    float gauss_s1_c0_m2 = expf(g1_factor * (e_c0_m2 * e_c0_m2));
    float gauss_s1_c0_m1 = expf(g1_factor * (e_c0_m1 * e_c0_m1));
    float gauss_s1_c0_z0 = 1.0f; //expf(g1_factor * (e_c0_z0 * e_c0_z0));
    float gauss_s1_c0_p1 = expf(g1_factor * (e_c0_p1 * e_c0_p1));
    float gauss_s1_c0_p2 = expf(g1_factor * (e_c0_p2 * e_c0_p2));
    float gauss_s1_c0_p3 = expf(g1_factor * (e_c0_p3 * e_c0_p3));
    float gauss_s1_c0_p4 = expf(g1_factor * (e_c0_p4 * e_c0_p4));

    float gauss_s1_c1_m4 = expf(g1_factor * (e_c1_m4 * e_c1_m4));
    float gauss_s1_c1_m3 = expf(g1_factor * (e_c1_m3 * e_c1_m3));
    float gauss_s1_c1_m2 = expf(g1_factor * (e_c1_m2 * e_c1_m2));
    float gauss_s1_c1_m1 = expf(g1_factor * (e_c1_m1 * e_c1_m1));
    float gauss_s1_c1_z0 = 1.0f; //expf(g1_factor * (e_c1_z0 * e_c1_z0));
    float gauss_s1_c1_p1 = expf(g1_factor * (e_c1_p1 * e_c1_p1));
    float gauss_s1_c1_p2 = expf(g1_factor * (e_c1_p2 * e_c1_p2));
    float gauss_s1_c1_p3 = expf(g1_factor * (e_c1_p3 * e_c1_p3));
    float gauss_s1_c1_p4 = expf(g1_factor * (e_c1_p4 * e_c1_p4));

    // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        for (int j=-8; j<8; j+=5) {

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            //int xs = (int) roundf(ipoint_x + (i+0.5f) * scale);
            //int ys = (int) roundf(ipoint_y + (j+0.5f) * scale);

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {

                //Get y coords of sample point
                int sample_y = (int) roundf(ipoint_y + l * scale);
                //float ys_sub_sample_y = (float) ys-sample_y;
                //float ys_sub_sample_y_squared = ys_sub_sample_y*ys_sub_sample_y;

                //Get y coords of sample point
                int sample_y_sub_int_scale = sample_y-int_scale;

                float gauss_s1_y = -1;
                if (j == -8) {
                    switch (gauss_index_l) {
                        case -4: gauss_s1_y = gauss_s1_c1_p4; break;
                        case -3: gauss_s1_y = gauss_s1_c1_p3; break;
                        case -2: gauss_s1_y = gauss_s1_c1_p2; break;
                        case -1: gauss_s1_y = gauss_s1_c1_p1; break;
                        case 0:  gauss_s1_y = gauss_s1_c1_z0; break;
                        case 1:  gauss_s1_y = gauss_s1_c1_m1; break;
                        case 2:  gauss_s1_y = gauss_s1_c1_m2; break;
                        case 3:  gauss_s1_y = gauss_s1_c1_m3; break;
                        case 4:  gauss_s1_y = gauss_s1_c1_m4; break;
                    };
                } else if (j == -3) {
                    switch (gauss_index_l) {
                        case -4: gauss_s1_y = gauss_s1_c0_p4; break;
                        case -3: gauss_s1_y = gauss_s1_c0_p3; break;
                        case -2: gauss_s1_y = gauss_s1_c0_p2; break;
                        case -1: gauss_s1_y = gauss_s1_c0_p1; break;
                        case 0:  gauss_s1_y = gauss_s1_c0_z0; break;
                        case 1:  gauss_s1_y = gauss_s1_c0_m1; break;
                        case 2:  gauss_s1_y = gauss_s1_c0_m2; break;
                        case 3:  gauss_s1_y = gauss_s1_c0_m3; break;
                        case 4:  gauss_s1_y = gauss_s1_c0_m4; break;
                    };
                } else if (j == 2) {
                    switch (gauss_index_l) {
                        case -4: gauss_s1_y = gauss_s1_c0_m4; break;
                        case -3: gauss_s1_y = gauss_s1_c0_m3; break;
                        case -2: gauss_s1_y = gauss_s1_c0_m2; break;
                        case -1: gauss_s1_y = gauss_s1_c0_m1; break;
                        case 0:  gauss_s1_y = gauss_s1_c0_z0; break;
                        case 1:  gauss_s1_y = gauss_s1_c0_p1; break;
                        case 2:  gauss_s1_y = gauss_s1_c0_p2; break;
                        case 3:  gauss_s1_y = gauss_s1_c0_p3; break;
                        case 4:  gauss_s1_y = gauss_s1_c0_p4; break;
                    };
                } else if (j == 7) {
                    switch (gauss_index_l) {
                        case -4: gauss_s1_y = gauss_s1_c1_m4; break;
                        case -3: gauss_s1_y = gauss_s1_c1_m3; break;
                        case -2: gauss_s1_y = gauss_s1_c1_m2; break;
                        case -1: gauss_s1_y = gauss_s1_c1_m1; break;
                        case 0:  gauss_s1_y = gauss_s1_c1_z0; break;
                        case 1:  gauss_s1_y = gauss_s1_c1_p1; break;
                        case 2:  gauss_s1_y = gauss_s1_c1_p2; break;
                        case 3:  gauss_s1_y = gauss_s1_c1_p3; break;
                        case 4:  gauss_s1_y = gauss_s1_c1_p4; break;
                    };
                }

                int gauss_index_k = -4;
                for (int k = i-4; k < i + 5; ++k, ++gauss_index_k) {
                
                    //Get x coords of sample point
                    int sample_x = (int) roundf(ipoint_x + k * scale);
                    //float xs_sub_sample_x = (float) xs-sample_x;
                    //float xs_sub_sample_x_squared = xs_sub_sample_x*xs_sub_sample_x;

                    //Get x coords of sample point
                    int sample_x_sub_int_scale = sample_x-int_scale;

                    float gauss_s1_x = -1;
                    if (i == -8) {
                        switch (gauss_index_k) {
                            case -4: gauss_s1_x = gauss_s1_c1_p4; break;
                            case -3: gauss_s1_x = gauss_s1_c1_p3; break;
                            case -2: gauss_s1_x = gauss_s1_c1_p2; break;
                            case -1: gauss_s1_x = gauss_s1_c1_p1; break;
                            case 0:  gauss_s1_x = gauss_s1_c1_z0; break;
                            case 1:  gauss_s1_x = gauss_s1_c1_m1; break;
                            case 2:  gauss_s1_x = gauss_s1_c1_m2; break;
                            case 3:  gauss_s1_x = gauss_s1_c1_m3; break;
                            case 4:  gauss_s1_x = gauss_s1_c1_m4; break;
                        };
                    } else if (i == -3) {
                        switch (gauss_index_k) {
                            case -4: gauss_s1_x = gauss_s1_c0_p4; break;
                            case -3: gauss_s1_x = gauss_s1_c0_p3; break;
                            case -2: gauss_s1_x = gauss_s1_c0_p2; break;
                            case -1: gauss_s1_x = gauss_s1_c0_p1; break;
                            case 0:  gauss_s1_x = gauss_s1_c0_z0; break;
                            case 1:  gauss_s1_x = gauss_s1_c0_m1; break;
                            case 2:  gauss_s1_x = gauss_s1_c0_m2; break;
                            case 3:  gauss_s1_x = gauss_s1_c0_m3; break;
                            case 4:  gauss_s1_x = gauss_s1_c0_m4; break;
                        };
                    } else if (i == 2) {
                        switch (gauss_index_k) {
                            case -4: gauss_s1_x = gauss_s1_c0_m4; break;
                            case -3: gauss_s1_x = gauss_s1_c0_m3; break;
                            case -2: gauss_s1_x = gauss_s1_c0_m2; break;
                            case -1: gauss_s1_x = gauss_s1_c0_m1; break;
                            case 0:  gauss_s1_x = gauss_s1_c0_z0; break;
                            case 1:  gauss_s1_x = gauss_s1_c0_p1; break;
                            case 2:  gauss_s1_x = gauss_s1_c0_p2; break;
                            case 3:  gauss_s1_x = gauss_s1_c0_p3; break;
                            case 4:  gauss_s1_x = gauss_s1_c0_p4; break;
                        };
                    } else if (i == 7) {
                        switch (gauss_index_k) {
                            case -4: gauss_s1_x = gauss_s1_c1_m4; break;
                            case -3: gauss_s1_x = gauss_s1_c1_m3; break;
                            case -2: gauss_s1_x = gauss_s1_c1_m2; break;
                            case -1: gauss_s1_x = gauss_s1_c1_m1; break;
                            case 0:  gauss_s1_x = gauss_s1_c1_z0; break;
                            case 1:  gauss_s1_x = gauss_s1_c1_p1; break;
                            case 2:  gauss_s1_x = gauss_s1_c1_p2; break;
                            case 3:  gauss_s1_x = gauss_s1_c1_p3; break;
                            case 4:  gauss_s1_x = gauss_s1_c1_p4; break;
                        };
                    }

/*
                    //Get the gaussian weighted x and y responses
                    float gauss_s1_x_real = expf(g1_factor * (xs_sub_sample_x_squared));
                    float gauss_s1_y_real = expf(g1_factor * (ys_sub_sample_y_squared));

                    if (gauss_s1_x != gauss_s1_x_real || gauss_s1_y != gauss_s1_y_real) {
                        std::cout << "gauss_s1_x: " << gauss_s1_x << std::endl;
                        std::cout << "gauss_s1_y: " << gauss_s1_y << std::endl;
                        std::cout << "gauss_s1_x_real: " << gauss_s1_x_real << std::endl;
                        std::cout << "gauss_s1_y_real: " << gauss_s1_y_real << std::endl;
                    }
*/

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = 0.0f;
                    float ry = 0.0f;
                    haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &rx, &ry);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            //float gauss_s2 = expf(g2_factor * (cx_squared + cy_squared)); 
            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            //float gauss_s2 = gauss_s2_precomputed[gauss_s2_index++];
            float gauss_s2;
            switch (gauss_s2_index) {
                case 0:  gauss_s2 = 0.026022f; break;
                case 1:  gauss_s2 = 0.040585f; break;
                case 2:  gauss_s2 = 0.040585f; break;
                case 3:  gauss_s2 = 0.026022f; break;
                case 4:  gauss_s2 = 0.040585f; break;
                case 5:  gauss_s2 = 0.063297f; break;
                case 6:  gauss_s2 = 0.063297f; break;
                case 7:  gauss_s2 = 0.040585f; break;
                case 8:  gauss_s2 = 0.040585f; break;
                case 9:  gauss_s2 = 0.063297f; break;
                case 10: gauss_s2 = 0.063297f; break;
                case 11: gauss_s2 = 0.040585f; break;
                case 12: gauss_s2 = 0.026022f; break;
                case 13: gauss_s2 = 0.040585f; break;
                case 14: gauss_s2 = 0.040585f; break;
                case 15: gauss_s2 = 0.026022f; break;
            };
            gauss_s2_index++;

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


void get_msurf_descriptors_precompute_gauss_case(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_precompute_gauss_case(iimage, &interest_points->at(i));
	}
}

float gauss_s1_c0[9] __attribute__((aligned(64)));
float gauss_s1_c1[9] __attribute__((aligned(64)));

extern const float gauss_s2_arr[16] __attribute__((aligned(64))) = {0.026022f, 0.040585f, 0.040585f, 0.026022f, 
                                        0.040585f, 0.063297f, 0.063297f, 0.040585f, 
                                        0.040585f, 0.063297f, 0.063297f, 0.040585f, 
                                        0.026022f, 0.040585f, 0.040585f, 0.026022f};

void get_msurf_descriptor_precompute_gauss_array(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - equal to get_msurf_descriptor_gauss_compute_once_case but with arrays insteat of case statements
    */

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale = (int) roundf(scale);
    // int int_scale_mul_2 =  2 * int_scale;
    float scale_squared = scale*scale;

    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    //float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    //r(2.5*s) + r(1.5*s) == -(r(2.5*s) - r(6.5*s))

    float s0  = roundf( 0.5 * scale);
    float s1  = roundf( 1.5 * scale);
    float s2  = roundf( 2.5 * scale);
    float s3  = roundf( 3.5 * scale);
    float s4  = roundf( 4.5 * scale);
    float s5  = roundf( 5.5 * scale);
    float s6  = roundf( 6.5 * scale);
    float s7  = roundf( 7.5 * scale);
    float s8  = roundf( 8.5 * scale);
    float s9  = roundf( 9.5 * scale);
    float s10 = roundf(10.5 * scale);
    float s11 = roundf(11.5 * scale);

    float e_c0_m4 = s2 + s1; // CAREFUL HERE!
    float e_c0_m3 = s2 + s0; // CAREFUL HERE!
    float e_c0_m2 = s2 - s0;
    float e_c0_m1 = s2 - s1;
    //float e_c0_z0 = s2 - s2;
    float e_c0_p1 = s2 - s3;
    float e_c0_p2 = s2 - s4;
    float e_c0_p3 = s2 - s5;
    float e_c0_p4 = s2 - s6;

    float e_c1_m4 = s7 - s3;
    float e_c1_m3 = s7 - s4;
    float e_c1_m2 = s7 - s5;
    float e_c1_m1 = s7 - s6;
    //float e_c1_z0 = s7 - s7;
    float e_c1_p1 = s7 - s8;
    float e_c1_p2 = s7 - s9;
    float e_c1_p3 = s7 - s10;
    float e_c1_p4 = s7 - s11;

    gauss_s1_c0[0] =  expf(g1_factor * (e_c0_m4 * e_c0_m4));
    gauss_s1_c0[1] =  expf(g1_factor * (e_c0_m3 * e_c0_m3));
    gauss_s1_c0[2] =  expf(g1_factor * (e_c0_m2 * e_c0_m2));
    gauss_s1_c0[3] =  expf(g1_factor * (e_c0_m1 * e_c0_m1));
    gauss_s1_c0[4] =  1.0f; //expf(g1_factor * (e_c0_z0 * e_c0_z0));
    gauss_s1_c0[5] =  expf(g1_factor * (e_c0_p1 * e_c0_p1));
    gauss_s1_c0[6] =  expf(g1_factor * (e_c0_p2 * e_c0_p2));
    gauss_s1_c0[7] =  expf(g1_factor * (e_c0_p3 * e_c0_p3));
    gauss_s1_c0[8] =  expf(g1_factor * (e_c0_p4 * e_c0_p4));

    gauss_s1_c1[0] =  expf(g1_factor * (e_c1_m4 * e_c1_m4));
    gauss_s1_c1[1] =  expf(g1_factor * (e_c1_m3 * e_c1_m3));
    gauss_s1_c1[2] =  expf(g1_factor * (e_c1_m2 * e_c1_m2));
    gauss_s1_c1[3] =  expf(g1_factor * (e_c1_m1 * e_c1_m1));
    gauss_s1_c1[4] =  1.0f; //expf(g1_factor * (e_c1_z0 * e_c1_z0));
    gauss_s1_c1[5] =  expf(g1_factor * (e_c1_p1 * e_c1_p1));
    gauss_s1_c1[6] =  expf(g1_factor * (e_c1_p2 * e_c1_p2));
    gauss_s1_c1[7] =  expf(g1_factor * (e_c1_p3 * e_c1_p3));
    gauss_s1_c1[8] =  expf(g1_factor * (e_c1_p4 * e_c1_p4));

    // calculate descriptor for this interest point

    for (int i=-8; i<8; i+=5) {

        float *i_arr;
        if (i == -8 || i == 7) {
            i_arr = gauss_s1_c1;
        } else {
            i_arr = gauss_s1_c0;
        } 

        int gauss_index_k_start = (i<0?8:0);
        int gauss_index_k_inc = (i<0?-1:1);

        for (int j=-8; j<8; j+=5) {

            float *j_arr;
            if (j == -8 || j == 7) {
                j_arr = gauss_s1_c1;
            } else {
                j_arr = gauss_s1_c0;
            } 

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int gauss_index_l = (j<0?8:0);
            int gauss_index_l_inc = (j<0?-1:1);
            for (int l = (j+8); l < (j + 17); ++l, gauss_index_l+=gauss_index_l_inc) {

                int sample_y = (int) roundf(ipoint_y + (l-12) * scale);
                int sample_y_sub_int_scale = sample_y-int_scale;

                float gauss_s1_y = j_arr[gauss_index_l];

                int gauss_index_k = gauss_index_k_start;
                for (int k = i+8; k < i + 17; ++k, gauss_index_k+=gauss_index_k_inc) {

                    int sample_x = (int) roundf(ipoint_x + (k-12) * scale);
                    int sample_x_sub_int_scale = sample_x-int_scale;
                    
                    float gauss_s1_x = i_arr[gauss_index_k];

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = 0.0f;
                    float ry = 0.0f;
                    haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &rx, &ry);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_arr[gauss_s2_index];
            gauss_s2_index++;

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

    /* Code pre optimization
    for (int i=-8; i<8; i+=5) {

        for (int j=-8; j<8; j+=5) {

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            //int xs = (int) roundf(ipoint_x + (i+0.5f) * scale);
            //int ys = (int) roundf(ipoint_y + (j+0.5f) * scale);

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {

                //Get y coords of sample point
                int sample_y = (int) roundf(ipoint_y + l * scale);
                //float ys_sub_sample_y = (float) ys-sample_y;
                //float ys_sub_sample_y_squared = ys_sub_sample_y*ys_sub_sample_y;

                //Get y coords of sample point
                int sample_y_sub_int_scale = sample_y-int_scale;

                float gauss_s1_y = -1;
                if (j == -8 ) {
                    gauss_s1_y = gauss_s1_c1[8-(gauss_index_l+4)];
                } else if (j == -3) {
                    gauss_s1_y = gauss_s1_c0[8-(gauss_index_l+4)];
                } else if (j == 2) {
                    gauss_s1_y = gauss_s1_c0[gauss_index_l+4];
                } else if (j == 7) {
                    gauss_s1_y = gauss_s1_c1[gauss_index_l+4];
                }

                int gauss_index_k = -4;
                for (int k = i-4; k < i + 5; ++k, ++gauss_index_k) {
                
                    //Get x coords of sample point
                    int sample_x = (int) roundf(ipoint_x + k * scale);
                    //float xs_sub_sample_x = (float) xs-sample_x;
                    //float xs_sub_sample_x_squared = xs_sub_sample_x*xs_sub_sample_x;

                    //Get x coords of sample point
                    int sample_x_sub_int_scale = sample_x-int_scale;

                    float gauss_s1_x = -1;
                    if (i == -8 ) {
                        gauss_s1_x = gauss_s1_c1[8-(gauss_index_k+4)];
                    } else if (i == -3) {
                        gauss_s1_x = gauss_s1_c0[8-(gauss_index_k+4)];
                    } else if (i == 2) {
                        gauss_s1_x = gauss_s1_c0[gauss_index_k+4];
                    } else if (i == 7) {
                        gauss_s1_x = gauss_s1_c1[gauss_index_k+4];
                    }

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = 0.0f;
                    float ry = 0.0f;
                    haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &rx, &ry);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            //float gauss_s2 = expf(g2_factor * (cx_squared + cy_squared)); 
            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            //float gauss_s2 = gauss_s2_precomputed[gauss_s2_index++];
            float gauss_s2 = gauss_s2_arr[gauss_s2_index];
            gauss_s2_index++;

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
    */
    // rescale to unit vector
    // NOTE: using sqrtf() for floats
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}


void get_msurf_descriptors_precompute_gauss_array(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_precompute_gauss_array(iimage, &interest_points->at(i));
	}
}



// float haarResponseX[24*24];
// float haarResponseY[24*24];

float * haarResponseX = (float*) aligned_malloc(24*24*sizeof(float), 64);
float * haarResponseY = (float*) aligned_malloc(24*24*sizeof(float), 64);

void get_msurf_descriptor_pecompute_haar(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - all of get_msurf_descriptor_gauss_compute_once_array
        - precompute haar wavelet responses and store them in two 24x24 arrays
    */

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale = (int) roundf(scale);
    // int int_scale_mul_2 =  2 * int_scale;
    float scale_squared = scale*scale;

    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    //float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

    float ipoint_x_sub_int_scale = ipoint_x-int_scale;
    float ipoint_y_sub_int_scale = ipoint_y-int_scale;

    int width = iimage->width;
    int height = iimage->height;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    //r(2.5*s) + r(1.5*s) == -(r(2.5*s) - r(6.5*s))

    // check if we ever hit a boundary
    if (((int) roundf(ipoint_x - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_y - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_x + 11*scale)) + int_scale > width 
        || ((int) roundf(ipoint_y + 11*scale)) + int_scale > height) 
    { // some outside

        for (int l=-12, l_count=0; l<12; ++l, ++l_count) {
            // int sample_y = (int) roundf(ipoint_y + l * scale);
            int sample_y_sub_int_scale = (int) roundf(ipoint_y_sub_int_scale + l * scale);

            for (int k=-12, k_count=0; k<12; ++k, k_count++) {

                //Get x coords of sample point
                // int sample_x = (int) roundf(ipoint_x + k * scale);
                int sample_x_sub_int_scale = (int) roundf(ipoint_x_sub_int_scale + k * scale);;

                // float rx = 0.0f;
                // float ry = 0.0f;
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &haarResponseX[l_count*24+k_count], &haarResponseY[l_count*24+k_count]);

                // haarResponseX[(l+12)*24+(k+12)] = rx;
                // haarResponseY[(l+12)*24+(k+12)] = ry;
            }

        }

    } else {

        for (int l=-12, l_count=0; l<12; ++l, ++l_count) {
            // int sample_y = (int) roundf(ipoint_y + l * scale);
            int sample_y_sub_int_scale = (int) roundf(ipoint_y_sub_int_scale + l * scale);

            for (int k=-12, k_count=0; k<12; ++k, ++k_count) {

                //Get x coords of sample point
                // int sample_x = (int) roundf(ipoint_x + k * scale);
                int sample_x_sub_int_scale = (int) roundf(ipoint_x_sub_int_scale + k * scale);;

                // float rx = 0.0f;
                // float ry = 0.0f;
                haarXY_unconditional(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &haarResponseX[l_count*24+k_count], &haarResponseY[l_count*24+k_count]);

                // haarResponseX[(l+12)*24+(k+12)] = rx;
                // haarResponseY[(l+12)*24+(k+12)] = ry;
            }

        }
        
    }

    float s0  = roundf( 0.5 * scale);
    float s1  = roundf( 1.5 * scale);
    float s2  = roundf( 2.5 * scale);
    float s3  = roundf( 3.5 * scale);
    float s4  = roundf( 4.5 * scale);
    float s5  = roundf( 5.5 * scale);
    float s6  = roundf( 6.5 * scale);
    float s7  = roundf( 7.5 * scale);
    float s8  = roundf( 8.5 * scale);
    float s9  = roundf( 9.5 * scale);
    float s10 = roundf(10.5 * scale);
    float s11 = roundf(11.5 * scale);

    float e_c0_m4 = s2 + s1; // CAREFUL HERE!
    float e_c0_m3 = s2 + s0; // CAREFUL HERE!
    float e_c0_m2 = s2 - s0;
    float e_c0_m1 = s2 - s1;
    //float e_c0_z0 = s2 - s2;
    float e_c0_p1 = s2 - s3;
    float e_c0_p2 = s2 - s4;
    float e_c0_p3 = s2 - s5;
    float e_c0_p4 = s2 - s6;

    float e_c1_m4 = s7 - s3;
    float e_c1_m3 = s7 - s4;
    float e_c1_m2 = s7 - s5;
    float e_c1_m1 = s7 - s6;
    //float e_c1_z0 = s7 - s7;
    float e_c1_p1 = s7 - s8;
    float e_c1_p2 = s7 - s9;
    float e_c1_p3 = s7 - s10;
    float e_c1_p4 = s7 - s11;

    gauss_s1_c0[0] =  expf(g1_factor * (e_c0_m4 * e_c0_m4));
    gauss_s1_c0[1] =  expf(g1_factor * (e_c0_m3 * e_c0_m3));
    gauss_s1_c0[2] =  expf(g1_factor * (e_c0_m2 * e_c0_m2));
    gauss_s1_c0[3] =  expf(g1_factor * (e_c0_m1 * e_c0_m1));
    gauss_s1_c0[4] =  1.0f; //expf(g1_factor * (e_c0_z0 * e_c0_z0));
    gauss_s1_c0[5] =  expf(g1_factor * (e_c0_p1 * e_c0_p1));
    gauss_s1_c0[6] =  expf(g1_factor * (e_c0_p2 * e_c0_p2));
    gauss_s1_c0[7] =  expf(g1_factor * (e_c0_p3 * e_c0_p3));
    gauss_s1_c0[8] =  expf(g1_factor * (e_c0_p4 * e_c0_p4));

    gauss_s1_c1[0] =  expf(g1_factor * (e_c1_m4 * e_c1_m4));
    gauss_s1_c1[1] =  expf(g1_factor * (e_c1_m3 * e_c1_m3));
    gauss_s1_c1[2] =  expf(g1_factor * (e_c1_m2 * e_c1_m2));
    gauss_s1_c1[3] =  expf(g1_factor * (e_c1_m1 * e_c1_m1));
    gauss_s1_c1[4] =  1.0f; //expf(g1_factor * (e_c1_z0 * e_c1_z0));
    gauss_s1_c1[5] =  expf(g1_factor * (e_c1_p1 * e_c1_p1));
    gauss_s1_c1[6] =  expf(g1_factor * (e_c1_p2 * e_c1_p2));
    gauss_s1_c1[7] =  expf(g1_factor * (e_c1_p3 * e_c1_p3));
    gauss_s1_c1[8] =  expf(g1_factor * (e_c1_p4 * e_c1_p4));

    // calculate descriptor for this interest point
    /*
    for (int i=-8; i<8; i+=5) {

        for (int j=-8; j<8; j+=5) {

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            //int xs = (int) roundf(ipoint_x + (i+0.5f) * scale);
            //int ys = (int) roundf(ipoint_y + (j+0.5f) * scale);

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {

                //Get y coords of sample point
                // int sample_y = (int) roundf(ipoint_y + l * scale);
                //float ys_sub_sample_y = (float) ys-sample_y;
                //float ys_sub_sample_y_squared = ys_sub_sample_y*ys_sub_sample_y;

                //Get y coords of sample point
                // int sample_y_sub_int_scale = sample_y-int_scale;

                float gauss_s1_y = -1;
                if (j == -8 ) {
                    gauss_s1_y = gauss_s1_c1[8-(gauss_index_l+4)];
                } else if (j == -3) {
                    gauss_s1_y = gauss_s1_c0[8-(gauss_index_l+4)];
                } else if (j == 2) {
                    gauss_s1_y = gauss_s1_c0[gauss_index_l+4];
                } else if (j == 7) {
                    gauss_s1_y = gauss_s1_c1[gauss_index_l+4];
                }

                int gauss_index_k = -4;
                for (int k = i-4; k < i + 5; ++k, ++gauss_index_k) {
                
                    //Get x coords of sample point
                    // int sample_x = (int) roundf(ipoint_x + k * scale);

                    //Get x coords of sample point
                    // int sample_x_sub_int_scale = sample_x-int_scale;

                    float gauss_s1_x = -1;
                    if (i == -8 ) {
                        gauss_s1_x = gauss_s1_c1[8-(gauss_index_k+4)];
                    } else if (i == -3) {
                        gauss_s1_x = gauss_s1_c0[8-(gauss_index_k+4)];
                    } else if (i == 2) {
                        gauss_s1_x = gauss_s1_c0[gauss_index_k+4];
                    } else if (i == 7) {
                        gauss_s1_x = gauss_s1_c1[gauss_index_k+4];
                    }

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = haarResponseX[(l+12)*24+(k+12)];
                    float ry = haarResponseY[(l+12)*24+(k+12)];
                    // haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &rx, &ry);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_arr[gauss_s2_index];
            gauss_s2_index++;

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
    */
    for (int i=-8; i<8; i+=5) {

        float *i_arr;
        if (i == -8 || i == 7) {
            i_arr = gauss_s1_c1;
        } else {
            i_arr = gauss_s1_c0;
        } 

        int gauss_index_k_start = (i<0?8:0);
        int gauss_index_k_inc = (i<0?-1:1);

        for (int j=-8; j<8; j+=5) {

            float *j_arr;
            if (j == -8 || j == 7) {
                j_arr = gauss_s1_c1;
            } else {
                j_arr = gauss_s1_c0;
            } 

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int gauss_index_l = (j<0?8:0);
            int gauss_index_l_inc = (j<0?-1:1);
            for (int l = (j+8)*24; l < (j + 17)*24; l+=24, gauss_index_l+=gauss_index_l_inc) {
                float gauss_s1_y = j_arr[gauss_index_l];

                int gauss_index_k = gauss_index_k_start;
                for (int k = i+8; k < i + 17; ++k, gauss_index_k+=gauss_index_k_inc) {
                    
                    float gauss_s1_x = i_arr[gauss_index_k];

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = haarResponseX[l+k];
                    float ry = haarResponseY[l+k];
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_arr[gauss_s2_index];
            gauss_s2_index++;

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


void get_msurf_descriptors_pecompute_haar(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_pecompute_haar(iimage, &interest_points->at(i));
	}
}


void get_msurf_descriptor_rounding(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - all of get_msurf_descriptors_pecompute_haar
        - simplify rounding as (int) (x+0.5) if x>=0 and (int) (x-0.5) oterwise where return of roundf is int
    */

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale = (int) roundf(scale);
    // int int_scale_mul_2 =  2 * int_scale;
    float scale_squared = scale*scale;

    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    //float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

    float ipoint_x_sub_int_scale = ipoint_x-int_scale;
    float ipoint_y_sub_int_scale = ipoint_y-int_scale;

    float ipoint_x_sub_int_scale_add_05 = ipoint_x-int_scale + 0.5;
    float ipoint_y_sub_int_scale_add_05 = ipoint_y-int_scale + 0.5;

    int width = iimage->width;
    int height = iimage->height;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    //r(2.5*s) + r(1.5*s) == -(r(2.5*s) - r(6.5*s))

    // check if we ever hit a boundary
    if (((int) roundf(ipoint_x - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_y - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_x + 11*scale)) + int_scale > width 
        || ((int) roundf(ipoint_y + 11*scale)) + int_scale > height) 
    { // some outside

        for (int l=-12, l_count=0; l<12; ++l, ++l_count) {
            // int sample_y = (int) roundf(ipoint_y + l * scale);
            float ipoint_y_sub_int_scale_add_l_mul_scale = ipoint_y_sub_int_scale + l * scale;
            int sample_y_sub_int_scale = (int) (ipoint_y_sub_int_scale_add_l_mul_scale + (ipoint_y_sub_int_scale_add_l_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; ++k, k_count++) {

                //Get x coords of sample point
                // int sample_x = (int) roundf(ipoint_x + k * scale);
                float ipoint_x_sub_int_scale_add_k_mul_scale = ipoint_x_sub_int_scale + k * scale;
                int sample_x_sub_int_scale = (int) (ipoint_x_sub_int_scale_add_k_mul_scale + (ipoint_x_sub_int_scale_add_k_mul_scale>=0 ? 0.5 : -0.5));

                // float rx = 0.0f;
                // float ry = 0.0f;
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &haarResponseX[l_count*24+k_count], &haarResponseY[l_count*24+k_count]);

                // haarResponseX[(l+12)*24+(k+12)] = rx;
                // haarResponseY[(l+12)*24+(k+12)] = ry;
            }

        }

    } else {

        for (int l=-12, l_count=0; l<12; ++l, ++l_count) {
            // int sample_y = (int) roundf(ipoint_y + l * scale);
            int sample_y_sub_int_scale = (int)(ipoint_y_sub_int_scale_add_05 + l * scale);

            for (int k=-12, k_count=0; k<12; ++k, ++k_count) {

                //Get x coords of sample point
                // int sample_x = (int) roundf(ipoint_x + k * scale);
                int sample_x_sub_int_scale = (int)(ipoint_x_sub_int_scale_add_05 + k * scale);

                // float rx = 0.0f;
                // float ry = 0.0f;
                haarXY_unconditional(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &haarResponseX[l_count*24+k_count], &haarResponseY[l_count*24+k_count]);

                // haarResponseX[(l+12)*24+(k+12)] = rx;
                // haarResponseY[(l+12)*24+(k+12)] = ry;
            }

        }
        
    }

    float s0  = roundf( 0.5 * scale);
    float s1  = roundf( 1.5 * scale);
    float s2  = roundf( 2.5 * scale);
    float s3  = roundf( 3.5 * scale);
    float s4  = roundf( 4.5 * scale);
    float s5  = roundf( 5.5 * scale);
    float s6  = roundf( 6.5 * scale);
    float s7  = roundf( 7.5 * scale);
    float s8  = roundf( 8.5 * scale);
    float s9  = roundf( 9.5 * scale);
    float s10 = roundf(10.5 * scale);
    float s11 = roundf(11.5 * scale);

    float e_c0_m4 = s2 + s1; // CAREFUL HERE!
    float e_c0_m3 = s2 + s0; // CAREFUL HERE!
    float e_c0_m2 = s2 - s0;
    float e_c0_m1 = s2 - s1;
    //float e_c0_z0 = s2 - s2;
    float e_c0_p1 = s2 - s3;
    float e_c0_p2 = s2 - s4;
    float e_c0_p3 = s2 - s5;
    float e_c0_p4 = s2 - s6;

    float e_c1_m4 = s7 - s3;
    float e_c1_m3 = s7 - s4;
    float e_c1_m2 = s7 - s5;
    float e_c1_m1 = s7 - s6;
    //float e_c1_z0 = s7 - s7;
    float e_c1_p1 = s7 - s8;
    float e_c1_p2 = s7 - s9;
    float e_c1_p3 = s7 - s10;
    float e_c1_p4 = s7 - s11;

    gauss_s1_c0[0] =  expf(g1_factor * (e_c0_m4 * e_c0_m4));
    gauss_s1_c0[1] =  expf(g1_factor * (e_c0_m3 * e_c0_m3));
    gauss_s1_c0[2] =  expf(g1_factor * (e_c0_m2 * e_c0_m2));
    gauss_s1_c0[3] =  expf(g1_factor * (e_c0_m1 * e_c0_m1));
    gauss_s1_c0[4] =  1.0f; //expf(g1_factor * (e_c0_z0 * e_c0_z0));
    gauss_s1_c0[5] =  expf(g1_factor * (e_c0_p1 * e_c0_p1));
    gauss_s1_c0[6] =  expf(g1_factor * (e_c0_p2 * e_c0_p2));
    gauss_s1_c0[7] =  expf(g1_factor * (e_c0_p3 * e_c0_p3));
    gauss_s1_c0[8] =  expf(g1_factor * (e_c0_p4 * e_c0_p4));

    gauss_s1_c1[0] =  expf(g1_factor * (e_c1_m4 * e_c1_m4));
    gauss_s1_c1[1] =  expf(g1_factor * (e_c1_m3 * e_c1_m3));
    gauss_s1_c1[2] =  expf(g1_factor * (e_c1_m2 * e_c1_m2));
    gauss_s1_c1[3] =  expf(g1_factor * (e_c1_m1 * e_c1_m1));
    gauss_s1_c1[4] =  1.0f; //expf(g1_factor * (e_c1_z0 * e_c1_z0));
    gauss_s1_c1[5] =  expf(g1_factor * (e_c1_p1 * e_c1_p1));
    gauss_s1_c1[6] =  expf(g1_factor * (e_c1_p2 * e_c1_p2));
    gauss_s1_c1[7] =  expf(g1_factor * (e_c1_p3 * e_c1_p3));
    gauss_s1_c1[8] =  expf(g1_factor * (e_c1_p4 * e_c1_p4));

    // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        float *i_arr;
        if (i == -8 || i == 7) {
            i_arr = gauss_s1_c1;
        } else {
            i_arr = gauss_s1_c0;
        } 

        int gauss_index_k_start = (i<0?8:0);
        int gauss_index_k_inc = (i<0?-1:1);

        for (int j=-8; j<8; j+=5) {

            float *j_arr;
            if (j == -8 || j == 7) {
                j_arr = gauss_s1_c1;
            } else {
                j_arr = gauss_s1_c0;
            } 

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int gauss_index_l = (j<0?8:0);
            int gauss_index_l_inc = (j<0?-1:1);
            for (int l = (j+8)*24; l < (j + 17)*24; l+=24, gauss_index_l+=gauss_index_l_inc) {
                float gauss_s1_y = j_arr[gauss_index_l];

                int gauss_index_k = gauss_index_k_start;
                for (int k = i+8; k < i + 17; ++k, gauss_index_k+=gauss_index_k_inc) {
                    
                    float gauss_s1_x = i_arr[gauss_index_k];

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = haarResponseX[l+k];
                    float ry = haarResponseY[l+k];
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_arr[gauss_s2_index];
            gauss_s2_index++;

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


void get_msurf_descriptors_rounding(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_rounding(iimage, &interest_points->at(i));
	}
}


void get_msurf_descriptor_rounding_unconditional(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    applied optimizations:
        - all of get_msurf_descriptor_gauss_pecompute_haar_rounding
        - replaced case statements with static arrays
    */

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale = (int) roundf(scale);
    // int int_scale_mul_2 =  2 * int_scale;
    float scale_squared = scale*scale;

    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    //float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

    float ipoint_x_sub_int_scale = ipoint_x-int_scale;
    float ipoint_y_sub_int_scale = ipoint_y-int_scale;

    float ipoint_x_sub_int_scale_add_05 = ipoint_x-int_scale + 0.5;
    float ipoint_y_sub_int_scale_add_05 = ipoint_y-int_scale + 0.5;
    
    int width = iimage->width;
    int height = iimage->height;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    // check if we ever hit a boundary

    for (int l=-12, l_count=0; l<12; ++l, ++l_count) {
        // int sample_y = (int) roundf(ipoint_y + l * scale);
        int sample_y_sub_int_scale = (int) roundf(ipoint_y_sub_int_scale + l * scale);// + (ipoint_y_sub_int_scale_add_l_mul_scale>=0 ? 0.5 : -0.5));

        for (int k=-12, k_count=0; k<12; ++k, k_count++) {
            int sample_x_sub_int_scale = (int) roundf(ipoint_x_sub_int_scale + k * scale);// + (ipoint_x_sub_int_scale_add_k_mul_scale>=0 ? 0.5 : -0.5));

            haarXY_unconditional(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &haarResponseX[l_count*24+k_count], &haarResponseY[l_count*24+k_count]);

        }

    }
        

    float s0  = roundf( 0.5 * scale);
    float s1  = roundf( 1.5 * scale);
    float s2  = roundf( 2.5 * scale);
    float s3  = roundf( 3.5 * scale);
    float s4  = roundf( 4.5 * scale);
    float s5  = roundf( 5.5 * scale);
    float s6  = roundf( 6.5 * scale);
    float s7  = roundf( 7.5 * scale);
    float s8  = roundf( 8.5 * scale);
    float s9  = roundf( 9.5 * scale);
    float s10 = roundf(10.5 * scale);
    float s11 = roundf(11.5 * scale);

    float e_c0_m4 = s2 + s1; // CAREFUL HERE!
    float e_c0_m3 = s2 + s0; // CAREFUL HERE!
    float e_c0_m2 = s2 - s0;
    float e_c0_m1 = s2 - s1;
    //float e_c0_z0 = s2 - s2;
    float e_c0_p1 = s2 - s3;
    float e_c0_p2 = s2 - s4;
    float e_c0_p3 = s2 - s5;
    float e_c0_p4 = s2 - s6;

    float e_c1_m4 = s7 - s3;
    float e_c1_m3 = s7 - s4;
    float e_c1_m2 = s7 - s5;
    float e_c1_m1 = s7 - s6;
    //float e_c1_z0 = s7 - s7;
    float e_c1_p1 = s7 - s8;
    float e_c1_p2 = s7 - s9;
    float e_c1_p3 = s7 - s10;
    float e_c1_p4 = s7 - s11;

    gauss_s1_c0[0] =  expf(g1_factor * (e_c0_m4 * e_c0_m4));
    gauss_s1_c0[1] =  expf(g1_factor * (e_c0_m3 * e_c0_m3));
    gauss_s1_c0[2] =  expf(g1_factor * (e_c0_m2 * e_c0_m2));
    gauss_s1_c0[3] =  expf(g1_factor * (e_c0_m1 * e_c0_m1));
    gauss_s1_c0[4] =  1.0f; //expf(g1_factor * (e_c0_z0 * e_c0_z0));
    gauss_s1_c0[5] =  expf(g1_factor * (e_c0_p1 * e_c0_p1));
    gauss_s1_c0[6] =  expf(g1_factor * (e_c0_p2 * e_c0_p2));
    gauss_s1_c0[7] =  expf(g1_factor * (e_c0_p3 * e_c0_p3));
    gauss_s1_c0[8] =  expf(g1_factor * (e_c0_p4 * e_c0_p4));

    gauss_s1_c1[0] =  expf(g1_factor * (e_c1_m4 * e_c1_m4));
    gauss_s1_c1[1] =  expf(g1_factor * (e_c1_m3 * e_c1_m3));
    gauss_s1_c1[2] =  expf(g1_factor * (e_c1_m2 * e_c1_m2));
    gauss_s1_c1[3] =  expf(g1_factor * (e_c1_m1 * e_c1_m1));
    gauss_s1_c1[4] =  1.0f; //expf(g1_factor * (e_c1_z0 * e_c1_z0));
    gauss_s1_c1[5] =  expf(g1_factor * (e_c1_p1 * e_c1_p1));
    gauss_s1_c1[6] =  expf(g1_factor * (e_c1_p2 * e_c1_p2));
    gauss_s1_c1[7] =  expf(g1_factor * (e_c1_p3 * e_c1_p3));
    gauss_s1_c1[8] =  expf(g1_factor * (e_c1_p4 * e_c1_p4));

    // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        float *i_arr;
        if (i == -8 || i == 7) {
            i_arr = gauss_s1_c1;
        } else {
            i_arr = gauss_s1_c0;
        } 

        int gauss_index_k_start = (i<0?8:0);
        int gauss_index_k_inc = (i<0?-1:1);

        for (int j=-8; j<8; j+=5) {

            float *j_arr;
            if (j == -8 || j == 7) {
                j_arr = gauss_s1_c1;
            } else {
                j_arr = gauss_s1_c0;
            } 

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int gauss_index_l = (j<0?8:0);
            int gauss_index_l_inc = (j<0?-1:1);
            for (int l = (j+8)*24; l < (j + 17)*24; l+=24, gauss_index_l+=gauss_index_l_inc) {
                float gauss_s1_y = j_arr[gauss_index_l];

                int gauss_index_k = gauss_index_k_start;
                for (int k = i+8; k < i + 17; ++k, gauss_index_k+=gauss_index_k_inc) {
                    
                    float gauss_s1_x = i_arr[gauss_index_k];

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = haarResponseX[l+k];
                    float ry = haarResponseY[l+k];
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_arr[gauss_s2_index];
            gauss_s2_index++;

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


void get_msurf_descriptors_rounding_unconditional(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_rounding_unconditional(iimage, &interest_points->at(i));
	}
}

 
void get_msurf_descriptor_gauss_pecompute_haar_unroll(struct integral_image* iimage, struct interest_point* ipoint) {
    /*
    Obsolete due to autotuning
    applied optimizations:
        - all of get_msurf_descriptors_gauss_pecompute_haar
        - unrolling for the precomputation
    */

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale = (int) roundf(scale);
    // int int_scale_mul_2 =  2 * int_scale;
    float scale_squared = scale*scale;

    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    //float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

    float ipoint_x_sub_int_scale = ipoint_x-int_scale;
    float ipoint_y_sub_int_scale = ipoint_y-int_scale;

    int width = iimage->width;
    int height = iimage->height;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    //r(2.5*s) + r(1.5*s) == -(r(2.5*s) - r(6.5*s))

    // check if we ever hit a boundary
    if (((int) roundf(ipoint_x - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_y - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_x + 11*scale)) + int_scale > width 
        || ((int) roundf(ipoint_y + 11*scale)) + int_scale > height) 
    { // some outside

        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            // int sample_y = (int) roundf(ipoint_y + l * scale);
            int l0 = l;
            int l1 = l+1;
            int l2 = l+2;
            int l3 = l+3;
            int l_count0 = l_count;
            int l_count1 = l_count+1;
            int l_count2 = l_count+2;
            int l_count3 = l_count+3;


            int sample_y_sub_int_scale0 = (int) roundf(ipoint_y_sub_int_scale + l0 * scale);
            int sample_y_sub_int_scale1 = (int) roundf(ipoint_y_sub_int_scale + l1 * scale);
            int sample_y_sub_int_scale2 = (int) roundf(ipoint_y_sub_int_scale + l2 * scale);
            int sample_y_sub_int_scale3 = (int) roundf(ipoint_y_sub_int_scale + l3 * scale);

            for (int k=-12, k_count=0; k<12; ++k, ++k_count) {

                //Get x coords of sample point
                // int sample_x = (int) roundf(ipoint_x + k * scale);
                int sample_x_sub_int_scale = (int) roundf(ipoint_x_sub_int_scale + k * scale);;

                // float rx = 0.0f;
                // float ry = 0.0f;
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale, int_scale, 
                                            &haarResponseX[l_count0*24+k_count], &haarResponseY[l_count0*24+k_count]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale, int_scale, 
                                            &haarResponseX[l_count1*24+k_count], &haarResponseY[l_count1*24+k_count]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale2, sample_x_sub_int_scale, int_scale, 
                                            &haarResponseX[l_count2*24+k_count], &haarResponseY[l_count2*24+k_count]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale3, sample_x_sub_int_scale, int_scale, 
                                            &haarResponseX[l_count3*24+k_count], &haarResponseY[l_count3*24+k_count]);

                // haarResponseX[(l+12)*24+(k+12)] = rx;
                // haarResponseY[(l+12)*24+(k+12)] = ry;
            }

        }

    } else {
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            // int sample_y = (int) roundf(ipoint_y + l * scale);
            int l0 = l;
            int l1 = l+1;
            int l2 = l+2;
            int l3 = l+3;
            int l_count0 = l_count;
            int l_count1 = l_count+1;
            int l_count2 = l_count+2;
            int l_count3 = l_count+3;


            int sample_y_sub_int_scale0 = (int) roundf(ipoint_y_sub_int_scale + l0 * scale);
            int sample_y_sub_int_scale1 = (int) roundf(ipoint_y_sub_int_scale + l1 * scale);
            int sample_y_sub_int_scale2 = (int) roundf(ipoint_y_sub_int_scale + l2 * scale);
            int sample_y_sub_int_scale3 = (int) roundf(ipoint_y_sub_int_scale + l3 * scale);

            for (int k=-12, k_count=0; k<12; ++k, k_count++) {

                //Get x coords of sample point
                // int sample_x = (int) roundf(ipoint_x + k * scale);
                int sample_x_sub_int_scale = (int) roundf(ipoint_x_sub_int_scale + k * scale);;

                // float rx = 0.0f;
                // float ry = 0.0f;
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale, int_scale, 
                                            &haarResponseX[l_count0*24+k_count], &haarResponseY[l_count0*24+k_count]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale, int_scale, 
                                            &haarResponseX[l_count1*24+k_count], &haarResponseY[l_count1*24+k_count]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale2, sample_x_sub_int_scale, int_scale, 
                                            &haarResponseX[l_count2*24+k_count], &haarResponseY[l_count2*24+k_count]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale3, sample_x_sub_int_scale, int_scale, 
                                            &haarResponseX[l_count3*24+k_count], &haarResponseY[l_count3*24+k_count]);

                // haarResponseX[(l+12)*24+(k+12)] = rx;
                // haarResponseY[(l+12)*24+(k+12)] = ry;
            }

        }

    }

    float s0  = roundf( 0.5 * scale);
    float s1  = roundf( 1.5 * scale);
    float s2  = roundf( 2.5 * scale);
    float s3  = roundf( 3.5 * scale);
    float s4  = roundf( 4.5 * scale);
    float s5  = roundf( 5.5 * scale);
    float s6  = roundf( 6.5 * scale);
    float s7  = roundf( 7.5 * scale);
    float s8  = roundf( 8.5 * scale);
    float s9  = roundf( 9.5 * scale);
    float s10 = roundf(10.5 * scale);
    float s11 = roundf(11.5 * scale);

    float e_c0_m4 = s2 + s1; // CAREFUL HERE!
    float e_c0_m3 = s2 + s0; // CAREFUL HERE!
    float e_c0_m2 = s2 - s0;
    float e_c0_m1 = s2 - s1;
    //float e_c0_z0 = s2 - s2;
    float e_c0_p1 = s2 - s3;
    float e_c0_p2 = s2 - s4;
    float e_c0_p3 = s2 - s5;
    float e_c0_p4 = s2 - s6;

    float e_c1_m4 = s7 - s3;
    float e_c1_m3 = s7 - s4;
    float e_c1_m2 = s7 - s5;
    float e_c1_m1 = s7 - s6;
    //float e_c1_z0 = s7 - s7;
    float e_c1_p1 = s7 - s8;
    float e_c1_p2 = s7 - s9;
    float e_c1_p3 = s7 - s10;
    float e_c1_p4 = s7 - s11;

    float gauss_s1_c0_m4 = expf(g1_factor * (e_c0_m4 * e_c0_m4));
    float gauss_s1_c0_m3 = expf(g1_factor * (e_c0_m3 * e_c0_m3));
    float gauss_s1_c0_m2 = expf(g1_factor * (e_c0_m2 * e_c0_m2));
    float gauss_s1_c0_m1 = expf(g1_factor * (e_c0_m1 * e_c0_m1));
    float gauss_s1_c0_z0 = 1.0f; //expf(g1_factor * (e_c0_z0 * e_c0_z0));
    float gauss_s1_c0_p1 = expf(g1_factor * (e_c0_p1 * e_c0_p1));
    float gauss_s1_c0_p2 = expf(g1_factor * (e_c0_p2 * e_c0_p2));
    float gauss_s1_c0_p3 = expf(g1_factor * (e_c0_p3 * e_c0_p3));
    float gauss_s1_c0_p4 = expf(g1_factor * (e_c0_p4 * e_c0_p4));

    float gauss_s1_c1_m4 = expf(g1_factor * (e_c1_m4 * e_c1_m4));
    float gauss_s1_c1_m3 = expf(g1_factor * (e_c1_m3 * e_c1_m3));
    float gauss_s1_c1_m2 = expf(g1_factor * (e_c1_m2 * e_c1_m2));
    float gauss_s1_c1_m1 = expf(g1_factor * (e_c1_m1 * e_c1_m1));
    float gauss_s1_c1_z0 = 1.0f; //expf(g1_factor * (e_c1_z0 * e_c1_z0));
    float gauss_s1_c1_p1 = expf(g1_factor * (e_c1_p1 * e_c1_p1));
    float gauss_s1_c1_p2 = expf(g1_factor * (e_c1_p2 * e_c1_p2));
    float gauss_s1_c1_p3 = expf(g1_factor * (e_c1_p3 * e_c1_p3));
    float gauss_s1_c1_p4 = expf(g1_factor * (e_c1_p4 * e_c1_p4));

    // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        for (int j=-8; j<8; j+=5) {

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            //int xs = (int) roundf(ipoint_x + (i+0.5f) * scale);
            //int ys = (int) roundf(ipoint_y + (j+0.5f) * scale);

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {

                //Get y coords of sample point
                // int sample_y = (int) roundf(ipoint_y + l * scale);
                //float ys_sub_sample_y = (float) ys-sample_y;
                //float ys_sub_sample_y_squared = ys_sub_sample_y*ys_sub_sample_y;

                //Get y coords of sample point
                // int sample_y_sub_int_scale = sample_y-int_scale;

                float gauss_s1_y = -1;
                if (j == -8) {
                    switch (gauss_index_l) {
                        case -4: gauss_s1_y = gauss_s1_c1_p4; break;
                        case -3: gauss_s1_y = gauss_s1_c1_p3; break;
                        case -2: gauss_s1_y = gauss_s1_c1_p2; break;
                        case -1: gauss_s1_y = gauss_s1_c1_p1; break;
                        case 0:  gauss_s1_y = gauss_s1_c1_z0; break;
                        case 1:  gauss_s1_y = gauss_s1_c1_m1; break;
                        case 2:  gauss_s1_y = gauss_s1_c1_m2; break;
                        case 3:  gauss_s1_y = gauss_s1_c1_m3; break;
                        case 4:  gauss_s1_y = gauss_s1_c1_m4; break;
                    };
                } else if (j == -3) {
                    switch (gauss_index_l) {
                        case -4: gauss_s1_y = gauss_s1_c0_p4; break;
                        case -3: gauss_s1_y = gauss_s1_c0_p3; break;
                        case -2: gauss_s1_y = gauss_s1_c0_p2; break;
                        case -1: gauss_s1_y = gauss_s1_c0_p1; break;
                        case 0:  gauss_s1_y = gauss_s1_c0_z0; break;
                        case 1:  gauss_s1_y = gauss_s1_c0_m1; break;
                        case 2:  gauss_s1_y = gauss_s1_c0_m2; break;
                        case 3:  gauss_s1_y = gauss_s1_c0_m3; break;
                        case 4:  gauss_s1_y = gauss_s1_c0_m4; break;
                    };
                } else if (j == 2) {
                    switch (gauss_index_l) {
                        case -4: gauss_s1_y = gauss_s1_c0_m4; break;
                        case -3: gauss_s1_y = gauss_s1_c0_m3; break;
                        case -2: gauss_s1_y = gauss_s1_c0_m2; break;
                        case -1: gauss_s1_y = gauss_s1_c0_m1; break;
                        case 0:  gauss_s1_y = gauss_s1_c0_z0; break;
                        case 1:  gauss_s1_y = gauss_s1_c0_p1; break;
                        case 2:  gauss_s1_y = gauss_s1_c0_p2; break;
                        case 3:  gauss_s1_y = gauss_s1_c0_p3; break;
                        case 4:  gauss_s1_y = gauss_s1_c0_p4; break;
                    };
                } else if (j == 7) {
                    switch (gauss_index_l) {
                        case -4: gauss_s1_y = gauss_s1_c1_m4; break;
                        case -3: gauss_s1_y = gauss_s1_c1_m3; break;
                        case -2: gauss_s1_y = gauss_s1_c1_m2; break;
                        case -1: gauss_s1_y = gauss_s1_c1_m1; break;
                        case 0:  gauss_s1_y = gauss_s1_c1_z0; break;
                        case 1:  gauss_s1_y = gauss_s1_c1_p1; break;
                        case 2:  gauss_s1_y = gauss_s1_c1_p2; break;
                        case 3:  gauss_s1_y = gauss_s1_c1_p3; break;
                        case 4:  gauss_s1_y = gauss_s1_c1_p4; break;
                    };
                }

                int gauss_index_k = -4;
                for (int k = i-4; k < i + 5; ++k, ++gauss_index_k) {
                
                    //Get x coords of sample point
                    // int sample_x = (int) roundf(ipoint_x + k * scale);

                    //Get x coords of sample point
                    // int sample_x_sub_int_scale = sample_x-int_scale;

                    float gauss_s1_x = -1;
                    if (i == -8) {
                        switch (gauss_index_k) {
                            case -4: gauss_s1_x = gauss_s1_c1_p4; break;
                            case -3: gauss_s1_x = gauss_s1_c1_p3; break;
                            case -2: gauss_s1_x = gauss_s1_c1_p2; break;
                            case -1: gauss_s1_x = gauss_s1_c1_p1; break;
                            case 0:  gauss_s1_x = gauss_s1_c1_z0; break;
                            case 1:  gauss_s1_x = gauss_s1_c1_m1; break;
                            case 2:  gauss_s1_x = gauss_s1_c1_m2; break;
                            case 3:  gauss_s1_x = gauss_s1_c1_m3; break;
                            case 4:  gauss_s1_x = gauss_s1_c1_m4; break;
                        };
                    } else if (i == -3) {
                        switch (gauss_index_k) {
                            case -4: gauss_s1_x = gauss_s1_c0_p4; break;
                            case -3: gauss_s1_x = gauss_s1_c0_p3; break;
                            case -2: gauss_s1_x = gauss_s1_c0_p2; break;
                            case -1: gauss_s1_x = gauss_s1_c0_p1; break;
                            case 0:  gauss_s1_x = gauss_s1_c0_z0; break;
                            case 1:  gauss_s1_x = gauss_s1_c0_m1; break;
                            case 2:  gauss_s1_x = gauss_s1_c0_m2; break;
                            case 3:  gauss_s1_x = gauss_s1_c0_m3; break;
                            case 4:  gauss_s1_x = gauss_s1_c0_m4; break;
                        };
                    } else if (i == 2) {
                        switch (gauss_index_k) {
                            case -4: gauss_s1_x = gauss_s1_c0_m4; break;
                            case -3: gauss_s1_x = gauss_s1_c0_m3; break;
                            case -2: gauss_s1_x = gauss_s1_c0_m2; break;
                            case -1: gauss_s1_x = gauss_s1_c0_m1; break;
                            case 0:  gauss_s1_x = gauss_s1_c0_z0; break;
                            case 1:  gauss_s1_x = gauss_s1_c0_p1; break;
                            case 2:  gauss_s1_x = gauss_s1_c0_p2; break;
                            case 3:  gauss_s1_x = gauss_s1_c0_p3; break;
                            case 4:  gauss_s1_x = gauss_s1_c0_p4; break;
                        };
                    } else if (i == 7) {
                        switch (gauss_index_k) {
                            case -4: gauss_s1_x = gauss_s1_c1_m4; break;
                            case -3: gauss_s1_x = gauss_s1_c1_m3; break;
                            case -2: gauss_s1_x = gauss_s1_c1_m2; break;
                            case -1: gauss_s1_x = gauss_s1_c1_m1; break;
                            case 0:  gauss_s1_x = gauss_s1_c1_z0; break;
                            case 1:  gauss_s1_x = gauss_s1_c1_p1; break;
                            case 2:  gauss_s1_x = gauss_s1_c1_p2; break;
                            case 3:  gauss_s1_x = gauss_s1_c1_p3; break;
                            case 4:  gauss_s1_x = gauss_s1_c1_p4; break;
                        };
                    }

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = haarResponseX[(l+12)*24+(k+12)];
                    float ry = haarResponseY[(l+12)*24+(k+12)];
                    // haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &rx, &ry);
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2;
            switch (gauss_s2_index) {
                case 0:  gauss_s2 = 0.026022f; break;
                case 1:  gauss_s2 = 0.040585f; break;
                case 2:  gauss_s2 = 0.040585f; break;
                case 3:  gauss_s2 = 0.026022f; break;
                case 4:  gauss_s2 = 0.040585f; break;
                case 5:  gauss_s2 = 0.063297f; break;
                case 6:  gauss_s2 = 0.063297f; break;
                case 7:  gauss_s2 = 0.040585f; break;
                case 8:  gauss_s2 = 0.040585f; break;
                case 9:  gauss_s2 = 0.063297f; break;
                case 10: gauss_s2 = 0.063297f; break;
                case 11: gauss_s2 = 0.040585f; break;
                case 12: gauss_s2 = 0.026022f; break;
                case 13: gauss_s2 = 0.040585f; break;
                case 14: gauss_s2 = 0.040585f; break;
                case 15: gauss_s2 = 0.026022f; break;
            };
            gauss_s2_index++;

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

void get_msurf_descriptors_gauss_pecompute_haar_unroll(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_gauss_pecompute_haar_unroll(iimage, &interest_points->at(i));
	}
}


void get_msurf_descriptor_rounding_unroll_2_24_True_winner(struct integral_image* iimage, struct interest_point* ipoint) {

    float scale = ipoint->scale;
    int int_scale = (int) roundf(scale);
    float scale_squared = scale*scale;
    float g1_factor = -0.08f / (scale_squared); 

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

    float ipoint_x_sub_int_scale = ipoint_x-int_scale;
    float ipoint_y_sub_int_scale = ipoint_y-int_scale;

    float ipoint_x_sub_int_scale_add_05 = ipoint_x-int_scale + 0.5;
    float ipoint_y_sub_int_scale_add_05 = ipoint_y-int_scale + 0.5;
    
    int width = iimage->width;
    int height = iimage->height;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    // check if we ever hit a boundary
    if (((int) roundf(ipoint_x - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_y - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_x + 11*scale)) + int_scale > width 
        || ((int) roundf(ipoint_y + 11*scale)) + int_scale > height) 
    {
        for (int l=-12, l_count=0; l<12; l+=2, l_count+=2) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=24, k_count+=24) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k4 = k + 4;
                int k5 = k + 5;
                int k6 = k + 6;
                int k7 = k + 7;
                int k8 = k + 8;
                int k9 = k + 9;
                int k10 = k + 10;
                int k11 = k + 11;
                int k12 = k + 12;
                int k13 = k + 13;
                int k14 = k + 14;
                int k15 = k + 15;
                int k16 = k + 16;
                int k17 = k + 17;
                int k18 = k + 18;
                int k19 = k + 19;
                int k20 = k + 20;
                int k21 = k + 21;
                int k22 = k + 22;
                int k23 = k + 23;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                int k_count4 = k_count + 4;
                int k_count5 = k_count + 5;
                int k_count6 = k_count + 6;
                int k_count7 = k_count + 7;
                int k_count8 = k_count + 8;
                int k_count9 = k_count + 9;
                int k_count10 = k_count + 10;
                int k_count11 = k_count + 11;
                int k_count12 = k_count + 12;
                int k_count13 = k_count + 13;
                int k_count14 = k_count + 14;
                int k_count15 = k_count + 15;
                int k_count16 = k_count + 16;
                int k_count17 = k_count + 17;
                int k_count18 = k_count + 18;
                int k_count19 = k_count + 19;
                int k_count20 = k_count + 20;
                int k_count21 = k_count + 21;
                int k_count22 = k_count + 22;
                int k_count23 = k_count + 23;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                float ipoint_x_sub_int_scale_add_k3_mul_scale = ipoint_x_sub_int_scale + k3 * scale;
                float ipoint_x_sub_int_scale_add_k4_mul_scale = ipoint_x_sub_int_scale + k4 * scale;
                float ipoint_x_sub_int_scale_add_k5_mul_scale = ipoint_x_sub_int_scale + k5 * scale;
                float ipoint_x_sub_int_scale_add_k6_mul_scale = ipoint_x_sub_int_scale + k6 * scale;
                float ipoint_x_sub_int_scale_add_k7_mul_scale = ipoint_x_sub_int_scale + k7 * scale;
                float ipoint_x_sub_int_scale_add_k8_mul_scale = ipoint_x_sub_int_scale + k8 * scale;
                float ipoint_x_sub_int_scale_add_k9_mul_scale = ipoint_x_sub_int_scale + k9 * scale;
                float ipoint_x_sub_int_scale_add_k10_mul_scale = ipoint_x_sub_int_scale + k10 * scale;
                float ipoint_x_sub_int_scale_add_k11_mul_scale = ipoint_x_sub_int_scale + k11 * scale;
                float ipoint_x_sub_int_scale_add_k12_mul_scale = ipoint_x_sub_int_scale + k12 * scale;
                float ipoint_x_sub_int_scale_add_k13_mul_scale = ipoint_x_sub_int_scale + k13 * scale;
                float ipoint_x_sub_int_scale_add_k14_mul_scale = ipoint_x_sub_int_scale + k14 * scale;
                float ipoint_x_sub_int_scale_add_k15_mul_scale = ipoint_x_sub_int_scale + k15 * scale;
                float ipoint_x_sub_int_scale_add_k16_mul_scale = ipoint_x_sub_int_scale + k16 * scale;
                float ipoint_x_sub_int_scale_add_k17_mul_scale = ipoint_x_sub_int_scale + k17 * scale;
                float ipoint_x_sub_int_scale_add_k18_mul_scale = ipoint_x_sub_int_scale + k18 * scale;
                float ipoint_x_sub_int_scale_add_k19_mul_scale = ipoint_x_sub_int_scale + k19 * scale;
                float ipoint_x_sub_int_scale_add_k20_mul_scale = ipoint_x_sub_int_scale + k20 * scale;
                float ipoint_x_sub_int_scale_add_k21_mul_scale = ipoint_x_sub_int_scale + k21 * scale;
                float ipoint_x_sub_int_scale_add_k22_mul_scale = ipoint_x_sub_int_scale + k22 * scale;
                float ipoint_x_sub_int_scale_add_k23_mul_scale = ipoint_x_sub_int_scale + k23 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_k3_mul_scale + (ipoint_x_sub_int_scale_add_k3_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale4 = (int) (ipoint_x_sub_int_scale_add_k4_mul_scale + (ipoint_x_sub_int_scale_add_k4_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale5 = (int) (ipoint_x_sub_int_scale_add_k5_mul_scale + (ipoint_x_sub_int_scale_add_k5_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale6 = (int) (ipoint_x_sub_int_scale_add_k6_mul_scale + (ipoint_x_sub_int_scale_add_k6_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale7 = (int) (ipoint_x_sub_int_scale_add_k7_mul_scale + (ipoint_x_sub_int_scale_add_k7_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale8 = (int) (ipoint_x_sub_int_scale_add_k8_mul_scale + (ipoint_x_sub_int_scale_add_k8_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale9 = (int) (ipoint_x_sub_int_scale_add_k9_mul_scale + (ipoint_x_sub_int_scale_add_k9_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale10 = (int) (ipoint_x_sub_int_scale_add_k10_mul_scale + (ipoint_x_sub_int_scale_add_k10_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale11 = (int) (ipoint_x_sub_int_scale_add_k11_mul_scale + (ipoint_x_sub_int_scale_add_k11_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale12 = (int) (ipoint_x_sub_int_scale_add_k12_mul_scale + (ipoint_x_sub_int_scale_add_k12_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale13 = (int) (ipoint_x_sub_int_scale_add_k13_mul_scale + (ipoint_x_sub_int_scale_add_k13_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale14 = (int) (ipoint_x_sub_int_scale_add_k14_mul_scale + (ipoint_x_sub_int_scale_add_k14_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale15 = (int) (ipoint_x_sub_int_scale_add_k15_mul_scale + (ipoint_x_sub_int_scale_add_k15_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale16 = (int) (ipoint_x_sub_int_scale_add_k16_mul_scale + (ipoint_x_sub_int_scale_add_k16_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale17 = (int) (ipoint_x_sub_int_scale_add_k17_mul_scale + (ipoint_x_sub_int_scale_add_k17_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale18 = (int) (ipoint_x_sub_int_scale_add_k18_mul_scale + (ipoint_x_sub_int_scale_add_k18_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale19 = (int) (ipoint_x_sub_int_scale_add_k19_mul_scale + (ipoint_x_sub_int_scale_add_k19_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale20 = (int) (ipoint_x_sub_int_scale_add_k20_mul_scale + (ipoint_x_sub_int_scale_add_k20_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale21 = (int) (ipoint_x_sub_int_scale_add_k21_mul_scale + (ipoint_x_sub_int_scale_add_k21_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale22 = (int) (ipoint_x_sub_int_scale_add_k22_mul_scale + (ipoint_x_sub_int_scale_add_k22_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale23 = (int) (ipoint_x_sub_int_scale_add_k23_mul_scale + (ipoint_x_sub_int_scale_add_k23_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale4, int_scale, &haarResponseX[l_count0*24+k_count4], &haarResponseY[l_count0*24+k_count4]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale4, int_scale, &haarResponseX[l_count1*24+k_count4], &haarResponseY[l_count1*24+k_count4]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale5, int_scale, &haarResponseX[l_count0*24+k_count5], &haarResponseY[l_count0*24+k_count5]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale5, int_scale, &haarResponseX[l_count1*24+k_count5], &haarResponseY[l_count1*24+k_count5]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale6, int_scale, &haarResponseX[l_count0*24+k_count6], &haarResponseY[l_count0*24+k_count6]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale6, int_scale, &haarResponseX[l_count1*24+k_count6], &haarResponseY[l_count1*24+k_count6]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale7, int_scale, &haarResponseX[l_count0*24+k_count7], &haarResponseY[l_count0*24+k_count7]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale7, int_scale, &haarResponseX[l_count1*24+k_count7], &haarResponseY[l_count1*24+k_count7]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale8, int_scale, &haarResponseX[l_count0*24+k_count8], &haarResponseY[l_count0*24+k_count8]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale8, int_scale, &haarResponseX[l_count1*24+k_count8], &haarResponseY[l_count1*24+k_count8]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale9, int_scale, &haarResponseX[l_count0*24+k_count9], &haarResponseY[l_count0*24+k_count9]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale9, int_scale, &haarResponseX[l_count1*24+k_count9], &haarResponseY[l_count1*24+k_count9]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale10, int_scale, &haarResponseX[l_count0*24+k_count10], &haarResponseY[l_count0*24+k_count10]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale10, int_scale, &haarResponseX[l_count1*24+k_count10], &haarResponseY[l_count1*24+k_count10]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale11, int_scale, &haarResponseX[l_count0*24+k_count11], &haarResponseY[l_count0*24+k_count11]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale11, int_scale, &haarResponseX[l_count1*24+k_count11], &haarResponseY[l_count1*24+k_count11]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale12, int_scale, &haarResponseX[l_count0*24+k_count12], &haarResponseY[l_count0*24+k_count12]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale12, int_scale, &haarResponseX[l_count1*24+k_count12], &haarResponseY[l_count1*24+k_count12]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale13, int_scale, &haarResponseX[l_count0*24+k_count13], &haarResponseY[l_count0*24+k_count13]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale13, int_scale, &haarResponseX[l_count1*24+k_count13], &haarResponseY[l_count1*24+k_count13]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale14, int_scale, &haarResponseX[l_count0*24+k_count14], &haarResponseY[l_count0*24+k_count14]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale14, int_scale, &haarResponseX[l_count1*24+k_count14], &haarResponseY[l_count1*24+k_count14]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale15, int_scale, &haarResponseX[l_count0*24+k_count15], &haarResponseY[l_count0*24+k_count15]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale15, int_scale, &haarResponseX[l_count1*24+k_count15], &haarResponseY[l_count1*24+k_count15]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale16, int_scale, &haarResponseX[l_count0*24+k_count16], &haarResponseY[l_count0*24+k_count16]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale16, int_scale, &haarResponseX[l_count1*24+k_count16], &haarResponseY[l_count1*24+k_count16]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale17, int_scale, &haarResponseX[l_count0*24+k_count17], &haarResponseY[l_count0*24+k_count17]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale17, int_scale, &haarResponseX[l_count1*24+k_count17], &haarResponseY[l_count1*24+k_count17]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale18, int_scale, &haarResponseX[l_count0*24+k_count18], &haarResponseY[l_count0*24+k_count18]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale18, int_scale, &haarResponseX[l_count1*24+k_count18], &haarResponseY[l_count1*24+k_count18]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale19, int_scale, &haarResponseX[l_count0*24+k_count19], &haarResponseY[l_count0*24+k_count19]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale19, int_scale, &haarResponseX[l_count1*24+k_count19], &haarResponseY[l_count1*24+k_count19]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale20, int_scale, &haarResponseX[l_count0*24+k_count20], &haarResponseY[l_count0*24+k_count20]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale20, int_scale, &haarResponseX[l_count1*24+k_count20], &haarResponseY[l_count1*24+k_count20]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale21, int_scale, &haarResponseX[l_count0*24+k_count21], &haarResponseY[l_count0*24+k_count21]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale21, int_scale, &haarResponseX[l_count1*24+k_count21], &haarResponseY[l_count1*24+k_count21]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale22, int_scale, &haarResponseX[l_count0*24+k_count22], &haarResponseY[l_count0*24+k_count22]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale22, int_scale, &haarResponseX[l_count1*24+k_count22], &haarResponseY[l_count1*24+k_count22]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale23, int_scale, &haarResponseX[l_count0*24+k_count23], &haarResponseY[l_count0*24+k_count23]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale23, int_scale, &haarResponseX[l_count1*24+k_count23], &haarResponseY[l_count1*24+k_count23]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=2, l_count+=2) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);

            for (int k=-12, k_count=0; k<12; k+=24, k_count+=24) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k4 = k + 4;
                int k5 = k + 5;
                int k6 = k + 6;
                int k7 = k + 7;
                int k8 = k + 8;
                int k9 = k + 9;
                int k10 = k + 10;
                int k11 = k + 11;
                int k12 = k + 12;
                int k13 = k + 13;
                int k14 = k + 14;
                int k15 = k + 15;
                int k16 = k + 16;
                int k17 = k + 17;
                int k18 = k + 18;
                int k19 = k + 19;
                int k20 = k + 20;
                int k21 = k + 21;
                int k22 = k + 22;
                int k23 = k + 23;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                int k_count4 = k_count + 4;
                int k_count5 = k_count + 5;
                int k_count6 = k_count + 6;
                int k_count7 = k_count + 7;
                int k_count8 = k_count + 8;
                int k_count9 = k_count + 9;
                int k_count10 = k_count + 10;
                int k_count11 = k_count + 11;
                int k_count12 = k_count + 12;
                int k_count13 = k_count + 13;
                int k_count14 = k_count + 14;
                int k_count15 = k_count + 15;
                int k_count16 = k_count + 16;
                int k_count17 = k_count + 17;
                int k_count18 = k_count + 18;
                int k_count19 = k_count + 19;
                int k_count20 = k_count + 20;
                int k_count21 = k_count + 21;
                int k_count22 = k_count + 22;
                int k_count23 = k_count + 23;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_05 + k3 * scale);
                int sample_x_sub_int_scale4 = (int) (ipoint_x_sub_int_scale_add_05 + k4 * scale);
                int sample_x_sub_int_scale5 = (int) (ipoint_x_sub_int_scale_add_05 + k5 * scale);
                int sample_x_sub_int_scale6 = (int) (ipoint_x_sub_int_scale_add_05 + k6 * scale);
                int sample_x_sub_int_scale7 = (int) (ipoint_x_sub_int_scale_add_05 + k7 * scale);
                int sample_x_sub_int_scale8 = (int) (ipoint_x_sub_int_scale_add_05 + k8 * scale);
                int sample_x_sub_int_scale9 = (int) (ipoint_x_sub_int_scale_add_05 + k9 * scale);
                int sample_x_sub_int_scale10 = (int) (ipoint_x_sub_int_scale_add_05 + k10 * scale);
                int sample_x_sub_int_scale11 = (int) (ipoint_x_sub_int_scale_add_05 + k11 * scale);
                int sample_x_sub_int_scale12 = (int) (ipoint_x_sub_int_scale_add_05 + k12 * scale);
                int sample_x_sub_int_scale13 = (int) (ipoint_x_sub_int_scale_add_05 + k13 * scale);
                int sample_x_sub_int_scale14 = (int) (ipoint_x_sub_int_scale_add_05 + k14 * scale);
                int sample_x_sub_int_scale15 = (int) (ipoint_x_sub_int_scale_add_05 + k15 * scale);
                int sample_x_sub_int_scale16 = (int) (ipoint_x_sub_int_scale_add_05 + k16 * scale);
                int sample_x_sub_int_scale17 = (int) (ipoint_x_sub_int_scale_add_05 + k17 * scale);
                int sample_x_sub_int_scale18 = (int) (ipoint_x_sub_int_scale_add_05 + k18 * scale);
                int sample_x_sub_int_scale19 = (int) (ipoint_x_sub_int_scale_add_05 + k19 * scale);
                int sample_x_sub_int_scale20 = (int) (ipoint_x_sub_int_scale_add_05 + k20 * scale);
                int sample_x_sub_int_scale21 = (int) (ipoint_x_sub_int_scale_add_05 + k21 * scale);
                int sample_x_sub_int_scale22 = (int) (ipoint_x_sub_int_scale_add_05 + k22 * scale);
                int sample_x_sub_int_scale23 = (int) (ipoint_x_sub_int_scale_add_05 + k23 * scale);

                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale4, int_scale, &haarResponseX[l_count0*24+k_count4], &haarResponseY[l_count0*24+k_count4]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale4, int_scale, &haarResponseX[l_count1*24+k_count4], &haarResponseY[l_count1*24+k_count4]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale5, int_scale, &haarResponseX[l_count0*24+k_count5], &haarResponseY[l_count0*24+k_count5]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale5, int_scale, &haarResponseX[l_count1*24+k_count5], &haarResponseY[l_count1*24+k_count5]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale6, int_scale, &haarResponseX[l_count0*24+k_count6], &haarResponseY[l_count0*24+k_count6]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale6, int_scale, &haarResponseX[l_count1*24+k_count6], &haarResponseY[l_count1*24+k_count6]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale7, int_scale, &haarResponseX[l_count0*24+k_count7], &haarResponseY[l_count0*24+k_count7]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale7, int_scale, &haarResponseX[l_count1*24+k_count7], &haarResponseY[l_count1*24+k_count7]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale8, int_scale, &haarResponseX[l_count0*24+k_count8], &haarResponseY[l_count0*24+k_count8]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale8, int_scale, &haarResponseX[l_count1*24+k_count8], &haarResponseY[l_count1*24+k_count8]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale9, int_scale, &haarResponseX[l_count0*24+k_count9], &haarResponseY[l_count0*24+k_count9]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale9, int_scale, &haarResponseX[l_count1*24+k_count9], &haarResponseY[l_count1*24+k_count9]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale10, int_scale, &haarResponseX[l_count0*24+k_count10], &haarResponseY[l_count0*24+k_count10]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale10, int_scale, &haarResponseX[l_count1*24+k_count10], &haarResponseY[l_count1*24+k_count10]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale11, int_scale, &haarResponseX[l_count0*24+k_count11], &haarResponseY[l_count0*24+k_count11]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale11, int_scale, &haarResponseX[l_count1*24+k_count11], &haarResponseY[l_count1*24+k_count11]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale12, int_scale, &haarResponseX[l_count0*24+k_count12], &haarResponseY[l_count0*24+k_count12]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale12, int_scale, &haarResponseX[l_count1*24+k_count12], &haarResponseY[l_count1*24+k_count12]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale13, int_scale, &haarResponseX[l_count0*24+k_count13], &haarResponseY[l_count0*24+k_count13]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale13, int_scale, &haarResponseX[l_count1*24+k_count13], &haarResponseY[l_count1*24+k_count13]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale14, int_scale, &haarResponseX[l_count0*24+k_count14], &haarResponseY[l_count0*24+k_count14]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale14, int_scale, &haarResponseX[l_count1*24+k_count14], &haarResponseY[l_count1*24+k_count14]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale15, int_scale, &haarResponseX[l_count0*24+k_count15], &haarResponseY[l_count0*24+k_count15]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale15, int_scale, &haarResponseX[l_count1*24+k_count15], &haarResponseY[l_count1*24+k_count15]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale16, int_scale, &haarResponseX[l_count0*24+k_count16], &haarResponseY[l_count0*24+k_count16]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale16, int_scale, &haarResponseX[l_count1*24+k_count16], &haarResponseY[l_count1*24+k_count16]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale17, int_scale, &haarResponseX[l_count0*24+k_count17], &haarResponseY[l_count0*24+k_count17]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale17, int_scale, &haarResponseX[l_count1*24+k_count17], &haarResponseY[l_count1*24+k_count17]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale18, int_scale, &haarResponseX[l_count0*24+k_count18], &haarResponseY[l_count0*24+k_count18]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale18, int_scale, &haarResponseX[l_count1*24+k_count18], &haarResponseY[l_count1*24+k_count18]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale19, int_scale, &haarResponseX[l_count0*24+k_count19], &haarResponseY[l_count0*24+k_count19]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale19, int_scale, &haarResponseX[l_count1*24+k_count19], &haarResponseY[l_count1*24+k_count19]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale20, int_scale, &haarResponseX[l_count0*24+k_count20], &haarResponseY[l_count0*24+k_count20]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale20, int_scale, &haarResponseX[l_count1*24+k_count20], &haarResponseY[l_count1*24+k_count20]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale21, int_scale, &haarResponseX[l_count0*24+k_count21], &haarResponseY[l_count0*24+k_count21]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale21, int_scale, &haarResponseX[l_count1*24+k_count21], &haarResponseY[l_count1*24+k_count21]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale22, int_scale, &haarResponseX[l_count0*24+k_count22], &haarResponseY[l_count0*24+k_count22]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale22, int_scale, &haarResponseX[l_count1*24+k_count22], &haarResponseY[l_count1*24+k_count22]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale23, int_scale, &haarResponseX[l_count0*24+k_count23], &haarResponseY[l_count0*24+k_count23]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale23, int_scale, &haarResponseX[l_count1*24+k_count23], &haarResponseY[l_count1*24+k_count23]);
            }
        }
    }

    float s0  = roundf( 0.5 * scale);
    float s1  = roundf( 1.5 * scale);
    float s2  = roundf( 2.5 * scale);
    float s3  = roundf( 3.5 * scale);
    float s4  = roundf( 4.5 * scale);
    float s5  = roundf( 5.5 * scale);
    float s6  = roundf( 6.5 * scale);
    float s7  = roundf( 7.5 * scale);
    float s8  = roundf( 8.5 * scale);
    float s9  = roundf( 9.5 * scale);
    float s10 = roundf(10.5 * scale);
    float s11 = roundf(11.5 * scale);

    float e_c0_m4 = s2 + s1; // CAREFUL HERE!
    float e_c0_m3 = s2 + s0; // CAREFUL HERE!
    float e_c0_m2 = s2 - s0;
    float e_c0_m1 = s2 - s1;
    //float e_c0_z0 = s2 - s2;
    float e_c0_p1 = s2 - s3;
    float e_c0_p2 = s2 - s4;
    float e_c0_p3 = s2 - s5;
    float e_c0_p4 = s2 - s6;

    float e_c1_m4 = s7 - s3;
    float e_c1_m3 = s7 - s4;
    float e_c1_m2 = s7 - s5;
    float e_c1_m1 = s7 - s6;
    //float e_c1_z0 = s7 - s7;
    float e_c1_p1 = s7 - s8;
    float e_c1_p2 = s7 - s9;
    float e_c1_p3 = s7 - s10;
    float e_c1_p4 = s7 - s11;

    gauss_s1_c0[0] =  expf(g1_factor * (e_c0_m4 * e_c0_m4));
    gauss_s1_c0[1] =  expf(g1_factor * (e_c0_m3 * e_c0_m3));
    gauss_s1_c0[2] =  expf(g1_factor * (e_c0_m2 * e_c0_m2));
    gauss_s1_c0[3] =  expf(g1_factor * (e_c0_m1 * e_c0_m1));
    gauss_s1_c0[4] =  1.0f; //expf(g1_factor * (e_c0_z0 * e_c0_z0));
    gauss_s1_c0[5] =  expf(g1_factor * (e_c0_p1 * e_c0_p1));
    gauss_s1_c0[6] =  expf(g1_factor * (e_c0_p2 * e_c0_p2));
    gauss_s1_c0[7] =  expf(g1_factor * (e_c0_p3 * e_c0_p3));
    gauss_s1_c0[8] =  expf(g1_factor * (e_c0_p4 * e_c0_p4));

    gauss_s1_c1[0] =  expf(g1_factor * (e_c1_m4 * e_c1_m4));
    gauss_s1_c1[1] =  expf(g1_factor * (e_c1_m3 * e_c1_m3));
    gauss_s1_c1[2] =  expf(g1_factor * (e_c1_m2 * e_c1_m2));
    gauss_s1_c1[3] =  expf(g1_factor * (e_c1_m1 * e_c1_m1));
    gauss_s1_c1[4] =  1.0f; //expf(g1_factor * (e_c1_z0 * e_c1_z0));
    gauss_s1_c1[5] =  expf(g1_factor * (e_c1_p1 * e_c1_p1));
    gauss_s1_c1[6] =  expf(g1_factor * (e_c1_p2 * e_c1_p2));
    gauss_s1_c1[7] =  expf(g1_factor * (e_c1_p3 * e_c1_p3));
    gauss_s1_c1[8] =  expf(g1_factor * (e_c1_p4 * e_c1_p4));
    
        // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        float *i_arr;
        if (i == -8 || i == 7) {
            i_arr = gauss_s1_c1;
        } else {
            i_arr = gauss_s1_c0;
        } 

        int gauss_index_k_start = (i<0?8:0);
        int gauss_index_k_inc = (i<0?-1:1);

        for (int j=-8; j<8; j+=5) {

            float *j_arr;
            if (j == -8 || j == 7) {
                j_arr = gauss_s1_c1;
            } else {
                j_arr = gauss_s1_c0;
            } 

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int gauss_index_l = (j<0?8:0);
            int gauss_index_l_inc = (j<0?-1:1);
            for (int l = (j+8)*24; l < (j + 17)*24; l+=24, gauss_index_l+=gauss_index_l_inc) {
                float gauss_s1_y = j_arr[gauss_index_l];

                int gauss_index_k = gauss_index_k_start;
                for (int k = i+8; k < i + 17; ++k, gauss_index_k+=gauss_index_k_inc) {
                    
                    float gauss_s1_x = i_arr[gauss_index_k];

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = haarResponseX[l+k];
                    float ry = haarResponseY[l+k];
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_arr[gauss_s2_index];
            gauss_s2_index++;

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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_rounding_unroll_2_24_True_winner(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_rounding_unroll_2_24_True_winner(iimage, &interest_points->at(i));
    }
}

void get_msurf_descriptor_rounding_unroll_2_24_True_winner_unconditional(struct integral_image* iimage, struct interest_point* ipoint) {

    float scale = ipoint->scale;
    int int_scale = (int) roundf(scale);
    float scale_squared = scale*scale;
    float g1_factor = -0.08f / (scale_squared); 

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

    float ipoint_x_sub_int_scale = ipoint_x-int_scale;
    float ipoint_y_sub_int_scale = ipoint_y-int_scale;

    float ipoint_x_sub_int_scale_add_05 = ipoint_x-int_scale + 0.5;
    float ipoint_y_sub_int_scale_add_05 = ipoint_y-int_scale + 0.5;
    
    int width = iimage->width;
    int height = iimage->height;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    // check if we ever hit a boundary
    // if (((int) roundf(ipoint_x - 12*scale)) - int_scale <= 0 
    //     || ((int) roundf(ipoint_y - 12*scale)) - int_scale <= 0 
    //     || ((int) roundf(ipoint_x + 11*scale)) + int_scale > width 
    //     || ((int) roundf(ipoint_y + 11*scale)) + int_scale > height) 
    // { // still need this for rounding to work -.-
        for (int l=-12, l_count=0; l<12; l+=2, l_count+=2) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            int sample_y_sub_int_scale0 = (int) roundf(ipoint_y_sub_int_scale_add_l0_mul_scale);// + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) roundf(ipoint_y_sub_int_scale_add_l1_mul_scale);// + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=24, k_count+=24) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k4 = k + 4;
                int k5 = k + 5;
                int k6 = k + 6;
                int k7 = k + 7;
                int k8 = k + 8;
                int k9 = k + 9;
                int k10 = k + 10;
                int k11 = k + 11;
                int k12 = k + 12;
                int k13 = k + 13;
                int k14 = k + 14;
                int k15 = k + 15;
                int k16 = k + 16;
                int k17 = k + 17;
                int k18 = k + 18;
                int k19 = k + 19;
                int k20 = k + 20;
                int k21 = k + 21;
                int k22 = k + 22;
                int k23 = k + 23;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                int k_count4 = k_count + 4;
                int k_count5 = k_count + 5;
                int k_count6 = k_count + 6;
                int k_count7 = k_count + 7;
                int k_count8 = k_count + 8;
                int k_count9 = k_count + 9;
                int k_count10 = k_count + 10;
                int k_count11 = k_count + 11;
                int k_count12 = k_count + 12;
                int k_count13 = k_count + 13;
                int k_count14 = k_count + 14;
                int k_count15 = k_count + 15;
                int k_count16 = k_count + 16;
                int k_count17 = k_count + 17;
                int k_count18 = k_count + 18;
                int k_count19 = k_count + 19;
                int k_count20 = k_count + 20;
                int k_count21 = k_count + 21;
                int k_count22 = k_count + 22;
                int k_count23 = k_count + 23;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                float ipoint_x_sub_int_scale_add_k3_mul_scale = ipoint_x_sub_int_scale + k3 * scale;
                float ipoint_x_sub_int_scale_add_k4_mul_scale = ipoint_x_sub_int_scale + k4 * scale;
                float ipoint_x_sub_int_scale_add_k5_mul_scale = ipoint_x_sub_int_scale + k5 * scale;
                float ipoint_x_sub_int_scale_add_k6_mul_scale = ipoint_x_sub_int_scale + k6 * scale;
                float ipoint_x_sub_int_scale_add_k7_mul_scale = ipoint_x_sub_int_scale + k7 * scale;
                float ipoint_x_sub_int_scale_add_k8_mul_scale = ipoint_x_sub_int_scale + k8 * scale;
                float ipoint_x_sub_int_scale_add_k9_mul_scale = ipoint_x_sub_int_scale + k9 * scale;
                float ipoint_x_sub_int_scale_add_k10_mul_scale = ipoint_x_sub_int_scale + k10 * scale;
                float ipoint_x_sub_int_scale_add_k11_mul_scale = ipoint_x_sub_int_scale + k11 * scale;
                float ipoint_x_sub_int_scale_add_k12_mul_scale = ipoint_x_sub_int_scale + k12 * scale;
                float ipoint_x_sub_int_scale_add_k13_mul_scale = ipoint_x_sub_int_scale + k13 * scale;
                float ipoint_x_sub_int_scale_add_k14_mul_scale = ipoint_x_sub_int_scale + k14 * scale;
                float ipoint_x_sub_int_scale_add_k15_mul_scale = ipoint_x_sub_int_scale + k15 * scale;
                float ipoint_x_sub_int_scale_add_k16_mul_scale = ipoint_x_sub_int_scale + k16 * scale;
                float ipoint_x_sub_int_scale_add_k17_mul_scale = ipoint_x_sub_int_scale + k17 * scale;
                float ipoint_x_sub_int_scale_add_k18_mul_scale = ipoint_x_sub_int_scale + k18 * scale;
                float ipoint_x_sub_int_scale_add_k19_mul_scale = ipoint_x_sub_int_scale + k19 * scale;
                float ipoint_x_sub_int_scale_add_k20_mul_scale = ipoint_x_sub_int_scale + k20 * scale;
                float ipoint_x_sub_int_scale_add_k21_mul_scale = ipoint_x_sub_int_scale + k21 * scale;
                float ipoint_x_sub_int_scale_add_k22_mul_scale = ipoint_x_sub_int_scale + k22 * scale;
                float ipoint_x_sub_int_scale_add_k23_mul_scale = ipoint_x_sub_int_scale + k23 * scale;
                int sample_x_sub_int_scale1 = (int) roundf(ipoint_x_sub_int_scale_add_k1_mul_scale);// + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale0 = (int) roundf(ipoint_x_sub_int_scale_add_k0_mul_scale);// + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) roundf(ipoint_x_sub_int_scale_add_k2_mul_scale);// + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale3 = (int) roundf(ipoint_x_sub_int_scale_add_k3_mul_scale);// + (ipoint_x_sub_int_scale_add_k3_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale4 = (int) roundf(ipoint_x_sub_int_scale_add_k4_mul_scale);// + (ipoint_x_sub_int_scale_add_k4_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale5 = (int) roundf(ipoint_x_sub_int_scale_add_k5_mul_scale);// + (ipoint_x_sub_int_scale_add_k5_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale6 = (int) roundf(ipoint_x_sub_int_scale_add_k6_mul_scale);// + (ipoint_x_sub_int_scale_add_k6_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale7 = (int) roundf(ipoint_x_sub_int_scale_add_k7_mul_scale);// + (ipoint_x_sub_int_scale_add_k7_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale8 = (int) roundf(ipoint_x_sub_int_scale_add_k8_mul_scale);// + (ipoint_x_sub_int_scale_add_k8_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale9 = (int) roundf(ipoint_x_sub_int_scale_add_k9_mul_scale);// + (ipoint_x_sub_int_scale_add_k9_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale10 = (int) roundf(ipoint_x_sub_int_scale_add_k10_mul_scale);// + (ipoint_x_sub_int_scale_add_k10_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale11 = (int) roundf(ipoint_x_sub_int_scale_add_k11_mul_scale);// + (ipoint_x_sub_int_scale_add_k11_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale12 = (int) roundf(ipoint_x_sub_int_scale_add_k12_mul_scale);// + (ipoint_x_sub_int_scale_add_k12_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale13 = (int) roundf(ipoint_x_sub_int_scale_add_k13_mul_scale);// + (ipoint_x_sub_int_scale_add_k13_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale14 = (int) roundf(ipoint_x_sub_int_scale_add_k14_mul_scale);// + (ipoint_x_sub_int_scale_add_k14_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale15 = (int) roundf(ipoint_x_sub_int_scale_add_k15_mul_scale);// + (ipoint_x_sub_int_scale_add_k15_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale16 = (int) roundf(ipoint_x_sub_int_scale_add_k16_mul_scale);// + (ipoint_x_sub_int_scale_add_k16_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale17 = (int) roundf(ipoint_x_sub_int_scale_add_k17_mul_scale);// + (ipoint_x_sub_int_scale_add_k17_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale18 = (int) roundf(ipoint_x_sub_int_scale_add_k18_mul_scale);// + (ipoint_x_sub_int_scale_add_k18_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale19 = (int) roundf(ipoint_x_sub_int_scale_add_k19_mul_scale);// + (ipoint_x_sub_int_scale_add_k19_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale20 = (int) roundf(ipoint_x_sub_int_scale_add_k20_mul_scale);// + (ipoint_x_sub_int_scale_add_k20_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale21 = (int) roundf(ipoint_x_sub_int_scale_add_k21_mul_scale);// + (ipoint_x_sub_int_scale_add_k21_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale22 = (int) roundf(ipoint_x_sub_int_scale_add_k22_mul_scale);// + (ipoint_x_sub_int_scale_add_k22_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale23 = (int) roundf(ipoint_x_sub_int_scale_add_k23_mul_scale);// + (ipoint_x_sub_int_scale_add_k23_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale4, int_scale, &haarResponseX[l_count0*24+k_count4], &haarResponseY[l_count0*24+k_count4]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale4, int_scale, &haarResponseX[l_count1*24+k_count4], &haarResponseY[l_count1*24+k_count4]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale5, int_scale, &haarResponseX[l_count0*24+k_count5], &haarResponseY[l_count0*24+k_count5]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale5, int_scale, &haarResponseX[l_count1*24+k_count5], &haarResponseY[l_count1*24+k_count5]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale6, int_scale, &haarResponseX[l_count0*24+k_count6], &haarResponseY[l_count0*24+k_count6]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale6, int_scale, &haarResponseX[l_count1*24+k_count6], &haarResponseY[l_count1*24+k_count6]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale7, int_scale, &haarResponseX[l_count0*24+k_count7], &haarResponseY[l_count0*24+k_count7]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale7, int_scale, &haarResponseX[l_count1*24+k_count7], &haarResponseY[l_count1*24+k_count7]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale8, int_scale, &haarResponseX[l_count0*24+k_count8], &haarResponseY[l_count0*24+k_count8]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale8, int_scale, &haarResponseX[l_count1*24+k_count8], &haarResponseY[l_count1*24+k_count8]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale9, int_scale, &haarResponseX[l_count0*24+k_count9], &haarResponseY[l_count0*24+k_count9]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale9, int_scale, &haarResponseX[l_count1*24+k_count9], &haarResponseY[l_count1*24+k_count9]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale10, int_scale, &haarResponseX[l_count0*24+k_count10], &haarResponseY[l_count0*24+k_count10]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale10, int_scale, &haarResponseX[l_count1*24+k_count10], &haarResponseY[l_count1*24+k_count10]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale11, int_scale, &haarResponseX[l_count0*24+k_count11], &haarResponseY[l_count0*24+k_count11]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale11, int_scale, &haarResponseX[l_count1*24+k_count11], &haarResponseY[l_count1*24+k_count11]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale12, int_scale, &haarResponseX[l_count0*24+k_count12], &haarResponseY[l_count0*24+k_count12]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale12, int_scale, &haarResponseX[l_count1*24+k_count12], &haarResponseY[l_count1*24+k_count12]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale13, int_scale, &haarResponseX[l_count0*24+k_count13], &haarResponseY[l_count0*24+k_count13]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale13, int_scale, &haarResponseX[l_count1*24+k_count13], &haarResponseY[l_count1*24+k_count13]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale14, int_scale, &haarResponseX[l_count0*24+k_count14], &haarResponseY[l_count0*24+k_count14]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale14, int_scale, &haarResponseX[l_count1*24+k_count14], &haarResponseY[l_count1*24+k_count14]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale15, int_scale, &haarResponseX[l_count0*24+k_count15], &haarResponseY[l_count0*24+k_count15]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale15, int_scale, &haarResponseX[l_count1*24+k_count15], &haarResponseY[l_count1*24+k_count15]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale16, int_scale, &haarResponseX[l_count0*24+k_count16], &haarResponseY[l_count0*24+k_count16]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale16, int_scale, &haarResponseX[l_count1*24+k_count16], &haarResponseY[l_count1*24+k_count16]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale17, int_scale, &haarResponseX[l_count0*24+k_count17], &haarResponseY[l_count0*24+k_count17]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale17, int_scale, &haarResponseX[l_count1*24+k_count17], &haarResponseY[l_count1*24+k_count17]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale18, int_scale, &haarResponseX[l_count0*24+k_count18], &haarResponseY[l_count0*24+k_count18]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale18, int_scale, &haarResponseX[l_count1*24+k_count18], &haarResponseY[l_count1*24+k_count18]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale19, int_scale, &haarResponseX[l_count0*24+k_count19], &haarResponseY[l_count0*24+k_count19]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale19, int_scale, &haarResponseX[l_count1*24+k_count19], &haarResponseY[l_count1*24+k_count19]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale20, int_scale, &haarResponseX[l_count0*24+k_count20], &haarResponseY[l_count0*24+k_count20]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale20, int_scale, &haarResponseX[l_count1*24+k_count20], &haarResponseY[l_count1*24+k_count20]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale21, int_scale, &haarResponseX[l_count0*24+k_count21], &haarResponseY[l_count0*24+k_count21]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale22, int_scale, &haarResponseX[l_count1*24+k_count22], &haarResponseY[l_count1*24+k_count22]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale23, int_scale, &haarResponseX[l_count0*24+k_count23], &haarResponseY[l_count0*24+k_count23]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale23, int_scale, &haarResponseX[l_count1*24+k_count23], &haarResponseY[l_count1*24+k_count23]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale21, int_scale, &haarResponseX[l_count1*24+k_count21], &haarResponseY[l_count1*24+k_count21]);
                haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale22, int_scale, &haarResponseX[l_count0*24+k_count22], &haarResponseY[l_count0*24+k_count22]);
            }
        }
    // } else {
    //     for (int l=-12, l_count=0; l<12; l+=2, l_count+=2) {
    //         int l0 = l + 0;
    //         int l1 = l + 1;
    //         int l_count0 = l_count + 0;
    //         int l_count1 = l_count + 1;
    //         int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
    //         int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);

    //         for (int k=-12, k_count=0; k<12; k+=24, k_count+=24) {
    //             int k0 = k + 0;
    //             int k1 = k + 1;
    //             int k2 = k + 2;
    //             int k3 = k + 3;
    //             int k4 = k + 4;
    //             int k5 = k + 5;
    //             int k6 = k + 6;
    //             int k7 = k + 7;
    //             int k8 = k + 8;
    //             int k9 = k + 9;
    //             int k10 = k + 10;
    //             int k11 = k + 11;
    //             int k12 = k + 12;
    //             int k13 = k + 13;
    //             int k14 = k + 14;
    //             int k15 = k + 15;
    //             int k16 = k + 16;
    //             int k17 = k + 17;
    //             int k18 = k + 18;
    //             int k19 = k + 19;
    //             int k20 = k + 20;
    //             int k21 = k + 21;
    //             int k22 = k + 22;
    //             int k23 = k + 23;
    //             int k_count0 = k_count + 0;
    //             int k_count1 = k_count + 1;
    //             int k_count2 = k_count + 2;
    //             int k_count3 = k_count + 3;
    //             int k_count4 = k_count + 4;
    //             int k_count5 = k_count + 5;
    //             int k_count6 = k_count + 6;
    //             int k_count7 = k_count + 7;
    //             int k_count8 = k_count + 8;
    //             int k_count9 = k_count + 9;
    //             int k_count10 = k_count + 10;
    //             int k_count11 = k_count + 11;
    //             int k_count12 = k_count + 12;
    //             int k_count13 = k_count + 13;
    //             int k_count14 = k_count + 14;
    //             int k_count15 = k_count + 15;
    //             int k_count16 = k_count + 16;
    //             int k_count17 = k_count + 17;
    //             int k_count18 = k_count + 18;
    //             int k_count19 = k_count + 19;
    //             int k_count20 = k_count + 20;
    //             int k_count21 = k_count + 21;
    //             int k_count22 = k_count + 22;
    //             int k_count23 = k_count + 23;
    //             int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
    //             int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
    //             int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);
    //             int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_05 + k3 * scale);
    //             int sample_x_sub_int_scale4 = (int) (ipoint_x_sub_int_scale_add_05 + k4 * scale);
    //             int sample_x_sub_int_scale5 = (int) (ipoint_x_sub_int_scale_add_05 + k5 * scale);
    //             int sample_x_sub_int_scale6 = (int) (ipoint_x_sub_int_scale_add_05 + k6 * scale);
    //             int sample_x_sub_int_scale7 = (int) (ipoint_x_sub_int_scale_add_05 + k7 * scale);
    //             int sample_x_sub_int_scale8 = (int) (ipoint_x_sub_int_scale_add_05 + k8 * scale);
    //             int sample_x_sub_int_scale9 = (int) (ipoint_x_sub_int_scale_add_05 + k9 * scale);
    //             int sample_x_sub_int_scale10 = (int) (ipoint_x_sub_int_scale_add_05 + k10 * scale);
    //             int sample_x_sub_int_scale11 = (int) (ipoint_x_sub_int_scale_add_05 + k11 * scale);
    //             int sample_x_sub_int_scale12 = (int) (ipoint_x_sub_int_scale_add_05 + k12 * scale);
    //             int sample_x_sub_int_scale13 = (int) (ipoint_x_sub_int_scale_add_05 + k13 * scale);
    //             int sample_x_sub_int_scale14 = (int) (ipoint_x_sub_int_scale_add_05 + k14 * scale);
    //             int sample_x_sub_int_scale15 = (int) (ipoint_x_sub_int_scale_add_05 + k15 * scale);
    //             int sample_x_sub_int_scale16 = (int) (ipoint_x_sub_int_scale_add_05 + k16 * scale);
    //             int sample_x_sub_int_scale17 = (int) (ipoint_x_sub_int_scale_add_05 + k17 * scale);
    //             int sample_x_sub_int_scale18 = (int) (ipoint_x_sub_int_scale_add_05 + k18 * scale);
    //             int sample_x_sub_int_scale19 = (int) (ipoint_x_sub_int_scale_add_05 + k19 * scale);
    //             int sample_x_sub_int_scale20 = (int) (ipoint_x_sub_int_scale_add_05 + k20 * scale);
    //             int sample_x_sub_int_scale21 = (int) (ipoint_x_sub_int_scale_add_05 + k21 * scale);
    //             int sample_x_sub_int_scale22 = (int) (ipoint_x_sub_int_scale_add_05 + k22 * scale);
    //             int sample_x_sub_int_scale23 = (int) (ipoint_x_sub_int_scale_add_05 + k23 * scale);

    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale4, int_scale, &haarResponseX[l_count0*24+k_count4], &haarResponseY[l_count0*24+k_count4]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale4, int_scale, &haarResponseX[l_count1*24+k_count4], &haarResponseY[l_count1*24+k_count4]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale5, int_scale, &haarResponseX[l_count0*24+k_count5], &haarResponseY[l_count0*24+k_count5]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale5, int_scale, &haarResponseX[l_count1*24+k_count5], &haarResponseY[l_count1*24+k_count5]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale6, int_scale, &haarResponseX[l_count0*24+k_count6], &haarResponseY[l_count0*24+k_count6]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale6, int_scale, &haarResponseX[l_count1*24+k_count6], &haarResponseY[l_count1*24+k_count6]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale7, int_scale, &haarResponseX[l_count0*24+k_count7], &haarResponseY[l_count0*24+k_count7]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale7, int_scale, &haarResponseX[l_count1*24+k_count7], &haarResponseY[l_count1*24+k_count7]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale8, int_scale, &haarResponseX[l_count0*24+k_count8], &haarResponseY[l_count0*24+k_count8]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale8, int_scale, &haarResponseX[l_count1*24+k_count8], &haarResponseY[l_count1*24+k_count8]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale9, int_scale, &haarResponseX[l_count0*24+k_count9], &haarResponseY[l_count0*24+k_count9]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale9, int_scale, &haarResponseX[l_count1*24+k_count9], &haarResponseY[l_count1*24+k_count9]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale10, int_scale, &haarResponseX[l_count0*24+k_count10], &haarResponseY[l_count0*24+k_count10]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale10, int_scale, &haarResponseX[l_count1*24+k_count10], &haarResponseY[l_count1*24+k_count10]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale11, int_scale, &haarResponseX[l_count0*24+k_count11], &haarResponseY[l_count0*24+k_count11]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale11, int_scale, &haarResponseX[l_count1*24+k_count11], &haarResponseY[l_count1*24+k_count11]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale12, int_scale, &haarResponseX[l_count0*24+k_count12], &haarResponseY[l_count0*24+k_count12]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale12, int_scale, &haarResponseX[l_count1*24+k_count12], &haarResponseY[l_count1*24+k_count12]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale13, int_scale, &haarResponseX[l_count0*24+k_count13], &haarResponseY[l_count0*24+k_count13]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale13, int_scale, &haarResponseX[l_count1*24+k_count13], &haarResponseY[l_count1*24+k_count13]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale14, int_scale, &haarResponseX[l_count0*24+k_count14], &haarResponseY[l_count0*24+k_count14]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale14, int_scale, &haarResponseX[l_count1*24+k_count14], &haarResponseY[l_count1*24+k_count14]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale15, int_scale, &haarResponseX[l_count0*24+k_count15], &haarResponseY[l_count0*24+k_count15]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale15, int_scale, &haarResponseX[l_count1*24+k_count15], &haarResponseY[l_count1*24+k_count15]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale16, int_scale, &haarResponseX[l_count0*24+k_count16], &haarResponseY[l_count0*24+k_count16]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale16, int_scale, &haarResponseX[l_count1*24+k_count16], &haarResponseY[l_count1*24+k_count16]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale17, int_scale, &haarResponseX[l_count0*24+k_count17], &haarResponseY[l_count0*24+k_count17]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale17, int_scale, &haarResponseX[l_count1*24+k_count17], &haarResponseY[l_count1*24+k_count17]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale18, int_scale, &haarResponseX[l_count0*24+k_count18], &haarResponseY[l_count0*24+k_count18]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale18, int_scale, &haarResponseX[l_count1*24+k_count18], &haarResponseY[l_count1*24+k_count18]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale19, int_scale, &haarResponseX[l_count0*24+k_count19], &haarResponseY[l_count0*24+k_count19]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale19, int_scale, &haarResponseX[l_count1*24+k_count19], &haarResponseY[l_count1*24+k_count19]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale20, int_scale, &haarResponseX[l_count0*24+k_count20], &haarResponseY[l_count0*24+k_count20]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale20, int_scale, &haarResponseX[l_count1*24+k_count20], &haarResponseY[l_count1*24+k_count20]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale21, int_scale, &haarResponseX[l_count0*24+k_count21], &haarResponseY[l_count0*24+k_count21]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale21, int_scale, &haarResponseX[l_count1*24+k_count21], &haarResponseY[l_count1*24+k_count21]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale22, int_scale, &haarResponseX[l_count0*24+k_count22], &haarResponseY[l_count0*24+k_count22]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale22, int_scale, &haarResponseX[l_count1*24+k_count22], &haarResponseY[l_count1*24+k_count22]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale23, int_scale, &haarResponseX[l_count0*24+k_count23], &haarResponseY[l_count0*24+k_count23]);
    //             haarXY_unconditional(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale23, int_scale, &haarResponseX[l_count1*24+k_count23], &haarResponseY[l_count1*24+k_count23]);
    //         }
    //     }
    // }

    float s0  = roundf( 0.5 * scale);
    float s1  = roundf( 1.5 * scale);
    float s2  = roundf( 2.5 * scale);
    float s3  = roundf( 3.5 * scale);
    float s4  = roundf( 4.5 * scale);
    float s5  = roundf( 5.5 * scale);
    float s6  = roundf( 6.5 * scale);
    float s7  = roundf( 7.5 * scale);
    float s8  = roundf( 8.5 * scale);
    float s9  = roundf( 9.5 * scale);
    float s10 = roundf(10.5 * scale);
    float s11 = roundf(11.5 * scale);

    float e_c0_m4 = s2 + s1; // CAREFUL HERE!
    float e_c0_m3 = s2 + s0; // CAREFUL HERE!
    float e_c0_m2 = s2 - s0;
    float e_c0_m1 = s2 - s1;
    //float e_c0_z0 = s2 - s2;
    float e_c0_p1 = s2 - s3;
    float e_c0_p2 = s2 - s4;
    float e_c0_p3 = s2 - s5;
    float e_c0_p4 = s2 - s6;

    float e_c1_m4 = s7 - s3;
    float e_c1_m3 = s7 - s4;
    float e_c1_m2 = s7 - s5;
    float e_c1_m1 = s7 - s6;
    //float e_c1_z0 = s7 - s7;
    float e_c1_p1 = s7 - s8;
    float e_c1_p2 = s7 - s9;
    float e_c1_p3 = s7 - s10;
    float e_c1_p4 = s7 - s11;

    gauss_s1_c0[0] =  expf(g1_factor * (e_c0_m4 * e_c0_m4));
    gauss_s1_c0[1] =  expf(g1_factor * (e_c0_m3 * e_c0_m3));
    gauss_s1_c0[2] =  expf(g1_factor * (e_c0_m2 * e_c0_m2));
    gauss_s1_c0[3] =  expf(g1_factor * (e_c0_m1 * e_c0_m1));
    gauss_s1_c0[4] =  1.0f; //expf(g1_factor * (e_c0_z0 * e_c0_z0));
    gauss_s1_c0[5] =  expf(g1_factor * (e_c0_p1 * e_c0_p1));
    gauss_s1_c0[6] =  expf(g1_factor * (e_c0_p2 * e_c0_p2));
    gauss_s1_c0[7] =  expf(g1_factor * (e_c0_p3 * e_c0_p3));
    gauss_s1_c0[8] =  expf(g1_factor * (e_c0_p4 * e_c0_p4));

    gauss_s1_c1[0] =  expf(g1_factor * (e_c1_m4 * e_c1_m4));
    gauss_s1_c1[1] =  expf(g1_factor * (e_c1_m3 * e_c1_m3));
    gauss_s1_c1[2] =  expf(g1_factor * (e_c1_m2 * e_c1_m2));
    gauss_s1_c1[3] =  expf(g1_factor * (e_c1_m1 * e_c1_m1));
    gauss_s1_c1[4] =  1.0f; //expf(g1_factor * (e_c1_z0 * e_c1_z0));
    gauss_s1_c1[5] =  expf(g1_factor * (e_c1_p1 * e_c1_p1));
    gauss_s1_c1[6] =  expf(g1_factor * (e_c1_p2 * e_c1_p2));
    gauss_s1_c1[7] =  expf(g1_factor * (e_c1_p3 * e_c1_p3));
    gauss_s1_c1[8] =  expf(g1_factor * (e_c1_p4 * e_c1_p4));
    
        // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        float *i_arr;
        if (i == -8 || i == 7) {
            i_arr = gauss_s1_c1;
        } else {
            i_arr = gauss_s1_c0;
        } 

        int gauss_index_k_start = (i<0?8:0);
        int gauss_index_k_inc = (i<0?-1:1);

        for (int j=-8; j<8; j+=5) {

            float *j_arr;
            if (j == -8 || j == 7) {
                j_arr = gauss_s1_c1;
            } else {
                j_arr = gauss_s1_c0;
            } 

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int gauss_index_l = (j<0?8:0);
            int gauss_index_l_inc = (j<0?-1:1);
            for (int l = (j+8)*24; l < (j + 17)*24; l+=24, gauss_index_l+=gauss_index_l_inc) {
                float gauss_s1_y = j_arr[gauss_index_l];

                int gauss_index_k = gauss_index_k_start;
                for (int k = i+8; k < i + 17; ++k, gauss_index_k+=gauss_index_k_inc) {
                    
                    float gauss_s1_x = i_arr[gauss_index_k];

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = haarResponseX[l+k];
                    float ry = haarResponseY[l+k];
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_arr[gauss_s2_index];
            gauss_s2_index++;

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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_rounding_unroll_2_24_True_winner_unconditional(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_rounding_unroll_2_24_True_winner_unconditional(iimage, &interest_points->at(i));
    }
}

const __m256 kk0_inner_arrays_simd = _mm256_setr_ps(-12.0f, -11.0f, -10.0f,  -9.0f,  -8.0f,  -7.0f,  -6.0f,  -5.0f);
const __m256 kk1_inner_arrays_simd = _mm256_setr_ps( -4.0f,  -3.0f,  -2.0f,  -1.0f,   0.0f,   1.0f,   2.0f,   3.0f);
const __m256 kk2_inner_arrays_simd = _mm256_setr_ps(  4.0f,   5.0f,   6.0f,   7.0f,   8.0f,   9.0f,  10.0f,  11.0f);

       
void get_msurf_descriptor_simd(struct integral_image* iimage, struct interest_point* ipoint) {

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale = (int) roundf(scale);

    __m256i int_scale_vec = _mm256_set1_epi32(int_scale);
    __m256 scale_vec = _mm256_set1_ps(scale);

    // int int_scale_mul_2 =  2 * int_scale;
    float scale_squared = scale*scale;

    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    //float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

    float ipoint_x_sub_int_scale = ipoint_x-int_scale;
    float ipoint_y_sub_int_scale = ipoint_y-int_scale;

    float ipoint_x_sub_int_scale_add_05 = ipoint_x-int_scale + 0.5;
    float ipoint_y_sub_int_scale_add_05 = ipoint_y-int_scale + 0.5;
    
    int width = iimage->width;
    int height = iimage->height;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    //r(2.5*s) + r(1.5*s) == -(r(2.5*s) - r(6.5*s))

    // check if we ever hit a boundary
    if (((int) roundf(ipoint_x - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_y - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_x + 11*scale)) + int_scale > width 
        || ((int) roundf(ipoint_y + 11*scale)) + int_scale > height) 
    { // some outside
        for (int l=-12, l_count=0; l<12; ++l, ++l_count) {
            float ipoint_y_sub_int_scale_add_l_mul_scale = ipoint_y_sub_int_scale + l * scale;
            int sample_y_sub_int_scale = (int) (ipoint_y_sub_int_scale_add_l_mul_scale + (ipoint_y_sub_int_scale_add_l_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; ++k, k_count++) {

                //Get x coords of sample point
                float ipoint_x_sub_int_scale_add_k_mul_scale = ipoint_x_sub_int_scale + k * scale;
                int sample_x_sub_int_scale = (int) (ipoint_x_sub_int_scale_add_k_mul_scale + (ipoint_x_sub_int_scale_add_k_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &haarResponseX[l_count*24+k_count], &haarResponseY[l_count*24+k_count]);
            }

        }

    } else {

        for (int l=-12, l_count=0; l<12; ++l, ++l_count) {
            int sample_y_sub_int_scale = (int)(ipoint_y_sub_int_scale_add_05 + l * scale);

            __m256i sample_y_sub_int_scale_vec = _mm256_set1_epi32(sample_y_sub_int_scale);

            __m256 kks0 = _mm256_mul_ps(kk0_inner_arrays_simd, scale_vec);
            __m256 kks1 = _mm256_mul_ps(kk1_inner_arrays_simd, scale_vec);
            __m256 kks2 = _mm256_mul_ps(kk2_inner_arrays_simd, scale_vec);

            __m256 ipoint_x_sub_int_scale_add_05_vec = _mm256_set1_ps(ipoint_x_sub_int_scale_add_05);

            // USE CVTTPS_EPI32 FOR TRUNCATION!!!
            __m256i sample_col0 = _mm256_cvttps_epi32(_mm256_add_ps(ipoint_x_sub_int_scale_add_05_vec, kks0));
            __m256i sample_col1 = _mm256_cvttps_epi32(_mm256_add_ps(ipoint_x_sub_int_scale_add_05_vec, kks1));
            __m256i sample_col2 = _mm256_cvttps_epi32(_mm256_add_ps(ipoint_x_sub_int_scale_add_05_vec, kks2));
            
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec, sample_col0, int_scale_vec, 
                                            &haarResponseX[l_count*24+0], &haarResponseY[l_count*24+0]);
            
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec, sample_col1, int_scale_vec, 
                                            &haarResponseX[l_count*24+8], &haarResponseY[l_count*24+8]);
            
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec, sample_col2, int_scale_vec, 
                                            &haarResponseX[l_count*24+16], &haarResponseY[l_count*24+16]);

        }
        
    }

    float s0  = roundf( 0.5 * scale);
    float s1  = roundf( 1.5 * scale);
    float s2  = roundf( 2.5 * scale);
    float s3  = roundf( 3.5 * scale);
    float s4  = roundf( 4.5 * scale);
    float s5  = roundf( 5.5 * scale);
    float s6  = roundf( 6.5 * scale);
    float s7  = roundf( 7.5 * scale);
    float s8  = roundf( 8.5 * scale);
    float s9  = roundf( 9.5 * scale);
    float s10 = roundf(10.5 * scale);
    float s11 = roundf(11.5 * scale);

    float e_c0_m4 = s2 + s1; // CAREFUL HERE!
    float e_c0_m3 = s2 + s0; // CAREFUL HERE!
    float e_c0_m2 = s2 - s0;
    float e_c0_m1 = s2 - s1;
    //float e_c0_z0 = s2 - s2;
    float e_c0_p1 = s2 - s3;
    float e_c0_p2 = s2 - s4;
    float e_c0_p3 = s2 - s5;
    float e_c0_p4 = s2 - s6;

    float e_c1_m4 = s7 - s3;
    float e_c1_m3 = s7 - s4;
    float e_c1_m2 = s7 - s5;
    float e_c1_m1 = s7 - s6;
    //float e_c1_z0 = s7 - s7;
    float e_c1_p1 = s7 - s8;
    float e_c1_p2 = s7 - s9;
    float e_c1_p3 = s7 - s10;
    float e_c1_p4 = s7 - s11;

    gauss_s1_c0[0] =  expf(g1_factor * (e_c0_m4 * e_c0_m4));
    gauss_s1_c0[1] =  expf(g1_factor * (e_c0_m3 * e_c0_m3));
    gauss_s1_c0[2] =  expf(g1_factor * (e_c0_m2 * e_c0_m2));
    gauss_s1_c0[3] =  expf(g1_factor * (e_c0_m1 * e_c0_m1));
    gauss_s1_c0[4] =  1.0f; //expf(g1_factor * (e_c0_z0 * e_c0_z0));
    gauss_s1_c0[5] =  expf(g1_factor * (e_c0_p1 * e_c0_p1));
    gauss_s1_c0[6] =  expf(g1_factor * (e_c0_p2 * e_c0_p2));
    gauss_s1_c0[7] =  expf(g1_factor * (e_c0_p3 * e_c0_p3));
    gauss_s1_c0[8] =  expf(g1_factor * (e_c0_p4 * e_c0_p4));

    gauss_s1_c1[0] =  expf(g1_factor * (e_c1_m4 * e_c1_m4));
    gauss_s1_c1[1] =  expf(g1_factor * (e_c1_m3 * e_c1_m3));
    gauss_s1_c1[2] =  expf(g1_factor * (e_c1_m2 * e_c1_m2));
    gauss_s1_c1[3] =  expf(g1_factor * (e_c1_m1 * e_c1_m1));
    gauss_s1_c1[4] =  1.0f; //expf(g1_factor * (e_c1_z0 * e_c1_z0));
    gauss_s1_c1[5] =  expf(g1_factor * (e_c1_p1 * e_c1_p1));
    gauss_s1_c1[6] =  expf(g1_factor * (e_c1_p2 * e_c1_p2));
    gauss_s1_c1[7] =  expf(g1_factor * (e_c1_p3 * e_c1_p3));
    gauss_s1_c1[8] =  expf(g1_factor * (e_c1_p4 * e_c1_p4));
    
    // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        float *i_arr;
        if (i == -8 || i == 7) {
            i_arr = gauss_s1_c1;
        } else {
            i_arr = gauss_s1_c0;
        } 

        int gauss_index_k_start = (i<0?8:0);
        int gauss_index_k_inc = (i<0?-1:1);

        for (int j=-8; j<8; j+=5) {

            float *j_arr;
            if (j == -8 || j == 7) {
                j_arr = gauss_s1_c1;
            } else {
                j_arr = gauss_s1_c0;
            } 

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int gauss_index_l = (j<0?8:0);
            int gauss_index_l_inc = (j<0?-1:1);
            for (int l = (j+8)*24; l < (j + 17)*24; l+=24, gauss_index_l+=gauss_index_l_inc) {
                float gauss_s1_y = j_arr[gauss_index_l];

                int gauss_index_k = gauss_index_k_start;
                for (int k = i+8; k < i + 17; ++k, gauss_index_k+=gauss_index_k_inc) {
                    
                    float gauss_s1_x = i_arr[gauss_index_k];

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = haarResponseX[l+k];
                    float ry = haarResponseY[l+k];
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_arr[gauss_s2_index];
            gauss_s2_index++;

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


void get_msurf_descriptors_simd(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_simd(iimage, &interest_points->at(i));
	}
}


void get_msurf_descriptor_simd_2_24(struct integral_image* iimage, struct interest_point* ipoint) {

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale = (int) roundf(scale);

    __m256i int_scale_vec = _mm256_set1_epi32(int_scale);
    __m256 scale_vec = _mm256_set1_ps(scale);

    // int int_scale_mul_2 =  2 * int_scale;
    float scale_squared = scale*scale;

    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    //float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

    float ipoint_x_sub_int_scale = ipoint_x-int_scale;
    float ipoint_y_sub_int_scale = ipoint_y-int_scale;

    float ipoint_x_sub_int_scale_add_05 = ipoint_x-int_scale + 0.5;
    float ipoint_y_sub_int_scale_add_05 = ipoint_y-int_scale + 0.5;
    
    int width = iimage->width;
    int height = iimage->height;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    //r(2.5*s) + r(1.5*s) == -(r(2.5*s) - r(6.5*s))

    // check if we ever hit a boundary
    if (((int) roundf(ipoint_x - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_y - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_x + 11*scale)) + int_scale > width 
        || ((int) roundf(ipoint_y + 11*scale)) + int_scale > height) 
    { // some outside
        for (int l=-12, l_count=0; l<12; l+=2, l_count+=2) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=24, k_count+=24) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k4 = k + 4;
                int k5 = k + 5;
                int k6 = k + 6;
                int k7 = k + 7;
                int k8 = k + 8;
                int k9 = k + 9;
                int k10 = k + 10;
                int k11 = k + 11;
                int k12 = k + 12;
                int k13 = k + 13;
                int k14 = k + 14;
                int k15 = k + 15;
                int k16 = k + 16;
                int k17 = k + 17;
                int k18 = k + 18;
                int k19 = k + 19;
                int k20 = k + 20;
                int k21 = k + 21;
                int k22 = k + 22;
                int k23 = k + 23;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                int k_count4 = k_count + 4;
                int k_count5 = k_count + 5;
                int k_count6 = k_count + 6;
                int k_count7 = k_count + 7;
                int k_count8 = k_count + 8;
                int k_count9 = k_count + 9;
                int k_count10 = k_count + 10;
                int k_count11 = k_count + 11;
                int k_count12 = k_count + 12;
                int k_count13 = k_count + 13;
                int k_count14 = k_count + 14;
                int k_count15 = k_count + 15;
                int k_count16 = k_count + 16;
                int k_count17 = k_count + 17;
                int k_count18 = k_count + 18;
                int k_count19 = k_count + 19;
                int k_count20 = k_count + 20;
                int k_count21 = k_count + 21;
                int k_count22 = k_count + 22;
                int k_count23 = k_count + 23;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                float ipoint_x_sub_int_scale_add_k3_mul_scale = ipoint_x_sub_int_scale + k3 * scale;
                float ipoint_x_sub_int_scale_add_k4_mul_scale = ipoint_x_sub_int_scale + k4 * scale;
                float ipoint_x_sub_int_scale_add_k5_mul_scale = ipoint_x_sub_int_scale + k5 * scale;
                float ipoint_x_sub_int_scale_add_k6_mul_scale = ipoint_x_sub_int_scale + k6 * scale;
                float ipoint_x_sub_int_scale_add_k7_mul_scale = ipoint_x_sub_int_scale + k7 * scale;
                float ipoint_x_sub_int_scale_add_k8_mul_scale = ipoint_x_sub_int_scale + k8 * scale;
                float ipoint_x_sub_int_scale_add_k9_mul_scale = ipoint_x_sub_int_scale + k9 * scale;
                float ipoint_x_sub_int_scale_add_k10_mul_scale = ipoint_x_sub_int_scale + k10 * scale;
                float ipoint_x_sub_int_scale_add_k11_mul_scale = ipoint_x_sub_int_scale + k11 * scale;
                float ipoint_x_sub_int_scale_add_k12_mul_scale = ipoint_x_sub_int_scale + k12 * scale;
                float ipoint_x_sub_int_scale_add_k13_mul_scale = ipoint_x_sub_int_scale + k13 * scale;
                float ipoint_x_sub_int_scale_add_k14_mul_scale = ipoint_x_sub_int_scale + k14 * scale;
                float ipoint_x_sub_int_scale_add_k15_mul_scale = ipoint_x_sub_int_scale + k15 * scale;
                float ipoint_x_sub_int_scale_add_k16_mul_scale = ipoint_x_sub_int_scale + k16 * scale;
                float ipoint_x_sub_int_scale_add_k17_mul_scale = ipoint_x_sub_int_scale + k17 * scale;
                float ipoint_x_sub_int_scale_add_k18_mul_scale = ipoint_x_sub_int_scale + k18 * scale;
                float ipoint_x_sub_int_scale_add_k19_mul_scale = ipoint_x_sub_int_scale + k19 * scale;
                float ipoint_x_sub_int_scale_add_k20_mul_scale = ipoint_x_sub_int_scale + k20 * scale;
                float ipoint_x_sub_int_scale_add_k21_mul_scale = ipoint_x_sub_int_scale + k21 * scale;
                float ipoint_x_sub_int_scale_add_k22_mul_scale = ipoint_x_sub_int_scale + k22 * scale;
                float ipoint_x_sub_int_scale_add_k23_mul_scale = ipoint_x_sub_int_scale + k23 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_k3_mul_scale + (ipoint_x_sub_int_scale_add_k3_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale4 = (int) (ipoint_x_sub_int_scale_add_k4_mul_scale + (ipoint_x_sub_int_scale_add_k4_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale5 = (int) (ipoint_x_sub_int_scale_add_k5_mul_scale + (ipoint_x_sub_int_scale_add_k5_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale6 = (int) (ipoint_x_sub_int_scale_add_k6_mul_scale + (ipoint_x_sub_int_scale_add_k6_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale7 = (int) (ipoint_x_sub_int_scale_add_k7_mul_scale + (ipoint_x_sub_int_scale_add_k7_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale8 = (int) (ipoint_x_sub_int_scale_add_k8_mul_scale + (ipoint_x_sub_int_scale_add_k8_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale9 = (int) (ipoint_x_sub_int_scale_add_k9_mul_scale + (ipoint_x_sub_int_scale_add_k9_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale10 = (int) (ipoint_x_sub_int_scale_add_k10_mul_scale + (ipoint_x_sub_int_scale_add_k10_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale11 = (int) (ipoint_x_sub_int_scale_add_k11_mul_scale + (ipoint_x_sub_int_scale_add_k11_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale12 = (int) (ipoint_x_sub_int_scale_add_k12_mul_scale + (ipoint_x_sub_int_scale_add_k12_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale13 = (int) (ipoint_x_sub_int_scale_add_k13_mul_scale + (ipoint_x_sub_int_scale_add_k13_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale14 = (int) (ipoint_x_sub_int_scale_add_k14_mul_scale + (ipoint_x_sub_int_scale_add_k14_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale15 = (int) (ipoint_x_sub_int_scale_add_k15_mul_scale + (ipoint_x_sub_int_scale_add_k15_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale16 = (int) (ipoint_x_sub_int_scale_add_k16_mul_scale + (ipoint_x_sub_int_scale_add_k16_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale17 = (int) (ipoint_x_sub_int_scale_add_k17_mul_scale + (ipoint_x_sub_int_scale_add_k17_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale18 = (int) (ipoint_x_sub_int_scale_add_k18_mul_scale + (ipoint_x_sub_int_scale_add_k18_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale19 = (int) (ipoint_x_sub_int_scale_add_k19_mul_scale + (ipoint_x_sub_int_scale_add_k19_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale20 = (int) (ipoint_x_sub_int_scale_add_k20_mul_scale + (ipoint_x_sub_int_scale_add_k20_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale21 = (int) (ipoint_x_sub_int_scale_add_k21_mul_scale + (ipoint_x_sub_int_scale_add_k21_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale22 = (int) (ipoint_x_sub_int_scale_add_k22_mul_scale + (ipoint_x_sub_int_scale_add_k22_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale23 = (int) (ipoint_x_sub_int_scale_add_k23_mul_scale + (ipoint_x_sub_int_scale_add_k23_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale4, int_scale, &haarResponseX[l_count0*24+k_count4], &haarResponseY[l_count0*24+k_count4]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale4, int_scale, &haarResponseX[l_count1*24+k_count4], &haarResponseY[l_count1*24+k_count4]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale5, int_scale, &haarResponseX[l_count0*24+k_count5], &haarResponseY[l_count0*24+k_count5]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale5, int_scale, &haarResponseX[l_count1*24+k_count5], &haarResponseY[l_count1*24+k_count5]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale6, int_scale, &haarResponseX[l_count0*24+k_count6], &haarResponseY[l_count0*24+k_count6]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale6, int_scale, &haarResponseX[l_count1*24+k_count6], &haarResponseY[l_count1*24+k_count6]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale7, int_scale, &haarResponseX[l_count0*24+k_count7], &haarResponseY[l_count0*24+k_count7]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale7, int_scale, &haarResponseX[l_count1*24+k_count7], &haarResponseY[l_count1*24+k_count7]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale8, int_scale, &haarResponseX[l_count0*24+k_count8], &haarResponseY[l_count0*24+k_count8]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale8, int_scale, &haarResponseX[l_count1*24+k_count8], &haarResponseY[l_count1*24+k_count8]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale9, int_scale, &haarResponseX[l_count0*24+k_count9], &haarResponseY[l_count0*24+k_count9]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale9, int_scale, &haarResponseX[l_count1*24+k_count9], &haarResponseY[l_count1*24+k_count9]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale10, int_scale, &haarResponseX[l_count0*24+k_count10], &haarResponseY[l_count0*24+k_count10]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale10, int_scale, &haarResponseX[l_count1*24+k_count10], &haarResponseY[l_count1*24+k_count10]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale11, int_scale, &haarResponseX[l_count0*24+k_count11], &haarResponseY[l_count0*24+k_count11]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale11, int_scale, &haarResponseX[l_count1*24+k_count11], &haarResponseY[l_count1*24+k_count11]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale12, int_scale, &haarResponseX[l_count0*24+k_count12], &haarResponseY[l_count0*24+k_count12]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale12, int_scale, &haarResponseX[l_count1*24+k_count12], &haarResponseY[l_count1*24+k_count12]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale13, int_scale, &haarResponseX[l_count0*24+k_count13], &haarResponseY[l_count0*24+k_count13]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale13, int_scale, &haarResponseX[l_count1*24+k_count13], &haarResponseY[l_count1*24+k_count13]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale14, int_scale, &haarResponseX[l_count0*24+k_count14], &haarResponseY[l_count0*24+k_count14]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale14, int_scale, &haarResponseX[l_count1*24+k_count14], &haarResponseY[l_count1*24+k_count14]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale15, int_scale, &haarResponseX[l_count0*24+k_count15], &haarResponseY[l_count0*24+k_count15]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale15, int_scale, &haarResponseX[l_count1*24+k_count15], &haarResponseY[l_count1*24+k_count15]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale16, int_scale, &haarResponseX[l_count0*24+k_count16], &haarResponseY[l_count0*24+k_count16]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale16, int_scale, &haarResponseX[l_count1*24+k_count16], &haarResponseY[l_count1*24+k_count16]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale17, int_scale, &haarResponseX[l_count0*24+k_count17], &haarResponseY[l_count0*24+k_count17]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale17, int_scale, &haarResponseX[l_count1*24+k_count17], &haarResponseY[l_count1*24+k_count17]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale18, int_scale, &haarResponseX[l_count0*24+k_count18], &haarResponseY[l_count0*24+k_count18]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale18, int_scale, &haarResponseX[l_count1*24+k_count18], &haarResponseY[l_count1*24+k_count18]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale19, int_scale, &haarResponseX[l_count0*24+k_count19], &haarResponseY[l_count0*24+k_count19]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale19, int_scale, &haarResponseX[l_count1*24+k_count19], &haarResponseY[l_count1*24+k_count19]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale20, int_scale, &haarResponseX[l_count0*24+k_count20], &haarResponseY[l_count0*24+k_count20]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale20, int_scale, &haarResponseX[l_count1*24+k_count20], &haarResponseY[l_count1*24+k_count20]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale21, int_scale, &haarResponseX[l_count0*24+k_count21], &haarResponseY[l_count0*24+k_count21]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale21, int_scale, &haarResponseX[l_count1*24+k_count21], &haarResponseY[l_count1*24+k_count21]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale22, int_scale, &haarResponseX[l_count0*24+k_count22], &haarResponseY[l_count0*24+k_count22]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale22, int_scale, &haarResponseX[l_count1*24+k_count22], &haarResponseY[l_count1*24+k_count22]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale0, sample_x_sub_int_scale23, int_scale, &haarResponseX[l_count0*24+k_count23], &haarResponseY[l_count0*24+k_count23]);
                haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale1, sample_x_sub_int_scale23, int_scale, &haarResponseX[l_count1*24+k_count23], &haarResponseY[l_count1*24+k_count23]);
            }
        }

    } else {

        __m256 kks0 = _mm256_mul_ps(kk0_inner_arrays_simd, scale_vec);
        __m256 kks1 = _mm256_mul_ps(kk1_inner_arrays_simd, scale_vec);
        __m256 kks2 = _mm256_mul_ps(kk2_inner_arrays_simd, scale_vec);

        __m256 ipoint_x_sub_int_scale_add_05_vec = _mm256_set1_ps(ipoint_x_sub_int_scale_add_05);

        // USE CVTTPS_EPI32 FOR TRUNCATION!!!
        __m256i sample_col0 = _mm256_cvttps_epi32(_mm256_add_ps(ipoint_x_sub_int_scale_add_05_vec, kks0));
        __m256i sample_col1 = _mm256_cvttps_epi32(_mm256_add_ps(ipoint_x_sub_int_scale_add_05_vec, kks1));
        __m256i sample_col2 = _mm256_cvttps_epi32(_mm256_add_ps(ipoint_x_sub_int_scale_add_05_vec, kks2));

        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {

            int l_count1 = l_count + 1;

            int sample_y_sub_int_scale0 = (int)(ipoint_y_sub_int_scale_add_05 + l * scale);
            int sample_y_sub_int_scale1 = (int)(ipoint_y_sub_int_scale_add_05 + (l+1) * scale);  

            __m256i sample_y_sub_int_scale_vec0 = _mm256_set1_epi32(sample_y_sub_int_scale0);
            __m256i sample_y_sub_int_scale_vec1 = _mm256_set1_epi32(sample_y_sub_int_scale1);
            
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec0, sample_col0, int_scale_vec, &haarResponseX[l_count*24+0], &haarResponseY[l_count*24+0]);
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec1, sample_col0, int_scale_vec, &haarResponseX[l_count1*24+0], &haarResponseY[l_count1*24+0]);
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec0, sample_col1, int_scale_vec, &haarResponseX[l_count*24+8], &haarResponseY[l_count*24+8]);
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec1, sample_col1, int_scale_vec, &haarResponseX[l_count1*24+8], &haarResponseY[l_count1*24+8]);
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec0, sample_col2, int_scale_vec, &haarResponseX[l_count*24+16], &haarResponseY[l_count*24+16]);
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec1, sample_col2, int_scale_vec, &haarResponseX[l_count1*24+16], &haarResponseY[l_count1*24+16]);

        }
        
    }

    float s0  = roundf( 0.5 * scale);
    float s1  = roundf( 1.5 * scale);
    float s2  = roundf( 2.5 * scale);
    float s3  = roundf( 3.5 * scale);
    float s4  = roundf( 4.5 * scale);
    float s5  = roundf( 5.5 * scale);
    float s6  = roundf( 6.5 * scale);
    float s7  = roundf( 7.5 * scale);
    float s8  = roundf( 8.5 * scale);
    float s9  = roundf( 9.5 * scale);
    float s10 = roundf(10.5 * scale);
    float s11 = roundf(11.5 * scale);

    float e_c0_m4 = s2 + s1; // CAREFUL HERE!
    float e_c0_m3 = s2 + s0; // CAREFUL HERE!
    float e_c0_m2 = s2 - s0;
    float e_c0_m1 = s2 - s1;
    //float e_c0_z0 = s2 - s2;
    float e_c0_p1 = s2 - s3;
    float e_c0_p2 = s2 - s4;
    float e_c0_p3 = s2 - s5;
    float e_c0_p4 = s2 - s6;

    float e_c1_m4 = s7 - s3;
    float e_c1_m3 = s7 - s4;
    float e_c1_m2 = s7 - s5;
    float e_c1_m1 = s7 - s6;
    //float e_c1_z0 = s7 - s7;
    float e_c1_p1 = s7 - s8;
    float e_c1_p2 = s7 - s9;
    float e_c1_p3 = s7 - s10;
    float e_c1_p4 = s7 - s11;

    gauss_s1_c0[0] =  expf(g1_factor * (e_c0_m4 * e_c0_m4));
    gauss_s1_c0[1] =  expf(g1_factor * (e_c0_m3 * e_c0_m3));
    gauss_s1_c0[2] =  expf(g1_factor * (e_c0_m2 * e_c0_m2));
    gauss_s1_c0[3] =  expf(g1_factor * (e_c0_m1 * e_c0_m1));
    gauss_s1_c0[4] =  1.0f; //expf(g1_factor * (e_c0_z0 * e_c0_z0));
    gauss_s1_c0[5] =  expf(g1_factor * (e_c0_p1 * e_c0_p1));
    gauss_s1_c0[6] =  expf(g1_factor * (e_c0_p2 * e_c0_p2));
    gauss_s1_c0[7] =  expf(g1_factor * (e_c0_p3 * e_c0_p3));
    gauss_s1_c0[8] =  expf(g1_factor * (e_c0_p4 * e_c0_p4));

    gauss_s1_c1[0] =  expf(g1_factor * (e_c1_m4 * e_c1_m4));
    gauss_s1_c1[1] =  expf(g1_factor * (e_c1_m3 * e_c1_m3));
    gauss_s1_c1[2] =  expf(g1_factor * (e_c1_m2 * e_c1_m2));
    gauss_s1_c1[3] =  expf(g1_factor * (e_c1_m1 * e_c1_m1));
    gauss_s1_c1[4] =  1.0f; //expf(g1_factor * (e_c1_z0 * e_c1_z0));
    gauss_s1_c1[5] =  expf(g1_factor * (e_c1_p1 * e_c1_p1));
    gauss_s1_c1[6] =  expf(g1_factor * (e_c1_p2 * e_c1_p2));
    gauss_s1_c1[7] =  expf(g1_factor * (e_c1_p3 * e_c1_p3));
    gauss_s1_c1[8] =  expf(g1_factor * (e_c1_p4 * e_c1_p4));
    
    // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        float *i_arr;
        if (i == -8 || i == 7) {
            i_arr = gauss_s1_c1;
        } else {
            i_arr = gauss_s1_c0;
        } 

        int gauss_index_k_start = (i<0?8:0);
        int gauss_index_k_inc = (i<0?-1:1);

        for (int j=-8; j<8; j+=5) {

            float *j_arr;
            if (j == -8 || j == 7) {
                j_arr = gauss_s1_c1;
            } else {
                j_arr = gauss_s1_c0;
            } 

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int gauss_index_l = (j<0?8:0);
            int gauss_index_l_inc = (j<0?-1:1);
            for (int l = (j+8)*24; l < (j + 17)*24; l+=24, gauss_index_l+=gauss_index_l_inc) {
                float gauss_s1_y = j_arr[gauss_index_l];

                int gauss_index_k = gauss_index_k_start;
                for (int k = i+8; k < i + 17; ++k, gauss_index_k+=gauss_index_k_inc) {
                    
                    float gauss_s1_x = i_arr[gauss_index_k];

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = haarResponseX[l+k];
                    float ry = haarResponseY[l+k];
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_arr[gauss_s2_index];
            gauss_s2_index++;

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


void get_msurf_descriptors_simd_2_24(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_simd_2_24(iimage, &interest_points->at(i));
	}
}


void get_msurf_descriptor_simd_2_24_unconditional(struct integral_image* iimage, struct interest_point* ipoint) {

    float scale = ipoint->scale;
    // float scale_mul_25f = 2.5f*scale;
    int int_scale = (int) roundf(scale);

    __m256i int_scale_vec = _mm256_set1_epi32(int_scale);
    __m256 scale_vec = _mm256_set1_ps(scale);

    // int int_scale_mul_2 =  2 * int_scale;
    float scale_squared = scale*scale;

    float g1_factor = -0.08f / (scale_squared); // since 0.08f / (scale*scale) == 1.0f / (2.0f * 2.5f * scale * 2.5f * scale)
    //float g2_factor = -1.0f / 4.5f; // since 1.0f / 4.5f == 1.0f / (2.0f * 1.5f * 1.5f)

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

    float ipoint_x_sub_int_scale = ipoint_x-int_scale;
    float ipoint_y_sub_int_scale = ipoint_y-int_scale;

    float ipoint_x_sub_int_scale_add_05 = ipoint_x-int_scale + 0.5;
    float ipoint_y_sub_int_scale_add_05 = ipoint_y-int_scale + 0.5;
    
    int width = iimage->width;
    int height = iimage->height;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    //r(2.5*s) + r(1.5*s) == -(r(2.5*s) - r(6.5*s))

    // check if we ever hit a boundary
   
        __m256 kks0 = _mm256_mul_ps(kk0_inner_arrays_simd, scale_vec);
        __m256 kks1 = _mm256_mul_ps(kk1_inner_arrays_simd, scale_vec);
        __m256 kks2 = _mm256_mul_ps(kk2_inner_arrays_simd, scale_vec);

        __m256 ipoint_x_sub_int_scale_vec = _mm256_set1_ps((float) ipoint_x_sub_int_scale);

        // USE CVTTPS_EPI32 FOR TRUNCATION!!!
        __m256i sample_col1 = _mm256_cvtps_epi32(_mm256_add_ps(ipoint_x_sub_int_scale_vec, kks1));
        __m256i sample_col0 = _mm256_cvtps_epi32(_mm256_add_ps(ipoint_x_sub_int_scale_vec, kks0));
        __m256i sample_col2 = _mm256_cvtps_epi32(_mm256_add_ps(ipoint_x_sub_int_scale_vec, kks2));

        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {

            int l_count1 = l_count + 1;

            int sample_y_sub_int_scale0 = (int) roundf(ipoint_y_sub_int_scale + l * scale);
            int sample_y_sub_int_scale1 = (int) roundf(ipoint_y_sub_int_scale + (l+1) * scale);  

            __m256i sample_y_sub_int_scale_vec0 = _mm256_set1_epi32(sample_y_sub_int_scale0);
            __m256i sample_y_sub_int_scale_vec1 = _mm256_set1_epi32(sample_y_sub_int_scale1);
            
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec0, sample_col0, int_scale_vec, &haarResponseX[l_count*24+0], &haarResponseY[l_count*24+0]);
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec1, sample_col0, int_scale_vec, &haarResponseX[l_count1*24+0], &haarResponseY[l_count1*24+0]);
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec0, sample_col1, int_scale_vec, &haarResponseX[l_count*24+8], &haarResponseY[l_count*24+8]);
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec1, sample_col1, int_scale_vec, &haarResponseX[l_count1*24+8], &haarResponseY[l_count1*24+8]);
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec0, sample_col2, int_scale_vec, &haarResponseX[l_count*24+16], &haarResponseY[l_count*24+16]);
            haarXY_unconditional_vectorized(iimage, sample_y_sub_int_scale_vec1, sample_col2, int_scale_vec, &haarResponseX[l_count1*24+16], &haarResponseY[l_count1*24+16]);
        }



    float s0  = roundf( 0.5 * scale);
    float s1  = roundf( 1.5 * scale);
    float s2  = roundf( 2.5 * scale);
    float s3  = roundf( 3.5 * scale);
    float s4  = roundf( 4.5 * scale);
    float s5  = roundf( 5.5 * scale);
    float s6  = roundf( 6.5 * scale);
    float s7  = roundf( 7.5 * scale);
    float s8  = roundf( 8.5 * scale);
    float s9  = roundf( 9.5 * scale);
    float s10 = roundf(10.5 * scale);
    float s11 = roundf(11.5 * scale);

    float e_c0_m4 = s2 + s1; // CAREFUL HERE!
    float e_c0_m3 = s2 + s0; // CAREFUL HERE!
    float e_c0_m2 = s2 - s0;
    float e_c0_m1 = s2 - s1;
    //float e_c0_z0 = s2 - s2;
    float e_c0_p1 = s2 - s3;
    float e_c0_p2 = s2 - s4;
    float e_c0_p3 = s2 - s5;
    float e_c0_p4 = s2 - s6;

    float e_c1_m4 = s7 - s3;
    float e_c1_m3 = s7 - s4;
    float e_c1_m2 = s7 - s5;
    float e_c1_m1 = s7 - s6;
    //float e_c1_z0 = s7 - s7;
    float e_c1_p1 = s7 - s8;
    float e_c1_p2 = s7 - s9;
    float e_c1_p3 = s7 - s10;
    float e_c1_p4 = s7 - s11;

    gauss_s1_c0[0] =  expf(g1_factor * (e_c0_m4 * e_c0_m4));
    gauss_s1_c0[1] =  expf(g1_factor * (e_c0_m3 * e_c0_m3));
    gauss_s1_c0[2] =  expf(g1_factor * (e_c0_m2 * e_c0_m2));
    gauss_s1_c0[3] =  expf(g1_factor * (e_c0_m1 * e_c0_m1));
    gauss_s1_c0[4] =  1.0f; //expf(g1_factor * (e_c0_z0 * e_c0_z0));
    gauss_s1_c0[5] =  expf(g1_factor * (e_c0_p1 * e_c0_p1));
    gauss_s1_c0[6] =  expf(g1_factor * (e_c0_p2 * e_c0_p2));
    gauss_s1_c0[7] =  expf(g1_factor * (e_c0_p3 * e_c0_p3));
    gauss_s1_c0[8] =  expf(g1_factor * (e_c0_p4 * e_c0_p4));

    gauss_s1_c1[0] =  expf(g1_factor * (e_c1_m4 * e_c1_m4));
    gauss_s1_c1[1] =  expf(g1_factor * (e_c1_m3 * e_c1_m3));
    gauss_s1_c1[2] =  expf(g1_factor * (e_c1_m2 * e_c1_m2));
    gauss_s1_c1[3] =  expf(g1_factor * (e_c1_m1 * e_c1_m1));
    gauss_s1_c1[4] =  1.0f; //expf(g1_factor * (e_c1_z0 * e_c1_z0));
    gauss_s1_c1[5] =  expf(g1_factor * (e_c1_p1 * e_c1_p1));
    gauss_s1_c1[6] =  expf(g1_factor * (e_c1_p2 * e_c1_p2));
    gauss_s1_c1[7] =  expf(g1_factor * (e_c1_p3 * e_c1_p3));
    gauss_s1_c1[8] =  expf(g1_factor * (e_c1_p4 * e_c1_p4));
    
    // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        float *i_arr;
        if (i == -8 || i == 7) {
            i_arr = gauss_s1_c1;
        } else {
            i_arr = gauss_s1_c0;
        } 

        int gauss_index_k_start = (i<0?8:0);
        int gauss_index_k_inc = (i<0?-1:1);

        for (int j=-8; j<8; j+=5) {

            float *j_arr;
            if (j == -8 || j == 7) {
                j_arr = gauss_s1_c1;
            } else {
                j_arr = gauss_s1_c0;
            } 

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int gauss_index_l = (j<0?8:0);
            int gauss_index_l_inc = (j<0?-1:1);
            for (int l = (j+8)*24; l < (j + 17)*24; l+=24, gauss_index_l+=gauss_index_l_inc) {
                float gauss_s1_y = j_arr[gauss_index_l];

                int gauss_index_k = gauss_index_k_start;
                for (int k = i+8; k < i + 17; ++k, gauss_index_k+=gauss_index_k_inc) {
                    
                    float gauss_s1_x = i_arr[gauss_index_k];

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = haarResponseX[l+k];
                    float ry = haarResponseY[l+k];
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_arr[gauss_s2_index];
            gauss_s2_index++;

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


void get_msurf_descriptors_simd_2_24_unconditional(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_simd_2_24_unconditional(iimage, &interest_points->at(i));
	}
}