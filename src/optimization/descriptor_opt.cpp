#include "descriptor.h"
#include "descriptor_opt.h"
#include "integral_image.h"
#include "interest_point.h"

#define _USE_MATH_DEFINES

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>


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
                    haarXY_nocheck_boundaries(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &rx, &ry);
                    
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


void get_msurf_descriptor_gauss_compute_once_case(struct integral_image* iimage, struct interest_point* ipoint) {
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


void get_msurf_descriptors_gauss_compute_once_case(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_gauss_compute_once_case(iimage, &interest_points->at(i));
	}
}

float haarResponseX[24*24];
float haarResponseY[24*24];

void get_msurf_descriptor_gauss_pecompute_haar(struct integral_image* iimage, struct interest_point* ipoint) {
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
                haarXY_nocheck_boundaries(iimage, sample_y_sub_int_scale, sample_x_sub_int_scale, int_scale, &haarResponseX[l_count*24+k_count], &haarResponseY[l_count*24+k_count]);

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


void get_msurf_descriptors_gauss_pecompute_haar(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_gauss_pecompute_haar(iimage, &interest_points->at(i));
	}
}

