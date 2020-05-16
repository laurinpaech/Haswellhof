#include "descriptor.h"
#include "descriptor_opt.h"
#include "descriptor_opt_gen.h"
#include "integral_image.h"
#include "interest_point.h"

#define _USE_MATH_DEFINES

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>

extern float * haarResponseX;
extern float * haarResponseY;

void get_msurf_descriptor_haar_unroll_1_1_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_1_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_1_1_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_1_1_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_1_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_1_1_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_1_2_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_1_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_1_2_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_1_2_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_1_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_1_2_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_1_3_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_1_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_1_3_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_1_3_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_1_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_1_3_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_1_4_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                float ipoint_x_sub_int_scale_add_k3_mul_scale = ipoint_x_sub_int_scale + k3 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_k3_mul_scale + (ipoint_x_sub_int_scale_add_k3_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_05 + k3 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_1_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_1_4_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_1_4_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                float ipoint_x_sub_int_scale_add_k3_mul_scale = ipoint_x_sub_int_scale + k3 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_k3_mul_scale + (ipoint_x_sub_int_scale_add_k3_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=1, l_count+=1) {
            int l0 = l + 0;
            int l_count0 = l_count + 0;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_05 + k3 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_1_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_1_4_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_2_1_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
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

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_2_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_2_1_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_2_1_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
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

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_2_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_2_1_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_2_2_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
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

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_2_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_2_2_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_2_2_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
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

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_2_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_2_2_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_2_3_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
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

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_2_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_2_3_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_2_3_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
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

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_2_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_2_3_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_2_4_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                float ipoint_x_sub_int_scale_add_k3_mul_scale = ipoint_x_sub_int_scale + k3 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_k3_mul_scale + (ipoint_x_sub_int_scale_add_k3_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
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

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_05 + k3 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_2_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_2_4_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_2_4_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                float ipoint_x_sub_int_scale_add_k3_mul_scale = ipoint_x_sub_int_scale + k3 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_k3_mul_scale + (ipoint_x_sub_int_scale_add_k3_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
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

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_05 + k3 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_2_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_2_4_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_3_1_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_3_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_3_1_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_3_1_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_3_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_3_1_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_3_2_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_3_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_3_2_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_3_2_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_3_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_3_2_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_3_3_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_3_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_3_3_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_3_3_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_3_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_3_3_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_3_4_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                float ipoint_x_sub_int_scale_add_k3_mul_scale = ipoint_x_sub_int_scale + k3 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_k3_mul_scale + (ipoint_x_sub_int_scale_add_k3_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count2*24+k_count3], &haarResponseY[l_count2*24+k_count3]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_05 + k3 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count2*24+k_count3], &haarResponseY[l_count2*24+k_count3]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_3_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_3_4_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_3_4_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                float ipoint_x_sub_int_scale_add_k3_mul_scale = ipoint_x_sub_int_scale + k3 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_k3_mul_scale + (ipoint_x_sub_int_scale_add_k3_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count2*24+k_count3], &haarResponseY[l_count2*24+k_count3]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=3, l_count+=3) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_05 + k3 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count2*24+k_count3], &haarResponseY[l_count2*24+k_count3]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_3_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_3_4_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_4_1_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            float ipoint_y_sub_int_scale_add_l3_mul_scale = ipoint_y_sub_int_scale + l3 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_l3_mul_scale + (ipoint_y_sub_int_scale_add_l3_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_05 + l3 * scale);

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_4_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_4_1_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_4_1_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            float ipoint_y_sub_int_scale_add_l3_mul_scale = ipoint_y_sub_int_scale + l3 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_l3_mul_scale + (ipoint_y_sub_int_scale_add_l3_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_05 + l3 * scale);

            for (int k=-12, k_count=0; k<12; k+=1, k_count+=1) {
                int k0 = k + 0;
                int k_count0 = k_count + 0;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_4_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_4_1_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_4_2_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            float ipoint_y_sub_int_scale_add_l3_mul_scale = ipoint_y_sub_int_scale + l3 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_l3_mul_scale + (ipoint_y_sub_int_scale_add_l3_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count3*24+k_count1], &haarResponseY[l_count3*24+k_count1]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_05 + l3 * scale);

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count3*24+k_count1], &haarResponseY[l_count3*24+k_count1]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_4_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_4_2_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_4_2_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            float ipoint_y_sub_int_scale_add_l3_mul_scale = ipoint_y_sub_int_scale + l3 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_l3_mul_scale + (ipoint_y_sub_int_scale_add_l3_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count3*24+k_count1], &haarResponseY[l_count3*24+k_count1]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_05 + l3 * scale);

            for (int k=-12, k_count=0; k<12; k+=2, k_count+=2) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count3*24+k_count1], &haarResponseY[l_count3*24+k_count1]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_4_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_4_2_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_4_3_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            float ipoint_y_sub_int_scale_add_l3_mul_scale = ipoint_y_sub_int_scale + l3 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_l3_mul_scale + (ipoint_y_sub_int_scale_add_l3_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count3*24+k_count1], &haarResponseY[l_count3*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count3*24+k_count2], &haarResponseY[l_count3*24+k_count2]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_05 + l3 * scale);

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count3*24+k_count1], &haarResponseY[l_count3*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count3*24+k_count2], &haarResponseY[l_count3*24+k_count2]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_4_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_4_3_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_4_3_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            float ipoint_y_sub_int_scale_add_l3_mul_scale = ipoint_y_sub_int_scale + l3 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_l3_mul_scale + (ipoint_y_sub_int_scale_add_l3_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count3*24+k_count1], &haarResponseY[l_count3*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count3*24+k_count2], &haarResponseY[l_count3*24+k_count2]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_05 + l3 * scale);

            for (int k=-12, k_count=0; k<12; k+=3, k_count+=3) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count3*24+k_count1], &haarResponseY[l_count3*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count3*24+k_count2], &haarResponseY[l_count3*24+k_count2]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_4_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_4_3_True(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_4_4_False(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            float ipoint_y_sub_int_scale_add_l3_mul_scale = ipoint_y_sub_int_scale + l3 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_l3_mul_scale + (ipoint_y_sub_int_scale_add_l3_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                float ipoint_x_sub_int_scale_add_k3_mul_scale = ipoint_x_sub_int_scale + k3 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_k3_mul_scale + (ipoint_x_sub_int_scale_add_k3_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count2*24+k_count3], &haarResponseY[l_count2*24+k_count3]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count3*24+k_count1], &haarResponseY[l_count3*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count3*24+k_count2], &haarResponseY[l_count3*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count3*24+k_count3], &haarResponseY[l_count3*24+k_count3]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_05 + l3 * scale);

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_05 + k3 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count2*24+k_count3], &haarResponseY[l_count2*24+k_count3]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count3*24+k_count1], &haarResponseY[l_count3*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count3*24+k_count2], &haarResponseY[l_count3*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count3*24+k_count3], &haarResponseY[l_count3*24+k_count3]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_4_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_4_4_False(iimage, &interest_points->at(i));
    }
}
void get_msurf_descriptor_haar_unroll_4_4_True(struct integral_image* iimage, struct interest_point* ipoint) {

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
    
    float *data = (float *)iimage->data;
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
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            float ipoint_y_sub_int_scale_add_l0_mul_scale = ipoint_y_sub_int_scale + l0 * scale;
            float ipoint_y_sub_int_scale_add_l1_mul_scale = ipoint_y_sub_int_scale + l1 * scale;
            float ipoint_y_sub_int_scale_add_l2_mul_scale = ipoint_y_sub_int_scale + l2 * scale;
            float ipoint_y_sub_int_scale_add_l3_mul_scale = ipoint_y_sub_int_scale + l3 * scale;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_l0_mul_scale + (ipoint_y_sub_int_scale_add_l0_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_l1_mul_scale + (ipoint_y_sub_int_scale_add_l1_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_l2_mul_scale + (ipoint_y_sub_int_scale_add_l2_mul_scale>=0 ? 0.5 : -0.5));
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_l3_mul_scale + (ipoint_y_sub_int_scale_add_l3_mul_scale>=0 ? 0.5 : -0.5));

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                float ipoint_x_sub_int_scale_add_k0_mul_scale = ipoint_x_sub_int_scale + k0 * scale;
                float ipoint_x_sub_int_scale_add_k1_mul_scale = ipoint_x_sub_int_scale + k1 * scale;
                float ipoint_x_sub_int_scale_add_k2_mul_scale = ipoint_x_sub_int_scale + k2 * scale;
                float ipoint_x_sub_int_scale_add_k3_mul_scale = ipoint_x_sub_int_scale + k3 * scale;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_k0_mul_scale + (ipoint_x_sub_int_scale_add_k0_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_k1_mul_scale + (ipoint_x_sub_int_scale_add_k1_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_k2_mul_scale + (ipoint_x_sub_int_scale_add_k2_mul_scale>=0 ? 0.5 : -0.5));
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_k3_mul_scale + (ipoint_x_sub_int_scale_add_k3_mul_scale>=0 ? 0.5 : -0.5));

                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count3*24+k_count1], &haarResponseY[l_count3*24+k_count1]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count3*24+k_count2], &haarResponseY[l_count3*24+k_count2]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count2*24+k_count3], &haarResponseY[l_count2*24+k_count3]);
                haarXY_precheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count3*24+k_count3], &haarResponseY[l_count3*24+k_count3]);
            }
        }
    } else {
        for (int l=-12, l_count=0; l<12; l+=4, l_count+=4) {
            int l0 = l + 0;
            int l1 = l + 1;
            int l2 = l + 2;
            int l3 = l + 3;
            int l_count0 = l_count + 0;
            int l_count1 = l_count + 1;
            int l_count2 = l_count + 2;
            int l_count3 = l_count + 3;
            int sample_y_sub_int_scale0 = (int) (ipoint_y_sub_int_scale_add_05 + l0 * scale);
            int sample_y_sub_int_scale1 = (int) (ipoint_y_sub_int_scale_add_05 + l1 * scale);
            int sample_y_sub_int_scale2 = (int) (ipoint_y_sub_int_scale_add_05 + l2 * scale);
            int sample_y_sub_int_scale3 = (int) (ipoint_y_sub_int_scale_add_05 + l3 * scale);

            for (int k=-12, k_count=0; k<12; k+=4, k_count+=4) {
                int k0 = k + 0;
                int k1 = k + 1;
                int k2 = k + 2;
                int k3 = k + 3;
                int k_count0 = k_count + 0;
                int k_count1 = k_count + 1;
                int k_count2 = k_count + 2;
                int k_count3 = k_count + 3;
                int sample_x_sub_int_scale0 = (int) (ipoint_x_sub_int_scale_add_05 + k0 * scale);
                int sample_x_sub_int_scale1 = (int) (ipoint_x_sub_int_scale_add_05 + k1 * scale);
                int sample_x_sub_int_scale2 = (int) (ipoint_x_sub_int_scale_add_05 + k2 * scale);
                int sample_x_sub_int_scale3 = (int) (ipoint_x_sub_int_scale_add_05 + k3 * scale);

                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count0*24+k_count0], &haarResponseY[l_count0*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count1*24+k_count0], &haarResponseY[l_count1*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count2*24+k_count0], &haarResponseY[l_count2*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale0, int_scale, &haarResponseX[l_count3*24+k_count0], &haarResponseY[l_count3*24+k_count0]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count0*24+k_count1], &haarResponseY[l_count0*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count1*24+k_count1], &haarResponseY[l_count1*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count2*24+k_count1], &haarResponseY[l_count2*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale1, int_scale, &haarResponseX[l_count3*24+k_count1], &haarResponseY[l_count3*24+k_count1]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count0*24+k_count2], &haarResponseY[l_count0*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count1*24+k_count2], &haarResponseY[l_count1*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count2*24+k_count2], &haarResponseY[l_count2*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale2, int_scale, &haarResponseX[l_count3*24+k_count2], &haarResponseY[l_count3*24+k_count2]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale0, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count0*24+k_count3], &haarResponseY[l_count0*24+k_count3]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale1, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count1*24+k_count3], &haarResponseY[l_count1*24+k_count3]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale2, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count2*24+k_count3], &haarResponseY[l_count2*24+k_count3]);
                haarXY_nocheck_boundaries(data, height, width, sample_y_sub_int_scale3, sample_x_sub_int_scale3, int_scale, &haarResponseX[l_count3*24+k_count3], &haarResponseY[l_count3*24+k_count3]);
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

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
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
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}

void get_msurf_descriptors_haar_unroll_4_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (size_t i=0; i<interest_points->size(); ++i) {
        get_msurf_descriptor_haar_unroll_4_4_True(iimage, &interest_points->at(i));
    }
}
