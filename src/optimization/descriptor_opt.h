#pragma once 

#include "integral_image.h"
#include "interest_point.h"
#include "helper.h"
#include <immintrin.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#include "stdio.h"
#include <vector>

// Alignment only works for powers of two!
#define ALIGN(x,a)              __ALIGN_MASK(x,(typeof(x))(a)-1)
#define __ALIGN_MASK(x,mask)    (((x)+(mask))&~(mask))

// https://stackoverflow.com/questions/38088732/explanation-to-aligned-malloc-implementation
// Alignment only works for powers of two!
inline void *aligned_malloc(size_t required_bytes, size_t alignment) {
    void *p1;  // original block
    void **p2; // aligned block
    int offset = alignment - 1 + sizeof(void *);
    if ((p1 = (void *)malloc(required_bytes + offset)) == NULL) {
       return NULL;
    }
    p2 = (void **)(((size_t)(p1) + offset) & ~(alignment - 1));
    p2[-1] = p1;
    return p2;
}

inline void aligned_free(void *p) {
    free(((void **)p)[-1]);
}


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



void get_msurf_descriptor_improved(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_improved(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_improved_flip(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_improved_flip(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_improved_flip_flip(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_improved_flip_flip(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_inlined(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_inlined(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_gauss_s1_separable_test(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_gauss_s1_separable_test(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_gauss_s2_precomputed(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_gauss_s2_precomputed(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_inlinedHaarWavelets(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_inlinedHaarWavelets(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_inlinedHaarWavelets_precheck_boundaries(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_inlinedHaarWavelets_precheck_boundaries(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_gauss_compute_once_case(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_gauss_compute_once_case(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_gauss_pecompute_haar(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_gauss_pecompute_haar(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_gauss_pecompute_haar_rounding(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_gauss_pecompute_haar_rounding(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_arrays(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_arrays(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_arrays_unconditional(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_arrays_unconditional(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_gauss_pecompute_haar_unroll(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_gauss_pecompute_haar_unroll(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_24_True_winner(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptor2_haar_unroll_2_24_True_winner(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_24_True_winner_unconditional(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptor2_haar_unroll_2_24_True_winner_unconditional(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_arrays_simd(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_arrays_simd(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);


inline void haarXY(struct integral_image *iimage, int row, int col, int scale, float* haarX, float* haarY) {
    
    // TODO: (Sebastian) fix for iimage with padding
    int width = iimage->width;
    int height = iimage->height;

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

    float *data = iimage->data;
    int data_width = iimage->data_width;

    if (r0 >= 0 && c0 >= 0)  {
        r0c0 = data[r0 * data_width + c0];
    }
    if (r0 >= 0 && c1 >= 0) {
        r0c1 = data[r0 * data_width + c1];
    }
    if (r0 >= 0 && c2 >= 0) {
        r0c2 = data[r0 * data_width + c2];
    }
    if (r1 >= 0 && c0 >= 0) {
        r1c0 = data[r1 * data_width + c0];
    }
    if (r1 >= 0 && c2 >= 0) {
        r1c2 = data[r1 * data_width + c2];
    }
    if (r2 >= 0 && c0 >= 0) {
        r2c0 = data[r2 * data_width + c0];
    }
    if (r2 >= 0 && c1 >= 0) {
        r2c1 = data[r2 * data_width + c1];
    }
    if (r2 >= 0 && c2 >= 0) {
        r2c2 = data[r2 * data_width + c2];
    }


    // *haarX = (r0c1 - r0c2 - r2c1 + r2c2) - (r0c0 - r0c1 - r2c0 + r2c1);
    // *haarY = (r1c0 - r1c2 - r2c0 + r2c2) - (r0c0 - r0c2 - r1c0 + r1c2);

    float r2c2_sub_r0c0 = r2c2 - r0c0;
    float r2c0_sub_r0c2 = r2c0 - r0c2;

    *haarX = 2*(r0c1 - r2c1) + r2c2_sub_r0c0 + r2c0_sub_r0c2;
    *haarY = 2*(r1c0 - r1c2) + r2c2_sub_r0c0 - r2c0_sub_r0c2;

    // for FMA
    // *haarX = (2*r0c1 + r2c0_sub_r0c2) - (2*r2c1 - r2c2_sub_r0c0);
    // *haarY = (2*r1c0 - r2c0_sub_r0c2) - (2*r1c2 - r2c2_sub_r0c0);
}

inline void haarXY_precheck_boundaries(struct integral_image *iimage, int row, int col, int scale, float* haarX, float* haarY) {

    int width = iimage->width;
    int height = iimage->height;
    // (row,col) is upper left corner of haar wavelet filter
    if (row <= 0 
        || col <= 0 
        || (row + 2*scale) > height 
        || (col + 2*scale) > width) {

        haarXY(iimage, row, col, scale, haarX, haarY);

        return;
        // wavelet filters that can not be applied completely could also be skipped
        // the result will deviate from the base implementation
        // but this is how the original surf implementation is handling it
    }
    
    // subtracting by one for row/col because row/col is inclusive.
    int r0 = row - 1;         
    int c0 = col - 1;         
    int r1 = row + scale - 1;  
    int c1 = col + scale - 1;   
    int r2 = row + 2*scale - 1;
    int c2 = col + 2*scale - 1;

    float *data = iimage->data;
    int data_width = iimage->data_width;

    float r0c0 = data[r0 * data_width + c0];
    float r0c1 = data[r0 * data_width + c1];
    float r0c2 = data[r0 * data_width + c2];
    float r1c0 = data[r1 * data_width + c0];
    float r1c2 = data[r1 * data_width + c2];
    float r2c0 = data[r2 * data_width + c0];
    float r2c1 = data[r2 * data_width + c1];
    float r2c2 = data[r2 * data_width + c2];

    // *haarX = (r0c1 - r0c2 - r2c1 + r2c2) - (r0c0 - r0c1 - r2c0 + r2c1);
    // *haarY = (r1c0 - r1c2 - r2c0 + r2c2) - (r0c0 - r0c2 - r1c0 + r1c2);

    float r2c2_sub_r0c0 = r2c2 - r0c0;
    float r2c0_sub_r0c2 = r2c0 - r0c2;

    *haarX = 2*(r0c1 - r2c1) + r2c2_sub_r0c0 + r2c0_sub_r0c2;
    *haarY = 2*(r1c0 - r1c2) + r2c2_sub_r0c0 - r2c0_sub_r0c2;

    // for FMA
    // *haarX = (2*r0c1 + r2c0_sub_r0c2) - (2*r2c1 - r2c2_sub_r0c0);
    // *haarY = (2*r1c0 - r2c0_sub_r0c2) - (2*r1c2 - r2c2_sub_r0c0);

}

inline void haarXY_unconditional(struct integral_image *iimage, int row, int col, int scale, float* haarX, float* haarY) {
    
    float *data = iimage->data;
    int data_width = iimage->data_width;
    
    // subtracting by one for row/col because row/col is inclusive.
    int r0 = row - 1;         
    int c0 = col - 1;         
    int r1 = row + scale - 1;  
    int c1 = col + scale - 1;   
    int r2 = row + 2*scale - 1;
    int c2 = col + 2*scale - 1;

    float r0c0 = data[r0 * data_width + c0];
    float r0c1 = data[r0 * data_width + c1];
    float r0c2 = data[r0 * data_width + c2];
    float r1c0 = data[r1 * data_width + c0];
    float r1c2 = data[r1 * data_width + c2];
    float r2c0 = data[r2 * data_width + c0];
    float r2c1 = data[r2 * data_width + c1];
    float r2c2 = data[r2 * data_width + c2];

    // *haarX = (r0c1 - r0c2 - r2c1 + r2c2) - (r0c0 - r0c1 - r2c0 + r2c1);
    // *haarY = (r1c0 - r1c2 - r2c0 + r2c2) - (r0c0 - r0c2 - r1c0 + r1c2);

    float r2c2_sub_r0c0 = r2c2 - r0c0;
    float r2c0_sub_r0c2 = r2c0 - r0c2;

    *haarX = 2*(r0c1 - r2c1) + r2c2_sub_r0c0 + r2c0_sub_r0c2;
    *haarY = 2*(r1c0 - r1c2) + r2c2_sub_r0c0 - r2c0_sub_r0c2;

    // for FMA
    // *haarX = (2*r0c1 + r2c0_sub_r0c2) - (2*r2c1 - r2c2_sub_r0c0);
    // *haarY = (2*r1c0 - r2c0_sub_r0c2) - (2*r1c2 - r2c2_sub_r0c0);

}



inline void haarXY_unconditional_vectorized(struct integral_image *iimage, __m256i row, __m256i col, __m256i scale, float* haarX, float* haarY) {

    // float *data = iimage->data;
    // int data_width = iimage->data_width;

    float *data = iimage->data;

    printf("data width: %d%", iimage->data_width);

    __m256 twos = _mm256_set1_ps(2.0f);
    __m256i ones = _mm256_set1_epi32(1);
    int data_width_int = iimage->data_width;
    __m256i data_width = _mm256_set1_epi32(data_width_int);

    __m256i r0 = _mm256_sub_epi32(row, ones);
    __m256i c0 = _mm256_sub_epi32(col, ones);
    __m256i r1 = _mm256_add_epi32(r0, scale);
    __m256i c1 = _mm256_add_epi32(c0, scale);
    __m256i r2 = _mm256_add_epi32(r1, scale);
    __m256i c2 = _mm256_add_epi32(c1, scale);

    __m256i r0w = _mm256_mul_epi32(r0, data_width);
    __m256i r1w = _mm256_mul_epi32(r1, data_width);
    __m256i r2w = _mm256_mul_epi32(r2, data_width);
    
    // subtracting by one for row/col because row/col is inclusive.wenn
    // int r0 = row - 1;         
    // int c0 = col - 1;         
    // int r1 = row + scale - 1;  
    // int c1 = col + scale - 1;   
    // int r2 = row + 2*scale - 1;
    // int c2 = col + 2*scale - 1;

    __m256i r0c0_idx = _mm256_add_epi32(r0w, c0);
    __m256i r0c1_idx = _mm256_add_epi32(r0w, c1);
    __m256i r0c2_idx = _mm256_add_epi32(r0w, c2);
    __m256i r1c0_idx = _mm256_add_epi32(r1w, c0);
    __m256i r1c2_idx = _mm256_add_epi32(r1w, c2);
    __m256i r2c0_idx = _mm256_add_epi32(r2w, c0);
    __m256i r2c1_idx = _mm256_add_epi32(r2w, c1);
    __m256i r2c2_idx = _mm256_add_epi32(r2w, c2);
    

    __m256 r0c0 = _mm256_i32gather_ps(data, r0c0_idx, 4);
    __m256 r0c1 = _mm256_i32gather_ps(data, r0c1_idx, 4);
    __m256 r0c2 = _mm256_i32gather_ps(data, r0c2_idx, 4);
    __m256 r1c0 = _mm256_i32gather_ps(data, r1c0_idx, 4);
    __m256 r1c2 = _mm256_i32gather_ps(data, r1c2_idx, 4);
    __m256 r2c0 = _mm256_i32gather_ps(data, r2c0_idx, 4);
    __m256 r2c1 = _mm256_i32gather_ps(data, r2c1_idx, 4);
    __m256 r2c2 = _mm256_i32gather_ps(data, r2c2_idx, 4);

    // *haarX = (r0c1 - r0c2 - r2c1 + r2c2) - (r0c0 - r0c1 - r2c0 + r2c1);
    // *haarY = (r1c0 - r1c2 - r2c0 + r2c2) - (r0c0 - r0c2 - r1c0 + r1c2);

    // float r2c2_sub_r0c0 = r2c2 - r0c0;
    // float r2c0_sub_r0c2 = r2c0 - r0c2;

    // *haarX = 2*(r0c1 - r2c1) + r2c2_sub_r0c0 + r2c0_sub_r0c2;
    // *haarY = 2*(r1c0 - r1c2) + r2c2_sub_r0c0 - r2c0_sub_r0c2;

    __m256 r2c2_sub_r0c0 = _mm256_sub_ps(r2c2, r0c0);
    __m256 r2c0_sub_r0c2 = _mm256_sub_ps(r2c0, r0c2);
    __m256 r0c1_sub_r2c1 = _mm256_sub_ps(r0c1, r2c1);
    __m256 r1c0_sub_r1c2 = _mm256_sub_ps(r1c0, r1c2);
    __m256 r0c1_sub_r2c1_2 = _mm256_mul_ps(r0c1_sub_r2c1, twos);
    __m256 r1c0_sub_r1c2_2 = _mm256_mul_ps(r1c0_sub_r1c2, twos);

    __m256 haarX_val = _mm256_add_ps(r0c1_sub_r2c1_2, _mm256_add_ps(r2c2_sub_r0c0, r2c0_sub_r0c2));
    __m256 haarY_val = _mm256_add_ps(r1c0_sub_r1c2_2, _mm256_sub_ps(r2c2_sub_r0c0, r2c0_sub_r0c2));

    _mm256_storeu_ps(haarX, r2c2);
    _mm256_storeu_ps(haarY, r2c2);

}
