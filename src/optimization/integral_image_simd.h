#pragma once

#include "integral_image.h"

#include <stdint.h>
#include <immintrin.h>


// https://stackoverflow.com/questions/34066228/how-to-perform-uint32-float-conversion-with-sse?rq=1
inline __m128 _my_mm_cvtepu32_ps(const __m128i v)
{
    __m128i v2 = _mm_srli_epi32(v, 1);     // v2 = v / 2
    __m128i v1 = _mm_sub_epi32(v, v2);     // v1 = v - (v / 2)
    __m128 v2f = _mm_cvtepi32_ps(v2);
    __m128 v1f = _mm_cvtepi32_ps(v1);
    return _mm_add_ps(v2f, v1f); 
}

inline __m256 _my_mm256_cvtepu32_ps(const __m256i v)
{
    __m256i v2 = _mm256_srli_epi32(v, 1);     // v2 = v / 2
    __m256i v1 = _mm256_sub_epi32(v, v2);     // v1 = v - (v / 2)
    __m256 v2f = _mm256_cvtepi32_ps(v2);
    __m256 v1f = _mm256_cvtepi32_ps(v1);
    return _mm256_add_ps(v2f, v1f); 
}

void compute_padded_integral_img_new(float *gray_image, struct integral_image *iimage);

void compute_padded_integral_img_faster_alg(float *gray_image, struct integral_image * iimage);

// Computes the integral image given an integer image as input and returning the float integral image
void compute_integral_img_int(uint8_t *gray_image, struct integral_image *iimage);

// Computes the integral image given an integer image as input and returning the float integral image
void compute_integral_img_simd_int(uint8_t *gray_image, struct integral_image *iimage);

// Computes the integral image given an integer image as input and returning the float integral image
// Converting previous row to  (1.0f/255.0f) * float cast on the fly
void compute_integral_img_simd_early_cast_int(uint8_t *gray_image, struct integral_image *iimage);

void compute_padded_integral_img_int(uint8_t *gray_image, struct integral_image *iimage);

void compute_padded_integral_img_simd_early_cast_int(uint8_t *gray_image, struct integral_image *iimage);

const __m256 box_zeros_vec = _mm256_setzero_ps();
const __m256i box_ones_vec = _mm256_set1_epi32(1);

inline __m256 box_integral_unconditional_simd(struct integral_image *iimage, __m256i r, __m256i c, __m256i rs, __m256i cs) {
    
    float *data = iimage->data;

    __m256i data_width = _mm256_set1_epi32(iimage->data_width);

    // TODO: Maybe resolve data dependency here
    __m256i r0 = _mm256_sub_epi32(r, box_ones_vec);
    __m256i r1 = _mm256_add_epi32(r0, rs);
    __m256i c0 = _mm256_sub_epi32(c, box_ones_vec);
    __m256i c1 = _mm256_add_epi32(c0, cs);

    __m256i r0w = _mm256_mul_epi32(r0, data_width);
    __m256i r1w = _mm256_mul_epi32(r1, data_width);

    __m256i A_idx = _mm256_add_epi32(r0w, c0);
    __m256i B_idx = _mm256_add_epi32(r0w, c1);
    __m256i C_idx = _mm256_add_epi32(r1w, c0);
    __m256i D_idx = _mm256_add_epi32(r1w, c1);

    __m256 A = _mm256_i32gather_ps(data, A_idx, sizeof(float));
    __m256 B = _mm256_i32gather_ps(data, B_idx, sizeof(float));
    __m256 C = _mm256_i32gather_ps(data, C_idx, sizeof(float));
    __m256 D = _mm256_i32gather_ps(data, D_idx, sizeof(float));

    __m256 res = _mm256_sub_ps(_mm256_sub_ps(A, B), _mm256_add_ps(C, D));

    return _mm256_max_ps(res, box_zeros_vec);

}
