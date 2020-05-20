#include "integral_image_simd.h"

#include <stdlib.h>
#include <immintrin.h>

void compute_padded_integral_img_new(float *gray_image, struct integral_image *iimage) {

    float *iimage_padded_data = iimage->padded_data;

    float *iimage_data = iimage->data;
   
    int data_width = iimage->data_width;
    int data_height = iimage->data_height;
    int width = iimage->width;
    int height = iimage->height;

    int border = (data_width - width) / 2; 
    //int border = PADDING_SIZE;

    // Block layout of padded integral image:
    // AAA
    // BIC
    // BDE

    // BLOCK A
    for (int j = 0; j < border; ++j) {
        for (int i = 0; i < data_width; ++i) {
            iimage_padded_data[j * data_width + i] = 0.0f;
        }
    }

    // BLOCK B
    for (int j = border; j < data_height; ++j) {
        for (int i = 0; i < border; ++i) {
            iimage_padded_data[j * data_width + i] = 0.0f;
        }
    }

    // BLOCK I
    {
        float row_sum = 0.0f;

        /* sum up the first row */
        for (int i = 0; i < width; ++i) {
            /* previous rows are 0 */
            row_sum += gray_image[i];
            iimage_data[i] = row_sum;
        }

        /* sum all remaining rows*/
        for (int j = 1; j < height; ++j) {
            row_sum = 0.0f;
            for (int i = 0; i < width; ++i) {
                row_sum += gray_image[j * width + i];
                /*add sum of current row until current idx to sum of all previous rows until current index */
                iimage_data[j * data_width + i] = row_sum + iimage_data[(j - 1) * data_width + i];
            }
        }
    
    }

    // BLOCK C
    for (int j = border; j < border + height; ++j) {
        float last_element_in_row = iimage_padded_data[j * data_width + border + width - 1]; // TODO: Check if this is correct
        for (int i = width + border; i < data_width; ++i) {
            iimage_padded_data[j * data_width + i] = last_element_in_row;
        }
    }

    // BLOCK D
    for (int j = border + height; j < data_height; ++j) {
        for (int i = border; i < width + border; ++i) {
            iimage_padded_data[j * data_width + i] = iimage_padded_data[(j-1) * data_width + i];
        }
    }

    // BLOCK E
    int index_last_element = (border + height - 1) * data_width + border + width - 1;
    float max_value = iimage_padded_data[index_last_element];
    for (int j = border + height; j < data_height; ++j) {
        for (int i = border + width; i < data_width; ++i) {
            iimage_padded_data[j * data_width + i] = max_value;
        }
    }

}

void compute_padded_integral_img_faster_alg(float *gray_image, struct integral_image * iimage) {
    
    float *iimage_padded_data = iimage->padded_data;

    int data_width = iimage->data_width;
    int data_height = iimage->data_height;
    int width = iimage->width;
    int height = iimage->height;

    int border = (data_width - width) / 2; 
    //int border = PADDING_SIZE;

    // Block layout of padded integral image:
    // AAA
    // BIC
    // BDE

    __m256 zeros_float = _mm256_setzero_ps();

    // BLOCK A
    for (int j = 0; j < border; ++j) {
        int i = 0;
        for (; i < data_width - 7; i+=8) {
            _mm256_storeu_ps(iimage_padded_data + j * data_width + i, zeros_float);
        }
        for (; i < data_width; ++i) {
            iimage_padded_data[j * data_width + i] = 0.0f;
        }
    }

    // BLOCK B
    for (int j = border; j < data_height; ++j) {
        int i = 0;
        for (; i < border - 7; i+=8) {
            _mm256_storeu_ps(iimage_padded_data + j * data_width + i, zeros_float);
        }
        for (; i < border; ++i) {
            iimage_padded_data[j * data_width + i] = 0.0f;
        }
    }

    // BLOCK I
    {
        float sum = 0.0f;
        float sum1 = 0.0f;
        
        int data_width = iimage->data_width;
        int width = iimage->width;
        int width_limit = width - 2;
        int height_limit = iimage->height - 2;
        float *iimage_data = iimage->data;

        int i = 0;
        // first row extra, since we don't have 0 padding
        for (i = 0; i < width_limit; i += 2) {
            sum += gray_image[i];
            iimage_data[i] = sum;
            sum += gray_image[i + 1];
            iimage_data[i + 1] = sum;

            sum1 += gray_image[width + i];
            iimage_data[data_width + i] = iimage_data[i] + sum1;
            sum1 += gray_image[width + i + 1];
            iimage_data[data_width + i + 1] = iimage_data[i + 1] + sum1;
        }

        // Handle last element of an image with odd width extra.
        // If previous stride = 2, then there should only be one element covered in this loop.
        // This might change if we increase the stride.
        for (; i < width; ++i) {
            sum += gray_image[i];
            iimage_data[i] = sum;
            // printf("row 0, col %i sum: %f\n", i, sum);
            sum1 += gray_image[width + i];
            // printf("row 1, col %i sum: %f\n", i, sum1);
            iimage_data[data_width + i] = iimage_data[i] + sum1;
        }

        // Handle all elements after the first row.
        i = 2;
        for (; i < height_limit; i += 2) {
            sum = 0.0f;
            sum1 = 0.0f;
            int j = 0;
            for (; j < width_limit; j += 2) {
                sum += gray_image[i * width + j];
                iimage_data[i * data_width + j] = iimage_data[(i - 1) * data_width + j] + sum;
                sum += gray_image[i * width + j + 1];
                iimage_data[i * data_width + j + 1] = iimage_data[(i - 1) * data_width + j + 1] + sum;

                sum1 += gray_image[(i + 1) * width + j];
                iimage_data[(i + 1) * data_width + j] = iimage_data[i * data_width + j] + sum1;
                sum1 += gray_image[(i + 1) * width + j + 1];
                iimage_data[(i + 1) * data_width + j + 1] = iimage_data[i * data_width + j + 1] + sum1;
            }
            // printf("j: %i\n", j);

            // Handle last element of an image with odd width extra.
            // If previous stride = 2, then there should only be one element covered in this loop.
            // This might change if we increase the stride.
            for (; j < width; ++j) {
                // printf("j: %i; sum1: %f\n", j, sum1);
                sum += gray_image[i * width + j];
                iimage_data[i * data_width + j] = iimage_data[(i - 1) * data_width + j] + sum;

                sum1 += gray_image[(i + 1) * width + j];
                // printf("sum1: %f\n", sum1);
                iimage_data[(i + 1) * data_width + j] = iimage_data[i * data_width + j] + sum1;
            }
        }

        // Handle last element of an image with odd height extra.
        // If previous stride = 2, then there should only be one row covered in this loop.
        // This might change if we increase the stride.
        for (; i < iimage->height ; ++i) {
            sum = 0.0f;
            int j = 0;
            for (; j < width_limit; j += 2) {
                // printf("j: %i; sum1: %f\n", j, sum);

                sum += gray_image[i * width + j];
                // printf("j: %i; sum1: %f\n", j, sum);
                iimage_data[i * data_width + j] = iimage_data[(i - 1) * data_width + j] + sum;
                sum += gray_image[i * width + j + 1];
                // printf("j: %i; sum1: %f\n", j, sum);

                iimage_data[i * data_width + j + 1] = iimage_data[(i - 1) * data_width + j + 1] + sum;
            }
            // printf("j: %i\n", j);

            for (; j < width; ++j) {
                sum += gray_image[i * width + j];
                // printf("j: %i; sum1: %f\n", j, sum);

                iimage_data[i * data_width + j] = iimage_data[(i - 1) * data_width + j] + sum;
            }
        }
    }

    // BLOCK C
    for (int j = border; j < border + height; ++j) {
        float last_element_in_row_float = iimage_padded_data[j * data_width + border + width - 1];
        __m256 last_element_in_row_vec = _mm256_set1_ps(last_element_in_row_float);
        int i = width + border;
        for (; i < data_width - 7; i+=8) {
            _mm256_storeu_ps(iimage_padded_data + j * data_width + i, last_element_in_row_vec);
        }
        for (; i < data_width; ++i) {
            iimage_padded_data[j * data_width + i] = last_element_in_row_float;
        }
    }

    // BLOCK D   
    for (int j = border + height; j < data_height; ++j) {
        int i = border;
        for (; i < width + border - 7; i+=8) {
            __m256 prev_elements_in_row = _mm256_loadu_ps(iimage_padded_data + (j-1) * data_width + i);
            _mm256_storeu_ps(iimage_padded_data + j * data_width + i, prev_elements_in_row);
        }
        for (; i < width + border; ++i) {
            iimage_padded_data[j * data_width + i] = iimage_padded_data[(j-1) * data_width + i];
        }
    }

    // BLOCK E
    int index_last_element = (border + height - 1) * data_width + border + width - 1;
    float max_value_float = iimage_padded_data[index_last_element];
    __m256 max_value_vec = _mm256_set1_ps(max_value_float);
    for (int j = border + height; j < data_height; ++j) {
        int i = border + width;
        for (; i < data_width - 7; i+=8) {
            _mm256_storeu_ps(iimage_padded_data + j * data_width + i, max_value_vec);
        }
        for (; i < data_width; ++i) {
            iimage_padded_data[j * data_width + i] = max_value_float;
        }
    }
    
}


void compute_integral_img_int(uint8_t *gray_image, struct integral_image *iimage) {

    uint32_t *iimage_idata = (uint32_t *) iimage->data;

    int data_width = iimage->data_width;
    int width = iimage->width;
    int height = iimage->height;

    uint32_t row_sum = 0;

    /* sum up the first row */
    for (size_t i = 0; i < width; ++i) {
        /* previous rows are 0 */
        row_sum += gray_image[i];
        iimage_idata[i] = row_sum;
    }

    /* sum all remaining rows*/
    for (size_t j = 1; j < height; ++j) {
        row_sum = 0;
        for (size_t i = 0; i < width; ++i) {
            row_sum += gray_image[j * width + i];
            /*add sum of current row until current idx to sum of all previous rows until current index */
            iimage_idata[j * data_width + i] = row_sum + iimage_idata[(j - 1) * data_width + i];
        }
    }

    
    float *iimage_data = iimage->data;

    for (size_t j = 0; j < height; ++j) {
        for (size_t i = 0; i < width; ++i) {
            //iimage_data[j * data_width + i] = (1.0f / 255.0f) * ((float) iimage_idata[j * data_width + i]);
            iimage_data[j * data_width + i] = (1.0f / 255.0f) * ((float) iimage_idata[j * data_width + i]);
            //iimage_data[j * data_width + i] = ((float) iimage_idata[j * data_width + i]) / 255.0f;
        }
    }
    

}

void compute_padded_integral_img_int(uint8_t *gray_image, struct integral_image *iimage) {

    float *iimage_padded_data = iimage->padded_data;
   
    int data_width = iimage->data_width;
    int data_height = iimage->data_height;
    int width = iimage->width;
    int height = iimage->height;

    int border = (data_width - width) / 2; 
    //int border = PADDING_SIZE;

    // Block layout of padded integral image:
    // AAA
    // BIC
    // BDE

    // BLOCK A
    for (int j = 0; j < border; ++j) {
        for (int i = 0; i < data_width; ++i) {
            iimage_padded_data[j * data_width + i] = 0.0f;
        }
    }

    // BLOCK B
    for (int j = border; j < data_height; ++j) {
        for (int i = 0; i < border; ++i) {
            iimage_padded_data[j * data_width + i] = 0.0f;
        }
    }

    // BLOCK I
    {
        uint32_t *iimage_idata = (uint32_t *) iimage->data;

        uint32_t row_sum = 0;

        /* sum up the first row */
        for (int i = 0; i < width; ++i) {
            /* previous rows are 0 */
            row_sum += gray_image[i];
            iimage_idata[i] = row_sum;
        }

        /* sum all remaining rows*/
        for (int j = 1; j < height; ++j) {
            row_sum = 0;
            for (int i = 0; i < width; ++i) {
                row_sum += gray_image[j * width + i];
                /*add sum of current row until current idx to sum of all previous rows until current index */
                iimage_idata[j * data_width + i] = row_sum + iimage_idata[(j - 1) * data_width + i];
            }
        }

        // Converting uint32_t to floats by casting them and multiplying with (1.0f/255.0f)
        float *iimage_data = iimage->data;
        for (int j = 0; j < height; ++j) {
            for (int i = 0; i < width; ++i) {
                //iimage_data[j * data_width + i] = (1.0f / 255.0f) * ((float) iimage_idata[j * data_width + i]);
                iimage_data[j * data_width + i] = (1.0f / 255.0f) * ((float) iimage_idata[j * data_width + i]);
                //iimage_data[j * data_width + i] = ((float) iimage_idata[j * data_width + i]) / 255.0f;
            }
        }
    
    }

    // BLOCK C
    for (int j = border; j < border + height; ++j) {
        float last_element_in_row = iimage_padded_data[j * data_width + border + width - 1]; // TODO: Check if this is correct
        for (int i = width + border; i < data_width; ++i) {
            iimage_padded_data[j * data_width + i] = last_element_in_row;
        }
    }

    // BLOCK D
    for (int j = border + height; j < data_height; ++j) {
        for (int i = border; i < width + border; ++i) {
            iimage_padded_data[j * data_width + i] = iimage_padded_data[(j-1) * data_width + i];
        }
    }

    // BLOCK E
    int index_last_element = (border + height - 1) * data_width + border + width - 1;
    float max_value = iimage_padded_data[index_last_element];
    for (int j = border + height; j < data_height; ++j) {
        for (int i = border + width; i < data_width; ++i) {
            iimage_padded_data[j * data_width + i] = max_value;
        }
    }

}

// Computes the integral image
void compute_integral_img_simd_int(uint8_t *gray_image, struct integral_image *iimage) {
    
    // Implementation adapted from:
    // https://stackoverflow.com/questions/46520275/how-to-speed-up-calculation-of-integral-image
    // https://github.com/ermig1979/Simd/blob/master/src/Simd/SimdAvx2Integral.cpp

    int data_width = iimage->data_width;
    int width = iimage->width;
    int height = iimage->height;

    __m256i MASK = _mm256_setr_epi64x(0x00000000000000FF, 0x000000000000FFFF, 0x0000000000FFFFFF, 0x00000000FFFFFFFF);
    __m256i PACK = _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7);
    //__m256i ZERO = _mm256_set1_epi32(0);
    __m256i ZERO = _mm256_setzero_si256();

    uint32_t *iimage_idata = (uint32_t *) iimage->data;

    size_t aligned_width = width / 4 * 4;

    {
        __m256i row_sums = ZERO;
        size_t i = 0;
        for(; i < aligned_width; i+=4) {

            // If (gray_image + i) has following uint8_t pixels: a b c d
            // a0000000 ab000000 abc00000 abcd0000
            __m256i img_rep_mask = _mm256_and_si256(_mm256_set1_epi32(*(uint32_t*)(gray_image + i)), MASK);
            
            // Summing together blocks of 8 uint8_t's and storing results
            // uint32_t eeee = a
            // uint32_t ffff = a + b
            // uint32_t gggg = a + b + c
            // uint32_t hhhh = a + b + c + d
            // eeee0000 ffff0000 gggg0000 hhhh0000
            __m256i img_sad_sum = _mm256_sad_epu8(img_rep_mask, ZERO);

            // Updating row_sum by adding previously computed img_sad_sum
            //   eeee0000 ffff0000 hhhh0000 hhhh0000
            // + iiii**** jjjj**** kkkk**** llll****
            // -------------------------------------
            //   mmmm**** nnnn**** oooo**** pppp****
            row_sums = _mm256_add_epi32(row_sums, img_sad_sum);

            // Permuting row_sum from first row to second and cutting off upper 128 bits
            // mmmm**** nnnn**** oooo**** pppp****
            // mmmmnnnn oooopppp ******** ********
            __m128i curr_row_sums = _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(row_sums, PACK));
            
            // Storing 128 bits (4*uint32_t) in curr_row_sums back to memory
            // 
            _mm_storeu_si128((__m128i*)(iimage_idata + i), curr_row_sums); // For unaligned
            //_mm_store_si128((__m128i *)(iimage_idata + i), curr_row_sums);
            
            // Updating row_sums to only include highest of four values (pppp)
            // mmmm**** nnnn**** oooo**** pppp****
            // pppp**** pppp**** pppp**** pppp****
            row_sums = _mm256_permute4x64_epi64(row_sums, 0xFF);

        }
        uint32_t row_sum = iimage_idata[i - 1];
        for (; i < width; ++i) {
            row_sum += gray_image[i];
            iimage_idata[i] = row_sum;
        }
    }

    for(size_t j = 1; j < height; ++j) {

        __m256i row_sums = ZERO;

        size_t i = 0;
        for(; i < aligned_width; i+=4) {

            // If (gray_image + i) has following uint8_t pixels: a b c d
            // a0000000 ab000000 abc00000 abcd0000
            __m256i img_rep_mask = _mm256_and_si256(_mm256_set1_epi32(*(uint32_t*)(gray_image + j * width + i)), MASK);
            
            // Summing together blocks of 8 uint8_t's and storing results
            // uint32_t eeee = a
            // uint32_t ffff = a + b
            // uint32_t gggg = a + b + c
            // uint32_t hhhh = a + b + c + d
            // eeee0000 ffff0000 gggg0000 hhhh0000
            __m256i img_sad_sum = _mm256_sad_epu8(img_rep_mask, ZERO);

            // Updating row_sum by adding previously computed img_sad_sum
            //   eeee0000 ffff0000 hhhh0000 hhhh0000
            // + iiii**** jjjj**** kkkk**** llll****
            // -------------------------------------
            //   mmmm**** nnnn**** oooo**** pppp****
            row_sums = _mm256_add_epi32(row_sums, img_sad_sum);

            // Permuting row_sum from first row to second and cutting off upper 128 bits
            // mmmm**** nnnn**** oooo**** pppp****
            // mmmmnnnn oooopppp ******** ********
            __m128i curr_row_sums = _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(row_sums, PACK));

            // Loading row_sum from previous row
            __m128i prev_row_sums = _mm_loadu_si128((__m128i *)(iimage_idata + (j-1) * data_width + i));
            //__m128i prev_row_sums = _mm_load_si128((__m128i *)(iimage_idata + (j-1) * data_width + i));

            // Adding curr_row_sums with prev_row_sums
            __m128i add_row_sums = _mm_add_epi32(curr_row_sums, prev_row_sums);

            // Storing 128 bits (4*uint32_t) in curr_row_sums back to memory
            _mm_storeu_si128((__m128i *)(iimage_idata + j * data_width + i), add_row_sums); // For unaligned
            //_mm_store_si128((__m128i *)(iimage_idata + j * data_width + i), add_row_sums);
            
            // Updating row_sums to only include highest of four values (pppp)
            // mmmm**** nnnn**** oooo**** pppp****
            // pppp**** pppp**** pppp**** pppp****
            row_sums = _mm256_permute4x64_epi64(row_sums, 0xFF);

        }

        uint32_t row_sum = iimage_idata[j * data_width + i - 1] - iimage_idata[(j-1) * data_width + i - 1];
        for (; i < width; ++i) {
            row_sum += gray_image[j * width + i];
            iimage_idata[j * data_width + i] = row_sum + iimage_idata[(j-1) * data_width + i];
        }

    }

    
    float *iimage_data = iimage->data;

    __m256 factor = _mm256_set1_ps(1.0f / 255.0f);

    for (size_t j = 0; j < height; ++j) {

        size_t i = 0;
        for (; i < width - 7; i+=8) {

            __m256i v_int = _mm256_loadu_si256((__m256i *)(iimage_idata + j * data_width + i));
            //__m256i v_int = _mm256_load_si256((__m256i *)(iimage_idata + j * data_width + i));

            // AVX512 only :(
            //__m256 v_float = _mm256_cvtepu32_ps(v_int);
            __m256 v_float = _mm256_cvtepi32_ps(v_int);
            //__m256 v_float = _my_mm256_cvtepu32_ps(v_int);

            // Multiplying converted float with (1.0f/255.0f) 
            __m256 v_iimage = _mm256_mul_ps(factor, v_float);

            // Storing value back to memory
            _mm256_storeu_ps(iimage_data + j * data_width + i, v_iimage);

        } 

        for (; i < width; ++i) {
            iimage_data[j * data_width + i] = (1.0f / 255.0f) * ((float) iimage_idata[j * data_width + i]);
        }

    }
    

}

// Computes the integral image
void compute_integral_img_simd_early_cast_int(uint8_t *gray_image, struct integral_image *iimage) {
    
    // Implementation adapted from:
    // https://stackoverflow.com/questions/46520275/how-to-speed-up-calculation-of-integral-image
    // https://github.com/ermig1979/Simd/blob/master/src/Simd/SimdAvx2Integral.cpp

    int data_width = iimage->data_width;
    int width = iimage->width;
    int height = iimage->height;

    __m256i MASK = _mm256_setr_epi64x(0x00000000000000FF, 0x000000000000FFFF, 0x0000000000FFFFFF, 0x00000000FFFFFFFF);
    __m256i PACK = _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7);
    //__m256i ZERO = _mm256_set1_epi32(0);
    __m256i ZERO = _mm256_setzero_si256();

    uint32_t *iimage_idata = (uint32_t *) iimage->data;

    size_t aligned_width = width / 4 * 4;

    {
        __m256i row_sums = ZERO;
        size_t i = 0;
        for(; i < aligned_width; i+=4) {

            // If (gray_image + i) has following uint8_t pixels: a b c d
            // a0000000 ab000000 abc00000 abcd0000
            __m256i img_rep_mask = _mm256_and_si256(_mm256_set1_epi32(*(uint32_t*)(gray_image + i)), MASK);
            
            // Summing together blocks of 8 uint8_t's and storing results
            // uint32_t eeee = a
            // uint32_t ffff = a + b
            // uint32_t gggg = a + b + c
            // uint32_t hhhh = a + b + c + d
            // eeee0000 ffff0000 gggg0000 hhhh0000
            __m256i img_sad_sum = _mm256_sad_epu8(img_rep_mask, ZERO);

            // Updating row_sum by adding previously computed img_sad_sum
            //   eeee0000 ffff0000 hhhh0000 hhhh0000
            // + iiii**** jjjj**** kkkk**** llll****
            // -------------------------------------
            //   mmmm**** nnnn**** oooo**** pppp****
            row_sums = _mm256_add_epi32(row_sums, img_sad_sum);

            // Permuting row_sum from first row to second and cutting off upper 128 bits
            // mmmm**** nnnn**** oooo**** pppp****
            // mmmmnnnn oooopppp ******** ********
            __m128i curr_row_sums = _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(row_sums, PACK));
            
            // Storing 128 bits (4*uint32_t) in curr_row_sums back to memory
            // 
            _mm_storeu_si128((__m128i*)(iimage_idata + i), curr_row_sums); // For unaligned
            //_mm_store_si128((__m128i *)(iimage_idata + i), curr_row_sums);
            
            // Updating row_sums to only include highest of four values (pppp)
            // mmmm**** nnnn**** oooo**** pppp****
            // pppp**** pppp**** pppp**** pppp****
            row_sums = _mm256_permute4x64_epi64(row_sums, 0xFF);

        }
        uint32_t row_sum = iimage_idata[i - 1];
        for (; i < width; ++i) {
            row_sum += gray_image[i];
            iimage_idata[i] = row_sum;
        }
    }

    float *iimage_data = iimage->data;
    __m128 factor = _mm_set1_ps(1.0f / 255.0f);

    for(size_t j = 1; j < height; ++j) {

        __m256i row_sums = ZERO;

        size_t i = 0;
        for(; i < aligned_width; i+=4) {

            // If (gray_image + i) has following uint8_t pixels: a b c d
            // a0000000 ab000000 abc00000 abcd0000
            __m256i img_rep_mask = _mm256_and_si256(_mm256_set1_epi32(*(uint32_t*)(gray_image + j * width + i)), MASK);
            
            // Summing together blocks of 8 uint8_t's and storing results
            // uint32_t eeee = a
            // uint32_t ffff = a + b
            // uint32_t gggg = a + b + c
            // uint32_t hhhh = a + b + c + d
            // eeee0000 ffff0000 gggg0000 hhhh0000
            __m256i img_sad_sum = _mm256_sad_epu8(img_rep_mask, ZERO);

            // Updating row_sum by adding previously computed img_sad_sum
            //   eeee0000 ffff0000 hhhh0000 hhhh0000
            // + iiii**** jjjj**** kkkk**** llll****
            // -------------------------------------
            //   mmmm**** nnnn**** oooo**** pppp****
            row_sums = _mm256_add_epi32(row_sums, img_sad_sum);

            // Permuting row_sum from first row to second and cutting off upper 128 bits
            // mmmm**** nnnn**** oooo**** pppp****
            // mmmmnnnn oooopppp ******** ********
            __m128i curr_row_sums = _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(row_sums, PACK));

            // Loading row_sum from previous row
            __m128i prev_row_sums = _mm_loadu_si128((__m128i *)(iimage_idata + (j-1) * data_width + i));
            //__m128i prev_row_sums = _mm_load_si128((__m128i *)(iimage_idata + (j-1) * data_width + i));

            // AVX512 only :(
            //__m128 prev_row_sums_cast = _mm_cvtepu32_ps(prev_row_sums);
            __m128 prev_row_sums_cast = _mm_cvtepi32_ps(prev_row_sums);
            //__128 prev_row_sums_cast = _my_mm_cvtepu32_ps(prev_row_sums);

            // Multiplying float version of previous row sum converted with (1.0f/255.0f) 
            __m128 prev_iimage_float = _mm_mul_ps(factor, prev_row_sums_cast);

            // Storing float of previous row sum back to memory as final integral image value
            _mm_storeu_ps(iimage_data +(j-1) * data_width + i, prev_iimage_float);

            // Adding curr_row_sums with prev_row_sums
            __m128i add_row_sums = _mm_add_epi32(curr_row_sums, prev_row_sums);

            // Storing 128 bits (4*uint32_t) in curr_row_sums back to memory
            _mm_storeu_si128((__m128i *)(iimage_idata + j * data_width + i), add_row_sums); // For unaligned
            //_mm_store_si128((__m128i *)(iimage_idata + j * data_width + i), add_row_sums);
            
            // Updating row_sums to only include highest of four values (pppp)
            // mmmm**** nnnn**** oooo**** pppp****
            // pppp**** pppp**** pppp**** pppp****
            row_sums = _mm256_permute4x64_epi64(row_sums, 0xFF);

        }

        uint32_t row_sum = iimage_idata[j * data_width + i - 1] - iimage_idata[(j-1) * data_width + i - 1];
        for (; i < width; ++i) {
            row_sum += gray_image[j * width + i];
            float prev_row_sum = iimage_idata[(j-1) * data_width + i];
            iimage_idata[j * data_width + i] = row_sum + prev_row_sum;
            iimage_data[(j-1) * data_width + i] = (1.0f / 255.0f) * ((float) prev_row_sum);
        }

    }

    // Converting last row to floats
    __m256 factor256 = _mm256_set1_ps(1.0f / 255.0f);
    size_t i = 0;
    for (; i < width - 7; i+=8) {

        __m256i v_int = _mm256_loadu_si256((__m256i *)(iimage_idata + (height - 1) * data_width + i));
        //__m256i v_int = _mm256_load_si256((__m256i *)(iimage_idata + j * data_width + i));

        // AVX512 only :(
        //__m256 v_float = _mm256_cvtepu32_ps(v_int);
        __m256 v_float = _mm256_cvtepi32_ps(v_int);
        //__m256 v_float = _my_mm256_cvtepu32_ps(v_int);

        // Multiplying converted float with (1.0f/255.0f) 
        __m256 v_iimage = _mm256_mul_ps(factor256, v_float);

        // Storing value back to memory
        _mm256_storeu_ps(iimage_data + (height - 1) * data_width + i, v_iimage);

    } 

    for (; i < width; ++i) {
        iimage_data[(height - 1) * data_width + i] = (1.0f / 255.0f) * ((float) iimage_idata[(height - 1) * data_width + i]);
    }

}


// Computes the integral image
void compute_padded_integral_img_simd_early_cast_int(uint8_t *gray_image, struct integral_image *iimage) {
    
    float *iimage_padded_data = iimage->padded_data;

    int data_width = iimage->data_width;
    int data_height = iimage->data_height;
    int width = iimage->width;
    int height = iimage->height;

    int border = (data_width - width) / 2; 
    //int border = PADDING_SIZE;

    // Block layout of padded integral image:
    // AAA
    // BIC
    // BDE

    __m256 zeros_float = _mm256_setzero_ps();

    // BLOCK A
    for (int j = 0; j < border; ++j) {
        int i = 0;
        for (; i < data_width - 7; i+=8) {
            _mm256_storeu_ps(iimage_padded_data + j * data_width + i, zeros_float);
        }
        for (; i < data_width; ++i) {
            iimage_padded_data[j * data_width + i] = 0.0f;
        }
    }

    // BLOCK B
    for (int j = border; j < data_height; ++j) {
        int i = 0;
        for (; i < border - 7; i+=8) {
            _mm256_storeu_ps(iimage_padded_data + j * data_width + i, zeros_float);
        }
        for (; i < border; ++i) {
            iimage_padded_data[j * data_width + i] = 0.0f;
        }
    }

    // BLOCK I
    {
        __m256i MASK = _mm256_setr_epi64x(0x00000000000000FF, 0x000000000000FFFF, 0x0000000000FFFFFF, 0x00000000FFFFFFFF);
        __m256i PACK = _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7);
        //__m256i ZERO = _mm256_set1_epi32(0);
        __m256i ZERO = _mm256_setzero_si256();
        
        // Implementation adapted from:
        // https://stackoverflow.com/questions/46520275/how-to-speed-up-calculation-of-integral-image
        // https://github.com/ermig1979/Simd/blob/master/src/Simd/SimdAvx2Integral.cpp

        uint32_t *iimage_idata = (uint32_t *) iimage->data;

        size_t aligned_width = width / 4 * 4;

        {
            __m256i row_sums = ZERO;
            size_t i = 0;
            for(; i < aligned_width; i+=4) {

                // If (gray_image + i) has following uint8_t pixels: a b c d
                // a0000000 ab000000 abc00000 abcd0000
                __m256i img_rep_mask = _mm256_and_si256(_mm256_set1_epi32(*(uint32_t*)(gray_image + i)), MASK);
                
                // Summing together blocks of 8 uint8_t's and storing results
                // uint32_t eeee = a
                // uint32_t ffff = a + b
                // uint32_t gggg = a + b + c
                // uint32_t hhhh = a + b + c + d
                // eeee0000 ffff0000 gggg0000 hhhh0000
                __m256i img_sad_sum = _mm256_sad_epu8(img_rep_mask, ZERO);

                // Updating row_sum by adding previously computed img_sad_sum
                //   eeee0000 ffff0000 hhhh0000 hhhh0000
                // + iiii**** jjjj**** kkkk**** llll****
                // -------------------------------------
                //   mmmm**** nnnn**** oooo**** pppp****
                row_sums = _mm256_add_epi32(row_sums, img_sad_sum);

                // Permuting row_sum from first row to second and cutting off upper 128 bits
                // mmmm**** nnnn**** oooo**** pppp****
                // mmmmnnnn oooopppp ******** ********
                __m128i curr_row_sums = _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(row_sums, PACK));
                
                // Storing 128 bits (4*uint32_t) in curr_row_sums back to memory
                // 
                _mm_storeu_si128((__m128i*)(iimage_idata + i), curr_row_sums); // For unaligned
                //_mm_store_si128((__m128i *)(iimage_idata + i), curr_row_sums);
                
                // Updating row_sums to only include highest of four values (pppp)
                // mmmm**** nnnn**** oooo**** pppp****
                // pppp**** pppp**** pppp**** pppp****
                row_sums = _mm256_permute4x64_epi64(row_sums, 0xFF);

            }
            uint32_t row_sum = iimage_idata[i - 1];
            for (; i < width; ++i) {
                row_sum += gray_image[i];
                iimage_idata[i] = row_sum;
            }
        }

        float *iimage_data = iimage->data;
        __m128 factor = _mm_set1_ps(1.0f / 255.0f);

        for(size_t j = 1; j < height; ++j) {

            __m256i row_sums = ZERO;

            size_t i = 0;
            for(; i < aligned_width; i+=4) {

                // If (gray_image + i) has following uint8_t pixels: a b c d
                // a0000000 ab000000 abc00000 abcd0000
                __m256i img_rep_mask = _mm256_and_si256(_mm256_set1_epi32(*(uint32_t*)(gray_image + j * width + i)), MASK);
                
                // Summing together blocks of 8 uint8_t's and storing results
                // uint32_t eeee = a
                // uint32_t ffff = a + b
                // uint32_t gggg = a + b + c
                // uint32_t hhhh = a + b + c + d
                // eeee0000 ffff0000 gggg0000 hhhh0000
                __m256i img_sad_sum = _mm256_sad_epu8(img_rep_mask, ZERO);

                // Updating row_sum by adding previously computed img_sad_sum
                //   eeee0000 ffff0000 hhhh0000 hhhh0000
                // + iiii**** jjjj**** kkkk**** llll****
                // -------------------------------------
                //   mmmm**** nnnn**** oooo**** pppp****
                row_sums = _mm256_add_epi32(row_sums, img_sad_sum);

                // Permuting row_sum from first row to second and cutting off upper 128 bits
                // mmmm**** nnnn**** oooo**** pppp****
                // mmmmnnnn oooopppp ******** ********
                __m128i curr_row_sums = _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(row_sums, PACK));

                // Loading row_sum from previous row
                __m128i prev_row_sums = _mm_loadu_si128((__m128i *)(iimage_idata + (j-1) * data_width + i));
                //__m128i prev_row_sums = _mm_load_si128((__m128i *)(iimage_idata + (j-1) * data_width + i));

                // AVX512 only :(
                //__m128 prev_row_sums_cast = _mm_cvtepu32_ps(prev_row_sums);
                __m128 prev_row_sums_cast = _mm_cvtepi32_ps(prev_row_sums);
                //__128 prev_row_sums_cast = _my_mm_cvtepu32_ps(prev_row_sums);

                // Multiplying float version of previous row sum converted with (1.0f/255.0f) 
                __m128 prev_iimage_float = _mm_mul_ps(factor, prev_row_sums_cast);

                // Storing float of previous row sum back to memory as final integral image value
                _mm_storeu_ps(iimage_data +(j-1) * data_width + i, prev_iimage_float);

                // Adding curr_row_sums with prev_row_sums
                __m128i add_row_sums = _mm_add_epi32(curr_row_sums, prev_row_sums);

                // Storing 128 bits (4*uint32_t) in curr_row_sums back to memory
                _mm_storeu_si128((__m128i *)(iimage_idata + j * data_width + i), add_row_sums); // For unaligned
                //_mm_store_si128((__m128i *)(iimage_idata + j * data_width + i), add_row_sums);
                
                // Updating row_sums to only include highest of four values (pppp)
                // mmmm**** nnnn**** oooo**** pppp****
                // pppp**** pppp**** pppp**** pppp****
                row_sums = _mm256_permute4x64_epi64(row_sums, 0xFF);

            }

            uint32_t row_sum = iimage_idata[j * data_width + i - 1] - iimage_idata[(j-1) * data_width + i - 1];
            for (; i < width; ++i) {
                row_sum += gray_image[j * width + i];
                float prev_row_sum = iimage_idata[(j-1) * data_width + i];
                iimage_idata[j * data_width + i] = row_sum + prev_row_sum;
                iimage_data[(j-1) * data_width + i] = (1.0f / 255.0f) * ((float) prev_row_sum);
            }

        }

        // Converting last row to floats
        __m256 factor256 = _mm256_set1_ps(1.0f / 255.0f);
        size_t i = 0;
        for (; i < width - 7; i+=8) { // TODO: Check if -7 correct here

            __m256i v_int = _mm256_loadu_si256((__m256i *)(iimage_idata + (height - 1) * data_width + i));
            //__m256i v_int = _mm256_load_si256((__m256i *)(iimage_idata + j * data_width + i));

            // AVX512 only :(
            //__m256 v_float = _mm256_cvtepu32_ps(v_int);
            __m256 v_float = _mm256_cvtepi32_ps(v_int);
            //__m256 v_float = _my_mm256_cvtepu32_ps(v_int);

            // Multiplying converted float with (1.0f/255.0f) 
            __m256 v_iimage = _mm256_mul_ps(factor256, v_float);

            // Storing value back to memory
            _mm256_storeu_ps(iimage_data + (height - 1) * data_width + i, v_iimage);

        } 

        for (; i < width; ++i) {
            iimage_data[(height - 1) * data_width + i] = (1.0f / 255.0f) * ((float) iimage_idata[(height - 1) * data_width + i]);
        }

    }

    // BLOCK C
    for (int j = border; j < border + height; ++j) {
        float last_element_in_row_float = iimage_padded_data[j * data_width + border + width - 1];
        __m256 last_element_in_row_vec = _mm256_set1_ps(last_element_in_row_float);
        int i = width + border;
        for (; i < data_width - 7; i+=8) {
            _mm256_storeu_ps(iimage_padded_data + j * data_width + i, last_element_in_row_vec);
        }
        for (; i < data_width; ++i) {
            iimage_padded_data[j * data_width + i] = last_element_in_row_float;
        }
    }

    // BLOCK D   
    for (int j = border + height; j < data_height; ++j) {
        int i = border;
        for (; i < width + border - 7; i+=8) {
            __m256 prev_elements_in_row = _mm256_loadu_ps(iimage_padded_data + (j-1) * data_width + i);
            _mm256_storeu_ps(iimage_padded_data + j * data_width + i, prev_elements_in_row);
        }
        for (; i < width + border; ++i) {
            iimage_padded_data[j * data_width + i] = iimage_padded_data[(j-1) * data_width + i];
        }
    }
 
// TODO: WRITE SECOND VERSION WITH COLMAJOR
/*  
    int i = border;
    for (; i < width + border - 7; i+=8) {
        __m256 prev_elements_in_row_vec = _mm256_loadu_ps(iimage_padded_data + (border + height - 1) * data_width + i);
        for (int j = border + height; j < data_height; ++j) {
            _mm256_storeu_ps(iimage_padded_data + j * data_width + i, prev_elements_in_row_vec);
        }
    }
    for (; i < width + border; ++i) {
        float prev_elements_in_row_float = iimage_padded_data[(border + height - 1) * data_width + i];
        for (int j = border + height; j < data_height; ++j) {
            iimage_padded_data[j * data_width + i] = prev_elements_in_row_float;
        }
    }
*/    

    // BLOCK E
    int index_last_element = (border + height - 1) * data_width + border + width - 1;
    float max_value_float = iimage_padded_data[index_last_element];
    __m256 max_value_vec = _mm256_set1_ps(max_value_float);
    for (int j = border + height; j < data_height; ++j) {
        int i = border + width;
        for (; i < data_width - 7; i+=8) {
            _mm256_storeu_ps(iimage_padded_data + j * data_width + i, max_value_vec);
        }
        for (; i < data_width; ++i) {
            iimage_padded_data[j * data_width + i] = max_value_float;
        }
    }
    
}
