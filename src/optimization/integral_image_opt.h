#pragma once

#include <stdio.h>

#include "integral_image.h"

// Returns the number of flops of the optimized function compute_integral_img_faster_alg.
inline long get_flops_compute_integral_img_faster_alg(int width, int height, int stride) {
    // First 2 rows = 3*width -> odd width + 3
    // Remaining rows (height - stride) * width + 2 -> odd width + 2*(height-2)
    // Odd height + 2*width and odd width + 2
    long num_flops = 3 * width + (width % stride * 3) + (height - stride) * width * 2 +
                     (width % stride) * (2 * (height - stride)) +
                     (height % stride) * (2 * width + (width % stride) * 2);
    return num_flops;
}

// An optimized function to compute the integral image.
// Parallelizes the additions which makes use of both addition ports.
// Computes two rows simultaneously.
void compute_integral_img_faster_alg(float *gray_image, struct integral_image *iimage);

// Creates the struct of the padded integral image.
// Has a larger size for the padding. the width is image_width + 2 * largest_lobe. The height is
// image_height + 2 * largest_border.
struct integral_image *create_padded_integral_img(int width, int height);

void compute_padded_integral_img(float *gray_image, struct integral_image *iimage);

inline float box_integral_improved(struct integral_image *iimage, int row, int col, int rows, int cols) {
    float *data = (float *)iimage->data;
    int data_width = iimage->data_width;
    int width = iimage->width;
    int height = iimage->height;
    float temp0, temp1, res;

    // subtracting by one for row/col because row/col is inclusive.
    int r0 = MIN(row, height) - 1;         // r - 3
    int c0 = MIN(col, width) - 1;          // c - b - 1
    int r1 = MIN(row + rows, height) - 1;  // r - 3 + 5
    int c1 = MIN(col + cols, width) - 1;   // c - b + filter_size - 1

    // Example with 9x9 filter at (0,0)
    // A: (-3, -5)
    // B: (-3, 8)
    // C: (2, -1)
    // D: (2, 4)

    float A = 0.0f;
    float B = 0.0f;
    float C = 0.0f;
    float D = 0.0f;

    if (r0 >= 0 && c0 >= 0) {
        A = data[r0 * data_width + c0];
    }
    if (r0 >= 0 && c1 >= 0) {
        B = data[r0 * data_width + c1];
    }
    if (r1 >= 0 && c0 >= 0) {
        C = data[r1 * data_width + c0];
    }
    if (r1 >= 0 && c1 >= 0) {
        D = data[r1 * data_width + c1];
    }

    // there was a floating point arithmetic bug in the original implementation
    // this fixes it and now fmaxf is not needed
    // use this for validation:
    temp0 = A - C;
    temp1 = D - B;
    res = temp0 + temp1;
    
    return res;
}

inline float box_integral_unconditional(struct integral_image *iimage, int row, int col, int rows, int cols) {

    float *data = iimage->data;
    int data_width = iimage->data_width;
    float temp0, temp1, res;

    // subtracting by one for row/col because row/col is inclusive.
    int r0 = row - 1;
    int c0 = col - 1;
    int r1 = row + rows - 1;
    int c1 = col + cols - 1;

    // r0 = x - lobe/2 -1
    // c0 = y - lobe
    // r1 = x + lobe/2 -1
    // c1 = y + lobe - 1

    float A = data[r0 * data_width + c0];
    float B = data[r0 * data_width + c1];
    float C = data[r1 * data_width + c0];
    float D = data[r1 * data_width + c1];

    // there was a floating point arithmetic bug in the original implementation
    // this fixes it and now fmaxf is not needed
    temp0 = A - C;
    temp1 = D - B;
    res = temp0 + temp1;

    // fmaxf instead of fmax (faster?)
    // return fmaxf(0.0f, A - B - C + D);
    return res;
}

inline float box_integral_unconditional_opt(struct integral_image *iimage, int r0, int c0, int r1, int c1) {

    float *data = iimage->data;
    int data_width = iimage->data_width;
    float temp0, temp1, res;

    float A = data[r0 * data_width + c0];
    float B = data[r0 * data_width + c1];
    float C = data[r1 * data_width + c0];
    float D = data[r1 * data_width + c1];

    // there was a floating point arithmetic bug in the original implementation
    // this fixes it and now fmaxf is not needed
    temp0 = A - C;
    temp1 = D - B;
    res = temp0 + temp1;

    // fmaxf instead of fmax (faster?)
    // return fmaxf(0.0f, A - B - C + D);
    return res;
}
