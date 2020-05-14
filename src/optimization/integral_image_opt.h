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
void compute_integral_img_faster_alg(float *gray_image, int width, int height, float *iimage_data);

// Creates the struct of the padded integral image.
// Has a larger size for the padding. the width is image_width + 2 * largest_lobe. The height is
// image_height + 2 * largest_border.
struct integral_image *create_padded_integral_img(int width, int height);

void compute_padded_integral_image(float *gray_image, int width, int height, float *iimage_data);

float box_integral_with_padding(struct integral_image *iimage, int row, int col, int rows, int cols, int print);