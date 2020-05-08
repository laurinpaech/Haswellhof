#pragma once

#include "integral_image.h"
#include <stdio.h>


// Returns the number of flops of the optimized function compute_integral_img_faster_alg.
static inline long get_flops_compute_integral_img_faster_alg( int width, int height, int stride){

    // First 2 rows = 3*width -> odd width + 3
    // Remaining rows (height - stride) * width + 2 -> odd width + 2*(height-2)
    // Odd height + 2*width and odd width + 2
    long num_flops = 3*width + (width%stride*3) + (height - stride)*width*2 + (width%stride)*(2*(height-stride)) + (height%stride) *(2*width+ (width%stride)*2);
    return num_flops;
}

// An optimized function to compute the integral image.
// Parallelizes the additions which makes use of both addition ports.
// Computes two rows simultaneously.
void compute_integral_img_faster_alg(float *gray_image, int width, int height, float *iimage_data);