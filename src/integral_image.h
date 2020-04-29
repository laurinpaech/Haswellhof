#pragma once

#include <stdlib.h>
#include <math.h>

struct integral_image {
    int width;
    int height;
    /* https://aticleworld.com/pointer-inside-a-structure/ */
    float *data;
};

struct integral_image* create_integral_img(float* gray_image, int width, int height);

static inline float box_integral(struct integral_image *iimage, int row, int col, int rows, int cols) {

    float *data = (float *) iimage->data;
    int width = iimage->width;
    int height = iimage->height;

    // subtracting by one for row/col because row/col is inclusive.
    int r0 = fmin(row, height) - 1;         // r - 3
    int c0 = fmin(col, width) - 1;          // c - b - 1
    int r1 = fmin(row + rows, height) - 1;  // r - 3 + 5
    int c1 = fmin(col + cols, width) - 1;   // c - b + filter_size - 1

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
        A = data[r0 * width + c0];
    }
    if (r0 >= 0 && c1 >= 0) {
        B = data[r0 * width + c1];
    }
    if (r1 >= 0 && c0 >= 0) {
        C = data[r1 * width + c0];
    }
    if (r1 >= 0 && c1 >= 0) {
        D = data[r1 * width + c1];
    }

    return fmax(0.0f, A - B - C + D);

}