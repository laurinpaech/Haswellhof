#pragma once

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#define MIN(a, b)               \
    ({                          \
        __typeof__(a) _a = (a); \
        __typeof__(b) _b = (b); \
        _a < _b ? _a : _b;      \
    })

struct integral_image {
    int width;
    int height;
    /* https://aticleworld.com/pointer-inside-a-structure/ */
    float *data;
};

// Creates the struct of the integral image with empty data
struct integral_image *create_integral_img(int width, int height);

// Computes the integral image
void compute_integral_img(float *gray_image, int width, int height, float *iimage_data);

inline float box_integral(struct integral_image *iimage, int row, int col, int rows, int cols) {
    float *data = (float *)iimage->data;
    int width = iimage->width;
    int height = iimage->height;

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

inline float box_integral_debug(struct integral_image *iimage, int row, int col, int rows, int cols, int print) {
    float *data = (float *)iimage->data;
    int width = iimage->width;
    int height = iimage->height;

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
    if(print == 1){
                printf(" row: %i, col: %i\n", row, col);

    printf("ORIGINAL: r0: %i, c0: %i, r1: %i, c1: %i, A: %f, B: %f, C: %f, D: %f\n", r0, c0, r1, c1, A, B, C, D);
    }
    return fmax(0.0f, A - B - C + D);
}
