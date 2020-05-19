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
    
    // Width and height of (original) image
    int width;
    int height;

    // Pointer to upper left corner (origin) of integral image
    float *data;

    // Width and height of (potentially) padded image (use this for indexing with data) 
    int data_width;
    int data_height;

    // Pointer to upper left corner (origin) of padding of integral image (same as data if no padding applies)
    float *padded_data;

};

// Creates the struct of the integral image with empty data
struct integral_image *create_integral_img(int width, int height);

// Computes the integral image
void compute_integral_img(float *gray_image, struct integral_image *iimage);

inline float box_integral(struct integral_image *iimage, int row, int col, int rows, int cols) {
    float *data = (float *)iimage->data;
    int data_width = iimage->data_width;
    int width = iimage->width;
    int height = iimage->height;
    float res, temp0, temp1;

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
    temp0 = A - C;
    temp1 = D - B;
    res = temp0 + temp1;

    // fmaxf instead of fmax (faster?)
    // return fmaxf(0.0f, A - B - C + D);
    return res;
}
