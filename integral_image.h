#pragma once

#include <stdlib.h>
#include <math.h>

struct integral_image {
    int width;
    int height;
    float *data;
};

float box_integral(struct integral_image *iimage, int row, int col, int rows, int cols) {

    float *data = (float *) iimage->data;
    int step = 1;

    // subtracting by one for row/col because row/col is inclusive.
    int r0 = fmin(row, iimage->height) - 1;
    int c0 = fmin(col, iimage->width) - 1;
    int r1 = fmin(row + rows, iimage->height) - 1;
    int c1 = fmin(col + cols, iimage->width) - 1;

    float A = 0.0f;
    float B = 0.0f; 
    float C = 0.0f; 
    float D = 0.0f;

    if (r0 >= 0 && c0 >= 0) {
        A = data[r0 * step + c0];
    }
    if (r0 >= 0 && c1 >= 0) {
        B = data[r0 * step + c1];
    }  
    if (r1 >= 0 && c0 >= 0) {
        C = data[r1 * step + c0];
    }
    if (r1 >= 0 && c1 >= 0) {
        D = data[r1 * step + c1];
    }

    return fmax(0.0f, A - B - C + D);

}
