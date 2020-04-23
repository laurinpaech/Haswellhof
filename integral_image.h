#pragma once

#include <stdlib.h>
#include <math.h>


struct integral_image {
    int width;
    int height;
    /* https://aticleworld.com/pointer-inside-a-structure/ */
    float *data;
};

struct integral_image* create_integral_img(float* gray_image, unsigned width, unsigned height);


float box_integral(struct integral_image *iimage, int row, int col, int rows, int cols);
