#include <stdio.h>

#include "integral_image_opt.h"


// An optimized function to compute the integral image.
// Parallelizes the additions which makes use of both addition ports.
void compute_integral_img_faster_alg(float *gray_image, int width, int height, float *iimage_data) {
    float sum = 0.0f;
    float sum1 = 0.0f;
    int width_limit = width - 2;
    int height_limit = height - 2;

    int i = 0;
    // first row extra, since we don't have 0 padding
    for (i = 0; i < width_limit; i += 2) {
        sum += gray_image[i];
        iimage_data[i] = sum;
        sum += gray_image[i + 1];
        iimage_data[i + 1] = sum;

        sum1 += gray_image[width + i];
        iimage_data[width + i] = iimage_data[i] + sum1;
        sum1 += gray_image[width + i + 1];
        iimage_data[width + i + 1] = iimage_data[i + 1] + sum1;
    }

    // Handle last element of an image with odd width extra. 
    // If previous stride = 2, then there should only be one element covered in this loop.
    // This might change if we increase the stride.
    for (; i < width; ++i) {
        sum += gray_image[i];
        iimage_data[i] = sum;
        //printf("row 0, col %i sum: %f\n", i, sum);
        sum1 += gray_image[width + i];
        //printf("row 1, col %i sum: %f\n", i, sum1);
        iimage_data[width + i] = iimage_data[i] +  sum1;
    }

    // Handle all elements after the first row.
    i = 2;
    for (; i < height_limit; i += 2) {
        sum = 0.0f;
        sum1 = 0.0f;
        int j = 0;
        for (; j < width_limit; j += 2) {
            sum += gray_image[i * width + j];
            iimage_data[i * width + j] = iimage_data[(i - 1) * width + j] + sum;
            sum += gray_image[i * width + j + 1];
            iimage_data[i * width + j + 1] = iimage_data[(i - 1) * width + j + 1] + sum;

            sum1 += gray_image[(i + 1) * width + j];
            iimage_data[(i + 1) * width + j] = iimage_data[i * width + j] + sum1;
            sum1 += gray_image[(i + 1) * width + j + 1];
            iimage_data[(i + 1) * width + j + 1] = iimage_data[i * width + j + 1] + sum1;
        }
        //printf("j: %i\n", j);

        // Handle last element of an image with odd width extra. 
        // If previous stride = 2, then there should only be one element covered in this loop.
        // This might change if we increase the stride.
        for (; j < width; ++j) {
            //printf("j: %i; sum1: %f\n", j, sum1);
            sum += gray_image[i * width + j];
            iimage_data[i * width + j] = iimage_data[(i - 1) * width + j] + sum;

            sum1 += gray_image[(i + 1) * width + j];
            //printf("sum1: %f\n", sum1);
            iimage_data[(i + 1) * width + j] = iimage_data[i * width + j] + sum1;
        }
    }

    // Handle last element of an image with odd height extra. 
    // If previous stride = 2, then there should only be one row covered in this loop.
    // This might change if we increase the stride.
    for (; i < height; ++i) {
        sum = 0.0f;
        int j = 0;
        for (; j < width_limit; j += 2) {
            //printf("j: %i; sum1: %f\n", j, sum);

            sum += gray_image[i * width + j];
           // printf("j: %i; sum1: %f\n", j, sum);
           iimage_data[i * width + j] = iimage_data[(i - 1) * width + j] + sum;
            sum += gray_image[i * width + j + 1];
            //printf("j: %i; sum1: %f\n", j, sum);

            iimage_data[i * width + j + 1] = iimage_data[(i - 1) * width + j + 1] + sum;
        }
        //printf("j: %i\n", j);

        for (; j < width; ++j) {
            sum += gray_image[i * width + j];
            //printf("j: %i; sum1: %f\n", j, sum);

            iimage_data[i * width + j] = iimage_data[(i - 1) * width + j] + sum;
        }
    }
}

// Computes the integral image
/**
void compute_integral_img_opt(float *gray_image, int width, int height, float *iimage_data) {
    float row1_col1 = 0.0f, row1_col2 = 0.0f, row1_col3 = 0.0f, row1_col4 = 0.0f;
    float row2_col1 = 0.0f, row2_col2 = 0.0f, row2_col3 = 0.0f, row2_col4 = 0.0f;
    float row3_col1 = 0.0f, row3_col2 = 0.0f, row3_col3 = 0.0f, row3_col4 = 0.0f;
    float row4_col1 = 0.0f, row4_col2 = 0.0f, row4_col3 = 0.0f, row4_col4 = 0.0f;

    float prev_row_sum1 = 0.0f, prev_row_sum2 = 0.0f, prev_row_sum3 = 0.0f, prev_row_sum4 = 0.0f;


    // compute first block
for (int j = 0; j < width; j+=4) {

        row1_col1 = prev_row_sum1 + gray_image[j];
        row1_col2 = row1_col1 + gray_image[j+1];
        row1_col3 = row1_col2 + gray_image[j+2];
        row1_col4 = row1_col3 + gray_image[j+3];


        row2_col1 = prev_row_sum2 + row1_col1 + gray_image[1*width+j];
        row2_col2 = row2_col1 + gray_image[(i+1)*width+j+1];
        row2_col3 = row2_col2 + gray_image[(i+1)*width+j+2];
        row2_col4 = row2_col3 + gray_image[(i+1)*width+j+3];


        row3_col1 = prev_row_sum3 + row2_col1 + gray_image[(i+2)*width+j];
        row3_col2 = row3_col1 + gray_image[(i+2)*width+j+1];
        row3_col3 = row3_col2 + gray_image[(i+2)*width+j+2];
        row3_col4 = row3_col3 + gray_image[(i+2)*width+j+3];


        row4_col1 = prev_row_sum4 + row3_col1 + gray_image[(i+3)*width+j];
        row4_col2 = row4_col1 + gray_image[(i+3)*width+j+1];
        row4_col3 = row4_col2 + gray_image[(i+3)*width+j+2];
        row4_col4 = row4_col3 + gray_image[(i+3)*width+j+3];

        prev_row_sum1 = prev_row_sum1 + row1_col4;
        prev_row_sum2 = prev_row_sum2 + row2_col4 - prev_row_sum1;
        prev_row_sum3 = prev_row_sum3 + row3_col4 - prev_row_sum2;
        prev_row_sum4 = prev_row_sum4 + row4_col4 - prev_row_sum3;

        iimage_data[i * width + j]

        }




    // sum up all rows
    for (int i = 0; i < height; i+=4) {
        for (int j = 0; j < width; j+=4) {

            row1_col1 = prev_row_sum1 + gray_image[i*width+j];
            row1_col2 = row1_col1 + gray_image[i*width+j+1];
            row1_col3 = row1_col2 + gray_image[i*width+j+2];
            row1_col4 = row1_col3 + gray_image[i*width+j+3];


            row2_col1 = prev_row_sum2 + row1_col1 + gray_image[(i+1)*width+j];
            row2_col2 = row2_col1 + gray_image[(i+1)*width+j+1];
            row2_col3 = row2_col2 + gray_image[(i+1)*width+j+2];
            row2_col4 = row2_col3 + gray_image[(i+1)*width+j+3];


            row3_col1 = prev_row_sum3 + row2_col1 + gray_image[(i+2)*width+j];
            row3_col2 = row3_col1 + gray_image[(i+2)*width+j+1];
            row3_col3 = row3_col2 + gray_image[(i+2)*width+j+2];
            row3_col4 = row3_col3 + gray_image[(i+2)*width+j+3];
penis

            row4_col1 = prev_row_sum4 + row3_col1 + gray_image[(i+3)*width+j];
            row4_col2 = row4_col1 + gray_image[(i+3)*width+j+1];
            row4_col3 = row4_col2 + gray_image[(i+3)*width+j+2];
            row4_col4 = row4_col3 + gray_image[(i+3)*width+j+3];

            prev_row_sum1 = prev_row_sum1 + row1_col4;
            prev_row_sum2 = prev_row_sum2 + row2_col4 - prev_row_sum1;
            prev_row_sum3 = prev_row_sum3 + row3_col4 - prev_row_sum2;
            prev_row_sum4 = prev_row_sum4 + row4_col4 - prev_row_sum3;

            iimage_data[i * width + j]

        }

        //TODO: handle remaining if row is not divisible by 4
    }
    */
