#include "integral_image_opt.h"

#include <stdio.h>

#include "fasthessian.h"

// An optimized function to compute the integral image.
// Parallelizes the additions which makes use of both addition ports.
// Computes two rows simultaneously.
void compute_integral_img_faster_alg(float *gray_image, struct integral_image * iimage) {
    
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

// Creates the struct of the padded integral image
// Has a larger size for the padding. the width is image_width + 2 * largest_border. The height is
// image_height + 2 * largest_border.
struct integral_image *create_padded_integral_img(int width, int height) {
    
    // Allocating memory to integral image that is returned
    struct integral_image *iimage = (struct integral_image *)malloc(sizeof(struct integral_image));

    // Border + 1 because A, B and C are always exclusive.
    int border = ((LARGEST_FILTER_SIZE - 1) / 2) + 1;

    // Getting data width  (width of padded image)
    // Border as padding to the left and right because of Dxx.
    int data_width = width + border * 2;

    // Getting data height (height of padded image)
    // Border as padding above and below the image because of Dyy.
    int data_height = height + border * 2;
    
    // Setting real image width and height
    iimage->width = width;
    iimage->height = height;
    
    // Setting data width and height (width and height of padded image)
    iimage->data_width = data_width;
    iimage->data_height = data_height;

    // Allocating data for storing values of integral image with padding
    iimage->padded_data = (float *)malloc(data_width * data_height * sizeof(float));
    
    // Setting data and padded data to be inner offset and origin of original image
    iimage->data = iimage->padded_data + (border * data_width) + border;

    return iimage;
    
}

void compute_padded_integral_img(float *gray_image, struct integral_image *iimage) {

    // TODO: (Sebastian) HACKY SOLUTION CLEANUP!

    // Getting width and height of original, unpadded image
    int width = iimage->width;
    int height = iimage->height;
    
    // Getting data width and data height, the width and height of the padded image
    int data_width = iimage->data_width;
    int data_height = iimage->data_height;

    // Getting pointer to upper left corner of padded image
    float *padded_data = iimage->padded_data;

    // The border and the lobe must be the same, since the filters are being turned dependant on Dxx/Dyy. The border is
    // always larger than the lobe. We have to make sure that the border is available as padding in all directions.
    // Border + 1 because A, B and C are always exclusive. We would need too many special cases if we'd work with
    // different upper and lower borders and left and right borders.
    int border = (data_width - width) / 2; 
    //int border = ((LARGEST_FILTER_SIZE - 1) / 2) + 1;

    // Pad the top part with 0.
    for (int i = 0; i < border; ++i) {
        for (int j = 0; j < data_width; ++j) {
            padded_data[i * data_width + j] = 0.0f;
        }
    }

    // Pad the remaining left part with 0.
    for (int i = border; i < data_height; ++i) {
        for (int j = 0; j < border; ++j) {
            padded_data[i * data_width + j] = 0.0f;
        }
    }

    float row_sum = 0.0f;
    float last_element_in_row = 0.0f;
    
    // Sum up the first row
    // The first element that is not padding.
    int ind = (border * data_width) + border;

    for (int i = 0; i < width; ++i, ++ind) {
        // Previous rows are 0
        row_sum += gray_image[i];
        padded_data[ind] = row_sum;
        last_element_in_row = row_sum;
    }
    // Pad the right side of the first row with the last element of the row.
    for (int j = width + border; j < data_width; ++j, ++ind) {
        padded_data[ind] = last_element_in_row;
    }

    // Sum all remaining rows of the original image
    for (int i = 1; i < height; ++i) {
        row_sum = 0.0f;
        ind += border;

        for (int j = 0; j < width; ++j, ++ind) {
            row_sum += gray_image[i * width + j];
            // Add sum of current row until current idx to sum of all previous rows until current index
            last_element_in_row = row_sum + padded_data[ind - data_width];
            padded_data[ind] = last_element_in_row;
        }

        // Pad the right side with the last element of each row.
        for (int j = width + border; j < data_width; ++j, ++ind) {
            padded_data[ind] = last_element_in_row;
        }
    }

    // Pad the middle lower part of the integral image with the elements of the last row.
    for (int i = height + border; i < data_height; ++i) {
        for (int j = border; j < data_width - border; ++j) {
            padded_data[i * data_width + j] = padded_data[(i - 1) * data_width + j];
        }
    }

    // Pad the right lower corner of the integral image with the max value of the integral image.
    int index_last_element = (height + border - 1) * data_width + width + border - 1;
    float max_value = padded_data[index_last_element];

    for (int i = border + height; i < data_height; ++i) {
        for (int j = border + width; j < data_width; ++j) {
            padded_data[i * data_width + j] = max_value;
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

        iimage_data[i * data_width + j]

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
