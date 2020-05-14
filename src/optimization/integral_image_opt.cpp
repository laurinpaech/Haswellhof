#include "integral_image_opt.h"

#include <stdio.h>

#include "fasthessian.h"

// An optimized function to compute the integral image.
// Parallelizes the additions which makes use of both addition ports.
// Computes two rows simultaneously.
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
        // printf("row 0, col %i sum: %f\n", i, sum);
        sum1 += gray_image[width + i];
        // printf("row 1, col %i sum: %f\n", i, sum1);
        iimage_data[width + i] = iimage_data[i] + sum1;
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
        // printf("j: %i\n", j);

        // Handle last element of an image with odd width extra.
        // If previous stride = 2, then there should only be one element covered in this loop.
        // This might change if we increase the stride.
        for (; j < width; ++j) {
            // printf("j: %i; sum1: %f\n", j, sum1);
            sum += gray_image[i * width + j];
            iimage_data[i * width + j] = iimage_data[(i - 1) * width + j] + sum;

            sum1 += gray_image[(i + 1) * width + j];
            // printf("sum1: %f\n", sum1);
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
            // printf("j: %i; sum1: %f\n", j, sum);

            sum += gray_image[i * width + j];
            // printf("j: %i; sum1: %f\n", j, sum);
            iimage_data[i * width + j] = iimage_data[(i - 1) * width + j] + sum;
            sum += gray_image[i * width + j + 1];
            // printf("j: %i; sum1: %f\n", j, sum);

            iimage_data[i * width + j + 1] = iimage_data[(i - 1) * width + j + 1] + sum;
        }
        // printf("j: %i\n", j);

        for (; j < width; ++j) {
            sum += gray_image[i * width + j];
            // printf("j: %i; sum1: %f\n", j, sum);

            iimage_data[i * width + j] = iimage_data[(i - 1) * width + j] + sum;
        }
    }
}

// Creates the struct of the padded integral image
// Has a larger size for the padding. the width is image_width + 2 * largest_lobe. The height is
// image_height + 2 * largest_border.
struct integral_image *create_padded_integral_img(int width, int height) {
    struct integral_image *iimage = (struct integral_image *)malloc(sizeof(struct integral_image));
    iimage->height = height;
    iimage->width = width;

    // Add the border above the image all values will be 0.
    // Add the border below the image, all values will be equivalent to the last row.

    // Add the border to the left side of the image, all values will be 0.
    // Add the border to the right side of the image, all values will be equivalent to the last column.

    // The right lower corner will contain only the max value of the integral image. (lowest right-most value)

    // Border + 1 because A, B and C are always exclusive.
    int border = ((LARGEST_FILTER_SIZE - 1) / 2) + 1;

    // Border as padding above and below the image because of Dyy.
    int padded_height = height + border * 2;
    // Border as padding to the left and right because of Dxx.
    int padded_width = width + border * 2;

    iimage->data = (float *)malloc(padded_width * padded_height * sizeof(float));

    return iimage;
}

void compute_padded_integral_image(float *gray_image, int original_image_width, int original_image_height,
                                   float *iimage_data) {
    // The border and the lobe must be the same, since the filters are being turned dependant on Dxx/Dyy. The border is
    // always larger than the lobe. We have to make sure that the border is available as padding in all directions.
    // Border + 1 because A, B and C are always exclusive. We would need too many special cases if we'd work with
    // different upper and lower borders and left and right borders.
    int border = ((LARGEST_FILTER_SIZE - 1) / 2) + 1;

    int padded_height = original_image_height + border * 2;

    // The width must also pad for the border, because of Dxx.
    int padded_width = original_image_width + border * 2;

    // Pad the top part with 0.
    for (int i = 0; i < border; i++) {
        for (int j = 0; j < padded_width; j++) {
            iimage_data[i * padded_width + j] = 0.0f;
        }
    }

    // Pad the remaining left part with 0.
    for (int i = border; i < padded_height; i++) {
        for (int j = 0; j < border; j++) {
            iimage_data[i * padded_width + j] = 0.0f;
        }
    }

    float row_sum = 0.0f;
    float last_element_in_row = 0.0f;
    // Sum up the first row
    // The first element that is not padding.
    int ind = (border * padded_width) + border;

    for (int i = 0; i < original_image_width; i++, ind++) {
        // Previous rows are 0
        row_sum += gray_image[i];
        iimage_data[ind] = row_sum;
        last_element_in_row = row_sum;
    }
    // Pad the right side with the last element of the row.
    for (int j = original_image_width + border; j < padded_width; j++, ind++) {
        iimage_data[ind] = last_element_in_row;
    }

    // Sum all remaining rows of the original image
    for (int i = 1; i < original_image_height; ++i) {
        row_sum = 0.0f;
        ind += border;

        for (int j = 0; j < original_image_width; ++j, ind++) {
            row_sum += gray_image[i * original_image_width + j];
            // Add sum of current row until current idx to sum of all previous rows until current index
            last_element_in_row = row_sum + iimage_data[ind - padded_width];
            iimage_data[ind] = last_element_in_row;
        }

        // Pad the right side with the last element of each row.
        for (int j = original_image_width + border; j < padded_width; j++, ind++) {
            iimage_data[ind] = last_element_in_row;
        }
    }

    // Pad the middle lower part of the integral image with the elements of the last row.
    for (int i = original_image_height + border; i < padded_height; i++) {
        for (int j = border; j < padded_width - border; j++) {
            iimage_data[i * padded_width + j] = iimage_data[(i - 1) * padded_width + j];
        }
    }

    // Pad the right lower corner of the integral image with the max value of the integral image.
    int index_last_element = (original_image_height + border - 1) * padded_width + original_image_width + border;
    float max_value = iimage_data[index_last_element];

    for (int i = border + original_image_height; i < padded_height; i++) {
        for (int j = border + original_image_width; j < padded_width; j++) {
            iimage_data[i * padded_width + j] = max_value;
        }
    }
}

float box_integral_with_padding(struct integral_image *iimage, int row, int col, int rows, int cols, int print) {
    float *data = (float *)iimage->data;
    int border = (LARGEST_FILTER_SIZE - 1) / 2 + 1;

    // The padded width must be + border because of Dxx. If we'd only account for Dyy lobe would be sufficient.
    int padded_width = iimage->width + (border * 2);

    // subtracting by one for row/col because row/col is inclusive.
    int r0 = row - 1;
    int c0 = col - 1;
    int r1 = row + rows - 1;
    int c1 = col + cols - 1;

    float A = data[r0 * padded_width + c0];
    float B = data[r0 * padded_width + c1];
    float C = data[r1 * padded_width + c0];
    float D = data[r1 * padded_width + c1];

    return fmax(0.0f, A - B - C + D);
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
