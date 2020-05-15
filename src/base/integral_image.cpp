#include "integral_image.h"

#include <stdio.h>

// Creates the struct of the integral image with empty data
struct integral_image *create_integral_img(int width, int height) {
    struct integral_image *iimage = (struct integral_image *)malloc(sizeof(struct integral_image));
    iimage->height = height;
    iimage->width = width;
    iimage->data = (float *)malloc(width * height * sizeof(float));

    return iimage;
}

// Computes the integral image
void compute_integral_img(float *gray_image, int width, int height, float *iimage_data) {
    float row_sum = 0.0f;

    int data_width = iimage->data_width;

    /* sum up the first row */

    for (int i = 0; i < width; i++) {
        /* previous rows are 0 */
        row_sum += gray_image[i];
        iimage_data[i] = row_sum;
    }

    /* sum all remaining rows*/
    for (int i = 1; i < height; ++i) {
        row_sum = 0.0f;
        for (int j = 0; j < width; ++j) {
            row_sum += gray_image[i * data_width + j];
            /*add sum of current row until current idx to sum of all previous rows until current index */
            iimage_data[i * data_width + j] = row_sum + iimage_data[(i - 1) * data_width + j];
        }
    }
}
