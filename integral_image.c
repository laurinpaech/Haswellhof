#include "integral_image.h"
#include <stdio.h>


struct integral_image* create_integral_img(float* gray_image, int width, int height){

    struct integral_image* integral_img = (struct integral_image *)malloc(sizeof(struct integral_image));
    integral_img->height = height;
    integral_img->width = width;
    float *data = (float*) malloc(width * height * sizeof(float));

    float row_sum = 0.0f;

    // TESTING
    // for (int i = 0; i < width; i++) {
    //     for (int j = 0; j < height; j++) {
    //         printf("%i, %i - %f\n", i, j, gray_image[i*width+j]);
    //     }
    // }

    /* sum up the first row */

    for(int i=0; i<width; i++)
    {
        /* previous rows are 0 */
        row_sum += gray_image[i];
        data[i] = row_sum;
    }

    /* sum all remaining rows*/
    for(int i=1; i<height; ++i)
    {
        row_sum = 0.0f;
        for(int j=0; j<width; ++j)
        {
            row_sum += gray_image[i*width+j];
            /*add sum of current row until current idx to sum of all previous rows until current index */
            data[i*width+j] = row_sum + data[(i-1)*width+j];

            // TESTING
            // printf("%i, %i - %f\n", i, j, data[i*width+j]);
            // TESTING
        }
    }

    integral_img->data = data;

    return integral_img;

}
