#include "integral_image.h"
#include <stdio.h>


struct integral_image* create_integral_img(float* gray_image, unsigned width, unsigned height){

    struct integral_image* integral_img = (struct integral_image *)malloc(sizeof(struct integral_image));
    integral_img-> height = height;
    integral_img-> width = width;
    float *data = (float*)malloc(width * height * sizeof(float));
   
    float row_sum = 0.0f;

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
        }
    }

    integral_img->data = data;

    return integral_img;
        
}

float box_integral(struct integral_image *iimage, int row, int col, int rows, int cols) {

    float *data = (float *) iimage->data;
    int width = iimage->width;
    int height = iimage->height;

    // subtracting by one for row/col because row/col is inclusive.
    int r0 = fmin(row, height) - 1;
    int c0 = fmin(col, width) - 1;
    int r1 = fmin(row + rows, height) - 1;
    int c1 = fmin(col + cols, width) - 1;

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
