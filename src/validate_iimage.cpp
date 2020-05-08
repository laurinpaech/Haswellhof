#include "validate_iimage.h"

#include <stdio.h>

bool evaluate_iimage(void (*original_function)(float *, int, int, float *),
                     void (*optimized_function)(float *, int, int, float *), int width, int height, float *image) {
    // Create integral image
    struct integral_image *original_iimage = create_integral_img(width, height);
    // Compute integral image
    original_function(image, original_iimage->width, original_iimage->height, original_iimage->data);

    // Create integral image
    struct integral_image *optimized_iimage = create_integral_img(width, height);
    // Compute integral image
    optimized_function(image, optimized_iimage->width, optimized_iimage->height, optimized_iimage->data);

    bool equal = are_equal(original_iimage->data, optimized_iimage->data, width, height);

    free(original_iimage->data);
    free(original_iimage);
    free(optimized_iimage->data);
    free(optimized_iimage);

    return equal;
}

bool evaluate_iimage_custom_array(void (*original_function)(float *, int, int, float *),
                                  void (*optimized_function)(float *, int, int, float *)) {
    int width = 4, height = 4;
    float *image = (float *)malloc(height * width * sizeof(float));
    int counter = 0;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            image[i * width + j] = counter;
            counter++;
        }
    }

    //float image[height][width] = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 }, { 7.0, 8.0, 9.0 } };


    printf("create integral image\n");
    // Create integral image
    struct integral_image *original_iimage = create_integral_img(width, height);
    // Compute integral image
    original_function(image, original_iimage->width, original_iimage->height, original_iimage->data);

    // Create integral image
    struct integral_image *optimized_iimage = create_integral_img(width, height);
    // Compute integral image
    optimized_function(image, optimized_iimage->width, optimized_iimage->height, optimized_iimage->data);

    print_debug(original_iimage->data, optimized_iimage->data, width, height);
    printf("before equal\n");
    bool equal = are_equal(original_iimage->data, optimized_iimage->data, width, height);
    printf("after equal\n");

    free(original_iimage->data);
    free(original_iimage);
    free(optimized_iimage->data);
    free(optimized_iimage);
    free(image);
    printf("after free");
    return equal;
}

bool are_equal(float *iimage1, float *iimage2, int width, int height) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (iimage1[i * width + j] != iimage2[i * width + j]) {
                return false;
            }
        }
    }
    return true;
}

void print_debug(float *iimage1, float *iimage2, int width, int height) {
    printf("print debug\n");
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            printf("i1: %f, i2: %f\n", iimage1[i * width + j], iimage2[i * width + j]);
        }
    }
}