#include "playground.h"

// NOTE: Fill this with whatever you want (preferably don't commit/push it though)

#include <iostream>

#include "benchmarking.h"
#include "descriptor.h"
#include "helper.h"
#include "integral_image.h"
#include "integral_image_opt.h"
#include "validation_integral_image.h"

void playground_function1(float *image, int width, int height) {
    // Create integral image
    struct integral_image *optimized_iimage = create_padded_integral_img(width, height);
    // Compute integral image

    compute_padded_integral_image(image, width, height, optimized_iimage->data);
    printf("PADDED IMAGE\n");
/**
    for (size_t i = 0; i < optimized_iimage->height; i++) {
        printf("row: %i \n", i);
        for (size_t j = 0; j < optimized_iimage->width; j++) {
            printf("%f ", optimized_iimage->data[i * optimized_iimage->width + j]);
        }
        printf("\n\n");
    }
*/
    free(optimized_iimage->data);
    free(optimized_iimage);
}

void playground_function2() {
    int width = 4, height = 4;

    // Create integral image
    struct integral_image *optimized_iimage = create_padded_integral_img(width, height);
    // Compute integral image

    float *image = (float *)malloc(height * width * sizeof(float));
    int counter = 0;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            image[i * width + j] = counter;
            counter++;
        }
    }

    compute_padded_integral_image(image, width, height, optimized_iimage->data);
    printf("PADDED IMAGE\n");

    for (size_t i = 0; i < optimized_iimage->height; i++) {
        printf("row: %i \n", i);
        for (size_t j = 0; j < optimized_iimage->width; j++) {
            printf("%f ", optimized_iimage->data[i * optimized_iimage->width + j]);
        }
        printf("\n\n");
    }

    free(optimized_iimage->data);
    free(optimized_iimage);
    free(image);
}
