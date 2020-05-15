
#include "validation_integral_image.h"
#include "helper.h"

#include <stdio.h>


// Easier for debugging.
// Creates an integral image based on a custom array, the base function and a list of other functions. 
// The results of the integral images are being compared.
// If one of the results difers from the base implementation false is being returned.
// Messages clarifying the equality of the results are being printed.
bool validate_integral_image_custom_matrix(void (*original_function)(float *, int, int, float *),
                                  const std::vector<void (*)(float *, int, int, float *)> &test_functions) {
    int width = 4, height = 4;
    float *image = (float *)malloc(height * width * sizeof(float));
    int counter = 0;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            image[i * width + j] = counter;
            counter++;
        }
    }

    bool all_functions_equal = true;

    // Create integral image
    struct integral_image *original_iimage = create_integral_img(width, height);
    // Compute integral image
    original_function(image, original_iimage->width, original_iimage->height, original_iimage->data);
    for (auto optimized_function : test_functions) {
        // Create integral image
        struct integral_image *optimized_iimage = create_integral_img(width, height);
        // Compute integral image
        optimized_function(image, optimized_iimage->width, optimized_iimage->height, optimized_iimage->data);

        //print_debug(original_iimage->data, optimized_iimage->data, width, height);
        if (!are_float_matrices_equal(original_iimage->data, optimized_iimage->data, width, height)) {
            all_functions_equal = false;
            printf("Error: The integral images are not equal.\n");
        }

        free(optimized_iimage->padded_data);
        free(optimized_iimage);
    }

    free(original_iimage->padded_data);
    free(original_iimage);
    free(image);

    return all_functions_equal;
}


// Prints all values of two given matrices.
void print_debug(float *iimage1, float *iimage2, int width, int height) {
    printf("print debug\n");
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            printf("i1: %f, i2: %f\n", iimage1[i * width + j], iimage2[i * width + j]);
        }
    }
}