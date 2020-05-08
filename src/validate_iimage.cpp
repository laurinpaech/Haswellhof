
#include "validate_iimage.h"

#include <stdio.h>

// Creates an integral image given an image, its corresponding height and width, the base function and a list of other functions. 
// The results of the integral images are being compared.
// If one of the results difers from the base implementation false is being returned.
// Messages clarifying the equality of the results are being printed.
bool validate_iimage(void (*original_function)(float *, int, int, float *),
                     std::vector<void (*)(float *, int, int, float *)> test_functions, int width, int height,
                     float *image) {
    // Create integral image
    struct integral_image *original_iimage = create_integral_img(width, height);
    // Compute integral image
    original_function(image, original_iimage->width, original_iimage->height, original_iimage->data);
    bool all_functions_equal = true;
    for (auto optimized_function : test_functions) {
        // Create integral image
        struct integral_image *optimized_iimage = create_integral_img(width, height);
        // Compute integral image
        optimized_function(image, optimized_iimage->width, optimized_iimage->height, optimized_iimage->data);

        bool equal = are_equal(original_iimage->data, optimized_iimage->data, width, height);

        if (equal == false) {
            all_functions_equal = false;
            printf("Error: The integral images are not equal.\n");
        } else {
            printf("The integral images are equal\n");
        }
        free(optimized_iimage->data);
        free(optimized_iimage);
    }

    free(original_iimage->data);
    free(original_iimage);

    return all_functions_equal;
}

// Easier for debugging.
// Creates an integral image based on a custom array, the base function and a list of other functions. 
// The results of the integral images are being compared.
// If one of the results difers from the base implementation false is being returned.
// Messages clarifying the equality of the results are being printed.
bool validate_iimage_custom_array(void (*original_function)(float *, int, int, float *),
                                  std::vector<void (*)(float *, int, int, float *)> test_functions) {
    int width = 4, height = 4;
    float *image = (float *)malloc(height * width * sizeof(float));
    int counter = 0;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            image[i * width + j] = counter;
            counter++;
        }
    }

    // float image[height][width] = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 }, { 7.0, 8.0, 9.0 } };
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
        bool equal = are_equal(original_iimage->data, optimized_iimage->data, width, height);
        if (equal == false) {
            all_functions_equal = false;
            printf("Error: The integral images are not equal.\n");
        } else {
            printf("The integral images are equal\n");
        }
        free(optimized_iimage->data);
        free(optimized_iimage);
    }

    free(original_iimage->data);
    free(original_iimage);
    free(image);

    return all_functions_equal;
}

// Compares two matrices and checks if the values are equal.
// Returns true if all values of the matrix are equal, false otherwise.
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

// Prints all values of two given matrices.
void print_debug(float *iimage1, float *iimage2, int width, int height) {
    printf("print debug\n");
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            printf("i1: %f, i2: %f\n", iimage1[i * width + j], iimage2[i * width + j]);
        }
    }
}