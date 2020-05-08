#pragma once

#include <vector>

#include "integral_image.h"
#include "integral_image_opt.h"


// Creates an integral image given an image, its corresponding height and width, the base function and a list of other functions. 
// The results of the integral images are being compared.
// If one of the results difers from the base implementation false is being returned.
// Messages clarifying the equality of the results are being printed.
bool validate_iimage(void (*original_function)(float *, int, int, float *), std::vector<void (*)(float *, int, int, float *)> test_functions, 
int width, int height, float* image);

// Easier for debugging.
// Creates an integral image based on a custom array, the base function and a list of other functions. 
// The results of the integral images are being compared.
// If one of the results difers from the base implementation false is being returned.
// Messages clarifying the equality of the results are being printed.
bool validate_iimage_custom_array(void (*original_function)(float *, int, int, float *),
                                  std::vector<void (*)(float *, int, int, float *)> test_functions);

// Compares two matrices and checks if the values are equal.
// Returns true if all values of the matrix are equal, false otherwise.
bool are_equal(float* iimage1, float* iimage2, int width, int height);

// Prints all values of two given matrices.
void print_debug(float *iimage1, float *iimage2, int width, int height);