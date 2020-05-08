#pragma once

#include <vector>

#include "integral_image.h"


// Easier for debugging.
// Creates an integral image based on a custom array, the base function and a list of other functions. 
// The results of the integral images are being compared.
// If one of the results difers from the base implementation false is being returned.
// Messages clarifying the equality of the results are being printed.
bool validate_integral_image_custom_matrix(void (*original_function)(float *, int, int, float *),
                                  const std::vector<void (*)(float *, int, int, float *)> &test_functions);


// Prints all values of two given matrices.
void print_debug(float *iimage1, float *iimage2, int width, int height);