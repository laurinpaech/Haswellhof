#pragma once

#include "integral_image.h"
#include "integral_image_opt.h"

bool evaluate_iimage(void (*original_function)(float *, int, int, float *), void (*optimized_function)(float *, int, int, float *), 
int width, int height, float* image);

bool evaluate_iimage_custom_array(void (*original_function)(float *, int, int, float *),
                                  void (*optimized_function)(float *, int, int, float *));

bool are_equal(float* iimage1, float* iimage2, int width, int height);

void print_debug(float *iimage1, float *iimage2, int width, int height);