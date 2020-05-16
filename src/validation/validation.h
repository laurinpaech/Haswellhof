#pragma once

#include "integral_image.h"
#include "interest_point.h"
#include "fasthessian.h"
#include "descriptor.h"
#include "descriptor_opt.h"
#include "fasthessian_opt.h"

#include <vector>

// Creates an integral image given an image, its corresponding height and width, the base function and a list of other functions. 
// The results of the integral images are being compared.
// If one of the results difers from the base implementation false is being returned.
// Messages clarifying the equality of the results are being printed.
bool validate_integral_image(void (*original_function)(float *, struct integral_image *), 
                             const std::vector<void (*)(float *, struct integral_image *)> &test_functions, 
                             int width, int height, float *image);

bool validate_compute_response_layer_custom_matrix(void (*original_function)(struct fasthessian *),
                                                   const std::vector<void (*)(struct fasthessian *)> &test_functions);

bool validate_compute_response_layers(void (*original_function)(struct fasthessian *),
                                      const std::vector<void (*)(struct fasthessian *)> &test_functions,
                                      struct integral_image *iimage);

bool validate_compute_response_layer_with_padding(
    void (*original_function)(struct response_layer *, struct integral_image *),
    const std::vector<void (*)(struct response_layer *, struct integral_image *)> &test_functions,
    float* original_image, int width, int height);

bool validate_get_interest_points(void (*original_function)(struct fasthessian *, std::vector<struct interest_point> *),
                                  const std::vector<void (*)(struct fasthessian *, std::vector<struct interest_point> *)> &test_functions,
                                  struct fasthessian *fh);

bool validate_get_msurf_descriptors(void (*original_function)(struct integral_image *, struct interest_point *),
                                    const std::vector<void (*)(struct integral_image *, struct interest_point *)> &test_functions,
                                    struct integral_image *iimage, const std::vector<struct interest_point> *interest_points);
