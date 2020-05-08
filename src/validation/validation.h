#pragma once

#include "integral_image.h"
#include "interest_point.h"
#include "fasthessian.h"
#include "descriptor.h"
#include "descriptor_opt.h"
#include "fasthessian_opt.h"

#include <vector>

bool validate_compute_response_layer(void (*original_function)(struct response_layer *, struct integral_image *), 
                                     const std::vector<void (*)(struct response_layer *, struct integral_image *)> &test_functions,
                                     struct response_layer* layer, struct integral_image* iimage);

bool validate_get_msurf_descriptors(void (*original_function)(struct integral_image *, struct interest_point *),
                                    const std::vector<void (*)(struct integral_image *, struct interest_point *)> &test_functions,
                                    struct integral_image *iimage, const std::vector<struct interest_point> *interest_points);
