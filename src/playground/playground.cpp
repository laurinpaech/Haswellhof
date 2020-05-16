#include "playground.h"

// NOTE: Fill this with whatever you want (preferably don't commit/push it though)

#include "benchmarking.h"
#include "descriptor.h"
#include "fasthessian_opt.h"
#include "fasthessian_opt_flat.h"
#include "helper.h"
#include "integral_image.h"
#include "integral_image_opt.h"
#include "validation.h"
#include "validation_integral_image.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


void playground_function1(struct integral_image *iimage) {
    // Fast-Hessian
    struct fasthessian *fh = create_fast_hessian(iimage);

    // Create octaves with response layers
    create_response_map(fh);

    // Compute responses for every layer
    compute_response_layers(fh);

    // Getting interest points with non-maximum supression
    std::vector<struct interest_point> interest_points;
    get_interest_points(fh, &interest_points);

    // ---------------------- Trying out fast hessian flat ------------------------- //

    struct fasthessian_flat fh_flat;
    create_fast_hessian_flat_and_response_map(iimage, &fh_flat);

    // Compute responses for every layer
    compute_response_layers_flat(&fh_flat);

    // Getting interest points with non-maximum supression
    std::vector<struct interest_point> interest_points_flat;
    get_interest_points_flat(&fh_flat, &interest_points_flat);

    std::cout << "interest_points.size: " << interest_points.size() << std::endl;
    std::cout << "interest_points_flat.size: " << interest_points_flat.size() << std::endl;

    assert(interest_points.size() == interest_points_flat.size());

    for (int i = 0; i < interest_points_flat.size(); ++i) {
        if (interest_points[i].x != interest_points_flat[i].x) {
            std::cout << "flat: x not equal" << std::endl;
        }
        if (interest_points[i].y != interest_points_flat[i].y) {
            std::cout << "flat: y not equal" << std::endl;
        }
        if (interest_points[i].scale != interest_points_flat[i].scale) {
            std::cout << "flat: scale not equal" << std::endl;
        }
        if (interest_points[i].orientation != interest_points_flat[i].orientation) {
            std::cout << "flat: orientation not equal" << std::endl;
        }
        if (interest_points[i].upright != interest_points_flat[i].upright) {
            std::cout << "flat: upright not equal" << std::endl;
        }
        if (interest_points[i].laplacian != interest_points_flat[i].laplacian) {
            std::cout << "flat: laplacian not equal" << std::endl;
        }
        // Not set yet
        // if (compare_arrays(interest_points[i].descriptor, interest_points_flat[i].descriptor, 65)) {
        //    std::cout << "flat: descriptors not equal" << std::endl;
        //}
    }

    aligned_free(fh_flat.response_map[0].response);

    for (int i = 0; i < NUM_LAYERS; ++i) {
        free(fh->response_map[i]->response);
        free(fh->response_map[i]->laplacian);
        free(fh->response_map[i]);
    }
    free(fh);
}
