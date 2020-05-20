#include "validation.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "helper.h"
#include "validation_integral_image.h"

//#define DEBUG_INFO
#define VALIDATION_PRECISION (1e-3)

// Creates an integral image given an image, its corresponding height and width, the base function and a list of other
// functions. The results of the integral images are being compared. If one of the results difers from the base
// implementation false is being returned. Messages clarifying the equality of the results are being printed.
bool validate_integral_image(void (*original_function)(float *, struct integral_image *),
                             const std::vector<void (*)(float *, struct integral_image *)> &test_functions, int width,
                             int height, float *image) {
    // Create integral image
    struct integral_image *original_iimage = create_integral_img(width, height);
    // Compute integral image
    original_function(image, original_iimage);
    bool all_functions_equal = true;
    for (auto optimized_function : test_functions) {
        // Create integral image
        struct integral_image *optimized_iimage = create_integral_img(width, height);
        // Compute integral image
        optimized_function(image, optimized_iimage);

        if (!are_float_matrices_equal(original_iimage->data, optimized_iimage->data, width, height)) {
            all_functions_equal = false;
            printf("Error: The integral images are not equal.\n");
        }

        free(optimized_iimage->padded_data);
        free(optimized_iimage);
    }

    free(original_iimage->padded_data);
    free(original_iimage);

    return all_functions_equal;
}

bool validate_compute_response_layer_custom_matrix(void (*original_function)(struct fasthessian *),
                                                   const std::vector<void (*)(struct fasthessian *)> &test_functions) {
    int width = 32, height = 32;
    float *image = (float *)malloc(height * width * sizeof(float));
    for (int i = 0; i < height; ++i) {
        int counter = 0;
        for (int j = 0; j < width; ++j) {
            image[i * width + j] = counter;
            counter++;
        }
    }

    // Create integral image
    struct integral_image *iimage = create_integral_img(width, height);

    // Compute integral image
    compute_integral_img(image, iimage);

    bool valid = validate_compute_response_layers(original_function, test_functions, iimage);

    free(iimage->padded_data);
    free(iimage);

    return valid;
}

bool validate_compute_response_layers(void (*original_function)(struct fasthessian *),
                                      const std::vector<void (*)(struct fasthessian *)> &test_functions,
                                      struct integral_image *iimage) {
    bool all_valid = true;

    // Fast-Hessian
    struct fasthessian *original_fh = create_fast_hessian(iimage);
    // Create octaves with response layers
    create_response_map(original_fh);

    // Compute responses for every layer
    original_function(original_fh);

    for (int j = 0; j < test_functions.size(); ++j) {
        // Fast-Hessian
        struct fasthessian *optimized_fh = create_fast_hessian(iimage);

        // Create octaves with response layers
        create_response_map(optimized_fh);

        // Compute responses for every layer
        test_functions[j](optimized_fh);

        if (original_fh->total_layers != optimized_fh->total_layers) {
            printf(
                "compute_response_layer() test function %d does not match original function - the number of layers "
                "differ.\n",
                j);
            all_valid = false;
            // If the numbers of layers doesn't match, don't do any further tests.
            continue;
        }

        // Compare each layer of each test function with the original
        for (int i = 0; i < original_fh->total_layers; i++) {
            struct response_layer *optimized_layer = optimized_fh->response_map[i];
            struct response_layer *original_layer = original_fh->response_map[i];
            if (original_layer->height != optimized_layer->height || original_layer->width != optimized_layer->width) {
                printf(
                    "compute_response_layer() test function %d does not match original function - the layer sizes "
                    "difer.\n",
                    j);
                all_valid = false;
                // If the sizes of layers don't match, don't do any further tests.
                continue;
            }
#ifdef DEBUG_INFO
            print_debug(original_layer->response, optimized_layer->response, original_layer->height,
                        original_layer->width);
#endif

            if (!are_float_matrices_equal(original_layer->response, optimized_layer->response, original_layer->width,
                original_layer->height) ||
                !are_bool_matrices_equal(original_layer->laplacian, optimized_layer->laplacian, original_layer->height,
                                         original_layer->width)) {
                printf("compute_response_layer() test function %d does not match original function for layer %i\n", j, i);
                all_valid = false;
            }
        }
        for (int i = 0; i < NUM_LAYERS; ++i) {
            free(optimized_fh->response_map[i]->response);
            free(optimized_fh->response_map[i]->laplacian);
            free(optimized_fh->response_map[i]);
        }
        free(optimized_fh);
    }

    for (int i = 0; i < NUM_LAYERS; ++i) {
        free(original_fh->response_map[i]->response);
        free(original_fh->response_map[i]->laplacian);
        free(original_fh->response_map[i]);
    }
    free(original_fh);

    return all_valid;
}

bool validate_compute_response_layer_with_padding(
    void (*original_function)(struct response_layer *, struct integral_image *),
    const std::vector<void (*)(struct response_layer *, struct integral_image *)> &test_functions,
    float *original_image, int width, int height) {
    bool all_valid = true;

    // Create integral image
    struct integral_image *original_integral_image = create_integral_img(width, height);
    // Compute integral image
    compute_integral_img(original_image, original_integral_image);

    // Fast-Hessian
    struct fasthessian *original_fh = create_fast_hessian(original_integral_image);
    // Create octaves with response layers
    create_response_map(original_fh);

    // Compute responses for every layer
    for (int i = 0; i < original_fh->total_layers; i++) {
#ifdef DEBUG_INFO
        printf("responselayer height: %i, width: %i\n", original_fh->response_map[i]->height, original_fh->response_map[i]->width);
#endif
        original_function(original_fh->response_map[i], original_integral_image);
    }

    // Create padded integral image
    struct integral_image *padded_integral_image = create_padded_integral_img(width, height);
    // Compute padded integral image
    compute_padded_integral_img(original_image, padded_integral_image);

    for (int j = 0; j < test_functions.size(); ++j) {
        // Fast-Hessian
        struct fasthessian *optimized_fh = create_fast_hessian(padded_integral_image);
        // Create octaves with response layers
        create_response_map(optimized_fh);
        // Compute responses for every layer
        for (int i = 0; i < optimized_fh->total_layers; i++) {
            test_functions[j](optimized_fh->response_map[i], padded_integral_image);
        }

        if (original_fh->total_layers != optimized_fh->total_layers) {
            printf(
                "compute_response_layer() test function %d does not match original function - the number of layers "
                "differ.\n",
                j);
            all_valid = false;
            // If the numbers of layers doesn't match, don't do any further tests.
            continue;
        }

        // Compare each layer of each test function with the original
        for (int i = 0; i < original_fh->total_layers; i++) {
            struct response_layer *optimized_layer = optimized_fh->response_map[i];
            struct response_layer *original_layer = original_fh->response_map[i];
            if (original_layer->height != optimized_layer->height || original_layer->width != optimized_layer->width) {
                printf(
                    "compute_response_layer() test function %d does not match original function - the layer sizes "
                    "differ.\n",
                    j);
                all_valid = false;
                // If the sizes of layers don't match, don't do any further tests.
                continue;
            }
#ifdef DEBUG_INFO
            print_debug(original_layer->response, optimized_layer->response, original_layer->height,
                        original_layer->width);
#endif

            if (!compare_arrays_close(original_layer->response, optimized_layer->response,
                                    (original_layer->height)*(original_layer->width), VALIDATION_PRECISION) ||
                !are_bool_matrices_equal(original_layer->laplacian, optimized_layer->laplacian, original_layer->height,
                                         original_layer->width)) {
                printf("compute_response_layer() test function %d does not match original function on layer %d/%d\n", j, i,original_fh->total_layers);
                all_valid = false;
            }
        }
        for (int i = 0; i < NUM_LAYERS; ++i) {
            free(optimized_fh->response_map[i]->response);
            free(optimized_fh->response_map[i]->laplacian);
            free(optimized_fh->response_map[i]);
        }
        free(optimized_fh);
    }

    for (int i = 0; i < NUM_LAYERS; ++i) {
        free(original_fh->response_map[i]->response);
        free(original_fh->response_map[i]->laplacian);
        free(original_fh->response_map[i]);
    }
    free(original_fh);
    free(original_integral_image->padded_data);
    free(original_integral_image);
    free(padded_integral_image->padded_data);
    free(padded_integral_image);

    return all_valid;
}

bool validate_get_interest_points(void (*original_function)(struct fasthessian *, std::vector<struct interest_point> *),
                                  const std::vector<void (*)(struct fasthessian *, std::vector<struct interest_point> *)> &test_functions,
                                  struct fasthessian *fh) {

    std::vector<struct interest_point> original_interest_points;
    original_function(fh, &original_interest_points);

    bool all_valid = true;

    for (int j = 0; j < test_functions.size(); ++j) {

        std::vector<struct interest_point> test_interest_points;
        test_functions[j](fh, &test_interest_points);

        bool valid = true;

        // Checking if original and test interest point vectors have same size
        if (original_interest_points.size() != test_interest_points.size()) {

            printf("ERROR: interest point test function %d does not match original function (size mismatch)\n", j);

            valid = false;

            // Updating flag indicating if all functions are valid
            all_valid &= valid;

            continue;

        }

        // Iterating through all interest points and comparing them
        for (int i = 0; i < original_interest_points.size(); ++i) {

            struct interest_point original_ipoint = original_interest_points[i];
            struct interest_point test_ipoint = test_interest_points[i];

            // Checking if all fields except descriptor arrays are identical
            if (original_ipoint.x != test_ipoint.x
                || original_ipoint.x != test_ipoint.x
                || original_ipoint.scale != test_ipoint.scale
                || original_ipoint.upright != test_ipoint.upright
                || original_ipoint.orientation != test_ipoint.orientation
                || original_ipoint.laplacian != test_ipoint.laplacian) {
                printf("ERROR: interest point test function %d does not match original function (interest point mismatch)\n", j);
                /*
                printf(" Interest point %d mismatch:\n", i);
                printf("  x: %d <-> %d\n", original_ipoint.x, test_ipoint.x);
                printf("  y: %d <-> %d\n", original_ipoint.y, test_ipoint.y);
                printf("  scale: %f <-> %f\n", original_ipoint.scale, test_ipoint.scale);
                printf("  upright: %d <-> %d\n", original_ipoint.upright ? 1 : 0, test_ipoint.upright ? 1 : 0);
                printf("  orientation: %f <-> %f\n", original_ipoint.orientation, test_ipoint.orientation);
                printf("  scale:  %f <-> %f\n", original_ipoint.laplacian ? 1 : 0, test_ipoint.laplacian  ? 1 : 0);
                */
                valid = false;
            }

        }

        // Updating flag indicating if all functions are valid
        all_valid &= valid;

    }

    return all_valid;

}

bool validate_get_msurf_descriptors(
    void (*original_function)(struct integral_image *, struct interest_point *),
    const std::vector<void (*)(struct integral_image *, struct interest_point *)> &test_functions,
    struct integral_image *iimage, const std::vector<struct interest_point> *interest_points) {
    // Copying over interest points to new vector
    std::vector<struct interest_point> original_interest_points = *interest_points;

    bool all_valid = true;

    for (int i = 0; i < interest_points->size(); ++i) {

        struct interest_point ref_ipoint = interest_points->at(i);
        original_function(iimage, &ref_ipoint);

        bool valid = true;
        for (int j = 0; j < test_functions.size(); ++j) {
            struct interest_point test_ipoint = interest_points->at(i);

            test_functions[j](iimage, &test_ipoint);

            bool valid = compare_arrays_close(ref_ipoint.descriptor, test_ipoint.descriptor, 64, VALIDATION_PRECISION);
            if (!valid) {
                printf("ERROR: MSURF descriptor test function %d does not match original function\n", j);
            }

            // Updating flag indicating if all functions are valid
            all_valid &= valid;
        }
    }

    return all_valid;
}
