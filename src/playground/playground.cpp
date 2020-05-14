#include "playground.h"

// NOTE: Fill this with whatever you want (preferably don't commit/push it though)

#include <iostream>

#include "benchmarking.h"
#include "descriptor.h"
#include "fasthessian_opt.h"
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
/*
    for (size_t i = 0; i < optimized_iimage->height; i++) {
        printf("row: %i \n", i);
        for (size_t j = 0; j < optimized_iimage->width; j++) {
            printf("%f ", optimized_iimage->data[i * optimized_iimage->width + j]);
        }
        printf("\n\n");
    }
    */

    // Create integral image
    struct integral_image *iimage = create_integral_img(width, height);
    // Compute integral image
    compute_integral_img(image, iimage->width, iimage->height, iimage->data);

    
        struct fasthessian *fh_original = create_fast_hessian(iimage);
        int start =6;

       int test_layer_size = fh_original->total_layers;
        // Create octaves with response layers
        create_response_map(fh_original);
        printf("ORIGINAL: \n");

        // Compute responses for every layer
            for (int i = start; i < test_layer_size; ++i) {
                    compute_response_layer_debug(fh_original->response_map[i], fh_original->iimage);
            }

        // Fast-Hessian
        struct fasthessian *optimized_fh = create_fast_hessian(optimized_iimage);

        // Create octaves with response layers
        create_response_map(optimized_fh);

        printf("Before compute response layer\n");
        // Compute responses for every layer
            printf("OPTIMIZED: \n");

        for (int i = start; i < test_layer_size; ++i) {
            compute_response_layer_with_padding(optimized_fh->response_map[i], optimized_fh->iimage);
        }

        if (fh_original->total_layers != optimized_fh->total_layers) {
            printf(
                "compute_response_layer() does not match original function - the number of layers "
                "differ.\n");
            // If the numbers of layers doesn't match, don't do any further tests.
        }

        // Compare each layer of each test function with the original
        for (int i = start; i < test_layer_size; i++) {
            struct response_layer *optimized_layer = optimized_fh->response_map[i];
            struct response_layer *original_layer = fh_original->response_map[i];
            printf("ORIGINAL size: %i, OPTIMIZED size: %i\n", original_layer->height, optimized_layer->height);
            if (original_layer->height != optimized_layer->height || original_layer->width != optimized_layer->width) {
                printf(
                    "compute_response_layer() does not match original function - the layer sizes "
                    "differ.\n");
                // If the sizes of layers don't match, don't do any further tests.
                continue;
            }
            // print_debug(original_layer->response, optimized_layer->response,
            // original_layer->height,original_layer->width);

            if (!are_float_matrices_equal(original_layer->response, optimized_layer->response, original_layer->height,
                                          original_layer->width) ||
                !are_bool_matrices_equal(original_layer->laplacian, optimized_layer->laplacian, original_layer->height,
                                         original_layer->width)) {
                printf("compute_response_layer() does not match original function\n");
            }
            printf("After compute response layer\n");
        }
                printf("successy\n");


        for (int i = 0; i < NUM_LAYERS; ++i) {
            free(optimized_fh->response_map[i]->response);
            free(optimized_fh->response_map[i]->laplacian);
            free(optimized_fh->response_map[i]);
            free(fh_original->response_map[i]->response);
            free(fh_original->response_map[i]->laplacian);
            free(fh_original->response_map[i]);
        }
        free(optimized_fh);
        free(fh_original);
        
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

    int padded_lobe = (LARGEST_FILTER_SIZE / 3) + 1;
    int padded_border = ((LARGEST_FILTER_SIZE - 1) / 2) + 1;
    int padded_width = optimized_iimage->width + 2 * padded_lobe;
    for (size_t i = 0; i < optimized_iimage->height + 2 * padded_border; i++) {
        printf("row: %i \n", i);
        for (size_t j = 0; j < padded_width; j++) {
            printf("%f ", optimized_iimage->data[i * padded_width + j]);
        }
        printf("\n\n");
    }

    free(optimized_iimage->data);
    free(optimized_iimage);
    free(image);
}
