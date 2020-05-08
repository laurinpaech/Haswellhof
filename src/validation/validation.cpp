#include "validation.h"
#include "helper.h"

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>



// Creates an integral image given an image, its corresponding height and width, the base function and a list of other functions. 
// The results of the integral images are being compared.
// If one of the results difers from the base implementation false is being returned.
// Messages clarifying the equality of the results are being printed.
bool validate_integral_image(void (*original_function)(float *, int, int, float *),
                     const std::vector<void (*)(float *, int, int, float *)> &test_functions, int width, int height,
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

        if (!are_matrices_equal(original_iimage->data, optimized_iimage->data, width, height)) {
            all_functions_equal = false;
            printf("Error: The integral images are not equal.\n");
        }

        free(optimized_iimage->data);
        free(optimized_iimage);
    }

    free(original_iimage->data);
    free(original_iimage);

    return all_functions_equal;
}

bool validate_compute_response_layer(void (*original_function)(struct response_layer *, struct integral_image *), 
                                     const std::vector<void (*)(struct response_layer *, struct integral_image *)> &test_functions,
                                     struct response_layer* layer, struct integral_image* iimage) {



}

bool validate_get_msurf_descriptors(void (*original_function)(struct integral_image *, struct interest_point *),
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

            bool valid = compare_arrays_close(ref_ipoint.descriptor, test_ipoint.descriptor, 64, EPSILON);
            if (!valid) {
                printf("get_msurf_descriptor() test function %d does not match original function\n", j);
            }

            // Updating flag indicating if all functions are valid
            all_valid &= valid;

        }

    } 

    return all_valid;

}

