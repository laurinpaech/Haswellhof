#include "validation.h"

#include "helper.h"

#include <stdbool.h>
#include <stdio.h>

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
