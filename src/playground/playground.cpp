#include "playground.h"

// NOTE: Fill this with whatever you want (preferably don't commit/push it though)

#include <iostream>

#include "benchmark_data_to_file.h"
#include "benchmarking.h"
#include "descriptor.h"
#include "helper.h"
#include "integral_image.h"

void playground_function1() {
    std::vector<void (*)(struct response_layer *, struct integral_image *)> test_functions;
    // test_functions.push_back(compute_response_layer_Dyy_leftcorner);
    // test_functions.push_back(compute_response_layer_Dyy_top);
    // test_functions.push_back(compute_response_layer_Dyy_top_mid);
    test_functions.push_back(super_sonic_Dyy);
    // bool valid = validate_compute_response_layer(compute_response_layer, test_functions, iimage);

    bool valid = validate_compute_response_layer_custom_matrix(compute_response_layer, test_functions);
    if (valid) {
        printf("COMPUTE RESPONSE LAYER VALIDATION:  \033[0;32mSUCCESS!\033[0m\n");
    } else {
        printf("COMPUTE RESPONSE LAYER VALIDATION:  \033[1;31mFAILED!\033[0m!\n");
    }
}

void playground_function2(struct integral_image *iimage) {
    std::vector<void (*)(struct response_layer *, struct integral_image *)> test_functions;
    test_functions.push_back(super_sonic_Dyy);
    // bool valid = validate_compute_response_layer(compute_response_layer, test_functions, iimage);

    bool valid = validate_compute_response_layer(compute_response_layer, test_functions, iimage);
    if (valid) {
        printf("COMPUTE RESPONSE LAYER VALIDATION:  \033[0;32mSUCCESS!\033[0m\n");
    } else {
        printf("COMPUTE RESPONSE LAYER VALIDATION:  \033[1;31mFAILED!\033[0m!\n");
    }
}

