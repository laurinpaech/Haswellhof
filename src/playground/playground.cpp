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
    test_functions.push_back(compute_response_layer_Dyy_leftcorner);
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
    test_functions.push_back(compute_response_layer_Dyy_leftcorner);
    // bool valid = validate_compute_response_layer(compute_response_layer, test_functions, iimage);

    bool valid = validate_compute_response_layer(compute_response_layer, test_functions, iimage);
    if (valid) {
        printf("COMPUTE RESPONSE LAYER VALIDATION:  \033[0;32mSUCCESS!\033[0m\n");
    } else {
        printf("COMPUTE RESPONSE LAYER VALIDATION:  \033[1;31mFAILED!\033[0m!\n");
    }
}

void playground_function3(struct integral_image *iimage, const char* image_name) {
    // Insert all compute_integral_img functions for benchmarking here
    std::vector<void (*)(struct response_layer *, struct integral_image *)> functions;
    functions.push_back(compute_response_layer);
    functions.push_back(compute_response_layer_Dyy_leftcorner);
    // functions.push_back(compute_integral_img_faster_alg);

    // TODO: (carla) how to take the image_name as argument
    struct benchmark_data default_data("never_gonna_give_you_up", iimage->width, iimage->height,
                                       "compute_response_layer", -1, 1 + iimage->height*iimage->width*13);
                                           std::vector<struct benchmark_data> data;
    data.push_back(default_data);

    struct benchmark_data default_data2("never_gonna_give_you_up", iimage->width, iimage->height,
                                       "compute_response_layer_Dyy_leftcorner", -1, 1 + iimage->height*iimage->width*13);

    data.push_back(default_data2);

    bench_compute_response_layer(functions, iimage, data);


    // Insert all respective benchmarking info for compute_integral_img here

    save_benchmark_data(data);

}
