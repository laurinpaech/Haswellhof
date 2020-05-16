#include "playground.h"

// NOTE: Fill this with whatever you want (preferably don't commit/push it though)

#include <iostream>

#include "benchmarking.h"
#include "descriptor.h"
#include "fasthessian_opt.h"
#include "helper.h"
#include "integral_image.h"
#include "integral_image_opt.h"
#include "validation.h"
#include "validation_integral_image.h"
#include "benchmark_data_to_file.h"

void playground_function1(float *image, int width, int height) {
    // Create integral image
    struct integral_image *iimage = create_integral_img(width, height);
    // Compute integral image
    compute_integral_img(image, iimage);

    std::vector<void (*)(struct response_layer *, struct integral_image *)> test_functions;
    test_functions.push_back(compute_response_layer_unconditional);

    bool valid =
        validate_compute_response_layer_with_padding(compute_response_layer, test_functions, image, width, height);

    if (valid) {
        printf("COMPUTE RESPONSE LAYER VALIDATION:    \033[0;32mSUCCESS!\033[0m\n");
    } else {
        printf("COMPUTE RESPONSE LAYER VALIDATION:    \033[1;31mFAILED!\033[0m\n");
    }
}

void playground_function3(float *image, int width, int height) {
    std::vector<struct benchmark_data> all_benchmark_data;

    // Create integral image
    struct integral_image *iimage = create_integral_img(width, height);
    // Compute integral image
    compute_integral_img(image, iimage);

    printf("compute_response_layer start\n");

    std::vector<void (*)(struct fasthessian *)> functions;
    functions.push_back(compute_response_layers);
    functions.push_back(compute_response_layers_at_once);

    struct benchmark_data default_data("never_gonna_give_you_up", width, height, "compute_response_layer", -1,
                                       (1 + height * width * 13));
    struct benchmark_data data1("never_gonna_give_you_up", width, height, "compute_response_layers_at_once", -1,
                                (1 + height * width * 13));

    std::vector<struct benchmark_data> data;
    data.push_back(default_data);
    data.push_back(data1);

    bench_compute_response_layer(functions, iimage, data);

    all_benchmark_data.insert(all_benchmark_data.end(), data.begin(), data.end());

    // Create integral image
    struct integral_image *padded_iimage = create_padded_integral_img(width, height);
    // Compute integral image
    compute_padded_integral_img(image, iimage);

    std::vector<void (*)(struct fasthessian *)> padded_functions;
    padded_functions.push_back(compute_response_layers_unconditional);

    struct benchmark_data padded_data("never_gonna_give_you_up", width, height, "compute_response_layers_unconditional", -1,
                                       (1 + height * width * 13));
    std::vector<struct benchmark_data> data_padded_functions;
    data_padded_functions.push_back(padded_data);
    bench_compute_response_layer(padded_functions, padded_iimage, data_padded_functions);
    all_benchmark_data.insert(all_benchmark_data.end(), data_padded_functions.begin(), data_padded_functions.end());

    printf("compute_response_layer end\n");
        save_benchmark_data(all_benchmark_data);
    
        free(iimage->padded_data);
        free(iimage);

        free(padded_iimage->padded_data);
        free(padded_iimage);       

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

    compute_padded_integral_img(image, optimized_iimage);
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

    free(optimized_iimage->padded_data);
    free(optimized_iimage);
    free(image);
}
