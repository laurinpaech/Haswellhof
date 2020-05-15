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
#include "validation.h"

void playground_function1(float *image, int width, int height) {
// Create integral image
	struct integral_image* iimage = create_integral_img(width, height);
	// Compute integral image
	compute_integral_img(image, iimage);

	std::vector<void (*)(struct response_layer *, struct integral_image *)> test_functions;
	test_functions.push_back(compute_response_layer_uncoditional);

	 bool valid =validate_compute_response_layer_with_padding(compute_response_layer,test_functions, image, width, height);
	
        if (valid) {
            printf("COMPUTE RESPONSE LAYER VALIDATION:    \033[0;32mSUCCESS!\033[0m\n");
        } else {
            printf("COMPUTE RESPONSE LAYER VALIDATION:    \033[1;31mFAILED!\033[0m\n");
        }
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
