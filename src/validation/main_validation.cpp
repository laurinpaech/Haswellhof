// has to be defined before stb includes
#define STB_IMAGE_IMPLEMENTATION


#include "stb_image.h"
#include "integral_image.h"
#include "fasthessian.h"
#include "interest_point.h"
#include "descriptor.h"

#include "integral_image_opt.h"
#include "integral_image_simd.h"
#include "descriptor_opt.h"
// #include "descriptor_opt_gen.h"
#include "fasthessian_opt.h"
// #include "fasthessian_opt_gen.h"

#include "validation.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <vector>

#define VALIDATE_INTEGRAL_IMAGE
#define VALIDATE_INTEGRAL_IMAGE_PADDED
#define VALIDATE_INTEGRAL_IMAGE_INT
#define VALIDATE_INTEGRAL_IMAGE_INT_PADDED
#define VALIDATE_COMPUTE_RESPONSE_LAYER
#define VALIDATE_COMPUTE_RESPONSE_LAYER_PADDED
#define VALIDATE_GET_INTEREST_POINTS
#define VALIDATE_GET_MSURF_DESCRIPTORS


int main(int argc, char const *argv[])
{
    if (argc != 3) {
        printf("Invaid argument count.\nUsage:\n\t./surf image_path descriptor_file_target_path\n");
        return 1;
    }
	int width, height, channels;

	// Load image
	stbi_ldr_to_hdr_gamma(1.0f);
	float* image = stbi_loadf(argv[1], &width, &height, &channels, STBI_grey);

    if (!image) {
		printf("Could not open or find image\n");
        return -1;
    }

#if defined(VALIDATE_INTEGRAL_IMAGE_INT) || defined(VALIDATE_INTEGRAL_IMAGE_INT_PADDED)
    int width_int, height_int, channels_int;

    // Load uint8_t version of image
    stbi_ldr_to_hdr_gamma(1.0f);
    uint8_t *image_int = stbi_load(argv[1], &width_int, &height_int, &channels_int, STBI_grey);

    if (!image_int) {
        printf("Could not open or find int image\n");
        return -1;
    }

    assert(width == width_int && height == height_int && channels == channels_int);
#endif

	// Create integral image
	struct integral_image* iimage = create_integral_img(width, height);

#ifdef VALIDATE_INTEGRAL_IMAGE
    {
        std::vector<void (*)(float *, struct integral_image *)> test_functions;
        test_functions.push_back(compute_integral_img_faster_alg);
        //test_functions.push_back(compute_padded_integral_img_new);
        //test_functions.push_back(compute_padded_integral_img);

        bool valid = validate_integral_image(compute_integral_img, test_functions, width, height, image, false);
        if (valid) {
            printf("INTEGRAL IMAGE VALIDATION:                 \033[0;32mSUCCESS!\033[0m\n");
        } else {
            printf("INTEGRAL IMAGE VALIDATION:                 \033[1;31mFAILED!\033[0m\n");
        }
    }
#endif

#ifdef VALIDATE_INTEGRAL_IMAGE_PADDED
    {
        std::vector<void (*)(float *, struct integral_image *)> test_functions;
        test_functions.push_back(compute_padded_integral_img_new);
        test_functions.push_back(compute_padded_integral_img_faster_alg);

        bool valid = validate_integral_image(compute_padded_integral_img, test_functions, width, height, image, true);
        if (valid) {
            printf("INTEGRAL IMAGE PADDED VALIDATION:          \033[0;32mSUCCESS!\033[0m\n");
        } else {
            printf("INTEGRAL IMAGE PADDED VALIDATION:          \033[1;31mFAILED!\033[0m\n");
        }
    }
#endif

#ifdef VALIDATE_INTEGRAL_IMAGE_INT
    {
        //if (width < 2895) {
        if (width <= 2048) {
            std::vector<void (*)(uint8_t *, struct integral_image *)> test_functions;
            test_functions.push_back(compute_integral_img_simd_int);
            test_functions.push_back(compute_integral_img_simd_early_cast_int);
            //test_functions.push_back(compute_integral_img_simd_original);

            bool valid = validate_integral_image_int(compute_integral_img_int, test_functions, width, height, image_int, false);
            if (valid) {
                printf("INTEGRAL IMAGE INT VALIDATION:             \033[0;32mSUCCESS!\033[0m\n");
            } else {
                printf("INTEGRAL IMAGE INT VALIDATION:             \033[1;31mFAILED!\033[0m\n");
            }
        }
    }
#endif

#ifdef VALIDATE_INTEGRAL_IMAGE_INT_PADDED
    {
        //if (width + 2 * PADDING_SIZE < 2895) {
        if (width + 2 * PADDING_SIZE <= 2048 + 2 * PADDING_SIZE) {
            std::vector<void (*)(uint8_t *, struct integral_image *)> test_functions;
            //test_functions.push_back(compute_padded_integral_img_int);
            test_functions.push_back(compute_padded_integral_img_simd_early_cast_int);
            //test_functions.push_back(compute_integral_img_simd_original);

            bool valid = validate_integral_image_int(compute_padded_integral_img_int, test_functions, width, height, image_int, true);
            if (valid) {
                printf("INTEGRAL IMAGE INT PADDED VALIDATION:      \033[0;32mSUCCESS!\033[0m\n");
            } else {
                printf("INTEGRAL IMAGE INT PADDED VALIDATION:      \033[1;31mFAILED!\033[0m\n");
            }
        }
    }
#endif

	// Compute integral image
	compute_integral_img(image, iimage);

	// Fast-Hessian
	struct fasthessian* fh = create_fast_hessian(iimage);

	// Create octaves with response layers
	create_response_map(fh);

	// Compute responses for every layer
	compute_response_layers(fh);

#ifdef VALIDATE_COMPUTE_RESPONSE_LAYER
    {
        std::vector<void (*)(struct fasthessian *)> test_functions;
        test_functions.push_back(compute_response_layers_precompute);
        test_functions.push_back(compute_response_layers_blocking);
        test_functions.push_back(compute_response_layers_at_once);
        test_functions.push_back(compute_response_layers_sonic_Dyy);

        // test_functions.push_back(compute_response_layers_blocking_3_3_False);
        // test_functions.push_back(compute_response_layers_blocking_3_7_False);
        // test_functions.push_back(compute_response_layers_blocking_7_3_False);
        // test_functions.push_back(compute_response_layers_blocking_7_7_False);

        bool valid = validate_compute_response_layers(compute_response_layers, test_functions, iimage);
        if (valid) {
            printf("COMPUTE RESPONSE LAYER VALIDATION:         \033[0;32mSUCCESS!\033[0m\n");
        } else {
            printf("COMPUTE RESPONSE LAYER VALIDATION:         \033[1;31mFAILED!\033[0m!\n");
        }
    }
#endif

#ifdef VALIDATE_COMPUTE_RESPONSE_LAYER_PADDED
    {
        std::vector<void (*)(struct response_layer *, struct integral_image *)> test_functions;
        test_functions.push_back(compute_response_layer_unconditional);
        test_functions.push_back(compute_response_layer_sonic_Dyy_unconditional);
        test_functions.push_back(compute_response_layer_sonic_Dyy_unconditional_opt);

        // test_functions.push_back(compute_response_layer_Dyy_laplacian_locality_uncond_opt_flops_invsqr); // is expected to fail for images < 128
        // test_functions.push_back(compute_response_layer_unconditional_strided); // is expected to fail on layers > 4

        bool valid = validate_compute_response_layer_with_padding(compute_response_layer, test_functions, image, width, height);

        if (valid) {
            printf("COMPUTE RESPONSE LAYER PADDED VALIDATION:  \033[0;32mSUCCESS!\033[0m\n");
        } else {
            printf("COMPUTE RESPONSE LAYER PADDED VALIDATION:  \033[1;31mFAILED!\033[0m\n");
        }
    }
#endif

// Getting interest points with non-maximum supression
    std::vector<struct interest_point> interest_points;
    get_interest_points(fh, &interest_points);

#ifdef VALIDATE_GET_INTEREST_POINTS
    {
        std::vector<void (*)(struct fasthessian *, std::vector<struct interest_point> *)> test_functions;
        //test_functions.push_back(get_interest_points);
        test_functions.push_back(get_interest_points_layers);

        bool valid = validate_get_interest_points(get_interest_points, test_functions, fh);
        if (valid) {
            printf("GET INTEREST POINTS VALIDATION:            \033[0;32mSUCCESS!\033[0m\n");
        } else {
            printf("GET INTEREST POINTS VALIDATION:            \033[1;31mFAILED!\033[0m\n");
        }
    }

#endif

    // Getting M-SURF descriptors for each interest point
	get_msurf_descriptors(iimage, &interest_points);

#ifdef VALIDATE_GET_MSURF_DESCRIPTORS
    {
        std::vector<void (*)(struct integral_image *, struct interest_point *)> test_functions;
        // test_functions.push_back(get_msurf_descriptor);
        test_functions.push_back(get_msurf_descriptor_improved);
        test_functions.push_back(get_msurf_descriptor_inlined);
        test_functions.push_back(get_msurf_descriptor_inlinedHaarWavelets);
        test_functions.push_back(get_msurf_descriptor_precompute_gauss_case);
        test_functions.push_back(get_msurf_descriptor_precompute_gauss_array);
        test_functions.push_back(get_msurf_descriptor_pecompute_haar);
        test_functions.push_back(get_msurf_descriptor_gauss_pecompute_haar_unroll);
        test_functions.push_back(get_msurf_descriptor_rounding);
        test_functions.push_back(get_msurf_descriptor_rounding_unroll_2_24_True_winner);
        test_functions.push_back(get_msurf_descriptor_simd);
        test_functions.push_back(get_msurf_descriptor_simd_2_24);

        // test_functions.push_back(get_msurf_descriptor_haar_unroll_4_1_False);
        // test_functions.push_back(get_msurf_descriptor_haar_unroll_1_4_True);


        bool valid = validate_get_msurf_descriptors(get_msurf_descriptor, test_functions, iimage, &interest_points);
        if (valid) {
            printf("MSURF DESCRIPTOR VALIDATION:               \033[0;32mSUCCESS!\033[0m\n");
        } else {
            printf("MSURF DESCRIPTOR VALIDATION:               \033[1;31mFAILED!\033[0m\n");
        }
    }
#endif


	// Write results to file
    FILE * fp = fopen(argv[2],"w");
    printf("%d %d %d\n", iimage->width, iimage->height, channels);
	for (size_t i=0; i<interest_points.size(); ++i) {
        fprintf(fp, "%f %f %f %d ", interest_points[i].x, interest_points[i].y, interest_points[i].scale, interest_points[i].laplacian);
        for(size_t j = 0; j < 64; j++) {
            fprintf(fp, "%f ", interest_points[i].descriptor[j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

	// Free memory
#if defined(VALIDATE_INTEGRAL_IMAGE_INT) || defined(VALIDATE_INTEGRAL_IMAGE_INT_PADDED)
    stbi_image_free(image_int);
#endif
	stbi_image_free(image); // possibly move this to create_integral_img
	free(iimage->padded_data);
	free(iimage);
	for (int i = 0; i < NUM_LAYERS; ++i) {
        free(fh->response_map[i]->response);
        free(fh->response_map[i]->laplacian);
		free(fh->response_map[i]);
	}
	free(fh);

    extern float* haarResponseX;
    extern float* haarResponseY;

    aligned_free(haarResponseX);
    aligned_free(haarResponseY);

	return 0;
}
