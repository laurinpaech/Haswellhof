// has to be defined before stb includes
#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"
#include "integral_image.h"
#include "fasthessian.h"
#include "interest_point.h"
#include "descriptor.h"

#include "integral_image_opt.h"
#include "integral_image_simd.h"
#include "fasthessian_opt.h"
#include "descriptor_opt.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>

#define CREATE_PADDED_INTEGRAL_IMAGE 0

//#define BASE_COMPUTE_INTEGRAL_IMAGE
//#define BASE_COMPUTE_RESPONSE_LAYERS
//#define BASE_GET_INTEREST_POINTS
//#define BASE_GET_MSURF_DESCRIPTORS

#define FULL_OPTIMIZATION_COMPUTE_INTEGRAL_IMAGE
#define FULL_OPTIMIZATION_COMPUTE_RESPONSE_LAYERS
#define FULL_OPTIMIZATION_GET_INTEREST_POINTS
#define FULL_OPTIMIZATION_GET_MSURF_DESCRIPTORS

#define FILE_OUTPUT

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

#if !CREATE_PADDED_INTEGRAL_IMAGE
    // Create integral image
    struct integral_image* iimage = create_integral_img(width, height);
#else
	// Create integral image
	struct integral_image* iimage = create_padded_integral_img(width, height);
#endif

#ifdef BASE_COMPUTE_INTEGRAL_IMAGE
	// Compute integral image
    compute_integral_img(image, iimage);
#endif

#ifdef FULL_OPTIMIZATION_COMPUTE_INTEGRAL_IMAGE
	// Compute integral image
	compute_integral_img_faster_alg(image, iimage);
#endif

    // Fast-Hessian
    struct fasthessian* fh = create_fast_hessian(iimage);

    // Create octaves with response layers
    create_response_map(fh);

#ifdef BASE_COMPUTE_RESPONSE_LAYERS
	// Compute responses for every layer
	compute_response_layers(fh);
#endif

#ifdef FULL_OPTIMIZATION_COMPUTE_RESPONSE_LAYERS
	// Compute responses for every layer
	compute_response_layers_sonic_Dyy(fh);
#endif

	std::vector<struct interest_point> interest_points;

#ifdef BASE_GET_INTEREST_POINTS
	// Getting interest points with non-maximum supression
    get_interest_points(fh, &interest_points);
#endif

#ifdef FULL_OPTIMIZATION_GET_INTEREST_POINTS
	// Getting interest points with non-maximum supression
	get_interest_points(fh, &interest_points);
#endif

#ifdef BASE_GET_MSURF_DESCRIPTORS
    // Getting M-SURF descriptors for each interest point
	get_msurf_descriptors(iimage, &interest_points);
#endif

#ifdef FULL_OPTIMIZATION_GET_MSURF_DESCRIPTORS
	// Getting M-SURF descriptors for each interest point
	get_msurf_descriptors_haar_unroll_2_24_True_winner(iimage, &interest_points);
#endif


#ifdef FILE_OUTPUT
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
#endif
    

    // Free memory
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
