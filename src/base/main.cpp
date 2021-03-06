// has to be defined before stb includes
#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"
#include "integral_image.h"
#include "fasthessian.h"
#include "interest_point.h"
#include "descriptor.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>

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

	// Create integral image
	struct integral_image* iimage = create_integral_img(width, height);
	// Compute integral image
	compute_integral_img(image, iimage);

	// Fast-Hessian
	struct fasthessian* fh = create_fast_hessian(iimage);

	// Create octaves with response layers
	create_response_map(fh);

	// Compute responses for every layer
	compute_response_layers(fh);

	// Getting interest points with non-maximum supression
	std::vector<struct interest_point> interest_points;
	get_interest_points(fh, &interest_points);

    // printf("interest points: %lu\n", interest_points.size());

	// Getting M-SURF descriptors for each interest point
	get_msurf_descriptors(iimage, &interest_points);

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
	stbi_image_free(image); // possibly move this to create_integral_img
	free(iimage->padded_data);
	free(iimage);
	for (int i = 0; i < NUM_LAYERS; ++i) {
		free(fh->response_map[i]->response);
		free(fh->response_map[i]->laplacian);
		free(fh->response_map[i]);
	}
	free(fh);

	return 0;
}
