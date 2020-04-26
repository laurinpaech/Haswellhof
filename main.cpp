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
	int width, height, channels;

	// Load image
	stbi_ldr_to_hdr_gamma(1.0f);
	float* image = stbi_loadf("images/img1.png", &width, &height, &channels, STBI_grey);

    if(!image) {
		printf("Could not open or find image\n");
        return -1;
    }

	// Calculate integral image
	struct integral_image* iimage = create_integral_img(image, width, height);

	// Fast-Hessian
	struct fasthessian* fh = create_fast_hessian(iimage);

	// Create octaves with response layers
	create_response_map(fh);

	// Compute responses for every layer
	for (size_t i = 0; i < fh->total_layers; i++) {
		compute_response_layer(fh->response_map[i], iimage);
	}

	// Getting interest points with non-maximum supression
	std::vector<struct interest_point> interest_points;
	get_interest_points(fh, &interest_points);

	// Descriptor stuff
    static float* GW = get_gaussian(3.3);
	for (int i=0; i<interest_points.size(); i++)
        get_descriptor(iimage, &interest_points[i], GW);

    free(GW);

	// Write results to file
    FILE * fp = fopen("desc.txt","w");
	for (int i=0; i<interest_points.size(); i++) {
        fprintf(fp, "%f %f ", interest_points[i].x, interest_points[i].y);
        // printf("%f ", BLOB_ORIENTATION); // TODO: how to get this value
        for(int j = 0; j < 64; j++) {
            fprintf(fp, "%f ", interest_points[i].descriptor[j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

	// Free memory
	stbi_image_free(image); // possibly move this to create_integral_img
	free(iimage->data);
	free(iimage);
	for (size_t i = 0; i < NUM_LAYERS; i++) {
		free(fh->response_map[i]);
	}
	free(fh);

	return 0;
}
