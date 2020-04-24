// has to be defined before stb includes
#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"
#include "integral_image.h"
#include "fasthessian.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char const *argv[])
{
	int width, height, channels;

	// Load image
	// stbi_ldr_to_hdr_gamma(1.0f)
	float* image = stbi_loadf("test.png", &width, &height, &channels, STBI_grey);

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
	for (size_t i = 0; i < fh->layers; i++) {
		compute_response_layer(fh->response_map[i], iimage);
	}

	// Non-maximum supression interest points
	// TODO

	// Descriptor stuff
	// for (ipoint in interest_points)
    //     get_descriptor(iimage, point);

	// Write results to file
    // dummy code:
    // FILE * fp = fopen("desc.txt","w");
	// for (ipoint in interest_points) {
    //     printf("%f %f ", ipoint->x, ipoint->y);
    //     // printf("%f ", BLOB_ORIENTATION); 
    //     for(int i = 0; i < 64; i++) {
    //         fprintf(fp, "%f ", ipoint->descriptor[i]);
    //     }
    //     printf("/n");
    // }
    // fclose(fp);

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
