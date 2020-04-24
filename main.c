// has to be defined before stb includes
#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"
#include "integral_image.h"
//#include "fasthessian.h"

#include <stdio.h>
#include <stdlib.h>

#include "tsc_x86.h"


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

	int start, end;
	double cycles = 0.;
	FILE *fptr_int_img;
	// Creates a file "fptr_int_img" 
    // with file acccess as write mode 
	fptr_int_img = fopen("benchmarking_files/integral_img.txt","w");

	// Calculate integral image
	start = start_tsc();
	struct integral_image* iimage = create_integral_img(image, width, height);
	end = stop_tsc(start);
	cycles = (double)end;

    fprintf(fptr_int_img,"%lf \n",cycles);

	// closes the file pointed by fptr_int_img 
    fclose(fptr_int_img); 

	// Fast-Hessian
	/*struct fasthessian* fh = create_fast_hessian(iimage);

	// Create octaves with response layers
	create_response_map(fh);

	// Compute responses for every layer
	for (size_t i = 0; i < fh->layers; i++) {
		compute_response_layer(fh->response_map[i], iimage);
	}*/

	// Non-maximum supression interest points
	// TODO

	// Detector stuff
	// TODO

	// Post-processing
	// TODO

	// Free memory
	stbi_image_free(image); // possibly move this to create_integral_img
	free(iimage->data);
	free(iimage);
	/*for (size_t i = 0; i < NUM_LAYERS; i++) {
		free(fh->response_map[i]);
	}
	free(fh);
	*/

	return 0;
}

