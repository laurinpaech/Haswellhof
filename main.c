#include <stdio.h>
#include <stdlib.h>

// has to be defined before stb includes
#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"
#include "integral_image.h"
#include "fasthessian.h"

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
	// TODO
	struct integral_image* iimage = Integral(image, width, height);

	// Fast-Hessian
	struct fasthessian* fh = createFastHessian(iimage);

	// Non-maximum supression interest points
	// TODO

	// Detector stuff
	// TODO

	// Post-processing
	// TODO

	stbi_image_free(image);
	free(iimage->data);
	free(iimage);

	return 0;
}
