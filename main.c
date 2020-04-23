#include <stdio.h>
#include <stdlib.h>

// has to be defined before stb includes
#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"

int main(int argc, char const *argv[])
{
	int width, height, channels;

	// image is gamma corrected by default:
	// stackoverflow.com/questions/59774647/loading-image-with-stb-image-library-as-float-outputs-wrong-values

	// load image
	// stbi_ldr_to_hdr_gamma(1.0f)
	float* image = stbi_loadf("test.png", &width, &height, &channels, STBI_grey);

    if(!image) {
		printf("Could not open or find image\n");
        return -1;
    }

	stbi_image_free(image);

	return 0;
}
