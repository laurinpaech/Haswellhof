#include "lodepng.h"
#include <stdio.h>
#include <stdlib.h>

/** https://github.com/lvandeve/lodepng */
unsigned char* decodeOneStep(const char* filename) {
	unsigned error;
	unsigned char* image = 0;
	unsigned width, height;

	/*error = lodepng_decode32_file(&image, &width, &height, filename);*/
	error = lodepng_decode_file(&image, &width, &height, filename, LCT_GREY, 8);

	if (error)  {
		printf("error %u: %s\n", error, lodepng_error_text(error));
	}

	printf("width: %d; height: %d\n%u %u %u %u\n", width, height, image[0], image[1], image[2], image[3]);
	/*convert image to grayscale*/
	return image;
}

int main(int argc, char const *argv[])
{
	const char* filename = argc > 1 ? argv[1] : "test.png";

	/* image in memory: RGBA RGBA RGBA RGBA RGBA */
	unsigned char* image = decodeOneStep(filename);

	/* use image here */

	free(image);

	return 0;
}
