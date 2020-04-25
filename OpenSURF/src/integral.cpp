/***********************************************************
*  --- OpenSURF ---                                       *
*  This library is distributed under the GNU GPL. Please   *
*  use the contact form at http://www.chrisevansdev.com    *
*  for more information.                                   *
*                                                          *
*  C. Evans, Research Into Robust Visual Features,         *
*  MSc University of Bristol, 2008.                        *
*                                                          *
************************************************************/
// has to be defined before stb includes
#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"
#include "utils.h"

#include "integral.h"

//! Computes the integral image of image img.  Assumes source image to be a
//! 32-bit floating point.  Returns IplImage of 32-bit float form.
struct integral_image* Integral()
{
  int height, width, channels;

  // Load image
  stbi_ldr_to_hdr_gamma(1.0f);
  float* data = stbi_loadf("imgs/test.png", &width, &height, &channels, STBI_grey);

  if(!data) {
      printf("Could not open or find image\n");
      return NULL;
  }

  struct integral_image* integral_img = (struct integral_image *) malloc(sizeof(struct integral_image));
  integral_img->height = height;
  integral_img->width = width;
  float *i_data = (float*) malloc(width * height * sizeof(float));

  int step = width;

  // TESTING
  // for (int i = 0; i < width; i++) {
  //     for (int j = 0; j < height; j++) {
  //         printf("%i, %i - %f\n", i, j, data[i*width+j]);
  //     }
  // }

  // first row only
  float rs = 0.0f;

  for(int j=0; j < width; j++)
  {
    rs += data[j];
    i_data[j] = rs;
  }

  // remaining cells are sum above and to the left
  for(int i=1; i < height; ++i)
  {
    rs = 0.0f;
    for(int j=0; j < width; ++j)
    {
      rs += data[i*step+j];
      i_data[i*step+j] = rs + i_data[(i-1)*step+j];
      // printf("%i, %i - %f\n", i, j, i_data[i*step+j]);
    }
  }

  integral_img->imageData = i_data;

  // release the gray image
  stbi_image_free(data);

  // return the integral image
  return integral_img;
}
