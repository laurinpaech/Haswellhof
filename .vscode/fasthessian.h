#pragma once

#include "response_layer.h"

#include <stdbool.h>

// default threshold of hessian response for non-maximum suppression
float threshold = 0.0004f;

// checking if (row, col) is maximum in 3x3x3 neighborhood
bool is_extremum(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom);

// interpolating maximum at (row, col) with quadratic in 3x3x3 neighborhood to get sub-pixel location
// and storing potential interest point if interpolation remains fairly 'close' to pixel location
void interpolate_extremum(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom);

// constructing hessian and negative gaussian to solve 3x3 linear system and get sub-pixel offsets 
void interpolate_step(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets[3]);
