#pragma once

#include "response_layer.h"

#include <stdlib.h>
#include <stdbool.h>

#define NUM_OCTAVES 4
#define NUM_LAYER 10

struct fasthessian {

    // Integral image
    struct integral_image* iimage;

    // Response stack of determinant of hessian values
    struct response_layer* response_map[NUM_LAYER];

    // Number of Octaves
    int octaves;

    // Number of layers per octave
    int layers;

    // Initial sampling step for Interest Point detection
    int step;

    // Threshold value for responses
    float thresh;

};

// default threshold of hessian response for non-maximum suppression
float threshold = 0.0004f;

// Create Fast-Hessian struct
void createFastHessian();

// Create octaves with response layers
void buildResponseMap();

// Compute responses for layer
void buildResponseLayer(struct response_layer* layer);

// checking if (row, col) is maximum in 3x3x3 neighborhood
bool is_extremum(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom);

// interpolating maximum at (row, col) with quadratic in 3x3x3 neighborhood to get sub-pixel location
// and storing potential interest point if interpolation remains fairly 'close' to pixel location
void interpolate_extremum(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom);

// constructing hessian and negative gaussian to solve 3x3 linear system and get sub-pixel offsets
void interpolate_step(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets[3]);
