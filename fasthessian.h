#pragma once

#include "response_layer.h"

#include <stdlib.h>
#include <stdbool.h>

// #define NUM_OCTAVES 4
// #define NUM_LAYER 10
// #define INITIAL_STEP 2
// #define THRESHOLD 0.0004f

static const int NUM_OCTAVES = 3;
static const int NUM_LAYER = 4;
static const int NUM_LAYERS = 8;
static const float THRESHOLD = 0.0004f;     // default threshold of hessian response for non-maximum suppression
static const int INITIAL_STEP = 2;

struct fasthessian {

    // Integral image
    struct integral_image* iimage;

    // Response stack of determinant of hessian values
    struct response_layer* response_map[NUM_LAYERS];

    // Number of Octaves
    int octaves;

    // Number of layers per octave
    int layers;

    // Initial sampling step for Interest Point detection
    int step;

    // Threshold value for hessian response in non-maximum suppression
    float thresh;

};

// Create Fast-Hessian struct
struct fasthessian* create_fast_hessian(struct integral_image *iimage);

// Create octaves with response layers
void create_response_map(struct fasthessian* fh);

// Compute responses for layer
void compute_response_layer(struct response_layer* layer, struct integral_image *iimage);

struct response_layer* initialise_response_layer(int filter_size, int width, int height, int init_step);

// checking if (row, col) is maximum in 3x3x3 neighborhood
bool is_extremum(struct fasthessian *fh, int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom);

// interpolating maximum at (row, col) with quadratic in 3x3x3 neighborhood to get sub-pixel location
// and storing potential interest point if interpolation remains fairly 'close' to pixel location
void interpolate_extremum(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom);

// constructing hessian and negative gaussian to solve 3x3 linear system and get sub-pixel offsets
void interpolate_step(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets[3]);

