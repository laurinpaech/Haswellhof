#pragma once

#include "fasthessian.h"

#include "helper.h"

struct fasthessian_flat {

    // Integral image
    struct integral_image* iimage;

    // Number of Octaves
    int octaves;

    // Number of layers per octave
    int layers;

    // Number of layers in total
    int total_layers;

    // Initial sampling step for Interest Point detection
    int step;

    // Threshold value for hessian response in non-maximum suppression
    float thresh;

    // Response stack of determinant of hessian values
    struct response_layer response_map[NUM_TOTAL_LAYERS];

};

void create_fast_hessian_flat_and_response_map(struct integral_image *iimage, struct fasthessian_flat *fh_flat);

void compute_response_layers_flat(struct fasthessian_flat* fh_flat);

void get_interest_points_flat(struct fasthessian_flat *fh_flat, std::vector<struct interest_point> *interest_points);

