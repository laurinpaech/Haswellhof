#pragma once

#include "fasthessian.h"

void get_interest_points_layers(struct fasthessian *fh, std::vector<struct interest_point> *interest_points);

void interpolate_step_gauss(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets[3]);
