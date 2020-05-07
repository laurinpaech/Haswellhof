#pragma once

#include "fasthessian.h"

void interpolate_step_gauss(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets[3]);

void compute_response_layers_at_once(struct fasthessian* fh, struct integral_image *iimage);
