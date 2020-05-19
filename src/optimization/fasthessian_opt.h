#pragma once

#include "fasthessian.h"


void super_sonic_Dyy(struct response_layer *layer, struct integral_image *iimage);

void blue_lines_Dyy(struct response_layer *layer, struct integral_image *iimage);

void double_blue_lines_Dyy(struct response_layer *layer, struct integral_image *iimage);

void get_interest_points_layers(struct fasthessian *fh, std::vector<struct interest_point> *interest_points);

void interpolate_step_gauss(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets[3]);

void compute_response_layer_Dyy_leftcorner(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layer_Dyy_top(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layer_Dyy_top_mid(struct response_layer* layer, struct integral_image* iimage);

void compute_response_map_Dyy(struct fasthessian *fh);

void compute_response_layer_Dyy(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_at_once(struct fasthessian* fh, struct integral_image *iimage);
