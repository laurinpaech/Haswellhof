#pragma once

#include "fasthessian.h"
#include "integral_image_opt.h"


void height_greater_border_width_greater_double_lobe_Dyy(struct response_layer *layer, struct integral_image *iimage);

void height_greater_border_width_less_double_lobe_Dyy(struct response_layer *layer, struct integral_image *iimage);

void get_interest_points_layers(struct fasthessian *fh, std::vector<struct interest_point> *interest_points);

void interpolate_step_gauss(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets[3]);

void compute_response_layer_Dyy_leftcorner(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layer_precompute(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_precompute(struct fasthessian* fh);

void compute_response_layer_blocking(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_blocking(struct fasthessian* fh);

void compute_response_layer_Dyy_top(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layer_Dyy_top_mid(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_Dyy(struct fasthessian *fh);

void compute_response_layer_Dyy(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layer_sonic_Dyy(struct response_layer *layer, struct integral_image *iimage);

void compute_response_layer_Dyy_laplacian(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_sonic_Dyy(struct fasthessian *fh);

void compute_response_layers_Dyy_laplacian(struct fasthessian *fh);

void compute_response_layers_Dyy_laplacian_localityloops(struct fasthessian *fh);

void compute_response_layer_Dyy_laplacian_localityloops(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_at_once(struct fasthessian* fh);

void compute_response_layers_unconditional(struct fasthessian* fh);

void compute_response_layer_unconditional(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_unconditional_strided(struct fasthessian* fh);

void compute_response_layer_unconditional_strided(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_sonic_Dyy_unconditional(struct fasthessian *fh);

void compute_response_layer_sonic_Dyy_unconditional(struct response_layer *layer, struct integral_image *iimage);

void compute_response_layer_Dyy_laplacian_localityloops_unconditional(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_Dyy_laplacian_localityloops_unconditional(struct fasthessian* fh);
