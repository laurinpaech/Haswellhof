#pragma once
#include "benchmarking.h"


void bench_interpolate_step(struct fasthessian *fh, std::vector<struct interest_point> *interest_points, struct benchmark_data* data);


//times the function create_integral_img from integral_image and returns the flops per cycle
void perf_interpolate_step(void (*function)(int, int, struct response_layer*, struct response_layer *, struct response_layer*, float[3]),int row, int col, struct response_layer *top, struct response_layer *middle, 
struct response_layer *bottom, struct benchmark_data* data);

