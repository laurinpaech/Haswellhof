#pragma once

#include <stdlib.h>

#include "integral_image.h"
#include "fasthessian.h"

//times the function create_integral_img from integral_image and returns the flops per cycle
double perf_test_integral_img(struct integral_image* (*function)(float*, int, int), float* gray_image, int width, int height,int flops);
//times the function create_integral_img from integral_image and saves the output to a file
void bench_integral_img(float* image, int width, int height);


//times the function compute_response_layer from fasthessian and returns the flops per cycle
double perf_test_compute_response_layer(void (*function)(struct response_layer*, struct integral_image*), struct response_layer* layer,struct integral_image* iimage, long flops);
//times the function compute_response_layer from fasthessian and saves the output to a file
void bench_compute_response_layer(struct response_layer* layer, struct integral_image* iimage, int width, int height);

//times the function get_interest_points from fasthessian and returns the flops per cycle
double perf_test_get_interest_points(void (*function)(struct fasthessian*, std::vector<struct interest_point>*), struct fasthessian *fh, long flops);
//times the function get_interest_points from fasthessian and saves the output to a file
void bench_get_interest_points(struct fasthessian *fh, int width, int height, int initial_step);


//times the function interpolate_step from fasthessian and returns the flops per cycle
double perf_test_interpolate_step(void (*function)(int, int, struct response_layer*, struct response_layer*, struct response_layer*, float), int row, int col, 
    struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets, int flops);
//times the function interpolate_step from fasthessian and saves the output to a file
void bench_interpolate_step(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets);



//saves a file containing the flops_per_cycle in the folder "benchmarking" under the specified file_name + the current datetime
void save_performance_file(double flops_per_cycle, const char *file_name);
//concatenates two given strings
char* concat(const char *s1, const char *s2);

