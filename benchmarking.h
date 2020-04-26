#pragma once

#include <stdlib.h>
#include "integral_image.h"
#include "fasthessian.h"

//times the function create_integral_img from integral_image and returns the flops per cycle
double perf_test_integral_img(struct integral_image* (*function)(float*, int, int), float* gray_image, int width, int height,int flops);
//times the function create_integral_img from integral_image and saves the output to a file
void bench_integral_img(float* image, int width, int height);


//times the function compute_response_layer from fasthessian and returns the flops per cycle
double perf_test_compute_response_layer(void (*function)(struct response_layer*, struct integral_image*), struct response_layer* layer,struct integral_image* iimage, int flops);
//times the function compute_response_layer from fasthessian and saves the output to a file
void bench_compute_response_layer(struct response_layer* layer, struct integral_image* iimage, int width, int height);


//saves a file containing the flops_per_cycle in the folder "benchmarking" under the specified file_name + the current datetime
void save_performance_file(double flops_per_cycle, const char *file_name);
//concatenates two given strings
char* concat(const char *s1, const char *s2);

