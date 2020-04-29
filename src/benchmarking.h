#pragma once

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "integral_image.h"
#include "fasthessian.h"
#include "descriptor.h"


struct benchmark_data {

    // Image name
    char image_name[256]; 

    // Image width
    int width;

    // Image height
    int height;

    // Name of fucntion that has been benchmarked
    char function_name[256]; 

    // Number of interest points
    int num_interest_points;

    // Number of flops -1 if no flops
    long num_flops;

    // Number of average cycles
    uint64_t avg_cycles;

    // Minimum amount of cycles needed
    uint64_t min_cycles;

    // Maximum amount of cycles needed
    uint64_t max_cycles;

    // Number of flops per cycle -1 if no flops 
    double flops_per_cycle;   

};

static inline struct benchmark_data* initialise_benchmark_data(char*image_name, int width, int height, char* function_name, int num_interest_points, long num_flops){
    struct benchmark_data* benchmark_data=(struct benchmark_data* )malloc(sizeof(struct benchmark_data));
    benchmark_data->width=width; 
    benchmark_data->height=height;
    benchmark_data->num_interest_points = num_interest_points;
    benchmark_data->num_flops = num_flops;
    strcpy(benchmark_data->image_name, image_name );
    strcpy(benchmark_data->function_name, function_name);
    return benchmark_data;
}


//times the function create_integral_img from integral_image and returns the flops per cycle
void perf_test_integral_img(struct integral_image* (*function)(float*, int, int), float* gray_image, struct benchmark_data* data);

//times the function compute_response_layer from fasthessian and returns the flops per cycle
void perf_test_compute_response_layer(void (*function)(struct response_layer*, struct integral_image*), struct response_layer* layer,struct integral_image* iimage, 
struct benchmark_data* data);

//times the function get_interest_points from fasthessian and returns the flops per cycle
void perf_test_get_interest_points(void (*function)(struct fasthessian*, std::vector<struct interest_point>*), struct fasthessian *fh, struct benchmark_data* data);

void bench_interpolate_step(struct fasthessian *fh, struct benchmark_data* data);


//times the function create_integral_img from integral_image and returns the flops per cycle
void perf_interpolate_step(void (*function)(int, int, struct response_layer*, struct response_layer *, struct response_layer*, float[3]),int row, int col, struct response_layer *top, struct response_layer *middle, 
struct response_layer *bottom, struct benchmark_data* data);

void bench_get_descriptor(struct integral_image* iimage, std::vector<struct interest_point> *interest_points, float* GW, struct benchmark_data* data);

//times the function create_integral_img from integral_image and returns the flops per cycle
void perf_get_descriptor(void (*function)(struct integral_image*, struct interest_point*, float*),struct integral_image* iimage, struct interest_point* ipoint, float* GW, struct benchmark_data* data);










