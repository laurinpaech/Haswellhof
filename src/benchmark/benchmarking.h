#pragma once

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "descriptor.h"
#include "fasthessian.h"
#include "helper.h"
#include "integral_image.h"

struct benchmark_data {
    benchmark_data()
        : width(-1),
          height(-1),
          num_interest_points(0),
          num_flops(-1),
          avg_cycles(0),
          min_cycles(0),
          max_cycles(0),
          flops_per_cycle(0.0) {}

    benchmark_data(char *i_name, int w, int h, char *f_name, int n_ipoints, long n_flops)
        : width(w),
          height(h),
          num_interest_points(n_ipoints),
          num_flops(n_flops),
          avg_cycles(0),
          min_cycles(0),
          max_cycles(0),
          flops_per_cycle(0.0) {
        // Copying image and function name strings into struct
        strncpy(image_name, i_name, 256);
        strncpy(function_name, f_name, 256);
    }

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

/*
inline void set_benchmark_data(char *image_name, int width, int height, char *function_name, int num_interest_points,
long num_flops, struct benchmark_data &data) {

    // Setting all values for benchmark data structure
    strcpy(data.image_name, image_name);
    data.width = width;
    data.height = height;
    strcpy(data.function_name, function_name);
    data.num_interest_points = num_interest_points;
    data.num_flops = num_flops;
    data.avg_cycles = 0;
    data.min_cycles = 0;
    data.max_cycles = 0;
    data.flops_per_cycle = 0.0;

}
*/

void bench_compute_integral_img(const std::vector<void (*)(float *, int, int, float *)> &functions, float *gray_image,
                                std::vector<struct benchmark_data> &data);

// Times the function compute_integral_img from integral_image. Stores the average, minimum and maximum number of cycles
// and the flops per cycle in benchmark_data. The number of flops for compute_integral_img must be set benchmark_data.
// The height and width of the image must be set in benchmark_data.
void perf_compute_integral_img(void (*function)(float *, int, int, float *), float *gray_image,
                               struct benchmark_data &data);

/*
void bench_compute_response_layer(const std::vector<void (*)(struct response_layer *, struct integral_image *)>
&functions, struct response_layer *layer, struct integral_image *iimage, std::vector<struct benchmark_data> &data);
*/

// Times the function compute_response_layer from fasthessian. Stores the average, minimum and maximum number of cycles
// and the flops per cycle in benchmark_data. The number of flops for compute_response_layer must be set benchmark_data.
void perf_compute_response_layer(void (*function)(struct response_layer *, struct integral_image *),
                                 struct response_layer *layer, struct integral_image *iimage,
                                 struct benchmark_data &data);

void perf_compute_response_all_layers(void (*function)(struct response_layer *, struct integral_image *),
                                      struct fasthessian *fh, struct benchmark_data &data);

void bench_compute_response_layer(
    const std::vector<void (*)(struct response_layer *, struct integral_image *)> &functions,
    struct integral_image *iimage, std::vector<struct benchmark_data> &data);

void bench_get_interest_points(
    const std::vector<void (*)(struct fasthessian *, std::vector<struct interest_point> *)> &functions,
    struct fasthessian *fh, std::vector<struct benchmark_data> &data);

// Times the function get_interest_points from fasthessian. Stores the average, minimum and maximum number of cycles and
// the flops per cycle in benchmark_data. The number of flops for get_interest_points must be set benchmark_data. The
// number of interest points that will be found, must be set in benchmark_data.
void perf_get_interest_points(void (*function)(struct fasthessian *, std::vector<struct interest_point> *),
                              struct fasthessian *fh, struct benchmark_data &data);

// Calls the timing function for get_interest_points. The number of times the timing should be conducted can be
// specified in this function. The number of flops for interpolate_step must be set in benchmark_data.
void bench_interpolate_step(const std::vector<void (*)(int, int, struct response_layer *, struct response_layer *,
                                                       struct response_layer *, float[3])> &functions,
                            struct fasthessian *fh, std::vector<struct benchmark_data> &data);

// Times the function interpolate_step from fasthessian. Stores the average, minimum and maximum number of cycles and
// the flops per cycle in benchmark_data. The number of flops for interpolate_step must be set in benchmark_data.
void perf_interpolate_step(void (*function)(int, int, struct response_layer *, struct response_layer *,
                                            struct response_layer *, float[3]),
                           int row, int col, struct response_layer *top, struct response_layer *middle,
                           struct response_layer *bottom, struct benchmark_data &data);

// Calls the timing function for get_descriptor. The number of times the timing should be conducted can be specified in
// this function. The number of flops for get_descriptor must be set in benchmark_data.
void bench_get_descriptor(
    const std::vector<void (*)(struct integral_image *, struct interest_point *, float *)> &functions,
    struct integral_image *iimage, std::vector<struct interest_point> *interest_points, float *GW,
    std::vector<struct benchmark_data> &data);

// Times the function get_descriptor from descriptor. Stores the average, minimum and maximum number of cycles and the
// flops per cycle in benchmark_data. The number of flops for get_descriptor must be set in benchmark_data.
void perf_get_descriptor(void (*function)(struct integral_image *, struct interest_point *, float *),
                         struct integral_image *iimage, struct interest_point *ipoint, float *GW,
                         struct benchmark_data &data);

// Calls the timing function for get_msurf_descriptor. The number of times the timing should be conducted can be
// specified in this function. The number of flops for get_msurf_descriptor must be set in benchmark_data.
void bench_get_msurf_descriptor(
    const std::vector<void (*)(struct integral_image *, struct interest_point *)> &functions,
    struct integral_image *iimage, std::vector<struct interest_point> *interest_points,
    std::vector<struct benchmark_data> &data);

// Times the function get_msurf_descriptor from descriptor. Stores the average, minimum and maximum number of cycles and
// the flops per cycle in benchmark_data. The number of flops for get_msurf_descriptor must be set in benchmark_data.
void perf_get_msurf_descriptor(void (*function)(struct integral_image *, struct interest_point *),
                               struct integral_image *iimage, struct interest_point *ipoint,
                               struct benchmark_data &data);
