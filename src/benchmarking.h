#pragma once

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "descriptor.h"
#include "fasthessian.h"
#include "integral_image.h"

// https://stackoverflow.com/questions/3437404/min-and-max-in-c
#define MIN(a, b)               \
    ({                          \
        __typeof__(a) _a = (a); \
        __typeof__(b) _b = (b); \
        _a < _b ? _a : _b;      \
    })

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

static inline struct benchmark_data *initialise_benchmark_data(char *image_name, int width, int height,
                                                               char *function_name, int num_interest_points,
                                                               long num_flops) {
    struct benchmark_data *benchmark_data = (struct benchmark_data *)calloc(1, sizeof(struct benchmark_data));
    benchmark_data->width = width;
    benchmark_data->height = height;
    benchmark_data->num_interest_points = num_interest_points;
    benchmark_data->num_flops = num_flops;
    strcpy(benchmark_data->image_name, image_name);
    strcpy(benchmark_data->function_name, function_name);
    return benchmark_data;
}

// Times the function compute_integral_img from integral_image. Stores the average, minimum and maximum number of cycles
// and the flops per cycle in benchmark_data. The number of flops for compute_integral_img must be set benchmark_data.
// The height and width of the image must be set in benchmark_data.
void perf_test_integral_img(void (*function)(float *, int, int, float *), float *gray_image,
                            struct benchmark_data *data);

// Times the function compute_response_layer from fasthessian. Stores the average, minimum and maximum number of cycles
// and the flops per cycle in benchmark_data. The number of flops for compute_response_layer must be set benchmark_data.
void perf_test_compute_response_layer(void (*function)(struct response_layer *, struct integral_image *),
                                      struct response_layer *layer, struct integral_image *iimage,
                                      struct benchmark_data *data);

// Times the function get_interest_points from fasthessian. Stores the average, minimum and maximum number of cycles and
// the flops per cycle in benchmark_data. The number of flops for get_interest_points must be set benchmark_data. The
// number of interest points that will be found, must be set in benchmark_data.
void perf_test_get_interest_points(void (*function)(struct fasthessian *, std::vector<struct interest_point> *),
                                   struct fasthessian *fh, struct benchmark_data *data);

// Calls the timing function for get_interest_points. The number of times the timing should be conducted can be
// specified in this function. The number of flops for interpolate_step must be set in benchmark_data.
void bench_interpolate_step(struct fasthessian *fh, struct benchmark_data *data);

// Times the function interpolate_step from fasthessian. Stores the average, minimum and maximum number of cycles and
// the flops per cycle in benchmark_data. The number of flops for interpolate_step must be set in benchmark_data.
void perf_interpolate_step(void (*function)(int, int, struct response_layer *, struct response_layer *,
                                            struct response_layer *, float[3]),
                           int row, int col, struct response_layer *top, struct response_layer *middle,
                           struct response_layer *bottom, struct benchmark_data *data);

// Calls the timing function for get_descriptor. The number of times the timing should be conducted can be specified in
// this function. The number of flops for get_descriptor must be set in benchmark_data.
void bench_get_descriptor(struct integral_image *iimage, std::vector<struct interest_point> *interest_points, float *GW,
                          struct benchmark_data *data);

// Times the function get_descriptor from descriptor. Stores the average, minimum and maximum number of cycles and the
// flops per cycle in benchmark_data. The number of flops for get_descriptor must be set in benchmark_data.
void perf_get_descriptor(void (*function)(struct integral_image *, struct interest_point *, float *),
                         struct integral_image *iimage, struct interest_point *ipoint, float *GW,
                         struct benchmark_data *data);
