#include "fasthessian_opt.h"
#include "benchmarking.h"

#include "integral_image_opt.h"

#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <vector>

#include "tsc_x86.h"

#define CYCLES_REQUIRED 1e8
#define REP 10
#define EPS (1e-3)
// Determines whether the warm up for the timing functions will be conducted.
#define WARM_UP


void bench_compute_integral_img(const std::vector<void (*)(float *, struct integral_image *)> &functions,
                                float *gray_image, std::vector<struct benchmark_data> &data, bool is_padded) {

    assert(functions.size() == data.size());

    for (int j = 0; j < functions.size(); ++j) {
        perf_compute_integral_img(functions[j], gray_image, data[j], is_padded);
    }

}

// Times the function compute_integral_img.
// Stores the average, minimum and maximum number of cycles and the flops per cycle in benchmark_data.
// The number of flops for compute_integral_img must be set in benchmark_data.
// The height and width of the image must be set in benchmark_data.
void perf_compute_integral_img(void (*function)(float *, struct integral_image *), float *gray_image,
                               struct benchmark_data &data, bool is_padded) {
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    uint64_t start, end;

    // Checking if integral image should be padded and initializing it
    struct integral_image *iimage;
    if (!is_padded) {
        iimage = create_integral_img(data.width, data.height);
    } else {
        iimage = create_padded_integral_img(data.width, data.height);
    }

#ifdef WARM_UP
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            (*function)(gray_image, iimage);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);
#endif

    std::vector<double> cycleslist;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            (*function)(gray_image, iimage);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cycleslist.push_back(cycles);
    }
    total_cycles /= REP;

    free(iimage->padded_data);
    free(iimage);

    cycles = total_cycles;  // cyclesList.front();
    double flops_per_cycle = round((100.0 * data.num_flops) / cycles) / 100.0;
    std::sort(cycleslist.begin(), cycleslist.end());
    data.avg_cycles = (uint64_t)cycles;
    data.min_cycles = (uint64_t)cycleslist.front();
    data.max_cycles = (uint64_t)cycleslist.back();
    data.flops_per_cycle = flops_per_cycle;
}

void bench_compute_integral_img_int(const std::vector<void (*)(uint8_t *, struct integral_image *)> &functions,
                                    uint8_t *gray_image, std::vector<struct benchmark_data> &data, bool is_padded) {
    assert(functions.size() == data.size());

    for (int j = 0; j < functions.size(); ++j) {
        perf_compute_integral_img_int(functions[j], gray_image, data[j], is_padded);
    }
}

// Times the function compute_integral_img.
// Stores the average, minimum and maximum number of cycles and the flops per cycle in benchmark_data.
// The number of flops for compute_integral_img must be set in benchmark_data.
// The height and width of the image must be set in benchmark_data.
void perf_compute_integral_img_int(void (*function)(uint8_t *, struct integral_image *), uint8_t *gray_image,
                                   struct benchmark_data &data, bool is_padded) {
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    uint64_t start, end;

    // Checking if integral image should be padded and initializing it
    struct integral_image *iimage;
    if (!is_padded) {
        iimage = create_integral_img(data.width, data.height);
    } else {
        iimage = create_padded_integral_img(data.width, data.height);
    }

#ifdef WARM_UP
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            (*function)(gray_image, iimage);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);
#endif

    std::vector<double> cycleslist;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            (*function)(gray_image, iimage);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cycleslist.push_back(cycles);
    }
    total_cycles /= REP;

    free(iimage->padded_data);
    free(iimage);

    cycles = total_cycles;  // cyclesList.front();
    double flops_per_cycle = round((100.0 * data.num_flops) / cycles) / 100.0;
    std::sort(cycleslist.begin(), cycleslist.end());
    data.avg_cycles = (uint64_t)cycles;
    data.min_cycles = (uint64_t)cycleslist.front();
    data.max_cycles = (uint64_t)cycleslist.back();
    data.flops_per_cycle = flops_per_cycle;
}

void bench_compute_response_layer(const std::vector<void (*)(struct fasthessian *)> &functions,
                                  struct integral_image *iimage, std::vector<struct benchmark_data> &data) {
    // TODO: (carla) do we need a new fast hessian for every function?
    struct fasthessian *fh = create_fast_hessian(iimage);
    // Create octaves with response layers
    create_response_map(fh);

    assert(functions.size() == data.size());

    for (int j = 0; j < functions.size(); ++j) {
        perf_compute_response_all_layers(functions[j], fh, data[j]);
    }
    for (int i = 0; i < NUM_LAYERS; ++i) {
        free(fh->response_map[i]->response);
        free(fh->response_map[i]->laplacian);
        free(fh->response_map[i]);
    }
    free(fh);
}

/* Todo: (valentin) can these be removed?
// Times the function compute_response_layer from fasthessian.
// Stores the average, minimum and maximum number of cycles and the flops per cycle in benchmark_data.
// The number of flops for compute_response_layer must be set in benchmark_data.
void perf_compute_response_layer(void (*function)(struct response_layer *, struct integral_image *),
                                 struct response_layer *layer, struct integral_image *iimage,
                                 struct benchmark_data &data) {
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    uint64_t start, end;

#ifdef WARM_UP
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            (*function)(layer, iimage);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);
#endif

    std::vector<double> cycleslist;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            (*function)(layer, iimage);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cycleslist.push_back(cycles);
    }
    total_cycles /= REP;

    cycles = total_cycles;  // cyclesList.front();
    double flops_per_cycle = round((100.0 * data.num_flops) / cycles) / 100.0;
    std::sort(cycleslist.begin(), cycleslist.end());
    data.avg_cycles = (uint64_t)cycles;
    data.min_cycles = (uint64_t)cycleslist.front();
    data.max_cycles = (uint64_t)cycleslist.back();
    data.flops_per_cycle = flops_per_cycle;
}
*/

// Times the function compute_response_layer from fasthessian.
// Stores the average, minimum and maximum number of cycles and the flops per cycle in benchmark_data.
// The number of flops for compute_response_layer must be set in benchmark_data.
void perf_compute_response_all_layers(void (*function)(struct fasthessian *), struct fasthessian *fh,
                                      struct benchmark_data &data) {
    double cycles = 0.;
    long num_runs = 5;
    double multiplier = 1;
    uint64_t start, end;

#ifdef WARM_UP
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            // TODO: (carla) should we do a warm up on more than one layer? feels like an overkill.
            (*function)(fh);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);
#endif

    std::vector<double> cycleslist;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            (*function)(fh);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cycleslist.push_back(cycles);
    }
    total_cycles /= REP;
    // total_cycles /= MAX(1, fh->total_layers);
    cycles = total_cycles;  // cyclesList.front();
    double flops_per_cycle = round((100.0 * data.num_flops) / cycles) / 100.0;
    std::sort(cycleslist.begin(), cycleslist.end());
    data.avg_cycles = (uint64_t)cycles;
    data.min_cycles = (uint64_t)cycleslist.front();  // / MAX(1, fh->total_layers);
    data.max_cycles = (uint64_t)cycleslist.back();   // / MAX(1, fh->total_layers);
    // TODO: (carla) do we have to divide by total_lyers here? the total_cycles have already been divided;
    data.flops_per_cycle = flops_per_cycle;
}

void bench_get_interest_points(
    const std::vector<void (*)(struct fasthessian *, std::vector<struct interest_point> *)> &functions,
    struct fasthessian *fh, std::vector<struct benchmark_data> &data) {
    assert(functions.size() == data.size());

    for (int j = 0; j < functions.size(); ++j) {

        perf_get_interest_points(functions[j], fh, data[j]);

    }

}

// Times the function get_interest_points from fasthessian.
// Stores the average, minimum and maximum number of cycles and the flops per cycle in benchmark_data.
// The number of flops for get_interest_points must be set in benchmark_data.
// The number of interest points that will be found, must be set in benchmark_data.
void perf_get_interest_points(void (*function)(struct fasthessian *, std::vector<struct interest_point> *),
                              struct fasthessian *fh, struct benchmark_data &data) {
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    uint64_t start, end;

#ifdef WARM_UP
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        std::vector<struct interest_point> dummy_interest_points;
        dummy_interest_points.reserve(num_runs * data.num_interest_points);

        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            (*function)(fh, &dummy_interest_points);
        }
        end = stop_tsc(start);
        dummy_interest_points.clear();

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);
#endif

    std::vector<double> cycleslist;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        std::vector<struct interest_point> dummy_interest_points;
        dummy_interest_points.reserve(num_runs * data.num_interest_points);

        start = start_tsc();

        for (size_t i = 0; i < num_runs; ++i) {
            (*function)(fh, &dummy_interest_points);
        }
        end = stop_tsc(start);

        dummy_interest_points.clear();

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cycleslist.push_back(cycles);
    }
    total_cycles /= REP;

    cycles = total_cycles;  // cyclesList.front();
    double flops_per_cycle = round((100.0 * data.num_flops) / cycles) / 100.0;
    std::sort(cycleslist.begin(), cycleslist.end());
    data.avg_cycles = (uint64_t)cycles;
    data.min_cycles = (uint64_t)cycleslist.front();
    data.max_cycles = (uint64_t)cycleslist.back();
    data.flops_per_cycle = flops_per_cycle;
}

// Calls the timing function for interpolate_step. The number of times the timing should be conducted can be
// specified in this function. 
// The number of flops for interpolate_step must be set in benchmark_data.
void bench_interpolate_step(const std::vector<void (*)(int, int, struct response_layer *, struct response_layer *, struct response_layer *, float[3])> &functions, 
                            struct fasthessian *fh, 
                            std::vector<struct benchmark_data> &data) {
    
    assert(functions.size() == data.size());

    // Iterating through all functions that should be benchmarked
    for (int j = 0; j < functions.size(); ++j) {

        int counter = 1;
        // Specifies how many times the benchmarikng for interpolate_step is called.
        // The average of the outcome will be taken
        int limit = 4;

        assert(fh != NULL);

        // filter index map
        const int filter_map[NUM_OCTAVES][NUM_LAYERS] = {
            {0, 1, 2, 3}, 
            {1, 3, 4, 5}, 
            {3, 5, 6, 7},
            //{5, 7, 8, 9},
            //{7, 9, 10, 11}
        };

        // getting response layers
        struct response_layer *bottom;
        struct response_layer *middle;
        struct response_layer *top;

        // iterating through all octaves and each layer of octave in window of three (top, middle, bottom)
        for (int o = 0; o < fh->octaves; ++o) {
            // TODO: (Sebastian) allow for fh->layers != 4 as well (note that fh->layers>=3 has to hold)
            for (int i = 0; i <= 1; ++i) {
                // assigning respective bottom, middle and top response layer
                bottom = fh->response_map[filter_map[o][i]];
                middle = fh->response_map[filter_map[o][i + 1]];
                top = fh->response_map[filter_map[o][i + 2]];

                // iterating over middle response layer at density of the most sparse layer (always top),
                // to find maxima accreoss scale and space
                for (int r = 0; r < top->height; ++r) {
                    for (int c = 0; c < top->width; ++c) {
                        // checking if current pixel position is local maximum in 3x3x3 maximum and above threshold
                        if (is_extremum(r, c, top, middle, bottom, fh->thresh)) {
                            perf_interpolate_step(functions[j], r, c, top, middle, bottom, data[j]);
                            counter++;
                            if (counter == limit) {
                                c = top->width;
                                r = top->height;
                                i = 2;
                                o = fh->octaves;
                            }
                        }
                    }
                }
            }
        }

        data[j].avg_cycles /= counter;
        data[j].max_cycles /= counter;
        data[j].min_cycles /= counter;
        data[j].flops_per_cycle /= counter;

    }

}

// Times the function interpolate_step from fasthessian. 
// Stores the average, minimum and maximum number of cycles and the flops per cycle in benchmark_data. 
// The number of flops for interpolate_step must be set in benchmark_data.
void perf_interpolate_step(void (*function)(int, int, struct response_layer *, struct response_layer *, struct response_layer *, float[3]),
                           int row, int col, struct response_layer *top, struct response_layer *middle,
                           struct response_layer *bottom, struct benchmark_data &data) {
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    uint64_t start, end;

#ifdef WARM_UP
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            float offsets[3];
            (*function)(row, col, top, middle, bottom, offsets);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);
#endif

    std::vector<double> cycleslist;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            float offsets[3];
            (*function)(row, col, top, middle, bottom, offsets);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cycleslist.push_back(cycles);
    }
    total_cycles /= REP;

    cycles = total_cycles;  // cyclesList.front();
    double flops_per_cycle = round((100.0 * data.num_flops) / cycles) / 100.0;
    std::sort(cycleslist.begin(), cycleslist.end());
    data.avg_cycles += (uint64_t)cycles;
    data.min_cycles += (uint64_t)cycleslist.front();
    data.max_cycles += (uint64_t)cycleslist.back();
    data.flops_per_cycle += flops_per_cycle;
}

// Calls the timing function for multiple versions of get_descriptor.
// The size of vector functions has to match the size of vector data.
// The number of times the timing should be conducted can be specified in this function.
// The number of flops for the version of get_descriptor must be set in the respective benchmark_data.
void bench_get_descriptor(
    const std::vector<void (*)(struct integral_image *, struct interest_point *, float *)> &functions,
    struct integral_image *iimage, std::vector<struct interest_point> *interest_points, float *GW,
    std::vector<struct benchmark_data> &data) {
    assert(functions.size() == data.size());

    // Iterating through all functions that should be benchmarked
    for (int j = 0; j < functions.size(); ++j) {
        // Specifies how many times the timing for get_descriptor will be called. The average will be taken.
        int counter = MIN(5, data[j].num_interest_points);
        for (int i = 0; i < counter; ++i) {
            struct interest_point ipoint = interest_points->at(i);
            perf_get_descriptor(functions[j], iimage, &ipoint, GW, data[j]);
        }

        // Take the average of the runs.
        if (counter != 0) {
            data[j].avg_cycles /= counter;
            data[j].max_cycles /= counter;
            data[j].min_cycles /= counter;
            data[j].flops_per_cycle /= counter;
        }
    }
}

// Times the function get_descriptor from descriptor.
// Stores the average, minimum and maximum number of cycles and the flops per cycle in benchmark_data.
// The number of flops for get_descriptor must be set in the benchmark_data.
void perf_get_descriptor(void (*function)(struct integral_image *, struct interest_point *, float *),
                         struct integral_image *iimage, struct interest_point *ipoint, float *GW,
                         struct benchmark_data &data) {
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    uint64_t start, end;

#ifdef WARM_UP
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            (*function)(iimage, ipoint, GW);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);
#endif

    std::vector<double> cycleslist;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            (*function)(iimage, ipoint, GW);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cycleslist.push_back(cycles);
    }
    total_cycles /= REP;

    cycles = total_cycles;  // cyclesList.front();
    double flops_per_cycle = round((100.0 * data.num_flops) / cycles) / 100.0;
    std::sort(cycleslist.begin(), cycleslist.end());
    data.avg_cycles += (uint64_t)cycles;
    data.min_cycles += (uint64_t)cycleslist.front();
    data.max_cycles += (uint64_t)cycleslist.back();
    data.flops_per_cycle += flops_per_cycle;
}

// Calls timing function for multiple versions of get_msurf_descriptor.
// The number of times the timing should be conducted can be specified in this function.
// The number of flops for the version of get_msurf_descriptor must be set in the respective benchmark_data.
void bench_get_msurf_descriptor(
    const std::vector<void (*)(struct integral_image *, struct interest_point *)> &functions,
    struct integral_image *iimage, std::vector<struct interest_point> *interest_points,
    std::vector<struct benchmark_data> &data) {
    assert(functions.size() == data.size());

    // Iterating through all functions that should be benchmarked
    for (int j = 0; j < functions.size(); ++j) {
        // Specifies how many times the timing for get_msurf_descriptor will be called. The average will be taken.
        int counter = MIN(5, data[j].num_interest_points);
        for (int i = 0; i < counter; ++i) {
            struct interest_point ipoint = interest_points->at(i);
            perf_get_msurf_descriptor(functions[j], iimage, &ipoint, data[j]);
        }

        // Take the average of the runs.
        if (counter != 0) {
            data[j].avg_cycles /= counter;
            data[j].max_cycles /= counter;
            data[j].min_cycles /= counter;
            data[j].flops_per_cycle /= counter;
        }
    }
}

// Times the function get_msurf_descriptor from descriptor.
// Stores the average, minimum and maximum number of cycles and the flops per cycle in benchmark_data.
// The number of flops for get_msurf_descriptor must be set in benchmark_data.
void perf_get_msurf_descriptor(void (*function)(struct integral_image *, struct interest_point *),
                               struct integral_image *iimage, struct interest_point *ipoint,
                               struct benchmark_data &data) {
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    uint64_t start, end;

#ifdef WARM_UP
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            (*function)(iimage, ipoint);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);
#endif

    std::vector<double> cycleslist;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            (*function)(iimage, ipoint);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cycleslist.push_back(cycles);
    }
    total_cycles /= REP;

    cycles = total_cycles;  // cyclesList.front();
    double flops_per_cycle = round((100.0 * data.num_flops) / cycles) / 100.0;
    std::sort(cycleslist.begin(), cycleslist.end());
    data.avg_cycles += (uint64_t)cycles;
    data.min_cycles += (uint64_t)cycleslist.front();
    data.max_cycles += (uint64_t)cycleslist.back();
    data.flops_per_cycle += flops_per_cycle;
}

// Calls timing function for multiple versions of get_msurf_descriptors. 
// The number of times the timing should be conducted can be specified in this function. 
// The number of flops for the version of get_msurf_descriptors must be set in the respective benchmark_data.
void bench_get_msurf_descriptors(const std::vector<void (*)(struct integral_image *, std::vector<struct interest_point> *)> &functions,
                                struct integral_image *iimage, std::vector<struct interest_point> *interest_points,
                                std::vector<struct benchmark_data> &data) {
    
    assert(functions.size() == data.size());

    // Iterating through all functions that should be benchmarked
    for (int j = 0; j < functions.size(); ++j) {
        perf_get_msurf_descriptors(functions[j], iimage, interest_points, data[j]);
    }
    
}

// Times the function get_msurf_descriptors from descriptor. 
// Stores the average, minimum and maximum number of cycles and the flops per cycle in benchmark_data. 
// The number of flops for get_msurf_descriptor must be set in benchmark_data.
void perf_get_msurf_descriptors(void (*function)(struct integral_image *, std::vector<struct interest_point> *),
                               struct integral_image *iimage, std::vector<struct interest_point> *ipoints,
                               struct benchmark_data &data) {
    double cycles = 0.;
    long num_runs = 1;
    double multiplier = 1;
    uint64_t start, end;

#ifdef WARM_UP
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            (*function)(iimage, ipoints);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);
#endif

    std::vector<double> cycleslist;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            (*function)(iimage, ipoints);
        }
        end = stop_tsc(start);

        // get avg cycles per keypoint
        cycles = ((double)end) / (num_runs*ipoints->size());
        total_cycles += cycles;

        cycleslist.push_back(cycles);
    }
    total_cycles /= REP;

    cycles = total_cycles;  // cyclesList.front();
    double flops_per_cycle = round((100.0 * data.num_flops) / cycles) / 100.0;
    std::sort(cycleslist.begin(), cycleslist.end());
    data.avg_cycles += (uint64_t)cycles;
    data.min_cycles += (uint64_t)cycleslist.front();
    data.max_cycles += (uint64_t)cycleslist.back();
    data.flops_per_cycle += flops_per_cycle;
}
