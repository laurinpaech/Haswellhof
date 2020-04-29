#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm> 


#include "tsc_x86.h"
#include "benchmarking.h"

#define CYCLES_REQUIRED 1e8
#define REP 50
#define EPS (1e-3)


/*
* 
* reports and returns the number of cycles required per iteration
*/
void perf_test_integral_img(struct integral_image* (*function)(float*, int, int), float* gray_image, struct benchmark_data* data)
{
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    int start, end;
    struct integral_image* iimage;

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            iimage = (*function)(gray_image, data->width, data->height);           
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);
        
    } while (multiplier > 2);

    std::vector<uint64_t>cycleslist;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            iimage = (*function)(gray_image, data->width, data->height);                   
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cycleslist.push_back(cycles);
    }
    total_cycles /= REP;

    free(iimage->data);
	free(iimage);

    cycles = total_cycles;//cyclesList.front();
    double flops_per_cycle = round((100.0 * data->num_flops) / cycles) / 100.0;
    std::sort(cycleslist.begin(), cycleslist.end());   
    data->avg_cycles = cycles;
    data->min_cycles = cycleslist.front();
    data->max_cycles = cycleslist.back();
    data->flops_per_cycle = flops_per_cycle;
}




void perf_test_compute_response_layer(void (*function)(struct response_layer*, struct integral_image*), struct response_layer* layer,struct integral_image* iimage, 
struct benchmark_data* data){

    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    int start, end;

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

    std::vector<uint64_t>cycleslist;

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

    cycles = total_cycles;//cyclesList.front();
    double flops_per_cycle = round((100.0 * data->num_flops) / cycles) / 100.0;
    std::sort(cycleslist.begin(), cycleslist.end());  
    data->avg_cycles = cycles;
    data->min_cycles = cycleslist.front();
    data->max_cycles = cycleslist.back();
    data->flops_per_cycle = flops_per_cycle;
}


//times the function get_interest_points from fasthessian and returns the flops per cycle
void perf_test_get_interest_points(void (*function)(struct fasthessian*, std::vector<struct interest_point>*), struct fasthessian *fh, struct benchmark_data* data){
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    int start, end;

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            std::vector<struct interest_point> interest_points;
            (*function)(fh, &interest_points);           
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);
        
    } while (multiplier > 2);

    std::vector<uint64_t>cycleslist;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            std::vector<struct interest_point> interest_points;
            (*function)(fh, &interest_points);           
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cycleslist.push_back(cycles);
    }
    total_cycles /= REP;

    cycles = total_cycles;//cyclesList.front();
    double flops_per_cycle = round((100.0 * data->num_flops) / cycles) / 100.0;
    std::sort(cycleslist.begin(), cycleslist.end());  
    data->avg_cycles = cycles;
    data->min_cycles = cycleslist.front();
    data->max_cycles = cycleslist.back();
    data->flops_per_cycle = flops_per_cycle;
}


void perf_interpolate_step(void (*function)(int, int, struct response_layer*, struct response_layer *, struct response_layer*, float[3]),int row, int col, struct response_layer *top, 
struct response_layer *middle, struct response_layer *bottom, struct benchmark_data* data){

    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    int start, end;

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

    std::vector<uint64_t>cycleslist;

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

    cycles = total_cycles;//cyclesList.front();
    double flops_per_cycle = round((100.0 * data->num_flops) / cycles) / 100.0;
    std::sort(cycleslist.begin(), cycleslist.end());  
    data->avg_cycles = cycles;
    data->min_cycles = cycleslist.front();
    data->max_cycles = cycleslist.back();
    data->flops_per_cycle = flops_per_cycle;
}



void bench_interpolate_step(struct fasthessian *fh, std::vector<struct interest_point> *interest_points, struct benchmark_data* data){


    int counter = 0;
    int limit = 3;

    assert(fh != NULL);
    assert(interest_points != NULL);

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
            middle = fh->response_map[filter_map[o][i+1]];
            top = fh->response_map[filter_map[o][i+2]];

            // iterating over middle response layer at density of the most sparse layer (always top),
            // to find maxima accreoss scale and space
            for (int r = 0; r < top->height; ++r) {
                for (int c = 0; c < top->width; ++c) {

                    // checking if current pixel position is local maximum in 3x3x3 maximum and above threshold
                    if (is_extremum(r, c, top, middle, bottom, fh->thresh)) {
                        
                        perf_interpolate_step(interpolate_step, r, c, top, middle, bottom, data);
                        counter++;
                        if(counter == limit){
                            c = top->width;
                            r = top->height;
                            i = 2;
                            o = fh->octaves;
                        }
                        printf("%i\n", counter);
                    }

                }
            }

        }
    }

    data->avg_cycles/=counter;
    data->max_cycles/=counter;
    data->min_cycles/=counter;
    data->flops_per_cycle/=counter;
}