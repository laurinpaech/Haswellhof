#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
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


//times the function interpolate_step from fasthessian and returns the flops per cycle
double perf_test_interpolate_step(void (*function)(int, int, struct response_layer*, struct response_layer*, struct response_layer*, float), int row, int col, 
    struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets, int flops){
        
    }



//times the function interpolate_step from fasthessian and saves the output to a file
void bench_interpolate_step(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets){
    int flops = 109;
    //printf("%li", flops);
    printf("bench_get_interest_points %li\n", flops);
	//double flops_per_cycle  = perf_test_interpolate_step(interpolate_step, row, col, top, middle, bottom, offsets, flops);
    //save_performance_file(flops_per_cycle, "/interpolate_step");
}

/**
void save_performance_file(double flops_per_cycle, const char *file_name){   
    // create the folder if it doesn't exist   
    char folder_name[] = "benchmarking_files";
    struct stat st = {0};
    if (stat(folder_name, &st) == -1) {
        mkdir(folder_name, 0700);
    }


	// Creates a file "fptr_int_img" 
    // with file acccess as write mode 
    FILE *fptr_int_img;

    // get the current time
    struct tm *timenow;
    time_t now = time(NULL);
    timenow = gmtime(&now);
    char current_time[30];
    strftime(current_time, sizeof(current_time), "_%Y-%m-%d_%H_%M.txt", timenow);

    // concatenate the path in the form  "benchmarking_files/integral_img_CURRENT_DATE"
    char* path_name = concat(folder_name, file_name);
    path_name = concat(path_name, current_time);

    // write to file
    fptr_int_img =fopen(path_name,"w");
    fprintf(fptr_int_img,"%lf \n",flops_per_cycle);

	// closes the file pointed by fptr_int_img 
    fclose(fptr_int_img); 
    free(path_name);
}

char* concat(const char *s1, const char *s2)
{
    const size_t len1 = strlen(s1);
    const size_t len2 = strlen(s2);
    char *result = (char *)malloc(len1 + len2 + 1); // +1 for the null-terminator
    // in real code you would check for errors in malloc here
    memcpy(result, s1, len1);
    memcpy(result + len1, s2, len2 + 1); // +1 to copy the null-terminator
    return result;
}
*/