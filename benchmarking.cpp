#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

#include "tsc_x86.h"
#include "benchmarking.h"

#define CYCLES_REQUIRED 1e8
#define REP 50
#define EPS (1e-3)

/*
* Checks the given function for validity. If valid, then computes and
* reports and returns the number of cycles required per iteration
*/
double perf_test_integral_img(struct integral_image* (*function)(float*, int, int), float* gray_image, int width, int height,int flops)
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
            iimage = (*function)(gray_image, width, height);           
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);
        
    } while (multiplier > 2);



    float *cyclesList = (float*)malloc(REP * sizeof(float));


    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            iimage = (*function)(gray_image, width, height);                   
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cyclesList[REP-j-1]=(cycles);
    }
    total_cycles /= REP;

    free(iimage->data);
	free(iimage);


    cycles = total_cycles;//cyclesList.front();
    free(cyclesList);
    return  round((100.0 * flops) / cycles) / 100.0;
}

// time the function create_integral_img and write it to a file
void bench_integral_img(float* image, int width, int height){
    int flops = width + 2*(height-1)*width;
	double flops_per_cycle  = perf_test_integral_img(create_integral_img, image, width, height, flops);
    save_performance_file(flops_per_cycle, "/integral_img");
}

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