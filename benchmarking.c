#include <stdio.h>

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
    return  round((100.0 * flops) / cycles) / 100.0;
}

void bench_integral_img(float* image, int width, int height){
    FILE *fptr_int_img;
	// Creates a file "fptr_int_img" 
    // with file acccess as write mode 
	fptr_int_img = fopen("benchmarking_files/integral_img.txt","w");
	int flops = width + 2*(height-1)*width;
	double x  = perf_test_integral_img(create_integral_img, image, width, height, flops);
    fprintf(fptr_int_img,"%lf \n",x);

	// closes the file pointed by fptr_int_img 
    fclose(fptr_int_img); 

}
