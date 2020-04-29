// has to be defined before stb includes
#define STB_IMAGE_IMPLEMENTATION

#define USE_MSURF 0

#include "stb_image.h"
#include "integral_image.h"
#include "fasthessian.h"
#include "interest_point.h"
#include "descriptor.h"
#include "benchmarking.h"
#include "benchmark_data_to_file.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>


const char* images[] = {
    "../images/sunflower/sunflower_32.jpg",
    "../images/sunflower/sunflower_64.jpg"
};
#define n_images (sizeof (images) / sizeof (const char *))
#define BENCHMARK_INTEGRAL_IMAGE
#define BENCHMARK_CREATE_RESPONSE_MAP
//#define BENCHMARK_INTEREST_POINTS
//#define BENCHMARK_INTERPOLATE_STEPS
#define BENCHMARK_GET_DESCRIPTORS

int main(int argc, char const *argv[])
{
    std::vector<struct benchmark_data*> all_benchmark_data;
    for (int i = 0; i < n_images; i++) {
        char* image_name = (char*) malloc(32*sizeof(char));
        strcpy(image_name, images[i]);

        int width, height, channels;

        // Load image
        stbi_ldr_to_hdr_gamma(1.0f);
        printf("%s\n", image_name);
        float* image = stbi_loadf(image_name, &width, &height, &channels, STBI_grey);
        if(!image) {
            printf("Could not open or find image\n");
            return -1;
        }

        // Calculate integral image
        struct integral_image* iimage = create_integral_img(image, width, height);

#ifdef BENCHMARK_INTEGRAL_IMAGE
        printf("create_integral_img start\n");
        struct benchmark_data* benchmark_integral_img=initialise_benchmark_data(image_name, width, height, "create_integral_img", -1, (width + 2*(height-1)*width));
        printf("benchmark data integral_img: %s\n",benchmark_integral_img->image_name);
        perf_test_integral_img(create_integral_img, image, benchmark_integral_img);
        all_benchmark_data.push_back(benchmark_integral_img);
        printf("create_integral_img done\n");
#endif

        // Fast-Hessian
        struct fasthessian* fh = create_fast_hessian(iimage);

        // Create octaves with response layers
        create_response_map(fh);

#ifdef BENCHMARK_CREATE_RESPONSE_MAP
        struct benchmark_data* benchmark_compute_response_layer=initialise_benchmark_data(image_name, width, height, "compute_response_layer", -1, (1 + height*width*13));
        perf_test_compute_response_layer(compute_response_layer, fh->response_map[0], iimage, benchmark_compute_response_layer);
        all_benchmark_data.push_back(benchmark_compute_response_layer);
#endif
        // Compute responses for every layer
        for (int i = 0; i < fh->total_layers; i++) {
            compute_response_layer(fh->response_map[i], iimage);
        }

        // Getting interest points with non-maximum supression
        std::vector<struct interest_point> interest_points;
        get_interest_points(fh, &interest_points);

#ifdef BENCHMARK_INTEREST_POINTS
        printf("get_interest_points start\n");       
        long flops = 109*interest_points.size();
        struct benchmark_data* benchmark_get_interest_points=initialise_benchmark_data(image_name, width, height, "get_interest_points", interest_points.size(), flops);
        perf_test_get_interest_points(get_interest_points, fh, benchmark_get_interest_points);
        all_benchmark_data.push_back(benchmark_get_interest_points);
        printf("get_interest_points end\n");
#endif 

#ifdef BENCHMARK_INTERPOLATE_STEPS
        printf("interpolate_step start\n");
        struct benchmark_data* benchmark_interpolate_step=initialise_benchmark_data(image_name, width, height, "interpolate_step", interest_points.size(), 109);
        bench_interpolate_step(fh, benchmark_interpolate_step);
        all_benchmark_data.push_back(benchmark_interpolate_step);
        printf("interpolate_step end\n");
#endif

    #if !USE_MSURF
        // Descriptor stuff
        float* GW = get_gaussian(3.3);

#ifdef BENCHMARK_GET_DESCRIPTORS
        printf("get_descriptor start");
        struct benchmark_data* benchmark_get_descriptor=initialise_benchmark_data(image_name, width, height, "get_descriptor", interest_points.size(), 5734);
        bench_get_descriptor(iimage, &interest_points, GW,benchmark_get_descriptor);
        all_benchmark_data.push_back(benchmark_get_descriptor);
        printf("get_descriptor end");
#endif
        for (size_t i=0; i<interest_points.size(); ++i)
            get_descriptor(iimage, &interest_points[i], GW);
        free(GW);
    #else
        // Alternative M-SURF descriptors as in OpenSURF
        for (size_t i=0; i<interest_points.size(); ++i)
            get_msurf_descriptor(iimage, &interest_points[i]);
    #endif


        // Free memory
        stbi_image_free(image); // possibly move this to create_integral_img
        free(iimage->data);

        free(iimage);
        for (size_t i = 0; i < NUM_LAYERS; i++) {
            free(fh->response_map[i]);
        }
        free(fh);
        free(image_name);
    }
    save_benchmark_data(all_benchmark_data);
    // free memory benchmarkdata
    // https://stackoverflow.com/questions/10464992/c-delete-vector-objects-free-memory
    std::vector<struct benchmark_data*>().swap(all_benchmark_data);
    printf("done\n");

	return 0;
}
