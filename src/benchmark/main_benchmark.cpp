// has to be defined before stb includes
#define STB_IMAGE_IMPLEMENTATION

#define USE_MSURF 1

#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "benchmark_data_to_file.h"
#include "benchmarking.h"
#include "descriptor.h"
#include "fasthessian.h"
#include "integral_image.h"
#include "interest_point.h"
#include "stb_image.h"
#include "descriptor_opt.h"

const char *images[] = {
    // "../images/sunflower/sunflower_32.jpg",  
    // "../images/sunflower/sunflower_64.jpg",
    // "../images/sunflower/sunflower_128.jpg", 
    // "../images/sunflower/sunflower_256.jpg"
    "../images/sunflower/sunflower_512.jpg",
    // "../images/sunflower/sunflower_1024.jpg",
    // "../images/sunflower/sunflower_2048.jpg"
    // "../images/sunflower/sunflower_4096.jpg"
};
#define n_images (sizeof(images) / sizeof(const char *))
#define BENCHMARK_INTEGRAL_IMAGE
#define BENCHMARK_CREATE_RESPONSE_MAP
#define BENCHMARK_INTEREST_POINTS
#define BENCHMARK_INTERPOLATE_STEPS
#define BENCHMARK_GET_DESCRIPTORS
#define BENCHMARK_GET_MSURF_DESCRIPTORS

int main(int argc, char const *argv[]) {
    std::vector<struct benchmark_data> all_benchmark_data;
    for (int i = 0; i < n_images; i++) {
        char *image_name = (char *)malloc(32 * sizeof(char));
        strcpy(image_name, images[i]);

        int width, height, channels;

        // Load image
        stbi_ldr_to_hdr_gamma(1.0f);
        printf("%s\n", image_name);
        float *image = stbi_loadf(image_name, &width, &height, &channels, STBI_grey);
        if (!image) {
            printf("Could not open or find image\n");
            return -1;
        }

        // Create integral image
        struct integral_image *iimage = create_integral_img(width, height);
        // Compute integral image
        compute_integral_img(image, iimage->width, iimage->height, iimage->data);



#ifdef BENCHMARK_INTEGRAL_IMAGE
        {
            printf("compute_integral_img start\n");

            // Insert all compute_integral_img functions for benchmarking here
            std::vector<void (*)(float *, int, int, float *)> functions;
            functions.push_back(compute_integral_img);
            //functions.push_back(compute_integral_img_faster_alg);

            struct benchmark_data default_data(image_name, width, height, "compute_integral_img", -1, (width + 2 * (height - 1) * width));

            // Insert all respective benchmarking info for compute_integral_img here
            std::vector<struct benchmark_data> data;
            data.push_back(default_data);
            //data.emplace_back(image_name, width, height, "compute_integral_img_faster_alg", -1, get_flops_compute_integral_img_faster_alg(width, height, 2));

            // Benchmarking all compute_integral_img functions and storing timing results in respective entries in data
            bench_compute_integral_img(functions, image, data);
            
            // Appending this data to all benchmarking data
            all_benchmark_data.insert(all_benchmark_data.end(), data.begin(), data.end());

            printf("compute_integral_img end\n");
        }
#endif

        // Fast-Hessian
        struct fasthessian *fh = create_fast_hessian(iimage);

        // Create octaves with response layers
        create_response_map(fh);
        
#ifdef BENCHMARK_CREATE_RESPONSE_MAP
        {
            printf("compute_response_layer start\n");
            struct benchmark_data data(image_name, width, height, "compute_response_layer", -1, (1 + height * width * 13));
            perf_compute_response_layer(compute_response_layer, fh->response_map[0], iimage, data);
            all_benchmark_data.push_back(data);
            printf("compute_response_layer end\n");
        }
#endif
        
        // Compute responses for every layer
	    compute_response_map(fh);

        // Getting interest points with non-maximum supression
        std::vector<struct interest_point> interest_points;
        get_interest_points(fh, &interest_points);

#ifdef BENCHMARK_INTEREST_POINTS
        {
            printf("get_interest_points start\n");
            
            // Insert all get_interest_points functions for benchmarking here
            std::vector<void (*)(struct fasthessian *, std::vector<struct interest_point> *)> functions;
            functions.push_back(get_interest_points);
            
            long flops = 109 * interest_points.size();
            struct benchmark_data default_data(image_name, width, height, "get_interest_points", interest_points.size(), flops);

            // Insert all respective benchmarking info for get_interest_points here
            std::vector<struct benchmark_data> data;
            data.push_back(default_data);

            // Benchmarking all get_interest_point functions and storing timing results in respective entries in data
            bench_get_interest_points(functions, fh, data);
            
            // Appending this data to all benchmarking data
            all_benchmark_data.insert(all_benchmark_data.end(), data.begin(), data.end());
            printf("get_interest_points end\n");
        }
#endif

#ifdef BENCHMARK_INTERPOLATE_STEPS
        {
            printf("interpolate_step start\n");
            
            // Insert all interpolate_step functions for benchmarking here
            std::vector<void (*)(int, int, struct response_layer *, struct response_layer *, struct response_layer *, float[3])> functions;
            functions.push_back(interpolate_step);

            struct benchmark_data default_data(image_name, width, height, "interpolate_step", interest_points.size(), 109);
            
            // Insert all respective benchmarking info for functions here
            std::vector<struct benchmark_data> data;
            data.push_back(default_data);

            // Benchmarking all interpolate_step functions and storing timing results in respective entries in data
            bench_interpolate_step(functions, fh, data);

            // Appending this data to all benchmarking data
            all_benchmark_data.insert(all_benchmark_data.end(), data.begin(), data.end());

            printf("interpolate_step end\n");
        }
#endif

#if !USE_MSURF
        // Descriptor stuff
        float *GW = get_gaussian(3.3);

#ifdef BENCHMARK_GET_DESCRIPTORS
        {
            printf("get_descriptor start\n");

            // Insert all interpolate_step functions for benchmarking here
            std::vector<void (*)(struct integral_image *, struct interest_point *, float *)> functions;
            functions.push_back(get_descriptor);

            struct benchmark_data default_data(image_name, width, height, "get_descriptor", interest_points.size(), 5734);

            // Insert all respective benchmarking info for functions here
            std::vector<struct benchmark_data> data;
            data.push_back(default_data);

            // Benchmarking all get_descriptor functions and storing timing results in respective entries in data
            bench_get_descriptor(functions, iimage, &interest_points, GW, data);

            // Appending this data to all benchmarking data
            all_benchmark_data.insert(all_benchmark_data.end(), data.begin(), data.end());

            printf("get_descriptor end\n");
        }
#endif
        for (size_t i = 0; i < interest_points.size(); ++i) {
            get_descriptor(iimage, &interest_points[i], GW);
        }
        free(GW);
#else

#ifdef BENCHMARK_GET_MSURF_DESCRIPTORS
        {
            printf("get_msurf_descriptor start\n");
            
            // Insert all interpolate_step functions for benchmarking here
            std::vector<void (*)(struct integral_image *, std::vector<struct interest_point> *)> functions;
            functions.push_back(get_msurf_descriptors);
            functions.push_back(get_msurf_descriptors_improved);
            functions.push_back(get_msurf_descriptors_improved_flip);
            functions.push_back(get_msurf_descriptors_improved_flip_flip);
            functions.push_back(get_msurf_descriptors_inlined);
            functions.push_back(get_msurf_descriptors_inlinedHaarWavelets);
            functions.push_back(get_msurf_descriptors_inlinedHaarWavelets_precheck_boundaries);
            
            // TODO: (Sebastian) find FLOPS count for get_msurf_descriptor
            struct benchmark_data default_data(image_name, width, height, "get_msurf_descriptors", interest_points.size(), -1);
            struct benchmark_data data1(image_name, width, height, "get_msurf_descriptors_improved", interest_points.size(), -1);   
            struct benchmark_data data2(image_name, width, height, "get_msurf_descriptors_improved_flip", interest_points.size(), -1);
            struct benchmark_data data22(image_name, width, height, "get_msurf_descriptors_improved_flip_flip", interest_points.size(), -1);   
            struct benchmark_data data3(image_name, width, height, "get_msurf_descriptors_inlined", interest_points.size(), -1);
            struct benchmark_data data4(image_name, width, height, "get_msurf_descriptors_inlinedHaarWavelets", interest_points.size(), -1);   
            struct benchmark_data data5(image_name, width, height, "get_msurf_descriptors_inlinedHaarWavelets_precheck_boundaries", interest_points.size(), -1);                     

            // Insert all respective benchmarking info for functions here
            std::vector<struct benchmark_data> data;
            data.push_back(default_data);
            data.push_back(data1);
            data.push_back(data2);
            data.push_back(data22);
            data.push_back(data3);
            data.push_back(data4);
            data.push_back(data5);

            // Benchmarking all get_msurf_descriptor functions and storing timing results in respective entries in data
            bench_get_msurf_descriptors(functions, iimage, &interest_points, data);

            // Appending this data to all benchmarking data
            all_benchmark_data.insert(all_benchmark_data.end(), data.begin(), data.end());

            printf("get_msurf_descriptor end\n");
        }
#endif

        // Alternative M-SURF descriptors as in OpenSURF
        for (size_t i = 0; i < interest_points.size(); ++i) {
            get_msurf_descriptor(iimage, &interest_points[i]);
        }
#endif

        // Free memory
        stbi_image_free(image);  // possibly move this to create_integral_img

        free(iimage->data);
        free(iimage);

        for (int i = 0; i < NUM_LAYERS; ++i) {
            free(fh->response_map[i]->response);
		    free(fh->response_map[i]->laplacian);
            free(fh->response_map[i]);
        }
        free(fh);
        
        free(image_name);
    }
    save_benchmark_data(all_benchmark_data);
    // free memory benchmarkdata
    // https://stackoverflow.com/questions/10464992/c-delete-vector-objects-free-memory
    //std::vector<struct benchmark_data *>().swap(all_benchmark_data);
    printf("Benchmarking done!\n");

    return 0;
}
