// has to be defined before stb includes
#define STB_IMAGE_IMPLEMENTATION

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
#include "descriptor_opt_gen.h"
#include "fasthessian_opt.h"

const char *images[] = {
    // "../images/sunflower/sunflower_32.jpg",  
    // "../images/sunflower/sunflower_64.jpg",
    // "../images/sunflower/sunflower_128.jpg", 
    // "../images/sunflower/sunflower_256.jpg",
    "../images/sunflower/sunflower_512.jpg",
    // "../images/sunflower/sunflower_1024.jpg",
    // "../images/sunflower/sunflower_2048.jpg",
    // "../images/sunflower/sunflower_4096.jpg",
};

#define n_images (sizeof(images) / sizeof(const char *))

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


        // Fast-Hessian
        struct fasthessian *fh = create_fast_hessian(iimage);

        // Create octaves with response layers
        create_response_map(fh);
        
        // Compute responses for every layer
        compute_response_layers(fh);

        // Getting interest points with non-maximum supression
        std::vector<struct interest_point> interest_points;
        get_interest_points(fh, &interest_points);

        {
            printf("get_msurf_descriptor start\n");
            
            // Insert all interpolate_step functions for benchmarking here
            std::vector<void (*)(struct integral_image *, std::vector<struct interest_point> *)> functions;
            // Insert all respective benchmarking info for functions here
            std::vector<struct benchmark_data> data;
            functions.push_back(get_msurf_descriptors_haar_unroll_1_1_False);
            struct benchmark_data data_1_1_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_1_False", interest_points.size(), -1);
            data.push_back(data_1_1_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_1_1_True);
            struct benchmark_data data_1_1_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_1_True", interest_points.size(), -1);
            data.push_back(data_1_1_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_1_2_False);
            struct benchmark_data data_1_2_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_2_False", interest_points.size(), -1);
            data.push_back(data_1_2_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_1_2_True);
            struct benchmark_data data_1_2_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_2_True", interest_points.size(), -1);
            data.push_back(data_1_2_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_1_3_False);
            struct benchmark_data data_1_3_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_3_False", interest_points.size(), -1);
            data.push_back(data_1_3_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_1_3_True);
            struct benchmark_data data_1_3_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_3_True", interest_points.size(), -1);
            data.push_back(data_1_3_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_1_4_False);
            struct benchmark_data data_1_4_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_4_False", interest_points.size(), -1);
            data.push_back(data_1_4_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_1_4_True);
            struct benchmark_data data_1_4_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_4_True", interest_points.size(), -1);
            data.push_back(data_1_4_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_2_1_False);
            struct benchmark_data data_2_1_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_2_1_False", interest_points.size(), -1);
            data.push_back(data_2_1_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_2_1_True);
            struct benchmark_data data_2_1_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_2_1_True", interest_points.size(), -1);
            data.push_back(data_2_1_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_2_2_False);
            struct benchmark_data data_2_2_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_2_2_False", interest_points.size(), -1);
            data.push_back(data_2_2_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_2_2_True);
            struct benchmark_data data_2_2_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_2_2_True", interest_points.size(), -1);
            data.push_back(data_2_2_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_2_3_False);
            struct benchmark_data data_2_3_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_2_3_False", interest_points.size(), -1);
            data.push_back(data_2_3_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_2_3_True);
            struct benchmark_data data_2_3_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_2_3_True", interest_points.size(), -1);
            data.push_back(data_2_3_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_2_4_False);
            struct benchmark_data data_2_4_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_2_4_False", interest_points.size(), -1);
            data.push_back(data_2_4_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_2_4_True);
            struct benchmark_data data_2_4_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_2_4_True", interest_points.size(), -1);
            data.push_back(data_2_4_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_3_1_False);
            struct benchmark_data data_3_1_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_3_1_False", interest_points.size(), -1);
            data.push_back(data_3_1_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_3_1_True);
            struct benchmark_data data_3_1_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_3_1_True", interest_points.size(), -1);
            data.push_back(data_3_1_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_3_2_False);
            struct benchmark_data data_3_2_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_3_2_False", interest_points.size(), -1);
            data.push_back(data_3_2_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_3_2_True);
            struct benchmark_data data_3_2_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_3_2_True", interest_points.size(), -1);
            data.push_back(data_3_2_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_3_3_False);
            struct benchmark_data data_3_3_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_3_3_False", interest_points.size(), -1);
            data.push_back(data_3_3_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_3_3_True);
            struct benchmark_data data_3_3_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_3_3_True", interest_points.size(), -1);
            data.push_back(data_3_3_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_3_4_False);
            struct benchmark_data data_3_4_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_3_4_False", interest_points.size(), -1);
            data.push_back(data_3_4_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_3_4_True);
            struct benchmark_data data_3_4_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_3_4_True", interest_points.size(), -1);
            data.push_back(data_3_4_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_4_1_False);
            struct benchmark_data data_4_1_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_1_False", interest_points.size(), -1);
            data.push_back(data_4_1_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_4_1_True);
            struct benchmark_data data_4_1_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_1_True", interest_points.size(), -1);
            data.push_back(data_4_1_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_4_2_False);
            struct benchmark_data data_4_2_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_2_False", interest_points.size(), -1);
            data.push_back(data_4_2_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_4_2_True);
            struct benchmark_data data_4_2_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_2_True", interest_points.size(), -1);
            data.push_back(data_4_2_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_4_3_False);
            struct benchmark_data data_4_3_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_3_False", interest_points.size(), -1);
            data.push_back(data_4_3_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_4_3_True);
            struct benchmark_data data_4_3_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_3_True", interest_points.size(), -1);
            data.push_back(data_4_3_True);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_4_4_False);
            struct benchmark_data data_4_4_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_4_False", interest_points.size(), -1);
            data.push_back(data_4_4_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unroll_4_4_True);
            struct benchmark_data data_4_4_True(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_4_True", interest_points.size(), -1);
            data.push_back(data_4_4_True);;
            
            // Benchmarking all get_msurf_descriptor functions and storing timing results in respective entries in data
            bench_get_msurf_descriptors(functions, iimage, &interest_points, data);

            // Appending this data to all benchmarking data
            all_benchmark_data.insert(all_benchmark_data.end(), data.begin(), data.end());

            printf("get_msurf_descriptor end\n");
        }

        // Getting M-SURF descriptors for each interest point
        get_msurf_descriptors(iimage, &interest_points);

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
    
    extern float* haarResponseX;
    extern float* haarResponseY;

    aligned_free(haarResponseX);
    aligned_free(haarResponseY);
    
    save_benchmark_data(all_benchmark_data);
    // free memory benchmarkdata
    // https://stackoverflow.com/questions/10464992/c-delete-vector-objects-free-memory
    //std::vector<struct benchmark_data *>().swap(all_benchmark_data);
    printf("Benchmarking done!\n");

    return 0;
}
