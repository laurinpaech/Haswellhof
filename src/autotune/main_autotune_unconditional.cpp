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
        compute_integral_img(image, iimage);


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
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_1_1_False);
            struct benchmark_data data_1_1_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_1_False", interest_points.size(), -1);
            data.push_back(data_1_1_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_1_4_False);
            struct benchmark_data data_1_4_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_4_False", interest_points.size(), -1);
            data.push_back(data_1_4_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_1_8_False);
            struct benchmark_data data_1_8_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_8_False", interest_points.size(), -1);
            data.push_back(data_1_8_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_1_12_False);
            struct benchmark_data data_1_12_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_12_False", interest_points.size(), -1);
            data.push_back(data_1_12_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_1_16_False);
            struct benchmark_data data_1_16_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_16_False", interest_points.size(), -1);
            data.push_back(data_1_16_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_1_20_False);
            struct benchmark_data data_1_20_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_20_False", interest_points.size(), -1);
            data.push_back(data_1_20_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_1_24_False);
            struct benchmark_data data_1_24_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_1_24_False", interest_points.size(), -1);
            data.push_back(data_1_24_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_4_1_False);
            struct benchmark_data data_4_1_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_1_False", interest_points.size(), -1);
            data.push_back(data_4_1_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_4_4_False);
            struct benchmark_data data_4_4_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_4_False", interest_points.size(), -1);
            data.push_back(data_4_4_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_4_8_False);
            struct benchmark_data data_4_8_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_8_False", interest_points.size(), -1);
            data.push_back(data_4_8_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_4_12_False);
            struct benchmark_data data_4_12_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_12_False", interest_points.size(), -1);
            data.push_back(data_4_12_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_4_16_False);
            struct benchmark_data data_4_16_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_16_False", interest_points.size(), -1);
            data.push_back(data_4_16_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_4_20_False);
            struct benchmark_data data_4_20_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_20_False", interest_points.size(), -1);
            data.push_back(data_4_20_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_4_24_False);
            struct benchmark_data data_4_24_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_4_24_False", interest_points.size(), -1);
            data.push_back(data_4_24_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_8_1_False);
            struct benchmark_data data_8_1_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_8_1_False", interest_points.size(), -1);
            data.push_back(data_8_1_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_8_4_False);
            struct benchmark_data data_8_4_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_8_4_False", interest_points.size(), -1);
            data.push_back(data_8_4_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_8_8_False);
            struct benchmark_data data_8_8_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_8_8_False", interest_points.size(), -1);
            data.push_back(data_8_8_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_8_12_False);
            struct benchmark_data data_8_12_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_8_12_False", interest_points.size(), -1);
            data.push_back(data_8_12_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_8_16_False);
            struct benchmark_data data_8_16_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_8_16_False", interest_points.size(), -1);
            data.push_back(data_8_16_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_8_20_False);
            struct benchmark_data data_8_20_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_8_20_False", interest_points.size(), -1);
            data.push_back(data_8_20_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_8_24_False);
            struct benchmark_data data_8_24_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_8_24_False", interest_points.size(), -1);
            data.push_back(data_8_24_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_12_1_False);
            struct benchmark_data data_12_1_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_12_1_False", interest_points.size(), -1);
            data.push_back(data_12_1_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_12_4_False);
            struct benchmark_data data_12_4_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_12_4_False", interest_points.size(), -1);
            data.push_back(data_12_4_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_12_8_False);
            struct benchmark_data data_12_8_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_12_8_False", interest_points.size(), -1);
            data.push_back(data_12_8_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_12_12_False);
            struct benchmark_data data_12_12_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_12_12_False", interest_points.size(), -1);
            data.push_back(data_12_12_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_12_16_False);
            struct benchmark_data data_12_16_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_12_16_False", interest_points.size(), -1);
            data.push_back(data_12_16_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_12_20_False);
            struct benchmark_data data_12_20_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_12_20_False", interest_points.size(), -1);
            data.push_back(data_12_20_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_12_24_False);
            struct benchmark_data data_12_24_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_12_24_False", interest_points.size(), -1);
            data.push_back(data_12_24_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_16_1_False);
            struct benchmark_data data_16_1_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_16_1_False", interest_points.size(), -1);
            data.push_back(data_16_1_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_16_4_False);
            struct benchmark_data data_16_4_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_16_4_False", interest_points.size(), -1);
            data.push_back(data_16_4_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_16_8_False);
            struct benchmark_data data_16_8_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_16_8_False", interest_points.size(), -1);
            data.push_back(data_16_8_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_16_12_False);
            struct benchmark_data data_16_12_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_16_12_False", interest_points.size(), -1);
            data.push_back(data_16_12_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_16_16_False);
            struct benchmark_data data_16_16_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_16_16_False", interest_points.size(), -1);
            data.push_back(data_16_16_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_16_20_False);
            struct benchmark_data data_16_20_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_16_20_False", interest_points.size(), -1);
            data.push_back(data_16_20_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_16_24_False);
            struct benchmark_data data_16_24_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_16_24_False", interest_points.size(), -1);
            data.push_back(data_16_24_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_20_1_False);
            struct benchmark_data data_20_1_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_20_1_False", interest_points.size(), -1);
            data.push_back(data_20_1_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_20_4_False);
            struct benchmark_data data_20_4_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_20_4_False", interest_points.size(), -1);
            data.push_back(data_20_4_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_20_8_False);
            struct benchmark_data data_20_8_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_20_8_False", interest_points.size(), -1);
            data.push_back(data_20_8_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_20_12_False);
            struct benchmark_data data_20_12_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_20_12_False", interest_points.size(), -1);
            data.push_back(data_20_12_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_20_16_False);
            struct benchmark_data data_20_16_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_20_16_False", interest_points.size(), -1);
            data.push_back(data_20_16_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_20_20_False);
            struct benchmark_data data_20_20_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_20_20_False", interest_points.size(), -1);
            data.push_back(data_20_20_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_20_24_False);
            struct benchmark_data data_20_24_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_20_24_False", interest_points.size(), -1);
            data.push_back(data_20_24_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_24_1_False);
            struct benchmark_data data_24_1_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_24_1_False", interest_points.size(), -1);
            data.push_back(data_24_1_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_24_4_False);
            struct benchmark_data data_24_4_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_24_4_False", interest_points.size(), -1);
            data.push_back(data_24_4_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_24_8_False);
            struct benchmark_data data_24_8_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_24_8_False", interest_points.size(), -1);
            data.push_back(data_24_8_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_24_12_False);
            struct benchmark_data data_24_12_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_24_12_False", interest_points.size(), -1);
            data.push_back(data_24_12_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_24_16_False);
            struct benchmark_data data_24_16_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_24_16_False", interest_points.size(), -1);
            data.push_back(data_24_16_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_24_20_False);
            struct benchmark_data data_24_20_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_24_20_False", interest_points.size(), -1);
            data.push_back(data_24_20_False);;
            
            functions.push_back(get_msurf_descriptors_haar_unconditional_unroll_24_24_False);
            struct benchmark_data data_24_24_False(image_name, width, height, "get_msurf_descriptor_haar_unroll_24_24_False", interest_points.size(), -1);
            data.push_back(data_24_24_False);;
            
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
