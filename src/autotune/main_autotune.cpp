// has to be defined before stb includes
#define STB_IMAGE_IMPLEMENTATION

#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "benchmark_data_to_file.h"
#include "benchmarking.h"
#include "descriptor.h"
#include "descriptor_opt.h"
#include "fasthessian.h"
#include "fasthessian_opt.h"
// #include "fasthessian_opt_gen.h"
#include "integral_image.h"
#include "interest_point.h"
#include "stb_image.h"

const char *images[] = {
    // "../images/sunflower/sunflower_32.jpg", 
    // "../images/sunflower/sunflower_64.jpg",
    //"../images/sunflower/sunflower_128.jpg",
    //"../images/sunflower/sunflower_256.jpg",
    "../images/sunflower/sunflower_512.jpg",
    //"../images/sunflower/sunflower_1024.jpg",
    //"../images/sunflower/sunflower_2048.jpg"
    //"../images/sunflower/sunflower_4096.jpg"
};
#define n_images (sizeof(images) / sizeof(const char *))
#define BENCHMARK_COMPUTE_RESPONSE_LAYERS
// BENCHMARK_COMPUTE_RESPONSE_LAYERS_PADDED only works with BENCHMARK_COMPUTE_RESPONSE_LAYERS enabled
#define BENCHMARK_COMPUTE_RESPONSE_LAYERS_PADDED

int main(int argc, char const *argv[]) {
    std::vector<struct benchmark_data> all_benchmark_data;
    for (int i = 0; i < n_images; i++) {
        char *image_name = (char *)malloc(512 * sizeof(char));
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

        // Create padded integral image
        struct integral_image *padded_iimage = create_padded_integral_img(width, height);
        // Compute padded integral image
        compute_padded_integral_img(image, padded_iimage);
        
        std::vector<void (*)(struct fasthessian *)> functions;
        std::vector<struct benchmark_data> data;
        /*
        functions.push_back(compute_response_layers_blocking_20_4_False);
        struct benchmark_data data_20_4_False(image_name, width, height, (char *)"compute_response_layers_blocking_20_4_False", -1, (1 + height * width * 13));
        data.push_back(data_20_4_False);;
        
        functions.push_back(compute_response_layers_blocking_20_8_False);
        struct benchmark_data data_20_8_False(image_name, width, height, (char *)"compute_response_layers_blocking_20_8_False", -1, (1 + height * width * 13));
        data.push_back(data_20_8_False);;
        
        functions.push_back(compute_response_layers_blocking_20_12_False);
        struct benchmark_data data_20_12_False(image_name, width, height, (char *)"compute_response_layers_blocking_20_12_False", -1, (1 + height * width * 13));
        data.push_back(data_20_12_False);;
        
        functions.push_back(compute_response_layers_blocking_20_16_False);
        struct benchmark_data data_20_16_False(image_name, width, height, (char *)"compute_response_layers_blocking_20_16_False", -1, (1 + height * width * 13));
        data.push_back(data_20_16_False);;
        
        functions.push_back(compute_response_layers_blocking_20_20_False);
        struct benchmark_data data_20_20_False(image_name, width, height, (char *)"compute_response_layers_blocking_20_20_False", -1, (1 + height * width * 13));
        data.push_back(data_20_20_False);;
        
        functions.push_back(compute_response_layers_blocking_20_24_False);
        struct benchmark_data data_20_24_False(image_name, width, height, (char *)"compute_response_layers_blocking_20_24_False", -1, (1 + height * width * 13));
        data.push_back(data_20_24_False);;
        
        functions.push_back(compute_response_layers_blocking_20_28_False);
        struct benchmark_data data_20_28_False(image_name, width, height, (char *)"compute_response_layers_blocking_20_28_False", -1, (1 + height * width * 13));
        data.push_back(data_20_28_False);;
        
        functions.push_back(compute_response_layers_blocking_20_32_False);
        struct benchmark_data data_20_32_False(image_name, width, height, (char *)"compute_response_layers_blocking_20_32_False", -1, (1 + height * width * 13));
        data.push_back(data_20_32_False);;
        */
        
        bench_compute_response_layer(functions, padded_iimage, data);

        all_benchmark_data.insert(all_benchmark_data.end(), data.begin(), data.end());
            
        // Compute responses for every layer
        compute_response_layers(fh);

        // Getting interest points with non-maximum supression
        std::vector<struct interest_point> interest_points;
        get_interest_points(fh, &interest_points);


        // Getting M-SURF descriptors for each interest point
        get_msurf_descriptors(iimage, &interest_points);

        // Free memory
        stbi_image_free(image);  // possibly move this to create_integral_img

        free(iimage->padded_data);
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
    // std::vector<struct benchmark_data *>().swap(all_benchmark_data);
    printf("Benchmarking done!\n");

    return 0;
}
