// has to be defined before stb includes
#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"
#include "integral_image.h"
#include "fasthessian.h"
#include "interest_point.h"
#include "descriptor.h"
#include "benchmarking.h"
#include "benchmark_interpolate_step.h"
#include "benchmark_data_to_file.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>

const char* images[] = {
    "img1",
    "img2",
};

#define n_images (sizeof (images) / sizeof (const char *))

int main(int argc, char const *argv[])
{
    char file_ending[16] = ".png";
    std::vector<struct benchmark_data*> all_benchmark_data;

    for (int i = 0; i < n_images; i++) {
        char image_name[32];
        strcpy(image_name, images[i]);

        int width, height, channels;

        // Load image
        stbi_ldr_to_hdr_gamma(1.0f);
        char* path_name = concat("../images/", image_name);
        path_name = concat(path_name, file_ending);
        printf("%s\n",path_name);
        float* image = stbi_loadf(path_name, &width, &height, &channels, STBI_grey);

        if(!image) {
            printf("Could not open or find image\n");
            return -1;
        }

        // Calculate integral image
        printf("integral start\n");
        struct integral_image* iimage = create_integral_img(image, width, height);

        struct benchmark_data* benchmark_integral_img=initialise_benchmark_data(image_name, width, height, "create_integral_img", -1, (width + 2*(height-1)*width));
        perf_test_integral_img(create_integral_img, image, benchmark_integral_img);
        all_benchmark_data.push_back(benchmark_integral_img);
        printf("integral done\n");
        
        // Fast-Hessian
        struct fasthessian* fh = create_fast_hessian(iimage);

        // Create octaves with response layers
        create_response_map(fh);
        
        printf("compute_response_layer start\n");
        /**
        struct benchmark_data* benchmark_compute_response_layer=initialise_benchmark_data(image_name, width, height, "compute_response_layer", -1, (1 + height*width*13));
        perf_test_compute_response_layer(compute_response_layer, fh->response_map[0], iimage, benchmark_compute_response_layer);
        all_benchmark_data.push_back(benchmark_compute_response_layer);
        **/
        
        // Compute responses for every layer
        for (size_t i = 0; i < fh->total_layers; i++) {
            compute_response_layer(fh->response_map[i], iimage);           
        }
        printf("compute_response_layer end\n");

        // Getting interest points with non-maximum supression
        std::vector<struct interest_point> interest_points;
        get_interest_points(fh, &interest_points);


        /**
        printf("get_interest_points start\n");       
        long flops = 109*interest_points.size();
        struct benchmark_data* benchmark_get_interest_points=initialise_benchmark_data(image_name, width, height, "get_interest_points", interest_points.size(), flops);
        perf_test_get_interest_points(get_interest_points, fh, benchmark_get_interest_points);
        all_benchmark_data.push_back(benchmark_get_interest_points);
        printf("get_interest_points end\n");
        **/

        printf("interpolate_step start\n");

        struct benchmark_data* benchmark_interpolate_step=initialise_benchmark_data(image_name, width, height, "interpolate_step", interest_points.size(), 109);
        bench_interpolate_step(fh, &interest_points, benchmark_interpolate_step);
        all_benchmark_data.push_back(benchmark_interpolate_step);
        printf("interpolate_step end\n");
        

    /**
        // Descriptor stuff
        static float* GW = get_gaussian(3.3);
        for (int i=0; i<interest_points.size(); i++)
            get_descriptor(iimage, &interest_points[i], GW);

        free(GW);

        // Write results to file
        FILE * fp = fopen("desc.txt","w");
        for (int i=0; i<interest_points.size(); i++) {
            fprintf(fp, "%f %f ", interest_points[i].x, interest_points[i].y);
            // printf("%f ", BLOB_ORIENTATION); // TODO: how to get this value
            for(int j = 0; j < 64; j++) {
                fprintf(fp, "%f ", interest_points[i].descriptor[j]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    */

        // Free memory
        stbi_image_free(image); // possibly move this to create_integral_img
        free(iimage->data);
        free(iimage);
        
        for (size_t i = 0; i < NUM_LAYERS; i++) {
            free(fh->response_map[i]);
        }
        free(fh);
        
    }
    printf("write tor file start\n");
    save_benchmark_data(all_benchmark_data);
    // free memory benchmarkdata
    // https://stackoverflow.com/questions/10464992/c-delete-vector-objects-free-memory
    std::vector<struct benchmark_data*>().swap(all_benchmark_data);

    return 0;
}
