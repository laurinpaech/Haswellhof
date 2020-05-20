#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os


# In[2]:


top = """
    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size/3;
    int border = (filter_size-1)/2;
    float inv_area = 1.f/(filter_size*filter_size);

    float inv_area_square = inv_area*inv_area;

    int lobe_div_2 = lobe / 2;
    int lobe_sub_1 = lobe - 1;
    int lobe_mul_2_sub_1 = 2*lobe - 1;

    int i;
"""

end = """
    for (; i < height; ++i) {
        int x = i*step;
        for (int j = 0; j < width; ++j) {            
            int y = j*step;

            float Dxx = box_integral_unconditional(iimage, x - lobe_sub_1, y - border, lobe_mul_2_sub_1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe_sub_1, y - lobe_div_2, lobe_mul_2_sub_1, lobe);
            float Dyy = box_integral_unconditional(iimage, x - border, y - lobe_sub_1, filter_size, lobe_mul_2_sub_1)
                    - 3 * box_integral_unconditional(iimage, x - lobe_div_2, y - lobe_sub_1, lobe, lobe_mul_2_sub_1);

            float Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            response[(i*width) + j] = inv_area_square*(Dxx * Dyy - 0.81f * Dxy * Dxy);

            laplacian[(i*width) + j] = (Dxx + Dyy >= 0 ? true : false);
        }
    }
}
"""


# In[3]:


h_include = """#pragma once
#include "fasthessian.h"
#include "integral_image_opt.h"

"""

cpp_include = """#include "fasthessian_opt.h"

"""


# In[4]:


def generate(outer, inner, flip=False):
    code = ""
    offset = ""

    code += f"void compute_response_layer_blocking_{outer}_{inner}_{flip}(struct response_layer* layer, struct integral_image* iimage) {{\n"
    offset += "    "
    code += top
    
    code += f"{offset}for (i = 0; i < height-({outer}-1); i+={outer}) {{\n"
    offset += "    "
    
    for i in range(outer):
        code += f"{offset}int i{i} = i+{i};\n"
    for i in range(outer):
        code += f"{offset}int x{i} = i{i}*step;\n"
    code += f"{offset}int j=0;\n"
    
    code += f"{offset}for (; j < width-({inner}-1); j+={inner}) {{\n"
    offset += "    "
    for j in range(inner):        
        code += f"{offset}int j{j} = j+{j};\n"
    for j in range(inner):        
        code += f"{offset}int y{j} = j{j}*step;\n"
    
    for i in range(outer):
        for j in range(inner):
            code += f"{offset}float Dxx{i}_{j} = box_integral_unconditional(iimage, x{i} - lobe_sub_1, y{j} - border, lobe_mul_2_sub_1, filter_size) - 3 * box_integral_unconditional(iimage, x{i} - lobe_sub_1, y{j} - lobe_div_2, lobe_mul_2_sub_1, lobe);\n"
    for i in range(outer):
        for j in range(inner):
            code += f"{offset}float Dyy{i}_{j} = box_integral_unconditional(iimage, x{i} - border, y{j} - lobe_sub_1, filter_size, lobe_mul_2_sub_1) - 3 * box_integral_unconditional(iimage, x{i} - lobe_div_2, y{j} - lobe_sub_1, lobe, lobe_mul_2_sub_1);\n"
    
    for i in range(outer):
        for j in range(inner):
            code += f"{offset}float Dxy{i}_{j} = box_integral_unconditional(iimage, x{i} - lobe, y{j} + 1, lobe, lobe) + box_integral_unconditional(iimage, x{i} + 1, y{j} - lobe, lobe, lobe) - box_integral_unconditional(iimage, x{i} - lobe, y{j} - lobe, lobe, lobe) - box_integral_unconditional(iimage, x{i} + 1, y{j} + 1, lobe, lobe);\n"
    
    for i in range(outer):
        for j in range(inner):
            code += f"{offset}response[(i{i}*width) + j{j}] = inv_area_square*(Dxx{i}_{j} * Dyy{i}_{j} - 0.81f * Dxy{i}_{j} * Dxy{i}_{j});\n"
            
    for i in range(outer):
        for j in range(inner):
            code += f"{offset}laplacian[(i{i}*width) + j{j}] = (Dxx{i}_{j} + Dyy{i}_{j} >= 0 ? true : false);\n"
            
    offset = offset[:-4]
    code += f"{offset}}}\n"

    code += f"{offset}for (; j < width; ++j) {{\n" 
    offset += "    "
    code += f"{offset}int y = j*step;\n"
    
    for i in range(outer):
        code += f"{offset}float Dxx{i} = box_integral_unconditional(iimage, x{i} - lobe_sub_1, y - border, lobe_mul_2_sub_1, filter_size) - 3 * box_integral_unconditional(iimage, x{i} - lobe_sub_1, y - lobe_div_2, lobe_mul_2_sub_1, lobe);\n"
        code += f"{offset}float Dyy{i} = box_integral_unconditional(iimage, x{i} - border, y - lobe_sub_1, filter_size, lobe_mul_2_sub_1) - 3 * box_integral_unconditional(iimage, x{i} - lobe_div_2, y - lobe_sub_1, lobe, lobe_mul_2_sub_1);\n"

        code += f"{offset}float Dxy{i} = box_integral_unconditional(iimage, x{i} - lobe, y + 1, lobe, lobe) + box_integral_unconditional(iimage, x{i} + 1, y - lobe, lobe, lobe) - box_integral_unconditional(iimage, x{i} - lobe, y - lobe, lobe, lobe) - box_integral_unconditional(iimage, x{i} + 1, y + 1, lobe, lobe);\n"
        
    for i in range(outer):
        code += f"{offset}response[(i{i}*width) + j] = inv_area_square*(Dxx{i} * Dyy{i} - 0.81f * Dxy{i} * Dxy{i});\n"

    for i in range(outer):
        code += f"{offset}laplacian[(i{i}*width) + j] = (Dxx{i} + Dyy{i} >= 0 ? true : false);\n"
    
    offset = offset[:-4]
    code += f"{offset}}}\n"
        
    offset = offset[:-4]
    code += f"{offset}}}\n"

    code += end
    offset = ''
    
    code += f"""\nvoid compute_response_layers_blocking_{outer}_{inner}_{flip}(struct fasthessian* fh) {{
    for (int i = 0; i < 1; ++i) {{//fh->total_layers
        compute_response_layer_blocking_{outer}_{inner}_{flip}(fh->response_map[i], fh->iimage);
    }}  

}}\n"""
    
    return code

def header(outer,inner,flip=False):
    code = f"void compute_response_layer_blocking_{outer}_{inner}_{flip}(struct response_layer* layer, struct integral_image* iimage);\n\n"
    code += f"void compute_response_layers_blocking_{outer}_{inner}_{flip}(struct fasthessian* fh);\n\n"
    return code 

# generate the benchmark code

def generate_benchmark_code(fn):
    benchmark_code = """// has to be defined before stb includes
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
#include "fasthessian_opt_gen.h"
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
        printf("%s\\n", image_name);
        float *image = stbi_loadf(image_name, &width, &height, &channels, STBI_grey);
        if (!image) {
            printf("Could not open or find image\\n");
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
"""

    for outer in outer_range:
        for inner in inner_range:
            for flip in FLIPS:

                benchmark_code += f"""
        functions.push_back(compute_response_layers_blocking_{outer}_{inner}_{flip});
        struct benchmark_data data_{outer}_{inner}_{flip}(image_name, width, height, (char *)"compute_response_layers_blocking_{outer}_{inner}_{flip}", -1, (1 + height * width * 13));
        data.push_back(data_{outer}_{inner}_{flip});;
        """
            
    benchmark_code += """
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
    printf("Benchmarking done!\\n");

    return 0;
}
"""

    with open(os.path.join(fp,fn), "w") as file:
        file.write(benchmark_code)


def generate_cpp_code(fn_cpp, fn_h):
    header_txt = h_include
    cpp_txt = cpp_include
    for outer in outer_range:
        for inner in inner_range:
            for flip in FLIPS:
                header_txt += header(outer,inner,flip)
                cpp_txt += generate(outer,inner,flip)

    with open(os.path.join(fp,fn_h), "w") as file:
        file.write(header_txt)

    with open(os.path.join(fp,fn_cpp), "w") as file:
        file.write(cpp_txt)
            


# In[5]:


# %cd ../build

os.chdir("../build")


# In[ ]:


from tqdm import tqdm
import subprocess

fp = "/Users/valentinwolf/Documents/Studium/ETH/ASL/team042/src/autotune/"
outer_range = [3,7] #+ list(range(4,25,4))#[1,4,8,12,16,20,24]
inner_range = list(range(4,33,4))#[1,4,8,12,16,20,24]
FLIPS = [False]#[True,False]

for i in tqdm(range(4,33,4)):
    outer_range = [i]
    generate_cpp_code('fasthessian_opt_gen.cpp', 'fasthessian_opt_gen.h')
    generate_benchmark_code('main_autotune.cpp')
#     !make -j7 autotune 
#     !./src/autotune

    subprocess.call(["make", "-j7"])
    subprocess.call(["./src/autotune"])

        


# for i in tqdm(range(33)):
#     outer_range = [i]
# 

# In[ ]:





# In[ ]:




