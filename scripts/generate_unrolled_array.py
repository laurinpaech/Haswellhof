#!/usr/bin/env python
# coding: utf-8

# In[2]:


top = """
    float scale = ipoint->scale;
    int int_scale = (int) roundf(scale);
    float scale_squared = scale*scale;
    float g1_factor = -0.08f / (scale_squared); 

    float ipoint_x = roundf(ipoint->x) + 0.5*scale;
    float ipoint_y = roundf(ipoint->y) + 0.5*scale;

    float ipoint_x_sub_int_scale = ipoint_x-int_scale;
    float ipoint_y_sub_int_scale = ipoint_y-int_scale;

    float ipoint_x_sub_int_scale_add_05 = ipoint_x-int_scale + 0.5;
    float ipoint_y_sub_int_scale_add_05 = ipoint_y-int_scale + 0.5;
    
    int width = iimage->width;
    int height = iimage->height;

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;

    // Initializing gauss_s2 index for precomputed array
    int gauss_s2_index = 0;

    // check if we ever hit a boundary
    if (((int) roundf(ipoint_x - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_y - 12*scale)) - int_scale <= 0 
        || ((int) roundf(ipoint_x + 11*scale)) + int_scale > width 
        || ((int) roundf(ipoint_y + 11*scale)) + int_scale > height) 
    {
"""

end = """
    float s0  = roundf( 0.5 * scale);
    float s1  = roundf( 1.5 * scale);
    float s2  = roundf( 2.5 * scale);
    float s3  = roundf( 3.5 * scale);
    float s4  = roundf( 4.5 * scale);
    float s5  = roundf( 5.5 * scale);
    float s6  = roundf( 6.5 * scale);
    float s7  = roundf( 7.5 * scale);
    float s8  = roundf( 8.5 * scale);
    float s9  = roundf( 9.5 * scale);
    float s10 = roundf(10.5 * scale);
    float s11 = roundf(11.5 * scale);

    float e_c0_m4 = s2 + s1; // CAREFUL HERE!
    float e_c0_m3 = s2 + s0; // CAREFUL HERE!
    float e_c0_m2 = s2 - s0;
    float e_c0_m1 = s2 - s1;
    //float e_c0_z0 = s2 - s2;
    float e_c0_p1 = s2 - s3;
    float e_c0_p2 = s2 - s4;
    float e_c0_p3 = s2 - s5;
    float e_c0_p4 = s2 - s6;

    float e_c1_m4 = s7 - s3;
    float e_c1_m3 = s7 - s4;
    float e_c1_m2 = s7 - s5;
    float e_c1_m1 = s7 - s6;
    //float e_c1_z0 = s7 - s7;
    float e_c1_p1 = s7 - s8;
    float e_c1_p2 = s7 - s9;
    float e_c1_p3 = s7 - s10;
    float e_c1_p4 = s7 - s11;

    gauss_s1_c0[0] =  expf(g1_factor * (e_c0_m4 * e_c0_m4));
    gauss_s1_c0[1] =  expf(g1_factor * (e_c0_m3 * e_c0_m3));
    gauss_s1_c0[2] =  expf(g1_factor * (e_c0_m2 * e_c0_m2));
    gauss_s1_c0[3] =  expf(g1_factor * (e_c0_m1 * e_c0_m1));
    gauss_s1_c0[4] =  1.0f; //expf(g1_factor * (e_c0_z0 * e_c0_z0));
    gauss_s1_c0[5] =  expf(g1_factor * (e_c0_p1 * e_c0_p1));
    gauss_s1_c0[6] =  expf(g1_factor * (e_c0_p2 * e_c0_p2));
    gauss_s1_c0[7] =  expf(g1_factor * (e_c0_p3 * e_c0_p3));
    gauss_s1_c0[8] =  expf(g1_factor * (e_c0_p4 * e_c0_p4));

    gauss_s1_c1[0] =  expf(g1_factor * (e_c1_m4 * e_c1_m4));
    gauss_s1_c1[1] =  expf(g1_factor * (e_c1_m3 * e_c1_m3));
    gauss_s1_c1[2] =  expf(g1_factor * (e_c1_m2 * e_c1_m2));
    gauss_s1_c1[3] =  expf(g1_factor * (e_c1_m1 * e_c1_m1));
    gauss_s1_c1[4] =  1.0f; //expf(g1_factor * (e_c1_z0 * e_c1_z0));
    gauss_s1_c1[5] =  expf(g1_factor * (e_c1_p1 * e_c1_p1));
    gauss_s1_c1[6] =  expf(g1_factor * (e_c1_p2 * e_c1_p2));
    gauss_s1_c1[7] =  expf(g1_factor * (e_c1_p3 * e_c1_p3));
    gauss_s1_c1[8] =  expf(g1_factor * (e_c1_p4 * e_c1_p4));
    
        // calculate descriptor for this interest point
    for (int i=-8; i<8; i+=5) {

        for (int j=-8; j<8; j+=5) {

            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            int gauss_index_l = -4;
            for (int l = j-4; l < j + 5; ++l, ++gauss_index_l) {
                float gauss_s1_y = -1;

                if (j == -8 ) {
                    gauss_s1_y = gauss_s1_c1[8-(gauss_index_l+4)];
                } else if (j == -3) {
                    gauss_s1_y = gauss_s1_c0[8-(gauss_index_l+4)];
                } else if (j == 2) {
                    gauss_s1_y = gauss_s1_c0[gauss_index_l+4];
                } else if (j == 7) {
                    gauss_s1_y = gauss_s1_c1[gauss_index_l+4];
                }

                int gauss_index_k = -4;
                for (int k = i-4; k < i + 5; ++k, ++gauss_index_k) {

                    float gauss_s1_x = -1;
                    if (i == -8 ) {
                        gauss_s1_x = gauss_s1_c1[8-(gauss_index_k+4)];
                    } else if (i == -3) {
                        gauss_s1_x = gauss_s1_c0[8-(gauss_index_k+4)];
                    } else if (i == 2) {
                        gauss_s1_x = gauss_s1_c0[gauss_index_k+4];
                    } else if (i == 7) {
                        gauss_s1_x = gauss_s1_c1[gauss_index_k+4];
                    }

                    float gauss_s1 = gauss_s1_x * gauss_s1_y;

                    float rx = haarResponseX[(l+12)*24+(k+12)];
                    float ry = haarResponseY[(l+12)*24+(k+12)];
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    float rrx = gauss_s1 * ry;
                    float rry = gauss_s1 * rx;

                    dx += rrx;
                    dy += rry;
                    mdx += fabsf(rrx);
                    mdy += fabsf(rry);
                }
            }

            // Precomputed 4x4 gauss_s2 with (x,y) = {-1.5, -0.5, 0.5, 1.5}^2 and sig = 1.5f
            float gauss_s2 = gauss_s2_arr[gauss_s2_index];
            gauss_s2_index++;

            // add the values to the descriptor vector
            float d1 = dx * gauss_s2;
            float d2 = dy * gauss_s2;
            float d3 = mdx * gauss_s2;
            float d4 = mdy * gauss_s2;

            descriptor[desc_idx] = d1;
            descriptor[desc_idx+1] = d2;
            descriptor[desc_idx+2] = d3;
            descriptor[desc_idx+3] = d4;

            // precompute for normaliztion
            sum_of_squares += (d1*d1 + d2*d2 + d3*d3 + d4*d4);

            desc_idx += 4;

        }
    }

    // rescale to unit vector
    float norm_factor = 1.0f / sqrtf(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }
}
"""


# In[3]:


h_include = """#include "integral_image.h"
#include "interest_point.h"
#include "helper.h"

#include <stdlib.h>
#include <vector>

"""

cpp_include = """#include "descriptor.h"
#include "descriptor_opt.h"
#include "descriptor_opt_gen.h"
#include "integral_image.h"
#include "interest_point.h"

#define _USE_MATH_DEFINES

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>

extern float * haarResponseX;
extern float * haarResponseY;

extern float gauss_s1_c0[9];
extern float gauss_s1_c1[9];
extern const float gauss_s2_arr[16];

"""


# In[4]:


def generate(outer, inner, flip=False):
    code = ""
    offset = ""

    code += f"void get_msurf_descriptor_haar_unroll_{outer}_{inner}_{flip}(struct integral_image* iimage, struct interest_point* ipoint) {{\n"
    code += top
    offset += "    "
    offset += "    "

    code += f"{offset}for (int l=-12, l_count=0; l<12; l+={outer}, l_count+={outer}) {{\n"
    offset += "    "

    for i in range(outer):
        code += f"{offset}int l{i} = l + {i};\n"
    for i in range(outer):
        code += f"{offset}int l_count{i} = l_count + {i};\n"
        
    for i in range(outer):
        code += f"{offset}float ipoint_y_sub_int_scale_add_l{i}_mul_scale = ipoint_y_sub_int_scale + l{i} * scale;\n"
        
    for i in range(outer):
        code += f"{offset}int sample_y_sub_int_scale{i} = (int) (ipoint_y_sub_int_scale_add_l{i}_mul_scale + (ipoint_y_sub_int_scale_add_l{i}_mul_scale>=0 ? 0.5 : -0.5));\n"
    code += "\n"
    code += f"{offset}for (int k=-12, k_count=0; k<12; k+={inner}, k_count+={inner}) {{\n"
    offset += "    "

    for i in range(inner):
        code += f"{offset}int k{i} = k + {i};\n"
    for i in range(inner):
        code += f"{offset}int k_count{i} = k_count + {i};\n"
    for i in range(inner):
        code += f"{offset}float ipoint_x_sub_int_scale_add_k{i}_mul_scale = ipoint_x_sub_int_scale + k{i} * scale;\n"
    for i in range(inner):
        code += f"{offset}int sample_x_sub_int_scale{i} = (int) (ipoint_x_sub_int_scale_add_k{i}_mul_scale + (ipoint_x_sub_int_scale_add_k{i}_mul_scale>=0 ? 0.5 : -0.5));\n"

    code += "\n"
    
    
    if flip:
        for j in range(inner):
            for i in range(outer):
                code += f"{offset}haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale{i}, sample_x_sub_int_scale{j}, int_scale, &haarResponseX[l_count{i}*24+k_count{j}], &haarResponseY[l_count{i}*24+k_count{j}]);\n"
    else:
        for i in range(outer):
            for j in range(inner):
                code += f"{offset}haarXY_precheck_boundaries(iimage, sample_y_sub_int_scale{i}, sample_x_sub_int_scale{j}, int_scale, &haarResponseX[l_count{i}*24+k_count{j}], &haarResponseY[l_count{i}*24+k_count{j}]);\n"

    offset = offset[:-4]
    code += f"{offset}}}\n"

    offset = offset[:-4]
    code += f"{offset}}}\n"

    offset = offset[:-4]

    code += f"{offset}}} else {{\n"
    offset += "    "

    code += f"{offset}for (int l=-12, l_count=0; l<12; l+={outer}, l_count+={outer}) {{\n"
    offset += "    "

    for i in range(outer):
        code += f"{offset}int l{i} = l + {i};\n"
    for i in range(outer):
        code += f"{offset}int l_count{i} = l_count + {i};\n"
    for i in range(outer):
        code += f"{offset}int sample_y_sub_int_scale{i} = (int) (ipoint_y_sub_int_scale_add_05 + l{i} * scale);\n"
    code += "\n"
    code += f"{offset}for (int k=-12, k_count=0; k<12; k+={inner}, k_count+={inner}) {{\n"
    offset += "    "

    for i in range(inner):
        code += f"{offset}int k{i} = k + {i};\n"
    for i in range(inner):
        code += f"{offset}int k_count{i} = k_count + {i};\n"
    for i in range(inner):
        code += f"{offset}int sample_x_sub_int_scale{i} = (int) (ipoint_x_sub_int_scale_add_05 + k{i} * scale);\n"

    code += "\n"

    if flip:
        for j in range(inner):
            for i in range(outer):
                code += f"{offset}haarXY_unconditional(iimage, sample_y_sub_int_scale{i}, sample_x_sub_int_scale{j}, int_scale, &haarResponseX[l_count{i}*24+k_count{j}], &haarResponseY[l_count{i}*24+k_count{j}]);\n"
    else:
        for i in range(outer):
            for j in range(inner):
                code += f"{offset}haarXY_unconditional(iimage, sample_y_sub_int_scale{i}, sample_x_sub_int_scale{j}, int_scale, &haarResponseX[l_count{i}*24+k_count{j}], &haarResponseY[l_count{i}*24+k_count{j}]);\n"

    offset = offset[:-4]
    code += f"{offset}}}\n"

    offset = offset[:-4]
    code += f"{offset}}}\n"

    offset = offset[:-4]

    code += f"{offset}}}\n"

    code += end
    offset = ''
    
    code += f"""\nvoid get_msurf_descriptors_haar_unroll_{outer}_{inner}_{flip}(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {{
    for (size_t i=0; i<interest_points->size(); ++i) {{
        get_msurf_descriptor_haar_unroll_{outer}_{inner}_{flip}(iimage, &interest_points->at(i));
    }}
}}\n"""
    
    return code

def header(outer,inner,flip=False):
    code = f"void get_msurf_descriptor_haar_unroll_{outer}_{inner}_{flip}(struct integral_image* iimage, struct interest_point* ipoint);\n\n"
    code += f"void get_msurf_descriptors_haar_unroll_{outer}_{inner}_{flip}(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);\n\n"
    return code 

# generate the benchmark code
def generate_benchmark_code(outer_range, inner_range, FLIPS):
    benchmark_code = """// has to be defined before stb includes
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
    "../images/sunflower/sunflower_64.jpg",
    "../images/sunflower/sunflower_128.jpg", 
    "../images/sunflower/sunflower_256.jpg",
    "../images/sunflower/sunflower_512.jpg",
    "../images/sunflower/sunflower_1024.jpg",
    // "../images/sunflower/sunflower_2048.jpg",
    // "../images/sunflower/sunflower_4096.jpg",
};

#define n_images (sizeof(images) / sizeof(const char *))

int main(int argc, char const *argv[]) {
    std::vector<struct benchmark_data> all_benchmark_data;
    initialize_folder_name();
    for (int i = 0; i < n_images; i++) {
        char *image_name = (char *)malloc(1024 * sizeof(char));
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
        
        // Compute responses for every layer
        compute_response_layers(fh);

        // Getting interest points with non-maximum supression
        std::vector<struct interest_point> interest_points;
        get_interest_points(fh, &interest_points);

        {
            printf("get_msurf_descriptor start\\n");
            
            // Insert all interpolate_step functions for benchmarking here
            std::vector<void (*)(struct integral_image *, std::vector<struct interest_point> *)> functions;
            // Insert all respective benchmarking info for functions here
            std::vector<struct benchmark_data> data;"""

    for outer in outer_range:
        for inner in inner_range:
            for flip in FLIPS:

                benchmark_code += f"""
            functions.push_back(get_msurf_descriptors_haar_unroll_{outer}_{inner}_{flip});
            struct benchmark_data data_{outer}_{inner}_{flip}(image_name, width, height, "get_msurf_descriptor_haar_unroll_{outer}_{inner}_{flip}", interest_points.size(), -1);
            data.push_back(data_{outer}_{inner}_{flip});;
            """

    benchmark_code += """
            // Benchmarking all get_msurf_descriptor functions and storing timing results in respective entries in data
            bench_get_msurf_descriptors(functions, iimage, &interest_points, data);

            // Appending this data to all benchmarking data
            all_benchmark_data.insert(all_benchmark_data.end(), data.begin(), data.end());

            printf("get_msurf_descriptor end\\n");
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
    
    // extern float* haarResponseX;
    // extern float* haarResponseY;

    // aligned_free(haarResponseX);
    // aligned_free(haarResponseY);
    
    save_benchmark_data(all_benchmark_data);
    // free memory benchmarkdata
    // https://stackoverflow.com/questions/10464992/c-delete-vector-objects-free-memory
    //std::vector<struct benchmark_data *>().swap(all_benchmark_data);
    printf("Benchmarking done!\\n");

    return 0;
}
"""
    return benchmark_code


# In[5]:


if __name__ == "__main__":
    import os
    import sys

    outer_range = [1, 2, 3, 4, 6, 8, 12, 24]
    print("generating all with outer unrolling=", outer_range)
    # assert(24%int(sys.argv[1]) == 0)
    inner_range = [1, 2, 3, 4, 6, 8, 12, 24]
    FLIPS = [True,False]

    header_txt = h_include
    cpp_txt = cpp_include
    for outer in outer_range:
        for inner in inner_range:
            for flip in FLIPS:
                header_txt += header(outer,inner,flip)
                cpp_txt += generate(outer,inner,flip)

    benchmark_txt = generate_benchmark_code(outer_range, inner_range, FLIPS)

    fp = "/Users/valentinwolf/Documents/Studium/ETH/ASL/team042/src/autotune/"
    with open(fp+'descriptor_opt_gen.h', "w") as file:
        file.write(header_txt)

    with open(fp+'descriptor_opt_gen.cpp', "w") as file:
        file.write(cpp_txt)

    with open(os.path.join(fp,'main_autotune.cpp'), "w") as file:
        file.write(benchmark_txt)

