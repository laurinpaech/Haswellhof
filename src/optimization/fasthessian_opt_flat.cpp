#include "fasthessian_opt_flat.h"

#include "helper.h"

#include <stdlib.h>
#include <assert.h>

#include <iostream>

void create_fast_hessian_flat_and_response_map(struct integral_image *iimage, struct fasthessian_flat *fh_flat) {

    assert(fh_flat != NULL);

    // Interest Image
    fh_flat->iimage = iimage;

    // Set Variables
    fh_flat->octaves = NUM_OCTAVES;
    fh_flat->layers = NUM_LAYERS;
    fh_flat->total_layers = NUM_TOTAL_LAYERS;
    fh_flat->step = INITIAL_STEP;
    fh_flat->thresh = THRESHOLD;

    int img_width = iimage->width;
    int img_height = iimage->height;
    int init_step = INITIAL_STEP;

    int w = img_width / init_step;
    int h = img_height / init_step;

    std::cout << "sizeof(bool): " << sizeof(bool) << std::endl;

    size_t step_o0_r = ALIGN(w * h * sizeof(float), 64);
    size_t step_o0_l = ALIGN(w * h * sizeof(bool), 64);
    size_t step_o0_rl = step_o0_r + step_o0_l;

    size_t step_o1_r = ALIGN(w/2 * h/2 * sizeof(float), 64);
    size_t step_o1_l = ALIGN(w/2 * h/2 * sizeof(bool), 64);
    size_t step_o1_rl = step_o1_r + step_o1_l;

    size_t step_o2_r = ALIGN(w/4 * h/4 * sizeof(float), 64);
    size_t step_o2_l = ALIGN(w/4 * h/4 * sizeof(bool), 64);
    size_t step_o2_rl = step_o2_r + step_o2_l;

/*
    const size_t offsets[] = {
        0*step_o0_r + 0*step_o0_l, 
        1*step_o0_r + 0*step_o0_l,
        1*step_o0_r + 1*step_o0_l,
        2*step_o0_r + 1*step_o0_l,
        2*step_o0_r + 2*step_o0_l,
        3*step_o0_r + 2*step_o0_l,
        3*step_o0_r + 3*step_o0_l,
        4*step_o0_r + 3*step_o0_l,
        4*step_o0_r + 4*step_o0_l + 0*step_o1_r + 0*step_o1_l,
        4*step_o0_r + 4*step_o0_l + 1*step_o1_r + 0*step_o1_l,
        4*step_o0_r + 4*step_o0_l + 1*step_o1_r + 1*step_o1_l,
        4*step_o0_r + 4*step_o0_l + 2*step_o1_r + 1*step_o1_l,
        4*step_o0_r + 4*step_o0_l + 2*step_o1_r + 2*step_o1_l + 0*step_o2_r + 0*step_o2_l,
        4*step_o0_r + 4*step_o0_l + 2*step_o1_r + 2*step_o1_l + 1*step_o2_r + 0*step_o2_l,
        4*step_o0_r + 4*step_o0_l + 2*step_o1_r + 2*step_o1_l + 1*step_o2_r + 1*step_o2_l,
        4*step_o0_r + 4*step_o0_l + 2*step_o1_r + 2*step_o1_l + 2*step_o2_r + 1*step_o2_l
    };
*/

    const size_t offsets[] = {
        0, 
        step_o0_r,
        step_o0_rl,
        step_o0_rl + step_o0_r,
        2*step_o0_rl,
        2*step_o0_rl + step_o0_r,
        3*step_o0_rl,
        3*step_o0_rl + step_o0_r,
        4*step_o0_rl,
        4*step_o0_rl + step_o1_r,
        4*step_o0_rl + step_o1_rl,
        4*step_o0_rl + step_o1_rl + step_o1_r,
        4*step_o0_rl + 2*step_o1_rl,
        4*step_o0_rl + 2*step_o1_rl + step_o2_r,
        4*step_o0_rl + 2*step_o1_rl + step_o2_rl,
        4*step_o0_rl + 2*step_o1_rl + step_o2_rl + step_o2_r
    };

    char *mem = (char *) aligned_malloc(4*step_o0_rl + 2*step_o1_rl + 2*step_o2_rl, 64);

    char *m0_r = mem + offsets[0];
    char *m0_l = mem + offsets[1];
    char *m1_r = mem + offsets[2];
    char *m1_l = mem + offsets[3];
    char *m2_r = mem + offsets[4];
    char *m2_l = mem + offsets[5];
    char *m3_r = mem + offsets[6];
    char *m3_l = mem + offsets[7];
    char *m4_r = mem + offsets[8];
    char *m4_l = mem + offsets[9];
    char *m5_r = mem + offsets[10];
    char *m5_l = mem + offsets[11];
    char *m6_r = mem + offsets[12];
    char *m6_l = mem + offsets[13];
    char *m7_r = mem + offsets[14];
    char *m7_l = mem + offsets[15];

 /*   
    // Computing number of total elements in response
    size_t num_elements = 4*(w * h) + 2*(w/2 * h/2) + 2*(w/4 * h/4);

    char *mem = (char *) aligned_malloc(num_elements * (sizeof(float) + sizeof(bool)));

    size_t offset = 0;
    char *m0_r = mem;           offset += (w * h) * sizeof(float);
    char *m0_l = mem + offset;  offset += (w * h) * sizeof(bool);

    char *m1_r = mem + offset;  offset += (w * h) * sizeof(float);
    char *m1_l = mem + offset;  offset += (w * h) * sizeof(bool);
    
    char *m2_r = mem + offset;  offset += (w * h) * sizeof(float);
    char *m2_l = mem + offset;  offset += (w * h) * sizeof(bool);

    char *m3_r = mem + offset;  offset += (w * h) * sizeof(float);
    char *m3_l = mem + offset;  offset += (w * h) * sizeof(bool);

    char *m4_r = mem + offset;  offset += (w/2 * h/2) * sizeof(float);
    char *m4_l = mem + offset;  offset += (w/2 * h/2) * sizeof(bool);

    char *m5_r = mem + offset;  offset += (w/2 * h/2) * sizeof(float);
    char *m5_l = mem + offset;  offset += (w/2 * h/2) * sizeof(bool);

    char *m6_r = mem + offset;  offset += (w/4 * h/4) * sizeof(float);
    char *m6_l = mem + offset;  offset += (w/4 * h/4) * sizeof(bool);

    char *m7_r = mem + offset;  offset += (w/4 * h/4) * sizeof(float);
    char *m7_l = mem + offset;  offset += (w/4 * h/4) * sizeof(bool);
*/

    {
        fh_flat->response_map[0].filter_size = 9;
        fh_flat->response_map[0].width = w;
        fh_flat->response_map[0].height = h;
        fh_flat->response_map[0].step = init_step;
        fh_flat->response_map[0].response = (float *) m0_r;
        fh_flat->response_map[0].laplacian = (bool *) m0_l;
    }

    {
        fh_flat->response_map[1].filter_size = 15;
        fh_flat->response_map[1].width = w;
        fh_flat->response_map[1].height = h;
        fh_flat->response_map[1].step = init_step;
        fh_flat->response_map[1].response = (float *) m1_r;
        fh_flat->response_map[1].laplacian = (bool *) m1_l;
    }

    {
        fh_flat->response_map[2].filter_size = 21;
        fh_flat->response_map[2].width = w;
        fh_flat->response_map[2].height = h;
        fh_flat->response_map[2].step = init_step;
        fh_flat->response_map[2].response = (float *) m2_r;
        fh_flat->response_map[2].laplacian = (bool *) m2_l;
    }

    {
        fh_flat->response_map[3].filter_size = 27;
        fh_flat->response_map[3].width = w;
        fh_flat->response_map[3].height = h;
        fh_flat->response_map[3].step = init_step;
        fh_flat->response_map[3].response = (float *) m3_r;
        fh_flat->response_map[3].laplacian = (bool *) m3_l;
    }

    {
        fh_flat->response_map[4].filter_size = 39;
        fh_flat->response_map[4].width = w / 2;
        fh_flat->response_map[4].height = h / 2;
        fh_flat->response_map[4].step = init_step * 2;
        fh_flat->response_map[4].response = (float *) m4_r;
        fh_flat->response_map[4].laplacian = (bool *) m4_l;
    }

    {
        fh_flat->response_map[5].filter_size = 51;
        fh_flat->response_map[5].width = w / 2;
        fh_flat->response_map[5].height = h / 2;
        fh_flat->response_map[5].step = init_step * 2;
        fh_flat->response_map[5].response = (float *) m5_r;
        fh_flat->response_map[5].laplacian = (bool *) m5_l;
    }

    {
        fh_flat->response_map[6].filter_size = 75;
        fh_flat->response_map[6].width = w / 4;
        fh_flat->response_map[6].height = h / 4;
        fh_flat->response_map[6].step = init_step * 4;
        fh_flat->response_map[6].response = (float *) m6_r;
        fh_flat->response_map[6].laplacian = (bool *) m6_l;
    }

    {
        fh_flat->response_map[7].filter_size = 99;
        fh_flat->response_map[7].width = w / 4;
        fh_flat->response_map[7].height = h / 4;
        fh_flat->response_map[7].step = init_step * 4;
        fh_flat->response_map[7].response = (float *) m7_r;
        fh_flat->response_map[7].laplacian = (bool *) m7_l;
    }

}

void compute_response_map_flat(struct fasthessian_flat* fh_flat) {
    
    int total_layers = fh_flat->total_layers;
    for (int i = 0; i < total_layers; ++i) {
		compute_response_layer(&fh_flat->response_map[i], fh_flat->iimage);
	} 

}

void get_interest_points_flat(struct fasthessian_flat *fh_flat, std::vector<struct interest_point> *interest_points) {

    assert(fh_flat != NULL);
    assert(interest_points != NULL);

    // filter index map
    const int filter_map[NUM_OCTAVES][NUM_LAYERS] = {
        {0, 1, 2, 3},
        {1, 3, 4, 5},
        {3, 5, 6, 7},
        //{5, 7, 8, 9},
        //{7, 9, 10, 11}
    };

    // getting response layers
    struct response_layer *bottom;
    struct response_layer *middle;
    struct response_layer *top;

    // iterating through all octaves and each layer of octave in window of three (top, middle, bottom)
    for (int o = 0; o < fh_flat->octaves; ++o) {

        // TODO: (Sebastian) allow for fh->layers != 4 as well (note that fh->layers>=3 has to hold)
        for (int i = 0; i <= 1; ++i) {

            // assigning respective bottom, middle and top response layer
            bottom = &fh_flat->response_map[filter_map[o][i]];
            middle = &fh_flat->response_map[filter_map[o][i+1]];
            top = &fh_flat->response_map[filter_map[o][i+2]];

            // iterating over middle response layer at density of the most sparse layer (always top),
            // to find maxima accreoss scale and space
            for (int r = 0; r < top->height; ++r) {
                for (int c = 0; c < top->width; ++c) {

                    // checking if current pixel position is local maximum in 3x3x3 maximum and above threshold
                    if (is_extremum(r, c, top, middle, bottom, fh_flat->thresh)) {

                        // sub-pixel interpolating local maxium and adding to resulting interest point vector
                        interpolate_extremum(r, c, top, middle, bottom, interest_points);

                    }

                }
            }

        }
    }

}
