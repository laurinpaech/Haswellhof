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

    //std::cout << "sizeof(bool): " << sizeof(bool) << std::endl;

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

void compute_response_layers_flat(struct fasthessian_flat* fh_flat) {
    
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

void compute_response_layers_at_once_flat(struct fasthessian_flat *fh_flat) {
    /* computes all 8 response layers at once, gives same results as base implementation
    valgrind reports no improvement for l1 misses, i.e. locality is not improved as expected
    guess due to different filter sizes always accesing different cachelines
    need to be smarter about computation order to gain benefit (i.e. same mem addresses are actually accessed at once)
    */
    float Dxx0, Dyy0, Dxy0;
    float Dxx1, Dyy1, Dxy1;
    float Dxx2, Dyy2, Dxy2;
    float Dxx3, Dyy3, Dxy3;
    float Dxx4, Dyy4, Dxy4;
    float Dxx5, Dyy5, Dxy5;
    float Dxx6, Dyy6, Dxy6;
    float Dxx7, Dyy7, Dxy7;

    int x, y;
    int step = fh_flat->step;
    struct integral_image* iimage = fh_flat->iimage;

    // step == 1: 0,1,2,3
    // step == 2: 4,5
    // step == 4: 6,7

    struct response_layer *layer0 = &fh_flat->response_map[0];  // step0 == step
    struct response_layer *layer1 = &fh_flat->response_map[1];  // step1 == step
    struct response_layer *layer2 = &fh_flat->response_map[2];  // step2 == step
    struct response_layer *layer3 = &fh_flat->response_map[3];  // step3 == step

    struct response_layer *layer4 = &fh_flat->response_map[4];  // step4 == 2*step
    struct response_layer *layer5 = &fh_flat->response_map[5];  // step5 == 2*step

    struct response_layer *layer6 = &fh_flat->response_map[6];  // step6 == 4*step
    struct response_layer *layer7 = &fh_flat->response_map[7];  // step7 == 4*step

    float *response0 = layer0->response;
    float *response1 = layer1->response;
    float *response2 = layer2->response;
    float *response3 = layer3->response;
    float *response4 = layer4->response;
    float *response5 = layer5->response;
    float *response6 = layer6->response;
    float *response7 = layer7->response;

    bool *laplacian0 = layer0->laplacian;
    bool *laplacian1 = layer1->laplacian;
    bool *laplacian2 = layer2->laplacian;
    bool *laplacian3 = layer3->laplacian;
    bool *laplacian4 = layer4->laplacian;
    bool *laplacian5 = layer5->laplacian;
    bool *laplacian6 = layer6->laplacian;
    bool *laplacian7 = layer7->laplacian;

    int filter_size0 = layer0->filter_size;
    int filter_size1 = layer1->filter_size;
    int filter_size2 = layer2->filter_size;
    int filter_size3 = layer3->filter_size;
    int filter_size4 = layer4->filter_size;
    int filter_size5 = layer5->filter_size;
    int filter_size6 = layer6->filter_size;
    int filter_size7 = layer7->filter_size;

    int height = layer0->height;
    int width = layer0->width;

    int lobe0 = filter_size0 / 3;
    int lobe1 = filter_size1 / 3;
    int lobe2 = filter_size2 / 3;
    int lobe3 = filter_size3 / 3;
    int lobe4 = filter_size4 / 3;
    int lobe5 = filter_size5 / 3;
    int lobe6 = filter_size6 / 3;
    int lobe7 = filter_size7 / 3;

    int border0 = (filter_size0 - 1) / 2;
    int border1 = (filter_size1 - 1) / 2;
    int border2 = (filter_size2 - 1) / 2;
    int border3 = (filter_size3 - 1) / 2;
    int border4 = (filter_size4 - 1) / 2;
    int border5 = (filter_size5 - 1) / 2;
    int border6 = (filter_size6 - 1) / 2;
    int border7 = (filter_size7 - 1) / 2;

    float inv_area0 = 1.f / (filter_size0 * filter_size0);
    float inv_area1 = 1.f / (filter_size1 * filter_size1);
    float inv_area2 = 1.f / (filter_size2 * filter_size2);
    float inv_area3 = 1.f / (filter_size3 * filter_size3);
    float inv_area4 = 1.f / (filter_size4 * filter_size4);
    float inv_area5 = 1.f / (filter_size5 * filter_size5);
    float inv_area6 = 1.f / (filter_size6 * filter_size6);
    float inv_area7 = 1.f / (filter_size7 * filter_size7);

    for (int i = 0, ind = 0, ind2 = 0, ind4 = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j, ind++) {
            // Image coordinates
            x = i * step;
            y = j * step;

            // layer0
            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx0 = box_integral(iimage, x - lobe0 + 1, y - border0, 2 * lobe0 - 1, filter_size0) -
                   3 * box_integral(iimage, x - lobe0 + 1, y - lobe0 / 2, 2 * lobe0 - 1, lobe0);
            Dyy0 = box_integral(iimage, x - border0, y - lobe0 + 1, filter_size0, 2 * lobe0 - 1) -
                   3 * box_integral(iimage, x - lobe0 / 2, y - lobe0 + 1, lobe0, 2 * lobe0 - 1);
            Dxy0 = box_integral(iimage, x - lobe0, y + 1, lobe0, lobe0) +
                   box_integral(iimage, x + 1, y - lobe0, lobe0, lobe0) -
                   box_integral(iimage, x - lobe0, y - lobe0, lobe0, lobe0) -
                   box_integral(iimage, x + 1, y + 1, lobe0, lobe0);

            // Normalize Responses with inverse area
            Dxx0 *= inv_area0;
            Dyy0 *= inv_area0;
            Dxy0 *= inv_area0;

            // Calculate Determinant
            response0[ind] = Dxx0 * Dyy0 - 0.81f * Dxy0 * Dxy0;

            // Calculate Laplacian
            laplacian0[ind] = (Dxx0 + Dyy0 >= 0 ? true : false);


            // layer1
            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx1 = box_integral(iimage, x - lobe1 + 1, y - border1, 2 * lobe1 - 1, filter_size1) -
                   3 * box_integral(iimage, x - lobe1 + 1, y - lobe1 / 2, 2 * lobe1 - 1, lobe1);
            Dyy1 = box_integral(iimage, x - border1, y - lobe1 + 1, filter_size1, 2 * lobe1 - 1) -
                   3 * box_integral(iimage, x - lobe1 / 2, y - lobe1 + 1, lobe1, 2 * lobe1 - 1);
            Dxy1 = box_integral(iimage, x - lobe1, y + 1, lobe1, lobe1) +
                   box_integral(iimage, x + 1, y - lobe1, lobe1, lobe1) -
                   box_integral(iimage, x - lobe1, y - lobe1, lobe1, lobe1) -
                   box_integral(iimage, x + 1, y + 1, lobe1, lobe1);

            // Normalize Responses with inverse area
            Dxx1 *= inv_area1;
            Dyy1 *= inv_area1;
            Dxy1 *= inv_area1;

            // Calculate Determinant
            response1[ind] = Dxx1 * Dyy1 - 0.81f * Dxy1 * Dxy1;

            // Calculate Laplacian
            laplacian1[ind] = (Dxx1 + Dyy1 >= 0 ? true : false);


            // layer2
            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx2 = box_integral(iimage, x - lobe2 + 1, y - border2, 2 * lobe2 - 1, filter_size2) -
                   3 * box_integral(iimage, x - lobe2 + 1, y - lobe2 / 2, 2 * lobe2 - 1, lobe2);
            Dyy2 = box_integral(iimage, x - border2, y - lobe2 + 1, filter_size2, 2 * lobe2 - 1) -
                   3 * box_integral(iimage, x - lobe2 / 2, y - lobe2 + 1, lobe2, 2 * lobe2 - 1);
            Dxy2 = box_integral(iimage, x - lobe2, y + 1, lobe2, lobe2) +
                   box_integral(iimage, x + 1, y - lobe2, lobe2, lobe2) -
                   box_integral(iimage, x - lobe2, y - lobe2, lobe2, lobe2) -
                   box_integral(iimage, x + 1, y + 1, lobe2, lobe2);

            // Normalize Responses with inverse area
            Dxx2 *= inv_area2;
            Dyy2 *= inv_area2;
            Dxy2 *= inv_area2;

            // Calculate Determinant
            response2[ind] = Dxx2 * Dyy2 - 0.81f * Dxy2 * Dxy2;

            // Calculate Laplacian
            laplacian2[ind] = (Dxx2 + Dyy2 >= 0 ? true : false);


            // layer3
            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx3 = box_integral(iimage, x - lobe3 + 1, y - border3, 2*lobe3 - 1, filter_size3)
                    - 3 * box_integral(iimage, x - lobe3 + 1, y - lobe3 / 2, 2*lobe3 - 1, lobe3);
            Dyy3 = box_integral(iimage, x - border3, y - lobe3 + 1, filter_size3, 2*lobe3 - 1)
                    - 3 * box_integral(iimage, x - lobe3 / 2, y - lobe3 + 1, lobe3, 2*lobe3 - 1);
            Dxy3 = box_integral(iimage, x - lobe3, y + 1, lobe3, lobe3)
                    + box_integral(iimage, x + 1, y - lobe3, lobe3, lobe3)
                    - box_integral(iimage, x - lobe3, y - lobe3, lobe3, lobe3)
                    - box_integral(iimage, x + 1, y + 1, lobe3, lobe3);

            // Normalize Responses with inverse area
            Dxx3 *= inv_area3;
            Dyy3 *= inv_area3;
            Dxy3 *= inv_area3;

            // Calculate Determinant
            response3[ind] = Dxx3 * Dyy3 - 0.81f * Dxy3 * Dxy3;

            // Calculate Laplacian
            laplacian3[ind] = (Dxx3 + Dyy3 >= 0 ? true : false);


            // 2*step
            if (i % 2 == 0 && j % 2 == 0) {
                // layer4
                // Calculate Dxx, Dyy, Dxy with Box Filter
                Dxx4 = box_integral(iimage, x - lobe4 + 1, y - border4, 2*lobe4 - 1, filter_size4)
                        - 3 * box_integral(iimage, x - lobe4 + 1, y - lobe4 / 2, 2*lobe4 - 1, lobe4);
                Dyy4 = box_integral(iimage, x - border4, y - lobe4 + 1, filter_size4, 2*lobe4 - 1)
                        - 3 * box_integral(iimage, x - lobe4 / 2, y - lobe4 + 1, lobe4, 2*lobe4 - 1);
                Dxy4 = box_integral(iimage, x - lobe4, y + 1, lobe4, lobe4)
                        + box_integral(iimage, x + 1, y - lobe4, lobe4, lobe4)
                        - box_integral(iimage, x - lobe4, y - lobe4, lobe4, lobe4)
                        - box_integral(iimage, x + 1, y + 1, lobe4, lobe4);

                // Normalize Responses with inverse area
                Dxx4 *= inv_area4;
                Dyy4 *= inv_area4;
                Dxy4 *= inv_area4;

                // Calculate Determinant
                response4[ind2] = Dxx4 * Dyy4 - 0.81f * Dxy4 * Dxy4;

                // Calculate Laplacian
                laplacian4[ind2] = (Dxx4 + Dyy4 >= 0 ? true : false);


                // layer5
                // Calculate Dxx, Dyy, Dxy with Box Filter
                Dxx5 = box_integral(iimage, x - lobe5 + 1, y - border5, 2*lobe5 - 1, filter_size5)
                        - 3 * box_integral(iimage, x - lobe5 + 1, y - lobe5 / 2, 2*lobe5 - 1, lobe5);
                Dyy5 = box_integral(iimage, x - border5, y - lobe5 + 1, filter_size5, 2*lobe5 - 1)
                        - 3 * box_integral(iimage, x - lobe5 / 2, y - lobe5 + 1, lobe5, 2*lobe5 - 1);
                Dxy5 = box_integral(iimage, x - lobe5, y + 1, lobe5, lobe5)
                        + box_integral(iimage, x + 1, y - lobe5, lobe5, lobe5)
                        - box_integral(iimage, x - lobe5, y - lobe5, lobe5, lobe5)
                        - box_integral(iimage, x + 1, y + 1, lobe5, lobe5);

                // Normalize Responses with inverse area
                Dxx5 *= inv_area5;
                Dyy5 *= inv_area5;
                Dxy5 *= inv_area5;

                // Calculate Determinant
                response5[ind2] = Dxx5 * Dyy5 - 0.81f * Dxy5 * Dxy5;

                // Calculate Laplacian
                laplacian5[ind2] = (Dxx5 + Dyy5 >= 0 ? true : false);


                ind2++;

                // 4*step
                if (i % 4 == 0 && j % 4 == 0) {
                    // layer6
                    // Calculate Dxx, Dyy, Dxy with Box Filter
                    Dxx6 = box_integral(iimage, x - lobe6 + 1, y - border6, 2*lobe6 - 1, filter_size6)
                            - 3 * box_integral(iimage, x - lobe6 + 1, y - lobe6 / 2, 2*lobe6 - 1, lobe6);
                    Dyy6 = box_integral(iimage, x - border6, y - lobe6 + 1, filter_size6, 2*lobe6 - 1)
                            - 3 * box_integral(iimage, x - lobe6 / 2, y - lobe6 + 1, lobe6, 2*lobe6 - 1);
                    Dxy6 = box_integral(iimage, x - lobe6, y + 1, lobe6, lobe6)
                            + box_integral(iimage, x + 1, y - lobe6, lobe6, lobe6)
                            - box_integral(iimage, x - lobe6, y - lobe6, lobe6, lobe6)
                            - box_integral(iimage, x + 1, y + 1, lobe6, lobe6);

                    // Normalize Responses with inverse area
                    Dxx6 *= inv_area6;
                    Dyy6 *= inv_area6;
                    Dxy6 *= inv_area6;

                    // Calculate Determinant
                    response6[ind4] = Dxx6 * Dyy6 - 0.81f * Dxy6 * Dxy6;

                    // Calculate Laplacian
                    laplacian6[ind4] = (Dxx6 + Dyy6 >= 0 ? true : false);


                    // layer7
                    // Calculate Dxx, Dyy, Dxy with Box Filter
                    Dxx7 = box_integral(iimage, x - lobe7 + 1, y - border7, 2*lobe7 - 1, filter_size7)
                            - 3 * box_integral(iimage, x - lobe7 + 1, y - lobe7 / 2, 2*lobe7 - 1, lobe7);
                    Dyy7 = box_integral(iimage, x - border7, y - lobe7 + 1, filter_size7, 2*lobe7 - 1)
                            - 3 * box_integral(iimage, x - lobe7 / 2, y - lobe7 + 1, lobe7, 2*lobe7 - 1);
                    Dxy7 = box_integral(iimage, x - lobe7, y + 1, lobe7, lobe7)
                            + box_integral(iimage, x + 1, y - lobe7, lobe7, lobe7)
                            - box_integral(iimage, x - lobe7, y - lobe7, lobe7, lobe7)
                            - box_integral(iimage, x + 1, y + 1, lobe7, lobe7);

                    // Normalize Responses with inverse area
                    Dxx7 *= inv_area7;
                    Dyy7 *= inv_area7;
                    Dxy7 *= inv_area7;

                    // Calculate Determinant
                    response7[ind4] = Dxx7 * Dyy7 - 0.81f * Dxy7 * Dxy7;

                    // Calculate Laplacian
                    laplacian7[ind4] = (Dxx7 + Dyy7 >= 0 ? true : false);


                    ind4++;
                }
            }

            // template for search and replace J to layer number
            // // layerJ
            // // Calculate Dxx, Dyy, Dxy with Box Filter
            // DxxJ = box_integral(iimage, x - lobeJ + 1, y - borderJ, 2*lobeJ - 1, filter_sizeJ)
            //         - 3 * box_integral(iimage, x - lobeJ + 1, y - lobeJ / 2, 2*lobeJ - 1, lobeJ);
            // DyyJ = box_integral(iimage, x - borderJ, y - lobeJ + 1, filter_sizeJ, 2*lobeJ - 1)
            //         - 3 * box_integral(iimage, x - lobeJ / 2, y - lobeJ + 1, lobeJ, 2*lobeJ - 1);
            // DxyJ = box_integral(iimage, x - lobeJ, y + 1, lobeJ, lobeJ)
            //         + box_integral(iimage, x + 1, y - lobeJ, lobeJ, lobeJ)
            //         - box_integral(iimage, x - lobeJ, y - lobeJ, lobeJ, lobeJ)
            //         - box_integral(iimage, x + 1, y + 1, lobeJ, lobeJ);

            // // Normalize Responses with inverse area
            // DxxJ *= inv_areaJ;
            // DyyJ *= inv_areaJ;
            // DxyJ *= inv_areaJ;

            // // Calculate Determinant
            // responseJ[ind] = DxxJ * DyyJ - 0.81f * DxyJ * DxyJ;

            // // Calculate Laplacian
            // laplacianJ[ind] = (DxxJ + DyyJ >= 0 ? true : false);
        }
    }

}
