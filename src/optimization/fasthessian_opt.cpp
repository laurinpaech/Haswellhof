#include "fasthessian_opt.h"


void get_interest_points_layers(struct fasthessian *fh, std::vector<struct interest_point> *interest_points) {

    assert(fh != NULL);
    assert(interest_points != NULL);

    // getting threshold
    float thresh = fh->thresh;

    // getting response layers
    struct response_layer *l0 = fh->response_map[0];
    struct response_layer *l1 = fh->response_map[1];
    struct response_layer *l2 = fh->response_map[2];
    struct response_layer *l3 = fh->response_map[3];
    struct response_layer *l4 = fh->response_map[4];
    struct response_layer *l5 = fh->response_map[5];
    struct response_layer *l6 = fh->response_map[6];
    struct response_layer *l7 = fh->response_map[7];

    // getting height and width for layers
    int height2 = l2->height;
    int width2 = l2->width;
    int height3 = l3->height;
    int width3 = l3->width;
    int height4 = l4->height;
    int width4 = l4->width;
    int height5 = l5->height;
    int width5 = l5->width;
    int height6 = l6->height;
    int width6 = l6->width;
    int height7 = l7->height;
    int width7 = l7->width;

    {
        // iterating over middle response layer at density of the most sparse layer (always top),
        // to find maxima accreoss scale and space
        for (int r = 0; r < height2; ++r) {
            for (int c = 0; c < width2; ++c) {

                // checking if current pixel position is local maximum in 3x3x3 maximum and above threshold
                if (is_extremum(r, c, l2, l1, l0, thresh)) {

                    // sub-pixel interpolating local maxium and adding to resulting interest point vector
                    interpolate_extremum(r, c, l2, l1, l0, interest_points);

                }

            }
        }
    }

    {
        // iterating over middle response layer at density of the most sparse layer (always top),
        // to find maxima accreoss scale and space
        for (int r = 0; r < height3; ++r) {
            for (int c = 0; c < width3; ++c) {

                // checking if current pixel position is local maximum in 3x3x3 maximum and above threshold
                if (is_extremum(r, c, l3, l2, l1, thresh)) {

                    // sub-pixel interpolating local maxium and adding to resulting interest point vector
                    interpolate_extremum(r, c, l3, l2, l1, interest_points);

                }

            }
        }
    }

    {
        // iterating over middle response layer at density of the most sparse layer (always top),
        // to find maxima accreoss scale and space
        for (int r = 0; r < height4; ++r) {
            for (int c = 0; c < width4; ++c) {

                // checking if current pixel position is local maximum in 3x3x3 maximum and above threshold
                if (is_extremum(r, c, l4, l3, l1, thresh)) {

                    // sub-pixel interpolating local maxium and adding to resulting interest point vector
                    interpolate_extremum(r, c, l4, l3, l1, interest_points);

                }

            }
        }
    }

    {
        // iterating over middle response layer at density of the most sparse layer (always top),
        // to find maxima accreoss scale and space
        for (int r = 0; r < height5; ++r) {
            for (int c = 0; c < width5; ++c) {

                // checking if current pixel position is local maximum in 3x3x3 maximum and above threshold
                if (is_extremum(r, c, l5, l4, l3, thresh)) {

                    // sub-pixel interpolating local maxium and adding to resulting interest point vector
                    interpolate_extremum(r, c, l5, l4, l3, interest_points);

                }

            }
        }
    }

    {
        // iterating over middle response layer at density of the most sparse layer (always top),
        // to find maxima accreoss scale and space
        for (int r = 0; r < height6; ++r) {
            for (int c = 0; c < width6; ++c) {

                // checking if current pixel position is local maximum in 3x3x3 maximum and above threshold
                if (is_extremum(r, c, l6, l5, l3, thresh)) {

                    // sub-pixel interpolating local maxium and adding to resulting interest point vector
                    interpolate_extremum(r, c, l6, l5, l3, interest_points);

                }

            }
        }
    }

    {
        // iterating over middle response layer at density of the most sparse layer (always top),
        // to find maxima accreoss scale and space
        for (int r = 0; r < height7; ++r) {
            for (int c = 0; c < width7; ++c) {

                // checking if current pixel position is local maximum in 3x3x3 maximum and above threshold
                if (is_extremum(r, c, l7, l6, l5, thresh)) {

                    // sub-pixel interpolating local maxium and adding to resulting interest point vector
                    interpolate_extremum(r, c, l7, l6, l5, interest_points);

                }

            }
        }
    }

}


void interpolate_step_gauss(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets[3]) {
    
}
