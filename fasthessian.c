#include "fasthessian.h"

#include "interest_point.h"
#include "helper.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>

struct fasthessian* create_fast_hessian(struct integral_image *iimage) {

    struct fasthessian* fh = (struct fasthessian*) malloc(sizeof(struct fasthessian));

    // Interest Image
    fh->iimage = iimage;

    // Set Variables
    fh->octaves = NUM_OCTAVES;
    fh->layers = NUM_LAYER;
    fh->step = INITIAL_STEP;
    fh->thresh = THRESHOLD;

    return fh;
}

void create_response_map(struct fasthessian* fh, struct integral_image *iimage) {
    int img_width = (fh->iimage)->width;
    int img_height = (fh->iimage)->height;
    int init_step = fh->step;

    int w = img_width / init_step;
    int h = (i_height / init_step);

    // Octave 1 - 9, 15, 21, 27
    struct response_layer* layer0 = (struct response_layer *) malloc(sizeof(struct response_layer));
    layer0->filter_size = 9;
    layer0->width = w;
    layer0->height = h;
    layer0->step = init_step;

    struct response_layer* layer1 = (struct response_layer *) malloc(sizeof(struct response_layer));
    layer1->filter_size = 15;
    layer1->width = w;
    layer1->height = h;
    layer1->step = init_step;

    struct response_layer* layer2 = (struct response_layer *) malloc(sizeof(struct response_layer));
    layer2->filter_size = 21;
    layer2->width = w;
    layer2->height = h;
    layer2->step = init_step;

    struct response_layer* layer3 = (struct response_layer *) malloc(sizeof(struct response_layer));
    layer3->filter_size = 27;
    layer3->width = w;
    layer3->height = h;
    layer3->step = init_step;

    fh->response_map[0] = layer0;
    fh->response_map[1] = layer1;
    fh->response_map[2] = layer2;
    fh->response_map[3] = layer3;

    // Octave 2 - 15, 27, 39, 51
    struct response_layer* layer4 = (struct response_layer *) malloc(sizeof(struct response_layer));
    layer4->filter_size = 39;
    layer4->width = w/2;
    layer4->height = h/2;
    layer4->step = init_step*2;

    struct response_layer* layer5 = (struct response_layer *) malloc(sizeof(struct response_layer));
    layer5->filter_size = 51;
    layer5->width = w/2;
    layer5->height = h/2;
    layer5->step = init_step*2;

    fh->response_map[4] = layer4;
    fh->response_map[5] = layer5;

    // Octave 3 - 27, 51, 75, 99
    struct response_layer* layer6 = (struct response_layer *) malloc(sizeof(struct response_layer));
    layer6->filter_size = 75;
    layer6->width = w/4;
    layer6->height = h/4;
    layer6->step = init_step*4;

    struct response_layer* layer7 = (struct response_layer *) malloc(sizeof(struct response_layer));
    layer7->filter_size = 99;
    layer7->width = w/4;
    layer7->height = h/4;
    layer7->step = init_step*4;

    fh->response_map[6] = layer6;
    fh->response_map[7] = layer7;
}

void compute_response_layer(struct response_layer* layer) {

}

bool is_extremum(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom) {

    assert(top != NULL && middle != NULL && bottom != NULL);

    // TODO: (Sebastian) Don't quite understand this with the steps in there...
    int layer_border = (top->filter_size + 1) / (2 * top->step)

    // checking if row or col are out of bounds
    if (row <= layer_border || top->height - layer_border <= row
        || col <= layer_border || top->width - layer_border <= col) {
        return false;
    }

    // getting candidate point in middle layer
    float candidate = middle->get_response(row, col, top);

    // checking if it passes the threshold
    if (candidate < threshold) {
        return false;
    }

    // iterating over 3x3x3 neighborhood and checking for local maxima
    for (int rr = -1; rr <= 1; ++rr) {
        for (int cc = -1; cc <= 1; ++cc) {
            // checking if any other response in the 3x3x3 neighborhood has a higher hessian response than the candidate
            if (candidate <= top->get_response(row+rr, col+cc)
                || ((rr != 0 || cc != 0) && candidate <= middle->get_response(row+rr, col+cc, top))
                || candidate <= bottom->get_response(row+rr, col+cc, top)) {
                return false;
            }
        }
    }

    return true;

}

void interpolate_extremum(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom) {

    // getting step distance between filters
    int filter_step = (middle->filter_size - bottom->filter_size);

    // checking if middle filter is mid way between top and bottom
    assert(filter_step > 0 && top->filter_size - middle->filter_size == middle->filter_size - bottom->filter_size);

    // getting sub-pixel offset of extremum to actual location
    float offsets[3];
    interpolate_step(row, col, top, middle, bottom, offsets);

    // getting sub-pixel offsets and assigning them to variables for clarity
    float dx = offsets[0];
    float dy = offsets[1];
    float ds = offsets[2];

    // checking if sub-pixel offsets are less than 0.5 and
    // would thus rounded sub-pixel location is the same as pixel location
    if (fabs(dx) < 0.5f && fabs(dy) < 0.5f && fabs(ds) < 0.5f) {

        // initializing interest point
        struct interest_point ip;
        ip.x = (float) ((col + dx) * top->step);
        ip.y = (float) ((row + dy) * top->step);
        ip.scale = (float) ((0.1333f) * (middle->filter_size + ds * filter_step));
        ip.laplacian = middle->get_laplacian(row, col, top);

        // TODO: (Sebastian) Only valid for U-SURF
        ip.orientation = 0.0;
        ip.upright = true;

        // TODO: (Sebastian) Add interest point to list
        // ...

    }

}

void interpolate_step(int row, int col,
                      struct response_layer *top, struct response_layer *middle, struct response_layer *bottom,
                      float offsets[3]) {

    assert(top != NULL && middle != NULL && bottom != NULL);
    assert(offsets != NULL);

    float hessian[9];
    {
        float v = middle->get_response(row, col, top);

        // computing second order partial derivatives in xy position, as well as scale direction
        float dxx = middle->get_response(row, col + 1, top) + middle->get_response(row, col - 1, top) - 2.0f * v;
        float dyy = middle->get_response(row + 1, col, top) + middle->get_response(row - 1, col, top) - 2.0f * v;
        float dss = top->get_response(row, col) + bottom->get_response(row, col, top) - 2.0f * v;

        float dxy = (middle->get_response(row + 1, col + 1, top) - middle->get_response(row + 1, col - 1, top) -
                     middle->get_response(row - 1, col + 1, top) + middle->get_response(row - 1, col - 1, top)) / 4.0f;
        float dxs = (top->get_response(row, col + 1) - top->get_response(row, col - 1) -
                     bottom->get_response(row, col + 1, top) + bottom->get_response(row, col - 1, top)) / 4.0f;
        float dys = (top->get_response(row + 1, col) - top->get_response(row - 1, col) -
                     bottom->get_response(row + 1, col, top) + bottom->get_response(row - 1, col, top)) / 4.0f;

        // constructing hessian 3x3 matrix:
        // dxx dxy dxs
        // dxy dyy dys
        // dxs dys dss
        hessian[0] = dxx;
        hessian[1] = dxy;
        hessian[2] = dxs;
        hessian[3] = dxy;
        hessian[4] = dyy;
        hessian[5] = dys;
        hessian[6] = dxs;
        hessian[7] = dys;
        hessian[8] = dss;
    }

    float neg_gradient[3];
    {
        float dx = (middle->get_response(row, col + 1, top) - middle->get_response(row, col - 1, top)) / 2.0f;
        float dy = (middle->get_response(row + 1, col, top) - middle->get_response(row - 1, col, top)) / 2.0f;
        float ds = (top->get_response(row, col) - bottom->get_response(row, col, top)) / 2.0f;

        // constructing negative gradient 3x1 vector
        neg_gradient[0] = -dx;
        neg_gradient[1] = -dy;
        neg_gradient[2] = -ds;
    }

    // solving linear system hessian * offsets = neg_gradient to get sub-pixel offsets
    solve_linear_3x3_system(hessian, neg_gradient, offsets);

}
