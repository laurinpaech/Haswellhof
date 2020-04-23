#pragma once

#include <assert.h>

struct response_layer {

    // Size of filter
    int filter_size;

    // Width and height of response layer
    // (smaller by factor of 1, 1/2, 1/4, ...  to original image depending on ocatave)
    int width;
    int height;

    // Sampling step for this layer
    int step;

    // Array with all hessian responses
    float *response;

    // Array with all flags marking positive laplacian (Dxx + Dyy >= 0)
    bool *laplacian;

};

inline bool get_laplacian(struct response_layer *rl, int row, int col) {
    assert(rl != NULL && rl->laplacian != NULL);
    return rl->laplacian[row * rl->width + col];
}

inline bool get_laplacian(struct response_layer *rl, int row, int col, struct response_layer *src) {
    assert(rl != NULL && rl->laplacian != NULL && src != NULL);
    int scale = rl->width / src->width;
    return rl->laplacian[(scale * row) * rl->width + (scale * col)];
}

inline float get_response(struct response_layer *rl, int row, int col) {
    assert(rl != NULL && rl->response != NULL);
    return rl->response[row * rl->width + col];
}

inline float get_response(struct response_layer *rl, int row, int col, struct response_layer *src) {
    assert(rl != NULL && rl->response != NULL && src != NULL);
    int scale = rl->width / src->width;
    return rl->response[(scale * row) * rl->width + (scale * col)];
}
