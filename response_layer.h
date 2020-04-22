#pragma once

#include <assert.h>

struct response_layer {
    
    // Size of filter
    int filter_size;

    // Width and height of response layer 
    // (smaller by factor of 1, 1/2, 1/4, ...  to original image depending on ocatave)
    int width;
    int height;

    // TODO: (Sebastian) decide if we want to have this initial samping steps in here or not
    int step = 2;

    // Array with all hessian responses
    float *response;

    // Array with all flags marking positive laplacian (Dxx + Dyy >= 0)
    bool *laplacian;

};
