#pragma once

#include <stdlib.h>
#include <stdbool.h> 

struct interest_point {
    
    // Detected (sub) pixel position
    float x;
    float y;

    // Detected scale
    float scale;

    // Detected orientation [radiant] (with 0.0 == upright)
    float orientation = 0.0f;

    // Flag indicating if U-SURF is used
    bool upright = true;

    // Flag indicating if laplacian Dxx + Dyy >= 0
    bool laplacian;

    // SURF descriptor
    float descriptor[64];

};
