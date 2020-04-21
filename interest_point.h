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

    // Flag for U-SURF
    bool upright = true;

    // SURF descriptor
    float descriptor[64];

};
