#include <stdlib.h>

#define NUM_OCTAVES 4
#define NUM_LAYER 10

struct fasthessian {

    // Integral image
    struct integral_image* iimage;

    // Response stack of determinant of hessian values
    struct response_layer* response_map[NUM_LAYER];

    // Number of Octaves
    int octaves;

    // Number of layers per octave
    int layers;

    // Initial sampling step for Interest Point detection
    int step;

    // Threshold value for responses
    float thresh;

};

// Create Fast-Hessian struct
void createFastHessian();

// Create octaves with response layers
void buildResponseMap();

// Compute responses for layer
void buildResponseLayer(struct response_layer* layer);
