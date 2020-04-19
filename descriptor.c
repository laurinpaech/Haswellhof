#include <math.h>
#include <stdlib.h>

void get_gaussian(float sigma, int size, float* dest) {
    // compute matrix of shape size*size 2d gaussian with variance sigma and mean (size/2, size/2) 
    // TODO
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            dest[i*size + j] = 1;
        }
    }
}

void get_descriptor(int keyp_x, int keyp_y, int scale, float* iimage, int height, int width, int PATCH_SIZE, *float descriptor) {

    // assert keyp_x, keyp_y is far enough out

    int top_right = keyp_y/2*width + keyp_x/2; // TODO: check if  indicescorrect

    // compute/load gaussian weighting matrix GW for PATCH_SIZE here 
    // due to resue over all patches we probably want to precompute it and pass it as a param 
    GW = (float*) malloc(PATCH_SIZE*PATCH_SIZE * sizeof(float));
    get_gaussian(3.3, PATCH_SIZE, GW);

    // store wavelet responses
    float dx[PATCH_SIZE][PATCH_SIZE];
    float dy[PATCH_SIZE][PATCH_SIZE];

    // compute haar wavelet responses
    for (int i=0, i<PATCH_SIZE, i++) { // x coordinate
        for (int j=0, j<PATCH_SIZE, j++) { // y coordinate
            float gw = GW[i*PATCH_SIZE + j];

            // compute needed corners of the 2sx2s patch
            // c1 c2 c3 
            // c4 ij c5
            // c6 c7 c8

            // TODO: check if  indicescorrect
            // int m = (j*width+i)*scale;
            int c1 = top_right + ((j)*scale-1)*width+((i)*scale-1);
            int c2 = top_right + ((j)*scale-1)*width+((i+1)*scale-1);
            int c3 = top_right + ((j)*scale-1)*width+((i+2)*scale-1);
            int c4 = top_right + ((j+1)*scale-1)*width+((i)*scale-1);
            int c6 = top_right + ((j+2)*scale-1)*width+((i)*scale-1);
            int c5 = top_right + ((j+1)*scale-1)*width+((i+2)*scale-1);
            int c7 = top_right + ((j+2)*scale-1)*width+((i+1)*scale-1);
            int c8 = top_right + ((j+2)*scale-1)*width+((i+2)*scale-1);

            dx[i][j] = gw*(iimages[c8] + iimages[c4] - iimages[c5] - iimages[c6] - (iimages[c5] + iimages[c1] - iimages[c3] - iimages[c4]));
            dy[i][j] = gw*(iimages[c8] + iimages[c2] - iimages[c3] - iimages[c7] - (iimages[c7] + iimages[c1] - iimages[c2] - iimages[c6]));

            // dx[i][j] = gw*(patch[(i+1)*PATCH_SIZE + j] - patch[i*PATCH_SIZE + j] + patch[(i+1)*PATCH_SIZE + j+1] - patch[i*PATCH_SIZE + j+1]);
            // dy[i][j] = gw*(patch[i*PATCH_SIZE + j+1] - patch[i*PATCH_SIZE + j] + patch[(i+1)*PATCH_SIZE + j+1] - patch[(i+1)*PATCH_SIZE + j]);
        }
    }

    // build descriptor
    int desc_idx = 0;
    float sum_of_squares = 0;
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) { // iterate over sub_patches 
            descriptor[desc_idx] = 0;
            descriptor[desc_idx+1] = 0;
            descriptor[desc_idx+2] = 0;
            descriptor[desc_idx+3] = 0;
            for (int k=i*5; k<i*5+5; k++) {
                for (int l=j*5; k<j*5+5; j++) { // iterate over pixels of 5x5 sub_patches
                    float x = dx[k][l]; 
                    float y = dy[k][l]; 

                    descriptor[desc_idx] += x; // sum(x)
                    descriptor[desc_idx+1] += y; // sum(y)
                    descriptor[desc_idx+2] += (float)fabs(x); // sum(abs(x))
                    descriptor[desc_idx+3] += (float)fabs(y); // sum(abs(y))
                }   
            }
            
            // precompute for normaliztion
            for (int m=0; m<4; m++) 
                sum_of_squares += descriptor[m]*descriptor[m];

            desc_idx += 4;
        }
    }

    // rescale to unit vector
    float norm_factor = 1./sqrt(sum_of_squares);
    for (int i=0; i<64; i++) 
        descriptor[i] *= norm_factor;

    return *descriptor;
}