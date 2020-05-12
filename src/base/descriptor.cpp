#include "descriptor.h"
#include "integral_image.h"
#include "interest_point.h"
#include "helper.h"

#define _USE_MATH_DEFINES

#include <math.h>
#include <stdlib.h>
#include <vector>

void get_msurf_descriptor(struct integral_image* iimage, struct interest_point* ipoint) {

    float scale = ipoint->scale;
    int ipoint_x = (int) round(ipoint->x);
    int ipoint_y = (int) round(ipoint->y);

    // build descriptor
    float* descriptor = ipoint->descriptor;
    int desc_idx = 0;
    float sum_of_squares = 0.0f;
    
    // subregion centers for the 4x4 gaussian weighting
    float cx = -0.5f;
    float cy = 0.0f;

    //float co = 1;
    //float si = 0;

    int i = -8;

    // calculate descriptor for this interest point
    while(i < 12) {

        int j = -8;
        i = i - 4;

        cx += 1.0f;
        cy = -0.5f;

        while(j < 12) {
            
            float dx = 0.0f;
            float dy = 0.0f; 
            float mdx = 0.0f; 
            float mdy = 0.0f;

            cy += 1.0f;

            j = j - 4;

            //int xs = (int) round(ipoint_x + ( -jx*scale*si + ix*scale*co));
            //int ys = (int) round(ipoint_y + ( jx*scale*co + ix*scale*si));
            // TODO: (Sebastian) I think this should be i + 4 and j + 4 (OpenSURF also wrong)
            //int xs = (int) round(ipoint_x + (i + 5) * scale);
            //int ys = (int) round(ipoint_y + (j + 5) * scale);
            int xs = (int) round(ipoint_x + (i + 4) * scale);
            int ys = (int) round(ipoint_y + (j + 4) * scale);

            for (int k = i; k < i + 9; ++k) {
                for (int l = j; l < j + 9; ++l) {

                    //Get coords of sample point on the rotated axis
                    //int sample_x = (int) round(ipoint_x + (-l*scale*si + k*scale*co));
                    //int sample_y = (int) round(ipoint_y + ( l*scale*co + k*scale*si));
                    int sample_x = (int) round(ipoint_x + k * scale);
                    int sample_y = (int) round(ipoint_y + l * scale);

                    float gauss_s1 = gaussian((float) xs-sample_x, (float) ys-sample_y, 2.5f * scale);
                                        
                    float rx = haarX(iimage, sample_y, sample_x, (int) 2.0 * round(scale));
                    float ry = haarY(iimage, sample_y, sample_x, (int) 2.0 * round(scale));
                    
                    //Get the gaussian weighted x and y responses on rotated axis
                    //float rrx = gauss_s1*(-rx*si + ry*co);
                    //float rry = gauss_s1*(rx*co + ry*si);
                    // TODO: (Sebastian) This seems weird...
                    float rrx = gauss_s1 * (ry);
                    float rry = gauss_s1 * (rx);

                    dx += rrx;
                    dy += rry;
                    mdx += (float) fabs(rrx);
                    mdy += (float) fabs(rry);

                }      
            }

            // TODO: (Sebastian) Precompute this...
            float gauss_s2 = gaussian(cx-2.0f, cy-2.0f, 1.5f);

            // add the values to the descriptor vector
            descriptor[desc_idx] = dx * gauss_s2;
            descriptor[desc_idx+1] = dy * gauss_s2;
            descriptor[desc_idx+2] = mdx * gauss_s2;
            descriptor[desc_idx+3] = mdy * gauss_s2;

            // precompute for normaliztion
            sum_of_squares += (dx*dx + dy*dy + mdx*mdx + mdy*mdy) * gauss_s2*gauss_s2;

            desc_idx += 4;

            j += 9;

        }

        i += 9;

    }

    // rescale to unit vector
    float norm_factor = 1.0f / sqrt(sum_of_squares);

    for (int i = 0; i < 64; ++i) {
        descriptor[i] *= norm_factor;
    }

}

void get_msurf_descriptors(struct integral_image* iimage, std::vector<struct interest_point> *interest_points) {
    for (int i = 0; i < interest_points->size(); ++i) {
        get_msurf_descriptor(iimage, &interest_points->at(i));
	}
}
