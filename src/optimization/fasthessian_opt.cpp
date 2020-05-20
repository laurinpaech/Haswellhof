#include "fasthessian_opt.h"

#include <stdio.h>
#include <stdlib.h>

// Only for images >= 128x128
void compute_response_layers_Dyy(struct fasthessian *fh) {
    for (int i = 0; i < fh->total_layers; ++i) {
        compute_response_layer_Dyy(fh->response_map[i], fh->iimage);
    }
}

// Only for images >= 128x128
void compute_response_layers_Dyy_laplacian(struct fasthessian *fh) {
    for (int i = 0; i < fh->total_layers; ++i) {
        compute_response_layer_Dyy_laplacian(fh->response_map[i], fh->iimage);
    }
}

// Only for images >= 128x128
void compute_response_layers_Dyy_laplacian_localityloops(struct fasthessian *fh) {
    for (int i = 0; i < fh->total_layers; ++i) {
        compute_response_layer_Dyy_laplacian_localityloops(fh->response_map[i], fh->iimage);
    }
}

void compute_response_layers_Dyy_leftcorner(struct fasthessian *fh) {
    for (int i = 0; i < fh->total_layers; ++i) {
        compute_response_layer_Dyy_leftcorner(fh->response_map[i], fh->iimage);
    }
}

void compute_response_layers_Dyy_top(struct fasthessian *fh) {
    for (int i = 0; i < fh->total_layers; ++i) {
        compute_response_layer_Dyy_top(fh->response_map[i], fh->iimage);
    }
}

void compute_response_layers_Dyy_top_mid(struct fasthessian *fh) {
    for (int i = 0; i < fh->total_layers; ++i) {
        compute_response_layer_Dyy_top_mid(fh->response_map[i], fh->iimage);
    }
}

void compute_response_layers_unconditional(struct fasthessian* fh){
    for (int i = 0; i < fh->total_layers; ++i) {
		compute_response_layer_unconditional(fh->response_map[i], fh->iimage);
	}
}

void compute_response_layers_sonic_Dyy(struct fasthessian *fh) {
    for (int i = 0; i < fh->total_layers; ++i) {
        compute_response_layer_sonic_Dyy(fh->response_map[i], fh->iimage);
    }
}

void compute_response_layers_sonic_Dyy_unconditional(struct fasthessian *fh) {
    for (int i = 0; i < fh->total_layers; ++i) {
        compute_response_layer_sonic_Dyy_unconditional(fh->response_map[i], fh->iimage);
    }
}

void compute_response_layers_sonic_Dyy_unconditional_opt(struct fasthessian *fh) {
    for (int i = 0; i < fh->total_layers; ++i) {
        compute_response_layer_sonic_Dyy_unconditional_opt(fh->response_map[i], fh->iimage);
    }
}

// void compute_response_layers_sonic_Dyy_unconditional_opt_naive(struct fasthessian *fh) {
//     for (int i = 0; i < fh->total_layers; ++i) {
//         compute_response_layer_sonic_Dyy_unconditional_opt_naive(fh->response_map[i], fh->iimage);
//     }
// }

/* Dyy coords
// whole box filter
r00 = x - border;
r01 = x + border;  // -> x-border+filtersize-1 = x-border+(border+border+1)-1
c00 = y - lobe;
c01 = y + lobe - 1; // -> y-lobe+1 + 2*lobe -1 -1

A =
B =
C =
D =

// neg part box filter
r10 = x - lobe / 2 - 1;
r11 = r10 + lobe;
c10 = y - lobe;
c11 = y + lobe - 1;

A = (r10, c10)
B = (r10, c11) = (x - lobe / 2 - 1, y + lobe - 1)
C = (r11, c10)
D = (r11, c11)
*/

void compute_response_layer_sonic_Dyy(struct response_layer *layer, struct integral_image *iimage) {

    int height = layer->height;
    int width = layer->width;

    int data_width = iimage->data_width;

    int iwidth = iimage->width;  // TODO: fix where needs to be fixed
    int iheight = iimage->height;

    float *data = (float *) iimage->data;
    int step = layer->step;
    int filter_size = layer->filter_size;
    int border = (filter_size - 1) / 2;
    int lobe = filter_size / 3;
    float inv_area = 1.f/(filter_size*filter_size);

    float Dxx, Dyy, Dxy, Dyy0, Dyy1, A, B, C, D;
    float A0, A1, B0, B1, C0, C1, D0, D1;
    int r10, r11, c10, c11, r00, r01, c00, c01;
    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int ind = 0;

    int i, j, k, k0;
    int x = 0;
    int y = 0;

    // 1. Case The filter is smaller than the image
    if (filter_size <= iheight) {
        // Split the image into 9 cases - corners, borders and middle part.
        compute_response_layer_Dyy_laplacian_localityloops_inlined(layer, iimage);

    } else {
        // 2. Case The filter is somewhat larger than the image
        if (iheight > border) {

            // 2.1. D is sometimes outside the image. Blue lines edition brings the sonic into Dyy
            // Idea: Do compute_response_layer_Dyy_leftcorner
            // but everytime all corners are outside, we just use row values above
            if (iwidth > 2 * lobe - 1) {
                height_greater_border_width_greater_double_lobe_Dyy_inlined(layer, iimage);
            } else {
                height_greater_border_width_less_double_lobe_Dyy_inlined(layer, iimage);
            }

        } else {
            // Case 2.2 D is always outside the image.
            // Half the filter height is longer than the image.
            if (iwidth <= lobe) {
                // Case 2.2a the filter is longer and wider than the image.

                // D is right bottom corner or to the right (possibly below image)
                // A, B, C = 0
                D0 = data[(iheight-1) * data_width + (iwidth-1)];

                // Differentiate 2 cases
                // for the negative part:
                // 1.
                // is negative / inner part of Dyy box filter completely too big
                // i.e. is inner D on the right bottom corner or even bigger (to the right or below)
                // SADLY THIS CASE NEVER HAPPENS. CRY :(

                // if (iheight <= lobe/2-1) {
                // D1 = D0;
                //
                // // combine Dyy
                // Dyy = D0 - 3 * D1;
                // Dyy *= inv_area;
                //
                // for (i = 0; i < height * step; i += step) {
                //
                //     for (j = 0; j < width * step; j += step) {
                //         // Image coordinates
                //         x = i;
                //         y = j;
                //
                //         // Calculate Dxx, Dyy, Dxy with Box Filter
                //         Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                //                 - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
                //         Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                //                 + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                //                 - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                //                 - box_integral(iimage, x + 1, y + 1, lobe, lobe);
                //
                //         // Normalize Responses with inverse area
                //         Dyy *= inv_area;
                //         Dxy *= inv_area;
                //
                //         // Calculate Determinant
                //         response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;
                //
                //         // Calculate Laplacian
                //         laplacian[ind] = Dxx + Dyy >= 0;
                //         ind += 1;
                //     }
                // }

                // 2. Case
                // Only for - image_size: 32, filter_size: 99
                // inner D is in last col but the inner part is not too big
                // col is irrelevant, B and D have always the same value
                // i.e. inner loop is irrelevant for coords and values

                // We divide it again in 3 parts.
                // 1. B outside, D inside
                // 2. B outside, D outside
                // 3. B inside, D outside

                // 1. Case: B outside, D inside
                for (x = 0; x < width*step-lobe/2; x += step) {
                    // negative part
                    r10 = x - lobe / 2 - 1;
                    r11 = r10 + lobe;  // TODO fixer

                    D1 = data[r11 * data_width + (iwidth-1)];

                    // Compute Dyy
                    Dyy = D0 - 3 * D1;
                    Dyy *= inv_area;

                    for (y = 0; y < width * step; y += step) {

                        // Calculate Dxx, Dyy, Dxy with Box Filter
                        Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                                - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
                        Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                                + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                                - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                                - box_integral(iimage, x + 1, y + 1, lobe, lobe);

                        // Normalize Responses with inverse area
                        Dxx *= inv_area;
                        Dxy *= inv_area;

                        // Calculate Determinant
                        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

                        // Calculate Laplacian
                        laplacian[ind] = Dxx + Dyy >= 0;
                        ind += 1;
                    }
                }

                // 2. Case: B outside, D outside
                Dyy = - 2 * D0;
                Dyy *= inv_area;

                for (; x < lobe/2+1; x += step) {

                    for (y = 0; y < width*step; y += step) {

                        // Calculate Dxx, Dyy, Dxy with Box Filter
                        Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                                - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
                        Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                                + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                                - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                                - box_integral(iimage, x + 1, y + 1, lobe, lobe);

                        // Normalize Responses with inverse area
                        Dxx *= inv_area;
                        Dxy *= inv_area;

                        // Calculate Determinant
                        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

                        // Calculate Laplacian
                        laplacian[ind] = Dxx + Dyy >= 0;
                        ind += 1;
                    }
                }

                // 3. Case: B inside, D outside
                for (; x < height * step; x += step) {

                    r10 = x - lobe / 2 - 1;

                    B = data[r10 * data_width + (iwidth-1)];
                    // D1 = D0 - B;
                    // Dyy = D0 - 3 * D1;
                    // Dyy *= inv_area;

                    // this can be simplified to:
                    // (small difference in original is this but error is <0.000001 eps)
                    Dyy = 3*B - 2*D0;
                    Dyy *= inv_area;

                    for (y = 0; y < width * step; y += step) {

                        // Calculate Dxx, Dyy, Dxy with Box Filter
                        Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                                - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
                        Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                                + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                                - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                                - box_integral(iimage, x + 1, y + 1, lobe, lobe);

                        // Normalize Responses with inverse area
                        Dxx *= inv_area;
                        Dxy *= inv_area;

                        // Calculate Determinant
                        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

                        // Calculate Laplacian
                        laplacian[ind] = Dxx + Dyy >= 0;
                        ind += 1;
                    }
                }

            } else {
                /* 2.2.b Half the filter is longer than the image, but narrower.
                // HERE: width > lobe and height <= border
                // Here only the whole filter is optimized, not negative part

                // Only happens for filter: 75 and image: 32

                // A, B = 0
                // D is outside (i.e. below) the image, but the columns change.
                // => D = [height-1, some_column]
                // C = 0 or [height-1, some_column]
                */

                // Create array for Dyy that has image width length
                // all rows (for big part) have same Dyy values
                float Dyy_arr[iwidth];  // stack is faster than heap

                // C = 0 and D = [height-1, some_column]
                // from y = 0 until D is (exclusive) in last column
                for (i = 0; i < width*step-lobe; i += step) {  // 0 - 7
                    // C = 0
                    // D = [height-1, i+lobe-1]
                    D = data[(iheight-1) * data_width + (i+lobe-1)];
                    Dyy_arr[i] = D;
                }

                // only bottom left corner value needed
                D = data[(iheight-1) * data_width + (iwidth-1)];

                // C = 0 and D = [height-1, width-1] (below or right of bottom corner)
                // C is still outside and D now too
                for (; i < lobe; i += step) {  // 7 - 25
                    Dyy_arr[i] = D;
                }

                // if y = lobe, then C = [height-1, 0]
                // C = [height-1, some_column] and D = [height-1, width-1] (below or right of bottom corner)
                for (; i < width*step; i += step) {  // 25 - 32
                    // C = [height-1, i-lobe]
                    // D = [height-1, width-1]
                    C = data[(iheight-1) * data_width + (i-lobe)];
                    Dyy_arr[i] = D - C;
                }

                // Use precomputation for faster compute
                for (x = 0; x < height*step; x += step) {

                    for (y = 0; y < width*step; y += step) {

                        // Calculate Dxx, Dyy, Dxy with Box Filter
                        Dyy = Dyy_arr[y] - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);
                        Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                                - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
                        Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                                + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                                - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                                - box_integral(iimage, x + 1, y + 1, lobe, lobe);

                        // Normalize Responses with inverse area
                        Dyy *= inv_area;
                        Dxx *= inv_area;
                        Dxy *= inv_area;

                        // Calculate Determinant
                        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

                        // Calculate Laplacian
                        laplacian[ind] = Dxx + Dyy >= 0;

                        // Increment index
                        ind += 1;
                    }
                }
            }
        }
    }
}

void compute_response_layer_sonic_Dyy_unconditional(struct response_layer *layer, struct integral_image *iimage) {

    int height = layer->height;
    int width = layer->width;

    int data_width = iimage->data_width;

    int iwidth = iimage->width;
    int iheight = iimage->height;

    float *data = (float *) iimage->data;
    int step = layer->step;
    int filter_size = layer->filter_size;
    int border = (filter_size - 1) / 2;
    int lobe = filter_size / 3;
    float inv_area = 1.f/(filter_size*filter_size);

    float Dxx, Dyy, Dxy, Dyy0, Dyy1, A, B, C, D;
    float A0, A1, B0, B1, C0, C1, D0, D1;
    int r10, r11, c10, c11, r00, r01, c00, c01;
    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int ind = 0;

    int k, k0, i, j;
    int x = 0;
    int y = 0;

    // 1. Case The filter is smaller than the image
    if (filter_size <= iheight) {
        // Split the image into 9 cases - corners, borders and middle part.
        compute_response_layer_Dyy_laplacian_locality_uncond_opt_flops_inlined(layer, iimage);

    } else {
        // 2. Case The filter is somewhat larger than the image
        if (iheight > border) {

            // 2.1. D is sometimes outside the image. Blue lines edition brings the sonic into Dyy
            // Idea: Do compute_response_layer_Dyy_leftcorner
            // but everytime all corners are outside, we just use row values above
            if (iwidth > 2 * lobe - 1) {
                height_greater_border_width_greater_double_lobe_Dyy_inlined(layer, iimage);
            } else {
                height_greater_border_width_less_double_lobe_Dyy_inlined(layer, iimage);
            }

        } else {
            // Case 2.2 D is always outside the image.
            // Half the filter height is longer than the image.
            if (iwidth <= lobe) {
                // Case 2.2a the filter is longer and wider than the image.

                // D is right bottom corner or to the right (possibly below image)
                // A, B, C = 0
                D0 = data[(iheight-1) * data_width + (iwidth-1)];

                // Differentiate 2 cases
                // for the negative part:
                // 1.
                // is negative / inner part of Dyy box filter completely too big
                // i.e. is inner D on the right bottom corner or even bigger (to the right or below)
                // SADLY THIS CASE NEVER HAPPENS. CRY :(

                // if (iheight <= lobe/2-1) {
                // D1 = D0;
                //
                // // combine Dyy
                // Dyy = D0 - 3 * D1;
                // Dyy *= inv_area;
                //
                // for (i = 0; i < height * step; i += step) {
                //
                //     for (j = 0; j < width * step; j += step) {
                //         // Image coordinates
                //         x = i;
                //         y = j;
                //
                //         // Calculate Dxx, Dyy, Dxy with Box Filter
                //         Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                //                 - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
                //         Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                //                 + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                //                 - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                //                 - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);
                //
                //         // Normalize Responses with inverse area
                //         Dyy *= inv_area;
                //         Dxy *= inv_area;
                //
                //         // Calculate Determinant
                //         response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;
                //
                //         // Calculate Laplacian
                //         laplacian[ind] = Dxx + Dyy >= 0;
                //         ind += 1;
                //     }
                // }

                // 2. Case
                // Only for - image_size: 32, filter_size: 99
                // inner D is in last col but the inner part is not too big
                // col is irrelevant, B and D have always the same value
                // i.e. inner loop is irrelevant for coords and values

                // We divide it again in 3 parts.
                // 1. B outside, D inside
                // 2. B outside, D outside
                // 3. B inside, D outside

                // 1. Case: B outside, D inside
                for (x = 0; x < width*step-lobe/2; x += step) {
                    // negative part
                    r10 = x - lobe / 2 - 1;
                    r11 = r10 + lobe;  // TODO: fix this

                    D1 = data[r11 * data_width + (iwidth-1)];

                    // Compute Dyy
                    Dyy = D0 - 3 * D1;
                    Dyy *= inv_area;

                    for (y = 0; y < width * step; y += step) {
                        // Calculate Dxx, Dyy, Dxy with Box Filter
                        Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                                - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
                        Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                                + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                                - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                                - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

                        // Normalize Responses with inverse area
                        Dxx *= inv_area;
                        Dxy *= inv_area;

                        // Calculate Determinant
                        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

                        // Calculate Laplacian
                        laplacian[ind] = Dxx + Dyy >= 0;
                        ind += 1;
                    }
                }

                // 2. Case: B outside, D outside
                Dyy = - 2 * D0;
                Dyy *= inv_area;

                for (; x < lobe/2+1; x += step) {

                    for (y = 0; y < width*step; y += step) {
                        // Calculate Dxx, Dyy, Dxy with Box Filter
                        Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                                - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
                        Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                                + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                                - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                                - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

                        // Normalize Responses with inverse area
                        Dxx *= inv_area;
                        Dxy *= inv_area;

                        // Calculate Determinant
                        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

                        // Calculate Laplacian
                        laplacian[ind] = Dxx + Dyy >= 0;
                        ind += 1;
                    }
                }

                // 3. Case: B inside, D outside
                for (; x < height * step; x += step) {
                    r10 = x - lobe / 2 - 1;

                    B = data[r10 * data_width + (iwidth-1)];
                    // D1 = D0 - B;
                    // Dyy = D0 - 3 * D1;
                    // Dyy *= inv_area;

                    // this can be simplified to:
                    // (small difference in original is this but error is <0.000001 eps)
                    Dyy = 3*B - 2*D0;
                    Dyy *= inv_area;

                    for (y = 0; y < width * step; y += step) {
                        // Calculate Dxx, Dyy, Dxy with Box Filter
                        Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                                - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
                        Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                                + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                                - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                                - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

                        // Normalize Responses with inverse area
                        Dxx *= inv_area;
                        Dxy *= inv_area;

                        // Calculate Determinant
                        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

                        // Calculate Laplacian
                        laplacian[ind] = Dxx + Dyy >= 0;
                        ind += 1;
                    }
                }

            } else {
                /* 2.2.b Half the filter is longer than the image, but narrower.
                // HERE: width > lobe and height <= border
                // Here only the whole filter is optimized, not negative part

                // Only happens for filter: 75 and image: 32

                // A, B = 0
                // D is outside (i.e. below) the image, but the columns change.
                // => D = [height-1, some_column]
                // C = 0 or [height-1, some_column]
                */

                // Create array for Dyy that has image width length
                // all rows (for big part) have same Dyy values
                float Dyy_arr[iwidth];  // stack is faster than heap

                // C = 0 and D = [height-1, some_column]
                // from y = 0 until D is (exclusive) in last column
                for (i = 0; i < width*step-lobe; i += step) {  // 0 - 7
                    // C = 0
                    // D = [height-1, i+lobe-1]
                    D = data[(iheight-1) * data_width + (i+lobe-1)];
                    Dyy_arr[i] = D;
                }

                // only bottom left corner value needed
                D = data[(iheight-1) * data_width + (iwidth-1)];

                // C = 0 and D = [height-1, width-1] (below or right of bottom corner)
                // C is still outside and D now too
                for (; i < lobe; i += step) {  // 7 - 25
                    Dyy_arr[i] = D;
                }

                // if y = lobe, then C = [height-1, 0]
                // C = [height-1, some_column] and D = [height-1, width-1] (below or right of bottom corner)
                for (; i < width*step; i += step) {  // 25 - 32
                    // C = [height-1, i-lobe]
                    // D = [height-1, width-1]
                    C = data[(iheight-1) * data_width + (i-lobe)];
                    Dyy_arr[i] = D - C;
                }

                // Use precomputation for faster compute
                for (x = 0; x < height * step; x += step) {

                    for (y = 0; y < width * step; y += step) {
                        // Calculate Dxx, Dyy, Dxy with Box Filter
                        Dyy = Dyy_arr[y] - 3 * box_integral_unconditional(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);
                        Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                                - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
                        Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                                + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                                - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                                - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

                        // Normalize Responses with inverse area
                        Dyy *= inv_area;
                        Dxx *= inv_area;
                        Dxy *= inv_area;

                        // Calculate Determinant
                        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

                        // Calculate Laplacian
                        laplacian[ind] = Dxx + Dyy >= 0;

                        // Increment index
                        ind += 1;
                    }
                }
            }
        }
    }
}

void compute_response_layer_sonic_Dyy_unconditional_opt(struct response_layer *layer, struct integral_image *iimage) {

    int height = layer->height;
    int width = layer->width;

    int data_width = iimage->data_width;

    int iwidth = iimage->width;
    int iheight = iimage->height;

    float *data = (float *) iimage->data;
    int step = layer->step;
    int filter_size = layer->filter_size;
    int border = (filter_size - 1) / 2;
    int lobe = filter_size / 3;
    float inv_area = 1.f/(filter_size*filter_size);

    float Dxx, Dyy, Dxy, Dyy0, Dyy1, A, B, C, D;
    float A0, A1, B0, B1, C0, C1, D0, D1;
    int r10, r11, c10, c11, r00, r01, c00, c01;
    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int ind = 0;

    int k, k0, i, j, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
    int x = 0;
    int y = 0;

    // 1. Case The filter is smaller than the image
    if (filter_size <= iheight) {
        // Split the image into 9 cases - corners, borders and middle part.
        compute_response_layer_Dyy_laplacian_locality_uncond_opt_flops_inlined(layer, iimage);

    } else {
        // 2. Case The filter is somewhat larger than the image
        if (iheight > border) {

            // 2.1. D is sometimes outside the image. Blue lines edition brings the sonic into Dyy
            // Idea: Do compute_response_layer_Dyy_leftcorner
            // but everytime all corners are outside, we just use row values above
            if (iwidth > 2 * lobe - 1) {
                height_greater_border_width_greater_double_lobe_Dyy_inlined(layer, iimage);
            } else {
                height_greater_border_width_less_double_lobe_Dyy_inlined(layer, iimage);
            }

        } else {
            // Case 2.2 D is always outside the image.
            // Half the filter height is longer than the image.
            if (iwidth <= lobe) {
                // Case 2.2a the filter is longer and wider than the image.

                // D is right bottom corner or to the right (possibly below image)
                // A, B, C = 0
                D0 = data[(iheight-1) * data_width + (iwidth-1)];

                // Differentiate 2 cases
                // for the negative part:
                // 1.
                // is negative / inner part of Dyy box filter completely too big
                // i.e. is inner D on the right bottom corner or even bigger (to the right or below)
                // SADLY THIS CASE NEVER HAPPENS. CRY :(

                // if (iheight <= lobe/2-1) {
                // D1 = D0;
                //
                // // combine Dyy
                // Dyy = D0 - 3 * D1;
                // Dyy *= inv_area;
                //
                // for (i = 0; i < height * step; i += step) {
                //
                //     for (j = 0; j < width * step; j += step) {
                //         // Image coordinates
                //         x = i;
                //         y = j;
                //
                //         // Calculate Dxx, Dyy, Dxy with Box Filter
                //         Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                //                 - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
                //         Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                //                 + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                //                 - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                //                 - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);
                //
                //         // Normalize Responses with inverse area
                //         Dyy *= inv_area;
                //         Dxy *= inv_area;
                //
                //         // Calculate Determinant
                //         response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;
                //
                //         // Calculate Laplacian
                //         laplacian[ind] = Dxx + Dyy >= 0;
                //         ind += 1;
                //     }
                // }

                // 2. Case
                // Only for - image_size: 32, filter_size: 99
                // inner D is in last col but the inner part is not too big
                // col is irrelevant, B and D have always the same value
                // i.e. inner loop is irrelevant for coords and values

                // We divide it again in 3 parts.
                // 1. B outside, D inside
                // 2. B outside, D outside
                // 3. B inside, D outside

                // 1. Case: B outside, D inside
                for (x = 0; x < width*step-lobe/2; x += step) {
                    // negative part
                    r10 = x - lobe / 2 - 1;
                    r11 = r10 + lobe;  // TODO: fix this

                    D1 = data[r11 * data_width + (iwidth-1)];

                    // Compute Dyy
                    Dyy = D0 - 3 * D1;
                    Dyy *= inv_area;

                    // precompute
                    t0 = x - lobe;
                    t1 = x + lobe - 1;
                    t2 = lobe/2 + 1;
                    t3 = x - 1;
                    t5 = x - lobe - 1;
                    t7 = x + lobe;

                    for (y = 0; y < width * step; y += step) {
                        // box_integral_unconditional precompute
                        t4 = y - 1;
                        t6 = y - lobe - 1;
                        t8 = y + lobe;

                        // Compute Dxx, Dxy
                        Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                            - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
                        Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                                + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                                - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                                - box_integral_unconditional_opt(iimage, x, y, t7, t8);

                        // Normalize Responses with inverse area
                        Dxx *= inv_area;
                        Dxy *= inv_area;

                        // Calculate Determinant
                        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

                        // Calculate Laplacian
                        laplacian[ind] = Dxx + Dyy >= 0;
                        ind += 1;
                    }
                }

                // 2. Case: B outside, D outside
                Dyy = - 2 * D0;
                Dyy *= inv_area;

                for (; x < lobe/2+1; x += step) {

                    // precompute
                    t0 = x - lobe;
                    t1 = x + lobe - 1;
                    t2 = lobe/2 + 1;
                    t3 = x - 1;
                    t5 = x - lobe - 1;
                    t7 = x + lobe;

                    for (y = 0; y < width*step; y += step) {
                        // box_integral_unconditional precompute
                        t4 = y - 1;
                        t6 = y - lobe - 1;
                        t8 = y + lobe;

                        // Compute Dxx, Dxy
                        Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                            - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
                        Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                                + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                                - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                                - box_integral_unconditional_opt(iimage, x, y, t7, t8);

                        // Normalize Responses with inverse area
                        Dxx *= inv_area;
                        Dxy *= inv_area;

                        // Calculate Determinant
                        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

                        // Calculate Laplacian
                        laplacian[ind] = Dxx + Dyy >= 0;
                        ind += 1;
                    }
                }

                // 3. Case: B inside, D outside
                for (; x < height * step; x += step) {
                    r10 = x - lobe / 2 - 1;

                    B = data[r10 * data_width + (iwidth-1)];
                    // D1 = D0 - B;
                    // Dyy = D0 - 3 * D1;
                    // Dyy *= inv_area;

                    // this can be simplified to:
                    // (small difference in original is this but error is <0.000001 eps)
                    Dyy = 3*B - 2*D0;
                    Dyy *= inv_area;

                    // precompute
                    t0 = x - lobe;
                    t1 = x + lobe - 1;
                    t2 = lobe/2 + 1;
                    t3 = x - 1;
                    t5 = x - lobe - 1;
                    t7 = x + lobe;

                    for (y = 0; y < width * step; y += step) {
                        // box_integral_unconditional precompute
                        t4 = y - 1;
                        t6 = y - lobe - 1;
                        t8 = y + lobe;

                        // Compute Dxx, Dxy
                        Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                            - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
                        Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                                + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                                - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                                - box_integral_unconditional_opt(iimage, x, y, t7, t8);

                        // Normalize Responses with inverse area
                        Dxx *= inv_area;
                        Dxy *= inv_area;

                        // Calculate Determinant
                        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

                        // Calculate Laplacian
                        laplacian[ind] = Dxx + Dyy >= 0;
                        ind += 1;
                    }
                }

            } else {
                /* 2.2.b Half the filter is longer than the image, but narrower.
                // HERE: width > lobe and height <= border
                // Here only the whole filter is optimized, not negative part

                // Only happens for filter: 75 and image: 32

                // A, B = 0
                // D is outside (i.e. below) the image, but the columns change.
                // => D = [height-1, some_column]
                // C = 0 or [height-1, some_column]
                */

                // Create array for Dyy that has image width length
                // all rows (for big part) have same Dyy values
                float Dyy_arr[iwidth];  // stack is faster than heap

                // C = 0 and D = [height-1, some_column]
                // from y = 0 until D is (exclusive) in last column
                for (i = 0; i < width*step-lobe; i += step) {  // 0 - 7
                    // C = 0
                    // D = [height-1, i+lobe-1]
                    D = data[(iheight-1) * data_width + (i+lobe-1)];
                    Dyy_arr[i] = D;
                }

                // only bottom left corner value needed
                D = data[(iheight-1) * data_width + (iwidth-1)];

                // C = 0 and D = [height-1, width-1] (below or right of bottom corner)
                // C is still outside and D now too
                for (; i < lobe; i += step) {  // 7 - 25
                    Dyy_arr[i] = D;
                }

                // if y = lobe, then C = [height-1, 0]
                // C = [height-1, some_column] and D = [height-1, width-1] (below or right of bottom corner)
                for (; i < width*step; i += step) {  // 25 - 32
                    // C = [height-1, i-lobe]
                    // D = [height-1, width-1]
                    C = data[(iheight-1) * data_width + (i-lobe)];
                    Dyy_arr[i] = D - C;
                }

                // Use precomputation for faster compute
                for (x = 0; x < height * step; x += step) {

                    // precompute
                    t0 = x - lobe;
                    t1 = x + lobe - 1;
                    t2 = lobe/2 + 1;
                    t3 = x - 1;
                    t5 = x - lobe - 1;
                    t7 = x + lobe;
                    t9 = x - lobe / 2 - 1;

                    for (y = 0; y < width * step; y += step) {
                        // box_integral_unconditional precompute
                        t4 = y - 1;
                        t6 = y - lobe - 1;
                        t8 = y + lobe;

                        // Calculate Dxx, Dyy, Dxy with Box Filter
                        Dyy = Dyy_arr[y] - 3 * box_integral_unconditional_opt(iimage, t9, t6 + 1, t9 + lobe , t8 - 2);

                        // Compute Dxx, Dxy
                        Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                            - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
                        Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                                + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                                - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                                - box_integral_unconditional_opt(iimage, x, y, t7, t8);

                        // Normalize Responses with inverse area
                        Dyy *= inv_area;
                        Dxx *= inv_area;
                        Dxy *= inv_area;

                        // Calculate Determinant
                        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

                        // Calculate Laplacian
                        laplacian[ind] = Dxx + Dyy >= 0;

                        // Increment index
                        ind += 1;
                    }
                }
            }
        }
    }
}

void compute_response_layer_sonic_Dyy_unconditional_opt_naive(struct response_layer *layer, struct integral_image *iimage) {

    int iheight = iimage->height;
    int filter_size = layer->filter_size;

    // 1. Case The filter is smaller than the image
    if (filter_size <= iheight) {
        // Split the image into 9 cases - corners, borders and middle part.
        compute_response_layer_Dyy_laplacian_locality_uncond_opt_flops_inlined(layer, iimage);

    } else {

        compute_response_layer_unconditional(layer, iimage);
    }
}

void height_greater_border_width_greater_double_lobe_Dyy(struct response_layer *layer, struct integral_image *iimage) {
    // Filter_size > height
    // 2. Case The filter is somewhat larger than the image
    /*  (height > border && (iwidth > 2 * lobe - 1))
          2.1. D is sometimes outside the image. Blue lines edition brings the sonic into Dyy
          Idea: Do compute_response_layer_Dyy_leftcorner
          but everytime all corners are outside, we just use row values above

          Top: Left - Middle (just one case) - Bottom
          Last Row Before Blue Lines: Left - Middle (just one case) - Bottom
          Blue Lines: Use previously stored values
          Bottom: Left - Middle (just one case) - Bottom

          */

    float Dxx, Dyy, Dxy;
    int x, y, k, k0, k1;
    int r00, r01, c00, c01, r10, r11, c10, c11;
    float Dyy0, Dyy1, A, B, C, D;

    float *response = layer->response;
    bool *laplacian = layer->laplacian;

    int data_width = iimage->data_width;
    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size / 3;
    int border = (filter_size - 1) / 2;
    float inv_area = 1.f / (filter_size * filter_size);

    int ind = 0;

    float *data = (float *)iimage->data;  // brauch hier keinen cast weil es eig float sein sollte
    int iheight = iimage->height;
    int iwidth = iimage->width;

    float dyy_row_before_blue[width];

    /****************
     *   TOP
     *****************/
    int i = 0;
    for (; i < height * step - border - 1; i += step) {
        int j = 0;
        // Top Left Corner
        for (; j < lobe; j += step) {
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;
            // A, B, C outside A, B, C = 0
            // D inside
            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0);
            ind += 1;
        }

        // Top Mid
        for (; j < width * step - lobe + 1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            // A, B outside A, B = 0
            // C, D inside
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }

        // Top Right
        for (; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            // A, B outside A, B = 0
            // C inside
            // D outside D = last element of row
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth - 1];
            Dyy0 = D - C;

            // neg part box filter
            Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0);
            ind += 1;
        }
    }

    /*********************************
     *    LAST ROW BEFORE BLUE LINE
     *********************************/
    int idx_store_row = ((height * step - border - 1 + step - 1) / step) * step;
    i = idx_store_row;
    int counter = 0;

    int j = 0;
    // Last Row Before Blue Line Left
    for (; j < lobe; j += step) {
        x = i;
        y = j;

        // Compute Dyy  
        // whole box filter
        c01 = y + lobe - 1;

        // A, B, C outside A, B, C = 0
        // D below D = last value of column.
        Dyy0 = data[(iheight - 1) * data_width + c01];

        // Store value for blue line.
        dyy_row_before_blue[counter] = Dyy0;

        // neg part box filter
        Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);
        // Compute Dxx, Dxy
        Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
              3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
        Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) + box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
              box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

        // Normalize Responses with inverse area
        Dxx *= inv_area;
        Dyy *= inv_area;
        Dxy *= inv_area;

        // Calculate Determinant
        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

        // Calculate Laplacian
        laplacian[ind] = (Dxx + Dyy >= 0);
        ind += 1;
        counter += 1;
    }

    // Last Row Before Blue Line Mid
    for (; j < width * step - lobe + 1; j += step) {
        // Image coordinates
        x = i;
        y = j;

        // Compute Dyy
        // whole box filter
        c00 = y - lobe;
        c01 = y + lobe - 1;

        // A, B outside A, B = 0
        // C and D below C, D last element of column.
        C = data[(iheight - 1) * data_width + c00];
        D = data[(iheight - 1) * data_width + c01];

        Dyy0 = D - C;

        // Store value for blue line.
        dyy_row_before_blue[counter] = Dyy0;

        // neg part box filter
        Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

        // Compute Dxx, Dxy
        Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
              3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
        Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) + box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
              box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

        // Normalize Responses with inverse area
        Dxx *= inv_area;
        Dyy *= inv_area;
        Dxy *= inv_area;

        // Calculate Determinant
        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

        // Calculate Laplacian
        laplacian[ind] = (Dxx + Dyy >= 0);
        ind += 1;
        counter += 1;
    }

    // Last Row Before Blue Line Right
    for (j; j < width * step; j += step) {
        // Image coordinates
        x = i;
        y = j;

        // whole filter
        c00 = y - lobe;

        // A, B outside A, B = 0
        // C below C = last element of column
        // D outside D = last element.
        C = data[(iheight - 1) * data_width + c00];
        D = data[(iheight - 1) * data_width + iwidth - 1];
        Dyy0 = D - C;

        // Store value for blue line.
        dyy_row_before_blue[counter] = Dyy0;

        // neg part box filter
        Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

        // Compute Dxx, Dxy
        Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
              3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
        Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) + box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
              box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

        // Normalize Responses with inverse area
        Dxx *= inv_area;
        Dyy *= inv_area;
        Dxy *= inv_area;

        // Calculate Determinant
        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

        // Calculate Laplacian
        laplacian[ind] = (Dxx + Dyy >= 0);
        ind += 1;
        counter += 1;
    }
    i += step;

    /*****************
     *    BLUE LINE
     *****************/

    for (i; i < (border + 1); i += step) {
        int counter = 0;
        for (int j = 0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            Dyy0 = dyy_row_before_blue[counter];

            // neg part box filter
            Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
            counter += 1;
        }
    }

    /****************
     *   BOTTOM
     *****************/

    // Bottom
    for (; i < height * step; i += step) {
        int j = 0;
        // Bottom Left
        for (; j < lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            r00 = x - border - 1;
            c01 = y + lobe - 1;

            // A, C outside A, C = 0
            // B inside
            // D below D = last element of column
            B = data[r00 * data_width + c01];
            D = data[(iheight - 1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0);
            ind += 1;
        }

        // Bottom Mid
        for (j; j < width * step - lobe + 1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            r00 = x - border - 1;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            // A, B inside
            // C, D below C, D last element of column
            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[(iheight - 1) * data_width + c00];
            D = data[(iheight - 1) * data_width + c01];

            Dyy0 = A - B - C + D;

            // neg part box filter
            Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }

        // Bottom Right
        for (j; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            r00 = x - border - 1;
            c00 = y - lobe;
            // A inside
            // B outside B = last element of row
            // C below C = last element of column
            // D outside D = last element
            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth - 1];
            C = data[(iheight - 1) * data_width + c00];
            D = data[(iheight - 1) * data_width + iwidth - 1];

            Dyy0 = A - B - C + D;

            // neg part box filter
            Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }
}

void height_greater_border_width_less_double_lobe_Dyy(struct response_layer *layer, struct integral_image *iimage) {
    // Filter_size > height
    // 2. Case The filter is somewhat larger than the image
    /*  (height > border && (iwidth < 2 * lobe))
          2.1. D is sometimes outside the image. Blue lines edition brings the sonic into Dyy
          Idea: Do compute_response_layer_Dyy_leftcorner
          but everytime all corners are outside, we just use row values above

          Top: Left - Middle (just one case) - Bottom
          Last Row Before Blue Lines: Left - Middle (just one case) - Bottom
          Blue Lines: Use previously stored values
          Bottom: Left - Middle (just one case) - Bottom

          */

    float Dxx, Dyy, Dxy;
    int x, y, k, k0, k1;
    int r00, r01, c00, c01, r10, r11, c10, c11;
    float Dyy0, Dyy1, A, B, C, D;

    float *response = layer->response;
    bool *laplacian = layer->laplacian;

    int data_width = iimage->data_width;
    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size / 3;
    int border = (filter_size - 1) / 2;
    float inv_area = 1.f / (filter_size * filter_size);

    int ind = 0;
    float *data = (float *)iimage->data;  // brauch hier keinen cast weil es eig float sein sollte
    int iheight = iimage->height;
    int iwidth = iimage->width;

    float dyy_row_before_blue[width];

    /****************
     *   TOP
     *****************/
    int i = 0;

    for (; i < height * step - border - 1; i += step) {
        int j = 0;
        // Top Left
        for (; j < width * step - lobe; j += step) {
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            // A, B, C outside A, B, C = 0
            // D inside
            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0);
            ind += 1;
        }

        // Top Mid

        // Image coordinates
        x = i;
        y = j;

        // whole filter
        r01 = x + border;

        // A, B, C outside A, B, C = 0
        // D outside, D = last element of row
        D = data[r01 * data_width + iwidth - 1];
        Dyy0 = D;

        // neg part box filter
        Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

        // Compute Dxx, Dxy
        Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
              3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
        Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) + box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
              box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

        // Normalize Responses with inverse area
        Dxx *= inv_area;
        Dyy *= inv_area;
        Dxy *= inv_area;

        // Calculate Determinant
        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

        // Calculate Laplacian
        laplacian[ind] = (Dxx + Dyy >= 0);
        ind += 1;
        j += step;

        // Top Right
        for (; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            // A, B outside A, B = 0
            // C inside
            // D outside D = last element of row
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth - 1];
            Dyy0 = D - C;

            // neg part box filter
            Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0);
            ind += 1;
        }
    }

    /*********************************
     *    LAST ROW BEFORE BLUE LINE
     *********************************/
    int idx_store_row = ((height * step - border - 1 + step - 1) / step) * step;
    i = idx_store_row;
    int counter = 0;

    int j = 0;
    // Last Row Before Blue Line Left
    for (; j < width * step - lobe; j += step) {
        x = i;
        y = j;

        // Compute Dyy  
        // whole box filter
        r01 = x + border;
        c01 = y + lobe - 1;

        // A, B, C outside A, B, C = 0
        // D below D = last value of column.
        Dyy0 = data[(iheight - 1) * data_width + c01];

        // Store value for blue line.
        dyy_row_before_blue[counter] = Dyy0;

        // neg part box filter
        Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);
        // Compute Dxx, Dxy
        Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
              3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
        Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) + box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
              box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

        // Normalize Responses with inverse area
        Dxx *= inv_area;
        Dyy *= inv_area;
        Dxy *= inv_area;

        // Calculate Determinant
        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

        // Calculate Laplacian
        laplacian[ind] = (Dxx + Dyy >= 0);
        ind += 1;
        counter += 1;
    }

    // Last Row Refore Blue Line Mid

    // Image coordinates
    x = i;
    y = j;

    // whole filter
    r01 = x + border;

    // A, B, C outside A, B, C = 0
    // D outside D = last element.
    D = data[(iheight - 1) * data_width + iwidth - 1];
    Dyy0 = D;
    // Store value for blue line.
    dyy_row_before_blue[counter] = Dyy0;

    // neg part box filter
    Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

    // Compute Dxx, Dxy
    Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
          3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
    Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) + box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
          box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

    // Normalize Responses with inverse area
    Dxx *= inv_area;
    Dyy *= inv_area;
    Dxy *= inv_area;

    // Calculate Determinant
    response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

    // Calculate Laplacian
    laplacian[ind] = (Dxx + Dyy >= 0);
    ind += 1;
    j += step;
    counter += 1;

    // Last Row Before Blue Line Right
    for (; j < width * step; j += step) {
        // Image coordinates
        x = i;
        y = j;

        // whole filter
        c00 = y - lobe;

        // A, B outside A, B = 0
        // C below C = last element of column
        // D outside D = last element.
        C = data[(iheight - 1) * data_width + c00];
        D = data[(iheight - 1) * data_width + iwidth - 1];
        Dyy0 = D - C;

        // Store value for blue line.
        dyy_row_before_blue[counter] = Dyy0;

        // neg part box filter
        Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

        // Compute Dxx, Dxy
        Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
              3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
        Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) + box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
              box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

        // Normalize Responses with inverse area
        Dxx *= inv_area;
        Dyy *= inv_area;
        Dxy *= inv_area;

        // Calculate Determinant
        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

        // Calculate Laplacian
        laplacian[ind] = (Dxx + Dyy >= 0);
        ind += 1;
        counter += 1;
    }
    i += step;

    /*****************
     *    BLUE LINE
     *****************/

    for (; i < (border + 1); i += step) {
        int counter = 0;
        for (int j = 0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            Dyy0 = dyy_row_before_blue[counter];

            // neg part box filter
            Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
            counter += 1;
        }
    }

    /****************
     *   BOTTOM
     *****************/

    // Bottom
    for(;i < height * step; i+=step) {
        int j = 0;
        // Bottom Left
        for (; j < width * step - lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            r00 = x - border - 1;
            c01 = y + lobe - 1;

            // A, C outside A, C = 0
            // B inside
            // D below D = last element of column
            B = data[r00 * data_width + c01];
            D = data[(iheight - 1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0);
            ind += 1;
        }

        // Bottom Mid

        // Image coordinates
        x = i;
        y = j;

        // whole filter
        r00 = x - border - 1;
        r01 = x + border;

        // A, C outside A, C = 0
        // B outside B = last element of row
        // D outside D = last element
        B = data[r00 * data_width + iwidth - 1];
        D = data[(iheight - 1) * data_width + iwidth - 1];
        Dyy0 = D - B;
        // Store value for blue line.
        dyy_row_before_blue[counter] = Dyy0;

        // neg part box filter
        Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

        // Compute Dxx, Dxy
        Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
              3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
        Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) + box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
              box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

        // Normalize Responses with inverse area
        Dxx *= inv_area;
        Dyy *= inv_area;
        Dxy *= inv_area;

        // Calculate Determinant
        response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

        // Calculate Laplacian
        laplacian[ind] = (Dxx + Dyy >= 0);
        ind += 1;
        j += step;

        // Bottom Right
        for (; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            r00 = x - border - 1;
            c00 = y - lobe;

            // A inside
            // B outside B = last element of row
            // C below C = last element of column
            // D outside D = last element
            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth - 1];
            C = data[(iheight - 1) * data_width + c00];
            D = data[(iheight - 1) * data_width + iwidth - 1];

            Dyy0 = A - B - C + D;

            // neg part box filter
            Dyy = Dyy0 - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }
}

// Only for images >= 128x128
void compute_response_layer_Dyy_leftcorner(struct response_layer *layer, struct integral_image *iimage) {
    float Dxx, Dyy, Dxy;
    int x, y;
    int r0, r1, c0, c1, r00, r01, c00, c01, r10, r11, c10, c11;
    float Dyy0, Dyy1;
    float A, B, C, D;
    int k;

    float *response = layer->response;
    bool *laplacian = layer->laplacian;

    int data_width = iimage->data_width;
    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size / 3;
    int border = (filter_size - 1) / 2;
    float inv_area = 1.f / (filter_size * filter_size);

    int ind = 0;  // oder alternativ (i+1)*j

    float *data = (float *) iimage->data;  // brauch hier keinen cast weil es eig float sein sollte
    int iheight = iimage->height;
    int iwidth = iimage->width;

    /*
        In our optimized function, i and j are raw coordinates in the integral
        / original image. In the original function these coordinates are
        coordinates of the layer. We changed it this way in the optimized function
        because its easier to about corner cases in raw coordinates.
    */

    // Top Left Corner - Case 1: B of neg part outside
    for (int i = 0; i < lobe / 2 + 1; i += step) {  // Inner B is outside, i.e. 0
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {  // c0 = col - 1 = (y - lobe + 1) -1 < 0
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];  // seg fault at: i: 12, j:0

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;  // -1 is already part of r10
            c11 = y + lobe - 1;

            D = data[r11 * data_width + c11];
            Dyy1 = D;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Left Corner - Case 2: B of neg part inside
    // initial value has to be rounded up to next bigger step
    k = (lobe / 2 + 1 + step - 1) / step * step;

    for (int i = k; i < border + 1; i += step) {  // Inner B is inside
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;   // -1 is already part of r10
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];
            Dyy1 = D - B;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Rest
    k = (lobe + step - 1) / step * step;

    for (int i = 0; i < border + 1; i += step) {
        ind = (i / step) * width + (k / step);

        for (int j = k; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dyy = box_integral(iimage, x - border, y - lobe + 1, filter_size, 2 * lobe - 1) -
                  3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // printf("OPTIMIZED: (%i, %i) - Dyy: %f, Dxx: %f, Dxy: %f\n\n", i, j, Dyy, Dxx, Dxy);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    k = (border + 1 + step - 1) / step * step;

    for (int i = k; i < height * step; i += step) {
        ind = (i / step) * width;

        for (int j = 0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dyy = box_integral(iimage, x - border, y - lobe + 1, filter_size, 2 * lobe - 1) -
                  3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }
}

// Only for images >= 128x128
void compute_response_layer_Dyy_top(struct response_layer* layer, struct integral_image* iimage) {
    float Dxx, Dyy, Dxy;
    int x, y, k, k0;
    int r00, r01, c00, c01, r10, r11, c10, c11;
    float Dyy0, Dyy1, A, B, C, D, temp0, temp1;

    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int data_width = iimage->data_width;
    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size / 3;
    int border = (filter_size - 1) / 2;
    float inv_area = 1.f / (filter_size*filter_size);

    int ind = 0;  // oder alternativ (i+1)*j

    float *data = (float*) iimage->data;  // brauch hier keinen cast weil es eig float sein sollte
    int iheight = iimage->height;
    int iwidth = iimage->width;

    // Top Left Corner - Case 1: B of neg part outside
    for (int i = 0; i < lobe / 2 + 1; i += step) {  // Inner B is outside, i.e. 0
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {  // c0 = col - 1 = (y - lobe + 1) -1 < 0
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;  // -1 is already part of r10
            c11 = y + lobe - 1;

            D = data[r11 * data_width + c11];
            Dyy1 = D;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Left Corner - Case 2: B of neg part inside
    // initial value has to be rounded up to next bigger step
    k = (lobe / 2 + 1 + step - 1) / step * step;

    for (int i = k; i < border + 1; i += step) {  // Inner B is inside
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;   // -1 is already part of r10
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];
            Dyy1 = D - B;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Mid - Case 1: B of neg part outside
    k0 = (lobe + step - 1) / step * step;  // initial y value to next multiple of step

    for (int i = 0; i < lobe / 2 + 1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step -lobe+1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            // A and B are outside (= 0), C and D inside
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            // A and B are outside (= 0), C and D inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;   // r10 und r11 kann man vereinfachen
            c10 = y - lobe;
            c11 = y + lobe - 1;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Mid - Case 2: B of neg part inside
    k = (lobe / 2 + 1 + step - 1) / step * step;  // round initial value to next highest step
    for (int i = k; i < border+1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step -lobe+1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Right - Case 1: A of neg part outside
    k0 = (width*step-lobe+1 + step - 1) / step * step;

    for (int i = 0; i < lobe/2+1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Right - Case 2: A of neg part inside
    k = (lobe / 2 + 1 + step - 1) / step * step;

    for (int i = k; i < border+1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dyy = box_integral(iimage, x - border, y - lobe + 1, filter_size, 2*lobe - 1)
                    - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2*lobe - 1);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Rest
    k = (border + 1 + step - 1) / step * step;

    for (int i = k; i < height * step; i += step) {
        ind = (i / step) * width;

        for (int j = 0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dyy = box_integral(iimage, x - border, y - lobe + 1, filter_size, 2 * lobe - 1) -
                  3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }
}

// Only for images >= 128x128
void compute_response_layer_Dyy_top_mid(struct response_layer* layer, struct integral_image* iimage) {
    float Dxx, Dyy, Dxy;
    int x, y, k, k0;
    int r00, r01, c00, c01, r10, r11, c10, c11;
    float Dyy0, Dyy1, A, B, C, D, temp0, temp1;

    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int data_width = iimage->data_width;
    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size/3;
    int border = (filter_size-1)/2;
    float inv_area = 1.f/(filter_size*filter_size);

    int ind = 0;  // oder alternativ (i+1)*j

    float *data = (float*) iimage->data;  // brauch hier keinen cast weil es eig float sein sollte
    int iheight = iimage->height;
    int iwidth = iimage->width;

    /****************
    *   TOP
    *****************/

    // Top Left Corner - Case 1: B of neg part outside
    for (int i = 0; i < lobe / 2 + 1; i += step) {  // Inner B is outside, i.e. 0
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {  // c0 = col - 1 = (y - lobe + 1) -1 < 0
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;  // -1 is already part of r10
            c11 = y + lobe - 1;

            D = data[r11 * data_width + c11];
            Dyy1 = D;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Left Corner - Case 2: B of neg part inside
    // initial value has to be rounded up to next bigger step
    k = (lobe / 2 + 1 + step - 1) / step * step;

    for (int i = k; i < border + 1; i += step) {  // Inner B is inside
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;   // -1 is already part of r10
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];
            Dyy1 = D - B;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Mid - Case 1: B of neg part outside
    k0 = (lobe + step - 1) / step * step;  // initial y value to next multiple of step

    for (int i = 0; i < lobe / 2 + 1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step -lobe+1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            // A and B are outside (= 0), C and D inside
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            // A and B are outside (= 0), C and D inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;   // r10 und r11 kann man vereinfachen
            c10 = y - lobe;
            c11 = y + lobe - 1;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Mid - Case 2: B of neg part inside
    k = (lobe / 2 + 1 + step - 1) / step * step;  // round initial value to next highest step
    for (int i = k; i < border+1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step -lobe+1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Right - Case 1: A of neg part outside
    k0 = (width*step-lobe+1 + step - 1) / step * step;

    for (int i = 0; i < lobe/2+1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Right - Case 2: A of neg part inside
    k = (lobe / 2 + 1 + step - 1) / step * step;

    for (int i = k; i < border+1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dyy = box_integral(iimage, x - border, y - lobe + 1, filter_size, 2*lobe - 1)
                    - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2*lobe - 1);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    /****************
    *   MID
    *****************/

    // Mid Left
    k = (border + 1 + step - 1) / step * step;

    for (int i = k; i < height * step - border; i += step) {
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside (= 0). B, D inside
            r00 = x - border - 1;
            r01 = x + border;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[r01 * data_width + c01];
            Dyy0 = D - B;

            // neg part box filter
            // A, C outside (= 0). B, D inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Mid Mid
    k0 = (lobe + step - 1) / step * step;

    for (int i = k; i < height * step - border; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step - lobe + 1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // All inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // All inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            // There is a very weird floating point arithmetic bug here
            // The original implementation is wrong too, we just try to do it
            // in the same way as them. Maybe fix this in the future
            // Depending on how you add/sub them they give different results

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dyy *= inv_area;
            Dxx *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Mid Right
    k0 = (width*step-lobe+1 + step - 1) / step * step;

    for (int i = k; i < height * step - border; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // B, D outside, A, C inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    /****************
    *   BOTTOM
    *****************/

    // Rest (Bottom)
    k = (height * step - border + step - 1) / step * step;

    for (int i = k; i < height * step; i += step) {
        ind = (i / step) * width;

        for (int j = 0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dyy = box_integral(iimage, x - border, y - lobe + 1, filter_size, 2 * lobe - 1) -
                  3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2 * lobe - 1);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }
}

// Only for images >= 128x128
void compute_response_layer_Dyy(struct response_layer* layer, struct integral_image* iimage) {
    float Dxx, Dyy, Dxy;
    int x, y, k, k0, k1;
    int r00, r01, c00, c01, r10, r11, c10, c11;
    float Dyy0, Dyy1, A, B, C, D, temp0, temp1;

    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int data_width = iimage->data_width;
    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size/3;
    int border = (filter_size-1)/2;
    float inv_area = 1.f/(filter_size*filter_size);

    int ind = 0;  // oder alternativ (i+1)*j

    float *data = (float*) iimage->data;  // brauch hier keinen cast weil es eig float sein sollte
    int iheight = iimage->height;
    int iwidth = iimage->width;

    /****************
    *   TOP
    *****************/

    // Top Left Corner - Case 1: B of neg part outside
    for (int i = 0; i < lobe / 2 + 1; i += step) {  // Inner B is outside, i.e. 0
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {  // c0 = col - 1 = (y - lobe + 1) -1 < 0
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;  // -1 is already part of r10
            c11 = y + lobe - 1;

            D = data[r11 * data_width + c11];
            Dyy1 = D;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Left Corner - Case 2: B of neg part inside
    // initial value has to be rounded up to next bigger step
    k = (lobe / 2 + 1 + step - 1) / step * step;

    for (int i = k; i < border + 1; i += step) {  // Inner B is inside
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;   // -1 is already part of r10
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];
            Dyy1 = D - B;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Mid - Case 1: B of neg part outside
    k0 = (lobe + step - 1) / step * step;  // initial y value to next multiple of step

    for (int i = 0; i < lobe / 2 + 1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step -lobe+1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            // A and B are outside (= 0), C and D inside
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            // A and B are outside (= 0), C and D inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;   // r10 und r11 kann man vereinfachen
            c10 = y - lobe;
            c11 = y + lobe - 1;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Mid - Case 2: B of neg part inside
    k = (lobe / 2 + 1 + step - 1) / step * step;  // round initial value to next highest step
    for (int i = k; i < border+1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step -lobe+1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Right - Case 1: A of neg part outside
    k0 = (width*step-lobe+1 + step - 1) / step * step;

    for (int i = 0; i < lobe/2+1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Top Right - Case 2: A of neg part inside
    k = (lobe / 2 + 1 + step - 1) / step * step;

    for (int i = k; i < border+1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dyy = box_integral(iimage, x - border, y - lobe + 1, filter_size, 2*lobe - 1)
                    - 3 * box_integral(iimage, x - lobe / 2, y - lobe + 1, lobe, 2*lobe - 1);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    /****************
    *   MID
    *****************/

    // Mid Left
    k = (border + 1 + step - 1) / step * step;

    for (int i = k; i < height * step - border; i += step) {
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside (= 0). B, D inside
            r00 = x - border - 1;
            r01 = x + border;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[r01 * data_width + c01];
            Dyy0 = D - B;

            // neg part box filter
            // A, C outside (= 0). B, D inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Mid Mid
    k0 = (lobe + step - 1) / step * step;

    for (int i = k; i < height * step - border; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step - lobe + 1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // All inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            /* Dyy coords
            // whole box filter
            r00 = x - border;
            r01 = x + border;  // -> x-border+filtersize-1 = x-border+(border+border+1)-1
            c00 = y - lobe;
            c01 = y + lobe - 1; // -> y-lobe+1 + 2*lobe -1 -1

            A =
            B =
            C =
            D =

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = (r10, c10)
            B = (r10, c11) = (x - lobe / 2 - 1, y + lobe - 1)
            C = (r11, c10)
            D = (r11, c11)
            */

            // neg part box filter
            // All inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            // There is a very weird floating point arithmetic bug here
            // The original implementation is wrong too, we just try to do it
            // in the same way as them. Maybe fix this in the future
            // Depending on how you add/sub them they give different results

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            // Dyy1 = A - B - C + D;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dyy *= inv_area;
            Dxx *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Mid Right
    k0 = (width*step-lobe+1 + step - 1) / step * step;

    for (int i = k; i < height * step - border; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // B, D outside, A, C inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    /****************
    *   BOTTOM
    *****************/

    // Bottom Left - Case 1: inner D inside
    k = (height * step - border + step - 1) / step * step;

    for (int i = k; i < height * step - lobe / 2; i += step) {
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            r00 = x - border - 1;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C outside, B, D inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Bottom Left - Case 2: inner D outside
    k = (height * step - lobe/2 + step - 1) / step * step;

    for (int i = k; i < height * step; i += step) {
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            r00 = x - border - 1;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C, D outside, B inside
            r10 = x - lobe / 2 - 1;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[(iheight - 1) * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }


    // Bottom Mid - Case 1: inner D inside
    k = (height * step - border + step - 1) / step * step;
    k0 = (lobe + step - 1) / step * step;

    for (int i = k; i < height * step - lobe / 2; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step - lobe + 1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            r00 = x - border - 1;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B, C, D inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Bottom Mid - Case 2: inner D outside
    k = (height * step - lobe/2 + step - 1) / step * step;

    for (int i = k; i < height * step; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step - lobe + 1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            r00 = x - border - 1;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B inside, C, D outside
            r10 = x - lobe / 2 - 1;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Bottom Right - Case 1: inner C inside
    k = (height * step - border + step - 1) / step * step;
    k0 = (width*step-lobe+1 + step - 1) / step * step;

    for (int i = k; i < height * step - lobe / 2; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            r00 = x - border - 1;
            c00 = y - lobe;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }

    // Bottom Right - Case 2: inner C outside
    k = (height * step - lobe / 2 + step - 1) / step * step;
    k0 = (width*step-lobe+1 + step - 1) / step * step;

    for (int i = k; i < height * step; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, C, D outside, A inside
            r10 = x - lobe / 2 - 1;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
            ind += 1;
        }
    }
}

// Only for images >= 128x128
void compute_response_layer_Dyy_laplacian(struct response_layer* layer, struct integral_image* iimage) {
    /*
        simple laplacian fix
    */
    float Dxx, Dyy, Dxy;
    int x, y, k, k0, k1;
    int r00, r01, c00, c01, r10, r11, c10, c11;
    float Dyy0, Dyy1, A, B, C, D, temp0, temp1;

    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int data_width = iimage->data_width;
    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size/3;
    int border = (filter_size-1)/2;
    float inv_area = 1.f/(filter_size*filter_size);

    int ind = 0;  // oder alternativ (i+1)*j

    float *data = (float*) iimage->data;  // brauch hier keinen cast weil es eig float sein sollte
    int iheight = iimage->height;
    int iwidth = iimage->width;

    /****************
    *   TOP
    *****************/

    // Top Left Corner - Case 1: B of neg part outside
    for (int i = 0; i < lobe / 2 + 1; i += step) {  // Inner B is outside, i.e. 0
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {  // c0 = col - 1 = (y - lobe + 1) -1 < 0
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;  // -1 is already part of r10
            c11 = y + lobe - 1;

            D = data[r11 * data_width + c11];
            Dyy1 = D;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    // Top Left Corner - Case 2: B of neg part inside
    // initial value has to be rounded up to next bigger step
    k = (lobe / 2 + 1 + step - 1) / step * step;

    for (int i = k; i < border + 1; i += step) {  // Inner B is inside
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;   // -1 is already part of r10
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];
            Dyy1 = D - B;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    // Top Mid - Case 1: B of neg part outside
    k0 = (lobe + step - 1) / step * step;  // initial y value to next multiple of step

    for (int i = 0; i < lobe / 2 + 1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step -lobe+1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            // A and B are outside (= 0), C and D inside
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            // A and B are outside (= 0), C and D inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;   // r10 und r11 kann man vereinfachen
            c10 = y - lobe;
            c11 = y + lobe - 1;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    // Top Mid - Case 2: B of neg part inside
    k = (lobe / 2 + 1 + step - 1) / step * step;  // round initial value to next highest step
    for (int i = k; i < border+1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step -lobe+1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    // Top Right - Case 1: A of neg part outside
    k0 = (width*step-lobe+1 + step - 1) / step * step;

    for (int i = 0; i < lobe/2+1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    // Top Right - Case 2: A of neg part inside
    k = (lobe / 2 + 1 + step - 1) / step * step;

    for (int i = k; i < border+1; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    /****************
    *   MID
    *****************/

    // Mid Left
    k = (border + 1 + step - 1) / step * step;

    for (int i = k; i < height * step - border; i += step) {
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside (= 0). B, D inside
            r00 = x - border - 1;
            r01 = x + border;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[r01 * data_width + c01];
            Dyy0 = D - B;

            // neg part box filter
            // A, C outside (= 0). B, D inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    // Mid Mid
    k0 = (lobe + step - 1) / step * step;

    for (int i = k; i < height * step - border; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step - lobe + 1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // All inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // All inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            // There is a very weird floating point arithmetic bug here
            // The original implementation is wrong too, we just try to do it
            // in the same way as them. Maybe fix this in the future
            // Depending on how you add/sub them they give different results

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            // Dyy1 = A - B - C + D;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dyy *= inv_area;
            Dxx *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    // Mid Right
    k0 = (width*step-lobe+1 + step - 1) / step * step;

    for (int i = k; i < height * step - border; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // B, D outside, A, C inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    /****************
    *   BOTTOM
    *****************/

    // Bottom Left - Case 1: inner D inside
    k = (height * step - border + step - 1) / step * step;

    for (int i = k; i < height * step - lobe / 2; i += step) {
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            r00 = x - border - 1;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C outside, B, D inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    // Bottom Left - Case 2: inner D outside
    k = (height * step - lobe/2 + step - 1) / step * step;

    for (int i = k; i < height * step; i += step) {
        ind = (i / step) * width;

        for (int j = 0; j < lobe; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            r00 = x - border - 1;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C, D outside, B inside
            r10 = x - lobe / 2 - 1;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[(iheight - 1) * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }


    // Bottom Mid - Case 1: inner D inside
    k = (height * step - border + step - 1) / step * step;
    k0 = (lobe + step - 1) / step * step;

    for (int i = k; i < height * step - lobe / 2; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step - lobe + 1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            r00 = x - border - 1;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B, C, D inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    // Bottom Mid - Case 2: inner D outside
    k = (height * step - lobe/2 + step - 1) / step * step;

    for (int i = k; i < height * step; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step - lobe + 1; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            r00 = x - border - 1;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B inside, C, D outside
            r10 = x - lobe / 2 - 1;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    // Bottom Right - Case 1: inner C inside
    k = (height * step - border + step - 1) / step * step;
    k0 = (width*step-lobe+1 + step - 1) / step * step;

    for (int i = k; i < height * step - lobe / 2; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            r00 = x - border - 1;
            c00 = y - lobe;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    // Bottom Right - Case 2: inner C outside
    k = (height * step - lobe / 2 + step - 1) / step * step;
    k0 = (width*step-lobe+1 + step - 1) / step * step;

    for (int i = k; i < height * step; i += step) {
        ind = (i / step) * width + (k0 / step);

        for (int j = k0; j < width * step; j += step) {
            // Image coordinates
            x = i;
            y = j;

            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, C, D outside, A inside
            r10 = x - lobe / 2 - 1;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }
}

// Only for images >= 128x128
void compute_response_layer_Dyy_laplacian_localityloops(struct response_layer* layer, struct integral_image* iimage) {
    float Dxx, Dyy, Dxy;
    int x, y, k, k0, k1;
    int r00, r01, c00, c01, r10, r11, c10, c11;
    float Dyy0, Dyy1, A, B, C, D, temp0, temp1;

    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int data_width = iimage->data_width;

    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size/3;
    int border = (filter_size-1)/2;
    float inv_area = 1.f/(filter_size*filter_size);

    int ind = 0;  // oder alternativ (i+1)*j

    float *data = (float*) iimage->data;  // brauch hier keinen cast weil es eig float sein sollte
    int iheight = iimage->height;
    int iwidth = iimage->width;

    int j, i;

    /****************
    *   TOP
    *****************/

    for (x = 0; x < lobe / 2 + 1; x += step) {  // Inner B is outside, i.e. 0

        // Top Left Corner - Case 1: B of neg part outside
        for (y = 0; y < lobe; y += step) {  // c0 = col - 1 = (y - lobe + 1) -1 < 0
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;  // -1 is already part of r10
            c11 = y + lobe - 1;

            D = data[r11 * data_width + c11];
            Dyy1 = D;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size)
                - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Mid - Case 1: B of neg part outside
        for (; y < width * step - lobe+1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy
            // whole box filter
            // A and B are outside (= 0), C and D inside
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            // A and B are outside (= 0), C and D inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;   // r10 und r11 kann man vereinfachen
            c10 = y - lobe;
            c11 = y + lobe - 1;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Right - Case 1: A of neg part outside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    for (; x < border + 1; x += step) {  // Inner B is inside

        // Top Left Corner - Case 2: B of neg part inside
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;   // -1 is already part of r10
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];
            Dyy1 = D - B;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Mid - Case 2: B of neg part inside
        for (; y < width * step -lobe+1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Right - Case 2: A of neg part inside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    /****************
    *   MID
    *****************/

    for (; x < height * step - border; x += step) {

        // Mid Left
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside (= 0). B, D inside
            r00 = x - border - 1;
            r01 = x + border;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[r01 * data_width + c01];
            Dyy0 = D - B;

            // neg part box filter
            // A, C outside (= 0). B, D inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Mid Mid
        for (; y < width * step - lobe + 1; y += step) {
            // Compute Dyy  
            // whole box filter
            // All inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // All inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            // There is a very weird floating point arithmetic bug here
            // The original implementation is wrong too, we just try to do it
            // in the same way as them. Maybe fix this in the future
            // Depending on how you add/sub them they give different results

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            // Dyy1 = A - B - C + D;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dyy *= inv_area;
            Dxx *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Mid Right
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // B, D outside, A, C inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    /****************
    *   BOTTOM
    *****************/

    for (; x < height * step - lobe / 2; x += step) {

        // Bottom Left - Case 1: inner D inside
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            r00 = x - border - 1;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C outside, B, D inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Mid - Case 1: inner D inside
        for (; y < width * step - lobe + 1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            r00 = x - border - 1;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B, C, D inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Right - Case 1: inner C inside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            r00 = x - border - 1;
            c00 = y - lobe;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    for (; x < height * step; x += step) {

        // Bottom Left - Case 2: inner D outside
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            r00 = x - border - 1;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C, D outside, B inside
            r10 = x - lobe / 2 - 1;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[(iheight - 1) * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Mid - Case 2: inner D outside
        for (; y < width * step - lobe + 1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            r00 = x - border - 1;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B inside, C, D outside
            r10 = x - lobe / 2 - 1;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Right - Case 2: inner C outside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, C, D outside, A inside
            r10 = x - lobe / 2 - 1;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe + 1, y - border, 2 * lobe - 1, filter_size) -
                  3 * box_integral(iimage, x - lobe + 1, y - lobe / 2, 2 * lobe - 1, lobe);
            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe) +
                  box_integral(iimage, x + 1, y - lobe, lobe, lobe) -
                  box_integral(iimage, x - lobe, y - lobe, lobe, lobe) - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

}

void compute_response_layer_unconditional(struct response_layer *layer, struct integral_image *iimage) {
    float Dxx, Dyy, Dxy;
    int x, y;

    float *response = layer->response;
    bool *laplacian = layer->laplacian;

    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size/3;
    int border = (filter_size-1)/2;
    float inv_area = 1.f/(filter_size*filter_size);
    for (int i = 0, ind = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j, ind++) {
            // Image coordinates
            x = i*step;
            y = j*step;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dyy = box_integral_unconditional(iimage, x - border, y - lobe + 1, filter_size, 2*lobe - 1)
                    - 3 * box_integral_unconditional(iimage, x - lobe / 2, y - lobe + 1, lobe, 2*lobe - 1);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            // response[ind] -= 0.81f * Dxy * Dxy;
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;
            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0);
        }
    }

}

// Only for images >= 128x128
void compute_response_layers_Dyy_laplacian_localityloops_unconditional(struct fasthessian* fh) {
    for (int i = 0; i < fh->total_layers; ++i) {
		compute_response_layer_Dyy_laplacian_localityloops_unconditional(fh->response_map[i], fh->iimage);
	}
}

// Only for images >= 128x128
void compute_response_layer_Dyy_laplacian_localityloops_unconditional(struct response_layer* layer, struct integral_image* iimage) {
    float Dxx, Dyy, Dxy;
    int x, y, k, k0, k1;
    int r00, r01, c00, c01, r10, r11, c10, c11;
    float Dyy0, Dyy1, A, B, C, D, temp0, temp1;

    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size/3;
    int border = (filter_size-1)/2;
    float inv_area = 1.f/(filter_size*filter_size);

    int ind = 0;  // oder alternativ (i+1)*j

    float *data = iimage->data;
    int iheight = iimage->height;
    int iwidth = iimage->width;
    int data_width = iimage->data_width;

    int j, i;

    /****************
    *   TOP
    *****************/

    for (x = 0; x < lobe / 2 + 1; x += step) {  // Inner B is outside, i.e. 0

        // Top Left Corner - Case 1: B of neg part outside
        for (y = 0; y < lobe; y += step) {  // c0 = col - 1 = (y - lobe + 1) -1 < 0
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;  // -1 is already part of r10
            c11 = y + lobe - 1;

            D = data[r11 * data_width + c11];
            Dyy1 = D;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Mid - Case 1: B of neg part outside
        for (; y < width * step - lobe+1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy
            // whole box filter
            // A and B are outside (= 0), C and D inside
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            // A and B are outside (= 0), C and D inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;   // r10 und r11 kann man vereinfachen
            c10 = y - lobe;
            c11 = y + lobe - 1;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Right - Case 1: A of neg part outside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    for (; x < border + 1; x += step) {  // Inner B is inside

        // Top Left Corner - Case 2: B of neg part inside
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;   // -1 is already part of r10
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];
            Dyy1 = D - B;

            Dyy = Dyy0 - 3 * Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Mid - Case 2: B of neg part inside
        for (; y < width * step -lobe+1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Right - Case 2: A of neg part inside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    /****************
    *   MID
    *****************/

    for (; x < height * step - border; x += step) {

        // Mid Left
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside (= 0). B, D inside
            r00 = x - border - 1;
            r01 = x + border;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[r01 * data_width + c01];
            Dyy0 = D - B;

            // neg part box filter
            // A, C outside (= 0). B, D inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Mid Mid
        for (; y < width * step - lobe + 1; y += step) {
            // Compute Dyy  
            // whole box filter
            // All inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // All inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            // There is a very weird floating point arithmetic bug here
            // The original implementation is wrong too, we just try to do it
            // in the same way as them. Maybe fix this in the future
            // Depending on how you add/sub them they give different results

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            // Dyy1 = A - B - C + D;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dyy *= inv_area;
            Dxx *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Mid Right
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // B, D outside, A, C inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    /****************
    *   BOTTOM
    *****************/

    for (; x < height * step - lobe / 2; x += step) {

        // Bottom Left - Case 1: inner D inside
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            r00 = x - border - 1;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C outside, B, D inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Mid - Case 1: inner D inside
        for (; y < width * step - lobe + 1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            r00 = x - border - 1;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B, C, D inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Right - Case 1: inner C inside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            r00 = x - border - 1;
            c00 = y - lobe;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    for (; x < height * step; x += step) {

        // Bottom Left - Case 2: inner D outside
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            r00 = x - border - 1;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C, D outside, B inside
            r10 = x - lobe / 2 - 1;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[(iheight - 1) * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Mid - Case 2: inner D outside
        for (; y < width * step - lobe + 1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            r00 = x - border - 1;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B inside, C, D outside
            r10 = x - lobe / 2 - 1;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Right - Case 2: inner C outside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, C, D outside, A inside
            r10 = x - lobe / 2 - 1;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }
}

// Only for images >= 128x128
void compute_response_layers_Dyy_laplacian_locality_uncond_opt(struct fasthessian* fh) {
    for (int i = 0; i < fh->total_layers; ++i) {
		compute_response_layer_Dyy_laplacian_locality_uncond_opt(fh->response_map[i], fh->iimage);
	}
}

// Only for images >= 128x128
void compute_response_layer_Dyy_laplacian_locality_uncond_opt(struct response_layer* layer, struct integral_image* iimage) {
    /*
        Optimized flop count and removed unnecessary computations
        and changed i to x and j to y
    */
    float Dxx, Dyy, Dxy;
    int x, y, k, k0, k1;
    int r00, r01, c00, c01, r10, r11, c10, c11;
    float Dyy0, Dyy1, A, B, C, D, temp0, temp1;

    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size/3;
    int border = (filter_size-1)/2;
    float inv_area = 1.f/(filter_size*filter_size);

    int ind = 0;  // oder alternativ (i+1)*j

    float *data = iimage->data;
    int iheight = iimage->height;
    int iwidth = iimage->width;
    int data_width = iimage->data_width;

    int j, i, t0, t1, t2, t3, t4, t5, t6, t7, t8;

    /****************
    *   TOP
    *****************/

    // Image coordinates x and y (row and column)

    for (x = 0; x < lobe / 2 + 1; x += step) {  // Inner B is outside, i.e. 0

        // precompute
        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Top Left Corner - Case 1: B of neg part outside
        for (y = 0; y < lobe; y += step) {  // c0 = col - 1 = (y - lobe + 1) -1 < 0
            // Compute Dyy  
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = x - lobe / 2 + lobe - 1;
            c11 = y + lobe - 1;

            D = data[r11 * data_width + c11];
            Dyy1 = D;

            Dyy = Dyy0 - 3 * Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Mid - Case 1: B of neg part outside
        for (; y < width * step - lobe+1; y += step) {
            // Compute Dyy
            // whole box filter
            // A and B are outside (= 0), C and D inside
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            // A and B are outside (= 0), C and D inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;   // r10 und r11 kann man vereinfachen
            c10 = y - lobe;
            c11 = y + lobe - 1;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Right - Case 1: A of neg part outside
        for (; y < width * step; y += step) {
            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;

            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    for (; x < border + 1; x += step) {  // Inner B is inside

        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Top Left Corner - Case 2: B of neg part inside
        for (y = 0; y < lobe; y += step) {
            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * data_width + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;   // -1 is already part of r10
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];
            Dyy1 = D - B;

            Dyy = Dyy0 - 3 * Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Mid - Case 2: B of neg part inside
        for (; y < width * step -lobe+1; y += step) {
            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Right - Case 2: A of neg part inside
        for (; y < width * step; y += step) {
            // whole filter
            r01 = x + border;
            c00 = y - lobe;

            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    /****************
    *   MID
    *****************/

    for (; x < height * step - border; x += step) {

        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Mid Left
        for (y = 0; y < lobe; y += step) {
            // Compute Dyy  
            // whole box filter
            // A, C outside (= 0). B, D inside
            r00 = x - border - 1;
            r01 = x + border;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[r01 * data_width + c01];
            Dyy0 = D - B;

            // neg part box filter
            // A, C outside (= 0). B, D inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Mid Mid
        for (; y < width * step - lobe + 1; y += step) {
            // Compute Dyy  
            // whole box filter
            // All inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // All inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dyy *= inv_area;
            Dxx *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Mid Right
        for (; y < width * step; y += step) {
            // Compute Dyy  
            // whole box filter
            // B, D outside, A, C inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[r01 * data_width + c00];
            D = data[r01 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    /****************
    *   BOTTOM
    *****************/

    for (; x < height * step - lobe / 2; x += step) {

        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Bottom Left - Case 1: inner D inside
        for (y = 0; y < lobe; y += step) {
            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            r00 = x - border - 1;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C outside, B, D inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[r11 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Mid - Case 1: inner D inside
        for (; y < width * step - lobe + 1; y += step) {
            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            r00 = x - border - 1;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B, C, D inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Right - Case 1: inner C inside
        for (; y < width * step; y += step) {
            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            r00 = x - border - 1;
            c00 = y - lobe;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            r10 = x - lobe / 2 - 1;
            r11 = x + lobe/2;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[r11 * data_width + c10];
            D = data[r11 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    for (; x < height * step; x += step) {

        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Bottom Left - Case 2: inner D outside
        for (y = 0; y < lobe; y += step) {
            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            r00 = x - border - 1;
            c01 = y + lobe - 1;

            B = data[r00 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C, D outside, B inside
            r10 = x - lobe / 2 - 1;
            c11 = y + lobe - 1;

            B = data[r10 * data_width + c11];
            D = data[(iheight - 1) * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Mid - Case 2: inner D outside
        for (; y < width * step - lobe + 1; y += step) {
            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            r00 = x - border - 1;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B inside, C, D outside
            r10 = x - lobe / 2 - 1;
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + c11];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Right - Case 2: inner C outside
        for (; y < width * step; y += step) {
            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            r00 = x - border - 1;
            r01 = x + border;
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00 * data_width + c00];
            B = data[r00 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, C, D outside, A inside
            r10 = x - lobe / 2 - 1;
            c10 = y - lobe;

            A = data[r10 * data_width + c10];
            B = data[r10 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }
}

// Only for images >= 128x128
void compute_response_layers_Dyy_laplacian_locality_uncond_opt_flops(struct fasthessian* fh) {
    for (int i = 0; i < fh->total_layers; ++i) {
		compute_response_layer_Dyy_laplacian_locality_uncond_opt_flops(fh->response_map[i], fh->iimage);
	}
}

// Only for images >= 128x128
void compute_response_layer_Dyy_laplacian_locality_uncond_opt_flops(struct response_layer* layer, struct integral_image* iimage) {
    /*
        Flops optimized even further
    */
    float Dxx, Dyy, Dxy;
    int x, y, k, k0, k1;
    int r00, r01, c00, c01, r10, r11, c10, c11;
    float Dyy0, Dyy1, A, B, C, D, temp0, temp1;

    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size/3;
    int border = (filter_size-1)/2;
    float inv_area = 1.f/(filter_size*filter_size);
    // float inv_area_squared = 1.f/(filter_size*filter_size) * 1.f/(filter_size*filter_size);

    int ind = 0;  // oder alternativ (i+1)*j

    float *data = iimage->data;
    int iheight = iimage->height;
    int iwidth = iimage->width;
    int data_width = iimage->data_width;

    int j, i, t0, t1, t2, t3, t4, t5, t6, t7, t8;
    int r01_T1, r11_T1, r01_T2, r10_T2, r11_T2;
    int r00_M1, r01_M1, r10_M1, r11_M1;
    int r00_B1, r10_B1, r11_B1, r00_B2, r10_B2;

    /****************
    *   TOP
    *****************/

    for (x = 0; x < lobe / 2 + 1; x += step) {  // Inner B is outside, i.e. 0

        // Box integral precompute
        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;  // TODO: get out of for loop
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Top 1 precompute
        r01_T1 = x + border;
        r11_T1 = x + lobe / 2;

        // Top Left Corner - Case 1: B of neg part outside
        for (y = 0; y < lobe; y += step) {  // c0 = col - 1 = (y - lobe + 1) -1 < 0
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            c01 = y + lobe - 1;  // TODO: compute only once get lobe-1 out of loop

            Dyy0 = data[r01_T1 * data_width + c01];

            // neg part box filter
            c11 = y + lobe - 1;

            D = data[r11_T1 * data_width + c11];
            Dyy1 = D;

            Dyy = Dyy0 - 3 * Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // response[ind] = inv_area * inv_area * Dxx * Dyy - 0.81f inv_area * inv_area* Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Mid - Case 1: B of neg part outside
        for (; y < width * step - lobe+1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy
            // whole box filter
            // A and B are outside (= 0), C and D inside
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01_T1 * data_width + c00];
            D = data[r01_T1 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            // A and B are outside (= 0), C and D inside
            c10 = y - lobe;
            c11 = y + lobe - 1;

            C = data[r11_T1 * data_width + c10];
            D = data[r11_T1 * data_width + c11];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Right - Case 1: A of neg part outside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // whole filter
            c00 = y - lobe;

            C = data[r01_T1 * data_width + c00];
            D = data[r01_T1 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            c10 = y - lobe;

            C = data[r11_T1 * data_width + c10];
            D = data[r11_T1 * data_width + iwidth-1];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    for (; x < border + 1; x += step) {  // Inner B is inside

        // Box integral precompute
        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Top 2 precompute
        r01_T2 = x + border;
        r10_T2 = x - lobe / 2 - 1;
        r11_T2 = x + lobe / 2;

        // Top Left Corner - Case 2: B of neg part inside
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy
            // whole box filter
            c01 = y + lobe - 1;

            Dyy0 = data[r01_T2 * data_width + c01];

            // neg part box filter
            c11 = y + lobe - 1;

            B = data[r10_T2 * data_width + c11];
            D = data[r11_T2 * data_width + c11];
            Dyy1 = D - B;

            Dyy = Dyy0 - 3 * Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Mid - Case 2: B of neg part inside
        for (; y < width * step -lobe+1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy
            // whole box filter
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01_T2 * data_width + c00];
            D = data[r01_T2 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10_T2 * data_width + c10];
            B = data[r10_T2 * data_width + c11];
            C = data[r11_T2 * data_width + c10];
            D = data[r11_T2 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Right - Case 2: A of neg part inside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // whole filter
            c00 = y - lobe;

            C = data[r01_T2 * data_width + c00];
            D = data[r01_T2 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            c10 = y - lobe;

            A = data[r10_T2 * data_width + c10];
            B = data[r10_T2 * data_width + iwidth-1];
            C = data[r11_T2 * data_width + c10];
            D = data[r11_T2 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    /****************
    *   MID
    *****************/

    for (; x < height * step - border; x += step) {

        // Box integral precompute
        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Mid 1 precompute
        r00_M1 = x - border - 1;
        r01_M1 = x + border;
        r10_M1 = x - lobe / 2 - 1;
        r11_M1 = x + lobe/2;

        // Mid Left
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside (= 0). B, D inside
            c01 = y + lobe - 1;

            B = data[r00_M1 * data_width + c01];
            D = data[r01_M1 * data_width + c01];
            Dyy0 = D - B;

            // neg part box filter
            // A, C outside (= 0). B, D inside
            c11 = y + lobe - 1;

            B = data[r10_M1 * data_width + c11];
            D = data[r11_M1 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Mid Mid
        for (; y < width * step - lobe + 1; y += step) {
            // Compute Dyy  
            // whole box filter
            // All inside
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00_M1 * data_width + c00];
            B = data[r00_M1 * data_width + c01];
            C = data[r01_M1 * data_width + c00];
            D = data[r01_M1 * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // All inside
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10_M1 * data_width + c10];
            B = data[r10_M1 * data_width + c11];
            C = data[r11_M1 * data_width + c10];
            D = data[r11_M1 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dyy *= inv_area;
            Dxx *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Mid Right
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // B, D outside, A, C inside
            c00 = y - lobe;

            A = data[r00_M1 * data_width + c00];
            B = data[r00_M1 * data_width + iwidth-1];
            C = data[r01_M1 * data_width + c00];
            D = data[r01_M1 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            c10 = y - lobe;

            A = data[r10_M1 * data_width + c10];
            B = data[r10_M1 * data_width + iwidth-1];
            C = data[r11_M1 * data_width + c10];
            D = data[r11_M1 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    /****************
    *   BOTTOM
    *****************/

    for (; x < height * step - lobe / 2; x += step) {

        // Box integral precompute
        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Bottom 1 precompute
        r00_B1 = x - border - 1;
        r10_B1 = x - lobe / 2 - 1;
        r11_B1 = x + lobe/2;

        // Bottom Left - Case 1: inner D inside
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            c01 = y + lobe - 1;

            B = data[r00_B1 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C outside, B, D inside
            c11 = y + lobe - 1;

            B = data[r10_B1 * data_width + c11];
            D = data[r11_B1 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Mid - Case 1: inner D inside
        for (; y < width * step - lobe + 1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00_B1 * data_width + c00];
            B = data[r00_B1 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B, C, D inside
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10_B1 * data_width + c10];
            B = data[r10_B1 * data_width + c11];
            C = data[r11_B1 * data_width + c10];
            D = data[r11_B1 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Right - Case 1: inner C inside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            c00 = y - lobe;

            A = data[r00_B1 * data_width + c00];
            B = data[r00_B1 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            c10 = y - lobe;

            A = data[r10_B1 * data_width + c10];
            B = data[r10_B1 * data_width + iwidth-1];
            C = data[r11_B1 * data_width + c10];
            D = data[r11_B1 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    for (; x < height * step; x += step) {

        // Box integral precompute
        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Bottom 2 precompute
        r00_B2 = x - border - 1;
        r10_B2 = x - lobe / 2 - 1;

        // Bottom Left - Case 2: inner D outside
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            c01 = y + lobe - 1;

            B = data[r00_B2 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C, D outside, B inside
            c11 = y + lobe - 1;

            B = data[r10_B2 * data_width + c11];
            D = data[(iheight - 1) * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Mid - Case 2: inner D outside
        for (; y < width * step - lobe + 1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00_B2 * data_width + c00];
            B = data[r00_B2 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B inside, C, D outside
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10_B2 * data_width + c10];
            B = data[r10_B2 * data_width + c11];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Right - Case 2: inner C outside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00_B2 * data_width + c00];
            B = data[r00_B2 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, C, D outside, A inside
            c10 = y - lobe;

            A = data[r10_B2 * data_width + c10];
            B = data[r10_B2 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;
            Dxy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy - 0.81f * Dxy * Dxy;

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }
}

// Only for images >= 128x128
void compute_response_layers_Dyy_laplacian_locality_uncond_opt_flops_invsqr(struct fasthessian* fh) {
    for (int i = 0; i < fh->total_layers; ++i) {
		compute_response_layer_Dyy_laplacian_locality_uncond_opt_flops_invsqr(fh->response_map[i], fh->iimage);
	}
}

// Only for images >= 128x128
void compute_response_layer_Dyy_laplacian_locality_uncond_opt_flops_invsqr(struct response_layer* layer, struct integral_image* iimage) {
    /*
        Flops noch weiter verbessert indem
    */
    float Dxx, Dyy, Dxy;
    int x, y, k, k0, k1;
    int r00, r01, c00, c01, r10, r11, c10, c11;
    float Dyy0, Dyy1, A, B, C, D, temp0, temp1;

    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size/3;
    int border = (filter_size-1)/2;
    // float inv_area = 1.f/(filter_size*filter_size);

    // TODO: write as 1/(filter_size*filter_size*filter_size*filter_size)
    float inv_area_squared = 1.f/(filter_size*filter_size) * 1.f/(filter_size*filter_size);

    int ind = 0;  // oder alternativ (i+1)*j

    float *data = iimage->data;
    int iheight = iimage->height;
    int iwidth = iimage->width;
    int data_width = iimage->data_width;

    int j, i, t0, t1, t2, t3, t4, t5, t6, t7, t8;
    int r01_T1, r11_T1, r01_T2, r10_T2, r11_T2;
    int r00_M1, r01_M1, r10_M1, r11_M1;
    int r00_B1, r10_B1, r11_B1, r00_B2, r10_B2;

    /****************
    *   TOP
    *****************/

    for (x = 0; x < lobe / 2 + 1; x += step) {  // Inner B is outside, i.e. 0

        // Box integral precompute
        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Top 1 precompute
        r01_T1 = x + border;
        r11_T1 = x + lobe / 2;

        // Top Left Corner - Case 1: B of neg part outside
        for (y = 0; y < lobe; y += step) {  // c0 = col - 1 = (y - lobe + 1) -1 < 0
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            c01 = y + lobe - 1;

            Dyy0 = data[r01_T1 * data_width + c01];

            // neg part box filter
            c11 = y + lobe - 1;

            D = data[r11_T1 * data_width + c11];
            Dyy1 = D;

            Dyy = Dyy0 - 3 * Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Mid - Case 1: B of neg part outside
        for (; y < width * step - lobe+1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy
            // whole box filter
            // A and B are outside (= 0), C and D inside
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01_T1 * data_width + c00];
            D = data[r01_T1 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            // A and B are outside (= 0), C and D inside
            c10 = y - lobe;
            c11 = y + lobe - 1;

            C = data[r11_T1 * data_width + c10];
            D = data[r11_T1 * data_width + c11];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Right - Case 1: A of neg part outside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // whole filter
            c00 = y - lobe;

            C = data[r01_T1 * data_width + c00];
            D = data[r01_T1 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            c10 = y - lobe;

            C = data[r11_T1 * data_width + c10];
            D = data[r11_T1 * data_width + iwidth-1];

            Dyy1 = D - C;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    for (; x < border + 1; x += step) {  // Inner B is inside

        // Box integral precompute
        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Top 2 precompute
        r01_T2 = x + border;
        r10_T2 = x - lobe / 2 - 1;
        r11_T2 = x + lobe / 2;

        // Top Left Corner - Case 2: B of neg part inside
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy
            // whole box filter
            c01 = y + lobe - 1;

            Dyy0 = data[r01_T2 * data_width + c01];

            // neg part box filter
            c11 = y + lobe - 1;

            B = data[r10_T2 * data_width + c11];
            D = data[r11_T2 * data_width + c11];
            Dyy1 = D - B;

            Dyy = Dyy0 - 3 * Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Mid - Case 2: B of neg part inside
        for (; y < width * step -lobe+1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy
            // whole box filter
            c00 = y - lobe;
            c01 = y + lobe - 1;

            C = data[r01_T2 * data_width + c00];
            D = data[r01_T2 * data_width + c01];
            Dyy0 = D - C;

            // neg part box filter
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10_T2 * data_width + c10];
            B = data[r10_T2 * data_width + c11];
            C = data[r11_T2 * data_width + c10];
            D = data[r11_T2 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Top Right - Case 2: A of neg part inside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // whole filter
            c00 = y - lobe;

            C = data[r01_T2 * data_width + c00];
            D = data[r01_T2 * data_width + iwidth-1];
            Dyy0 = D - C;

            // neg part box filter
            c10 = y - lobe;

            A = data[r10_T2 * data_width + c10];
            B = data[r10_T2 * data_width + iwidth-1];
            C = data[r11_T2 * data_width + c10];
            D = data[r11_T2 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    /****************
    *   MID
    *****************/

    for (; x < height * step - border; x += step) {

        // Box integral precompute
        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Mid 1 precompute
        r00_M1 = x - border - 1;
        r01_M1 = x + border;
        r10_M1 = x - lobe / 2 - 1;
        r11_M1 = x + lobe/2;

        // Mid Left
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside (= 0). B, D inside
            c01 = y + lobe - 1;

            B = data[r00_M1 * data_width + c01];
            D = data[r01_M1 * data_width + c01];
            Dyy0 = D - B;

            // neg part box filter
            // A, C outside (= 0). B, D inside
            c11 = y + lobe - 1;

            B = data[r10_M1 * data_width + c11];
            D = data[r11_M1 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Mid Mid
        for (; y < width * step - lobe + 1; y += step) {
            // Compute Dyy  
            // whole box filter
            // All inside
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00_M1 * data_width + c00];
            B = data[r00_M1 * data_width + c01];
            C = data[r01_M1 * data_width + c00];
            D = data[r01_M1 * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // All inside
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10_M1 * data_width + c10];
            B = data[r10_M1 * data_width + c11];
            C = data[r11_M1 * data_width + c10];
            D = data[r11_M1 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Mid Right
        for (; y < width * step; y += step) {
            // Compute Dyy  
            // whole box filter
            // B, D outside, A, C inside
            c00 = y - lobe;

            A = data[r00_M1 * data_width + c00];
            B = data[r00_M1 * data_width + iwidth-1];
            C = data[r01_M1 * data_width + c00];
            D = data[r01_M1 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            c10 = y - lobe;

            A = data[r10_M1 * data_width + c10];
            B = data[r10_M1 * data_width + iwidth-1];
            C = data[r11_M1 * data_width + c10];
            D = data[r11_M1 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

    }

    /****************
    *   BOTTOM
    *****************/

    for (; x < height * step - lobe / 2; x += step) {

        // Box integral precompute
        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Bottom 1 precompute
        r00_B1 = x - border - 1;
        r10_B1 = x - lobe / 2 - 1;
        r11_B1 = x + lobe/2;

        // Bottom Left - Case 1: inner D inside
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            c01 = y + lobe - 1;

            B = data[r00_B1 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C outside, B, D inside
            c11 = y + lobe - 1;

            B = data[r10_B1 * data_width + c11];
            D = data[r11_B1 * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Mid - Case 1: inner D inside
        for (; y < width * step - lobe + 1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00_B1 * data_width + c00];
            B = data[r00_B1 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B, C, D inside
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10_B1 * data_width + c10];
            B = data[r10_B1 * data_width + c11];
            C = data[r11_B1 * data_width + c10];
            D = data[r11_B1 * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Right - Case 1: inner C inside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            c00 = y - lobe;

            A = data[r00_B1 * data_width + c00];
            B = data[r00_B1 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, D outside, A, C inside
            c10 = y - lobe;

            A = data[r10_B1 * data_width + c10];
            B = data[r10_B1 * data_width + iwidth-1];
            C = data[r11_B1 * data_width + c10];
            D = data[r11_B1 * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }

    for (; x < height * step; x += step) {

        // Box integral precompute
        t0 = x - lobe;
        t1 = x + lobe - 1;
        t2 = lobe/2 + 1;
        t3 = x - 1;
        t5 = x - lobe - 1;
        t7 = x + lobe;

        // Bottom 2 precompute
        r00_B2 = x - border - 1;
        r10_B2 = x - lobe / 2 - 1;

        // Bottom Left - Case 2: inner D outside
        for (y = 0; y < lobe; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, C outside, B inside, D below
            c01 = y + lobe - 1;

            B = data[r00_B2 * data_width + c01];
            D = data[(iheight-1) * data_width + c01];

            Dyy0 = D - B;

            // neg part box filter
            // A, C, D outside, B inside
            c11 = y + lobe - 1;

            B = data[r10_B2 * data_width + c11];
            D = data[(iheight - 1) * data_width + c11];

            Dyy1 = D - B;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Mid - Case 2: inner D outside
        for (; y < width * step - lobe + 1; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // A, B inside, C, D outside
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00_B2 * data_width + c00];
            B = data[r00_B2 * data_width + c01];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + c01];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // A, B inside, C, D outside
            c10 = y - lobe;
            c11 = y + lobe - 1;

            A = data[r10_B2 * data_width + c10];
            B = data[r10_B2 * data_width + c11];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + c11];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }

        // Bottom Right - Case 2: inner C outside
        for (; y < width * step; y += step) {
            // Image coordinates
            // x = i;
            // y = j;

            // Compute Dyy  
            // whole box filter
            // B, C, D outside A inside
            c00 = y - lobe;
            c01 = y + lobe - 1;

            A = data[r00_B2 * data_width + c00];
            B = data[r00_B2 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c00];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy0 = temp0 + temp1;

            // neg part box filter
            // B, C, D outside, A inside
            c10 = y - lobe;

            A = data[r10_B2 * data_width + c10];
            B = data[r10_B2 * data_width + iwidth-1];
            C = data[(iheight-1) * data_width + c10];
            D = data[(iheight-1) * data_width + iwidth-1];

            temp0 = A - C;
            temp1 = D - B;
            Dyy1 = temp0 + temp1;

            Dyy = Dyy0 - 3*Dyy1;

            // box_integral_unconditional precompute
            t4 = y - 1;
            t6 = y - lobe - 1;
            t8 = y + lobe;

            // Compute Dxx, Dxy
            Dxx = box_integral_unconditional_opt(iimage, t0, y-border-1, t1, y+border)
                - 3 * box_integral_unconditional_opt(iimage, t0, y-t2, t1, t8-t2);
            Dxy = box_integral_unconditional_opt(iimage, t5, y, t3, t8)
                    + box_integral_unconditional_opt(iimage, x, t6, t7, t4)
                    - box_integral_unconditional_opt(iimage, t5, t6, t3, t4)
                    - box_integral_unconditional_opt(iimage, x, y, t7, t8);

            // Calculate Determinant
            response[ind] = inv_area_squared * (Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian
            laplacian[ind] = Dxx + Dyy >= 0;
            ind += 1;
        }
    }
}

void compute_response_layer_precompute(struct response_layer* layer, struct integral_image* iimage) {
    /*
    optimizations:
        - simplifying normalization and precomputing inv_area_square => -2 flops per loop iteration
        - precomputation of indices
    */

    float Dxx, Dyy, Dxy;
    int x, y;

    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size/3;
    int border = (filter_size-1)/2;
    float inv_area = 1.f/(filter_size*filter_size);

    float inv_area_square = inv_area*inv_area;

    int lobe_div_2 = lobe / 2;
    int lobe_sub_1 = lobe - 1;
    int lobe_mul_2_sub_1 = 2*lobe - 1;

    for (int i = 0, ind = 0; i < height; ++i) {
        x = i*step;
        for (int j = 0; j < width; ++j, ind++) {
            y = j*step;

            // // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx = box_integral(iimage, x - lobe_sub_1, y - border, lobe_mul_2_sub_1, filter_size)
                    - 3 * box_integral(iimage, x - lobe_sub_1, y - lobe_div_2, lobe_mul_2_sub_1, lobe);
            Dyy = box_integral(iimage, x - border, y - lobe_sub_1, filter_size, lobe_mul_2_sub_1)
                    - 3 * box_integral(iimage, x - lobe_div_2, y - lobe_sub_1, lobe, lobe_mul_2_sub_1);

            Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            // Calculate Determinant & normalize
            response[ind] = inv_area_square*(Dxx * Dyy - 0.81f * Dxy * Dxy);

            // Calculate Laplacian (scaling does not matter as rescaling both summands keeps the sign)
            laplacian[ind] = (Dxx + Dyy >= 0 ? true : false);
        }
    }

}

void compute_response_layers_precompute(struct fasthessian* fh) {

    for (int i = 0; i < fh->total_layers; ++i) {
		compute_response_layer_precompute(fh->response_map[i], fh->iimage);
	}

}

void compute_response_layer_blocking(struct response_layer* layer, struct integral_image* iimage) {
    float* response = layer->response;
    bool* laplacian = layer->laplacian;

    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size/3;
    int border = (filter_size-1)/2;
    float inv_area = 1.f/(filter_size*filter_size);

    float inv_area_square = inv_area*inv_area;

    int lobe_div_2 = lobe / 2;
    int lobe_sub_1 = lobe - 1;
    int lobe_mul_2_sub_1 = 2*lobe - 1;

    int i;
    for (i = 0; i < height-(1-1); i+=1) {

        int i0 = (i+0);

        int x0 = i0*step;
        int j=0;
        for (; j < width-(2-1); j+=2) {
            int j0 = (j+0);
            int j1= (j+1);

            int y0 = j0*step;
            int y1 = (j+1)*step;

            float Dxx0_0 = box_integral(iimage, x0 - lobe_sub_1, y0 - border, lobe_mul_2_sub_1, filter_size)
                    - 3 * box_integral(iimage, x0 - lobe_sub_1, y0 - lobe_div_2, lobe_mul_2_sub_1, lobe);
            float Dxx0_1 = box_integral(iimage, x0 - lobe_sub_1, y1 - border, lobe_mul_2_sub_1, filter_size)
                    - 3 * box_integral(iimage, x0 - lobe_sub_1, y1 - lobe_div_2, lobe_mul_2_sub_1, lobe);

            float Dyy0_0 = box_integral(iimage, x0 - border, y0 - lobe_sub_1, filter_size, lobe_mul_2_sub_1)
                    - 3 * box_integral(iimage, x0 - lobe_div_2, y0 - lobe_sub_1, lobe, lobe_mul_2_sub_1);
            float Dyy0_1 = box_integral(iimage, x0 - border, y1 - lobe_sub_1, filter_size, lobe_mul_2_sub_1)
                    - 3 * box_integral(iimage, x0 - lobe_div_2, y1 - lobe_sub_1, lobe, lobe_mul_2_sub_1);

            float Dxy0_0 = box_integral(iimage, x0 - lobe, y0 + 1, lobe, lobe)
                    + box_integral(iimage, x0 + 1, y0 - lobe, lobe, lobe)
                    - box_integral(iimage, x0 - lobe, y0 - lobe, lobe, lobe)
                    - box_integral(iimage, x0 + 1, y0 + 1, lobe, lobe);
            float Dxy0_1 = box_integral(iimage, x0 - lobe, y1 + 1, lobe, lobe)
                    + box_integral(iimage, x0 + 1, y1 - lobe, lobe, lobe)
                    - box_integral(iimage, x0 - lobe, y1 - lobe, lobe, lobe)
                    - box_integral(iimage, x0 + 1, y1 + 1, lobe, lobe);

            response[(i0*width) + j0] = inv_area_square*(Dxx0_0 * Dyy0_0 - 0.81f * Dxy0_0 * Dxy0_0);
            response[(i0*width) + j1] = inv_area_square*(Dxx0_1 * Dyy0_1 - 0.81f * Dxy0_1 * Dxy0_1);

            laplacian[(i0*width) + j0] = (Dxx0_0 + Dyy0_0 >= 0 ? true : false);
            laplacian[(i0*width) + j1] = (Dxx0_1 + Dyy0_1 >= 0 ? true : false);
        }

        for (; j < width; ++j) {
            int y = j*step;

            float Dxx = box_integral(iimage, x0 - lobe_sub_1, y - border, lobe_mul_2_sub_1, filter_size)
                    - 3 * box_integral(iimage, x0 - lobe_sub_1, y - lobe_div_2, lobe_mul_2_sub_1, lobe);
            float Dyy = box_integral(iimage, x0 - border, y - lobe_sub_1, filter_size, lobe_mul_2_sub_1)
                    - 3 * box_integral(iimage, x0 - lobe_div_2, y - lobe_sub_1, lobe, lobe_mul_2_sub_1);

            float Dxy = box_integral(iimage, x0 - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x0 + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x0 - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x0 + 1, y + 1, lobe, lobe);

            response[(i0*width) + j] = inv_area_square*(Dxx * Dyy - 0.81f * Dxy * Dxy);

            laplacian[(i0*width) + j] = (Dxx + Dyy >= 0 ? true : false);
        }
    }

    for (; i < height; ++i) {
        int x = i*step;
        for (int j = 0; j < width; ++j) {
            int y = j*step;

            float Dxx = box_integral(iimage, x - lobe_sub_1, y - border, lobe_mul_2_sub_1, filter_size)
                    - 3 * box_integral(iimage, x - lobe_sub_1, y - lobe_div_2, lobe_mul_2_sub_1, lobe);
            float Dyy = box_integral(iimage, x - border, y - lobe_sub_1, filter_size, lobe_mul_2_sub_1)
                    - 3 * box_integral(iimage, x - lobe_div_2, y - lobe_sub_1, lobe, lobe_mul_2_sub_1);

            float Dxy = box_integral(iimage, x - lobe, y + 1, lobe, lobe)
                    + box_integral(iimage, x + 1, y - lobe, lobe, lobe)
                    - box_integral(iimage, x - lobe, y - lobe, lobe, lobe)
                    - box_integral(iimage, x + 1, y + 1, lobe, lobe);

            response[(i*width) + j] = inv_area_square*(Dxx * Dyy - 0.81f * Dxy * Dxy);

            laplacian[(i*width) + j] = (Dxx + Dyy >= 0 ? true : false);
        }
    }

}

void compute_response_layers_blocking(struct fasthessian* fh) {

    for (int i = 0; i < fh->total_layers; ++i) {
		compute_response_layer_blocking(fh->response_map[i], fh->iimage);
	}

}

void compute_response_layers_at_once(struct fasthessian* fh) {
    /*
    TODO: (valentin) add index precomputation as in compute_response_layer_precompute

    optimizations:
        - precompute inv_area_squared as in compute_response_layer_precompute
        - computes all 8 response layers at once, gives same results as base implementation


    Results:
        - slower than before :(
        - valgrind reports no improvement for l1 misses, i.e. locality is not improved as expected
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
    int step = fh->step;
    struct integral_image* iimage = fh->iimage;

    // step == 1: 0,1,2,3
    // step == 2: 4,5
    // step == 4: 6,7

    struct response_layer *layer0 = fh->response_map[0];  // step0 == step
    struct response_layer *layer1 = fh->response_map[1];  // step1 == step
    struct response_layer *layer2 = fh->response_map[2];  // step2 == step
    struct response_layer *layer3 = fh->response_map[3];  // step3 == step

    struct response_layer *layer4 = fh->response_map[4];  // step4 == 2*step
    struct response_layer *layer5 = fh->response_map[5];  // step5 == 2*step

    struct response_layer *layer6 = fh->response_map[6];  // step6 == 4*step
    struct response_layer *layer7 = fh->response_map[7];  // step7 == 4*step

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

    int lobe0 = filter_size0/3;
    int lobe1 = filter_size1/3;
    int lobe2 = filter_size2/3;
    int lobe3 = filter_size3/3;
    int lobe4 = filter_size4/3;
    int lobe5 = filter_size5/3;
    int lobe6 = filter_size6/3;
    int lobe7 = filter_size7/3;

    int border0 = (filter_size0-1)/2;
    int border1 = (filter_size1-1)/2;
    int border2 = (filter_size2-1)/2;
    int border3 = (filter_size3-1)/2;
    int border4 = (filter_size4-1)/2;
    int border5 = (filter_size5-1)/2;
    int border6 = (filter_size6-1)/2;
    int border7 = (filter_size7-1)/2;

    float inv_area0 = 1.f/(filter_size0*filter_size0);
    float inv_area1 = 1.f/(filter_size1*filter_size1);
    float inv_area2 = 1.f/(filter_size2*filter_size2);
    float inv_area3 = 1.f/(filter_size3*filter_size3);
    float inv_area4 = 1.f/(filter_size4*filter_size4);
    float inv_area5 = 1.f/(filter_size5*filter_size5);
    float inv_area6 = 1.f/(filter_size6*filter_size6);
    float inv_area7 = 1.f/(filter_size7*filter_size7);

    float inv_area_squared0 = (inv_area0*inv_area0);
    float inv_area_squared1 = (inv_area1*inv_area1);
    float inv_area_squared2 = (inv_area2*inv_area2);
    float inv_area_squared3 = (inv_area3*inv_area3);
    float inv_area_squared4 = (inv_area4*inv_area4);
    float inv_area_squared5 = (inv_area5*inv_area5);
    float inv_area_squared6 = (inv_area6*inv_area6);
    float inv_area_squared7 = (inv_area7*inv_area7);


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
            // Dxx0 *= inv_area0;
            // Dyy0 *= inv_area0;
            // Dxy0 *= inv_area0;

            // Calculate Determinant
            response0[ind] = inv_area_squared0*(Dxx0 * Dyy0 - 0.81f * Dxy0 * Dxy0);

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
            // Dxx1 *= inv_area1;
            // Dyy1 *= inv_area1;
            // Dxy1 *= inv_area1;

            // Calculate Determinant
            response1[ind] = inv_area_squared1*(Dxx1 * Dyy1 - 0.81f * Dxy1 * Dxy1);

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
            // Dxx2 *= inv_area2;
            // Dyy2 *= inv_area2;
            // Dxy2 *= inv_area2;

            // Calculate Determinant
            response2[ind] = inv_area_squared2*(Dxx2 * Dyy2 - 0.81f * Dxy2 * Dxy2);

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
            // Dxx3 *= inv_area3;
            // Dyy3 *= inv_area3;
            // Dxy3 *= inv_area3;

            // Calculate Determinant
            response3[ind] = inv_area_squared3*(Dxx3 * Dyy3 - 0.81f * Dxy3 * Dxy3);

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
                // Dxx4 *= inv_area4;
                // Dyy4 *= inv_area4;
                // Dxy4 *= inv_area4;

                // Calculate Determinant
                response4[ind2] = inv_area_squared4*(Dxx4 * Dyy4 - 0.81f * Dxy4 * Dxy4);

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
                // Dxx5 *= inv_area5;
                // Dyy5 *= inv_area5;
                // Dxy5 *= inv_area5;

                // Calculate Determinant
                response5[ind2] = inv_area_squared5*(Dxx5 * Dyy5 - 0.81f * Dxy5 * Dxy5);

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
                    // Dxx6 *= inv_area6;
                    // Dyy6 *= inv_area6;
                    // Dxy6 *= inv_area6;

                    // Calculate Determinant
                    response6[ind4] = inv_area_squared6*(Dxx6 * Dyy6 - 0.81f * Dxy6 * Dxy6);

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
                    // Dxx7 *= inv_area7;
                    // Dyy7 *= inv_area7;
                    // Dxy7 *= inv_area7;

                    // Calculate Determinant
                    response7[ind4] = inv_area_squared7*(Dxx7 * Dyy7 - 0.81f * Dxy7 * Dxy7);

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

void compute_response_layer_unconditional_strided(struct response_layer *layer, struct integral_image *iimage) {
    // float Dxx, Dyy, Dxy;
    // int x, y;

    float *response = layer->response;
    bool *laplacian = layer->laplacian;

    int step = layer->step;
    int filter_size = layer->filter_size;
    int height = layer->height;
    int width = layer->width;

    int lobe = filter_size/3;
    int border = (filter_size-1)/2;
    float inv_area = 1.f/(filter_size*filter_size);

    int Dxy_stride = (lobe+1)/step;
    // printf("\nDxy_stride: lobe+1:%d step:%d Dxy_stride:%d\n", lobe+1, step, Dxy_stride);

    // int Dxx_col_stride = lobe/step;
    // printf("Dyy_row_stride/Dxx_col_stride: lobe:%d step:%d Dxx_col_stride:%d\n", lobe, step, Dxx_col_stride);

    int Dxx_col_stride = (2*lobe)/step;
    // printf("2*Dyy_row_stride/2*Dxx_col_stride: lobe:%d step:%d 2*Dxx_col_stride:%d\n", lobe, step, Dxx_col_stride);

    // int Dxx_row_stride = (2*lobe - 1)/step;
    // printf("Dyy_col_stride/Dxx_row_stride: 2*lobe-1:%d step:%d Dxx_row_stride:%d\n", 2*lobe-1, step, Dxx_row_stride);

    // int Dyy_col_stride = (2*lobe - 1)/step;
    int Dyy_row_stride = (2*lobe)/step;

    for (int i = 0, ind = 0; i < height; ++i) {
        int x = i*step;
        for (int j = 0; j < width; ++j, ind++) {
            // Image coordinates
            int y = j*step;

            // Calculate Dxx, Dyy, Dxy with Box Filter
            float Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                    - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
            float Dyy = box_integral_unconditional(iimage, x - border, y - lobe + 1, filter_size, 2*lobe - 1)
                    - 3 * box_integral_unconditional(iimage, x - lobe / 2, y - lobe + 1, lobe, 2*lobe - 1);

            // Normalize Responses with inverse area
            Dxx *= inv_area;
            Dyy *= inv_area;

            // Calculate Determinant
            response[ind] = Dxx * Dyy;
            // Calculate Laplacian
            laplacian[ind] = (Dxx + Dyy >= 0);
        }
    }
    /*
    {
        int i_unroll = 24;
    // printf("\n");
    // printf("Dyy_row_stride %d step %d \n", Dyy_row_stride, step);
    // for (int i_offset=0; i_offset<Dyy_row_stride-i_unroll+1; i_offset+=i_unroll) {
        int i_offset=0;
        for (; i_offset<Dyy_row_stride; i_offset+=1) {
            // printf("j_offset: %d - %d ", j_offset, j_offset+j_unroll-1);
            int i = i_offset;
            for (; i < height-(Dyy_row_stride*i_unroll)+1; i+=Dyy_row_stride*i_unroll) {

                for (int j = 0; j < width; j++) {
                    int y = j*step;
                    int x = i*step;
                    for (int ii = i; ii<i+(Dyy_row_stride*i_unroll); ii+=Dyy_row_stride) {
                        // Image coordinates
                        int x = ii*step;

                        // Calculate Dxx, Dyy, Dxy with Box Filter
                        float Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                                - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
                        float Dyy = box_integral_unconditional(iimage, x - border, y - lobe + 1, filter_size, 2*lobe - 1)
                                - 3 * box_integral_unconditional(iimage, x - lobe / 2, y - lobe + 1, lobe, 2*lobe - 1);

                        // Normalize Responses with inverse area
                        Dxx *= inv_area;
                        Dyy *= inv_area;

                        // Calculate Determinant
                        response[ii*width + j] = Dxx * Dyy;
                        // Calculate Laplacian
                        laplacian[ii*width + j] = (Dxx + Dyy >= 0);
                    }
                }
            }

            for (; i < height; i+=Dxy_stride) {
                int x = i*step;
                for (int j = 0; j < width; j++) {
                    // Image coordinates
                    int y = j*step;

                    // Calculate Dxx, Dyy, Dxy with Box Filter
                    float Dxx = box_integral_unconditional(iimage, x - lobe + 1, y - border, 2*lobe - 1, filter_size)
                            - 3 * box_integral_unconditional(iimage, x - lobe + 1, y - lobe / 2, 2*lobe - 1, lobe);
                    float Dyy = box_integral_unconditional(iimage, x - border, y - lobe + 1, filter_size, 2*lobe - 1)
                            - 3 * box_integral_unconditional(iimage, x - lobe / 2, y - lobe + 1, lobe, 2*lobe - 1);

                    // Normalize Responses with inverse area
                    Dxx *= inv_area;
                    Dyy *= inv_area;

                    // Calculate Determinant
                    response[i*width + j] = Dxx * Dyy;
                    // Calculate Laplacian
                    laplacian[i*width + j] = (Dxx + Dyy >= 0);

                }

            }
        }
    }
    */



    // strided Dxy
    {
        int i_unroll = 12;
    // printf("\n");
    // printf("Dxy_stride %d step %d \n", Dxy_stride, step);
    // for (int i_offset=0; i_offset<Dxy_stride-i_unroll+1; i_offset+=i_unroll) {
        int i_offset=0;
        /**/
        for (; i_offset<Dxy_stride; i_offset+=1) {
            // printf("j_offset: %d - %d ", j_offset, j_offset+j_unroll-1);
            int i = i_offset;
            for (; i < height-(Dxy_stride*i_unroll)+1; i+=Dxy_stride*i_unroll) {

                for (int j = 0; j < width; j++) {
                    int y = j*step;
                    int x = i*step;
                    float bottom_left = box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe);
                    float bottom_right = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe);
                    for (int ii = i; ii<i+(Dxy_stride*i_unroll); ii+=Dxy_stride) {
                        // Image coordinates
                        int x = ii*step;

                        float top_left = bottom_left;//box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe);
                        float top_right = bottom_right;//box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe);

                        bottom_left = box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe);
                        bottom_right = box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);


                        // Calculate Dxx, Dyy, Dxy with Box Filter
                        float Dxy = top_right //box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                                + bottom_left //box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                                - top_left //box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                                - bottom_right; //box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

                        // printf("i:%d j:%d ii:%d jj:%d New Filter: %d %d | %d %d | %d %d | %d %d\n", i, j, ii, jj, x - lobe, y + 1, x + 1, y - lobe, x - lobe, y - lobe, x + 1, y + 1);


                        // Normalize Responses with inverse area
                        Dxy *= inv_area;

                        // Calculate Determinant
                        response[ii*width + j] -= 0.81f * Dxy * Dxy;
                        // printf("ii:%d ",ii);
                    }
                }
            }
            // printf("i:%d \n",i);

            for (; i < height; i+=Dxy_stride) {
                int x = i*step;
                for (int j = 0; j < width; j++) {
                    // Image coordinates
                    int y = j*step;

                    // Calculate Dxx, Dyy, Dxy with Box Filter
                    float Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                            + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                            - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                            - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

                    // printf("i:%d j:%d ii:%d jj:%d New Filter: %d %d | %d %d | %d %d | %d %d\n", i, j, ii, jj, x - lobe, y + 1, x + 1, y - lobe, x - lobe, y - lobe, x + 1, y + 1);


                    // Normalize Responses with inverse area
                    Dxy *= inv_area;

                    // Calculate Determinant
                    response[i*width + j] -= 0.81f * Dxy * Dxy;

                }

            }
        }
    }
        /**/
        /*
        for (; i_offset<Dxy_stride; i_offset+=1) {
            // printf("j_offset: %d ", j_offset);
            for (int i = i_offset; i < height; i+=Dxy_stride) {
                int x = i*step;
                for (int j = 0; j < width; j+=1) {
                    // Image coordinates
                    int y = j*step;

                    // Calculate Dxx, Dyy, Dxy with Box Filter
                    float Dxy = box_integral_unconditional(iimage, x - lobe, y + 1, lobe, lobe)
                            + box_integral_unconditional(iimage, x + 1, y - lobe, lobe, lobe)
                            - box_integral_unconditional(iimage, x - lobe, y - lobe, lobe, lobe)
                            - box_integral_unconditional(iimage, x + 1, y + 1, lobe, lobe);

                    // printf("i:%d j:%d ii:%d jj:%d New Filter: %d %d | %d %d | %d %d | %d %d\n", i, j, ii, jj, x - lobe, y + 1, x + 1, y - lobe, x - lobe, y - lobe, x + 1, y + 1);


                    // Normalize Responses with inverse area
                    Dxy *= inv_area;

                    // Calculate Determinant
                    response[i*width + j] -= 0.81f * Dxy * Dxy;

                }

            }
        }
        */
    // }
}

void compute_response_layers_unconditional_strided(struct fasthessian* fh){

    for (int i = 0; i < 4; ++i) {
		compute_response_layer_unconditional_strided(fh->response_map[i], fh->iimage);
	}
    for (int i = 4; i < fh->total_layers; ++i) {
		compute_response_layer_unconditional(fh->response_map[i], fh->iimage);
	}
}

void get_interest_points_block(struct fasthessian *fh, std::vector<struct interest_point> *interest_points) {

    const int block_size = 16;

    assert(fh != NULL);
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
    for (int o = 0; o < fh->octaves; ++o) {

        // TODO: (Sebastian) allow for fh->layers != 4 as well (note that fh->layers>=3 has to hold)
        for (int i = 0; i <= 1; ++i) {

            // assigning respective bottom, middle and top response layer
            bottom = fh->response_map[filter_map[o][i]];
            middle = fh->response_map[filter_map[o][i+1]];
            top = fh->response_map[filter_map[o][i+2]];

            // iterating over middle response layer at density of the most sparse layer (always top),
            // to find maxima accreoss scale and space
            for (int r = 0; r < top->height; r+=block_size) {
                for (int c = 0; c < top->width; c+=block_size) {
                    int br_max = MIN(r + block_size, top->height);
                    for (int br = r; br < br_max; ++br) {
                        int bc_max = MIN(c + block_size, top->width);
                        for (int bc = c; bc < bc_max; ++bc) {
                            // checking if current pixel position is local maximum in 3x3x3 maximum and above threshold
                            if (is_extremum(br, bc, top, middle, bottom, fh->thresh)) {

                                // sub-pixel interpolating local maxium and adding to resulting interest point vector
                                interpolate_extremum(br, bc, top, middle, bottom, interest_points);

                            }
                        }
                    }
                }
            }

        }
    }
}
