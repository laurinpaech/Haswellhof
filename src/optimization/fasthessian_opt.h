#pragma once

#include "fasthessian.h"
#include "integral_image_opt.h"

void height_greater_border_width_greater_double_lobe_Dyy(struct response_layer *layer, struct integral_image *iimage);

void height_greater_border_width_less_double_lobe_Dyy(struct response_layer *layer, struct integral_image *iimage);

void get_interest_points_layers(struct fasthessian *fh, std::vector<struct interest_point> *interest_points);

void interpolate_step_gauss(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets[3]);

void compute_response_layer_Dyy_leftcorner(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layer_precompute(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_precompute(struct fasthessian* fh);

void compute_response_layer_blocking(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_blocking(struct fasthessian* fh);

void compute_response_layer_Dyy_top(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layer_Dyy_top_mid(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_Dyy(struct fasthessian *fh);

void compute_response_layer_Dyy(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layer_sonic_Dyy(struct response_layer *layer, struct integral_image *iimage);

void compute_response_layer_Dyy_laplacian(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_sonic_Dyy(struct fasthessian *fh);

void compute_response_layers_Dyy_laplacian(struct fasthessian *fh);

void compute_response_layers_Dyy_laplacian_localityloops(struct fasthessian *fh);

void compute_response_layer_Dyy_laplacian_localityloops(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_at_once(struct fasthessian* fh);

void compute_response_layers_unconditional(struct fasthessian* fh);

void compute_response_layer_unconditional(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_unconditional_strided(struct fasthessian* fh);

void compute_response_layer_unconditional_strided(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_sonic_Dyy_unconditional(struct fasthessian *fh);

void compute_response_layer_sonic_Dyy_unconditional(struct response_layer *layer, struct integral_image *iimage);

void compute_response_layer_Dyy_laplacian_localityloops_unconditional(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_Dyy_laplacian_localityloops_unconditional(struct fasthessian* fh);

void compute_response_layers_Dyy_laplacian_locality_uncond_opt(struct fasthessian* fh);

void compute_response_layer_Dyy_laplacian_locality_uncond_opt(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_Dyy_laplacian_locality_uncond_opt_flops(struct fasthessian* fh);

void compute_response_layer_Dyy_laplacian_locality_uncond_opt_flops(struct response_layer* layer, struct integral_image* iimage);

void compute_response_layers_Dyy_laplacian_locality_uncond_opt_flops_invsqr(struct fasthessian* fh);

void compute_response_layer_Dyy_laplacian_locality_uncond_opt_flops_invsqr(struct response_layer* layer, struct integral_image* iimage);

inline void height_greater_border_width_greater_double_lobe_Dyy_inlined(struct response_layer *layer, struct integral_image *iimage) {
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

inline void height_greater_border_width_less_double_lobe_Dyy_inlined(struct response_layer *layer, struct integral_image *iimage) {
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
