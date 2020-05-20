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

inline void compute_response_layer_Dyy_laplacian_localityloops_inlined(struct response_layer* layer, struct integral_image* iimage) {
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

inline void compute_response_layer_Dyy_laplacian_locality_uncond_opt_flops_inlined(struct response_layer* layer, struct integral_image* iimage) {
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
