#include "fasthessian_opt.h"


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
r11 = r10 + lobe - 1;
c10 = y - lobe;
c11 = y + lobe - 1;

A = (r10, c10)
B = (r10, c11) = (x - lobe / 2 - 1, y + lobe - 1)
C = (r11, c10)
D = (r11, c11)
*/

void compute_response_layer_Dyy_leftcorner(struct response_layer* layer, struct integral_image* iimage) {
    float Dxx, Dyy, Dxy;
    int x, y;
    int r0, r1, c0, c1, r00, r01, c00, c01, r10, r11, c10, c11;
    float Dyy0, Dyy1;
    float A, B, C, D;

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

    float *data = (float*) iimage->data;  // brauch hier keinen cast weil es eig float sein sollte
    int iheight = iimage->height;
    int iwidth = iimage->width;

    // Top Left Corner - Case 1: B of neg part outside
    for (int i = 0; i < lobe/2+1; i++) {    // Inner B is outside, i.e. 0
        for (int j = 0; j < lobe; j++) {    // c0 = col - 1 = (y - lobe + 1) -1 < 0
            // Image coordinates
            x = i*step;
            y = j*step;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * iwidth + c01];

            // neg part box filter
            r11 = r10 + lobe - 1;
            c11 = y + lobe - 1;

            D = data[r11 * iwidth + c11];
            Dyy1 = D;

            Dyy = Dyy0 - 3 * Dyy1;

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

    // Top Left Corner - Case 2: B of neg part inside
    for (int i = lobe/2+1; i < border+1; i++) {     // Inner B is inside
        for (int j = 0; j < lobe; j++) {
            // Image coordinates
            x = i*step;
            y = j*step;

            // Compute Dyy
            // whole box filter
            r01 = x + border;
            c01 = y + lobe - 1;

            Dyy0 = data[r01 * iwidth + c01];

            // neg part box filter
            r10 = x - lobe / 2 - 1;
            r11 = r10 + lobe;  // -1 is already part of r10
            c11 = y + lobe - 1;

            B = data[r10 * iwidth + c11];
            D = data[r11 * iwidth + c11];
            Dyy1 = D - B;

            Dyy = Dyy0 - 3 * Dyy1;

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

    // Rest
    for (int i = 0; i < height; ++i) {
        for (int j = lobe; j < width; ++j) {
            // Image coordinates
            x = i*step;
            y = j*step;

            // Calculate Dxx, Dyy, Dxy with Box Filter
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
    for (int i = border+1; i < height; ++i) {
        for (int j = 0; j < lobe; ++j) {
            // Image coordinates
            x = i*step;
            y = j*step;

            // Calculate Dxx, Dyy, Dxy with Box Filter
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

}

void compute_response_layers_at_once(struct fasthessian* fh, struct integral_image* iimage) {
    /* computes all 8 response layers at once, gives same results as base implementation
    valgrind reports no improvement for l1 misses, i.e. locality is not improved as expected
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

    // step == 1: 0,1,2,3
    // step == 2: 4,5
    // step == 4: 6,7

    struct response_layer* layer0 = fh->response_map[0]; // step0 == step
    struct response_layer* layer1 = fh->response_map[1]; // step1 == step
    struct response_layer* layer2 = fh->response_map[2]; // step2 == step
    struct response_layer* layer3 = fh->response_map[3]; // step3 == step

    struct response_layer* layer4 = fh->response_map[4]; // step4 == 2*step
    struct response_layer* layer5 = fh->response_map[5]; // step5 == 2*step

    struct response_layer* layer6 = fh->response_map[6]; // step6 == 4*step
    struct response_layer* layer7 = fh->response_map[7]; // step7 == 4*step

    float* response0 = layer0->response;
    float* response1 = layer1->response;
    float* response2 = layer2->response;
    float* response3 = layer3->response;
    float* response4 = layer4->response;
    float* response5 = layer5->response;
    float* response6 = layer6->response;
    float* response7 = layer7->response;

    bool* laplacian0 = layer0->laplacian;
    bool* laplacian1 = layer1->laplacian;
    bool* laplacian2 = layer2->laplacian;
    bool* laplacian3 = layer3->laplacian;
    bool* laplacian4 = layer4->laplacian;
    bool* laplacian5 = layer5->laplacian;
    bool* laplacian6 = layer6->laplacian;
    bool* laplacian7 = layer7->laplacian;

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


    for (int i = 0, ind = 0, ind2 = 0, ind4 = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j, ind++) {
            // Image coordinates
            x = i*step;
            y = j*step;

            // layer0
            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx0 = box_integral(iimage, x - lobe0 + 1, y - border0, 2*lobe0 - 1, filter_size0)
                    - 3 * box_integral(iimage, x - lobe0 + 1, y - lobe0 / 2, 2*lobe0 - 1, lobe0);
            Dyy0 = box_integral(iimage, x - border0, y - lobe0 + 1, filter_size0, 2*lobe0 - 1)
                    - 3 * box_integral(iimage, x - lobe0 / 2, y - lobe0 + 1, lobe0, 2*lobe0 - 1);
            Dxy0 = box_integral(iimage, x - lobe0, y + 1, lobe0, lobe0)
                    + box_integral(iimage, x + 1, y - lobe0, lobe0, lobe0)
                    - box_integral(iimage, x - lobe0, y - lobe0, lobe0, lobe0)
                    - box_integral(iimage, x + 1, y + 1, lobe0, lobe0);

            // Normalize Responses with inverse area
            Dxx0 *= inv_area0;
            Dyy0 *= inv_area0;
            Dxy0 *= inv_area0;

            // Calculate Determinant
            response0[ind] = Dxx0 * Dyy0 - 0.81f * Dxy0 * Dxy0;

            // Calculate Laplacian
            laplacian0[ind] = (Dxx0 + Dyy0 >= 0 ? true : false);


            // layer1
            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx1 = box_integral(iimage, x - lobe1 + 1, y - border1, 2*lobe1 - 1, filter_size1)
                    - 3 * box_integral(iimage, x - lobe1 + 1, y - lobe1 / 2, 2*lobe1 - 1, lobe1);
            Dyy1 = box_integral(iimage, x - border1, y - lobe1 + 1, filter_size1, 2*lobe1 - 1)
                    - 3 * box_integral(iimage, x - lobe1 / 2, y - lobe1 + 1, lobe1, 2*lobe1 - 1);
            Dxy1 = box_integral(iimage, x - lobe1, y + 1, lobe1, lobe1)
                    + box_integral(iimage, x + 1, y - lobe1, lobe1, lobe1)
                    - box_integral(iimage, x - lobe1, y - lobe1, lobe1, lobe1)
                    - box_integral(iimage, x + 1, y + 1, lobe1, lobe1);

            // Normalize Responses with inverse area
            Dxx1 *= inv_area1;
            Dyy1 *= inv_area1;
            Dxy1 *= inv_area1;

            // Calculate Determinant
            response1[ind] = Dxx1 * Dyy1 - 0.81f * Dxy1 * Dxy1;

            // Calculate Laplacian
            laplacian1[ind] = (Dxx1 + Dyy1 >= 0 ? true : false);


            // layer2
            // Calculate Dxx, Dyy, Dxy with Box Filter
            Dxx2 = box_integral(iimage, x - lobe2 + 1, y - border2, 2*lobe2 - 1, filter_size2)
                    - 3 * box_integral(iimage, x - lobe2 + 1, y - lobe2 / 2, 2*lobe2 - 1, lobe2);
            Dyy2 = box_integral(iimage, x - border2, y - lobe2 + 1, filter_size2, 2*lobe2 - 1)
                    - 3 * box_integral(iimage, x - lobe2 / 2, y - lobe2 + 1, lobe2, 2*lobe2 - 1);
            Dxy2 = box_integral(iimage, x - lobe2, y + 1, lobe2, lobe2)
                    + box_integral(iimage, x + 1, y - lobe2, lobe2, lobe2)
                    - box_integral(iimage, x - lobe2, y - lobe2, lobe2, lobe2)
                    - box_integral(iimage, x + 1, y + 1, lobe2, lobe2);

            // Normalize Responses with inverse area
            Dxx2 *= inv_area2;
            Dyy2 *= inv_area2;
            Dxy2 *= inv_area2;

            // Calculate Determinant
            response2[ind] = Dxx2 * Dyy2 - 0.81f * Dxy2 * Dxy2;

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
            Dxx3 *= inv_area3;
            Dyy3 *= inv_area3;
            Dxy3 *= inv_area3;

            // Calculate Determinant
            response3[ind] = Dxx3 * Dyy3 - 0.81f * Dxy3 * Dxy3;

            // Calculate Laplacian
            laplacian3[ind] = (Dxx3 + Dyy3 >= 0 ? true : false);


            // 2*step
            if (i%2 == 0 && j%2 == 0) {
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
                Dxx4 *= inv_area4;
                Dyy4 *= inv_area4;
                Dxy4 *= inv_area4;

                // Calculate Determinant
                response4[ind2] = Dxx4 * Dyy4 - 0.81f * Dxy4 * Dxy4;

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
                Dxx5 *= inv_area5;
                Dyy5 *= inv_area5;
                Dxy5 *= inv_area5;

                // Calculate Determinant
                response5[ind2] = Dxx5 * Dyy5 - 0.81f * Dxy5 * Dxy5;

                // Calculate Laplacian
                laplacian5[ind2] = (Dxx5 + Dyy5 >= 0 ? true : false);


                ind2++;

                // 4*step
                if (i%4 == 0 && j%4 == 0) {
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
                    Dxx6 *= inv_area6;
                    Dyy6 *= inv_area6;
                    Dxy6 *= inv_area6;

                    // Calculate Determinant
                    response6[ind4] = Dxx6 * Dyy6 - 0.81f * Dxy6 * Dxy6;

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
                    Dxx7 *= inv_area7;
                    Dyy7 *= inv_area7;
                    Dxy7 *= inv_area7;

                    // Calculate Determinant
                    response7[ind4] = Dxx7 * Dyy7 - 0.81f * Dxy7 * Dxy7;

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

void interpolate_step_gauss(int row, int col, struct response_layer *top, struct response_layer *middle, struct response_layer *bottom, float offsets[3]) {

}
