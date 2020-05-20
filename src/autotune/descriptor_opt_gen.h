#include "integral_image.h"
#include "interest_point.h"
#include "helper.h"

#include <stdlib.h>
#include <vector>

void get_msurf_descriptor_haar_unroll_1_1_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_1_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_2_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_2_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_3_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_3_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_4_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_4_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_6_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_6_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_6_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_6_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_8_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_8_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_8_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_8_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_12_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_12_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_12_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_12_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_24_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_24_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_24_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_24_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_1_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_1_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_2_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_2_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_3_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_3_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_4_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_4_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_6_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_6_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_6_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_6_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_8_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_8_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_8_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_8_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_12_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_12_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_12_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_12_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_24_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_24_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_24_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_24_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_1_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_1_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_2_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_2_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_3_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_3_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_4_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_4_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_6_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_6_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_6_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_6_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_8_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_8_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_8_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_8_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_12_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_12_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_12_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_12_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_24_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_24_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_24_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_24_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_1_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_1_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_2_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_2_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_3_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_3_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_4_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_4_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_6_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_6_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_6_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_6_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_8_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_8_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_8_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_8_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_12_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_12_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_12_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_12_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_24_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_24_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_24_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_24_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_1_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_1_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_2_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_2_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_3_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_3_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_4_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_4_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_6_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_6_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_6_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_6_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_8_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_8_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_8_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_8_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_12_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_12_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_12_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_12_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_24_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_24_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_6_24_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_6_24_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_1_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_1_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_2_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_2_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_3_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_3_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_4_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_4_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_6_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_6_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_6_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_6_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_8_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_8_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_8_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_8_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_12_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_12_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_12_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_12_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_24_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_24_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_8_24_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_8_24_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_1_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_1_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_2_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_2_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_3_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_3_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_4_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_4_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_6_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_6_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_6_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_6_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_8_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_8_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_8_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_8_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_12_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_12_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_12_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_12_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_24_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_24_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_12_24_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_12_24_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_1_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_1_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_2_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_2_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_3_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_3_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_4_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_4_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_6_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_6_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_6_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_6_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_8_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_8_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_8_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_8_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_12_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_12_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_12_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_12_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_24_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_24_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_24_24_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_24_24_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

