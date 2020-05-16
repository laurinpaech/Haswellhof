#include "integral_image.h"
#include "interest_point.h"
#include "helper.h"

#include <stdlib.h>
#include <vector>

void get_msurf_descriptor_haar_unroll_1_1_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_1_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_2_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_2_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_3_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_3_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_4_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_1_4_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_1_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_1_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_1_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_2_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_2_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_3_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_3_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_4_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_2_4_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_2_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_1_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_1_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_2_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_2_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_3_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_3_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_4_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_3_4_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_3_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_1_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_1_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_1_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_1_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_2_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_2_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_2_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_2_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_3_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_3_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_3_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_3_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_4_False(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_4_False(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

void get_msurf_descriptor_haar_unroll_4_4_True(struct integral_image* iimage, struct interest_point* ipoint);

void get_msurf_descriptors_haar_unroll_4_4_True(struct integral_image* iimage, std::vector<struct interest_point> *interest_points);

