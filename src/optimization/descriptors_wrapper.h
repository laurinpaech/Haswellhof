#pragma once

#include "integral_image.h"
#include "interest_point.h"
#include "descriptor_opt.h"

#include <vector>

#define DECLARE_GET_MSURF_DESCRIPTORS(get_msurf_descriptors_name) \
    void get_msurf_descriptors_name(struct integral_image *iimage, std::vector<struct interest_point> *interest_points);

#define DEFINE_GET_MSURF_DESCRIPTORS(get_msurf_descriptors_name, get_msurf_descriptor_name) \
    void get_msurf_descriptors_name(struct integral_image *iimage, std::vector<struct interest_point> *interest_points) { \
        for (int i = 0; i < interest_points->size(); ++i) { \
            get_msurf_descriptor_name(iimage, &interest_points->at(i)); \
        } \
    }


/*
DECLARE_GET_MSURF_DESCRIPTORS(get_msurf_descriptors_improved)

DECLARE_GET_MSURF_DESCRIPTORS(get_msurf_descriptors_improved_flip)

DECLARE_GET_MSURF_DESCRIPTORS(get_msurf_descriptors_improved_flip_flip)

DECLARE_GET_MSURF_DESCRIPTORS(get_msurf_descriptors_inlined)

DECLARE_GET_MSURF_DESCRIPTORS(get_msurf_descriptors_gauss_s1_separable_test)

DECLARE_GET_MSURF_DESCRIPTORS(get_msurf_descriptors_gauss_s2_precomputed)

DECLARE_GET_MSURF_DESCRIPTORS(get_msurf_descriptors_inlinedHaarWavelets)

DECLARE_GET_MSURF_DESCRIPTORS(get_msurf_descriptors_inlinedHaarWavelets_precheck_boundaries)

DECLARE_GET_MSURF_DESCRIPTORS(get_msurf_descriptors_gauss_compute_once_case)
*/
