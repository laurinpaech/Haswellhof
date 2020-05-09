#pragma once

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

// https://stackoverflow.com/questions/3437404/min-and-max-in-c
#define MIN(a, b)               \
    ({                          \
        __typeof__(a) _a = (a); \
        __typeof__(b) _b = (b); \
        _a < _b ? _a : _b;      \
    })

#define MAX(a, b)               \
    ({                          \
        __typeof__(a) _a = (a); \
        __typeof__(b) _b = (b); \
        _a > _b ? _a : _b;      \
    }) 

#define EPSILON (1e-3)

void solve_linear_3x3_system(float A[9], float b[3], float ret_x[3]);

inline bool compare_arrays(float a[], float b[], int n) {
    for (int i = 0; i < n; ++i) {
        if (a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

inline bool compare_arrays_close(float a[], float b[], int n, float epsilon = EPSILON) {
    for (int i = 0; i < n; ++i) {
        if (fabsf(a[i] - b[i]) > epsilon) {
            return false;
        }
    }
    return true;
}

// Compares two matrices and checks if the values are equal.
// Returns true if all values of the matrix are equal, false otherwise.
static inline bool are_matrices_equal(float* iimage1, float* iimage2, int width, int height){
     for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (iimage1[i * width + j] != iimage2[i * width + j]) {
                return false;
            }
        }
    }
    return true;
}
