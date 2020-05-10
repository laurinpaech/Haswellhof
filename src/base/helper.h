#pragma once

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

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

// Compares two matrices of floats and checks if the values are equal.
// Returns true if all values of the matrix are equal, false otherwise.
inline bool are_float_matrices_equal(float *matrix1, float *matrix2, int width, int height) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (matrix1[i * width + j] != matrix2[i * width + j]) {
                return false;
            }
        }
    }
    return true;
}

// Compares two matrices of booleans and checks if the values are equal.
// Returns true if all values of the matrix are equal, false otherwise.
inline bool are_bool_matrices_equal(bool *matrix1, bool *matrix2, int width, int height) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (matrix1[i * width + j] != matrix2[i * width + j]) {
                return false;
            }
        }
    }
    return true;
}

// Alignment only works for powers of two!
#define ALIGN(x,a)              __ALIGN_MASK(x,(typeof(x))(a)-1)
#define __ALIGN_MASK(x,mask)    (((x)+(mask))&~(mask))

// https://stackoverflow.com/questions/38088732/explanation-to-aligned-malloc-implementation
// Alignment only works for powers of two!
inline void *aligned_malloc(size_t required_bytes, size_t alignment) {
    void *p1;  // original block
    void **p2; // aligned block
    int offset = alignment - 1 + sizeof(void *);
    if ((p1 = (void *)malloc(required_bytes + offset)) == NULL) {
       return NULL;
    }
    p2 = (void **)(((size_t)(p1) + offset) & ~(alignment - 1));
    p2[-1] = p1;
    return p2;
}

inline void aligned_free(void *p) {
    free(((void **)p)[-1]);
}

void solve_linear_3x3_system(float A[9], float b[3], float ret_x[3]);
