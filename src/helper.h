#pragma once

#include <stdlib.h>
#include <assert.h>

void solve_linear_3x3_system(float A[9], float b[3], float ret_x[3]) {

    assert((A != NULL) && (b != NULL) && (ret_x != NULL));

    // computing determinant and inverse determinant of A
    float det = 0.0f;
    for (int i = 0; i < 3; ++i) {
       det += A[0*3 + i] * (A[1*3 + (i+1)%3] * A[2*3 + (i+2)%3] - A[1*3 + (i+2)%3] * A[2*3 + (i+1)%3]);
    }
    float det_inv = 1.0f / det;

    // computing inverse of 3x3 matrix A with inverse determinant
    float A_inv[9];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            A_inv[i*3 + j] = det_inv * ((A[((j+1)%3)*3 + (i+1)%3] * A[((j+2)%3)*3 + (i+2)%3]) - (A[((j+1)%3)*3 + (i+2)%3] * A[((j+2)%3)*3 + (i+1)%3]));
        }
    }

    // computing solution vector x = A_inv * b
    ret_x[0] = 0.0f;
    ret_x[1] = 0.0f;
    ret_x[2] = 0.0f;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            ret_x[i] += A_inv[i*3 + j] * b[j];
        }
    }

}
