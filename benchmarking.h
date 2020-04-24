#pragma once

#include <stdlib.h>
#include "integral_image.h"

double perf_test_integral_img(struct integral_image* (*function)(float*, int, int), float* gray_image, int width, int height,int flops);
void bench_integral_img(float* image, int width, int height);
