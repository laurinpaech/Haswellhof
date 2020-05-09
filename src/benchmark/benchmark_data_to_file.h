#pragma once

#include <vector>

#include "benchmarking.h"

// Saves a vector of different benchmarking data to their respective files. One file contains the timing data of
// multiple images and one function.
void save_benchmark_data(const std::vector<struct benchmark_data> &all_benchmark_data);

// Saves a file containing the information of benchmark_data
void save_performance_file(const struct benchmark_data &data, char *folder_name);

// Creates the folder with the specified name if it doesn't exist.
// For nested directories the previous folders must exist.
void create_folder(const char *folder_name);

// Concatenates two given strings
char *concat(const char *s1, const char *s2);