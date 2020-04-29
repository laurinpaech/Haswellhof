#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <inttypes.h>

#include "benchmark_data_to_file.h"

void save_benchmark_data(std::vector<struct benchmark_data *> all_benchmark_data)
{
    // get the current time
    struct tm *timenow;
    time_t now = time(NULL);
    timenow = gmtime(&now);
    char folder_name[32] = "../benchmarking_files";
    create_folder(folder_name);
    printf("created first folder\n");

    char current_time[32];
    strftime(current_time, sizeof(current_time), "/%Y_%m_%d_%H_%M", timenow);

    char *nested_folder = concat(folder_name, current_time);
    create_folder(nested_folder);
    printf("created second folder\n");

    // concatenate the path in the form  "benchmarking_files/CURRENT_DATE/function_name"
    for (int i = 0; i < all_benchmark_data.size(); i++)
    {
        save_performance_file(all_benchmark_data[i], nested_folder);
    }

    free(nested_folder);
}

void create_folder(const char *folder_name)
{
    // create the folder if it doesn't exist
    struct stat st = {0};
    if (stat(folder_name, &st) == -1)
    {
        mkdir(folder_name, 0700);
    }
}

// Saves a file containing the information of benchmark_data
void save_performance_file(struct benchmark_data *data, char *folder_name)
{
    char *path_name = concat(folder_name, "/");
    //path_name = concat(path_name, data->image_name);
    path_name = concat(path_name, data->function_name);

    printf("%s\n", path_name);
    path_name = concat(path_name, ".txt");
    printf("%s\n", path_name);

    // open the file for reading and writing - create it if necessary
    FILE *fp = fopen(path_name, "a+");
    printf("opened %s\n", path_name);
    fprintf(fp, "%s,%d,%d,%d,%ld,%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%lf\n", data->image_name, data->width, data->height, data->num_interest_points, data->num_flops, data->avg_cycles, data->min_cycles, data->max_cycles, data->flops_per_cycle);

    // closes the file pointed by fptr_int_img
    fclose(fp);
    free(path_name);
}

char *concat(const char *s1, const char *s2)
{
    const size_t len1 = strlen(s1);
    const size_t len2 = strlen(s2);
    char *result = (char *)malloc(len1 + len2 + 1); // +1 for the null-terminator
    // in real code you would check for errors in malloc here
    memcpy(result, s1, len1);
    memcpy(result + len1, s2, len2 + 1); // +1 to copy the null-terminator
    return result;
}
