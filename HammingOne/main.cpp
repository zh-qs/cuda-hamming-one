#include "structures/vector_table.h"
#include "structures/time_statistics.h"
#include "util/util.h"
#include "cuda/vector_operations.cuh"

#include "cpu/cpu_vector_operations.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <iostream>
#include <fstream>
#include <cstring>

#define READ_ARGUMENTS

// allocates vector_table table in host memory and fills details (size, bit count, ...) fields with values
void alloc_table(vector_table& table, int n, int length)
{
    const int bits_in_level = BITS_IN_LEVEL;

    table.n = n;

    if (length % bits_in_level == 0)
        table.l = length / bits_in_level;
    else
        table.l = length / bits_in_level + 1;
    table.bit_count = length;
    table.levels = new bits_t * [table.l];
    bits_t* bits = new bits_t[table.l * table.n];
    table.first_differences = new int[table.n];
    for (int i = 0; i < table.l; ++i)
        table.levels[i] = bits + i * table.n;
}

bool read_arguments(char *filename, vector_table &table)
{
    std::ifstream s(filename);
    int length;
    std::string str;
    const int bits_in_level = BITS_IN_LEVEL;
    if (s.good())
    {
        s >> table.n >> length;
        if (length % bits_in_level == 0)
            table.l = length / bits_in_level;
        else
            table.l = length / bits_in_level + 1;
        table.bit_count = length;
        table.levels = new bits_t*[table.l];
        bits_t* bits = new bits_t[table.l * table.n];
        table.first_differences = new int[table.n];
        for (int i = 0; i < table.l; ++i)
            table.levels[i] = bits + i * table.n;
        for (int i = 0; i < table.n; ++i)
        {
            s >> str;
            for (int j = 0; j < table.l; ++j)
            {
                table.levels[j][i] = bit_string_to_bits(str.substr(j * bits_in_level, bits_in_level));
            }
        }
        return true;
    }
    return false;
}

void delete_table(vector_table& table)
{
    delete[] table.levels[0];
    delete[] table.levels;
    delete[] table.first_differences;
}

// runs GPU algorithm and stores time measurements
int run_cuda_algorithm(vector_table &h_table, time_statistics &stats)
{
    clock_t begin, end;

    vector_table d_table, d_rev_table;

    begin = clock();
    alloc_and_copy_to_device(h_table, d_table);
    end = clock();
    stats.memcpy_h_to_d = end - begin;

    alloc_only(h_table, d_rev_table);

    int* d_indexes, * d_rev_indexes, * d_diffs, * d_rev_diffs, * d_res, * d_abstract_classes, * d_rev_abstract_classes, * d_clear_indexes;
    bits_t* d_level;
    unsigned long long* d_class_pairs;
    cudaError_t err = cudaMalloc(&d_class_pairs, ((7 + d_table.bit_count) * sizeof(int) + sizeof(unsigned long long) + sizeof(bits_t)) * d_table.n);

    if (err != cudaSuccess)
    {
        std::cout << "cudaMalloc error: " << cudaGetErrorString(err) << std::endl;
        delete_table(h_table);
        delete_device_table(d_table);
        delete_device_table(d_rev_table);
        return 1;
    }

    d_diffs = (int*)(d_class_pairs + d_table.n);
    d_rev_diffs = d_diffs + d_table.n;
    d_res = d_rev_diffs + d_table.n;
    d_abstract_classes = d_res + d_table.n * d_table.bit_count;
    d_rev_abstract_classes = d_abstract_classes + d_table.n;
    d_clear_indexes = d_rev_abstract_classes + d_table.n;
    d_indexes = d_clear_indexes + d_table.n;
    d_rev_indexes = d_indexes + d_table.n;
    d_level = (bits_t*)(d_rev_indexes + d_table.n);

    cudaMemset(d_res, 0, sizeof(int) * d_table.n * d_table.bit_count);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);

    reverse_vectors(d_table, d_rev_table);
  
    radix_sort(d_table, d_indexes, d_level);
    radix_sort(d_rev_table, d_rev_indexes, d_level);

    calculate_first_differences_with_previous_vector(d_table, d_indexes, d_diffs);
    calculate_first_differences_with_previous_vector(d_rev_table, d_rev_indexes, d_rev_diffs);

    for (int i = 0; i < d_table.bit_count; ++i)
    {
        calculate_abstract_classes_and_clear_indexes(d_table, d_rev_table, d_indexes, d_rev_indexes, d_diffs, d_rev_diffs, d_abstract_classes, d_rev_abstract_classes, d_class_pairs, i, d_clear_indexes);

        get_hamming_one_vector_pairs_to_array(d_table, d_clear_indexes, d_res, d_class_pairs, i);
    }

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&stats.cuda_algorithm, start, stop);

    int* h_res = new int[h_table.n * h_table.bit_count];

    begin = clock();
    cudaMemcpy(h_res, d_res, sizeof(int) * h_table.n * h_table.bit_count, cudaMemcpyDeviceToHost);
    end = clock();
    stats.memcpy_d_to_h = end - begin;

    begin = clock();
    int index;
    for (int i = 0; i < h_table.n * h_table.bit_count; ++i)
    {
        index = i / h_table.bit_count + 1;
        if (h_res[i] > 0)
        {
            printf("%d;%d\n", index, h_res[i]);
        }
    }
    end = clock();
    stats.data_print = end - begin;

    delete[] h_res;
    cudaFree(d_class_pairs);
    cudaFree(d_indexes);
    cudaFree(d_rev_indexes);
    delete_device_table(d_table);
    delete_device_table(d_rev_table);

    return 0;
}
// runs CPU algorithm and stores time measurements; print_results -- flag determining if CPU output should be printed on stdout
int run_cpu_algorithm(vector_table &table, time_statistics &stats, bool print_results = true)
{
    int* indexes, * rev_indexes, * diffs, * rev_diffs, *abstract_classes, *rev_abstract_classes, *clear_indexes;
    unsigned long long* class_pairs;

    vector_table rev_table;

    alloc_table(rev_table, table.n, table.bit_count);

    diffs = new int[table.n];
    rev_diffs = new int[table.n];
    abstract_classes = new int[table.n];
    rev_abstract_classes = new int[table.n];
    clear_indexes = new int[table.n];
    class_pairs = new unsigned long long[table.n];

    clock_t begin = clock();

    cpu_radix_sort_and_alloc_indexes(table, &indexes);

    cpu_reverse_vectors(table, rev_table);
    cpu_radix_sort_and_alloc_indexes(rev_table, &rev_indexes);

    cpu_calculate_first_differences_with_previous_vector(table, indexes, diffs);
    cpu_calculate_first_differences_with_previous_vector(rev_table, rev_indexes, rev_diffs);

    for (int i = 0; i < table.bit_count; ++i)
    {
        cpu_calculate_abstract_classes_and_clear_indexes(table, rev_table, indexes, rev_indexes, diffs, rev_diffs, abstract_classes, rev_abstract_classes, class_pairs, i, clear_indexes);

        cpu_print_found_hamming_one_vector_pairs(table, clear_indexes, class_pairs, i, print_results);
    }

    clock_t end = clock();
    stats.cpu_algorithm = end - begin;

    delete[] indexes;
    delete[] rev_indexes;
    delete[] diffs;
    delete[] rev_diffs;
    delete[] abstract_classes;
    delete[] rev_abstract_classes;
    delete[] clear_indexes;
    delete[] class_pairs;

    return 0;
}

int main(int argc, char **argv)
{
    const char csv_separator = ';';
    bool only_cuda_output = true;

    vector_table h_table;
    time_statistics stats;

#ifdef READ_ARGUMENTS

    // USAGE: prog_name input_file [csv_results [-pcpu]]
    // input_file -- file containing list of vectors of the same length (L) in format:
    //              N L
    //              vector1
    //              vector2
    //              ...
    //              vectorN
    // csv_results -- CSV file to store time measurements
    // -pcpu -- optional, if present, CPU results will also be printed

    if (argc < 2 || argc > 4)
    {
        std::cout << "USAGE: %s input_file [csv_results [-pcpu]]" << std::endl;
        return 1;
    }
    if (argc == 4)
    {
        if (strncmp(argv[3], "-pcpu", 5) == 0) only_cuda_output = false;
        else
        {
            std::cout << "USAGE: %s input_file [csv_results [-pcpu]]" << std::endl;
            return 1;
        }
    }

    clock_t begin = clock();
    if (!read_arguments(argv[1], h_table)) return 1;
    clock_t end = clock();
    stats.data_read = end - begin;
#else
    h_table.n = 3;
    h_table.bit_count = 4;
    h_table.l = 1;
    h_table.first_differences = new int[h_table.n];
    bits_t* bits = new bits_t[h_table.l * h_table.n];
    h_table.levels = new bits_t*[h_table.l];
    h_table.levels[0] = bits;

    h_table.levels[0][0] = 0b1011 << 28;
    h_table.levels[0][1] = 0b0011 << 28;
    h_table.levels[0][2] = 0b1000 << 28;
    h_table.table = bits;
#endif // READ_ARGUMENTS

    stats.n = h_table.n;
    stats.l = h_table.bit_count;

    if (!only_cuda_output) 
        std::cout << "Output from GPU:" << std::endl;

    int res = run_cuda_algorithm(h_table, stats);

    if (argc >= 3) // we want to compare gpu vs cpu and write results to file
    {
        if (!only_cuda_output)
            std::cout << std::endl << "Output from CPU:" << std::endl;

        run_cpu_algorithm(h_table, stats, !only_cuda_output);

        std::ofstream s(argv[2], std::ios_base::app);
        if (s.good())
        {
            s << argv[1] << csv_separator << stats.n << csv_separator << stats.l << csv_separator << stats.data_read << csv_separator << stats.memcpy_h_to_d << csv_separator << stats.cuda_algorithm << csv_separator << stats.memcpy_d_to_h << csv_separator << stats.cpu_algorithm << csv_separator << stats.data_print << std::endl;
            s.close();
        }
    }

    delete_table(h_table);

    return res;
}

