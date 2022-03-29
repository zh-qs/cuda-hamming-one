#ifndef VECTOR_OPERATIONS_CUH
#define VECTOR_OPERATIONS_CUH

#include "../structures/vector_table.h"

// allocates d_table on GPU with same parameters as h_table and copies contents of h_table into d_table
void alloc_and_copy_to_device(vector_table& h_table, vector_table& d_table);
// allocates d_table on GPU with same parameters as h_table
void alloc_only(vector_table& h_table, vector_table& d_table);
// copies d_table content to h_table
void copy_to_host(vector_table& h_table, vector_table& d_table);
// frees GPU resources allocated in d_table
void delete_device_table(vector_table& d_table);
// calculates ascending order of vector indexes and stores in array indexes
void radix_sort(vector_table& table, int* d_indexes, bits_t* d_level);
// stores in dest reversed bit order vectors from table
void reverse_vectors(vector_table& table, vector_table &dest);
// stores into array diffs value: diffs[0] = 0, diffs[i+1] = first high-order bit, which is different in table[permutations[i+1]] and table[permutations[i]]
void calculate_first_differences_with_previous_vector(vector_table& table, int* d_permutation, int* d_diffs);
// stores to clear_indexes values clear_indexes[i] = i; stores into abstract_classes abstract classes indexes of relation "i-1 first bits are the same", into rev_abstract_classes abstract classes of relation "n-i first bits are the same"; class_pairs[i]=abstract_classes[i]<<32 & rev_abstract_classes[i]
void calculate_abstract_classes_and_clear_indexes(vector_table& table, vector_table& rev_table, int* indexes, int* rev_indexes, int* diffs, int* rev_diffs, int* abstract_classes, int* rev_abstract_classes, unsigned long long* class_pairs, int bit, int* clear_indexes);
// finds duplicates in class_pairs and corresponding indexes in clear_indexes and (if print_results set) stores found pairs of duplicates in array result: if no result, result[i]=0, else, if found pair is (a,b), a < b, result[a*bit_count + (bit where a and b differ)] = b
void get_hamming_one_vector_pairs_to_array(vector_table& table, int* clear_indexes, int* result, unsigned long long* class_pairs, int bit);

// functions used to print array stored on a device; used mostly for debugging
void print_device_array(int* d_table, int n);
void print_device_array(unsigned long long* d_table, int n);

#endif // !VECTOR_OPERATIONS_CUH