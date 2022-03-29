#ifndef CPU_VECTOR_OPERATIONS_CUH
#define CPU_VECTOR_OPERATIONS_CUH

#include "../structures/vector_table.h"

// calculates ascending order of vector indexes and stores in inside-function-allocated array indexes
void cpu_radix_sort_and_alloc_indexes(vector_table& table, int** indexes);
// stores in dest reversed bit order vectors from table
void cpu_reverse_vectors(vector_table& table, vector_table& dest);
// stores into array diffs value: diffs[0] = 0, diffs[i+1] = first high-order bit, which is different in table[permutations[i+1]] and table[permutations[i]]
void cpu_calculate_first_differences_with_previous_vector(vector_table& table, int* permutation, int* diffs);
// stores to clear_indexes values clear_indexes[i] = i; stores into abstract_classes abstract classes indexes of relation "i-1 first bits are the same", into rev_abstract_classes abstract classes of relation "n-i first bits are the same"; class_pairs[i]=abstract_classes[i]<<32 & rev_abstract_classes[i]
void cpu_calculate_abstract_classes_and_clear_indexes(vector_table& table, vector_table& rev_table, int* indexes, int* rev_indexes, int* diffs, int* rev_diffs, int* abstract_classes, int* rev_abstract_classes, unsigned long long* class_pairs, int bit, int* clear_indexes);
// finds duplicates in class_pairs and corresponding indexes in clear_indexes and (if print_results set) prints found pairs of duplicates
void cpu_print_found_hamming_one_vector_pairs(vector_table& table, int* clear_indexes, unsigned long long* class_pairs, int bit, bool print_results = true);

#endif // !CPU_VECTOR_OPERATIONS_CUH