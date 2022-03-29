#ifndef UTIL_H
#define UTIL_H

#include <string>
#include "../structures/vector_table.h"

// converts string of size <= BITS_IN_LEVEL to integer type
bits_t bit_string_to_bits(std::string str);
// prints integer type in binary notation; mostly used for debugging
void print_bits(bits_t bits);
// prints vector_table table contents; mostly used for debugging
void print_table(vector_table& table);

#endif // !UTIL_H