#ifndef VECTOR_TABLE_H
#define VECTOR_TABLE_H

#define BITS_IN_BYTE 8
#define BITS_IN_LEVEL (BITS_IN_BYTE * sizeof(bits_t))

typedef unsigned int bits_t;

// main structure to contain vectors
// vector are aligned by levels: whole vector is represented by arrays: levels[0][i], level[1][i], ..., level[l-1][i].
struct vector_table
{
	int n;
	int l;
	int bit_count;
	bits_t** levels;
	bits_t* table;
	int* first_differences;
};

#endif // ! VECTOR_TABLE_H

