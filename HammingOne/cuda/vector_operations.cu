#include "vector_operations.cuh"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <thrust/scan.h>
#include <thrust/iterator/reverse_iterator.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

#include <cmath>

#include <climits>

__global__ void fill_up_down_arrays(vector_table table, int bit, int *down, int *up, int *bits, int *indexes, bool first)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int level = bit / BITS_IN_LEVEL;
	int bit_index = bit % BITS_IN_LEVEL;
	
	if (index < table.n)
	{
		if (first)
		{
			indexes[index] = index;
		}
		bits[index] = table.levels[level][indexes[index]] & (1 << bit_index);
	}
}

__global__ void fill_level(vector_table table, int level_index, bits_t* level, int* indexes)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;

	if (index < table.n)
	{
		level[index] = table.levels[level_index][indexes[index]];
	}
}

__global__ void init_indexes(int n, int* indexes)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;

	if (index < n)
	{
		indexes[index] = index;
	}
}

void radix_sort(vector_table& table, int* d_indexes, bits_t* d_level)
{
	int block_size = 1024;
	int blocks = ceil((float)table.n / block_size);
	init_indexes<<<blocks, block_size>>>(table.n, d_indexes);
	for (int i = table.l - 1; i >= 0; --i)
	{
		fill_level<<<blocks, block_size>>>(table, i, d_level, d_indexes);
		// thrust radix-sorts integers
		thrust::sort_by_key(thrust::device, d_level, d_level + table.n, d_indexes);
	}
}

__global__ void reverse_to(vector_table table, vector_table to)
{
	int index = threadIdx.x + blockDim.x * blockIdx.x;
	
	if (index < table.n)
	{
		int shift = table.l * BITS_IN_LEVEL - table.bit_count;

		bits_t rev = 0;
		for (int i = table.l - 1; i >= 0; --i)
		{
			bits_t bits = table.levels[i][index];
			for (int j = shift; j < BITS_IN_LEVEL; ++j)
			{
				rev <<= 1;
				rev |= (bits & (1 << j)) >> j;
				
			}
			if (i > 0)
			{
				bits = table.levels[i - 1][index];
				for (int j = 0; j < shift; ++j)
				{
					rev <<= 1;
					rev |= (bits & (1 << j)) >> j;
					
				}
			}
			else
			{
				rev <<= shift;
			}
			to.levels[table.l - i - 1][index] = rev;
		}
	}
}

void reverse_vectors(vector_table& table, vector_table& dest)
{
	int block_size = 1024;
	int blocks = ceil((float)table.n / block_size);
	reverse_to<<<blocks, block_size>>>(table, dest);
}

__global__ void first_differences(vector_table table, int* perm, int* diffs)
{
	int index = threadIdx.x + blockDim.x * blockIdx.x;

	if (index == 0) diffs[0] = -1;
	else if (index < table.n)
	{
		int num_zeros_at_beginning = 0;
		for (int i = 0; i < table.l; ++i)
		{
			bits_t bits = table.levels[i][perm[index]],
				bits_prev = table.levels[i][perm[index - 1]];
			int clz = __clz(bits_prev ^ bits);
			num_zeros_at_beginning += clz;
			if (clz < BITS_IN_LEVEL) break;
		}
		diffs[index] = num_zeros_at_beginning;
	}
}

void calculate_first_differences_with_previous_vector(vector_table& table, int* d_permutation, int *d_diffs)
{
	int block_size = 1024;
	int blocks = ceil((float)table.n / block_size);

	first_differences<<<blocks, block_size>>>(table, d_permutation, d_diffs);
}

void alloc_and_copy_to_device(vector_table& h_table, vector_table& d_table)
{
	d_table.n = h_table.n;
	d_table.l = h_table.l;
	d_table.bit_count = h_table.bit_count;
	cudaMalloc(&d_table.first_differences, sizeof(int) * d_table.n);
	bits_t* bits;
	cudaMalloc(&bits, sizeof(bits_t) * d_table.n * d_table.l);
	bits_t** levels = new bits_t*[d_table.l];
	cudaMalloc(&d_table.levels, sizeof(bits_t*) * d_table.l);
	for (int i = 0; i < d_table.l; ++i)
	{
		levels[i] = bits + i * d_table.n;
		
	}
	d_table.table = bits;
	cudaMemcpy(levels[0], h_table.levels[0], sizeof(bits_t) * d_table.n * d_table.l, cudaMemcpyHostToDevice);
	cudaMemcpy(d_table.levels, levels, sizeof(bits_t*) * d_table.l, cudaMemcpyHostToDevice);
}

void alloc_only(vector_table& h_table, vector_table& d_table)
{
	d_table.n = h_table.n;
	d_table.l = h_table.l;
	d_table.bit_count = h_table.bit_count;
	cudaMalloc(&d_table.first_differences, sizeof(int) * d_table.n);
	bits_t* bits;
	cudaMalloc(&bits, sizeof(bits_t) * d_table.n * d_table.l);
	bits_t** levels = new bits_t * [d_table.l];
	cudaMalloc(&d_table.levels, sizeof(bits_t*) * d_table.l);
	for (int i = 0; i < d_table.l; ++i)
	{
		levels[i] = bits + i * d_table.n;
	}
	d_table.table = bits;
	cudaMemcpy(d_table.levels, levels, sizeof(bits_t*) * d_table.l, cudaMemcpyHostToDevice);
}

void copy_to_host(vector_table& h_table, vector_table& d_table)
{
	bits_t* level;
	cudaMemcpy(&level, d_table.levels, sizeof(bits_t*), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_table.levels[0], level, sizeof(bits_t) * d_table.n * d_table.l, cudaMemcpyDeviceToHost);
}

void delete_device_table(vector_table& d_table)
{
	bits_t* level;
	cudaMemcpy(&level, d_table.levels, sizeof(bits_t*), cudaMemcpyDeviceToHost);
	cudaFree(level);
	cudaFree(d_table.levels);
}

__global__ void beginnings_of_classes(int n, int bit_count, int* diffs, int* rev_diffs, int* abstract_classes, int* rev_abstract_classes, int bit)
{
	int index = threadIdx.x + blockDim.x * blockIdx.x;

	if (index == 0) rev_abstract_classes[0] = abstract_classes[0] = 1;
	else if (index < n)
	{
		abstract_classes[index] = diffs[index] < bit ? 1 : 0;
		rev_abstract_classes[index] = rev_diffs[index] < bit_count - bit - 1 ? 1 : 0;
	}
}

__global__ void create_class_pairs_and_indexes(int n, int* indexes, int* rev_indexes, int* abstract_classes, int* rev_abstract_classes, unsigned long long* class_pairs, int* clear_indexes)
{
	int index = threadIdx.x + blockDim.x * blockIdx.x;

	if (index < n)
	{
		((int*)(&class_pairs[indexes[index]]))[0] = abstract_classes[index];
		((int*)(&class_pairs[rev_indexes[index]]))[1] = rev_abstract_classes[index];
		clear_indexes[index] = index;
	}
}

void calculate_abstract_classes_and_clear_indexes(vector_table& table, vector_table& rev_table, int* indexes, int* rev_indexes, int* diffs, int *rev_diffs, int* abstract_classes, int* rev_abstract_classes, unsigned long long *class_pairs, int bit, int *clear_indexes)
{
	int block_size = 1024;
	int blocks = ceil((float)table.n / block_size);

	beginnings_of_classes<<<blocks, block_size>>>(table.n, table.bit_count, diffs, rev_diffs, abstract_classes, rev_abstract_classes, bit);

	thrust::inclusive_scan(thrust::device, abstract_classes, abstract_classes + table.n, abstract_classes);
	thrust::inclusive_scan(thrust::device, rev_abstract_classes, rev_abstract_classes + table.n, rev_abstract_classes);

	create_class_pairs_and_indexes<<<blocks, block_size>>>(table.n, indexes, rev_indexes, abstract_classes, rev_abstract_classes, class_pairs, clear_indexes);
}

__global__ void write_results_to_array(int n, int* clear_indexes, int* result, unsigned long long* class_pairs, int bit, int bit_count)
{
	int index = threadIdx.x + blockDim.x * blockIdx.x;

	if (index < n && index > 0)
	{
		if (class_pairs[index] == class_pairs[index - 1])
		{
			int c1 = clear_indexes[index], c2 = clear_indexes[index - 1];
			int cmin = min(c1, c2), cmax = max(c1, c2);
			result[cmin * bit_count + (bit_count - bit - 1)] = cmax + 1; 
			// we number vectors in output from 1;
			// there can be also (bit) instead of (bit_count-bit-1), but this solution gives results in the same order as my brute-force program, checking is easier
		}
	}
}

void get_hamming_one_vector_pairs_to_array(vector_table& table, int* clear_indexes, int* result, unsigned long long* class_pairs, int bit)
{
	// thrust radix-sorts long long integers
	thrust::sort_by_key(thrust::device, class_pairs, class_pairs + table.n, clear_indexes);

	int block_size = 1024;
	int blocks = ceil((float)table.n / block_size);

	write_results_to_array<<<blocks, block_size>>>(table.n, clear_indexes, result, class_pairs, bit, table.bit_count);
}



void print_device_array(int* d_table, int n)
{
	int* h_table = new int[n];
	cudaMemcpy(h_table, d_table, n*sizeof(int), cudaMemcpyDeviceToHost);
	for (int i = 0; i < n; ++i)
	{
		printf("%d ", h_table[i]);
	}
	printf("\n");
	delete[] h_table;
}

void print_device_array(unsigned long long* d_table, int n)
{
	unsigned long long* h_table = new unsigned long long[n];
	cudaMemcpy(h_table, d_table, n * sizeof(unsigned long long), cudaMemcpyDeviceToHost);
	for (int i = 0; i < n; ++i)
	{
		printf("(%lld,%lld) ", h_table[i] >> 32, (h_table[i] << 32) >> 32 );
	}
	printf("\n");
	delete[] h_table;
}
 
