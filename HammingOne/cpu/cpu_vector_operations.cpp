#include "cpu_vector_operations.h"
#include <iostream>

#define BITS_RADIX_SORT 8

struct node
{
	int index;
	node* next = nullptr;

	node(int index) : index(index) { }
};

void cpu_radix_sort_and_alloc_indexes(vector_table& table, int** indexes)
{
	*indexes = new int[table.n];
	node* buckets_heads[1 << BITS_RADIX_SORT];
	node* buckets_tails[1 << BITS_RADIX_SORT];

	for (int i = 0; i < 1 << BITS_RADIX_SORT; ++i)
		buckets_heads[i] = nullptr;

	for (int i = 0; i < table.n; ++i)
	{
		(*indexes)[i] = i;
	}
	for (int j = table.l - 1; j >= 0; --j)
	{
		for (int k = 0; k < BITS_IN_LEVEL / BITS_RADIX_SORT; ++k)
		{
			for (int i = 0; i < table.n; ++i)
			{
				int l = (table.levels[j][(*indexes)[i]] >> (k * BITS_RADIX_SORT)) & 0xff;
				if (buckets_heads[l] == nullptr)
				{
					buckets_heads[l] = buckets_tails[l] = new node((*indexes)[i]);
				}
				else
				{
					buckets_tails[l]->next = new node((*indexes)[i]);
					buckets_tails[l] = buckets_tails[l]->next;
				}
			}
			int index = 0;
			for (int i = 0; i < (1 << BITS_RADIX_SORT); ++i)
			{
				while (buckets_heads[i])
				{
					(*indexes)[index++] = buckets_heads[i]->index;
					node* head = buckets_heads[i];
					buckets_heads[i] = buckets_heads[i]->next;
					delete head;
				}
			}
		}
	}

}

void cpu_reverse_vectors(vector_table& table, vector_table& dest)
{
	for (int i = 0; i < table.n; ++i)
	{
		int shift = table.l * BITS_IN_LEVEL - table.bit_count;

		bits_t rev = 0;
		for (int j = table.l - 1; j >= 0; --j)
		{
			bits_t bits = table.levels[j][i];
			for (int k = shift; k < BITS_IN_LEVEL; ++k)
			{
				rev <<= 1;
				rev |= (bits & (1 << k)) >> k;

			}
			if (j > 0)
			{
				bits = table.levels[j - 1][i];
				for (int k = 0; k < shift; ++k)
				{
					rev <<= 1;
					rev |= (bits & (1 << k)) >> k;

				}
			}
			else
			{
				rev <<= shift;
			}
			dest.levels[table.l - j - 1][i] = rev;
		}
	}
}

int number_of_high_order_zeros(bits_t bits)
{
	int res = 0;
	for (int i = BITS_IN_LEVEL - 1; i >= 0; --i)
	{
		if ((bits >> i) & 1) break;
		++res;
	}
	return res;
}

void cpu_calculate_first_differences_with_previous_vector(vector_table& table, int* permutation, int* diffs)
{
	diffs[0] = -1;
	for (int index = 1; index < table.n; ++index)
	{
		int num_zeros_at_beginning = 0;
		for (int i = 0; i < table.l; ++i)
		{
			bits_t bits = table.levels[i][permutation[index]],
				bits_prev = table.levels[i][permutation[index - 1]];
			int clz = number_of_high_order_zeros(bits_prev ^ bits);
			num_zeros_at_beginning += clz;
			if (clz < BITS_IN_LEVEL) break;
		}
		diffs[index] = num_zeros_at_beginning;
	}
}

void cpu_calculate_abstract_classes_and_clear_indexes(vector_table& table, vector_table& rev_table, int* indexes, int* rev_indexes, int* diffs, int* rev_diffs, int* abstract_classes, int* rev_abstract_classes, unsigned long long* class_pairs, int bit, int* clear_indexes)
{
	abstract_classes[0] = rev_abstract_classes[0] = 1;
	clear_indexes[0] = 0;

	((int*)(&class_pairs[indexes[0]]))[0] = 1;
	((int*)(&class_pairs[rev_indexes[0]]))[1] = 1;

	for (int i = 1; i < table.n; ++i)
	{
		abstract_classes[i] = abstract_classes[i - 1] + (diffs[i] < bit ? 1 : 0);
		rev_abstract_classes[i] = rev_abstract_classes[i - 1] + (rev_diffs[i] < table.bit_count - bit - 1 ? 1 : 0);
		clear_indexes[i] = i;

		((int*)(&class_pairs[indexes[i]]))[0] = abstract_classes[i];
		((int*)(&class_pairs[rev_indexes[i]]))[1] = rev_abstract_classes[i];
	}
}

// inspired by https://www.sanfoundry.com/c-program-implement-radix-sort/
void radix_sort_ull_by_key(unsigned long long* tab, int n, int* keys)
{
	unsigned long long* output = new unsigned long long[n];
	int* key_output = new int[n];
	for (int i = 0; i < sizeof(unsigned long long) * BITS_IN_BYTE; i += BITS_RADIX_SORT)
	{
		int count[1 << BITS_RADIX_SORT]{ 0 };

		for (int j = 0; j < n; ++j)
		{
			count[(tab[j] >> (i * BITS_RADIX_SORT)) & 0xff]++;
			output[j] = 0;
		}

		for (int j = 1; j < (1 << BITS_RADIX_SORT); ++j)
		{
			count[j] += count[j - 1];
		}

		for (int j = n - 1; j >= 0; j--) 
		{
			output[count[(tab[j] >> (i * BITS_RADIX_SORT)) & 0xff] - 1] = tab[j];
			key_output[count[(tab[j] >> (i * BITS_RADIX_SORT)) & 0xff] - 1] = keys[j];

			count[(tab[j] >> (i * BITS_RADIX_SORT)) & 0xff]--;
		}

		for (int j = 0; j < n; j++)
		{
			tab[j] = output[j];
			keys[j] = key_output[j];
		}
	}
	delete[] output;
	delete[] key_output;
}

void cpu_print_found_hamming_one_vector_pairs(vector_table& table, int* clear_indexes, unsigned long long* class_pairs, int bit, bool print_results)
{
	radix_sort_ull_by_key(class_pairs, table.n, clear_indexes);

	for (int i = 1; i < table.n; ++i)
	{
		if (class_pairs[i] == class_pairs[i - 1] && print_results)
		{
			if (clear_indexes[i] < clear_indexes[i - 1])
			{
				std::cout << clear_indexes[i] + 1 << ";" << clear_indexes[i - 1] + 1 << std::endl;
			}
			else
			{
				std::cout << clear_indexes[i - 1] + 1 << ";" << clear_indexes[i] + 1 << std::endl;
			}
		}
	}
}