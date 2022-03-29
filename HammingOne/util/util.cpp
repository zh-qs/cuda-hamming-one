#include "util.h"

bits_t bit_string_to_bits(std::string str)
{
	// we assume that str's size is <=BITS_IN_LEVEL (default 32)
	// if str's size is <BITS_IN_LEVEL, bits are written at the beginning of the number (i.e: "111" => 111000..00)

	bits_t bits = 0;
	for (int i = 0; i < str.size(); ++i)
	{
		bits |= (str[i] == '1' ? 1 : 0) << (BITS_IN_LEVEL - i - 1);
	}
	return bits;
}

const char* bit_rep[16] = {
	"0000", "0001", "0010", "0011",
	"0100", "0101", "0110", "0111",
	"1000", "1001", "1010", "1011",
	"1100", "1101", "1110", "1111",
};

void print_bits(bits_t bits)
{
	printf("%s%s%s%s%s%s%s%s", bit_rep[bits >> 28], bit_rep[(bits >> 24) & 0x0F], bit_rep[(bits >> 20) & 0x0F], bit_rep[(bits >> 16) & 0x0F], bit_rep[(bits >> 12) & 0x0F], bit_rep[(bits >> 8) & 0x0F], bit_rep[(bits >> 4) & 0x0F], bit_rep[(bits) & 0x0F]);
}

void print_table(vector_table& table)
{
	printf("n: %d, l: %d, bits: %d\n", table.n, table.l, table.bit_count);
	for (int i = 0; i < table.n; ++i)
	{
		for (int j = 0; j < table.l; ++j)
		{
			print_bits(table.levels[j][i]);
		}
		printf("\n");
	}
}