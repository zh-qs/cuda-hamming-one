#include <iostream>
#include <cstdlib>
#include <random>

using namespace std;

#define BITS_IN_BYTE 8

const char* bit_rep[16] = {
    "0000", "0001", "0010", "0011",
    "0100", "0101", "0110", "0111",
    "1000", "1001", "1010", "1011",
    "1100", "1101", "1110", "1111",
};

void print_bits(random_device::result_type bits)
{
    printf("%s%s%s%s%s%s%s%s", bit_rep[bits >> 28], bit_rep[(bits >> 24) & 0x0F], bit_rep[(bits >> 20) & 0x0F], bit_rep[(bits >> 16) & 0x0F], bit_rep[(bits >> 12) & 0x0F], bit_rep[(bits >> 8) & 0x0F], bit_rep[(bits >> 4) & 0x0F], bit_rep[(bits) & 0x0F]);
}

void print_bits(random_device::result_type bits, int count)
{
    for (int i = 0; i < count; ++i)
    {
        cout << ((bits >> i) & 1);
    }
}

int main(int argc, char ** argv)
{
    const int bits_in_rd = BITS_IN_BYTE * sizeof(random_device::result_type);
    if (argc != 3)
    {
        cout << "USAGE: " << argv[0] << " n l" << endl;
        return 1;
    }
    int n = atoi(argv[1]);
    int l = atoi(argv[2]);
    random_device rd;

    cout << n << " " << l << endl;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < l / bits_in_rd; ++j)
        {
            print_bits(rd());
        }
        print_bits(rd(), l % bits_in_rd);
        cout << endl;
    }
    return 0;
}
