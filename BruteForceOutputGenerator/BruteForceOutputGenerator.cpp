#include <iostream>
#include <fstream>
#include <string>

using namespace std;

bool areHammingOne(string s1, string s2, int l)
{
    int d = 0;
    for (int i = 0; i < l; ++i)
    {
        if (s1[i] != s2[i])
        {
            if (++d > 1) return false;
        }
    }
    return d == 1;
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        cout << "USAGE: %s input_file" << endl;
        return 1;
    }
    ifstream s(argv[1]);

    int n, l;
    string* strings;
    if (s.good())
    {
        s >> n >> l;
        strings = new string[n];
        for (int i = 0; i < n; ++i)
        {
            s >> strings[i];
        }
    }
    else return 1;

    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            if (areHammingOne(strings[i], strings[j], l))
                cout << i+1 << ";" << j+1 << endl;
        }
    }

    return 0;
}