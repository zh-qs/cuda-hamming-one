#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string>

// https://stackoverflow.com/questions/478898/how-do-i-execute-a-command-and-get-the-output-of-the-command-within-c-using-po
std::string exec(const char* cmd) {
    char buffer[128];
    std::string result = "";
    FILE* pipe = _popen(cmd, "r");
    if (!pipe) throw std::runtime_error("popen() failed!");
    try {
        while (fgets(buffer, sizeof buffer, pipe) != NULL) {
            result += buffer;
        }
    }
    catch (...) {
        _pclose(pipe);
        throw;
    }
    _pclose(pipe);
    return result;
}

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cout << "USAGE: %s test_file" << std::endl;
        return 2;
    }
    std::string ham = exec((std::string("HammingOne ") + std::string(argv[1])).c_str());
    std::string br = exec((std::string("BruteForceOutputGenerator ") + std::string(argv[1])).c_str());
    if (ham == br)
    {
        std::cout << "EQUAL!" << std::endl;
        return 0;
    }
    std::cout << "NOT EQUAL!" << std::endl;
    return 1;
}
