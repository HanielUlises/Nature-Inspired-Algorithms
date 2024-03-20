#include <string>

int binToDec(std::string line) {
    int decimalNumber = 0;
    int multiplier = 1;

    for (int i = line.length(); i > 0; i--) {
        int binaryNumber = line[i - 1] - '0';
        decimalNumber += binaryNumber * multiplier;
        multiplier *= 2;
    }

    return decimalNumber;
}