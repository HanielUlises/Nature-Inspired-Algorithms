#include <string>

void decToBin(int decimalNumber, std::string& binaryNumber) {
    binaryNumber = (decimalNumber == 0) ? "0" : "";

    while (decimalNumber > 0) {
        binaryNumber = std::to_string(decimalNumber % 2) + binaryNumber;
        decimalNumber /= 2;
    }
}

