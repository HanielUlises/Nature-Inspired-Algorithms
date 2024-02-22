#include "Solution.h"
#include "Utils.h"

#include <sstream>
#include <cmath>

Solution::Solution(int numberOfBits, int low, int high) :
    mNumberOfBits(numberOfBits),
    mLow(low),
    mHigh(high)
{
    bits.reserve(numberOfBits);
    for(int i = 0; i < numberOfBits; i++){
        bits.push_back(rand() % 2);
    }
}

std::string Solution::toString() {
    std::ostringstream stream;
    stream << "[";
    for (int bit : bits) {
        stream << bit << " ";
    }
    stream.seekp(-1, std::ios_base::end);
    stream << "]" << bitsToDouble() << " fitness = " << fitness ();
    return stream.str();
}

double Solution::bitsToDouble () {
    double prec = precision(mLow, mHigh, mNumberOfBits);
    double sum = 0.0;
    for (int i = 0; i < mNumberOfBits; i++) {
        if(bits[i] == 1) {
            sum += sum + pow(2,i) * bits[i];
        }
    }
    return mLow + sum * prec;
}

// f(x) = x + 2 * sin(x)
double Solution::fitness () {
    double x = bitsToDouble ();
    return x + 2 * sin(x);
}