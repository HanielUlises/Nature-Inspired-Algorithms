#include "Solution.h"
#include "Utils.h"

#include <sstream>
#include <cmath>

Solution::Solution(int numberOfBits, int low, int high) :
    numberOfBits{numberOfBits},
    low{low},
    high{high}
{
    bits.reserve(numberOfBits);
    for(int i = 0; i < numberOfBits; i++){
        bits.push_back(rand() % 2);
    }
}

Solution::Solution(int numberOfBits, int low, int high, std::vector<int> bits){
    this-> numberOfBits = numberOfBits;
    this-> low = low;
    this-> high = high;
    this-> bits = bits;
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
    double prec = precision(low, high, numberOfBits);
    double sum = 0.0;
    for (int i = 0; i < numberOfBits; i++) {
        if(bits[i] == 1) {
            sum += sum + pow(2,i) * bits[i];
        }
    }
    return low + sum * prec;
}

// f(x) = x + 2 * sin(x)
double Solution::fitness () {
    double x = bitsToDouble ();
    return x + 2 * sin(x);
}

std::vector<Solution> Solution::single_point_crossover (Solution other, double crossoverProbability){
    bool cross = randomProbability(crossoverProbability);

    if(cross){
        int crossPoint = rand() % numberOfBits;


        std::vector<int> bits1;
        std::vector<int> bits2;

        std::copy(this-> bits.begin(), this -> bits.begin() + crossPoint, std::back_inserter(bits1));
        std::copy(other.bits.begin() + crossPoint, other.bits.end(), std::back_inserter(bits1));

        std::copy(other.bits.begin(), other.bits.begin() + crossPoint, std::back_inserter(bits2));
        std::copy(this -> bits.begin() + crossPoint, this-> bits.end(), std::back_inserter(bits2));

        Solution child1{numberOfBits, low, high, bits1};
        Solution child2{numberOfBits, low, high, bits2};

        return {child1, child2};
    }else{
        return {*this, other};
    }
}