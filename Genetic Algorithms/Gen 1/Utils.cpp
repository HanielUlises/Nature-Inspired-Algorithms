#include "Utils.h"

#include <cstdlib>
#include <stdexcept>
#include <set>

double randomDouble (){
    return double(rand() % 1000)/1000;
}

bool randomProbability (double probability){
    if(probability < 0 || probability > 1){
        throw std::runtime_error("Probabilities are only between 0 and 1 \n");
    }

    double r = randomDouble();

    return r < probability ? true : false;
}

double precision (int low, int high, int numberOfBits){
    double prec = (double)(high - low)/(double)(pow(2,numberOfBits)-1);
    return prec;
}

std::set<int> randomDistinctNumbers (int upperLimit, int amountNumbers){
    std::set<int> numbers;
    while (numbers.size() < amountNumbers){
        numbers.insert(rand() % upperLimit);
    }

    return numbers;
}