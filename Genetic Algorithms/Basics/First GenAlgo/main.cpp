#include <cmath>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

#include "encode.h"
#include "decode.h"
#include "GeneticAlgorithm.h"

int countDecimalPlaces(double number) {
    // Like stod in C but inverted jiji
    std::string numberAsString = std::to_string(number);
    // Trim trailing zeros
    numberAsString.erase(numberAsString.find_last_not_of('0') + 1, std::string::npos); 
    auto decimalPos = numberAsString.find('.');

     // No decimal point found
    if (decimalPos == std::string::npos) 
        return 0;

    return numberAsString.length() - decimalPos - 1; 
}

int bitsNeeded(double x, double y, double& range) {
    // Determine the maximum number of decimal places between the two numbers
    int maxDecimals = std::max(countDecimalPlaces(x), countDecimalPlaces(y));
    
    // Formulae implemented
    // log2(upLim * 10^dec - lowLim * 10^dec)
    double scaledX = x * std::pow(10, maxDecimals);
    double scaledY = y * std::pow(10, maxDecimals);
    
    // Difference between the scaled numbers
    range = std::abs(scaledY - scaledX);
    return std::ceil(std::log2(range));
}

void conversion (double lowerLimit, double range){
    int iter = (int) range;
    int aux = 0;
    std::vector<int> result_it;
    std::string bin_string;

    for(int i = 0; i <= iter; i++){
        decToBin (i, bin_string);
        std::cout << lowerLimit << " : ";
        std::cout << bin_string << std::endl;
        std::cout << std::endl;
        lowerLimit += 0.01;

        aux = binToDec(bin_string);
        result_it.push_back(aux);
        aux = 0;
    }

}


int main() {
    /*
    
    double x, y;
    std::cout << "Enter the first number (x): ";
    std::cin >> x;
    std::cout << "Enter the second number (y): ";
    std::cin >> y;

    double range = 0;

    int totalBitsNeeded = bitsNeeded(x, y, range);
    conversion (x,range);
    std::cout << "Total number of bits needed: " << totalBitsNeeded << std::endl;
    */
   
    int population = 100;
    int num_gen = 500;
    double cross_rate = 0.9f;
    double mut_rate = 0.1f;

    GeneticAlgorithm newGen (population, num_gen, cross_rate, mut_rate);
    //option 1: Rosenbrock
    //option 2: Ackley
    newGen.run(1);


    return 0;
}
