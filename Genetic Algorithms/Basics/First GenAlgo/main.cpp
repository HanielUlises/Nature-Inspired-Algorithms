#include <cmath>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

#include "GeneticAlgorithm.h"

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
    int num_gen = 100;
    double cross_rate = 0.9f;
    double mut_rate = 0.1f;

    GeneticAlgorithm newGen (population, num_gen, cross_rate, mut_rate);
    //option 1: Rosenbrock
    //option 2: Ackley
    newGen.run2(1);


    return 0;
}
