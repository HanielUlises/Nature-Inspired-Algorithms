#include <iostream>

#include "EvolutionaryAlgorithm.h"

int main (){
    int size;
    std::cout << "Size: ";
    std::cin >> size;
    std::cout << std::endl;

    int population_size = 1000;
    int generations = 5000;
    //  Cross rate
    double cross_rate = 1.0f;
    // Mutation rate
    double mut_rate = 1.0f;

    EvolutionaryAlgorithm test (size, population_size, generations, cross_rate, mut_rate);
    test.solve();

    return 0;
}