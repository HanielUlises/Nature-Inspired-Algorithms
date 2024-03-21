#include <cmath>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

#include "GeneticAlgorithm.h"

int main() {
    // Individuals
    int population = 100;
    // Generations
    int num_gen = 100;
    // Crossover rate
    double cross_rate = 0.65f;
    // Mutation rate
    double mut_rate = 0.05f;

    GeneticAlgorithm newGen (population, num_gen, cross_rate, mut_rate);
    // option 1: Rosenbrock
    // option 2: Ackley
    
    // newGen.binary_performance(1);
    // newGen.real_performance(1);
    // newGen.binary_performance(2);
    newGen.real_performance(2);


    return 0;
}
