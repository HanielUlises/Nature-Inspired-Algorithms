#include "EvolutionaryAlgorithm.h"

int main (){
    int size = 4; // e.g. 4x4 square
    int population_size = 600;
    int generations = 5000;
    //  Cross rate
    double cross_rate = 1.0f;
    // Mutation rate
    double mut_rate = 1.0f;
    EvolutionaryAlgorithm test (size, population_size, generations, cross_rate, mut_rate);
    test.solve();

    return 0;
}