#include "EvolutionaryAlgorithm.h"

int main (){
    int size = 4; // e.g. 3x3 square
    int population_size = 100;
    int generations = 5000;

    EvolutionaryAlgorithm test (size, population_size, generations);
    test.solve();

    return 0;
}