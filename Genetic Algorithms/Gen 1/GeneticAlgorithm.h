#pragma once

#ifndef GENETIC_ALGORITHM
#define GENETIC_ALGORITHM

#include "Solution.h"

#include <iostream>
#include <ctime>
#include <set>

class GeneticAlgorithm{
    public:
    GeneticAlgorithm(int populationSize, int generation, int tournamentGroupSize);
    Solution perform (int numberOfBits, int low, int high);
    std::vector<Solution> tournamentWinners (std::vector<Solution>const & currentGeneration);
    std::vector<Solution> tournamentCrossover (std::vector<Solution>const & currentGeneration);

    private:
    int populationSize;
    int generations;
    int tournamentGroupSize;
};

#endif // GENETIC_ALGORITHM