#pragma once

#ifndef GENETIC_ALGORITHM
#define GENETIC_ALGORITHM

#include "Solution.h"

#include <iostream>
#include <ctime>
#include <set>

class GeneticAlgorithm{
    public:
    GeneticAlgorithm(int populationSize, int generation, int tournamentGroupSize, double crossoverProbability);
    Solution perform (int numberOfBits, int low, int high);

    private:
    Solution tournamentWinners (std::vector<Solution>const & currentGeneration);
    std::vector<Solution> tournamentCrossover (std::vector<Solution>const & currentGeneration);

    int populationSize;
    int generations;
    int tournamentGroupSize;
    double crossoverProbability;
};

#endif // GENETIC_ALGORITHM