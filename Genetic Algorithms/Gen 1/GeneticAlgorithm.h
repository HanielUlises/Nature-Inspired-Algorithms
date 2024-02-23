#pragma once

#ifndef GENETIC_ALGORITHM
#define GENETIC_ALGORITHM

#include "Solution.h"

#include <iostream>

class GeneticAlgorithm{
    public:
    GeneticAlgorithm(int populationSize);
    Solution perform (int numberOfBits, int low, int high);

    private:
    int populationSize;
};

#endif // GENETIC_ALGORITHM