#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <vector>
#include "Solution.h"

class Genetic_Algorithm {
private:
    std::vector<Solution> population;
    int populationSize;
    int numCities;
    double mutationRate;

    Solution tournamentSelection();
    void crossover(Solution& parent1, Solution& parent2, Solution& child);
    void mutate(Solution& sol);

public:
    Genetic_Algorithm(int popSize, int nCities, double mutRate);
    std::vector<Solution>& getPopulation();
    void evolve();
    Solution getBestSolution();
};

#endif