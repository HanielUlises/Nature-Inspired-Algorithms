#include "Genetic_Algorithm.h"
#include <cstdlib>

Genetic_Algorithm::Genetic_Algorithm(int popSize, int nCities, double mutRate) {
    populationSize = popSize;
    numCities = nCities;
    mutationRate = mutRate;
    population.resize(popSize);
    for (int i = 0; i < popSize; i++) {
        population[i] = Solution(nCities);
    }
}

std::vector<Solution>& Genetic_Algorithm::getPopulation() {
    return population;
}

Solution Genetic_Algorithm::tournamentSelection() {
    int tournamentSize = 5;
    Solution best = population[rand() % populationSize];
    for (int i = 1; i < tournamentSize; i++) {
        Solution contender = population[rand() % populationSize];
        if (contender.getFitness() > best.getFitness()) {
            best = contender;
        }
    }
    return best;
}

void Genetic_Algorithm::crossover(Solution& parent1, Solution& parent2, Solution& child) {
    std::vector<int>& p1 = parent1.getRoute();
    std::vector<int>& p2 = parent2.getRoute();
    std::vector<int>& c = child.getRoute();
    
    int start = rand() % numCities;
    int end = rand() % numCities;
    if (start > end) {
        int temp = start;
        start = end;
        end = temp;
    }

    std::vector<bool> used(numCities, false);
    for (int i = start; i <= end; i++) {
        c[i] = p1[i];
        used[p1[i]] = true;
    }

    int currentPos = (end + 1) % numCities;
    for (int i = 0; i < numCities; i++) {
        if (currentPos == start) {
            currentPos = end + 1;
        }
        if (!used[p2[i]]) {
            c[currentPos] = p2[i];
            used[p2[i]] = true;
            currentPos = (currentPos + 1) % numCities;
        }
    }
}

void Genetic_Algorithm::mutate(Solution& sol) {
    if ((double)rand() / RAND_MAX < mutationRate) {
        int i = rand() % numCities;
        int j = rand() % numCities;
        sol.swapCities(i, j);
    }
}

void Genetic_Algorithm::evolve() {
    std::vector<Solution> newPopulation(populationSize, Solution(numCities));
    newPopulation[0] = getBestSolution();
    
    for (int i = 1; i < populationSize; i++) {
        Solution parent1 = tournamentSelection();
        Solution parent2 = tournamentSelection();
        crossover(parent1, parent2, newPopulation[i]);
        mutate(newPopulation[i]);
    }
    
    population = newPopulation;
}

Solution Genetic_Algorithm::getBestSolution() {
    Solution best = population[0];
    for (int i = 1; i < populationSize; i++) {
        if (population[i].getFitness() > best.getFitness()) {
            best = population[i];
        }
    }
    return best;
}