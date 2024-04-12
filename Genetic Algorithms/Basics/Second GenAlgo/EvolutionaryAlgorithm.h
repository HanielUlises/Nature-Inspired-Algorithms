// EvolutionaryAlgorithm.h
#ifndef EVOLUTIONARY_ALGORITHM_H
#define EVOLUTIONARY_ALGORITHM_H

#include <vector>
struct IndividualWithFitness {
    std::vector<int> individual;
    double fitness;

    IndividualWithFitness(std::vector<int> ind, double fit)
        : individual(std::move(ind)), fitness(fit) {}
};

class EvolutionaryAlgorithm {
public:
    EvolutionaryAlgorithm(int size, int populationSize, int generations);
    void solve();

private:
    int size; // Size of the magic square
    int populationSize; // Number of individuals in the population
    int generations; // Number of generations to evolve
    std::vector<std::vector<int>> population; // Population of solutions

    void initializePopulation();
    std::vector<int> generateIndividual();
    int calculateFitness(const std::vector<int>& individual);
    std::vector<int> selection(std::vector<int>fitnessP);
    std::vector<int> tournamentSelection(std::vector<int>fitnessValues,int tournamentSize);
    std::vector<int> crossover(std::vector<int>& parent1, std::vector<int>& parent2);
    std::vector<int> mutation(std::vector<int> individual);
    std::vector<std::vector<int>> elitism(std::vector<std::vector<int>> population);
    bool isMagicSquare(const std::vector<int>& square);
    void printSolution(const std::vector<int>& solution);
};

#endif // EVOLUTIONARY_ALGORITHM_H
