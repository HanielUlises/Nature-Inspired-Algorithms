// EvolutionaryAlgorithm.h
#ifndef EVOLUTIONARY_ALGORITHM_H
#define EVOLUTIONARY_ALGORITHM_H

#include <vector>
#include <random>
struct IndividualWithFitness {
    std::vector<int> individual;
    double fitness;

    IndividualWithFitness(std::vector<int> ind, double fit)
        : individual(std::move(ind)), fitness(fit) {}
};

class EvolutionaryAlgorithm {
public:
    EvolutionaryAlgorithm(int size, int populationSize, int generations, double cross_rate, double mut_rate);
    void solve();

private:
    int size; // Size of the magic square
    int populationSize; // Number of individuals in the population
    int generations; // Number of generations to evolve
    double cross_rate;
    double mut_rate;
    std::vector<std::vector<int>> population; // Population of solutions
    std::uniform_real_distribution<double> disCross; // Cross register

    void initializePopulation();
    std::vector<int> generateIndividual();
    int calculateFitness(const std::vector<int>& individual);
    std::vector<int> selection(std::vector<int>fitnessP);
    std::vector<int> tournamentSelection(std::vector<int>fitnessValues,int tournamentSize);
    std::vector<int> crossoverOX(std::vector<int>& parent1, std::vector<int>& parent2);
    
    std::vector<std::vector<int>> crossoverPMX(const std::vector<int>& parent1,const  std::vector<int>& parent2);
    std::vector<int> mutationIn(std::vector<int>& individual);
    std::vector<int> mutationDes(std::vector<int> individual);
    std::vector<std::vector<int>> elitism(std::vector<std::vector<int>> population);
    bool isMagicSquare(const std::vector<int>& square);
    void printSolution(const std::vector<int>& solution);
};

#endif // EVOLUTIONARY_ALGORITHM_H
