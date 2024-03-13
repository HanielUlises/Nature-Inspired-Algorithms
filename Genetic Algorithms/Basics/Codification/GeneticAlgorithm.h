#pragma once
#include <vector>
#include <functional>

class GeneticAlgorithm {
public:
    GeneticAlgorithm(int populationSize, int numberOfGenerations, double crossoverRate, double mutationRate);
    virtual ~GeneticAlgorithm();

    void run();
    // Test functions
    double rosenbrockFunction(const std::vector<double>& individual);
    double ackleyFunction(const std::vector<double>& individual);

private:
    int populationSize;
    int numberOfGenerations;
    double crossoverRate;
    double mutationRate;
    std::vector<std::vector<double>> population;
    std::vector<double> fitnessValues;

    void initializePopulation();
    void evaluateFitness(std::function<double(const std::vector<double>&)> objectiveFunction);
    std::vector<int> selection();
    void crossover(std::vector<int>& selectedParents);
    void mutation();
    bool shouldStop(int currentGeneration);
};
