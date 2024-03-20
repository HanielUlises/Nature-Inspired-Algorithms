#pragma once

#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <vector>
#include <functional>

class GeneticAlgorithm {
public:
    GeneticAlgorithm(int populationSize, int numberOfGenerations, double crossoverRate, double mutationRate);
    virtual ~GeneticAlgorithm();

    void run(int option);
    void plotConvergenceGraph();
    // Test functions
    double rosenbrockFunction(const std::vector<double>& individual);
    double ackleyFunction(const std::vector<double>& individual);

    std::vector<double> bestFitnessHistory;
    std::vector<double> worstFitnessHistory;
    std::vector<double> averageFitnessHistory;

private:
    int populationSize;
    int numberOfGenerations;
    double crossoverRate;
    double mutationRate;
    std::vector<std::vector<double>> population;
    std::vector<double> fitnessValues;
    void elitismParents();
    void initializePopulation(int option);
    void evaluateFitness(std::function<double(const std::vector<double>&)> objectiveFunction);
    std::vector<int> selection();
    std::vector<std::vector<double>> crossover(std::vector<int>& selectedParents);
    std::vector<std::vector<double>> mutation(std::vector<std::vector<double>> hijos);
    bool shouldStop(int currentGeneration);
};

#endif