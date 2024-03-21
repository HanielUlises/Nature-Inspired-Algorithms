#pragma once

#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <vector>
#include <functional>
#include <string>

struct IndividualWithFitness {
    std::vector<std::string> individual;
    double fitness;

    IndividualWithFitness(std::vector<std::string> ind, double fit)
        : individual(std::move(ind)), fitness(fit) {}
};

class GeneticAlgorithm {
public:
    GeneticAlgorithm(int populationSize, int numberOfGenerations, double crossoverRate, double mutationRate);
    virtual ~GeneticAlgorithm();

    void real_performance(int option);
    void binary_performance(int option);
    void plotConvergenceGraph(std::string function);
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
    std::vector<std::vector<std::string>> populationBinary;
    std::vector<double> fitnessValues;

    void elitismParents();
    void elitismParentsBinary();

    void initializePopulation(int option);
    void initializePopulationBinary(int option);
    void evaluateFitness(std::function<double(const std::vector<double>&)> objectiveFunction, int option,  double liminf, double limsup);

    std::vector<int> selection();
    std::vector<int> tournamentSelection(int tournamentSize);

    std::vector<std::vector<double>> crossover(std::vector<int>& selectedParents);
    std::vector<std::vector<std::string>> crossover_binary(std::vector<int>& selectedParents);

    std::vector<std::vector<double>> mutation(std::vector<std::vector<double>> hijos);
    std::vector<std::vector<std::string>> mutation_binary(std::vector<std::vector<std::string>> hijos, int currentGeneration);
    bool shouldStop(int currentGeneration, const int option);

    double adaptiveMutationRate(int currentGeneration);
};

#endif