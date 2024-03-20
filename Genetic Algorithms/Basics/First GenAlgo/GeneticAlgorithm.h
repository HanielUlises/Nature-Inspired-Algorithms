#pragma once

#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <vector>
#include <functional>
#include <string>

class GeneticAlgorithm {
public:
    GeneticAlgorithm(int populationSize, int numberOfGenerations, double crossoverRate, double mutationRate);
    virtual ~GeneticAlgorithm();

    void run(int option);
    void run2(int option);
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
    std::vector<std::vector<std::string>> populationBinary;
    std::vector<double> fitnessValues;
    void elitismParents();
    void elitismParents2();
    void initializePopulation(int option);
    void initializePopulationBinary(int option);
    void evaluateFitness(std::function<double(const std::vector<double>&)> objectiveFunction, int option,  double liminf, double limsup);
    std::vector<int> selection();
    std::vector<std::vector<double>> crossover(std::vector<int>& selectedParents);
    std::vector<std::vector<std::string>> crossoverB(std::vector<int>& selectedParents);
    std::vector<std::vector<double>> mutation(std::vector<std::vector<double>> hijos);
    std::vector<std::vector<std::string>> mutation2(std::vector<std::vector<std::string>> hijos);
    bool shouldStop(int currentGeneration);
};

#endif