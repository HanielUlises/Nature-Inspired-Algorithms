#include "GeneticAlgorithm.h"
#include <algorithm>
#include <cmath>
#include <random>

// Random number generator initialization
std::random_device rd;
std::mt19937 gen(rd());

GeneticAlgorithm::GeneticAlgorithm(int populationSize, int numberOfGenerations, double crossoverRate, double mutationRate)
    : populationSize(populationSize), numberOfGenerations(numberOfGenerations),
      crossoverRate(crossoverRate), mutationRate(mutationRate) {
}

GeneticAlgorithm::~GeneticAlgorithm() {
}

void GeneticAlgorithm::run() {
    initializePopulation();
        for (int gen = 0; gen < numberOfGenerations; ++gen) {
        // Functions test
        evaluateFitness(std::bind(&GeneticAlgorithm::rosenbrockFunction, this, std::placeholders::_1));
        // evaluateFitness(std::bind(&GeneticAlgorithm::ackleyFunction, this, std::placeholders::_1));
        auto selectedParents = selection();
        crossover(selectedParents);
        mutation();
        if (shouldStop(gen)) break;
    }
}

void GeneticAlgorithm::initializePopulation() {
}

void GeneticAlgorithm::evaluateFitness(std::function<double(const std::vector<double>&)> objectiveFunction) {
    for (auto& individual : population) {
        double fitness = objectiveFunction(individual);
        fitnessValues.push_back(fitness);
    }
}

std::vector<int> GeneticAlgorithm::selection() {
    return std::vector<int>(); // Placeholder
}

void GeneticAlgorithm::crossover(std::vector<int>& selectedParents) {
}

void GeneticAlgorithm::mutation() {
}

bool GeneticAlgorithm::shouldStop(int currentGeneration) {
    return false; // Placeholder
}

double GeneticAlgorithm::rosenbrockFunction(const std::vector<double>& individual) {
    double rosenbrock = 0.0;
    for (size_t i = 0; i < individual.size() - 1; ++i) {
        rosenbrock += 100 * std::pow((individual[i+1] - std::pow(individual[i], 2)), 2) + std::pow((1 - individual[i]), 2);
    }
    return rosenbrock;
}

double GeneticAlgorithm::ackleyFunction(const std::vector<double>& individual) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (auto x_i : individual) {
        sum1 += std::pow(x_i, 2);
        sum2 += std::cos(2 * M_PI * x_i);
    }
    double ackley = -20 * std::exp(-0.2 * std::sqrt(sum1 / individual.size())) - std::exp(sum2 / individual.size()) + 20 + M_E;
    return ackley;
}
