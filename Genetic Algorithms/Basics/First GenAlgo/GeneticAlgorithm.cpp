#include "GeneticAlgorithm.h"
#include <algorithm>
#include <cmath>
#include <random>

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
    population.clear();

    // CHECK
    std::uniform_real_distribution<double> dis(-10.0, 10.0);

    int numGenes = 100;

    for (int i = 0; i < populationSize; ++i) {
        std::vector<double> individual;
        for (int j = 0; j < numGenes; ++j) {
            double gene = dis(gen);
            individual.push_back(gene);
        }
        population.push_back(individual);
    }
}


void GeneticAlgorithm::evaluateFitness(std::function<double(const std::vector<double>&)> objectiveFunction) {
    for (auto& individual : population) {
        double fitness = objectiveFunction(individual);
        fitnessValues.push_back(fitness);
    }
}

std::vector<int> GeneticAlgorithm::selection() {
    std::vector<int> selectedParents;

    double totalFitness = 0.0;
    for (double fitness : fitnessValues) {
        totalFitness += fitness;
    }

    // Normalizar los valores de ajuste (fitness) para obtener probabilidades de selección
    std::vector<double> probabilities;
    for (double fitness : fitnessValues) {
        probabilities.push_back(fitness / totalFitness);
    }

    // Fathers selection without replacement
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    while (selectedParents.size() < populationSize) {
        // Random number [0, 1)
        double r = dis(gen);
        double cumulativeProbability = 0.0;
        for (int i = 0; i < populationSize; ++i) {
            cumulativeProbability += probabilities[i];
            if (cumulativeProbability > r) {
                selectedParents.push_back(i);
                break;
            }
        }
    }

    return selectedParents;
}


// Single point crossover
void GeneticAlgorithm::crossover(std::vector<int>& selectedParents) {
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < selectedParents.size(); i += 2) {
        if (dis(gen) < crossoverRate) {
            // Randomly select crossover point
            int crossoverPoint = std::uniform_int_distribution<>(1, population[0].size() - 1)(gen);

            for (int j = crossoverPoint; j < population[0].size(); ++j) {
                std::swap(population[selectedParents[i]][j], population[selectedParents[i + 1]][j]);
            }
        }
    }
}

void GeneticAlgorithm::mutation() {
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    for (auto& individual : population) {
        for (double& gene : individual) {
            // Mutation at a given rate
            if (dis(gen) < mutationRate) {
                // Random change within the gene
                // Mean 0
                // Standard deviation 0.1
                double mutationChange = std::normal_distribution<double>(0.0, 0.1)(gen);
            }
        }
    }
}

bool GeneticAlgorithm::shouldStop(int currentGeneration) {
    double bestFitness = *std::min_element(fitnessValues.begin(), fitnessValues.end());
    double worstFitness = *std::max_element(fitnessValues.begin(), fitnessValues.end());
    double optimalGlobalFitness = 0.0; // Hay que ajustar esto;

    // |f(⃗xbest)−f(⃗x∗)| ≤ ε
    if (std::abs(bestFitness - optimalGlobalFitness) <= 0.001) {
        return true;
    }

    // f(⃗xworst)−f(⃗xbest)| ≤ ε
    if (std::abs(worstFitness - bestFitness) <= 0.001) {
        return true;
    }

    return false;
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
