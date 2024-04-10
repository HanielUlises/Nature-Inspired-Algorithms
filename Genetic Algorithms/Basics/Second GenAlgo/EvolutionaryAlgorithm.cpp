// EvolutionaryAlgorithm.cpp
#include "EvolutionaryAlgorithm.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>
#include <ctime>

EvolutionaryAlgorithm::EvolutionaryAlgorithm(int size, int populationSize, int generations)
    : size(size), populationSize(populationSize), generations(generations) {
    std::srand(static_cast<unsigned int>(std::time(nullptr))); // Seed for random number generation
}

void EvolutionaryAlgorithm::initializePopulation() {
    population.clear();
    for (int i = 0; i < populationSize; ++i) {
        population.push_back(generateIndividual());
    }
}

std::vector<int> EvolutionaryAlgorithm::generateIndividual() {
    std::vector<int> individual(size * size);
    std::iota(individual.begin(), individual.end(), 1);
    std::random_shuffle(individual.begin(), individual.end()); // Shuffle to generate a random permutation
    return individual;
}

int EvolutionaryAlgorithm::calculateFitness(const std::vector<int>& individual) {
    int magicConstant = size * (size * size + 1) / 2;
    int fitness = 0;

    // Check rows and columns
    for (int i = 0; i < size; ++i) {
        int rowSum = 0, colSum = 0;
        for (int j = 0; j < size; ++j) {
            rowSum += individual[i * size + j];
            colSum += individual[j * size + i];
        }
        fitness += abs(magicConstant - rowSum) + abs(magicConstant - colSum);
    }

    // Check diagonals
    int diagSum1 = 0, diagSum2 = 0;
    for (int i = 0; i < size; ++i) {
        diagSum1 += individual[i * size + i];
        diagSum2 += individual[(size - i - 1) * size + i];
    }
    fitness += abs(magicConstant - diagSum1) + abs(magicConstant - diagSum2);

    return fitness;
}

void EvolutionaryAlgorithm::selection() {
    // Placeholder for selection logic
    // Typically, you might use tournament selection, roulette wheel selection, or rank selection
}

void EvolutionaryAlgorithm::crossover(std::vector<int>& parent1, std::vector<int>& parent2) {
    // Placeholder for crossover logic
    // One common approach is a single-point crossover
}

void EvolutionaryAlgorithm::mutation(std::vector<int>& individual) {
    // Placeholder for mutation logic
    // A simple mutation might involve swapping two elements
}

bool EvolutionaryAlgorithm::isMagicSquare(const std::vector<int>& square) {
    // This method utilizes calculateFitness and checks if the fitness is 0, indicating a perfect magic square
    return calculateFitness(square) == 0;
}

void EvolutionaryAlgorithm::solve() {
    initializePopulation();
    for (int gen = 0; gen < generations; ++gen) {
        // Selection
        selection();

        // Crossover and Mutation
        for (int i = 0; i < populationSize; i += 2) {
            crossover(population[i], population[i+1]);
            mutation(population[i]);
            mutation(population[i+1]);
        }

        // Optionally, implement elitism or other mechanisms to preserve the best individual(s)

        // Find and print solution if one exists in the current generation
        for (const auto& individual : population) {
            if (isMagicSquare(individual)) {
                std::cout << "Magic square found in generation " << gen << ":\n";
                printSolution(individual);
                return;
            }
        }
    }
    std::cout << "Solution not found in " << generations << " generations.\n";
}

void EvolutionaryAlgorithm::printSolution(const std::vector<int>& solution) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            std::cout << solution[i * size + j] << ' ';
        }
        std::cout << std::endl;
    }
}
