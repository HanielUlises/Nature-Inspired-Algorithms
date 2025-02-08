#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <random>
#include <functional>

namespace NatureInspired {

// Generates a random real number within a range
float generateRandomFloat(float lower_bound, float upper_bound);

// Initializes a vector of random floats of given size and range
std::vector<float> initializePopulation(size_t size, float lower_bound, float upper_bound);

// Applies a fitness function to a population and returns the evaluated scores
std::vector<float> evaluateFitness(const std::vector<float>& population, std::function<float(float)> fitness_function);

// Selects individuals from a population based on their fitness scores (roulette-wheel selection)
std::vector<float> selectParents(const std::vector<float>& population, const std::vector<float>& fitness_scores);

// Performs crossover between two individuals
float crossover(float parent1, float parent2);

// Mutates an individual by a given mutation rate
float mutate(float individual, float mutation_rate, float lower_bound, float upper_bound);

}

#endif