#include "utils.h"
#include <numeric>

namespace NatureInspired {

float generateRandomFloat(float lower_bound, float upper_bound) {
    static std::random_device rd;
    static std::mt19937 generator(rd());
    std::uniform_real_distribution<float> distribution(lower_bound, upper_bound);
    return distribution(generator);
}

std::vector<float> initializePopulation(size_t size, float lower_bound, float upper_bound) {
    std::vector<float> population(size);
    for (auto& individual : population) {
        individual = generateRandomFloat(lower_bound, upper_bound);
    }
    return population;
}

std::vector<float> evaluateFitness(const std::vector<float>& population, std::function<float(float)> fitness_function) {
    std::vector<float> fitness_scores(population.size());
    for (size_t i = 0; i < population.size(); ++i) {
        fitness_scores[i] = fitness_function(population[i]);
    }
    return fitness_scores;
}

std::vector<float> selectParents(const std::vector<float>& population, const std::vector<float>& fitness_scores) {
    std::vector<float> selected_parents;
    float total_fitness = std::accumulate(fitness_scores.begin(), fitness_scores.end(), 0.0f);
    if (total_fitness == 0.0f) return population; 

    std::vector<float> probabilities(fitness_scores.size());
    for (size_t i = 0; i < fitness_scores.size(); ++i) {
        probabilities[i] = fitness_scores[i] / total_fitness;
    }

    static std::random_device rd;
    static std::mt19937 generator(rd());
    std::discrete_distribution<size_t> distribution(probabilities.begin(), probabilities.end());

    for (size_t i = 0; i < population.size(); ++i) {
        selected_parents.push_back(population[distribution(generator)]);
    }

    return selected_parents;
}

float crossover(float parent1, float parent2) {
    return (parent1 + parent2) / 2.0f;
}

float mutate(float individual, float mutation_rate, float lower_bound, float upper_bound) {
    if (generateRandomFloat(0.0f, 1.0f) < mutation_rate) {
        return generateRandomFloat(lower_bound, upper_bound);
    }
    return individual;
}

} 
