#include "Evolutionary_Strategy.h"
#include <iostream>
#include <cmath>
#include <random>

EvolutionaryStrategy::EvolutionaryStrategy(int population_size, int dimension, int generation_max, double lower_bound, double upper_bound, std::function<double(const std::vector<double>&)> objFunc, bool plusStrategy)
    : population_size_(population_size), dimension_(dimension), generation_max_(generation_max), lower_bound_(lower_bound), upper_bound_(upper_bound), plusStrategy_(plusStrategy), objectiveFunction_(objFunc), generator_(std::random_device{}()), distribution_(lower_bound, upper_bound) {
    initializePopulation();
}

void EvolutionaryStrategy::initializePopulation() {
    population_.resize(population_size_);
    for (auto& individual : population_) {
        individual.resize(dimension_);
        for (auto& x : individual) {
            x = distribution_(generator_);
        }
    }
}

double EvolutionaryStrategy::runEvolution() {
    for (int generation = 0; generation < generation_max_; ++generation) {
        auto parents = selectParents();
        generateOffspring(parents);
        evaluatePopulation();
        if (plusStrategy_) {
            competeForSurvival();
        }
    }
    // Return the best fitness value from the last generation
    return objectiveFunction_(population_.front());
}

void EvolutionaryStrategy::evaluatePopulation() {
    std::sort(population_.begin(), population_.end(), [this](const std::vector<double>& a, const std::vector<double>& b) {
        return objectiveFunction_(a) < objectiveFunction_(b);
    });
}

std::vector<std::vector<double>> EvolutionaryStrategy::selectParents() {
    std::vector<std::vector<double>> parents(population_size_);
    std::copy_n(population_.begin(), population_size_, parents.begin());
    return parents;
}

void EvolutionaryStrategy::generateOffspring(const std::vector<std::vector<double>>& parents) {
    population_.clear();
    for (int i = 0; i < population_size_; ++i) {
        std::vector<double> offspring = parents[i % population_size_]; // Create offspring from a parent
        mutateOffspring(offspring);
        population_.push_back(offspring);
    }
}

void EvolutionaryStrategy::mutateOffspring(std::vector<double>& offspring) {
    for (auto& gene : offspring) {
        gene += distribution_(generator_);
    }
}

void EvolutionaryStrategy::competeForSurvival() {
    if (population_.size() > population_size_) {
        std::sort(population_.begin(), population_.end(), [this](const std::vector<double>& a, const std::vector<double>& b) {
            return objectiveFunction_(a) < objectiveFunction_(b);
        });
        population_.resize(population_size_);
    }
}
