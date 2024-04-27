#include "Evolutionary_Strategy.h"
#include <iostream>
#include <cmath>

EvolutionaryStrategy::EvolutionaryStrategy(int mu, int lambda, std::function<double(const std::vector<double>&)> objFunc, bool plusStrategy)
    : mu_(mu), lambda_(lambda), plusStrategy_(plusStrategy), objectiveFunction_(objFunc), generator_(std::random_device{}()), distribution_(0.0, 1.0) {
    initializePopulation();
}

void EvolutionaryStrategy::initializePopulation() {
    population_.resize(mu_);
    for (auto& individual : population_) {
        individual.resize(10);  // Assuming each individual has 10 dimensions
        for (auto& x : individual) {
            x = distribution_(generator_);  // Initialize with random values from a normal distribution
        }
    }
}

double EvolutionaryStrategy::runEvolution() {
    for (int generation = 0; generation < 100; ++generation) {
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
    std::vector<std::vector<double>> parents(mu_);
    std::copy_n(population_.begin(), mu_, parents.begin());
    return parents;
}

void EvolutionaryStrategy::generateOffspring(const std::vector<std::vector<double>>& parents) {
    population_.clear();
    for (int i = 0; i < lambda_; ++i) {
        std::vector<double> offspring = parents[i % mu_]; // Create offspring from a parent
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
    if (population_.size() > mu_) {
        std::sort(population_.begin(), population_.end(), [this](const std::vector<double>& a, const std::vector<double>& b) {
            return objectiveFunction_(a) < objectiveFunction_(b);
        });
        population_.resize(mu_);
    }
}
