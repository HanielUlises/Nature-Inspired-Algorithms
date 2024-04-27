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
        // Individuals with dimension of 10
        individual.resize(10); 
        for (auto& x : individual) {
            // Initializes with random values from a normal distribution
            x = distribution_(generator_); 
        }
    }
}

void EvolutionaryStrategy::runEvolution() {
    std::vector<std::vector<double>> allCandidates;
    for (int generation = 0; generation < 100; ++generation) {
        auto parents = selectParents();
        generateOffspring(parents);
        if (plusStrategy_) {
            allCandidates = population_;
            allCandidates.insert(allCandidates.end(), parents.begin(), parents.end());
            population_ = allCandidates;
        }
        evaluatePopulation();
        if (plusStrategy_) {
            competeForSurvival();
        }
        std::cout << "Generation " << generation << ": Best Fitness = " << objectiveFunction_(population_.front()) << std::endl;
    }
}

void EvolutionaryStrategy::competeForSurvival() {
    // Resize the population to keep only the best mu individuals
    if (population_.size() > mu_) {
        std::sort(population_.begin(), population_.end(), [this](const std::vector<double>& a, const std::vector<double>& b) {
            return objectiveFunction_(a) < objectiveFunction_(b);
        });
        population_.resize(mu_);
    }
}

void EvolutionaryStrategy::evaluatePopulation() {
    // Sort the population based on the fitness values determined by the objective function
    std::sort(population_.begin(), population_.end(), [this](const std::vector<double>& a, const std::vector<double>& b) {
        return objectiveFunction_(a) < objectiveFunction_(b);
    });
}

std::vector<std::vector<double>> EvolutionaryStrategy::selectParents() {
    // Select the best mu individuals as parents
    std::vector<std::vector<double>> parents(mu_);
    std::copy_n(population_.begin(), mu_, parents.begin());
    return parents;
}

void EvolutionaryStrategy::generateOffspring(const std::vector<std::vector<double>>& parents) {
    population_.clear();
    for (int i = 0; i < lambda_; ++i) {
        int parent_index = i % mu_;
        std::vector<double> offspring = parents[parent_index]; // Create offspring from a parent
        mutateOffspring(offspring);                           // Mutate the offspring
        population_.push_back(offspring);                     // Add mutated offspring to the population
    }
}

void EvolutionaryStrategy::mutateOffspring(std::vector<double>& offspring) {
    for (auto& gene : offspring) {
        gene += distribution_(generator_);  // Apply a Gaussian mutation to each gene
    }
}
