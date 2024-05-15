#include "Evolutionary_Strategy.h"
#include <iostream>
#include <cmath>
#include <random>

EvolutionaryStrategy::EvolutionaryStrategy(int population_size, int dimension, int generation_max, double lower_bound, double upper_bound, std::function<double(const std::vector<double>&)> objFunc, bool plusStrategy)
    : population_size_(population_size), dimension_(dimension), generation_max_(generation_max),
      generation_(0), lower_bound_(lower_bound), upper_bound_(upper_bound), plusStrategy_(plusStrategy),
      objectiveFunction_(objFunc), generator_(std::random_device{}()), distribution_(0.0, 1.0),
      mutation_strength_(0.1), successfulMutations_(0), totalMutations_(0) {
    initializePopulation();
}

void EvolutionaryStrategy::initializePopulation() {
    population_.resize(population_size_);
    for (auto& individual : population_) {
        individual.resize(dimension_);
        for (auto& x : individual) {
            x = lower_bound_ + (upper_bound_ - lower_bound_) * ((double)rand() / RAND_MAX);
        }
    }
}

double EvolutionaryStrategy::runEvolution() {
    std::vector<std::vector<double>> newPopulation;
    for (generation_ = 0; generation_ < generation_max_; ++generation_) {
        auto parents = selectParents();
        generateOffspring(parents);
        if (plusStrategy_) {
            newPopulation.insert(newPopulation.end(), parents.begin(), parents.end());
        }
        newPopulation.insert(newPopulation.end(), population_.begin(), population_.end());
        elitism(newPopulation);
        population_ = std::move(newPopulation);
        newPopulation.clear();
        evaluatePopulation();
        adaptMutationStrength();
    }
    return objectiveFunction_(population_.front());
}


double EvolutionaryStrategy::calculateSuccessRate() {
    if (totalMutations_ == 0) return 0.0; 
    double successRate = static_cast<double>(successfulMutations_) / totalMutations_;
    successfulMutations_ = 0;  // Reset after calculation
    totalMutations_ = 0;       // Reset after calculation
    return successRate;
}

double EvolutionaryStrategy::sharingFunction(double distance, double sigmaShare) {
    if (distance < sigmaShare) {
        return 1 - (distance / sigmaShare);
    }
    return 0;
}

double EvolutionaryStrategy::euclideanDistance(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return std::sqrt(sum);
}


void EvolutionaryStrategy::evaluatePopulation() {
    std::vector<double> fitnesses(population_size_);
    for (int i = 0; i < population_size_; ++i) {
        fitnesses[i] = objectiveFunction_(population_[i]);
    }

    // Apply fitness sharing
    std::vector<double> sharedFitnesses(population_size_, 0.0);
    double sigmaShare = 0.5;  // This value might need tuning based on the problem domain

    for (int i = 0; i < population_size_; ++i) {
        double nicheCount = 0.0;
        for (int j = 0; j < population_size_; ++j) {
            double distance = euclideanDistance(population_[i], population_[j]);
            nicheCount += sharingFunction(distance, sigmaShare);
        }
        sharedFitnesses[i] = fitnesses[i] / nicheCount;  // Adjust fitness by the niche count
    }

    // Sort the population based on shared fitness
    std::vector<int> indices(population_size_);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&sharedFitnesses](int i1, int i2) {
        return sharedFitnesses[i1] < sharedFitnesses[i2];
    });
    std::vector<std::vector<double>> newPopulation(population_size_);
    for (int i = 0; i < population_size_; ++i) {
        newPopulation[i] = population_[indices[i]];
    }
    population_ = std::move(newPopulation);
}

std::vector<std::vector<double>> EvolutionaryStrategy::selectParents() {
    std::vector<std::vector<double>> parents;
    int tournament_size = 5; 

    for (int i = 0; i < population_size_; ++i) {
        int best_index = rand() % population_.size();
        for (int j = 1; j < tournament_size; ++j) {
            int idx = rand() % population_.size();
            if (objectiveFunction_(population_[idx]) < objectiveFunction_(population_[best_index])) {
                best_index = idx;
            }
        }
        parents.push_back(population_[best_index]);
    }
    return parents;
}


void EvolutionaryStrategy::generateOffspring(const std::vector<std::vector<double>>& parents) {
    std::vector<std::vector<double>> offspring;
    recombineParents(offspring, parents);
    for (int i = 0; i < offspring.size(); i++) {
        mutateOffspring(offspring[i], parents[i % parents.size()]);
    }
    population_ = std::move(offspring);
}


void EvolutionaryStrategy::recombineParents(std::vector<std::vector<double>>& offspring, const std::vector<std::vector<double>>& parents) {
    // BLX-Alpha
    double alpha = 0.5;  // Blend factor
    for (int i = 0; i < population_size_; i += 2) {
        int p1 = i % parents.size();
        int p2 = (i + 1) % parents.size();

        std::vector<double> child1(dimension_), child2(dimension_);
        for (int d = 0; d < dimension_; ++d) {
            double range = std::abs(parents[p1][d] - parents[p2][d]);
            double min_val = std::min(parents[p1][d], parents[p2][d]);
            child1[d] = min_val - alpha * range + (1 + 2 * alpha) * range * ((double)rand() / RAND_MAX);
            child2[d] = min_val - alpha * range + (1 + 2 * alpha) * range * ((double)rand() / RAND_MAX);
        }

        offspring.push_back(child1);
        offspring.push_back(child2);
    }
}

void EvolutionaryStrategy::mutateOffspring(std::vector<double>& offspring, const std::vector<double>& parent) {
    double originalFitness = objectiveFunction_(parent);
    for (auto& gene : offspring) {
        gene += mutation_strength_ * distribution_(generator_);
    }
    double newFitness = objectiveFunction_(offspring);
    totalMutations_++;
    if (newFitness < originalFitness) {
        successfulMutations_++;
    }
}

void EvolutionaryStrategy::elitism(std::vector<std::vector<double>>& newPopulation) {
    std::sort(newPopulation.begin(), newPopulation.end(), [this](const std::vector<double>& a, const std::vector<double>& b) {
        return objectiveFunction_(a) < objectiveFunction_(b);
    });
    // Select the top mu individuals to form the next generation
    newPopulation.resize(population_size_);
}

void EvolutionaryStrategy::adaptMutationStrength() {
    double successRate = calculateSuccessRate();
    double targetRate = 0.20;
    if (successRate > targetRate) {
        mutation_strength_ *= 1.2;
    } else {
        mutation_strength_ /= 1.2;
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