#include <vector>
#include <cmath>
#include <iostream>
#include <numeric> 
#include <random> 

void f1(double x1, double x2, std::vector<double>& z){
    double goal = -(x1 + x2) + 400;
    z.push_back(goal);
}

void f2(double x1, double x2, std::vector<double>& z){
    double goal = -x1 * std::sin(std::sqrt(std::abs(x1))) - x2 * std::sin(std::sqrt(std::abs(x2))) + 837.9657745448674f;
    z.push_back(goal);
}

// Mean of goal set z
double mean (std::vector<double> z){
    if (z.empty()) return 0;
    double sum = 0;
    for(int i = 0; i < z.size(); i++){
        sum += z[i];
    }
    return sum / z.size();
}


void expected_val (double z, double mean_z, std::vector<double>& exp_z){
    double result = 0;
    result = z / mean_z;
    exp_z.push_back(result);
}

// Roulette selection without replacement
template<typename T>
T& rouletteSelection(std::vector<T>& population, const std::vector<double>& fitness) {
    static std::mt19937 rng(std::random_device{}());

    // Calculate the total fitness
    double totalFitness = std::accumulate(fitness.begin(), fitness.end(), 0.0);

    // Generate a random number between 0 and totalFitness
    std::uniform_real_distribution<double> dist(0.0, totalFitness);
    double randomSelection = dist(rng);

    // Iterate through the population to find the selected individual!!
    double accumulated = 0.0;
    for (size_t i = 0; i < population.size(); ++i) {
        accumulated += fitness[i];
        if (accumulated >= randomSelection) {
            return population[i];
        }
    }

    // In case the roulette does not select properly, return a default (first) element
    // This should not happen, but it's good practice tho
    return population.front();
}

// Roulette selection with replacement :D
template<typename T>
std::vector<size_t> rouletteSelectionWithReplacement(std::vector<T>& population, const std::vector<double>& fitness, size_t numSelections) {
    std::vector<size_t> selectedIndices;
    static std::mt19937 rng(std::random_device{}());
    double totalFitness = std::accumulate(fitness.begin(), fitness.end(), 0.0);

    std::uniform_real_distribution<double> dist(0.0, totalFitness);

    for (size_t n = 0; n < numSelections; ++n) {
        double randomSelection = dist(rng);
        double accumulated = 0.0;
        for (size_t i = 0; i < population.size(); ++i) {
            accumulated += fitness[i];
            if (accumulated >= randomSelection) {
                selectedIndices.push_back(i);
                break;
            }
        }
    }

    return selectedIndices;
}