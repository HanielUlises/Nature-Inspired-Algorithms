#include "Evolutionary_Strategy.h"
#include "Differential_Evolution.h"
#include "ObjectiveFunctions.h"
#include <numeric>
#include <cmath>
#include <iostream>

int main() {
    // Parameters for Differential Evolution
    int population_size = 50;
    int dimension = 10;
    int generation_max = 1000;
    double lower_bound = -5.12;
    double upper_bound = 5.12;
    double F = 0.6;
    double Cr = 0.9;

    // Parameters for Evolutionary Strategy
    int mu = 20;
    int lambda = 30;

    // Set of objective functions
    std::vector<std::function<double(const std::vector<double>&)>> objectiveFunctions = {
        rosenbrockFunction, ackleyFunction, GriewankFunction, RastriginFunction
    };
    std::vector<std::string> functionNames = {"Rosenbrock", "Ackley", "Griewank", "Rastrigin"};

    // Differential Evolution Tests
    // ----------------------------
    std::cout << "Differential Evolution Tests" << std::endl;
    std::vector<double> deBestHistory;
    for (size_t funcIndex = 0; funcIndex < objectiveFunctions.size(); ++funcIndex) {
        deBestHistory.clear();
        std::cout << "Testing " << functionNames[funcIndex] << " Function" << std::endl;

        for (size_t i = 0; i < 20; i++) {
            DifferentialEvolution de(population_size, dimension, generation_max, lower_bound, upper_bound, F, Cr);
            double bestResult = de.runEvolution(best, expo);
            deBestHistory.push_back(bestResult);
        }

        double deAverage = std::accumulate(deBestHistory.begin(), deBestHistory.end(), 0.0) / deBestHistory.size();
        double deVariation = 0.0;
        for (double best : deBestHistory) {
            deVariation += std::pow((best - deAverage), 2);
        }
        double deStdDev = std::sqrt(deVariation / deBestHistory.size());
        std::cout << "Average Best: " << deAverage << ", Standard Deviation: " << deStdDev << std::endl;
    }

    // Evolutionary Strategy Tests
    // ---------------------------
    std::cout << "\nEvolutionary Strategy Tests" << std::endl;
    for (bool plusStrategy : {false, true}) {
        std::cout << "Testing with strategy: " << (plusStrategy ? "μ+λ" : "μ,λ") << std::endl;

        for (size_t funcIndex = 0; funcIndex < objectiveFunctions.size(); ++funcIndex) {
            std::vector<double> esBestHistory;
            std::cout << "Using " << functionNames[funcIndex] << " function." << std::endl;

            for (size_t i = 0; i < 20; i++) {
                EvolutionaryStrategy es(mu, lambda, objectiveFunctions[funcIndex], plusStrategy);
                es.runEvolution();
                double bestR = es.getBestFitness();  // Assume getBestFitness() returns the best fitness from the last generation
                esBestHistory.push_back(bestR);
            }

            double average = std::accumulate(esBestHistory.begin(), esBestHistory.end(), 0.0) / esBestHistory.size();
            double variation = 0.0;
            for (double best : esBestHistory) {
                variation += (std::pow((best - average), 2) / (esBestHistory.size() - 1));
            }

            double stddev = std::sqrt(variation);
            std::cout << "Average Best Fitness: " << average << std::endl;
            std::cout << "Standard Deviation: " << stddev << std::endl;
        }
    }

    return 0;
}
