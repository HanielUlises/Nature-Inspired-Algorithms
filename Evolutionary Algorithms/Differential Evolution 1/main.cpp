#include "Evolutionary_Strategy.h"
#include "Differential_Evolution.h"
#include "../../Utils/ObjectiveFunctions.h"
#include <numeric>
#include <cmath>
#include <iostream>

int main() {
    // Constants for selection and recombination strategies
    const int rand = 0;
    const int best = 1; // "best" strategy for mutation
    const int bin = 0;  // "binomial" crossover
    const int expo = 1; // "exponential" crossover

    // Parameters 
    int population_size = 50;
    int dimension = 10;
    int generation_max = 1000;
    double lower_bound = -5.12;
    double upper_bound = 5.12;
    double F = 0.6;
    double Cr = 0.9;

    // Set of objective functions
    std::vector<std::function<double(const std::vector<double>&)>> objectiveFunctions = {
        rosenbrockFunction, ackleyFunction, GriewankFunction, RastriginFunction
    };
    std::vector<std::string> functionNames = {"Rosenbrock", "Ackley", "Griewank", "Rastrigin"};

    // Differential Evolution Tests
    /*
    std::cout << "Differential Evolution Tests" << std::endl;
    for (size_t funcIndex = 0; funcIndex < objectiveFunctions.size(); ++funcIndex) {
        std::cout << "============================" << std::endl;
        std::vector<double> deBestHistory;
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
    */
    // Evolutionary strategies
    std::cout << "\nEvolutionary Strategy Tests" << std::endl;
    for (bool plusStrategy : {true, false}) {
        std::cout << "============================" << std::endl;
        std::cout << "Testing with strategy: " << (plusStrategy ? "μ+λ" : "μ,λ") << std::endl;

        for (size_t funcIndex = 0; funcIndex < objectiveFunctions.size(); ++funcIndex) {
            
            std::cout << "============================" << std::endl;
            std::cout << "Using " << functionNames[funcIndex] << " function." << std::endl;
            std::vector<double> esBestHistory;

            for (size_t i = 0; i < 20; i++) {
                EvolutionaryStrategy es(population_size, dimension, generation_max, lower_bound, upper_bound, objectiveFunctions[funcIndex], plusStrategy);
                double bestResult = es.runEvolution();
                esBestHistory.push_back(bestResult);
            }

            double average = std::accumulate(esBestHistory.begin(), esBestHistory.end(), 0.0) / esBestHistory.size();
            double variation = std::accumulate(esBestHistory.begin(), esBestHistory.end(), 0.0,
                [average](double acc, double x) { return acc + (x - average) * (x - average); });
            double stddev = std::sqrt(variation / esBestHistory.size());

            std::cout << "Average Best Fitness: " << average << std::endl;
            std::cout << "Standard Deviation: " << stddev << std::endl;
        }
    }

    return 0;
}
