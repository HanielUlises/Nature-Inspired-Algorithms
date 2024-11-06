#include "PSO.h"
#include <iostream>
#include "LinearRegression.h"
#include "../Utils/ObjectiveFunctions.h"

int main() {
    int dimensions = 10;  // Dimensions
    int swarm_size = 100;  // Particles
    int max_iterations = 60;
    // Start and end values for adaptive omega
    double omega_start = 0.9, omega_end = 0.4; 
    double phi_p = 1.0, phi_g = 2.0;
    double threshold = 1000.0; 

    // Lista de funciones objetivo y sus nombres correspondientes
    std::vector<std::function<double(const std::vector<double>&)>> objectiveFunctions = {
        rosenbrockFunction, ackleyFunction, GriewankFunction, RastriginFunction
    };
    std::vector<std::string> functionNames = {"Rosenbrock", "Ackley", "Griewank", "Rastrigin"};

    for (size_t i = 0; i < objectiveFunctions.size(); ++i) {
        std::cout << "Optimizing with " << functionNames[i] << " function...\n";
        std::vector<std::vector<double>> history_ejecutions; // Move this inside the loop
        for (size_t j = 0; j < 20; j++){
            PSO pso(swarm_size, dimensions, objectiveFunctions[i], max_iterations, threshold,i);
            pso.optimize(max_iterations, omega_start, omega_end, phi_p, phi_g);
            pso.printResults();
            std::cout << "\n"; 
            std::vector<double> ejecution = pso.history_global_best_score;
            history_ejecutions.push_back(ejecution);
        }
        plotting(history_ejecutions);
    }
    
    dimensions = 2;  // Dimensions
    swarm_size = 5000;  // Particles
    max_iterations = 100;
    omega_start = 0.9, omega_end = 0.4;
    PSO pso(swarm_size, dimensions, objectiveFunctions[0], max_iterations, threshold,3); // example with the first function
    pso.minimizeError(max_iterations, omega_start, phi_p, phi_g);

    return 0;
}
