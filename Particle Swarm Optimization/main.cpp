#include "PSO.h"
#include <iostream>
#include "../Utils/ObjectiveFunctions.h"

int main() {
    int dimensions = 10;  // Número de dimensiones
    int swarm_size = 30;  // Número de partículas
    int max_iterations = 1000;
    double omega = 0.5, phi_p = 1.0, phi_g = 2.0;

    // Lista de funciones objetivo y sus nombres correspondientes
    std::vector<std::function<double(const std::vector<double>&)>> objectiveFunctions = {
        rosenbrockFunction, ackleyFunction, GriewankFunction, RastriginFunction
    };
    std::vector<std::string> functionNames = {"Rosenbrock", "Ackley", "Griewank", "Rastrigin"};

    for (size_t i = 0; i < objectiveFunctions.size(); ++i) {
        std::cout << "Optimizing with " << functionNames[i] << " function...\n";

        PSO pso(swarm_size, dimensions, objectiveFunctions[i]);
        pso.optimize(max_iterations, omega, phi_p, phi_g);
        pso.printResults();

        std::cout << "\n"; 
    }

    return 0;
}