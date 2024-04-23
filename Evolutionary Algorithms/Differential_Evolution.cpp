#include "differential_evolution.h"

#include <iostream>

DifferentialEvolution::DifferentialEvolution(int population_size, int dimension)
    : population_size_(population_size), dimension_(dimension) {}

void DifferentialEvolution::runEvolution(const std::string& strategy) {
    std::cout << "Running Differential Evolution with strategy: " << strategy << std::endl;
}
