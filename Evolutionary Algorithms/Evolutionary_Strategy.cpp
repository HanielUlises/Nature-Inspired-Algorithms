#include "evolutionary_strategy.h"

#include <iostream>

EvolutionaryStrategy::EvolutionaryStrategy(int mu, int lambda)
    : mu_(mu), lambda_(lambda) {}

void EvolutionaryStrategy::runEvolution() {
    std::cout << "Running evolutionary strategy with mu=" << mu_ << " and lambda=" << lambda_ << std::endl;
}
