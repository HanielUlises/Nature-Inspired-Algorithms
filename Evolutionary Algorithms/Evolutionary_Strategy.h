#ifndef EVOLUTIONARY_STRATEGY_H
#define EVOLUTIONARY_STRATEGY_H

#include <vector>
#include <random>
#include <algorithm>
#include <functional>

class EvolutionaryStrategy {
public:
    EvolutionaryStrategy(int mu, int lambda, std::function<double(const std::vector<double>&)> objFunc, bool plusStrategy = false);
    void runEvolution();

private:
    int mu_;
    int lambda_;
    bool plusStrategy_;
    std::vector<std::vector<double>> population_;
    std::default_random_engine generator_;
    std::normal_distribution<double> distribution_;
    std::function<double(const std::vector<double>&)> objectiveFunction_;  // Objective function

    void initializePopulation();
    void evaluatePopulation();
    std::vector<std::vector<double>> selectParents();
    void generateOffspring(const std::vector<std::vector<double>>& parents);
    void mutateOffspring(std::vector<double>& offspring);
    void competeForSurvival(); // Function to handle μ+λ competition
};

#endif // EVOLUTIONARY_STRATEGY_H
