#ifndef EVOLUTIONARY_STRATEGY_H
#define EVOLUTIONARY_STRATEGY_H

#include <vector>
#include <functional>
#include <random>
#include <numeric>

class EvolutionaryStrategy {
public:
    EvolutionaryStrategy(int population_size, int dimension, int generation_max, double lower_bound, double upper_bound, std::function<double(const std::vector<double>&)> objFunc, bool plusStrategy);

    double runEvolution();

private:
    int population_size_;
    int dimension_;
    int generation_max_;
    int generation_;
    double lower_bound_;
    double upper_bound_;
    bool plusStrategy_;

    std::vector<std::vector<double>> population_;
    std::function<double(const std::vector<double>&)> objectiveFunction_;

    std::default_random_engine generator_;
    std::uniform_real_distribution<double> distribution_;

    double mutation_strength_;
    int successfulMutations_;
    int totalMutations_;

    void initializePopulation();
    void evaluatePopulation();
    void adaptMutationStrength();
    std::vector<std::vector<double>> selectParents();
    void generateOffspring(const std::vector<std::vector<double>>& parents);
    void recombineParents(std::vector<std::vector<double>>& offspring, const std::vector<std::vector<double>>& parents);
    void mutateOffspring(std::vector<double>& offspring, const std::vector<double>& parent);
    void elitism(std::vector<std::vector<double>>& newPopulation);
    double calculateSuccessRate();
    void competeForSurvival();

    double euclideanDistance(const std::vector<double>& a, const std::vector<double>& b);
    double sharingFunction(double distance, double sigmaShare);
};

#endif // EVOLUTIONARY_STRATEGY_H
