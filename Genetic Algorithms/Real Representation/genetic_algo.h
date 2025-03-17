#ifndef GENETIC_ALGO_H
#define GENETIC_ALGO_H

#include <vector>
#include <random>
#include <fstream>

#include "PM/polynomial_mut.h"
#include "SBX/sbx.h"
#include "../../Utils/ObjectiveFunctions.h"

class GeneticAlgorithmReal {
private:
    int population_size;
    int num_generations;
    int num_genes;
    double crossover_prob;
    double mutation_prob;
    std::vector<double> lower_bound;
    std::vector<double> upper_bound;

    std::ofstream output_file;

    std::vector<std::vector<double>> population;
    std::vector<double> fitness_values;

    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;

public:
    GeneticAlgorithmReal(int pop_size, int num_genes, int num_generations, double crossover_prob, 
                        double mutation_prob, std::vector<double> lower_bound, std::vector<double> upper_bound);
    
    ~GeneticAlgorithmReal();

    void initialize_population();
    void evaluate_fitness();
    int tournament_selection();
    void crossover();
    void mutate();
    void select_best_individual();
    void apply_elitism();
    void operator()();
};

#endif // GENETIC_ALGO_H