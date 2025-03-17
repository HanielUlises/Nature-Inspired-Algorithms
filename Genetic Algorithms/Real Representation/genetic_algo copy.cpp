#include "genetic_algo.h"
#include <iostream>
#include <algorithm>

GeneticAlgorithmReal::GeneticAlgorithmReal(int pop_size, int num_genes, int num_generations, 
                                         double crossover_prob, double mutation_prob, 
                                         std::vector<double> lower_bound, std::vector<double> upper_bound)
    : population_size(pop_size), num_genes(num_genes), num_generations(num_generations),
      crossover_prob(crossover_prob), mutation_prob(mutation_prob),
      lower_bound(lower_bound), upper_bound(upper_bound),
      gen(rd()), dis(0.0, 1.0) {

    population.resize(population_size, std::vector<double>(num_genes));
    fitness_values.resize(population_size);
    initialize_population();
    
    output_file.open("results.txt");
    if (!output_file) {
        std::cerr << "Error opening file.\n";
        exit(1);
    }
}

GeneticAlgorithmReal::~GeneticAlgorithmReal() {
    if (output_file.is_open()) {
        output_file.close();
    }
}

void GeneticAlgorithmReal::initialize_population() {
    for (auto &individual : population) {
        for (int i = 0; i < num_genes; ++i) {
            individual[i] = lower_bound[i] + dis(gen) * (upper_bound[i] - lower_bound[i]);
        }
    }
}

void GeneticAlgorithmReal::evaluate_fitness() {
    for (int i = 0; i < population_size; ++i) {
        // fitness_values[i] = Langermann(population[i]);
        // fitness_values[i] = dropWave(population[i]);
    }
}

int GeneticAlgorithmReal::tournament_selection() {
    std::uniform_int_distribution<> dis_int(0, population_size - 1);
    int idx1 = dis_int(gen);
    int idx2 = dis_int(gen);
    return (fitness_values[idx1] < fitness_values[idx2]) ? idx1 : idx2;
}

void GeneticAlgorithmReal::crossover() {
    std::vector<std::vector<double>> new_population = population;

    for (int i = 0; i < population_size / 2; ++i) {
        int parent1_idx = tournament_selection();
        int parent2_idx = tournament_selection();

        if (dis(gen) <= crossover_prob) {
            auto [offspring1, offspring2] = sbx(population[parent1_idx], population[parent2_idx], 
                                               crossover_prob, lower_bound, upper_bound);
            new_population[parent1_idx] = offspring1;
            new_population[parent2_idx] = offspring2;
        }
    }

    population = new_population;
}

void GeneticAlgorithmReal::mutate() {
    for (int i = 0; i < population_size; ++i) {
        if (dis(gen) <= mutation_prob) {
            population[i] = polynomial_mutation(population[i], mutation_prob, 20);
        }
    }
}

void GeneticAlgorithmReal::select_best_individual() {
    int best_idx = std::min_element(fitness_values.begin(), fitness_values.end()) - fitness_values.begin();
    
    std::cout << "Best individual: ";
    for (double v : population[best_idx]) {
        std::cout << v << " ";
        output_file << v << " ";
    }
    
    std::cout << " | Fitness: " << fitness_values[best_idx] << "\n";
    output_file << " | Fitness: " << fitness_values[best_idx] << "\n";
}

void GeneticAlgorithmReal::apply_elitism() {
    int best_idx = std::min_element(fitness_values.begin(), fitness_values.end()) - fitness_values.begin();
    int worst_idx = std::max_element(fitness_values.begin(), fitness_values.end()) - fitness_values.begin();
    
    population[worst_idx] = population[best_idx];
    fitness_values[worst_idx] = fitness_values[best_idx];
}

void GeneticAlgorithmReal::operator()() {
    for (int generation = 0; generation < num_generations; ++generation) {
        std::cout << "Generation " << generation + 1 << std::endl;
        output_file << "Generation: " << generation + 1 << " | ";

        evaluate_fitness();
        select_best_individual();

        crossover();
        mutate();
        evaluate_fitness();

        apply_elitism();
    }
}

// int main() {
//     int population_size = 50;
//     int num_genes = 2;
//     int num_generations = 200;
//     double crossover_prob = 0.9;
//     double mutation_prob = 0.03;


//     // // For Langermann
//     // std::vector<double> lower_bound = {0.0, 0.0}; 
//     // std::vector<double> upper_bound = {10.0, 10.0}; 

//     // For drop wave
//     std::vector<double> lower_bound = {-5.12,-5,12};
//     std::vector<double> upper_bound = {5.12, 5.12}; 

//     GeneticAlgorithmReal ga(population_size, num_genes, num_generations, crossover_prob, mutation_prob, lower_bound, upper_bound);
//     ga();

//     return 0;
// }
