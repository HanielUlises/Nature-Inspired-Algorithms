#include <vector>
#include <pair>

#include "PM/polynomial_mut.h"
#include "SBX/sbx.h"
#include "../Utils/ObjectiveFunctions.h"

class GeneticAlgorithmReal {
private:
    int population_size;
    int num_generations;
    int num_genes;
    double crossover_prob;
    double mutation_prob;
    double lower_bound;
    double upper_bound;

    std::vector<std::vector<double>> population;
    std::vector<double> fitness_values;

public:
    GeneticAlgorithmReal(int pop_size, int num_genes, int num_generations, double crossover_prob, double mutation_prob, double lower_bound, double upper_bound) 
        : population_size(pop_size), num_genes(num_genes), num_generations(num_generations), crossover_prob(crossover_prob), mutation_prob(mutation_prob), lower_bound(lower_bound), upper_bound(upper_bound) {
        // Initialize population
        population.resize(population_size, std::vector<double>(num_genes));
        fitness_values.resize(population_size);
        initialize_population();
    }

    // Initialize population with random values within bounds
    void initialize_population() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(lower_bound, upper_bound);

        for (int i = 0; i < population_size; ++i) {
            for (int j = 0; j < num_genes; ++j) {
                population[i][j] = dis(gen);
            }
        }
    }

    void evaluate_fitness() {
        for (int i = 0; i < population_size; ++i) {
            fitness_values[i] = 0;
            for (int j = 0; j < num_genes; ++j) {
                fitness_values[i] += 
            }
        }
    }

    // Tournament selection (2 individuals randomly selected, the one with the best fitness is returned)
    int tournament_selection() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, population_size - 1);

        int idx1 = dis(gen);
        int idx2 = dis(gen);

        return (fitness_values[idx1] < fitness_values[idx2]) ? idx1 : idx2;
    }

    // Crossover (Simulated Binary Crossover - SBX)
    void crossover() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        for (int i = 0; i < population_size / 2; ++i) {
            int parent1_idx = tournament_selection();
            int parent2_idx = tournament_selection();

            if (dis(gen) <= crossover_prob) {
                auto [offspring1, offspring2] = sbx(population[parent1_idx], population[parent2_idx]);

                population[parent1_idx] = offspring1;
                population[parent2_idx] = offspring2;
            }
        }
    }

    // Mutation (Polynomial Mutation)
    void mutate() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        for (int i = 0; i < population_size; ++i) {
            if (dis(gen) <= mutation_prob) {
                population[i] = polynomial_mutation(population[i], dis(gen), 20);
            }
        }
    }

    // Select the best individual based on fitness
    void select_best_individual() {
        int best_idx = std::min_element(fitness_values.begin(), fitness_values.end()) - fitness_values.begin();
        std::cout << "Best individual: ";
        for (double v : population[best_idx]) {
            std::cout << v << " ";
        }
        std::cout << "\n";
    }

    void operator() (){
        for (int generation = 0; generation < num_generations; ++generation) {
            std::cout << "Generation " << generation + 1 << std::endl;

            evaluate_fitness();
            select_best_individual();
            crossover();
            mutate();
        }
    }
};

int main() {
    int population_size = 10;
    int num_genes = 2;
    int num_generations = 50;
    double crossover_prob = 0.9;
    double mutation_prob = 0.1;
    double lower_bound = -5.0;
    double upper_bound = 5.0;

    GeneticAlgorithmReal ga(population_size, num_genes, num_generations, crossover_prob, mutation_prob, lower_bound, upper_bound);
    ga();

    return 0;
}
