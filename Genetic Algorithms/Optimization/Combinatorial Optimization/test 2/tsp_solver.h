#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include <vector>
#include <string>

struct City {
    std::string name;
    double time_window_start;
    double time_window_end;
};

class TspSolver {
public:
    TspSolver(const std::vector<std::vector<double>>& distances,
              const std::vector<City>& cities,
              int population_size,
              int generations,
              double crossover_rate,
              double mutation_rate,
              double random_insert_prob,
              double penalty_factor);

    void solve();
    double get_best_distance() const;
    std::vector<int> get_best_tour() const;

private:
    std::vector<std::vector<double>> distances_;
    std::vector<City> cities_;
    int population_size_;
    int generations_;
    double crossover_rate_;
    double mutation_rate_;
    double random_insert_prob_;
    double penalty_factor_;

    std::vector<std::vector<int>> population_;
    std::vector<double> fitness_;
    std::vector<int> best_tour_;
    double best_distance_;

    void initialize_population();
    void evaluate_population();
    void remove_abrupts(std::vector<int>& tour);
    double calculate_fitness(const std::vector<int>& tour);
    double calculate_distance(const std::vector<int>& tour);
    std::vector<int> cycle_crossover(const std::vector<int>& parent1, const std::vector<int>& parent2);
    void next_generation();
};

#endif