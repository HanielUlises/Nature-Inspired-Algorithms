#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include <string>
#include <vector>

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
    void initialize_population();
    void remove_abrupts(std::vector<int>& tour, int m = 3);
    void evaluate_population();
    double calculate_fitness(const std::vector<int>& tour, int start_city);
    std::vector<int> cycle_crossover(const std::vector<int>& parent1, const std::vector<int>& parent2, int start_city);
    void next_generation();
    std::vector<int> _get_sorted_nearest_neighbors(int city, const std::vector<std::vector<double>>& dist_matrix, int m);
    std::vector<std::vector<int>> _generate_candidate_insertions(int city, const std::vector<int>& base_tour, const std::vector<int>& neighbors);
    std::vector<int> tournament_selection(int tournament_size);

    const std::vector<std::vector<double>>& distances_;
    const std::vector<City>& cities_;
    int population_size_;
    int generations_;
    double crossover_rate_;
    double mutation_rate_;
    double random_insert_prob_;
    double penalty_factor_;
    std::vector<std::vector<int>> population_;
    std::vector<int> start_cities_;
    std::vector<double> fitness_;
    double best_distance_;
    std::vector<int> best_tour_;
    std::vector<std::vector<int>> nearest_neighbors_cache_;
};

#endif