#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include <vector>
#include <string>
#include <random>

struct City {
    std::string name;
    double window_start;
    double window_end;
};

class TspSolver {
public:
    TspSolver(const std::vector<City>& cities, const std::vector<std::vector<double>>& distances);
    
    void run(int population_size, int generations, double mutation_rate, double random_injection_prob);
    
    std::vector<int> get_best_route() const;
    double get_best_route_time() const;

private:
    std::vector<City> cities_;
    std::vector<std::vector<double>> distances_;
    std::vector<std::vector<int>> population_;
    std::vector<int> best_route_;
    double best_route_time_;
    mutable std::mt19937 gen_;
    
    double evaluate_route(const std::vector<int>& route) const;
    std::vector<int> generate_random_route() const;
    void mutate(std::vector<int>& route);
    
    std::vector<int> cycle_crossover(const std::vector<int>& parent1, const std::vector<int>& parent2);
    
    std::vector<int> remove_abrupts(const std::vector<int>& individual);
    
    void initialize_population(int population_size);
    void apply_remove_abrupts_to_population();
    
    void next_generation(double mutation_rate, double random_injection_prob);
    
    int select_parent_index(int max_index);
    
    void select_survivors();
};

#endif