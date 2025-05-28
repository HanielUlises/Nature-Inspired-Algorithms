#include "tsp_solver.h"
#include <algorithm>
#include <random>
#include <numeric>
#include <limits>
#include <utility>
#include <ctime>

TspSolver::TspSolver(const std::vector<City>& cities, const std::vector<std::vector<double>>& distances)
    : cities_(cities), distances_(distances), best_route_time_(std::numeric_limits<double>::infinity()) {
    std::random_device rd;
    gen_.seed(rd());
}

double TspSolver::evaluate_route(const std::vector<int>& route) const {
    if (route.empty() || route.size() != cities_.size()) {
        return std::numeric_limits<double>::infinity();
    }
    
    double time = 0.0;
    if (route[0] != 0) {
        return std::numeric_limits<double>::infinity();
    }
    
    for (size_t i = 0; i < route.size() - 1; ++i) {
        int from = route[i];
        int to = route[i + 1];
        
        if (from < 0 || from >= static_cast<int>(cities_.size()) || 
            to < 0 || to >= static_cast<int>(cities_.size())) {
            return std::numeric_limits<double>::infinity();
        }
        
        time += distances_[from][to];
        
        if (time < cities_[to].window_start) {
            time = cities_[to].window_start;
        }
        
        if (time > cities_[to].window_end) {
            return std::numeric_limits<double>::infinity();
        }
    }
    
    time += distances_[route.back()][route[0]];
    
    return time;
}

std::vector<int> TspSolver::generate_random_route() const {
    std::vector<int> route(cities_.size());
    std::iota(route.begin(), route.end(), 0);
    std::shuffle(route.begin() + 1, route.end(), gen_);
    return route;
}

void TspSolver::mutate(std::vector<int>& route) {
    if (route.size() <= 2) return; 
    
    std::uniform_int_distribution<> dist(1, static_cast<int>(route.size()) - 1);
    int i = dist(gen_);
    int j = dist(gen_);
    std::swap(route[i], route[j]);
}

std::vector<int> TspSolver::cycle_crossover(const std::vector<int>& parent1, const std::vector<int>& parent2) {
    int n = static_cast<int>(parent1.size());
    std::vector<int> child(n, -1);

    std::vector<bool> visited(n, false);
    int index = 0;
    while (!visited[index]) {
        child[index] = parent1[index];
        visited[index] = true;
        int val = parent2[index];
        auto it = std::find(parent1.begin(), parent1.end(), val);
        if (it != parent1.end()) {
            index = std::distance(parent1.begin(), it);
        } else {
            break;
        }
    }

    for (int i = 0; i < n; ++i) {
        if (child[i] == -1) child[i] = parent2[i];
    }

    return child;
}

std::vector<int> TspSolver::remove_abrupts(const std::vector<int>& individual) {
    std::vector<int> new_individual = individual;
    int n = static_cast<int>(new_individual.size());

    if (n < 3) return new_individual; 

    for (int i = 1; i < n - 1; ++i) {
        int prev = new_individual[i - 1];
        int curr = new_individual[i];
        int next = new_individual[i + 1];

        double dist_prev_curr = distances_[prev][curr];
        double dist_curr_next = distances_[curr][next];
        double dist_prev_next = distances_[prev][next];

        if (dist_prev_curr + dist_curr_next > dist_prev_next) {
            std::swap(new_individual[i], new_individual[i + 1]);
        }
    }
    return new_individual;
}

void TspSolver::initialize_population(int population_size) {
    population_.clear();
    population_.reserve(population_size);
    for (int i = 0; i < population_size; ++i) {
        population_.push_back(generate_random_route());
    }
}

void TspSolver::apply_remove_abrupts_to_population() {
    for (auto& individual : population_) {
        individual = remove_abrupts(individual);
    }
}

int TspSolver::select_parent_index(int max_index) {
    std::uniform_int_distribution<> dist(0, max_index);
    return dist(gen_);
}

void TspSolver::select_survivors() {
    std::vector<std::pair<double, std::vector<int>>> fitness_individuals;
    
    for (const auto& individual : population_) {
        double fitness = evaluate_route(individual);
        fitness_individuals.push_back({fitness, individual});
    }
    
    std::sort(fitness_individuals.begin(), fitness_individuals.end());
    
    population_.clear();
    for (size_t i = 0; i < std::min(fitness_individuals.size(), population_.capacity()); ++i) {
        population_.push_back(fitness_individuals[i].second);
    }
}

void TspSolver::next_generation(double mutation_rate, double random_injection_prob) {
    std::vector<std::vector<int>> new_population = population_;
    
    std::uniform_real_distribution<> mut_dist(0.0, 1.0);
    std::uniform_real_distribution<> inj_dist(0.0, 1.0);
    
    while (new_population.size() < population_.capacity()) {
        int parent1_idx = select_parent_index(static_cast<int>(population_.size()) - 1);
        int parent2_idx = select_parent_index(static_cast<int>(population_.size()) - 1);
        const std::vector<int>& parent1 = population_[parent1_idx];
        const std::vector<int>& parent2 = population_[parent2_idx];

        std::vector<int> child = cycle_crossover(parent1, parent2);

        if (mut_dist(gen_) < mutation_rate) {
            mutate(child);
        }

        child = remove_abrupts(child);
        new_population.push_back(child);
    }
    
    if (inj_dist(gen_) < random_injection_prob && !new_population.empty()) {
        std::uniform_int_distribution<> replace_dist(0, static_cast<int>(new_population.size()) - 1);
        int replace_idx = replace_dist(gen_);
        new_population[replace_idx] = generate_random_route();
        new_population[replace_idx] = remove_abrupts(new_population[replace_idx]);
    }
    
    population_ = std::move(new_population);
    select_survivors();
}

void TspSolver::run(int population_size, int generations, double mutation_rate, double random_injection_prob) {
    initialize_population(population_size);
    apply_remove_abrupts_to_population();

    for (const auto& individual : population_) {
        double fitness = evaluate_route(individual);
        if (fitness < best_route_time_) {
            best_route_time_ = fitness;
            best_route_ = individual;
        }
    }

    for (int gen = 0; gen < generations; ++gen) {
        next_generation(mutation_rate, random_injection_prob);

        for (const auto& individual : population_) {
            double fitness = evaluate_route(individual);
            if (fitness < best_route_time_) {
                best_route_time_ = fitness;
                best_route_ = individual;
            }
        }
    }
}

std::vector<int> TspSolver::get_best_route() const {
    if (best_route_.empty()) {
        std::vector<int> default_route(cities_.size());
        std::iota(default_route.begin(), default_route.end(), 0);
        return default_route;
    }
    return best_route_;
}

double TspSolver::get_best_route_time() const {
    return best_route_time_;
}