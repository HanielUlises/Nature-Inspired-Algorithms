#include "tsp_solver.h"
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>
#include <iostream>

// Constructor
TspSolver::TspSolver(const std::vector<std::vector<double>>& distances,
                     const std::vector<City>& cities,
                     int population_size,
                     int generations,
                     double crossover_rate,
                     double mutation_rate,
                     double random_insert_prob,
                     double penalty_factor)
    : distances_(distances), cities_(cities), population_size_(population_size),
      generations_(generations), crossover_rate_(crossover_rate), mutation_rate_(mutation_rate),
      random_insert_prob_(random_insert_prob), penalty_factor_(penalty_factor),
      best_distance_(std::numeric_limits<double>::max()) {}

void TspSolver::initialize_population() {
    population_.clear();
    std::vector<int> base_tour(cities_.size() - 1);
    std::iota(base_tour.begin(), base_tour.end(), 1);
    std::random_device rd;
    std::mt19937 g(rd());

    for (int i = 0; i < population_size_; ++i) {
        std::shuffle(base_tour.begin(), base_tour.end(), g);
        std::vector<int> tour = {0};
        tour.insert(tour.end(), base_tour.begin(), base_tour.end());
        population_.push_back(tour);
    }
}

// Remove abrupt changes by swapping cities that break continuity
void TspSolver::remove_abrupts(std::vector<int>& tour) {
    for (size_t i = 1; i + 1 < tour.size(); ++i) {
        double d1 = distances_[tour[i - 1]][tour[i]];
        double d2 = distances_[tour[i]][tour[i + 1]];
        double d_swap = distances_[tour[i - 1]][tour[i + 1]];

        if (d1 + d2 > d_swap * 1.5) {
            std::swap(tour[i], tour[i + 1]);
        }
    }
}

// Evaluate fitness of all tours
void TspSolver::evaluate_population() {
    fitness_.resize(population_.size());
    for (size_t i = 0; i < population_.size(); ++i) {
        fitness_[i] = calculate_fitness(population_[i]);
        if (fitness_[i] < best_distance_) {
            best_distance_ = fitness_[i];
            best_tour_ = population_[i];
        }
    }
}

// Compute fitness of a single tour
double TspSolver::calculate_fitness(const std::vector<int>& tour) {
    double total_distance = 0.0;
    double time = 0.0;
    for (size_t i = 0; i + 1 < tour.size(); ++i) {
        int from = tour[i];
        int to = tour[i + 1];
        double dist = distances_[from][to];
        time += dist;
        if (time < cities_[to].time_window_start || time > cities_[to].time_window_end) {
            total_distance += penalty_factor_ * dist;
        } else {
            total_distance += dist;
        }
    }
    return total_distance;
}

// Distance without penalty (used in abrupt removal)
double TspSolver::calculate_distance(const std::vector<int>& tour) {
    double distance = 0.0;
    for (size_t i = 0; i + 1 < tour.size(); ++i) {
        distance += distances_[tour[i]][tour[i + 1]];
    }
    return distance;
}

// Cycle Crossover (CX) with New York fixed at index 0
std::vector<int> TspSolver::cycle_crossover(const std::vector<int>& parent1, const std::vector<int>& parent2) {
    std::vector<int> child(parent1.size(), -1);
    child[0] = 0;

    std::vector<int> p1_sub(parent1.begin() + 1, parent1.end());
    std::vector<int> p2_sub(parent2.begin() + 1, parent2.end());
    std::vector<int> sub_child(p1_sub.size(), -1);
    std::vector<bool> visited(p1_sub.size(), false);

    int start = 0;
    bool from_parent1 = true;

    while (!visited[start]) {
        int current = start;
        do {
            sub_child[current] = from_parent1 ? p1_sub[current] : p2_sub[current];
            visited[current] = true;
            int next_val = p2_sub[current];
            current = std::find(p1_sub.begin(), p1_sub.end(), next_val) - p1_sub.begin();
        } while (current != start);
        from_parent1 = !from_parent1;
        start = std::find(visited.begin(), visited.end(), false) - visited.begin();
        if (start >= (int)p1_sub.size()) break;
    }

    std::copy(sub_child.begin(), sub_child.end(), child.begin() + 1);
    return child;
}

// Proceed to next generation
void TspSolver::next_generation() {
    std::vector<std::vector<int>> new_population;
    std::random_device rd;
    std::mt19937 g(rd());

    while ((int)new_population.size() < population_size_) {
        std::uniform_int_distribution<int> dist(0, population_size_ - 1);
        const auto& parent1 = population_[dist(g)];
        const auto& parent2 = population_[dist(g)];
        auto child = cycle_crossover(parent1, parent2);

        std::uniform_real_distribution<double> prob(0.0, 1.0);
        if (prob(g) < mutation_rate_) {
            std::shuffle(child.begin() + 1, child.end(), g); 
        }

        remove_abrupts(child);
        new_population.push_back(child);
    }

    population_ = new_population;
    evaluate_population();

    std::uniform_real_distribution<double> rand_prob(0.0, 1.0);
    if (rand_prob(g) < random_insert_prob_) {
        std::vector<int> random_tour(cities_.size() - 1);
        std::iota(random_tour.begin(), random_tour.end(), 1);
        std::shuffle(random_tour.begin(), random_tour.end(), g);
        std::vector<int> full_tour = {0};
        full_tour.insert(full_tour.end(), random_tour.begin(), random_tour.end());
        population_[rd() % population_size_] = full_tour;
    }
}

// Solve the problem
void TspSolver::solve() {
    initialize_population();
    for (auto& tour : population_) {
        remove_abrupts(tour);
    }
    evaluate_population();

    for (int i = 0; i < generations_; ++i) {
        next_generation();
    }
}

double TspSolver::get_best_distance() const {
    return best_distance_;
}

std::vector<int> TspSolver::get_best_tour() const {
    return best_tour_;
}
