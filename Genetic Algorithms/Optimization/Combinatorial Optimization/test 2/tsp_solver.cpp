#include "tsp_solver.h"
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>
#include <vector>
#include <cmath>

// Constructor: Initialize nearest neighbors cache
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
      best_distance_(std::numeric_limits<double>::max()) {
    nearest_neighbors_cache_.resize(cities.size());
    for (size_t i = 0; i < cities.size(); ++i) {
        nearest_neighbors_cache_[i] = _get_sorted_nearest_neighbors(i, distances_, 3);
    }
}

// Paso 1: Generar población inicial aleatoria con depósitos aleatorios
void TspSolver::initialize_population() {
    population_.clear();
    start_cities_.clear();
    std::vector<int> all_cities(cities_.size());
    std::iota(all_cities.begin(), all_cities.end(), 0);
    std::random_device rd;
    std::mt19937 g(rd());

    for (int i = 0; i < population_size_; ++i) {
        std::shuffle(all_cities.begin(), all_cities.end(), g);
        int start_city = all_cities[0];
        std::vector<int> tour = {start_city};
        for (size_t j = 1; j < all_cities.size(); ++j) {
            tour.push_back(all_cities[j]);
        }
        population_.push_back(tour);
        start_cities_.push_back(start_city);
    }
}

// Helper: Obtener los m vecinos más cercanos de una ciudad
std::vector<int> TspSolver::_get_sorted_nearest_neighbors(int city, const std::vector<std::vector<double>>& dist_matrix, int m) {
    std::vector<std::pair<double, int>> distance_pairs;
    for (size_t i = 0; i < dist_matrix.size(); ++i) {
        if (static_cast<int>(i) != city) {
            distance_pairs.emplace_back(dist_matrix[city][i], i);
        }
    }
    std::sort(distance_pairs.begin(), distance_pairs.end());
    std::vector<int> neighbors;
    for (size_t i = 0; i < std::min(static_cast<size_t>(m), distance_pairs.size()); ++i) {
        neighbors.push_back(distance_pairs[i].second);
    }
    return neighbors;
}

// Helper: Generar tours candidatos insertando una ciudad cerca de sus vecinos
std::vector<std::vector<int>> TspSolver::_generate_candidate_insertions(int city, const std::vector<int>& base_tour, const std::vector<int>& neighbors) {
    std::vector<std::vector<int>> candidate_tours;
    for (int neighbor : neighbors) {
        auto it = std::find(base_tour.begin(), base_tour.end(), neighbor);
        if (it == base_tour.end()) continue;
        size_t pos = it - base_tour.begin();

        std::vector<int> tour_before = base_tour;
        tour_before.insert(tour_before.begin() + pos, city);
        candidate_tours.push_back(tour_before);

        if (pos + 1 < tour_before.size()) {
            std::vector<int> tour_after = base_tour;
            tour_after.insert(tour_after.begin() + pos + 1, city);
            candidate_tours.push_back(tour_after);
        }
    }
    return candidate_tours;
}

// Heurística de "Remoción de Abruptos"
void TspSolver::remove_abrupts(std::vector<int>& tour, int m) {
    if (tour.empty() || tour.size() == 1) return;

    int start_city = tour[0];
    std::vector<int> current_best_tour = tour;
    double current_best_cost = calculate_fitness(current_best_tour, start_city);

    std::vector<int> cities_to_optimize = current_best_tour;
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(cities_to_optimize.begin() + 1, cities_to_optimize.end(), g);
    size_t max_cities = (cities_to_optimize.size() + 1) / 2;
    cities_to_optimize.resize(std::min(max_cities, cities_to_optimize.size()));

    for (int city_to_move : cities_to_optimize) {
        if (city_to_move == start_city) continue;
        if (std::find(current_best_tour.begin(), current_best_tour.end(), city_to_move) == current_best_tour.end()) {
            continue;
        }

        const std::vector<int>& nearest_neighbors = nearest_neighbors_cache_[city_to_move];

        std::vector<int> temp_tour_base;
        for (int c : current_best_tour) {
            if (c != city_to_move) {
                temp_tour_base.push_back(c);
            }
        }

        std::vector<std::vector<int>> candidate_tours = _generate_candidate_insertions(city_to_move, temp_tour_base, nearest_neighbors);

        for (const auto& candidate_tour : candidate_tours) {
            if (candidate_tour[0] != start_city) continue;
            double cost_candidate = calculate_fitness(candidate_tour, start_city);
            if (cost_candidate < current_best_cost) {
                current_best_tour = candidate_tour;
                current_best_cost = cost_candidate;
            }
        }
    }

    tour = current_best_tour;
}

// Evaluar la aptitud de la población
void TspSolver::evaluate_population() {
    fitness_.resize(population_.size());
    for (size_t i = 0; i < population_.size(); ++i) {
        fitness_[i] = calculate_fitness(population_[i], start_cities_[i]);
        if (fitness_[i] < best_distance_) {
            best_distance_ = fitness_[i];
            best_tour_ = population_[i];
        }
    }
}

// Calcular la aptitud (VFO)
double TspSolver::calculate_fitness(const std::vector<int>& tour, int start_city) {
    if (tour.size() < 2) return 0.0;

    double total_distance = 0.0;
    double total_penalty = 0.0;
    double time = 0.0;

    for (size_t i = 0; i + 1 < tour.size(); ++i) {
        int from = tour[i];
        int to = tour[i + 1];
        double dist = distances_[from][to];
        time += dist;
        double e_i = cities_[to].time_window_start;
        double l_i = cities_[to].time_window_end;
        if (time < e_i) {
            total_penalty += (e_i - time) * (e_i - time);
            time = e_i;
        } else if (time > l_i) {
            total_penalty += (time - l_i) * (time - l_i);
        }
        total_distance += dist;
    }

    int last = tour.back();
    double dist_back = distances_[last][start_city];
    time += dist_back;
    double e_i = cities_[start_city].time_window_start;
    double l_i = cities_[start_city].time_window_end;
    if (time < e_i) {
        total_penalty += (e_i - time) * (e_i - time);
    } else if (time > l_i) {
        total_penalty += (time - l_i) * (time - l_i);
    }
    total_distance += dist_back;

    return total_distance + penalty_factor_ * total_penalty;
}

// Operador de cruce: Cycle Crossover
std::vector<int> TspSolver::cycle_crossover(const std::vector<int>& parent1, const std::vector<int>& parent2, int start_city) {
    std::vector<int> child(parent1.size(), -1);
    child[0] = start_city;

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
        if (start >= static_cast<int>(p1_sub.size())) break;
    }

    std::copy(sub_child.begin(), sub_child.end(), child.begin() + 1);
    return child;
}

// Selección por torneo
std::vector<int> TspSolver::tournament_selection(int tournament_size) {
    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_int_distribution<int> dist(0, population_size_ - 1);

    std::vector<std::pair<int, double>> candidates;
    for (int i = 0; i < tournament_size; ++i) {
        int idx = dist(g);
        candidates.emplace_back(idx, fitness_[idx]);
    }
    auto best = std::min_element(candidates.begin(), candidates.end(),
                                 [](const auto& a, const auto& b) { return a.second < b.second; });
    return population_[best->first];
}

// Generar la próxima generación
void TspSolver::next_generation() {
    std::vector<std::vector<int>> new_population;
    std::vector<int> new_start_cities;
    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_real_distribution<double> prob(0.0, 1.0);
    std::uniform_int_distribution<int> dist(1, population_.size() - 1);

    // Elitismo: Transferir el mejor individuo
    auto best_idx = std::distance(fitness_.begin(), std::min_element(fitness_.begin(), fitness_.end()));
    new_population.push_back(population_[best_idx]);
    new_start_cities.push_back(start_cities_[best_idx]);

    // Generar hijos
    std::vector<std::vector<int>> children;
    std::vector<int> child_start_cities;
    while (children.size() < population_size_ - 1) { // -1 para el elite
        // Selección por torneo
        auto parent1 = tournament_selection(3);
        auto parent2 = tournament_selection(3);
        int idx1 = std::find_if(population_.begin(), population_.end(),
                                [&parent1](const auto& p) { return p == parent1; }) - population_.begin();
        int start_city = start_cities_[idx1];

        std::vector<int> child = parent1; // Por defecto, copiar parent1
        if (prob(g) < crossover_rate_) {
            child = cycle_crossover(parent1, parent2, start_city);
        }

        if (prob(g) < mutation_rate_) {
            if (child.size() > 2) {
                int i = dist(g);
                int j = dist(g);
                while (i == j) j = dist(g);
                std::swap(child[i], child[j]);
            }
        }

        remove_abrupts(child, 3);
        children.push_back(child);
        child_start_cities.push_back(start_city);
    }

    // Combinar elite y hijos
    new_population.insert(new_population.end(), children.begin(), children.end());
    new_start_cities.insert(new_start_cities.end(), child_start_cities.begin(), child_start_cities.end());

    // Actualizar población y evaluar
    population_ = new_population;
    start_cities_ = new_start_cities;
    evaluate_population();

    // Inyección aleatoria: Reemplazar 5% de la población
    int num_random = std::ceil(population_size_ * random_insert_prob_);
    std::uniform_int_distribution<int> pop_dist(0, population_size_ - 1);
    for (int i = 0; i < num_random; ++i) {
        std::vector<int> all_cities(cities_.size());
        std::iota(all_cities.begin(), all_cities.end(), 0);
        std::shuffle(all_cities.begin(), all_cities.end(), g);
        int start_city = all_cities[0];
        std::vector<int> random_tour = {start_city};
        for (size_t j = 1; j < all_cities.size(); ++j) {
            random_tour.push_back(all_cities[j]);
        }
        remove_abrupts(random_tour, 3);
        int idx = pop_dist(g);
        population_[idx] = random_tour;
        start_cities_[idx] = start_city;
    }

    // Reemplazo: Seleccionar los mejores TSP_POP_SIZE individuos
    std::vector<std::pair<std::vector<int>, double>> pop_with_fitness;
    for (size_t i = 0; i < population_.size(); ++i) {
        pop_with_fitness.emplace_back(population_[i], calculate_fitness(population_[i], start_cities_[i]));
    }
    std::sort(pop_with_fitness.begin(), pop_with_fitness.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    population_.clear();
    start_cities_.clear();
    for (size_t i = 0; i < static_cast<size_t>(population_size_); ++i) {
        population_.push_back(pop_with_fitness[i].first);
        start_cities_.push_back(pop_with_fitness[i].first[0]);
    }

    evaluate_population();
}

// Resolver el TSP
void TspSolver::solve() {
    initialize_population();
    for (size_t i = 0; i < population_.size(); ++i) {
        remove_abrupts(population_[i], 3);
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