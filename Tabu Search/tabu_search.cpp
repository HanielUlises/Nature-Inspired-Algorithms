#include "tabu_search.h"
#include <algorithm>
#include <random>
#include <limits>
#include <cmath>

Container::Container(int id, double capacity) : id(id), capacity(capacity), current_weight(0) {}

bool Container::addObject(const Object& obj) {
    if (current_weight + obj.weight <= capacity) {
        objects.push_back(obj);
        current_weight += obj.weight;
        return true;
    }
    return false;
}

bool Container::removeObject(const Object& obj) {
    auto it = std::find(objects.begin(), objects.end(), obj);
    if (it != objects.end()) {
        current_weight -= obj.weight;
        objects.erase(it);
        return true;
    }
    return false;
}

bool Object::operator==(const Object& other) const {
    return id == other.id && weight == other.weight;
}

bool Container::operator==(const Container& other) const {
    return id == other.id &&
           capacity == other.capacity &&
           current_weight == other.current_weight &&
           objects == other.objects;
}

TabuSearch::TabuSearch(const std::vector<Object>& objects, double container_capacity, int max_iterations, int tabu_list_size)
    : objects(objects), container_capacity(container_capacity), max_iterations(max_iterations), tabu_list_size(tabu_list_size) {}

void TabuSearch::solve() {
    generateInitialSolution();
    best_solution = tabu_list.front();

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        auto neighbour = getNeighbourSolution(tabu_list.back());
        if (!isTabu(neighbour) || aspirationCriterion(neighbour)) {
            tabu_list.push_back(neighbour);
            evaluateSolution(neighbour);
            if (neighbour.size() < best_solution.size()) {
                best_solution = neighbour;
            }
            updateTabuList(neighbour);
        }

        if (iteration % 10 == 0) { // Diversification every 10 iterations
            diversify();
        }
    }
}

void TabuSearch::printSolution() {
    std::cout << "Best solution found:\n";
    for (const auto& container : best_solution) {
        std::cout << "Container " << container.id << " (current weight: " << container.current_weight << "): ";
        for (const auto& obj : container.objects) {
            std::cout << "Object " << obj.id << " (weight: " << obj.weight << ") ";
        }
        std::cout << std::endl;
    }
    std::cout << "Total containers used: " << best_solution.size() << std::endl;
}

void TabuSearch::generateInitialSolution() {
    std::vector<Container> initial_solution;
    initial_solution.emplace_back(1, container_capacity);

    for (const auto& obj : objects) {
        bool placed = false;
        for (auto& container : initial_solution) {
            if (container.addObject(obj)) {
                placed = true;
                break;
            }
        }
        if (!placed) {
            Container new_container(initial_solution.size() + 1, container_capacity);
            new_container.addObject(obj);
            initial_solution.push_back(new_container);
        }
    }
    tabu_list.push_back(initial_solution);
}

void TabuSearch::evaluateSolution(std::vector<Container>& solution) {
    double solution_cost = calculateSolutionCost(solution);
    std::cout << solution_cost << std::endl;
}

std::vector<Container> TabuSearch::getNeighbourSolution(const std::vector<Container>& current_solution) {
    std::vector<Container> neighbour = current_solution;
    std::random_device rd;
    std::mt19937 gen(rd());

    // Edge case where there's only one container
    if (neighbour.size() < 2) return neighbour;

    // Select two random containers
    std::uniform_int_distribution<> dis(0, neighbour.size() - 1);
    int i = dis(gen);
    int j = dis(gen);

    // These random containers must be different
    while (i == j) {
        j = dis(gen);
    }

    if (!neighbour[i].objects.empty()) {
        // Move a random object from container i to container j
        std::uniform_int_distribution<> obj_dis(0, neighbour[i].objects.size() - 1);
        int obj_idx = obj_dis(gen);
        Object obj = neighbour[i].objects[obj_idx];

        if (neighbour[j].addObject(obj)) {
            neighbour[i].removeObject(obj);

            if (neighbour[i].objects.empty()) {
                neighbour.erase(neighbour.begin() + i);
            }

            object_move_history.push_back({obj.id, neighbour[j].id});
        }
    }
    return neighbour;
}

bool TabuSearch::isTabu(const std::vector<Container>& solution) {
    for (const auto& tabu_solution : tabu_list) {
        if (tabu_solution == solution) {
            return true;
        }
    }
    return false;
}

void TabuSearch::updateTabuList(const std::vector<Container>& solution) {
    if (tabu_list.size() >= static_cast<size_t>(tabu_list_size)) {
        tabu_list.pop_front();
    }
    tabu_list.push_back(solution);
}

bool TabuSearch::aspirationCriterion(const std::vector<Container>& solution) {
    // Allow moves that improve the best-known solution
    return solution.size() < best_solution.size();
}

void TabuSearch::diversify() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, objects.size() - 1);

    int random_obj_idx = dis(gen);
    Object random_obj = objects[random_obj_idx];

    std::vector<Container> diversified_solution = best_solution;
    for (auto& container : diversified_solution) {
        container.removeObject(random_obj);
    }

    bool placed = false;
    for (auto& container : diversified_solution) {
        if (container.addObject(random_obj)) {
            placed = true;
            break;
        }
    }

    if (!placed) {
        Container new_container(diversified_solution.size() + 1, container_capacity);
        new_container.addObject(random_obj);
        diversified_solution.push_back(new_container);
    }

    tabu_list.push_back(diversified_solution);
}

double TabuSearch::calculateSolutionCost(const std::vector<Container>& solution) {
    double weight_variance = 0.0;
    double mean_weight = std::accumulate(solution.begin(), solution.end(), 0.0,
                                         [](double sum, const Container& container) { return sum + container.current_weight; }) / solution.size();

    for (const auto& container : solution) {
        weight_variance += std::pow(container.current_weight - mean_weight, 2);
    }

    weight_variance /= solution.size();

    double penalty_factor = 1.0;  // Adjust this factor to balance container usage and weight balance
    return solution.size() + penalty_factor * weight_variance;
}