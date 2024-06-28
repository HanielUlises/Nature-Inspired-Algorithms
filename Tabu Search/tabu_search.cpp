#include "tabu_search.h"
#include <algorithm>
#include <random>
#include <limits>
#include <cmath>

// Constructor de Container
Container::Container(int id, double capacity) : id(id), capacity(capacity), current_weight(0) {}

bool Container::addObject(const Object& obj) {
    if (current_weight + obj.weight <= capacity) {
        objects.push_back(obj);
        current_weight += obj.weight;
        return true;
    }
    return false;
}

// Constructor de TabuSearch
TabuSearch::TabuSearch(const std::vector<Object>& objects, double container_capacity, int max_iterations, int tabu_list_size)
    : objects(objects), container_capacity(container_capacity), max_iterations(max_iterations), tabu_list_size(tabu_list_size) {}

void TabuSearch::solve() {
    generateInitialSolution();
    best_solution = solutions[0];

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        auto neighbour = getNeighbourSolution(solutions.back());
        if (!isTabu(neighbour)) {
            solutions.push_back(neighbour);
            if (neighbour.size() < best_solution.size()) {
                best_solution = neighbour;
            }
            updateTabuList(neighbour);
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
    solutions.push_back(initial_solution);
}

void TabuSearch::evaluateSolution(std::vector<Container>& solution) {
    // TODO: implementation
    // No-op for this implementation as we're primarily interested in the number of containers used
}

std::vector<Container> TabuSearch::getNeighbourSolution(const std::vector<Container>& current_solution) {
    std::vector<Container> neighbour = current_solution;
    std::random_device rd;
    std::mt19937 gen(rd());

    // Select two random containers
    std::uniform_int_distribution<> dis(0, neighbour.size() - 1);
    int i = dis(gen);
    int j = dis(gen);

    // Ensure they are different
    while (i == j) {
        j = dis(gen);
    }

    if (!neighbour[i].objects.empty()) {
        // Move a random object from container i to container j
        std::uniform_int_distribution<> obj_dis(0, neighbour[i].objects.size() - 1);
        int obj_idx = obj_dis(gen);
        Object obj = neighbour[i].objects[obj_idx];

        if (neighbour[j].addObject(obj)) {
            neighbour[i].objects.erase(neighbour[i].objects.begin() + obj_idx);
            neighbour[i].current_weight -= obj.weight;

            if (neighbour[i].objects.empty()) {
                neighbour.erase(neighbour.begin() + i);
            }
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
    tabu_list.push_back(solution);
    if (tabu_list.size() > tabu_list_size) {
        tabu_list.erase(tabu_list.begin());
    }
}
