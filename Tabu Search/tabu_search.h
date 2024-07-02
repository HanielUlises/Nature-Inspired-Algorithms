#ifndef TABU_SEARCH_H
#define TABU_SEARCH_H

#include <vector>
#include <deque>
#include <unordered_set>
#include <numeric>
#include <iostream>

struct Object {
    int id;
    double weight;

    bool operator==(const Object& other) const;
};

struct Container {
    int id;
    double capacity;
    std::vector<Object> objects;
    double current_weight;

    Container(int id, double capacity);
    bool addObject(const Object& obj);
    bool removeObject(const Object& obj);

    bool operator==(const Container& other) const;
};

class TabuSearch {
public:
    TabuSearch(const std::vector<Object>& objects, double container_capacity, int max_iterations, int tabu_list_size);
    void solve();
    void printSolution();
    void plotConvergenceGraph();

private:
    std::vector<Object> objects;
    double container_capacity;
    int max_iterations;
    int tabu_list_size;
    double best_evaluation;
    std::vector<Container> best_solution;
    std::vector<std::vector<Container>> solutions;
    std::deque<std::vector<Container>> tabu_list;
    std::unordered_set<size_t> tabu_set;
    std::vector<std::pair<int, int>> object_move_history;
    std::vector<double> evaluation_history;

    void generateInitialSolution();
    double evaluateSolution(std::vector<Container>& solution);
    std::vector<Container> getNeighbourSolution(const std::vector<Container>& current_solution);
    bool isTabu(const std::vector<Container>& solution);
    void updateTabuList(const std::vector<Container>& solution);
    bool aspirationCriterion(const std::vector<Container>& solution);
    void diversify();
    double calculateSolutionCost(const std::vector<Container>& solution);
    size_t hashSolution(const std::vector<Container>& solution) const;
    
};

#endif // TABU_SEARCH_H


