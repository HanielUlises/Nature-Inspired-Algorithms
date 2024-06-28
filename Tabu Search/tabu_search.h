#ifndef TABU_SEARCH_H
#define TABU_SEARCH_H

#include <vector>
#include <iostream>

struct Object {
    int id;
    double weight;
};

struct Container {
    int id;
    double capacity;
    std::vector<Object> objects;
    double current_weight;

    Container(int id, double capacity);
    bool addObject(const Object& obj);
};

class TabuSearch {
public:
    TabuSearch(const std::vector<Object>& objects, double container_capacity, int max_iterations, int tabu_list_size);
    void solve();
    void printSolution();

private:
    std::vector<Object> objects;
    double container_capacity;
    int max_iterations;
    int tabu_list_size;
    std::vector<Container> best_solution;
    std::vector<std::vector<Container>> solutions;
    std::vector<std::vector<Container>> tabu_list;

    void generateInitialSolution();
    void evaluateSolution(std::vector<Container>& solution);
    std::vector<Container> getNeighbourSolution(const std::vector<Container>& current_solution);
    bool isTabu(const std::vector<Container>& solution);
    void updateTabuList(const std::vector<Container>& solution);
};

#endif // TABU_SEARCH_H
