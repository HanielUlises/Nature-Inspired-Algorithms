#ifndef ACO_H
#define ACO_H

#include <vector>
#include <functional>

class ACO {
public:
    
    int num_ants;
    int max_iterations;
    std::vector<std::vector<int>> distance;
    std::vector<std::vector<int>> flows;
    std::vector<std::vector<double>> pheromone;
    double evaporation; 
    int alpha; 
    int beta; 
    int q;

    ACO(int num_ants, int max_iterations, std::vector<std::vector<int>> distance,std::vector<std::vector<int>> flows, double evaporation, int alpha, int beta, int q);
    int selectNext_city(int actual, std::vector<int> busy) ;
    std::vector<std::vector<int>> constructRoute(int num_city);
    std::vector<int> objectiveFunction(std::vector<std::vector<int>> routes);
    void updatePheromone(std::vector<std::vector<int>> routes, std::vector<int> distances);
    std::vector<int> optimize();
    void printSolution(std::vector<int> route);
};



#endif