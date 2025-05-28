#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>

class Solution {
private:
    std::vector<int> route;
    double fitness;

public:
    Solution(int numCities);
    std::vector<int>& getRoute();
    void setFitness(double fit);
    double getFitness() const;
    void swapCities(int i, int j);
};

#endif