#pragma once

#include <vector>
#include <string>

class Solution{
    public:
    Solution(int numberOfBits, int low, int high);
    Solution(int numberOfBits, int low, int high, std::vector<int> bits);

    std::string toString();
    double bitsToDouble ();
    // Restriction  function to be optimized
    double fitness();
    std::vector<Solution> single_point_crossover (Solution other, double crossoverProbability);

    private:
    int numberOfBits, low, high;
    std::vector<int> bits; // vector of chromosomes

};