#include "GeneticAlgorithm.h"

GeneticAlgorithm::GeneticAlgorithm (int populationSize){
    this -> populationSize = populationSize;
}

Solution GeneticAlgorithm::perform (int numberOfBits, int low, int high){
    Solution best (numberOfBits, low, high);

    // 1st generation of the population
    std::vector<Solution> currentGen;

    for(int i = 0; i < populationSize; i++){
        currentGen.push_back(Solution(numberOfBits, low, high));
    }

    for(Solution s: currentGen){
        std::cout << s.toString() << std::endl;
    }

    return best;
}