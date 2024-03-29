#include "GeneticAlgorithm.h"
#include "Utils.h"

GeneticAlgorithm::GeneticAlgorithm (int populationSize, int generations ,int tournamentGroupSize, double crossoverProbability){
    this -> populationSize = populationSize;
    this -> generations = generations;
    this -> tournamentGroupSize = tournamentGroupSize;
    this -> crossoverProbability = crossoverProbability;
}

Solution GeneticAlgorithm::perform (int numberOfBits, int low, int high){
    Solution best (numberOfBits, low, high);
    srand(time(NULL));

    // Vector of individuals generations
    std::vector<Solution> currentGen;

    for(int i = 0; i < populationSize; i++){
        currentGen.push_back(Solution(numberOfBits, low, high));
    }

    std::cout << "First generation" << std::endl;
    for(Solution s: currentGen){
        std::cout << s.toString() << std::endl;
    }

    std::cout << "Crossed solutions" << std::endl;

    for(int i = 0; i < generations; i++){
        std::vector<Solution> crossedSolutions = tournamentCrossover (currentGen);
    }

    return best;
}

Solution GeneticAlgorithm::tournamentWinners (std::vector<Solution> const& currentGeneration){
    std::vector<Solution> candidates;

    std::set<int> randomNumbers = randomDistinctNumbers (populationSize, tournamentGroupSize);

    for(int rn : randomNumbers){
        candidates.push_back(currentGeneration[rn]);
    }

    Solution best_candidate = candidates[0];
        double max_fitness = candidates[0].fitness();

        for(Solution s: candidates){
            double fitness = s.fitness();
            if(fitness < max_fitness){
                max_fitness = fitness;
                best_candidate = s;
            }
        }
    return best_candidate;
}

std::vector<Solution> GeneticAlgorithm::tournamentCrossover (std::vector<Solution> const& currentGeneration){
    std::vector<Solution> newSolutions;

    while (newSolutions.size() < populationSize){
        Solution winner_1 = tournamentWinners(currentGeneration);
        Solution winner_2 = tournamentWinners(currentGeneration);

        // Crossover selection
        // Children
        std::vector<Solution> chosen = winner_1.single_point_crossover(winner_2, crossoverProbability);
        newSolutions.push_back(chosen[0]);
        newSolutions.push_back(chosen[1]);
    }

    return newSolutions;
}