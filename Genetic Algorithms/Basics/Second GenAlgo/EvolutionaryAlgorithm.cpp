// EvolutionaryAlgorithm.cpp
#include "EvolutionaryAlgorithm.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>
#include <ctime>
std::random_device rd;
std::mt19937 gen(rd());

EvolutionaryAlgorithm::EvolutionaryAlgorithm(int size, int populationSize, int generations)
    : size(size), populationSize(populationSize), generations(generations) {
    std::srand(static_cast<unsigned int>(std::time(nullptr))); // Seed for random number generation
}

void EvolutionaryAlgorithm::initializePopulation() {
    population.clear();
    for (int i = 0; i < populationSize; ++i) {
        population.push_back(generateIndividual());
    }
}

std::vector<int> EvolutionaryAlgorithm::generateIndividual() {
    std::vector<int> individual(size * size);
    std::iota(individual.begin(), individual.end(), 1);
    
    std::mt19937 rng(std::random_device{}());
    std::shuffle(individual.begin(), individual.end(), rng);
    
    return individual;
}

int EvolutionaryAlgorithm::calculateFitness(const std::vector<int>& individual) {
    int magicConstant = size * (size * size + 1) / 2;
    int fitness = 0;

    // Check rows and columns
    for (int i = 0; i < size; ++i) {
        int rowSum = 0, colSum = 0;
        for (int j = 0; j < size; ++j) {
            rowSum += individual[i * size + j];
            colSum += individual[j * size + i];
        }
        fitness += abs(magicConstant - rowSum) + abs(magicConstant - colSum);
    }

    // Check diagonals
    int diagSum1 = 0, diagSum2 = 0;
    for (int i = 0; i < size; ++i) {
        diagSum1 += individual[i * size + i];
        diagSum2 += individual[(size - i - 1) * size + i];
    }
    fitness += abs(magicConstant - diagSum1) + abs(magicConstant - diagSum2);

    return fitness;
}

std::vector<int> EvolutionaryAlgorithm::selection(std::vector<int>fitnessP) {
    std::vector<int> selectedParents;
    double totalFitness = accumulate(fitnessP.begin(), fitnessP.end(), 0.0);

    std::vector<double> probabilities;
    for (double fitness : fitnessP) {
        probabilities.push_back(fitness / totalFitness);
    }

    std::uniform_real_distribution<double> dis(0.0, 1.0);
    for (int i = 0; i < populationSize; ++i) {
        double r = dis(gen);
        double cumulativeProbability = 0.0;
        for (int j = 0; j < probabilities.size(); ++j) {
            cumulativeProbability += probabilities[j];
            if (cumulativeProbability > r) {
                selectedParents.push_back(j);
                break;
            }
        }
    }
    return selectedParents;
}

std::vector<int> EvolutionaryAlgorithm::tournamentSelection(std::vector<int>&fitnessValues,int tournamentSize) {
    std::vector<int> selectedParents;
    std::uniform_int_distribution<int> dis(0, populationSize - 1);

    for (int i = 0; i < populationSize; ++i) {
        double bestFitness = std::numeric_limits<double>::max();
        int bestIndividual = -1;
        
        // Run a tournament
        for (int j = 0; j < tournamentSize; ++j) {
            int contenderIndex = dis(gen);
            if (fitnessValues[contenderIndex] < bestFitness) {
                bestFitness = fitnessValues[contenderIndex];
                bestIndividual = contenderIndex;
            }
        }

        selectedParents.push_back(bestIndividual);
    }
    return selectedParents;
}
bool isPair(int num){
    if( num % 2 ==0 )return true;
    else return false;
}

int notBusy(std::vector<int> busynum, int numA){
    int boole;
    for(auto num:busynum){
        if(num==numA){
            boole=1;
            break;
        }else{ boole=0;}
    }
    return boole;
}
std::vector<int> EvolutionaryAlgorithm::crossover(std::vector<int>& parent1, std::vector<int>& parent2) {
    int sizeParent = parent1.size();
    std::vector<int> child1;
    std::vector<int> child2;
    std::vector<int> busynum;
    int sizeSubString = 0;
    if(isPair(sizeParent))sizeSubString = (sizeParent-2)/2;
    else sizeSubString = ((sizeParent-1)/2);

    for(size_t i=0; i<sizeSubString+2; i++){
        if(i<2){
            child1.push_back(0);
        }else {
            child1.push_back(parent1[i]);
            busynum.push_back(parent1[i]);
        }
    }
    int j = sizeSubString+2;
    for (size_t i = sizeSubString+2; i < sizeParent; i++){
        
        if(j<sizeParent){
            int allel = parent2[j];
            bool boole=notBusy(busynum,allel);
            if(boole==0){
                child1.push_back(allel);
                busynum.push_back(allel);
            }else{
                i--;
            }
        }else{
            int allel = parent2[j-sizeParent];
            bool boole=notBusy(busynum,allel);
            if(boole==0){
                child1.push_back(allel);
                busynum.push_back(allel);
            }else{
                i--;
            }
        }
        j++;

    }
    int k=0;
    for(size_t i=0; i<2; i++){
        int allel = parent2[k];
        bool boole=notBusy(busynum,allel);
        if(boole==0){
            child1[i]=allel;
        }else{
            i--;
        }
        k++;
    }
    
    return child1;
}

void EvolutionaryAlgorithm::mutation(std::vector<int>& individual) {
    // Placeholder for mutation logic
}

bool EvolutionaryAlgorithm::isMagicSquare(const std::vector<int>& square) {
    // This method utilizes calculateFitness and checks if the fitness is 0, indicating a perfect magic square
    return calculateFitness(square) == 0;
}

void EvolutionaryAlgorithm::solve() {
    initializePopulation();
    
    for (int gen = 0; gen < generations; ++gen) {
        // Selection
        std::vector<int> fitnessP;
        for(const auto& individual : population){
            int fitness = calculateFitness (individual);
            fitnessP.push_back(fitness);
        }
        //tournament selection
        //int tournamentSize = 2 + gen * (10 - 2) / generations;
        //selectionTournament();

        //stochastic selection
        auto selectedParents = selection(fitnessP);

        // Crossover and Mutation
        std::vector<std::vector<int>> children;
        for (int i = 0; i < populationSize; i += 2) {
            children.push_back(crossover(population[selectedParents[i]], population[selectedParents[i+1]]));
            children.push_back(crossover(population[selectedParents[i+1]], population[selectedParents[i]]));
            break;
            mutation(population[i]);
            mutation(population[i+1]);
        }
        break;

        // Find and print solution if one exists in the current generation
        for (const auto& individual : population) {
            if (isMagicSquare(individual)) {
                std::cout << "Magic square found in generation " << gen << ":\n";
                printSolution(individual);
                return;
            }
        }
    }
    std::cout << "Solution not found in " << generations << " generations.\n";
}

void EvolutionaryAlgorithm::printSolution(const std::vector<int>& solution) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            std::cout << solution[i * size + j] << ' ';
        }
        std::cout << std::endl;
    }
}
