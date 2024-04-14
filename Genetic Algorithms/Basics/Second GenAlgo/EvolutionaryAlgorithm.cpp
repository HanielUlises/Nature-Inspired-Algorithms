// EvolutionaryAlgorithm.cpp
#include "EvolutionaryAlgorithm.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>
#include <ctime>
#include "ploteado.h"
std::random_device rd;
std::mt19937 gen(rd());

EvolutionaryAlgorithm::EvolutionaryAlgorithm(int size, int populationSize, int generations, double cross_rate, double  mut_rate)
    : size(size), populationSize(populationSize), generations(generations), cross_rate(cross_rate),mut_rate(mut_rate){
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
    // n * (n²+1)/2
    int magicConstant = (size * (size * size + 1)) / 2;
    int fitness = 0;

    // Check rows and columns
    int cont=0;
    for (int i = 0; i < size; ++i) {
        int rowSum = 0, colSum = 0;
        for (int j = 0; j < size; ++j) {
                rowSum += individual[i * size + j];
                colSum += individual[j * size + i];
        }
        fitness += abs(magicConstant - rowSum) + abs(magicConstant - colSum);
        cont++;
    }

    // Check diagonals
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

std::vector<int> EvolutionaryAlgorithm::tournamentSelection(std::vector<int>fitnessValues,int tournamentSize) {
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

std::vector<int> EvolutionaryAlgorithm::crossoverOX(std::vector<int>& parent1, std::vector<int>& parent2) {
    int sizeParent = parent1.size();
    std::vector<int> child1;
    std::vector<int> child2;
    std::vector<int> busynum;
    int sizeSubString = 0;
    if(isPair(sizeParent))sizeSubString = (sizeParent-2)/2;
    else sizeSubString = ((sizeParent-1)/2);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

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

int chooseCrossingPoint(int size, int busy){
    std::uniform_int_distribution dis(0, size-1);
    bool keep=true;
    int point;
    while(keep){
        point = dis(gen);
        if (point!=busy)break;
    }
    return point;
}

int position(std::vector<int> child, int numA){
    int cont=0;
    for(auto num:child){
        if(num==numA){
            return cont;
            break;
        }else{ cont++;}
    }
    return cont;
}

bool isBusyPMX(std::vector<int> busynum, int numA){
    bool boole;
    for(auto num:busynum){
        if(num==numA){
            boole=true;
            return boole;
            break;
        }else{boole=false;}
    }
    return boole;
}

std::vector<std::vector<int>> EvolutionaryAlgorithm::crossoverPMX(const std::vector<int>& parent1, const std::vector<int>& parent2) {
    std::vector<std::vector<int>> children;
    int size = parent1.size();
    std::vector<int> child1(size), child2(size);
    std::unordered_map<int, int> mapping1, mapping2;

    int point1 = chooseCrossingPoint(size, size);
    int point2 = chooseCrossingPoint(size, point1);

    if (point2 < point1) std::swap(point1, point2);

    // Create the mapping and construct the middle segment of the children
    for (int i = point1; i <= point2; ++i) {
        child1[i] = parent2[i];
        child2[i] = parent1[i];
        mapping1[parent2[i]] = parent1[i];
        mapping2[parent1[i]] = parent2[i];
    }

    // Fill in the rest of the slots in the children
    for (int i = 0; i < size; ++i) {
        if (i < point1 || i > point2) {
            int nextParent1 = parent1[i], nextParent2 = parent2[i];

            while (mapping1.count(nextParent1)) nextParent1 = mapping1[nextParent1];
            while (mapping2.count(nextParent2)) nextParent2 = mapping2[nextParent2];

            child1[i] = nextParent1;
            child2[i] = nextParent2;
        }
    }

    children.push_back(child1);
    children.push_back(child2);
    return children;
}

std::vector<int> EvolutionaryAlgorithm::mutationIn(std::vector<int>& individual) {
    double mutationChance = std::uniform_real_distribution<double>(0.0, 1.0)(gen);
    if (mutationChance < mut_rate) {
        int i = std::uniform_int_distribution<int>(0, individual.size() - 1)(gen);
        int j;
        do {
            j = std::uniform_int_distribution<int>(0, individual.size() - 1)(gen);
        } while (i == j);
        std::swap(individual[i], individual[j]);
    }
    return individual;
}


std::vector<int> EvolutionaryAlgorithm::mutationDes(std::vector<int> individual) {

    int sizeChild = individual.size();
    std::uniform_real_distribution<double> dis3(0.0, 1.0);
    std::uniform_int_distribution dis(1, sizeChild);
    std::uniform_int_distribution dis2(1, sizeChild-1);
    int num_mutation = dis(gen);
    int i,j;


    for(size_t k=0; k<num_mutation; k++){
        if(dis3(gen)<mut_rate){
            i = dis2(gen);
            j = dis2(gen);
            if(i<j){
                int aux = individual[j]; //elemento q voy a mover
                int dif = j-i;
                for(size_t w=i; w<=j;w++){
                    int aux2=individual[w];
                    individual[w]=aux;
                    aux=aux2;
                }
            }else if(i>j){
                int aux = individual[i]; //elemento q voy a mover
                int dif = i-j;
                for(size_t w=j; w<=i;w++){
                    int aux2=individual[w];
                    individual[w]=aux;
                    aux=aux2;
                }
            }else{
                k--;
            }
            
        }
    }
    return individual;

}

std::vector<std::vector<int>> EvolutionaryAlgorithm::elitism(std::vector<std::vector<int>> population){
    std::vector<int> fitnessP;
    for(const auto& individual : population){
        int fitness = calculateFitness (individual);
        fitnessP.push_back(fitness);
    }
        // This vector  will hold individuals along with their fitness
    std::vector<IndividualWithFitness> sortedPopulation;
    sortedPopulation.reserve(population.size());

    // Vector with individuals and their corresponding fitness
    for (size_t i = 0; i < population.size(); ++i) {
        sortedPopulation.emplace_back(population[i], fitnessP[i]);
    }

    // Sorting the combined vector based on fitness in descending order
    // NNE
    std::sort(sortedPopulation.begin(), sortedPopulation.end(), [](const IndividualWithFitness& a, const IndividualWithFitness& b) {
        return a.fitness < b.fitness;
    });

    population.clear();
    fitnessP.clear();

    // Putting the top individuals back into the population until the desired population size is reached
    for (size_t i = 0; i < populationSize && i < sortedPopulation.size(); ++i) {
        population.push_back(sortedPopulation[i].individual);
        fitnessP.push_back(sortedPopulation[i].fitness);
    }
    return population;
}

bool EvolutionaryAlgorithm::isMagicSquare(const std::vector<int>& square) {
    // This method utilizes calculateFitness and checks if the fitness is 0, indicating a perfect magic square
    return calculateFitness(square) == 0;
}

int EvolutionaryAlgorithm::maximizeSuccessCount(const std::vector<int>& square) {
    int successCount = 0;
    int magicConstant = size * (size * size + 1) / 2;

    // Row and column sum
    for (int i = 0; i < size; ++i) {
        int rowSum = 0, colSum = 0;
        for (int j = 0; j < size; ++j) {
            rowSum += square[i * size + j];
            colSum += square[j * size + i];
        }
        if (rowSum == magicConstant) ++successCount;
        if (colSum == magicConstant) ++successCount;
    }

    // Diagonal sum
    int diagSum1 = 0, diagSum2 = 0;
    for (int i = 0; i < size; ++i) {
        diagSum1 += square[i * size + i];
        diagSum2 += square[i * size + (size - i - 1)];
    }
    if (diagSum1 == magicConstant) ++successCount;
    if (diagSum2 == magicConstant) ++successCount;

    return successCount;
}

int EvolutionaryAlgorithm::minimizeMagicConstantError(const std::vector<int>& square) {
    int error = 0;
    int magicConstant = size * (size * size + 1) / 2;

    // Note: Formula based error
    // Row and column error
    for (int i = 0; i < size; ++i) {
        int rowSum = 0, colSum = 0;
        for (int j = 0; j < size; ++j) {
            rowSum += square[i * size + j];
            colSum += square[j * size + i];
        }
        error += std::abs(magicConstant - rowSum);
        error += std::abs(magicConstant - colSum);
    }

    // Diagonal error
    int diagSum1 = 0, diagSum2 = 0;
    for (int i = 0; i < size; ++i) {
        diagSum1 += square[i * size + i];
        diagSum2 += square[i * size + (size - i - 1)];
    }
    error += std::abs(magicConstant - diagSum1);
    error += std::abs(magicConstant - diagSum2);

    return error;
}

void EvolutionaryAlgorithm::solve() {
    initializePopulation();

    std::vector<double> averageFitnessHistory;
    std::vector<double> bestFitnessHistory;
    std::vector<double> worstFitnessHistory;
    
    for (int gener = 0; gener < generations; ++gener) {
        std::vector<IndividualWithFitness> evaluatedPopulation;
        std::vector<int> fitnessP;
        for(const auto& individual : population){
            int fitness = calculateFitness(individual);
            fitnessP.push_back(fitness);
            if (fitness == 0) { // Check immediately when fitness is calculated
                std::cout << "Magic square found in generation " << gener << ":\n";
                printSolution(individual);
                plotConvergenceGraph(averageFitnessHistory, bestFitnessHistory, worstFitnessHistory);
                return;
            }
        }
        
        int tournamentSize = 2;
        auto selectedParents = tournamentSelection(fitnessP, tournamentSize);

        std::vector<std::vector<int>> children;
        for (int i = 0; i < populationSize; i += 2) {
            if (i+1 < selectedParents.size()) {
                auto childrenaux = crossoverPMX(population[selectedParents[i]], population[selectedParents[i+1]]);
                for (auto& child : childrenaux) {
                    child = mutationIn(child);
                    children.push_back(child);
                }
            }
        }

        for(auto& child : children){
            population.push_back(child);
        }
        
        population = elitism(population);

        // Recalculate fitness for stats and check for magic squares again (in case mutations created one)
        fitnessP.clear();
        for(const auto& individual : population){
            int fitness = calculateFitness(individual);
            fitnessP.push_back(fitness);
            if (fitness == 0) {
                std::cout << "Magic square found in generation " << gener << ":\n";
                printSolution(individual);
                plotConvergenceGraph(averageFitnessHistory, bestFitnessHistory, worstFitnessHistory);
                return;
            }
        }

        // Calculate statistics for the generation
        double averageFitness = std::accumulate(fitnessP.begin(), fitnessP.end(), 0.0) / fitnessP.size();
        double bestFitness = *std::min_element(fitnessP.begin(), fitnessP.end());
        double worstFitness = *std::max_element(fitnessP.begin(), fitnessP.end());

        std::cout << "Best: " << bestFitness << "|" << " Average: " << averageFitness  << "|" << " Worst: " << worstFitness << std::endl;

        averageFitnessHistory.push_back(averageFitness);
        bestFitnessHistory.push_back(bestFitness);
        worstFitnessHistory.push_back(worstFitness);
    }
    
    std::cout << "Solution not found in " << generations << " generations.\n";
    plotConvergenceGraph(averageFitnessHistory, bestFitnessHistory, worstFitnessHistory);
}

void EvolutionaryAlgorithm::printSolution(const std::vector<int>& solution) {
    
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            std::cout << solution[i * size + j] << ' ';
        }
        std::cout << std::endl;
    }
}
