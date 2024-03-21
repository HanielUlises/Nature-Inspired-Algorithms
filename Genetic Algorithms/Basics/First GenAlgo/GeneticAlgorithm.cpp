#include "GeneticAlgorithm.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <iostream>
#include <numeric>
#include "encode.h"
#include "decode.h"

// #define M_PI 3.14159265358979323846
// #define M_E 2.71828182845904523536

std::random_device rd;
std::mt19937 gen(rd());

GeneticAlgorithm::GeneticAlgorithm(int populationSize, int numberOfGenerations, double crossoverRate, double mutationRate)
    : populationSize(populationSize), numberOfGenerations(numberOfGenerations),
      crossoverRate(crossoverRate), mutationRate(mutationRate) {
}

GeneticAlgorithm::~GeneticAlgorithm() {
}

// Perfomance of the genetic algorithm on ℝ_0.01
void GeneticAlgorithm::real_performance(int option) {
    std::vector<std::vector<double>> child;
    initializePopulation(option);
    std::function<double(const std::vector<double>&)> objectiveFunction;
    if (option==1){
        objectiveFunction=std::bind(&GeneticAlgorithm::rosenbrockFunction, this, std::placeholders::_1);
    }else if (option==2){
        objectiveFunction=std::bind(&GeneticAlgorithm::ackleyFunction, this, std::placeholders::_1);    
    }else{
        std::cout<<"cuak"<<std::endl;
    }

    for (int i = 0; i < numberOfGenerations; ++i) {
        
        // Functions test
        evaluateFitness(objectiveFunction,1,0.0f,0.0f);
        auto selectedParents = selection();
        //se muere
        child=crossover(selectedParents);
        child=mutation(child); //child mutados
        
        for (size_t i = 0; i < child.size(); i++){
            population.push_back(child[i]);
        }
        evaluateFitness(objectiveFunction,1,0.0f,0.0f);
        elitismParents();

        double bestFitness = *std::min_element(fitnessValues.begin(), fitnessValues.end());
        double worstFitness = *std::max_element(fitnessValues.begin(), fitnessValues.end());
        double averageFitness = accumulate(fitnessValues.begin(), fitnessValues.end(), 0.0) / fitnessValues.size();

        bestFitnessHistory.push_back(bestFitness);
        worstFitnessHistory.push_back(worstFitness);
        averageFitnessHistory.push_back(averageFitness);

        if (shouldStop(i, option)) break;
        std::cout<<i<<std::endl;
    }
    if (option==1){
        plotConvergenceGraph("Rosenbrock function");
    }else if (option==2){
        plotConvergenceGraph("Ackley function");
    }
}

// Binary codification for the genetic algorithm
void GeneticAlgorithm::binary_performance(const int option) {
    std::vector<std::vector<std::string>> child;
    initializePopulationBinary(option);
    std::function<double(const std::vector<double>&)> objectiveFunction;
    
    double liminf;
    double limsup;

    // Bind the appropriate objective function and set limits based on the option chosen
    if (option==1){
        objectiveFunction=std::bind(&GeneticAlgorithm::rosenbrockFunction, this, std::placeholders::_1);
        liminf=-2.048;
        limsup=2.048;
    }else if (option==2){
        objectiveFunction=std::bind(&GeneticAlgorithm::ackleyFunction, this, std::placeholders::_1);    
        liminf=-32.768;
        limsup=32.768;
    }else{
        std::cout<<"Not an actual option"<<std::endl;
    }

    for (int i = 0; i < numberOfGenerations; ++i) {
        // Evaluate fitness for the current population
        evaluateFitness(objectiveFunction,2,liminf,limsup);    
        
        // Selection process, here using tournament selection to increase selective pressure
        int tournamentSize = 2 + i * (10 - 2) / numberOfGenerations; // Gradually increase tournament size
        auto selectedParents = tournamentSelection(tournamentSize);
        
        // Crossover among selected parents
        child = crossover_binary(selectedParents);
        
        // Mutate the offspring produced by crossover
        // Mutate based on dynamic mutation rate
        child = mutation_binary(child, i); 
        
        // Add children to the current population
        for (size_t j = 0; j < child.size(); ++j){
            populationBinary.push_back(child[j]);
        }
        
        // Evaluate the fitness of the new population, including offspring
        evaluateFitness(objectiveFunction,2,liminf,limsup);
        
        // Apply elitism to select the best individuals for the next generation
        elitismParentsBinary();

        // Record fitness statistics for the current generation
        double bestFitness = *std::min_element(fitnessValues.begin(), fitnessValues.end());
        double worstFitness = *std::max_element(fitnessValues.begin(), fitnessValues.end());
        double averageFitness = accumulate(fitnessValues.begin(), fitnessValues.end(), 0.0) / fitnessValues.size();

        bestFitnessHistory.push_back(bestFitness);
        worstFitnessHistory.push_back(worstFitness);
        averageFitnessHistory.push_back(averageFitness);

        // Check if the stopping criteria have been met
        if (shouldStop(i, option)) break;

        // Output current generation number
        std::cout << "Generation: " << i << std::endl;
    }

    // Plot the convergence graph at the end of the run
    if (option==1){
        plotConvergenceGraph("Rosenbrock function (binary)");
    }else if (option==2){
        plotConvergenceGraph("Ackley function (binary)");
    }
}

int countDecimalPlaces(double number) {
    // Like stod in C but inverted jiji
    std::string numberAsString = std::to_string(number);
    // Trim trailing zeros
    numberAsString.erase(numberAsString.find_last_not_of('0') + 1, std::string::npos); 
    auto decimalPos = numberAsString.find('.');

     // No decimal point found
    if (decimalPos == std::string::npos) 
        return 0;

    return numberAsString.length() - decimalPos - 1; 
}

int bitsNeeded(double x, double y, double& range) {
    // Determine the maximum number of decimal places between the two numbers
    int maxDecimals = std::max(countDecimalPlaces(x), countDecimalPlaces(y));
    
    // Formulae implemented
    // log2(upLim * 10^dec - lowLim * 10^dec)
    double scaledX = x * std::pow(10, maxDecimals);
    double scaledY = y * std::pow(10, maxDecimals);
    
    // Difference between the scaled numbers
    range = std::abs(scaledY - scaledX);
    return std::ceil(std::log2(range));
}

std::vector<double> decodeAllele(std::vector<std::string> bin_string, double lim,double limS){
    std::vector<double> allel;
    int aux;
    for (std::string indivi : bin_string){
        aux = binToDec(indivi);
        int maxDecimals = countDecimalPlaces(lim);
        
        // Formulae implemented
        // log2(upLim * 10^dec - lowLim * 10^dec)
        double scaledX = lim * std::pow(10, maxDecimals);
        
        double real = scaledX +aux;
        double realReduce = real * std::pow(10, -3);
        if(realReduce>limS){
            realReduce=limS;
        }
        allel.push_back(realReduce);
    }
    return allel;
}

void GeneticAlgorithm::initializePopulationBinary(int option){
    populationBinary.clear();
    if (option==1){
        double range = 0;
        int totalBitsNeeded = bitsNeeded(-2.048, 2.048, range);

        // Conversion within(-2.048,range);
        std::cout << "Total number of bits needed: " << totalBitsNeeded << std::endl;
    
        std::uniform_int_distribution<> dis(0, 1);
        int numGenes = 10;

        for (int i = 0; i < populationSize; ++i) { 
            std::vector<std::string> allel;
            for (int j = 0; j < numGenes; ++j) {
                std::string individual;
                for (int k = 0; k < totalBitsNeeded; k++){
                    int gene = dis(gen);
                    std::string s = std::to_string(gene);
                    individual+=s;
                }
                allel.push_back(individual);
            }
            populationBinary.push_back(allel);
        }
    }else if (option==2){
        double range = 0;
        int totalBitsNeeded = bitsNeeded(-32.768, 32.768, range);

        // Conversion within(-32.768,range);
        std::cout << "Total number of bits needed: " << totalBitsNeeded << std::endl;
    
        std::uniform_int_distribution<> dis(0, 1);
        int numGenes = 10;

        for (int i = 0; i < populationSize; ++i) {
            std::vector<std::string> allel;
            for (int j = 0; j < numGenes; ++j) {
                std::string individual;
                for (int k = 0; k < totalBitsNeeded; k++){
                    int gene = dis(gen);
                    std::string s = std::to_string(gene);
                    individual+=s;
                }
                allel.push_back(individual);
            }
            populationBinary.push_back(allel);
        }
    }else{
        std::cout<<"Esta raro nah"<<std::endl;
    }
}

void GeneticAlgorithm::initializePopulation(int option) {
    population.clear();

    if (option==1){
        std::uniform_real_distribution<double> dis(-2.048, 2.048);
        int numGenes = 10;

        for (int i = 0; i < populationSize; ++i) {
            std::vector<double> individual;
            for (int j = 0; j < numGenes; ++j) {
                double gene = dis(gen);
                individual.push_back(gene);
            }
            population.push_back(individual);
        }
    }else if (option==2){
        std::uniform_real_distribution<double> dis(-32.768, 32.768);
        int numGenes = 10;

        for (int i = 0; i < populationSize; ++i) {
            std::vector<double> individual;
            for (int j = 0; j < numGenes; ++j) {
                double gene = dis(gen);
                individual.push_back(gene);
            }
            population.push_back(individual);
        }
    }else{
        std::cout<<"Esta raro na"<<std::endl;
    }

}

void GeneticAlgorithm::elitismParents() {

    std::vector<double> temp;
    double temp2; 
    for(int i=0;i<population.size();i++){
        for(int q=i+1;q<population.size();q++){
            if(fitnessValues[i]<fitnessValues[q]){
                temp=population[q]; 
                population[q]=population[i]; 
                population[i]=temp;
                temp2=fitnessValues[q];
                fitnessValues[q]=fitnessValues[i];
                fitnessValues[i]=temp2;
            } 
        }
    }
    while (population.size()>100)
    {
       population.erase(population.begin());
    }

    population.resize(populationSize);
    fitnessValues.resize(populationSize);
}


void GeneticAlgorithm::elitismParentsBinary() {
    // This vector  will hold individuals along with their fitness
    std::vector<IndividualWithFitness> sortedPopulation;
    sortedPopulation.reserve(populationBinary.size());

    // Vector with individuals and their corresponding fitness
    for (size_t i = 0; i < populationBinary.size(); ++i) {
        sortedPopulation.emplace_back(populationBinary[i], fitnessValues[i]);
    }

    // Sorting the combined vector based on fitness in descending order
    std::sort(sortedPopulation.begin(), sortedPopulation.end(), [](const IndividualWithFitness& a, const IndividualWithFitness& b) {
        return a.fitness < b.fitness;
    });

    populationBinary.clear();
    fitnessValues.clear();

    // Putting the top individuals back into the population until the desired population size is reached
    for (size_t i = 0; i < populationSize && i < sortedPopulation.size(); ++i) {
        populationBinary.push_back(sortedPopulation[i].individual);
        fitnessValues.push_back(sortedPopulation[i].fitness);
    }
}

void GeneticAlgorithm::evaluateFitness(std::function<double(const std::vector<double>&)> objectiveFunction, int option, double liminf,double limsup) {
    switch (option){
    case 1:
        fitnessValues.clear();
        for (auto& individual : population) {
            double fitness = objectiveFunction(individual);
            fitnessValues.push_back(fitness);
        }
        break;
    case 2:
        fitnessValues.clear();
        for (auto& bina : populationBinary) {
            std::vector<double> individual = decodeAllele(bina, liminf, limsup);
            double fitness = objectiveFunction(individual);
            fitnessValues.push_back(fitness);
        }
        break;
    default:
        break;
    }
}

std::vector<int> GeneticAlgorithm::selection() {
    std::vector<int> selectedParents;
    double totalFitness = accumulate(fitnessValues.begin(), fitnessValues.end(), 0.0);

    std::vector<double> probabilities;
    for (double fitness : fitnessValues) {
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

std::vector<int> GeneticAlgorithm::tournamentSelection(int tournamentSize) {
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

// Two point crossover
std::vector<std::vector<double>> GeneticAlgorithm::crossover(std::vector<int>& selectedParents) {
    std::vector<std::vector<double>> children;
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    std::uniform_int_distribution<> pointDist(1, population[0].size() - 2); // Ensure valid points

    for (size_t i = 0; i < selectedParents.size() - 1; i += 2) {
        if (dis(gen) < crossoverRate) {
            std::vector<double>& parent1 = population[selectedParents[i]];
            std::vector<double>& parent2 = population[selectedParents[i + 1]];

            // Generate two points for crossover
            int point1 = pointDist(gen);
            int point2 = pointDist(gen);
            // Ensure point1 < point2
            if (point1 > point2) std::swap(point1, point2);

            std::vector<double> child1 = parent1;
            std::vector<double> child2 = parent2;

            // Perform the crossover
            for (int j = point1; j <= point2; ++j) {
                child1[j] = parent2[j];
                child2[j] = parent1[j];
            }

            children.push_back(child1);
            children.push_back(child2);
        }
    }
    return children;
}

std::vector<std::vector<std::string>> GeneticAlgorithm::crossover_binary(std::vector<int>& selectedParents) {
    std::vector<std::vector<std::string>> children;
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    std::uniform_int_distribution<> pointDist(1, populationBinary[0].size() - 2); // Ensure valid points

    for (size_t i = 0; i < selectedParents.size() - 1; i += 2) {
        if (dis(gen) < crossoverRate) {
            std::vector<std::string>& parent1 = populationBinary[selectedParents[i]];
            std::vector<std::string>& parent2 = populationBinary[selectedParents[i + 1]];

            int point1 = pointDist(gen);
            int point2 = pointDist(gen);
            if (point1 > point2) std::swap(point1, point2);

            std::vector<std::string> child1 = parent1;
            std::vector<std::string> child2 = parent2;

            for (int j = point1; j <= point2; ++j) {
                child1[j] = parent2[j];
                child2[j] = parent1[j];
            }

            children.push_back(child1);
            children.push_back(child2);
        }
    }
    return children;
}

std::vector<std::vector<double>> GeneticAlgorithm::mutation(std::vector<std::vector<double>> child) {
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    for (auto& individual : child) {
        for (double& gene : individual) {
            // Mutation at a given rate
            if (dis(gen) < mutationRate) {
                // Random change within the gene
                // Mean 0
                // Standard deviation 0.1
                double mutationChange = std::normal_distribution<double>(0.0, 0.1)(gen);
                gene += mutationChange;
            }
        }
    }
    return child;
}

double GeneticAlgorithm::adaptiveMutationRate(int currentGeneration) {
    // Start with a high mutation rate and decay exponentially
    double decayRate = 0.05; // Experiment with this rate
    return mutationRate * exp(-decayRate * currentGeneration);
}


std::vector<std::vector<std::string>> GeneticAlgorithm::mutation_binary(std::vector<std::vector<std::string>> children, int currentGeneration) {
    // Decreasing the mutation rate as the number of generations increases
   double dynamicMutationRate = mutationRate / (1.0 + (double)currentGeneration / numberOfGenerations);
    
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    for (auto& individual : children) {
        for (std::string& gene : individual) {
            for (char& bit : gene) {
                if (dis(gen) < dynamicMutationRate) {
                    // Perform bit flip mutation
                    bit = (bit == '0') ? '1' : '0';
                }
            }
        }
    }
    return children;
}


bool GeneticAlgorithm::shouldStop(int currentGeneration, const int option) {
    double bestFitness = *std::min_element(fitnessValues.begin(), fitnessValues.end());
    double worstFitness = *std::max_element(fitnessValues.begin(), fitnessValues.end());

    std::cout << bestFitness << std::endl;
    std::cout << worstFitness << std::endl;

    // |f(⃗xbest)−f(⃗x∗)| ≤ ε

    // Criterion for the Rosenbrock function, where the optimal point x* is a vector of 1's
    if (option == 1) {
        for (const auto& individual : population) {
            double fitness = rosenbrockFunction(individual);
            // Fitness close to 0 means we're close to the vector of 1's.
            if (std::abs(fitness) <= 0.001) {
                return true; 
            }
        }
    }

    // Criterion for the Ackley function, where the optimal point x* is 0.
    if (option == 2) {
        if (std::abs(bestFitness) <= 0.001) {
            return true;
        }
    }

    // f(⃗xworst)−f(⃗xbest)| ≤ ε
    if (std::abs(worstFitness - bestFitness) <= 0.001) {
        return true;
    }
    return false;
}

double GeneticAlgorithm::rosenbrockFunction(const std::vector<double>& individual) {
    double rosenbrock = 0.0;
    for (size_t i = 0; i < individual.size() - 1; ++i) {
        rosenbrock += 100 * std::pow((individual[i+1] - std::pow(individual[i], 2)), 2) + std::pow((1 - individual[i]), 2);
    }
    return rosenbrock;
}

double GeneticAlgorithm::ackleyFunction(const std::vector<double>& individual) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (auto x_i : individual) {
        sum1 += std::pow(x_i, 2);
        sum2 += std::cos(2 * M_PI * x_i);
    }
    double ackley = -20 * std::exp(-0.2 * std::sqrt(sum1 /  individual.size())) - std::exp(sum2 / individual.size()) + 20 + M_E;
    return ackley;
}