#include "GeneticAlgorithm.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <iostream>

std::random_device rd;
std::mt19937 gen(rd());

GeneticAlgorithm::GeneticAlgorithm(int populationSize, int numberOfGenerations, double crossoverRate, double mutationRate)
    : populationSize(populationSize), numberOfGenerations(numberOfGenerations),
      crossoverRate(crossoverRate), mutationRate(mutationRate) {
}

GeneticAlgorithm::~GeneticAlgorithm() {
}

void GeneticAlgorithm::run() {
    std::vector<std::vector<double>> child;
    initializePopulation();

    for (int i = 0; i < numberOfGenerations; ++i) {
        
        // Functions test
        evaluateFitness(std::bind(&GeneticAlgorithm::rosenbrockFunction, this, std::placeholders::_1));    
        // evaluateFitness(std::bind(&GeneticAlgorithm::ackleyFunction, this, std::placeholders::_1));
        auto selectedParents = selection();
        //se muere
        child=crossover(selectedParents);
        child=mutation(child); //child mutados
        
        for (size_t i = 0; i < child.size(); i++){
            population.push_back(child[i]);
        }
        evaluateFitness(std::bind(&GeneticAlgorithm::rosenbrockFunction, this, std::placeholders::_1));
        elitismParents();

        // Analysis of fitness
        double bestFitness = *std::min_element(fitnessValues.begin(), fitnessValues.end());
        double worstFitness = *std::max_element(fitnessValues.begin(), fitnessValues.end());
        double averageFitness = std::accumulate(fitnessValues.begin(), fitnessValues.end(), 0.0) / fitnessValues.size();

        bestFitnessHistory.push_back(bestFitness);
        worstFitnessHistory.push_back(worstFitness);
        averageFitnessHistory.push_back(averageFitness);

        if (shouldStop(i)) break;
        std::cout<<i<<std::endl;
    }
    plotConvergenceGraph();
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

void GeneticAlgorithm::initializePopulation() {
    population.clear();

    // CHECK
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
}

void GeneticAlgorithm::evaluateFitness(std::function<double(const std::vector<double>&)> objectiveFunction) {
    fitnessValues.clear();
    for (auto& individual : population) {
        double fitness = objectiveFunction(individual);
        fitnessValues.push_back(fitness);
    }
}

std::vector<int> GeneticAlgorithm::selection() {
    std::vector<int> selectedParents;
    double totalFitness = std::accumulate(fitnessValues.begin(), fitnessValues.end(), 0.0);

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

// Two point crossover
std::vector<std::vector<double>> GeneticAlgorithm::crossover(std::vector<int>& selectedParents) {
    std::vector<std::vector<double>> children;
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    std::uniform_int_distribution<> pointDist(1, population[0].size() - 2); 

    for (size_t i = 0; i < selectedParents.size() - 1; i += 2) {
        if (dis(gen) < crossoverRate) {
            std::vector<double>& parent1 = population[selectedParents[i]];
            std::vector<double>& parent2 = population[selectedParents[i + 1]];

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

bool GeneticAlgorithm::shouldStop(int currentGeneration) {
    double bestFitness = *std::min_element(fitnessValues.begin(), fitnessValues.end());
    double worstFitness = *std::max_element(fitnessValues.begin(), fitnessValues.end());
    double optimalGlobalFitness = 0.0; // Hay que ajustar esto;

    std::cout << bestFitness << std::endl;
    std::cout << worstFitness << std::endl;
    // |f(⃗xbest)−f(⃗x∗)| ≤ ε
    if (std::abs(bestFitness - optimalGlobalFitness) <= 0.001) {
        return true;
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