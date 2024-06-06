#include <iostream>
#include "ACO.h"
#include <vector>
#include <random>
#include <cmath>
#include <numeric>

std::random_device rd;
std::mt19937 gen(rd());

ACO::ACO(int num_ants, int max_iterations, std::vector<std::vector<int>> distance,std::vector<std::vector<int>> flows, double evaporation, int alpha, int beta, int q) 
    : num_ants(num_ants), max_iterations(max_iterations), distance(distance), flows(flows), evaporation(evaporation), alpha(alpha), beta(beta), q(q) {
    for(size_t i = 0; i < distance.size(); ++i){
        std::vector<double> phe;
        for(size_t j = 0; j < distance[0].size(); ++j){
            double a = (double) 1 / distance.size();
            phe.push_back(a);
        }
        pheromone.push_back(phe);
    }
}

int ACO::selectNext_city(int actual, std::vector<int> not_busy){
    std::vector<double> probabilities;
    for (size_t i = 0; i < not_busy.size(); i++){
        auto phe = std::pow(pheromone[actual][not_busy[i]], alpha);
        double visibility = std::pow((1.0 / flows[actual][not_busy[i]]), beta);
        double probability = phe * visibility;
        probabilities.push_back(probability);

    }
    double sum_probabilities = accumulate(probabilities.begin(), probabilities.end(), 0.0);
    for (size_t i = 0; i < not_busy.size(); i++){
        probabilities[i]=probabilities[i]/sum_probabilities;
    }
    double max_probability = *std::max_element(probabilities.begin(), probabilities.end());
    int index=0;
    for(const auto& prob:probabilities){
        if(prob==max_probability) return index;
        else index++;
    }
    
}

std::vector<std::vector<int>> ACO::constructRoute(int num_city){
    std::uniform_int_distribution<> dis(0, num_city-1);
    std::vector<std::vector<int>> routes;

    for (size_t i = 0; i < num_ants; i++){
        int init = dis(gen);
        std::vector<int> route; 
        std::vector<int> not_busy(num_city);
        std::iota(not_busy.begin(), not_busy.end(), 0);
        route.push_back(init);
        auto ref = std::find(not_busy.begin(), not_busy.end(), init); 
        not_busy.erase(ref);
        
        int actual = init;
        while (not_busy.size()!=0){
            auto siguiente = selectNext_city(actual, not_busy);
            route.push_back(not_busy[siguiente]);
            actual = not_busy[siguiente];
            auto ref = std::find(not_busy.begin(), not_busy.end(), not_busy[siguiente]); 
            not_busy.erase(ref);
        }
        routes.push_back(route);
    }
    return routes;
}

std::vector<int> ACO::objectiveFunction(std::vector<std::vector<int>> routes){
    std::vector<int> fitness;
    for (const auto& route:routes){
        double fit = 0.0f;
        for (size_t i=0; i<route.size();++i){
            for(size_t j=0; j<route.size();++j){
                fit += distance[i][j] * flows[route[i]][route[j]];
            }   
        }
        fitness.push_back(fit);
    }
    return fitness;
}

void ACO::updatePheromone(std::vector<std::vector<int>> routes, std::vector<int> fitness){
    for(const auto& phe : pheromone){
        for (auto p:phe){
            p *= (1.0f - evaporation);
        }
    }
    for (size_t i = 0; i<routes.size(); ++i){
        std::vector<int> route = routes[i];
        int k=0;
        for (size_t j = 0; j < route.size()-1; j++){
            pheromone[route[j]][route[j+1]] += (double) q / fitness[i];
            k++;
        }
        pheromone[route[route.size()-1]][route[0]] += (double) q / fitness[i];
    }
    
}

std::vector<int> ACO::optimize(){
    std::vector<int> best_route;
    std::vector<std::vector<int>> routes;
    std::vector<int> fitness;
    double best_fit = 1000000000000000;
    int num_city=distance[0].size();

    for (size_t i = 0; i < max_iterations; i++){

        routes = constructRoute (num_city);
        fitness = objectiveFunction(routes);
        updatePheromone(routes,fitness);
        int min_fit = *std::min_element(fitness.begin(), fitness.end());
        if (min_fit<best_fit){
            best_fit = min_fit;
            int index=0;
            for(const auto& fit:fitness){
                if(fit==best_fit) break;
                else index++;
            }
            best_route = routes[index];
        }
        
    }
    std::cout<<"\nBest_route: ";
    for (const auto& r:best_route){
        std::cout<<r<<" | ";
    }
    std::cout<<std::endl;
    return best_route;
    
}

void ACO::printSolution(std::vector<int> route){
    double fit = 0.0f;
    for (size_t i=0; i<route.size();++i){
        for(size_t j=0; j<route.size();++j){
            fit += distance[i][j] * flows[route[i]][route[j]];
        }   
    }
    
    std::cout<<"Evaluation with combination : "<<"| Num Ants = "<<num_ants;
    std::cout<<" | Alpha = "<<alpha;
    std::cout<<" | Beta = "<<beta<<std::endl;
    std::cout<<"evaluations result: "<<fit<<std::endl;

}