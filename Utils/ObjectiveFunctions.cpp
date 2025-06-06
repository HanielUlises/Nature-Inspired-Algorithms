#include "ObjectiveFunctions.h"
#include <iostream>
#include <algorithm>
#include <cmath>

double rosenbrockFunction(const std::vector<double>& individual) {
    double rosenbrock = 0.0;
    for (size_t i = 0; i < individual.size() - 1; ++i) {
        rosenbrock += 100 * std::pow((individual[i+1] - std::pow(individual[i], 2)), 2) + std::pow((1 - individual[i]), 2);
    }
    return rosenbrock;
}

double ackleyFunction(const std::vector<double>& individual) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (auto x_i : individual) {
        sum1 += std::pow(x_i, 2);
        sum2 += std::cos(2 * M_PI * x_i);
    }
    double ackley = -20 * std::exp(-0.2 * std::sqrt(sum1 /  individual.size())) - std::exp(sum2 / individual.size()) + 20 + M_E;
    return ackley;
}
double GriewankFunction(const std::vector<double>& individual){
    double sum = 0.0f;
    double mult = 0.0f;
    for (size_t i = 0; i < individual.size(); i++){
        auto x_i=individual[i];
        sum += std::pow(x_i, 2)/4000.0f;
        mult *= std::cos(x_i / sqrt(i+1));
    }
    
    double griewank = (sum)-(mult)+1;

    return griewank;

}
double RastriginFunction(const std::vector<double>& individual){
    double sum = 0.0;
    for (auto x_i : individual) {
        sum += (std::pow(x_i, 2))-(10*(std::cos(2*M_PI*x_i)));
    }
    double rastrigin = sum+(10*individual.size());
    return rastrigin;
}


double Langermann(const std::vector<double> &individual) {
    std::vector<double> a = {3.0, 5.0, 2.0, 1.0, 7.0};
    std::vector<double> b = {5.0, 2.0, 1.0, 4.0, 9.0};
    std::vector<double> c = {1.0, 2.0, 5.0, 2.0, 3.0};

    int size_vec = a.size();
    double sum = 0.0;

    for (int i = 0; i < size_vec; i++) {
        double term = std::pow(individual[0] - a[i], 2) + std::pow(individual[1] - b[i], 2);
        sum += c[i] * std::cos(M_PI * term) / std::exp(term / M_PI);
    }
    return -sum; 
}

double dropWave(const std::vector<double> &individual){
    return - 1 + std::cos(12*std::sqrt(std::pow(individual[0],2) + std::pow(individual[1],2))/
            (0.5)* (std::pow(individual[0],2) + std::pow(individual[1],2)) + 2);
}
