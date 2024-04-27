#include "Evolutionary_Strategy.h"
#include "Differential_Evolution.h"
#include <numeric>
#include <cmath>
#include <iostream>

int main() {
    int mu = 20;
    int lambda = 30;
    
    EvolutionaryStrategy es(mu, lambda);
    es.runEvolution();

    int population_size = 50;
    int dimension = 10;
    int generation_max=1000;
    double lower_bound=-5.12;
    double upper_bound=5.12;
    double F=0.6;
    double Cr=0.9;

    const int bin=0;
    const int expo=1;
    const int rand=0;
    const int best=1;

    std::vector<double> bestHistory;
    for (size_t i = 0; i < 20; i++){
        DifferentialEvolution de(population_size, dimension,generation_max,lower_bound,upper_bound,F,Cr);
        double bestR;
        //bestR=de.runEvolution(rand,bin);//de.runEvolution("rand/1/bin");
        //bestR=de.runEvolution(rand,expo);//de.runEvolution("rand/1/exp");
        //bestR=de.runEvolution(best,bin);//de.runEvolution("best/1/bin");
        bestR=de.runEvolution(best,expo);
        bestHistory.push_back(bestR);
    }
    double average = accumulate(bestHistory.begin(), bestHistory.end(), 0.0) / bestHistory.size();
    double variation=0.0f;
    for (double bestx:bestHistory){
        variation=variation+(std::pow((bestx-average),2)/(bestHistory.size()-1));
    }

    double desv=sqrt(variation);
    std::cout<<"Promedio: "<<average<<std::endl;
    std::cout<<"Desviacion estandar"<<desv<<std::endl;
    
    return 0;
}
