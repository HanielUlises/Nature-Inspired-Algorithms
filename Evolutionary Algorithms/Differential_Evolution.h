#ifndef DIFFERENTIAL_EVOLUTION_H
#define DIFFERENTIAL_EVOLUTION_H

#include <vector>
#include <string>

struct single {
    std::vector<double> values;
    double fitness;
};

class DifferentialEvolution {
public:
    DifferentialEvolution(int population_size, int dimension, int generation_max, double lower_bound, double upper_bound, 
    double F, double Cr);
    void runEvolution(const std::string& strategy);
    void initializePopulation();
    std::vector<int> selection_r();

private:    
    int population_size_;
    int dimension_;
    int generation_max_;
    double lower_bound_; 
    double upper_bound_; 
    double F_; 
    double Cr_;
    std::vector<single> population;
    single bestx;
};


#endif // DIFFERENTIAL_EVOLUTION_H
