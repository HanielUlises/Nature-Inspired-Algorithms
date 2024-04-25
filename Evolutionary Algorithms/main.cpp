#include "Evolutionary_Strategy.h"
#include "Differential_Evolution.h"

int main() {
    int mu = 20;
    int lambda = 30;
    
    EvolutionaryStrategy es(mu, lambda);
    es.runEvolution();

    int population_size = 50;
    int dimension = 10;
    int generation_max=100;
    double lower_bound=-10;
    double upper_bound=10;
    double F=0.6;
    double Cr=0.9;
    
    DifferentialEvolution de(population_size, dimension,generation_max,lower_bound,upper_bound,F,Cr);
    
    de.runEvolution("rand/1/bin");
    de.runEvolution("rand/1/exp");
    de.runEvolution("best/1/bin");
    de.runEvolution("best/1/exp");
    return 0;
}
