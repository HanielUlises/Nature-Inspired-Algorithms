#include "Evolutionary_Strategy.h"
#include "Differential_Evolution.h"

int main() {
    int mu = 20;
    int lambda = 30;
    
    EvolutionaryStrategy es(mu, lambda);
    es.runEvolution();

    int population_size = 50;
    int dimension = 10;
    int generation_max=500;
    double lower_bound=-10;
    double upper_bound=10;
    double F=0.6;
    double Cr=0.9;

    const int bin=0;
    const int expo=1;
    const int rand=0;
    const int best=1;
    
    DifferentialEvolution de(population_size, dimension,generation_max,lower_bound,upper_bound,F,Cr);
    
    //de.runEvolution(rand,bin);//de.runEvolution("rand/1/bin");
    //de.runEvolution(best,bin);//de.runEvolution("best/1/bin");
    de.runEvolution(rand,expo);//de.runEvolution("rand/1/exp");
    //de.runEvolution("best/1/exp");
    return 0;
}
