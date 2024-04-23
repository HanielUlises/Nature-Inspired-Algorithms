#include "Evolutionary_Strategy.h"
#include "Differential_Evolution.h"

int main() {
    int mu = 20;
    int lambda = 30;
    
    EvolutionaryStrategy es(mu, lambda);
    es.runEvolution();

    int population_size = 50;
    int dimension = 10;
    
    DifferentialEvolution de(population_size, dimension);
    
    de.runEvolution("rand/1/bin");
    de.runEvolution("rand/1/exp");
    de.runEvolution("best/1/bin");
    de.runEvolution("best/1/exp");
    return 0;
}
