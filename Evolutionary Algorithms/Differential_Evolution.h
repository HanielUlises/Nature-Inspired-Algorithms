#ifndef DIFFERENTIAL_EVOLUTION_H
#define DIFFERENTIAL_EVOLUTION_H

#include <vector>
#include <string>

class DifferentialEvolution {
public:
    DifferentialEvolution(int population_size, int dimension);
    void runEvolution(const std::string& strategy);

private:
    int population_size_;
    int dimension_;
};

#endif // DIFFERENTIAL_EVOLUTION_H
