#ifndef EVOLUTIONARY_STRATEGY_H
#define EVOLUTIONARY_STRATEGY_H

#include <vector>

class EvolutionaryStrategy {
public:
    EvolutionaryStrategy(int mu, int lambda);
    void runEvolution();

private:
    int mu_;
    int lambda_;
};

#endif // EVOLUTIONARY_STRATEGY_H
