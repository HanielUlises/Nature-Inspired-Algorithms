#ifndef CYCLECROSSOVER_H
#define CYCLECROSSOVER_H

#include <vector>

class CycleCrossover {
public:
    static std::pair<std::vector<int>, std::vector<int>> crossover(
        const std::vector<int>& parent1, const std::vector<int>& parent2
    );
};

#endif // CYCLECROSSOVER_H
