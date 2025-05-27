#include "CycleCrossover.h"
#include <vector>
#include <utility>
#include <algorithm>

std::pair<std::vector<int>, std::vector<int>> CycleCrossover::crossover(
    const std::vector<int>& parent1, const std::vector<int>& parent2
) {
    size_t size = parent1.size();
    std::vector<int> child1(size, -1);
    std::vector<int> child2(size, -1);

    size_t i = 0;
    bool inverse = false;

    while (true) {
        while (child1[i] == -1) {
            child1[i] = inverse ? parent2[i] : parent1[i];
            child2[i] = inverse ? parent1[i] : parent2[i];

            int value_to_find = inverse ? parent1[i] : parent2[i];
            i = std::distance(parent1.begin(), std::find(parent1.begin(), parent1.end(), value_to_find));
        }

        if (std::all_of(child1.begin(), child1.end(), [](int val) { return val != -1; })) {
            break;
        }

        i = std::distance(child1.begin(), std::find(child1.begin(), child1.end(), -1));
        inverse = !inverse;
    }

    return { child1, child2 };
}
