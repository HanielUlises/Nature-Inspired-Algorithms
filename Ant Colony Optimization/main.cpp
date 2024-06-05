#include <cmath>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

#include "ACO.h"

int main() {
    // Individuals
    std::vector<std::vector<int>> distance = {
    {0, 12, 3, 23, 1},
    {12,  0,  9, 18, 3},
    {3, 9, 0, 89, 56},
    {23, 18, 89, 0, 87},
    {1, 3, 56, 87, 0}
    };
    ACO ant(10, 10, distance, 0.5, 1, 5, 1);
    ant.optimize();


    return 0;
}
