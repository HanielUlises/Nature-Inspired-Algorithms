#include "Solution.h"
#include "GeneticAlgorithm.h"

int main (){
    int population = 100;
    int bits = 8;
    int low = 0, high = 4;
    GeneticAlgorithm test(population);
    Solution test_sol = test.perform(bits, low, high);

    return 1;
}