#include "Solution.h"
#include "GeneticAlgorithm.h"

int main (){
    int population = 100, generations = 10;
    int bits = 8;
    int low = 0, high = 4;
    int t_gp_size = 3;

    double crossoverProbability;

    GeneticAlgorithm test(population, generations, t_gp_size, crossoverProbability);
    Solution test_sol = test.perform(bits, low, high);

    return 1;
}