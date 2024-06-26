#include <iostream>
#include <vector>
#include <fstream>
#include "ACO.h"

int main() {
    // Matriz de distancias y flujos
    std::vector<std::vector<int>> distance = {
        {0, 21, 95, 82, 56, 41, 6, 25, 10, 4, 63, 6},
        {21, 0, 44, 40, 75, 79, 0, 89, 35, 9, 1, 85},
        {95, 44, 0, 84, 12, 0, 26, 91, 11, 35, 82, 26},
        {82, 40, 84, 0, 69, 56, 86, 45, 91, 59, 18, 76},
        {56, 75, 12, 69, 0, 39, 18, 57, 36, 61, 36, 21},
        {41, 79, 0, 56, 39, 0, 71, 11, 29, 82, 82, 6},
        {6, 0, 26, 86, 18, 71, 0, 71, 8, 77, 74, 30},
        {25, 89, 91, 45, 57, 11, 71, 0, 89, 76, 76, 40},
        {10, 35, 11, 91, 36, 29, 8, 89, 0, 93, 56, 1},
        {4, 9, 35, 59, 61, 82, 77, 76, 93, 0, 50, 4},
        {63, 1, 82, 18, 36, 82, 74, 76, 56, 50, 0, 36},
        {6, 85, 26, 76, 21, 6, 30, 40, 1, 4, 36, 0}
    };

    std::vector<std::vector<int>> flows = {
        {0, 27, 85, 2, 1, 15, 11, 35, 11, 20, 21, 61},
        {27, 0, 80, 58, 21, 76, 72, 44, 85, 94, 90, 51},
        {85, 80, 0, 3, 48, 29, 90, 66, 41, 15, 83, 96},
        {2, 58, 3, 0, 74, 45, 65, 40, 54, 83, 14, 71},
        {1, 21, 48, 74, 0, 77, 36, 53, 37, 26, 87, 76},
        {15, 76, 29, 45, 77, 0, 91, 13, 29, 11, 77, 32},
        {11, 72, 90, 65, 36, 91, 0, 87, 67, 94, 79, 2},
        {35, 44, 66, 40, 53, 13, 87, 0, 10, 99, 56, 70},
        {11, 85, 41, 54, 37, 29, 67, 10, 0, 99, 60, 4},
        {20, 94, 15, 83, 26, 11, 94, 99, 99, 0, 56, 2},
        {21, 90, 83, 14, 87, 77, 79, 56, 60, 56, 0, 60},
        {61, 51, 96, 71, 76, 32, 2, 70, 4, 2, 60, 0}
    };

    // Parámetros para la optimización
    std::vector<int> num_ants = {10, 25, 50, 100};
    std::vector<int> alpha = {50, 100, 200};
    std::vector<int> beta = {2, 10};

    std::ofstream resultsFile;
    resultsFile.open("ACO_Results.csv");
    resultsFile << "Num_Ants,Alpha,Beta,Fit\n";

    for (auto nants : num_ants) {
        for (auto a : alpha) {
            for (auto b : beta) {
                ACO ant(nants, 1000, distance, flows, 0.5, a, b, 1);
                std::vector<int> solution = ant.optimize();
                double fit = ant.evaluateSolution(solution);
                resultsFile << nants << "," << a << "," << b << "," << fit << "\n";
                std::cout << "Simulation completed for: Num Ants = " << nants << ", Alpha = " << a << ", Beta = " << b << ", Fitness = " << fit << std::endl;
            }
        }
    }

    resultsFile.close();
    return 0;
}
