#include "function_test_1.h"

int main (){
    // Parte uno
    std::vector<std::pair<double, double>> points_p1 = {
        {-5.5, 2.3}, {1.2, -1.1}, {6.0, 3.0}, {1.0, 0.9}, {-0.9, 0.8}, {2, 4.0}
    };

    std::vector<double> z, expected_z;

    std::cout << "Points 1 tested" << std::endl;
    for (const auto& point : points_p1) {
        std::cout << "(" << point.first << ", " << point.second << ")" << std::endl;
        f1(point.first, point.second, z);
    }

    // Mean and expectancy values
    double mean_z = mean(z);
    for (double val : z) {
        expected_val(val, mean_z, expected_z);
    }

    // Display expected values
    std::cout << std::endl << "Expected values" << std::endl;
    for (double i : expected_z) {
        std::cout << i << std::endl;
    }

    // Select a point using roulette selection
    auto& selectedPoint = rouletteSelection(points_p1, expected_z);
    std::cout << std::endl << "Selected Point: (" << selectedPoint.first << ", " << selectedPoint.second << ")" << std::endl;
}